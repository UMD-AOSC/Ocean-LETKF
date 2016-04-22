#!/bin/ksh --login
set -e
#module load mpt
#module load intel
#module load netcdf
#module load nco
#module load wgrib

#INPUTS:
# EXP_DATA
# executable
# YYYYMMDDHH
# gsupdir 
# days
# months
# rtype
# REGRID

# START:
echo "GODAS solo run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/godas_solo
mkdir -p ${workdir}
cd ${workdir}

# setup directory structure
gINPUT=$workdir/INPUT
gRESTART=$workdir/RESTART
mkdir -p $gINPUT
mkdir -p $gRESTART

OBSDIR=$gsupdir/PTMP_SAL
VARDIR=$gsupdir/VARr
SBCDIR=$gsupdir/PREP_SBC

# STEVE: allow godas to use difference surface forcing fields:
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
if [ "$MEM2" -eq "00" ]; then
  # STEVE: for the 00 member, assume this is a 'control' run and use the mean forcing fields only:
  USE_MFLX=1
  USE_EFLX=0
  FLUXDIR=$FLXDIR
else
  USE_MFLX=0
  USE_EFLX=1
  FLUXDIR=$FLXDIR2
fi
echo "USE_MFLX = $USE_MFLX" >> godas_solo.out
echo "USE_EFLX = $USE_EFLX" >> godas_solo.out

#STEVE: active code:
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the date for the next analysis cycle
date=/bin/date
inc=$days #5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

# Update the date for the previous analysis cycle
pinc=5
pinc_units='days ago'
PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
PN=$IN
PS=$IS

#########################################################
### Processs inputs

#STEVE:
# Access model forecast (this is applied to the mean field in the hybrid)
MODEL=${EXP_DATA}/${YYYYMMDDHH}/model

# RESTART directory constructed from MODEL
RESTART=${EXP_DATA}/${YYYYMMDDHH}/model/$MEM3/RESTART
#STEVE: general input files for all timesteps are in this directory:
INPUT=$INPUT_INIT

#########################################################
# GODAS solo code:
# Minimal runscript: gdsSolo

#npes=$PBS_NP        # number of processors
                     # Note: If you change npes you may need to change
                     # the layout in the corresponding namelist
#===========================================================================

# Check the time_stamp
  echo "${IY}-${IM}-${ID} is the current date"
  echo "$IY	$IM	$ID	$IH	$IN	$IS" > time_stamp.out

#
#Check the existance of essential input files
#
   if [ ! -f $REGRID/grid_spec.nc ]; then
     echo "ERROR: required input file does not exist $REGRID/grid_spec.nc "
     exit 1
   fi
   if [ ! -f $RESTART/ocean_temp_salt.res.nc ]; then
     echo "ERROR: required input file does not exist $RESTART/ocean_temp_salt.res.nc "
     exit 1
   fi
   if [ ! -f $RESTART/$rtype.res ]; then
     echo "ERROR: required input file does not exist $RESTART/$rtype.res "
     exit 1
   fi

#-------------------------------------------

  ############################
  # Prepare observation data
  ############################
  # arcdir (directory with profile observations)
  # here (current directory)
  # expdir ()
  # d1, d2 (range of dates to look forward and backward for observations)
  # gridDIR (directory with grid_spec.nc file)
  # scriptDir (directory of support C script)

  #STEVE: I copied the obs that were converted to LETKF format to this directory (arcdir):
  #STEVE: should the timestamp be for the beginning or end (at the analysis time) of the cycle?
  cat > $OBSDIR/cDG.header.csh <<EOF
  set arcdir=$GOBSDIR
  set d1=0
  set d2=4
  set gridDir=$REGRID
  set scriptDir=$OBSDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $OBSDIR
  echo "In directory: "
  pwd
  rm -f tmpa.mom
  rm -f sala.mom

  $OBSDIR/cDG.csh

  #STEVE: observation files:
  mv $OBSDIR/tmpa.mom $gINPUT/
  mv $OBSDIR/sala.mom $gINPUT/

  ############################
  # Prepare background error covariance
  ############################
  # RESTART (RESTART directory from last ocean model run)
  # gridDir (directory of grid_spec.nc file)
  # scriptDir (directory of support C script)

  cat > $VARDIR/mEN4.header.csh <<EOF
  set RESTART=$RESTART
  set gridDir=$REGRID
  set scriptDir=$VARDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $VARDIR
  echo "In directory: "
  pwd
  rm -f tvv.mom
  rm -f svv.mom

  $VARDIR/mEN4.csh

  #STEVE: background error covariance files:
  mv $VARDIR/tvv.mom $gINPUT/
  mv $VARDIR/svv.mom $gINPUT/

  ############################
  # Prepare surface boundary conditions (SST, SSS)
  ############################
  # SBC_DIR (surface forcing from the atmos model)
  # SST_DIR (sst reanalysis data)
  # SSS_DIR (sss climatology data, salt12.nc)
  # L_YR (last year of data)
  # gridDir
  # scriptDir

  cp /sw/xe6/wgrib/1.8.1.0b/sles11.1_gnu4.3.4/bin/wgrib ${workdir}
  cat > $SBCDIR/mPS.header.csh <<EOF
  set SBC_DIR=$FLUXDIR
  set SST_DIR=$SSTDIR
  set SSS_DIR=$EXP_DATA/INPUT
  @ L_YR = 2013
  set gridDir=$REGRID
  set scriptDir=$SBCDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $SBCDIR
  echo "In directory: "
  pwd

  $SBCDIR/mPS.csh

  mv $SBCDIR/*_sfc_restore.nc $gINPUT/

  ############################
  # Prep the input files for godas_solo:
  ############################
  cd $gINPUT
  ln -f $RESTART/ocean_temp_salt.res.nc $gINPUT/
  ln -f $RESTART/$rtype.res $gINPUT/
  ln -f $INPUT_INIT/grid_spec.nc $gINPUT/

  ############################
  cd $workdir

  datatable=$gsupdir/data_table            # path to the data override table.
  if [ ! -f $namelist ]; then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  fi
  ln -fs $datatable $workdir

  namelist=$gsupdir/namelist   # path to namelist file w/o coupler_nml
  if [ ! -f $datatable ]; then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  fi
  ln -fs $namelist $workdir/

  ############################
  # set input.nml
  ############################

 cat > input.nml <<EOF
 &godas_solo_nml
        months = $months,
        days   = $days,
        hours=0
        minutes=0
        seconds=0
        calendar = 'julian',
/
EOF
  cat $namelist >> input.nml
  ln -fs $datatable  $gINPUT/data_table

# START RUN GODAS SOLO ################################################

# mpiexec_mpt -np $PBS_NP $executable
  aprun -n $PBS_NP $executable

# END   RUN GODS SOLO ################################################

#STEVE:
## Add g4p1 ocean correction to the input mean field
#cp $gRESTART/ocean_cor.res.nc $gRESTART/ocean_cor2.res.nc
#echo "Calling:"
#echo "ncrename -O -v tcor,temp -v scor,salt $gRESTART/ocean_cor2.res.nc"
#ncrename -O -v tcor,temp -v scor,salt $gRESTART/ocean_cor2.res.nc

#echo "Calling:"
#echo "ncflint -w 1,1 $RESTART/ocean_temp_salt.res.nc $gRESTART/ocean_cor.res.nc $gRESTART/ocean_temp_salt.res.nc"
#ncflint -w 1,1 $RESTART/ocean_temp_salt.res.nc $gRESTART/ocean_cor2.res.nc $gRESTART/ocean_temp_salt.res.nc
##STEVE: !!!!! Changed godas to output this file: $gRESTART/ocean_temp_salt.res.nc

#Write out a file indicating the processing is complete.
echo "This is the GODAS SOLO run assimilation for cycle ${YYYYMMDDHH}" > godas_solo.out

exit 0

