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
AOBDIR=$gsupdir/PREP_ALT
VARDIR=$gsupdir/VARr
SBCDIR=$gsupdir/PREP_SBC

# STEVE: allow godas to use different surface forcing fields:
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
RESTART=${EXP_DATA}/${YYYYMMDDHH}/model/$MEM2/RESTART
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
#  cat > $OBSDIR/cDG.header.csh <<EOF
#  set arcdir_t=$GOBSDIR_t
#  set arcdir_s=$GOBSDIR_s
#  set d1=$d1
#  set d2=$d2
#  set gridDir=$REGRID
#  set scriptDir=$OBSDIR
#  set yr=$IY
#  set mo=$IM
#  set dy=$ID
#EOF

# cd $OBSDIR
  cd $gINPUT #STEVE: this is on the lustre filesystem, $OBSDIR is not
  echo "In directory: "
  pwd
  rm -f tmpa.mom
  rm -f sala.mom

  # cDG.csh (a) extracts obs profiles from yearly tar files for TMP and SAL (via c program: extDysP4nc.c)
  #         (b) reads daily profile files, combines them and writes them in NEW M4 assimilation format (via c pograms: cmbDLstTh4nc.c and cmbDLstPs4nc.c)
# $OBSDIR/cDG.csh
  # Do obs setup and Create the observation data files here instead:
  # MARK
  inc_units=days
  D1Y=`$date -d "$IY-$IM-$ID $d1 $inc_units" +%Y`
  D1M=`$date -d "$IY-$IM-$ID $d1 $inc_units" +%m`
  D1D=`$date -d "$IY-$IM-$ID $d1 $inc_units" +%d`
  D1H=`$date -d "$IY-$IM-$ID $d1 $inc_units" +%H`
  D1N=$IN
  D1S=$IS

  #STEVE: Handle the case at the beginning of the experiment period:
  if [ $D1Y -lt 1991 ]; then
    D1Y=$IY
    D1M=$IM
    D1D=$ID
    D1H=$IH
    D1N=$IN
    D1S=$IS
  fi

  D2Y=`$date -d "$IY-$IM-$ID $d2 $inc_units" +%Y`
  D2M=`$date -d "$IY-$IM-$ID $d2 $inc_units" +%m`
  D2D=`$date -d "$IY-$IM-$ID $d2 $inc_units" +%d`
  D2H=`$date -d "$IY-$IM-$ID $d2 $inc_units" +%H`
  D2N=$IN
  D2S=$IS

  sinc=1
  sinc_units=day
  YYYY=$D1Y
  MM=$D1M
  DD=$D1D
  while test $YYYY$MM$DD -lt $D2Y$D2M$D2D
  do
    NY=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%Y`
    NM=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%m`
    ND=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%d`
    # Start 1 day later, end $days days later

    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_t ..."
    ln -f $GOBSDIR_t/${YYYY}${MM}${DD}tmp.nc .
    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_s ..."
    ln -f $GOBSDIR_s/${YYYY}${MM}${DD}sal.nc .

    YYYY=$NY
    MM=$NM
    DD=$ND
  done

  echo "Creating list of *.nc obs files for tmp and sal..."
  ls ????????tmp.nc > tmpFiles
  ls ????????sal.nc > salFiles

  # STEVE: NOTE! These are specifically for the synthesized observations.
  #              The original versions are needed when going back to real observations
  #              to convert in situ temps to potential temp (same is needed for LETKF) 
  #              and to compute the estimated standard deviation
  ln -f $INPUT_INIT/grid_spec.nc .
  cp -f $OBSDIR/cmbDLstTh4nc_synth .
  cp -f $OBSDIR/cmbDLstPs4nc_synth . 
  nlv=40
  echo "./cmbDLstTh4nc_synth -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k $nlv"
  ./cmbDLstTh4nc_synth -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k $nlv 
  echo "./cmbDLstPs4nc_synth -f salFiles -o sala.mom -g grid_spec.nc -k $nlv"
  ./cmbDLstPs4nc_synth -f salFiles -o sala.mom -g grid_spec.nc -k $nlv

  if [ ! -f tmpa.mom ]; then
    echo "tmpa.mom does not exist! Exiting godas_solo.ksh..."
    exit 1
  fi
  if [ ! -f sala.mom ]; then
    echo "sala.mom does not exist! Exiting godas_solo.ksh..."
    exit 2
  fi

  #STEVE: observation files:
# mv $OBSDIR/tmpa.mom $gINPUT/
# mv $OBSDIR/sala.mom $gINPUT/
  rm ????????tmp.nc
  rm ????????sal.nc

  #STEVE: do the same to set up altimetry if being used:
  if [ $USE_ALTIMETRY -eq "1" ]; then

    echo "Create and combine weekly altimetry files:"
    echo " "
    echo "Begin"

    #(Skip the use of cWG.csh and just implement it here:)
    latmin=-40
    latmax=50
    nlv=40
    dte=$YYYY$MM$DD

    # Remove existing files, if any:
    rm -f runSWA4
    rm -f txj1j2_cTbl
    rm -f spltWkAlt4
    rm -f cmbWksAlt4
    rm -f alta.mom

#   ln -f $INPUT_INIT/grid_spec.nc .
    cp -f $AOBDIR/runSWA4 .
    cp -f $AOBDIR/txj1j2_cTbl .
    cp -f $AOBDIR/spltWkAlt4 .
    nweeks=5 # number of weeks to combine
    echo "./runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv"
#   ./runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv
    aprun -n 1 $gINPUT/runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv

#   if [ ! -f runSWA4.out ]; then
    if [ ! -f split.ksh ]; then
      exit 5
    fi

    chmod 755 split.ksh
    ./split.ksh
    
    #STEVE: had to edit runSWA4.c to output comands to a file: runSWA4.out
#   while read line
#   do
#     # This should evaluate spltWkAlt4 for each week file:
#     echo "Evaluating line: $line"
#     eval "$line"
#   done < runSWA4.out

    ls -1 swssh.?????? > cWGLst

    cp -f $AOBDIR/cmbWksAlt4 .
    echo "./cmbWksAlt4 -f cWGLst -o alta.mom"
    ./cmbWksAlt4 -f cWGLst -o alta.mom

    if [ ! -f alta.mom ]; then
      echo "alta.mom does not exist! Exiting godas_solo.ksh..."
      exit 3
    fi
    
    # Copy godas support files for altimetry:
    # (1) The climatology
    cp $gsupdir/aEtaCds9399.i3e $gINPUT/
    # (2) The coefficients for the linear transformation
    cp $gsupdir/cdnzTS.m4 $gINPUT/

#   mv alta.mom $gINPUT/
    rm cWGLst
    rm swssh.??????
  fi

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

  echo "Preparing background error covariance..."
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

  echo "Preparing surface boundary conditions..."

  cp /sw/xe6/wgrib/1.8.1.0b/sles11.1_gnu4.3.4/bin/wgrib ${workdir}
  cat > $SBCDIR/mPS.header.csh <<EOF
  set SBC_DIR=$FLUXDIR
  set SST_DIR=$SSTDIR
  set SSS_DIR=$INPUT_INIT/../
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

  # Prepare surface restore files:
  # salt12.nc and grid_spec.nc must be in $SBCDIR: 
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

  namelist=$gsupdir/$namelist_in   # path to namelist file w/o coupler_nml
  if [ ! -f $namelist ]; then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  fi
  cp -f $namelist $workdir/

  datatable=$INPUT_INIT/data_table            # path to the data override table.
  if [ ! -f $datatable ]; then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  fi
  cp -f $datatable $workdir

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
  cp -f $datatable  $gINPUT/data_table

# START RUN GODAS SOLO ################################################
  echo "Running godas_solo..."
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

