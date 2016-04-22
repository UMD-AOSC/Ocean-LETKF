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
echo "GODAS 3DVar run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/g4p1
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
echo "USE_MFLX = $USE_MFLX" >> g4p1.out
echo "USE_EFLX = $USE_EFLX" >> g4p1.out

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
pinc=$inc
pinc_units='days ago'
PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
PN=$IN
PS=$IS
workdir0=${EXP_DATA}/$PY$PM$PD$PH/g4p1

#########################################################
### Processs inputs

#STEVE:
# Compute ensemble mean of model inputs and set up in
HYBRID=${EXP_DATA}/${YYYYMMDDHH}/hybrid
mkdir -p $HYBRID

# RESTART directory constructed from LETKF analysis mean
RESTART=${EXP_DATA}/${YYYYMMDDHH}/hybrid/letkfmean
mkdir -p $RESTART
#STEVE: general input files for all timesteps are in this directory:
INPUT=$INPUT_INIT

# use nco to average RESTART files and put in RESTART/
cp /sw/xe6/nco/4.0.8/sles11.1_netcdf4.1.3_gnu4.6.2/bin/ncea .
rm -f $RESTART/ocean_temp_salt.res.nc
ncea ${EXP_DATA}/${YYYYMMDDHH}/letkf/anal???.ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc

#########################################################
# GODAS solo code:
# Minimal runscript: gdsSolo

#npes=$PBS_NP         # number of processors
                     # Note: If you change npes you may need to change
                     # the layout in the corresponding namelist
#===========================================================================

# Check the time_stamp
  echo "${IY}-${IM}-${ID} is the current date"
  echo "$IY	$IM	$ID	$IH	$IN	$IS" > time_stamp.out

#
#Check the existance of essential input files
#
#  if [ ! -f $REGRID/grid_spec.nc ]; then
   if [ ! -f $INPUT_INIT/grid_spec.nc ]; then
     echo "ERROR: required input file does not exist $INPUT_INIT/grid_spec.nc "
     exit 1
   fi
   if [ ! -f $RESTART/ocean_temp_salt.res.nc ]; then
     echo "ERROR: required input file does not exist $RESTART/ocean_temp_salt.res.nc "
     exit 1
   fi

   #STEVE: I'm not sure if this file must be in this location for the sub-shell scripts
   #MEM2 was 01
   cp -r ${EXP_DATA}/${YYYYMMDDHH}/model/$MEM2/RESTART/$rtype.res $RESTART/
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

  cat > $OBSDIR/cDG.header.csh <<EOF
  set arcdir_t=$GOBSDIR_t
  set arcdir_s=$GOBSDIR_s
  set d1=1
  set d2=5
  set gridDir=$REGRID
  set scriptDir=$OBSDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

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
  sinc=1
  sinc_units=day
  YYYY=$IY
  MM=$IM
  DD=$ID
  while test $YYYY$MM$DD -lt $TY$TM$TD
  do
    NY=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%Y`
    NM=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%m`
    ND=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%d`
    # Start 1 day later, end $days days later
    YYYY=$NY
    MM=$NM
    DD=$ND

    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_t ..."
    ln -f $GOBSDIR_t/${YYYY}${MM}${DD}tmp.nc .
    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_s ..."
    ln -f $GOBSDIR_sm/${YYYY}${MM}${DD}sal.nc .
  done

  echo "Creating list of *.nc obs files for tmp and sal..."
  ls ????????tmp.nc > tmpFiles
  ls ????????sal.nc > salFiles

  # STEVE: NOTE! These are specifically for the synthesized observations.
  #              The original versions are needed when going back to real observations
  #              to convert in situ temps to potential temp (same is needed for LETKF) 
  #              and to compute the estimated standard deviation
  # ADDITIONAL NOTE: for real observations, use script: g4p1_robs.ksh instead!
  rm -f grid_spec.nc
  ln -f $INPUT_INIT/grid_spec.nc .
  cp -f $OBSDIR/cmbDLstTh4nc_synth .
  cp -f $OBSDIR/cmbDLstPs4nc_synth . 

  echo "Creating list of *.nc obs files for tmp and sal..."
  ls ????????tmp.nc > tmpFiles
  ls ????????sal.nc > salFiles

  echo "./cmbDLstTh4nc -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 40"
# ./cmbDLstTh4nc_synth -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 30
  ./cmbDLstTh4nc_synth -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 40
  if [ ! -f tmpa.mom ]; then
    echo "tmpa.mom does not exist! Exiting godas_solo.ksh..."
    exit 1
  fi

  echo "./cmbDLstPs4nc -f salFiles -o sala.mom -g grid_spec.nc -k 40"
# ./cmbDLstPs4nc_synth -f salFiles -o sala.mom -g grid_spec.nc -k 30
  ./cmbDLstPs4nc_synth -f salFiles -o sala.mom -g grid_spec.nc -k 40
  if [ ! -f sala.mom ]; then
    echo "sala.mom does not exist! Exiting godas_solo.ksh..."
    exit 2
  fi 
 
  #STEVE: observation files:
# mv $OBSDIR/tmpa.mom $gINPUT/
# mv $OBSDIR/sala.mom $gINPUT/

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
  # Prep the input files for g4p1:
  ############################
  cd $gINPUT
  ln -f $RESTART/ocean_temp_salt.res.nc $gINPUT/
  ln -f ${EXP_DATA}/${YYYYMMDDHH}/model/$MEM2/RESTART/$rtype.res $gINPUT/
  ln -f $INPUT_INIT/grid_spec.nc $gINPUT/

  ############################
  cd $workdir

  datatable=$gsupdir/data_table            # path to the data override table.
  if [ ! -f $namelist ]; then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  fi
  cp -f $datatable $workdir

  namelist=$gsupdir/namelist   # path to namelist file w/o coupler_nml
  if [ ! -f $datatable ]; then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  fi
  cp -f $namelist $workdir/

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
  rm -f wgrib

# START RUN G4P1 ################################################

  aprun -n $PBS_NP $executable

# END   RUN G4P1 ################################################

#########################################
# Do the Hybrid-LETKF processing
#########################################

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

#STEVE:
# Next, combine the GODAS result (gRESTART/ocean_cor.res.nc) with the 
# LETKF mean analysis and recenter the ensemble around the new analysis
#
cp /sw/xe6/nco/4.0.8/sles11.1_netcdf4.1.3_gnu4.6.2/bin/ncflint .
#alpha=0.5 #1.0 #0.5
#NOTE: alpha is now an input value
one_minus_alpha=`echo "scale=2; 1 - $alpha" | bc -l`
if [[ $alpha -gt 0 && $alpha -lt 1 ]]; then
  echo "alpha   = $alpha"
  echo "1-alpha = $one_minus_alpha"
  echo "Calling:"
  echo "ncflint -w $alpha,$one_minus_alpha $gRESTART/ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc"
  ncflint -w $alpha,$one_minus_alpha $gRESTART/ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc
elif [[ $alpha -eq 1 ]]; then
  echo "alpha   = $alpha"
  echo "1-alpha = $one_minus_alpha"
  echo "Calling:"
  echo "ln -f $gRESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc"
  ln -f $gRESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc
else
  echo "alpha is not in range 0 < alpha <= 1"
  echo "alpha   = $alpha"
  echo "Exiting..."
  exit 1
fi

##################################################################################
#STEVE:
# Finally, recenter all of the letkf ocean_temp_salt.res.nc analysis fields with the new mean field

# (1) Create mean correction field {(new mean) - (old mean)}
LETKF=${EXP_DATA}/${YYYYMMDDHH}/letkf
echo "Calling:"
echo "ncflint -w 1,-1 $HYBRID/ocean_temp_salt.res.nc $HYBRID/letkfmean/ocean_temp_salt.res.nc $HYBRID/hybrid_cor.res.nc"
rm -f $HYBRID/hybrid_cor.res.nc
ncflint -w 1,-1 $HYBRID/ocean_temp_salt.res.nc $HYBRID/letkfmean/ocean_temp_salt.res.nc $HYBRID/hybrid_cor.res.nc

# (2) Add to each member
MEM=1
while test $MEM -le $MEMBERS
do
  MEM3=`printf %.3d ${MEM}` 
  afile=anal${MEM3}.ocean_temp_salt.res.nc
  # Add the correction to this member
  echo "Calling:" 
  echo "ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile"
  rm -f $HYBRID/$afile
# ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile &
  ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile
  MEM=`expr $MEM + 1`
done
rm -f ncflint

#STEVE: wait for all of the background jobs to finish
#wait

# Generate netcdf background ensemble mean for viewing and testing output
cp /sw/xe6/nco/4.0.8/sles11.1_netcdf4.1.3_gnu4.6.2/bin/ncea .
rm -f $HYBRID/gues05_temp_salt.res.nc
echo "Calling:"
echo "ncea $LETKF/gs05???.ocean_temp_salt.res.nc $HYBRID/gues05_temp_salt.res.n"
ncea $LETKF/gs05???.ocean_temp_salt.res.nc $HYBRID/gues05_temp_salt.res.nc
rm -f ncea

#Write out a file indicating the processing is complete.
echo "This is the G4P1 run and Hybrid-LETKF assimilation for cycle ${YYYYMMDDHH}" > g4p1.out

exit 0

