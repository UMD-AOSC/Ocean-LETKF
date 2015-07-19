#!/bin/ksh --login
#===============================================================================
# SCRIPT:
# model_prep.ksh
#
# PURPOSE:
# This script prepares the model input files prior to running the model.
# A separate instance can be called for each ensemble member.
#
# MODULES USED:
#  (e.g. on Gaea)
#  module swap PrgEnv-pgi PrgEnv-intel
#  module load netcdf
#
#  (If accessible, use nco, otherwise, copy nco executables from another source)
#  module load nco
#
# INPUTS:
#  YYYYMMDDHH    :: string containing 4-digit year, 2-digit month, 2-digit day, 2-digit hour
#  MEMBERID      :: Ensemble member number
#  EXP_DATA      :: directory containing experiment output data
#  days          :: forecast length (in integer days)
#  INPUT_INIT    :: Directory containing static model input files (e.g. climatologies, grid dimensions, etc.)
#  USE_ALTIMETRY :: Flag for assimilating altimetry (1==true,0==false)
#  
#  
#
#===============================================================================

set -e
wgrib=/sw/xe6/wgrib/1.8.1.0b/sles11.1_gnu4.3.4/bin/wgrib

echo "================================================================================"
echo "Model preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_prep/${MEMBERID}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
mkdir -p ${workdir2}
cd ${workdir}

echo "This is the model preparation step for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > model_prep.out

# STEVE: for the 00 member, assume this is a 'control' run and use the mean forcing fields only:
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
if [ "$MEM2" -eq "00" ]; then
  USE_MFLX=1
  USE_EFLX=0
else
  USE_MFLX=0
  USE_EFLX=1
fi
echo "USE_MFLX = $USE_MFLX" >> model_prep.out
echo "USE_EFLX = $USE_EFLX" >> model_prep.out

#DO_SFCFLUXES=1  #1 #Input via xml script
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

mkdir -p ${workdir2}/INPUT
mkdir -p ${workdir2}/RESTART

# If this is the first timestep, then we'll have to start up the model from a fixed INPUT directory
isfirst=0
#if [ "$IY" -eq "1985" ]; then
#  isfirst=1
#fi
echo "IY=$IY." >> model_prep.out
echo "isfirst=$isfirst." >> model_prep.out

if [ "$isfirst" -eq "1" ]; then
  # Copy INPUT files, identical for all members
  ln -f $INPUT_INIT/* ${workdir2}/INPUT/
else
  # Update the date for the previous analysis cycle
  date=/bin/date
  pinc=$days
  pinc_units='days ago'
  PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
  PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
  PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
  PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
  PN=$IN
  PS=$IS

  #STEVE: copy INPUT/RESTART files to here
  #STEVE: getting it started with identical initial conditions:
  ln -f $INPUT_INIT/* ${workdir2}/INPUT/
# ln -fs $EXP_DATA/INPUT/chl.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/grid_spec.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/namelist ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/ncar_precip_clim.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/RUNOFF.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/salt12.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/SNOW.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/sst_ice_clim.nc ${workdir2}/INPUT  
# ln -fs $EXP_DATA/INPUT/input.nml ${workdir2}/INPUT
  rm -f ${workdir2}/INPUT/input.nml
  cp -f $SBCDIR/input.nml ${workdir2}/INPUT/

  #STEVE: copy last timestep's RESTART files to here
  workdir0=${EXP_DATA}/$PY$PM$PD$PH/model/${MEMBERID}
  for file in `ls -d ${workdir0}/RESTART/ocean_*res*`; 
  do
    echo "linking $file to ${workdir2}/INPUT ..."
    ln -f $file ${workdir2}/INPUT
  done

  #STEVE: Just for safety, remove ocean_temp_salt.res.nc, ocean_velocity.res.nc and ocean_sbc.res.nc
  #       so they don't accidentally get used as the restart fiels
  rm -f ${workdir2}/INPUT/ocean_temp_salt.res.nc
  rm -f ${workdir2}/INPUT/ocean_velocity.res.nc
  rm -f ${workdir2}/INPUT/ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    rm -f ${workdir2}/INPUT/ocean_barotropic.res.nc
  fi

  #STEVE: copy restart files from last timestep (need to set this up manually for first timestep)
  #STEVE: if an analysis exists, replace the primary restart files with the analysis
  ##################################################################################
  # TEMPERATURE AND SALINITY Restart Files:
  ##################################################################################
  echo "datype=$datype, setting working directory for LETKF..."
  echo "Using Standard-LETKF Analysis from previous step as initial conditions..."
  workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/letkf

  echo "datype=$datype, running for LETKF..."
  if [ -f $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc ]; then
    ln -f $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc ${workdir2}/INPUT/ocean_temp_salt.res.nc
  else
    echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc"
    exit 1
  fi

  echo "Using Standard-LETKF Analysis from previous step as initial conditions for ocean_velocity and ocean_sbc ..."
  workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/letkf

  ##################################################################################
  # VELOCITY Restart Files:
  ##################################################################################
  if [ -f $workdir_analysis/anal${MEM3}.ocean_velocity.res.nc ]; then
    ln -f $workdir_analysis/anal${MEM3}.ocean_velocity.res.nc ${workdir2}/INPUT/ocean_velocity.res.nc
  else
    echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_velocity.res.nc"
    exit 1
  fi

  ##################################################################################
  # SBC Restart Files:
  ##################################################################################
  if [ -f $workdir_analysis/anal${MEM3}.ocean_sbc.res.nc ]; then
    ln -f $workdir_analysis/anal${MEM3}.ocean_sbc.res.nc ${workdir2}/INPUT/ocean_sbc.res.nc
  else
    echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_sbc.res.nc"
    exit 1
  fi

  ##################################################################################
  # eta Restart Files:
  ##################################################################################
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    if [ -f $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc ]; then
      ln -f $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc ${workdir2}/INPUT/ocean_barotropic.res.nc
    else
      echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc"
      exit 1
    fi
  fi

fi

#Check for timestamp restart file before continuing...
#It should have been in previous RESTART directory
if [ ! -f ${workdir2}/INPUT/$rtype.res ]; then
  echo "LINK ERROR: ${workdir2}/INPUT/$rtype.res does not exist. Unfortunately."
  exit 1
fi

#############################################
# set up the surface forcing fluxes  #STEVE: need t60 version for R2 and 20CRv2:
#############################################
mkDlySBCnc4=mkDlySBCnc4
mDS=mDS_R2CR.bash

# Use only the mean fluxes
if [ "$USE_MFLX" -eq "1" ]; then
  # OVERRIDE: !!!!*!*!*!*!*!*!
  mkDlySBCnc4=mkDlySBCnc4
  mDS=mDS.bash
fi

cd ${workdir}
# /usr/bin/perl -i.bak -pe "s/L_YR=\d+/L_YR=$IY/" ./mDS.sh
cp $SBCDIR/$mkDlySBCnc4 .
cp $SBCDIR/mkDlySss4nci .
cp $SBCDIR/mkDlySst4i .
ln -f $INPUT_INIT/grid_spec.nc .    # Model grid definitions
ln -f $INPUT_INIT/../salt12.nc .    # Salinity climatology file

#if [ "$isfirst" -eq "1" ]; then
#  # change first day so it doesn't need data from 1984:
#  ID=02
#fi

echo "$IY  $IM  $ID  $IH  $IN  $IS" >> ${workdir}/time_stamp.out
echo "In: $PWD"
echo "Running:"
cp $wgrib ${workdir}
if [ $USE_EFLX -eq 1 ]; then
  echo "${SBCDIR}/$mDS ${workdir2} $FLXDIR $FLXDIR2 $SSTDIR $days"
  ${SBCDIR}/$mDS ${workdir2} $FLXDIR $FLXDIR2/$MEM2 $SSTDIR $days
  #STEVE: need to format each ensemble member exactly as R2 (to make my life easier)
elif [ $USE_MFLX -eq 1 ]; then
  echo "${SBCDIR}/$mDS ${workdir2} $FLXDIR $days"
  ${SBCDIR}/$mDS ${workdir2} $FLXDIR $SSTDIR $days
fi

# link to each model working directory
for file in `ls -d ${workdir}/*_sfc_restore.nc`; do
  ln -f $file ${workdir2}/INPUT
done

#STEVE: link the surface flux files to each model INPUT directory
echo "Linking surface forcings for RA2_daily_*.nc, member: $MEM2"
for file in `ls -d ${workdir}/RA2_daily_*.nc`; do
  ln -f $file ${workdir2}/INPUT
done

exit 0
