#!/bin/ksh --login
# STEVE: test these to replace copying wgrib and nco exe's to the local directory below:

#source $MODULESHOME/init/ksh
#module load wgrib/1.8.1.0b
##module load wgrib

set -ex

echo "Model preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_prep/${MEMBERID}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
mkdir -p ${workdir2}
cd ${workdir}
export PATH=$PATH:`pwd` #STEVE: need to access wgrib and nco executables copied to here, since modules are not available.

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
  #STEVE: getting the members started with identical initial conditions:
  ln -f $INPUT_INIT/* ${workdir2}/INPUT/
#STEVE: this links these files:
# ln -fs $INPUT_INIT/chl.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/grid_spec.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/namelist ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/ncar_precip_clim.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/RUNOFF.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/salt12.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/SNOW.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/sst_ice_clim.nc ${workdir2}/INPUT  
# ln -fs $INPUT_INIT/input.nml ${workdir2}/INPUT

  # Check for compressed files (possible if re-running this timestep):
  workdir0=${EXP_DATA}/$PY$PM$PD$PH/model/${MEMBERID}
  for file in `ls -d ${workdir0}/RESTART/ocean_*res*.gz`; 
  do
    echo "unzipping $file ..."
    gunzip $file
  done

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

  #STEVE: copy restart files from last timestep (need to set this up manually for first timestep)
  #STEVE: if an analysis exists, replace the primary restart files with the analysis
  ##################################################################################
  # TEMPERATURE AND SALINITY Restart Files:
  ##################################################################################
  if [[ "$datype" == "HYBRID" && ! -d ${EXP_DATA}/$PY$PM$PD$PH/go ]]; then
    echo "datype=$datype, setting working directory for HYBRID..."
    #STEVE: if using the hybrid, then get the recentered data
    echo "Using Hybrid-LETKF Analysis from previous step as initial conditions for ocean_temp_salt..."
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/hybrid
  elif [ "$datype" == "GODAS" ]; then
    echo "datype=$datype, setting working directory for GODAS..."
    echo "Using GODAS-SOLO Analysis from previous step as initial conditions for ocean_temp_salt..."
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/godas_solo
  else
    echo "datype=$datype, setting working directory for LETKF..."
    echo "Using Standard-LETKF Analysis from previous step as initial conditions..."
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/letkf
  fi

  if [ "$datype" == "GODAS" ]; then # GODAS GODAS GODAS !!!
    echo "datype=$datype, running for GODAS..."
    if [ -f $workdir_analysis/RESTART/ocean_temp_salt.res.nc ]; then
      ln -f $workdir_analysis/RESTART/ocean_temp_salt.res.nc ${workdir2}/INPUT/ocean_temp_salt.res.nc
    else
     echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/RESTART/ocean_temp_salt.res.nc"
      exit 1
    fi

    echo "Using Model forecast from previous step as initial conditions for ocean_velocity and ocean_sbc ..."
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/model/$MEM2

    # VELOCITY Restart Files:
    srcfile=$workdir_analysis/RESTART/ocean_velocity.res.nc
    if [ -f $srcfile ]; then
      ln -f $srcfile ${workdir2}/INPUT/ocean_velocity.res.nc
    else
      echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $srcfile"
      exit 1
    fi

    # SBC Restart Files:
    srcfile=$workdir_analysis/RESTART/ocean_sbc.res.nc
    if [ -f $srcfile ]; then
      ln -f $srcfile ${workdir2}/INPUT/ocean_sbc.res.nc
    else
      echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $srcfile"
      exit 1
    fi

    # Barotropic Restart Files:
    if [ $USE_ALTIMETRY -eq "1" ]; then
      srcfile=$workdir_analysis/RESTART/ocean_barotropic.res.nc
      if [ -f $srcfile ]; then
        ln -f $srcfile ${workdir2}/INPUT/ocean_barotropic.res.nc
      else
        echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $srcfile"
        exit 1
      fi
    fi

  else # LETKF LETKF LETKF !!! -or- HYBRID HYBRID HYBRID !!!
    echo "datype=$datype, running for LETKF..."

    # Check if the HYBRID temp_salt analysis files have been archived
    if [ -f $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc.gz ]; then
      echo "Unzipping..."
      gunzip $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc.gz
    fi

    # Continue with linking files...
    if [ -f $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc ]; then
      ln -f $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc ${workdir2}/INPUT/ocean_temp_salt.res.nc
    else
      echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_temp_salt.res.nc"
      exit 1
    fi

    echo "Using Standard-LETKF Analysis from previous step as initial conditions for ocean_velocity and ocean_sbc ..."
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/letkf

    # Check if the LETKF files for non-temp_salt data have been archived
    if [ -f $workdir_analysis/anal${MEM3}.ocean_velocity.res.nc.gz ]; then
      echo "Unzipping..."
      gunzip $workdir_analysis/anal${MEM3}.ocean_velocity.res.nc.gz
    fi
    if [ -f $workdir_analysis/anal${MEM3}.ocean_sbc.res.nc.gz ]; then
      echo "Unzipping..."
      gunzip $workdir_analysis/anal${MEM3}.ocean_sbc.res.nc.gz
    fi
    if [ -f $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc.gz ]; then
      echo "Unzipping..."
      gunzip $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc.gz
    fi

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
    # Barotropic Restart Files:
    ##################################################################################
    if [ $USE_ALTIMETRY -eq "1" ]; then
      if [ -f $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc ]; then
        ln -f $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc ${workdir2}/INPUT/ocean_barotropic.res.nc
      else
        echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_analysis/anal${MEM3}.ocean_barotropic.res.nc"
        exit 1
      fi
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

# OVERRIDE: !!!!*!*!*!*!*!*!
# Use only the mean fluxes
if [ "$USE_MFLX" -eq "1" ]; then
  mkDlySBCnc4=mkDlySBCnc4
  mDS=mDS.bash
fi

cd ${workdir}
# /usr/bin/perl -i.bak -pe "s/L_YR=\d+/L_YR=$IY/" ./mDS.sh
cp $SBCDIR/$mkDlySBCnc4 .
cp $SBCDIR/mkDlySss4nci .
cp $SBCDIR/mkDlySst4i .
ln -f $INPUT_INIT/grid_spec.nc .
ln -f $INPUT_INIT/../salt12.nc .

#if [ "$isfirst" -eq "1" ]; then
#  # change first day so it doesn't need data from 1984:
#  ID=02
#fi

echo "$IY  $IM  $ID  $IH  $IN  $IS" >> ${workdir}/time_stamp.out
echo "In: $PWD"
echo "Running:"
cp /sw/xe6/wgrib/1.8.1.0b/sles11.1_gnu4.3.4/bin/wgrib ${workdir}
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
