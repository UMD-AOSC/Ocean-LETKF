#!/bin/ksh --login
#module load mpt
#module load intel
#module load netcdf #/4.1.3-intel
#module load nco

set -e

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
    workdir_analysis=${EXP_DATA}/$PY$PM$PD$PH/model/00

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

# Establish previous model_prep working directory that contained RA2 wind files from previous step
workdir_prev=${EXP_DATA}/$PY$PM$PD$PH/model_prep/$MEM2

#########################################################################
# Modify the surface forcing with analyzed data, if it exists.
#
# STEVE: This implements a quasi-analysis of the surface forcing fields.
#        There is an effective bias correction to the surface forcing
#        fields enacted in the process.
#########################################################################

if [ -f ${workdir_analysis}/SFA_${MEM3}_daily_pres.nc ]; then
  echo "Surface analysis exists in LETKF output, and DO_SFCFLUXES=$DO_SFCFLUXES."
else
  echo "Surface analysis does not exists in LETKF output. Switching DO_SFCFLUXES=0."
  DO_SFCFLUXES=0
fi

if [ "$DO_SFCFLUXES" -eq 1 ]; then
  halflife=6 # (in cycles; 6 cycles ~ 30 days)
  hlcoeff=0.89089871814 # 2^(-1/halflife)
  one_minus_hlc=0.10910128186
  time_smoothing=0.05  #(gamma in paper, at the moment)

  # To start, get the original RA2 sfc files from the previous analysis cycle
  # I.E., NOT the previously adjusted ones created for the previous analysis cycle
  for file in `ls -d ${workdir_prev}/RA2_daily_*.nc`; do
    oldRA2file=`basename $file`
    newRA2file=`echo $oldRA2file | sed s/RA2_daily_/RA2_prev_daily_/`          #STEVE: FIX: now using the original RA2_ files instead of the SFC_ files
    ln -f $file ${workdir}/$newRA2file   #STEVE: this is the original SFC RA2 that should be used for the increment for bias correction
  done

  cd ${workdir}

  # FIRST, create the analysis increments
  echo "Doing Part 1: create SFI files..."
  if [ -f ${workdir_analysis}/SFA_${MEM3}_daily_pres.nc ]; then
    ln -f ${workdir_analysis}/SFA_${MEM3}_*.nc ${workdir}       #STEVE: this is the SFC 'analysis' generated by LETKF
    for file in `ls -d SFA_${MEM3}_*.nc`; do
    #STEVE: also, use nco function to compute perturbations
    file0=`echo $file | sed s/SFA_${MEM3}_/RA2_prev_/`          #STEVE: use the RA2_ files from the previous cycle as baseline
    file1=`echo $file | sed s/SFA_${MEM3}_/SFI_${MEM3}_/`       #STEVE: create analysis increment for this cycle
    #STEVE: need to either separate file2 or duplicate file1:
    ti=5
    #STEVE: get data for this specific time from background surface forcing
    ncks -d TIME,$ti $file0 $file0.$ti
    #STEVE: subtract at time t=5 days from start
    echo "adding $file minus $file0.$ti to get $file1..."
    #STEVE: SFA_ - RA2_prev = SFI_
    ncflint -w 1,-1 $file $file0.$ti $file1
    done
    rm *.nc.5
  else
    echo "There are no surface analysis fields from LETKF. Skipping creation of rolling SFR_*.nc"
    echo "(This should only occur on the first analysis cycle.)"
  fi

  #SECOND, if the rolling increment exists (i.e. every cycle after the first analysis cycle), 
  # then apply a weighted average of these two increments for the next timestep's sfc forcing
  # Choose the weighted average to give a 'halflife' to the new increment of about 30 days.
  echo "Doing Part 2: create SFR_new files..."
  if [ -f ${workdir_prev}/SFR_${MEM3}_daily_pres.nc ]; then
    echo "Doing Part 2a: create SFR_new files..."
    for file in `ls -d ${workdir_prev}/SFR_${MEM3}_*.nc`; do
      oldSFRfile=`basename $file`
      file2=`echo $oldSFRfile | sed s/SFR_${MEM3}_/SFI_${MEM3}_/`
      file3=$oldSFRfile #`echo $oldSFRfile | sed s/SFR_${MEM3}_/SFR_new_/`
      #STEVE: i.e. halflife_coefficient * old_rolling_increment + (1-hlc) * new_increment = new rolling increment
      ncflint -w $hlcoeff,$one_minus_hlc $file $file2 $file3
    done
  else
   echo "Doing Part 2b: copy SFI to SFR file..."
    echo "(because ${workdir_prev}/SFR_${MEM3}_daily_pres.nc does not exist)"
    #STEVE: this must mean it's the second analysis cycle, so just copy the new increment to the rolling increment file
    if [ -f SFI_${MEM3}_daily_pres.nc ]; then
      for file in `ls -d SFI_*.nc`; do
        file2=`echo $file | sed s/SFI_${MEM3}_/SFR_${MEM3}_/`
        ln -fs $file $file2
      done

      #STEVE: FIX: correct the SFR_new files to have the correct timestamp, for reading into grads
    else
      #STEVE: this must mean it's the first analysis cycle, so just move along...
      echo "There are no surface increment fields generated by model_prep. Skipping initial creation of SFR_*.nc"
      echo "(This should only occur on the first analysis cycle.)"
    fi
  fi

  # THIRD, add rolling analysis increments to each new surface field
  echo "Doing Part 3: create SFU files..."
  if [ -f SFR_${MEM3}_daily_pres.nc ]; then
    echo "Doing Part 3a: create SFU"
    for file in `ls -d SFR_${MEM3}_*.nc`; do
      file2=`echo $file | sed s/SFR_${MEM3}_/RA2_/`
      file3=`echo $file | sed s/SFR_${MEM3}_/SFU_/`    #Create surface update increment
      #STEVE: need to either separate file2 or duplicate file1:
      ti=0
      tf=6
      while test $ti -le $tf
      do
        #STEVE: get each new sfc forcing time, one by one
        ncks -d TIME,$ti $file2 $file2.$ti
        #STEVE: add a small percentage of the analysis increment to the new forcing fields,
        #       creating a simple recursive filter.
        ncflint -w 1,$time_smoothing $file2.$ti $file $file3.$ti
        ti=`expr $ti + 1`
      done
      ncrcat $file3.0 $file3.1 $file3.2 $file3.3 $file3.4 $file3.5 $file3.6 $file3
      #Fix time scaling applied to time coordinate by ncflint:
      ncks -A -v TIME $file2 $file3
      rm *.nc.[0-6]

      #STEVE: link the updated (SFU) sfc forcing to the RA2_ file in the model INPUT directory for this member.
      ln -fs $workdir/$file3 ${workdir2}/INPUT/$file2
    done
  else
    echo "Part 3b: the file SFR_${MEM3}_daily_pres.nc does not exist. (a test to see if SFR_${MEM3}_*.nc files are available)"
  fi
else
  #STEVE: otherwise, link the original surface flux files to each model INPUT directory
  echo "Using original (unaltered) surface forcings for RA2_daily_*.nc"
  for file in `ls -d ${workdir}/RA2_daily_*.nc`; do
    ln -fs $file ${workdir2}/INPUT
  done
fi

# STEVE: prevent error with too many hard links to grid_spec.nc:
rm -f ${workdir}/grid_spec.nc

exit 0
