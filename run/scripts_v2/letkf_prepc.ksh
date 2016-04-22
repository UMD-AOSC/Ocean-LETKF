#!/bin/ksh --login
# This script prepares background and observation data for letkf 
# Tripolar background grid is converted to spherical if necessary
#
# This 'compact' version runs all processing for each members (not each member & timeslot)
#
# Written by Dr. Stephen G. Penny
#
set -ex

# Load the modules
source $MODULESHOME/init/ksh
#module use /sw/eslogin-c3/modulefiles
module load nco
ncatted=ncatted

echo "LETKF preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir_fcst=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/RESTART
workdir_hist=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/history
workdir_inpt=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/INPUT
workdir2=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir2}

MEM3=`printf %.3d ${MEMBERID}`

echo "USE_ALTIMETRY == $USE_ALTIMETRY"
if [ $USE_ALTIMETRY -eq "1" ]; then
  echo "Assimilating Altimetry!"
fi

#STEVE: active code:
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00
YYYYMMDD=$IY$IM$ID

FCST=0
while test $FCST -lt $NSLOTS
do
  FCST=`expr $FCST + 1`
  ISLOT2=`printf %.2d ${FCST}`
  echo "Posting forecast ${ISLOT2}"
  workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf_prep/${MEMBERID}/${ISLOT2}
  mkdir -p ${workdir}
  cd ${workdir}
  export PATH=$PATH:`pwd` #STEVE: need to access nco executables copied to here, since modules are not available.
  echo "This is the LETKF preparation step for member ${MEMBERID}, for cycle ${YYYYMMDDHH}" > $workdir/letkf_prep_${MEMBERID}.out

  # Update the ISLOT date
  date=/bin/date
  #sinc=`expr $FCST - 1` #STEVE: use this if I want the dates offset
  sinc=1
  sinc_units=days
  NY=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%Y`
  NM=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%m`
  ND=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%d`
  NH=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%H`
  IY=$NY
  IM=$NM
  ID=$ND
  IH=$NH
  echo "processing cycle as: $IY$IM$ID$IH"

  # FIRST, link the background model data

  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ $USE_ALTIMETRY -eq "1" ]; then
    ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_barotropic.res.nc gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi
  if [ "$USE_MLD" -eq "1" ]; then
    echo "Using mixed layer depth localization, USE_MLD = $USE_MLD"
    echo "ln -f ${workdir_hist}/${YYYYMMDD}.ocean_TS.nc $workdir2/gs${ISLOT2}$MEM3.ocean_TS.nc"
    ln -f ${workdir_hist}/${YYYYMMDD}.ocean_TS.nc $workdir2/gs${ISLOT2}$MEM3.ocean_TS.nc
  else
    echo "Not using mixed layer depth localization, USE_MLD = $USE_MLD"
  fi

  #STEVE: add 'fill value' to netcdf files
  nca=/sw/xe6/nco/4.0.8/sles11.1_netcdf4.2.0_gnu4.7.0/bin/ncatted
# nca3=/sw/eslogin-c3/nco/4.5.2/sles11.3_gnu5.1.0/bin/ncatted
# if [ -f $nca ]; then
    ncatted=$nca
# elif [ -f $nca3 ]; then
#   ncatted=$nca3
# fi
  $ncatted -O -a _FillValue,temp,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  $ncatted -O -a _FillValue,salt,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  $ncatted -O -a _FillValue,u,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  $ncatted -O -a _FillValue,v,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  $ncatted -O -a _FillValue,sea_lev,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ $USE_ALTIMETRY -eq "1" ]; then
    $ncatted -O -a _FillValue,eta_t,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi

  # FIRST (a), Link background files to the letkf working directory:

  #STEVE: need a template for the output analysis files:
  #       MUST COPY, NOT LINK
  if [ "${ISLOT2}" -eq "${ATIME}" ]; then #STEVE: trying to only do it once per member...
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/anal$MEM3.ocean_temp_salt.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/anal$MEM3.ocean_velocity.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/anal$MEM3.ocean_sbc.res.nc
    if [ $USE_ALTIMETRY -eq "1" ]; then
      cp ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc       ${workdir2}/anal$MEM3.ocean_barotropic.res.nc
    fi
    #STEVE: may want to zero these out for peace of mind...
  fi

  # SECOND (b), Run separate obs operators on each dataset:
  echo "Processing obs for timestep: $IY$IM$ID$IH to obs${ISLOT2}${MEM3}.dat"
  ln -f ${INPUT_INIT}/grid_spec.nc .
  obsoutfile=${workdir2}/obs${ISLOT2}${MEM3}.dat
  rm -f $obsoutfile
  touch $obsoutfile

# exec_command="aprun "
  exec_command="./"
  #STEVE: if the directory was provided for a certain obs type, then apply the observation operator to it.
  if [ -d "${OBSDIR_T}" ]; then
    # copying data necessary for potential to in situ temperature conversion:
    if [ -f ${PT2IS_DATA_DIR}/$PT2IS_DATA_FILE ]; then
      cp ${PT2IS_DATA_DIR}/$PT2IS_DATA_FILE . 
    else
      echo "SOURCE: ${PT2IS_DATA_DIR}/$PT2IS_DATA_FILE"
      echo "The necessary file, $PT2IS_DATA_FILE, required for potential temperature to in situ conversion is missing."
      exit 82
    fi
    # copying executable to local working directory:
    cp ${ODIR}/$OBSOPexe_T .
    echo "Running $OBSOPexe_T"
    # link observation data from repository:
    ln -fs ${OBSDIR_T}/${IY}${IM}${ID}${OBS_T_SUFFIX} obsin_t.nc
    # Check for synthetic observations, use for in situ to potential temperature conversion:
    if [ "TEST$OBSDIR_SM" != "TEST" ]; then
      ln -fs $OBSDIR_SM/${IY}${IM}${ID}${OBS_S_SUFFIX} obsin2_s.nc
      obsin2opt='-obsin2 obsin2_s.nc'
    fi
    # run observation operator:
    ${exec_command}$OBSOPexe_T -obsin obsin_t.nc $obsin2opt -obsout obsout_t.dat -gues gs${ISLOT2}$MEM3 $CONV_OPT
    cp $obsoutfile scratch.dat
    cat scratch.dat obsout_t.dat > $obsoutfile 
    # cleanup:
    rm -f $OBSOPexe_T
    rm -f $PT2IS_DATA_FILE
    rm -f scratch.dat
    rm -f obsin_t.nc
    rm -f obsout_t.dat
  else
    echo "Skipping Temperature Profile observation data."
  fi
  if [ -d "${OBSDIR_S}" ]; then
    cp ${ODIR}/$OBSOPexe_S .
    echo "Running $OBSOPexe_S"
    ln -fs ${OBSDIR_S}/${IY}${IM}${ID}${OBS_S_SUFFIX} obsin_s.nc
    ${exec_command}${OBSOPexe_S} -obsin obsin_s.nc -obsout obsout_s.dat -gues gs${ISLOT2}$MEM3
    rm -f $OBSOPexe_S
    cp $obsoutfile scratch.dat
    cat scratch.dat obsout_s.dat > $obsoutfile 
    # cleanup:
    rm -f $OBSOPexe_S
    rm -f scratch.dat
    rm -f obsin_s.nc
    rm -f obsout_s.dat
  else
    echo "Skipping Salinity Profile observation data."
  fi
  if [ -d "${OBSDIR_SST}" ]; then
    cp ${ODIR}/$OBSOPexe_SST .
    echo "Running $OBSOPexe_SST"
    JDAY=`$date -d "$IY-$IM-$ID" +%j`
    echo "Searching directory: $OBSDIR_SST"
    echo "Searching directory: ${OBSDIR_SST}"
    OBSIN=`ls ${OBSDIR_SST}/*${IY}${JDAY}${OBS_SST_SUFFIX}`  #STEVE: there should only be one date with this suffix.
#########
    echo "TEST trailing dir error:"
    echo "$OBSIN"
    echo "${OBSIN}"
    echo "(both should be identical.)"
##########
    ln -fs ${OBSIN} obsin_sst.nc
    ${exec_command}$OBSOPexe_SST -obsin obsin_sst.nc -obsout obsout_sst.dat -gues gs${ISLOT2}$MEM3 -superob ".true."
#   ${exec_command}$OBSOPexe_SST -obsin obsin_sst.nc -obsout obsout_sst.dat -gues gs${ISLOT2}$MEM3 -thin $ThinSST
    rm -f $OBSOPexe_SST
    cp $obsoutfile scratch.dat
    cat scratch.dat obsout_sst.dat > $obsoutfile 
    # cleanup:
    rm -f $OBSOPexe_SST
    rm -f scratch.dat
    rm -f obsin_sst.nc
    rm -f obsout_sst.dat
  else
    echo "Skipping Sea Surface Temperature observation data."
  fi
  
  # Link background files to the letkf working directory:
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc ${workdir2}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi

  #STEVE: the hard-link limit is running out (65000 max), so best to delete unnecessary links
  rm -f grid_spec.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    rm -f $altimetry_climatology_file
  fi

done #DONE FCST loop

exit 0
