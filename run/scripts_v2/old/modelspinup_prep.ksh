#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Model preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
echo "at time: ${NNSS}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_prep_${NNSS}/${MEMBERID}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model_${NNSS}/${MEMBERID}
mkdir -p ${workdir}
mkdir -p ${workdir2}
cd ${workdir}

echo "This is the model spinup preparation step for member ${MEMBERID} for cycle ${YYYYMMDDHH}_${NNSS}" > model_prep.out

#STEVE: active code:
USE_MFLX=0
USE_EFLX=1
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

echo "USE_MFLX = $USE_MFLX"
echo "USE_EFLX = $USE_EFLX"

# t190 to t254 cutoff date (t254:: 2/14/2012)
t254_start=2012021500
if [ ${YYYYMMDDHH} -lt $t254_start ]; then
  FLUXDIR=$FLXDIR
  mkDlyPSBCnc4=mkDlyPSBCnc4_t190.x
  mPS=mPS_t190.csh
else
  FLUXDIR=$FLXDIR2
  mkDlyPSBCnc4=mkDlyPSBCnc4_t254.x
  mPS=mPS_t254.csh
fi

# OVERRIDE: !!!!*!*!*!*!*!*!
if [ $USE_MFLX -eq 1 ]; then
  FLUXDIR=/scratch2/portfolios/NCEPDEV/climate/noscrub/David.Behringer/SBC/R2/DAILYnc
  mkDlyPSBCnc4=mkDlyPSBCnc4.x
  mPS=mPS.csh
fi

echo "Using fluxes from: $FLUXDIR"

mkdir -p ${workdir2}/INPUT
mkdir -p ${workdir2}/RESTART

#STEVE: copy INPUT/RESTART files to here
#STEVE: getting it started with identical initial conditions:
#cp $INPUT/$MEM2/INPUT/* ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/cfc.bc.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/chl.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/grid_spec.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/namelist ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/ncar_precip_clim.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/RUNOFF.nc ${workdir2}/INPUT  
#ln -fs $EXP_DATA/INPUT/sala.mom ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/salt12.i3e ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/salt12.i3eS ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/salt12.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/SNOW.nc ${workdir2}/INPUT  
ln -fs $EXP_DATA/INPUT/sst_ice_clim.nc ${workdir2}/INPUT  
#ln -fs $EXP_DATA/INPUT/svv.mom ${workdir2}/INPUT  
#ln -fs $EXP_DATA/INPUT/tmpa.mom ${workdir2}/INPUT  
#ln -fs $EXP_DATA/INPUT/tvv.mom ${workdir2}/INPUT  

#STEVE: copy restart files from last timestep (need to set this up manually for first timestep)
workdir0=${EXP_DATA}/${YYYYMMDDHH}/model_${PNNSS}/$MEM2
for file in `ls -d ${workdir0}/RESTART/[a-z]*`; do
  echo "linking $file to ${workdir2}/INPUT ..."
  ln -fs $file ${workdir2}/INPUT
done

#STEVE: replace the primary restart files with the final perturbations re-centered at the original model mean
#STEVE: NOTE! This must be done in a separate script that is called once only
#workdir_prev=${EXP_DATA}/${YYYYMMDDHH}/model_${PNNSS}/${MEMBERID}/RESTART
#if [ -f $workdir_prev/ocean_temp_salt.res.nc ]; then
#  ln -fs $workdir_prev/ocean_temp_salt.res.nc ${workdir2}/INPUT/ocean_temp_salt.res.nc
#else
#  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/ocean_temp_salt.res.nc"
#  exit 1
#fi
#if [ -f $workdir_prev/ocean_velocity.res.nc ]; then
#  ln -fs $workdir_prev/ocean_velocity.res.nc ${workdir2}/INPUT/ocean_velocity.res.nc
#else
#  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/ocean_velocity.res.nc"
#  exit 1
#fi
#if [ -f $workdir_prev/ocean_sbc.res.nc ]; then
#  ln -fs $workdir_prev/ocean_sbc.res.nc ${workdir2}/INPUT/ocean_sbc.res.nc
#else
#  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/ocean_sbc.res.nc"
#  exit 1
#fi

## Get coupler files from previous input directory so it has the right times...

workdir_base=${EXP_DATA}/${YYYYMMDDHH}/model_00/$MEM2
if [ -f $workdir_base/INPUT/coupler.res ]; then
  ln -fs $workdir_base/INPUT/coupler.res ${workdir2}/INPUT/coupler.res
else
  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir0/INPUT/coupler.res"
  exit 1
fi

#if [ -f $workdir0/INPUT/coupler.intermediate.res ]; then
#  ln -fs $workdir0/INPUT/coupler.intermediate.res ${workdir2}/INPUT/coupler.intermediate.res
#else
#  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir0/INPUT/coupler.intermediate.res"
#  exit 1
#fi

#Check before continuing...
if [ ! -f ${workdir2}/INPUT/coupler.res ]; then
  echo "LINK ERROR: ${workdir2}/INPUT/coupler.res does not exist. Unfortunately."
  exit 1
fi

# set up the surface forcing fluxes
  cd ${workdir}
# /usr/bin/perl -i.bak -pe "s/L_YR=\d+/L_YR=$IY/" ./mPS.sh
  cp $SBCDIR/time_stamp.out ${workdir}
  ln -fs $SBCDIR/$mkDlyPSBCnc4 .
  ln -fs $SBCDIR/mkDlySss4nci .
  ln -fs $SBCDIR/mkDlySst4i .
  echo "$IY  $IM  $ID  $IH  $IN  $IS" >> ${workdir}/time_stamp.out
  echo "In: $PWD"
  echo "Running:"
if [ $USE_EFLX -eq 1 ]; then
  echo "${SBCDIR}/mPS.sh ${workdir2} $FLUXDIR/$MEM2"
  ${SBCDIR}/$mPS ${workdir2} $FLUXDIR/$MEM2 
elif [ $USE_MFLX -eq 1 ]; then
  echo "${SBCDIR}/mPS.sh ${workdir2} $FLUXDIR"
  ${SBCDIR}/$mPS ${workdir2} $FLUXDIR
fi
  # link to each model working directory
  for file in `ls -d ${workdir}/RA2_daily_*.nc`; do
    ln -fs $file ${workdir2}/INPUT
  done
  for file in `ls -d ${workdir}/*_sfc_restore.nc`; do
    ln -fs $file ${workdir2}/INPUT
  done

exit 0
