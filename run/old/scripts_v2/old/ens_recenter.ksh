#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Ensemble Recentering step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
echo "at time: ${NNSS}"
workdir0=${EXP_DATA}/${YYYYMMDDHH}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model_${NNSS}
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_${NNSS}/meansprd
mkdir -p ${workdir}

# Change to model directory
cd ${workdir0}/model_${NNSS}

echo "This is the model spinup preparation step for member ${MEMBERID} for cycle ${YYYYMMDDHH}_${NNSS}" > model_prep.out

#STEVE: active code:
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

#STEVE: create prep output netcdf files
workdir_prev=${EXP_DATA}/${YYYYMMDDHH}/model_${PNNSS}/meansprd
if [ -f $workdir_prev/mean.ocean_temp_salt.res.nc ]; then
  cp $workdir_prev/mean.ocean_temp_salt.res.nc ${workdir0}/model_${NNSS}/
  cp $workdir_prev/std.ocean_temp_salt.res.nc ${workdir0}/model_${NNSS}/
else
  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/mean.ocean_temp_salt.res.nc"
  exit 1
fi
if [ -f $workdir_prev/mean.ocean_velocity.res.nc ]; then
  cp $workdir_prev/mean.ocean_velocity.res.nc ${workdir0}/model_${NNSS}/
  cp $workdir_prev/std.ocean_velocity.res.nc ${workdir0}/model_${NNSS}/
else
  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/mean.ocean_velocity.res.nc"
  exit 1
fi
if [ -f $workdir_prev/mean.ocean_sbc.res.nc ]; then
  cp $workdir_prev/mean.ocean_sbc.res.nc ${workdir0}/model_${NNSS}/
  cp $workdir_prev/std.ocean_sbc.res.nc ${workdir0}/model_${NNSS}/
else
  echo "ERROR: ANALYSIS FILE DOES NOT EXIST: $workdir_prev/mean.ocean_sbc.res.nc"
  exit 1
fi

# Do the recentering:
${workdir0}/msmr.x -basedir ${workdir0}/model_00/meansprd -indir ${workdir0}/model_${PNNSS} -outdir ${workdir0}/model_${NNSS} -gsdir ${workdir0}/model_00/01/INPUT

# Store the mean and std netcdf files
mv $workdir0/model_${NNSS}/mean.*.nc $workdir
mv $workdir0/model_${NNSS}/std.*.nc  $workdir
cp $workdir_prev/*.ctl $workdir

exit 0
