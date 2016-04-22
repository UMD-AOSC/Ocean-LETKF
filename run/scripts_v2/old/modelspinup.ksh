#!/bin/ksh --login
set -e

module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Model step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
echo "at time: ${NNSS}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_${NNSS}/${MEMBERID}
mkdir -p ${workdir}
cd ${workdir}

#for fcst in ${FCSTS}; do
#  echo "This is the model output for forecast ${fcst} for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > model_${fcst}.out
#done

#STEVE: active code:
MEM3=`printf %.3d ${MEMBERID}`
MEM2=`printf %.2d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

#STEVE: need mpiexec to run across all procs
echo $0
if [ ! -f ${workdir}/INPUT/ocean_temp_salt.res.nc ]; then
  echo "Input file does not exist... but it should: ${workdir}/INPUT/ocean_temp_salt.res.nc"
  exit 1
fi

echo "Running Member::$MEM3"
ln -fs $MOM4dir/$MOM4exe .

cd ${workdir}

echo "Running model ensemble members..."

#STEVE: Submit model run
export MPI_VERBOSE=1
export MPI_DISPLAY_SETTINGS=1
export MPI_BUFS_PER_PROC=128
export MPI_BUFS_PER_HOST=128
export MPI_IB_RAILS=2
export MPI_GROUP_MAX=128
cp $MOM4run ${workdir} #STEVE: mostly for running in debugging purposes
$MOM4run ${TMPDIR}/model_${NNSS} ${MEMBERID} mom4p1_$mtype $PBS_NP

# Output a file containing
echo "$MOM4run ${TMPDIR}/model_${NNSS} ${MEMBERID} mom4p1_$mtype $PBS_NP" > model_run_command.txt

exit 0
