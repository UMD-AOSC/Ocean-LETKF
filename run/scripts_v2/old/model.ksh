#!/bin/ksh --login
set -e

module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Model step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
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

# Update the date for the next analysis
date=/bin/date
ainc=4
ainc_units=days
AY=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%Y`
AM=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%m`
AD=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%d`
AH=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%H`
AN=$IN
AS=$IS

# Start date of next analysis cycle
inc=5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

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
$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_coupled $PBS_NP

# Output a file containing
echo "$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_coupled $PBS_NP" > model_run_command.txt

exit 0
