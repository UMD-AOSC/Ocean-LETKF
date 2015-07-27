#!/bin/ksh --login
set -e

#module load intel
#module load netcdf
#module load nco

echo "Model step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
cd ${workdir}

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

# Start date of next analysis cycle
date=/bin/date
inc=$days #5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

#STEVE: in case this is a re-run, cleanup intermediate files:
rm -f ${workdir}/ocean_TS.nc*
rm -f ${workdir}/ocean_UV.nc*

#STEVE: need mpiexec to run across all procs
echo $0
if [ ! -f ${workdir}/INPUT/ocean_temp_salt.res.nc ]; then
  echo "Input file does not exist... but it should: ${workdir}/INPUT/ocean_temp_salt.res.nc"
  exit 1
fi

echo "Running Member::$MEM3"
cd ${workdir}
echo "Running model ensemble members..."

#STEVE: Submit model run
cp $MOM4run ${workdir} #STEVE: mostly for running in debugging purposes
# (Calls aprun within this script)
$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_$mtype $PBS_NP

# Output a file containing
echo "$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_$mtype $PBS_NP" > model_run_command.txt
echo "model run for ${MEMBERID} is finished." > ${workdir}/model.out

exit 0
