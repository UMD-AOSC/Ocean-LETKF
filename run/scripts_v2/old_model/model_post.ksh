#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Model Post step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
echo "Posting forecast ${FCST}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_post/${MEMBERID}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
cd ${workdir}

echo "This is the Model post output for forecast ${FCST} for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > post_${FCST}.out

#STEVE: active code:
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
ISLOTL=`printf %.2d $FCST`

IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the date for the next analysis cycle
date=/bin/date
inc=5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

# Update the ISLOT date
date=/bin/date
sinc=`expr $FCST - 1`
sinc_units=days
NY=`$date -d "$TY-$TM-$TD $sinc $sinc_units" +%Y`
NM=`$date -d "$TY-$TM-$TD $sinc $sinc_units" +%m`
ND=`$date -d "$TY-$TM-$TD $sinc $sinc_units" +%d`
NH=`$date -d "$TY-$TM-$TD $sinc $sinc_units" +%H`
IY=$NY
IM=$NM
ID=$ND
IH=$NH
echo "processing cycle as: $IY$IM$ID$IH"

# Update the "Next" date:
inc=5
inc_units=days
NY=`$date -d "$TY-$TM-$TD $inc $inc_units" +%Y`
NM=`$date -d "$TY-$TM-$TD $inc $inc_units" +%m`
ND=`$date -d "$TY-$TM-$TD $inc $inc_units" +%d`
NH=`$date -d "$TY-$TM-$TD $inc $inc_units" +%H`
TY=$NY
TM=$NM
TD=$ND
TH=$NH


mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc $OUTPUT/gues/$MEM3
mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  $OUTPUT/gues/$MEM3
mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       $OUTPUT/gues/$MEM3

#STEVE: copy files from model output to new gues files
#if [ -f ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc ]; then
#  mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc $OUTPUT/gues/$MEM3
#else
#  echo "${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc does not exist after model run. EXITING..."
#  exit 1
#fi
#if [ -f ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc ]; then
#  mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  $OUTPUT/gues/$MEM3
#else
#  echo "${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc does not exist after model run. EXITING..."
#  exit 1
#fi
##if [ $USE_SFC ]; then
#  if [ -f ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc ]; then
#    mv ${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc     $OUTPUT/gues/$MEM3
#  else
#    echo "${workdir2}/RESTART/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc does not exist after model run. EXITING..."
#    exit 1
#  fi
##fi

exit
