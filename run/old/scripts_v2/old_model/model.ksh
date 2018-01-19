#!/bin/ksh --login
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
echo "Running Member::$MEM3"
ln -fs $MOM4dir/$MOM4exe .

# Make sure to use correct date for model run
cd ${workdir}/INPUT
# Get perl script to fix coupler restart files to the analysis time
#(STEVE: this shouldn't be necessary, but just in case it gets off from where it is supposed to be...)
cp $CDIR/fix_coupler.*.pl .

#STEVE: make sure the couple.res and coupler.intermediate.res have the correct starting times: (i.e. the last day of the last 5-day analysis cycle)
#use perl script to edit ymdhns of 3rd line of coupler.res
/usr/bin/perl fix_coupler.res.pl $AY $AM $AD $AH $AN $AS
echo "Converting to date $AY $AM $AD $AH $AN $AS"
mv coupler.res coupler.res.old
mv coupler.res.new coupler.res
cat coupler.res

#use perl script to edit ymdhns of 1st line of coupler.intermediate.res
/usr/bin/perl fix_coupler.intermediate.res.pl $AY $AM $AD $AH $AN $AS
echo "Converting to date $AY $AM $AD $AH $AN $AS"
mv coupler.intermediate.res coupler.intermediate.res.old
mv coupler.intermediate.res.new coupler.intermediate.res
cat coupler.intermediate.res

cd ${workdir}

echo "Running model ensemble members..."

#STEVE: Submit model run
cp $MOM4run ${workdir} #STEVE: mostly for running in debugging purposes
$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_coupled $PBS_NP

rsrtOUT=$OUTPUT/rsrt/$TY$TM$TD/$MEM3
mkdir -p $rsrtOUT
#STEVE: Move necessary files for doing a run from restart
#       These first two are important for tracking model time:
#rm -f ${workdir}/RESTART/RA2_daily_*.nc #STEVE: just in case they end up here, don't need the sfc forcing files
#if [ -f ${workdir}/RESTART/coupler.res ]; then
  mv ${workdir}/RESTART/coupler.res              $rsrtOUT/
#fi
#if [ -f ${workdir}/RESTART/coupler.intermediate.res ]; then
  mv ${workdir}/RESTART/coupler.intermediate.res $rsrtOUT/
#fi

mv ${workdir}/RESTART/[a-z]* $rsrtOUT/

#STEVE: also move the dated files for safer restarting:
cp ${workdir}/RESTART/$TY$TM$TD.* $rsrtOUT/
#if [ ! -f ${workdir}/RESTART/ocean_temp_salt.res.nc ]; then
 # rm -f ${workdir}/RESTART/[1-2]*
#fi

exit
