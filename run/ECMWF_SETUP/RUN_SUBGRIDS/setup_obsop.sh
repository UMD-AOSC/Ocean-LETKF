#!/bin/bash
set -ex
date=/bin/date

source params.sh

function run_obsop {

# local EXEDIR=$1
# local EXE_tprof=$2
# local EXE_sprof=$3
# local ISLOT2=$4
# local MEM3=$5

  # Process temperature profiles
  args="-obsin $obsinfile -obsout obs${ISLOT2}${MEM3}.tprof.dat -obs_name $GLOBAL_TPROF_OBS_NAME -hxb_name $GLOBAL_TPROF_Hxb_NAME"
  $EXEDIR/$EXE_tprof $args

  # Process salinity profiles
  args="-obsin $obsinfile -obsout obs${ISLOT2}${MEM3}.sprof.dat -obs_name $GLOBAL_SPROF_OBS_NAME -hxb_name $GLOBAL_SPROF_Hxb_NAME"
  $EXEDIR/$EXE_sprof $args

  # cat files together and move to working directory for letkf:
  `cat obs${ISLOT2}${MEM3}.tprof.dat obs${ISLOT2}${MEM3}.sprof.dat > obs${ISLOT2}${MEM3}.dat`

  # Remove individual files:
  rm -f obs${ISLOT2}${MEM3}.tprof.dat
  rm -f obs${ISLOT2}${MEM3}.sprof.dat

}

# Set up the observation/innovation .dat files by running the letkf obsop_tprof.x / obsop_sprof.x routines.

MEM=$GLOBAL_MEM_START
MEMBERS=$GLOBAL_MEMBERS

IYYYY=$GLOBAL_IYYYY
IMM=$GLOBAL_IMM
IDD=$GLOBAL_IDD

inc=$GLOBAL_inc
inc_units=$GLOBAL_inc_units

FILE_PREFIX=$GLOBAL_OBS_FILE_PREFIX #0001_nrt_5_2_
FILE_SUFFIX=$GLOBAL_OBS_FILE_SUFFIX #_profb_01_fdbk.nc

ISLOT=1  # Use if there are multiple observation time bins
ISLOT2=`printf "%02d" $ISLOT`

#SCRATCH=/scratch/rd/dasp
SCRATCH=$GLOBAL_SCRATCH #/lus/snx11064/scratch/rd/dasp
EXPNAME=$GLOBAL_EXPNAME #CERA20C_50mem

LETKF_ROOT=$GLOBAL_LETKF_ROOT
name=$GLOBAL_LETKF_NAME
EXEDIR=$LETKF_ROOT/build/build_obsop/${name}.build
EXE_tprof=obsop.$name.tprof.x
EXE_sprof=obsop.$name.sprof.x

GRDDIR=$GLOBAL_GRDDIR
GRDFILE=$GLOBAL_GRDFILE

DSTDIR=$SCRATCH/$EXPNAME/WORK/obs
mkdir -p $DSTDIR
cp -f $GRDDIR/$GRDFILE $DSTDIR/
cd $DSTDIR

YYYY=$IYYYY
MM=$IMM
DD=$IDD
  
NY=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%Y`
NM=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%m`
ND=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%d`

echo "For reference, the current date is:"
echo "$YYYY-$MM-$DD"
echo "The next date is:"
echo "$NY-$NM-$ND"

DATE_RANGE= #${YYYY}${MM}${DD}_${NY}${NM}${ND}   #STEVE: current files don't use a date range

while [ $MEM -lt $MEMBERS ]; do

  MEM1=`printf "%01d" $MEM`
  MEM2=`printf "%02d" $MEM`

  FILE=sorted_${MEM2}_${FILE_PREFIX}${DATE_RANGE}${FILE_SUFFIX}
  obsinfile=$FILE
  echo "The input observation/innovation file is: $obsinfile"

  MEM=`expr $MEM + 1`
  MEM3=`printf "%03d" $MEM`

  run_obsop & #

# echo "setup_obsop.sh:: temporarily ending early..."
# exit 1

done

time wait
