#!/bin/bash
date=/bin/date

# Set up the observation/innovation .dat files by running the letkf obsop_tprof.x / obsop_sprof.x routines.

MEM=0
MEMBERS=5

IYYYY=2016
IMM=12
IDD=26

inc=6
inc_units='days'

ISLOT=1  # Use if there are multiple observation time bins
ISLOT2=`printf "%02d" $ISLOT`

EXEDIR=/home/rd/dasp/perm/Ocean-LETKF/build/build_obsop/ECMWF_nemo.nemovarHx0.build

GRDDIR=/home/rd/dasp/perm/DATA/STATIC
GRDFILE=mesh_mask.nc

name=ECMWF_nemo.nemovarHx0
EXE_tprof=obsop.$name.tprof.x
EXE_sprof=obsop.$name.sprof.x

FILE_PREFIX=0001_nrt_5_2_
#FILE_SUFFIX=_profbqc_01_fdbk.nc
FILE_SUFFIX=_profb_01_fdbk.nc

DSTDIR=/home/rd/dasp/scratch/WORK
mkdir -p $DSTDIR
cp -f $GRDDIR/$GRDFILE $DSTDIR/
cd $DSTDIR

YYYY=$IYYYY
MM=$IMM
DD=$IDD
while [ $MEM -lt $MEMBERS ]; do
  
  NY=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%Y`
  NM=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%m`
  ND=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%d`

  echo "For reference, the current date is:"
  echo "$YYYY-$MM-$DD"
  echo "The next date is:"
  echo "$NY-$NM-$ND"

  MEM1=`printf "%01d" $MEM`
  MEM2=`printf "%02d" $MEM`

  FILE=sorted_${MEM2}_${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
  obsinfile=$FILE
  echo "The input observation/innovation file is: $obsinfile"

  MEM=`expr $MEM + 1`
  MEM3=`printf "%03d" $MEM`

  # Process temperature profiles
  args="-obsin $obsinfile -obsout obs${ISLOT2}${MEM3}.tprof.dat"
  $EXEDIR/$EXE_tprof $args

  # Process salinity profiles
  args="-obsin $obsinfile -obsout obs${ISLOT2}${MEM3}.sprof.dat"
  $EXEDIR/$EXE_sprof $args

  # cat files together and move to working directory for letkf:
  `cat obs${ISLOT2}${MEM3}.tprof.dat obs${ISLOT2}${MEM3}.sprof.dat > $DSTDIR/obs${ISLOT2}${MEM3}.dat`

  # Remove individual files:
  rm -f obs${ISLOT2}${MEM3}.tprof.dat
  rm -f obs${ISLOT2}${MEM3}.sprof.dat

# echo "setup_obsop.sh:: temporarily ending early..."
# exit 1

done
