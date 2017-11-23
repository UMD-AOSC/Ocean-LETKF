#!/bin/bash
date=/bin/date
set -ex

IYYYY=2017
IMM=01
IDD=01
IHH=00

inc=5
inc_units='days'

SCRATCH=/lus/snx11064/scratch/rd/dasp
#SRCDIR=/fws4/lb/work/rd/dipb/gu5w/$IYYYY$IMM$IDD$IHH #2010010100
SRCDIR=/fws10/lb/work/rd/dihz/gucq/$IYYYY$IMM$IDD$IHH #2010010100
EXPNAME=TEST_5memOceanSubgrids

NY=`$date -d "$IYYYY-$IMM-$IDD $inc $inc_units" +%Y`
NM=`$date -d "$IYYYY-$IMM-$IDD $inc $inc_units" +%m`
ND=`$date -d "$IYYYY-$IMM-$IDD $inc $inc_units" +%d`

DSTDIR=$SCRATCH/$EXPNAME/DATA/FCST
#FILE_PREFIX=assim_background_state_DI_
#e.g. gucq_20170106_000000_restart_0479.nc
FILE_PREFIX=gucq_${NY}${NM}${ND}_000000_restart_
FILE_SUFFIX=.nc

MEM=0
MEMBERS=5 # exclusive

while test $MEM -lt $MEMBERS
do
  MEM2=`printf "%02d" $MEM`
  mkdir -p $DSTDIR/$MEM2

  FILES=$SRCDIR/opa${MEM}/restarts/${FILE_PREFIX}????${FILE_SUFFIX}

  echo "Copying files: $FILES"
  echo "Copying to: $DSTDIR/$MEM2/"
  cp -f $FILES $DSTDIR/$MEM2/
  
  MEM=`expr $MEM + 1`
done
