#!/bin/bash
date=/bin/date
set -ex

IYYYY=2017
IMM=01
IDD=01
IHH=00

SCRATCH=/lus/snx11064/scratch/rd/dasp
#SRCDIR=/fws4/lb/work/rd/dipb/gu5w/$IYYYY$IMM$IDD$IHH #2010010100
SRCDIR=/fws10/lb/work/rd/dihz/gucq/$IYYYY$IMM$IDD$IHH #2010010100
EXPNAME=TEST_5memOceanSubgrids

DSTDIR=$SCRATCH/$EXPNAME/DATA/OBS
FILE_PREFIX=profb
FILE_SUFFIX=_01_fdbk_00.nc

MEM=0
MEMBERS=5

while test $MEM -lt $MEMBERS
do
  MEM2=`printf "%02d" $MEM`
  mkdir -p $DSTDIR/$MEM2

  FILES=$SRCDIR/opa${MEM}/${FILE_PREFIX}*${FILE_SUFFIX}

  echo "Copying files: $FILES"
  echo "Copying to: $DSTDIR/$MEM2/"
  cp -f $FILES $DSTDIR/$MEM2/ &
  
  MEM=`expr $MEM + 1`
done

time wait
