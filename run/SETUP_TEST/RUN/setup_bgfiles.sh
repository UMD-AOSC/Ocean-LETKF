#!/bin/bash
date=/bin/date

# Set up the forecast files as background, and add template analysis files

MEM=0
MEMBERS=5

ISLOT=1  # Use if there are multiple observation time bins

SRCDIR=/home/rd/dasp/perm/DATA/FCST
FILE_SUFFIX=.restart.nc
DSTDIR=/home/rd/dasp/scratch/WORK
mkdir -p $DSTDIR

while [ $MEM -lt $MEMBERS ]; do
  
  MEM2=`printf "%02d" $MEM`
  ISLOT2=`printf "%02d" $ISLOT`

  # Copy forecast files to here:
  MEM=`expr $MEM + 1`
  MEM3=`printf "%03d" $MEM`
  bfile=gs${ISLOT2}${MEM3}${FILE_SUFFIX}
  afile=anal${MEM3}${FILE_SUFFIX}
  echo "Creating file: $bfile"
  cp -f $SRCDIR/$MEM2/nvars.nc $DSTDIR/$bfile
  echo "Creating file: $afile"
  cp -f $SRCDIR/$MEM2/nvars.nc $DSTDIR/$afile

done
