#!/bin/bash
set -ex
date=/bin/date

source params.sh

# Set up the forecast files as background, and add template analysis files
MEMBERS=$GLOBAL_MEMBERS
SUBGRIDS=$GLOBAL_num_subgrids

ISLOT=1  # Use if there are multiple observation time bins
FILE_SUFFIX=$GLOBAL_LETKF_FILE_SUFFIX

#SCRATCH=/scratch/rd/dasp
SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME

SRCDIR=$SCRATCH/$EXPNAME/DATA/FCST
DSTDIR0=$SCRATCH/$EXPNAME/WORK
mkdir -p $DSTDIR0

# Create a subdirectory in the working directory for each subgrid analysis
SG=$GLOBAL_subgrid_start
while [ $SG -lt $SUBGRIDS ]; do
  SG4=`printf "%04d" $SG`

  DSTDIR=$DSTDIR0/$SG4
  mkdir -p $DSTDIR
  echo "Processing all $SUBGRIDS subgrids..."

  MEM=$GLOBAL_MEM_START
  while [ $MEM -lt $MEMBERS ]; do
  
    MEM2=`printf "%02d" $MEM`
    MEMp1=`expr $MEM + 1`
    MEM3=`printf "%03d" $MEMp1`
    ISLOT2=`printf "%02d" $ISLOT`

    echo "Processing member $MEM2 of subgrid $SG4..."

    # Copy forecast files to here:
    bfile=gs${ISLOT2}${MEM3}${FILE_SUFFIX}
    afile=anal${MEM3}${FILE_SUFFIX}

    echo "Creating file: $bfile"
    cp -f $SRCDIR/$MEM2/${GLOBAL_RESTART_FILE_PREFIX}${SG4}${GLOBAL_RESTART_FILE_SUFFIX} $DSTDIR/$bfile &

    echo "Creating file: $afile"
    cp -f $SRCDIR/$MEM2/${GLOBAL_RESTART_FILE_PREFIX}${SG4}${GLOBAL_RESTART_FILE_SUFFIX} $DSTDIR/$afile &

    MEM=`expr $MEM + 1`
  done

  time wait

  SG=`expr $SG + 1`
done

