#!/bin/bash
set -ex
date=/bin/date

source params.sh

# Link the observations to the subgrid working directories

MEM=$GLOBAL_MEM_START
MEMBERS=$GLOBAL_MEMBERS
SUBGRIDS=$GLOBAL_num_subgrids

ISLOT=1  # Use if there are multiple observation time bins

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME

SRCDIR=$SCRATCH/$EXPNAME/WORK/obs
DSTDIR0=$SCRATCH/$EXPNAME/WORK
mkdir -p $DSTDIR0


# Create a subdirectory in the working directory for each subgrid analysis
SG=$GLOBAL_subgrid_start
while [ $SG -lt $SUBGRIDS ]; do
  SG4=`printf "%04d" $SG`
  echo "Linking to subgrid number $SG out of $SUBGRIDS..."

  DSTDIR=$DSTDIR0/$SG4
  mkdir -p $DSTDIR
  cd $DSTDIR

  # Link observation files to here:
  echo "Creating file: $bfile"
  ln -f $SRCDIR/obs*.dat .

  SG=`expr $SG + 1`
done

