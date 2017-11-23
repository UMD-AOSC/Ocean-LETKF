#!/bin/bash
set -ex

# Get initial start time
start0=`date +%s`

module load cdo

source params.sh

# Run sample subgrid:
SG=200
SG4=`printf "%04d" $SG`

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME/WORK/$SG4
cd $DSTDIR

# Compute forecast ensemble mean and spread
cdo ensavg gs01???.restart.nc bmean.restart.nc &
cdo ensstd gs01???.restart.nc bsprd.restart.nc &

# Compute analysis ensemble mean and spread
cdo ensavg anal???.restart.nc amean.restart.nc &
cdo ensstd anal???.restart.nc asprd.restart.nc &

time wait
end=`date +%s`
runtime0=$((end-start0))
echo "runtime :: $runtime0"
