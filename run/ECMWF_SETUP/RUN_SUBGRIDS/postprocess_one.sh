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
cdo -O ensavg gs01???.restart.nc bmean.restart.nc &
cdo -O ensstd gs01???.restart.nc bsprd.restart.nc &

# Compute analysis ensemble mean and spread
cdo -O ensavg anal???.restart.nc amean.restart.nc &
cdo -O ensstd anal???.restart.nc asprd.restart.nc &

time wait

cdo -O sub amean.restart.nc bmean.restart.nc ainc.restart.nc

end=`date +%s`
runtime0=$((end-start0))
echo "runtime :: $runtime0"
