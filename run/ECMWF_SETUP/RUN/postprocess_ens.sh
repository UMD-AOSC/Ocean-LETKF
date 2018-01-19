#!/bin/bash
set -ex

module load cdo

DSTDIR=/home/rd/dasp/scratch/WORK
cd $DSTDIR

# Get initial start time
start0=`date +%s`

# Compute forecast ensemble mean and spread
start=`date +%s`
cdo ensavg gs01???.restart.nc bmean.restart.nc &
cdo ensstd gs01???.restart.nc bsprd.restart.nc &
time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (guess) :: $runtime / $runtime0"

# Compute analysis ensemble mean and spread
start=`date +%s`
cdo ensavg anal???.restart.nc amean.restart.nc &
cdo ensstd anal???.restart.nc asprd.restart.nc &
time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (analysis) :: $runtime / $runtime0"

# Sample timing:
#
#+ cdo ensavg gs01001.restart.nc gs01002.restart.nc gs01003.restart.nc gs01004.restart.nc gs01005.restart.nc bmean.restart.nc
#cdo ensavg: Processed 2215784410 values from 25 variables over 5 timesteps ( 110.81s )
#+ cdo ensstd gs01001.restart.nc gs01002.restart.nc gs01003.restart.nc gs01004.restart.nc gs01005.restart.nc bsprd.restart.nc
#cdo ensstd: Processed 2215784410 values from 25 variables over 5 timesteps ( 182.50s )
#++ date +%s
#+ end=1510230506
#+ runtime0=555
#+ runtime=555
#+ echo 'runtime (guess) :: 555 / 555'
#runtime (guess) :: 555 / 555
#++ date +%s
#+ start=1510230506
#+ cdo ensavg anal001.restart.nc anal002.restart.nc anal003.restart.nc anal004.restart.nc anal005.restart.nc amean.restart.nc
#cdo ensavg: Processed 2215784410 values from 25 variables over 5 timesteps ( 117.63s )
#+ cdo ensstd anal001.restart.nc anal002.restart.nc anal003.restart.nc anal004.restart.nc anal005.restart.nc asprd.restart.nc
#cdo ensstd: Processed 2215784410 values from 25 variables over 5 timesteps ( 122.81s )
#++ date +%s
#+ end=1510230909
#+ runtime0=958
#+ runtime=403
#+ echo 'runtime (analysis) :: 403 / 958'
#runtime (analysis) :: 403 / 958

