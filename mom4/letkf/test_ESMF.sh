#!/bin/sh --login
#PBS -r y                                                              #This job is restartable
#PBS -S /bin/sh                                                        #Do not change this - it keeps your job from issuing a false alarm
#PBS -E                                                                #Do not change this - it gives your job more and more useful Moab environment variables
# -- Request 120 cores
#PBS -l size=1
# 
# -- Specify a maximum wallclock of 4 hours
#PBS -l walltime=0:02:00
#
# -- Specify under which account a job should run
#PBS -A cpo_orr
#
# -- Set the name of the job, or moab will default to STDIN
#PBS -N cpo_ESMF_Regrid
#
# -- Set the queue: debug, batch, novel, bigmem, ldtn
#PBS -q debug
# 
# -- Set the partition (for Gaea)
##PBS -l partition=es
#PBS -l partition=c1:c2
#
# -- Set this as working directory
#PBS -d .
#
# -- Send output and error to same file
#PBS -j oe
#
#set -exv
set -x

source $MODULESHOME/init/sh 
module load esmf
module load cray-netcdf

src_file='basinmask_01.nc'
dst_file='grid_spec.nc'
#method='patch'
method='nearestdtos'
# [bilinear|patch|nearestdtos|neareststod|conserve]
src_type='GRIDSPEC'
dst_type='GRIDSPEC'
weight='regrid_basin_weights.nc'
lon_s='lon'
lat_s='lat'
lon_d='grid_x_T'
lat_d='grid_y_T'

echo "ESMF_BINDIR=$ESMF_BINDIR"

echo "TEST 1"
/sw/eslogin-c3/esmf/6.2.0/sles11.3_intel15.0.2/bin/ESMF_RegridWeightGen -s $src_file -d $dst_file -m $method -w $weight --src_type $src_type --dst_type $dst_type --src_coordinates $lon_s,$lat_s --dst_coordinates $lon_d,$lat_d

echo "TEST 2"
ESMF_RegridWeightGen -s $src_file -d $dst_file -m $method -w $weight --src_type $src_type --dst_type $dst_type --src_coordinates $lon_s,$lat_s --dst_coordinates $lon_d,$lat_d

echo "TEST 3"
$ESMF_BINDIR/ESMF_RegridWeightGen -s $src_file -d $dst_file -m $method -w $weight --src_type $src_type --dst_type $dst_type --src_coordinates $lon_s,$lat_s --dst_coordinates $lon_d,$lat_d

