#!/bin/bash

#PBS -N OcnLETKF_OneSubgrid
#PBS -q np
#PBS -l EC_total_tasks=72
#PBS -l EC_hyperthreads=2
#PBS -l EC_threads_per_task=1

#PBS -l EC_nodes=1
##PBS -l EC_predicted_walltime=600

#PBS -l EC_max_threads_per_node=72
#PBS -l EC_tasks_per_node=72

#PBS -l EC_job_tmpdir=DEFAULT
#PBS -l EC_threads_per_numa_node=36
#PBS -l EC_tmpdir_mem=0
#PBS -l EC_billing_account=ecrdasdm

##PBS -l select=1:vntype=cray_login:EC_accept_from_queue=np:ncpus=0:mem=300MB+1:vntype=cray_compute:EC_accept_from_queue=np:mem=120GB

set -ex

# Run sample subgrid:
SG=200
SG4=`printf "%04d" $SG`

CDIR=/home/rd/dasp/perm/DATA/RUN_SUBGRIDS
source $CDIR/params.sh

RUNDIR=$GLOBAL_RUNDIR
name=$GLOBAL_LETKF_NAME
EXEDIR=$GLOBAL_LETKF_ROOT/build/build_letkf/$name.build
EXEFILE=letkf.$name.x

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME/WORK/$SG4
cd $DSTDIR

# Copy necessary files to here:
cp $RUNDIR/input.nml .    # input namelist
cp $EXEDIR/$EXEFILE .     # executable file

# Run Ocean-LETKF DA system:
echo "Running Ocean-LETKF..."
aprun -N $EC_tasks_per_node -n $EC_total_tasks -j $EC_hyperthreads ./$EXEFILE
echo "Run complete. Exiting..."

