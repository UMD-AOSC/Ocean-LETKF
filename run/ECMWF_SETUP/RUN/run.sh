#!/bin/bash

#PBS -N Ocean-LETKF_test
#PBS -q np
#PBS -l EC_total_tasks=288
#PBS -l EC_hyperthreads=2
#PBS -l EC_threads_per_task=1

#PBS -l EC_nodes=2
#PBS -l EC_predicted_walltime=600

#PBS -l EC_max_threads_per_node=72
#PBS -l EC_tasks_per_node=72

#PBS -l EC_job_tmpdir=DEFAULT
#PBS -l EC_threads_per_numa_node=36
#PBS -l EC_tmpdir_mem=0
#PBS -l EC_billing_account=ecrdasdm

##PBS -l select=1:vntype=cray_login:EC_accept_from_queue=np:ncpus=0:mem=300MB+1:vntype=cray_compute:EC_accept_from_queue=np:mem=120GB

set -ex

source params.sh

RUNDIR=$GLOBAL_RUNDIR
name=$GLOBAL_LETKF_NAME
EXEDIR=$GLOBAL_LETKF_ROOT/build/build_letkf/$name.build
EXEFILE=letkf.$name.x

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME/WORK
cd $DSTDIR

# Copy necessary files to here:
cp $RUNDIR/input.nml .    # input namelist
cp $EXEDIR/$EXEFILE .     # executable file

# Ensure that analysis template files match background files
MEMBERS=$GLOBAL_MEMBERS
FILE_SUFFIX=$RESTART_FILE_SUFFIX
ISLOT=1
ISLOT2=`printf "%02d" $ISLOT`
MEM=1
while [ $MEM -le $MEMBERS ]; do
  MEM2=`printf "%02d" $MEM`
  MEM3=`printf "%03d" $MEM`

  # Copy forecast files to here:
  bfile=gs${ISLOT2}${MEM3}${FILE_SUFFIX}
  afile=anal${MEM3}${FILE_SUFFIX}

  echo "Creating file: $afile"
  cp -f $bfile $DSTDIR/$afile &

  MEM=`expr $MEM + 1`
done
time wait

# Run Ocean-LETKF DA system:
echo "Running Ocean-LETKF..."
aprun -N $EC_tasks_per_node -n $EC_total_tasks -j $EC_hyperthreads ./$EXEFILE
echo "Run complete. Exiting..."

exit 0
