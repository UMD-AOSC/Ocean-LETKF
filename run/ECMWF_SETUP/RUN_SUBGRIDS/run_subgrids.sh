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

CDIR=/home/rd/dasp/perm/DATA/RUN_SUBGRIDS
source $CDIR/params.sh

RUNDIR=$GLOBAL_RUNDIR
name=$GLOBAL_LETKF_NAME
EXEDIR=$GLOBAL_LETKF_ROOT/build/build_letkf/$name.build
EXEFILE=letkf.$name.x

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME/WORK
cd $DSTDIR

MEMBERS=$GLOBAL_MEMBERS
FILE_SUFFIX=$GLOBAL_LETKF_FILE_SUFFIX
ISLOT=1
ISLOT2=`printf "%02d" $ISLOT`

SUBGRIDS=$GLOBAL_num_subgrids
SG=$GLOBAL_subgrid_start
while [ $SG -le $SUBGRIDS ]; do
  SG4=`printf "%04d" $SG`

  # Change to working directory specific to this subgrid
  cd $DSTDIR/$SG4
   
# # Ensure that analysis template files match background files
# MEM=1
# while [ $MEM -le $MEMBERS ]; do
#   MEM2=`printf "%02d" $MEM`
#   MEM3=`printf "%03d" $MEM`

#   # Copy forecast files to here:
#   bfile=gs${ISLOT2}${MEM3}${FILE_SUFFIX}
#   afile=anal${MEM3}${FILE_SUFFIX}

#   echo "Creating file: $afile"
#   cp -f $bfile $DSTDIR/$afile &

#   MEM=`expr $MEM + 1`
# done
# time wait

  # Run Ocean-LETKF DA system:
  if [ -f no_ocean.skip ]; then
    echo "-------------------------------------------------------------"
    echo "Skipping subgrid $SG4. No ocean points in this region."
    echo "-------------------------------------------------------------"
    continue
  else
    # Copy necessary files to here:
    cp $RUNDIR/input.nml .    # input namelist
    cp $EXEDIR/$EXEFILE .     # executable file

    echo "-------------------------------------------------------------"
    echo "Running Ocean-LETKF for subgrid $SG4..."
    echo "..."
    aprun -N $EC_tasks_per_node -n $EC_total_tasks -j $EC_hyperthreads ./$EXEFILE
    echo "..."
    echo "Run complete for subgrid $SG4. Continuing..."
    echo "-------------------------------------------------------------"
  fi

done

exit 0
