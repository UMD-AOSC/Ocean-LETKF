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

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME
cd $DSTDIR

MEM_START=`expr $GLOBAL_MEM_START + 1`
MEM_END=`expr $GLOBAL_MEMBERS + 1`
NTILES=$GLOBAL_NTILES

SUBGRIDS_X=$(($GLOBAL_NX / GLOBAL_SG_NX))
SUBGRIDS_Y=$(($GLOBAL_NY / GLOBAL_SG_NY))
SUBGRIDS=$(($SUBGRIDS_X * SUBGRIDS_Y))
echo "Estimated total subgrids: $SUBGRIDS"
echo "Input total subgrids via params.sh: $GLOBAL_NTILES"

# Make the tile list:
tiledef_file='tiledef.txt'
rm -f $tiledef_file
for tile in [0-9][0-9][0-9][0-9]
do
  istart=$(awk -F'[ ,]' '/istart = / { print $4 } ' $tile/input.nml)
  iend=$(awk -F'[ ,]' '/iend = / { print $4 } ' $tile/input.nml)
  jstart=$(awk -F'[ ,]' '/jstart = / { print $4 } ' $tile/input.nml)
  jend=$(awk -F'[ ,]' '/jend = / { print $4 } ' $tile/input.nml)
  echo -e "${istart[$ij]} \t ${iend[$ij]} \t ${jstart[$ij]} \t ${jend[$ij]}" >> $tiledef_file
done

# Merge the output tiles into a single global analysis file for each ensemble member
for var in t s u v
do
  for ((mem=$MEM_START;mem<=$MEM_END;mem++))
  do
    MEM3=`printf %0.3d $mem`
    filetype="anal$MEM3.$var.dat"

    #---------------------
    # Create the filelist
    #---------------------
    filelist=`ls -1 ????/$filetype`
    filelist_file=filelist.txt #.$MEM3.$var.txt
    echo "This file list is:"
    rm -f $filelist_file
    for file in $filelist
    do
      echo $file >> $filelist_file
    done

    #---------------------
    # Create the grid specification file
    #---------------------

    # Run fortran program to read filelist and tilelist and output a global anaylsis
    bgfile="gs01$MEM3.$var.dat"
    infiles=$filelist_file
    outfile=anal${MEM3}.${var}.dat
#   echo "length of filelist = ${#filelist[@]}"
    $GLOBAL_MERGE_EXE -f $bgfile -flist $infiles -o $outfile -kmax 32 -inc .true. -ntiles 5 #$NTILES

    exit 1
  done
done

echo "Merge complete. Exiting..."

exit 0
