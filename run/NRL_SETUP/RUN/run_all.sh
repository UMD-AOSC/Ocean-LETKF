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
DSTDIR=$SCRATCH/$EXPNAME
cd $DSTDIR

SUBGRIDS_X=$(($GLOBAL_NX / GLOBAL_SG_NX))
SUBGRIDS_Y=$(($GLOBAL_NY / GLOBAL_SG_NY))
SUBGRIDS=$(($SUBGRIDS_X * SUBGRIDS_Y))
echo "Estimated total subgrids: $SUBGRIDS"
#exit 1

IDX=0
SG=0
SGX=1
while [ $SGX -le $GLOBAL_NX ]; do
  SGY=1
  while [ $SGY -le $GLOBAL_NY ]; do
    SG=`expr $SG + 1`
    SG4=`printf "%04d" $SG`
    
    #---------------------------------------------------------------------------
    # Change to working directory specific to this subgrid
    #---------------------------------------------------------------------------
    mkdir -p $DSTDIR/$SG4
    cd $DSTDIR/$SG4

    #---------------------------------------------------------------------------
    # Copy necessary files to here:
    #---------------------------------------------------------------------------
    cp $RUNDIR/input.nml .    # input namelist
    cp $EXEDIR/$EXEFILE .     # executable file

    # Specify the grid range
    SGXF=`expr $SGX + $GLOBAL_SG_NX - 1`
    SGYF=`expr $SGY + $GLOBAL_SG_NY - 1`

    # Keep a separate file for the record to verify the input.nml subgrid definitions
    echo "X0 = $SGX, Y0 = $SGY" > subgrid_def.txt
    echo "Xf = $SGXF, Yf = $SGYF" >> subgrid_def.txt

    # Fix the grid range in the input.nml file
    sed -i "s/^.*istart.*$/ istart = $SGX,/"  input.nml
    sed -i "s/^.*iend.*$/ iend = $SGXF,/"     input.nml
    sed -i "s/^.*jstart.*$/ jstart = $SGY,/"  input.nml
    sed -i "s/^.*jend.*$/ jend = $SGYF,/"     input.nml

#   echo "Exiting early on purpose (1)..."
#   exit 1

    #---------------------------------------------------------------------------
    # Check if it is subgrid with ocean points
    #---------------------------------------------------------------------------
    if [ -f no_ocean.skip ]; then # If there are no ocean points
      echo "-------------------------------------------------------------"
      echo "Skipping subgrid $SG4. No ocean points in this region."
      echo "-------------------------------------------------------------"
    else
      echo "-------------------------------------------------------------"
      echo "Processing LETKF analysis for subgrid ${SG4} ..."
      echo "-------------------------------------------------------------"
      IDX=`expr $IDX + 1`

      #---------------------------------------------------------------------------
      # Copy all links from parent directory and link to analysis files
      #---------------------------------------------------------------------------
      cp -d $DSTDIR/gs*.dat .
      cp -d $DSTDIR/gl??.dat .
      cp -d $DSTDIR/kmt.dat .
      cp $DSTDIR/gridnl .
      cp $DSTDIR/oanl .
      ln -fs $DSTDIR/obs?????.dat .
      if [ "$GLOBAL_DO_WRITE_TILE" == "1" ]; then
        echo "letkf executable will write out this tile's analysis field to a new file."
        rm -f anal???.*.dat
      else
        ln -fs $DSTDIR/anal*.dat .
      fi

      #---------------------------------------------------------------------------
      # Run Ocean-LETKF DA system:
      #---------------------------------------------------------------------------
      echo "-------------------------------------------------------------"
      echo "Running Ocean-LETKF for subgrid $SG4..."
      echo "..."
      #aprun -N $EC_tasks_per_node -n $EC_total_tasks -j $EC_hyperthreads ./$EXEFILE
      mpirun -n 4 ./$EXEFILE
#     mpirun -n 4 ./$EXEFILE &
      echo "..."
      echo "Run complete for subgrid $SG4. Continuing..."
      echo "-------------------------------------------------------------"

#     echo "Exiting early on purpose (1)..."
#     exit 1

      # Wait for every 4 LETKF runs to finish
#     n=`expr $IDX % 4`
#     if [ "$n" == "0" ]; then
#       time wait 
#     fi
    fi

#   echo "Exiting early on purpose (2)..."
#   exit 2
#   read -n 1 -p "Press any key to continue:" dummy

    SGY=`expr $SGYF + 1`
  done
  SGX=`expr $SGXF + 1`
done

echo "Run complete. Exiting..."

exit 0
