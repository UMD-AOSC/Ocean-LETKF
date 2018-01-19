#!/bin/bash
set -ex

source params.sh

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME
#---------------------
# Copy the executable to destination directory and change directories
#---------------------
cp mt.x $DSTDIR/
cd $DSTDIR

MEM_START=`expr $GLOBAL_MEM_START + 1`
MEM_END=`expr $GLOBAL_MEMBERS + 1`
NTILES=$GLOBAL_NTILES

# Get initial start time
start0=`date +%s`

#-------------------------------------------------------------------------------
# Read in all tiles and compute a new global analysis for each member
#-------------------------------------------------------------------------------
start=`date +%s`

#---------------------
# get the tile dimensions
#---------------------

#istart=($(awk -F'[ ,]' '/istart = / { print $4 } ' ????/input.nml))
#iend=($(awk -F'[ ,]' '/iend = / { print $4 } ' ????/input.nml))
#jstart=($(awk -F'[ ,]' '/jstart = / { print $4 } ' ????/input.nml))
#jend=($(awk -F'[ ,]' '/jend = / { print $4 } ' ????/input.nml))

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

#---------------------
# Create the input_mt.nml file
#---------------------

#------------------------------------------
# Cycle through variables and ensemble members:
#------------------------------------------
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
    infile=$filelist_file
    outfile=anal${MEM3}.${var}.dat
    ./mt.x -f $infile -o $outfile -ntiles $GLOBAL_NTILES

    exit 1
  done
done

time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (guess) :: $runtime / $runtime0"

#-------------------------------------------------------------------------------
# Read global analyses and forecasts and compute global increments for each member
#-------------------------------------------------------------------------------
start=`date +%s`

time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (guess) :: $runtime / $runtime0"

#-------------------------------------------------------------------------------
# Compute forecast ensemble mean and spread
#-------------------------------------------------------------------------------
start=`date +%s`


time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (guess) :: $runtime / $runtime0"

#-------------------------------------------------------------------------------
# Compute analysis ensemble mean and spread
#-------------------------------------------------------------------------------
start=`date +%s`


time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (analysis) :: $runtime / $runtime0"

#-------------------------------------------------------------------------------
# Compute mean analysis increment
#-------------------------------------------------------------------------------
start=`date +%s`


time wait
end=`date +%s`
runtime0=$((end-start0))
runtime=$((end-start))
echo "runtime (analysis) :: $runtime / $runtime0"

