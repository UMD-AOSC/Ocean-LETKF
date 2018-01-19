#!/bin/bash

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
    cp -d $DSTDIR/kmt.dat .

    # Specify the grid range
    SGXF=`expr $SGX + $GLOBAL_SG_NX - 1`
    SGYF=`expr $SGY + $GLOBAL_SG_NY - 1`

    # Fix the grid range in the input.nml file
    sed -i "s/^.*istart.*$/ istart = $SGX,/"  input.nml
    sed -i "s/^.*iend.*$/ iend = $SGXF,/"     input.nml
    sed -i "s/^.*jstart.*$/ jstart = $SGY,/"  input.nml
    sed -i "s/^.*jend.*$/ jend = $SGYF,/"     input.nml

    #---------------------------------------------------------------------------
    # Check if it is subgrid with ocean points
    #---------------------------------------------------------------------------

    # Open land sea mask and check for ocean points in this subgrid
    rm -f no_ocean.skip
    ${GLOBAL_LSCHECK_EXE} -f kmt.dat -is ${SGX} -ie ${SGXF} -js ${SGY} -je ${SGYF} -debug .false.

    if [ -f no_ocean.skip ]; then # If there are no ocean points
      echo "-------------------------------------------------------------"
      echo "Skipping subgrid $SG4. No ocean points in this region."
      echo "-------------------------------------------------------------"
    else
      echo "-------------------------------------------------------------"
      echo "Will process analysis for subgrid ${SG4} ..."
      echo "-------------------------------------------------------------"

      #---------------------------------------------------------------------------
      # Copy all links from parent directory and link to analysis files
      #---------------------------------------------------------------------------
      cp -d $DSTDIR/gs*.dat .
      cp -d $DSTDIR/gl??.dat .
      cp $DSTDIR/gridnl .
      cp $DSTDIR/oanl .
      ln -fs $DSTDIR/obs?????.dat .
      if [ "$GLOBAL_DO_WRITE_TILE" == "1" ]; then
        echo "letkf executable will write out this tile's analysis field to a new file."
        rm -f anal???.*.dat
      elif [ "$GLOBAL_DO_ANALYSIS_TEMPLATE" == "1" ]; then
        echo "Using global analysis template file as template for local tile."
        ln -fs $DSTDIR/anal*.dat .
      fi
    fi

    SGY=`expr $SGYF + 1`
  done
  SGX=`expr $SGXF + 1`
done

echo "Run complete. Exiting..."

exit 0
