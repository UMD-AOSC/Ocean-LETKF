#!/bin/bash

MEM=0
MEMBERS=5
IYYYY=2016
IMM=12
IDD=31 #27

#SCRATCH=/scratch/rd/dasp
SCRATCH=/lus/snx11064/scratch/rd/dasp/
SRCDIR0=/emos/OCEA/5/0001
EXPNAME=OCEA5_5mem

TMPDIR=$SCRATCH/$EXPNAME/DATA_ORIG/FCST
DSTDIR0=$SCRATCH/$EXPNAME/DATA/FCST
mkdir -p $DSTDIR0
FILE_PREFIX=0001_
FILE_SUFFIX=_000000_restart.nc
YYYY=$IYYYY
MM=$IMM
DD=$IDD
while [ $MEM -lt $MEMBERS ]; do

  MEM1=`printf "%01d" $MEM`
  MEM2=`printf "%02d" $MEM`

  echo "creating temporary directory: $TMPDIR/$MEM2"
  mkdir -p $TMPDIR/$MEM2
  cd $TMPDIR/$MEM2
  DSTDIR=$DSTDIR0/$MEM2
  mkdir -p $DSTDIR

  SRCDIR=$SRCDIR0/opa${MEM1}/restart/${YYYY}
  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}${FILE_SUFFIX}

  # Copy the data file:
  #ecp ec:/emos/OCEA/5/0001/opa3/restart/2016/0001_nrt_20161228_000000_restart.nc.gz .
  echo "getting file from ec: ..."
  if [ ! -f ${FILE}.gz ]; then
    ecp ec:${SRCDIR}/${FILE}.gz .
  fi

  # Decompress the data file:
  echo "Unzipping file ..."
  if [ -f ${FILE}.gz ]; then
    rm -f ${FILE}
    gunzip ${FILE}.gz
  fi

  # Get the grid file info:
# echo "Extracting grid information..."
# ncks -a -v nav_lon,nav_lat,nav_lev ${FILE} $DSTDIR/grid_spec.nc

  # Get the state vector data:
  #echo "Extracting b-type variables ..."
  #ncks -a -v tb,sb,ub,vb,sshb ${FILE} $DSTDIR/bvars.nc
  echo "Extracting n-type variables ..."
  ncks -a -v tn,sn,un,vn,sshn ${FILE} $DSTDIR/nvars.nc

  echo "Removing file ${FILE} ..."
  if [ -f ${FILE} ]; then
    rm ${FILE}
  fi

  MEM=`expr $MEM + 1`

done
