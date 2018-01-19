#!/bin/bash
set -ex

# Date and directory settings:
MEM=0
MEMBERS=50
IYYYY=2005
IMM=10
IDD=01
IHH=18
INN=00
ISS=00

#SCRATCH=/scratch/rd/dasp
SCRATCH=/lus/snx11064/scratch/rd/dasp
SRCDIR0=/ERAS/cera20c/3253/an/restart
EXPNAME=CERA20C_50mem

TMPDIR=$SCRATCH/$EXPNAME/DATA_ORIG/FCST
DSTDIR0=$SCRATCH/$EXPNAME/DATA/FCST
FILE_PREFIX=3253_
FILE_SUFFIX=_${IHH}${INN}${ISS}_restart_2.nc

# Get files:
mkdir -p $DSTDIR0
YYYY=$IYYYY
MM=$IMM
DD=$IDD
HH=$IHH
NN=$INN
SS=$ISS
while [ $MEM -le $MEMBERS ]; do

  MEM2=`printf "%02d" $MEM`

  echo "creating temporary directory: $TMPDIR/$MEM2"
  mkdir -p $TMPDIR/$MEM2
  cd $TMPDIR/$MEM2
  DSTDIR=$DSTDIR0/$MEM2
  mkdir -p $DSTDIR

  # Specify source file format:
  if [ "$MEM" -eq "0" ]; then
    SRCDIR=$SRCDIR0/control/${YYYY}
  else
    SRCDIR=$SRCDIR0/enda${MEM2}/${YYYY}
  fi
  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}${FILE_SUFFIX}

  # Copy the data file:
  #ecp ec:/ERAS/cera20c/3253/an/restart/enda01/2005/3253_20051001_180000_restart_2.nc.gz
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

  # Get the state vector data:
  echo "Extracting n-type variables ..."
  ncks -a -v tn,sn,un,vn,sshn ${FILE} $DSTDIR/nvars.nc

  echo "Removing file ${FILE} ..."
  if [ -f ${FILE} ]; then
    rm ${FILE}
  fi

  MEM=`expr $MEM + 1`

done
