#!/bin/bash
date=/bin/date
set -ex

function get_data {

  # Copy the observation/innovation data file:
  echo "getting file from ec: ..."
  if [ ! -f ${FILE} ]; then
    ecp ec:${SRCDIR}/${FILE} .
  fi

  # Decompress the observation/innovation data file:
  echo "Unzipping file ..."
  if [ -f ${FILE} ]; then
    gunzip ${FILE}
  fi

  # Untar specific files needed 
  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
  tar -xf $FILE_tar $FILE 
  mv $FILE $DSTDIR/

  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_qc}
  tar -xf $FILE_tar $FILE 
  mv $FILE $DSTDIR/

  rm -f $FILE_tar

}

# Setup observation innovations computed from NEMOVAR for use in LETKF

MEM=0
MEMBERS=50

IYYYY=2005
IMM=10
IDD=01
IHH=18
INN=00
ISS=00

EYYYY=2005
EMM=10
EDD=01
EHH=18
ENN=00
ESS=00

inc=1
inc_units='days'

# els ec:/ERAS/cera20c/3253/an/fdbk/

#SCRATCH=/scratch/rd/dasp
SCRATCH=/lus/snx11064/scratch/rd/dasp
SRCDIR0=/ERAS/cera20c/3253/an/fdbk
EXPNAME=CERA20C_50mem

TMPDIR=$SCRATCH/$EXPNAME/DATA_ORIG/OBS
DSTDIR0=$SCRATCH/$EXPNAME/DATA/OBS
FILE_PREFIX=3253_
FILE_SUFFIX_gz=_fdbk.tar.gz
FILE_SUFFIX_tar=_fdbk.tar
FILE_SUFFIX=_profb_01_fdbk.nc
FILE_SUFFIX_qc=_profbqc_01_fdbk.nc
YYYY=$IYYYY
MM=$IMM
DD=$IDD
HH=$IHH
NN=$INN
SS=$ISS

while [ $MEM -le $MEMBERS ]; do

  NY=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%Y`
  NM=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%m`
  ND=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%d`

  echo "For reference, the current date is:"
  echo "$YYYY-$MM-$DD"
  echo "The next date is:"
  echo "$NY-$NM-$ND"

  MEM1=`printf "%01d" $MEM`
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
  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_gz}
  FILE_tar=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_tar}

  # Get NEMO feedback file for observation and model equivalent data:
  get_data &

  MEM=`expr $MEM + 1`

done

time wait

# YYYY=$NY
# MM=$NM
# DD=$ND
