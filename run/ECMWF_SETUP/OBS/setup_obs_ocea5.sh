#!/bin/bash
date=/bin/date

# Setup observation innovations computed from NEMOVAR for use in LETKF

MEM=0
MEMBERS=5

IYYYY=2016
IMM=12
IDD=26

EYYYY=2017
EMM=01
EDD=01

inc=6
inc_units='days'

# ec:/emos/OCEA/5/0001/opa0/fdbk/2016/0001_nrt_5_2_20161226_20170101_fdbk.tar.gz

TMPDIR=/scratch/rd/dasp
SRCDIR0=/emos/OCEA/5/0001
DSTDIR0=/home/rd/dasp/perm/DATA/OBS
FILE_PREFIX=0001_nrt_5_2_
FILE_SUFFIX_gz=_fdbk.tar.gz
FILE_SUFFIX_tar=_fdbk.tar
FILE_SUFFIX=_profb_01_fdbk.nc
FILE_SUFFIX_qc=_profbqc_01_fdbk.nc
YYYY=$IYYYY
MM=$IMM
DD=$IDD

while [ $MEM -lt $MEMBERS ]; do

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

  SRCDIR=$SRCDIR0/opa${MEM1}/fdbk/${YYYY}
  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_gz}

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
  FILE_tar=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_tar}

  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
  tar -xf $FILE_tar $FILE 
  mv $FILE $DSTDIR/

  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX_qc}
  tar -xf $FILE_tar $FILE 
  mv $FILE $DSTDIR/

  rm -f $FILE_tar

  MEM=`expr $MEM + 1`

done

# YYYY=$NY
# MM=$NM
# DD=$ND
