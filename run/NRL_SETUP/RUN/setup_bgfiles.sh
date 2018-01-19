#!/bin/bash
set -ex
date=/bin/date

source params.sh

# Set up the forecast files as background, and add template analysis files

MEM=$GLOBAL_MEM_START
MEMBERS=$GLOBAL_MEMBERS

nlon=$GLOBAL_NX
nlat=$GLOBAL_NY

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
DSTDIR=$SCRATCH/$EXPNAME

SRCDIR=$GLOBAL_FCST_SRCDIR
FILE_SUFFIX=$GLOBAL_RSRT_FILE_SUFFIX

YYYY=$GLOBAL_IYYYY
MM=$GLOBAL_IMM
DD=$GLOBAL_IDD
HH=$GLOBAL_IHH
ID1=$GLOBAL_ID1
ID2=$GLOBAL_ID2
ID3=$GLOBAL_ID3

TFILE=seatmp_pre_${ID1}_${ID2}_1o${nlon}x${nlat}_${YYYY}${MM}${DD}${HH}_${ID3}_fcstfld
SFILE=salint_pre_${ID1}_${ID2}_1o${nlon}x${nlat}_${YYYY}${MM}${DD}${HH}_${ID3}_fcstfld
UFILE=uucurr_pre_${ID1}_${ID2}_1o${nlon}x${nlat}_${YYYY}${MM}${DD}${HH}_${ID3}_fcstfld
VFILE=vvcurr_pre_${ID1}_${ID2}_1o${nlon}x${nlat}_${YYYY}${MM}${DD}${HH}_${ID3}_fcstfld
SSHFILE=seahgt_sfc_${ID1}_000000_1o${nlon}x${nlat}_${YYYY}${MM}${DD}${HH}_${ID3}_fcstfld

# Change to working directory
mkdir -p $DSTDIR
cd $DSTDIR

# Get NCODA namelists with grid definitions
cp $SRCDIR/${EXPNAME}000/analysism/gridnl .
cp $SRCDIR/${EXPNAME}000/analysism/oanl .

ISLOT=1
ISLOT2=`printf %0.2d $ISLOT`

while [ "$MEM" -le "$MEMBERS" ];
do
  MEM3=`printf %0.3d $MEM`
  MEMp1=`expr $MEM + 1`
  MEM3p1=`printf %0.3d $MEMp1`
  echo "Member: $MEM3"

  # Get the background restart files for the LETKF DA
  ln -fs $SRCDIR/$EXPNAME${MEM3}/restartm/$TFILE   gs${ISLOT2}${MEM3p1}.t${FILE_SUFFIX}
  ln -fs $SRCDIR/$EXPNAME${MEM3}/restartm/$SFILE   gs${ISLOT2}${MEM3p1}.s${FILE_SUFFIX}
  ln -fs $SRCDIR/$EXPNAME${MEM3}/restartm/$UFILE   gs${ISLOT2}${MEM3p1}.u${FILE_SUFFIX}
  ln -fs $SRCDIR/$EXPNAME${MEM3}/restartm/$VFILE   gs${ISLOT2}${MEM3p1}.v${FILE_SUFFIX}
  ln -fs $SRCDIR/$EXPNAME${MEM3}/restartm/$SSHFILE gs${ISLOT2}${MEM3p1}.ssh${FILE_SUFFIX}

  # Set up restart template files to be overwritten by the LETKF analysis
  if [ "$GLOBAL_DO_ANALYSIS_TEMPLATE" == "1"]; then
    for var in t s u v ssh
    do
      cp gs${ISLOT2}${MEM3p1}.${var}.dat anal${MEM3p1}.${var}.dat
    done
  fi

  MEM=$MEMp1
done
