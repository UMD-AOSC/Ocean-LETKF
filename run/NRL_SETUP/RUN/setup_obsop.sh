#!/bin/bash
set -ex
date=/bin/date

source params.sh

function run_obsop {

# local EXEDIR=$1
# local EXE_tprof=$2
# local EXE_sprof=$3
# local ISLOT2=$4
# local MEM3=$5

  # Process T/S profiles
  obsinfile=MVO_prp.dat
  ensinfile=MVO_ENS_obs.dat
  args="-obsin $obsinfile -ensin $ensinfile -obsout obs${ISLOT2}${MEM3}.${type1}.dat -mem $MEM"
  $EXEDIR/$EXE_3d $args

  # Process sea surface height data
  if [ "$DO_SSH" == "1" ]; then
    obsinfile=SSH_prp.dat
    ensinfile=SSH_ENS_obs.dat
    args="-obsin $obsinfile -ensin $ensinfile -obsout obs${ISLOT2}${MEM3}.${type2}.dat -mem $MEM -superob .true."
    $EXEDIR/$EXE_ssh $args
  else
    rm -f obs${ISLOT2}${MEM3}.${type2}.dat
    touch obs${ISLOT2}${MEM3}.${type2}.dat
  fi

  # cat files together and move to working directory for letkf:
  `cat obs${ISLOT2}${MEM3}.${type1}.dat obs${ISLOT2}${MEM3}.${type2}.dat > $DSTDIR/obs${ISLOT2}${MEM3}.dat`

  # Remove individual files:
  rm -f obs${ISLOT2}${MEM3}.${type1}.dat
  rm -f obs${ISLOT2}${MEM3}.${type2}.dat

}

# Set up the observation/innovation .dat files by running the letkf obsop_*.x routines

MEM=$GLOBAL_MEM_START
MEMBERS=$GLOBAL_MEMBERS

#IYYYY=$GLOBAL_IYYYY
#IMM=$GLOBAL_IMM
#IDD=$GLOBAL_IDD
#IHH=$GLOBAL_IHH

# Step forward number of forecast hours to get observations at the appropriate datetime
IYYYY=`$date -d "$GLOBAL_IYYYY-$GLOBAL_IMM-$GLOBAL_IDD $GLOBAL_IHH $GLOBAL_FCST_HRS4 hours" +%Y`
IMM=`$date -d "$GLOBAL_IYYYY-$GLOBAL_IMM-$GLOBAL_IDD $GLOBAL_IHH $GLOBAL_FCST_HRS4 hours" +%m`
IDD=`$date -d "$GLOBAL_IYYYY-$GLOBAL_IMM-$GLOBAL_IDD $GLOBAL_IHH $GLOBAL_FCST_HRS4 hours" +%d`
IHH=`$date -d "$GLOBAL_IYYYY-$GLOBAL_IMM-$GLOBAL_IDD $GLOBAL_IHH $GLOBAL_FCST_HRS4 hours" +%H`

#echo "$IYYYY-$IMM-$IDD $IHH"
#exit 1

inc=$GLOBAL_inc
inc_units=$GLOBAL_inc_units

SRCDIR=$GLOBAL_OBS_SRCDIR

FILE_PREFIX=$GLOBAL_OBS_FILE_PREFIX #0001_nrt_5_2_
FILE_SUFFIX=$GLOBAL_OBS_FILE_SUFFIX #_profb_01_fdbk.nc

ISLOT=1  # Use if there are multiple observation time bins
ISLOT2=`printf "%02d" $ISLOT`

SCRATCH=$GLOBAL_SCRATCH #/lus/snx11064/scratch/rd/dasp
EXPNAME=$GLOBAL_EXPNAME #CERA20C_50mem

LETKF_ROOT=$GLOBAL_LETKF_ROOT
name=$GLOBAL_LETKF_NAME
EXEDIR=$LETKF_ROOT/build/build_obsop/${name}.build
type1=ncoda_3d
type2=ncoda_ssh
EXE_3d=obsop.$name.${type1}.x
EXE_ssh=obsop.$name.${type2}.x

GRDDIR=$GLOBAL_GRDDIR
GRDFILE=$GLOBAL_GRDFILE

DSTDIR=$SCRATCH/$EXPNAME
mkdir -p $DSTDIR
ln -fs $GRDDIR/restartm/$GRDFILE_LON $DSTDIR/glon.dat
ln -fs $GRDDIR/restartm/$GRDFILE_LAT $DSTDIR/glat.dat
ln -fs $GRDDIR/restartm/$GRDFILE_LEV $DSTDIR/glev.dat
ln -fs $GRDDIR/restartm/$GRDFILE_KMT $DSTDIR/kmt.dat
cd $DSTDIR

YYYY=$IYYYY
MM=$IMM
DD=$IDD
HH=$IHH
while [ $MEM -le $MEMBERS ]; do
  
  NY=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%Y`
  NM=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%m`
  ND=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%d`

  echo "For reference, the current date is:"
  echo "$YYYY-$MM-$DD :: $HH"
  echo "The next date is:"
  echo "$NY-$NM-$ND :: $HH"

  MEM=`expr $MEM + 1`
  MEM3=`printf "%03d" $MEM`

  # Prepare the T/S profiles
  file1=${GLOBAL_OBS_FILE_PREFIX}MVO_prp${GLOBAL_OBS_FILE_SUFFIX}${YYYY}${MM}${DD}${HH}
  file2=${GLOBAL_OBS_FILE_PREFIX}MVO_ENS_obs${GLOBAL_OBS_FILE_SUFFIX}${YYYY}${MM}${DD}${HH}
  ln -fs $SRCDIR/$file1 MVO_prp.dat
  ln -fs $SRCDIR/$file2 MVO_ENS_obs.dat

  # Prepare the satellite altimetry
  if [ "$DO_SSH"=="1" ]; then
    file1=${GLOBAL_OBS_FILE_PREFIX}SSH_prp${GLOBAL_OBS_FILE_SUFFIX}${YYYY}${MM}${DD}${HH}
    file2=${GLOBAL_OBS_FILE_PREFIX}SSH_ENS_obs${GLOBAL_OBS_FILE_SUFFIX}${YYYY}${MM}${DD}${HH}
    ln -fs $SRCDIR/$file1 SSH_prp.dat
    ln -fs $SRCDIR/$file2 SSH_ENS_obs.dat
  fi

  run_obsop # use locally defined function (above) to permit parallel execution

# echo "setup_obsop.sh:: temporarily ending early..."
# exit 1

done

time wait
