#!/bin/bash
set -ex
date=/bin/date

source params.sh

# Set up the observation/innovation .dat files by running the letkf obsop_tprof.x / obsop_sprof.x routines.

MEM=$GLOBAL_MEM_START
MEMBERS=$GLOBAL_MEMBERS

IYYYY=$GLOBAL_IYYYY
IMM=$GLOBAL_IMM
IDD=$GLOBAL_IDD

inc=$GLOBAL_inc
inc_units=$GLOBAL_inc_units

ISLOT=1  # Use if there are multiple observation time bins
ISLOT2=`printf "%02d" $ISLOT`

EXEDIR=/home/rd/dasp/perm/DATA/RUN

SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
QC_REF_FILE=$GLOBAL_QC_REF

OBSDIR=$SCRATCH/$EXPNAME/DATA/OBS

name=ECMWF_nemo.nemovarHx0
EXE_tprof=obsop.$name.tprof.x
EXE_sprof=obsop.$name.sprof.x

FILE_PREFIX=$OBS_FILE_PREFIX #0001_nrt_5_2_
#FILE_SUFFIX=_profbqc_01_fdbk.nc
FILE_SUFFIX=$OBS_FILE_SUFFIX #_profb_01_fdbk.nc

OBS_FILELIST='observation_infile_list.txt'
DO_QC=$GLOBAL_DO_QC
if [ $DO_QC == "1" ]; then
  # Use a feedback qc file to filter out NEMOVAR qc'd obs
  # then reorder obs in all members so they align properly by index
  ALIGN_OBS_EXE="reduce_obs.py $QC_REF_FILE"
else
  # Use a regular feedback file prior to NEMOVAR qc
  # Reorder observations in all members so they align properly by index
  ALIGN_OBS_EXE=align_obs.py
fi

DSTDIR=$SCRATCH/$EXPNAME/WORK
mkdir -p $DSTDIR
cd $DSTDIR
rm -f $DSTDIR/$OBS_FILELIST
cp -f $EXEDIR/$ALIGN_OBS_EXE .

#-------------------------------------------------------------------------------
# Get original observations
#-------------------------------------------------------------------------------
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

  MEM2=`printf "%02d" $MEM`
  echo "Preparing to sort observations for member: $MEM2"

  FILE=${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
  obsinfile=$OBSDIR/$MEM2/$FILE

  echo $obsinfile >> $OBS_FILELIST

  MEM=`expr $MEM + 1`
done

#-------------------------------------------------------------------------------
# Run python script to align observations:
#-------------------------------------------------------------------------------
echo "Sorting observations (python script) ..."
python $ALIGN_OBS_EXE

#-------------------------------------------------------------------------------
# Replace the control member observations to all members:
#-------------------------------------------------------------------------------
SOURCE_FILE=sorted_00_${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
MEM=`expr $GLOBAL_MEM_START + 1`
MEM=`echo $(($MEM>0?$MEM:1))`  # i.e. start at 1 unless the desired starting member is > 0 (not the control)

# Get original observations
YYYY=$IYYYY
MM=$IMM
DD=$IDD
while [ $MEM -lt $MEMBERS ]; do
  
  MEM2=`printf "%02d" $MEM`
  echo "Replacing observations with control for member: $MEM2"

  DEST_FILE=sorted_${MEM2}_${FILE_PREFIX}${YYYY}${MM}${DD}_${NY}${NM}${ND}${FILE_SUFFIX}
  ncks -A -v POTM_OBS,PSAL_OBS,LONGITUDE,LATITUDE $SOURCE_FILE $DEST_FILE &

  echo $obsinfile >> $OBS_FILELIST

  MEM=`expr $MEM + 1`
done

time wait

# Next, run the obsop_XXX.x on each member ...
