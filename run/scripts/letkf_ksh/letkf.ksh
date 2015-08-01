#!/bin/ksh --login
#===============================================================================
# SCRIPT:
# letkf_prepc.ksh
#
# PURPOSE:
# This script prepares background and observation data for letkf 
# This 'compact' version runs all processing for each member 
# (as opposed to a separate instance not each member & timeslot)
#
# MODULES USED:
#  (e.g. on Gaea)
#  module swap PrgEnv-pgi PrgEnv-intel
#  module load netcdf
#
# INPUTS:
#  YYYYMMDDHH    :: string containing 4-digit year, 2-digit month, 2-digit day, 2-digit hour
#  MEMBERID      :: Ensemble member number
#  EXP_DATA      :: directory containing experiment output data
#  NSLOTS        :: number of timeslots to use for 4D-LETKF (e.g. "5" for 5 days)
#  days          :: forecast length (in integer days)
#  OBSOPexe      :: execuatable for LETKF observation operator
#  OBSDIR1       :: Observation directory, primary
#  OBSDIR5       :: Observation directory to use for analysis time only
#  INPUT_INIT    :: Directory containing static model input files
#  LDIR          :: directory of letkf executable and obsoperator executable
#  USE_ALTIMETRY :: flag to assimilate altimetry data (1==true,0==false)
#  altimetry_climatology_file       :: model eta climatology, for assimilating AVISO altimetry  
# 
#===============================================================================
# Author      :: Stephen G. Penny
# Institution :: University of Maryland (UMD) 
#                Department of Atmospheric and Oceanic Science (AOSC), and
#                National Centers for Environmental Prediction (NCEP)
#                National Oceanograpic and Atmospheric Administration (NOAA)
# Email       :: Steve.Penny@noaa.gov
#===============================================================================

set -e
set -x
set -v

echo "LETKF run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir}
cd ${workdir}

#STEVE: active code:
USE_INFLADJ=0
USE_ADAPOBS=0
#DO_SFCFLUXES=1   #This is input via the xml script
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the date for the next analysis cycle
date=/bin/date
inc=5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

# Update the analysis time and date
ainc=$inc #`expr $inc - 1`
ainc_units=days
AY=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%Y`
AM=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%m`
AD=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%d`
AH=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%H`
AN=$IN
AS=$IS

# Update the date for the previous analysis cycle
pinc=$inc
pinc_units='days ago'
PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
PN=$IN
PS=$IS
workdir0=${EXP_DATA}/$PY$PM$PD$PH/letkf

### Processs inputs

# Copy the inflation (copy instead of link, because letkf rewrites over the file)
if [ "$USE_INFLADJ" -eq "1" -a -f $workdir0/infl_redux.grd ]; then
  ln -f $workdir0/infl_redux.grd infl_mul.grd  
else
# if test -f $workdir0/infl_mul.grd
  if test -f $workdir0/infl_out.grd
  then
#   cp $workdir0/infl_mul.grd infl_mul.grd
    cp $workdir0/infl_out.grd infl_mul.grd
  else
    echo "WARNING: There is no inflation file from the previous timestep. If this is not the first analysis cycle, there is a possible ERROR."
  fi
fi

# Also, if there is adaptive obervation error, copy that file as well...
if test -f $workdir0/adapt_oer.grd
then
  ln -f $workdir0/adapt_oer.grd adapt_inp.grd
fi

# (This is needed to read and write mom4p1 netcdf files in letkf)
ln -f $INPUT_INIT/grid_spec.nc grid_spec.nc
if [ "$USE_ALTIMETRY" -eq "1" ]; then
  if [ -f $altimetry_climatology_file ]; then
    ln -f $altimetry_climatology_file aEtaCds9399.nc
  else
    echo "The aEtaCds9399.nc is not present: $altimetry_climatology_file"
  fi
else
  echo "NOT using ALTIMETRY"
  echo "USE_ALTIMETRY = $USE_ALTIMETRY"
  echo "altimetry_climatology_file = $altimetry_climatology_file"
fi

# START RUN LETKF ################################################
cp $LDIR/$LETKFexe .
echo "Running LETKF..."
aprun -n $PBS_NP $LETKFexe
echo "This is the LETKF run for cycle ${YYYYMMDDHH}" > ${workdir}/letkf.out
# END   RUN LETKF ################################################

### process outputs
#if test -f infl_mul.grd
if test -f infl_out.grd
then
  if [ "$USE_INFLADJ" -eq "1" -a ! -f infl_redux.grd ]; then  #STEVE: if inflation relaxation is needed, it is applied here
    echo "Using inflation adjust..."
    cp $LDIR/$INFLadj .
    aprun -n 1 $INFLadj
  else
    echo "Not using inflation adjust..."
  fi
else
  echo "WARNING: no inflation output file infl_out.grd from letkf."
fi

#STEVE: If being used, this should be made another job so that it can be run in parallel with the model
#       prior to the next analysis cycle with letkf
#
# Get the files for adaptive obs error
if test -f adap_obserr.grd
then
  if [ "$USE_ADAPOBS" -eq "1" ]; then
    echo "Using adaptive observations..."
    cp $LDIR/$ADAPTexe .
    aprun -n $PBS_NP $ADAPTexe
  else
    echo "Using prescribed observation error. (Not adaptive observation error)"
  fi
fi

#STEVE: avoid issue with too many hard links being created:
rm -f $workdir/aEtaCds9399.nc
rm -f $workdir/grid_spec.nc

exit 0

