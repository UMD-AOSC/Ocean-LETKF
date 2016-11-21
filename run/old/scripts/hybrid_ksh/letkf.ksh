#!/bin/ksh --login
set -e

#module load intel
#module load netcdf/4.1.3-intel
#module load nco

echo "LETKF run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir}
cd ${workdir}
export PATH=$PATH:`pwd` #STEVE: need to access nco executables copied to here, since modules are not available.

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
# This is for assimilating the surface forcing fields
#cp $LDIR/coeff_m2s.nc coeff_m2s.nc
#cp $LDIR/coeff_s2m.nc coeff_s2m.nc
cp $REGRID/coeff_m2s.nc coeff_m2s.nc
cp $REGRID/coeff_s2m.nc coeff_s2m.nc

# If an input namelist exists, then copy it here to be read in at runtime:
if [ -f $LDIR/input.nml ]; then
  cp $LDIR/input.nml .
fi
if [ -f $INPUT_INIT/../$altimetry_climatology_file ]; then
  ln -f $INPUT_INIT/../$altimetry_climatology_file .
fi

# START RUN LETKF ################################################
cp $LDIR/$LETKFexe .
echo "Running LETKF..."
aprun -n $PBS_NP $LETKFexe
echo "This is the LETKF run for cycle ${YYYYMMDDHH}" > ${workdir}/letkf.out
# END   RUN LETKF ################################################

# DO_SFCFLUXES
#On the next analysis cycle in the model_prep stage, add the perturbations to the SFC fields for each member.
if [ $DO_SFCFLUXES -eq "1" ]; then
  cp /sw/xe6/nco/4.0.8/sles11.1_netcdf4.2.0_gnu4.7.0/bin/ncatted .
  for file in `ls -d ${workdir}/SFA_*.nc`; do
    ncatted -O -a units,TIME,c,c,"days since $AY-$AM-$AD 12:00:00" $file
  done
fi

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

exit 0

