#!/bin/ksh --login
set -ex

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

# (This is needed to read and write mom4p1 netcdf files in letkf)
ln -f $INPUT_INIT/grid_spec.nc grid_spec.nc

# If an input namelist exists, then copy it here to be read in at runtime:
if [ -f $INPUT_INIT/../letkf/input.nml ]; then
  cp $INPUT_INIT/../letkf/input.nml .
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

rm -f grid_spec.nc

#echo "letkf.ksh :: Forcing (intentional) ERROR on exit for further experimentation..."
#exit 1 #STEVE: force error on exit to enable repeating experiment with different versions of letkf
exit 0
