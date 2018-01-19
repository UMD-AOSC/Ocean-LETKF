#!/bin/ksh
#set -e
date=/bin/date

# Calling each member independently...
echo "I am member $MEMBERID"
echo "INPUTS:"
echo "days := $days"
echo "YYYYMMDDHH := $YYYYMMDDHH"
echo "USE_ALTIMETRY := $USE_ALTIMETRY"
echo "altimetry_climatology_file := $altimetry_climatology_file"
echo "EXP_DATA := $EXP_DATA"
echo "UTIL := $UTIL"
echo "INIT := $INIT"
echo "INPUT_INIT := $INPUT_INIT"
#STEVE: override so only profiles are computed:
OBSDIR1=/lustre/f1/unswept/Steve.Penny/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep
OBSDIR2=$OBSDIR1
echo "MEMBERS := $MEMBERS"
LDIR=/autofs/na1_home1/Steve.Penny/letkf/mom4/letkf_GAEA
OBSOPexe=obsope.${MEMBERS}
echo "OBSOPexe := $OBSOPexe"
echo "OBSDIR1 := $OBSDIR1"
echo "OBSDIR5 := $OBSDIR5"
echo "ATIME := $ATIME"
echo "LDIR := $LDIR"
echo "NSLOTS := $NSLOTS"

################################################################################
# Setup times
################################################################################
inc=5
inc_units=days
ainc=$inc
YMDH=$YYYYMMDDHH
YYYY=${YMDH:0:4}
MM=${YMDH:4:2}
DD=${YMDH:6:2}
HH=${YMDH:8:2}
NN=00
SS=00
NY=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%Y`
NM=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%m`
ND=`$date -d "$YYYY-$MM-$DD $inc $inc_units" +%d`
NH=`$date -d "$YYYY-$MM-$DD $HH:$NN:$SS $inc $inc_units" +%H`
NYMD=$NY$NM$ND
echo "NYMD := $NYMD"
AY=`$date -d "$YYYY-$MM-$DD $ainc $inc_units" +%Y`
AM=`$date -d "$YYYY-$MM-$DD $ainc $inc_units" +%m`
AD=`$date -d "$YYYY-$MM-$DD $ainc $inc_units" +%d`
AH=`$date -d "$YYYY-$MM-$DD $HH:$NN:$SS $ainc $inc_units" +%H`
AYMD=$AY$AM$AD
echo "AYMD := $AYMD"

CUTDATE=2003070700 #1999020100
GDIR=$INPUT_INIT

################################################################################
# Set up working directory
################################################################################
MEM=$MEMBERID
PDIR=${EXP_DATA}/$YYYYMMDDHH/post
mkdir -p $PDIR
MEM2=`printf %.2d ${MEM}`
MEM3=`printf %.3d ${MEM}`
echo "MEM2 := $MEM2"
echo "MEM3 := $MEM3"
letkf=${EXP_DATA}/$YYYYMMDDHH/letkf
hybrid=${EXP_DATA}/$YYYYMMDDHH/hybrid
model=${EXP_DATA}/$YYYYMMDDHH/model

################################################################################
# Compute the oma increments for the hybrid
################################################################################
workdir=${PDIR}/oma_hyb/${MEM2}
mkdir -p $workdir
cd $workdir

OUTPUT=$EXP_DATA/STORE/anal/oma_hyb_ts
mkdir -p $OUTPUT/${MEM2}

ln -fs $GDIR/grid_spec.nc .
ln -fs $INIT/$altimetry_climatology_file .
cp $LDIR/$OBSOPexe .

if [ ${AYMD}00 -lt $CUTDATE ]; then
  OBSDIR=$OBSDIR1
else
  OBSDIR=$OBSDIR2
fi
echo "Processing obs: ${OBSDIR}/$AYMD.dat to obsin.dat"
ln -f $OBSDIR/$AYMD.dat obsin.dat

ln -fs ${hybrid}/anal$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
ln -fs ${letkf}/anal$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
ln -fs ${letkf}/anal$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
# $OBSOPexe -obsin $OBSDIR/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
./$OBSOPexe
mv obsout.dat $OUTPUT/${MEM2}/$AYMD.dat
rm -f grid_spec.nc
rm -f $altimetry_climatology_file

################################################################################
# Compute the oma increments for letkf
################################################################################
workdir=${PDIR}/oma_prehyb/${MEM2}
mkdir -p $workdir
cd $workdir

OUTPUT=$EXP_DATA/STORE/anal/oma_prehyb_ts
mkdir -p $OUTPUT/${MEM2}

ln -fs $GDIR/grid_spec.nc .
ln -fs $INIT/$altimetry_climatology_file .
cp $LDIR/$OBSOPexe .

if [ ${AYMD}00 -lt $CUTDATE ]; then
  OBSDIR=$OBSDIR1
else
  OBSDIR=$OBSDIR2
fi
echo "Processing obs: ${OBSDIR}/$AYMD.dat to obsin.dat"
ln -f $OBSDIR/$AYMD.dat obsin.dat

ln -fs ${letkf}/anal$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
ln -fs ${letkf}/anal$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
ln -fs ${letkf}/anal$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
# $OBSOPexe -obsin $OBSDIR1/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
./$OBSOPexe
mv obsout.dat $OUTPUT/${MEM2}/$AYMD.dat
rm -f grid_spec.nc
rm -f $altimetry_climatology_file

################################################################################
# Compute the omf increments for only T and S profiles, if using altimetry
################################################################################
workdir=${PDIR}/omf/${MEM2}
mkdir -p $workdir
cd $workdir

OUTPUT=$EXP_DATA/STORE/fcst/omf_ts
mkdir -p $OUTPUT/${MEM2}

ln -fs $GDIR/grid_spec.nc .
ln -fs $INIT/$altimetry_climatology_file .
cp $LDIR/$OBSOPexe .

if [ ${SYMD}00 -le $CUTDATE ]; then
  OBSDIR=$OBSDIR1
else
  OBSDIR=$OBSDIR2
fi

ISLOT=1
while test $ISLOT -le $NSLOTS
do
  ISLOT2=`printf %.2d ${ISLOT}`
  SY=`$date -d "$YYYY-$MM-$DD $ISLOT $inc_units" +%Y`
  SM=`$date -d "$YYYY-$MM-$DD $ISLOT $inc_units" +%m`
  SD=`$date -d "$YYYY-$MM-$DD $ISLOT $inc_units" +%d`
  SYMD=$SY$SM$SD
  echo "================================================="
  echo "Observation date: $SYMD"
  echo "================================================="
  echo "Processing obs: ${OBSDIR}/$SYMD.dat to obsin.dat"
  ln -f $OBSDIR/$SYMD.dat obsin.dat

  ln -fs ${letkf}/gs$ISLOT2$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
  ln -fs ${letkf}/gs$ISLOT2$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
  ln -fs ${letkf}/gs$ISLOT2$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
  # $OBSOPexe -obsin $OBSDIR/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
  ./$OBSOPexe
  mv obsout.dat $OUTPUT/${MEM2}/$SYMD.dat
  echo "***********************************"
  echo "mv obsout.dat $OUTPUT/${MEM2}/$SYMD.dat"
  ls -l $OUTPUT/${MEM2}/$SYMD.dat
  echo "***********************************"
  rm -f obsin.dat

  ISLOT=`expr $ISLOT + 1`
done

rm -f grid_spec.nc
rm -f $altimetry_climatology_file

################################################################################
# Erase extra files
# Delete most files if this is not the last of the month.
################################################################################
TDIR=${EXP_DATA}/$YYYY$MM$DD$HH
ERASE_DIR=$TDIR/model
if [ $MM -eq $NM ]; then
  echo "Running date: $YMDH, deleting..." #, analysis time: $AYMD"
  rm -rf $TDIR/model_prep/$MEM2
  rm -f $TDIR/model/${MEM2}/INPUT/*
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_bih_friction.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_con_temp.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_density.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_frazil.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_neutralB.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_neutral.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_residency.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_sigma_transport.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_solo.intermediate.res
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_solo.res
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_thickness.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_velocity_advection.res.nc
  rm -f $TDIR/model/${MEM2}/RESTART/ocean_tracer.res
  rm -rf $TDIR/letkf_prep/$MEM2
  rm -f $TDIR/letkf/gs??${MEM3}.ocean_*.nc
  rm -f $TDIR/letkf/grid_spec.nc
  rm -f $TDIR/letkf/aEtaCds9399.nc
  rm -f $TDIR/letkf/coeff_*.nc
  rm -f $TDIR/g4p1/tmpZ.????
  rm -f $TDIR/g4p1/INPUT/grid_spec.nc
  rm -rf $TDIR/g4p1/INPUT/work_*
  rm -f $TDIR/g4p1/INPUT/temp_sfc_restore.nc
  rm -f $TDIR/g4p1/INPUT/salt_sfc_restore.nc
  rm -f $TDIR/g4p1/INPUT/aEtaCds9399.i3e
  rm -f $TDIR/g4p1/INPUT/salt12.nc
else
  echo "Keeping date: $YMDH, skipping delete..."
  echo "(keeping last date of every month)"
  echo "(but clearing gs\*.nc files from $TDIR/letkf directory)"
  rm -f $TDIR/letkf/gs??${MEM3}.ocean_*.nc
  rm -f $TDIR/letkf/grid_spec.nc
  rm -f $TDIR/letkf/aEtaCds9399.nc
  rm -f $TDIR/letkf/coeff_*.nc
fi

################################################################################
# Erase extra files and gzip the important netcdf data files
################################################################################

# Before zipping, we need to delete the model INPUT directories of the next timestep
model_next=${EXP_DATA}/$NY$NM$ND$NH/model/$MEM2
rm -f $model_next/INPUT/*

# Zip up all of the analysis files
# NOTE: these can only be zipped after the next timestep INPUT files are removed
echo "gzip letkf analysis files..."
gzip $letkf/anal${MEM3}.ocean_*.res.nc
echo "done zipping letkf analysis files."
echo "gzip hybrid analysis files..."
gzip $hybrid/anal${MEM3}.ocean_*.res.nc
echo "done zipping hybrid analysis files."

# Zip up all of the model/RESTART files
echo "gzip model RESTART/ocean_*.nc files..."
gzip $model/$MEM2/RESTART/ocean_*.nc
echo "done zipping model RESTART files."

# Remove the daily restart files
echo "Removing daily model RESTART files..."
rm -f $TDIR/model/${MEM2}/RESTART/${YYYY}*
rm -f $TDIR/model/${MEM2}/RESTART/${NY}*
echo "done removing model RESTART files."

# Before zipping, we need to delete the post processing directories that are linked
echo "Remove post processing directories and contents for member $MEM2 ..."
rm -rf $PDIR/*/$MEM2

################################################################################
# Erase extra files and gzip the last year if on a new year for this cycle
# Check if this is a new year
################################################################################

exit 0
