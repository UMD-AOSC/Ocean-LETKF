#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "Model preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model_prep/${MEMBERID}
workdir2=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
mkdir -p ${workdir2}
cd ${workdir}

echo "This is the model preparation step for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > model_prep.out

#STEVE: active code:
USE_EFLX=1
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the date for the next analysis cycle
date=/bin/date
ainc=4
ainc_units=days
AY=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%Y`
AM=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%m`
AD=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%d`
AH=`$date -d "$IY-$IM-$ID $ainc $ainc_units" +%H`
AN=$IN
AS=$IS

# set up the ENSEMBLE surface forcing fluxes
if [ $USE_EFLX ]; then
  cp $FLXDIR/mkDlyPSBC4 ${workdir}
  cp $FLXDIR/mPS.sh ${workdir}
  echo "$IY$IM$ID" > run_date
  /usr/bin/perl -i.bak -pe "s/L_YR=\d+/L_YR=$IY/" ./mPS.sh
  ./mPS.sh $IY
fi

mkdir -p ${workdir2}/INPUT
mkdir -p ${workdir2}/RESTART

#STEVE: copy INPUT/RESTART files to here
rsrtOUT=$OUTPUT/rsrt/$IY$IM$ID
cp $OUTPUT/rsrt/INPUT/* ${workdir2}/INPUT   #STEVE: temporary fix...
cp $rsrtOUT/$MEM3/* ${workdir2}/INPUT
#   cp $OUTPUT/gues/$MEM3/$TY$TM$TD.$TH$TN$TS.* ${workdir2}
cp $rsrtOUT/$MEM3/input.nml ${workdir2}
#cp $OUTPUT/rsrt/$MEM3/data_table ${workdir2}
cp $OUTPUT/rsrt/INPUT/data_table.orig ${workdir2}/INPUT/data_table  #STEVE: the RA2 files are causing the model to crash
cp $OUTPUT/rsrt/INPUT/data_table.orig ${workdir2}/data_table        #STEVE: the RA2 files are causing the model to crash
cp $rsrtOUT/$MEM3/diag_table ${workdir2}
cp $rsrtOUT/$MEM3/field_table ${workdir2}

# Need tripolar grid_spec.nc (for sure)
cp $NATURE/grid_spec_om3_core3.nc ${workdir2}/INPUT/grid_spec.nc

# SFC Fluxes
if [ $USE_EFLX ]; then
  for file in `ls -d ${workdir}/RA2_daily_*.nc`; do
    ln -fs $file ${workdir2}/INPUT
  done
fi
#STEVE: now copy over new files from analysis to be used as background
#STEVE: This is for LETKF-BASIC and a partial requirement for LETKF-RIP.
#STEVE: For LETKF-IAU, instead of this, we will copy a correctors file.
#NOTE: sometimes, if the script halted, these files have already been moved, so just move on in that case

exit 0
#STEVE: not doing IAU for now

if [ $USE_IAU ]; then
  #STEVE: Use letkf-calculated means to calculate the analysis increment:
  #STEVE: NOTE - should be compared at FIVE-DAY-MARK - is that what this is??? check that letkf calc's mean at analysis time
  #STEVE: should change these dates to match the analysis time $AY$AM$AD$AH.grd
  ln -fs $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$INN$ISS.ocean_temp_salt.res.nc anal.ocean_temp_salt.res.nc
  ln -fs $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$INN$ISS.ocean_velocity.res.nc  anal.ocean_velocity.res.nc
  ln -fs $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$INN$ISS.ocean_sbc.res.nc       anal.ocean_sbc.res.nc
  ln -fs $OUTPUT/gues/$MEM3/$AY$AM$AD.$AH$AN$AS.ocean_temp_salt.res.nc   gues.ocean_temp_salt.res.nc
  ln -fs $OUTPUT/gues/$MEM3/$AY$AM$AD.$AH$AN$AS.ocean_velocity.res.nc    gues.ocean_velocity.res.nc
  ln -fs $OUTPUT/gues/$MEM3/$AY$AM$AD.$AH$AN$AS.ocean_sbc.res.nc         gues.ocean_sbc.res.nc

  # Get grid_spec.nc file for kmt data
  if [ ! -f grid_spec.nc ]; then
    ln -fs $NATURE/INPUT/grid_spec.nc .
  fi

  # Copy template netCDF files:
  #cp -f $NATURE/TEMPLATE/*_increment.nc .

  # Generate increment files to be used by mom4p1
  $MOM4/letkf/grd2cor.x -$inc_units $inc -afile anal -bfile gues
  #STEVE: temp_increment.nc file should now exist:
  if [ -f "temp_increment.nc" ]; then
    echo "temp_increment.nc file generated for $MEM3"
  else
    echo "error, no temp_increment.nc file was generated"
    exit 1
  fi
  if [ -f "salt_increment.nc" ]; then
    echo "salt_increment.nc file generated for $MEM3"
  else
    echo "error, no salt_increment.nc file was generated"
    exit 1
  fi
  if [ -f "u_increment.nc" ]; then
    echo "u_increment.nc file generated for $MEM3"
  else
    echo "error, no u_increment.nc file was generated"
    exit 1
  fi
  if [ -f "v_increment.nc" ]; then
    echo "v_increment.nc file generated for $MEM3"
  else
    echo "error, no v_increment.nc file was generated"
    exit 1
  fi
  if [ -f "eta_increment.nc" ]; then
    echo "eta_increment.nc file generated for $MEM3"
  else
    echo "error, no eta_increment.nc file was generated"
    exit 1
  fi

  rm anal.ocean_*.res.nc
  rm gues.ocean_*.res.nc
  # ID correctors at time of the beginning of the IAU run
  mv temp_increment.nc $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.temp_increment.nc
  mv salt_increment.nc $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.salt_increment.nc
  mv u_increment.nc    $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.u_increment.nc
  mv v_increment.nc    $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.v_increment.nc
  mv eta_increment.nc  $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.eta_increment.nc

  echo "ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.temp_increment.nc ${workdir2}/INPUT/temp_increment.nc"
  ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.temp_increment.nc       ${workdir2}/INPUT/temp_increment.nc

  echo "ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.salt_increment.nc ${workdir2}/INPUT/salt_increment.nc"
  ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.salt_increment.nc       ${workdir2}/INPUT/salt_increment.nc

  echo "ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.u_increment.nc    ${workdir2}/INPUT/u_increment.nc"
  ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.u_increment.nc          ${workdir2}/INPUT/u_increment.nc

  echo "ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.v_increment.nc    ${workdir2}/INPUT/v_increment.nc"
  ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.v_increment.nc          ${workdir2}/INPUT/v_increment.nc

  echo "ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.eta_increment.nc  ${workdir2}/INPUT/eta_increment.nc"
  ln -fs $OUTPUT/increment/$MEM3/$IY$IM$ID$IH.eta_increment.nc        ${workdir2}/INPUT/eta_increment.nc
fi
 

exit
