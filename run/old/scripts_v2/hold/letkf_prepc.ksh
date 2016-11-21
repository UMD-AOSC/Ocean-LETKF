#!/bin/ksh --login
# This script prepares background and observation data for letkf 
# Tripolar background grid is converted to spherical if necessary
#
# This 'compact' version runs all processing for each members (not each member & timeslot)
#
# Written by Dr. Stephen G. Penny
#
set -e
#module load mpt
#module load intel
#module load netcdf
#module load nco

echo "LETKF preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir_fcst=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/RESTART
workdir_inpt=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/INPUT
workdir2=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir2}

MEM3=`printf %.3d ${MEMBERID}`

#STEVE: active code:
#DO_SFCFLUXES=1  #This is input via the xml script
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

FCST=0
while test $FCST -lt $NSLOTS
do
  FCST=`expr $FCST + 1`
  ISLOT2=`printf %.2d ${FCST}`
  echo "Posting forecast ${ISLOT2}"
  workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf_prep/${MEMBERID}/${ISLOT2}
  mkdir -p ${workdir}
  cd ${workdir}
  echo "This is the LETKF preparation step for member ${MEMBERID}, for cycle ${YYYYMMDDHH}" > $workdir/letkf_prep_${MEMBERID}.out

  # Update the ISLOT date
  date=/bin/date
  #sinc=`expr $FCST - 1` #STEVE: use this if I want the dates offset
  sinc=1
  sinc_units=days
  NY=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%Y`
  NM=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%m`
  ND=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%d`
  NH=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%H`
  IY=$NY
  IM=$NM
  ID=$ND
  IH=$NH
  echo "processing cycle as: $IY$IM$ID$IH"

  # FIRST, link the background model data

  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       gs${ISLOT2}$MEM3.ocean_sbc.res.nc

  #STEVE: add 'fill value' to netcdf files
  cp /sw/xe6/nco/4.0.8/sles11.1_netcdf4.2.0_gnu4.7.0/bin/ncatted .
  ncatted -O -a _FillValue,temp,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ncatted -O -a _FillValue,salt,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ncatted -O -a _FillValue,u,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ncatted -O -a _FillValue,v,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ncatted -O -a _FillValue,sea_lev,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_sbc.res.nc

  # FIRST (a), Link background files to the letkf working directory:
  #STEVE: this should NOT be necessary with new external obs operator for letkf:
# for file in `ls -d ${workdir}/gs${ISLOT2}$MEM3.*`; do
#   echo "linking $file to ${workdir2}..."
#   ln -f $file ${workdir2}/
# done

  #STEVE: need a template for the output analysis files:
  #       MUST COPY, NOT LINK
  if [ "${ISLOT2}" -eq "${ATIME}" ]; then #STEVE: trying to only do it once per member...
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/anal$MEM3.ocean_temp_salt.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/anal$MEM3.ocean_velocity.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/anal$MEM3.ocean_sbc.res.nc
    #STEVE: may want to zero these out for peace of mind...
  fi

  # SECOND (b), link the observations
  #Use different observation source collection for different slots.
  # (i.e. only use surface obs at time of analysis, use profiles through analysis cycle window)
  echo "For ISLOT2=$ISLOT2, and ATIME=${ATIME}"
  if [ "$ISLOT2" -eq "${ATIME}" ]; then
    OBSDIR=$OBSDIR5
  else
    OBSDIR=$OBSDIR1
  fi
  echo "Processing obs: ${OBSDIR}/$IY$IM$ID$IH.dat to obs${ISLOT2}${MEM3}.dat"
  ln -f $INPUT_INIT/grid_spec.nc .
  cp $LDIR/$OBSOPexe .
  ln -f $OBSDIR/$IY$IM$ID.dat obsin.dat
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
# $OBSOPexe -obsin $OBSDIR/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
  $OBSOPexe 
  ln -f obsout.dat ${workdir2}/obs${ISLOT2}${MEM3}.dat
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc
# ls ${workdir2}

#ln -s $OBSDIR/$IYYYY$IMM$IDD$IHH.dat obsin.dat
#ln -s $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd gs01${MEM3}.grd
#ln -fs $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd gues.grd
#./$OBSOPE -obsin $OBSDIR/$IYYYY$IMM$IDD$IHH.dat -gues $OUTPUT/gues/${MEM3}/$IYYYY$IMM$IDD$IHH.grd -obsout obs${ISLOT2}${MEM3}.dat > obsope.log
#mv obsout.dat obs01${MEM3}.dat

done #DONE FCST loop

#STEVE: DO_SFCFLUXES
if [ "${DO_SFCFLUXES}" -eq "1" ]; then
  #STEVE: for now, this is hard-coded into LETKF. Later, change this to a namelist option
  #       so that it can be turned off if desired to shorten runtime.
    ln -f ${workdir_inpt}/RA2_daily_dlw.nc    ${workdir2}/SFC_${MEM3}_daily_dlw.nc
    ln -f ${workdir_inpt}/RA2_daily_dsw.nc    ${workdir2}/SFC_${MEM3}_daily_dsw.nc
    ln -f ${workdir_inpt}/RA2_daily_LONGWV.nc ${workdir2}/SFC_${MEM3}_daily_LONGWV.nc
    ln -f ${workdir_inpt}/RA2_daily_PRATE.nc  ${workdir2}/SFC_${MEM3}_daily_PRATE.nc
    ln -f ${workdir_inpt}/RA2_daily_pres.nc   ${workdir2}/SFC_${MEM3}_daily_pres.nc
    ln -f ${workdir_inpt}/RA2_daily_q2m.nc    ${workdir2}/SFC_${MEM3}_daily_q2m.nc
    ln -f ${workdir_inpt}/RA2_daily_QFLUX.nc  ${workdir2}/SFC_${MEM3}_daily_QFLUX.nc
    ln -f ${workdir_inpt}/RA2_daily_SHRTWV.nc ${workdir2}/SFC_${MEM3}_daily_SHRTWV.nc
    ln -f ${workdir_inpt}/RA2_daily_t2m.nc    ${workdir2}/SFC_${MEM3}_daily_t2m.nc
    ln -f ${workdir_inpt}/RA2_daily_TAUX.nc   ${workdir2}/SFC_${MEM3}_daily_TAUX.nc
    ln -f ${workdir_inpt}/RA2_daily_TAUY.nc   ${workdir2}/SFC_${MEM3}_daily_TAUY.nc
    ln -f ${workdir_inpt}/RA2_daily_TFLUX.nc  ${workdir2}/SFC_${MEM3}_daily_TFLUX.nc
    ln -f ${workdir_inpt}/RA2_daily_U10.nc    ${workdir2}/SFC_${MEM3}_daily_U10.nc
    ln -f ${workdir_inpt}/RA2_daily_V10.nc    ${workdir2}/SFC_${MEM3}_daily_V10.nc
fi

exit 0
