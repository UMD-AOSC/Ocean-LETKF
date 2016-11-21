#!/bin/ksh --login
# This script prepares background and observation data for letkf 
#
# Written by Dr. Stephen G. Penny
#
set -e
#module load intel
#module load netcdf
#module load nco

echo "LETKF preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
echo "Posting forecast ${FCST}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf_prep/${MEMBERID}/${FCST}
mkdir -p ${workdir}
cd ${workdir}
workdir_fcst=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/RESTART
workdir_inpt=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/INPUT

echo "This is the LETKF preparation step for member ${MEMBERID}, forecast ${FCST}, for cycle ${YYYYMMDDHH}" > $workdir/letkf_prep_${FCST}.out

#STEVE: active code:
#DO_SFCFLUXES=1  #This is input via the xml script
USE_SFC=1
USE_TRI2SPH=0
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the ISLOT date
date=/bin/date
#sinc=`expr $FCST - 1` #STEVE: use this if I want the dates offset
sinc=$FCST
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

# FIRST, link the (0-360 degree spherical) topography for letkf to read grid information

MEM3=`printf %.3d ${MEMBERID}`
ISLOTL=`printf %.2d ${FCST}`

ln -fs $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc gs${ISLOTL}$MEM3.ocean_temp_salt.res.nc
ln -fs $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  gs${ISLOTL}$MEM3.ocean_velocity.res.nc
# if [ $USE_SFC ]; then
  ln -fs $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       gs${ISLOTL}$MEM3.ocean_sbc.res.nc
# fi

#STEVE: add 'fill value' to netcdf files
ncatted -O -a _FillValue,temp,o,f,-1.e+34 gs${ISLOTL}$MEM3.ocean_temp_salt.res.nc
ncatted -O -a _FillValue,salt,o,f,-1.e+34 gs${ISLOTL}$MEM3.ocean_temp_salt.res.nc
ncatted -O -a _FillValue,u,o,f,-1.e+34 gs${ISLOTL}$MEM3.ocean_velocity.res.nc
ncatted -O -a _FillValue,v,o,f,-1.e+34 gs${ISLOTL}$MEM3.ocean_velocity.res.nc
ncatted -O -a _FillValue,sea_lev,o,f,-1.e+34 gs${ISLOTL}$MEM3.ocean_sbc.res.nc

# FIRST (a), Link background files to the letkf working directory:
workdir2=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir2}
for file in `ls -d ${workdir}/gs${ISLOTL}$MEM3.*`; do
  echo "linking $file to ${workdir2}..."
  ln -fs $file ${workdir2}/
done

#STEVE: need a template for the output files:
if [ "${ISLOTL}" -eq "${ATIME}" ]; then #STEVE: trying to only do it once per member...
  cp ${workdir}/gs${ISLOTL}$MEM3.ocean_temp_salt.res.nc ${workdir2}/anal$MEM3.ocean_temp_salt.res.nc
  cp ${workdir}/gs${ISLOTL}$MEM3.ocean_velocity.res.nc  ${workdir2}/anal$MEM3.ocean_velocity.res.nc
#if [ $USE_SFC ]; then
  cp ${workdir}/gs${ISLOTL}$MEM3.ocean_sbc.res.nc     ${workdir2}/anal$MEM3.ocean_sbc.res.nc
#fi
fi

#STEVE: DO_SFCFLUXES
#if [ "${DO_SFCFLUXES}" -eq "1" ]; then
# for file in `ls -d ${workdir_inpt}/RA2_*.nc`; do
#   echo "linking $file to ${workdir2}..."
#   ln -fs $file ${workdir2}/${MEM3}.${file}
# done
    ln -fs ${workdir_inpt}/RA2_daily_dlw.nc    ${workdir2}/SFC_${MEM3}_daily_dlw.nc
    ln -fs ${workdir_inpt}/RA2_daily_dsw.nc    ${workdir2}/SFC_${MEM3}_daily_dsw.nc
    ln -fs ${workdir_inpt}/RA2_daily_LONGWV.nc ${workdir2}/SFC_${MEM3}_daily_LONGWV.nc
    ln -fs ${workdir_inpt}/RA2_daily_PRATE.nc  ${workdir2}/SFC_${MEM3}_daily_PRATE.nc
    ln -fs ${workdir_inpt}/RA2_daily_pres.nc   ${workdir2}/SFC_${MEM3}_daily_pres.nc
    ln -fs ${workdir_inpt}/RA2_daily_q2m.nc    ${workdir2}/SFC_${MEM3}_daily_q2m.nc
    ln -fs ${workdir_inpt}/RA2_daily_QFLUX.nc  ${workdir2}/SFC_${MEM3}_daily_QFLUX.nc
    ln -fs ${workdir_inpt}/RA2_daily_SHRTWV.nc ${workdir2}/SFC_${MEM3}_daily_SHRTWV.nc
    ln -fs ${workdir_inpt}/RA2_daily_t2m.nc    ${workdir2}/SFC_${MEM3}_daily_t2m.nc
    ln -fs ${workdir_inpt}/RA2_daily_TAUX.nc   ${workdir2}/SFC_${MEM3}_daily_TAUX.nc
    ln -fs ${workdir_inpt}/RA2_daily_TAUY.nc   ${workdir2}/SFC_${MEM3}_daily_TAUY.nc
    ln -fs ${workdir_inpt}/RA2_daily_TFLUX.nc  ${workdir2}/SFC_${MEM3}_daily_TFLUX.nc
    ln -fs ${workdir_inpt}/RA2_daily_U10.nc    ${workdir2}/SFC_${MEM3}_daily_U10.nc
    ln -fs ${workdir_inpt}/RA2_daily_V10.nc    ${workdir2}/SFC_${MEM3}_daily_V10.nc
#fi

# SECOND (b), process and link the observations
if [ ${MEMBERID} -lt "02" ]; then #STEVE: trying to only do it once per slot...
 #Use different observation source collection for different slots.
 # (i.e. only use surface obs at time of analysis, use profiles through analysis cycle window)
echo "For ISLOTL=$ISLOTL, and ATIME=${ATIME}"
 if [ "$ISLOTL" -eq "${ATIME}" ]; then
  if [ -s ${OBSDIR5}/$IY$IM$ID$IH.dat ]; then
    echo "Linking obs: ${OBSDIR5}/$IY$IM$ID$IH.dat to obs${ISLOTL}.dat"
    ln -fs $OBSDIR5/$IY$IM$ID$IH.dat ${workdir2}/obs${ISLOTL}.dat
  else
    echo "Linking obs: ${OBSDIR5}/$IY$IM$ID.dat to obs${ISLOTL}.dat"
    ln -fs $OBSDIR5/$IY$IM$ID.dat ${workdir2}/obs${ISLOTL}.dat
  fi
 else
  if [ -s ${OBSDIR1}/$IY$IM$ID$IH.dat ]; then
    echo "Linking obs: ${OBSDIR1}/$IY$IM$ID$IH.dat to obs${ISLOTL}.dat"
    ln -fs $OBSDIR1/$IY$IM$ID$IH.dat ${workdir2}/obs${ISLOTL}.dat
  else
    echo "Linking obs: ${OBSDIR1}/$IY$IM$ID.dat to obs${ISLOTL}.dat"
    ln -fs $OBSDIR1/$IY$IM$ID.dat ${workdir2}/obs${ISLOTL}.dat
  fi
 fi
fi

#ln -s $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd gs01${MEM3}.grd
#ln -fs $OUTPUT/gues/$MEM/$IYYYY$IMM$IDD$IHH.grd gues.grd
#./$OBSOPE > obsope.log
#mv obsout.dat obs01${MEM3}.dat


exit 0
