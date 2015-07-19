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
ncatted=/sw/xe6/nco/4.0.8/sles11.1_netcdf4.2.0_gnu4.7.0/bin/ncatted

echo "LETKF preparation step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir_fcst=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/RESTART
workdir_inpt=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}/INPUT
workdir2=${EXP_DATA}/${YYYYMMDDHH}/letkf
workdir_alt=${EXP_DATA}/${YYYYMMDDHH}/model/??/RESTART
mkdir -p ${workdir2}

MEM3=`printf %.3d ${MEMBERID}`

#STEVE: active code:
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

  #-----------------------------------------------------------------------------
  # FIRST, link the background model data
  #-----------------------------------------------------------------------------

  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc     gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc      gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc           gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f $workdir_fcst/$IY$IM$ID.$IH$IN$IS.ocean_barotropic.res.nc  gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi

  #STEVE: add 'fill value' to netcdf files for identification of missing values
  cp $ncatted .
  ncatted -O -a _FillValue,temp,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ncatted -O -a _FillValue,salt,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ncatted -O -a _FillValue,u,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ncatted -O -a _FillValue,v,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ncatted -O -a _FillValue,sea_lev,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ncatted -O -a _FillValue,eta_t,o,f,-1.e+34 gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi

  #STEVE: need a template for the output analysis files:
  #       MUST COPY, NOT LINK. The files WILL be overwritten.
  if [ "${ISLOT2}" -eq "${ATIME}" ]; then #STEVE: trying to only do it once per member...
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/anal$MEM3.ocean_temp_salt.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/anal$MEM3.ocean_velocity.res.nc
    cp ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/anal$MEM3.ocean_sbc.res.nc
    if [ "$USE_ALTIMETRY" -eq "1" ]; then
      cp ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc  ${workdir2}/anal$MEM3.ocean_barotropic.res.nc
    fi
    #STEVE: may want to zero these out for peace of mind...
  fi

  #-----------------------------------------------------------------------------
  # SECOND, link the observations
  #-----------------------------------------------------------------------------

  # The next conditional gives the ability to use different observation source 
  # collection for different slots. (e.g. only use surface obs at time of 
  # analysis, use profiles through analysis cycle window)
  echo "For ISLOT2=$ISLOT2, and ATIME=${ATIME}"
  if [ "$ISLOT2" -eq "${ATIME}" ]; then
    OBSDIR=$OBSDIR5
  else
    OBSDIR=$OBSDIR1
  fi

  echo "Processing obs: ${OBSDIR}/$IY$IM$ID.dat to obs${ISLOT2}${MEM3}.dat"
  ln -f $INPUT_INIT/grid_spec.nc .
  cp $LDIR/$OBSOPexe .
  ln -f $OBSDIR/$IY$IM$ID.dat obsin.dat
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc gues.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  gues.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       gues.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc       gues.ocean_barotropic.res.nc
    echo "Linking: $altimetry_climatology_file to here..."
    if [ -f "$altimetry_climatology_file" ]; then
      ln -f $altimetry_climatology_file .
    else
      echo "$altimetry_climatology_file does not exist..."
      echo "Exiting..."
      exit 1
    fi
  fi

  ####################################################################################################################################
  # Running LETKF Obs Operator executable to generate the observation innovations for each member at each timestep:
  ####################################################################################################################################
# $OBSOPexe -obsin $OBSDIR/$IY$IM$ID.dat -gues gs${ISLOT2}$MEM3 -obsout ${workdir2}/obs${ISLOT2}${MEM3}.dat > obsope.log
  #STEVE: (perhaps make parallel and call with aprun)
  $OBSOPexe 

  if [ -f "obsout.dat" ]; then
    ln -f obsout.dat ${workdir2}/obs${ISLOT2}${MEM3}.dat
  else
    echo "output obs2 formatted file not created by $OBSOPexe."
    pwd
    ls
    echo "Exiting..."
    exit 2
  fi
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc ${workdir2}/gs${ISLOT2}$MEM3.ocean_temp_salt.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc  ${workdir2}/gs${ISLOT2}$MEM3.ocean_velocity.res.nc
  ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc       ${workdir2}/gs${ISLOT2}$MEM3.ocean_sbc.res.nc
  if [ "$USE_ALTIMETRY" -eq "1" ]; then
    ln -f ${workdir}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc       ${workdir2}/gs${ISLOT2}$MEM3.ocean_barotropic.res.nc
  fi

  #STEVE: the hard-link limit is running out (65000 max), so best to delete unnecessary links
  rm -f grid_spec.nc
  rm -f ncatted

done #DONE FCST loop

exit 0
