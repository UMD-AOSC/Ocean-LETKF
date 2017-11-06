#!/bin/bash
# 
set -e
CWD=`pwd`
#cd ..
#wdir=`pwd`
#cd $CWD

#STEVE: mpi command for Gaea system:
MPIRUN= #"aprun -n 1 "

expdir=$1
SBC_DIR=$2
SBC_DIR2=$3
SST_DIR=$4
sfcdays=$5 #`expr $5 + 1`
echo "Using SBC source directory1: $SBC_DIR"
echo "Using SBC source directory2: $SBC_DIR2"
echo "Using SST source directory1: $SST_DIR"
echo "Using sfcdays = $sfcdays"
echo "REQUIRES: grid_spec.nc and salt12.nc in working directory"

touch RA2_daily_dum.nc
touch dum_sfc_restore.nc
rm RA2_daily_*.nc
rm *_sfc_restore.nc

L_YR=2013

#cp ../time_stamp.out .
ymd=(`tail -1 time_stamp.out`)
yr=${ymd[0]}
mo=${ymd[1]}
dy=${ymd[2]}
#STEVE: fix for case where 08 and 09 are treated as octal instead of decimal:
mo=$(( 10#$mo ))
dy=$(( 10#$dy ))
dte=`printf "%4.4d%2.2d%2.2d" $yr $mo $dy`

echo "Making TFLUX.nc : sensible heat W/m^2"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 TFLUX -d $dte -n $sfcdays -o RA2_daily_TFLUX.nc -y $L_YR &

echo "Making SHRTWV.nc : net shortwave W/m^2"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR2}/R2CR SHRTWV -d $dte -n $sfcdays -o RA2_daily_SHRTWV.nc -y $L_YR &

echo "Making LONGWV.nc : net longwave W/m^2"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 LONGWV -d $dte -n $sfcdays -o RA2_daily_LONGWV.nc -y $L_YR &

echo "Making PRATE.nc : precipatation rate kg/m^2/sec"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR2}/R2CR PRATE -d $dte -n $sfcdays -o RA2_daily_PRATE.nc -y $L_YR &

echo "Making QFLUX.nc : evaporation rate kg/m^2/sec"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 QFLUX -d $dte -n $sfcdays -o RA2_daily_QFLUX.nc -y $L_YR &

echo "Making TAUX.nc : X-momentum flux Pa=N/m^2"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR2}/R2CR TAUX -d $dte -n $sfcdays -o RA2_daily_TAUX.nc -y $L_YR &

echo "Making TAUY.nc : Y-momentum flux Pa=N/m^2"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR2}/R2CR TAUY -d $dte -n $sfcdays -o RA2_daily_TAUY.nc -y $L_YR &

#echo "Making U10.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 U10 -d $dte -n $sfcdays -o RA2_daily_U10.nc -y $L_YR &
#echo "Making V10.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 V10 -d $dte -n $sfcdays -o RA2_daily_V10.nc -y $L_YR &

#echo "Making DSW.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 dsw -d $dte -n $sfcdays -o RA2_daily_dsw.nc -y $L_YR &
#echo "Making DLW.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 dlw -d $dte -n $sfcdays -o RA2_daily_dlw.nc -y $L_YR &

#echo "Making T2m.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 t2m -d $dte -n $sfcdays -o RA2_daily_t2m.nc -y $L_YR &
#echo "Making Q2m.nc"
#$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 q2m -d $dte -n $sfcdays -o RA2_daily_q2m.nc -y $L_YR &

echo "Making pres.nc"
$MPIRUN $CWD/mkDlySBCnc4 -f ${SBC_DIR}/RA2 pres -d $dte -n $sfcdays -o RA2_daily_pres.nc -y $L_YR &

# check for leap year and possibly adjust length of sfcdays

#remainder=`expr $yr % 4`
#if ( $remainder == 0 ) then
#  if ( $mo == 2 && $dy == 25 ) then
#    sfcdays=6
#  fi
#fi

#ln -sf ${wdir}/INPUT/grid_spec.nc grid_spec.nc

echo "Making temp_sfc_restore.nc"
$MPIRUN $CWD/mkDlySst4i -p ${SST_DIR} -g grid_spec.nc -y 65 -d $dte -n $sfcdays -o temp_sfc_restore.nc &

#ln -sf ${wdir}/INPUT/salt12.nc salt12.nc

echo "Making salt_sfc_restore.nc"
$MPIRUN $CWD/mkDlySss4nci -f salt12.nc -g grid_spec.nc -y 65 -d $dte -n $sfcdays -o salt_sfc_restore.nc &

time wait

exit 0
