#!/bin/sh --login
#SBATCH -n 20      
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH -J letkf_dr_mom4p1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

IYYYYMMDDHH=1980010100
EYYYYMMDDHH=1980011100
YYYYMMDDHH=$IYYYYMMDDHH

echo "Start running the experiment..."

root=/lustre/lysun/models/Ocean-LETKF-test/mom4
root_run=${root}/run

ENS_NUM=5
days=5
isfirst=1

# sh pre_run.sh # generate obsop.001 obsop.DRIFTERS.001 letkf.DRIFTERS.005

while [ "${YYYYMMDDHH}" -lt "${EYYYYMMDDHH}" ]
do
  TYYYY=${YYYYMMDDHH:0:4}
  TMM=${YYYYMMDDHH:4:2}
  TDD=${YYYYMMDDHH:6:2}
  THH=${YYYYMMDDHH:8:2}
  TNN=00
  TSS=00 
  # sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst>  
  sh model_mom4p1.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM} ${isfirst}
  
  wait

  # letkf.sh <root> <TIME> <days> <ENS_NUM>
  sh letkf.sh ${root} ${YYYYMMDDHH} ${days} ${ENS_NUM}

  isfirst=0

  date=/bin/date
  sinc=$days
  sinc_units=days
  NYYYY=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%Y`
  NMM=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%m`
  NDD=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%d`
  NHH=`$date -d "$TYYYY-$TMM-$TDD $sinc $sinc_units" +%H`
  NNN=$TNN
  NSS=$TSS

  YYYYMMDDHH=${NYYYY}${NMM}${NDD}${NHH}

done

echo "NATURAL END."

exit 0
























