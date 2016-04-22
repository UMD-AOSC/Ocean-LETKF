#!/bin/ksh
YYYYMMDDHH=2009010100
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00
echo "IY = $IY"
echo "IM = $IM"
echo "ID = $ID"
echo "IH = $IH"

FCST=5
sinc=`expr $FCST - 1`

if [ "01" -lt "02" ]; then
  echo "sinc = $sinc"
fi

workdir=tmp/2009010100/model_prep/01
for file in `ls -d ${workdir}/RA2_daily_*.nc`; do
    echo "file = $file"
done
