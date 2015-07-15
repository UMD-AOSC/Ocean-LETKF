#!/bin/bash
set -e
date=date
CWD=`pwd`

INDIR1=~/lf1u/OBS/historical/ARGO_raw
INDIR2=~/lf1u/OBS/historical/ARGO_raw
OUTDIR=~/lf1u/OBS/historical/letkf_fmt/PROFS_rawnc_TS
#ZIPDIR=$CWD/ARCS
IY=2005
IM=01
ID=01
IH=00

EY=2005
EM=01
ED=01
EH=00

inc=1
inc_units='days'

mkdir -p $OUTDIR

# Unzip and untar profiles for first year's run
#rm -f $INDIR/*.nc
#cd $INDIR
#tar -xvf $ZIPDIR/${IY}tmp.tar
#tar -xvf $ZIPDIR/${IY}sal.tar

# MAIN LOOP
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
  cd $CWD
  TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
  TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
  TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
  TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`

  echo "Computing for $IY$IM$ID"
  ./p2l_rawnc.x -y $IY -m $IM -d $ID -indir1 $INDIR1 -indir2 $INDIR2 -outdir $OUTDIR > $OUTDIR/$IY$IM$ID.out

  # Untar profiles for subsequent years' runs
# if [ $IY -ne $TY ]; then
#   rm -f $INDIR/*.nc
#   cd $INDIR
#   tar -xvf $ZIPDIR/${TY}tmp.tar
#   tar -xvf $ZIPDIR/${TY}sal.tar
# fi

  IY=$TY
  IM=$TM
  ID=$TD
  IH=$TH

done
