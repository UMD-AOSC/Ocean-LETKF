#!/bin/bash
set -e
date=date
CWD=`pwd`

arcdir_t=~/lf1u/OBS/historical/TMP_profs
arcdir_so=~/lf1u/OBS/historical/SAL_profs_O #observed only
arcdir_sm=~/lf1u/OBS/historical/SAL_profs_M #observed and synthetic, for potential temperature computation

scriptDir=/autofs/na1_home1/Steve.Penny/letkf/mom4/obs/NCEP_PROF/PTMP_SAL/SRC
oaTexe=cmbDLstTh4nc_obsa
oaSexe=cmbDLstPs4nc_obsa

INDIR1=~/lf1u/OBS/historical/TMP_profs_tmpa_deep
INDIR2=~/lf1u/OBS/historical/SAL_profs_O_sala_deep
OUTDIR=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep
mkdir -p $INDIR1
mkdir -p $INDIR2
mkdir -p $OUTDIR
#ZIPDIR=$CWD/ARCS

gridDir=$CWD

IY=1991
IM=01
ID=01
IH=00

EY=1998
EM=12
ED=31
EH=00

#IY=1999
#IM=01
#ID=01
#IH=00

#EY=2012
#EM=01
#ED=01
#EH=00

inc=1
inc_units='days'

mkdir -p $OUTDIR

# Unzip and untar profiles for first year's run
#rm -f $INDIR/*.nc
#cd $INDIR
#tar -xvf $ZIPDIR/${IY}tmp.tar
#tar -xvf $ZIPDIR/${IY}sal.tar

# MAIN LOOP
cd $OUTDIR
ln -fs $gridDir/grid_spec.nc grid_spec.nc
while test $IY$IM$ID$IH -le $EY$EM$ED$EH
do
  cd $CWD
  TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
  TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
  TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
  TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`

  echo "==================================================== "
  echo "  Computing for $IY$IM$ID"
  echo "==================================================== "

  if [ 1 -eq 1 ]; then
    # Create the csh-script prep files:
    rm -f tmpFiles
    echo ${IY}${IM}${ID}tmp.nc > tmpFiles

    rm -f salFiles
    echo ${IY}${IM}${ID}sal.nc > salFiles

    rm -f ????????tmp.nc
    rm -f ????????sal.nc
    ln -fs $arcdir_t/${IY}${IM}${ID}tmp.nc .
    ln -fs $arcdir_sm/${IY}${IM}${ID}sal.nc .

    # Run the netcdf to obsa-format conversion:
    # Compute temperature observations & errors
    rm -f tmpa.mom
    $scriptDir/$oaTexe -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 40
    mv tmpa.mom $INDIR1/${IY}${IM}${ID}.tmpa.mom

    rm -f salFiles
    echo ${IY}${IM}${ID}sal.nc > salFiles
  
    rm -f ????????sal.nc
    if [ -f $arcdir_so/${IY}${IM}${ID}sal.nc ]; then
      # Run the netcdf to obsa-format conversion:
      # Compute salinity observations & errors
      ln -fs $arcdir_so/${IY}${IM}${ID}sal.nc .
      rm -f sala.mom
      $scriptDir/$oaSexe -f salFiles -o sala.mom -g grid_spec.nc  -k 40
      mv sala.mom $INDIR2/${IY}${IM}${ID}.sala.mom
    else
      echo "Does not exist: $arcdir_so/${IY}${IM}${ID}sal.nc"
#     exit 1
    fi
  fi

  # Run the obsa-format to letkf-format conversion:
  echo "Running p2l_obsa.x ..."
  $CWD/p2l_obsa.x -y $IY -m $IM -d $ID -indir1 $INDIR1 -indir2 $INDIR2 -outdir $OUTDIR > $OUTDIR/$IY$IM$ID.out

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
