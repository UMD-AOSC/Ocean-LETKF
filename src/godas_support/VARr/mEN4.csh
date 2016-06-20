#!/bin/csh 
#

#STEVE: This file creates the background error covariance matrix to be used by godas

#STEVE:
# INPUTS:
# RESTART (RESTART directory from last ocean model run)
# gridDir (directory of grid_spec.nc file)
# scriptDir (directory of support C script)
# yr
# mo
# dy
source mEN4.header.csh

#set RESTART = $expdir/RESTART
#set gridDir = $expdir/INPUT

rm tvv.mom
rm svv.mom

#STEVE: remove leading zeros, or else csh interprets this as an octal number
set yr = `echo $yr | sed 's/0*//'`
set mo = `echo $mo | sed 's/0*//'`
set dy = `echo $dy | sed 's/0*//'`
set dte = `printf "%4.4d%2.2d%2.2d" $yr $mo $dy`

echo $dte

if ( ! -f $RESTART/ocean_temp_salt.res.nc ) then
  echo "Missing ocean_temp_salt.res.nc file: $RESTART/ocean_temp_salt.res.nc"
  echo "Cannot compute background VAR"
  exit 99
else
  ln -sf $RESTART/ocean_temp_salt.res.nc TS.nc
endif

if ( ! -f TS.nc ) then
  echo "$0"
  pwd
  echo "Missing TS.nc file"
  echo "Cannot compute background VAR"
  exit 99
endif

if ( ! -f $gridDir/grid_spec.nc ) then
  echo "Missing grid_spec.nc file: $gridDir/grid_spec.nc"
  echo "Cannot compute background VAR"
  exit 99
else
  ln -sf $gridDir/grid_spec.nc gSpec.nc
endif

if ( ! -f gSpec.nc ) then
  echo "$0"
  pwd
  echo "Missing gSpec.nc file"
  echo "Cannot compute background VAR"
  exit 99
endif

#$scriptDir/mkEvNc4r -f TS.nc -gr gSpec.nc -d $dte -o tvv.mom -p -Sv -so svv.mom -k 30 -SqRt 1 -sc 1.0
$scriptDir/mkEvNc4r -f TS.nc -gr gSpec.nc -d $dte -o tvv.mom -p -Sv -so svv.mom -k 40 -SqRt 1 -sc 1.0

rm TS.nc
rm gSpec.nc
