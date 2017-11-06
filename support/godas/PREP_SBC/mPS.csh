#!/bin/csh
#
#source ~/mom4p1/bin/environs.ia64
module load wgrib
#export PATH=$PATH:`pwd` #STEVE: need to access wgrib and nco executables copied to here, since modules are not available.
set CWD = `pwd`
#setenv PATH $PATH : $CWD
set path = ( $path $CWD )

#STEVE:
# INPUTS:
#
# SST_DIR (sst reanalysis data)
# SSS_DIR (sss climatology data)
# L_YR (last year of data)
# gridDir
# scriptDir
# yr
# mo
# dy

source mPS.header.csh

#set SST_DIR=${dataDir}/SST2/DAILY
#set SSS_DIR=${wdir}/INPUT
#@ L_YR = 2013
#set gridDir=${wdir}/INPUT

touch dum_sfc_restore.nc
rm *_sfc_restore.nc

#STEVE: remove leading zeros, or else csh interprets this as an octal number
set yr = `echo $yr | sed 's/0*//'`
set mo = `echo $mo | sed 's/0*//'`
set dy = `echo $dy | sed 's/0*//'`

set dte = `printf "%4.4d%2.2d%2.2d" $yr $mo $dy`

# check for leap year and possibly adjust length of days

@ days = 10

@ remainder = $yr % 4

if ( $remainder == 0 ) then
   if ( $mo == 2 && $dy >= 20 ) then
     @ days = 11
   endif
endif

ln -sf $gridDir/grid_spec.nc grid_spec.nc

echo "Making temp_sfc_restore.nc"
$scriptDir/mkDlySst4i -p ${SST_DIR} -g grid_spec.nc -y 65 -d $dte -n $days -o temp_sfc_restore.nc

ln -sf $SSS_DIR/salt12.nc salt12.nc

echo "Making salt_sfc_restore.nc"
$scriptDir/mkDlySss4nci -f salt12.nc -g grid_spec.nc -y 65 -d $dte -n $days -o salt_sfc_restore.nc

