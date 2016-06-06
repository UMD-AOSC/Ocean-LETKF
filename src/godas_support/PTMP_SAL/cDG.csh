#!/bin/csh
#

#STEVE:
# INPUTS:
#
# arcdir (directory with profile observations)
# d1, d2 (range of dates to look forward and backward for observations)
# gridDir (directory with grid_spec.nc file)
# scriptDir (directory of support C script)
# yr
# mo
# dy

source cDG.header.csh

echo "Using inputs:"
echo "arcdir_t :: $arcdir_t"
echo "arcdir_s :: $arcdir_s"
echo "d1, d2 :: $d1, $d2"
echo "gridDir :: $gridDir"
echo "scriptDir :: $scriptDir"

#STEVE: remove leading zeros, or else csh interprets this as an octal number
set yr = `echo $yr | sed 's/0*//'`
set mo = `echo $mo | sed 's/0*//'`
set dy = `echo $dy | sed 's/0*//'`

echo "yr = $yr"
echo "mo = $mo"
echo "dy = $dy"

#set arcdir = /climate/save/wx24db/G_PROFILES/AUTOnc
#set arcdir = /scratch2/portfolios/NCEPDEV/climate/noscrub/David.Behringer/ARCS/DAILYnc
#set arcdir = /home/Steve.Penny/godas4p1/gdsSolo/OBS

#set gridDir = $expdir/INPUT
#set INPUT = $expdir/INPUT

set dte = `printf "%4.4d%2.2d%2.2d" $yr $mo $dy`
echo "dte = $dte"

#@ d1 = -5
#@ d2 = 5

set dt1 = `./incDate $dte $d1`
set dt2 = `./incDate $dte $d2`

echo "======================="
echo "From cDG.csh...        "
echo "GODAS OBS date1 :: $dt1"
echo "GODAS OBS date2 :: $dt2"
echo "======================="

set yr1 = `echo $dt1 | cut -c -4`
set mn1 = `echo $dt1 | cut -c 5-6`
set dy1 = `echo $dt1 | cut -c 7-8`

echo $yr1 $mn1 $dy1

set yr2 = `echo $dt2 | cut -c -4`
set mn2 = `echo $dt2 | cut -c 5-6`
set dy2 = `echo $dt2 | cut -c 7-8`

echo $yr2 $mn2 $dy2

set opt_t = O
$scriptDir/extDysP4nc -p $arcdir_t -t $opt_t -d $dte -b $d1 $d2

ls ????????tmp.nc > tmpFiles

#set opt_t = M
set opt_t = O
$scriptDir/extDysP4nc -p $arcdir_s -t $opt_t -d $dte -b $d1 $d2 -s

ls ????????sal.nc > salFiles

ln -f $gridDir/grid_spec.nc grid_spec.nc

# Create the observation data files
set nlev = 40
$scriptDir/cmbDLstTh4nc -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k $nlev
$scriptDir/cmbDLstPs4nc -f salFiles -o sala.mom -g grid_spec.nc -k $nlev

# Clear out the temporary files
rm ????????tmp.nc
rm ????????sal.nc
rm tmpFiles
rm salFiles
