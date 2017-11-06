#!/bin/csh -f
#

source cWG.header.csh

echo "Using inputs:"
echo "altdir :: $altdir"
echo "yr :: $yr"
echo "mo :: $mo"
echo "dy :: $dy"
echo "outdir :: $outdir"
echo "NOTE: grid_spec.nc must be in this directory!"
#echo "nlv :: $nlv"

#STEVE: remove leading zeros, or else csh interprets this as an octal number
set yr = `echo $yr | sed 's/0*//'`
set mo = `echo $mo | sed 's/0*//'`
set dy = `echo $dy | sed 's/0*//'`

echo "yr = $yr"
echo "mo = $mo"
echo "dy = $dy"

#set altdir = /lustre/fs/scratch/ncep/David.Behringer/ptmp/AVISO/CYCLE
@ nlv = 40

#set ymd = (`tail -1 time_stamp.out`)
#@ yr = $ymd[1]
#@ mn = $ymd[2]
#@ dy = $ymd[3]
set dte = `printf "%4.4d%2.2d%2.2d" $yr $mn $dy`

echo "Create and combine weekly altimetry files:"
echo " "
echo "Begin"

#ln -sf ../../INPUT/grid_spec.nc grid_spec.nc

runSWA4 -d $dte -n 5 -f txj1j2_cTbl -p $altdir -m grid_spec.nc -lt -40 50 -k $nlv

ls -1 swssh.?????? > cWGLst

cmbWksAlt4 -f cWGLst -o alta.mom

mv alta.mom $outdir

rm cWGLst
rm swssh.??????
#rm grid_spec.nc

