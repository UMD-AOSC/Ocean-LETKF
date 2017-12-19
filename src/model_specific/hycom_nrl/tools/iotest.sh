#!/bin/bash

SRCDIR=/u/prob/scratch/relo/lig1/analysis/restart

YYYY=2017
MM=12
DD=01
HH=00

NX4=0187
NY4=0121

# grids are 187x121
# The naming convention puts the lowest level (2978m) in the filename, 
# but there are 46 levels in this analysis
GRIDFILE_lon=grdlon_sfc_000000_000000_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00000000_datafld
GRIDFILE_lat=grdlat_sfc_000000_000000_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00000000_datafld
GRIDFILE_lev=grdlvl_pre_000000_000000_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00000000_datafld
GRIDFILE_kmt=maskls_pre_000000_000000_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00000000_datafld

ln -fs $SRCDIR/$GRIDFILE_lon ./glon.dat
ln -fs $SRCDIR/$GRIDFILE_lat ./glat.dat
ln -fs $SRCDIR/$GRIDFILE_lev ./glev.dat
ln -fs $SRCDIR/$GRIDFILE_kmt ./kmt.dat

TFILE=seatmp_pre_000000_002978_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00240000_fcstfld
SFILE=salint_pre_000000_002978_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00240000_fcstfld
UFILE=uucurr_pre_000000_002978_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00240000_fcstfld
VFILE=vvcurr_pre_000000_002978_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00240000_fcstfld

SSHFILE=seahgt_sfc_000000_000000_1o${NX4}x${NY4}_${YYYY}${MM}${DD}${HH}_00240000_fcstfld

ftype='gues'
ln -fs $SRCDIR/$TFILE ${ftype}.t.dat
ln -fs $SRCDIR/$SFILE ${ftype}.s.dat
ln -fs $SRCDIR/$UFILE ${ftype}.u.dat
ln -fs $SRCDIR/$VFILE ${ftype}.v.dat
ln -fs $SRCDIR/$SSHFILE ${ftype}.ssh.dat

# Set up the analysis file template:
ftype1='gues'
ftype2='anal'
for var in t s u v ssh
do
  cp -f ${ftype1}.${var}.dat ${ftype2}.${var}.dat
done

./iotest.x
