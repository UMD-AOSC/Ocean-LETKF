#!/bin/sh
F90=ftn
#OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep/20020101.dat
OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep/20030712.dat
ln -s $OBS fort.3
$F90 obsdump.f90
./a.out
rm a.out
rm fort.3
