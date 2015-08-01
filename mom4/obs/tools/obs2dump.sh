#!/bin/sh
F90=ftn
OBS=/lustre/f1/unswept/Steve.Penny/OUTPUT/tmp_letkf_mom6/2005010100/letkf/obs01001.dat
ln -s $OBS fort.3
$F90 obs2dump.f90
./a.out
rm a.out
rm fort.3
