#!/bin/sh
F90=ftn
LDIR=/autofs/na1_home1/Steve.Penny/letkf/mom4/letkf_GAEA
OBSB=/lustre/f1/unswept/Steve.Penny/OUTPUT/tmp_letkf_ts/STORE/fcst/omf/01/19910111.dat
OBSA=/lustre/f1/unswept/Steve.Penny/OUTPUT/tmp_letkf_ts/STORE/anal/oma/01/19910111.dat
ln -sf $OBSB fort.3
ln -sf $OBSA fort.4
$F90 $LDIR/obs2_verify.f90 -o o2v.x
./o2v.x
rm o2v.x
rm fort.3
rm fort.4
