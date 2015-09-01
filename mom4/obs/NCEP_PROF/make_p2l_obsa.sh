#!/bin/bash --login

set -exv

CONFIGDIR=../../..
source $CONFIGDIR/config/machine.sh
source $CONFIGDIR/config/$MACHINE.fortran.sh
source $CONFIGDIR/config/$MACHINE.netcdf.sh

PGM=p2l_obsa.x

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INCLUDE_NETCDF -c prof2letkf_obsa.f90
$F90 $OMP $F90OPT $INLINE -o ${PGM} *.o $LIB_NETCDF

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
