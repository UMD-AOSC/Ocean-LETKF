#!/bin/bash
set -ex

source ../../../config/machine.sh
source ../../../config/$MACHINE.fortran.sh
source ../../../config/$MACHINE.netcdf.sh

CDIR=../../../common
LDIR=../../letkf

PGM=p2l_rawnc.x
#F90OPT='-ip -mcmodel=medium -shared-intel -fast -O2' #'-O3 -xSSE4.2 -ip -mcmodel=medium -shared-intel'

rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c $CDIR/SFMT.f90
$F90 $OMP $F90OPT $INLINE -c $CDIR/common.f90
$F90 $OMP $F90OPT $INLINE -c $LDIR/params_obs.f90
$F90 $OMP $F90OPT $INLINE -c $LDIR/params_model.f90
##$F90 $OMP $F90OPT $INLINE -c calendar.f90
$F90 $OMP $F90OPT $INCLUDE_NETCDF -c prof2letkf_rawnc.f90
$F90 $OMP $F90OPT $INLINE -o ${PGM} *.o $LIB_NETCDF

rm -f *.mod
rm -f *.o

echo "NORMAL END"
