#!/bin/sh
set -ex

source ../../config/machine.sh
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh

MEM=056
PGM=obsop.$MEM
#F90OPT='-ftz -ip -ipo -O2 -parallel -i_dynamic -what -fpp -fno-alias -stack_temps -safe_cray_ptr -fast'

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT -c params_model.f90
$F90 $OMP $F90OPT -c vars_model.f90
$F90 $OMP $F90OPT -c params_letkf.f90
$F90 $OMP $F90OPT $INLINE -c common_mom4.f90
$F90 $OMP $F90OPT -c params_obs.f90
$F90 $OMP $F90OPT -c vars_obs.f90
$F90 $OMP $F90OPT -c common_obs_mom4.f90
$F90 $OMP $F90OPT -c obsop.f90
$F90 $OMP $F90OPT -o ${PGM} *.o

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
