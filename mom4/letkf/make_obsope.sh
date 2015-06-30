#!/bin/sh
set -ex
MEM=056
#MEM=028
#MEM=002
#MEM=004
PGM=obsope_alt.$MEM
#PGM=obsope.$MEM
F90=ftn #pgf90
OMP=
#F90OPT='-ftz -ip -ipo -O2 -parallel -i_dynamic -assume byterecl -i4 -r8 -what -fpp -fno-alias -stack_temps -safe_cray_ptr -fast'
F90OPT='-ftz -ip -ipo -O2 -parallel -i_dynamic -what -fpp -fno-alias -stack_temps -safe_cray_ptr -fast'
INLINE= #"-Minline"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE -c SFMT.f90
$F90 $OMP $F90OPT $INLINE -c common.f90
$F90 $OMP $F90OPT $INLINE -c common_mom4.f90
$F90 $OMP $F90OPT -c common_obs_mom4.f90
$F90 $OMP $F90OPT -c obsope.f90
$F90 $OMP $F90OPT -o ${PGM} *.o

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
