#!/bin/bash --login

set -exv

# sh make_obsop_drifters.sh $MEM
CONFIGDIR=../../..
source $CONFIGDIR/config/machine.sh
source $CONFIGDIR/config/$MACHINE.fortran.sh
source $CONFIGDIR/config/$MACHINE.netcdf.sh

PGM=p2l_synth.x

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c SFMT.f90 
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c common.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c params_model.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c vars_model.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c params_letkf.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c common_mom4.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c params_obs.f90
$F90 $OMP $F90OPT $INLINE $NETCDF_INC -c common_obs_mom4.f90

$F90 $OMP $F90OPT $NETCDF_INC -c prof2letkf.f90
$F90 $OMP $F90OPT $INLINE -o ${PGM} *.o $NETCDF_LIB

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
