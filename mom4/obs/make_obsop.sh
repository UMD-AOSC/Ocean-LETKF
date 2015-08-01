#!/bin/sh
set -exv

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

F90_OPT=''

$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG common_mom4.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG common_obs_mom4.f90
#--
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_oceanographic_toolbox.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_pot_to_insitu.f90
#--
$F90 $OMP $F90OPT $F90_OBJECT_FLAG obsop.f90
$F90 $OMP $F90OPT -o ${PGM} *.o $NETCDF_LIB

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh
#--
#--

echo "NORMAL END"
