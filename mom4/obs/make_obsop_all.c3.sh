#!/bin/sh
set -exv

source $MODULESHOME/init/sh
module use /sw/eslogin-c3/modulefiles

#source ../../config/machine.sh
MACHINE="GAEA3"
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

name=TESTc3
PGM=obsop.$name

$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG common_mom4.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG common_obs_mom4.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG compute_profile_error.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_argo.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_avhrr_pathfinder.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_aviso_adt.f90
#--
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_oceanographic_toolbox.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_pot_to_insitu.f90
#--
#$F90 $OMP $F90_OPT $F90_OBJECT_FLAG obsop_tprof.f90
#$F90 $OMP $F90_OPT $F90_OBJECT_FLAG obsop_sprof.f90
#$F90 $OMP $F90_OPT $F90_OBJECT_FLAG obsop_adt.f90
#$F90 $OMP $F90_OPT $F90_OBJECT_FLAG obsop_sst.f90
$F90 $OMP $F90_OPT obsop_tprof.f90 -o ${PGM}.tprof.x *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_sprof.f90 -o ${PGM}.sprof.x *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_adt.f90   -o ${PGM}.adt.x   *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_sst.f90   -o ${PGM}.sst.x   *.o $NETCDF_LIB

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh
#--
#--

echo "NORMAL END"
