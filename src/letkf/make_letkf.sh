#!/bin/bash
set -ex

source ../../config/machine.sh
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh
source ../../config/$MACHINE.mpi.sh

# Ensemble size
# STEVE: figure out how to read from params_letkf.f90 and put here (e.g. with awk/perl/etc.)
#        -> grep and sed seem to work ok:
#        (Set the ensemble size in params_letkf.f90, it will read it in here)
#MEM=`grep nbv= params_letkf.f90 | sed -r 's/INTEGER,PARAMETER :: nbv=([0-9]+)/\1/'`
#MEM=`grep nbv= params_letkf.f90 | sed -r 's/INTEGER :: nbv=([0-9]+)/\1/'`
#echo "MEM=$MEM"
#MEM3=`printf %.3d ${MEM}`

# Experiment name
#name=TEST5
#name=TEST5_dyn
name=TEST6
# Executable for letkf
PGM=letkf.$name.x

OMP=
PWD=`pwd`
#BLAS=1 #0: no blas 1: using blas
#STEVE: ask about blas on zeus
sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

cat netlib.f > netlib2.f
if test $BLAS -eq 1
then
  LBLAS="-L${CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_intel -lsci_intel_mp"
else
  cat netlibblas.f >> netlib2.f
  LBLAS=""
fi

F90_FPP='-fpp' # Fortran preprocessor
# If dynamic model grid specification is desired (e.g. via input.nml namelist):
#FPP_FLAG='-DDYNAMIC'
FPP_FLAG=""

$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_mpi.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG common_mtx.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG netlib2.f
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $FPP_FLAG $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $FPP_FLAG $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG common_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_obs_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $NETCDF_INC common_mpi_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG vars_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG kdtree.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_local.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_local.o letkf_tools.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $FPP_FLAG $F90_OBJECT_FLAG letkf.f90

# Generate executable:
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE -o ${PGM} *.o $MPI_LIB $NETCDF_LIB $LBLAS

#STEVE: keep a record of the build:
mkdir -p CONFIG_$PGM
cp *.f90 CONFIG_$PGM/

rm -f *.mod
rm -f *.o
rm -f netlib2.f
sh ulnkcommon.sh

echo "STEVE: don't forget - in phys2ijk, obs above model level 1 are set to model level 1"
echo "NORMAL END"
