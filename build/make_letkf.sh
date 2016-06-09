#!/bin/bash
set -ex

CDIR=`pwd`

# Machine configuration
CONFIGDIR=../config
source $CONFIGDIR/machine.sh
source $CONFIGDIR/$MACHINE.fortran.sh
source $CONFIGDIR/$MACHINE.netcdf.sh
source $CONFIGDIR/$MACHINE.mpi.sh
source $CONFIGDIR/$MACHINE.modules.sh

# Model name
model=mom4 #mom4, mom6, hycom, roms

# Experiment name
name=test_$model
#name=TESTc3

# Build directory
BDIR=$CDIR/letkf_build/$name.build
mkdir -p $BDIR
cd $BDIR

# Executable for letkf
PGM=letkf.$name.x

OMP=
#BLAS=1 #0: no blas 1: using blas
sh $CDIR/ulnkcommon.sh
sh $CDIR/lnkcommon.sh $model $CDIR/../
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

$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_mpi.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG common_mtx.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG netlib2.f
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG common_$model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG common_obs_$model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $NETCDF_INC common_mpi_$model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG vars_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG kdtree.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_local.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG letkf_local.o letkf_tools.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_FPP $F90_OBJECT_FLAG letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE -o ${PGM} *.o $MPI_LIB $NETCDF_LIB $LBLAS

#STEVE: keep a record of the build:
mkdir -p CONFIG_$PGM
cp *.f90 CONFIG_$PGM/

rm -f *.mod
rm -f *.o
rm -f netlib2.f
sh $CDIR/ulnkcommon.sh

echo "STEVE: min temp is set to -4 ºC and max salt is set to 50.0 psu in params_model.f90"
echo "STEVE: don't forget - in phys2ijk, obs above model level 1 are set to model level 1"
echo "NORMAL END"
