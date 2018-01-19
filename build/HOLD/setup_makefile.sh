#!/bin/bash
set -exv

root=../..
model=mom4

MAIN=letkf.x

source ../../config/machine.sh
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh
source ../../config/$MACHINE.mpi.sh

# Source directories, other than the local directory:
C1=../../common
C2=../common
O1=../obs
L1=../letkf

#SPATH="*.f90 $root/$model/common/*.f90 $root/$model/obs/*.f90 $root/common/*.f90"

# Ensemble size
# STEVE: figure out how to read from params_letkf.f90 and put here (e.g. with awk/perl/etc.)
#        -> grep and sed seem to work ok:
#        (Set the ensemble size in params_letkf.f90, it will read it in here)
MEM=`grep nbv= $L1/params_letkf.f90 | sed -r 's/INTEGER,PARAMETER :: nbv=([0-9]+)/\1/'`
echo "MEM=$MEM"
MEM3=`printf %.3d ${MEM}`

# Experiment name
name=MAKE
# Executable for letkf
PGM=letkf.$name.$MEM3

if [ ! -f netlib2.f ]; then
  cat $C1/netlib.f > netlib2.f
fi
if test $BLAS -eq 1
then
  LBLAS="-L${CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_intel -lsci_intel_mp"
else
  cat netlibblas.f >> netlib2.f
  LBLAS=""
fi

CC=icc
TEMPLATE=template
cat <<EOM >$TEMPLATE
FC=$F90
CC=$CC
LD=$F90
LDFLAGS=$F90_OPT $NETCDF_LIB $LBLAS
EOM

FILELIST=filelist
cat <<EOM >$FILELIST
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG $C1/SFMT.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG $C1/common.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $C1/common_mpi.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG $C1/common_mtx.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG netlib2.f
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $F90_OBJECT_FLAG $L1/params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $C1/common_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/params_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/vars_model.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG $C2/common_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $O1/params_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $O1/vars_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $C2/common_obs_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $NETCDF_INC $C2/common_mpi_mom4.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/letkf_obs.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/vars_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/letkf_local.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/letkf_tools.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG $L1/letkf.f90
EOM
#$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE -o ${PGM} *.o $MPI_LIB $NETCDF_LIB $LBLAS

echo "Running mkmf ..."
echo "NETCDF directory: $NETCDF_DIR"
echo "NETCDF include files: $NETCDF_INC"
echo "NETCDF library files: $NETCDF_LIB"
echo ""

#perl mkmf -t $TEMPLATE -o $NETCDF_INC -p ${PGM} $FILELIST
perl mkmf -t $TEMPLATE -p ${PGM} $FILELIST
