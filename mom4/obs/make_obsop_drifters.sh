#!/bin/sh
#        /usr/bin/modulecmd tcsh !*
#        
#        . /etc/profile
#        if your are a sh/bash/ksh user, and
#
#        source /etc/csh.cshrc
#        in case your script is csh derived.
set -exv

# sh make_obsop_drifters.sh $MEM
source ../../config/machine.sh
source ../../config/$MACHINE.netcdf.sh

MEM=$1
if [ ${MEM} -lt 100 ]; then
  MEM=0${MEM}
fi

if [ ${MEM} -lt 10 ]; then
  MEM=0${MEM}
fi 

name=DRIFTERS
PGM=obsop.$name.$MEM

#STEVE: put the directory for netcdf here:
NETCDF_ROOT=/cell_root/software/netcdf/4.3.2/intel/2013.1.039/openmpi/1.8.1/hdf5/1.8.13/hdf4/4.2.10/sys
NETCDFF_ROOT=/cell_root/software/netcdf-fortran/4.4.1/netcdf/4.3.2/intel/2013.1.039/openmpi/1.8.1/sys
LIB_NETCDF="-L$NETCDFF_ROOT/lib -lnetcdff -Wl,-rpath,$NETCDFF_ROOT/lib"
INC_NETCDF="-I$NETCDFF_ROOT/include -I$NETCDF_ROOT/include"
INC=$INC_NETCDF

F90=mpif90
F90s=mpif90 #STEVE: in case we need a different compiler for serial runs
OMP=
PWD=`pwd`
F90OPT='-O3'
# explanation of -mcmodel=medium and -shared-intel: http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
INLINE= #"-Q -qinline"
DEBUG_OPT= #'-g -qfullpath -v -C -qsigtrap=xl__trcedump' # -qflttrap=en:nanq -qsigtrap'
BLAS=0 #0: no blas 1: using blas

INCLUDE_MPI= 
LIB_MPI="-lmpi"
sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

cat netlib.f > netlib2.f
if test $BLAS -eq 1
then
  LBLAS="-L/usr/lib -lblas"
else
  #cat netlibblas.f >> netlib2.f90
  cat netlibblas.f >> netlib2.f
  rm netlibblas.f
  LBLAS=""
fi

OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG SFMT.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG common.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG params_model.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG params_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG vars_model.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG common_mpi.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG common_mtx.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG netlib2.f
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG common_letkf.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $INC_NETCDF $OBJECT_FLAG common_mom4.f90 $LIB_NETCDF
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG $INC_NETCDF common_obs_mom4.f90 $LIB_NETCDF		
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG $INC_NETCDF common_mpi_mom4.f90 $LIB_NETCDF
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_drifters_local.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_local.f90
#$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_local.o letkf_tools.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG $INC_NETCDF letkf_drifters_local.o letkf_drifters_tools.f90 $LIB_NETCDF
#--
#$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_oceanographic_toolbox.f90
#$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_pot_to_insitu.f90
#--
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG $INC_NETCDF letkf_drifters_tools.o obsop_drifters.f90 $LIB_NETCDF
$F90 $OMP $F90_OPT $INC_NETCDF -o ${PGM} *.o $LIB_NETCDF
#$F90 $OMP $F90OPT $DEBUG_OPT -o ${PGM} *.o $LIB_NETCDF $LBLAS $LIB_MPI
#$F90 $OMP $F90OPT $DEBUG_OPT $INLINE -o ${PGM} *.o $LIB_MPI $LIB_NETCDF $LBLAS

echo "finish keeping the record."

rm -rf *.mod
rm -rf *.o
rm -rf netlib2.f
sh ulnkcommon.sh

echo "NORMAL END"
