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
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh

# Ensemble size
# STEVE: figure out how to read from params_letkf.f90 and put here (e.g. with awk/perl/etc.)
#        -> grep and sed seem to work ok:
#        (Set the ensemble size in params_letkf.f90, it will read it in here)
MEM=`grep nbv= params_letkf.f90 | sed -r 's/INTEGER,PARAMETER :: nbv=([0-9]+)/\1/'`
echo "MEM=$MEM"
MEM3=`printf %.3d ${MEM}`

name=DRIFTERS
PGM=obsop.$name.$MEM3

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

$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG common_mpi.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG common_mtx.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $F90_OBJECT_FLAG netlib2.f
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG common_letkf.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $NETCDF_INC $F90_OBJECT_FLAG common_mom4.f90 #$NETCDF_LIB
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG $NETCDF_INC common_obs_mom4.f90 #$NETCDF_LIB		
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG $NETCDF_INC common_mpi_mom4.f90 #$NETCDF_LIB
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG letkf_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG letkf_drifters_local.f90
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG letkf_local.f90
#$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG letkf_local.o letkf_tools.f90
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG $NETCDF_INC letkf_drifters_local.o letkf_drifters_tools.f90 # $NETCDF_LIB
#--
#$F90 $OMP $F90_OPT $F90_DEBUG $F90_F90_OBJECT_FLAG gsw_oceanographic_toolbox.f90
#$F90 $OMP $F90_OPT $F90_DEBUG $F90_F90_OBJECT_FLAG gsw_pot_to_insitu.f90
#--
$F90 $OMP $F90OPT $DEBUG_OPT $F90_OBJECT_FLAG $NETCDF_INC letkf_drifters_tools.o obsop_drifters.f90  # $NETCDF_LIB
$F90 $OMP $F90_OPT $NETCDF_INC -o ${PGM} *.o $NETCDF_LIB

echo "finish keeping the record."

rm -rf *.mod
rm -rf *.o
rm -rf netlib2.f
sh ulnkcommon.sh

echo "NORMAL END"
