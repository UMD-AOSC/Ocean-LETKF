#!/bin/sh
#        /usr/bin/modulecmd tcsh !*
#        
#        . /etc/profile
#        if your are a sh/bash/ksh user, and
#
#        source /etc/csh.cshrc
#        in case your script is csh derived.
set -ex

# Load modules
#. /etc/profile
#module load intel
#module load mpt
#module load netcdf/4.1.3-intel

#MEM=008
MEM=056
#MEM=028
#MEM=002
#MEM=004
#name=CPO_DEBUG
#name=CPO_SOLO
#name=day1analysis
name=CPO_ALT
#name=DRIFTERS
#name=NCEP_SFCFLUX
#name=NCEP_TEST
PGM=letkf.$name.$MEM
AGM=aoerl.$name.$MEM

F90=ftn #ifort
F90s=ftn #ifort #STEVE: in case we need a different compiler for serial runs
OMP=
PWD=`pwd`
F90OPT='-O3 -parallel -what'
#STEVE: -mcmodel=medium needed for large model grid sizes (e.g. higher than 1 degree resolution of om3_core3)
# explanation of -mcmodel=medium and -shared-intel: http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
INLINE= #"-Q -qinline"
DEBUG_OPT= #'-g -qfullpath -v -C -qsigtrap=xl__trcedump' # -qflttrap=en:nanq -qsigtrap'
BLAS=1 #0: no blas 1: using blas
#STEVE: ask about blas on zeus
INCLUDE_NETCDF="-I${NETCDF}/include"
LIB_NETCDF="-L${NETCDF}/lib -lnetcdf -lnetcdff"
INCLUDE_MPI= 
LIB_MPI= #"-lmpi"
sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

#cat netlib.f > netlib2.f90
cat netlib.f > netlib2.f
if test $BLAS -eq 1
then
  LBLAS="-L${CRAY_LIBSCI_PREFIX_DIR}/lib -lsci_intel -lsci_intel_mp"
else
  cat netlibblas.f >> netlib2.f
  LBLAS=""
fi

IEEE='-fltconsistency' #'-Kieee' #'-fltconsistency'
OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG SFMT.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG common.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $IEEE $OBJECT_FLAG common.o isa.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG common_mpi.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG common_mtx.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $OBJECT_FLAG netlib2.f
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG common_letkf.f90
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE $INCLUDE_NETCDF $OBJECT_FLAG common_mom4.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG common_obs_mom4.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG $INCLUDE_NETCDF common_mpi_mom4.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_obs.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_local.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf_local.o letkf_tools.f90
$F90 $OMP $F90OPT $DEBUG_OPT $OBJECT_FLAG letkf.f90
rm isa.o
$F90 $OMP $F90OPT $DEBUG_OPT $INLINE -o ${PGM} *.o $LIB_MPI $LIB_NETCDF $LBLAS
#$F90 $OMP $F90OPT $DEBUG_OPT $INLINE -o ${PGM} $OFILES $LIB_MPI $LIB_NETCDF $LBLAS

OFILES='SFMT.o netlib2.o common.o common_mtx.o common_mom4.o common_mpi.o common_mpi_mom4.o common_obs_mom4.o common_letkf.o letkf_obs.o letkf_local.o letkf_tools.o'

#STEVE: keep a record of the build:
mkdir -p CONFIG_$PGM
cp *.f90 CONFIG_$PGM/

rm -f *.mod
rm -f *.o
rm -f netlib2.f
sh ulnkcommon.sh

echo "STEVE: min temp is set to -4 ºC and max salt is set to 50.0 psu incommon_mom4:: write_grd4"
echo "STEVE: don't forget - in phys2ijk, obs above model level 1 are set to model level 1"
echo "STEVE: reset line 1028 to .true. in ../common/common_mom4.f90 :: write_grd4"
echo "NORMAL END"
