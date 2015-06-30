#!/bin/bash --login
set -ex

PGM=p2l_synth.x
F90=ftn
OMP=
#F90OPT='-byteswapio -tp amd64 -fast -O3'
F90OPT='-ip -mcmodel=medium -shared-intel -fast -O2' #'-O3 -xSSE4.2 -ip -mcmodel=medium -shared-intel'
INLINE=   #'-Q -qipa=inline'
NETCDF=/opt/cray/netcdf/4.2.0/intel/120
INCLUDE_NETCDF="-I${NETCDF}/include"
LIB_NETCDF="-L${NETCDF}/lib -lnetcdf -lnetcdff"
#LIB_NETCDF="-L${NETCDF}/lib -lnetcdf"

sh ulnkcommon.sh
sh lnkcommon.sh
rm -f *.mod
rm -f *.o

#$F90 $OMP $F90OPT $INLINE -c SFMT.f90
#$F90 $OMP $F90OPT $INLINE -c common.f90
#$F90 $OMP $F90OPT $INLINE -Kieee -c common.o isa.f90
##$F90 $OMP $F90OPT $INLINE -c common_obs.f90
#$F90 $OMP $F90OPT $INLINE $INCLUDE_NETCDF -c isa.o common_mom4.f90
#$F90 $OMP $F90OPT $INLINE -c common_obs_mom4.f90
##$F90 $OMP $F90OPT $INLINE -c calendar.f90
##$F90 $OMP $F90OPT $INCLUDE_NETCDF -c common.o common_obs_mom4.o prof2letkf.f90
$F90 $OMP $F90OPT $INCLUDE_NETCDF -c prof2letkf_Synth.f90
$F90 $OMP $F90OPT $INLINE -o ${PGM} *.o $LIB_NETCDF

rm -f *.mod
rm -f *.o
sh ulnkcommon.sh

echo "NORMAL END"
