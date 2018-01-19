#!/bin/bash
set -ex 

exe=mt.x

F90_ROOT=/common/pgi/15.7
F90_BIN=$F90_ROOT/linux86-64/15.7/bin
F90_LIB=$F90_ROOT/linux86-64/15.7/lib
LM_LICENSE_FILE=$F90_ROOT/license.dat
MPI_F90_ROOT=/common/openmpi/el7/pgi15.7/1.10.2
MPI_F90_BIN=$MPI_F90_ROOT/bin

F90=$MPI_F90_BIN/mpif90
F90_OPT='-fast -byteswapio -Mextend -tp=nehalem,sandybridge -r4 -i4' # -lpgas-dmapp'
#F90_DEBUG='-C -g -gopt -Mbounds -Mchkfpstk -Mchkptr -Mchkstk -Mcoff -Mdwarf1 -Mdwarf2 -Mdwarf3 -Melf -Mnodwarf -Mpgicoff -traceback'
F90_DEBUG='-Mbounds -traceback'
#F90_OPT='-fast -byteswapio -Mextend -tp=nehalem,sandybridge -r4 -i8' # -lpgas-dmapp'


$F90 $F90_OPT $F90_DEBUG merge_tiles.f90 -o $exe
