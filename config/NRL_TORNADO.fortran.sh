#!/bin/bash

F90_ROOT=/common/pgi/15.7
F90_BIN=$F90_ROOT/linux86-64/15.7/bin
F90_LIB=$F90_ROOT/linux86-64/15.7/lib
LM_LICENSE_FILE=$F90_ROOT/license.dat
MPI_F90_ROOT=/common/openmpi/el7/pgi15.7/1.10.2
MPI_F90_BIN=$MPI_F90_ROOT/bin

F90=$MPI_F90_BIN/mpif90
F90_OPT='-fast -byteswapio -Mextend -tp=nehalem,sandybridge -r4 -i4' # -lpgas-dmapp'
F90_INLINE= #"-Q -qinline"
F90_DEBUG='-Mbounds -traceback'
#F90_DEBUG=
F90_IEEE='-Kieee' #'-Kieee' #'-fltconsistency'
F90_OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g
BLAS=1
# If BLAS=1, then must specify machine's BLAS location and compiler options:
BLAS_DIR=$F90_LIB
BLAS_LIBS="-llapack -lblas"
# For preprocessor to be applied to fortran compile:
F90_FPP='-Mpreprocess'   # pgi
#F90_FPP='-fpp'          # intel
#F90_FPP='-eZ'           # cray

