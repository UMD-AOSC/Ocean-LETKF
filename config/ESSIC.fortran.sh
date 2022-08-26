#!/bin/bash

F90=mpif90
F90s="gfortran"   #used to compile BUFRLIBS
#F90_OPT='-O2 -ffree-form -ffree-line-length-none'
#STEVE: -mcmodel=medium needed for large model grid sizes (e.g. higher than 1 degree resolution of om3_core3)
# explanation of -mcmodel=medium and -shared-intel: http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
#F90_OPT='-O2 -ffree-line-length-none -frecord-marker=4 -finit-local-zero' #CDA: for gfortran version < 10
F90_OPT='-O2 -ffree-line-length-none -frecord-marker=4 -finit-local-zero -fbacktrace -fcheck=bounds' #CDA: for gfortran version < 10
F90_OPT="$F90_OPT -fallow-argument-mismatch" #CDA: for gfortran version >= 10
F90_INLINE=
F90_DEBUG=
F90_IEEE= #'-Kieee' #'-fltconsistency'
F90_OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g
BLAS=0

F90_FPP='-cpp'  # for gfortran

