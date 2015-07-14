#!/bin/bash

F90=mpif90
F90s=mpif90
F90_OPT='-O3'
#STEVE: -mcmodel=medium needed for large model grid sizes (e.g. higher than 1 degree resolution of om3_core3)
# explanation of -mcmodel=medium and -shared-intel: http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
F90_INLINE=
F90_DEBUG=
F90_IEEE= #'-Kieee' #'-fltconsistency'
F90_OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g

