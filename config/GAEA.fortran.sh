#!/bin/bash

F90=ftn
F90_OPT='-O3 -parallel -what'
#STEVE: -mcmodel=medium needed for large model grid sizes (e.g. higher than 1 degree resolution of om3_core3)
# explanation of -mcmodel=medium and -shared-intel: http://software.intel.com/en-us/forums/showthread.php?t=43717#18089
F90_INLINE= #"-Q -qinline"
F90_DEBUG= #'-g -qfullpath -v -C -qsigtrap=xl__trcedump' # -qflttrap=en:nanq -qsigtrap'
F90_IEEE='-fltconsistency' #'-Kieee' #'-fltconsistency'
F90_OBJECT_FLAG='-c' #STEVE: for some reason, mpxlf doesn't use -c, but rather -g
BLAS=1

#F90_OPT='-ip -mcmodel=medium -shared-intel -fast -O2'
#F90_OPT='-ftz -ip -ipo -O2 -parallel -i_dynamic -what -fpp -fno-alias -stack_temps -safe_cray_ptr -fast'
#F90_OPT='-byteswapio -tp amd64 -fast -O3'
#F90_OPT='-ip -mcmodel=medium -shared-intel -fast -O2' #'-O3 -xSSE4.2 -ip -mcmodel=medium -shared-intel'
