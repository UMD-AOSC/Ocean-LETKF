#!/bin/bash

F90=ftn
F90OPT='O3 -parallel -what'
#F90OPT='-ip -mcmodel=medium -shared-intel -fast -O2'
#F90OPT='-ftz -ip -ipo -O2 -parallel -i_dynamic -what -fpp -fno-alias -stack_temps -safe_cray_ptr -fast'
#F90OPT='-byteswapio -tp amd64 -fast -O3'
#F90OPT='-ip -mcmodel=medium -shared-intel -fast -O2' #'-O3 -xSSE4.2 -ip -mcmodel=medium -shared-intel'
