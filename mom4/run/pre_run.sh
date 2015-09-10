#!/bin/sh

# This code is to generate the executive file that is necessary for the whole program.

# sh pre_run.sh <root> <ENS_NUM>

root=$1
ENS_NUM=$2

###########################################################################
#   Generate OBSOP
###########################################################################
cd ${root}/obs
sh make_obsop.sh 1 		# generate obsop.001
sh make_obsop_drifters.sh 1     # generate obsop.DRIFTERS.001

###########################################################################
#   Generate LETKF.DRIFTERS
###########################################################################
cd ${root}/letkf_drifters
sh make_drifters.DT2.sh $ENS_NUM   # generate letkf.DRIFTERS.$ENS_NUM

exit 0


