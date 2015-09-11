#!/bin/bash

root=../..
model=mom4

MAIN=letkf.f90

source ../../config/machine.sh
source ../../config/$MACHINE.fortran.sh
source ../../config/$MACHINE.netcdf.sh
source ../../config/$MACHINE.mpi.sh

export FMKMF_F90="$F90 $NETCDF_INC"
export FMKMF_LINKOPTS="$NETCDF_LIB"
export FMKMF_SPATH=".:$root/$model/common:$root/$model/obs:$root/common"
#export FMKMF_SFTAG=

echo "FMKMF_F90 = $FMKMF_F90"
echo "FMKMF_LINKOPTS = $FMKMF_LINKOPTS"
echo "FMKMF_SPATH = $FMKMF_SPATH"

echo "Running fmkmf ..."
echo "NETCDF directory: $NETCDF_DIR"
echo "NETCDF include files: $NETCDF_INC"
echo ""

#perl fmkmf.pl $MAIN > Makefile
perl fmkmf.pl $MAIN > Makefile
