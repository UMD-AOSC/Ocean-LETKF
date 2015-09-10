#!/bin/sh
# Special compile script just for UMD's DT2

f90file=create_drifters.f90
exefile=cd.x

ifort -I$NETCDF_FORTRAN_ROOT/include $f90file -o $exefile -L$NETCDF_FORTRAN_ROOT/lib -lnetcdff -Wl,-rpath,$NETCDF_FORTRAN_ROOT/lib

# Generate "drifters_inp.nc"
./$exefile



