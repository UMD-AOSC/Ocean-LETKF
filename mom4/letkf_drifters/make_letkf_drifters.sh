#!/bin/bash
# Special compile script just for UMD's DT2

set f90file = letkf_drifters.f90
set exefile = out.sh
set 

ifort -I$NETCDF_FORTRAN_ROOT/include $f90file common*.f90 letkf_drifters_local.f90 letkf_obs.f90 params_*.f90 -o $exefile -L$NETCDF_FORTRAN_ROOT/lib -lnetcdff -Wl,-rpath,$NETCDF_FORTRAN_ROOT/lib
