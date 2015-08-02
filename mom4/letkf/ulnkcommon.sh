#!/bin/sh
# for deleting links of common sources
set -e

flist='SFMT.f90 common.f90 common_mpi.f90 common_mtx.f90 common_letkf.f90 netlib.f netlibblas.f common_mom4.f90 common_mpi_mom4.f90 common_obs_mom4.f90 params_obs.f90 vars_obs.f90 letkf_drifters_tools.f90' 

for f in $flist
do
  # For safety, make sure each file is a symbolic link before deleting:
  if test -h "$f"; then
    rm -f $f
  fi
done

exit 0
