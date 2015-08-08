#!/bin/sh
# for deleting links of common sources
set -e

flist='SFMT.f90 common.f90 common_mpi.f90 common_mtx.f90 common_letkf.f90 netlib.f netlibblas.f common_mom4.f90 common_mpi_mom4.f90 common_obs_mom4.f90 letkf_local.f90 letkf_obs.f90 letkf_tools.f90 obsop.f90 params_obs.f90 vars_obs.f90 params_model.f90 vars_model.f90'

for f in $flist
do
  # For safety, make sure each file is a symbolic link before deleting:
  if test -h "$f"; then
    rm -f $f
  fi
done

exit 0

rm -f SFMT.f90
rm -f common.f90
rm -f common_mpi.f90
rm -f common_mtx.f90
rm -f common_letkf.f90
rm -f netlib.f
rm -f netlibblas.f

rm -f common_mom4.f90
rm -f common_mpi_mom4.f90
rm -f common_obs_mom4.f90

rm -f letkf.f90
rm -f letkf_local.f90
rm -f letkf_obs.f90
rm -f letkf_tools.f90
rm -f obsop.f90
