#!/bin/sh
# for deleting links of common sources
set -e

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
