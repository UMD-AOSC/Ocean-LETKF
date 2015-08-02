#!/bin/sh
# for making link of common sources
set -e
COMMONDIR=../../common
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/netlibblas.f ./

ln -fs ../common/common_mom4.f90 ./
ln -fs ../common/common_mpi_mom4.f90 ./
ln -fs ../common/common_obs_mom4.f90 ./

ln -fs ../obs/params_obs.f90
ln -fs ../obs/vars_obs.f90

# for drifters (DRIFTERS)
#ln -fs ../letkf_drifters/letkf_drifters_tools.f90 ./


