#!/bin/sh
# for making link of common sources
set -e
COMMONDIR=../common_all
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/netlibblas.f ./
ln -fs $COMMONDIR/kdtree.f90 ./

model='mom4'
MODELDIR=../model_specific/$model
ln -fs $MODELDIR/common_$model.f90 ./
ln -fs $MODELDIR/common_mpi_$model.f90 ./
ln -fs $MODELDIR/common_obs_$model.f90 ./
ln -fs $MODELDIR/params_model.f90 ./

OBSDIR=../obs
ln -fs $OBSDIR/params_obs.f90

# for drifters (DRIFTERS)
#ln -fs ../letkf_drifters/letkf_drifters_tools.f90 ./
