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

model='mom4'
MODELDIR=../model_specific/$model
ln -fs $MODELDIR/common_$model.f90 ./
ln -fs $MODELDIR/common_mpi_$model.f90 ./
ln -fs $MODELDIR/common_obs_$model.f90 ./

LETKFDIR=../letkf
ln -fs $LETKFDIR/params_letkf.f90 ./
ln -fs $LETKFDIR/params_model.f90 ./
ln -fs $LETKFDIR/vars_model.f90 ./
ln -fs $LETKFDIR/letkf_obs.f90 ./
ln -fs $LETKFDIR/letkf_local.f90 ./
ln -fs $LETKFDIR/letkf_tools.f90 ./
ln -fs $LETKFDIR/vars_obs.f90 ./

#ln -fs ../letkf_drifters/letkf_drifters_local.f90 ./
#ln -fs ../letkf_drifters/letkf_drifters_tools.f90 ./

ln -fs gsw_fortran_v3_03/gsw_oceanographic_toolbox.f90 ./
ln -fs gsw_fortran_v3_03/gsw_data_v3_0.dat ./
