#!/bin/sh
# for making link of common sources
set -e

model=$1
root=$2

COMMONDIR=$root/src/common_all
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/netlibblas.f ./

MODELDIR=$root/src/model_specific/$model
ln -fs $MODELDIR/common_$model.f90 ./
ln -fs $MODELDIR/common_mpi_$model.f90 ./
ln -fs $MODELDIR/common_obs_$model.f90 ./

LETKFDIR=$root/src/letkf
ln -fs $LETKFDIR/params_letkf.f90 ./
ln -fs $LETKFDIR/params_model.f90 ./
ln -fs $LETKFDIR/vars_model.f90 ./
ln -fs $LETKFDIR/letkf_obs.f90 ./
ln -fs $LETKFDIR/letkf_local.f90 ./
ln -fs $LETKFDIR/letkf_tools.f90 ./
ln -fs $LETKFDIR/vars_obs.f90 ./

OBSDIR=$root/src/obs
ln -fs $OBSDIR/params_obs.f90 ./
ln -fs $OBSDIR/compute_profile_error.f90 ./
ln -fs $OBSDIR/read_argo.f90 ./
ln -fs $OBSDIR/read_avhrr_pathfinder.f90 ./
ln -fs $OBSDIR/read_aviso_adt.f90 ./
ln -fs $OBSDIR/obsop_tprof.f90 ./
ln -fs $OBSDIR/obsop_sprof.f90 ./
ln -fs $OBSDIR/obsop_sst.f90 ./
ln -fs $OBSDIR/obsop_adt.f90 ./
ln -fs $OBSDIR/gsw_pot_to_insitu.f90 ./

# Software package for equation of state computation from TEOS-2010
GSWDIR=$root/src/obs/gsw_fortran_v3_03
ln -fs $GSWDIR/gsw_oceanographic_toolbox.f90 ./
ln -fs $GSWDIR/gsw_data_v3_0.dat ./
