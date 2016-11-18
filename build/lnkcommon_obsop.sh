#!/bin/sh
# for making link of common sources
set -e

model=$1
root=$2

COMMONDIR=$root/src/common_all
cp $COMMONDIR/SFMT.f90 ./
cp $COMMONDIR/common.f90 ./
cp $COMMONDIR/common_mpi.f90 ./
cp $COMMONDIR/common_mtx.f90 ./
cp $COMMONDIR/common_letkf.f90 ./
cp $COMMONDIR/netlib.f ./
cp $COMMONDIR/netlibblas.f ./
cp $COMMONDIR/kdtree.f90 ./

MODELDIR=$root/src/model_specific/$model
cp $MODELDIR/common_$model.f90 ./
cp $MODELDIR/common_mpi_$model.f90 ./
cp $MODELDIR/common_obs_$model.f90 ./
cp $MODELDIR/params_model.f90 ./
cp $MODELDIR/vars_model.f90 ./
if [ -f $MODELDIR/${model}_io.f90 ]; then
  cp $MODELDIR/${model}_io.f90 ./
fi

LETKFDIR=$root/src/letkf
cp $LETKFDIR/params_letkf.f90 ./
cp $LETKFDIR/letkf_obs.f90 ./
cp $LETKFDIR/letkf_local.f90 ./
cp $LETKFDIR/letkf_tools.f90 ./
cp $LETKFDIR/vars_obs.f90 ./

OBSDIR=$root/src/obs
cp $OBSDIR/params_obs.f90 ./
cp $OBSDIR/compute_profile_error.f90 ./
cp $OBSDIR/read_argo.f90 ./
cp $OBSDIR/read_avhrr_pathfinder.f90 ./
cp $OBSDIR/read_aviso_adt.f90 ./
cp $OBSDIR/obsop_tprof.f90 ./
cp $OBSDIR/obsop_sprof.f90 ./
cp $OBSDIR/obsop_sst.f90 ./
cp $OBSDIR/obsop_adt.f90 ./
cp $OBSDIR/gsw_pot_to_insitu.f90 ./
cp $OBSDIR/read_ice_txt.f90 ./
cp $OBSDIR/obsop_icefrac.f90 ./

# Software package for equation of state computation from TEOS-2010
GSWDIR=$root/src/obs/gsw_fortran_v3_03
cp $GSWDIR/gsw_oceanographic_toolbox.f90 ./
cp $GSWDIR/gsw_data_v3_0.dat ./
