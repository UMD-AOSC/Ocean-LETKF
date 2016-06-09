#!/bin/sh
# for making link of common sources
set -e

model=$1 #input as 1st script argument, e.g. 'mom4'
root=$2  #input as 2nd script argument

# Link common source files
COMMONDIR=$root/src/common_all
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_mtx.f90 ./
ln -fs $COMMONDIR/common_letkf.f90 ./
ln -fs $COMMONDIR/netlib.f ./
ln -fs $COMMONDIR/netlibblas.f ./
ln -fs $COMMONDIR/kdtree.f90 ./

# Link model-specific source files
MODELDIR=$root/src/model_specific/$model
ln -fs $MODELDIR/common_$model.f90 ./
ln -fs $MODELDIR/common_mpi_$model.f90 ./
ln -fs $MODELDIR/common_obs_$model.f90 ./
ln -fs $MODELDIR/params_model.f90 ./

# Link letkf source files
LDIR=$root/src/letkf/
ln -fs $LDIR/letkf.f90 ./
ln -fs $LDIR/letkf_tools.f90 ./
ln -fs $LDIR/letkf_local.f90 ./
ln -fs $LDIR/letkf_obs.f90 ./
ln -fs $LDIR/params_letkf.f90 ./
ln -fs $LDIR/vars_letkf.f90 ./
ln -fs $LDIR/vars_obs.f90 ./
ln -fs $LDIR/vars_model.f90 ./

OBSDIR=$root/src/obs
ln -fs $OBSDIR/params_obs.f90

