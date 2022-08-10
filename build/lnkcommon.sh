#!/bin/sh
# for making link of common sources
set -e

model=$1 #input as 1st script argument, e.g. 'mom4'
root=$2  #input as 2nd script argument

# Link common source files
COMMONDIR=$root/src/common_all
cp $COMMONDIR/SFMT.f90 ./
cp $COMMONDIR/common.f90 ./
cp $COMMONDIR/common_mpi.f90 ./
cp $COMMONDIR/common_mtx.f90 ./
cp $COMMONDIR/common_letkf.f90 ./
cp $COMMONDIR/netlib.f ./
cp $COMMONDIR/netlibblas.f ./
cp $COMMONDIR/kdtree.f90 ./

# Link model-specific source files
MODELDIR=$root/src/model_specific/$model
cp $MODELDIR/common_$model.f90 ./
cp $MODELDIR/common_mpi_$model.f90 ./
cp $MODELDIR/common_obs_$model.f90 ./
cp $MODELDIR/params_model.f90 ./
cp $MODELDIR/vars_model.f90 ./
if [ -f $MODELDIR/${model}_io.f90 ]; then
  cp $MODELDIR/${model}_io.f90 ./
fi
cp $MODELDIR/input_nml_${model}.f90 ./
cp $MODELDIR/common_debug_${model}.f90 ./

# Link letkf source files
LDIR=$root/src/letkf/
cp $LDIR/letkf.f90 ./
cp $LDIR/letkf_tools.f90 ./
cp $LDIR/letkf_local.f90 ./
cp $LDIR/letkf_obs.f90 ./
cp $LDIR/params_letkf.f90 ./
cp $LDIR/vars_letkf.f90 ./
cp $LDIR/vars_obs.f90 ./

#if [ ${model} == "hycom" ] || [ ${model} == "hycom_nrl" ]; then
if [ ${model} == "hycom" ]; then
   cp $MODELDIR/mod_xc.F ./
   cp $MODELDIR/mod_za.F ./
   cp $MODELDIR/mod_ppsw.F ./
   cp $MODELDIR/wtime.F ./
   cp $MODELDIR/hycom_intrp.f ./
   cp $MODELDIR/layer2z.f ./
else
  echo "in lnkcommon.sh: model = $model"
fi

if [ ${model} == "hycom_nrl" ]; then
   echo "in lnkcommon.sh: model = $model: copy additional files"
   cp $MODELDIR/tools/test_io_subgrid.f90 .
   cp $MODELDIR/tools/test_read_bin.f90 .
fi 

OBSDIR=$root/src/obs
cp $OBSDIR/params_obs.f90 ./
