#!/bin/bash

#OBS_FILE_PREFIX=0001_nrt_5_2_
# (Usually preceeds date)
OBS_FILE_PREFIX=
#OBS_FILE_SUFFIX=_profb_01_fdbk.nc
OBS_FILE_SUFFIX=profb_01_fdbk_00.nc
OBS_FILE_SUFFIX_qc=profbqc_01_fdbk_00.nc

RESTART_FILE_PREFIX=assim_background_state_DI_
# Here, preceeds the 4-digit zero-padded subgrid id number ????
RESTART_FILE_SUFFIX=.nc

GLOBAL_LETKF_FILE_SUFFIX=.restart.nc

GLOBAL_SRC=/home/rd/dasp/perm
GLOBAL_SCRATCH=/lus/snx11064/scratch/rd/dasp
GLOBAL_EXPNAME=TEST_preMPPNCCOMBINE

GLOBAL_MEM_START=0
GLOBAL_MEMBERS=20 # non-inclusive
GLOBAL_subgrid_start=0
GLOBAL_num_subgrids=360 # non-inclusive

GLOBAL_IYYYY=2010
GLOBAL_IMM=10
GLOBAL_IDD=01

GLOBAL_EYYYY=2010
GLOBAL_EMM=10
GLOBAL_EDD=02

GLOBAL_inc=1
GLOBAL_inc_units='days'

# Use a qc feedback file to QC the original feedback file by filtering out qc'd obs from NEMO
# 1 = true, 0 = false
GLOBAL_DO_QC=1
# QC reference file (control member) to use for quality control
#GLOBAL_QC_REF=$GLOBAL_SCRATCH/$GLOBAL_EXPNAME/DATA/OBS/00/${OBS_FILE_PREFIX}${GLOBAL_IYYYY}${GLOBAL_IMM}${GLOBAL_IDD}_${GLOBAL_EYYYY}${GLOBAL_EMM}${GLOBAL_EDD}${OBS_FILE_SUFFIX_qc}
GLOBAL_QC_REF=$GLOBAL_SCRATCH/$GLOBAL_EXPNAME/DATA/OBS/00/${OBS_FILE_PREFIX}${OBS_FILE_SUFFIX_qc}

GLOBAL_LETKF_ROOT=/home/rd/dasp/perm/Ocean-LETKF
GLOBAL_LETKF_NAME=ECMWF_nemo.subgrid

GLOBAL_GRDDIR=/home/rd/dasp/perm/DATA/STATIC
GLOBAL_GRDFILE=mesh_mask.nc

GLOBAL_RUNDIR=/home/rd/dasp/perm/DATA/RUN

GLOBAL_TPROF_OBS_NAME=POTM_OBS
GLOBAL_TPROF_Hxb_NAME=POTM_Hx
GLOBAL_SPROF_OBS_NAME=PSAL_OBS
GLOBAL_SPROF_Hxb_NAME=PSAL_Hx

