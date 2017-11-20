#!/bin/bash

#OBS_FILE_PREFIX=0001_nrt_5_2_
OBS_FILE_PREFIX=3253_

#OBS_FILE_SUFFIX=_profbqc_01_fdbk.nc
OBS_FILE_SUFFIX=_profb_01_fdbk.nc
OBS_FILE_SUFFIX_qc=_profbqc_01_fdbk.nc

RESTART_FILE_SUFFIX=.restart.nc

GLOBAL_SCRATCH=/lus/snx11064/scratch/rd/dasp
GLOBAL_EXPNAME=CERA20C_50mem

GLOBAL_MEM_START=0
GLOBAL_MEMBERS=51

GLOBAL_IYYYY=2005
GLOBAL_IMM=10
GLOBAL_IDD=01

GLOBAL_EYYYY=2005
GLOBAL_EMM=10
GLOBAL_EDD=02

GLOBAL_inc=1
GLOBAL_inc_units='days'

# Use a qc feedback file to QC the original feedback file by filtering out qc'd obs from NEMO
# 1 = true, 0 = false
GLOBAL_DO_QC=1
# QC reference file (control member) to use for quality control
GLOBAL_QC_REF=$GLOBAL_SCRATCH/$GLOBAL_EXPNAME/DATA/OBS/00/${OBS_FILE_PREFIX}${GLOBAL_IYYYY}${GLOBAL_IMM}${GLOBAL_IDD}_${GLOBAL_EYYYY}${GLOBAL_EMM}${GLOBAL_EDD}${OBS_FILE_SUFFIX_qc}

GLOBAL_LETKF_ROOT=/home/rd/dasp/perm/Ocean-LETKF
GLOBAL_LETKF_NAME=ECMWF_nemo.nemovarHx0

GLOBAL_GRDDIR=/home/rd/dasp/perm/DATA/STATIC
GLOBAL_GRDFILE=mesh_mask.nc

GLOBAL_RUNDIR=/home/rd/dasp/perm/DATA/RUN
