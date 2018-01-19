#!/bin/bash

# (Below) Probably needs to be changed for each new experiment --------------------------------------

# Experiment working directory name:
GLOBAL_EXPNAME=TEST_5memOceanSubgrids

# Specify number of members and subgrids:
GLOBAL_MEM_START=0
GLOBAL_MEMBERS=5        # non-inclusive
GLOBAL_subgrid_start=0
GLOBAL_num_subgrids=480 # non-inclusive

# Specify dates:
GLOBAL_IYYYY=2017
GLOBAL_IMM=01
GLOBAL_IDD=01

GLOBAL_EYYYY=2017
GLOBAL_EMM=01
GLOBAL_EDD=06

GLOBAL_inc=5
GLOBAL_inc_units='days'

# Specify obs data file naming formats:
# (Usually preceeds date)
#OBS_FILE_PREFIX=0001_nrt_5_2_
#OBS_FILE_SUFFIX=_profb_01_fdbk.nc
GLOBAL_OBS_FILE_PREFIX=
#GLOBAL_OBS_FILE_SUFFIX=profb_01_fdbk_00.nc
GLOBAL_OBS_FILE_SUFFIX=profbqc_01_fdbk_00.nc    #STEVE: skipping preprocessing because obs locatiosn aren't perturbed in this run
GLOBAL_OBS_FILE_SUFFIX_qc=profbqc_01_fdbk_00.nc

# Specify restart file naming formats:
GLOBAL_RESTART_FILE_PREFIX=gucq_${GLOBAL_EYYYY}${GLOBAL_EMM}${GLOBAL_EDD}_000000_restart_
# Here, preceeds the 4-digit zero-padded subgrid id number ????
GLOBAL_RESTART_FILE_SUFFIX=.nc

# (Below) Probably only needs to be changed once per new setup --------------------------------------

GLOBAL_PERM=/home/rd/dasp/perm
GLOBAL_SCRATCH=/lus/snx11064/scratch/rd/dasp

GLOBAL_LETKF_ROOT=$GLOBAL_PERM/Ocean-LETKF
GLOBAL_LETKF_NAME=ECMWF_nemo.subgrid

GLOBAL_GRDDIR=$GLOBAL_PERM/DATA/STATIC
GLOBAL_GRDFILE=mesh_mask.nc

GLOBAL_RUNDIR=$GLOBAL_PERM/DATA/RUN_SUBGRIDS

# (Below) Probably doesn't need to be changed --------------------------------------

GLOBAL_LETKF_FILE_SUFFIX=.restart.nc

# Use a qc feedback file to QC the original feedback file by filtering out qc'd obs from NEMO
# 1 = true, 0 = false
GLOBAL_DO_QC=0 #1 (STEVE: the qc preprocessing doesn't work if the observation indices exceed the number of obs in the control file)

# QC reference file (control member) to use for quality control
#GLOBAL_QC_REF=$GLOBAL_SCRATCH/$GLOBAL_EXPNAME/DATA/OBS/00/${OBS_FILE_PREFIX}${GLOBAL_IYYYY}${GLOBAL_IMM}${GLOBAL_IDD}_${GLOBAL_EYYYY}${GLOBAL_EMM}${GLOBAL_EDD}${OBS_FILE_SUFFIX_qc}
GLOBAL_QC_REF=$GLOBAL_SCRATCH/$GLOBAL_EXPNAME/DATA/OBS/00/${GLOBAL_OBS_FILE_PREFIX}${GLOBAL_OBS_FILE_SUFFIX_qc}

# Specify names of obs data and model equivalent in the feedback files being used for this experiment
GLOBAL_TPROF_OBS_NAME=POTM_OBS
GLOBAL_TPROF_Hxb_NAME=POTM_Hx
GLOBAL_SPROF_OBS_NAME=PSAL_OBS
GLOBAL_SPROF_Hxb_NAME=PSAL_Hx
# NOTE: these may need to be changed manually in the python scripts
