#!/bin/bash

# (Below) Probably needs to be changed for each new experiment --------------------------------------

GLOBAL_EXPNAME=expt_75.1m #000
#GLOBAL_SCRATCH=/home/spenny/Research/TESTCASE
GLOBAL_SCRATCH=/u/gens/penny/TESTCASE

GLOBAL_OBS_FILE_PREFIX=coda.

GLOBAL_OBS_FILE_SUFFIX=.nest01.
GLOBAL_OBS_FILE_SUFFIX_qc=

GLOBAL_NX=4500
GLOBAL_NY=3298
GLOBAL_SG_NX=250        # (Some) factors of 4500 :: 1500, 900, 500, 300, 250, 225, 180, 150  (note - 500 caused seg fault)
GLOBAL_SG_NY=194        # Factors of 3298 :: 1649, 194, 97, 34, 17, 2
GLOBAL_NTILES=308       # ($GLOBAL_NX/$GLOBAL_SG_NX) * ($GLOBAL_NY/$GLOBAL_SG_NY)
                        
GLOBAL_MEM_START=0
GLOBAL_MEMBERS=4 #1 #5 # inclusive

# Forecast date
GLOBAL_IYYYY=2014
GLOBAL_IMM=02
GLOBAL_IDD=02 #1 #2 #01
GLOBAL_IHH=12

# Forecast end date
GLOBAL_EYYYY=2014
GLOBAL_EMM=02
GLOBAL_EDD=02 #1 #2 #01
GLOBAL_EHH=12

# Grid data date
GLOBAL_GYYYY=2014
GLOBAL_GMM=02
GLOBAL_GDD=02
GLOBAL_GHH=12

GLOBAL_ID1=000000
GLOBAL_ID2=002500
GLOBAL_ID3=00240000

GLOBAL_inc=1
GLOBAL_inc_units='days'

DO_SSH=0  # 1=yes, 0=no

# (Below) Probably only needs to be changed once per new setup --------------------------------------

GLOBAL_FCST_SRCDIR=/u/scrh/smedstad/GLBocn0.08
GLOBAL_OBS_SRCDIR=/u/scrh/smedstad/GLBocn0.08/work2

GLOBAL_LETKF_ROOT=/home/spenny/Research/Ocean-LETKF
GLOBAL_LETKF_NAME=NRL_TORNADO_hycom_nrl

GLOBAL_GRDDIR=$GLOBAL_FCST_SRCDIR/${GLOBAL_EXPNAME}000
GLOBAL_GRDFILE=maskls_pre_000000_000000_1o${GLOBAL_NX}x${GLOBAL_NY}_${GLOBAL_GYYYY}${GLOBAL_GMM}${GLOBAL_GDD}${GLOBAL_GHH}_00000000_datafld

GLOBAL_RUNDIR=/home/spenny/Research/Ocean-LETKF/run/NRL_SETUP/RUN

GRDFILE_LAT=grdlat_sfc_000000_000000_1o${GLOBAL_NX}x${GLOBAL_NY}_${GLOBAL_GYYYY}${GLOBAL_GMM}${GLOBAL_GDD}${GLOBAL_GHH}_00000000_datafld
GRDFILE_LON=grdlon_sfc_000000_000000_1o${GLOBAL_NX}x${GLOBAL_NY}_${GLOBAL_GYYYY}${GLOBAL_GMM}${GLOBAL_GDD}${GLOBAL_GHH}_00000000_datafld
GRDFILE_LEV=grdlvl_pre_000000_000000_1o${GLOBAL_NX}x${GLOBAL_NY}_${GLOBAL_GYYYY}${GLOBAL_GMM}${GLOBAL_GDD}${GLOBAL_GHH}_00000000_datafld
GRDFILE_KMT=maskls_pre_000000_000000_1o${GLOBAL_NX}x${GLOBAL_NY}_${GLOBAL_GYYYY}${GLOBAL_GMM}${GLOBAL_GDD}${GLOBAL_GHH}_00000000_datafld

GLOBAL_MERGE_EXE=${GLOBAL_RUNDIR}/mt.x
GLOBAL_LSCHECK_EXE=${GLOBAL_RUNDIR}/cfo.x

# (Below) Probably doesn't need to be changed --------------------------------------

GLOBAL_RSRT_FILE_SUFFIX=.dat

# Use a qc feedback file to QC the original feedback file by filtering out qc'd obs from NCODA
# 1 = true, 0 = false
GLOBAL_DO_QC=0
# QC reference file (control member) to use for quality control
GLOBAL_QC_REF=

GLOBAL_DO_WRITE_TILE=1

