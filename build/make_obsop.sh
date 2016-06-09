#!/bin/sh --login
#PBS -r y                                                              #This job is restartable
#PBS -S /bin/sh                                                        #Do not change this - it keeps your job from issuing a false alarm
#PBS -E                                                                #Do not change this - it gives your job more and more useful Moab environment variables
# -- Request 120 cores
#PBS -l size=1
# 
# -- Specify a maximum wallclock of 4 hours
#PBS -l walltime=0:01:00
#
# -- Specify under which account a job should run
#PBS -A cpo_orr
#
# -- Set the name of the job, or moab will default to STDIN
#PBS -N cpo_compile_obs
#
# -- Set the queue: debug, batch, novel, bigmem
#PBS -q ldtn
# 
# -- Set the partition (for Gaea)
#PBS -l partition=es
#
# -- Set this as working directory
#PBS -d .
#
# -- Send output and error to same file
#PBS -j oe
#
set -exv

CDIR=`pwd`

#MACHINE="GAEA3"
CONFIGDIR=../config
source $CONFIGDIR/machine.sh
source $CONFIGDIR/$MACHINE.fortran.sh
source $CONFIGDIR/$MACHINE.netcdf.sh
source $CONFIGDIR/$MACHINE.modules_ldtn.sh

# Model name:
model=mom4

# Experiment name:
name=test_$model

# Build directory
BDIR=$CDIR/obsop_build/$name.build
mkdir -p $BDIR
cd $BDIR

sh $CDIR/ulnkcommon.sh
sh $CDIR/lnkcommon_obsop.sh $model $CDIR/../
rm -f *.mod
rm -f *.o

name=test_$model
#name=TESTc3
PGM=obsop.$name

$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG SFMT.f90
$F90 $OMP $F90_OPT $INLINE $F90_OBJECT_FLAG common.f90
$F90 $OMP $F90_OPT $F90_FPP $F90_OBJECT_FLAG params_model.f90
$F90 $OMP $F90_OPT $F90_FPP $F90_OBJECT_FLAG vars_model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_letkf.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_INLINE $NETCDF_INC $F90_OBJECT_FLAG common_$model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG params_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG vars_obs.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG common_obs_$model.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG compute_profile_error.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_argo.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_avhrr_pathfinder.f90
$F90 $OMP $F90_OPT $F90_OBJECT_FLAG read_aviso_adt.f90
#--
# Equation of state for converting from model to obs space is in here:
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_oceanographic_toolbox.f90
$F90 $OMP $F90_OPT $F90_DEBUG $F90_OBJECT_FLAG gsw_pot_to_insitu.f90
#--
$F90 $OMP $F90_OPT obsop_tprof.f90 -o ${PGM}.tprof.x *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_sprof.f90 -o ${PGM}.sprof.x *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_adt.f90   -o ${PGM}.adt.x   *.o $NETCDF_LIB
$F90 $OMP $F90_OPT obsop_sst.f90   -o ${PGM}.sst.x   *.o $NETCDF_LIB

rm -f *.mod
rm -f *.o
sh $CDIR/ulnkcommon.sh
#--
#--

echo "NORMAL END"
