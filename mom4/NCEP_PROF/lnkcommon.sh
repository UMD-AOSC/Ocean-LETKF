#!/bin/sh
# for making link of common sources
set -e
COMMONDIR=/gpfs/v/cdev/save/wx23sgp/letkf-penny/common
MOM4DIR=/gpfs/v/cdev/save/wx23sgp/letkf-penny/mom4
ln -fs $COMMONDIR/SFMT.f90 ./
ln -fs $COMMONDIR/common.f90 ./
#ln -fs $COMMONDIR/common_mpi.f90 ./
ln -fs $COMMONDIR/common_obs.f90 ./
#ln -fs $COMMONDIR/common_mtx.f90 ./
#ln -fs $COMMONDIR/common_letkf.f90 ./
#ln -fs $COMMONDIR/netlib.f ./

ln -fs $MOM4DIR/common/common_mom4.f90 ./
#ln -fs ../../common/common_mpi_mom4.f90 ./
ln -fs $MOM4DIR/common/common_obs_mom4.f90 ./

ln -fs $MOM4DIR/letkf/isa.f90 ./
