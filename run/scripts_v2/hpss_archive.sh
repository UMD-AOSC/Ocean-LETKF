#!/bin/bash --login
#PBS -d /lustre/f1/unswept/Steve.Penny/
#PBS -l partition=es,size=1,walltime=02:00:00
#PBS -q rdtn
#PBS -l qos=windfall
#PBS -S /bin/sh                                                        #Do not change this - it keeps your job from issuing a false alarm
#PBS -E                                                                #Do not change this - it gives your job more and more useful Moab environment variables
#PBS -A cpo_orr
#PBS -N archive
#PBS -j oe

set -e

module load hsi
#date=/bin/date
date=/lustre/f1/unswept/Steve.Penny/date
DO_GZIP=1

echo "HPSS Archiving Step"
echo "processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"

echo "Checking date function..."
`ls $date`
`which date`

#HPSSDIR=/NCEPDEV/hpssuser/g01/wx23sgp
#OUTDIR=/lustre/f1/unswept/Steve.Penny/OUTPUT
#TMPDIR=/lustre/f1/Steve.Penny
#EXPNAME=tmp_hybrid_ts
#SDIR=$OUTDIR/$EXPNAME/STORE

# STEVE: for the 00 member, assume this is a 'control' run and use the mean forcing fields only:
MEM=$((10#$MEMBERID))
MEM2=`printf %.2d ${MEM}`
MEM3=`printf %.3d ${MEM}`

IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

mkdir -p $TMPDIR/work_archive
TEMP=$TMPDIR/work_archive/$MEM2
mkdir -p $TEMP

inc=$days
inc_units='days'

hsi mkdir -p $EXPNAME
hsi ls

HDIR=$HPSSDIR/$EXPNAME

######################
# MAIN ###############
######################
# MEM=1
# while test $MEM -le $MEMBERS
# do
#   MEM2=`printf %.2d ${MEM}`

    echo "RUNNING... $IY$IM$ID :: MEMBER: $MEM"

    # FIRST, store the analysis restart files
    ADIR=$SDIR/anal/$MEM2
#   HDIR=$HPSSDIR/$EXPNAME/anal/$MEM2
    HDIR=$EXPNAME/anal/$MEM2
    hsi "mkdir -p $HDIR"
    hsi "ls -P $HDIR"

    if [ "$DO_GZIP" -eq "1" ]; then
      echo "tar/gzip'ing: $TEMP/$IY$IM$ID.ocean.res.tgz ..."
      cd $ADIR
      tar cvzf $TEMP/$IY$IM$ID.ocean.res.tgz $IY$IM$ID.000000.ocean_temp_salt.res.nc  $IY$IM$ID.000000.ocean_velocity.res.nc $IY$IM$ID.000000.ocean_sbc.res.nc 
      cd $TEMP
      hsi "cput $IY$IM$ID.ocean.res.tgz $HDIR/$IY$IM$ID.ocean.res.tgz" 
      rm -f $TEMP/$IY$IM$ID.ocean.res.tgz
    fi

    hsi "ls -P $HDIR"

    # SECOND, store the forecast / background restart files
    BDIR=$SDIR/fcst/$MEM2
#   HDIR=$HPSSDIR/$EXPNAME/fcst/$MEM2
    HDIR=$EXPNAME/fcst/$MEM2
    hsi "mkdir -p $HDIR"
    hsi "ls -P $HDIR"

    if [ "$DO_GZIP" -eq "1" ]; then
      echo "tar/gzip'ing: $TEMP/$IY$IM$ID.ocean.res.tgz ..."
      cd $BDIR
      tar cvzf $TEMP/$IY$IM$ID.ocean.res.tgz $IY$IM$ID.000000.ocean_temp_salt.res.nc  $IY$IM$ID.000000.ocean_velocity.res.nc $IY$IM$ID.000000.ocean_sbc.res.nc 
      tar cvzf $TEMP/$IY$IM$ID.ocean.res.tgz $IY$IM$ID.000000.ocean_barotropic.res.nc $IY$IM$ID.000000.ocean_bih_friction.res.nc $IY$IM$ID.000000.ocean_con_temp.res.nc $IY$IM$ID.000000.ocean_density.res.nc $IY$IM$ID.000000.ocean_frazil.res.nc $IY$IM$ID.000000.ocean_neutralB.res.nc $IY$IM$ID.000000.ocean_neutral.res.nc $IY$IM$ID.000000.ocean_residency.res.nc $IY$IM$ID.000000.ocean_sbc.res.nc $IY$IM$ID.000000.ocean_sigma_transport.res.nc $IY$IM$ID.000000.ocean_temp_salt.res.nc $IY$IM$ID.000000.ocean_thickness.res.nc $IY$IM$ID.000000.ocean_velocity_advection.res.nc $IY$IM$ID.000000.ocean_velocity.res.nc $IY$IM$ID.000000.ocean_solo.intermediate.res $IY$IM$ID.000000.ocean_solo.res
      cd $TEMP
      hsi "cput $IY$IM$ID.ocean.res.tgz $HDIR/$IY$IM$ID.ocean.res.tgz" 
      rm -f $TEMP/$IY$IM$ID.ocean.res.tgz
    fi

    hsi "ls -P $HDIR"

#   MEM=`expr $MEM + 1`
# done

echo "The archiving step is complete for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > $SDIR/archive_log/archive_${YYYYMMDDHH}.${MEM2}.out

