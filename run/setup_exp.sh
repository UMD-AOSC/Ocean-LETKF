#!/bin/bash
set -ex

name=m2o_sst_letkf
#name=m2o_adt_letkf
#name=o2m_adt
YMDH=2006122800
DO_RESTART=1
DO_ANALYSIS=1

#---

eroot=/ncrc/home1/Steve.Penny/lf1u/OUTPUT
basename=tmp_robs_hybrid_ts
MEMBERS=56

expname=${basename}__$name
EXPSRC=$eroot/$basename
EXPDST=$eroot/$expname

# Create the new experiment directory:
mkdir -p $EXPDST
if [ ! -f tmp_$name ]; then
  ln -fs $EXPDST tmp_$name
fi

# Copy RESTART files for each member:
if [ "$DO_RESTART" -eq "1" ]; then
  MEM=1
  while test $MEM -le $MEMBERS
  do
    MEM2=`printf %.2d $MEM`

    SRCDIR=$EXPSRC/$YMDH/model/$MEM2/RESTART
    DSTDIR=$EXPDST/$YMDH/model/$MEM2/RESTART
    mkdir -p $DSTDIR

#   filelist=(
#             ocean_temp_salt.res.nc
#             ocean_velocity.res.nc
#             ocean_sbc.res.nc
#             ocean_barotropic.res.nc
#            )

#   for file in ${filelist[@]}
    cd $SRCDIR
    for file in *.gz
    do
      #${file%.gz}
      #${file##*/}
      [ -e "$file" ] || continue
      gunzip -c $SRCDIR/${file%.gz} > $DSTDIR/${file%.gz}
    done

    for file in *.nc
    do
      [ -e "$file" ] || continue
      ln -f $SRCDIR/$file $DSTDIR/
    done

    for file in *.res
    do
      [ -e "$file" ] || continue
      cp $file $DSTDIR/
    done

    MEM=`expr $MEM + 1`
  done
fi

# Copy LETKF and hybrid analysis files:
if [ "$DO_ANALYSIS" -eq "1" ]; then

  # Get letkf files:
  SRCDIR=$EXPSRC/$YMDH/letkf
  DSTDIR=$EXPDST/$YMDH/letkf
  mkdir -p $DSTDIR
  cd $SRCDIR
  for file in anal*.nc.gz
  do
    [ -e "$file" ] || continue
    gunzip -c $file >  $DSTDIR/${file%.gz}
  done
  for file in anal*.nc
  do
    [ -e "$file" ] || continue
    cp $file $DSTDIR/$file
  done

  # Get hybrid files:
  SRCDIR=$EXPSRC/$YMDH/hybrid
  DSTDIR=$EXPDST/$YMDH/hybrid
  mkdir -p $DSTDIR
  cd $SRCDIR
  for file in anal*.nc.gz
  do
    [ -e "$file" ] || continue
    gunzip -c $file >  $DSTDIR/${file%.gz}
  done
fi

# Get INIT directory, and change permissions to update input.nml files:
SRCDIR=$EXPSRC/INIT
DSTDIR=$EXPDST/INIT
rm -rf $DSTDIR
mkdir -p $DSTDIR
cp -rl $SRCDIR/* $DSTDIR/
mv $DSTDIR/INPUT/input.nml $DSTDIR/INPUT/input.nml.orig
mv $DSTDIR/letkf/input.nml $DSTDIR/letkf/input.nml.orig
#mkdir -p $DSTDIR/INPUT
#mkdir -p $DSTDIR/letkf
chmod -w $DSTDIR/INPUT/*
chmod +w $DSTDIR/INPUT
chmod +w $DSTDIR/letkf
cp $DSTDIR/INPUT/input.nml.orig $DSTDIR/INPUT/input.nml
cp $DSTDIR/letkf/input.nml.orig $DSTDIR/letkf/input.nml
chmod +w $DSTDIR/INPUT/input.nml
chmod +w $DSTDIR/letkf/input.nml

# Get all grads scripts and support files
cp $EXPSRC/*.ctl $EXPDST/
cp $EXPSRC/*.gs  $EXPDST/

# Set up log directory
mkdir -p $EXPDST/log_$name
if [ ! -f log_$name ]; then
  ln -fs $EXPDST/log_$name log_$name
fi
mkdir -p $EXPDST/log_$name/workflow
mkdir -p $EXPDST/log_$name/model_prep
mkdir -p $EXPDST/log_$name/model
mkdir -p $EXPDST/log_$name/letkf_prep
mkdir -p $EXPDST/log_$name/letkf

# Setup to run on first cycle:
mkdir -p $EXPDST/$YMDH/go

# Set mppnccombine executable:
cp -f $EXPSRC/INIT/mppnccombine.ftn $EXPDST/INIT/mppnccombine.ftn

echo " FINAL STEPS: setup submit_job or cronjob that calls rocoto "
echo "              specify experiment parameters in appropriate xml/run.xml script "
echo "              specify experiment parameters in INIT/letkf/input.nml "
