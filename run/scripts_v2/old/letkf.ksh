#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
#module load nco

echo "LETKF run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf
mkdir -p ${workdir}
cd ${workdir}

#STEVE: active code:
USE_INFLADJ=1
USE_ADAPOBS=0
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the date for the next analysis cycle
date=/bin/date
inc=5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

# Update the date for the previous analysis cycle
pinc=5
pinc_units='days ago'
PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
PN=$IN
PS=$IS
workdir0=${EXP_DATA}/$PY$PM$PD$PH/letkf

### Processs inputs

# Copy the inflation (copy instead of link, because letkf rewrites over the file)
#if test -f $OUTPUT/infl_mul/$IY$IM$ID$IH.grd
if [ $USE_INFLADJ -a -f $workdir0/infl_redux.grd ]; then
  cp $workdir0/infl_redux.grd infl_mul.grd  
else
  if test -f $workdir0/infl_mul.grd
  then
    cp $workdir0/infl_mul.grd infl_mul.grd
  fi
fi

# Also, if there is adaptive obervation error, copy that file as well...
if test -f $workdir0/adapt_oer.grd
then
  cp $workdir0/adapt_oer.grd adapt_inp.grd
fi

# link the (0-360 degree spherical) topography for letkf to read grid information
# (This is needed to read and write mom4p1 netcdf files in letkf)
ln -fs ../../../regrid/grid_spec.nc grid_spec.nc

# START RUN LETKF ################################################
export MPI_VERBOSE=1
export MPI_DISPLAY_SETTINGS=1
export MPI_BUFS_PER_PROC=128
export MPI_BUFS_PER_HOST=128
export MPI_IB_RAILS=2
export MPI_GROUP_MAX=128
mpiexec_mpt -np $PBS_NP $LETKFexe
echo "This is the LETKF run for cycle ${YYYYMMDDHH}" > letkf.out
# END   RUN LETKF ################################################


### process outputs
if [ -s $TMPDIR/letkf/NOUT-000 ]; then
  mv NOUT-000 $OUTPUT/log/$IY$IM$ID$IH.log
fi
if test -f infl_mul.grd
then
  if [ $USE_INFLADJ -a ! -f infl_redux.grd ]; then  #STEVE: if inflation relaxation is needed, it is applied here
    echo "Using inflation adjust..."
    $INFLadj
#   cp -f infl_redux.grd $OUTPUT/infl_mul/$TY$TM$TD$TH.grd
#   cp -f infl_mul.grd $OUTPUT/infl_out/$TY$TM$TD$TH.grd
# else
#   cp infl_mul.grd $OUTPUT/infl_mul/$TY$TM$TD$TH.grd
  fi
fi

# Get the files for adaptive obs
if test -f adap_obserr.grd
then
  if [ $USE_ADAPOBS ]; then
    echo "Using adaptive observations..."
#   $LDIR/adapobs.x
#   cp -f adap_newerr.grd $OUTPUT/adap_obs/$TY$TM$TD$TH.grd
#   ADAPTexe=./aoerl.x
    export MPI_VERBOSE=1
    export MPI_DISPLAY_SETTINGS=1
    export MPI_BUFS_PER_PROC=128
    export MPI_BUFS_PER_HOST=128
    export MPI_IB_RAILS=2
    export MPI_GROUP_MAX=128
    mpiexec_mpt -np $PBS_NP $ADAPTexe
  else
    echo "Using prescribed observation error."
  fi
fi

exit 0

#STEVE: should change these dates to match the analysis time $AY$AM$AD$AH.grd
if [ -f gues_me.grd ]; then
  cp gues_me.grd $OUTPUT/gues/mean/$IY$IM$ID$IH.grd
  cp gues_sp.grd $OUTPUT/gues/sprd/$IY$IM$ID$IH.grd
  cp anal_me.grd $OUTPUT/anal/mean/$IY$IM$ID$IH.grd
  cp anal_sp.grd $OUTPUT/anal/sprd/$IY$IM$ID$IH.grd
fi
if [ -f obs.grd ]; then
  cp obs.grd     $OUTPUT/obs/value/$IY$IM$ID$IH.grd
  cp obscnt.grd  $OUTPUT/obs/count/$IY$IM$ID$IH.grd
  cp obserr.grd  $OUTPUT/obs/error/$IY$IM$ID$IH.grd
fi

#STEVE: get the observation increment files
if [ -f gues_rms.grd ]; then
  cp gues_rms.grd $OUTPUT/gues/oinc_rms/$IY$IM$ID$IH.grd
  cp gues_dep.grd $OUTPUT/gues/oinc_dep/$IY$IM$ID$IH.grd
#   mv gues_lev.out $OUTPUT/gues/oinc_lev/$IY$IM$ID$IH.grd
  cp anal_rms.grd $OUTPUT/anal/oinc_rms/$IY$IM$ID$IH.grd
  cp anal_dep.grd $OUTPUT/anal/oinc_dep/$IY$IM$ID$IH.grd
#   mv anal_lev.out $OUTPUT/anal/oinc_lev/$IY$IM$ID$IH.grd
fi

