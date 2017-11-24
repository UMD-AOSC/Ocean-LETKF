#!/bin/bash
set -ex

# Get initial start time
start0=`date +%s`

module load cdo

source params.sh

MPPCOMB_DIR=/perm/rd/nemc/ifs_build_nemo_git_craypat/nemofcm_build/ecmwf-cray-ftn-opt_default/build-tools/bin
MPPCOMB_EXE=mppcomb.exe

ftypes=(amean asprd bmean bsprd ainc)

# Run postprocessing on all subgrids:
SCRATCH=$GLOBAL_SCRATCH
EXPNAME=$GLOBAL_EXPNAME
WORKDIR=$SCRATCH/$EXPNAME/WORK
OUTDIR=$WORKDIR/output

mkdir -p $OUTDIR
cd $OUTDIR

# For each filetype, combine the netcdf files to form a global gridfile in the OUTPUT directory:
for f in ${ftypes[@]}; do

  filename=${f}${GLOBAL_LETKF_FILE_SUFFIX}
  infiles="$WORKDIR/????/$filename"
  outfile="$OUTDIR/$filename"

  $MPPCOMB_DIR/$MPPCOMB_EXE $outfile $infiles

done

end=`date +%s`
runtime0=$((end-start0))
echo "runtime :: $runtime0"
