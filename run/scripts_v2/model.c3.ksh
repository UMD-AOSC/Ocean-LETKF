#!/bin/ksh --login
set -ex

#module load intel
#module load netcdf
#module load nco

echo "Model step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/model/${MEMBERID}
mkdir -p ${workdir}
cd ${workdir}

#STEVE: active code:
MEM3=`printf %.3d ${MEMBERID}`
MEM2=`printf %.2d ${MEMBERID}`
TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Start date of next analysis cycle
date=/bin/date
inc=$days #5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

#STEVE: need mpiexec to run across all procs
echo $0
if [ ! -f ${workdir}/INPUT/ocean_temp_salt.res.nc ]; then
  echo "Input file does not exist... but it should: ${workdir}/INPUT/ocean_temp_salt.res.nc"
  exit 1
fi

echo "Running Member::$MEM3"
cd ${workdir}
echo "Running model ensemble members..."

#STEVE: clear out RESTART directory in case this is rerunning a failed attempt
rm -f $workdir/RESTART/*
rm -f $workdir/ocean_??.nc
rm -f $workdir/ocean_??.nc.*

if [ ! -f ${workdir}/INPUT/grid_spec.nc ]; then
  echo "ERROR: Input file does not exist... but it should: ${workdir}/INPUT/grid_spec.nc"
  echo "Linking $INPUT_INIT/grid_spec.nc to ${workdir}/INPUT/grid_spec.nc"
  ln -f $INPUT_INIT/grid_spec.nc ${workdir}/INPUT/
fi

#STEVE: Submit model run
#STEVE: make sure there are no intermediate files if the run restarted:
rm -f ${workdir}/ocean_??.nc.????

###############################################################
# RUN THE MODEL
###############################################################
##STEVE: copy the model executable, mostly for running in debugging purposes:
#cp $MOM4run ${workdir} 
# (note: Calls aprun within this script, a machine-dependent mpi command)
#$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_$mtype $PBS_NP

## REPLACE WITH:
inputDataDir=$workdir/INPUT                                  # This is path to the directory that contains the input data for this experiment.
diagtable=$inputDataDir/diag_table  # path to diagnositics table
datatable=$inputDataDir/data_table  # path to the data override table.
fieldtable=$inputDataDir/field_table # path to the field table
namelist=$inputDataDir/input.nml   # path to namelist file

if [ ! -e $namelist ]; then
  echo "ERROR: required input file does not exist $namelist"
  exit 1
fi
if [ ! -e $datatable ]; then
  echo "ERROR: required input file does not exist $datatable"
  exit 1
fi
if [ ! -e $diagtable ]; then
  echo "ERROR: required input file does not exist $diagtable"
  exit 1
fi
if [ ! -e $fieldtable ]; then
  echo "ERROR: required input file does not exist $fieldtable"
  exit 1
fi

cp -f $namelist   $workdir/input.nml
cp -f $datatable  $workdir/data_table
cp -f $diagtable  $workdir/diag_table
cp -f $fieldtable $workdir/field_table

platform=ftn                                                 # A unique identifier for your platform
etype=mom4p1_$mtype
executable=$mroot/exec_$platform/$etype/fms_$etype.x   # executable that was created after compilation

cd $workdir
cp $executable .
echo "Calling: aprun -n $PBS_NP $workdir/fms_$etype.x > $workdir/fmt.out"
aprun -n $PBS_NP $workdir/fms_$etype.x > $workdir/fmt.out

##----------------------------------------------------------------------------------------------
# Cleanup model output files:
##----------------------------------------------------------------------------------------------

##----------------------------------------------------------------------------------------------
## get a tar restart file
cp $workdir/input.nml $workdir/RESTART/
cp $workdir/*_table $workdir/RESTART/

cd $workdir
mkdir -p history
mkdir -p ascii

##----------------------------------------------------------------------------------------------
## rename ascii files with the date
for out in *.out
do
  mv $out ascii/$YYYYMMDDHH.$out
done

# STEVE: prevent error with too many hard links to grid_spec.nc:
rm -f ${workdir}/INPUT/grid_spec.nc

# Output a file containing
echo "$MOM4run ${TMPDIR}/model ${MEMBERID} mom4p1_$mtype $PBS_NP" > model_run_command.txt
echo "model run for ${MEMBERID} is finished." > ${workdir}/model.out

exit 0
