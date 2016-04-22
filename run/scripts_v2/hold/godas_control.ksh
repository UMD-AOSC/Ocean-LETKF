#!/bin/ksh --login
#module load mpt
#module load intel
#module load netcdf/4.1.3-intel
#module load nco
#module load wgrib

#INPUTS:
# EXP_DATA
# executable
# YYYYMMDDHH
# gsupdir 
# days
# months
# rtype
# REGRID

# Mean SFC forcing:
FLUXDIR=/scratch2/portfolios/NCEPDEV/climate/noscrub/Steve.Penny/SFLUX/DAILYnc/00
OBSDIR=$gsupdir/PTMP_SAL
VARDIR=$gsupdir/VARr
SBCDIR=$gsupdir/PREP_SBC

echo "GODAS CONTROL run step"
echo "Processing cycle: ${YYYYMMDDHH}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/godas_control
mkdir -p ${workdir}
cd ${workdir}

# setup directory structure
gINPUT=$workdir/INPUT
gRESTART=$workdir/RESTART
mkdir -p $gINPUT
mkdir -p $gRESTART

LETKF=${EXP_DATA}/${YYYYMMDDHH}/letkf

# RESTART directory constructed from background mean
RESTART=$workdir/bgmean
mkdir -p $RESTART

#STEVE: active code:
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
PYMDH=$PY$PM$PD$PH
workdir0=${EXP_DATA}/$PYMDH/godas_control

#########################################################
### Processs inputs

#STEVE:
# Compute ensemble mean of model inputs and set up in
# INPUT directory of GODAS CONTROL (done below into gINPUT)
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_temp_salt.res.nc $INPUT/ocean_temp_salt.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_velocity.res.nc $INPUT/ocean_velocity.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_sbc.res.nc $INPUT/ocean_sbc.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_barotropic.res.nc $INPUT/ocean_barotropic.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_bih_friction.res.nc $INPUT/ocean_bih_friction.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_con_temp.res.nc $INPUT/ocean_con_temp.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_density.res.nc $INPUT/ocean_density.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_frazil.res.nc $INPUT/ocean_frazil.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_neutralB.res.nc $INPUT/ocean_neutralB.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_neutral.res.nc $INPUT/ocean_neutral.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_residency.res.nc $INPUT/ocean_residency.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_sigma_transport.res.nc $INPUT/ocean_sigma_transport.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_thickness.res.nc $INPUT/ocean_thickness.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_tracer.res.nc $INPUT/ocean_tracer.res.nc &
#ncea ${EXP_DATA}/${YMDH}/model/??/INPUT/ocean_velocity_advection.res.nc $INPUT/ocean_velocity_advection.res.nc &

#wait

# Copy the rest of the non-varying data
cp ${EXP_DATA}/${PYMDH}/model/01/RESTART/ocean_solo.intermediate.res $gINPUT/ocean_solo.intermediate.res
cp ${EXP_DATA}/${PYMDH}/model/01/RESTART/ocean_solo.res $gINPUT/ocean_solo.res

#########################################################
# GODAS solo code:
# Minimal runscript: gdsSolo

#npes=$PBS_NP         # number of processors
                     # Note: If you change npes you may need to change
                     # the layout in the corresponding namelist
#===========================================================================

# Check the time_stamp
  echo "${IY}-${IM}-${ID} is the current date"
  echo "$IY	$IM	$ID	$IH	$IN	$IS" > time_stamp.out

#
#Check the existance of essential input files
#
   if [ ! -f $REGRID/grid_spec.nc ]; then
     echo "ERROR: required input file does not exist $REGRID/grid_spec.nc "
     exit 1
   fi

   # Generate background ensemble mean for INPUTTING INTO GODAS (previoualy viewing and testing output)
   ncea $LETKF/gs05???.ocean_temp_salt.res.nc $RESTART/gues05_temp_salt.res.nc
   cp $RESTART/gues05_temp_salt.res.nc $gINPUT/ocean_temp_salt.res.nc
   if [ ! -f $gINPUT/ocean_temp_salt.res.nc ]; then
     echo "ERROR: required input file does not exist $gINPUT/ocean_temp_salt.res.nc "
     exit 1
   fi

   #STEVE: I'm not sure if this file must be in this location for the sub-shell scripts
   cp -r ${EXP_DATA}/${YYYYMMDDHH}/model/01/RESTART/$rtype.res $gINPUT/$rtype.res
   if [ ! -f $gINPUT/$rtype.res ]; then
     echo "ERROR: required input file does not exist $gINPUT/$rtype.res "
     exit 1
   fi

#-------------------------------------------

  ############################
  # Prepare observation data
  ############################
  # arcdir (directory with profile observations)
  # here (current directory)
  # expdir ()
  # d1, d2 (range of dates to look forward and backward for observations)
  # gridDIR (directory with grid_spec.nc file)
  # scriptDir (directory of support C script)

  cat > $OBSDIR/cDG.header.csh <<EOF
  set arcdir=$GOBSDIR
  set d1=0
  set d2=4
  set gridDir=$REGRID
  set scriptDir=$OBSDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $OBSDIR
  echo "In directory: "
  pwd
  rm -f tmpa.mom
  rm -f sala.mom

  $OBSDIR/cDG.csh

  #STEVE: observation files:
  mv $OBSDIR/tmpa.mom $gINPUT/
  mv $OBSDIR/sala.mom $gINPUT/

  ############################
  # Prepare background error covariance
  ############################
  # RESTART (RESTART directory containing ocean_temp_salt.res.nc file, e.g. from last ocean model run)
  # gridDir (directory of grid_spec.nc file)
  # scriptDir (directory of support C script)

  cat > $VARDIR/mEN4.header.csh <<EOF
  set RESTART=$gINPUT
  set gridDir=$REGRID
  set scriptDir=$VARDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $VARDIR
  echo "In directory: "
  pwd
  rm -f tvv.mom
  rm -f svv.mom

  $VARDIR/mEN4.csh

  #STEVE: background error covariance files:
  mv $VARDIR/tvv.mom $gINPUT/
  mv $VARDIR/svv.mom $gINPUT/

  ############################
  # Prepare surface boundary conditions (SST, SSS)
  ############################
  # SBC_DIR (surface forcing from the atmos model)
  # SST_DIR (sst reanalysis data)
  # SSS_DIR (sss climatology data, salt12.nc)
  # L_YR (last year of data)
  # gridDir
  # scriptDir

  dataDir=/scratch2/portfolios/NCEPDEV/climate/noscrub/David.Behringer/SBC
  cat > $SBCDIR/mPS.header.csh <<EOF
  set SBC_DIR=$FLUXDIR
  set SST_DIR=$dataDir/SST2/DAILY
  set SSS_DIR=$EXP_DATA/INPUT
  @ L_YR = 2013
  set gridDir=$REGRID
  set scriptDir=$SBCDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

  cd $SBCDIR
  echo "In directory: "
  pwd

  $SBCDIR/mPS.csh

  mv $SBCDIR/*_sfc_restore.nc $gINPUT/

  ############################
  # Prep the input files for godas_control:
  ############################
  cd $gINPUT
# cp $RESTART/ocean_temp_salt.res.nc $gINPUT/

  cp -r ${EXP_DATA}/${YYYYMMDDHH}/model/01/RESTART/$rtype.res $gINPUT/
  ln -fs $REGRID/grid_spec.nc $gINPUT/

  ############################
  cd $workdir

  datatable=$gsupdir/data_table            # path to the data override table.
  if [ ! -f $namelist ]; then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  fi
  ln -fs $datatable $workdir

  namelist=$gsupdir/namelist   # path to namelist file w/o coupler_nml
  if [ ! -f $datatable ]; then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  fi
  ln -fs $namelist $workdir/

  ############################
  # set input.nml
  ############################

 cat > input.nml <<EOF
 &godas_solo_nml
        months = $months,
        days   = $days,
        hours=0
        minutes=0
        seconds=0
        calendar = 'julian',
/
EOF
  cat $namelist >> input.nml
  ln -fs $datatable  $gINPUT/data_table

# START RUN G4P1 ################################################

 export MPI_VERBOSE=1
 export MPI_DISPLAY_SETTINGS=1
 #export MPI_BUFS_PER_PROC=128
 #export MPI_BUFS_PER_HOST=128
 #export MPI_IB_RAILS=2
 #export MPI_GROUP_MAX=128
  mpiexec_mpt -np $PBS_NP $executable

# END   RUN G4P1 ################################################

#Write out a file indicating the processing is complete.
echo "This is the GODAS Control run assimilation for cycle ${YYYYMMDDHH}" > $workdir/gc.out

exit 0

