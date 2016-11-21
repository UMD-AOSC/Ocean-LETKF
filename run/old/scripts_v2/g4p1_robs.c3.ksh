#!/bin/ksh --login
set -ex

#source $MODULESHOME/init/ksh
#module use /sw/eslogin-c3/modulefiles
#module load nco
#module load wgrib

# Location of wgrib:
wgrib=/sw/xe6/wgrib/1.8.1.0b/sles11.1_gnu4.3.4/bin/wgrib
# Location of ncflint:
#ncflint=/sw/xe6/nco/4.0.8/sles11.1_netcdf4.1.3_gnu4.6.2/bin/ncflint
ncflint=/sw/eslogin-c3/nco/4.5.2/sles11.3_gnu5.1.0/bin/ncflint
# Location of ncea:
#ncea=/sw/xe6/nco/4.0.8/sles11.1_netcdf4.1.3_gnu4.6.2/bin/ncea
ncea=/sw/eslogin-c3/nco/4.5.2/sles11.3_gnu5.1.0/bin/ncea

#module load intel
#module load netcdf
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

# START:
echo "GODAS 3DVar run step" > g4p1.out
echo "Processing cycle: ${YYYYMMDDHH}" >> g4p1.out
workdir=${EXP_DATA}/${YYYYMMDDHH}/g4p1
mkdir -p ${workdir}
cd ${workdir}
export PATH=$PATH:`pwd` #STEVE: need to access wgrib and nco executables copied to here, since modules are not available.

# setup directory structure
gINPUT=$workdir/INPUT
gRESTART=$workdir/RESTART
mkdir -p $gINPUT
mkdir -p $gRESTART

OBSDIR=$gsupdir/OBSA #PTMP_SAL
AOBDIR=$gsupdir/PREP_ALT
VARDIR=$gsupdir/VARr
SBCDIR=$gsupdir/PREP_SBC

# STEVE: allow godas to use different surface forcing fields:
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`

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
inc=$days #5
inc_units=days
TY=`$date -d "$IY-$IM-$ID $inc $inc_units" +%Y`
TM=`$date -d "$IY-$IM-$ID $inc $inc_units" +%m`
TD=`$date -d "$IY-$IM-$ID $inc $inc_units" +%d`
TH=`$date -d "$IY-$IM-$ID $inc $inc_units" +%H`
TN=$IN
TS=$IS

# Update the date for the previous analysis cycle
pinc=$inc
pinc_units='days ago'
PY=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%Y`
PM=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%m`
PD=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%d`
PH=`$date -d "$IY-$IM-$ID $pinc $pinc_units" +%H`
PN=$IN
PS=$IS
workdir0=${EXP_DATA}/$PY$PM$PD$PH/g4p1

#########################################################
### Processs inputs

#STEVE:
# Compute ensemble mean of model inputs and set up in
HYBRID=${EXP_DATA}/${YYYYMMDDHH}/hybrid
mkdir -p $HYBRID

# RESTART directory constructed from LETKF analysis mean
RESTART=${EXP_DATA}/${YYYYMMDDHH}/hybrid/letkfmean
mkdir -p $RESTART
#STEVE: general input files for all timesteps are in this directory:
INPUT=$INPUT_INIT

# use nco to average RESTART files and put in RESTART/
cp $ncea .
rm -f $RESTART/ocean_temp_salt.res.nc
ncea ${EXP_DATA}/${YYYYMMDDHH}/letkf/anal???.ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc

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
   if [ ! -f $INPUT_INIT/grid_spec.nc ]; then
     echo "ERROR: required input file does not exist $INPUT_INIT/grid_spec.nc "
     exit 1
   fi
   if [ ! -f $RESTART/ocean_temp_salt.res.nc ]; then
     echo "ERROR: required input file does not exist $RESTART/ocean_temp_salt.res.nc "
     exit 1
   fi

   #STEVE: I'm not sure if this file must be in this location for the sub-shell scripts
   #MEM2 was 01
   cp -r ${EXP_DATA}/${YYYYMMDDHH}/model/$MEM2/RESTART/$rtype.res $RESTART/
   if [ ! -f $RESTART/$rtype.res ]; then
     echo "ERROR: required input file does not exist $RESTART/$rtype.res "
     exit 1
   fi

#-------------------------------------------

  ############################
  # Prepare observation data
  ############################

  cd $gINPUT #STEVE: this is on the lustre filesystem, $OBSDIR is not
  echo "In directory: "
  pwd
  rm -f tmpa.mom
  rm -f sala.mom
  rm -f *tmp.nc
  rm -f *sal.nc

  # cDG.csh (a) extracts obs profiles from yearly tar files for TMP and SAL (via c program: extDysP4nc.c)
  #         (b) reads daily profile files, combines them and writes them in NEW M4 assimilation format (via c pograms: cmbDLstTh4nc.c and cmbDLstPs4nc.c)
# $OBSDIR/cDG.csh
  # Do obs setup and Create the observation data files here instead:
  mkdir -p $gINPUT/work_tmpa
  cd $gINPUT/work_tmpa
  sinc=1
  sinc_units=day
  YYYY=$IY
  MM=$IM
  DD=$ID
  while test $YYYY$MM$DD -lt $TY$TM$TD
  do
    NY=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%Y`
    NM=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%m`
    ND=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%d`
    # Start 1 day later, end $days days later
    YYYY=$NY
    MM=$NM
    DD=$ND

    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_t ..."
    ln -f $GOBSDIR_t/${YYYY}${MM}${DD}tmp.nc .
    #STEVE: merged salinity observations must exist for computation of potential temperature:
    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_sm ..."
    ln -f $GOBSDIR_sm/${YYYY}${MM}${DD}sal.nc .
  done

  echo "Creating list of *.nc obs files for tmp and sal..."
  rm -f tmpFiles
  rm -f salFiles
  ls ????????tmp.nc > tmpFiles
  ls ????????sal.nc > salFiles

  # STEVE: 
  #              The original versions are needed when going back to real observations
  #              to convert in situ temps to potential temp (same is needed for LETKF) 
  #              and to compute the estimated standard deviation
  rm -f grid_spec.nc
  ln -f $INPUT_INIT/grid_spec.nc .
  cp -f $OBSDIR/cmbDLstTh4nc .
  echo "./cmbDLstTh4nc -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 40" > g4p1.out
# ./cmbDLstTh4nc -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 30 >> g4p1.out
  ./cmbDLstTh4nc -f tmpFiles salFiles -o tmpa.mom -g grid_spec.nc -k 40 >> g4p1.out
  ln -f tmpa.mom $gINPUT/

  if [ ! -f $gINPUT/tmpa.mom ]; then
    echo "$gINPUT/tmpa.mom does not exist! Exiting g4p1_robs.ksh..."
    exit 5701
  else
    echo "Finished preparing Potential Temperature observations."
#   chmod 0444 $gINPUT/tmpa.mom
    pwd
  fi
# rm -f *tmp.nc
# rm -f *sal.nc

  # Since we need the merged salinities above, we may want the actual observed salinities below:
  mkdir -p $gINPUT/work_sala
  cd $gINPUT/work_sala
  sinc=1
  sinc_units=day
  YYYY=$IY
  MM=$IM
  DD=$ID
  while test $YYYY$MM$DD -lt $TY$TM$TD
  do
    NY=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%Y`
    NM=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%m`
    ND=`$date -d "$YYYY-$MM-$DD $sinc $sinc_units" +%d`
    # Start 1 day later, end $days days later
    YYYY=$NY
    MM=$NM
    DD=$ND

    echo "Hard linking obs for $YYYY$MM$DD from $GOBSDIR_s ..."
    if [ -f $GOBSDIR_s/${YYYY}${MM}${DD}sal.nc ]; then
      ln -f $GOBSDIR_s/${YYYY}${MM}${DD}sal.nc .
    else
      echo "WARNING!!! Salinity file does not exist for date: $YYYY$MM$DD"
      echo "WARNING!!! Make sure there is actually no observation file!"
    fi
  done

  echo "Creating list of *.nc obs files for sal..."
  rm -f salFiles
  ls ????????sal.nc > salFiles
  rm -f grid_spec.nc
  ln -f $INPUT_INIT/grid_spec.nc .
  cp -f $OBSDIR/cmbDLstPs4nc . 
  echo "./cmbDLstPs4nc -f salFiles -o sala.mom -g grid_spec.nc -k 40" >> g4p1.out
# ./cmbDLstPs4nc -f salFiles -o sala.mom -g grid_spec.nc -k 30 >> g4p1.out
  ./cmbDLstPs4nc -f salFiles -o sala.mom -g grid_spec.nc -k 40 >> g4p1.out
  ln -f sala.mom $gINPUT/

  if [ ! -f $gINPUT/sala.mom ]; then
    echo "$gINPUT/sala.mom does not exist! Exiting g4p1_robs.ksh..."
    exit 2191
  else
    echo "Finished preparing Observed Salinity observations."
#   chmod 0444 $gINPUT/sala.mom
    pwd
  fi 

  #STEVE: do the same to set up altimetry if being used:
  if [[ $USE_ALTIMETRY -eq "1" && "TEST$GOBSDIR_a" != "TEST" ]]; then

    echo "Create and combine weekly altimetry files:"
    echo " "
    echo "Begin"
    mkdir -p $gINPUT/work_alta
    cd $gINPUT/work_alta

    #(Skip the use of cWG.csh and just implement it here:)
    latmin=-40
    latmax=50
    nlv=40
    dte=$YYYY$MM$DD

    # Remove existing files, if any:
    rm -f runSWA4
    rm -f txj1j2_cTbl
    rm -f spltWkAlt4
    rm -f cmbWksAlt4
    rm -f alta.mom

    ln -fs $INPUT_INIT/grid_spec.nc .
    cp -f $AOBDIR/runSWA4 .
    cp -f $AOBDIR/txj1j2_cTbl .
    cp -f $AOBDIR/spltWkAlt4 .
    nweeks=5 # number of weeks to combine
#   ./runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv
    echo "aprun -n 1 $gINPUT/work_alta/runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv"
    aprun -n 1 $gINPUT/work_alta/runSWA4 -d $dte -n $nweeks -f txj1j2_cTbl -p $GOBSDIR_a -m grid_spec.nc -lt $latmin $latmax -k $nlv

#   if [ ! -f runSWA4.out ]; then
    if [ ! -f split.ksh ]; then
      exit 5
    fi

    chmod 755 split.ksh
    ./split.ksh

    #STEVE: had to edit runSWA4.c to output comands to a file: runSWA4.out
#   while read line
#   do
#     # This should evaluate spltWkAlt4 for each week file:
#     echo "Evaluating line: $line"
#     eval "$line"
#   done < runSWA4.out

    ls -1 swssh.?????? > cWGLst

    cp -f $AOBDIR/cmbWksAlt4 .
    echo "./cmbWksAlt4 -f cWGLst -o alta.mom"
    ./cmbWksAlt4 -f cWGLst -o alta.mom

    if [ ! -f alta.mom ]; then
      echo "alta.mom does not exist! Exiting godas_solo.ksh..."
      exit 3
    fi

    # Copy godas support files for altimetry:
    # (1) The climatology
    cp $gsupdir/aEtaCds9399.i3e $gINPUT/
    # (2) The coefficients for the linear transformation
    cp $gsupdir/cdnzTS.m4 $gINPUT/

#   mv alta.mom $gINPUT/
    rm cWGLst
    rm swssh.??????
  fi
 
  ############################
  # Prepare background error covariance
  ############################
  # RESTART (RESTART directory from last ocean model run)
  # gridDir (directory of grid_spec.nc file)
  # scriptDir (directory of support C script)

  mkdir -p $gINPUT/work_var
  cd $gINPUT/work_var

# cat > $VARDIR/mEN4.header.csh <<EOF
  cat > mEN4.header.csh <<EOF
  set RESTART=$RESTART
  set gridDir=$INPUT_INIT
  set scriptDir=$VARDIR
  set yr=$TY
  set mo=$TM
  set dy=$TD
EOF
# set yr=$IY
# set mo=$IM
# set dy=$ID
  
  echo "Preparing background error covariance..."
# cd $VARDIR
  echo "In directory: "
  pwd
  rm -f tvv.mom
  rm -f svv.mom

  $VARDIR/mEN4.csh

  #STEVE: background error covariance files:
  if [ -f tvv.mom ]; then
    echo "tvv.mom created succesfully."
    ln -f tvv.mom $gINPUT/
  else
    echo "tvv.mom doens't exist. Exiting..."
    pwd
    exit 54
  fi
  if [ -f svv.mom ]; then
    echo "svv.mom created succesfully."
    ln -f svv.mom $gINPUT/
  else
    echo "svv.mom doens't exist. Exiting..."
    pwd
    exit 24
  fi

  ############################
  # Prepare surface boundary conditions (SST, SSS)
  ############################
  # SST_DIR (sst reanalysis data)
  # SSS_DIR (sss climatology data, salt12.nc)
  # L_YR (last year of data)
  # gridDir
  # scriptDir

  echo "Preparing surface boundary conditions..."

  mkdir -p $gINPUT/work_sfc_rst
  cd $gINPUT/work_sfc_rst
  export PATH=$PATH:`pwd` #STEVE: need to access wgrib and nco executables copied to here, since modules are not available.
  cp $wgrib $gINPUT/work_sfc_rst/
# set SBC_DIR=$FLUXDIR
# cat > $SBCDIR/mPS.header.csh <<EOF
  cat > mPS.header.csh <<EOF
  set SST_DIR=$SSTDIR
  set SSS_DIR=$INPUT_INIT/../
  @ L_YR = 2013
  set gridDir=$INPUT_INIT
  set scriptDir=$SBCDIR
  set yr=$IY
  set mo=$IM
  set dy=$ID
EOF

# cd $SBCDIR
  echo "In directory: "
  pwd

  # Prepare surface restore files:
  # salt12.nc and grid_spec.nc must be in $SBCDIR: 
  ln -f $INPUT_INIT/../salt12.nc $gINPUT/
  ln -f $INPUT_INIT/grid_spec.nc $gINPUT/
  $SBCDIR/mPS.csh

# mv $SBCDIR/*_sfc_restore.nc $gINPUT/
  if [ -f temp_sfc_restore.nc ]; then
    echo "temp_sfc_restore.nc exists." 
    ln -f temp_sfc_restore.nc $gINPUT/
  else
    echo "temp_sfc_restore.nc does not exist. Exiting..."
    exit 53
  fi

  if [ -f salt_sfc_restore.nc ]; then
    echo "salt_sfc_restore.nc exists." 
    ln -f salt_sfc_restore.nc $gINPUT/
  else
    echo "salt_sfc_restore.nc does not exist. Exiting..."
    exit 21
  fi

  ############################
  # Prep the input files for g4p1:
  ############################
  ln -f $RESTART/ocean_temp_salt.res.nc $gINPUT/
  ln -f ${EXP_DATA}/${YYYYMMDDHH}/model/$MEM2/RESTART/$rtype.res $gINPUT/
  ln -f $INPUT_INIT/grid_spec.nc $gINPUT/

  ############################
  cd $workdir
# cp $wgrib ${workdir}/

  datatable=$gsupdir/data_table            # path to the data override table.
  if [ ! -f $namelist ]; then
    echo "ERROR: required input file does not exist $namelist "
    exit 1
  fi
  cp -f $datatable $workdir

  namelist=$gsupdir/namelist   # path to namelist file w/o coupler_nml
  if [ ! -f $datatable ]; then
    echo "ERROR: required input file does not exist $datatable "
    exit 1
  fi
  cp -f $namelist $workdir/

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
  cp -f $datatable  $gINPUT/data_table
# rm -f wgrib

# START RUN G4P1 ################################################

  cd $workdir
  aprun -n $PBS_NP $executable

# END   RUN G4P1 ################################################

#########################################
# Do the Hybrid-LETKF processing
#########################################

#STEVE:
## Add g4p1 ocean correction to the input mean field
#cp $gRESTART/ocean_cor.res.nc $gRESTART/ocean_cor2.res.nc
#echo "Calling:"
#echo "ncrename -O -v tcor,temp -v scor,salt $gRESTART/ocean_cor2.res.nc"
#ncrename -O -v tcor,temp -v scor,salt $gRESTART/ocean_cor2.res.nc

#echo "Calling:"
#echo "ncflint -w 1,1 $RESTART/ocean_temp_salt.res.nc $gRESTART/ocean_cor.res.nc $gRESTART/ocean_temp_salt.res.nc"
#ncflint -w 1,1 $RESTART/ocean_temp_salt.res.nc $gRESTART/ocean_cor2.res.nc $gRESTART/ocean_temp_salt.res.nc
##STEVE: !!!!! Changed godas to output this file: $gRESTART/ocean_temp_salt.res.nc

#STEVE:
# Next, combine the GODAS result (gRESTART/ocean_cor.res.nc) with the 
# LETKF mean analysis and recenter the ensemble around the new analysis
#
cp $ncflint .
# Remove target file for new hybrid solution:
rm -f $HYBRID/ocean_temp_salt.res.nc
#alpha=0.5 #1.0 #0.5
#NOTE: alpha is now an input value
one_minus_alpha=`echo "scale=2; 1 - $alpha" | bc -l`
if [[ $alpha -gt 0 && $alpha -lt 1 ]]; then
  echo "alpha   = $alpha"
  echo "1-alpha = $one_minus_alpha"
  echo "Calling:"
  echo "ncflint -w $alpha,$one_minus_alpha $gRESTART/ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc"
  ncflint -w $alpha,$one_minus_alpha $gRESTART/ocean_temp_salt.res.nc $RESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc
elif [[ $alpha -eq 1 ]]; then
  echo "alpha   = $alpha"
  echo "1-alpha = $one_minus_alpha"
  echo "Calling:"
  echo "ln -f $gRESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc"
  ln -f $gRESTART/ocean_temp_salt.res.nc $HYBRID/ocean_temp_salt.res.nc
else
  echo "alpha is not in range 0 < alpha <= 1"
  echo "alpha   = $alpha"
  echo "Exiting..."
  exit 1
fi

##################################################################################
#STEVE:
# Finally, recenter all of the letkf ocean_temp_salt.res.nc analysis fields with the new mean field

# (1) Create mean correction field {(new mean) - (old mean)}
LETKF=${EXP_DATA}/${YYYYMMDDHH}/letkf
echo "Calling:"
echo "ncflint -w 1,-1 $HYBRID/ocean_temp_salt.res.nc $HYBRID/letkfmean/ocean_temp_salt.res.nc $HYBRID/hybrid_cor.res.nc"
rm -f $HYBRID/hybrid_cor.res.nc
ncflint -w 1,-1 $HYBRID/ocean_temp_salt.res.nc $HYBRID/letkfmean/ocean_temp_salt.res.nc $HYBRID/hybrid_cor.res.nc

# (2) Add to each member
MEM=1
while test $MEM -le $MEMBERS
do
  MEM3=`printf %.3d ${MEM}` 
  afile=anal${MEM3}.ocean_temp_salt.res.nc
  # Add the correction to this member
  echo "Calling:" 
  echo "ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile"
  rm -f $HYBRID/$afile
# ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile &
  ncflint -w 1,1 $HYBRID/hybrid_cor.res.nc $LETKF/$afile $HYBRID/$afile
  MEM=`expr $MEM + 1`
done
rm -f ncflint

#STEVE: wait for all of the background jobs to finish
#wait

# Generate netcdf background ensemble mean for viewing and testing output
cp $ncea .
rm -f $HYBRID/gues05_temp_salt.res.nc
ncea $LETKF/gs05???.ocean_temp_salt.res.nc $HYBRID/gues05_temp_salt.res.nc
rm -f ncea

#Write out a file indicating the processing is complete.
echo "This is the G4P1 run and Hybrid-LETKF assimilation for cycle ${YYYYMMDDHH}" >> g4p1.out

echo "REMOVING DAILY MODEL RESTART FILES!!!"
rm -f ${EXP_DATA}/${YYYYMMDDHH}/model/??/RESTART/${YYYYMMDDHH}*


exit 0

