#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst>
set -e

root=$1
root_run=${root}/run
EXP_DIR=${root}/OUTPUT
INPUT_INI=${root}/INPUT_INIT

YYYYMMDDHH=$2
days=$3
ENS_NUM=$4
isfirst=$5

IYYYY=${YYYYMMDDHH:0:4}
IMM=${YYYYMMDDHH:4:2}
IDD=${YYYYMMDDHH:6:2}
IHH=${YYYYMMDDHH:8:2}
INN=00
ISS=00
MEM3=1


echo "==================================================================="
echo "Running MOM4p1_drifters"
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  
  while [ "${MEM3}" -le "$ENS_NUM" ]
  do
    if [ "${MEM3}" -lt 100 ]; then
      MEM3=0${MEM3}
    fi
    if [ "${MEM3}" -lt 10 ]; then 
      MEM3=0${MEM3}
    fi

    echo "I am member ${MEM3}"
    workdir=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/model/${MEM3}
    mkdir -p ${workdir}
    cd ${workdir}

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART
    mkdir -p ${workdir}/DRIFTERS

    if [ "$isfirst" -eq "1" ]; then
      cp ${INPUT_INI}/${MEM3}/INPUT/input.nml ${workdir}/INPUT/
      cp ${INPUT_INI}/${MEM3}/INPUT/drifters_inp.nc ${workdir}/INPUT/
      cp ${INPUT_INI}/${MEM3}/INPUT/ocean_temp_salt.res.nc ${workdir}/INPUT/
         
      cp ${INPUT_INI}/${MEM3}/INPUT/*_table ${workdir}/INPUT/
    else
      # Update the date for the previous analysis cycle
      date=/bin/date
      pinc=$days
      pinc_units='days ago'
      PYYYY=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%Y`
      PMM=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%m`
      PDD=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%d`
      PHH=`$date -d "$IYYYY-$IMM-$IDD $pinc $pinc_units" +%H`
      PNN=$INN
      PSS=$ISS

      # LINK all the analytic files in letkf of previous step or from the restart file from previous step. ****FINISH Later
      # First link the files in PREVIOUS INPUT folder:
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/input.nml ${workdir}/INPUT/input.nml
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/drifters_inp.nc ${workdir}/INPUT/drifters_inp.nc
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc         
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/model/${MEM3}/RESTART/*_table ${workdir}/INPUT/

      # Next, check whether the previous step has LETKF:
      echo "Checking whether the Standard-LETKF Analysis from previous step exists. If yes, link as initial conditions......"
      workdir_analysis=${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/letkf
      if [ -e ${workdir_analysis} ]; then

        ##### TEMPERATURE and SALINITY file:
        if [ -f ${workdir_analysis}/anal${MEM3}.ocean_temp_salt.res.nc ]; then
          echo "The anal${MEM3}.ocean_temp_salt.res.nc exists... and start linking...."
          rm -rf ${workdir}/INPUT/ocean_temp_salt.res.nc
          ln -f ${workdir_analysis}/anal${MEM3}.ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc
        else
          echo "ANALYSIS FILE DOES NOT EXIST: anal${MEM3}.ocean_temp_salt.res.nc..."
        fi
       
        #### DRIFTERS file:
        if [ -f ${workdir_analysis}/anal${MEM3}.drifters_inp.nc ]; then
          echo "The anal${MEM3}.drifters_inp.nc exists..."
          rm -rf ${workdir}/INPUT/drifters_inp.nc
          ln -f ${workdir_analysis}/anal${MEM3}.drifters_inp.nc ${workdir}/INPUT/drifters_inp.nc
        else
          echo "ANALYSIS FILE DOES NOT EXISTT: anal${MEM3}.drifters_inp.nc..."
        fi
      fi
    fi # END IF on isfirst

    sed -i "335s/.*/            days = $days/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "336s/.*/            date_init = ${IYYYY},${IMM},${IDD},${IHH},${INN},${ISS}/" ${workdir}/INPUT/input.nml

    cp ${INPUT_INI}/${MEM3}/INPUT/grid_spec.nc ${workdir}/INPUT
    cp ${INPUT_INI}/${MEM3}/INPUT/salt_sfc_restore.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/ssw_atten_depth.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/sw_flux.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/tau.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/temp_sfc_restore.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/t_flux.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/${MEM3}/INPUT/water_flux.nc ${workdir}/INPUT/
    
    
    cp ${root_run}/mom4p1_solo_*.csh ${workdir}

    csh mom4p1_solo_compile.csh 
    # sbatch mom4p1_solo_run.csh <full path to working directory> $MEM3
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}
 
    cp ${root_run}/drifters_code/* ${workdir}/DRIFTERS/
        cp ${root_run}/drifters_code/* ${workdir}/DRIFTERS/

    cd ${workdir}/DRIFTERS/
    sh make_rd.sh # LUYU: If you are going to change the number of processor running the model, you need to modify the corresponding number in make_rd.sh
    sh make_cd.sh
    
    if [ -e ${workdir}/RESTART/drifters_inp.nc ]; then
      rm -rf ${workdir}/RESTART/drifters_inp.nc
      cp drifters_inp.nc ${workdir}/RESTART/
    fi
    
    cd ${workdir} 
   
    MEM3=`expr $MEM3 + 1`
    
  done # END the loop of ENS_NUM

echo 'Naturual Finish.'
exit 0
