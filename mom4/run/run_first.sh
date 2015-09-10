#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 02:00:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# This code is for generating initial conditions for different ensembles.

set -e

root=/lustre/lysun/models/Ocean-LETKF-test/mom4
root_run=${root}/run
EXP_DIR=${root}/INPUT_INIT
INPUT_INI=${root}/INPUT_INIT/INPUT_INIT

YYYYMMDDHH=1980010100
days=20
ENS_NUM=31
isfirst=1

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
  
  while [ "${MEM3}" -le "$ENS_NUM" ]
  do
    if [ "${MEM3}" -lt 100 ]; then
      MEM3=0${MEM3}
    fi
    if [ "${MEM3}" -lt 10 ]; then 
      MEM3=0${MEM3}
    fi
    
    if [ -e ${EXP_DIR}/${MEM3} ]; then
      rm -rf ${EXP_DIR}/${MEM3}
    fi

    mkdir -p ${EXP_DIR}/${MEM3}

    echo "I am member ${MEM3}"
    workdir=${EXP_DIR}/${MEM3}
    cd ${workdir}

    mkdir -p ${workdir}/INPUT
    mkdir -p ${workdir}/RESTART
    mkdir -p ${workdir}/DRIFTERS

    if [ "$isfirst" -eq "1" ]; then
      ln ${INPUT_INI}/input.nml ${workdir}/INPUT/
      ln ${INPUT_INI}/ocean_temp_salt.res.nc ${workdir}/INPUT/        
      ln ${INPUT_INI}/*_table ${workdir}/INPUT/
    else
      PMEM3=`expr $MEM3 - 1`
      if [ "${PMEM3}" -lt 100 ]; then
        PMEM3=0${PMEM3}
      fi
      if [ "${PMEM3}" -lt 10 ]; then 
        PMEM3=0${PMEM3}
      fi

      # LINK all the analytic files in letkf of previous step or from the restart file from previous step. ****FINISH Later
      # First link the files in PREVIOUS INPUT folder:
      ln ${EXP_DIR}/${PMEM3}/RESTART/input.nml ${workdir}/INPUT/input.nml
      ln ${EXP_DIR}/${PMEM3}/RESTART/ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc         
      ln ${EXP_DIR}/${PMEM3}/RESTART/*_table ${workdir}/INPUT/

    fi # END IF on isfirst

    sed -i "335s/.*/            days = $days/" ${workdir}/INPUT/input.nml # Don't forget to add "-i" to save
    sed -i "336s/.*/            date_init = ${IYYYY},${IMM},${IDD},${IHH},${INN},${ISS}/" ${workdir}/INPUT/input.nml

    ln ${INPUT_INI}/drifters_inp.nc ${workdir}/INPUT/

    cp ${INPUT_INI}/grid_spec.nc ${workdir}/INPUT
    cp ${INPUT_INI}/salt_sfc_restore.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/ssw_atten_depth.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/sw_flux.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/tau.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/temp_sfc_restore.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/t_flux.nc ${workdir}/INPUT/
    cp ${INPUT_INI}/water_flux.nc ${workdir}/INPUT/
    
    
    cp ${root_run}/mom4p1_solo_*.csh ${workdir}

    csh mom4p1_solo_compile.csh 
    # sbatch mom4p1_solo_run.csh <full path to working directory> $MEM3
    csh mom4p1_solo_run.csh ${workdir} ${MEM3}
   
    MEM3=`expr $MEM3 + 1`
    isfirst=0
    
  done # END the loop of ENS_NUM

echo 'Naturual Finish.'
exit 0
