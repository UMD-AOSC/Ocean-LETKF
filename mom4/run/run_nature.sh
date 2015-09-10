#!/bin/sh
#SBATCH -n 20         
#SBATCH -t 06:30:00  
#SBATCH -A aosc-hi 
#SBATCH --mail-type=ALL
#SBATCH --mail-user=lysun@umd.edu

# sh run_mom4p1 <root> <TIME> <days> <ENS_NUM> <isfirst>
set -e

root=/lustre/lysun/models/Ocean-LETKF-test/mom4
root_run=${root}/run
EXP_DIR=${root}/OBS
INPUT_INI=${root}/INPUT_INIT

IYYYYMMDDHH=1980010100
EYYYYMMDDHH=1980030100
days=1
ENS_NUM=1 # This is a nature run, so only one ensemble.
isfirst=1

IYYYY=${IYYYYMMDDHH:0:4}
IMM=${IYYYYMMDDHH:4:2}
IDD=${IYYYYMMDDHH:6:2}
IHH=${IYYYYMMDDHH:8:2}
INN=00
ISS=00
MEM3=031 # ONLY have one ensemble, so set up MEM3 with 0('s)



while [ "${IYYYY}${IMM}${IDD}${IHH}" -le "${EYYYYMMDDHH}" ]
do
echo "==================================================================="
echo "Running MOM4p1_drifters for NATURE RUN."
echo "processing cycle: ${IYYYY}${IMM}${IDD}${IHH}"

  if [ -e ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH} ]; then
    rm -rf ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}
  fi

  mkdir ${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}

    echo "I am member ${MEM3}"
    workdir=${EXP_DIR}/${IYYYY}${IMM}${IDD}${IHH}/${MEM3}
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
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/input.nml ${workdir}/INPUT/input.nml
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/drifters_inp.nc ${workdir}/INPUT/drifters_inp.nc
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/ocean_temp_salt.res.nc ${workdir}/INPUT/ocean_temp_salt.res.nc         
      ln ${EXP_DIR}/${PYYYY}${PMM}${PDD}${PHH}/${MEM3}/RESTART/*_table ${workdir}/INPUT/

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

    cd ${workdir}/DRIFTERS/
    sh make_rd.sh # LUYU: If you are going to change the number of processor running the model, you need to modify the corresponding number in make_rd.sh
    sh make_cd.sh
    
    if [ -e ${workdir}/RESTART/drifters_inp.nc ]; then
      rm -rf ${workdir}/RESTART/drifters_inp.nc
      cp drifters_inp.nc ${workdir}/RESTART/
    fi

    cd ${workdir}  

isfirst=0

# Update the time
date=/bin/date
sinc=$days
sinc_units='days'
TYYYY=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%Y`
TMM=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%m`
TDD=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%d`
THH=`$date -d "$IYYYY-$IMM-$IDD $sinc $sinc_units" +%H`
TNN=$INN
TSS=$ISS

IYYYY=${TYYYY}
IMM=${TMM}
IDD=${TDD}
IHH=${THH}
INN=${TNN}
ISS=${TSS}

done # END the loop of time

echo 'Naturual Finish.'
exit 0
