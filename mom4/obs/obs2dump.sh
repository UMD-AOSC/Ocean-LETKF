#!/bin/sh
F90=ftn
#EXP=mk3p51
#OBS=../DATA/$EXP/obs/1980010200.dat
#OBS=../DATA/obs/profs_ts/20080102.dat
#OBS=NCEP_PROF/PROFILES_LETKF/20100101.dat
#OBS=obsout.dat
#OBS=/ncrc/home1/Steve.Penny/LETKF/tmp_alt/1992092700/letkf/obs01001.dat
#OBS=obsin.dat
#OBS=~/HYBRID/tmp_robs2/1991010100/letkf/obs01001.dat
#OBS=~/HYBRID/tmp_robs2/1991010100/letkf/obs02001.dat
#OBS=~/HYBRID/tmp_robs2/1991010100/letkf/obs03001.dat
#OBS=~/HYBRID/tmp_robs2/1991010100/letkf/obs04001.dat
#OBS=~/HYBRID/tmp_robs2/1991010100/letkf/obs05001.dat

OBS=/lustre/f1/unswept/Steve.Penny/OUTPUT/tmp_letkf_mom6/2005010100/letkf/obs01001.dat

ln -s $OBS fort.3
$F90 obs2dump.f90
./a.out
rm a.out
rm fort.3
