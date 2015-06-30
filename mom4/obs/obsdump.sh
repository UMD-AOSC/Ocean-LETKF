#!/bin/sh
F90=ftn
#EXP=mk3p51
#OBS=../DATA/$EXP/obs/1980010200.dat
#OBS=../DATA/obs/profs_ts/20080102.dat
#OBS=NCEP_PROF/PROFILES_LETKF/20100101.dat
#OBS=obsin.dat #~/lf1u/OBS/synthetic/letkf_fmt/PROFS_gerr_TS/19910106.dat
#OBS=~/lf1u/OBS/synthetic/letkf_fmt/PROFS_gerr_TS_ALT/19920924.dat
#OBS=obsin.dat
#OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS/19910101.dat
#OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS/19910102.dat
#OBS=~/lf1u/OBS/synthetic/letkf_fmt/PROFS_gerr_TS/19910102.dat
#OBS=obs01001.dat
#OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep/20020101.dat
OBS=~/lf1u/OBS/historical/letkf_fmt/PROFS_gerr_TS_deep/20030712.dat
ln -s $OBS fort.3
$F90 obsdump.f90
./a.out
rm a.out
rm fort.3

