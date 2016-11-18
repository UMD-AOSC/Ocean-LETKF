#!/bin/csh
ln -s  /scratch3/NCEPDEV/marine/noscrub/Bin.Li/hycom/saved_results/tarv_109l/041_archv.2009_335_01.a input.a
ln -s  /scratch3/NCEPDEV/marine/noscrub/Bin.Li/hycom/saved_results/tarv_109l/041_archv.2009_335_01.b input.b
ln -s ../from_archive2bin/data.bin .
rm -f output.a output.b
../src/write_archive>check.txt

