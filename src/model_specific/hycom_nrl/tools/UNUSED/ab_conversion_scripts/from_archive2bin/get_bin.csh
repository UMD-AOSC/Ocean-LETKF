#!/bin/csh
ln -s  /scratch3/NCEPDEV/marine/noscrub/Bin.Li/hycom/saved_results/tarv_109l/041_archv.2009_335_01.a input.a
ln -s  /scratch3/NCEPDEV/marine/noscrub/Bin.Li/hycom/saved_results/tarv_109l/041_archv.2009_335_01.b input.b
../src/read_archive > output.txt
