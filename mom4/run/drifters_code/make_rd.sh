#!/bin/sh
# Special compile script just for UMD's DT2

f90file=read_drifters.f90
exefile=rd.x

a=0
while [ $a -lt 10 ]
do
  file1="drifters_out.nc.0000 $a"
  file2=drifters_out.nc.00000$a
  if [ -f "$file1" ]; then
    mv "$file1" $file2
  fi
  a=`expr $a + 1`
done

ifort -I$NETCDF_FORTRAN_INCDIR $f90file -o $exefile -L$NETCDF_FORTRAN_LIBDIR -lnetcdff -Wl,-rpath,$NETCDF_FORTRAN_LIBDIR

txtfile=drifters_inp.txt
outfile=drifters_inpc.txt

if [ -f $txtfile ]; then
  rm -rf $txtfile
fi
./rd.x -np 19 -o $txtfile > rd.out # "79" depends on how many processors you have used.


if [ -f $outfile ]; then
  rm -rf $outfile
fi

head -1 $txtfile > a.out
num_time=`awk '{print $6}' a.out`
rm -rf a.out

head -2 $txtfile > grd.drifters_inp.txt

head -2 $txtfile > $outfile
awk '$8=='`expr $num_time - 1`'' $txtfile >> $outfile

 