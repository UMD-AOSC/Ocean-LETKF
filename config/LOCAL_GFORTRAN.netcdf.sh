#!/bin/bash

#NETCDF_DIR=/opt/local
#NETCDFF_DIR=/opt/local
#NETCDF_LIB="-L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk -arch x86_64 -lnetcdf -lnetcdf -lm"
#NETCDF_INC="-I/opt/local/include"

#NETCDF_DIR=/opt/local
#NETCDFF_DIR=/opt/local
NETCDF_LIB=`nf-config --flibs`
NETCDF_INC=`nf-config --fflags`
echo "[$0] NETCDF_LIB=$NETCDF_LIB"
echo "[$0] NETCDF_INC=$NETCDF_INC"
