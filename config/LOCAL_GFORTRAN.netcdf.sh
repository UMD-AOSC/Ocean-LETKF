#!/bin/bash

NETCDF_DIR=/opt/local
NETCDFF_DIR=/opt/local
NETCDF_LIB="-L/opt/local/lib -lnetcdff -L/opt/local/lib -Wl,-headerpad_max_install_names -Wl,-syslibroot,/Library/Developer/CommandLineTools/SDKs/MacOSX10.15.sdk -arch x86_64 -lnetcdf -lnetcdf -lm"
NETCDF_INC="-I/opt/local/include"
