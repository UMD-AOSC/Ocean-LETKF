#!/bin/bash

NETCDF_DIR=/apps/hdf4/4.2.7-intel/
NETCDFF_DIR=/apps/netcdf/4.3.0-intel/
NETCDF_LIB="-L$NETCDFF_DIR/lib -lnetcdff -Wl,-rpath,$NETCDFF_DIR/lib"
NETCDF_INC="-I$NETCDFF_DIR/include -I$NETCDF_DIR/include"
