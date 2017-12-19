#!/bin/bash

MPI_DIR=/common/openmpi/el7/pgi15.7/1.10.2
MPI_INC="-I$MPI_DIR/include"
MPI_LIB="-L$MPI_DIR/lib"

OMP="$MPI_INC"
