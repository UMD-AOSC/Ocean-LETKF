#!/bin/bash

# only H5_LIB and H5_INC will be used by other scripts

if command -v h5pfc &> /dev/null
then
   echo "use parallel hdf5 wrapper: h5pfc to set hdf5 envs"
   H5_LIB=`h5pfc -show TESTSRC | awk -F"TESTSRC " '{print $2}'`
   H5_INCDIR=`h5pfc -show TESTSRC | awk -F"TESTSRC" '{print $1}' | awk -F"-I" '{print $2}'`
   H5_INC="-I${H5_INCDIR}"
else
   echo "use serial hdf5 wrapper: h5fc to set hdf5 envs"
   H5_LIB=`h5fc -show TESTSRC | awk -F"TESTSRC " '{print $2}'`
   H5_INCDIR=`h5fc -show TESTSRC | awk -F"TESTSRC" '{print $1}' | awk -F"-I" '{print $2}'`
   H5_INC="-I$H5_INCDIR"
fi

echo "[$0] H5_LIB=$H5_LIB"
echo "[$0] H5_INC=$H5_INC"

