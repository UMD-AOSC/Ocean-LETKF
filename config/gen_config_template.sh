#!/bin/bash

#
# generate template for machine config files
#

if [ $# -ne 1 ]; then
   echo "$0: Generate template for machine config files"
   echo "usage: "
   echo "  $0 CONFIG_FILE_PREFIX"
   exit 1
fi

machine=$1

if [ -e flist_LOCAL_GFORTRAN ]; then
    rm -f flist_LOCAL_GFORTRAN
fi

ls LOCAL_GFORTRAN.*.sh > flist_LOCAL_GFORTRAN
while read fname_template; do
    fend=$(echo $fname_template|cut -d "." -f2-)
    fname_out="${machine}.${fend}"
    echo "$fname_out <---------- $fname_template"
    cp $fname_template $fname_out
done<flist_LOCAL_GFORTRAN

if [ -e flist_LOCAL_GFORTRAN ]; then
    rm -f flist_LOCAL_GFORTRAN
fi
