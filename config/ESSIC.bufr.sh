#!/bin/bash

# BUFR lib is only used for model=hycom

BUILD_BUFRLIB="yes"  # "yes" or "no"

# variables to set if BUILD_BUFRLIB = "no"
BUFR_DIR=/Users/cda14/Desktop/bufr_util/src/bufr_10.2.3
BUFR_LIB="-L$BUFR_DIR/ -lbufr"

# variables to set if BUILD_BUFRLIB = "yes"
CC=gcc

#--------------------------------------------------------------------
# build BUFRLIB under support:
# uses F90s and F90_OPT in MACHINE.fortran.sh
# (generally you should not need to modify below if you use 
#  Intel or GNU compiler. Otherwise you need to modify below )
#
if [ "$BUILD_BUFRLIB" = "yes" ]; then
   BUILD_DIR=$PWD
   SUPPORT_DIR="${BUILD_DIR}/../support"
   BUFR_DIR_NAME="bufr_10.2.3_LE"
   BUFR_PATH="$SUPPORT_DIR/$BUFR_DIR_NAME"
   BUFR_URL="https://github.com/cd10kfsu/${BUFR_DIR_NAME}.git"
   if [ -e $BUFR_PATH/libbufr.a ]; then
       echo "bufrlib found at: ${SUPPORT_DIR}/libbufr.a"
   else
       echo "bufrlib not found at: $SUPPORT_DIR"
       if [ -d $BUFR_PATH ]; then
          echo "bufrlib directory found at ($BUFR_PATH) but no libbufs.a: removing this directory"
          rm -rf $BUFR_PATH
       fi
       # clone bufrlibs under support/
       cd $SUPPORT_DIR && git clone $BUFR_URL
       cd $BUFR_DIR_NAME

       # setting fortran compiler flags
       FC=${F90s}
       if [ "$($FC --version|grep -i 'intel\|ifort'|wc -l)0" -gt 0 ]; then
          # Intel compiler
          FC_BUFR_OPT="$F90_OPT -DUNDERSCORE"
       elif [ "$($FC --version|grep -i 'gnu\|gcc\|gfortran'|wc -l)0" -gt 0 ]; then
          # GNU compiler
          FC_BUFR_OPT="$F90_OPT -DUNDERSCORE -fno-second-underscore"
       else
          echo "[error] unrecognized Fortran compiler. Exit..."
          exit 1
       fi

       # setting C compiler flags
       if [ "$($CC --version|grep -i 'intel\|icc'|wc -l)0" -gt 0 ]; then
          # Intel compiler
          CC_BUFR_OPT="-DUNDERSCORE"
       elif [ "$($CC --version|grep -i 'gnu\|gcc'|wc -l)0" -gt 0 ]; then
          # GNU compiler
          CC_BUFR_OPT="-std=c90 -DUNDERSCORE"
       else
          echo "[error] unrecognized C compiler. Exit..."
          exit 1
       fi

       # print config
       echo "======================================="
       echo "config for building BUFRLIBS"
       echo "CC=$CC"
       echo "CC_BUFR_OPT=$CC_BUFR_OPT"
       echo "FC=$FC"
       echo "FC_BUFR_OPT=$FC_BUFR_OPT"

       # build libs
       $CC $CC_BUFR_OPT -c *.c
       $FC $FC_BUFR_OPT -c *.f
       ar crv libbufr.a *.o
   fi

   cd $BUILD_DIR

   BUFR_DIR="$BUFR_PATH"
   BUFR_LIB="-L$BUFR_DIR/ -lbufr" 

   echo "BUFR_DIR=$BUFR_DIR"
   echo "BUFR_LIB=$BUFR_LIB"

fi

