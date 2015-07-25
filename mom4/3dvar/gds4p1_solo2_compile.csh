#!/bin/csh -f
# Minimal compile script for godas solo meaning stand-alone analysis

# set echo
set platform      = ftn            # A unique identifier for your platform
                                   # It could be one of: ia64, ibm, ifc, nec, pgi, sgi, t3e, t90 
                                   # This corresponds to the mkmf templates in $root/bin dir.
#DAVE
set LIB_NETCDF    = "/opt/cray/netcdf/4.2.0/intel/120/lib/libnetcdf.a /opt/cray/netcdf/4.2.0/intel/120/lib/libnetcdff.a /opt/cray/netcdf/4.2.0/intel/120/lib/libnetcdf_intel.a /opt/cray/netcdf/4.2.0/intel/120/lib/libnetcdff_intel.a"
set INC_NETCDF    = "-I/opt/cray/netcdf/4.2.0/intel/120/include"
#
# User does not need to change anything below!
#
set type          = gds4p1_solo2                      # Name of the experiment
set root          = /autofs/na1_home1/Steve.Penny/godas4p1 #$cwd:h                            # The directory you created when you checkout
set fix_dir       = $root/fix
set code_dir      = $root/      #godas4p1                    # source code directory
set executable    = $root/exec_$platform/$type/fms_$type.x      # executable created after compilation
set pathnames     = $fix_dir/path_names_$type         # path to file containing list of source paths
set mppnccombine  = $root/bin/mppnccombine.$platform  # path to executable mppnccombine
set mkmfTemplate  = $root/bin/mkmf.template.$platform # path to template for your platform
set mkmf          = $root/bin/mkmf                    # path to executable mkmf
set MPI           = 1              # 1 if you want to use MPI, 0 otherwise
if($MPI) set cppDefs  = ( "-Duse_netCDF -Duse_netCDF3 -Duse_libMPI -DENABLE_GDS" )

#
# Users must ensure the correct environment file exists for their platform.
#
source $root/bin/environs.$platform  # environment variables and loadable modules

# setup directory structure
  if ( ! -d $executable:h )    mkdir -p $executable:h

#
# compile mppnccombine.c, needed only if $npes > 1
#NOTE: On some platforms you may need to specify the location for netcdf.h and libnetcdf.a
#      by modifying the following -I and -L
  if ( $MPI && ! -f $mppnccombine ) then
    cc -O -o $mppnccombine -I/usr/local/include $INC_NETCDF $code_dir/postprocessing/mppnccombine/mppnccombine.c $LIB_NETCDF
  endif

# set srcList = (  m4p1/src/mom4p1/drivers m4p1/src/mom4p1/ocean_diag m4p1/src/mom4p1/ocean_param/sources m4p1/src/mom4p1/ocean_param/mixing m4p1/src/mom4p1/ocean_param/gotm-4.0/include m4p1/src/mom4p1/ocean_param/gotm-4.0/turbulence m4p1/src/mom4p1/ocean_param/gotm-4.0/util m4p1/src/mom4p1/ocean_tracers  m4p1/src/ice_param m4p1/src/land_null )

# compile the model code and create executable
set makeFile      = Make_$type
  cd $executable:h
  $mkmf -f -m $makeFile -a $code_dir -t $mkmfTemplate -p $executable:t -c "$cppDefs" $pathnames $root/include $code_dir/shared/include $code_dir/shared/mpp/include 
#/usr/local/include

  make -f $makeFile

