#!/bin/ksh --login
module load mpt
module load intel
module load netcdf/4.1.3-intel
module load nco

echo "LETKF Post-processing step"
echo "Processing cycle: ${YYYYMMDDHH}"
echo "I am member ${MEMBERID}"
workdir=${EXP_DATA}/${YYYYMMDDHH}/letkf_post/${MEMBERID}
mkdir -p ${workdir}
cd ${workdir}

echo "This is the LETKF post processing for member ${MEMBERID} for cycle ${YYYYMMDDHH}" > post.out

#STEVE: active code:
USE_IAU= #(False)
USE_TRI2SPH=0
MEM2=`printf %.2d ${MEMBERID}`
MEM3=`printf %.3d ${MEMBERID}`

TMPDIR=${EXP_DATA}/${YYYYMMDDHH}
IY=${YYYYMMDDHH:0:4}
IM=${YYYYMMDDHH:4:2}
ID=${YYYYMMDDHH:6:2}
IH=${YYYYMMDDHH:8:2}
IN=00
IS=00

# Update the ISLOT date
#date=/bin/date
#sinc=`expr $FCST - 1`
#sinc_units=days
#NY=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%Y`
#NM=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%m`
#ND=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%d`
#NH=`$date -d "$IY-$IM-$ID $sinc $sinc_units" +%H`
#IY=$NY
#IM=$NM
#ID=$ND
#IH=$NH
#echo "processing cycle as: $IY$IM$ID$IH"

#if [ $USE_IAU ]; then
#  #STEVE: if using IAU, use the original gues and apply IAU increments over the model run.
#  cp $OUTPUT/gues/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc $OUTPUT/anal/$MEM3/
#  cp $OUTPUT/gues/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  $OUTPUT/anal/$MEM3/
#  cp $OUTPUT/gues/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       $OUTPUT/anal/$MEM3/
#  if [ $USE_TRI2SPH ]; then  #STEVE: actually, this is SPH2TRI:
#    echo "ERROR: must add code to convert spherical increments to tripolar grid."
#    exit 1
#  fi
#else
if [ "$USE_TRI2SPH" -eq 1 ]; then  #STEVE: actually, this is SPH2TRI:
  src_grid="grid_spec_p5deg_sph_stripped.nc"
  dst_grid="grid_spec_mom4_stripped.nc"
  ln -fs $REGRID/$dst_grid .
  ln -fs $REGRID/$src_grid .
  ln -fs $TMPDIR/letkf/anal$MEM3.ocean_temp_salt.res.nc .
  ln -fs $TMPDIR/letkf/anal$MEM3.ocean_velocity.res.nc .
  ln -fs $TMPDIR/letkf/anal$MEM3.ocean_sbc.res.nc .

  # Must create an interpolated version in spherical grid format
  # First, set up the namelist
  # Next, run the regrid utility (postprocessing: tripolar to spherical)
  # For temp/salt
  cat >input.nml <<!
&regrid_3d_nml
  src_file = "anal$MEM3.ocean_temp_salt.res.nc",
  numfields = 2
  src_field_name = 'temp', 'salt'
  dest_grid = '$dst_grid',
  dest_file = "$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc",
  ntimes_saved = 1,
  timelevel_saved = 1,
  num_nbrs = 4 /
&fms_nml
  domains_stack_size = 395850 /
!
    mpiexec_mpt -np $PBS_NP $REGRID/regrid_3d_sph2tri.x

  # For u/v
  cat >input.nml <<!
&regrid_nml
  src_data = "anal$MEM3.ocean_velocity.res.nc",
  src_grid = '$src_grid',
  dst_grid = '$dst_grid',
  dst_data = "$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc",
  num_flds = 2
  fld_name = 'u','v'
  fld_pos  =  'C','C'
  vector_fld = .true., .true.
  use_source_vertical_grid = .false.
  apply_mask = .true.
  debug      = .false. /
!
    mpiexec_mpt -np $PBS_NP $REGRID/regrid_sph2tri.x

  # For sbc
#  if [ $USE_SFC ]; then
#
#         num_flds       = 6
#         fld_name       = 't_surf','s_surf','u_surf','v_surf','sea_lev','frazil'
#         fld_pos        =  'T','T','C','C','T','T'
#         vector_fld     = .false.,.false.,.true.,.true.,.false.,.false.

    cat >input.nml <<!
&regrid_nml
  src_data       = "anal$MEM3.ocean_sbc.res.nc",
  src_grid = '$src_grid',
  dst_grid = '$dst_grid',
  dst_data       = "$IY$IM$ID.$IH$IN$IS.ocean_sbc_ts.res.nc",
  num_flds       = 2
  fld_name       = 't_surf','s_surf'
  fld_pos        =  'T','T'
  vector_fld     = .false.
  debug          = .false. /
!
    mpiexec_mpt -np $PBS_NP $REGRID/regrid_sph2tri.x

    cat >input.nml <<!
&regrid_nml
  src_data       = "anal$MEM3.ocean_sbc.res.nc",
  src_grid = '$src_grid',
  dst_grid = '$dst_grid',
  dst_data       = "$IY$IM$ID.$IH$IN$IS.ocean_sbc_uv.res.nc",
  num_flds       = 2
  fld_name       = 'u_surf','v_surf'
  fld_pos        =  'C','C'
  vector_fld     = .true.,.true.
  debug          = .false. /
!
    mpiexec_mpt -np $PBS_NP $REGRID/regrid_sph2tri.x

    cat >input.nml <<!
&regrid_nml
  src_data       = "anal$MEM3.ocean_sbc.res.nc",
  src_grid = '$src_grid',
  dst_grid = '$dst_grid',
  dst_data       = "$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc",
  num_flds       = 2
  fld_name       = 'sea_lev','frazil'
  fld_pos        =  'T','T'
  vector_fld     = .false.
  debug          = .false. /
!
    mpiexec_mpt -np $PBS_NP $REGRID/regrid_sph2tri.x

    # Use nco command "ncks -h -A" to combine (append) these files
    ncks -h -A $IY$IM$ID.$IH$IN$IS.ocean_sbc_uv.res.nc $IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc
    ncks -h -A $IY$IM$ID.$IH$IN$IS.ocean_sbc_ts.res.nc $IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc
    # and remove the extra files:
#   rm -f $IY$IM$ID.$IH$IN$IS.ocean_sbc_uv.res.nc
#   rm -f $IY$IM$ID.$IH$IN$IS.ocean_sbc_ts.res.nc

    # Link the output files:
    ln -fs $IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc ${workdir}/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc
    ln -fs $IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc  ${workdir}/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc
    ln -fs $IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc       ${workdir}/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc

#  fi
  else
    if [ -f $TMPDIR/letkf/anal$MEM3.ocean_temp_salt.res.nc ]; then
      cp $TMPDIR/letkf/anal$MEM3.ocean_temp_salt.res.nc $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_temp_salt.res.nc
    fi
    if [ -f $TMPDIR/letkf/anal$MEM3.ocean_velocity.res.nc ]; then
      cp $TMPDIR/letkf/anal$MEM3.ocean_velocity.res.nc  $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_velocity.res.nc
    fi
#    if [ $USE_SFC ]; then
#      if [ -f $TMPDIR/letkf/anal$MEM3.ocean_sbc.res.nc ]; then
        cp $TMPDIR/letkf/anal$MEM3.ocean_sbc.res.nc     $OUTPUT/anal/$MEM3/$IY$IM$ID.$IH$IN$IS.ocean_sbc.res.nc
#      fi
#    fi
  fi
#fi

    #STEVE: Move the files to experiment output directory
#   cp ${workdir}/$IY$IM$ID.$IH$IN$IS.ocean_*.res.nc $OUTPUT/anal/$MEM3/
    #STEVE: NOTICE - I'm leaving these here for now.

exit 0
