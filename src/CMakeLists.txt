#cmake_minimum_required (VERSION 3.18)

include(GNUInstallDirs)

# This is for selecting the MOM6 ocean model
## We default to MOM6 and look if something is passed in...
set (DEFAULT_MODEL "mom6")
if (NOT MODEL)
  set (MODEL ${DEFAULT_MODEL})
endif ()
set (MODEL ${MODEL} CACHE STRING "Building ocean LETKF with a set model")

## ... and make sure it is allowed
set (ALLOWED_MODEL "cice5" "hycom" "hycom_nrl" "mom4" "mom6" "nemo" "roms" "sis")
if (NOT MODEL IN_LIST ALLOWED_MODEL)
  message (FATAL_ERROR "The only allowed values for MODEL are ${ALLOWED_MODEL}")
endif ()

set (SRCS
# from common_all/
  common_all/SFMT.f90
  common_all/common.f90
  common_all/common_mpi.f90
  common_all/common_mtx.f90
  common_all/common_letkf.f90
  common_all/kdtree.f90
  common_all/netlib.f
  common_all/netlibblas.f
# from ../support/io/
  ../support/io/m_ncio.f90
  ../support/io/nc_rdatt.f90.inc
  ../support/io/nc_rdvar.f90.inc
  ../support/io/m_h5io.f90
  ../support/io/h5_rdatt.f90.inc
  ../support/io/h5_rdvar.f90.inc
#--
# from letkf/
  letkf/params_letkf.f90
  letkf/letkf_local.f90
  letkf/letkf_obs.f90
  letkf/letkf_tools.f90
  letkf/vars_letkf.f90
  letkf/vars_obs.f90
#--
# from model_specific/${MODEL}/
  model_specific/${MODEL}/vars_model.f90
  model_specific/${MODEL}/params_model.f90
  model_specific/${MODEL}/input_nml_${MODEL}.f90
  model_specific/${MODEL}/common_obs_${MODEL}.f90
  model_specific/${MODEL}/common_mpi_${MODEL}.f90
  model_specific/${MODEL}/common_${MODEL}.f90
  model_specific/${MODEL}/common_debug_${MODEL}.f90
#--
# from obs/
  obs/w3movdat_full.f
  obs/params_obs.f90
  obs/compute_profile_error.f90
  obs/read_argo.f90
  obs/read_avhrr_pathfinder.f90
  obs/read_geostationary.f90
  obs/read_viirs.f90
  obs/read_smap.f90
  obs/read_aviso_adt.f90
  obs/read_ice_txt.f90
  obs/gsw_pot_to_insitu.f90
# obs/
# for EOS from gsw-- this must be replaced - later, to be "in sync" with above choice of ${MODEL}
  obs/gsw_fortran_v3_03/gsw_oceanographic_toolbox.f90
  )

add_library(OCN.letkf ${SRCS})
if (${SOLO_BUILD})
    set (OCNLETKF_LD_LIBS NetCDF::NetCDF_Fortran MPI::MPI_Fortran HDF5::HDF5)
else()
    set (OCNLETKF_LD_LIBS NetCDF::NetCDF_Fortran MPI::MPI_Fortran HDF5::HDF5 Threads::Threads)
endif()
target_link_libraries(OCN.letkf PUBLIC ${OCNLETKF_LD_LIBS})
set_target_properties (OCN.letkf PROPERTIES Fortran_PREPROCESS ON) # Need to preprocess .f90 files for ifdef on next line
target_compile_definitions(OCN.letkf PRIVATE DYNAMIC MOM6)
install(TARGETS OCN.letkf)

## Programs
# letkf
add_executable(OCN.letkf_${MODEL}.x letkf/letkf.f90)
target_link_libraries(OCN.letkf_${MODEL}.x PRIVATE OCN.letkf)
install(TARGETS OCN.letkf_${MODEL}.x)

# different observation operators or H(.)
set(operators
  adt
  sst
  sst_geostationary
  sst_viirs
  sss
  tprof
  sprof
  icefrac
  )
foreach(operator ${operators})
  add_executable(OCN.obsOp_${MODEL}.${operator}.x obs/obsop_${operator}.f90)
  target_link_libraries(OCN.obsOp_${MODEL}.${operator}.x PRIVATE OCN.letkf NetCDF::NetCDF_Fortran HDF5::HDF5)
  set_target_properties (OCN.obsOp_${MODEL}.${operator}.x  PROPERTIES Fortran_PREPROCESS ON) 
  target_compile_definitions(OCN.obsOp_${MODEL}.${operator}.x PRIVATE DYNAMIC)
  install(TARGETS OCN.obsOp_${MODEL}.${operator}.x)
endforeach()

# different observation decoders
set(decoders
  sst_viirs
  sss_smap
  )
foreach(decoder ${decoders})
  add_executable(OCN.decoder.${decoder}.x obs/decoder_${decoder}.f90)
  target_link_libraries(OCN.decoder.${decoder}.x PRIVATE OCN.letkf NetCDF::NetCDF_Fortran HDF5::HDF5)
  set_target_properties (OCN.decoder.${decoder}.x  PROPERTIES Fortran_PREPROCESS ON) 
  target_compile_definitions(OCN.decoder.${decoder}.x PRIVATE DYNAMIC)
  install(TARGETS OCN.decoder.${decoder}.x)
endforeach()

# superob utils
add_executable(OCN.superob.x obs/superob.f90)
target_link_libraries(OCN.superob.x PRIVATE OCN.letkf NetCDF::NetCDF_Fortran HDF5::HDF5)
set_target_properties (OCN.superob.x  PROPERTIES Fortran_PREPROCESS ON) 
target_compile_definitions(OCN.superob.x PRIVATE DYNAMIC)
install(TARGETS OCN.superob.x)



