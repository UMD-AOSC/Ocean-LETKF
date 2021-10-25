esma_set_this (OVERRIDE OCN.letkf)

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
  obs/params_obs.f90
  obs/compute_profile_error.f90
  obs/read_argo.f90
  obs/read_avhrr_pathfinder.f90
  obs/read_aviso_adt.f90
  obs/read_ice_txt.f90
  obs/gsw_pot_to_insitu.f90
# obs/
# for EOS from gsw-- this must be replaced - later, to be "in sync" with above choice of ${MODEL}
  obs/gsw_fortran_v3_03/gsw_oceanographic_toolbox.f90
  )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

# MAT There is a bug in common_mom4.f90 that cannot handle bounds-checking
#     For now build with release flags even in debug
set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE}")

esma_add_library (${this}
  SRCS ${SRCS}
  INCLUDES ${INC_NETCDF}
  DEPENDENCIES ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES}
  )

set_target_properties (${this} PROPERTIES COMPILE_OPTIONS ${PP})
target_compile_definitions(${this} PRIVATE DYNAMIC)

ecbuild_add_executable (
   TARGET OCN.letkf_${MODEL}.x 
   SOURCES letkf/letkf.f90 
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )

# different observation operators or H(.)
#adt
 ecbuild_add_executable (
    TARGET OCN.obsOp_${MODEL}.adt.x 
    SOURCES obs/obsop_adt.f90
    LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
    )
#sst
#ecbuild_add_executable (
#   TARGET OCN.obsOp_${MODEL}.sst.x 
#   SOURCES obs/obsop_sst.f90
#   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
#   )
#t-profile
 ecbuild_add_executable (
    TARGET OCN.obsOp_${MODEL}.tprof.x 
    SOURCES obs/obsop_tprof.f90
    LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
    )
#s-profile
 ecbuild_add_executable (
    TARGET OCN.obsOp_${MODEL}.sprof.x 
    SOURCES obs/obsop_sprof.f90
    LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
    )
#sea-ice concentration
 ecbuild_add_executable (
    TARGET OCN.obsOp_${MODEL}.icefrac.x 
    SOURCES obs/obsop_icefrac.f90
    LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
    )
