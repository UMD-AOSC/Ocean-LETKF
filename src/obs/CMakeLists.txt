# SA: esma needs to go out!
esma_set_this ()


set (SRCS
  compute_profile_error.f90 gsw_pot_to_insitu.f90 params_obs.f90
  obsop_adt.f90
  obsop_bufr_adt.f90
  obsop_bufr_sprof.f90
  obsop_bufr_sprof_hycom.f90
  obsop_bufr_sst.f90
  obsop_bufr_tprof.f90
  obsop_bufr_tprof_hycom.f90
  obsop_ecmwf_sprof.f90
  obsop_ecmwf_tprof.f90
  obsop_icefrac.f90
  obsop_ncoda_3d.f90*
  obsop_ncoda_ssh.f90*
  obsop_sla.f90
  obsop_sprof.f90
  obsop_sst.f90
  obsop_tprof.f90
  read_argo.f90 read_avhrr_pathfinder.f90 read_aviso_adt.f90 read_bufr.f90 read_bufr_hycom.f90 
  read_ecmwf_fdbk.f90 read_ice_txt.f90
#SA: not like this, only for _now_
  gsw_fortran_v3_03/gsw_check_functions.f90
  gsw_fortran_v3_03/gsw_oceanographic_toolbox.f90
 )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

#SA: is this flag needed?
set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE}")

esma_add_library (${this}
  SRCS ${SRCS}
  INCLUDES ${INC_NETCDF}
  DEPENDENCIES ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES}
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

#SA: depends on common_all, model_specific/mom6/, letkf and ?? 
#    doing this for obsop_adt.f90, but has to be repeated in some _nice_ way for other obsop_xx.f90 as well.
ecbuild_add_executable (
   TARGET obsop_adt.x
   SOURCES obsop_adt.f90
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )
