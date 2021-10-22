# SA: esma needs to go out!
esma_set_this ()


set (SRCS
  common_debug_mom6.f90 common_mom6.f90 common_mpi_mom6.f90 common_obs_mom6.f90 
  gdas_diagfile_io.f90 input_nml_mom6.f90 params_model.f90 vars_model.f90
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

# SA: this gives us a library that we use elsewhere.
