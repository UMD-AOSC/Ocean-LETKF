# SA: esma needs to go out!
esma_set_this ()


set (SRCS
  letkf.f90 letkf_local.f90 letkf_obs.f90 letkf_tools.f90 params_letkf.f90
  vars_letkf.f90 vars_obs.f90
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

#SA: depends on ..?
ecbuild_add_executable (
   TARGET ocean_letkf.x 
   SOURCES letkf.f90 
   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
   )
