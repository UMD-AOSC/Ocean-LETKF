# SA: esma needs to go out!
esma_set_this ()


set (SRCS
   common.f90 common_letkf.f90 common_mpi.f90 common_mtx.f90 common_obs.f90
   kdtree.f90 netlibblas.f netlib.f SFMT.f90
  )

if (CMAKE_Fortran_COMPILER_ID MATCHES Intel)
   set( CMAKE_Fortran_FLAGS_RELEASE "${FOPT3} ${BIG_ENDIAN} ${BYTERECLEN} ${EXTENDED_SOURCE} ${ALIGNCOM}")
endif ()

#SA: following should is N/A here, since we do not use common_mom4.f90 -- but do we still have the same bug??
# MAT There is a bug in common_mom4.f90 that cannot handle bounds-checking
#     For now build with release flags even in debug
set( CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_RELEASE}")

esma_add_library (${this}
  SRCS ${SRCS}
  INCLUDES ${INC_NETCDF}
# SA: if we need mom6, keep it as following:
  DEPENDENCIES ${NETCDF_LIBRARIES} mom6 ${MPI_Fortran_LIBRARIES}
  )

target_compile_definitions (${this} PUBLIC _LAPACK_ gmao_intf HAVE_ESMF)

#SA: no exec from this dir, only a library that is used by others.
#ecbuild_add_executable (
#   TARGET oceanda.x 
#   SOURCES letkf.f90 
#   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
#   )
#
#ecbuild_add_executable (
#   TARGET oceanobsop.x 
#   SOURCES obsop.f90 
#   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
#   )
#
#ecbuild_add_executable (
#   TARGET oceanobs_nc2bin.x 
#   SOURCES oceanobs_nc2bin.f90 
#   LIBS ${this} ${NETCDF_LIBRARIES} ${MPI_Fortran_LIBRARIES} Threads::Threads
#   )
