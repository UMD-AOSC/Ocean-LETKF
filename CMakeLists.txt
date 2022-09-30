cmake_minimum_required(VERSION 3.18)
option(SOLO_BUILD "build ocean-letkf with local libs" OFF)

if (SOLO_BUILD) # build ocean-letkf with local libs. Most users should use this 

    project(ocean-letkf LANGUAGES Fortran)

    message(STATUS "OCEAN-LETKF: SOLO_BUILD=${SOLO_BUILD}, search local libs: MPI, NetCDF")

    # find MPI libs
    find_package(MPI COMPONENTS Fortran REQUIRED)

    # find HDF5 libs
    find_package(HDF5 COMPONENTS Fortran REQUIRED)

    # find NetCDF Fortran libs
    ## download additional CMakeModules for NetCDF if FindNetCDF.cmake is not found.
    ## this is done only once
    set(EXTERNAL_CMAKE_MODULE_URL "https://github.com/NOAA-EMC/CMakeModules.git")
    set(EXTERNAL_CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/support/CMakeModules/Modules)
    if (NOT EXISTS ${EXTERNAL_CMAKE_MODULE_PATH}/FindNetCDF.cmake)
       message(WARNING "OCEAN-LETKF: File (${EXTERNAL_CMAKE_MODULE_PATH}/FindNetCDF.cmake) not found.")
       message(STATUS "Pull additional CMakeModules from ${EXTERNAL_CMAKE_MODULE_URL}")

       find_package(Git REQUIRED)
       execute_process(COMMAND git clone -b release/public-v1 ${EXTERNAL_CMAKE_MODULE_URL}
                       WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/support) 
    endif()
    list(APPEND CMAKE_MODULE_PATH ${EXTERNAL_CMAKE_MODULE_PATH};${CMAKE_MODULE_PATH})
    find_package(NetCDF COMPONENTS Fortran REQUIRED)

    # compiler-dependent flags   
    if (CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        set(CMAKE_Fortran_FLAGS "-O2 -ffree-line-length-none -frecord-marker=4 -finit-local-zero -fbacktrace")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -ffree-line-length-none -frecord-marker=4 -finit-local-zero -fbacktrace -fcheck=mem,pointer,bounds -fno-realloc-lhs")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffree-line-length-none -frecord-marker=4 -finit-local-zero")
        ## additional flags for gfortran version >= 10.0
        if (${CMAKE_Fortran_COMPILER_VERSION} VERSION_GREATER_EQUAL 10) 
           set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -fallow-argument-mismatch -fallow-invalid-boz")
        endif ()
    elseif (CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
        set(CMAKE_Fortran_FLAGS "-O2 -g -traceback")
        set(CMAKE_Fortran_FLAGS_DEBUG "-O0 -traceback -check bounds")
        set(CMAKE_Fortran_FLAGS_RELEASE "-O3")
    endif()

else() # build as an application with GEOS-ESM infra

    message(STATUS "OCEAN-LETKF: SOLO_BUILD=${SOLO_BUILD}, build with GEOS-ESM infra")

endif()

add_subdirectory(src)
