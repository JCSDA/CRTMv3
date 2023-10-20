# CRTM
# Copyright and License: see LICENSE

cmake_minimum_required( VERSION 3.12 )
project( crtm VERSION 3.1.0 LANGUAGES Fortran )

option(OPENMP "Build crtm with OpenMP support" ON)

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

## Configuration options
include( ${PROJECT_NAME}_compiler_flags )
include(GNUInstallDirs)

## Dependencies
if(OPENMP)
  find_package( OpenMP COMPONENTS Fortran )

	# Set OMP_NUM_THREADS as a cached variable
	# probably need to have a consistent method for selecting this
	set(OMP_NUM_THREADS "8" CACHE STRING "Number of threads for OpenMP")

	# Export OMP_NUM_THREADS to the environment
	set(ENV{OMP_NUM_THREADS} ${OMP_NUM_THREADS})
	
endif()

if(DEFINED ENV{NETCDF_PATH})
  list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF_PATH})
elseif(DEFINED ENV{NETCDF})
  list(APPEND CMAKE_PREFIX_PATH $ENV{NETCDF})
endif()
find_package(NetCDF REQUIRED Fortran)


## Sources
add_subdirectory( src )
add_subdirectory( test )
include(CTest)

# Installation
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-config
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})

install(EXPORT ${PROJECT_NAME}-config
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Summary
# todo select status message based on build type selected
message(STATUS "Configuration summary")
message(STATUS "Project name : ${PROJECT_NAME}")
message(STATUS "Project version : ${PROJECT_VERSION}")
message(STATUS "Fortran compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler flags : ${CMAKE_Fortran_FLAGS}")
message(STATUS "Fortran compiler flags for release : ${CMAKE_Fortran_FLAGS_RELEASE}")
message(STATUS "Fortran compiler flags for debug : ${CMAKE_Fortran_FLAGS_DEBUG}")

set(crtm_LIBRARIES crtm)

