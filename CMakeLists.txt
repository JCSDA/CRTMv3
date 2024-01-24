cmake_minimum_required( VERSION 3.12 )

project(crtm VERSION 3.1.0 LANGUAGES Fortran)

option(OPENMP "Build crtm with OpenMP support" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})
#if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
set(CMAKE_INSTALL_PREFIX ${CMAKE_CURRENT_BINARY_DIR} CACHE PATH "default install path" FORCE)
#endif()
message("CMAKE_INSTALL_PREFIX set to ${CMAKE_INSTALL_PREFIX}")

## Configuration options
include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

## Dependencies
if(OPENMP)
  find_package(OpenMP COMPONENTS Fortran)

	# Check if OMP_NUM_THREADS is set in the environment
	if(DEFINED ENV{OMP_NUM_THREADS})
		set(OMP_NUM_THREADS $ENV{OMP_NUM_THREADS})
	else()
		# Set a default value if not set
		set(OMP_NUM_THREADS "1")
	endif()
	
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
add_subdirectory(src)

# Generate and configure the project config file
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/crtm-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/crtm-config.cmake"
  @ONLY
)

# Installation
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-config
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
#add_dependencies(${PROJECT_NAME} install)

# Generate the version file for the config package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/crtm-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install the generated or configured files
install(FILES "${CMAKE_CURRENT_BINARY_DIR}/crtm-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/crtm-config-version.cmake"
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
#add_dependencies(${PROJECT_NAME} install)



# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}-config
        FILE "${PROJECT_NAME}-targets.cmake"
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


export(EXPORT ${PROJECT_NAME}-config
        FILE "${PROJECT_NAME}-targets.cmake"
        NAMESPACE ${PROJECT_NAME}::)

add_subdirectory(test)
include(CTest)


# Summary
message(STATUS "Configuration summary")
message(STATUS "Project name : ${PROJECT_NAME}")
message(STATUS "Project version : ${PROJECT_VERSION}")
message(STATUS "Fortran compiler : ${CMAKE_Fortran_COMPILER}")
message(STATUS "Fortran compiler flags : ${CMAKE_Fortran_FLAGS}")
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

string(TOUPPER "${CMAKE_BUILD_TYPE}" CMAKE_BUILD_TYPE_UPPER)
message(STATUS "Build type : ${CMAKE_BUILD_TYPE}")
if(CMAKE_BUILD_TYPE_UPPER STREQUAL "RELEASE")
    # Commands for Release build
		message(STATUS "Fortran compiler flags for release : ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
    # Commands for Debug build
		message(STATUS "Fortran compiler flags for debug : ${CMAKE_Fortran_FLAGS_DEBUG}")
else()
    # Commands for other build types
		message(STATUS "No build type specified, defaulting to Release.")
endif()




set(crtm_LIBRARIES crtm)
