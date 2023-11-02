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
		set(OMP_NUM_THREADS "4")
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
add_subdirectory(test)
include(CTest)

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
        DESTINATION .)
#add_dependencies(${PROJECT_NAME} install)

# Install the export set for use with the install-tree
set(BUILD_DIR ${CMAKE_BINARY_DIR})
install(EXPORT ${PROJECT_NAME}-config
        FILE ${PROJECT_NAME}-targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${BUILD_DIR})
#add_dependencies(${PROJECT_NAME} install)

#add_custom_target(default COMMAND ${CMAKE_MAKE_PROGRAM} install)
#add_custom_target(my_all COMMAND ${CMAKE_MAKE_PROGRAM} crtm && ${CMAKE_MAKE_PROGRAM} install)
#add_custom_target(my_all ALL COMMAND ${CMAKE_MAKE_PROGRAM} install)

#add_custom_target(install_all
#                  COMMAND ${CMAKE_COMMAND} --build . --target all
#                  COMMAND ${CMAKE_COMMAND} --build . --target install
#                  WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

#add_custom_target(run_install ALL
#                  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target install
#                  DEPENDS ${PROJECT_NAME})

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
