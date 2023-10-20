project(crtm VERSION 3.1.0 LANGUAGES Fortran)

option(OPENMP "Build crtm with OpenMP support" ON)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

## Configuration options
include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

## Dependencies
if(OPENMP)
  find_package(OpenMP COMPONENTS Fortran)
  set(OMP_NUM_THREADS "8" CACHE STRING "Number of threads for OpenMP")
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
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})

# Install the export set for use with the install-tree
install(EXPORT ${PROJECT_NAME}-config
        FILE ${PROJECT_NAME}-targets.cmake
        NAMESPACE ${PROJECT_NAME}::
        DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME})


install( EXPORT ${PROJECT_NAME}-targets
         DESTINATION "${INSTALL_CMAKE_DIR}" )

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
