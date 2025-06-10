cmake_minimum_required(VERSION 3.12)

project(crtm VERSION 3.1.1 LANGUAGES Fortran)

option(OPENMP "Build crtm with OpenMP support" ON)
option(FIX_FILE_PATH "Path to fix files (default: fix/)" OFF)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake)
set(CMAKE_DIRECTORY_LABELS ${PROJECT_NAME})

# Include GNUInstallDirs to get the standard installation directory variables
include(GNUInstallDirs)

# Set custom installation paths for crtm only (MODIFIED SECTION)
if(CMAKE_INSTALL_PREFIX_INITIALIZED_TO_DEFAULT)
    set(CRTM_INSTALL_PREFIX "${CMAKE_BINARY_DIR}/lib" CACHE PATH "crtm install path" FORCE)
else()
    set(CRTM_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX})
endif()

message("CRTM_INSTALL_PREFIX set to ${CRTM_INSTALL_PREFIX}")

# Set output directories (unchanged from original)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)

# Function to get the current git commit hash (unchanged from original)
function(get_git_hash result)
    execute_process(
        COMMAND git rev-parse --short HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_HASH
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    set(${result} ${GIT_HASH} PARENT_SCOPE)
endfunction()

# Get the git hash (unchanged from original)
get_git_hash(PROJECT_HASH)

# Define ANSI color code for green (unchanged from original)
string(ASCII 27 ESCAPE)
set(RESET "${ESCAPE}[0m")
set(GREEN "${ESCAPE}[32m")

# Output the project name, version, and hash in green text (unchanged from original)
message(STATUS "${GREEN}[${PROJECT_NAME}] (${PROJECT_VERSION}) [${PROJECT_HASH}]${RESET}")

# Output the various install/output paths (unchanged from original)
message("CMAKE_INSTALL_PREFIX set to ${CMAKE_INSTALL_PREFIX}")
message("CMAKE_INSTALL_LIBDIR set to ${CMAKE_INSTALL_LIBDIR}")
message("CMAKE_LIBRARY_OUTPUT_DIRECTORY set to ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")

## Configuration options (unchanged from original)
include(${PROJECT_NAME}_compiler_flags)
include(GNUInstallDirs)

## Dependencies (unchanged from original)
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

## Sources (unchanged from original)
add_subdirectory(src)

# Generate and configure the project config file (unchanged from original)
configure_file(
  "${CMAKE_CURRENT_SOURCE_DIR}/crtm-config.cmake.in"
  "${CMAKE_CURRENT_BINARY_DIR}/crtm-config.cmake"
  @ONLY
)

# Installation for crtm (MODIFIED SECTION)
install(TARGETS ${PROJECT_NAME}
        EXPORT ${PROJECT_NAME}-config
        ARCHIVE DESTINATION ${CRTM_INSTALL_PREFIX}/lib
        LIBRARY DESTINATION ${CRTM_INSTALL_PREFIX}/lib)

# Generate the version file for the config package (unchanged from original)
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/crtm-config-version.cmake"
  VERSION ${PROJECT_VERSION}
  COMPATIBILITY AnyNewerVersion
)

# Install the generated or configured files (unchanged from original)
install(FILES 
              "${CMAKE_CURRENT_BINARY_DIR}/crtm-config.cmake"
              "${CMAKE_CURRENT_BINARY_DIR}/crtm-config-version.cmake"
        DESTINATION ${CRTM_INSTALL_PREFIX}/cmake/${PROJECT_NAME})

# Install the export set for use with the install-tree (unchanged, with modified destination)
install(EXPORT ${PROJECT_NAME}-config
	FILE "${PROJECT_NAME}-targets.cmake"
	DESTINATION ${CRTM_INSTALL_PREFIX}/cmake/${PROJECT_NAME})

# For build-tree linkage (unchanged from original)
export(EXPORT ${PROJECT_NAME}-config
	FILE "${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}-targets.cmake"
)

# Additional settings for test and build type (unchanged from original)
set(CRTM_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

include(CTest)
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# Summary (unchanged from original)
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
    message(STATUS "Fortran compiler flags for release : ${CMAKE_Fortran_FLAGS_RELEASE}")
elseif(CMAKE_BUILD_TYPE_UPPER STREQUAL "DEBUG")
    message(STATUS "Fortran compiler flags for debug : ${CMAKE_Fortran_FLAGS_DEBUG}")
else()
    message(STATUS "No build type specified, defaulting to Release.")
endif()

#export the CRTM_SOURCE_DIR for use in the test framework (unchanged from original)

