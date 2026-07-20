# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required (VERSION 3.12)
project("CRTM_Tests" VERSION 3.2.0 LANGUAGES Fortran C)

enable_testing ()


message (STATUS "Building tests for CRTM v${PROJECT_VERSION}.") 

list( APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake )
set( CMAKE_DIRECTORY_LABELS ${PROJECT_NAME} )

# macro to create a symlink from src to dst
function(CREATE_SYMLINK src dst)
    foreach (FILENAME ${ARGN})
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${FILENAME} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK)

# macro to create a symlink from src to dst with just filename
function(CREATE_SYMLINK_FILENAME src dst)
    foreach (FILENAME ${ARGN})
        get_filename_component(filename ${FILENAME} NAME )
        execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
            ${src}/${FILENAME}
            ${dst}/${filename} )
        endforeach(FILENAME)
endfunction(CREATE_SYMLINK_FILENAME)

function(ADD_CRTM_UNIT_TEST target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE crtm)
  add_test(NAME test_${target}
           COMMAND $<TARGET_FILE:${target}>)
endfunction()

# Create Data directory for test input config and symlink all files

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/adjoint)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/tangent_linear)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/k_matrix)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/forward)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/unit)


if ( DEFINED ENV{CRTM_TEST_ROOT})
  set(CRTM_TEST_ROOT "$ENV{CRTM_TEST_ROOT}")
else()
  set(CRTM_TEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
endif()

if(NOT FIX_FILE_PATH)
  if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES})
    set(FIX_FILE_PATH "$ENV{LOCAL_PATH_JEDI_TESTFILES}/fix")
  else()
    set(FIX_FILE_PATH ${CMAKE_SOURCE_DIR}/fix)
  endif()
endif()

# If local path to testfiles is defined don't download 
IF(EXISTS  ${FIX_FILE_PATH})
  get_filename_component(CRTM_COEFFS_PATH "${FIX_FILE_PATH}/.." ABSOLUTE)
  message(STATUS "Using existing fix files: ${CRTM_COEFFS_PATH}/fix")
  set( CRTM_COEFFS_BRANCH "" )

else() # Download CRTM coefficients
  set( CRTM_COEFFS_BRANCH_PREFIX "" )  #preserves the structure of the paths that have been used in jedi previously vs the local path above
  set( CRTM_COEFFS_BRANCH "fix_REL-3.2.0.0" )

  set(CRTM_COEFFS_PATH ${CMAKE_SOURCE_DIR}/test-data-release)
  file(MAKE_DIRECTORY ${CRTM_COEFFS_PATH})

  set(DOWNLOAD_BASE_URL "https://bin.ssec.wisc.edu/pub/s4/CRTM/")
  set(test_files_dirname ${CRTM_COEFFS_BRANCH}.tgz)
  set(checksum "3dcef94c129efb78c85cdf542fca55ae") # MD5SUM of fix_REL-3.2.0.0.tgz
  
  # Check if the CRTM binary tarball is already present otherwise download it.
  if (DEFINED ENV{CRTM_BINARY_FILES_TARBALL})
    file(MD5 "$ENV{CRTM_BINARY_FILES_TARBALL}" computed_checksum )
    if (NOT "${computed_checksum}" STREQUAL "${checksum}")
      message(FATAL_ERROR "Environment provided CRTM_BINARY_FILES_TARBALL=$ENV{CRTM_BINARY_FILES_TARBALL} had checksum '${computed_checksum}' not matching expected '${checksum}'")
    endif()
    message(STATUS "Using CRTM coeffs tarball file CRTM_BINARY_FILES_TARBALL=$ENV{CRTM_BINARY_FILES_TARBALL}")
    set(DOWNLOAD_DEST "$ENV{CRTM_BINARY_FILES_TARBALL}")
  else()
    set(DOWNLOAD_URL "${DOWNLOAD_BASE_URL}/${test_files_dirname}")
    set(DOWNLOAD_DEST "${CRTM_COEFFS_PATH}/${test_files_dirname}")
    message(STATUS "Downloading CRTM coeffs files from: ${DOWNLOAD_URL} to ${DOWNLOAD_DEST}")
    file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_DEST}
      EXPECTED_MD5 ${checksum}
      SHOW_PROGRESS)
  endif()

  # -- UNTAR -- 
  
  # Define the directory to untar into
  set(UNTAR_DEST "${CMAKE_SOURCE_DIR}/test-data-release")
    message(STATUS "Checking if ${UNTAR_DEST}/${CRTM_COEFFS_BRANCH} already exists...")
  # Only untar if the CHECK_PATH does not exist
  if(NOT EXISTS ${UNTAR_DEST}/${CRTM_COEFFS_BRANCH}/)
    message(STATUS "Untarring the downloaded file (~2 minutes) to ${UNTAR_DEST}")
    # Use the built-in tar feature to work around an issue on Orion with the
    # system tar command. See: https://github.com/JCSDA/spack-stack/issues/1355
    execute_process(COMMAND ${CMAKE_COMMAND} -E chdir ${UNTAR_DEST}
                    ${CMAKE_COMMAND} -E tar xzf ${DOWNLOAD_DEST}
                    RESULT_VARIABLE result
                    OUTPUT_QUIET ERROR_QUIET)
    if(NOT result EQUAL 0)
      message(FATAL_ERROR "Failed to untar the file.")
    endif()
  else()
    message(STATUS "Files have already been untarred. Skipping.")
  endif()
endif()

  set (CRTM_TESTFILES_PATH ${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH}/fix)
  message(STATUS "Downloaded CRTM files are placed at ${CRTM_TESTFILES_PATH} ")
  set_property (GLOBAL PROPERTY CRTM_TESTFILES_PATH ${CRTM_TESTFILES_PATH})

# Create list of sensor ids for testing
list( APPEND Simple_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
        amsua_n19
)

list( APPEND ScatteringSwitch_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)

list( APPEND User_Emissivity_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)

list( APPEND SOI_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18 
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)

list( APPEND VerticalCoordinates_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18 
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)

list( APPEND SSU_Sensor_Ids
        ssu_n06
        ssu_n14
)

list( APPEND ClearSky_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)

# Create list of sensor ids for testing
list( APPEND AOD_Sensor_Ids
        cris-fsr_n21
        v.abi_g18
        cris399_npp
        v.abi_gr
        abi_g18
        airs_aqua
)

# Create list of sensor ids for testing
list( APPEND Zeeman_Sensor_Ids
      ssmis_f20
      ssmis_f19
      ssmis_f18
      ssmis_f17
      ssmis_f16
)


list( APPEND ChannelSubset_Sensor_Ids
      iasi_metop-b
)

list( APPEND Aircraft_Sensor_Ids
      cris-fsr_n21
      crisB1_npp
)

list( APPEND Downwelling_Radiance_Sensor_Ids
        atms_n21
        cris-fsr_n21
        v.abi_g18
        atms_npp
        cris399_npp
        v.abi_gr
        abi_g18
        modis_aqua
)


list( APPEND OMPoverChannels_Sensor_Ids
        atms_n21
)

# OMP_Speedup uses cris-fsr_n21 (~431 channels) so the channel-level OMP
# parallelism in CRTM_Forward has enough work to demonstrate clear speedup.
list( APPEND OMP_Speedup_Sensor_Ids
        cris-fsr_n21
)

# OMP_Consistency: Forward + K_Matrix must be bit-identical across thread counts.
# A spread of paths: MW radiometer, IR imager, and a hyperspectral IR sounder
# (the latter exercises the per-channel NLTE/Zeeman predictor handling).
list( APPEND OMP_Consistency_Sensor_Ids
        atms_n21
        abi_g18
        cris399_npp
)

# ChannelSubset_OMP: like OMP_Consistency, but applies several channel subsets
# first (sparse / front-loaded / split) so the channel-thread chunking and the
# inactive-channel bookkeeping in CRTM_Forward / _K_Matrix get exercised under a
# thread-count sweep. iasi_metop-b gives a large channel count for meaningful
# chunk sizes, cris399_npp adds the NLTE per-channel path, and atms_n21 covers
# the small-channel-count clamps. (Follow-up coverage for JCSDA/CRTMv3#164.)
list( APPEND ChannelSubset_OMP_Sensor_Ids
        iasi_metop-b
        cris399_npp
        atms_n21
)

list (APPEND common_tests
     Simple
     AOD
     Zeeman
     ChannelSubset
     ClearSky
     Aircraft
     Downwelling_Radiance
     ScatteringSwitch
     User_Emissivity
     SOI
     SSU
     VerticalCoordinates
)

list (APPEND omp_tests
     OMPoverChannels
     OMP_Speedup
     OMP_Consistency
     ChannelSubset_OMP
)


list (APPEND regression_types
     forward
     k_matrix
)

list (APPEND TLAD_types
     adjoint
     tangent_linear
)

list (APPEND TLAD_tests
     Simple
     ClearSky
)


#####################################################################

add_executable(check_tropics mains/application/check_tropics.f90)
target_link_libraries(check_tropics PRIVATE crtm)

add_executable(test_check_crtm mains/application/check_crtm.F90)
target_link_libraries(test_check_crtm PRIVATE crtm)
add_test(NAME test_check_crtm 
         COMMAND test_check_crtm)

add_executable(Unit_TL_TEST mains/unit/Unit_Test/test_TL.f90)
target_link_libraries(Unit_TL_TEST PRIVATE crtm)
add_test(NAME test_Unit_TL_TEST
         COMMAND $<TARGET_FILE:Unit_TL_TEST>)

add_executable(Unit_Aerosol_Bypass mains/unit/Unit_Test/test_Aerosol_Bypass.f90)
target_link_libraries(Unit_Aerosol_Bypass PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass
         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass>)

add_executable(Unit_AerosolScatter_TL mains/unit/Unit_Test/test_AerosolScatter_TL.f90)
target_link_libraries(Unit_AerosolScatter_TL PRIVATE crtm)
add_test(NAME test_Unit_AerosolScatter_TL
         COMMAND $<TARGET_FILE:Unit_AerosolScatter_TL>)

add_executable(Unit_AerosolScatter_AD mains/unit/Unit_Test/test_AerosolScatter_AD.f90)
target_link_libraries(Unit_AerosolScatter_AD PRIVATE crtm)
add_test(NAME test_Unit_AerosolScatter_AD
         COMMAND $<TARGET_FILE:Unit_AerosolScatter_AD>)

add_executable(Unit_AerosolScatter_K mains/unit/Unit_Test/test_AerosolScatter_K.f90)
target_link_libraries(Unit_AerosolScatter_K PRIVATE crtm)
add_test(NAME test_Unit_AerosolScatter_K
         COMMAND $<TARGET_FILE:Unit_AerosolScatter_K>)

add_executable(Unit_Aerosol_Bypass_TL mains/unit/Unit_Test/test_Aerosol_Bypass_TL.f90)
target_link_libraries(Unit_Aerosol_Bypass_TL PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_TL
	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_TL>)

add_executable(Unit_Aerosol_Bypass_adjoint mains/unit/Unit_Test/test_Aerosol_Bypass_adjoint.f90)
target_link_libraries(Unit_Aerosol_Bypass_adjoint PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_adjoint
	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_adjoint>)

add_executable(Unit_Aerosol_Bypass_k_matrix mains/unit/Unit_Test/test_Aerosol_Bypass_k_matrix.f90)
target_link_libraries(Unit_Aerosol_Bypass_k_matrix PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_k_matrix
         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_k_matrix>)

# MW land surface Jacobian validation (issue #281): analytic LAI / Vegetation_Fraction
# Jacobians (K-matrix adjoint and tangent-linear) vs finite differences.
add_executable(Unit_Land_Jacobian mains/unit/Unit_Test/test_Land_Jacobian.f90)
target_link_libraries(Unit_Land_Jacobian PRIVATE crtm)
add_test(NAME test_Unit_Land_Jacobian
         COMMAND $<TARGET_FILE:Unit_Land_Jacobian>)

# Legacy Fastem1 MW-water SST (Water_Temperature) Jacobian vs finite differences.
add_executable(Unit_Fastem1_SST_Jacobian mains/unit/Unit_Test/test_Fastem1_SST_Jacobian.f90)
target_link_libraries(Unit_Fastem1_SST_Jacobian PRIVATE crtm)
add_test(NAME test_Unit_Fastem1_SST_Jacobian
         COMMAND $<TARGET_FILE:Unit_Fastem1_SST_Jacobian>)

# Pins the v3.2.0 DDA-ARTS ICE_CLOUD scattering change (ICE_CLOUD now scatters via
# the full branch / IconCloudIce habit). Registered only when a DDA-ARTS CloudCoeff
# is present in the staged fix tree.
set(_CC_DDA_FIX "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2024.nc")
add_executable(test_DDA_ICE_CLOUD_Forward mains/unit/Unit_Test/test_DDA_ICE_CLOUD_Forward.f90)
target_link_libraries(test_DDA_ICE_CLOUD_Forward PRIVATE crtm)
if(EXISTS "${_CC_DDA_FIX}")
  add_test(NAME test_DDA_ICE_CLOUD_Forward
           COMMAND $<TARGET_FILE:test_DDA_ICE_CLOUD_Forward>)
else()
  message(STATUS "DDA-ARTS CloudCoeff not found at ${_CC_DDA_FIX}; test_DDA_ICE_CLOUD_Forward not registered.")
endif()
unset(_CC_DDA_FIX)

# Issue #238: coefficient file paths must survive a long File_Path without being
# truncated by fixed-length buffers. Initializes CRTM through a ~300-character
# symlinked path that exceeds the old 80/128/256 caps.
add_executable(Unit_Long_Path_Init mains/unit/Unit_Test/test_Long_Path_Init.f90)
target_link_libraries(Unit_Long_Path_Init PRIVATE crtm)
add_test(NAME test_Unit_Long_Path_Init
         COMMAND $<TARGET_FILE:Unit_Long_Path_Init>)

# General TL-vs-FD and AD consistency check across sensor types (from #280).
# Parameterized by sensor id and mode: "fd" checks the tangent-linear against a
# central finite difference of the forward (per channel), "ad" checks the TL-AD
# dot-product adjoint identity. Perturbs atmospheric temperature. Covers the
# three main sensor types (all present in the common suite): atms_n21 (MW
# sounder), cris-fsr_n21 (IR sounder), v.abi_g18 (IR imager).
add_executable(test_FD_consistency mains/unit/Unit_Test/test_FD_consistency.f90)
target_link_libraries(test_FD_consistency PRIVATE crtm)
foreach(_fdc_sensor atms_n21 cris-fsr_n21 v.abi_g18)
  add_test(NAME test_FD_consistency_${_fdc_sensor} COMMAND test_FD_consistency ${_fdc_sensor} fd)
  add_test(NAME test_AD_consistency_${_fdc_sensor} COMMAND test_FD_consistency ${_fdc_sensor} ad)
endforeach()

# Surface downwelling radiance (RTSolution%Down_Radiance) TL/AD/K correctness:
# baseline-independent finite-difference, adjoint dot-product, and K-vs-AD checks
# for a TOA-radiance control plus the surface downwelling output.
add_executable(Unit_Downwelling_TLADK mains/unit/Unit_Test/test_Downwelling_TLADK.f90)
target_link_libraries(Unit_Downwelling_TLADK PRIVATE crtm)
add_test(NAME test_Unit_Downwelling_TLADK
         COMMAND $<TARGET_FILE:Unit_Downwelling_TLADK>)

# Multi-sensor single-call consistency: one CRTM_Forward/TL/K_Matrix call with
# ChannelInfo(1:2) must match the per-sensor calls bit-for-bit (guards the
# cumulative 'ln' channel-offset bookkeeping in the OpenMP sensor/channel loops,
# which no single-sensor test can see).
add_executable(Unit_MultiSensor_SingleCall mains/unit/Unit_Test/test_MultiSensor_SingleCall.f90)
target_link_libraries(Unit_MultiSensor_SingleCall PRIVATE crtm)
add_test(NAME test_Unit_MultiSensor_SingleCall
         COMMAND $<TARGET_FILE:Unit_MultiSensor_SingleCall>)

# Performance micro-benchmark for the opt-in level-resolved radiance profiles
# (Forward + K_Matrix timing, overcast ADA). Not a ctest; run manually.
add_executable(bench_Profile_Perf mains/unit/Unit_Test/bench_Profile_Perf.f90)
target_link_libraries(bench_Profile_Perf PRIVATE crtm)

# Surface netCDF I/O round-trip (REL-3.2.0 baseline-format conversion):
# write a rank-2 Surface(L x M) with NetCDF=.TRUE., read it back, and verify
# every field round-trips exactly.
add_executable(Unit_Surface_netCDF_io mains/unit/Unit_Test/test_Surface_netCDF_io.f90)
target_link_libraries(Unit_Surface_netCDF_io PRIVATE crtm)
add_test(NAME test_Unit_Surface_netCDF_io
         COMMAND $<TARGET_FILE:Unit_Surface_netCDF_io>)

# Atmosphere netCDF I/O round-trip (REL-3.2.0 baseline-format conversion):
# write a rank-2 Atmosphere(L x M) with clouds+aerosols and NetCDF=.TRUE.,
# read it back, and verify every serialized field round-trips exactly.
add_executable(Unit_Atmosphere_netCDF_io mains/unit/Unit_Test/test_Atmosphere_netCDF_io.f90)
target_link_libraries(Unit_Atmosphere_netCDF_io PRIVATE crtm)
add_test(NAME test_Unit_Atmosphere_netCDF_io
         COMMAND $<TARGET_FILE:Unit_Atmosphere_netCDF_io>)

#SpcCoeff utilities
list (APPEND SCoeff_Utils
  SpcCoeff_Edit
  SpcCoeff_Inspect
  SpcCoeff_BIN2NC
  SpcCoeff_NC2BIN
)
foreach(testtype IN LISTS SCoeff_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/SpcCoeff/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#TauCoeff ODAS utilities
list (APPEND ODAS_Utils
  ODAS_BIN2NC
  ODAS_NC2BIN
)
foreach(testtype IN LISTS ODAS_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/ODAS/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#TauCoeff ODPS utilities
list (APPEND ODPS_Utils
  ODPS_BIN2NC
  ODPS_NC2BIN
)
foreach(testtype IN LISTS ODPS_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/ODPS/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#TauCoeff ODSSU utilities (SSU container, BIN -> NC only)
list (APPEND ODSSU_Utils
  ODSSUBIN2NC
)
foreach(testtype IN LISTS ODSSU_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/ODSSU/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#TauCoeff general utilities
list (APPEND TCoeff_Utils
  TauCoeff_Inspect
)
foreach(testtype IN LISTS TCoeff_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#CloudCoeff utilities 
list (APPEND CCoef_Utils
  CloudCoeff_NC2BIN
  CloudCoeff_BIN2NC
  CloudCoeff_Inspect
)
foreach(testtype IN LISTS CCoef_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/CloudCoeff/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

#AerosolCoeff utilities
list (APPEND ACoeff_Utils
  AerosolCoeff_Inspect
)
foreach(testtype IN LISTS ACoeff_Utils)
	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/AerosolCoeff/${testtype}/${testtype}.f90)
	target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

list (APPEND SEcategoryCoef_Utils
  SEcategory_NC2BIN
  SEcategory_BIN2NC
)
foreach(testtype IN LISTS SEcategoryCoef_Utils)
        add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/EmisCoeff/SEcategory/${testtype}/${testtype}.f90)
        target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

list (APPEND IRsnowCoef_Utils
#  IRsnowCoeff_BIN2NC  !** doesn't exist yet (BTJ)
  IRsnowCoeff_NC2BIN
)
foreach(testtype IN LISTS IRsnowCoef_Utils)
        add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/EmisCoeff/IR_Snow/${testtype}/${testtype}.f90)
        target_link_libraries(${testtype} PRIVATE crtm)
endforeach()

list (APPEND IRwaterCoef_Utils
  IRwaterCoeff_BIN2NC
  IRwaterCoeff_NC2BIN
)
foreach(testtype IN LISTS IRwaterCoef_Utils)
        add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/EmisCoeff/IR_Water/${testtype}/${testtype}.f90)
        target_link_libraries(${testtype} PRIVATE crtm)
endforeach()



#---------------------------------------------------------------------------------
#unit tests

#first upper level Unit_Test

# REL-3.2.0: previously-binary IO tests now exercise NetCDF (the new default).
# The dedicated _NC sibling tests have been removed because they are now
# duplicative; if a binary regression test is needed in the future, add
# it back as an explicit case.
list(APPEND io_tests
  Spc_IO mains/unit/input_output/test_SpcCoeff/test_spc_io.f90
  AerosolCoeff_IO mains/unit/input_output/test_AerosolCoeff/test_aerosol_coeff_io.f90
  CloudCoeff_IO mains/unit/input_output/test_CloudCoeff/test_cloud_coeff_io.f90
  EmisCoeff_IO mains/unit/input_output/test_EmisCoeff/test_emis_coeff_io.f90
)

list(LENGTH io_tests io_tests_len)
math(EXPR io_tests_max "${io_tests_len} - 1")
foreach(io_test_index RANGE 0 ${io_tests_max} 2)
  list(GET io_tests ${io_test_index} io_test_name)
  math(EXPR io_test_source_index "${io_test_index} + 1")
  list(GET io_tests ${io_test_source_index} io_test_source)
  ADD_CRTM_UNIT_TEST(${io_test_name} ${io_test_source})
endforeach()

#(PS) Hypsometric Equation test
add_executable(hypsometric_eq mains/unit/Unit_Test/test_Hypsometric.f90)
target_link_libraries(hypsometric_eq PRIVATE crtm)
add_test(NAME test_hypsometric_eq
         COMMAND $<TARGET_FILE:hypsometric_eq>)

# test_MWwater_io: removed in REL-3.2.0. The deprecated CRTM_MWwaterCoeff_Load
# function it exercised reads a binary-only EmisCoeff file; with the LUT-based
# FASTEM path retired in favor of CRTM_MWwaterCoeff_Load_FASTEM (no file IO),
# this comparison test no longer has an active counterpart to compare against.

# test_MWSurfEM
add_executable(test_MWSurfEM mains/unit/Unit_Test/test_MWSurfEM.f90)
target_link_libraries(test_MWSurfEM PRIVATE crtm)
add_test(NAME test_MWSurfEM
         COMMAND $<TARGET_FILE:test_MWSurfEM>)

# test_CONST_MIXED_Polarization: exercises the CONST_MIXED_POLARIZATION (=13)
# surface-optics V/H mixing in CRTM_Compute_SfcOptics. Pol type 13 is used only
# by TMS (TROPICS / tomorrow.io) sensors, none of which are in the common
# regression suite, so this is the dedicated coverage for that code path. It
# needs only the TMS SpcCoeff and the (already-staged) FASTEM6 MW-water
# EmisCoeff; CRTM_Compute_SfcOptics reads no TauCoeff. Registered only when the
# TMS SpcCoeff is present in the staged fix tree.
set(_TMS_SPC "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/SpcCoeff/netCDF/tms_tropics-01.SpcCoeff.nc")
add_executable(test_CONST_MIXED_Polarization mains/unit/Unit_Test/test_CONST_MIXED_Polarization.f90)
target_link_libraries(test_CONST_MIXED_Polarization PRIVATE crtm)
if(EXISTS "${_TMS_SPC}")
  set(TMS_COEFFS_PRESENT TRUE)
  add_test(NAME test_CONST_MIXED_Polarization
           COMMAND $<TARGET_FILE:test_CONST_MIXED_Polarization>)
else()
  set(TMS_COEFFS_PRESENT FALSE)
  message(STATUS "TMS SpcCoeff (tms_tropics-01.SpcCoeff.nc) not found; test_CONST_MIXED_Polarization will not be registered.")
endif()
unset(_TMS_SPC)

# A full forward run for TROPICS (used by the >=200 GHz PARMIO delta sweep
# below) additionally needs the TMS TauCoeff, which CONST_MIXED does not.
set(_TMS_TAU "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/TauCoeff/ODPS/netCDF/tms_tropics-01.TauCoeff.nc")
set(_TMS_SPC "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/SpcCoeff/netCDF/tms_tropics-01.SpcCoeff.nc")
if(EXISTS "${_TMS_SPC}" AND EXISTS "${_TMS_TAU}")
  set(TMS_TAU_PRESENT TRUE)
else()
  set(TMS_TAU_PRESENT FALSE)
endif()
unset(_TMS_TAU)
unset(_TMS_SPC)

# test_TELSEM2_MWland
# Stage the TELSEM2 atlas under two names in the shared testinput:
#   * TELSEM2.MWland.test.nc       -- explicit-opt-in path (MWlandCoeff_File).
#   * TELSEM2.MWland.EmisCoeff.nc  -- the DEFAULT drop-in name, staged ON PURPOSE
#     so the test can prove the opt-in gate: with the default-named atlas present
#     on the coefficient path, CRTM_Init must still use NESDIS_LandEM unless the
#     caller opts in (Use_MWland_Atlas=.TRUE. or MWlandCoeff_File). Because the
#     atlas is now opt-in rather than presence-activated, staging the default name
#     is inert for every other land test -- and any regression of the gate would
#     surface as those tests changing. Registered only when the atlas is present
#     in the staged fix tree.
set(_TELSEM2_FIX "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/EmisCoeff/MW_Land/netCDF/TELSEM2.MWland.EmisCoeff.nc")
add_executable(test_TELSEM2_MWland mains/unit/Unit_Test/test_TELSEM2_MWland.f90)
target_link_libraries(test_TELSEM2_MWland PRIVATE crtm)
if(EXISTS "${_TELSEM2_FIX}")
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                   "${_TELSEM2_FIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/TELSEM2.MWland.test.nc" )
  execute_process( COMMAND ${CMAKE_COMMAND} -E create_symlink
                   "${_TELSEM2_FIX}"
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/TELSEM2.MWland.EmisCoeff.nc" )
  add_test(NAME test_TELSEM2_MWland
           COMMAND $<TARGET_FILE:test_TELSEM2_MWland>)
else()
  message(STATUS "TELSEM2 atlas not found at ${_TELSEM2_FIX}; test_TELSEM2_MWland will not be registered.")
endif()
unset(_TELSEM2_FIX)

# Check for canonical PARMIO LUT and AWS coefficients inside the staged
# test_data tree. Tests are registered only when the files are present at
# their canonical paths; nothing outside build/test_data/** is referenced.
set(_PARMIO_FIX "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix")

if(EXISTS "${_PARMIO_FIX}/EmisCoeff/MW_Water/netCDF/PARMIO.MWwater.EmisCoeff.nc")
  set(PARMIO_LUT_PRESENT TRUE)
else()
  message(STATUS "PARMIO LUT not found at ${_PARMIO_FIX}/EmisCoeff/MW_Water/netCDF/PARMIO.MWwater.EmisCoeff.nc; PARMIO regression tests will not be registered.")
  set(PARMIO_LUT_PRESENT FALSE)
endif()

if(EXISTS "${_PARMIO_FIX}/SpcCoeff/netCDF/mwr_aws.SpcCoeff.nc"
   AND EXISTS "${_PARMIO_FIX}/TauCoeff/ODPS/netCDF/mwr_aws.TauCoeff.nc")
  set(AWS_COEFFS_PRESENT TRUE)
else()
  message(STATUS "AWS coefficients (mwr_aws.{Spc,Tau}Coeff.nc) not found in ${_PARMIO_FIX}; AWS PARMIO test will not be registered.")
  set(AWS_COEFFS_PRESENT FALSE)
endif()
unset(_PARMIO_FIX)

# test_CloudCoeff_Exp_Forward: end-to-end coverage of the experimental
# ('CRTM-Exp') cloud-optics scheme. Runs the mwr_aws MW sensor with
# Cloud_Model='CRTM-Exp' and the complete 6-habit experimental LUT
# (CloudCoeff_Exp_Full6.nc), asserting a physical graupel-scattering TB
# depression that grows with water content. The legacy cloud schemes have no
# such coverage. Registered only when both the AWS coeffs and the experimental
# LUT are present in the staged fix tree.
set(_CC_EXP_FIX "${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix/CloudCoeff/netCDF/CloudCoeff_Exp_Full6.nc")
if(EXISTS "${_CC_EXP_FIX}")
  set(CLOUDCOEFF_EXP_PRESENT TRUE)
else()
  set(CLOUDCOEFF_EXP_PRESENT FALSE)
  message(STATUS "Experimental LUT (CloudCoeff_Exp_Full6.nc) not found at ${_CC_EXP_FIX}; test_CloudCoeff_Exp_Forward will not be registered.")
endif()
unset(_CC_EXP_FIX)
add_executable(test_CloudCoeff_Exp_Forward mains/unit/Unit_Test/test_CloudCoeff_Exp_Forward.f90)
target_link_libraries(test_CloudCoeff_Exp_Forward PRIVATE crtm)
if(AWS_COEFFS_PRESENT AND CLOUDCOEFF_EXP_PRESENT)
  add_test(NAME test_CloudCoeff_Exp_Forward
           COMMAND $<TARGET_FILE:test_CloudCoeff_Exp_Forward>)
else()
  message(STATUS "test_CloudCoeff_Exp_Forward not registered (needs mwr_aws coeffs + CloudCoeff_Exp_Full6.nc).")
endif()

# test_VectorRT_TLADK: TL/AD/K correctness for the vector-RT (n_Stokes>1)
# cloud-scattering path. The n_Stokes>1 ADA branch needs a >=6-phase-element
# LUT (the CRTM-Exp scheme), so it is unreachable by every other test; this is
# its only Jacobian coverage. Verifies TL-vs-central-FD on both Stokes
# components (d/dWater_Content for the polarized phase chain, d/dTemperature
# for the thermal source), the full-Stokes adjoint dot-product, K-vs-AD, and
# an n_Stokes=1 scalar control on the same scene. Same gating as the Exp
# forward test.
add_executable(test_VectorRT_TLADK mains/unit/Unit_Test/test_VectorRT_TLADK.f90)
target_link_libraries(test_VectorRT_TLADK PRIVATE crtm)
if(AWS_COEFFS_PRESENT AND CLOUDCOEFF_EXP_PRESENT)
  add_test(NAME test_VectorRT_TLADK
           COMMAND $<TARGET_FILE:test_VectorRT_TLADK>)
else()
  message(STATUS "test_VectorRT_TLADK not registered (needs mwr_aws coeffs + CloudCoeff_Exp_Full6.nc).")
endif()

# test_MW_O3_TLAD: TL/AD/K parity for the MW scene-ozone ODPS component
# (GROUP_MW_O3, Group_Index=7). Verifies TL-vs-central-FD for O3/H2O/T column
# perturbations, the adjoint dot-product over T+H2O+O3, and K-vs-AD on the
# most O3-sensitive channel. Self-adapting: against a 2-component group-3
# TauCoeff (the stock testinput link) it asserts the O3 response is
# identically zero (backward-compatibility control); against a group-7 file
# it exercises the ozone predictor TL/AD/K blocks end-to-end.
add_executable(test_MW_O3_TLAD mains/unit/Unit_Test/test_MW_O3_TLAD.f90)
target_link_libraries(test_MW_O3_TLAD PRIVATE crtm)
if(AWS_COEFFS_PRESENT)
  add_test(NAME test_MW_O3_TLAD
           COMMAND $<TARGET_FILE:test_MW_O3_TLAD>)
else()
  message(STATUS "test_MW_O3_TLAD not registered (needs mwr_aws coeffs).")
endif()

# test_Grazing_SfcOptics: regression test for the catastrophic-reflectivity
# guard in the MW-water surface optics (CRTM_FastemX and CRTM_PARMIO). Drives
# both backends directly at grazing zenith angles + 325 GHz (where the reflection
# correction extrapolates to ~1e35 without the guard) and asserts reflectivity
# stays bounded and the guard fires. Uses the AWS coeffs (mwr_aws init loads the
# FASTEM and PARMIO surface coefficients). Registered only when both the AWS
# coeffs and the PARMIO LUT are present.
add_executable(test_Grazing_SfcOptics mains/unit/Unit_Test/test_Grazing_SfcOptics.f90)
target_link_libraries(test_Grazing_SfcOptics PRIVATE crtm)
if(AWS_COEFFS_PRESENT AND PARMIO_LUT_PRESENT)
  add_test(NAME test_Grazing_SfcOptics
           COMMAND $<TARGET_FILE:test_Grazing_SfcOptics>)
else()
  message(STATUS "test_Grazing_SfcOptics not registered (needs mwr_aws coeffs + PARMIO LUT).")
endif()

set(_PARMIO_LUT_TESTINPUT "${CMAKE_CURRENT_BINARY_DIR}/testinput/PARMIO.MWwater.EmisCoeff.nc")

# test_PARMIO_TLAD
add_executable(test_PARMIO_TLAD mains/regression/parmio_tlad/test_PARMIO_TLAD.f90)
target_link_libraries(test_PARMIO_TLAD PRIVATE crtm)
if(PARMIO_LUT_PRESENT)
  add_test(NAME test_PARMIO_TLAD
           COMMAND $<TARGET_FILE:test_PARMIO_TLAD> "${_PARMIO_LUT_TESTINPUT}")
endif()

# test_PARMIO_M_Group_Probe (build-only diagnostic; not registered with CTest)
add_executable(test_PARMIO_M_Group_Probe mains/regression/parmio_tlad/test_PARMIO_M_Group_Probe.f90)
target_link_libraries(test_PARMIO_M_Group_Probe PRIVATE crtm)

# test_PARMIO_TLAD_RefValues (build-only helper to refresh test_PARMIO_TLAD's expected_e)
add_executable(test_PARMIO_TLAD_RefValues mains/regression/parmio_tlad/test_PARMIO_TLAD_RefValues.f90)
target_link_libraries(test_PARMIO_TLAD_RefValues PRIVATE crtm)

# test_PARMIO_FASTEM_VH_Sweep (build-only V/H emissivity comparison at 89/166/183/325 GHz)
add_executable(test_PARMIO_FASTEM_VH_Sweep mains/regression/parmio_tlad/test_PARMIO_FASTEM_VH_Sweep.f90)
target_link_libraries(test_PARMIO_FASTEM_VH_Sweep PRIVATE crtm)

# test_PARMIO_RC_Residual: build-only diagnostic. The driver compares CRTM's
# atmosphere-on TB (Kirchhoff (1-e)*Mod sky-reflection via FASTEM RCCoeff)
# against PARMIO standalone Tb.f's atmosphere-on output (specular Rvv0 sky
# reflection). Those are different atmospheric conventions, so the residual
# is structural rather than a defect on either side. Not registered as a
# CTest gate; runnable manually for diagnosis.
add_executable(test_PARMIO_RC_Residual mains/regression/parmio_tlad/test_PARMIO_RC_Residual.f90)
target_link_libraries(test_PARMIO_RC_Residual PRIVATE crtm)

add_executable(test_PARMIO_FASTEM_DeltaSweep mains/regression/parmio_tlad/test_PARMIO_FASTEM_DeltaSweep.f90)
target_link_libraries(test_PARMIO_FASTEM_DeltaSweep PRIVATE crtm)

if(PARMIO_LUT_PRESENT)
  add_test(NAME test_PARMIO_FASTEM_DeltaSweep
           COMMAND $<TARGET_FILE:test_PARMIO_FASTEM_DeltaSweep>
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                   "${_PARMIO_LUT_TESTINPUT}")
  set_tests_properties(test_PARMIO_FASTEM_DeltaSweep
                       PROPERTIES
                         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

  # Same target, AWS sensor invocation. Reaches 325 GHz, well above ATMS's 183 GHz.
  if(AWS_COEFFS_PRESENT)
    add_test(NAME test_PARMIO_FASTEM_DeltaSweep_AWS
             COMMAND $<TARGET_FILE:test_PARMIO_FASTEM_DeltaSweep>
                     "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                     "${_PARMIO_LUT_TESTINPUT}"
                     "mwr_aws")
    set_tests_properties(test_PARMIO_FASTEM_DeltaSweep_AWS
                         PROPERTIES
                           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()

  # Same target, TROPICS sensor invocation. Its 204.8 GHz channel is the
  # stronger PARMIO witness: unlike the AWS 325 GHz channels (which sit on the
  # H2O line and are surface-blind over moist atmospheres, so PARMIO-vs-FASTEM
  # delta_tb ~ 0.1 K), 204.8 GHz is window-like and PARMIO moves brightness
  # temperature by ~1-3 K -- a clear, regression-detectable surface signal.
  if(TMS_TAU_PRESENT)
    add_test(NAME test_PARMIO_FASTEM_DeltaSweep_TMS
             COMMAND $<TARGET_FILE:test_PARMIO_FASTEM_DeltaSweep>
                     "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                     "${_PARMIO_LUT_TESTINPUT}"
                     "tms_tropics-01")
    set_tests_properties(test_PARMIO_FASTEM_DeltaSweep_TMS
                         PROPERTIES
                           WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")
  endif()
endif()
unset(_PARMIO_LUT_TESTINPUT)

add_executable(test_PARMIO_AWS1_ObsSmoke mains/regression/parmio_tlad/test_PARMIO_AWS1_ObsSmoke.f90)
target_link_libraries(test_PARMIO_AWS1_ObsSmoke PRIVATE crtm)

add_executable(test_PARMIO_GMI_ObsSpace mains/regression/parmio_tlad/test_PARMIO_GMI_ObsSpace.f90)
target_link_libraries(test_PARMIO_GMI_ObsSpace PRIVATE crtm)

#=================================================================================
#forward and k_matrix regression tests
foreach(regtype IN LISTS regression_types)
  string(COMPARE EQUAL ${regtype} "k_matrix" isregtype)
  foreach(testtype IN LISTS common_tests)
    string(COMPARE EQUAL ${testtype} "Aircraft" istesttype)
    if (isregtype AND istesttype)
       continue() #skip Aircraft type for k_Matrix
    endif()
    string(COMPARE EQUAL ${testtype} "Downwelling_Radiance" istesttype)
    if (isregtype AND istesttype)
       continue() #skip Downwelling_Radiance type for k_Matrix
    endif()

    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.f90)
    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)

    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
    endforeach()
  endforeach()
endforeach()

#---------------------------------------------------------------------------------
#TLAD Regression tests
foreach(regtype IN LISTS TLAD_types)
  foreach(testtype IN LISTS TLAD_tests)
    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.f90)
    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)

    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
    endforeach()
  endforeach()
endforeach()

#---------------------------------------------------------------------------------
# >=200 GHz PARMIO default-path stored-reference regression (issue #311)
#
# No common-suite sensor reaches the 200 GHz PARMIO dispatch threshold, so the
# default MW-water ocean-emissivity path for mwr_aws (4 channels at ~325 GHz) is
# otherwise covered only by self-consistency drivers (test_PARMIO_TLAD, the
# FASTEM delta sweep) with no stored-reference truth file. Register the ClearSky
# driver for mwr_aws across forward/k_matrix/adjoint/tangent_linear so a
# numerical drift in the PARMIO default path is caught by the same
# RTSolution_Compare truth-file mechanism as every other sensor. ClearSky
# (n_Clouds = n_Aerosols = 0) isolates the surface/PARMIO contribution. Guarded
# on AWS_COEFFS_PRESENT AND PARMIO_LUT_PRESENT: without the LUT the dispatch
# falls back to FASTEM and this would not exercise PARMIO. The
# test_{regtype}_test_ClearSky executables are already defined by the loops
# above (ClearSky is in both common_tests and TLAD_tests).
if(AWS_COEFFS_PRESENT AND PARMIO_LUT_PRESENT)
  foreach(regtype IN LISTS regression_types TLAD_types)
    add_test(NAME test_${regtype}_ClearSky_mwr_aws
             COMMAND $<TARGET_FILE:test_${regtype}_test_ClearSky> "mwr_aws")
  endforeach()
endif()

#=================================================================================
#OpenMP regression tests
foreach(regtype IN LISTS regression_types)
  string(COMPARE EQUAL ${regtype} "k_matrix" isregtype)
  if (isregtype)
    continue() #skip all k_matrix tests
  endif()
  foreach(testtype IN LISTS omp_tests)
    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.F90)
    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)

    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
        # The speedup test asserts a wall-clock threading speedup (>= 1.20x);
        # under `ctest -j` the host is oversubscribed by sibling tests and the
        # measurement flakes. Run it with the machine to itself.
        if(testtype STREQUAL "OMP_Speedup")
          set_tests_properties(test_${regtype}_${testtype}_${sensor_id}
                               PROPERTIES RUN_SERIAL TRUE)
        endif()
    endforeach()
  endforeach()
endforeach()


#---------------------------------------------------------------------------------
# Add test_Active_Sensor executable
add_executable(test_Active_Sensor mains/unit/Unit_Test/test_active_sensor.f90)
target_link_libraries(test_Active_Sensor crtm)

# Add test for test_Active_Sensor
add_test(NAME test_Active_Sensor COMMAND test_Active_Sensor)

# Add test_AD_Active_Sensor executable
add_executable(test_AD_Active_Sensor mains/unit/Unit_Test/test_AD_Active_Sensor.f90)
target_link_libraries(test_AD_Active_Sensor crtm)

# Add test for test_AD_Active_Sensor
add_test(NAME test_AD_Active_Sensor COMMAND test_AD_Active_Sensor)

# Add test_TL_convergence_active_sensor executable
add_executable(test_TL_convergence_active_sensor mains/unit/Unit_Test/test_TL_convergence_active_sensor.f90)
target_link_libraries(test_TL_convergence_active_sensor crtm)

# Add test for test_TL_convergence_active_sensor
add_test(NAME test_TL_convergence_active_sensor COMMAND test_TL_convergence_active_sensor)


#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
# REL-3.2.0: NetCDF is the canonical LUT format; every coefficient symlink
# below points at a .nc file (SSU TauCoeff was converted via ODSSUBIN2NC and
# the Zeeman SSMIS zssmis_*.TauCoeff.nc set is complete). The
# fix_REL-3.2.0.0 tarball ships NO .bin coefficient files, so none are
# staged and the per-loader netCDF->Binary fallback (Resolve_Coeff_Format /
# SpcCoeff / TauCoeff) is not exercised by this suite: the binary
# coefficient READ path is untested at the suite level for this release.
AerosolCoeff/netCDF/AerosolCoeff.nc
AerosolCoeff/netCDF/AerosolCoeff.GOCART-GEOS5.nc
AerosolCoeff/netCDF/AerosolCoeff.GOCART-GEOS5.BRC.kb.v2.nc
AerosolCoeff/netCDF/AerosolCoeff.CMAQ.nc
CloudCoeff/netCDF/CloudCoeff.nc
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2022.nc
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2024.nc
BeCoeff/netCDF/BeCoeff.nc
EmisCoeff/MW_Water/netCDF/FASTEM6.MWwater.EmisCoeff.nc
EmisCoeff/MW_Water/netCDF/FASTEM4.MWwater.EmisCoeff.nc
EmisCoeff/IR_Ice/SEcategory/netCDF/NPOESS.IRice.EmisCoeff.nc
EmisCoeff/IR_Land/SEcategory/netCDF/NPOESS.IRland.EmisCoeff.nc
EmisCoeff/IR_Snow/SEcategory/netCDF/NPOESS.IRsnow.EmisCoeff.nc
EmisCoeff/IR_Snow/Nalli/netCDF/Nalli.IRsnow.EmisCoeff.nc
EmisCoeff/IR_Snow/Nalli/netCDF/Nalli2.IRsnow.EmisCoeff.nc
EmisCoeff/VIS_Ice/SEcategory/netCDF/NPOESS.VISice.EmisCoeff.nc
EmisCoeff/VIS_Land/SEcategory/netCDF/NPOESS.VISland.EmisCoeff.nc
EmisCoeff/VIS_Snow/SEcategory/netCDF/NPOESS.VISsnow.EmisCoeff.nc
EmisCoeff/VIS_Water/SEcategory/netCDF/NPOESS.VISwater.EmisCoeff.nc
EmisCoeff/IR_Water/netCDF/Nalli.IRwater.EmisCoeff.nc
EmisCoeff/IR_Water/netCDF/Nalli2.IRwater.EmisCoeff.nc
EmisCoeff/IR_Land/SEcategory/netCDF/USGS.IRland.EmisCoeff.nc
EmisCoeff/VIS_Land/SEcategory/netCDF/USGS.VISland.EmisCoeff.nc
EmisCoeff/VIS_Snow/SNICAR/netCDF/SNICAR.VISsnow.EmisCoeff.nc
SpcCoeff/netCDF/hirs4_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4_metop-a.TauCoeff.nc
SpcCoeff/netCDF/amsua_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n19.TauCoeff.nc
SpcCoeff/netCDF/amsua_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_metop-a.TauCoeff.nc
SpcCoeff/netCDF/amsua_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_aqua.TauCoeff.nc
SpcCoeff/netCDF/gmi_gpm.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gmi_gpm.TauCoeff.nc
SpcCoeff/netCDF/seviri_m08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/seviri_m08.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_n21.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_npp.TauCoeff.nc
SpcCoeff/netCDF/iasi_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_metop-a.TauCoeff.nc
SpcCoeff/netCDF/iasi_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_metop-b.TauCoeff.nc
SpcCoeff/netCDF/mhs_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_n19.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g15.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g15.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g15.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g15.TauCoeff.nc
SpcCoeff/netCDF/airs_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airs_aqua.TauCoeff.nc
SpcCoeff/netCDF/modis_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modis_aqua.TauCoeff.nc
SpcCoeff/netCDF/atms_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_n21.TauCoeff.nc
# NLTECoeff / ACCoeff sibling files for sensors that need them.
# SpcCoeff_netCDF_ReadFile loads <sensor>.{NLTE,AC}Coeff.nc from the same
# directory as the SpcCoeff file when present; without these, NLTE-active
# IR channels (e.g. cris399_npp 4.3 um band) and AMSU/MHS antenna
# corrections silently default to zero.
NLTECoeff/netCDF/cris399_npp.NLTECoeff.nc
NLTECoeff/netCDF/iasi_metop-b.NLTECoeff.nc
NLTECoeff/netCDF/airs_aqua.NLTECoeff.nc
# cris-fsr_n21 is exercised by many regression/OMP tests above; its NLTECoeff
# lives only under the canonical NLTECoeff/netCDF/ tree (no SpcCoeff/netCDF/
# copy). CREATE_SYMLINK_FILENAME stages by basename, so the testinput/ symlink
# still lands next to cris-fsr_n21.SpcCoeff.nc for the flat-directory sibling
# lookup the test drivers rely on (File_Path=./testinput/).
NLTECoeff/netCDF/cris-fsr_n21.NLTECoeff.nc
ACCoeff/netCDF/amsua_n19.ACCoeff.nc
SpcCoeff/netCDF/v.viirs-m_j2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/v.viirs-m_j2.TauCoeff.nc
SpcCoeff/netCDF/v.abi_g18.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_g18.TauCoeff.nc
SpcCoeff/netCDF/abi_g18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_g18.TauCoeff.nc
SpcCoeff/netCDF/cris399_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris399_npp.TauCoeff.nc
SpcCoeff/netCDF/crisB1_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB1_npp.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-m_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-m_npp.TauCoeff.nc
SpcCoeff/netCDF/v.abi_gr.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_gr.TauCoeff.nc
TauCoeff/ODPS/netCDF/zssmis_f20.TauCoeff.nc
TauCoeff/ODPS/netCDF/zssmis_f19.TauCoeff.nc
TauCoeff/ODPS/netCDF/zssmis_f18.TauCoeff.nc
TauCoeff/ODPS/netCDF/zssmis_f17.TauCoeff.nc
TauCoeff/ODPS/netCDF/zssmis_f16.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f20.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f19.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f18.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f17.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f16.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f20.SpcCoeff.nc
SpcCoeff/netCDF/ssmis_f18.SpcCoeff.nc
SpcCoeff/netCDF/ssmis_f19.SpcCoeff.nc
SpcCoeff/netCDF/ssmis_f16.SpcCoeff.nc
SpcCoeff/netCDF/ssmis_f17.SpcCoeff.nc
# SSU TauCoeff: ODSSU container in netCDF, generated from the legacy
# REL-3.1.2.0 binaries via src/.../ODSSU/ODSSUBIN2NC.
TauCoeff/ODPS/netCDF/ssu_n06.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n07.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n08.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n09.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n11.TauCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n14.TauCoeff.nc
SpcCoeff/netCDF/ssu_n06.SpcCoeff.nc
SpcCoeff/netCDF/ssu_n07.SpcCoeff.nc
SpcCoeff/netCDF/ssu_n08.SpcCoeff.nc
SpcCoeff/netCDF/ssu_n09.SpcCoeff.nc
SpcCoeff/netCDF/ssu_n11.SpcCoeff.nc
SpcCoeff/netCDF/ssu_n14.SpcCoeff.nc
)

# Locally-staged coefficients (PARMIO LUT, AWS) live in test_data/.../fix/
# under canonical paths (see PARMIO/AWS staging block above) and are added to
# the symlink list only when available.
if(PARMIO_LUT_PRESENT)
  list(APPEND crtm_test_input
       EmisCoeff/MW_Water/netCDF/PARMIO.MWwater.EmisCoeff.nc)
endif()
if(AWS_COEFFS_PRESENT)
  list(APPEND crtm_test_input
       SpcCoeff/netCDF/mwr_aws.SpcCoeff.nc
       TauCoeff/ODPS/netCDF/mwr_aws.TauCoeff.nc)
endif()
# Experimental cloud-optics LUT for test_CloudCoeff_Exp_Forward; staged from the
# tarball's fix tree so testinput/CloudCoeff_Exp_Full6.nc resolves there.
if(CLOUDCOEFF_EXP_PRESENT)
  list(APPEND crtm_test_input
       CloudCoeff/netCDF/CloudCoeff_Exp_Full6.nc)
endif()
# TMS SpcCoeff for test_CONST_MIXED_Polarization (polarization type 13).
if(TMS_COEFFS_PRESENT)
  list(APPEND crtm_test_input
       SpcCoeff/netCDF/tms_tropics-01.SpcCoeff.nc)
endif()
# TMS TauCoeff for the >=200 GHz PARMIO delta sweep (full forward run).
if(TMS_TAU_PRESENT)
  list(APPEND crtm_test_input
       TauCoeff/ODPS/netCDF/tms_tropics-01.TauCoeff.nc)
endif()

# Symlink all CRTM files
# Version 3
CREATE_SYMLINK_FILENAME( ${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput
                         ${crtm_test_input} )
