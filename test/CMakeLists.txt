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
  set( CRTM_COEFFS_BRANCH "fix_REL-3.2.0.0" )  #Note switch to 3.2.0 when tarball ready

  set(CRTM_COEFFS_PATH ${CMAKE_BINARY_DIR}/test_data/${PROJECT_VERSION})
  file(MAKE_DIRECTORY ${CRTM_COEFFS_PATH})

  set(DOWNLOAD_BASE_URL "https://bin.ssec.wisc.edu/pub/s4/CRTM/")
  set(test_files_dirname ${CRTM_COEFFS_BRANCH}.tgz)
  set(checksum "5777242387228359869325e1a0505f85") # MD5SUM of fix_REL-3.2.0.0.tgz
  
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
  set(UNTAR_DEST "${CMAKE_BINARY_DIR}/test_data/${PROJECT_VERSION}")
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

set(_PARMIO_LUT_TESTINPUT "${CMAKE_CURRENT_BINARY_DIR}/testinput/PARMIO.MWwater.EmisCoeff.nc")

# test_PARMIO_TLAD
add_executable(test_PARMIO_TLAD mains/regression/parmio_tlad/test_PARMIO_TLAD.f90)
target_link_libraries(test_PARMIO_TLAD PRIVATE crtm)
if(PARMIO_LUT_PRESENT)
  add_test(NAME test_PARMIO_TLAD
           COMMAND $<TARGET_FILE:test_PARMIO_TLAD> "${_PARMIO_LUT_TESTINPUT}")
endif()

add_executable(test_PARMIO_Dispatcher mains/regression/parmio_tlad/test_PARMIO_Dispatcher.f90)
target_link_libraries(test_PARMIO_Dispatcher PRIVATE crtm)
if(PARMIO_LUT_PRESENT)
  add_test(NAME test_PARMIO_Dispatcher
           COMMAND $<TARGET_FILE:test_PARMIO_Dispatcher>
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                   "${_PARMIO_LUT_TESTINPUT}")
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

# Phase-4 lifecycle integration test: exercises CRTM_Init(PARMIOCoeff_File=...)
# plus dispatch routing and silent FASTEM fallback. Filename is passed as a
# basename relative to the testinput/ directory; CRTM_Init handles prepending.
add_executable(test_PARMIO_Lifecycle mains/regression/parmio_tlad/test_PARMIO_Lifecycle.f90)
target_link_libraries(test_PARMIO_Lifecycle PRIVATE crtm)

if(PARMIO_LUT_PRESENT)
  add_test(NAME test_PARMIO_FASTEM_DeltaSweep
           COMMAND $<TARGET_FILE:test_PARMIO_FASTEM_DeltaSweep>
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                   "${_PARMIO_LUT_TESTINPUT}")
  set_tests_properties(test_PARMIO_FASTEM_DeltaSweep
                       PROPERTIES
                         WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}")

  add_test(NAME test_PARMIO_Lifecycle
           COMMAND $<TARGET_FILE:test_PARMIO_Lifecycle>
                   "${CMAKE_CURRENT_BINARY_DIR}/testinput/"
                   "PARMIO.MWwater.EmisCoeff.nc")
  set_tests_properties(test_PARMIO_Lifecycle
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
endif()
unset(_PARMIO_LUT_TESTINPUT)

add_executable(test_PARMIO_AWS1_ObsSmoke mains/regression/parmio_tlad/test_PARMIO_AWS1_ObsSmoke.f90)
target_link_libraries(test_PARMIO_AWS1_ObsSmoke PRIVATE crtm)

add_executable(test_PARMIO_GMI_ObsSpace mains/regression/parmio_tlad/test_PARMIO_GMI_ObsSpace.f90)
target_link_libraries(test_PARMIO_GMI_ObsSpace PRIVATE crtm)

add_executable(test_PARMIO_ATMS_ClearOcean_ObsSpace mains/regression/parmio_tlad/test_PARMIO_ATMS_ClearOcean_ObsSpace.f90)
target_link_libraries(test_PARMIO_ATMS_ClearOcean_ObsSpace PRIVATE crtm)

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
# REL-3.2.0: NetCDF is the canonical LUT format. Symlinks point at .nc
# files in the netCDF/ subtrees of the test data. Two exceptions retain
# .bin symlinks because no NetCDF analog is published yet:
#   - SSU TauCoeff (no TauCoeff/ODPS/netCDF/ssu_n*.TauCoeff.nc on disk)
#   - Zeeman SSMIS TauCoeff (loaded by a path that hardcodes .bin)
# CRTM_Init's per-loader fallback (Resolve_Coeff_Format / SpcCoeff /
# TauCoeff) silently switches to .bin when an .nc file is missing, so
# missing entries here are recoverable for sensors that have a binary copy.
AerosolCoeff/netCDF/AerosolCoeff.nc
AerosolCoeff/netCDF/AerosolCoeff.GOCART-GEOS5.nc
AerosolCoeff/netCDF/AerosolCoeff.GOCART-GEOS5.BRC.kb.v2.nc
AerosolCoeff/netCDF/AerosolCoeff.CMAQ.nc
CloudCoeff/netCDF/CloudCoeff.nc
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2022.nc
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
SpcCoeff/netCDF/cris399_npp.NLTECoeff.nc
SpcCoeff/netCDF/iasi_metop-b.NLTECoeff.nc
SpcCoeff/netCDF/airs_aqua.NLTECoeff.nc
# cris-fsr_n21 is exercised by many regression/OMP tests above; its NLTECoeff
# lives only under the canonical NLTECoeff/netCDF/ tree (no SpcCoeff/netCDF/
# copy). CREATE_SYMLINK_FILENAME stages by basename, so the testinput/ symlink
# still lands next to cris-fsr_n21.SpcCoeff.nc for the flat-directory sibling
# lookup the test drivers rely on (File_Path=./testinput/).
NLTECoeff/netCDF/cris-fsr_n21.NLTECoeff.nc
SpcCoeff/netCDF/amsua_n19.ACCoeff.nc
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

# Symlink all CRTM files
# Version 3
CREATE_SYMLINK_FILENAME( ${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput
                         ${crtm_test_input} )
