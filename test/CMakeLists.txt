# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required (VERSION 3.12)
project("CRTM_Tests" VERSION 3.1.0 LANGUAGES Fortran C)

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

# Create Data directory for test input config and symlink all files

file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/testinput)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/adjoint)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/tangent_linear)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/k_matrix)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/forward)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/results/unit)


if( DEFINED ENV{LOCAL_PATH_JEDI_TESTFILES})
  set(LOCAL_PATH_JEDI_TESTFILES "$ENV{LOCAL_PATH_JEDI_TESTFILES}")
endif()

if ( DEFINED ENV{CRTM_TEST_ROOT})
  set(CRTM_TEST_ROOT "$ENV{CRTM_TEST_ROOT}")
else()
  set(CRTM_TEST_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
endif()



# If local path to testfiles is defined don't download 
IF(EXISTS  ${CMAKE_SOURCE_DIR}/fix)
  set( CRTM_COEFFS_PATH ${CMAKE_SOURCE_DIR}/fix )
  message(STATUS "use LOCAL_PATH_JEDI_TESTFILES: ${CRTM_COEFFS_PATH}")
  message("Using existing local fix directory instead of downloading.")
  set( CRTM_COEFFS_BRANCH "" )

else() # Download CRTM coefficients
  set( CRTM_COEFFS_BRANCH_PREFIX "" )  #preserves the structure of the paths that have been used in jedi previously vs the local path above
  set( CRTM_COEFFS_BRANCH "fix_REL-3.1.0.1_01252024" ) #this version is the CRTM version, not the jedi / skylab version (again, following prior installations -- I'm not tied to this in any way)
  
  set(CRTM_COEFFS_PATH ${CMAKE_BINARY_DIR}/test_data/${PROJECT_VERSION})
  file(MAKE_DIRECTORY ${CRTM_COEFFS_PATH})
  
  set(DOWNLOAD_BASE_URL "https://bin.ssec.wisc.edu/pub/s4/CRTM/")
  set(test_files_dirname ${CRTM_COEFFS_BRANCH}.tgz) 
  set(checksum "b37602ea463be55451de1e61ff80b2ba") # MD5SUM of fix_REL-3.1.0.1_01252024.tgz
  
  set(DOWNLOAD_URL "${DOWNLOAD_BASE_URL}/${test_files_dirname}")
  set(DOWNLOAD_DEST "${CRTM_COEFFS_PATH}/${test_files_dirname}")
  
  # Download during CMake configure step
  message(STATUS "Downloading CRTM coeffs files from: ${DOWNLOAD_URL} to ${DOWNLOAD_DEST}")
  file(DOWNLOAD ${DOWNLOAD_URL} ${DOWNLOAD_DEST}
    EXPECTED_MD5 ${checksum}
    SHOW_PROGRESS)

  # -- UNTAR -- 
  
  # Define the directory to untar into
  set(UNTAR_DEST "${CMAKE_BINARY_DIR}/test_data/${PROJECT_VERSION}")
	message(STATUS "Checking if ${UNTAR_DEST}/${CRTM_COEFFS_BRANCH} already exists...")
  # Only untar if the CHECK_PATH does not exist
  if(NOT EXISTS ${UNTAR_DEST}/${CRTM_COEFFS_BRANCH}/)
		message(STATUS "Untarring the downloaded file (~2 minutes) to ${UNTAR_DEST}")
    execute_process(COMMAND tar -xzf ${DOWNLOAD_DEST} -C ${UNTAR_DEST}
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

list( APPEND OMPoverChannels_Sensor_Ids
        atms_n21
)

list (APPEND common_tests
     Simple
     AOD
     Zeeman
     ChannelSubset
     ClearSky
     Aircraft
     ScatteringSwitch
     SOI
     SSU
     VerticalCoordinates
)

list (APPEND omp_tests
     OMPoverChannels
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



add_executable(test_check_crtm mains/application/check_crtm.F90)
target_link_libraries(test_check_crtm PRIVATE crtm)
add_test(NAME test_check_crtm 
         COMMAND test_check_crtm)
set_tests_properties(test_check_crtm PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(test_check_crtm_random  mains/application/check_crtm_random_profiles.F90)
target_link_libraries(test_check_crtm_random PRIVATE crtm)
add_test(NAME test_check_crtm_random
         COMMAND test_check_crtm_random)
set_tests_properties(test_check_crtm_random PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")


add_executable(Unit_TL_TEST mains/unit/Unit_Test/test_TL.f90)
target_link_libraries(Unit_TL_TEST PRIVATE crtm)
add_test(NAME test_Unit_TL_TEST
         COMMAND $<TARGET_FILE:Unit_TL_TEST>)
set_tests_properties(test_Unit_TL_TEST PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(Unit_Aerosol_Bypass mains/unit/Unit_Test/test_Aerosol_Bypass.f90)
target_link_libraries(Unit_Aerosol_Bypass PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass
         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass>)
set_tests_properties(test_Unit_Aerosol_Bypass PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(Unit_Aerosol_Bypass_TL mains/unit/Unit_Test/test_Aerosol_Bypass_TL.f90)
target_link_libraries(Unit_Aerosol_Bypass_TL PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_TL
	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_TL>)
set_tests_properties(test_Unit_Aerosol_Bypass_TL PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(Unit_Aerosol_Bypass_adjoint mains/unit/Unit_Test/test_Aerosol_Bypass_adjoint.f90)
target_link_libraries(Unit_Aerosol_Bypass_adjoint PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_adjoint
	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_adjoint>)
set_tests_properties(test_Unit_Aerosol_Bypass_adjoint  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(Unit_Aerosol_Bypass_k_matrix mains/unit/Unit_Test/test_Aerosol_Bypass_k_matrix.f90)
target_link_libraries(Unit_Aerosol_Bypass_k_matrix PRIVATE crtm)
add_test(NAME test_Unit_Aerosol_Bypass_k_matrix
         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_k_matrix>)
set_tests_properties(test_Unit_Aerosol_Bypass_k_matrix  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

add_executable(SpcCoeff_Inspect mains/application/SpcCoeff_Inspect.f90)
target_link_libraries(SpcCoeff_Inspect PRIVATE crtm)
add_test(NAME run_SpcCoeff_Inspect
         COMMAND $<TARGET_FILE:SpcCoeff_Inspect> "testinput/atms_n21.SpcCoeff.bin")
set_tests_properties(run_SpcCoeff_Inspect PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")


#---------------------------------------------------------------------------------
#unit tests

#first upper level Unit_Test

#(PS) SpcCoeff Version 3 I/O test
add_executable(Spc_IO mains/unit/input_output/test_SpcCoeff/test_spc_io.f90)
target_link_libraries(Spc_IO PRIVATE crtm)
add_test(NAME test_Spc_IO
         COMMAND $<TARGET_FILE:Spc_IO>)
set_tests_properties(test_Spc_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(PS) SpcCoeff netCDF I/O test
add_executable(Spc_IO_NC mains/unit/input_output/test_SpcCoeff_NC/test_spc_io_nc.f90)
target_link_libraries(Spc_IO_NC PRIVATE crtm)
add_test(NAME test_Spc_IO_NC
         COMMAND $<TARGET_FILE:Spc_IO_NC>)
set_tests_properties(test_Spc_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(PS) TauCoeff netCDF I/O test
add_executable(TauCoeff_IO_NC mains/unit/input_output/test_TauCoeff_NC/test_taucoeff_io_nc.f90)
target_link_libraries(TauCoeff_IO_NC PRIVATE crtm)
add_test(NAME test_TauCoeff_IO_NC
         COMMAND $<TARGET_FILE:TauCoeff_IO_NC>)
set_tests_properties(test_TauCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(CD) AerosolCoeff I/O test
add_executable(AerosolCoeff_IO mains/unit/input_output/test_AerosolCoeff/test_aerosol_coeff_io.f90)
target_link_libraries(AerosolCoeff_IO PRIVATE crtm)
add_test(NAME test_AerosolCoeff_IO
         COMMAND $<TARGET_FILE:AerosolCoeff_IO>)
set_tests_properties(test_AerosolCoeff_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(CD) AerosolCoeff netCDF I/O test
add_executable(AerosolCoeff_IO_NC mains/unit/input_output/test_AerosolCoeff_NC/test_aerosol_coeff_io_nc.f90)
target_link_libraries(AerosolCoeff_IO_NC PRIVATE crtm)
add_test(NAME test_AerosolCoeff_IO_NC
         COMMAND $<TARGET_FILE:AerosolCoeff_IO_NC>)
set_tests_properties(test_AerosolCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(CD) CloudCoeff Binary I/O test
add_executable(CloudCoeff_IO mains/unit/input_output/test_CloudCoeff/test_cloud_coeff_io.f90)
target_link_libraries(CloudCoeff_IO PRIVATE crtm)
add_test(NAME test_CloudCoeff_IO
         COMMAND $<TARGET_FILE:CloudCoeff_IO>)
set_tests_properties(test_CloudCoeff_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(CD) CloudCoeff netCDF I/O test
add_executable(CloudCoeff_IO_NC mains/unit/input_output/test_CloudCoeff_NC/test_cloud_coeff_io_nc.f90)
target_link_libraries(CloudCoeff_IO_NC PRIVATE crtm)
add_test(NAME test_CloudCoeff_IO_NC
         COMMAND $<TARGET_FILE:CloudCoeff_IO_NC>)
set_tests_properties(test_CloudCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#(PS) Hypsometric Equation test
add_executable(hypsometric_eq mains/unit/Unit_Test/test_Hypsometric.f90)
target_link_libraries(hypsometric_eq PRIVATE crtm)
add_test(NAME test_hypsometric_eq
         COMMAND $<TARGET_FILE:hypsometric_eq>)
set_tests_properties(test_hypsometric_eq PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")

#=================================================================================
#forward and k_matrix regression tests
foreach(regtype IN LISTS regression_types)
  string(COMPARE EQUAL ${regtype} "k_matrix" isregtype)
  foreach(testtype IN LISTS common_tests)
    string(COMPARE EQUAL ${testtype} "Aircraft" istesttype)
    if (isregtype AND istesttype)
       continue() #skip Aircraft type for k_Matrix
    endif()

    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.f90)
    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)

    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
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
        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
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
        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
    endforeach()
  endforeach()
endforeach()


#---------------------------------------------------------------------------------
# Add test_Active_Sensor executable
add_executable(test_Active_Sensor mains/unit/Unit_Test/test_active_sensor.f90)
target_link_libraries(test_Active_Sensor crtm)

# Add test for test_Active_Sensor
add_test(NAME test_Active_Sensor COMMAND test_Active_Sensor)
set_tests_properties(test_Active_Sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})

# Add test_AD_Active_Sensor executable
add_executable(test_AD_Active_Sensor mains/unit/Unit_Test/test_AD_Active_Sensor.f90)
target_link_libraries(test_AD_Active_Sensor crtm)

# Add test for test_AD_Active_Sensor
add_test(NAME test_AD_Active_Sensor COMMAND test_AD_Active_Sensor)
set_tests_properties(test_AD_Active_Sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})

# Add test_TL_convergence_active_sensor executable
add_executable(test_TL_convergence_active_sensor mains/unit/Unit_Test/test_TL_convergence_active_sensor.f90)
target_link_libraries(test_TL_convergence_active_sensor crtm)

# Add test for test_TL_convergence_active_sensor
add_test(NAME test_TL_convergence_active_sensor COMMAND test_TL_convergence_active_sensor)
set_tests_properties(test_TL_convergence_active_sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})


#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
AerosolCoeff/netCDF/AerosolCoeff.nc4
CloudCoeff/Little_Endian/CloudCoeff.bin
CloudCoeff/netCDF/CloudCoeff.nc4
CloudCoeff/netCDF/CloudCoeff_DDA_ARTS.nc4
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2022.nc4
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/IR_Ice/SEcategory/Little_Endian/NPOESS.IRice.EmisCoeff.bin
EmisCoeff/IR_Ice/SEcategory/netCDF/NPOESS.IRice.EmisCoeff.nc4
EmisCoeff/IR_Land/SEcategory/Little_Endian/NPOESS.IRland.EmisCoeff.bin
EmisCoeff/IR_Land/SEcategory/netCDF/NPOESS.IRland.EmisCoeff.nc4
EmisCoeff/IR_Snow/SEcategory/Little_Endian/NPOESS.IRsnow.EmisCoeff.bin
EmisCoeff/IR_Snow/SEcategory/netCDF/NPOESS.IRsnow.EmisCoeff.nc4
EmisCoeff/IR_Snow/Nalli/Little_Endian/Nalli.IRsnow.EmisCoeff.bin
EmisCoeff/IR_Snow/Nalli/netCDF/Nalli.IRsnow.EmisCoeff.nc4
EmisCoeff/IR_Snow/Nalli/netCDF/Nalli2.IRsnow.EmisCoeff.nc4
EmisCoeff/VIS_Ice/SEcategory/Little_Endian/NPOESS.VISice.EmisCoeff.bin
EmisCoeff/VIS_Ice/SEcategory/netCDF/NPOESS.VISice.EmisCoeff.nc4
EmisCoeff/VIS_Land/SEcategory/Little_Endian/NPOESS.VISland.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/netCDF/NPOESS.VISland.EmisCoeff.nc4
EmisCoeff/VIS_Snow/SEcategory/Little_Endian/NPOESS.VISsnow.EmisCoeff.bin
EmisCoeff/VIS_Snow/SEcategory/netCDF/NPOESS.VISsnow.EmisCoeff.nc4
EmisCoeff/VIS_Water/SEcategory/Little_Endian/NPOESS.VISwater.EmisCoeff.bin
EmisCoeff/VIS_Water/SEcategory/netCDF/NPOESS.VISwater.EmisCoeff.nc4
EmisCoeff/IR_Water/Little_Endian/Nalli.IRwater.EmisCoeff.bin
EmisCoeff/IR_Water/Little_Endian/Nalli2.IRwater.EmisCoeff.bin
EmisCoeff/IR_Water/netCDF/Nalli.IRwater.EmisCoeff.nc4
EmisCoeff/IR_Water/netCDF/Nalli2.IRwater.EmisCoeff.nc4
EmisCoeff/IR_Land/SEcategory/Little_Endian/USGS.IRland.EmisCoeff.bin
EmisCoeff/VIS_Land/SEcategory/Little_Endian/USGS.VISland.EmisCoeff.bin
SpcCoeff/Little_Endian/hirs4_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/hirs4_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_n19.TauCoeff.bin
SpcCoeff/Little_Endian/amsua_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/amsua_metop-a.TauCoeff.bin
SpcCoeff/netCDF/amsua_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_aqua.TauCoeff.nc
SpcCoeff/Little_Endian/gmi_gpm.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/gmi_gpm.TauCoeff.bin
SpcCoeff/Little_Endian/seviri_m08.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/seviri_m08.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n21.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_npp.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-a.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-a.TauCoeff.bin
SpcCoeff/Little_Endian/iasi_metop-b.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/iasi_metop-b.TauCoeff.bin
SpcCoeff/Little_Endian/mhs_n19.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/mhs_n19.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD1_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD1_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD2_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD2_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD3_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD3_g15.TauCoeff.bin
SpcCoeff/Little_Endian/sndrD4_g15.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/sndrD4_g15.TauCoeff.bin
SpcCoeff/Little_Endian/airs_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/airs_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/modis_aqua.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/modis_aqua.TauCoeff.bin
SpcCoeff/Little_Endian/cris-fsr_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris-fsr_n21.TauCoeff.bin
SpcCoeff/Little_Endian/atms_n21.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_n21.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_j2.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/v.viirs-m_j2.TauCoeff.bin
SpcCoeff/Little_Endian/v.abi_g18.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.abi_g18.TauCoeff.bin
SpcCoeff/Little_Endian/abi_g18.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/abi_g18.TauCoeff.bin
SpcCoeff/Little_Endian/cris399_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/cris399_npp.TauCoeff.bin
SpcCoeff/Little_Endian/crisB1_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/crisB1_npp.TauCoeff.bin
SpcCoeff/Little_Endian/atms_npp.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/atms_npp.TauCoeff.bin
SpcCoeff/Little_Endian/v.viirs-m_npp.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.viirs-m_npp.TauCoeff.bin
SpcCoeff/Little_Endian/v.abi_gr.SpcCoeff.bin
TauCoeff/ODAS/Little_Endian/v.abi_gr.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/zssmis_f20.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/zssmis_f19.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/zssmis_f18.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/zssmis_f17.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/zssmis_f16.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f20.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f19.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f18.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f17.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssmis_f16.TauCoeff.bin
SpcCoeff/Little_Endian/ssmis_f20.SpcCoeff.bin
SpcCoeff/Little_Endian/ssmis_f18.SpcCoeff.bin
SpcCoeff/Little_Endian/ssmis_f19.SpcCoeff.bin
SpcCoeff/Little_Endian/ssmis_f16.SpcCoeff.bin
SpcCoeff/Little_Endian/ssmis_f17.SpcCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n06.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n07.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n08.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n09.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n11.TauCoeff.bin
TauCoeff/ODPS/Little_Endian/ssu_n14.TauCoeff.bin
SpcCoeff/Little_Endian/ssu_n06.SpcCoeff.bin
SpcCoeff/Little_Endian/ssu_n07.SpcCoeff.bin
SpcCoeff/Little_Endian/ssu_n08.SpcCoeff.bin
SpcCoeff/Little_Endian/ssu_n09.SpcCoeff.bin
SpcCoeff/Little_Endian/ssu_n11.SpcCoeff.bin
SpcCoeff/Little_Endian/ssu_n14.SpcCoeff.bin
)

# Symlink all CRTM files
# Version 3
CREATE_SYMLINK_FILENAME( ${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput
                         ${crtm_test_input} )
