# (C) Copyright 2017-2018 UCAR.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

cmake_minimum_required (VERSION 3.12)
project("CRTM_Tests" VERSION 3.1.1 LANGUAGES Fortran C)

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
  set( CRTM_COEFFS_BRANCH "fix_REL-3.1.1.3" ) 
  
  set(CRTM_COEFFS_PATH ${CMAKE_BINARY_DIR}/test_data/${PROJECT_VERSION})
  file(MAKE_DIRECTORY ${CRTM_COEFFS_PATH})
  
  set(DOWNLOAD_BASE_URL "https://bin.ssec.wisc.edu/pub/s4/CRTM/")
  set(test_files_dirname ${CRTM_COEFFS_BRANCH}.tgz) 
  set(checksum "400d1bd73ef527ec336ef3fff78e74fd") # MD5SUM of fix_REL-3.1.1.3.tgz
  
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

list (APPEND check_crtm_Sensor_Ids
  iasiB3_metop-a
  amsua_aqua
  viirs-iD09_j1
  viirs-iD17_j1
  viirs-mD13_j1
  v.viirs-iD29_j1
  iasi_g13
  modisD05S_terra
  v.viirs-mD13_npp
  sndrD1_g13
  v.sndr_g11
  micromas_cs03
  metimage_metop-sg-a1
  imgr_g13
  modisD06_aqua
  mersiD02_fy3a
  v.ahi_himawari8
  v.imgrD1_g12
  ssmi_f11
  tms_tropics-04
  viirs-mD01_npp
  hirs3-UWS_n17
  cris-fsr_npp
  v.viirs-mD11_j1
  v.modisD02_aqua
  v.viirs-iD06_npp
  v.imgrD8_g13
  v.imgrD2_g14
  v.viirs-iD15_j1
  iasi_metop-a
  mersiD17_fy3a
  amsua_n16
  hirs4_metop-a
  v.sndrD4_g12
  hirs4-UWS_n18
  mws_metop-sg-a1
  v.viirs-iD27_npp
  amsua_n19
  mersiD29_fy3a
  v.viirs-iD18_npp
  viirs-iD14_npp
  aatsr_envisat
  sndrD2_g15
  modisD07S_aqua
  mersiD03_fy3a
  v.sndrD3_g13
  mwri_fy3b
  modisD07_aqua
  modisS_terra
  mviriBKUP_m06
  abi-81K_g17
  viirs-mD11_npp
  avhrr2_n14
  v.modisD02_terra
  v.imgr_g14
  v.viirs-iD30_j1
  v.sndrD1_g08
  v.sndrD2_g11
  viirs-mD14_j1
  iasiB1_metop-b
  modisD03_terra
  v.seviri_m08
  v.imgrD5_g15
  avhrr3_n16
  v.viirs-mD03_npp
  viirs-iD10_j1
  airsM1b_aqua
  avhrr3_n19
  mersiD28_fy3a
  v.modisD09_terra
  modisD09S_terra
  sndr_g09
  modisD09Shift1_terra
  hirs2-UWS_n10
  msu_n14
  v.imgrD6_g13
  v.viirs-iD08_npp
  viirs-iD04_npp
  modisD03Shift2_aqua
  modisD08Shift2_aqua
  avhrr2_tirosn
  v.mi-lD01_coms
  v.mi-mD08_coms
  viirs-iD25_npp
  v.viirs-iD29_npp
  abi_g19
  hirs4-UWS_metop-b
  abi_g16
  v.viirs-iD12_j1
  atms_n21
  v.viirs-mD16_j1
  mersiD16_fy3a
  v.viirs-mD08_j1
  modisD09Shift1_aqua
  modisD08_terra
  viirs-iD32_j1
  v.imgrD7_g11
  sndr_insat-3ds
  modisD02Shift1_aqua
  cris-fsr_n20
  v.modisD03_aqua
  v.viirs-iD16_npp
  tms_tomorrow-s05_v4
  mersiD15_fy3a
  modisD09Shift2_aqua
  modisD02Shift2_aqua
  vtprS4_itos
  v.viirs-iD26_npp
  tropics_designed_v1
  iasi_metop-c
  geostorm_proposed.boxcar
  sndrD3_g15
  v.viirs-iD19_npp
  viirs-iD15_npp
  v.imgrD8_g11
  modisD03Shift1_aqua
  v.viirs-mD01_j1
  v.iras_fy3a
  viirs-iD25_j1
  mi-lD04_coms
  modisD08Shift1_aqua
  v.abi_g16
  mi-l_coms
  v.sndrD4_g10
  v.mi-lD02_coms
  imgrD1S2_g13
  modisD01Shift1_terra
  v.viirs-iD05_j1
  amsre_aqua
  v.imgrD3_g14
  hirs2_n11
  v.avhrr3_metop-a
  v.viirs-iD07_npp
  cris-fsrB1_npp
  modisD07Shift2_terra
  ssu_n07
  v.imgr_mt2
  modisD05_terra
  aster_terra
  ssu_n08
  v.viirs-m_j1
  airsM1a_aqua
  iasiB3_metop-c
  vhrr_kalpana1
  modisD10_terra
  micromas_cs01
  sndrD1_g11
  v.sndr_g13
  seviri_m11
  viirs-iD19_j1
  viirs-iD07_j1
  imgr_g11
  modisD01S_aqua
  hirs3-UWS_n15
  modisD02S_terra
  viirs-i_n21
  v.modisD04_terra
  v.viirs-mD12_npp
  ssmi_f13
  v.viirs-iD27_j1
  tms_tropics-06
  cris399_npp
  modisD04_aqua
  viirs-mD03_j1
  hirs2-UWS_n12
  atsr2_ers2
  v.viirs-iD02_j1
  mi-m_coms
  v.mi-lD03_coms
  ivissr_fy2e
  mwhs_fy3b
  viirs-iD24_npp
  v.viirs-iD28_npp
  v.sndrD1_g15
  viirs-i_j2
  v.viirs-iD17_npp
  v.imgrD6_g11
  viirs-iD22_j1
  mviriNOM_m03
  v.viirs-mD06_j1
  v.modisD01_aqua
  modisD10_aqua
  imgrD1_g14
  mersiD14_fy3a
  gxs_geoxo_lw
  tmi_trmm
  avhrr2_n06
  tms_tomorrow-s02_v4
  v.viirs-iD09_npp
  v.imgrD7_g13
  viirs-iD05_npp
  modisD03Shift2_terra
  avhrr2_n09
  atms-ng_v1
  gxs_geoxo_mw
  fci_mtg-i1
  cris399_n20
  mviriBKUP_m04
  ssmis_f18.v1-boxcar
  modisD05Shift1_terra
  viirs-mD04_j1
  v.sndrD3_g11
  v.viirs-iD20_j1
  modisD05_aqua
  v.viirs-mD02_npp
  mersiD01_fy3a
  v.imgrD4_g15
  avhrr3JM_n17
  u.omps-npAllFOV_j2
  v.sndrD2_g13
  avhrr3JM_n18
  sndr_g14
  msu_n06
  sndrD4_g14
  tempest-D_cubesat
  msu_n09
  viirs-mD10_npp
  cris-fsrB1_n20
  ssmis_f20
  mersiD39_fy3a
  viirs-iD06_npp
  v.imgrD4_g13
  modisD07S_terra
  avhrr2_n10
  imgrD2_g14
  sndrD4_g12
  v.viirs-iD14_npp
  viirs-iD18_npp
  modisD04Shift2_terra
  modisS_aqua
  viirs-iD12_j1
  avhrr3_metop-c
  v.imgrD5_g11
  v.viirs-iD32_j1
  viirs-mD08_j1
  viirs-mD16_j1
  mersiD07_fy3a
  modisD03_aqua
  sndr_g12
  v.imgrD1_mt2
  v.sndrD2_g15
  iasi316_metop-b
  viirs-iD27_npp
  mersiD12_fy3a
  mviriNOM_m05
  msu_n10
  viirs-mD13_npp
  cris-fsrB2_n20
  vtprS3_itos
  v.modisD07_aqua
  modisD02Shift1_terra
  ivissr_fy2c
  hirs2-UWS_n14
  v.sndrD1_g13
  ssmt2_f14
  v.viirs-mD14_j1
  v.viirs-mD01_npp
  viirs-iD30_j1
  iasiB2_metop-b
  v.mi-lD05_coms
  mi-lD03_coms
  v.imgrD7_g15
  iasi300_metop-a
  vas_g05
  amsua_metop-c
  imgr_mt1r
  v.sndrD3_g08
  v.viirs-iD10_j1
  hirs2_n07
  v.sndr_g15
  viirs-mD11_j1
  v.ivissr_fy2d
  hirs2_n08
  v.viirs-iD04_npp
  viirs-iD08_npp
  ssmis_f19
  modisD02_aqua
  ssmis_f16
  mersiD06_fy3a
  v.sndrD4_g09
  modisD01Shift1_aqua
  viirs-m_npp
  viirs-iD15_j1
  modisD09_terra
  v.sndr_insat-3ds
  viirs-iD16_npp
  mersiD38_fy3a
  modisD06Shift1_terra
  mhs_metop-a
  ssmi_f15
  v.modisD08_terra
  viirs-iD29_npp
  v.viirs-iD25_npp
  modisD10Shift2_terra
  saphir_meghat
  modisD02_terra
  v.mi-lD04_coms
  iasi616_metop-c
  sndrD3_g13
  mi-lD02_coms
  v.viirs-iD17_j1
  ssmis_f18.v2-ff_passband
  v.viirs-iD09_j1
  v.viirs-mD11_npp
  hirs3_n15
  v.viirs-i_npp
  mhs_n18
  viirs-iD29_j1
  v.modis_terra
  v.viirs-mD13_j1
  v.avhrr2_n14
  imgr_g08
  viirs-mD03_npp
  vtprS2_itos
  v.modisD06_aqua
  atms_npp.nasa.T20_Vlow_LOpri
  sndrD2_g11
  modisD06S_aqua
  v.imgrD3_g12
  seviri_m08
  mersiD13_fy3a
  ahi_himawari8
  v.modisD03_terra
  cris-fsrB2_npp
  v.avhrr3_n16
  airs324_aqua
  ssu_n11
  v.avhrr3_n19
  masc_cubesat
  v.imgrD6_g15
  mwts_fy3b
  mviriNOM_m07
  msu_n12
  viirs-m_j1
  v.modisD10_terra
  viirs-iD20_j1
  v.viirs-mD04_j1
  iasi-ng_metop-sg-a1
  v.modisD05_aqua
  v.sndrD1_g11
  vtprS1_itos
  v.sndrD2_g08
  mersiD10_fy3a
  v.modisD05_terra
  modisD08Shift2_terra
  iasi300_metop-c
  vas_g07
  cris-fsrB3_n20
  viirs-mD12_npp
  modisD04_terra
  v.mi-lD07_coms
  mi-lD01_coms
  amsua_metop-a
  v.mi-l_coms
  avhrr2_n12
  iras_fy3a
  v.imgrD4_g11
  ssmis_f19.v2-ff_passband
  viirs-iD02_j1
  v.imgr_g12
  gmi_gpm
  v.viirs-iD15_npp
  viirs-iD19_npp
  gmi_gpm.TNOM
  viirs-iD26_npp
  v.sndrD3_g15
  v.metimage_metop-sg-a1
  viirs-mD06_j1
  v.viirs-iD22_j1
  modisD01_aqua
  v.modisD10_aqua
  avhrr3_metop-a
  viirs-iD07_npp
  v.imgrD5_g13
  modisD01Shift2_aqua
  virr_fy3a
  hirs4_n19
  sndrD4_g10
  mersiD05_fy3a
  airsM9_aqua
  sndr_g10
  hirs2-UWS_n06
  hirs2-UWS_n09
  hirs3_n17
  v.imgrD2_g12
  ssmi_f08
  giirsB1_fsr_fy4a
  v.imgrD8_g15
  ami_gk2
  sndrD3_g11
  iasi616_metop-a
  viirs-mD02_npp
  v.mi-lD06_coms
  cris-fsrB3_npp
  u.omps-tcAllFOV_j2
  v.viirs-mD10_npp
  v.viirs-iD07_j1
  v.viirs-iD19_j1
  v.sndr_g08
  sndrD2_g13
  windsat_coriolis
  v.sndrD4_g14
  mersiD11_fy3a
  v.modisD04_aqua
  v.viirs-mD03_j1
  viirs-iD27_j1
  mersiD04_fy3a
  airsM8_aqua
  hsb_aqua
  viirs-iD17_npp
  micromas_cs05
  sndrD1_g15
  v.viirs-i_j2
  msu_tirosn
  v.ivissr_fy2f
  dpr_gpm
  amsub_n16
  viirs-iD28_npp
  v.viirs-iD24_npp
  cowvr_ors6
  modisShift2_terra
  v.viirs-iD05_npp
  viirs-iD09_npp
  mhs_metop-c
  tms_tropics-02
  v.imgrD1_g14
  v.viirs-iD25_j1
  viirs-mD01_j1
  ssmis_f19.v1-boxcar
  v.mi-m_coms
  viirs-iD05_j1
  imgr_g15
  viirs-iD01_npp
  v.imgrD3_g13
  ahi_himawari9
  seviri_m09
  airsM4d_aqua
  modisD10Shift1_terra
  sndrD2_g10
  viirs-iD21_j1
  v.viirs-mD05_j1
  v.avhrr3_n18
  modisD06Shift2_terra
  v.viirs-iD01_j1
  modisD05Shift2_aqua
  mersiD08_fy3a
  crisB1_n20
  v.avhrr3_n17
  cris-fsr431_n20
  v.viirs-iD32_npp
  v.imgrD2_g11
  mersiD36_fy3a
  tms_tomorrow-s01_v4
  sndrD3_g12
  modisD01S_terra
  iasi616_metop-b
  v.viirs-iD13_npp
  atms_npp.nasa.T-10_Vnom_LOpri
  imgr_insat-3ds
  imgr_g09
  viirs-iD20_npp
  mhs_n19
  modisD04Shift1_aqua
  modisD10S_aqua
  tms_tropics-01
  ssmi_f14
  atms_n21-srf
  viirs-mD14_npp
  viirs-iD03_j1
  mersiD40_fy3a
  miras_smos
  v.mi-mD03_coms
  mersiD23_fy3a
  modisD05S_aqua
  ssmis_f20.v2-ff_passband
  mwi_metop-sg-a1
  viirs-mD07_j1
  v.viirs-iD23_j1
  hirs2_n09
  v.sndr_g14
  hirs2_n06
  v.viirs-i_j1
  v.modisD08_aqua
  v.viirs-mD06_npp
  v.ivissr_fy2e
  v.sndrD4_g08
  amsub_n15
  ssmis_f17
  cris374_n20
  modisD01_terra
  ssmis_f18
  v.imgrD7_g14
  vas_g04
  iasiB2_metop-c
  v.viirs-iD03_npp
  mersi_fy3a
  mersiD37_fy3a
  hirs2_tirosn
  tms_tomorrow-s06_v4
  v.sndrD3_g09
  modisD04Shift1_terra
  viirs-iD30_npp
  v.viirs-i_n21
  imgrD1_g13
  amsua_metop-b
  cris-fsr431_npp
  mviriNOM_m04
  v.viirs-iD06_j1
  mwts_fy3a
  v.viirs-iD18_j1
  msu_n11
  crisB1_npp
  viirs-m_j2
  viirs-iD11_npp
  v.aster_terra
  v.viirs-iD22_npp
  ssmt2_f15
  v.sndrD1_g12
  v.viirs-mD02_j1
  airsM5_aqua
  mersiD09_fy3a
  viirs-iD26_j1
  v.modisD09_aqua
  imgr_mt2
  v.viirs-mD16_npp
  avhrr3_metop-b
  cris374_npp
  modisShift1_aqua
  sndrD4_g13
  viirs-m_n21
  v.sndrD2_g14
  iasi316_metop-c
  sndr_g13
  ssmt1_f13
  atms_j2
  airs_g13
  avhrr2_n11
  iras_fy3b
  v.imgrD4_g12
  mi-mD04_coms
  v.mi-mD02_coms
  v.viirs-mD08_npp
  v.imgr_g11
  viirs-mD04_npp
  v.viirs-iD24_j1
  modisD02Shift2_terra
  viirs-iD04_j1
  imgrD2_g15
  mviriBKUP_m03
  airs_aqua
  mersiD22_fy3a
  v.viirs-mD07_npp
  v.imgrD1_g15
  mersiD21_fy3a
  mhs_metop-b
  tms_tropics-03
  v.mi-mD01_coms
  v.mi-lD08_coms
  amsr2_gcom-w1
  imgr_g14
  atms_npp.nasa.T50_Vnom_LOpri
  sndrD1_g14
  micromas_cs04
  slstr_sentinel3a
  viirs-iD13_j1
  viirs-mD15_npp
  atms_n20-srf
  modisShift2_aqua
  viirs-mD09_j1
  giirsB2_fsr_fy4a
  amsub_n17
  modis_aqua
  ssmis_f16.v1-boxcar
  v.viirs-iD12_npp
  v.sndr_g09
  sndrD2_g12
  v.imgrD3_g11
  airsM6_aqua
  modisD06S_terra
  hirs2_n14
  hamsr_grip
  v.avhrr3_n15
  ssu_pseudo
  viirs-iD21_npp
  v.sndrD4_g15
  v.viirs-mD15_j1
  modisShift1_terra
  sndrD3_g10
  gmi_gpm.TLO
  viirs-iD31_j1
  v.imgrD8_g14
  hirs3_n16
  v.imgrD2_g13
  v.viirs-iD11_j1
  mersiD34_fy3a
  viirs-mD10_j1
  hirs4_n18
  sndrD4_g11
  v.viirs-mD09_npp
  viirs-mD05_npp
  v.imgrD5_g12
  modisD07_terra
  airs281_aqua
  sndr_g11
  hirs2-UWS_n07
  viirs-iD14_j1
  modisD08Shift1_terra
  iasi316_metop-a
  v.modisD06_terra
  v.imgr_g13
  mersiD20_fy3a
  modisD03S_aqua
  v.sndrD3_g14
  viirs-iD10_npp
  iasiB2_metop-a
  v.viirs-iD16_j1
  v.viirs-iD08_j1
  modisD04Shift2_aqua
  iasi300_metop-b
  vas_g06
  atms_npp.nasa.T20_Vhi_LOpri
  mersiD35_fy3a
  viirs-iD28_j1
  v.viirs-mD12_j1
  v.viirs-iD23_npp
  airsM7_aqua
  modisD05Shift1_aqua
  modisD08S_aqua
  mviriNOM_m06
  v.viirs-iD02_npp
  v.imgrD6_g14
  ssu_tirosn
  v.seviri_m10
  viirs-iD31_npp
  v.sndrD2_g09
  v.sndrD1_g10
  atms_n20
  airsM2b_aqua
  tms_tomorrow-s04_v4
  modisD09_aqua
  abi_g17
  abi_g18
  ssmis_f16.v3-eff_passband
  cris-fsr_n21
  modisD10Shift1_aqua
  vissrDetB_gms5
  airsM4a_aqua
  v.viirs-mD14_npp
  hirs2-UWS_n11
  ivissr_fy2f
  v.modisD01_terra
  sndr_g08
  mwhs_fy3a
  viirs-i_j1
  viirs-iD24_j1
  viirs-mD06_npp
  modisD03S_terra
  v.viirs-iD04_j1
  mersiD33_fy3a
  v.imgrD6_g12
  modisD04S_aqua
  v.seviri_m09
  cris-fsr_j2
  viirs-iD32_npp
  v.sndrD2_g10
  v.sndrD1_g09
  iasiB1_metop-c
  mersiD26_fy3a
  avhrr3_n18
  v.mi-mD06_coms
  avhrr3_n17
  v.viirs-iD01_npp
  v.imgrD5_g14
  modisD07Shift1_aqua
  viirs-iD18_j1
  viirs-iD06_j1
  mviriBKUP_m07
  madras_meghat
  v.sndrD3_g12
  mi-l_coms.v2
  v.viirs-iD20_npp
  v.imgr_g15
  viirs-iD13_npp
  modisD09Shift2_terra
  crisB3_npp
  modisD06Shift2_aqua
  v.imgr_insat-3ds
  v.viirs-iD26_j1
  viirs-mD02_j1
  mersiD18_fy3a
  v.viirs-iD03_j1
  atms_npp.nasa.T20_Vnom_LOpri
  amsua_n18
  atms_j2-srf
  v.iras_fy3b
  ssu_n14
  amsua_n17
  v.sndrD4_g13
  hirs4-UWS_n19
  ssmis_f17.v2-ff_passband
  sndrD2_g14
  mersiD32_fy3a
  viirs-mD16_npp
  viirs-iD23_j1
  hirs2_n12
  v.avhrr3_metop-b
  v.viirs-mD07_j1
  radiometer_smap
  modisD10S_terra
  modisD08_aqua
  v.viirs-mD04_npp
  viirs-mD08_npp
  tms_tomorrow-s03_v4
  airsM12_aqua
  v.imgrD2_g15
  v.imgrD8_g12
  atms_npp
  v.viirs-iD30_npp
  v.ahi_himawari9
  imgr_g12
  viirs-mD05_j1
  crisB3_n20
  v.viirs-iD21_j1
  hirs3-UWS_n16
  mersiD19_fy3a
  v.abi_gr
  amsua_metop-a_v2
  viirs-iD01_j1
  tms_tropics-05
  ssmi_f10
  v.imgrD1_g13
  viirs-iD03_npp
  geostorm_proposed.gaussian
  mersiD27_fy3a
  viirs-iD22_npp
  v.viirs-m_j2
  v.viirs-dnb_j1
  micromas_cs02
  v.sndr_g10
  sndrD1_g12
  v.viirs-iD11_npp
  v.mi-mD07_coms
  mi-mD01_coms
  ssmis_f16.v2-ff_passband
  sndr_g15
  modisD03Shift1_terra
  v.viirs-iD21_npp
  ssmt1_f15
  v.modis_aqua
  v.sndrD2_g12
  mwi_wsf-m1
  iasiB1_metop-a
  gmi_gpm.THI
  v.mi-mD04_coms
  viirs-iD08_j1
  mi-mD02_coms
  viirs-iD16_j1
  viirs-mD12_j1
  msu_n08
  mersiD24_fy3a
  v.viirs-iD28_j1
  avhrr3_n15
  sndrD4_g15
  crisB2_npp
  msu_n07
  viirs-iD12_npp
  viirs-i_npp
  v.sndrD3_g10
  mwri_fy3a
  imgrD2_g13
  mviriBKUP_m05
  v.imgrD4_g14
  modis_terra
  cris_n20
  eon_mw.v1
  airsM4c_aqua
  v.viirs-mD10_j1
  ssmis_f20.v1-boxcar
  hirs4-UWS_metop-a
  imgrD1_g15
  abi_gr
  airsM11_aqua
  modisD04S_terra
  avhrr2_n08
  airsM3_aqua
  v.imgrD7_g12
  v.viirs-iD14_j1
  viirs-mD07_npp
  avhrr2_n07
  mersiD31_fy3a
  v.sndrD1_g14
  ivissr_fy2d
  ssmis_f17.v1-boxcar
  v.viirs-mD15_npp
  v.viirs-m_npp
  modisD05Shift2_terra
  modisD06Shift1_aqua
  modisD08S_terra
  viirs-iD23_npp
  imgr_g10
  modisD01Shift2_terra
  modisD09S_aqua
  v.imgrD1_g11
  tms_tropics-07
  v.viirs-iD10_npp
  v.viirs-iD31_j1
  v.mi-mD05_coms
  iasiB3_metop-b
  mi-mD03_coms
  crisB2_n20
  viirs-mD15_j1
  ssu_n09
  modisD07Shift2_aqua
  v.viirs-iD31_npp
  ssu_n06
  seviri_m10
  viirs-iD02_npp
  viirs-iD11_j1
  mersiD25_fy3a
  v.sndr_g12
  sndrD1_g10
  v.ivissr_fy2c
  cpr_cloudsat
  micromas_cs00
  hirs4_metop-b
  amsua_n15
  v.sndrD4_g11
  v.abi_g17
  v.modisD07_terra
  v.abi_g18
  mersiD30_fy3a
  modisD02S_aqua
  v.imgrD3_g15
  hirs2_n10
  atsr1_ers1
  v.viirs-mD05_npp
  viirs-mD09_npp
  airsM10_aqua
  vissrDetA_gms5
  modisD10Shift2_aqua
  v.viirs-iD13_j1
  airsM4b_aqua
  v.viirs-mD09_j1
  airsM2a_aqua
  modisD06_terra
  sndrD3_g14
  modisD07Shift1_terra
  iasi_metop-b
  cris_npp  
)

#####################################################################

#add_executable(check_tropics mains/application/check_tropics.f90)
#target_link_libraries(check_tropics PRIVATE crtm)
#
add_executable(test_check_crtm mains/application/check_crtm.F90)
target_link_libraries(test_check_crtm PRIVATE crtm)

foreach(sensor_id IN LISTS check_crtm_Sensor_Ids)
  add_test(NAME test_check_crtm_${sensor_id}
    COMMAND test_check_crtm "${sensor_id}")
  set_tests_properties(test_check_crtm_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
endforeach()


#add_executable(test_check_crtm_random  mains/application/check_crtm_random_profiles.F90)
#target_link_libraries(test_check_crtm_random PRIVATE crtm)
#add_test(NAME test_check_crtm_random
#         COMMAND test_check_crtm_random)
#set_tests_properties(test_check_crtm_random PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
#
#add_executable(Unit_TL_TEST mains/unit/Unit_Test/test_TL.f90)
#target_link_libraries(Unit_TL_TEST PRIVATE crtm)
#add_test(NAME test_Unit_TL_TEST
#         COMMAND $<TARGET_FILE:Unit_TL_TEST>)
#set_tests_properties(test_Unit_TL_TEST PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
#add_executable(Unit_Aerosol_Bypass mains/unit/Unit_Test/test_Aerosol_Bypass.f90)
#target_link_libraries(Unit_Aerosol_Bypass PRIVATE crtm)
#add_test(NAME test_Unit_Aerosol_Bypass
#         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass>)
#set_tests_properties(test_Unit_Aerosol_Bypass PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
#add_executable(Unit_Aerosol_Bypass_TL mains/unit/Unit_Test/test_Aerosol_Bypass_TL.f90)
#target_link_libraries(Unit_Aerosol_Bypass_TL PRIVATE crtm)
#add_test(NAME test_Unit_Aerosol_Bypass_TL
#	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_TL>)
#set_tests_properties(test_Unit_Aerosol_Bypass_TL PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
#add_executable(Unit_Aerosol_Bypass_adjoint mains/unit/Unit_Test/test_Aerosol_Bypass_adjoint.f90)
#target_link_libraries(Unit_Aerosol_Bypass_adjoint PRIVATE crtm)
#add_test(NAME test_Unit_Aerosol_Bypass_adjoint
#	 COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_adjoint>)
#set_tests_properties(test_Unit_Aerosol_Bypass_adjoint  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
#add_executable(Unit_Aerosol_Bypass_k_matrix mains/unit/Unit_Test/test_Aerosol_Bypass_k_matrix.f90)
#target_link_libraries(Unit_Aerosol_Bypass_k_matrix PRIVATE crtm)
#add_test(NAME test_Unit_Aerosol_Bypass_k_matrix
#         COMMAND $<TARGET_FILE:Unit_Aerosol_Bypass_k_matrix>)
#set_tests_properties(test_Unit_Aerosol_Bypass_k_matrix  PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##SpcCoeff utilities
#list (APPEND SCoeff_Utils
#	SpcCoeff_Edit
#	SpcCoeff_Inspect
#	SpcCoeff_BIN2NC
#	SpcCoeff_NC2BIN
#)
#foreach(testtype IN LISTS SCoeff_Utils)
#	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/SpcCoeff/${testtype}/${testtype}.f90)
#	target_link_libraries(${testtype} PRIVATE crtm)
#endforeach()
#
##TauCoeff ODAS utilities
#list (APPEND ODAS_Utils
#	ODAS_BIN2NC
#	ODAS_NC2BIN
#)
#foreach(testtype IN LISTS ODAS_Utils)
#	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/ODAS/${testtype}/${testtype}.f90)
#	target_link_libraries(${testtype} PRIVATE crtm)
#endforeach()
#
##TauCoeff ODPS utilities
#list (APPEND ODPS_Utils
#	ODPS_BIN2NC
#	ODPS_NC2BIN
#)
#foreach(testtype IN LISTS ODPS_Utils)
#	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/TauCoeff/ODPS/${testtype}/${testtype}.f90)
#	target_link_libraries(${testtype} PRIVATE crtm)
#endforeach()
#
##CloudCoeff utilities 
#list (APPEND CCoeff_Utils
#  CloudCoeff_NC2BIN
#	CloudCoeff_BIN2NC
#	CloudCoeff_Inspect
#)
#foreach(testtype IN LISTS CCoeff_Utils)
#	add_executable(${testtype} ${CRTM_SOURCE_DIR}/src/Coefficients/CloudCoeff/${testtype}/${testtype}.f90)
#	target_link_libraries(${testtype} PRIVATE crtm)
#endforeach()
#
#

#---------------------------------------------------------------------------------
#unit tests

#first upper level Unit_Test

#(PS) SpcCoeff Version 3 I/O test
#add_executable(Spc_IO mains/unit/input_output/test_SpcCoeff/test_spc_io.f90)
#target_link_libraries(Spc_IO PRIVATE crtm)
#add_test(NAME test_Spc_IO
#         COMMAND $<TARGET_FILE:Spc_IO>)
#set_tests_properties(test_Spc_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(PS) SpcCoeff netCDF I/O test
#add_executable(Spc_IO_NC mains/unit/input_output/test_SpcCoeff_NC/test_spc_io_nc.f90)
#target_link_libraries(Spc_IO_NC PRIVATE crtm)
#add_test(NAME test_Spc_IO_NC
#         COMMAND $<TARGET_FILE:Spc_IO_NC>)
#set_tests_properties(test_Spc_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(PS) TauCoeff netCDF I/O test
#add_executable(TauCoeff_IO_NC mains/unit/input_output/test_TauCoeff_NC/test_taucoeff_io_nc.f90)
#target_link_libraries(TauCoeff_IO_NC PRIVATE crtm)
#add_test(NAME test_TauCoeff_IO_NC
#         COMMAND $<TARGET_FILE:TauCoeff_IO_NC>)
#set_tests_properties(test_TauCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(CD) AerosolCoeff I/O test
#add_executable(AerosolCoeff_IO mains/unit/input_output/test_AerosolCoeff/test_aerosol_coeff_io.f90)
#target_link_libraries(AerosolCoeff_IO PRIVATE crtm)
#add_test(NAME test_AerosolCoeff_IO
#         COMMAND $<TARGET_FILE:AerosolCoeff_IO>)
#set_tests_properties(test_AerosolCoeff_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(CD) AerosolCoeff netCDF I/O test
#add_executable(AerosolCoeff_IO_NC mains/unit/input_output/test_AerosolCoeff_NC/test_aerosol_coeff_io_nc.f90)
#target_link_libraries(AerosolCoeff_IO_NC PRIVATE crtm)
#add_test(NAME test_AerosolCoeff_IO_NC
#         COMMAND $<TARGET_FILE:AerosolCoeff_IO_NC>)
#set_tests_properties(test_AerosolCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(CD) CloudCoeff Binary I/O test
#add_executable(CloudCoeff_IO mains/unit/input_output/test_CloudCoeff/test_cloud_coeff_io.f90)
#target_link_libraries(CloudCoeff_IO PRIVATE crtm)
#add_test(NAME test_CloudCoeff_IO
#         COMMAND $<TARGET_FILE:CloudCoeff_IO>)
#set_tests_properties(test_CloudCoeff_IO PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(CD) CloudCoeff netCDF I/O test
#add_executable(CloudCoeff_IO_NC mains/unit/input_output/test_CloudCoeff_NC/test_cloud_coeff_io_nc.f90)
#target_link_libraries(CloudCoeff_IO_NC PRIVATE crtm)
#add_test(NAME test_CloudCoeff_IO_NC
#         COMMAND $<TARGET_FILE:CloudCoeff_IO_NC>)
#set_tests_properties(test_CloudCoeff_IO_NC PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##(PS) Hypsometric Equation test
#add_executable(hypsometric_eq mains/unit/Unit_Test/test_Hypsometric.f90)
#target_link_libraries(hypsometric_eq PRIVATE crtm)
#add_test(NAME test_hypsometric_eq
#         COMMAND $<TARGET_FILE:hypsometric_eq>)
#set_tests_properties(test_hypsometric_eq PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
## test_MWwater_io
#add_executable(MWwater_IO mains/unit/input_output/test_MWwater/test_MWwater_io.f90)
#target_link_libraries(MWwater_IO PRIVATE crtm)
#add_test(NAME test_MWwater_io
#         COMMAND $<TARGET_FILE:MWwater_IO>)
#set_tests_properties(test_MWwater_io PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
## test_MWSurfEM
#add_executable(test_MWSurfEM mains/unit/Unit_Test/test_MWSurfEM.f90)
#target_link_libraries(test_MWSurfEM PRIVATE crtm)
#add_test(NAME test_MWSurfEM
#         COMMAND $<TARGET_FILE:test_MWSurfEM>)
#set_tests_properties(test_MWSurfEM PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#
##=================================================================================
##forward and k_matrix regression tests
#foreach(regtype IN LISTS regression_types)
#  string(COMPARE EQUAL ${regtype} "k_matrix" isregtype)
#  foreach(testtype IN LISTS common_tests)
#    string(COMPARE EQUAL ${testtype} "Aircraft" istesttype)
#    if (isregtype AND istesttype)
#       continue() #skip Aircraft type for k_Matrix
#    endif()
#    string(COMPARE EQUAL ${testtype} "Downwelling_Radiance" istesttype)
#    if (isregtype AND istesttype)
#       continue() #skip Downwelling_Radiance type for k_Matrix
#    endif()
#
#    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.f90)
#    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)
#
#    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
#        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
#                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
#        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#    endforeach()
#  endforeach()
#endforeach()
#
##---------------------------------------------------------------------------------
##TLAD Regression tests
#foreach(regtype IN LISTS TLAD_types)
#  foreach(testtype IN LISTS TLAD_tests)
#    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.f90)
#    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)
#
#    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
#        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
#                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
#        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#    endforeach()
#  endforeach()
#endforeach()
#
##=================================================================================
##OpenMP regression tests
#foreach(regtype IN LISTS regression_types)
#  string(COMPARE EQUAL ${regtype} "k_matrix" isregtype)
#  if (isregtype)
#    continue() #skip all k_matrix tests
#  endif()
#  foreach(testtype IN LISTS omp_tests)
#    add_executable(test_${regtype}_test_${testtype} mains/regression/${regtype}/test_${testtype}/test_${testtype}.F90)
#    target_link_libraries(test_${regtype}_test_${testtype} PRIVATE crtm)
#
#    foreach(sensor_id IN LISTS ${testtype}_Sensor_Ids)
#        add_test(NAME test_${regtype}_${testtype}_${sensor_id}
#                 COMMAND $<TARGET_FILE:test_${regtype}_test_${testtype}> "${sensor_id}")
#        set_tests_properties(test_${regtype}_${testtype}_${sensor_id} PROPERTIES ENVIRONMENT "OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS}")
#    endforeach()
#  endforeach()
#endforeach()
#
#
##---------------------------------------------------------------------------------
## Add test_Active_Sensor executable
#add_executable(test_Active_Sensor mains/unit/Unit_Test/test_active_sensor.f90)
#target_link_libraries(test_Active_Sensor crtm)
#
## Add test for test_Active_Sensor
#add_test(NAME test_Active_Sensor COMMAND test_Active_Sensor)
#set_tests_properties(test_Active_Sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})
#
## Add test_AD_Active_Sensor executable
#add_executable(test_AD_Active_Sensor mains/unit/Unit_Test/test_AD_Active_Sensor.f90)
#target_link_libraries(test_AD_Active_Sensor crtm)
#
## Add test for test_AD_Active_Sensor
#add_test(NAME test_AD_Active_Sensor COMMAND test_AD_Active_Sensor)
#set_tests_properties(test_AD_Active_Sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})
#
## Add test_TL_convergence_active_sensor executable
#add_executable(test_TL_convergence_active_sensor mains/unit/Unit_Test/test_TL_convergence_active_sensor.f90)
#target_link_libraries(test_TL_convergence_active_sensor crtm)
#
## Add test for test_TL_convergence_active_sensor
#add_test(NAME test_TL_convergence_active_sensor COMMAND test_TL_convergence_active_sensor)
#set_tests_properties(test_TL_convergence_active_sensor PROPERTIES ENVIRONMENT OMP_NUM_THREADS=$ENV{OMP_NUM_THREADS})


#####################################################################
# Files for CRTM tests
#####################################################################

list( APPEND crtm_test_input
AerosolCoeff/Little_Endian/AerosolCoeff.bin
AerosolCoeff/netCDF/AerosolCoeff.nc4
CloudCoeff/Little_Endian/CloudCoeff.bin
CloudCoeff/netCDF/CloudCoeff.nc4
CloudCoeff/netCDF/CloudCoeff_DDA_Moradi_2022.nc4
EmisCoeff/MW_Water/Little_Endian/FASTEM6.MWwater.EmisCoeff.bin
EmisCoeff/MW_Water/Little_Endian/FASTEM4.MWwater.EmisCoeff.bin
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
SpcCoeff/netCDF/iasiB3_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB3_metop-a.TauCoeff.nc
SpcCoeff/netCDF/amsua_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD09_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD09_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD17_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD17_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD13_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD13_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD29_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD29_j1.TauCoeff.nc
SpcCoeff/netCDF/iasi_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_g13.TauCoeff.nc
SpcCoeff/netCDF/modisD05S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05S_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD13_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD13_npp.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g13.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g11.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs03.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs03.TauCoeff.nc
SpcCoeff/netCDF/metimage_metop-sg-a1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/metimage_metop-sg-a1.TauCoeff.nc
SpcCoeff/netCDF/imgr_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g13.TauCoeff.nc
SpcCoeff/netCDF/modisD06_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06_aqua.TauCoeff.nc
SpcCoeff/netCDF/mersiD02_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD02_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.ahi_himawari8.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ahi_himawari8.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_g12.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f11.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-04.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-04.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD01_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD01_npp.TauCoeff.nc
SpcCoeff/netCDF/hirs3-UWS_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3-UWS_n17.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD11_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD11_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD02_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD02_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD06_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD06_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD8_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD8_g13.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD2_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD2_g14.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD15_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD15_j1.TauCoeff.nc
SpcCoeff/netCDF/iasi_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_metop-a.TauCoeff.nc
SpcCoeff/netCDF/mersiD17_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD17_fy3a.TauCoeff.nc
SpcCoeff/netCDF/amsua_n16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n16.TauCoeff.nc
SpcCoeff/netCDF/hirs4_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4_metop-a.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g12.TauCoeff.nc
SpcCoeff/netCDF/hirs4-UWS_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4-UWS_n18.TauCoeff.nc
SpcCoeff/netCDF/mws_metop-sg-a1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mws_metop-sg-a1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD27_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD27_npp.TauCoeff.nc
SpcCoeff/netCDF/amsua_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n19.TauCoeff.nc
SpcCoeff/netCDF/mersiD29_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD29_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD18_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD18_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD14_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD14_npp.TauCoeff.nc
SpcCoeff/netCDF/aatsr_envisat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/aatsr_envisat.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g15.TauCoeff.nc
SpcCoeff/netCDF/modisD07S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07S_aqua.TauCoeff.nc
SpcCoeff/netCDF/mersiD03_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD03_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g13.TauCoeff.nc
SpcCoeff/netCDF/mwri_fy3b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwri_fy3b.TauCoeff.nc
SpcCoeff/netCDF/modisD07_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisS_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisS_terra.TauCoeff.nc
SpcCoeff/netCDF/mviriBKUP_m06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriBKUP_m06.TauCoeff.nc
SpcCoeff/netCDF/abi-81K_g17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi-81K_g17.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD11_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD11_npp.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n14.TauCoeff.nc
SpcCoeff/netCDF/v.modisD02_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD02_terra.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_g14.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD30_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD30_j1.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g08.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g11.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD14_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD14_j1.TauCoeff.nc
SpcCoeff/netCDF/iasiB1_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB1_metop-b.TauCoeff.nc
SpcCoeff/netCDF/modisD03_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03_terra.TauCoeff.nc
SpcCoeff/netCDF/v.seviri_m08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.seviri_m08.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD5_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD5_g15.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_n16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_n16.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD03_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD03_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD10_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD10_j1.TauCoeff.nc
SpcCoeff/netCDF/airsM1b_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM1b_aqua.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_n19.TauCoeff.nc
SpcCoeff/netCDF/mersiD28_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD28_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.modisD09_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD09_terra.TauCoeff.nc
SpcCoeff/netCDF/modisD09S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09S_terra.TauCoeff.nc
SpcCoeff/netCDF/sndr_g09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g09.TauCoeff.nc
SpcCoeff/netCDF/modisD09Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n10.TauCoeff.nc
SpcCoeff/netCDF/msu_n14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n14.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD6_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD6_g13.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD08_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD08_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD04_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD04_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD03Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD08Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_tirosn.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_tirosn.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD01_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD01_coms.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD08_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD08_coms.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD25_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD25_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD29_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD29_npp.TauCoeff.nc
SpcCoeff/netCDF/abi_g19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_g19.TauCoeff.nc
SpcCoeff/netCDF/hirs4-UWS_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4-UWS_metop-b.TauCoeff.nc
SpcCoeff/netCDF/abi_g16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_g16.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD12_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD12_j1.TauCoeff.nc
SpcCoeff/netCDF/atms_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_n21.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD16_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD16_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD16_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD16_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD08_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD08_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD09Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD08_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD32_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD32_j1.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD7_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD7_g11.TauCoeff.nc
SpcCoeff/netCDF/sndr_insat-3ds.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_insat-3ds.TauCoeff.nc
SpcCoeff/netCDF/modisD02Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_n20.TauCoeff.nc
SpcCoeff/netCDF/v.modisD03_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD03_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD16_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD16_npp.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s05_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s05_v4.TauCoeff.nc
SpcCoeff/netCDF/mersiD15_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD15_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD09Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD02Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/vtprS4_itos.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vtprS4_itos.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD26_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD26_npp.TauCoeff.nc
SpcCoeff/netCDF/tropics_designed_v1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tropics_designed_v1.TauCoeff.nc
SpcCoeff/netCDF/iasi_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_metop-c.TauCoeff.nc
SpcCoeff/netCDF/geostorm_proposed.boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/geostorm_proposed.boxcar.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g15.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD19_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD19_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD15_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD15_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD8_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD8_g11.TauCoeff.nc
SpcCoeff/netCDF/modisD03Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD01_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD01_j1.TauCoeff.nc
SpcCoeff/netCDF/v.iras_fy3a.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.iras_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD25_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD25_j1.TauCoeff.nc
SpcCoeff/netCDF/mi-lD04_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-lD04_coms.TauCoeff.nc
SpcCoeff/netCDF/modisD08Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.abi_g16.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_g16.TauCoeff.nc
SpcCoeff/netCDF/mi-l_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-l_coms.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g10.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD02_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD02_coms.TauCoeff.nc
SpcCoeff/netCDF/imgrD1S2_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD1S2_g13.TauCoeff.nc
SpcCoeff/netCDF/modisD01Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD05_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD05_j1.TauCoeff.nc
SpcCoeff/netCDF/amsre_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsre_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD3_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD3_g14.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n11.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_metop-a.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_metop-a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD07_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD07_npp.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB1_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB1_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD07Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/ssu_n07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n07.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_mt2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_mt2.TauCoeff.nc
SpcCoeff/netCDF/modisD05_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05_terra.TauCoeff.nc
SpcCoeff/netCDF/aster_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/aster_terra.TauCoeff.nc
SpcCoeff/netCDF/ssu_n08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n08.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-m_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-m_j1.TauCoeff.nc
SpcCoeff/netCDF/airsM1a_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM1a_aqua.TauCoeff.nc
SpcCoeff/netCDF/iasiB3_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB3_metop-c.TauCoeff.nc
SpcCoeff/netCDF/vhrr_kalpana1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vhrr_kalpana1.TauCoeff.nc
SpcCoeff/netCDF/modisD10_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10_terra.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs01.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs01.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g11.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g13.TauCoeff.nc
SpcCoeff/netCDF/seviri_m11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/seviri_m11.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD19_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD19_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD07_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD07_j1.TauCoeff.nc
SpcCoeff/netCDF/imgr_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g11.TauCoeff.nc
SpcCoeff/netCDF/modisD01S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01S_aqua.TauCoeff.nc
SpcCoeff/netCDF/hirs3-UWS_n15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3-UWS_n15.TauCoeff.nc
SpcCoeff/netCDF/modisD02S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02S_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-i_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-i_n21.TauCoeff.nc
SpcCoeff/netCDF/v.modisD04_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD04_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD12_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD12_npp.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f13.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD27_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD27_j1.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-06.TauCoeff.nc
SpcCoeff/netCDF/cris399_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris399_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD04_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD03_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD03_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n12.TauCoeff.nc
SpcCoeff/netCDF/atsr2_ers2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atsr2_ers2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD02_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD02_j1.TauCoeff.nc
SpcCoeff/netCDF/mi-m_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-m_coms.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD03_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD03_coms.TauCoeff.nc
SpcCoeff/netCDF/ivissr_fy2e.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ivissr_fy2e.TauCoeff.nc
SpcCoeff/netCDF/mwhs_fy3b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwhs_fy3b.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD24_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD24_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD28_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD28_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g15.TauCoeff.nc
SpcCoeff/netCDF/viirs-i_j2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-i_j2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD17_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD17_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD6_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD6_g11.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD22_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD22_j1.TauCoeff.nc
SpcCoeff/netCDF/mviriNOM_m03.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriNOM_m03.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD06_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD06_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD01_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD01_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD10_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10_aqua.TauCoeff.nc
SpcCoeff/netCDF/imgrD1_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD1_g14.TauCoeff.nc
SpcCoeff/netCDF/mersiD14_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD14_fy3a.TauCoeff.nc
SpcCoeff/netCDF/gxs_geoxo_lw.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gxs_geoxo_lw.TauCoeff.nc
SpcCoeff/netCDF/tmi_trmm.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tmi_trmm.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n06.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s02_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s02_v4.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD09_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD09_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD7_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD7_g13.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD05_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD05_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD03Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n09.TauCoeff.nc
SpcCoeff/netCDF/atms-ng_v1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms-ng_v1.TauCoeff.nc
SpcCoeff/netCDF/gxs_geoxo_mw.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gxs_geoxo_mw.TauCoeff.nc
SpcCoeff/netCDF/fci_mtg-i1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/fci_mtg-i1.TauCoeff.nc
SpcCoeff/netCDF/cris399_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris399_n20.TauCoeff.nc
SpcCoeff/netCDF/mviriBKUP_m04.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriBKUP_m04.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f18.v1-boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f18.v1-boxcar.TauCoeff.nc
SpcCoeff/netCDF/modisD05Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD04_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD04_j1.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g11.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD20_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD20_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD05_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD02_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD02_npp.TauCoeff.nc
SpcCoeff/netCDF/mersiD01_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD01_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD4_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD4_g15.TauCoeff.nc
SpcCoeff/netCDF/avhrr3JM_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3JM_n17.TauCoeff.nc
SpcCoeff/netCDF/u.omps-npAllFOV_j2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/u.omps-npAllFOV_j2.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g13.TauCoeff.nc
SpcCoeff/netCDF/avhrr3JM_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3JM_n18.TauCoeff.nc
SpcCoeff/netCDF/sndr_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g14.TauCoeff.nc
SpcCoeff/netCDF/msu_n06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n06.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g14.TauCoeff.nc
SpcCoeff/netCDF/tempest-D_cubesat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tempest-D_cubesat.TauCoeff.nc
SpcCoeff/netCDF/msu_n09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n09.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD10_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD10_npp.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB1_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB1_n20.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f20.TauCoeff.nc
SpcCoeff/netCDF/mersiD39_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD39_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD06_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD06_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD4_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD4_g13.TauCoeff.nc
SpcCoeff/netCDF/modisD07S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07S_terra.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n10.TauCoeff.nc
SpcCoeff/netCDF/imgrD2_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD2_g14.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g12.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD14_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD14_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD18_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD18_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD04Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/modisS_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisS_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD12_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD12_j1.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_metop-c.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD5_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD5_g11.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD32_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD32_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD08_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD08_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD16_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD16_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD07_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD07_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD03_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03_aqua.TauCoeff.nc
SpcCoeff/netCDF/sndr_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g12.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_mt2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_mt2.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g15.TauCoeff.nc
SpcCoeff/netCDF/iasi316_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi316_metop-b.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD27_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD27_npp.TauCoeff.nc
SpcCoeff/netCDF/mersiD12_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD12_fy3a.TauCoeff.nc
SpcCoeff/netCDF/mviriNOM_m05.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriNOM_m05.TauCoeff.nc
SpcCoeff/netCDF/msu_n10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n10.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD13_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD13_npp.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB2_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB2_n20.TauCoeff.nc
SpcCoeff/netCDF/vtprS3_itos.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vtprS3_itos.TauCoeff.nc
SpcCoeff/netCDF/v.modisD07_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD07_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD02Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/ivissr_fy2c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ivissr_fy2c.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n14.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g13.TauCoeff.nc
SpcCoeff/netCDF/ssmt2_f14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmt2_f14.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD14_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD14_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD01_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD01_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD30_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD30_j1.TauCoeff.nc
SpcCoeff/netCDF/iasiB2_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB2_metop-b.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD05_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD05_coms.TauCoeff.nc
SpcCoeff/netCDF/mi-lD03_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-lD03_coms.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD7_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD7_g15.TauCoeff.nc
SpcCoeff/netCDF/iasi300_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi300_metop-a.TauCoeff.nc
SpcCoeff/netCDF/vas_g05.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vas_g05.TauCoeff.nc
SpcCoeff/netCDF/amsua_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_metop-c.TauCoeff.nc
SpcCoeff/netCDF/imgr_mt1r.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_mt1r.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g08.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD10_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD10_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n07.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g15.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD11_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD11_j1.TauCoeff.nc
SpcCoeff/netCDF/v.ivissr_fy2d.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ivissr_fy2d.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n08.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD04_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD04_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD08_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD08_npp.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f19.TauCoeff.nc
SpcCoeff/netCDF/modisD02_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02_aqua.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f16.TauCoeff.nc
SpcCoeff/netCDF/mersiD06_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD06_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g09.TauCoeff.nc
SpcCoeff/netCDF/modisD01Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-m_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-m_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD15_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD15_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD09_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09_terra.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_insat-3ds.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_insat-3ds.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD16_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD16_npp.TauCoeff.nc
SpcCoeff/netCDF/mersiD38_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD38_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD06Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/mhs_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_metop-a.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f15.TauCoeff.nc
SpcCoeff/netCDF/v.modisD08_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD08_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD29_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD29_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD25_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD25_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD10Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/saphir_meghat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/saphir_meghat.TauCoeff.nc
SpcCoeff/netCDF/modisD02_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02_terra.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD04_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD04_coms.TauCoeff.nc
SpcCoeff/netCDF/iasi616_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi616_metop-c.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g13.TauCoeff.nc
SpcCoeff/netCDF/mi-lD02_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-lD02_coms.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD17_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD17_j1.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f18.v2-ff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f18.v2-ff_passband.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD09_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD09_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD11_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD11_npp.TauCoeff.nc
SpcCoeff/netCDF/hirs3_n15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3_n15.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-i_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-i_npp.TauCoeff.nc
SpcCoeff/netCDF/mhs_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_n18.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD29_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD29_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modis_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modis_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD13_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD13_j1.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr2_n14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr2_n14.TauCoeff.nc
SpcCoeff/netCDF/imgr_g08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g08.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD03_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD03_npp.TauCoeff.nc
SpcCoeff/netCDF/vtprS2_itos.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vtprS2_itos.TauCoeff.nc
SpcCoeff/netCDF/v.modisD06_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD06_aqua.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.nasa.T20_Vlow_LOpri.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.nasa.T20_Vlow_LOpri.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g11.TauCoeff.nc
SpcCoeff/netCDF/modisD06S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06S_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD3_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD3_g12.TauCoeff.nc
SpcCoeff/netCDF/seviri_m08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/seviri_m08.TauCoeff.nc
SpcCoeff/netCDF/mersiD13_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD13_fy3a.TauCoeff.nc
SpcCoeff/netCDF/ahi_himawari8.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ahi_himawari8.TauCoeff.nc
SpcCoeff/netCDF/v.modisD03_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD03_terra.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB2_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB2_npp.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_n16.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_n16.TauCoeff.nc
SpcCoeff/netCDF/airs324_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airs324_aqua.TauCoeff.nc
SpcCoeff/netCDF/ssu_n11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n11.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_n19.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_n19.TauCoeff.nc
SpcCoeff/netCDF/masc_cubesat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/masc_cubesat.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD6_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD6_g15.TauCoeff.nc
SpcCoeff/netCDF/mwts_fy3b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwts_fy3b.TauCoeff.nc
SpcCoeff/netCDF/mviriNOM_m07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriNOM_m07.TauCoeff.nc
SpcCoeff/netCDF/msu_n12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n12.TauCoeff.nc
SpcCoeff/netCDF/viirs-m_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-m_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD10_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD10_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD20_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD20_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD04_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD04_j1.TauCoeff.nc
SpcCoeff/netCDF/iasi-ng_metop-sg-a1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi-ng_metop-sg-a1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD05_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD05_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g11.TauCoeff.nc
SpcCoeff/netCDF/vtprS1_itos.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vtprS1_itos.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g08.TauCoeff.nc
SpcCoeff/netCDF/mersiD10_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD10_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.modisD05_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD05_terra.TauCoeff.nc
SpcCoeff/netCDF/modisD08Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/iasi300_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi300_metop-c.TauCoeff.nc
SpcCoeff/netCDF/vas_g07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vas_g07.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB3_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB3_n20.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD12_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD12_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD04_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04_terra.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD07_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD07_coms.TauCoeff.nc
SpcCoeff/netCDF/mi-lD01_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-lD01_coms.TauCoeff.nc
SpcCoeff/netCDF/amsua_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_metop-a.TauCoeff.nc
SpcCoeff/netCDF/v.mi-l_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-l_coms.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n12.TauCoeff.nc
SpcCoeff/netCDF/iras_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iras_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD4_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD4_g11.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f19.v2-ff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f19.v2-ff_passband.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD02_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD02_j1.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_g12.TauCoeff.nc
SpcCoeff/netCDF/gmi_gpm.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gmi_gpm.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD15_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD15_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD19_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD19_npp.TauCoeff.nc
SpcCoeff/netCDF/gmi_gpm.TNOM.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gmi_gpm.TNOM.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD26_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD26_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g15.TauCoeff.nc
SpcCoeff/netCDF/v.metimage_metop-sg-a1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.metimage_metop-sg-a1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD06_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD06_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD22_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD22_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD01_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.modisD10_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD10_aqua.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_metop-a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD07_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD07_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD5_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD5_g13.TauCoeff.nc
SpcCoeff/netCDF/modisD01Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/virr_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/virr_fy3a.TauCoeff.nc
SpcCoeff/netCDF/hirs4_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4_n19.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g10.TauCoeff.nc
SpcCoeff/netCDF/mersiD05_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD05_fy3a.TauCoeff.nc
SpcCoeff/netCDF/airsM9_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM9_aqua.TauCoeff.nc
SpcCoeff/netCDF/sndr_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g10.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n06.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n09.TauCoeff.nc
SpcCoeff/netCDF/hirs3_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3_n17.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD2_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD2_g12.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f08.TauCoeff.nc
SpcCoeff/netCDF/giirsB1_fsr_fy4a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/giirsB1_fsr_fy4a.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD8_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD8_g15.TauCoeff.nc
SpcCoeff/netCDF/ami_gk2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ami_gk2.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g11.TauCoeff.nc
SpcCoeff/netCDF/iasi616_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi616_metop-a.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD02_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD02_npp.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD06_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD06_coms.TauCoeff.nc
SpcCoeff/netCDF/cris-fsrB3_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsrB3_npp.TauCoeff.nc
SpcCoeff/netCDF/u.omps-tcAllFOV_j2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/u.omps-tcAllFOV_j2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD10_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD10_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD07_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD07_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD19_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD19_j1.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g08.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g13.TauCoeff.nc
SpcCoeff/netCDF/windsat_coriolis.SpcCoeff.nc
TauCoeff/ODPS/netCDF/windsat_coriolis.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g14.TauCoeff.nc
SpcCoeff/netCDF/mersiD11_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD11_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.modisD04_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD04_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD03_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD03_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD27_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD27_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD04_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD04_fy3a.TauCoeff.nc
SpcCoeff/netCDF/airsM8_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM8_aqua.TauCoeff.nc
SpcCoeff/netCDF/hsb_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hsb_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD17_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD17_npp.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs05.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs05.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g15.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-i_j2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-i_j2.TauCoeff.nc
SpcCoeff/netCDF/msu_tirosn.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_tirosn.TauCoeff.nc
SpcCoeff/netCDF/v.ivissr_fy2f.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ivissr_fy2f.TauCoeff.nc
SpcCoeff/netCDF/dpr_gpm.SpcCoeff.nc
TauCoeff/ODPS/netCDF/dpr_gpm.TauCoeff.nc
SpcCoeff/netCDF/amsub_n16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsub_n16.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD28_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD28_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD24_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD24_npp.TauCoeff.nc
SpcCoeff/netCDF/cowvr_ors6.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cowvr_ors6.TauCoeff.nc
SpcCoeff/netCDF/modisShift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisShift2_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD05_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD05_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD09_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD09_npp.TauCoeff.nc
SpcCoeff/netCDF/mhs_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_metop-c.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-02.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-02.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_g14.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD25_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD25_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD01_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD01_j1.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f19.v1-boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f19.v1-boxcar.TauCoeff.nc
SpcCoeff/netCDF/v.mi-m_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-m_coms.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD05_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD05_j1.TauCoeff.nc
SpcCoeff/netCDF/imgr_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g15.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD01_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD01_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD3_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD3_g13.TauCoeff.nc
SpcCoeff/netCDF/ahi_himawari9.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ahi_himawari9.TauCoeff.nc
SpcCoeff/netCDF/seviri_m09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/seviri_m09.TauCoeff.nc
SpcCoeff/netCDF/airsM4d_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM4d_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD10Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g10.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD21_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD21_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD05_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD05_j1.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_n18.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_n18.TauCoeff.nc
SpcCoeff/netCDF/modisD06Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD01_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD01_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD05Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/mersiD08_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD08_fy3a.TauCoeff.nc
SpcCoeff/netCDF/crisB1_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB1_n20.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_n17.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_n17.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr431_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr431_n20.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD32_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD32_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD2_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD2_g11.TauCoeff.nc
SpcCoeff/netCDF/mersiD36_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD36_fy3a.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s01_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s01_v4.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g12.TauCoeff.nc
SpcCoeff/netCDF/modisD01S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01S_terra.TauCoeff.nc
SpcCoeff/netCDF/iasi616_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi616_metop-b.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD13_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD13_npp.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.nasa.T-10_Vnom_LOpri.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.nasa.T-10_Vnom_LOpri.TauCoeff.nc
SpcCoeff/netCDF/imgr_insat-3ds.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_insat-3ds.TauCoeff.nc
SpcCoeff/netCDF/imgr_g09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g09.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD20_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD20_npp.TauCoeff.nc
SpcCoeff/netCDF/mhs_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_n19.TauCoeff.nc
SpcCoeff/netCDF/modisD04Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD10S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10S_aqua.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-01.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-01.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f14.TauCoeff.nc
SpcCoeff/netCDF/atms_n21-srf.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_n21-srf.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD14_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD14_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD03_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD03_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD40_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD40_fy3a.TauCoeff.nc
SpcCoeff/netCDF/miras_smos.SpcCoeff.nc
TauCoeff/ODPS/netCDF/miras_smos.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD03_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD03_coms.TauCoeff.nc
SpcCoeff/netCDF/mersiD23_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD23_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD05S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05S_aqua.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f20.v2-ff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f20.v2-ff_passband.TauCoeff.nc
SpcCoeff/netCDF/mwi_metop-sg-a1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwi_metop-sg-a1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD07_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD07_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD23_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD23_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n09.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g14.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n06.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-i_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-i_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD08_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD08_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD06_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD06_npp.TauCoeff.nc
SpcCoeff/netCDF/v.ivissr_fy2e.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ivissr_fy2e.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g08.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g08.TauCoeff.nc
SpcCoeff/netCDF/amsub_n15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsub_n15.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f17.TauCoeff.nc
SpcCoeff/netCDF/cris374_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris374_n20.TauCoeff.nc
SpcCoeff/netCDF/modisD01_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01_terra.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f18.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD7_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD7_g14.TauCoeff.nc
SpcCoeff/netCDF/vas_g04.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vas_g04.TauCoeff.nc
SpcCoeff/netCDF/iasiB2_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB2_metop-c.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD03_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD03_npp.TauCoeff.nc
SpcCoeff/netCDF/mersi_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersi_fy3a.TauCoeff.nc
SpcCoeff/netCDF/mersiD37_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD37_fy3a.TauCoeff.nc
SpcCoeff/netCDF/hirs2_tirosn.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_tirosn.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s06_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s06_v4.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g09.TauCoeff.nc
SpcCoeff/netCDF/modisD04Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD30_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD30_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-i_n21.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-i_n21.TauCoeff.nc
SpcCoeff/netCDF/imgrD1_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD1_g13.TauCoeff.nc
SpcCoeff/netCDF/amsua_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_metop-b.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr431_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr431_npp.TauCoeff.nc
SpcCoeff/netCDF/mviriNOM_m04.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriNOM_m04.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD06_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD06_j1.TauCoeff.nc
SpcCoeff/netCDF/mwts_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwts_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD18_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD18_j1.TauCoeff.nc
SpcCoeff/netCDF/msu_n11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n11.TauCoeff.nc
SpcCoeff/netCDF/crisB1_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB1_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-m_j2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-m_j2.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD11_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD11_npp.TauCoeff.nc
SpcCoeff/netCDF/v.aster_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.aster_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD22_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD22_npp.TauCoeff.nc
SpcCoeff/netCDF/ssmt2_f15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmt2_f15.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g12.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD02_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD02_j1.TauCoeff.nc
SpcCoeff/netCDF/airsM5_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM5_aqua.TauCoeff.nc
SpcCoeff/netCDF/mersiD09_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD09_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD26_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD26_j1.TauCoeff.nc
SpcCoeff/netCDF/v.modisD09_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD09_aqua.TauCoeff.nc
SpcCoeff/netCDF/imgr_mt2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_mt2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD16_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD16_npp.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_metop-b.TauCoeff.nc
SpcCoeff/netCDF/cris374_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris374_npp.TauCoeff.nc
SpcCoeff/netCDF/modisShift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisShift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g13.TauCoeff.nc
SpcCoeff/netCDF/viirs-m_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-m_n21.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g14.TauCoeff.nc
SpcCoeff/netCDF/iasi316_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi316_metop-c.TauCoeff.nc
SpcCoeff/netCDF/sndr_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g13.TauCoeff.nc
SpcCoeff/netCDF/ssmt1_f13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmt1_f13.TauCoeff.nc
SpcCoeff/netCDF/atms_j2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_j2.TauCoeff.nc
SpcCoeff/netCDF/airs_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airs_g13.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n11.TauCoeff.nc
SpcCoeff/netCDF/iras_fy3b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iras_fy3b.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD4_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD4_g12.TauCoeff.nc
SpcCoeff/netCDF/mi-mD04_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-mD04_coms.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD02_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD02_coms.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD08_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD08_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_g11.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD04_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD04_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD24_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD24_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD02Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD04_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD04_j1.TauCoeff.nc
SpcCoeff/netCDF/imgrD2_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD2_g15.TauCoeff.nc
SpcCoeff/netCDF/mviriBKUP_m03.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriBKUP_m03.TauCoeff.nc
SpcCoeff/netCDF/airs_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airs_aqua.TauCoeff.nc
SpcCoeff/netCDF/mersiD22_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD22_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD07_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD07_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_g15.TauCoeff.nc
SpcCoeff/netCDF/mersiD21_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD21_fy3a.TauCoeff.nc
SpcCoeff/netCDF/mhs_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mhs_metop-b.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-03.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-03.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD01_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD01_coms.TauCoeff.nc
SpcCoeff/netCDF/v.mi-lD08_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-lD08_coms.TauCoeff.nc
SpcCoeff/netCDF/amsr2_gcom-w1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsr2_gcom-w1.TauCoeff.nc
SpcCoeff/netCDF/imgr_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g14.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.nasa.T50_Vnom_LOpri.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.nasa.T50_Vnom_LOpri.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g14.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs04.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs04.TauCoeff.nc
SpcCoeff/netCDF/slstr_sentinel3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/slstr_sentinel3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD13_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD13_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD15_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD15_npp.TauCoeff.nc
SpcCoeff/netCDF/atms_n20-srf.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_n20-srf.TauCoeff.nc
SpcCoeff/netCDF/modisShift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisShift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD09_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD09_j1.TauCoeff.nc
SpcCoeff/netCDF/giirsB2_fsr_fy4a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/giirsB2_fsr_fy4a.TauCoeff.nc
SpcCoeff/netCDF/amsub_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsub_n17.TauCoeff.nc
SpcCoeff/netCDF/modis_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modis_aqua.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f16.v1-boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f16.v1-boxcar.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD12_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD12_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g09.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g12.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD3_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD3_g11.TauCoeff.nc
SpcCoeff/netCDF/airsM6_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM6_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD06S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06S_terra.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n14.TauCoeff.nc
SpcCoeff/netCDF/hamsr_grip.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hamsr_grip.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_n15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_n15.TauCoeff.nc
SpcCoeff/netCDF/ssu_pseudo.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_pseudo.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD21_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD21_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g15.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD15_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD15_j1.TauCoeff.nc
SpcCoeff/netCDF/modisShift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisShift1_terra.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g10.TauCoeff.nc
SpcCoeff/netCDF/gmi_gpm.TLO.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gmi_gpm.TLO.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD31_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD31_j1.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD8_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD8_g14.TauCoeff.nc
SpcCoeff/netCDF/hirs3_n16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3_n16.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD2_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD2_g13.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD11_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD11_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD34_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD34_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD10_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD10_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs4_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4_n18.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g11.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD09_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD09_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD05_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD05_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD5_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD5_g12.TauCoeff.nc
SpcCoeff/netCDF/modisD07_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07_terra.TauCoeff.nc
SpcCoeff/netCDF/airs281_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airs281_aqua.TauCoeff.nc
SpcCoeff/netCDF/sndr_g11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g11.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n07.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD14_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD14_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD08Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/iasi316_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi316_metop-a.TauCoeff.nc
SpcCoeff/netCDF/v.modisD06_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD06_terra.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_g13.TauCoeff.nc
SpcCoeff/netCDF/mersiD20_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD20_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD03S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03S_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g14.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD10_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD10_npp.TauCoeff.nc
SpcCoeff/netCDF/iasiB2_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB2_metop-a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD16_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD16_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD08_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD08_j1.TauCoeff.nc
SpcCoeff/netCDF/modisD04Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/iasi300_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi300_metop-b.TauCoeff.nc
SpcCoeff/netCDF/vas_g06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vas_g06.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.nasa.T20_Vhi_LOpri.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.nasa.T20_Vhi_LOpri.TauCoeff.nc
SpcCoeff/netCDF/mersiD35_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD35_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD28_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD28_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD12_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD12_j1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD23_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD23_npp.TauCoeff.nc
SpcCoeff/netCDF/airsM7_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM7_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD05Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD08S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08S_aqua.TauCoeff.nc
SpcCoeff/netCDF/mviriNOM_m06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriNOM_m06.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD02_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD02_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD6_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD6_g14.TauCoeff.nc
SpcCoeff/netCDF/ssu_tirosn.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_tirosn.TauCoeff.nc
SpcCoeff/netCDF/v.seviri_m10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.seviri_m10.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD31_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD31_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g09.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g10.TauCoeff.nc
SpcCoeff/netCDF/atms_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_n20.TauCoeff.nc
SpcCoeff/netCDF/airsM2b_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM2b_aqua.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s04_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s04_v4.TauCoeff.nc
SpcCoeff/netCDF/modisD09_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09_aqua.TauCoeff.nc
SpcCoeff/netCDF/abi_g17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_g17.TauCoeff.nc
SpcCoeff/netCDF/abi_g18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_g18.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f16.v3-eff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f16.v3-eff_passband.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_n21.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_n21.TauCoeff.nc
SpcCoeff/netCDF/modisD10Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/vissrDetB_gms5.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vissrDetB_gms5.TauCoeff.nc
SpcCoeff/netCDF/airsM4a_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM4a_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD14_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD14_npp.TauCoeff.nc
SpcCoeff/netCDF/hirs2-UWS_n11.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2-UWS_n11.TauCoeff.nc
SpcCoeff/netCDF/ivissr_fy2f.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ivissr_fy2f.TauCoeff.nc
SpcCoeff/netCDF/v.modisD01_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD01_terra.TauCoeff.nc
SpcCoeff/netCDF/sndr_g08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g08.TauCoeff.nc
SpcCoeff/netCDF/mwhs_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwhs_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-i_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-i_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD24_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD24_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD06_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD06_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD03S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03S_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD04_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD04_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD33_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD33_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD6_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD6_g12.TauCoeff.nc
SpcCoeff/netCDF/modisD04S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04S_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.seviri_m09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.seviri_m09.TauCoeff.nc
SpcCoeff/netCDF/cris-fsr_j2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris-fsr_j2.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD32_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD32_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g10.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g09.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g09.TauCoeff.nc
SpcCoeff/netCDF/iasiB1_metop-c.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB1_metop-c.TauCoeff.nc
SpcCoeff/netCDF/mersiD26_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD26_fy3a.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_n18.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD06_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD06_coms.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_n17.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD01_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD01_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD5_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD5_g14.TauCoeff.nc
SpcCoeff/netCDF/modisD07Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD18_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD18_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD06_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD06_j1.TauCoeff.nc
SpcCoeff/netCDF/mviriBKUP_m07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriBKUP_m07.TauCoeff.nc
SpcCoeff/netCDF/madras_meghat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/madras_meghat.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g12.TauCoeff.nc
SpcCoeff/netCDF/mi-l_coms.v2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-l_coms.v2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD20_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD20_npp.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_g15.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD13_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD13_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD09Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/crisB3_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB3_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD06Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgr_insat-3ds.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgr_insat-3ds.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD26_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD26_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD02_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD02_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD18_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD18_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD03_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD03_j1.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.nasa.T20_Vnom_LOpri.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.nasa.T20_Vnom_LOpri.TauCoeff.nc
SpcCoeff/netCDF/amsua_n18.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n18.TauCoeff.nc
SpcCoeff/netCDF/atms_j2-srf.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_j2-srf.TauCoeff.nc
SpcCoeff/netCDF/v.iras_fy3b.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.iras_fy3b.TauCoeff.nc
SpcCoeff/netCDF/ssu_n14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n14.TauCoeff.nc
SpcCoeff/netCDF/amsua_n17.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n17.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g13.TauCoeff.nc
SpcCoeff/netCDF/hirs4-UWS_n19.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4-UWS_n19.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f17.v2-ff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f17.v2-ff_passband.TauCoeff.nc
SpcCoeff/netCDF/sndrD2_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD2_g14.TauCoeff.nc
SpcCoeff/netCDF/mersiD32_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD32_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD16_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD16_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD23_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD23_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n12.TauCoeff.nc
SpcCoeff/netCDF/v.avhrr3_metop-b.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.avhrr3_metop-b.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD07_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD07_j1.TauCoeff.nc
SpcCoeff/netCDF/radiometer_smap.SpcCoeff.nc
TauCoeff/ODPS/netCDF/radiometer_smap.TauCoeff.nc
SpcCoeff/netCDF/modisD10S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10S_terra.TauCoeff.nc
SpcCoeff/netCDF/modisD08_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD04_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD04_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD08_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD08_npp.TauCoeff.nc
SpcCoeff/netCDF/tms_tomorrow-s03_v4.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tomorrow-s03_v4.TauCoeff.nc
SpcCoeff/netCDF/airsM12_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM12_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD2_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD2_g15.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD8_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD8_g12.TauCoeff.nc
SpcCoeff/netCDF/atms_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atms_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD30_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD30_npp.TauCoeff.nc
SpcCoeff/netCDF/v.ahi_himawari9.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ahi_himawari9.TauCoeff.nc
SpcCoeff/netCDF/imgr_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g12.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD05_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD05_j1.TauCoeff.nc
SpcCoeff/netCDF/crisB3_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB3_n20.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD21_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD21_j1.TauCoeff.nc
SpcCoeff/netCDF/hirs3-UWS_n16.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs3-UWS_n16.TauCoeff.nc
SpcCoeff/netCDF/mersiD19_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD19_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.abi_gr.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_gr.TauCoeff.nc
SpcCoeff/netCDF/amsua_metop-a_v2.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_metop-a_v2.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD01_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD01_j1.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-05.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-05.TauCoeff.nc
SpcCoeff/netCDF/ssmi_f10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmi_f10.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_g13.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_g13.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD03_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD03_npp.TauCoeff.nc
SpcCoeff/netCDF/geostorm_proposed.gaussian.SpcCoeff.nc
TauCoeff/ODPS/netCDF/geostorm_proposed.gaussian.TauCoeff.nc
SpcCoeff/netCDF/mersiD27_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD27_fy3a.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD22_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD22_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-m_j2.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-m_j2.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-dnb_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-dnb_j1.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs02.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs02.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g10.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g12.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g12.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD11_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD11_npp.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD07_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD07_coms.TauCoeff.nc
SpcCoeff/netCDF/mi-mD01_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-mD01_coms.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f16.v2-ff_passband.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f16.v2-ff_passband.TauCoeff.nc
SpcCoeff/netCDF/sndr_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndr_g15.TauCoeff.nc
SpcCoeff/netCDF/modisD03Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD03Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD21_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD21_npp.TauCoeff.nc
SpcCoeff/netCDF/ssmt1_f15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmt1_f15.TauCoeff.nc
SpcCoeff/netCDF/v.modis_aqua.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modis_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD2_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD2_g12.TauCoeff.nc
SpcCoeff/netCDF/mwi_wsf-m1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwi_wsf-m1.TauCoeff.nc
SpcCoeff/netCDF/iasiB1_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB1_metop-a.TauCoeff.nc
SpcCoeff/netCDF/gmi_gpm.THI.SpcCoeff.nc
TauCoeff/ODPS/netCDF/gmi_gpm.THI.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD04_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD04_coms.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD08_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD08_j1.TauCoeff.nc
SpcCoeff/netCDF/mi-mD02_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-mD02_coms.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD16_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD16_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD12_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD12_j1.TauCoeff.nc
SpcCoeff/netCDF/msu_n08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n08.TauCoeff.nc
SpcCoeff/netCDF/mersiD24_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD24_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD28_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD28_j1.TauCoeff.nc
SpcCoeff/netCDF/avhrr3_n15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr3_n15.TauCoeff.nc
SpcCoeff/netCDF/sndrD4_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD4_g15.TauCoeff.nc
SpcCoeff/netCDF/crisB2_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB2_npp.TauCoeff.nc
SpcCoeff/netCDF/msu_n07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/msu_n07.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD12_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD12_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-i_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-i_npp.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD3_g10.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD3_g10.TauCoeff.nc
SpcCoeff/netCDF/mwri_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mwri_fy3a.TauCoeff.nc
SpcCoeff/netCDF/imgrD2_g13.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD2_g13.TauCoeff.nc
SpcCoeff/netCDF/mviriBKUP_m05.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mviriBKUP_m05.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD4_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD4_g14.TauCoeff.nc
SpcCoeff/netCDF/modis_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modis_terra.TauCoeff.nc
SpcCoeff/netCDF/cris_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris_n20.TauCoeff.nc
SpcCoeff/netCDF/eon_mw.v1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/eon_mw.v1.TauCoeff.nc
SpcCoeff/netCDF/airsM4c_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM4c_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD10_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD10_j1.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f20.v1-boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f20.v1-boxcar.TauCoeff.nc
SpcCoeff/netCDF/hirs4-UWS_metop-a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4-UWS_metop-a.TauCoeff.nc
SpcCoeff/netCDF/imgrD1_g15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgrD1_g15.TauCoeff.nc
SpcCoeff/netCDF/abi_gr.SpcCoeff.nc
TauCoeff/ODPS/netCDF/abi_gr.TauCoeff.nc
SpcCoeff/netCDF/airsM11_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM11_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD04S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD04S_terra.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n08.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n08.TauCoeff.nc
SpcCoeff/netCDF/airsM3_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM3_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD7_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD7_g12.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD14_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD14_j1.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD07_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD07_npp.TauCoeff.nc
SpcCoeff/netCDF/avhrr2_n07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/avhrr2_n07.TauCoeff.nc
SpcCoeff/netCDF/mersiD31_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD31_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD1_g14.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD1_g14.TauCoeff.nc
SpcCoeff/netCDF/ivissr_fy2d.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ivissr_fy2d.TauCoeff.nc
SpcCoeff/netCDF/ssmis_f17.v1-boxcar.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssmis_f17.v1-boxcar.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD15_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD15_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-m_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-m_npp.TauCoeff.nc
SpcCoeff/netCDF/modisD05Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD05Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/modisD06Shift1_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06Shift1_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD08S_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD08S_terra.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD23_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD23_npp.TauCoeff.nc
SpcCoeff/netCDF/imgr_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/imgr_g10.TauCoeff.nc
SpcCoeff/netCDF/modisD01Shift2_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD01Shift2_terra.TauCoeff.nc
SpcCoeff/netCDF/modisD09S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD09S_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD1_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD1_g11.TauCoeff.nc
SpcCoeff/netCDF/tms_tropics-07.SpcCoeff.nc
TauCoeff/ODPS/netCDF/tms_tropics-07.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD10_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD10_npp.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD31_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD31_j1.TauCoeff.nc
SpcCoeff/netCDF/v.mi-mD05_coms.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.mi-mD05_coms.TauCoeff.nc
SpcCoeff/netCDF/iasiB3_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasiB3_metop-b.TauCoeff.nc
SpcCoeff/netCDF/mi-mD03_coms.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mi-mD03_coms.TauCoeff.nc
SpcCoeff/netCDF/crisB2_n20.SpcCoeff.nc
TauCoeff/ODPS/netCDF/crisB2_n20.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD15_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD15_j1.TauCoeff.nc
SpcCoeff/netCDF/ssu_n09.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n09.TauCoeff.nc
SpcCoeff/netCDF/modisD07Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD31_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD31_npp.TauCoeff.nc
SpcCoeff/netCDF/ssu_n06.SpcCoeff.nc
TauCoeff/ODPS/netCDF/ssu_n06.TauCoeff.nc
SpcCoeff/netCDF/seviri_m10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/seviri_m10.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD02_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD02_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-iD11_j1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-iD11_j1.TauCoeff.nc
SpcCoeff/netCDF/mersiD25_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD25_fy3a.TauCoeff.nc
SpcCoeff/netCDF/v.sndr_g12.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndr_g12.TauCoeff.nc
SpcCoeff/netCDF/sndrD1_g10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD1_g10.TauCoeff.nc
SpcCoeff/netCDF/v.ivissr_fy2c.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.ivissr_fy2c.TauCoeff.nc
SpcCoeff/netCDF/cpr_cloudsat.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cpr_cloudsat.TauCoeff.nc
SpcCoeff/netCDF/micromas_cs00.SpcCoeff.nc
TauCoeff/ODPS/netCDF/micromas_cs00.TauCoeff.nc
SpcCoeff/netCDF/hirs4_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs4_metop-b.TauCoeff.nc
SpcCoeff/netCDF/amsua_n15.SpcCoeff.nc
TauCoeff/ODPS/netCDF/amsua_n15.TauCoeff.nc
SpcCoeff/netCDF/v.sndrD4_g11.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.sndrD4_g11.TauCoeff.nc
SpcCoeff/netCDF/v.abi_g17.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_g17.TauCoeff.nc
SpcCoeff/netCDF/v.modisD07_terra.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.modisD07_terra.TauCoeff.nc
SpcCoeff/netCDF/v.abi_g18.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.abi_g18.TauCoeff.nc
SpcCoeff/netCDF/mersiD30_fy3a.SpcCoeff.nc
TauCoeff/ODPS/netCDF/mersiD30_fy3a.TauCoeff.nc
SpcCoeff/netCDF/modisD02S_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD02S_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.imgrD3_g15.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.imgrD3_g15.TauCoeff.nc
SpcCoeff/netCDF/hirs2_n10.SpcCoeff.nc
TauCoeff/ODPS/netCDF/hirs2_n10.TauCoeff.nc
SpcCoeff/netCDF/atsr1_ers1.SpcCoeff.nc
TauCoeff/ODPS/netCDF/atsr1_ers1.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD05_npp.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD05_npp.TauCoeff.nc
SpcCoeff/netCDF/viirs-mD09_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/viirs-mD09_npp.TauCoeff.nc
SpcCoeff/netCDF/airsM10_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM10_aqua.TauCoeff.nc
SpcCoeff/netCDF/vissrDetA_gms5.SpcCoeff.nc
TauCoeff/ODPS/netCDF/vissrDetA_gms5.TauCoeff.nc
SpcCoeff/netCDF/modisD10Shift2_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD10Shift2_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-iD13_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-iD13_j1.TauCoeff.nc
SpcCoeff/netCDF/airsM4b_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM4b_aqua.TauCoeff.nc
SpcCoeff/netCDF/v.viirs-mD09_j1.SpcCoeff.nc
TauCoeff/ODAS/netCDF/v.viirs-mD09_j1.TauCoeff.nc
SpcCoeff/netCDF/airsM2a_aqua.SpcCoeff.nc
TauCoeff/ODPS/netCDF/airsM2a_aqua.TauCoeff.nc
SpcCoeff/netCDF/modisD06_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD06_terra.TauCoeff.nc
SpcCoeff/netCDF/sndrD3_g14.SpcCoeff.nc
TauCoeff/ODPS/netCDF/sndrD3_g14.TauCoeff.nc
SpcCoeff/netCDF/modisD07Shift1_terra.SpcCoeff.nc
TauCoeff/ODPS/netCDF/modisD07Shift1_terra.TauCoeff.nc
SpcCoeff/netCDF/iasi_metop-b.SpcCoeff.nc
TauCoeff/ODPS/netCDF/iasi_metop-b.TauCoeff.nc
SpcCoeff/netCDF/cris_npp.SpcCoeff.nc
TauCoeff/ODPS/netCDF/cris_npp.TauCoeff.nc
)

# Symlink all CRTM files
# Version 3
CREATE_SYMLINK_FILENAME( ${CRTM_COEFFS_PATH}/${CRTM_COEFFS_BRANCH_PREFIX}/${CRTM_COEFFS_BRANCH}/fix
                         ${CMAKE_CURRENT_BINARY_DIR}/testinput
                         ${crtm_test_input} )
