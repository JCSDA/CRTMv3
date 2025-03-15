!
! test_MWSurfEM.f90
!
! Program to provide a (relatively) simple example of how
! to test the CRTM MW surface emissivity functions.
!
! Copyright Patrick Stegmann, 2025
!
PROGRAM test_MWSurfEM

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE UnitTest_Define, ONLY: UnitTest_type,   &
                             UnitTest_Init,   &
                             UnitTest_Setup,  &
                             UnitTest_Assert, &
                             UnitTest_Passed
  USE Type_Kinds,               ONLY: fp
  USE CRTM_Surface_Define, ONLY: CRTM_Surface_type
  USE SpcCoeff_Define, ONLY: SpcCoeff_type
  USE CRTM_Geometry_Define, ONLY: CRTM_Geometry_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_SetValue, &
                                      CRTM_GeometryInfo_Destroy
  USE CRTM_SpcCoeff
  USE Message_Handler, ONLY: SUCCESS, Display_Message
  !USE CRTM_MW_Water_SfcOptics, ONLY: iVar_type, &
  !                                   Compute_MW_Water_SfcOptics
  USE CRTM_SfcOptics_Define, ONLY: CRTM_SfcOptics_type, &
                                   CRTM_SfcOptics_Create, &
                                   CRTM_SfcOptics_Destroy
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_SfcOptics, ONLY: iVar_type, &
                            CRTM_Compute_SfcOptics

  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_MWSurfEM'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/unit/'
  LOGICAL,      PARAMETER :: Quiet = .TRUE.

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS EXAMPLE ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  ! ...but only ONE Sensor at a time
  INTEGER, PARAMETER :: N_SENSORS = 1

  ! Test GeometryInfo angles. The test scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  ! ============================================================================


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256), DIMENSION(1) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: SensorIndex
  INTEGER :: ChannelIndex
  INTEGER :: ii
  INTEGER :: testresult
  TYPE(UnitTest_type) :: emTest
  LOGICAL :: testPassed


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_GeometryInfo_type), DIMENSION(N_PROFILES) :: GeometryInfo(N_PROFILES)
  TYPE(iVar_type)                         :: iVar
  ! Define the FORWARD variables
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_SfcOptics_type)               :: SfcOptics

  ! ============================================================================
  ! Initialize Unit test:
  CALL UnitTest_Init(emTest, .TRUE.)
  CALL UnitTest_Setup(emTest, 'test_MWSurfEM', Program_Name, .TRUE.)

  ! Get sensor id from user
  ! -----------------------
  SensorIndex = 1 ! Number of sensors
  !WRITE( *,'(/5x,"Enter sensor id [hirs4_n18, amsua_metop-a, or mhs_n18]: ")',ADVANCE='NO' )
  !READ( *,'(a)' ) Sensor_Id
  Sensor_Id = (/'atms_npp'/)
  Sensor_Id = ADJUSTL(Sensor_Id(SensorIndex))
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id(SensorIndex))

  ! Load the instrument spectral coefficients
  ! -----------------------
  Error_Status = 3
  Error_Status = CRTM_SpcCoeff_Load( &
                                    Sensor_ID         = Sensor_Id        , &
                                    File_Path         = COEFFICIENTS_PATH, &
                                    netCDF            = .FALSE.           , &
                                    Quiet             = Quiet          )
  CALL UnitTest_Assert(emTest, (Error_Status==SUCCESS) )

  Error_Status = CRTM_MWwaterCoeff_Load_FASTEM( &
                        'FASTEM6',             &
                        Quiet             = Quiet            )
  CALL UnitTest_Assert(emTest, (Error_Status==SUCCESS) )
  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  ChannelIndex = 1 ! Channel number
  n_Channels = SC(SensorIndex)%n_Channels
  ! ============================================================================
  WRITE(*,*) "Channels: ", n_Channels


  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !

  ! 4a. Surface input
  ! --------------------------------
  CALL Load_Sfc_Data()


  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_GeometryInfo_SetValue( GeometryInfo, &
                                   Source_Azimuth_Angle = 0._fp, &
                                   Sensor_Azimuth_Angle = 0._fp, &
                                   Sensor_Scan_Angle   = SCAN_ANGLE, & 
                                   Sensor_Zenith_Angle = ZENITH_ANGLE )
  ! ============================================================================

  CALL CRTM_SfcOptics_Create(SfcOptics, &
                              1, & ! n_Angles
                              1)   ! n_Stokes

  !Error_Status = Compute_MW_Water_SfcOptics( &
  !  Sfc(1)     , &  ! Input
  !  GeometryInfo(1), &  ! Input
  !  SensorIndex , &  ! Input
  !  ChannelIndex, &  ! Input
  !  SfcOptics   , &  ! Output
  !  iVar        )

  ChannelLoop: DO ChannelIndex = 1, n_Channels
    Error_Status = CRTM_Compute_SfcOptics( &
                                          Sfc(1)     , &  ! Input
                                          GeometryInfo(1), &  ! Input
                                          SensorIndex , &  ! Input
                                          ChannelIndex, &  ! Input
                                          SfcOptics   , &  ! Output
                                          iVar        )
    WRITE(*,*) "Channel: ", ChannelIndex, " Emissivity: ", SfcOptics%Emissivity
  END DO ChannelLoop

  CALL UnitTest_Assert(emTest, (Error_Status==SUCCESS) )

  ! 5. Cleanup
  CALL CRTM_SfcOptics_Destroy(SfcOptics)
  CALL CRTM_GeometryInfo_Destroy(GeometryInfo)

  testPassed = UnitTest_Passed(emTest)
  IF(testPassed) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

CONTAINS

  SUBROUTINE Load_Sfc_Data()
    
    
    ! 4a.0 Surface type definitions for default SfcOptics definitions
    !      For IR and VIS, this is the NPOESS reflectivities.
    ! ---------------------------------------------------------------
    INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
    INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
    INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
    INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics



    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc(1)%Land_Coverage     = 0.1_fp
    sfc(1)%Land_Type         = TUNDRA_SURFACE_TYPE
    sfc(1)%Land_Temperature  = 272.0_fp
    sfc(1)%Lai               = 0.17_fp
    sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
    sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc(1)%Water_Coverage    = 0.5_fp
    sfc(1)%Water_Type        = SEA_WATER_TYPE
    sfc(1)%Water_Temperature = 275.0_fp
    ! ...Snow coverage characteristics
    sfc(1)%Snow_Coverage    = 0.25_fp
    sfc(1)%Snow_Type        = FRESH_SNOW_TYPE
    sfc(1)%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc(1)%Ice_Coverage    = 0.15_fp
    sfc(1)%Ice_Type        = FRESH_ICE_TYPE
    sfc(1)%Ice_Temperature = 269.0_fp



    ! 4a.2 Profile #2
    ! ---------------
    ! Surface data
    sfc(2)%Land_Coverage    = 1.0_fp
    sfc(2)%Land_Type        = SCRUB_SURFACE_TYPE
    sfc(2)%Land_Temperature = 318.0_fp
    sfc(2)%Lai              = 0.65_fp
    sfc(2)%Soil_Type        = COARSE_SOIL_TYPE
    sfc(2)%Vegetation_Type  = BARE_SOIL_VEGETATION_TYPE

  END SUBROUTINE Load_Sfc_Data

END PROGRAM test_MWSurfEM
