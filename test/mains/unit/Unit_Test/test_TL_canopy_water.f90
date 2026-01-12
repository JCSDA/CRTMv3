!
! test_TL_canopy_water.f90
!
! Program to test the CRTM tangent-linear response to canopy water content.
!
PROGRAM test_TL_canopy_water

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  USE CRTM_Module
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_TL_canopy_water'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS   = 1
  INTEGER, PARAMETER :: TEST_CHANNEL = 1
  INTEGER, PARAMETER :: TEST_PROFILE = 2
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  REAL(fp), PARAMETER :: TOLERANCE = 0.2_fp


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  REAL(fp) :: Perturbation
  REAL(fp) :: Ratio

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_Perturb(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:)
  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_TL(N_PROFILES)


  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Program to test canopy water content TL response.', &
    'CRTM Version: '//TRIM(Version) )

  Sensor_Id = 'atms_n21'
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)

  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo  , &
                            File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), &
            RTSolution_Perturb( n_Channels, N_PROFILES ), &
            RTSolution_TL( n_Channels, N_PROFILES ), &
            STAT = Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Create( Atmosphere_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_TL)) ) THEN
    Message = 'Error allocating CRTM Atmosphere_TL structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()

  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  Sfc(:)%Land_Coverage  = 1.0_fp
  Sfc(:)%Water_Coverage = 0.0_fp
  Sfc(:)%Snow_Coverage  = 0.0_fp
  Sfc(:)%Ice_Coverage   = 0.0_fp
  Sfc(:)%Canopy_Water_Content = 0.15_fp

  CALL CRTM_Atmosphere_Zero( Atmosphere_TL )
  CALL CRTM_Surface_Zero( Surface_TL )
  Perturbation = 0.01_fp
  Surface_TL(TEST_PROFILE)%Canopy_Water_Content = Perturbation

  Error_Status = CRTM_Tangent_Linear( Atm , &
                                      Sfc , &
                                      Atmosphere_TL , &
                                      Surface_TL , &
                                      Geometry , &
                                      ChannelInfo , &
                                      RTSolution , &
                                      RTSolution_TL  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Tangent-linear Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Error_Status = CRTM_Forward( Atm         , &
                               Sfc         , &
                               Geometry    , &
                               ChannelInfo , &
                               RTSolution  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Sfc(TEST_PROFILE)%Canopy_Water_Content = Sfc(TEST_PROFILE)%Canopy_Water_Content + Perturbation
  Error_Status = CRTM_Forward( Atm         , &
                               Sfc         , &
                               Geometry    , &
                               ChannelInfo , &
                               RTSolution_Perturb )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in perturbed CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  IF ( ABS(RTSolution_TL(TEST_CHANNEL,TEST_PROFILE)%Radiance) <= 0.0_fp ) THEN
    WRITE(*,*) 'TL radiance:', RTSolution_TL(TEST_CHANNEL,TEST_PROFILE)%Radiance
    Message = 'TL radiance is zero for canopy water perturbation'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Ratio = ( RTSolution_Perturb(TEST_CHANNEL,TEST_PROFILE)%Radiance - &
            RTSolution(TEST_CHANNEL,TEST_PROFILE)%Radiance ) / &
            RTSolution_TL(TEST_CHANNEL,TEST_PROFILE)%Radiance
  WRITE(*,*) 'Base radiance:', RTSolution(TEST_CHANNEL,TEST_PROFILE)%Radiance
  WRITE(*,*) 'Perturb radiance:', RTSolution_Perturb(TEST_CHANNEL,TEST_PROFILE)%Radiance
  WRITE(*,*) 'TL radiance:', RTSolution_TL(TEST_CHANNEL,TEST_PROFILE)%Radiance
  WRITE(*,*) 'Ratio: ', Ratio

  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Destroy(Atmosphere_TL)
  CALL CRTM_Atmosphere_Destroy(Atm)

  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_Perturb, &
              STAT = Allocate_Status )

  IF ( ABS(1.0_fp - Ratio) < TOLERANCE ) THEN
    STOP 0
  END IF

  WRITE(*,*) 'FAIL abs(1 - Ratio)=', ABS(1.0_fp - Ratio), ' TOL=', TOLERANCE
  STOP 1

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_TL_canopy_water
