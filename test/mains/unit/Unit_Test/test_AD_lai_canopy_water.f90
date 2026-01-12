!
! test_AD_lai_canopy_water.f90
!
! Program to test the CRTM adjoint response for LAI and canopy water content.
!
PROGRAM test_AD_lai_canopy_water

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_AD_lai_canopy_water'
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
  REAL(fp), PARAMETER :: TOLERANCE = 1.0e-5_fp


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  REAL(fp) :: Perturb_Lai
  REAL(fp) :: Perturb_Canopy
  REAL(fp) :: Dot_TL
  REAL(fp) :: Dot_AD
  REAL(fp) :: Dot_Diff
  REAL(fp) :: Dot_Norm

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_AD(:,:)
  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_TL(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_AD(N_PROFILES)


  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Program to test LAI and canopy water AD response.', &
    'CRTM Version: '//TRIM(Version) )

  CALL GET_COMMAND_ARGUMENT(1, Sensor_Id)
  IF ( LEN_TRIM(Sensor_Id) == 0 ) Sensor_Id = 'atms_n21'
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
            RTSolution_TL( n_Channels, N_PROFILES ), &
            RTSolution_AD( n_Channels, N_PROFILES ), &
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

  CALL CRTM_Atmosphere_Create( Atmosphere_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_AD)) ) THEN
    Message = 'Error allocating CRTM Atmosphere_AD structure'
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
  Sfc(:)%Lai = 2.0_fp

  CALL CRTM_Atmosphere_Zero( Atmosphere_TL )
  CALL CRTM_Surface_Zero( Surface_TL )
  Perturb_Lai = 0.01_fp
  Perturb_Canopy = -0.02_fp
  Surface_TL(TEST_PROFILE)%Lai = Perturb_Lai
  Surface_TL(TEST_PROFILE)%Canopy_Water_Content = Perturb_Canopy

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

  CALL CRTM_Atmosphere_Zero( Atmosphere_AD )
  CALL CRTM_Surface_Zero( Surface_AD )
  RTSolution_AD%Radiance = ZERO
  RTSolution_AD%Brightness_Temperature = ZERO
  RTSolution_AD(TEST_CHANNEL,TEST_PROFILE)%Radiance = ONE

  Error_Status = CRTM_Adjoint( Atm , &
                               Sfc , &
                               RTSolution_AD, &
                               Geometry, &
                               ChannelInfo, &
                               Atmosphere_AD, &
                               Surface_AD, &
                               RTSolution  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Adjoint Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Dot_TL = RTSolution_TL(TEST_CHANNEL,TEST_PROFILE)%Radiance
  Dot_AD = (Surface_AD(TEST_PROFILE)%Lai * Surface_TL(TEST_PROFILE)%Lai) + &
           (Surface_AD(TEST_PROFILE)%Canopy_Water_Content * Surface_TL(TEST_PROFILE)%Canopy_Water_Content)
  Dot_Diff = ABS(Dot_TL - Dot_AD)
  Dot_Norm = MAX(ABS(Dot_TL), ONE)

  WRITE(*,*) 'TL dot:', Dot_TL
  WRITE(*,*) 'AD dot:', Dot_AD
  WRITE(*,*) 'Diff:', Dot_Diff

  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Destroy(Atmosphere_AD)
  CALL CRTM_Atmosphere_Destroy(Atmosphere_TL)
  CALL CRTM_Atmosphere_Destroy(Atm)

  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_AD, &
              STAT = Allocate_Status )

  IF ( ABS(Dot_TL) <= 0.0_fp ) THEN
    Message = 'TL dot product is zero for LAI/canopy perturbations'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  IF ( (Dot_Diff / Dot_Norm) < TOLERANCE ) THEN
    STOP 0
  END IF

  WRITE(*,*) 'FAIL diff/norm=', Dot_Diff / Dot_Norm, ' TOL=', TOLERANCE
  STOP 1

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_AD_lai_canopy_water
