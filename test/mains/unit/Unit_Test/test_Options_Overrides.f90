!
! test_Options_Overrides.f90
!
! Unit test to verify option overrides affect Forward/TL/AD/K outputs
! only when the Options argument is present.
!
PROGRAM test_Options_Overrides

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  USE CRTM_Module
  IMPLICIT NONE
  ! ============================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Options_Overrides'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS   = 1
  INTEGER, PARAMETER :: N_STREAMS   = 8
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  REAL(fp), PARAMETER :: DIFF_TOL     = 1.0e-6_fp

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: m
  REAL(fp) :: r_noopt, r_opt
  INTEGER :: n_full_opt

  REAL(fp), ALLOCATABLE :: Emissivity(:)

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_noopt(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_opt(:,:)
  TYPE(CRTM_Options_type)                 :: Opt(N_PROFILES)

  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_TL(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:)

  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_AD(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_AD(:,:)

  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type) , ALLOCATABLE   :: Surface_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)

  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, '',  &
    'CRTM Version: '//TRIM(Version) )

  Sensor_Id = 'atms_n21'
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)

  ! ============================================================================
  ! 1. **** INITIALIZE THE CRTM ****
  !
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ! ============================================================================
  ! 2. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), &
            RTSolution_noopt( n_Channels, N_PROFILES ), &
            RTSolution_opt( n_Channels, N_PROFILES ), &
            RTSolution_TL( n_Channels, N_PROFILES ), &
            RTSolution_AD( n_Channels, N_PROFILES ), &
            RTSolution_K( n_Channels, N_PROFILES ), &
            Atmosphere_K( n_Channels, N_PROFILES ), &
            Surface_K( n_Channels, N_PROFILES ), &
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
  CALL CRTM_Atmosphere_Create( Atmosphere_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  ! ============================================================================
  ! 3. **** ASSIGN INPUT DATA ****
  !
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()

  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! ============================================================================
  ! 4. **** SET UP OPTIONS ****
  !
  ALLOCATE( Emissivity( n_Channels ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating Emissivity'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Emissivity(:) = 0.5_fp

  CALL CRTM_Options_Create( Opt, n_Channels )
  IF ( .NOT. ALL(CRTM_Options_Associated( Opt )) ) THEN
    Message = 'Error allocating Options structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  DO m = 1, N_PROFILES
    CALL CRTM_Options_SetValue( Opt(m), n_Streams = N_STREAMS )
    CALL CRTM_Options_SetEmissivity( Opt(m), Emissivity )
  END DO

  ! ============================================================================
  ! 5. **** FORWARD WITHOUT OPTIONS ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model (no options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  r_noopt = RTSolution(1,1)%Radiance
  RTSolution_noopt = RTSolution

  ! ============================================================================
  ! 6. **** FORWARD WITH OPTIONS ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution , &
                               Options = Opt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model (options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  r_opt = RTSolution(1,1)%Radiance
  RTSolution_opt = RTSolution
  n_full_opt = RTSolution(1,1)%n_Full_Streams
  IF ( n_full_opt /= (Opt(1)%n_Streams + 2) ) THEN
    WRITE(*,'(a,2i6)') 'n_Full_Streams mismatch: ', n_full_opt, Opt(1)%n_Streams + 2
    STOP 1
  END IF
  IF ( ABS(r_opt - r_noopt) <= DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'Radiance not affected by options: ', r_opt, r_noopt
    STOP 1
  END IF

  ! ============================================================================
  ! 7. **** TL WITH/WITHOUT OPTIONS ****
  !
  CALL CRTM_Atmosphere_Zero( Atmosphere_TL )
  CALL CRTM_Surface_Zero( Surface_TL )
  CALL CRTM_RTSolution_Zero( RTSolution_TL )
  Error_Status = CRTM_Tangent_Linear( Atm , &
                                      Sfc , &
                                      Atmosphere_TL , &
                                      Surface_TL , &
                                      Geometry , &
                                      ChannelInfo , &
                                      RTSolution_noopt , &
                                      RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Tangent-Linear Model (no options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( ABS(RTSolution_noopt(1,1)%Radiance - r_noopt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'TL forward mismatch (no options): ', RTSolution_noopt(1,1)%Radiance, r_noopt
    STOP 1
  END IF

  CALL CRTM_RTSolution_Zero( RTSolution_TL )
  Error_Status = CRTM_Tangent_Linear( Atm , &
                                      Sfc , &
                                      Atmosphere_TL , &
                                      Surface_TL , &
                                      Geometry , &
                                      ChannelInfo , &
                                      RTSolution_opt , &
                                      RTSolution_TL , &
                                      Options = Opt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Tangent-Linear Model (options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( RTSolution(1,1)%n_Full_Streams /= (Opt(1)%n_Streams + 2) ) THEN
    WRITE(*,'(a,2i6)') 'TL n_Full_Streams mismatch: ', RTSolution(1,1)%n_Full_Streams, Opt(1)%n_Streams + 2
    STOP 1
  END IF
  IF ( ABS(RTSolution_opt(1,1)%Radiance - r_opt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'TL forward mismatch (options): ', RTSolution_opt(1,1)%Radiance, r_opt
    STOP 1
  END IF

  ! ============================================================================
  ! 8. **** AD WITH/WITHOUT OPTIONS ****
  !
  CALL CRTM_Atmosphere_Zero( Atmosphere_AD )
  CALL CRTM_Surface_Zero( Surface_AD )
  CALL CRTM_RTSolution_Zero( RTSolution_AD )
  Error_Status = CRTM_Adjoint( Atm          , &
                               Sfc          , &
                               RTSolution_AD, &
                               Geometry     , &
                               ChannelInfo  , &
                               Atmosphere_AD, &
                               Surface_AD   , &
                               RTSolution_noopt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Adjoint Model (no options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( ABS(RTSolution_noopt(1,1)%Radiance - r_noopt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'AD forward mismatch (no options): ', RTSolution_noopt(1,1)%Radiance, r_noopt
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Zero( Atmosphere_AD )
  CALL CRTM_Surface_Zero( Surface_AD )
  CALL CRTM_RTSolution_Zero( RTSolution_AD )
  Error_Status = CRTM_Adjoint( Atm          , &
                               Sfc          , &
                               RTSolution_AD, &
                               Geometry     , &
                               ChannelInfo  , &
                               Atmosphere_AD, &
                               Surface_AD   , &
                               RTSolution_opt , &
                               Options = Opt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Adjoint Model (options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( RTSolution(1,1)%n_Full_Streams /= (Opt(1)%n_Streams + 2) ) THEN
    WRITE(*,'(a,2i6)') 'AD n_Full_Streams mismatch: ', RTSolution(1,1)%n_Full_Streams, Opt(1)%n_Streams + 2
    STOP 1
  END IF
  IF ( ABS(RTSolution_opt(1,1)%Radiance - r_opt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'AD forward mismatch (options): ', RTSolution_opt(1,1)%Radiance, r_opt
    STOP 1
  END IF

  ! ============================================================================
  ! 9. **** K-MATRIX WITH/WITHOUT OPTIONS ****
  !
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Surface_Zero( Surface_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  Error_Status = CRTM_K_Matrix( Atm         , &
                                Sfc         , &
                                RTSolution_K, &
                                Geometry    , &
                                ChannelInfo , &
                                Atmosphere_K, &
                                Surface_K   , &
                                RTSolution_noopt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM K-Matrix Model (no options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( ABS(RTSolution_noopt(1,1)%Radiance - r_noopt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'K forward mismatch (no options): ', RTSolution_noopt(1,1)%Radiance, r_noopt
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Surface_Zero( Surface_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  Error_Status = CRTM_K_Matrix( Atm         , &
                                Sfc         , &
                                RTSolution_K, &
                                Geometry    , &
                                ChannelInfo , &
                                Atmosphere_K, &
                                Surface_K   , &
                                RTSolution_opt , &
                                Options = Opt )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM K-Matrix Model (options)'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  IF ( RTSolution(1,1)%n_Full_Streams /= (Opt(1)%n_Streams + 2) ) THEN
    WRITE(*,'(a,2i6)') 'K n_Full_Streams mismatch: ', RTSolution(1,1)%n_Full_Streams, Opt(1)%n_Streams + 2
    STOP 1
  END IF
  IF ( ABS(RTSolution_opt(1,1)%Radiance - r_opt) > DIFF_TOL ) THEN
    WRITE(*,'(a,2es14.6)') 'K forward mismatch (options): ', RTSolution_opt(1,1)%Radiance, r_opt
    STOP 1
  END IF

  STOP 0

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_Options_Overrides
