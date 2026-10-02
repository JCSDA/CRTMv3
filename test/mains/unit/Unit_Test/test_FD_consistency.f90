!
! test_FD_consistency
!
! Finite-difference consistency check for CRTM tangent-linear output.
!

PROGRAM test_FD_consistency

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  USE CRTM_Module
  IMPLICIT NONE
  ! ============================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_FD_consistency'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0
  INTEGER, PARAMETER :: N_SENSORS   = 1

  ! Test GeometryInfo angles. The test scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! Finite-difference settings
  INTEGER,  PARAMETER :: PROFILE_IDX = 1
  INTEGER,  PARAMETER :: LAYER_IDX   = 50
  REAL(fp), PARAMETER :: DELTA_T     = 0.01_fp
  REAL(fp), PARAMETER :: RTOL        = 1.0e-3_fp
  REAL(fp), PARAMETER :: ATOL        = 1.0e-6_fp

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, k, m
  INTEGER :: n_fail, exit_status, channel_max
  INTEGER :: arg_count
  REAL(fp) :: fd, tl, diff, tol, max_diff
  REAL(fp) :: lhs, rhs, rel
  REAL(fp), PARAMETER :: AD_RTOL = 1.0e-3_fp
  REAL(fp), PARAMETER :: AD_ATOL = 1.0e-6_fp
  LOGICAL :: run_fd, run_ad
  CHARACTER(16) :: mode

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm_TL(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm_Prtb(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc_AD(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_Prtb(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_AD(:,:)

  ! First, make sure the right number of inputs have been provided
  arg_count = COMMAND_ARGUMENT_COUNT()
  IF ( arg_count < 1 .OR. arg_count > 2 ) THEN
    WRITE(*,*) TRIM(PROGRAM_NAME)//': ERROR, one or two command-line arguments required, returning'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, Sensor_Id)
  mode = 'fd'
  IF ( arg_count == 2 ) THEN
    CALL GET_COMMAND_ARGUMENT(2, mode)
    mode = ADJUSTL(mode)
  END IF
  run_fd = (TRIM(mode) == 'fd' .OR. TRIM(mode) == 'both')
  run_ad = (TRIM(mode) == 'ad' .OR. TRIM(mode) == 'both')
  IF ( .NOT. run_fd .AND. .NOT. run_ad ) THEN
    WRITE(*,*) TRIM(PROGRAM_NAME)//': ERROR, invalid mode (fd, ad, or both expected)'
    STOP 1
  END IF

  ! Program header
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Finite-difference consistency check for CRTM TL output.', &
    'CRTM Version: '//TRIM(Version) )

  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)

  ! ============================================================================
  ! 1. **** INITIALIZE THE CRTM ****
  !
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
            RTSolution_TL( n_Channels, N_PROFILES ), &
            RTSolution_Prtb( n_Channels, N_PROFILES ), &
            RTSolution_AD( n_Channels, N_PROFILES ), &
            STAT = Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_Prtb, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AD, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_TL)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_Prtb)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_AD)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_Prtb, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atm_TL)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atm_Prtb)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atm_AD)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Atm%Add_Extra_Layers = .FALSE.
  Atm_TL%Add_Extra_Layers = .FALSE.
  Atm_Prtb%Add_Extra_Layers = .FALSE.
  Atm_AD%Add_Extra_Layers = .FALSE.

  ! ============================================================================
  ! 3. **** ASSIGN INPUT DATA ****
  !
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()

  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! ============================================================================
  ! 4. **** SET PERTURBATIONS ****
  !
  Atm_Prtb = Atm
  Atm_TL = Atm
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Sfc_TL = Sfc
  CALL CRTM_Surface_Zero( Sfc_TL )

  Atm_TL(PROFILE_IDX)%Temperature(LAYER_IDX) = DELTA_T
  Atm_Prtb(PROFILE_IDX)%Temperature(LAYER_IDX) = &
    Atm(PROFILE_IDX)%Temperature(LAYER_IDX) + DELTA_T

  ! ============================================================================
  ! 5. **** CALL THE CRTM MODELS ****
  !
  Error_Status = CRTM_Tangent_Linear( Atm          , &
                                      Sfc          , &
                                      Atm_TL       , &
                                      Sfc_TL       , &
                                      Geometry     , &
                                      ChannelInfo  , &
                                      RTSolution   , &
                                      RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Tangent-Linear Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  IF ( run_fd ) THEN
    Error_Status = CRTM_Forward( Atm_Prtb     , &
                                 Sfc          , &
                                 Geometry     , &
                                 ChannelInfo  , &
                                 RTSolution_Prtb )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in perturbed CRTM Forward Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  ! ============================================================================
  ! 6. **** COMPARE FD WITH TL ****
  !
  exit_status = 0
  IF ( run_fd ) THEN
    WRITE(*,'(/5x,"FD vs TL per-channel differences:")')
    WRITE(*,'(5x,"chan  dR(FD)            dR(TL)            |FD-TL|          tol")')
    n_fail = 0
    max_diff = -1.0_fp
    channel_max = -1
    DO l = 1, n_Channels
      fd = RTSolution_Prtb(l,PROFILE_IDX)%Radiance - &
           RTSolution(l,PROFILE_IDX)%Radiance
      tl = RTSolution_TL(l,PROFILE_IDX)%Radiance
      diff = ABS(fd - tl)
      tol = ATOL + RTOL * MAX(ABS(fd), ABS(tl))
      WRITE(*,'(5x,i4,1x,4es16.6)') RTSolution(l,PROFILE_IDX)%Sensor_Channel, &
        fd, tl, diff, tol
      IF ( diff > tol ) THEN
        n_fail = n_fail + 1
        IF ( diff > max_diff ) THEN
          max_diff = diff
          channel_max = RTSolution(l,PROFILE_IDX)%Sensor_Channel
        END IF
      END IF
    END DO

    IF ( n_fail > 0 ) THEN
      WRITE(*,'(/5x,"FD vs TL check failed for ",i0," channels.")') n_fail
      WRITE(*,'(5x,"Largest |FD-TL| at channel ",i0,": ",es12.4)') channel_max, max_diff
      exit_status = 1
    ELSE
      WRITE(*,'(/5x,"FD vs TL check passed.")')
    END IF
  END IF

  ! ============================================================================
  ! 7. **** TL-AD DOT-PRODUCT CHECK ****
  !
  IF ( run_ad ) THEN
    Atm_AD = Atm
    Sfc_AD = Sfc
    CALL CRTM_Atmosphere_Zero( Atm_AD )
    CALL CRTM_Surface_Zero( Sfc_AD )
    CALL CRTM_RTSolution_Zero( RTSolution_AD )
    RTSolution_AD%Radiance = 1.0_fp

    Error_Status = CRTM_Adjoint( Atm          , &
                                 Sfc          , &
                                 RTSolution_AD, &
                                 Geometry     , &
                                 ChannelInfo  , &
                                 Atm_AD       , &
                                 Sfc_AD       , &
                                 RTSolution    )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Adjoint Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF

    lhs = 0.0_fp
    DO m = 1, N_PROFILES
      DO l = 1, n_Channels
        lhs = lhs + RTSolution_TL(l,m)%Radiance * RTSolution_AD(l,m)%Radiance
      END DO
    END DO

    rhs = 0.0_fp
    DO m = 1, N_PROFILES
      DO k = 1, N_LAYERS
        rhs = rhs + Atm_TL(m)%Temperature(k) * Atm_AD(m)%Temperature(k)
      END DO
    END DO

    rel = ABS(lhs - rhs) / MAX(ABS(lhs), ABS(rhs), AD_ATOL)
    WRITE(*,'(/5x,"TL-AD dot-product check:")')
    WRITE(*,'(5x,"<TL,AD_in> = ",es16.6)') lhs
    WRITE(*,'(5x,"<TL_in,AD_out> = ",es16.6)') rhs
    WRITE(*,'(5x,"rel diff = ",es12.4," (tol ",es12.4,")")') rel, AD_RTOL
    IF ( rel > AD_RTOL ) THEN
      WRITE(*,'(5x,"TL-AD check failed.")')
      exit_status = 1
    ELSE
      WRITE(*,'(5x,"TL-AD check passed.")')
    END IF
  END IF

  ! ============================================================================
  ! 7. **** CLEAN UP ****
  !
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    exit_status = 1
  END IF

  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_Prtb )
  CALL CRTM_Atmosphere_Destroy( Atm_AD )
  DEALLOCATE( RTSolution, RTSolution_TL, RTSolution_Prtb, RTSolution_AD, STAT=Allocate_Status )

  STOP exit_status

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_FD_consistency
