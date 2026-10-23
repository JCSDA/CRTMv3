!
! test_OMP_Thread_Policy
!
! Guards two properties of how CRTM decides to use OpenMP threads. Both are
! about the single-profile call, which is what GSI issues (crtm_interface.f90
! passes an atmosphere array of dimension(1)) and which is the case that used
! to be handled worst.
!
!   A. CRTM must not leave the OpenMP runtime reconfigured behind the caller's
!      back. CRTM raises max-active-levels to run its nested channel loop; that
!      setting is global and outlives the call, so a host doing its own
!      threading would find its nesting policy silently replaced. This check is
!      exact and carries no timing.
!
!   B. Threading must never be dramatically slower than not threading. Channel
!      threading gives every thread its own AtmOptics/SfcOptics/RTV/scatter
!      scratch, sized by layers and stream count rather than by the channels the
!      thread receives, so splitting a small sensor across many threads once
!      cost far more than it saved: one ATMS profile on 16 threads measured
!      0.03x, roughly 30 times slower than one thread.
!
! Part B is a timing check and so is bounded loosely on purpose. It asserts only
! that the threaded path is not more than SLOWDOWN_LIMIT times slower than the
! serial path, against a regression that was ~30x. That leaves room for a busy
! machine while still catching the defect. Registered RUN_SERIAL for the same
! reason. If this test fails on timing alone, suspect the channel-thread gate
! (MIN_CHANNELS_PER_CHANNEL_THREAD in CRTM_Parameters) before suspecting the host.
!

PROGRAM test_OMP_Thread_Policy

  USE CRTM_Module
  USE OMP_LIB
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_OMP_Thread_Policy'
  CHARACTER(*), PARAMETER :: SENSOR_ID    = 'atms_n21'
  CHARACTER(*), PARAMETER :: COEFF_PATH   = './testinput/'

  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  INTEGER,  PARAMETER :: N_PROFILES  = 1     ! the GSI case
  INTEGER,  PARAMETER :: N_TRIALS    = 3
  INTEGER,  PARAMETER :: N_CALLS     = 400   ! ~0.1 s per trial; above clock noise
  REAL(fp), PARAMETER :: SLOWDOWN_LIMIT = 3.0_fp

  INTEGER  :: err, n_Channels, k, m, itrial, icall
  INTEGER  :: lev_before, lev_after_fwd, lev_after_k, max_threads
  INTEGER  :: c0, c1, crate
  REAL(fp) :: t_serial, t_parallel, t_trial, ratio, f, p_top, p_sfc
  INTEGER  :: n_fail

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm(:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc(:)
  TYPE(CRTM_Geometry_type),   ALLOCATABLE :: Geo(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  n_fail = 0

  WRITE(*,'(/5x,a)') '**********************************************************'
  WRITE(*,'(5x,a)')  '                 test_OMP_Thread_Policy'
  WRITE(*,'(5x,a)')  '**********************************************************'

  err = CRTM_Init( (/ SENSOR_ID /), ChannelInfo, &
                   File_Path         = COEFF_PATH, &
                   Load_CloudCoeff   = .FALSE., &
                   Load_AerosolCoeff = .FALSE., &
                   Quiet             = .TRUE. )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( Atm(N_PROFILES), Sfc(N_PROFILES), Geo(N_PROFILES), &
            RTSolution(n_Channels, N_PROFILES),   &
            RTSolution_K(n_Channels, N_PROFILES), &
            Atm_K(n_Channels, N_PROFILES), Sfc_K(n_Channels, N_PROFILES) )
  CALL CRTM_Atmosphere_Create( Atm,   N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_Atmosphere_Create( Atm_K, N_LAYERS, N_ABSORBERS, 0, 0 )
  CALL CRTM_RTSolution_Create( RTSolution,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K, N_LAYERS )

  ! A plausible mid-latitude sounding. Absolute realism does not matter here,
  ! only that the profile is valid and exercises the full layer count.
  p_top = 0.1_fp; p_sfc = 1013.0_fp
  DO m = 1, N_PROFILES
    Atm(m)%Climatology    = US_STANDARD_ATMOSPHERE
    Atm(m)%Absorber_Id    = (/ H2O_ID, O3_ID /)
    Atm(m)%Absorber_Units = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    Atm(m)%Level_Pressure(0) = p_top
    DO k = 1, N_LAYERS
      f = REAL(k,fp) / REAL(N_LAYERS,fp)
      Atm(m)%Level_Pressure(k) = p_top * EXP( f * LOG(p_sfc/p_top) )
      Atm(m)%Pressure(k)    = 0.5_fp*(Atm(m)%Level_Pressure(k-1)+Atm(m)%Level_Pressure(k))
      Atm(m)%Temperature(k) = 215.0_fp + 75.0_fp*f
      Atm(m)%Absorber(k,1)  = MAX( 1.0e-2_fp, 12.0_fp * f**3 )
      Atm(m)%Absorber(k,2)  = MAX( 1.0e-2_fp, 8.0_fp * (1.0_fp - f)**2 )
    END DO
    Sfc(m)%Water_Coverage    = 1.0_fp
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 6.25_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geo(m), Sensor_Zenith_Angle = 30.0_fp, &
                                         Sensor_Scan_Angle   = 26.37_fp )
  END DO

  max_threads = OMP_GET_MAX_THREADS()
  WRITE(*,'(/5x,a,a)')  'Sensor                : ', SENSOR_ID
  WRITE(*,'(5x,a,i0)')  'Channels              : ', n_Channels
  WRITE(*,'(5x,a,i0)')  'Profiles per call     : ', N_PROFILES
  WRITE(*,'(5x,a,i0)')  'Threads available     : ', max_threads

  ! ------------------------------------------------------------------
  ! A. The caller's nesting policy must survive a CRTM call
  ! ------------------------------------------------------------------
  CALL OMP_SET_MAX_ACTIVE_LEVELS(1)
  lev_before = OMP_GET_MAX_ACTIVE_LEVELS()

  err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTSolution )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
  END IF
  lev_after_fwd = OMP_GET_MAX_ACTIVE_LEVELS()

  CALL OMP_SET_MAX_ACTIVE_LEVELS(1)
  CALL CRTM_Atmosphere_Zero( Atm_K ); CALL CRTM_Surface_Zero( Sfc_K )
  CALL CRTM_RTSolution_Zero( RTSolution_K )
  RTSolution_K%Brightness_Temperature = ONE
  err = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geo, ChannelInfo, &
                       Atm_K, Sfc_K, RTSolution )
  IF ( err /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_K_Matrix failed', FAILURE ); STOP 1
  END IF
  lev_after_k = OMP_GET_MAX_ACTIVE_LEVELS()

  WRITE(*,'(/5x,a)')    '--- A. caller OpenMP nesting policy ---'
  WRITE(*,'(5x,a,i0)')  'set by caller         : ', lev_before
  WRITE(*,'(5x,a,i0)')  'after CRTM_Forward    : ', lev_after_fwd
  WRITE(*,'(5x,a,i0)')  'after CRTM_K_Matrix   : ', lev_after_k
  IF ( lev_after_fwd /= lev_before .OR. lev_after_k /= lev_before ) THEN
    n_fail = n_fail + 1
    CALL Display_Message( PROGRAM_NAME, &
         'CRTM changed the caller max-active-levels and did not restore it', FAILURE )
  END IF

  ! ------------------------------------------------------------------
  ! B. Threading must not be pathologically slower than serial
  ! ------------------------------------------------------------------
  IF ( max_threads > 1 ) THEN

    CALL OMP_SET_NUM_THREADS(1)
    err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTSolution )   ! warm up
    t_serial = HUGE(t_serial)
    DO itrial = 1, N_TRIALS
      CALL SYSTEM_CLOCK(c0, crate)
      DO icall = 1, N_CALLS
        err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTSolution )
        IF ( err /= SUCCESS ) THEN
          CALL Display_Message( PROGRAM_NAME, 'serial CRTM_Forward failed', FAILURE ); STOP 1
        END IF
      END DO
      CALL SYSTEM_CLOCK(c1)
      t_trial = REAL(c1-c0,fp)/REAL(crate,fp)
      t_serial = MIN(t_serial, t_trial)
    END DO

    CALL OMP_SET_NUM_THREADS(max_threads)
    err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTSolution )   ! warm up
    t_parallel = HUGE(t_parallel)
    DO itrial = 1, N_TRIALS
      CALL SYSTEM_CLOCK(c0, crate)
      DO icall = 1, N_CALLS
        err = CRTM_Forward( Atm, Sfc, Geo, ChannelInfo, RTSolution )
        IF ( err /= SUCCESS ) THEN
          CALL Display_Message( PROGRAM_NAME, 'threaded CRTM_Forward failed', FAILURE ); STOP 1
        END IF
      END DO
      CALL SYSTEM_CLOCK(c1)
      t_trial = REAL(c1-c0,fp)/REAL(crate,fp)
      t_parallel = MIN(t_parallel, t_trial)
    END DO

    ratio = t_parallel / MAX(t_serial, TINY(t_serial))
    WRITE(*,'(/5x,a)')      '--- B. threaded vs serial, small sensor ---'
    WRITE(*,'(5x,a,f10.4,a)') 'serial   wall (best) : ', t_serial,   ' s'
    WRITE(*,'(5x,a,f10.4,a)') 'threaded wall (best) : ', t_parallel, ' s'
    WRITE(*,'(5x,a,f10.3,a)') 'slowdown             : ', ratio, ' x'
    WRITE(*,'(5x,a,f10.3,a)') 'limit                : ', SLOWDOWN_LIMIT, ' x'
    IF ( ratio > SLOWDOWN_LIMIT ) THEN
      n_fail = n_fail + 1
      CALL Display_Message( PROGRAM_NAME, &
           'threading a small sensor is far slower than not threading it', FAILURE )
    END IF
  ELSE
    WRITE(*,'(/5x,a)') '--- B. skipped, only one thread available ---'
  END IF

  err = CRTM_Destroy( ChannelInfo )

  WRITE(*,'(/5x,a)') '**********************************************************'
  IF ( n_fail == 0 ) THEN
    WRITE(*,'(5x,a/)') 'PASS: OpenMP thread policy is sane.'
    STOP 0
  ELSE
    WRITE(*,'(5x,a,i0,a/)') 'FAIL: ', n_fail, ' OpenMP thread-policy check(s) failed.'
    STOP 1
  END IF

END PROGRAM test_OMP_Thread_Policy
