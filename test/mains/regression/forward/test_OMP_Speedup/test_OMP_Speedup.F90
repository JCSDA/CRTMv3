!
! test_OMP_Speedup
!
! Times CRTM_Forward in serial (1 thread) and parallel (max threads)
! configurations on the same input, then asserts that the parallel run
! achieves a meaningful speedup. Demonstrates that the OpenMP build is
! both functional and effective end-to-end.
!
! Behavior when OpenMP is not enabled at compile time, or when only a
! single hardware thread is available, the test prints an explanatory
! message and exits successfully (treated as a no-op).
!

PROGRAM test_OMP_Speedup

  USE CRTM_Module
#ifdef _OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_OMP_Speedup'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 1
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  INTEGER,  PARAMETER :: N_REPEATS   = 5
  ! Each phase is timed N_TRIALS times and the fastest trial is kept. A single
  ! measurement is not robust: on this class of host the observed nvfortran
  ! speedup ranged 1.11x to 2.25x across repeat runs of the identical binary on
  ! an idle machine, straddling MIN_SPEEDUP. Taking the best trial suppresses
  ! transient interference (scheduler, thermal, page-cache) without inflating
  ! the result, since the fastest run is the one least perturbed.
  INTEGER,  PARAMETER :: N_TRIALS    = 3

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! Pass threshold. Conservative on purpose: portable CI machines and laptops
  ! vary widely in OMP scaling. We just want to confirm meaningful speedup.
  REAL(fp), PARAMETER :: MIN_SPEEDUP = 1.20_fp

  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status, Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: irep

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

#ifdef _OPENMP
  REAL(fp) :: t0, t_trial, t_serial, t_parallel, speedup
  INTEGER  :: max_threads, itrial
#endif

  ! --- Argument parsing ---
  IF ( COMMAND_ARGUMENT_COUNT() /= 1 ) THEN
    WRITE(*,*) PROGRAM_NAME//': ERROR, requires one argument: <Sensor_Id>'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, Sensor_Id)
  Sensor_Id = ADJUSTL(Sensor_Id)

  CALL CRTM_Version(Version)
  CALL Program_Message( PROGRAM_NAME, &
    'OpenMP forward-model speedup test.', &
    'CRTM Version: '//TRIM(Version) )
  WRITE( *,'(/5x,"Sensor: ",a)' ) TRIM(Sensor_Id)

  ! --- OpenMP availability gates ---
#ifndef _OPENMP
  WRITE(*,'(/5x,a)') 'CRTM was built without OpenMP (_OPENMP undefined).'
  WRITE(*,'(5x,a)')  'Speedup test is a no-op in this configuration.'
  STOP 0
#else
  max_threads = OMP_GET_MAX_THREADS()
  IF ( max_threads <= 1 ) THEN
    WRITE(*,'(/5x,a,i0,a)') 'OMP_GET_MAX_THREADS() returned ', max_threads, &
      ' — cannot demonstrate speedup. Skipping (PASS).'
    STOP 0
  END IF
#endif

  ! --- Initialize CRTM ---
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            File_Path = COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE )
    STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ! --- Allocate ---
  ALLOCATE( RTSolution(n_Channels, N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating RTSolution', FAILURE )
    STOP 1
  END IF
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating Atm', FAILURE )
    STOP 1
  END IF

  ! --- Populate inputs ---
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! --- Warmup call (untimed) to amortize first-call setup costs (page faults,
  !     coefficient lookup tables warming, etc.) so they don't show up as
  !     fake "serial overhead" on the first timed run.
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Warmup CRTM_Forward failed', FAILURE )
    STOP 1
  END IF

#ifdef _OPENMP
  ! --- Serial timing (1 thread), best of N_TRIALS ---
  CALL OMP_SET_NUM_THREADS(1)
  t_serial = HUGE(t_serial)
  DO itrial = 1, N_TRIALS
    t0 = OMP_GET_WTIME()
    DO irep = 1, N_REPEATS
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'Serial CRTM_Forward failed', FAILURE )
        STOP 1
      END IF
    END DO
    t_trial  = OMP_GET_WTIME() - t0
    t_serial = MIN( t_serial, t_trial )
  END DO

  ! --- Parallel timing (max threads), best of N_TRIALS ---
  CALL OMP_SET_NUM_THREADS(max_threads)
  t_parallel = HUGE(t_parallel)
  DO itrial = 1, N_TRIALS
    t0 = OMP_GET_WTIME()
    DO irep = 1, N_REPEATS
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'Parallel CRTM_Forward failed', FAILURE )
        STOP 1
      END IF
    END DO
    t_trial    = OMP_GET_WTIME() - t0
    t_parallel = MIN( t_parallel, t_trial )
  END DO

  speedup = t_serial / t_parallel

  WRITE(*,'(/5x,a)') '======================================================'
  WRITE(*,'(5x,a,a)')      'OpenMP speedup report — ', TRIM(Sensor_Id)
  WRITE(*,'(5x,a)') '======================================================'
  WRITE(*,'(5x,a,i0)')     'Threads (parallel run):   ', max_threads
  WRITE(*,'(5x,a,i0)')     'Profiles per call:        ', N_PROFILES
  WRITE(*,'(5x,a,i0)')     'Channels per call:        ', n_Channels
  WRITE(*,'(5x,a,i0)')     'Forward calls per phase:  ', N_REPEATS
  WRITE(*,'(5x,a,i0)')     'Timed trials per phase:   ', N_TRIALS
  WRITE(*,'(5x,a,f10.4,a)') 'Serial wall time (best):  ', t_serial,    ' s'
  WRITE(*,'(5x,a,f10.4,a)') 'Parallel wall time (best):', t_parallel,  ' s'
  WRITE(*,'(5x,a,f10.3,a)') 'Speedup (serial/parallel):', speedup,     ' x'
  WRITE(*,'(5x,a,f10.3,a)') 'Pass threshold:           ', MIN_SPEEDUP, ' x'
  WRITE(*,'(5x,a)') '======================================================'

  IF ( speedup < MIN_SPEEDUP ) THEN
    WRITE(*,'(/5x,a)') 'FAIL: parallel speedup below threshold.'
    STOP 1
  END IF
  WRITE(*,'(/5x,a)') 'PASS: OpenMP speedup demonstrated.'
#endif

  ! --- Cleanup ---
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error destroying CRTM', FAILURE )
    STOP 1
  END IF
  CALL CRTM_Atmosphere_Destroy(Atm)
  DEALLOCATE(RTSolution, STAT=Allocate_Status)

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_OMP_Speedup
