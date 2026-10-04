!
! test_OMP_Consistency
!
! Thread-safety regression test (JCSDA/CRTMv3#111). For a given sensor it runs
! CRTM_Forward, CRTM_Tangent_Linear and CRTM_K_Matrix on the same input at OMP_NUM_THREADS = 1
! (the serial reference) and then again at an increasing sweep of thread counts
! (2, 4, 8, ... up to the number available on the host), asserting that every
! result is BIT-IDENTICAL to the serial run.
!
! Channel-thread parallelism in CRTM_Forward / _Tangent_Linear / _K_Matrix must not change any
! per-channel value, so the parallel and serial outputs are required to be
! exactly equal -- not merely "close". This is the kind of check that surfaces
! the races fixed for #111 (unindexed RTV broadcasts, off-by-one / overshoot in
! the channel-chunk math -> out-of-bounds writes, per-channel NLTE/Zeeman
! predictor contamination -> wrong Jacobians, shared Error_Status writes).
!
! No reference data files are needed -- the invariant is parallel == serial.
!
! No-op (treated as PASS) when CRTM is built without OpenMP, or when only a
! single hardware thread is available.
!

PROGRAM test_OMP_Consistency

  USE CRTM_Module
#ifdef _OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_OMP_Consistency'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 1
  INTEGER,  PARAMETER :: N_SENSORS   = 1

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status, Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: l, isweep, nthr, n_mismatch

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  ! Forward
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_ref(:,:)
  ! Tangent-Linear
  TYPE(CRTM_Atmosphere_type)              :: Atmosphere_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Surface_TL(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TL(:,:), RTSolution_TL_ref(:,:)
  ! K-Matrix
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:), Atmosphere_K_ref(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: Surface_K(:,:)   , Surface_K_ref(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:), RTSolution_K_ref(:,:)

#ifdef _OPENMP
  INTEGER, PARAMETER :: MAX_SWEEP = 8
  INTEGER :: sweep(MAX_SWEEP), n_sweep, n_threads_avail, cand
  INTEGER :: i
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
    'OpenMP / serial consistency test for CRTM_Forward, CRTM_Tangent_Linear and CRTM_K_Matrix.', &
    'CRTM Version: '//TRIM(Version) )
  WRITE( *,'(/5x,"Sensor: ",a)' ) TRIM(Sensor_Id)

#ifndef _OPENMP
  WRITE(*,'(/5x,a)') 'CRTM was built without OpenMP (_OPENMP undefined).'
  WRITE(*,'(5x,a)')  'Consistency test is a no-op in this configuration (PASS).'
  STOP 0
#else
  ! How much parallelism does this host actually offer?  Query BEFORE CRTM_Init,
  ! which coerces the thread count to 1 when OMP_NUM_THREADS is unset/empty.
  n_threads_avail = OMP_GET_MAX_THREADS()
  IF ( n_threads_avail <= 1 ) THEN
    WRITE(*,'(/5x,a,i0,a)') 'OMP_GET_MAX_THREADS() = ', n_threads_avail, &
      ' -- no parallelism available, nothing to compare. Skipping (PASS).'
    STOP 0
  END IF

  ! Build the thread-count sweep: 1, then 2,4,8,... up to n_threads_avail,
  ! and n_threads_avail itself (deduplicated).
  sweep(1) = 1
  n_sweep  = 1
  cand     = 2
  DO WHILE ( cand < n_threads_avail .AND. n_sweep < MAX_SWEEP-1 )
    n_sweep        = n_sweep + 1
    sweep(n_sweep) = cand
    cand           = cand * 2
  END DO
  IF ( sweep(n_sweep) /= n_threads_avail ) THEN
    n_sweep        = n_sweep + 1
    sweep(n_sweep) = n_threads_avail
  END IF
  WRITE(*,'(/5x,a,8(i0,1x))') 'Thread-count sweep: ', sweep(1:n_sweep)
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
  WRITE(*,'(5x,a,i0,a,i0,a)') 'Channels: ', n_Channels, '   Profiles: ', N_PROFILES, ''

  ! --- Allocate ---
  ALLOCATE( RTSolution      (n_Channels, N_PROFILES), &
            RTSolution_ref  (n_Channels, N_PROFILES), &
            RTSolution_TL   (n_Channels, N_PROFILES), &
            RTSolution_TL_ref(n_Channels, N_PROFILES), &
            Atmosphere_K    (n_Channels, N_PROFILES), &
            Atmosphere_K_ref(n_Channels, N_PROFILES), &
            Surface_K       (n_Channels, N_PROFILES), &
            Surface_K_ref   (n_Channels, N_PROFILES), &
            RTSolution_K    (n_Channels, N_PROFILES), &
            RTSolution_K_ref(n_Channels, N_PROFILES), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating result arrays', FAILURE )
    STOP 1
  END IF
  CALL CRTM_Atmosphere_Create( Atm,           N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_K)) .OR. &
       ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_TL)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating Atmosphere structures', FAILURE )
    STOP 1
  END IF

  ! --- Populate inputs ---
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! Tangent-linear perturbation: zero, then +0.5 K on temperature (matches the
  ! standard tangent_linear regression setup). The actual perturbation values
  ! are immaterial here -- the invariant under test is parallel == serial.
  Atmosphere_TL = Atm
  CALL CRTM_Atmosphere_Zero( Atmosphere_TL )
  DO l = 1, N_PROFILES
    Atmosphere_TL(l)%Temperature = 0.5_fp
  END DO
  Surface_TL = Sfc
  CALL CRTM_Surface_Zero( Surface_TL )

#ifdef _OPENMP
  ! ===================  reference run @ 1 thread  ===================
  CALL OMP_SET_NUM_THREADS(1)

  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_ref )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Serial CRTM_Forward failed', FAILURE )
    STOP 1
  END IF

  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atmosphere_TL, Surface_TL, &
                                      Geometry, ChannelInfo, RTSolution, RTSolution_TL )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Serial CRTM_Tangent_Linear failed', FAILURE )
    STOP 1
  END IF
  RTSolution_TL_ref = RTSolution_TL

  CALL Init_K_Inputs( Atmosphere_K, Surface_K, RTSolution_K )
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atmosphere_K, Surface_K, RTSolution )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Serial CRTM_K_Matrix failed', FAILURE )
    STOP 1
  END IF
  Atmosphere_K_ref = Atmosphere_K
  Surface_K_ref    = Surface_K
  RTSolution_K_ref = RTSolution_K

  ! ===================  parallel sweep  ===================
  n_mismatch = 0
  DO isweep = 2, n_sweep
    nthr = sweep(isweep)
    CALL OMP_SET_NUM_THREADS(nthr)

    ! --- Forward ---
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE(Message,'("CRTM_Forward failed at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF
    IF ( .NOT. ALL(RTSolution == RTSolution_ref) ) THEN
      WRITE(Message,'("Forward result differs from the 1-thread run at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      n_mismatch = n_mismatch + 1
    END IF

    ! --- Tangent-Linear ---
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atmosphere_TL, Surface_TL, &
                                        Geometry, ChannelInfo, RTSolution, RTSolution_TL )
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE(Message,'("CRTM_Tangent_Linear failed at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF
    IF ( .NOT. ALL(RTSolution_TL == RTSolution_TL_ref) ) THEN
      WRITE(Message,'("Tangent-Linear RTSolution_TL differs from the 1-thread run at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      n_mismatch = n_mismatch + 1
    END IF

    ! --- K-Matrix ---
    CALL Init_K_Inputs( Atmosphere_K, Surface_K, RTSolution_K )
    Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                  Atmosphere_K, Surface_K, RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE(Message,'("CRTM_K_Matrix failed at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF
    IF ( .NOT. ALL(RTSolution_K == RTSolution_K_ref) ) THEN
      WRITE(Message,'("K-Matrix RTSolution_K differs from the 1-thread run at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      n_mismatch = n_mismatch + 1
    END IF
    IF ( .NOT. ALL(Atmosphere_K == Atmosphere_K_ref) ) THEN
      WRITE(Message,'("K-Matrix Atmosphere_K differs from the 1-thread run at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      n_mismatch = n_mismatch + 1
    END IF
    IF ( .NOT. ALL(Surface_K == Surface_K_ref) ) THEN
      WRITE(Message,'("K-Matrix Surface_K differs from the 1-thread run at OMP_NUM_THREADS=",i0)') nthr
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      n_mismatch = n_mismatch + 1
    END IF

    WRITE(*,'(5x,"OMP_NUM_THREADS=",i0,": Forward, Tangent-Linear & K-Matrix bit-identical to serial.")') nthr
  END DO

  IF ( n_mismatch > 0 ) THEN
    WRITE(*,'(/5x,"FAIL: ",i0," parallel result(s) differed from the serial run.")') n_mismatch
    STOP 1
  END IF
  WRITE(*,'(/5x,"PASS: Forward, Tangent-Linear & K-Matrix are thread-count invariant.")')
#endif

  ! --- Cleanup ---
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error destroying CRTM', FAILURE )
    STOP 1
  END IF
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K_ref )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_TL )
  DEALLOCATE( RTSolution, RTSolution_ref, RTSolution_TL, RTSolution_TL_ref, &
              Atmosphere_K, Atmosphere_K_ref, &
              Surface_K, Surface_K_ref, RTSolution_K, RTSolution_K_ref, &
              STAT=Allocate_Status )

CONTAINS

  ! Reset the K-matrix inputs the way the regression k_matrix tests do:
  ! zero the Jacobian accumulators, and unit-perturb the brightness temperature
  ! (radiance for visible) so the call produces the BT/radiance Jacobian.
  SUBROUTINE Init_K_Inputs( Atm_K, Sfc_K, RTSol_K )
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: Atm_K(:,:)
    TYPE(CRTM_Surface_type)   , INTENT(IN OUT) :: Sfc_K(:,:)
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTSol_K(:,:)
    CALL CRTM_Atmosphere_Zero( Atm_K )
    CALL CRTM_Surface_Zero( Sfc_K )
    IF ( ChannelInfo(1)%Sensor_Type == INFRARED_SENSOR .OR. &
         ChannelInfo(1)%Sensor_Type == MICROWAVE_SENSOR ) THEN
      RTSol_K(:,:)%Radiance               = ZERO
      RTSol_K(:,:)%Brightness_Temperature = ONE
    ELSE
      RTSol_K(:,:)%Radiance               = ONE
      RTSol_K(:,:)%Brightness_Temperature = ZERO
    END IF
  END SUBROUTINE Init_K_Inputs

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_OMP_Consistency
