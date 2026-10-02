!
! test_ChannelSubset_OMP
!
! Thread-safety regression test for CRTM channel subsetting under OpenMP
! (follow-up to JCSDA/CRTMv3#164, exercised the way #111's test_OMP_Consistency
! exercises the full-channel path).
!
! For a given sensor it applies several CRTM_ChannelInfo_Subset patterns and,
! for each, runs CRTM_Forward and CRTM_K_Matrix at OMP_NUM_THREADS = 1 (the
! serial reference) then again at an increasing sweep of thread counts
! (2, 4, 8, ... up to the number available on the host), asserting that every
! result is BIT-IDENTICAL to the serial run.
!
! The subset patterns are chosen to stress the channel-chunking / inactive-
! channel bookkeeping in CRTM_Forward / _K_Matrix:
!   * "sparse"  - channels spread roughly uniformly across the full range, so
!                 every channel-thread chunk holds a mix of active/inactive
!                 channels (the prefix-sum offsets must all be non-trivial);
!   * "front"   - only the first handful of channels, so the leading chunk is
!                 fully active and every trailing chunk is entirely inactive
!                 (a thread that processes zero channels);
!   * "split"   - a few channels near the start and a few near the end, so the
!                 leading chunk is partly active, the middle chunks are empty,
!                 and the trailing chunk picks up again -- this is the case
!                 that historically produced the off-by-one in the "ln" output
!                 index when split across threads.
!
! No reference data files are needed -- the invariant is parallel == serial.
!
! No-op (treated as PASS) when CRTM is built without OpenMP, or when only a
! single hardware thread is available.
!

PROGRAM test_ChannelSubset_OMP

  USE CRTM_Module
#ifdef _OPENMP
  USE OMP_LIB
#endif
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_ChannelSubset_OMP'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 1
  INTEGER,  PARAMETER :: N_SENSORS   = 1

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  INTEGER,  PARAMETER :: N_SUBSET_MODES = 3
  CHARACTER(8), PARAMETER :: MODE_NAME(N_SUBSET_MODES) = (/ 'sparse  ', &
                                                           'front   ', &
                                                           'split   ' /)

  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status, Allocate_Status
  INTEGER :: n_full, n_sub
  INTEGER :: imode, n_mismatch_total

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)

  INTEGER, ALLOCATABLE :: subset(:)        ! channel numbers requested this mode
  INTEGER, ALLOCATABLE :: full_channels(:) ! all of this sensor's channel numbers

#ifdef _OPENMP
  INTEGER, PARAMETER :: MAX_SWEEP = 8
  INTEGER :: sweep(MAX_SWEEP), n_sweep, n_threads_avail, cand
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
    'OpenMP / serial consistency test for CRTM channel subsetting.', &
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
  n_full = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  WRITE(*,'(5x,a,i0,a,i0)') 'Full channel count: ', n_full, '   Profiles: ', N_PROFILES
  ! Capture this sensor's actual channel numbers (NOT assumed to be 1..n_full --
  ! "subset" SpcCoeff files such as cris399_npp carry a sparse channel list).
  ALLOCATE( full_channels(n_full), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating channel-number array', FAILURE )
    STOP 1
  END IF
  full_channels = CRTM_ChannelInfo_Channels( ChannelInfo(1) )

  ! --- Static inputs (shared across all subset modes) ---
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating Atmosphere structures', FAILURE )
    STOP 1
  END IF
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  n_mismatch_total = 0

  Subset_Mode_Loop: DO imode = 1, N_SUBSET_MODES

    ! Restore the full channel set, then apply this mode's subset.
    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Reset = .TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error resetting ChannelInfo subset', FAILURE )
      STOP 1
    END IF

    CALL Build_Subset( imode, full_channels, subset )
    Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), Channel_Subset = subset )
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE(Message,'("Error applying """,a,""" channel subset")') TRIM(MODE_NAME(imode))
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF

    n_sub = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
    ! CRTM_ChannelInfo_Subset must turn on exactly the channels requested.
    IF ( n_sub /= SIZE(subset) ) THEN
      WRITE(Message,'("Subset """,a,""": requested ",i0," channels but ",i0," are active")') &
        TRIM(MODE_NAME(imode)), SIZE(subset), n_sub
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF
    IF ( ANY( CRTM_ChannelInfo_Channels(ChannelInfo(1)) /= subset ) ) THEN
      WRITE(Message,'("Subset """,a,""": active channel list does not match the request")') &
        TRIM(MODE_NAME(imode))
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE )
      STOP 1
    END IF

    WRITE(*,'(/5x,"--- Subset mode """,a,""": ",i0," of ",i0," channels ---")') &
      TRIM(MODE_NAME(imode)), n_sub, n_full

#ifdef _OPENMP
    CALL Run_Sweep( n_sub, n_mismatch_total )
#else
    WRITE(*,'(5x,a)') '(no OpenMP - subset applied but no thread sweep performed)'
#endif

    DEALLOCATE( subset )

  END DO Subset_Mode_Loop

#ifdef _OPENMP
  IF ( n_mismatch_total > 0 ) THEN
    WRITE(*,'(/5x,"FAIL: ",i0," parallel result(s) differed from the serial run.")') n_mismatch_total
    STOP 1
  END IF
  WRITE(*,'(/5x,"PASS: subset Forward & K-Matrix are thread-count invariant for all modes.")')
#endif

  ! --- Cleanup ---
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error destroying CRTM', FAILURE )
    STOP 1
  END IF
  CALL CRTM_Atmosphere_Destroy( Atm )
  DEALLOCATE( full_channels, STAT=Allocate_Status )

CONTAINS

  ! Build the requested channel-number list for a given subset mode. Picks
  ! positions within the sensor's channel list (ch = all of this sensor's
  ! channel numbers, in storage order) so the result is valid for both
  ! contiguous-numbered sensors and sparse "subset" SpcCoeff files. The picked
  ! positions are what drives the channel-thread chunking in CRTM_Forward, so
  ! the modes are described in terms of positions, not channel numbers.
  ! Allocates "list" to the exact size needed.
  SUBROUTINE Build_Subset( mode, ch, list )
    INTEGER,              INTENT(IN)  :: mode
    INTEGER,              INTENT(IN)  :: ch(:)
    INTEGER, ALLOCATABLE, INTENT(OUT) :: list(:)
    INTEGER :: nfull, step, n, i, nf, nt

    nfull = SIZE(ch)

    SELECT CASE ( mode )

    CASE ( 1 )  ! "sparse": ~13 positions spread across the whole range
      step = MAX( 1, nfull / 13 )
      n    = ( nfull - 1 ) / step + 1
      ALLOCATE( list(n) )
      DO i = 1, n
        list(i) = ch( 1 + (i-1)*step )
      END DO

    CASE ( 2 )  ! "front": just the first handful of positions
      n = MIN( nfull, 7 )
      ALLOCATE( list(n) )
      DO i = 1, n
        list(i) = ch(i)
      END DO

    CASE DEFAULT  ! "split": a few at the front and a few at the tail
      nf = MIN( nfull, 4 )                 ! positions at the front
      nt = MIN( nfull - nf, 3 )            ! positions at the tail (no overlap)
      ALLOCATE( list(nf + nt) )
      DO i = 1, nf
        list(i) = ch(i)
      END DO
      DO i = 1, nt
        list(nf + i) = ch( nfull - nt + i )
      END DO

    END SELECT
  END SUBROUTINE Build_Subset

#ifdef _OPENMP
  ! Run the serial reference + parallel thread sweep for the currently-active
  ! channel subset (n_sub channels). Increments n_mismatch by the number of
  ! parallel results that differed from the serial reference.
  SUBROUTINE Run_Sweep( n_sub, n_mismatch )
    INTEGER, INTENT(IN)    :: n_sub
    INTEGER, INTENT(INOUT) :: n_mismatch
    INTEGER :: isweep, nthr
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_ref(:,:)
    TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:), Atmosphere_K_ref(:,:)
    TYPE(CRTM_Surface_type)   , ALLOCATABLE :: Surface_K(:,:)   , Surface_K_ref(:,:)
    TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:), RTSolution_K_ref(:,:)

    ALLOCATE( RTSolution      (n_sub, N_PROFILES), &
              RTSolution_ref  (n_sub, N_PROFILES), &
              Atmosphere_K    (n_sub, N_PROFILES), &
              Atmosphere_K_ref(n_sub, N_PROFILES), &
              Surface_K       (n_sub, N_PROFILES), &
              Surface_K_ref   (n_sub, N_PROFILES), &
              RTSolution_K    (n_sub, N_PROFILES), &
              RTSolution_K_ref(n_sub, N_PROFILES), &
              STAT = Allocate_Status )
    IF ( Allocate_Status /= 0 ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error allocating result arrays', FAILURE )
      STOP 1
    END IF
    CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_K)) ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error allocating Atmosphere_K structures', FAILURE )
      STOP 1
    END IF

    ! --- reference run @ 1 thread ---
    CALL OMP_SET_NUM_THREADS(1)
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution_ref )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Serial CRTM_Forward failed', FAILURE )
      STOP 1
    END IF
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

    ! --- parallel sweep ---
    DO isweep = 2, n_sweep
      nthr = sweep(isweep)
      CALL OMP_SET_NUM_THREADS(nthr)

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

      WRITE(*,'(7x,"OMP_NUM_THREADS=",i0,": Forward & K-Matrix bit-identical to serial.")') nthr
    END DO

    CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
    CALL CRTM_Atmosphere_Destroy( Atmosphere_K_ref )
    DEALLOCATE( RTSolution, RTSolution_ref, Atmosphere_K, Atmosphere_K_ref, &
                Surface_K, Surface_K_ref, RTSolution_K, RTSolution_K_ref, &
                STAT = Allocate_Status )
  END SUBROUTINE Run_Sweep
#endif

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

END PROGRAM test_ChannelSubset_OMP
