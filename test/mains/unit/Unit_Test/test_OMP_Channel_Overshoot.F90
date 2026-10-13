!
! test_OMP_Channel_Overshoot
!
! Regression test for JCSDA/CRTMv3#383: the channel-thread loop in the
! Forward, Tangent-Linear and K-Matrix modules split a sensor's channels into
! chunks of CEILING(n_Channels/n_channel_threads) but did not clamp the end of
! the non-last chunks. Whenever CEILING(n/t)*(t-1) >= n, the threads past the
! last channel read ChannelInfo%Process_Channel and Channel_Index beyond their
! bounds, and ODPS or CRTM_Compute_nStreams could segfault.
!
! A 10-channel sensor with 2 profiles at 16 threads gives 8 channel threads
! with chunks of 2, so threads 6 to 8 start past channel 10. 19 threads gives
! 9 channel threads. The test runs FWD, TL and K at 1 thread as the reference,
! then at 16 and 19 threads, and requires bit-identical results. The thread
! count is set here with OMP_SET_NUM_THREADS, so the test does not depend on
! the OMP_NUM_THREADS value ctest passes in.
!
PROGRAM test_OMP_Channel_Overshoot

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  USE CRTM_Module
  USE UnitTest_Define, ONLY: UnitTest_type
!$ USE OMP_LIB
  IMPLICIT NONE
  ! ============================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_OMP_Channel_Overshoot'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Profile dimensions (matching the Load_*_Data include files); clear sky
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0

  ! 10-channel sensor: CEILING(10/t)*(t-1) >= 10 for t = 8 and 9 channel threads
  CHARACTER(*), PARAMETER :: SENSOR_ID(1) = (/ 'abi_g18' /)

  ! Thread counts compared against the single-thread reference
  INTEGER, PARAMETER :: N_THREAD_CASES = 2
  INTEGER, PARAMETER :: THREAD_CASES(N_THREAD_CASES) = (/ 16, 19 /)

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  INTEGER :: Error_Status
  INTEGER :: n_Channels, i, m

  TYPE(UnitTest_type) :: test

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)

  ! Single-thread reference outputs
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_ref(:,:), rtsTL_ref(:,:), rtsKout_ref(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atmK_ref(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfcK_ref(:,:)

  ! Multi-thread outputs
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:), rtsTL(:,:), rtsKout(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atmK(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfcK(:,:)

  ! ============================================================================
  ! 1. **** INITIALISE ****
  CALL test%Init(.TRUE.)
  CALL test%Setup(PROGRAM_NAME, PROGRAM_NAME, .TRUE.)

  Error_Status = CRTM_Init( SENSOR_ID, ChannelInfo, File_Path=COEFFICIENTS_PATH, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE )
    STOP 1
  END IF
  n_Channels = CRTM_ChannelInfo_n_Channels(ChannelInfo(1))
  WRITE( *,'(5x,"Sensor: ",a," (",i0," ch), ",i0," profiles")' ) &
    TRIM(SENSOR_ID(1)), n_Channels, N_PROFILES

  ! ============================================================================
  ! 2. **** INPUT DATA ****
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating Atmosphere', FAILURE ); STOP 1
  END IF
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! TL perturbation: +0.5 K at every layer, everything else zero
  Atm_TL = Atm
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Sfc_TL = Sfc
  CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    Atm_TL(m)%Temperature = 0.5_fp
  END DO

  ALLOCATE( rts_ref(n_Channels,N_PROFILES), rtsTL_ref(n_Channels,N_PROFILES), &
            rtsKout_ref(n_Channels,N_PROFILES), atmK_ref(n_Channels,N_PROFILES), &
            sfcK_ref(n_Channels,N_PROFILES), &
            rts(n_Channels,N_PROFILES), rtsTL(n_Channels,N_PROFILES), &
            rtsKout(n_Channels,N_PROFILES), atmK(n_Channels,N_PROFILES), &
            sfcK(n_Channels,N_PROFILES) )
  CALL CRTM_Atmosphere_Create( atmK_ref, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( atmK,     N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  ! ============================================================================
  ! 3. **** SINGLE-THREAD REFERENCE ****
  CALL Set_Threads( 1 )
  CALL Run_All( rts_ref, rtsTL_ref, rtsKout_ref, atmK_ref, sfcK_ref )

  ! ============================================================================
  ! 4. **** MULTI-THREAD RUNS, BIT-IDENTICAL TO THE REFERENCE ****
  DO i = 1, N_THREAD_CASES
    CALL Set_Threads( THREAD_CASES(i) )
    CALL Run_All( rts, rtsTL, rtsKout, atmK, sfcK )

    CALL Assert_Equal( 'FWD', ALL(CRTM_RTSolution_Compare(rts, rts_ref)) )
    CALL Assert_Equal( 'TL',  ALL(CRTM_RTSolution_Compare(rtsTL, rtsTL_ref)) )
    CALL Assert_Equal( 'K-rts', ALL(CRTM_RTSolution_Compare(rtsKout, rtsKout_ref)) )
    CALL Assert_Equal( 'K Atmosphere Jacobians', ALL(CRTM_Atmosphere_Compare(atmK, atmK_ref)) )
    CALL Assert_Equal( 'K Surface Jacobians', ALL(CRTM_Surface_Compare(sfcK, sfcK_ref)) )
  END DO

  ! ============================================================================
  ! 5. **** REPORT & CLEAN UP ****
  CALL test%Report()

  Error_Status = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( atmK_ref )
  CALL CRTM_Atmosphere_Destroy( atmK )

  IF ( test%n_Failed() == 0 ) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

CONTAINS

  ! Set the OpenMP thread count CRTM sees (a no-op in a non-OpenMP build).
  SUBROUTINE Set_Threads( n )
    INTEGER, INTENT(IN) :: n
!$  CALL OMP_SET_NUM_THREADS( n )
    WRITE( *,'(5x,"Running with ",i0," thread(s)")' ) n
  END SUBROUTINE Set_Threads

  ! Run FWD, TL and K on the profiles, asserting each call succeeds.
  SUBROUTINE Run_All( rts_out, rtsTL_out, rtsKout_out, atmK_out, sfcK_out )
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: rts_out(:,:), rtsTL_out(:,:), rtsKout_out(:,:)
    TYPE(CRTM_Atmosphere_type), INTENT(IN OUT) :: atmK_out(:,:)
    TYPE(CRTM_Surface_type),    INTENT(IN OUT) :: sfcK_out(:,:)
    TYPE(CRTM_RTSolution_type) :: rtsK_in(n_Channels,N_PROFILES)

    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, rts_out )
    CALL test%Assert( Error_Status == SUCCESS )
    IF ( Error_Status /= SUCCESS ) &
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM_Forward', FAILURE )

    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                        ChannelInfo, rts_out, rtsTL_out )
    CALL test%Assert( Error_Status == SUCCESS )
    IF ( Error_Status /= SUCCESS ) &
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM_Tangent_Linear', FAILURE )

    CALL CRTM_Atmosphere_Zero( atmK_out )
    CALL CRTM_Surface_Zero( sfcK_out )
    rtsK_in%Radiance               = ZERO
    rtsK_in%Brightness_Temperature = ONE
    Error_Status = CRTM_K_Matrix( Atm, Sfc, rtsK_in, Geometry, ChannelInfo, &
                                  atmK_out, sfcK_out, rtsKout_out )
    CALL test%Assert( Error_Status == SUCCESS )
    IF ( Error_Status /= SUCCESS ) &
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM_K_Matrix', FAILURE )
  END SUBROUTINE Run_All

  ! Assert one comparison against the single-thread reference.
  SUBROUTINE Assert_Equal( tag, ok )
    CHARACTER(*), INTENT(IN) :: tag
    LOGICAL,      INTENT(IN) :: ok
    CALL test%Assert( ok )
    IF ( .NOT. ok ) THEN
      WRITE( Message,'(a," results at ",i0," threads differ from the single-thread run")' ) &
        tag, THREAD_CASES(i)
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
  END SUBROUTINE Assert_Equal

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_OMP_Channel_Overshoot
