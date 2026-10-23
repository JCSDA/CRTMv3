!
! test_MultiSensor_SingleCall
!
! Consistency test for a single CRTM call covering MULTIPLE sensors
! (ChannelInfo(1:2)) against the equivalent per-sensor calls.
!
! The RTSolution (and K-matrix) arrays are indexed by a cumulative channel
! counter 'ln' that must carry the previous sensors' channel count across the
! Sensor_Loop in the OpenMP channel-thread loops of the Forward, Tangent-Linear
! and K-Matrix modules. A bug that resets 'ln' per sensor makes sensor 2
! overwrite sensor 1's outputs, leaving the tail of the arrays untouched --
! invisible to every other test because they all pass ChannelInfo one sensor
! at a time. This test pins the combined call bit-for-bit to the per-sensor
! calls for FWD, TL and K.
!
PROGRAM test_MultiSensor_SingleCall

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  USE CRTM_Module
  USE UnitTest_Define, ONLY: UnitTest_type
  IMPLICIT NONE
  ! ============================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_MultiSensor_SingleCall'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Profile dimensions (matching the Load_*_Data include files)
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 1
  INTEGER, PARAMETER :: N_AEROSOLS  = 1

  ! Two microwave sensors processed in ONE call
  INTEGER, PARAMETER :: N_SENSORS = 2
  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = &
    (/ 'amsua_metop-a', 'mhs_n19      ' /)

  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  INTEGER :: Error_Status, Allocate_Status
  INTEGER :: n_Channels, n1, n2
  INTEGER :: n, l, m, l0

  TYPE(UnitTest_type) :: test

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)

  ! Combined-call outputs (all sensors at once)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_all(:,:), rtsTL_all(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rtsK_all(:,:), rtsKout_all(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atmK_all(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfcK_all(:,:)

  ! Per-sensor-call outputs
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_one(:,:), rtsTL_one(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rtsK_one(:,:), rtsKout_one(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atmK_one(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: sfcK_one(:,:)

  ! ============================================================================
  ! 1. **** INITIALISE ****
  CALL test%Init(.TRUE.)
  CALL test%Setup(PROGRAM_NAME, PROGRAM_NAME, .TRUE.)

  Error_Status = CRTM_Init( SENSOR_ID, ChannelInfo, File_Path=COEFFICIENTS_PATH, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE )
    STOP 1
  END IF

  n1 = CRTM_ChannelInfo_n_Channels(ChannelInfo(1))
  n2 = CRTM_ChannelInfo_n_Channels(ChannelInfo(2))
  n_Channels = n1 + n2
  WRITE( *,'(5x,"Sensors: ",a," (",i0," ch) + ",a," (",i0," ch)")' ) &
    TRIM(SENSOR_ID(1)), n1, TRIM(SENSOR_ID(2)), n2

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

  ! ============================================================================
  ! 3. **** FORWARD: combined vs per-sensor ****
  ALLOCATE( rts_all(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) STOP 1
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, rts_all )
  CALL test%Assert( Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in combined CRTM_Forward', FAILURE ); STOP 1
  END IF

  ! The direct symptom of the 'ln' reset bug: sensor 2's block must carry
  ! sensor 2's identity (the bug leaves sensor 1's channels there untouched).
  CALL test%Assert( TRIM(rts_all(n1+1,1)%Sensor_Id) == TRIM(SENSOR_ID(2)) )
  IF ( TRIM(rts_all(n1+1,1)%Sensor_Id) /= TRIM(SENSOR_ID(2)) ) THEN
    WRITE( Message,'("RTSolution(",i0,") Sensor_Id is ",a," -- expected ",a)' ) &
      n1+1, TRIM(rts_all(n1+1,1)%Sensor_Id), TRIM(SENSOR_ID(2))
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  END IF

  l0 = 0
  DO n = 1, N_SENSORS
    ALLOCATE( rts_one(CRTM_ChannelInfo_n_Channels(ChannelInfo(n)),N_PROFILES) )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo(n:n), rts_one )
    CALL test%Assert( Error_Status == SUCCESS )
    CALL Assert_RTS_Block_Equal( 'FWD', n, rts_all, rts_one, l0 )
    l0 = l0 + SIZE(rts_one,DIM=1)
    DEALLOCATE( rts_one )
  END DO

  ! ============================================================================
  ! 4. **** TANGENT-LINEAR: combined vs per-sensor ****
  ALLOCATE( rtsTL_all(n_Channels,N_PROFILES) )
  Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                      ChannelInfo, rts_all, rtsTL_all )
  CALL test%Assert( Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in combined CRTM_Tangent_Linear', FAILURE ); STOP 1
  END IF

  l0 = 0
  DO n = 1, N_SENSORS
    ALLOCATE( rts_one(CRTM_ChannelInfo_n_Channels(ChannelInfo(n)),N_PROFILES), &
              rtsTL_one(CRTM_ChannelInfo_n_Channels(ChannelInfo(n)),N_PROFILES) )
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                        ChannelInfo(n:n), rts_one, rtsTL_one )
    CALL test%Assert( Error_Status == SUCCESS )
    CALL Assert_RTS_Block_Equal( 'TL', n, rtsTL_all, rtsTL_one, l0 )
    l0 = l0 + SIZE(rtsTL_one,DIM=1)
    DEALLOCATE( rts_one, rtsTL_one )
  END DO

  ! ============================================================================
  ! 5. **** K-MATRIX: combined vs per-sensor ****
  ALLOCATE( atmK_all(n_Channels,N_PROFILES), sfcK_all(n_Channels,N_PROFILES), &
            rtsK_all(n_Channels,N_PROFILES), rtsKout_all(n_Channels,N_PROFILES) )
  CALL CRTM_Atmosphere_Create( atmK_all, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Zero( atmK_all )
  CALL CRTM_Surface_Zero( sfcK_all )
  rtsK_all%Radiance               = ZERO
  rtsK_all%Brightness_Temperature = ONE

  Error_Status = CRTM_K_Matrix( Atm, Sfc, rtsK_all, Geometry, ChannelInfo, &
                                atmK_all, sfcK_all, rtsKout_all )
  CALL test%Assert( Error_Status == SUCCESS )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in combined CRTM_K_Matrix', FAILURE ); STOP 1
  END IF

  l0 = 0
  DO n = 1, N_SENSORS
    l = CRTM_ChannelInfo_n_Channels(ChannelInfo(n))
    ALLOCATE( atmK_one(l,N_PROFILES), sfcK_one(l,N_PROFILES), &
              rtsK_one(l,N_PROFILES), rtsKout_one(l,N_PROFILES) )
    CALL CRTM_Atmosphere_Create( atmK_one, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    CALL CRTM_Atmosphere_Zero( atmK_one )
    CALL CRTM_Surface_Zero( sfcK_one )
    rtsK_one%Radiance               = ZERO
    rtsK_one%Brightness_Temperature = ONE

    Error_Status = CRTM_K_Matrix( Atm, Sfc, rtsK_one, Geometry, ChannelInfo(n:n), &
                                  atmK_one, sfcK_one, rtsKout_one )
    CALL test%Assert( Error_Status == SUCCESS )

    CALL Assert_RTS_Block_Equal( 'K-rts', n, rtsKout_all, rtsKout_one, l0 )
    CALL test%Assert( ALL(CRTM_Atmosphere_Compare( atmK_all(l0+1:l0+l,:), atmK_one )) )
    IF ( .NOT. ALL(CRTM_Atmosphere_Compare( atmK_all(l0+1:l0+l,:), atmK_one )) ) THEN
      WRITE( Message,'("K Atmosphere Jacobians differ for sensor ",i0," (",a,")")' ) &
        n, TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    CALL test%Assert( ALL(CRTM_Surface_Compare( sfcK_all(l0+1:l0+l,:), sfcK_one )) )
    IF ( .NOT. ALL(CRTM_Surface_Compare( sfcK_all(l0+1:l0+l,:), sfcK_one )) ) THEN
      WRITE( Message,'("K Surface Jacobians differ for sensor ",i0," (",a,")")' ) &
        n, TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF

    l0 = l0 + l
    CALL CRTM_Atmosphere_Destroy( atmK_one )
    DEALLOCATE( atmK_one, sfcK_one, rtsK_one, rtsKout_one )
  END DO

  ! ============================================================================
  ! 6. **** REPORT & CLEAN UP ****
  CALL test%Report()

  Error_Status = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( atmK_all )
  DEALLOCATE( rts_all, rtsTL_all, atmK_all, sfcK_all, rtsK_all, rtsKout_all, &
              STAT=Allocate_Status )

  IF ( test%n_Failed() == 0 ) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

CONTAINS

  ! Compare one sensor's block of the combined-call RTSolution array against
  ! the per-sensor call, reporting the first differing channel/profile.
  SUBROUTINE Assert_RTS_Block_Equal( tag, sensor, all_rts, one_rts, offset )
    CHARACTER(*),               INTENT(IN) :: tag
    INTEGER,                    INTENT(IN) :: sensor
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: all_rts(:,:), one_rts(:,:)
    INTEGER,                    INTENT(IN) :: offset
    LOGICAL :: ok(SIZE(one_rts,DIM=1),SIZE(one_rts,DIM=2))
    INTEGER :: il, im
    ok = CRTM_RTSolution_Compare( all_rts(offset+1:offset+SIZE(one_rts,DIM=1),:), one_rts )
    CALL test%Assert( ALL(ok) )
    IF ( .NOT. ALL(ok) ) THEN
      DO im = 1, SIZE(ok,DIM=2)
        DO il = 1, SIZE(ok,DIM=1)
          IF ( .NOT. ok(il,im) ) THEN
            WRITE( Message,'(a," results differ: sensor ",i0," (",a,"), channel ",i0, &
                   &", profile ",i0," (combined index ",i0,")")' ) &
              tag, sensor, TRIM(SENSOR_ID(sensor)), il, im, offset+il
            CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
            RETURN
          END IF
        END DO
      END DO
    END IF
  END SUBROUTINE Assert_RTS_Block_Equal

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_MultiSensor_SingleCall
