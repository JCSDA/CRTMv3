!
! test_Fastem1_SST_Jacobian
!
! Validates the sea-surface-temperature (Water_Temperature) Jacobian on the
! legacy Fastem1 MW-water path (Options%Use_Old_MWSSEM=.TRUE., frequency >= 20 GHz).
!
! Fastem1 returns only wind-speed derivatives, so d(emissivity)/d(SST) used to be
! silently zero -> Surface_K%Water_Temperature carried only the skin-emission term
! and disagreed with a finite difference of the forward. The forward now caches a
! central-difference d(emissivity)/d(Water_Temperature), so the analytic K-matrix
! SST Jacobian must match a central finite difference of the forward.
!
! amsua_n19 over ocean; all channels are >= 20 GHz, so all use Fastem1 here.
! Exit status: STOP 0 = success, STOP 1 = failure.
!
PROGRAM test_Fastem1_SST_Jacobian

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_Fastem1_SST_Jacobian'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR_ID         = 'amsua_n19'
  INTEGER,  PARAMETER :: N_PROFILES  = 2   ! matches Load_Atm_Data.inc
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  ! Ocean base state
  REAL(fp), PARAMETER :: TS0 = 290.0_fp, WIND0 = 5.0_fp, SAL0 = 33.0_fp
  REAL(fp), PARAMETER :: DTS = 1.0e-2_fp          ! central-FD SST perturbation (K)
  REAL(fp), PARAMETER :: TOL_REL = 5.0e-3_fp, TOL_ABS = 1.0e-4_fp
  REAL(fp), PARAMETER :: ACTIVE = 1.0e-2_fp       ! surface-sensitive threshold (K/K)

  CHARACTER(256) :: Message, Version
  INTEGER :: Error_Status, Alloc_Status, n_Channels, l, m, n_active
  LOGICAL :: failed
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS_p(:,:), RTS_m(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Surface_K(:,:)
  REAL(fp), ALLOCATABLE :: ad_ts(:,:), fd_ts(:,:)

  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Validate the legacy Fastem1 SST (Water_Temperature) Jacobian vs finite differences.', &
    'CRTM Version: '//TRIM(Version) )

  Error_Status = CRTM_Init( (/SENSOR_ID/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error initializing CRTM', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_K(n_Channels,N_PROFILES), &
            Atmosphere_K(n_Channels,N_PROFILES), Surface_K(n_Channels,N_PROFILES), &
            RTS_p(n_Channels,N_PROFILES), RTS_m(n_Channels,N_PROFILES), &
            ad_ts(n_Channels,N_PROFILES), fd_ts(n_Channels,N_PROFILES), STAT=Alloc_Status )
  IF ( Alloc_Status /= 0 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error allocating arrays', FAILURE ); STOP 1
  END IF
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_RTSolution_Create( RTSolution,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_p,        N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_m,        N_LAYERS )

  CALL Load_Atm_Data()
  CALL Load_Ocean_Surface()
  CALL CRTM_Geometry_SetValue( Geometry, Sensor_Zenith_Angle=ZENITH_ANGLE, Sensor_Scan_Angle=SCAN_ANGLE )
  ! Force the legacy Fastem1 path
  Options%Use_Old_MWSSEM = .TRUE.

  ! Analytic SST Jacobian (K-matrix / adjoint path)
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  CALL CRTM_Surface_Zero( Surface_K )
  RTSolution_K%Radiance = ZERO
  RTSolution_K%Brightness_Temperature = ONE
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                Atmosphere_K, Surface_K, RTSolution, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Error in CRTM K_Matrix', FAILURE ); STOP 1
  END IF
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      ad_ts(l,m) = Surface_K(l,m)%Water_Temperature
    END DO
  END DO

  ! Central finite difference of the forward
  Sfc%Water_Temperature = TS0 + DTS
  CALL Run_Forward( RTS_p, 'Ts+' )
  Sfc%Water_Temperature = TS0 - DTS
  CALL Run_Forward( RTS_m, 'Ts-' )
  fd_ts = (RTS_p%Brightness_Temperature - RTS_m%Brightness_Temperature)/(2.0_fp*DTS)
  Sfc%Water_Temperature = TS0

  failed = .FALSE.; n_active = 0
  WRITE(*,'(/5x,a)') 'Fastem1 SST Jacobian dTb/dWater_Temperature (K/K). Columns: AD / FD'
  WRITE(*,'(5x,a)')  '  m  ch        AD           FD'
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      WRITE(*,'(5x,i3,i4,2x,2es14.5)') m, RTSolution(l,m)%Sensor_Channel, ad_ts(l,m), fd_ts(l,m)
      IF ( ABS(fd_ts(l,m)) > ACTIVE ) n_active = n_active + 1
      CALL Check( 'AD dTb/dTs', m, RTSolution(l,m)%Sensor_Channel, ad_ts(l,m), fd_ts(l,m), failed )
    END DO
  END DO
  WRITE(*,'(/5x,"Surface-sensitive channels (|FD|>",es8.1,"): ",i0)') ACTIVE, n_active
  IF ( n_active < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'No SST-sensitive channels exercised -- test not meaningful', FAILURE )
    failed = .TRUE.
  END IF

  Error_Status = CRTM_Destroy( ChannelInfo )
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atmosphere_K )
  DEALLOCATE( RTSolution, RTSolution_K, Atmosphere_K, Surface_K, RTS_p, RTS_m, ad_ts, fd_ts )

  IF ( failed ) THEN
    CALL Display_Message( PROGRAM_NAME, 'FAILED: Fastem1 SST Jacobian disagrees with finite differences', FAILURE )
    STOP 1
  ELSE
    CALL Display_Message( PROGRAM_NAME, 'PASSED: Fastem1 SST Jacobian matches finite differences', INFORMATION )
    STOP 0
  END IF

CONTAINS

  SUBROUTINE Run_Forward( RTS, label )
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTS(:,:)
    CHARACTER(*),               INTENT(IN)     :: label
    INTEGER :: stat
    stat = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
    IF ( stat /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Error in CRTM Forward ('//label//')', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_Forward

  SUBROUTINE Check( name, m, ch, analytic, fd, failed )
    CHARACTER(*), INTENT(IN)     :: name
    INTEGER,      INTENT(IN)     :: m, ch
    REAL(fp),     INTENT(IN)     :: analytic, fd
    LOGICAL,      INTENT(IN OUT) :: failed
    REAL(fp) :: tol
    tol = TOL_ABS + TOL_REL*ABS(fd)
    IF ( ABS(analytic - fd) > tol ) THEN
      WRITE(Message,'(a," mismatch: profile ",i0," channel ",i0,": analytic=",es13.5,&
            &" FD=",es13.5," |diff|=",es11.3," tol=",es11.3)') &
            TRIM(name), m, ch, analytic, fd, ABS(analytic-fd), tol
      CALL Display_Message( PROGRAM_NAME, TRIM(Message), FAILURE ); failed = .TRUE.
    END IF
  END SUBROUTINE Check

  SUBROUTINE Load_Ocean_Surface()
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm)%Water_Coverage    = 1.0_fp
      Sfc(mm)%Land_Coverage     = 0.0_fp
      Sfc(mm)%Snow_Coverage     = 0.0_fp
      Sfc(mm)%Ice_Coverage      = 0.0_fp
      Sfc(mm)%Water_Type        = 1          ! SEA_WATER
      Sfc(mm)%Water_Temperature = TS0
      Sfc(mm)%Wind_Speed        = WIND0
      Sfc(mm)%Salinity          = SAL0
    END DO
  END SUBROUTINE Load_Ocean_Surface

  INCLUDE 'Load_Atm_Data.inc'

END PROGRAM test_Fastem1_SST_Jacobian
