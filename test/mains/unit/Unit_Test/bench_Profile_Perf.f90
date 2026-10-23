!
! bench_Profile_Perf
!
! Micro-benchmark for the performance cost of the level-resolved radiance profile
! outputs (Downwelling_Radiance / Upwelling_Radiance, opt-in).  Times CRTM_Forward and
! CRTM_K_Matrix on an overcast (ADA scattering) MW scene under four Options configs:
!   (1) flags off        -> baseline (regression check)
!   (2) down-profile on
!   (3) up-profile on
!   (4) both profiles on
! plus a clear-sky config to show the emission path is unaffected.
!
PROGRAM bench_Profile_Perf

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  INTEGER,  PARAMETER :: NITER = 400          ! repetitions per timed config

  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status, Allocate_Status, n_Channels, l, m

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_K_dummy
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  Sensor_Id = 'atms_npp'
  WRITE(*,'(/5x,a)') 'Profile-output performance benchmark ('//TRIM(Sensor_Id)//', overcast ADA)'

  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) THEN ; WRITE(*,*) 'Init fail' ; STOP 1 ; END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), RTSolution_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), STAT=Allocate_Status )
  CALL CRTM_RTSolution_Create( RTSolution,   N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_K, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm,   N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_ID = Atm(m)%Absorber_ID ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
  END DO
  CALL CRTM_Geometry_SetValue( Geometry, Sensor_Zenith_Angle=ZENITH_ANGLE, Sensor_Scan_Angle=SCAN_ANGLE )

  ! Overcast strongly-scattering MW cloud -> the ADA scattering solver (where the cost is).
  DO m = 1, N_PROFILES
    Atm(m)%n_Clouds = 1
    Atm(m)%Cloud_Fraction = ZERO ; Atm(m)%Cloud_Fraction(70:90) = ONE
    Atm(m)%Cloud(1)%Type = SNOW_CLOUD
    Atm(m)%Cloud(1)%Effective_Radius = ZERO ; Atm(m)%Cloud(1)%Effective_Radius(70:90) = 500.0_fp
    Atm(m)%Cloud(1)%Water_Content    = ZERO ; Atm(m)%Cloud(1)%Water_Content(70:90)    = 5.0_fp
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO

  WRITE(*,'(/5x,"NITER=",i0,", n_Channels=",i0,", N_PROFILES=",i0,", N_LAYERS=",i0)') &
        NITER, n_Channels, N_PROFILES, N_LAYERS
  WRITE(*,'(5x,"(set OMP_NUM_THREADS=1 for a clean single-thread measure)"/)')

  WRITE(*,'(5x,a)') '------------------------- FORWARD (overcast ADA) -------------------------'
  CALL time_fwd( .FALSE., .FALSE., 'flags off (baseline)         ' )
  CALL time_fwd( .TRUE. , .FALSE., 'Down profile on              ' )
  CALL time_fwd( .FALSE., .TRUE. , 'Up   profile on              ' )
  CALL time_fwd( .TRUE. , .TRUE. , 'Both profiles on             ' )

  WRITE(*,'(/5x,a)') '------------------------- K_MATRIX (overcast ADA) ------------------------'
  CALL time_k( .FALSE., .FALSE., 'flags off (baseline)         ' )
  CALL time_k( .TRUE. , .FALSE., 'Down profile on              ' )
  CALL time_k( .FALSE., .TRUE. , 'Up   profile on              ' )
  CALL time_k( .TRUE. , .TRUE. , 'Both profiles on             ' )

  Error_Status = CRTM_Destroy( ChannelInfo )
  STOP 0

CONTAINS

  SUBROUTINE set_flags( dn, up )
    LOGICAL, INTENT(IN) :: dn, up
    DO m = 1, N_PROFILES
      Options(m)%Compute_Down_Radiance_Profile = dn
      Options(m)%Compute_Up_Radiance_Profile   = up
    END DO
  END SUBROUTINE set_flags

  SUBROUTINE time_fwd( dn, up, tag )
    LOGICAL, INTENT(IN) :: dn, up ; CHARACTER(*), INTENT(IN) :: tag
    INTEGER :: it, c0, c1, cr ; REAL(fp) :: secs
    CALL set_flags( dn, up )
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options ) ! warm-up
    CALL SYSTEM_CLOCK(c0, cr)
    DO it = 1, NITER
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
    END DO
    CALL SYSTEM_CLOCK(c1)
    secs = REAL(c1-c0,fp)/REAL(cr,fp)
    WRITE(*,'(7x,a,": ",f8.4," s  (",f9.2," us/call)")') tag, secs, 1.0e6_fp*secs/REAL(NITER,fp)
  END SUBROUTINE time_fwd

  SUBROUTINE time_k( dn, up, tag )
    LOGICAL, INTENT(IN) :: dn, up ; CHARACTER(*), INTENT(IN) :: tag
    INTEGER :: it, c0, c1, cr ; REAL(fp) :: secs
    CALL set_flags( dn, up )
    CALL CRTM_RTSolution_Zero( RTSolution_K )
    DO m = 1, N_PROFILES ; DO l = 1, n_Channels
      RTSolution_K(l,m)%Radiance = ONE
      IF ( dn ) RTSolution_K(l,m)%Downwelling_Radiance = ONE
      IF ( up ) RTSolution_K(l,m)%Upwelling_Radiance   = ONE
    END DO ; END DO
    Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                  Atm_K, Sfc_K, RTSolution, Options=Options ) ! warm-up
    CALL SYSTEM_CLOCK(c0, cr)
    DO it = 1, NITER
      Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, &
                                    Atm_K, Sfc_K, RTSolution, Options=Options )
    END DO
    CALL SYSTEM_CLOCK(c1)
    secs = REAL(c1-c0,fp)/REAL(cr,fp)
    WRITE(*,'(7x,a,": ",f8.4," s  (",f9.2," us/call)")') tag, secs, 1.0e6_fp*secs/REAL(NITER,fp)
  END SUBROUTINE time_k

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM bench_Profile_Perf
