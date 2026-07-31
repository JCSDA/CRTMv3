!
! dump_scalar_fullprec
!
! Prints scalar-path (n_Stokes = 1) forward radiances and brightness
! temperatures at full double precision, for a bit-identity comparison between
! two builds. Not a registered test: it asserts nothing. It exists so that
! "the default scalar path is unchanged" can be a measurement rather than an
! argument.
!
! The regression suite compares against its stored references at
! DEFAULT_N_SIGFIG, which is SP_N_SIGFIG, roughly six significant figures. That
! is not sensitive enough to detect a change in the last bits, so a green suite
! is not by itself proof of bit-identity.
!
! Usage: build in both trees, run both, diff the output. Any difference at all
! is a change to the scalar path.
!
! The scene deliberately exercises the code paths touched by the polarimetric
! work while running scalar: microwave over ocean (the coverage aggregation),
! clear sky (the emission solver and its Stokes completion, which must return
! immediately at n_Stokes = 1), and a cloudy column (the scattering solver and
! the azimuthal Fourier accumulation).
!

PROGRAM dump_scalar_fullprec

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'dump_scalar_fullprec'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'amsua_n19'

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86

  INTEGER :: Error_Status, Allocate_Status, n_Channels, l, m, icase
  REAL(fp) :: wc, cf

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, &
                            File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_ECMWF84_Atm_Data()

  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 12.0_fp
    Sfc(m)%Wind_Direction    = 100.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle  = ZENITH, &
                                              Sensor_Azimuth_Angle = 40.0_fp )
    Options(m)%n_Stokes = 1
  END DO

  ! Three scenes: clear, overcast, fractional
  DO icase = 1, 3
    SELECT CASE ( icase )
      CASE (1) ; wc = ZERO      ; cf = ZERO
      CASE (2) ; wc = 1.0_fp    ; cf = ONE
      CASE (3) ; wc = 1.0_fp    ; cf = 0.5_fp
    END SELECT
    DO m = 1, N_PROFILES
      Atm(m)%n_Clouds                           = N_CLOUDS
      Atm(m)%Cloud_Fraction                     = ZERO
      Atm(m)%Cloud(1)%Type                      = SNOW_CLOUD
      Atm(m)%Cloud(1)%Effective_Radius          = ZERO
      Atm(m)%Cloud(1)%Water_Content             = ZERO
      Atm(m)%Cloud_Fraction(KC1:KC2)            = cf
      Atm(m)%Cloud(1)%Effective_Radius(KC1:KC2) = 500.0_fp
      Atm(m)%Cloud(1)%Water_Content(KC1:KC2)    = wc
    END DO

    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
    END IF

    DO m = 1, N_PROFILES
      DO l = 1, n_Channels
        WRITE(*,'("case",i2," prof",i2," ch",i3,"  R=",ES26.17E3,"  Tb=",ES26.17E3)') &
              icase, m, l, RTSolution(l,m)%Radiance, RTSolution(l,m)%Brightness_Temperature
      END DO
    END DO
  END DO

  Error_Status = CRTM_Destroy( ChannelInfo )

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM dump_scalar_fullprec
