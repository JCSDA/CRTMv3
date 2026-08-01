!
! dump_scalar_fullprec
!
! Prints scalar-path (n_Stokes = 1) results at full double precision, for a
! bit-identity comparison between two builds. Not a registered test: it asserts
! nothing. It exists so that "the default scalar path is unchanged" can be a
! measurement rather than an argument.
!
! Why it is needed
! ---------------
! The regression suite compares against its stored references at
! DEFAULT_N_SIGFIG, which is SP_N_SIGFIG, roughly six significant figures. That
! cannot detect a change in the last bits, so a green suite is not by itself
! proof of bit-identity.
!
! Coverage, chosen to hit everything the polarimetric work touched that a
! scalar user can reach
! ---------------------------------------------------------------------------
!   Entry points : Forward, Tangent_Linear, Adjoint, K_Matrix. The azimuthal
!                  Fourier accumulation was refactored in all four, so covering
!                  Forward alone would leave three quarters of that change
!                  unmeasured.
!   Sensors      : a microwave sounder and a VISIBLE imager. The visible sensor
!                  is the important one: it is the only class for which
!                  n_Azi > 0, so it is the only case in which the accumulation
!                  weight COS(mth_Azi*dphi) is evaluated at a non-zero
!                  argument. The geometry below deliberately sets the sensor
!                  and source azimuths apart so dphi is not zero.
!   Surfaces     : ocean and land. The microwave coverage aggregation was
!                  changed at the water sites only, so both need checking.
!   Cloud states : clear, overcast and fractional. The fractional case is the
!                  one that exercises RTV_Clear, whose n_Stokes plumbing was
!                  changed in the tangent linear, adjoint and K-matrix modules.
!
! Usage: build in both trees, run both, diff the output. Any difference at all
! is a change to the scalar path.
!

PROGRAM dump_scalar_fullprec

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'dump_scalar_fullprec'
  CHARACTER(*), PARAMETER :: PATH = './testinput/'
  ! Microwave sounder plus a visible imager: the visible one is what drives
  ! n_Azi > 0 and so the only non-trivial evaluation of the cosine weight.
  CHARACTER(*), PARAMETER :: SENSORS(2) = (/ 'amsua_n19', 'v.abi_g18' /)

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86
  INTEGER,  PARAMETER :: KP  = 82           ! probed layer

  REAL(fp), PARAMETER :: ZENITH     = 40.0_fp
  REAL(fp), PARAMETER :: SENSOR_AZI = 60.0_fp
  REAL(fp), PARAMETER :: SOURCE_ZEN = 45.0_fp
  REAL(fp), PARAMETER :: SOURCE_AZI = 30.0_fp   ! != SENSOR_AZI on purpose

  INTEGER :: Error_Status, Allocate_Status, n_Channels, l, m, isfc, icase
  REAL(fp) :: wc, cf

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(SIZE(SENSORS))
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES), Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:), RTS_TL(:,:), RTS_AD(:,:), RTS_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  Error_Status = CRTM_Init( SENSORS, ChannelInfo, File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTS(n_Channels,N_PROFILES), RTS_TL(n_Channels,N_PROFILES), &
            RTS_AD(n_Channels,N_PROFILES), RTS_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTS,    N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_TL, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_AD, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTS_K,  N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm,    N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_AD, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  CALL CRTM_Atmosphere_Create( Atm_K,  N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )

  CALL Load_ECMWF84_Atm_Data()

  DO m = 1, N_PROFILES
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_TL(m)%Cloud(1)%Type = SNOW_CLOUD
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Cloud(1)%Type = SNOW_CLOUD
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
      Atm_K(l,m)%Cloud(1)%Type = SNOW_CLOUD
    END DO
    CALL CRTM_Geometry_SetValue( Geometry(m), &
           Sensor_Zenith_Angle  = ZENITH,     Sensor_Azimuth_Angle = SENSOR_AZI, &
           Source_Zenith_Angle  = SOURCE_ZEN, Source_Azimuth_Angle = SOURCE_AZI )
    Options(m)%n_Stokes = 1
  END DO

  ! isfc = 1 ocean, 2 land ; icase = 1 clear, 2 overcast, 3 fractional
  DO isfc = 1, 2
    DO m = 1, N_PROFILES
      Sfc(m) = CRTM_Surface_type()
      IF ( isfc == 1 ) THEN
        Sfc(m)%Water_Coverage    = ONE
        Sfc(m)%Water_Type        = 1
        Sfc(m)%Water_Temperature = 290.0_fp
        Sfc(m)%Wind_Speed        = 12.0_fp
        Sfc(m)%Wind_Direction    = 100.0_fp
        Sfc(m)%Salinity          = 33.0_fp
      ELSE
        Sfc(m)%Land_Coverage     = ONE
        Sfc(m)%Land_Type         = 1
        Sfc(m)%Land_Temperature  = 288.0_fp
        Sfc(m)%Soil_Moisture_Content = 0.2_fp
        Sfc(m)%Vegetation_Fraction   = 0.4_fp
        Sfc(m)%LAI                   = 2.0_fp
      END IF
    END DO

    DO icase = 1, 3
      SELECT CASE ( icase )
        CASE (1) ; wc = ZERO   ; cf = ZERO
        CASE (2) ; wc = 1.0_fp ; cf = ONE
        CASE (3) ; wc = 1.0_fp ; cf = 0.5_fp
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

      ! ---- Forward ----
      Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'FWD failed', FAILURE ); STOP 1
      END IF
      DO m = 1, N_PROFILES
        DO l = 1, n_Channels
          WRITE(*,'("FWD s",i1," c",i1," p",i1," ch",i4,2(1x,ES26.17E3))') &
                isfc, icase, m, l, RTS(l,m)%Radiance, RTS(l,m)%Brightness_Temperature
        END DO
      END DO

      ! ---- Tangent linear ----
      CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
      DO m = 1, N_PROFILES
        Atm_TL(m)%Temperature(KP)            = ONE
        Atm_TL(m)%Cloud(1)%Water_Content(KP) = 0.1_fp
      END DO
      Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                          ChannelInfo, RTS, RTS_TL, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'TL failed', FAILURE ); STOP 1
      END IF
      DO m = 1, N_PROFILES
        DO l = 1, n_Channels
          WRITE(*,'("TL  s",i1," c",i1," p",i1," ch",i4,1x,ES26.17E3)') &
                isfc, icase, m, l, RTS_TL(l,m)%Radiance
        END DO
      END DO

      ! ---- Adjoint ----
      CALL CRTM_RTSolution_Zero( RTS_AD )
      DO m = 1, N_PROFILES
        DO l = 1, n_Channels
          RTS_AD(l,m)%Radiance = ONE
        END DO
      END DO
      CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
      Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                                   Atm_AD, Sfc_AD, RTS, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'AD failed', FAILURE ); STOP 1
      END IF
      DO m = 1, N_PROFILES
        WRITE(*,'("AD  s",i1," c",i1," p",i1,3(1x,ES26.17E3))') &
              isfc, icase, m, Atm_AD(m)%Temperature(KP), &
              Atm_AD(m)%Cloud(1)%Water_Content(KP), Sfc_AD(m)%Water_Temperature
      END DO

      ! ---- K matrix ----
      CALL CRTM_RTSolution_Zero( RTS_K )
      DO m = 1, N_PROFILES
        DO l = 1, n_Channels
          RTS_K(l,m)%Radiance = ONE
        END DO
      END DO
      CALL CRTM_Atmosphere_Zero( Atm_K ) ; CALL CRTM_Surface_Zero( Sfc_K )
      Error_Status = CRTM_K_Matrix( Atm, Sfc, RTS_K, Geometry, ChannelInfo, &
                                    Atm_K, Sfc_K, RTS, Options=Options )
      IF ( Error_Status /= SUCCESS ) THEN
        CALL Display_Message( PROGRAM_NAME, 'K failed', FAILURE ); STOP 1
      END IF
      DO m = 1, N_PROFILES
        DO l = 1, n_Channels
          WRITE(*,'("K   s",i1," c",i1," p",i1," ch",i4,2(1x,ES26.17E3))') &
                isfc, icase, m, l, Atm_K(l,m)%Temperature(KP), Sfc_K(l,m)%Water_Temperature
        END DO
      END DO

    END DO
  END DO

  Error_Status = CRTM_Destroy( ChannelInfo )

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM dump_scalar_fullprec
