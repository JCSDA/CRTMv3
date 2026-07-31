!
! test_VectorRT_StokesOutput
!
! Proves that the third and fourth Stokes components survive the azimuthal
! Fourier accumulation and reach RTSolution%Stokes, and that the polarimetric
! reference frame established at the surface is preserved end to end by the
! vector solver.
!
! Background
! ----------
! Emergent radiances are written into RTSolution%Stokes in Assign_Common_Output
! (Common_RTSolution.f90) as an azimuthal Fourier accumulation:
!
!     Stokes(1:2)          += Radiance(1:2) * COS( mth_Azi * dphi )
!     Stokes(3:n_Stokes)   += Radiance(3:n) * SIN( mth_Azi * dphi )
!
! That cosine/sine split is the standard convention for a solar problem, where
! the m = 0 Fourier coefficients of U and V vanish identically. But CRTM sets
! n_Azi > 0 only for visible channels (CRTM_Forward_Module.f90:993 versus
! :1011), and the coupled polarimetric surface branch exists only for
! microwave. For every configuration in which n_Stokes > 1 is meaningful,
! mth_Azi is therefore always 0, SIN(0) is 0, and components 3 and 4 are
! multiplied by zero on their way out. The solver's U and V were discarded at
! the last step regardless of what it computed.
!
! In that single m = 0 solve the azimuth dependence is not carried by the
! Fourier series at all: it is carried by the surface, which is evaluated at
! the actual relative wind azimuth. The correct accumulation weight for m = 0
! is therefore unity, exactly as it already is for components 1 and 2.
!
! Changing the m = 0 weight cannot perturb a solar or visible run. At m = 0 the
! generalized spherical function T_l^m (RTV%Pminus) vanishes identically:
! Gl2n (CRTM_Utility.f90:1295) drops its n argument when MF = 0, in both the
! seed and the recursion, so Pminus = (Gl2n(-2) - Gl2n(2))/2 is exactly zero.
! Every phase-matrix block carrying a Pminus factor, which is all of (1,3),
! (3,1), (2,3), (3,2), (2,4) and (4,2), vanishes with it, leaving the m = 0
! phase matrix block diagonal in {I,Q} and {U,V}. The visible and infrared
! surface sets component 1 only and the thermal source is intensity only, so
! at m = 0 the U and V sources are zero and so is their solution: the quantity
! whose weight changes is identically zero on those paths.
!
! What this test does
! -------------------
! It runs one overcast snow column over ocean at n_Stokes = 4, as two profiles
! that are identical in every respect except the relative wind azimuth, which
! is +phi for the first and -phi for the second. That is a reflection of the
! scene through the vertical plane containing the view direction, under which
! a Stokes vector referred to the meridional frame must transform as
!
!     (I, Q, U, V)  ->  (I, Q, -U, -V) .
!
! It then asserts
!
!   (a) some channel has a non-zero third and fourth Stokes component;
!   (b) Stokes 1 and 2 are even under the reflection;
!   (c) Stokes 3 and 4 are odd under the reflection.
!
! Assertion (a) fails against the unfixed code at exactly zero. Assertions (b)
! and (c) then confirm that the solver transports the surface's polarimetric
! signal without corrupting the frame it is referred to; they hold exactly
! because the m = 0 phase matrix is block diagonal, so negating the U and V
! surface source negates the U and V solution and leaves I and Q untouched.
!
! Gating
! ------
! The n_Stokes > 1 scattering branch needs a cloud lookup table with at least
! six phase elements, so this uses the experimental CRTM-Exp scheme exactly as
! test_VectorRT_TLADK does. It also loads FASTEM4, because the FASTEM6 default
! has no third or fourth Stokes azimuth model and would leave the surface with
! no polarimetric signal to transport (see test_VectorRT_SurfaceFrame).
!

PROGRAM test_VectorRT_StokesOutput

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_StokesOutput'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'
  CHARACTER(*), PARAMETER :: LUT    = 'CloudCoeff_Exp_Full6.nc'
  CHARACTER(*), PARAMETER :: MWWATER_SCHEME = 'FASTEM4'

  ! Column setup: same overcast snow band as test_VectorRT_TLADK, which keeps
  ! the optical depth away from the thin-cloud conditioning problem at
  ! sub-millimetre frequencies.
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86
  REAL(fp), PARAMETER :: REFF_S = 500.0_fp
  REAL(fp), PARAMETER :: WC_S   = 1.0_fp

  ! Relative azimuth built as WIND_DIR - SENSOR_AZI, so the mirrored profile is
  ! the exact negation of the direct one and the symmetry checks are exact.
  REAL(fp), PARAMETER :: WIND_DIR     = 100.0_fp
  REAL(fp), PARAMETER :: SENSOR_AZI_P =  40.0_fp   ! -> +60
  REAL(fp), PARAMETER :: SENSOR_AZI_M = 160.0_fp   ! -> -60
  REAL(fp), PARAMETER :: WIND_SPEED   =  12.0_fp

  ! Non-degeneracy floor in radiance units, far above round-off.
  REAL(fp), PARAMETER :: SIGNAL_FLOOR = 1.0e-12_fp
  ! Symmetry tolerance, relative to the intensity scale of the scene.
  REAL(fp), PARAMETER :: TOL_REL = 1.0e-10_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels
  INTEGER :: l, m
  LOGICAL :: ok_signal, ok_even, ok_odd, all_ok
  REAL(fp) :: d_even, d_odd, scale, max_U, max_V

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT Stokes output and frame-preservation verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  ! --------------------------------------------------------------------------
  ! Initialize with the experimental cloud optics (>= 6 phase elements) and the
  ! FASTEM4 microwave water scheme (carries third/fourth Stokes azimuth terms)
  ! --------------------------------------------------------------------------
  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo,                &
                            Cloud_Model         = 'CRTM-Exp',         &
                            CloudCoeff_File     = LUT,                &
                            CloudCoeff_Format   = 'netCDF',           &
                            MWwaterCoeff_Scheme = MWWATER_SCHEME,     &
                            File_Path           = PATH,               &
                            Quiet               = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init (Cloud_Model=CRTM-Exp) failed', FAILURE )
    STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'no channels loaded for '//SENSOR, FAILURE )
    STOP 1
  END IF

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )

  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Atmosphere_Create failed', FAILURE )
    STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()          ! fills Atm(1)
  Atm(2) = Atm(1)
  DO m = 1, N_PROFILES
    Atm(m)%n_Clouds                           = 1
    Atm(m)%Cloud_Fraction                     = ZERO
    Atm(m)%Cloud_Fraction(KC1:KC2)            = ONE     ! overcast: isolate the solver
    Atm(m)%Cloud(1)%Type                      = SNOW_CLOUD
    Atm(m)%Cloud(1)%Effective_Radius          = ZERO
    Atm(m)%Cloud(1)%Water_Content             = ZERO
    Atm(m)%Cloud(1)%Effective_Radius(KC1:KC2) = REFF_S
    Atm(m)%Cloud(1)%Water_Content(KC1:KC2)    = WC_S
  END DO

  ! Identical ocean scenes; the two profiles differ only by the sign of the
  ! relative wind azimuth, which is the mirror reflection through the view plane.
  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = WIND_SPEED
    Sfc(m)%Wind_Direction    = WIND_DIR
    Sfc(m)%Salinity          = 33.0_fp
  END DO
  CALL CRTM_Geometry_SetValue( Geometry(1), Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = SENSOR_AZI_P )
  CALL CRTM_Geometry_SetValue( Geometry(2), Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = SENSOR_AZI_M )

  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 4
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO

  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
  END IF

  ! ------------------------------------------------
  ! Assertions
  ! ------------------------------------------------
  max_U  = ZERO
  max_V  = ZERO
  d_even = ZERO
  d_odd  = ZERO
  scale  = ZERO
  DO l = 1, n_Channels
    max_U = MAX( max_U, ABS(RTSolution(l,1)%Stokes(3)) )
    max_V = MAX( max_V, ABS(RTSolution(l,1)%Stokes(4)) )
    scale = MAX( scale, ABS(RTSolution(l,1)%Stokes(1)) )
    ! (b) Stokes 1,2 even under the reflection
    d_even = MAX( d_even, ABS(RTSolution(l,1)%Stokes(1) - RTSolution(l,2)%Stokes(1)) )
    d_even = MAX( d_even, ABS(RTSolution(l,1)%Stokes(2) - RTSolution(l,2)%Stokes(2)) )
    ! (c) Stokes 3,4 odd under the reflection
    d_odd  = MAX( d_odd,  ABS(RTSolution(l,1)%Stokes(3) + RTSolution(l,2)%Stokes(3)) )
    d_odd  = MAX( d_odd,  ABS(RTSolution(l,1)%Stokes(4) + RTSolution(l,2)%Stokes(4)) )
  END DO

  ok_signal = ( max_U > SIGNAL_FLOOR .AND. max_V > SIGNAL_FLOOR )
  ok_even   = ( d_even < TOL_REL*MAX(scale,ONE) )
  ok_odd    = ( d_odd  < TOL_REL*MAX(scale,ONE) )

  WRITE(*,'(5x,a,i0,a)') 'sensor '//SENSOR//', ', n_Channels, &
        ' channels, n_Stokes = 4, overcast snow over ocean'
  WRITE(*,'(5x,a,f6.2,a)') 'relative wind azimuth +/-', WIND_DIR-SENSOR_AZI_P, ' deg'
  WRITE(*,'(/5x,a)') 'ch      Stokes I        Stokes Q        Stokes U        Stokes V'
  DO l = 1, MIN(n_Channels,8)
    WRITE(*,'(5x,i2,4(2x,es14.6))') l, RTSolution(l,1)%Stokes(1), RTSolution(l,1)%Stokes(2), &
                                       RTSolution(l,1)%Stokes(3), RTSolution(l,1)%Stokes(4)
  END DO

  WRITE(*,'(/5x,a,es12.4)')      'intensity scale                = ', scale
  WRITE(*,'(5x,a,es12.4,a,l1)')  'max |Stokes U|                 = ', max_U,  '   pass = ', ok_signal
  WRITE(*,'(5x,a,es12.4)')       'max |Stokes V|                 = ', max_V
  WRITE(*,'(5x,a,es12.4,a,l1)')  'Stokes I,Q even under mirror   = ', d_even, '   pass = ', ok_even
  WRITE(*,'(5x,a,es12.4,a,l1)')  'Stokes U,V odd  under mirror   = ', d_odd,  '   pass = ', ok_odd

  all_ok = ok_signal .AND. ok_even .AND. ok_odd

  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a)') 'PASS: U and V reach RTSolution%Stokes, and the solver'
    WRITE(*,'(5x,a/)') '      preserves the meridional Stokes frame.'
    STOP 0
  ELSE
    IF ( .NOT. ok_signal ) THEN
      WRITE(*,'(/5x,a)') 'FAIL: RTSolution%Stokes(3:4) are zero on every channel, so the'
      WRITE(*,'(5x,a)')  '      polarimetric output is annihilated before it is reported.'
    ELSE
      WRITE(*,'(/5x,a)') 'FAIL: the solver does not preserve the surface Stokes frame.'
    END IF
    WRITE(*,'(a)') ''
    STOP 1
  END IF

CONTAINS

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_StokesOutput
