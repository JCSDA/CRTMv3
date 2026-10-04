!
! test_VectorRT_Physics
!
! Physical invariants of the emergent Stokes vector, asserted without any cloud
! lookup table and without any external reference radiance.
!
! Why this test is possible at all
! -------------------------------
! The polarimetric lookup tables are known to be unsuitable for full-Stokes
! work, so any check that depends on their content cannot separate a code
! defect from a data defect. This test avoids the question entirely by running
! clear sky. With no scattering the atmosphere is polarization neutral, so
! there is no cloud optics in the problem: the emergent Stokes vector is fixed
! by the surface model, the gas absorption and the radiative transfer itself.
! Every assertion below is then a statement about the machinery.
!
! What is asserted
! ----------------
! 1. TRUNCATION LADDER. Without scattering the Stokes components do not couple
!    to one another anywhere: each is a surface boundary value transported
!    upward with no source. Running the same scene at n_Stokes = 2, 3 and 4
!    must therefore return bit-identical I and Q, and bit-identical U between
!    3 and 4. Anything that leaks between components, or any dimension-
!    dependent indexing error, breaks this. It also exercises n_Stokes = 3,
!    which nothing else does: the solver guards U with n_Stokes > 2 and V with
!    n_Stokes == 4, so the three-component truncation is a distinct path.
!
! 2. ODD-HARMONIC DEGENERACY. The surface third and fourth Stokes components
!    are built from sin(m*phi) harmonics of the relative wind azimuth. At
!    phi = 0 and phi = 180 degrees every one of those vanishes, so U and V must
!    come back exactly zero. This is the check that catches a U or V that is
!    not actually the azimuthal signal but some other quantity leaking into
!    those slots: such a leak would be non-zero here.
!
! 3. POLARIZATION BOUND. Physically realisable radiation satisfies
!    I^2 >= Q^2 + U^2 + V^2. This is asserted on the emergent radiance at a
!    relative azimuth where the polarized components are genuinely non-zero.
!
! 4. POSITIVITY. Stokes I is a total intensity and must be positive.
!
! FASTEM4 is loaded because the FASTEM6 default carries no third or fourth
! Stokes azimuth model, which would make assertions 2 and 3 vacuous.
!

PROGRAM test_VectorRT_Physics

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_Physics'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'amsua_n19'
  CHARACTER(*), PARAMETER :: MWWATER_SCHEME = 'FASTEM4'

  INTEGER,  PARAMETER :: N_PROFILES  = 2     ! the ECMWF84 loader fills atm(1) and atm(2)
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  REAL(fp), PARAMETER :: WIND_DIR    = 100.0_fp
  REAL(fp), PARAMETER :: AIRCRAFT_P  = 300.0_fp   ! hPa

  ! Machine-precision: the ladder and the degeneracy are exact statements.
  REAL(fp), PARAMETER :: TOL         = 1.0e-14_fp
  REAL(fp), PARAMETER :: SIGNAL_FLOOR = 1.0e-12_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, Allocate_Status, n_Channels, l, m
  LOGICAL :: ok_ladder, ok_deg, ok_bound, ok_pos, ok_signal, all_ok
  LOGICAL :: ok_air_bound, ok_air_diff
  REAL(fp) :: d_ladder, d_deg, worst_bound, min_I, max_pol
  REAL(fp) :: SA(4,64,2), air_bound, air_diff
  REAL(fp) :: S2(4,64,2), S3(4,64,2), S4(4,64,2)
  REAL(fp) :: pol2

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT emergent Stokes physical invariants (clear sky, no cloud LUT)'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo,             &
                            File_Path           = PATH,            &
                            MWwaterCoeff_Scheme = MWWATER_SCHEME,  &
                            Quiet               = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels > 64 ) THEN
    CALL Display_Message( PROGRAM_NAME, 'raise the 64-channel scratch bound', FAILURE ); STOP 1
  END IF

  ALLOCATE( RTSolution(n_Channels,N_PROFILES), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Atmosphere_Create failed', FAILURE ); STOP 1
  END IF

  CALL Load_ECMWF84_Atm_Data()          ! fills Atm(1) and Atm(2)

  DO m = 1, N_PROFILES
    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 12.0_fp
    Sfc(m)%Wind_Direction    = WIND_DIR
    Sfc(m)%Salinity          = 33.0_fp
  END DO

  ! ------------------------------------------------------------------
  ! 1 and 3 and 4: relative azimuth 60 degrees, where U and V are real
  ! ------------------------------------------------------------------
  CALL set_azimuth( WIND_DIR - 60.0_fp )
  CALL run( 2, S2 )
  CALL run( 3, S3 )
  CALL run( 4, S4 )

  ! Truncation ladder
  d_ladder = ZERO
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      d_ladder = MAX( d_ladder, ABS(S2(1,l,m)-S3(1,l,m)), ABS(S2(2,l,m)-S3(2,l,m)) )
      d_ladder = MAX( d_ladder, ABS(S3(1,l,m)-S4(1,l,m)), ABS(S3(2,l,m)-S4(2,l,m)) )
      d_ladder = MAX( d_ladder, ABS(S3(3,l,m)-S4(3,l,m)) )
    END DO
  END DO
  ok_ladder = ( d_ladder < TOL )

  ! Polarization bound and positivity, on the full four-component run
  ! Start below any achievable value so the reported number is the true worst
  ! margin, not the initialiser. The margin is negative when the bound holds.
  worst_bound = -HUGE(ONE)
  min_I       = HUGE(ONE)
  max_pol     = ZERO
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      pol2 = S4(2,l,m)**2 + S4(3,l,m)**2 + S4(4,l,m)**2
      ! positive when the bound is violated
      worst_bound = MAX( worst_bound, pol2 - S4(1,l,m)**2 )
      min_I       = MIN( min_I, S4(1,l,m) )
      max_pol     = MAX( max_pol, SQRT(pol2) )
    END DO
  END DO
  ok_bound = ( worst_bound <= ZERO )
  ok_pos   = ( min_I > ZERO )
  ! Guard against the bound holding only because U and V are zero
  ok_signal = ( max_pol > SIGNAL_FLOOR )

  ! ------------------------------------------------------------------
  ! 2: relative azimuth 0 and 180, where every odd harmonic vanishes
  ! ------------------------------------------------------------------
  d_deg = ZERO
  CALL set_azimuth( WIND_DIR )              ! relative azimuth 0
  CALL run( 4, S4 )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      d_deg = MAX( d_deg, ABS(S4(3,l,m)), ABS(S4(4,l,m)) )
    END DO
  END DO
  CALL set_azimuth( WIND_DIR - 180.0_fp )   ! relative azimuth 180
  CALL run( 4, S4 )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      d_deg = MAX( d_deg, ABS(S4(3,l,m)), ABS(S4(4,l,m)) )
    END DO
  END DO
  ok_deg = ( d_deg < TOL )

  ! ------------------------------------------------------------------
  ! 5: aircraft observer. CRTM_Emission_Stokes transports the polarized
  !    components to the observer level, not unconditionally to the top of
  !    the atmosphere, and that branch is otherwise never executed.
  ! ------------------------------------------------------------------
  CALL set_azimuth( WIND_DIR - 60.0_fp )
  CALL run( 4, S4 )                          ! reference: top of atmosphere
  DO m = 1, N_PROFILES
    Options(m)%Aircraft_Pressure = AIRCRAFT_P
  END DO
  CALL run( 4, SA )
  DO m = 1, N_PROFILES
    Options(m)%Aircraft_Pressure = ZERO      ! restore
  END DO
  air_bound = -HUGE(ONE)
  air_diff  = ZERO
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      air_bound = MAX( air_bound, SA(2,l,m)**2 + SA(3,l,m)**2 + SA(4,l,m)**2 - SA(1,l,m)**2 )
      ! The aircraft view must not simply reproduce the top-of-atmosphere one,
      ! otherwise the observer level is being ignored and the check is vacuous.
      air_diff = MAX( air_diff, ABS(SA(1,l,m) - S4(1,l,m)) )
    END DO
  END DO
  ok_air_bound = ( air_bound <= ZERO )
  ok_air_diff  = ( air_diff > SIGNAL_FLOOR )

  WRITE(*,'(5x,a,i0,a)') 'sensor '//SENSOR//', ', n_Channels, ' channels, clear sky over ocean'
  WRITE(*,'(/5x,a,es12.4,a,l1)') 'truncation ladder n_Stokes 2/3/4 = ', d_ladder, '   pass = ', ok_ladder
  WRITE(*,'(5x,a,es12.4,a,l1)')  'U,V zero at rel azimuth 0 and 180= ', d_deg,    '   pass = ', ok_deg
  WRITE(*,'(5x,a,es12.4,a,l1)')  'max (Q2+U2+V2 - I2), <0 is good  = ', worst_bound, '   pass = ', ok_bound
  WRITE(*,'(5x,a,es12.4,a,l1)')  'min Stokes I                     = ', min_I,    '   pass = ', ok_pos
  WRITE(*,'(5x,a,es12.4,a,l1)')  'max sqrt(Q2+U2+V2) (must be > 0) = ', max_pol,  '   pass = ', ok_signal

  WRITE(*,'(5x,a,es12.4,a,l1)')  'aircraft obs: Q2+U2+V2 - I2      = ', air_bound, '   pass = ', ok_air_bound
  WRITE(*,'(5x,a,es12.4,a,l1)')  'aircraft obs: differs from TOA   = ', air_diff,  '   pass = ', ok_air_diff

  all_ok = ok_ladder .AND. ok_deg .AND. ok_bound .AND. ok_pos .AND. ok_signal &
           .AND. ok_air_bound .AND. ok_air_diff

  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: emergent Stokes vector satisfies the physical invariants'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: emergent Stokes vector violates a physical invariant'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE set_azimuth( sensor_azi )
    REAL(fp), INTENT(IN) :: sensor_azi
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      CALL CRTM_Geometry_SetValue( Geometry(mm), Sensor_Zenith_Angle  = ZENITH, &
                                                 Sensor_Azimuth_Angle = sensor_azi )
    END DO
  END SUBROUTINE set_azimuth

  SUBROUTINE run( ns, S )
    INTEGER,  INTENT(IN)  :: ns
    REAL(fp), INTENT(OUT) :: S(4,64,2)
    INTEGER :: mm, ll, kk
    S = ZERO
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTSolution, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE(*,'(5x,a,i0,a)') 'CRTM_Forward failed at n_Stokes = ', ns, ''
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
    END IF
    DO mm = 1, N_PROFILES
      DO ll = 1, n_Channels
        DO kk = 1, ns
          S(kk,ll,mm) = RTSolution(ll,mm)%Stokes(kk)
        END DO
      END DO
    END DO
  END SUBROUTINE run

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_Physics
