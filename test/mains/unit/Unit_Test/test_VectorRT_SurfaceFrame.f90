!
! test_VectorRT_SurfaceFrame
!
! Pins the polarimetric reference frame of the microwave surface optics, and
! proves that the third and fourth Stokes components computed by the surface
! model survive the coverage aggregation into the vector solver input.
!
! The frame question
! ------------------
! A polarimetric Stokes vector is meaningless without the plane it is
! referred to. The vector radiative transfer solver carries (I,Q,U,V) in the
! meridional frame of each quadrature direction: the phase matrix is assembled
! from the generalized spherical functions R_l^m and T_l^m
! (Common_RTSolution.f90:2017-2019, RTV%Pplus and RTV%Pminus), which is the
! standard meridional-frame azimuthal Fourier expansion with the
! scattering-plane rotations folded in analytically. Consistent with that,
! no rotation matrix exists anywhere in src/RTSolution or src/SfcOptics.
!
! So if the surface models referred their (V,H,U,V) vector to any plane other
! than the meridional plane of the view direction, a rotation would be
! required at the handoff and its absence would be a defect. Critically, such
! a defect is invisible at nadir, where the meridional plane is degenerate,
! and invisible at zero relative azimuth, which is where the surface has no
! polarimetric signal at all.
!
! This test settles the question by measuring the symmetry that defines the
! frame. Reflect the scene through the vertical plane containing the view
! direction. The viewing geometry is unchanged and the relative wind azimuth
! maps phi -> -phi. A Stokes vector referred to a frame lying in that mirror
! plane must transform as
!
!     (I, Q, U, V)  ->  (I, Q, -U, -V)
!
! because I and Q are defined by intensities along axes that the reflection
! maps to themselves, while U and V change handedness. Observing exactly that
! even/odd split is what identifies the reference plane as the view plane,
! which for a plane-parallel atmosphere is the meridional plane. Any other
! reference plane would mix the components and break the split.
!
! What this test does
! -------------------
! Over open ocean at 45 degrees it calls CRTM_Compute_SfcOptics three times on
! the same scene:
!
!   1. scalar (n_Stokes=1) at relative azimuth +phi. The scalar branch writes
!      only component 1, so components 3 and 4 still hold what the surface
!      model itself produced. This is the reference: the raw FASTEM U and V.
!   2. vector (n_Stokes=4) at relative azimuth +phi.
!   3. vector (n_Stokes=4) at relative azimuth -phi.
!
! and asserts
!
!   (a) the vector path receives the surface model's U and V unchanged,
!       U_vector == U_raw and V_vector == V_raw;
!   (b) I and Q are even under phi -> -phi;
!   (c) U and V are odd under phi -> -phi;
!   (d) U and V are not identically zero, so (b) and (c) are not vacuous.
!
! Assertion (a) is what fails against the unfixed code, and it fails at
! exactly zero rather than by a tolerance margin: the microwave coverage
! aggregation in CRTM_SfcOptics copied only components 1 and 2 out of the
! surface model, into an array that was zero-initialised, so the solver was
! handed U = V = 0 no matter what FASTEM computed.
!
! Coefficient note
! ----------------
! This test loads FASTEM4 explicitly rather than taking the CRTM default.
! The default is FASTEM6, whose azimuth model (Kazumori,
! Azimuth_Emissivity_F6_Module.f90:187-188) parameterises the vertical and
! horizontal components only and returns the third and fourth Stokes
! components as identically zero. FASTEM4 and FASTEM5 use
! Azimuth_Emissivity_Module.f90:139-142, which carries all four. PARMIO, used
! at and above 200 GHz, also carries all four
! (PARMIO_Azimuth_Module.f90:89-91). A polarimetric run therefore has a real
! surface U and V only on the FASTEM4/5 or PARMIO backends, never on the
! shipped default.
!
! Like test_VectorRT_SurfaceBasis, this needs no cloud lookup table and no
! reference radiances. It is a statement about the surface handoff and its
! reference frame, which the self-consistency tests (TL against finite
! difference, adjoint dot product, K against AD) structurally cannot check.
!

PROGRAM test_VectorRT_SurfaceFrame

  ! -----------------
  ! Environment setup
  ! -----------------
  USE CRTM_Module
  USE CRTM_SfcOptics_Define   , ONLY: CRTM_SfcOptics_type      , &
                                      CRTM_SfcOptics_Create    , &
                                      CRTM_SfcOptics_Destroy   , &
                                      CRTM_SfcOptics_Associated
  USE CRTM_SfcOptics          , ONLY: CRTM_Compute_SfcOptics, iVar_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_SetValue
  USE CRTM_GeometryInfo       , ONLY: CRTM_GeometryInfo_Compute
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_SurfaceFrame'
  ! Two sensors: amsua_n19 channel 1 (23.8 GHz) exercises the FASTEM backend,
  ! mwr_aws channel 16 (325 GHz) is above PARMIO_FREQ_THRESHOLD so the MW-water
  ! dispatcher routes it to PARMIO instead. PARMIO carries its own independent
  ! four-Stokes azimuth model, and nothing exercised it: it is the default
  ! backend at and above 200 GHz, but every sensor that reaches that far sits
  ! on a water-vapour line, so its polarimetric surface never survives to the
  ! top of the atmosphere and a full radiative-transfer test cannot see it.
  ! Checking it here, at the surface interface, is the only way to reach it.
  CHARACTER(*), PARAMETER :: SENSORS(2)   = (/ 'amsua_n19', 'mwr_aws  ' /)
  CHARACTER(*), PARAMETER :: PATH         = 'testinput/'
  ! FASTEM4 carries the third and fourth Stokes azimuth harmonics; the FASTEM6
  ! default does not (see header). Selection is by scheme name, not filename:
  ! the file-based MWwaterCoeff load in CRTM_LifeCycle.f90 is commented out, so
  ! MWwaterCoeff_File does not choose the model. MWwaterCoeff_Scheme does.
  CHARACTER(*), PARAMETER :: MWWATER_SCHEME = 'FASTEM4'

  ! Scene: open ocean, brisk wind so the azimuthal signal is well above noise.
  REAL(fp), PARAMETER :: WIND_SPEED  = 12.0_fp     ! m/s
  REAL(fp), PARAMETER :: WATER_TEMP  = 285.0_fp    ! K
  REAL(fp), PARAMETER :: SALINITY    = 33.0_fp     ! ppmv
  REAL(fp), PARAMETER :: ZENITH      = 45.0_fp     ! deg, away from nadir so the
                                                   ! meridional plane is well defined
  ! Relative azimuth is built as WIND_DIR - SENSOR_AZI so that the mirrored
  ! case is the exact negation of the direct case, which makes the even/odd
  ! assertions exact rather than approximate.
  REAL(fp), PARAMETER :: WIND_DIR    = 100.0_fp    ! deg
  REAL(fp), PARAMETER :: SENSOR_AZI_P =  40.0_fp   ! -> relative azimuth = +60
  REAL(fp), PARAMETER :: SENSOR_AZI_M = 160.0_fp   ! -> relative azimuth = -60

  INTEGER , PARAMETER :: N_ANGLES    = 1
  INTEGER , PARAMETER :: CHANNEL     = 1           ! 23.8 GHz

  ! The identities are exact algebra plus an exactly-negated trigonometric
  ! argument, so only round-off slack is needed.
  REAL(fp), PARAMETER :: TOL = 1.0e-12_fp
  ! Non-degeneracy floor: far above round-off, far below any real signal.
  REAL(fp), PARAMETER :: SIGNAL_FLOOR = 1.0e-8_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status, n_Channels
  LOGICAL :: ok_pass_U, ok_pass_V, ok_even, ok_odd, ok_signal, all_ok
  LOGICAL :: ok_fastem, ok_parmio
  REAL(fp) :: eU_raw, eV_raw
  REAL(fp) :: eI_p, eQ_p, eU_p, eV_p
  REAL(fp) :: eI_m, eQ_m, eU_m, eV_m
  REAL(fp) :: d_passU, d_passV, d_even, d_odd

  TYPE(CRTM_ChannelInfo_type)  :: ChannelInfo(2)
  TYPE(CRTM_Surface_type)      :: Sfc
  TYPE(CRTM_GeometryInfo_type) :: gInfo_p, gInfo_m
  TYPE(CRTM_SfcOptics_type)    :: SfcOptics_s, SfcOptics_p, SfcOptics_m
  TYPE(iVar_type)              :: iVar

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT surface polarimetric frame verification'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  ! --------------
  ! Initialize
  ! --------------
  Error_Status = CRTM_Init( SENSORS, ChannelInfo,                  &
                            File_Path           = PATH,            &
                            MWwaterCoeff_Scheme = MWWATER_SCHEME,  &
                            Quiet               = .TRUE.           )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ! ------------------------
  ! Ocean surface + geometry
  ! ------------------------
  Sfc%Water_Coverage    = ONE
  Sfc%Land_Coverage     = ZERO
  Sfc%Snow_Coverage     = ZERO
  Sfc%Ice_Coverage      = ZERO
  Sfc%Water_Type        = 1                ! sea water
  Sfc%Water_Temperature = WATER_TEMP
  Sfc%Wind_Speed        = WIND_SPEED
  Sfc%Wind_Direction    = WIND_DIR
  Sfc%Salinity          = SALINITY

  CALL CRTM_GeometryInfo_SetValue( gInfo_p, Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = SENSOR_AZI_P )
  CALL CRTM_GeometryInfo_Compute( gInfo_p )
  CALL CRTM_GeometryInfo_SetValue( gInfo_m, Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = SENSOR_AZI_M )
  CALL CRTM_GeometryInfo_Compute( gInfo_m )

  CALL CRTM_SfcOptics_Create( SfcOptics_s, N_ANGLES, MAX_N_STOKES )
  CALL CRTM_SfcOptics_Create( SfcOptics_p, N_ANGLES, MAX_N_STOKES )
  CALL CRTM_SfcOptics_Create( SfcOptics_m, N_ANGLES, MAX_N_STOKES )
  IF ( .NOT. CRTM_SfcOptics_Associated(SfcOptics_s) .OR. &
       .NOT. CRTM_SfcOptics_Associated(SfcOptics_p) .OR. &
       .NOT. CRTM_SfcOptics_Associated(SfcOptics_m) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'SfcOptics_Create failed', FAILURE ); STOP 1
  END IF

  SfcOptics_s%Angle(1)      = ZENITH
  SfcOptics_s%Weight(1)     = ONE
  SfcOptics_s%Index_Sat_Ang = 1
  SfcOptics_s%n_Angles      = N_ANGLES
  SfcOptics_p = SfcOptics_s
  SfcOptics_m = SfcOptics_s

  ! Both microwave-water backends. amsua_n19 channel 1 is 23.8 GHz, below the
  ! PARMIO frequency threshold, so it uses FASTEM. mwr_aws channel 16 is
  ! 325 GHz, above it, so the dispatcher routes it to PARMIO, which has its own
  ! independent four-Stokes azimuth model.
  CALL check_backend( 1,  1, 'FASTEM (amsua_n19 ch1, 23.8 GHz)', ok_fastem )
  CALL check_backend( 2, 16, 'PARMIO (mwr_aws ch16, 325 GHz)  ', ok_parmio )
  all_ok = ok_fastem .AND. ok_parmio

  CALL CRTM_SfcOptics_Destroy( SfcOptics_s )
  CALL CRTM_SfcOptics_Destroy( SfcOptics_p )
  CALL CRTM_SfcOptics_Destroy( SfcOptics_m )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a)') 'PASS: surface Stokes frame is the view (meridional) plane,'
    WRITE(*,'(5x,a/)') '      and U, V reach the vector solver input.'
    STOP 0
  ELSE
    IF ( .NOT. ok_signal ) THEN
      WRITE(*,'(/5x,a)') 'FAIL: the surface third/fourth Stokes components are zero at the'
      WRITE(*,'(5x,a)')  '      solver input, so no polarimetric surface signal exists.'
    ELSE
      WRITE(*,'(/5x,a)') 'FAIL: surface polarimetric frame or handoff is not as asserted.'
    END IF
    WRITE(*,'(a)') ''
    STOP 1
  END IF

CONTAINS

  SUBROUTINE check_backend( sidx, chan, label, ok )
    INTEGER,      INTENT(IN)  :: sidx, chan
    CHARACTER(*), INTENT(IN)  :: label
    LOGICAL,      INTENT(OUT) :: ok

  ! ------------------------------------------------------------------
  ! 1. Raw surface-model U and V, read through the scalar path.
  !    The n_Stokes==1 branch writes component 1 only, so components 3
  !    and 4 still hold exactly what Compute_MW_Water_SfcOptics wrote.
  ! ------------------------------------------------------------------
  SfcOptics_s%n_Stokes = 1
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo_p, sidx, chan, SfcOptics_s, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (scalar) failed', FAILURE ); STOP 1
  END IF
  eU_raw = SfcOptics_s%Emissivity(1,3)
  eV_raw = SfcOptics_s%Emissivity(1,4)

  ! ------------------------------------------------
  ! 2. Vector path at relative azimuth +phi
  ! ------------------------------------------------
  SfcOptics_p%n_Stokes = 4
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo_p, sidx, chan, SfcOptics_p, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (vector +phi) failed', FAILURE ); STOP 1
  END IF
  eI_p = SfcOptics_p%Emissivity(1,1)
  eQ_p = SfcOptics_p%Emissivity(1,2)
  eU_p = SfcOptics_p%Emissivity(1,3)
  eV_p = SfcOptics_p%Emissivity(1,4)

  ! ------------------------------------------------
  ! 3. Vector path at relative azimuth -phi
  ! ------------------------------------------------
  SfcOptics_m%n_Stokes = 4
  Error_Status = CRTM_Compute_SfcOptics( Sfc, gInfo_m, sidx, chan, SfcOptics_m, iVar )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics (vector -phi) failed', FAILURE ); STOP 1
  END IF
  eI_m = SfcOptics_m%Emissivity(1,1)
  eQ_m = SfcOptics_m%Emissivity(1,2)
  eU_m = SfcOptics_m%Emissivity(1,3)
  eV_m = SfcOptics_m%Emissivity(1,4)

  ! ------------------------------------------------
  ! Assertions
  ! ------------------------------------------------
  ! (a) the surface model's U and V reach the solver input untouched
  d_passU = ABS( eU_p - eU_raw )
  d_passV = ABS( eV_p - eV_raw )
  ! (b) I and Q are even under the mirror reflection
  d_even  = MAX( ABS(eI_p - eI_m), ABS(eQ_p - eQ_m) )
  ! (c) U and V are odd under the mirror reflection
  d_odd   = MAX( ABS(eU_p + eU_m), ABS(eV_p + eV_m) )

  ok_pass_U = ( d_passU < TOL )
  ok_pass_V = ( d_passV < TOL )
  ok_even   = ( d_even  < TOL )
  ok_odd    = ( d_odd   < TOL )
  ! (d) non-degeneracy
  ok_signal = ( ABS(eU_p) > SIGNAL_FLOOR .AND. ABS(eV_p) > SIGNAL_FLOOR )

  WRITE(*,'(/5x,a)') '--- backend: '//label//' ---'
    WRITE(*,'(5x,a,i0,a,f6.2,a,f6.2,a)') 'Channel ', chan, ' at ', ZENITH, &
        ' deg, relative azimuth +/-', WIND_DIR-SENSOR_AZI_P, ' deg over ocean'
  WRITE(*,'(5x,a,es14.6,a,es14.6)') 'surface model  U = ', eU_raw, '   V = ', eV_raw
  WRITE(*,'(5x,a,es14.6,a,es14.6)') 'solver input   U = ', eU_p,   '   V = ', eV_p
  WRITE(*,'(5x,a,es14.6,a,es14.6)') 'mirrored       U = ', eU_m,   '   V = ', eV_m
  WRITE(*,'(5x,a,es14.6,a,es14.6)') 'mirrored       I = ', eI_m,   '   Q = ', eQ_m
  WRITE(*,'(5x,a,es14.6,a,es14.6)') 'direct         I = ', eI_p,   '   Q = ', eQ_p

  WRITE(*,'(/5x,a,es12.4,a,l1)') 'U reaches solver   |diff| = ', d_passU, '   pass = ', ok_pass_U
  WRITE(*,'(5x,a,es12.4,a,l1)')  'V reaches solver   |diff| = ', d_passV, '   pass = ', ok_pass_V
  WRITE(*,'(5x,a,es12.4,a,l1)')  'I,Q even in azimuth       = ', d_even,  '   pass = ', ok_even
  WRITE(*,'(5x,a,es12.4,a,l1)')  'U,V odd  in azimuth       = ', d_odd,   '   pass = ', ok_odd
  WRITE(*,'(5x,a,l1)')           'U,V above signal floor      ................  pass = ', ok_signal

  ok = ok_pass_U .AND. ok_pass_V .AND. ok_even .AND. ok_odd .AND. ok_signal


  END SUBROUTINE check_backend

END PROGRAM test_VectorRT_SurfaceFrame
