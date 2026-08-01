!
! test_VectorRT_StokesSign
!
! Pins the adopted sign convention of the third and fourth Stokes components
! of the microwave water surface, independently for the FASTEM and PARMIO
! backends.
!
! Why this test exists
! --------------------
! The azimuthal emissivity expansion puts V and H on cosine harmonics and the
! third and fourth Stokes components on sine harmonics (Liu et al., FASTEM-4
! validation, NWPSAF-MO-VS-045, equations 2a-2d; implemented in
! Azimuth_Emissivity_Module.f90:139-142 and PARMIO_Azimuth_Module.f90:88-91).
! Cosine is even in the relative azimuth and sine is odd. A global sign error
! in U therefore cancels out of I and Q entirely and is invisible to every
! test that checks them.
!
! It is invisible to the polarimetric tests too. test_VectorRT_SurfaceFrame
! asserts that U and V4 reach the solver unchanged, that I and Q are even
! under phi -> -phi, and that U and V4 are odd. Negating U globally preserves
! all three: pass-through still holds because the reference is negated with
! it, oddness is preserved under negation, and the magnitude floor is
! unchanged. The self-consistency instruments are no help either, since
! TL against finite difference, the adjoint dot product and K against AD all
! compare the model to itself.
!
! So nothing in the suite fails if U's sign flips. That is what this test
! fixes.
!
! What it asserts
! ---------------
! With the relative azimuth built as phi = Wind_Direction - Sensor_Azimuth
! (the convention defined in CRTM_MW_Water_SfcOptics.f90), at phi = +90 the
! first harmonic is at its extremum and the second vanishes, so U reduces to
! the first-harmonic amplitude alone and its sign is read directly:
!
!   (a) U and V4 vanish at phi = 0 and phi = 180. A sine expansion carries no
!       constant term, so a leak here means a cosine term has been mixed in.
!   (b) U(+90) = -U(-90) and V4(+90) = -V4(-90), the oddness that identifies
!       the components as the ones changing handedness under reflection.
!   (c) the SIGN of U(+90) and V4(+90) matches the adopted convention, per
!       backend. This is the assertion the rest of the suite cannot make.
!   (d) both are above a non-degeneracy floor, so (a) and (b) cannot pass on
!       zeros.
!
! Scope and honesty about what this proves
! ----------------------------------------
! This test pins the convention AS ADOPTED in
! docs/design/polarimetric_conventions.md. It is not evidence that the sign
! is correct against nature.
!
! The FASTEM-4 report defines the harmonic form but never defines the origin
! or sense of its relative azimuth, and no accessible RTTOV or NWP SAF
! document states it either. Whether CRTM's phi origin matches the one the
! coefficients were regressed under remains open, and closing it needs an
! external reference: an RTTOV run at nonzero wind direction, or WindSat's
! published upwind/downwind harmonic amplitudes.
!
! What this test does guarantee is that the convention cannot change
! silently. If the external check later shows the adopted sign is wrong, this
! test is the thing that has to be edited deliberately to change it, and the
! edit is reviewable.
!

PROGRAM test_VectorRT_StokesSign

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

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_StokesSign'
  ! amsua_n19 channel 1 (23.8 GHz) routes to FASTEM; mwr_aws channel 16
  ! (325 GHz) is above PARMIO_FREQ_THRESHOLD and routes to PARMIO, which
  ! carries its own independently fitted four-Stokes azimuth model.
  CHARACTER(*), PARAMETER :: SENSORS(2) = (/ 'amsua_n19', 'mwr_aws  ' /)
  CHARACTER(*), PARAMETER :: PATH       = 'testinput/'
  ! FASTEM6, the CRTM default, has no third or fourth Stokes azimuth model
  ! and returns both as identically zero. FASTEM4 carries all four.
  ! Selection is by scheme name, not filename.
  CHARACTER(*), PARAMETER :: MWWATER_SCHEME = 'FASTEM4'

  ! Scene: open ocean, brisk wind so the azimuthal signal is well clear of
  ! round-off.
  REAL(fp), PARAMETER :: WIND_SPEED = 12.0_fp     ! m/s
  REAL(fp), PARAMETER :: WATER_TEMP = 285.0_fp    ! K
  REAL(fp), PARAMETER :: SALINITY   = 33.0_fp     ! ppmv
  REAL(fp), PARAMETER :: ZENITH     = 45.0_fp     ! deg

  ! phi = WIND_DIR - SENSOR_AZI. At phi = +/-90 the first harmonic is at its
  ! extremum and sin(2 phi) = 0, so U is the first-harmonic amplitude alone.
  REAL(fp), PARAMETER :: WIND_DIR    =  90.0_fp
  REAL(fp), PARAMETER :: AZI_PLUS90  =   0.0_fp   ! -> phi = +90
  REAL(fp), PARAMETER :: AZI_MINUS90 = 180.0_fp   ! -> phi = -90
  REAL(fp), PARAMETER :: AZI_ZERO    =  90.0_fp   ! -> phi =   0
  REAL(fp), PARAMETER :: AZI_ONE80   = 270.0_fp   ! -> phi = -180

  ! ------------------------------------------------------------------
  ! THE ADOPTED CONVENTION.
  ! Sign of the third and fourth Stokes surface emissivity at phi = +90,
  ! per backend. Changing either value changes CRTM's polarimetric sign
  ! convention and must be a deliberate, documented decision. See
  ! docs/design/polarimetric_conventions.md.
  ! ------------------------------------------------------------------
  REAL(fp), PARAMETER :: SIGN_U_FASTEM  = -1.0_fp
  REAL(fp), PARAMETER :: SIGN_V4_FASTEM = -1.0_fp
  REAL(fp), PARAMETER :: SIGN_U_PARMIO  = -1.0_fp
  REAL(fp), PARAMETER :: SIGN_V4_PARMIO = -1.0_fp

  INTEGER , PARAMETER :: N_ANGLES = 1

  ! phi = 0 and phi = 180 give sin(m phi) = 0 exactly for m = 1 and 2, so the
  ! null assertions need only round-off slack.
  REAL(fp), PARAMETER :: TOL          = 1.0e-12_fp
  REAL(fp), PARAMETER :: SIGNAL_FLOOR = 1.0e-8_fp

  CHARACTER(256) :: Version
  INTEGER :: Error_Status
  LOGICAL :: ok_fastem, ok_parmio, all_ok

  TYPE(CRTM_ChannelInfo_type)  :: ChannelInfo(2)
  TYPE(CRTM_Surface_type)      :: Sfc
  TYPE(CRTM_SfcOptics_type)    :: SfcOptics_ref

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT third/fourth Stokes sign convention'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( SENSORS, ChannelInfo,                 &
                            File_Path           = PATH,           &
                            MWwaterCoeff_Scheme = MWWATER_SCHEME, &
                            Quiet               = .TRUE.          )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF

  Sfc%Water_Coverage    = ONE
  Sfc%Land_Coverage     = ZERO
  Sfc%Snow_Coverage     = ZERO
  Sfc%Ice_Coverage      = ZERO
  Sfc%Water_Type        = 1
  Sfc%Water_Temperature = WATER_TEMP
  Sfc%Wind_Speed        = WIND_SPEED
  Sfc%Wind_Direction    = WIND_DIR
  Sfc%Salinity          = SALINITY

  CALL CRTM_SfcOptics_Create( SfcOptics_ref, N_ANGLES, MAX_N_STOKES )
  IF ( .NOT. CRTM_SfcOptics_Associated(SfcOptics_ref) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'SfcOptics_Create failed', FAILURE ); STOP 1
  END IF
  SfcOptics_ref%Angle(1)      = ZENITH
  SfcOptics_ref%Weight(1)     = ONE
  SfcOptics_ref%Index_Sat_Ang = 1
  SfcOptics_ref%n_Angles      = N_ANGLES

  CALL check_backend( 1,  1, 'FASTEM (amsua_n19 ch1, 23.8 GHz)', &
                      SIGN_U_FASTEM, SIGN_V4_FASTEM, ok_fastem )
  CALL check_backend( 2, 16, 'PARMIO (mwr_aws ch16, 325 GHz)  ', &
                      SIGN_U_PARMIO, SIGN_V4_PARMIO, ok_parmio )

  all_ok = ok_fastem .AND. ok_parmio

  WRITE(*,'(/5x,a)') '=================================================='
  IF ( all_ok ) THEN
    WRITE(*,'(5x,a)') 'RESULT: PASS - adopted Stokes sign convention held'
  ELSE
    WRITE(*,'(5x,a)') 'RESULT: FAIL - Stokes sign convention violated'
  END IF
  WRITE(*,'(5x,a/)') '=================================================='

  CALL CRTM_SfcOptics_Destroy( SfcOptics_ref )
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    STOP 0
  ELSE
    STOP 1
  END IF

CONTAINS

  ! Evaluate the surface optics at one relative azimuth and return U and V4.
  SUBROUTINE eval_at( sidx, chan, sensor_azi, eU, eV4 )
    INTEGER     , INTENT(IN)  :: sidx, chan
    REAL(fp)    , INTENT(IN)  :: sensor_azi
    REAL(fp)    , INTENT(OUT) :: eU, eV4
    TYPE(CRTM_GeometryInfo_type) :: gInfo
    TYPE(CRTM_SfcOptics_type)    :: SfcOptics
    TYPE(iVar_type)              :: iVar
    INTEGER                      :: err

    CALL CRTM_GeometryInfo_SetValue( gInfo, Sensor_Zenith_Angle  = ZENITH, &
                                            Sensor_Azimuth_Angle = sensor_azi )
    CALL CRTM_GeometryInfo_Compute( gInfo )

    SfcOptics = SfcOptics_ref
    SfcOptics%n_Stokes = 4
    err = CRTM_Compute_SfcOptics( Sfc, gInfo, sidx, chan, SfcOptics, iVar )
    IF ( err /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'Compute_SfcOptics failed', FAILURE ); STOP 1
    END IF
    eU  = SfcOptics%Emissivity(1,3)
    eV4 = SfcOptics%Emissivity(1,4)
  END SUBROUTINE eval_at


  SUBROUTINE check_backend( sidx, chan, label, want_U, want_V4, ok )
    INTEGER     , INTENT(IN)  :: sidx, chan
    CHARACTER(*), INTENT(IN)  :: label
    REAL(fp)    , INTENT(IN)  :: want_U, want_V4
    LOGICAL     , INTENT(OUT) :: ok

    REAL(fp) :: eU_p, eV4_p, eU_m, eV4_m, eU_0, eV4_0, eU_180, eV4_180
    LOGICAL  :: ok_null, ok_odd, ok_sign, ok_signal

    CALL eval_at( sidx, chan, AZI_PLUS90 , eU_p  , eV4_p   )
    CALL eval_at( sidx, chan, AZI_MINUS90, eU_m  , eV4_m   )
    CALL eval_at( sidx, chan, AZI_ZERO   , eU_0  , eV4_0   )
    CALL eval_at( sidx, chan, AZI_ONE80  , eU_180, eV4_180 )

    ! (a) no constant term: a pure sine expansion vanishes at 0 and 180
    ok_null = ( ABS(eU_0)   < TOL ) .AND. ( ABS(eV4_0)   < TOL ) .AND. &
              ( ABS(eU_180) < TOL ) .AND. ( ABS(eV4_180) < TOL )
    ! (b) odd under phi -> -phi
    ok_odd  = ( ABS(eU_p + eU_m) < TOL ) .AND. ( ABS(eV4_p + eV4_m) < TOL )
    ! (c) the adopted sign
    ok_sign = ( SIGN(ONE, eU_p)  == SIGN(ONE, want_U) ) .AND. &
              ( SIGN(ONE, eV4_p) == SIGN(ONE, want_V4) )
    ! (d) not vacuous
    ok_signal = ( ABS(eU_p) > SIGNAL_FLOOR ) .AND. ( ABS(eV4_p) > SIGNAL_FLOOR )

    WRITE(*,'(/5x,a)') '--- backend: '//label//' ---'
    WRITE(*,'(5x,a,f6.2,a,f6.2,a)') 'ocean, zenith ', ZENITH, &
        ' deg, wind ', WIND_SPEED, ' m/s'
    WRITE(*,'(5x,a,es14.6,a,es14.6)') 'phi = +90   U = ', eU_p  , '   V4 = ', eV4_p
    WRITE(*,'(5x,a,es14.6,a,es14.6)') 'phi = -90   U = ', eU_m  , '   V4 = ', eV4_m
    WRITE(*,'(5x,a,es14.6,a,es14.6)') 'phi =   0   U = ', eU_0  , '   V4 = ', eV4_0
    WRITE(*,'(5x,a,es14.6,a,es14.6)') 'phi = 180   U = ', eU_180, '   V4 = ', eV4_180
    WRITE(*,'(5x,a,f5.1,a,f5.1)')     'adopted sign at +90:  U = ', want_U, &
                                      '   V4 = ', want_V4

    WRITE(*,'(/5x,a,l1)') 'vanish at phi = 0 and 180 ..................  pass = ', ok_null
    WRITE(*,'(5x,a,l1)')  'odd under phi -> -phi .....................  pass = ', ok_odd
    WRITE(*,'(5x,a,l1)')  'sign matches adopted convention ...........  pass = ', ok_sign
    WRITE(*,'(5x,a,l1)')  'above non-degeneracy floor ................  pass = ', ok_signal

    ok = ok_null .AND. ok_odd .AND. ok_sign .AND. ok_signal

  END SUBROUTINE check_backend

END PROGRAM test_VectorRT_StokesSign
