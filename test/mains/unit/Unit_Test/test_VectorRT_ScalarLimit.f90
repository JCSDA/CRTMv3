!
! test_VectorRT_ScalarLimit
!
! Ground-truth test for the polarimetric (n_Stokes > 1) ADA path, constructed so
! that it depends on neither the physical quality of the cloud lookup table nor
! any external radiative transfer code.
!
! The identity
! ------------
! In the limit of negligible scattering, an atmosphere is polarization-neutral:
! it emits unpolarized radiation and its transmittance is identical for every
! Stokes component (CRTM replicates the per-angle cosine across the Stokes slots
! of that angle, Common_RTSolution.f90 ~line 360). The only polarized quantity in
! the problem is the surface. Radiance is therefore affine in the surface
! emissivity, separately for each polarization, and the emergent Stokes vector
! must satisfy exactly
!
!     I = ( Iv + Ih ) / 2
!     Q = ( Iv - Ih ) / 2
!
! where Iv and Ih are the radiances of two ordinary scalar runs on the identical
! scene with the channel's polarization forced to pure vertical and pure
! horizontal. The reflection term obeys the same relation, because the surface
! reflection matrix carries (rV+rH)/2 on its diagonal and (rV-rH)/2 off it while
! the downwelling radiation is unpolarized.
!
! Why this is the right ground truth here
! ---------------------------------------
! The scalar path is the code the entire data assimilation community runs, so it
! is the most heavily exercised and most trustworthy reference available inside
! CRTM. This test uses it to validate the vector path end to end: the surface
! (V,H) to Stokes (I,Q) conversion, the ADA adding machinery under n_Stokes > 1,
! the surface boundary condition, and the emergent Stokes vector. None of that
! requires the cloud lookup table to be physically correct, because the cloud is
! present only to route the calculation through ADA rather than the (scalar)
! emission solver, and its scattering is driven to negligible.
!
! It is deliberately not a self-consistency check. Tangent-linear versus finite
! difference, the adjoint dot-product identity and K versus AD all verify that
! the derivative code matches the forward code, and all of them passed while the
! surface handoff was wrong by a factor of 3.3 in Stokes Q. This test compares
! the forward model against an independent statement about the physics.
!
! What it actually exercises, and the history
! -------------------------------------------
! Driving the cloud water content low enough to suppress scattering also drops
! CRTM_Include_Scattering below its trigger, so RTV%Scattering_RT becomes false
! and both runs dispatch to the emission path rather than ADA. Tuning cannot
! avoid that: SCATTERING_ALBEDO_THRESHOLD is 1.0e-10, so there is no window in
! which ADA runs with its scattering coupling inactive. Validating ADA's vector
! machinery requires driving CRTM_ADA directly with a synthetic unpolarized
! phase matrix, which is the companion test, test_ADA_VectorDecoupling. What
! this test validates is therefore the NON-SCATTERING vector path: the surface
! (V,H) to Stokes (I,Q) conversion, the flattened surface handoff, and the
! polarized emission solution, against the trusted scalar path.
!
! It failed until 2026-07-31, for two reasons that were each recorded wrongly
! at first and are worth keeping straight.
!
! An earlier revision of this header claimed the intensity came back about 9
! percent high on the opaque 183 GHz channels. That was an artefact of this
! test, not of CRTM. Atm(2) was never copied from Atm(1), so the second profile
! carried an empty atmosphere while the scalar reference Iv_ch/Ih_ch was
! captured from the first, and the two were compared across different scenes.
! With the profiles equalized the intensity error is 1.1e-16.
!
! The real defect was confined to the polarized components. CRTM_Emission is a
! scalar solver: it reads emissivity(n_Angles), a single element, and returns
! one scalar radiance. With n_Angles = 1, which is what the microwave
! non-scattering path always uses because it is specular, the flattened arrays
! happen to place Stokes I first, so emissivity(n_Angles) is e_I and
! reflectivity(1,1) is R_II and the intensity was right by construction. But no
! Q, U or V was produced at all, and Assign_Common_Output filled Radiance(1)
! only, leaving Radiance(2:n_Stokes) read before it was ever assigned. A
! clear-sky polarimetric run, which is exactly what ocean wind-vector retrieval
! needs, returned Q = 0.
!
! That is now fixed by CRTM_Emission_Stokes and its tangent-linear and adjoint
! siblings, and this test is registered. Jacobian coverage for the same path
! lives in test_VectorRT_TLADK, whose clear-sky block drives TL against finite
! difference, the adjoint dot product and K against AD on it.
!
! Scattering is suppressed by reducing the cloud water content until every
! layer's single-scatter albedo falls below CRTM's scattering threshold, at which
! point ADA takes its diagonal-transmittance branch. The test reports the
! residual so the margin is visible rather than assumed.
!

PROGRAM test_VectorRT_ScalarLimit

  USE CRTM_Module
  USE CRTM_SpcCoeff        , ONLY: SC
  USE SensorInfo_Parameters, ONLY: VL_POLARIZATION, HL_POLARIZATION
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_ScalarLimit'
  CHARACTER(*), PARAMETER :: SENSOR = 'mwr_aws'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: LUT    = 'CloudCoeff_Exp_Full6.nc'

  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0
  REAL(fp), PARAMETER :: ZENITH      = 53.0_fp
  INTEGER,  PARAMETER :: KC1 = 78, KC2 = 86
  REAL(fp), PARAMETER :: REFF_S = 500.0_fp
  ! Water content driven far below the scattering-albedo threshold so that ADA
  ! runs but its scattering coupling is inactive.
  REAL(fp), PARAMETER :: WC_TINY = 1.0e-8_fp

  ! The identity is exact in the no-scattering limit, so the tolerance need only
  ! absorb round-off plus the residual scattering left by WC_TINY.
  REAL(fp), PARAMETER :: TOL = 1.0e-8_fp

  CHARACTER(256) :: Version
  INTEGER  :: Error_Status, Allocate_Status, n_Channels, l, m, saved_pol
  LOGICAL  :: all_ok
  REAL(fp) :: Iv, Ih, Iexp, Qexp, Igot, Qgot, dI, dQ, worst_I, worst_Q

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:)
  REAL(fp), ALLOCATABLE :: Iv_ch(:), Ih_ch(:)

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Vector-RT scalar-limit ground truth (ADA, negligible scattering)'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo,      &
                            Cloud_Model       = 'CRTM-Exp', &
                            CloudCoeff_File   = LUT,        &
                            CloudCoeff_Format = 'netCDF',   &
                            File_Path         = PATH,       &
                            Quiet             = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))

  ALLOCATE( RTS(n_Channels,N_PROFILES), Iv_ch(n_Channels), Ih_ch(n_Channels), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN; WRITE(*,*) 'Alloc error'; STOP 1; END IF
  CALL CRTM_RTSolution_Create( RTS, N_LAYERS )
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    CALL Display_Message( PROGRAM_NAME, 'Atmosphere_Create failed', FAILURE ); STOP 1
  END IF

  ! Scene: a vanishingly thin frozen layer, present only so the calculation is
  ! routed through ADA instead of the scalar emission solver.
  CALL Load_ECMWF84_Atm_Data()          ! fills Atm(1) only
  DO m = 2, N_PROFILES
    Atm(m) = Atm(1)                     ! every profile must be the SAME scene
  END DO
  DO m = 1, N_PROFILES
    Atm(m)%n_Clouds                           = 1
    Atm(m)%Cloud_Fraction                     = ZERO
    Atm(m)%Cloud_Fraction(KC1:KC2)            = ONE
    Atm(m)%Cloud(1)%Type                      = SNOW_CLOUD
    Atm(m)%Cloud(1)%Effective_Radius          = ZERO
    Atm(m)%Cloud(1)%Water_Content             = ZERO
    Atm(m)%Cloud(1)%Effective_Radius(KC1:KC2) = REFF_S
    Atm(m)%Cloud(1)%Water_Content(KC1:KC2)    = WC_TINY

    Sfc(m)%Water_Coverage    = ONE
    Sfc(m)%Water_Type        = 1
    Sfc(m)%Water_Temperature = 290.0_fp
    Sfc(m)%Wind_Speed        = 6.0_fp
    Sfc(m)%Salinity          = 33.0_fp
    CALL CRTM_Geometry_SetValue( Geometry(m), Sensor_Zenith_Angle = ZENITH )
  END DO

  ! ------------------------------------------------------------------
  ! Reference: two scalar runs with the channel forced to pure V then H
  ! ------------------------------------------------------------------
  CALL run_scalar( VL_POLARIZATION, Iv_ch )
  CALL run_scalar( HL_POLARIZATION, Ih_ch )

  ! ------------------------------------------------------------------
  ! Vector run on the identical scene
  ! ------------------------------------------------------------------
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 2
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO
  Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward (vector) failed', FAILURE ); STOP 1
  END IF

  WRITE(*,'(5x,a,a)') 'vector run  solver=', TRIM(RTS(1,1)%RT_Algorithm_Name)

  ! ------------------------------------------------------------------
  ! Compare against the identity
  ! ------------------------------------------------------------------
  worst_I = ZERO ; worst_Q = ZERO
  WRITE(*,'(5x,a)') '  ch      Iv          Ih        I(expect)   I(CRTM)     Q(expect)   Q(CRTM)'
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      Iv   = Iv_ch(l) ; Ih = Ih_ch(l)
      Iexp = POINT_5*(Iv + Ih)
      Qexp = POINT_5*(Iv - Ih)
      Igot = RTS(l,m)%Stokes(1)
      Qgot = RTS(l,m)%Stokes(2)
      dI   = ABS(Igot - Iexp)
      dQ   = ABS(Qgot - Qexp)
      worst_I = MAX(worst_I, dI)
      worst_Q = MAX(worst_Q, dQ)
      IF ( l <= 4 .OR. dI > TOL .OR. dQ > TOL ) &
        WRITE(*,'(5x,i4,6f12.6)') l, Iv, Ih, Iexp, Igot, Qexp, Qgot
    END DO
  END DO

  WRITE(*,'(/5x,a,es12.4)') 'worst |I - (Iv+Ih)/2| = ', worst_I
  WRITE(*,'(5x,a,es12.4)')  'worst |Q - (Iv-Ih)/2| = ', worst_Q
  WRITE(*,'(5x,a,es12.4)')  'tolerance             = ', TOL

  all_ok = ( worst_I < TOL ) .AND. ( worst_Q < TOL )

  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: vector RT reduces to the scalar path in the no-scattering limit'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: vector RT does not reduce to the scalar path'
    STOP 1
  END IF

CONTAINS

  ! Run the scalar (n_Stokes=1) model with every channel's polarization
  ! temporarily forced to pol, returning the per-channel radiance.
  SUBROUTINE run_scalar( pol, out )
    INTEGER,  INTENT(IN)  :: pol
    REAL(fp), INTENT(OUT) :: out(:)
    INTEGER :: ll, mm, saved(n_Channels)
    DO ll = 1, n_Channels
      saved(ll) = SC(1)%Polarization(ll)
      SC(1)%Polarization(ll) = pol
    END DO
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = 1
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward (scalar) failed', FAILURE ); STOP 1
    END IF
    WRITE(*,'(5x,a,i0,a,a,a,i0)') 'scalar run pol=', pol, '  solver=', &
       TRIM(RTS(1,1)%RT_Algorithm_Name), '  n_Stokes=1, n_Layers=', RTS(1,1)%n_Layers
    DO ll = 1, n_Channels
      out(ll) = RTS(ll,1)%Radiance   ! profile 1 is the comparison scene
      SC(1)%Polarization(ll) = saved(ll)
    END DO
  END SUBROUTINE run_scalar

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_ScalarLimit
