!
! test_VectorRT_VIS_Rayleigh
!
! First test of the VISIBLE vector radiative transfer path.
!
! Background
! ----------
! Until 2026-08-03 every test touching n_Stokes ran a microwave sensor, so
! the visible vector path had never been executed by anything in the suite.
! It was also documented as unsupported: the roadmap and the release notes
! both described polarimetric support as microwave-only, inferred from the
! surface models rather than measured.
!
! That was wrong. CRTM_MoleculeScatter fills the polarized Rayleigh
! expansion coefficients under an n_Stokes gate, with matching tangent
! linear and adjoint, and visible/ultraviolet is the only sensor class for
! which the azimuthal Fourier machinery (RTV%n_Azi) is switched on at all.
! A clear-sky visible run at n_Stokes = 4 therefore produces a large, real
! and correctly structured polarization signal, entirely molecular: the
! visible surface writes the first Stokes component only, so nothing
! polarized comes from it.
!
! This test exists so that capability cannot regress unnoticed, and so the
! claim "the visible produces real polarized radiance" is a measurement
! that reruns rather than a paragraph in a document.
!
! What is asserted
! ----------------
!   1. U vanishes at relative azimuth 0 and 180 and is large at 90. Every
!      odd azimuthal harmonic must vanish in the principal plane, so this
!      distinguishes Rayleigh polarization from something leaking into the
!      U slot. It is the single most diagnostic check here.
!   2. Q is largest in magnitude at relative azimuth 0 and is NEGATIVE
!      there. At Delta-phi = 0 the scattering plane is the meridian plane,
!      and Rayleigh scattering polarizes perpendicular to the scattering
!      plane, which is Q = I_v - I_h < 0. A sign error flips this.
!   3. The degree of polarization falls with wavelength, following the
!      lambda^-4 Rayleigh optical depth. Checked between the 0.47 um and
!      0.86 um bands.
!   4. V is zero to machine precision. Rayleigh scattering by molecules
!      produces no circular polarization.
!   5. The polarization bound I^2 >= Q^2+U^2+V^2 holds and is non-vacuous.
!   6. The vector and scalar intensities DIFFER. This is the scalar
!      approximation error in a Rayleigh atmosphere, and asserting it is
!      nonzero guards against the vector path silently degenerating into
!      the scalar one, which would make every other check here vacuous.
!   7. Tangent linear against central finite differences on a polarized
!      component, the adjoint dot-product identity over all four Stokes
!      components, and K against AD. None of the visible vector Jacobian
!      chain had ever been differentiated before this.
!
! Clear sky throughout: no cloud lookup table is involved, so a failure
! cannot be blamed on coefficient quality.
!
! Exit: STOP 0 if every check passes, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 03-Aug-2026
!
PROGRAM test_VectorRT_VIS_Rayleigh

  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_VectorRT_VIS_Rayleigh'
  CHARACTER(*), PARAMETER :: PATH   = './testinput/'
  CHARACTER(*), PARAMETER :: SENSOR = 'v.abi_g18'

  INTEGER,  PARAMETER :: N_PROFILES  = 2     ! the ECMWF84 loader fills atm(1) and atm(2)
  INTEGER,  PARAMETER :: N_LAYERS    = 100
  INTEGER,  PARAMETER :: N_ABSORBERS = 6
  INTEGER,  PARAMETER :: N_CLOUDS    = 0
  INTEGER,  PARAMETER :: N_AEROSOLS  = 0

  ! v.abi_g18 carries ABI bands 1 to 6. Index 1 is 0.47 um (blue), index 3
  ! is 0.86 um (near infrared). Rayleigh scattering is far stronger at the
  ! blue end, which is what check 3 measures.
  INTEGER,  PARAMETER :: CH_BLUE = 1
  INTEGER,  PARAMETER :: CH_NIR  = 3

  INTEGER,  PARAMETER :: N_AZI = 3
  REAL(fp), PARAMETER :: DPHI(N_AZI) = (/ 0.0_fp, 90.0_fp, 180.0_fp /)
  INTEGER,  PARAMETER :: I_PRINCIPAL = 1, I_PERP = 2, I_ANTI = 3

  REAL(fp), PARAMETER :: ZENITH   = 45.0_fp
  REAL(fp), PARAMETER :: SOLAR    = 30.0_fp
  REAL(fp), PARAMETER :: SENS_AZI = 30.0_fp

  INTEGER,  PARAMETER :: IDX_O3 = 3
  REAL(fp), PARAMETER :: FRAC_O3 = 1.0e-5_fp   ! relative FD step on the ozone profile

  REAL(fp), PARAMETER :: TOL_ZERO  = 1.0e-14_fp  ! U in the principal plane, relative to its peak
  REAL(fp), PARAMETER :: TOL_FD    = 1.0e-5_fp
  ! Two adjoint tiers, set from the measured ladder documented at the
  ! dot-product block below. Cosine tier (I,Q) closes like the scalar path;
  ! the sine tier (U,V) has no m = 0 anchor and runs about three digits
  ! looser. Do not collapse these into one number.
  REAL(fp), PARAMETER :: TOL_ADJ_COS = 1.0e-11_fp
  REAL(fp), PARAMETER :: TOL_ADJ_SIN = 1.0e-8_fp
  REAL(fp), PARAMETER :: TOL_K     = 1.0e-11_fp

  CHARACTER(256) :: Version
  INTEGER  :: Error_Status, Allocate_Status, n_Channels, l, m, ia, k
  LOGICAL  :: all_ok
  REAL(fp) :: LHS, RHS, rel, worst_bound, maxV, pol_blue, pol_nir, dI_max
  REAL(fp) :: tl_val, fd_val, upk, umax, absum, rel_scalar
  REAL(fp), ALLOCATABLE :: S(:,:,:), Rs(:,:), Sp(:), Sm(:), o3_ref(:,:)

  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)  :: Atm_TL(N_PROFILES), Atm_AD(N_PROFILES)
  TYPE(CRTM_Surface_type)     :: Sfc_TL(N_PROFILES), Sfc_AD(N_PROFILES)
  TYPE(CRTM_Options_type)     :: Options(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTS(:,:), RTS_TL(:,:), RTS_AD(:,:), RTS_K(:,:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type),    ALLOCATABLE :: Sfc_K(:,:)

  all_ok = .TRUE.

  CALL CRTM_Version(Version)
  WRITE(*,'(/5x,a)') 'Visible vector RT: Rayleigh polarization, clear sky'
  WRITE(*,'(5x,a/)') 'CRTM Version: '//TRIM(Version)

  Error_Status = CRTM_Init( (/ SENSOR /), ChannelInfo, File_Path = PATH, Quiet = .TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, 'CRTM_Init failed', FAILURE ); STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  IF ( n_Channels < CH_NIR ) THEN
    CALL Display_Message( PROGRAM_NAME, 'unexpected channel count', FAILURE ); STOP 1
  END IF

  ALLOCATE( RTS(n_Channels,N_PROFILES), RTS_TL(n_Channels,N_PROFILES), &
            RTS_AD(n_Channels,N_PROFILES), RTS_K(n_Channels,N_PROFILES), &
            Atm_K(n_Channels,N_PROFILES), Sfc_K(n_Channels,N_PROFILES), &
            S(4,n_Channels,N_AZI), Rs(n_Channels,N_AZI), Sp(4), Sm(4), &
            o3_ref(N_LAYERS,N_PROFILES), STAT=Allocate_Status )
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
    o3_ref(:,m) = Atm(m)%Absorber(:,IDX_O3)
    Atm_TL(m)%Climatology = Atm(m)%Climatology
    Atm_TL(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_TL(m)%Absorber_Units = Atm(m)%Absorber_Units
    Atm_AD(m)%Climatology = Atm(m)%Climatology
    Atm_AD(m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_AD(m)%Absorber_Units = Atm(m)%Absorber_Units
    DO l = 1, n_Channels
      Atm_K(l,m)%Climatology = Atm(m)%Climatology
      Atm_K(l,m)%Absorber_Id = Atm(m)%Absorber_Id ; Atm_K(l,m)%Absorber_Units = Atm(m)%Absorber_Units
    END DO
  END DO

  ! ==================================================================
  ! Forward: relative-azimuth sweep, vector against scalar
  ! ==================================================================
  DO ia = 1, N_AZI
    CALL Set_Geometry( DPHI(ia) )
    CALL Run_Forward( 1 )
    DO l = 1, n_Channels ; Rs(l,ia) = RTS(l,1)%Radiance ; END DO
    CALL Run_Forward( 4 )
    DO l = 1, n_Channels ; S(1:4,l,ia) = RTS(l,1)%Stokes(1:4) ; END DO
  END DO

  WRITE(*,'(7x,a)') 'blue channel (0.47 um) by relative azimuth:'
  DO ia = 1, N_AZI
    WRITE(*,'(9x,f6.1,a,4(es22.15,1x))') DPHI(ia), ' deg  I,Q,U,V = ', S(1:4,CH_BLUE,ia)
  END DO

  ! 1. U vanishes in the principal plane, and is large out of it
  umax = MAXVAL(ABS(S(3,1:n_Channels,I_PERP)))
  upk  = MAX( MAXVAL(ABS(S(3,1:n_Channels,I_PRINCIPAL))), &
              MAXVAL(ABS(S(3,1:n_Channels,I_ANTI))) )
  CALL judge( 'U vanishes at relative azimuth 0 and 180 (odd harmonics)', &
              upk < TOL_ZERO * umax )
  CALL judge( 'U is large at relative azimuth 90 (non-vacuous)', umax > 1.0e-3_fp )
  WRITE(*,'(7x,a,es12.4,a,es12.4)') 'max |U| in principal plane = ', upk, &
        '   at 90 deg = ', umax

  ! 2. Q peaks at Delta-phi = 0 and is negative there
  CALL judge( 'Q is largest in magnitude in the principal plane', &
              ABS(S(2,CH_BLUE,I_PRINCIPAL)) > ABS(S(2,CH_BLUE,I_PERP)) )
  CALL judge( 'Q is negative in the principal plane (perpendicular to scattering plane)', &
              S(2,CH_BLUE,I_PRINCIPAL) < ZERO )

  ! 3. Degree of polarization falls with wavelength
  pol_blue = ABS(S(2,CH_BLUE,I_PRINCIPAL))/S(1,CH_BLUE,I_PRINCIPAL)
  pol_nir  = ABS(S(2,CH_NIR ,I_PRINCIPAL))/S(1,CH_NIR ,I_PRINCIPAL)
  CALL judge( 'degree of polarization falls from 0.47 um to 0.86 um', pol_blue > pol_nir )
  WRITE(*,'(7x,a,f8.4,a,f8.4)') '|Q|/I at 0.47 um = ', pol_blue, &
        '   at 0.86 um = ', pol_nir

  ! 4. V is zero to machine precision
  maxV = MAXVAL(ABS(S(4,1:n_Channels,:)))
  CALL judge( 'V is zero to machine precision (Rayleigh is not circular)', &
              maxV < 1.0e-12_fp * MAXVAL(ABS(S(1,1:n_Channels,:))) )
  WRITE(*,'(7x,a,es12.4)') 'max |V| = ', maxV

  ! 5. Polarization bound
  worst_bound = MINVAL( S(1,1:n_Channels,:)**2 - S(2,1:n_Channels,:)**2 &
                      - S(3,1:n_Channels,:)**2 - S(4,1:n_Channels,:)**2 )
  CALL judge( 'polarization bound I^2 >= Q^2+U^2+V^2 holds', worst_bound >= ZERO )

  ! 6. The vector intensity is NOT the scalar intensity. Without this the
  !    whole test could pass on a path that had silently gone scalar.
  dI_max = MAXVAL(ABS(S(1,1:n_Channels,:) - Rs(1:n_Channels,:)))
  CALL judge( 'vector and scalar intensities differ (scalar-approximation error)', &
              dI_max > ZERO )
  WRITE(*,'(7x,a,es12.4,a,es12.4)') 'max |I(vector) - I(scalar)| = ', dI_max, &
        '   relative = ', dI_max/MAXVAL(ABS(Rs(1:n_Channels,:)))

  ! ==================================================================
  ! 7. Jacobians at relative azimuth 90, where all of Q and U are alive
  ! ==================================================================
  CALL Set_Geometry( 90.0_fp )

  ! ...tangent linear against central finite differences, on U
  CALL CRTM_Atmosphere_Zero( Atm_TL ) ; CALL CRTM_Surface_Zero( Sfc_TL )
  DO m = 1, N_PROFILES
    Atm_TL(m)%Absorber(:,IDX_O3) = o3_ref(:,m)
  END DO
  CALL Run_TL( 4 )
  tl_val = RTS_TL(CH_BLUE,1)%Stokes(3)

  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_O3) = (ONE + FRAC_O3)*o3_ref(:,m)
  END DO
  CALL Run_Forward( 4 ) ; Sp(1:4) = RTS(CH_BLUE,1)%Stokes(1:4)
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_O3) = (ONE - FRAC_O3)*o3_ref(:,m)
  END DO
  CALL Run_Forward( 4 ) ; Sm(1:4) = RTS(CH_BLUE,1)%Stokes(1:4)
  DO m = 1, N_PROFILES
    Atm(m)%Absorber(:,IDX_O3) = o3_ref(:,m)
  END DO
  fd_val = (Sp(3) - Sm(3))/(TWO*FRAC_O3)
  rel = ABS(tl_val - fd_val)/MAX(ABS(fd_val),1.0e-300_fp)
  CALL judge( 'TL dU/d(ozone) against central finite differences', &
              rel < TOL_FD .AND. ABS(fd_val) > ZERO )
  WRITE(*,'(7x,a,es22.15,a,es12.4)') 'dU/dO3  TL = ', tl_val, '   rel = ', rel

  ! ...adjoint dot-product identity, run at n_Stokes = 2 AND 4 with different
  !    tolerances, which is a measured distinction rather than a convenience.
  !
  !    Measured residual ladder on this scene, ozone control vector:
  !       n_Stokes = 2  ->  1.94e-13
  !       n_Stokes = 3  ->  9.71e-11
  !       n_Stokes = 4  ->  1.87e-10
  !    and the n_Stokes = 1 scalar path, seeded through %Radiance, sits at
  !    2.41e-13. So I and Q close to the same precision as the scalar path,
  !    and the degradation appears exactly when U enters.
  !
  !    That has a structural explanation rather than being a defect to chase.
  !    U and V ride SINE azimuthal harmonics, so their m = 0 term is
  !    identically zero and they are assembled entirely from the higher
  !    Fourier components. I and Q ride cosines and take most of their value
  !    from the well-conditioned m = 0 solve. The higher-harmonic adding-
  !    doubling solves are more poorly conditioned, so the U and V adjoint
  !    inherits roughly three digits less accuracy. It is not the phase-matrix
  !    positivity clamp: that was instrumented on this exact scene and fires
  !    0 times in 19200 assembled blocks.
  !
  !    Both tiers are therefore asserted. The cosine tier is held tight, so a
  !    genuine non-transpose in I or Q cannot hide behind the loose tier, and
  !    the U/V tier is held at a level that still catches a real break while
  !    tolerating the conditioning. If the tight tier ever fails, that is a
  !    defect; if only the loose tier fails, suspect conditioning first and
  !    re-measure the ladder before changing any code.
  CALL Run_Forward( 2 )
  CALL Run_TL( 2 )
  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTS_AD )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_AD(l,m)%Stokes(1:2) = RTS_TL(l,m)%Stokes(1:2)
      LHS = LHS + SUM( RTS_TL(l,m)%Stokes(1:2)**2 )
    END DO
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 2
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL judge( 'CRTM_Adjoint call (cosine tier)', .FALSE. )
  ELSE
    RHS = ZERO
    DO m = 1, N_PROFILES
      DO k = 1, N_LAYERS
        RHS = RHS + o3_ref(k,m)*Atm_AD(m)%Absorber(k,IDX_O3)
      END DO
    END DO
    rel_scalar = ABS(LHS - RHS)/MAX(ABS(LHS),1.0e-300_fp)
    CALL judge( 'adjoint identity, cosine tier (I and Q), tight', &
                rel_scalar < TOL_ADJ_COS .AND. LHS > ZERO )
    WRITE(*,'(7x,a,es12.4)') 'n_Stokes=2 adjoint residual = ', rel_scalar
  END IF

  ! ...full four-component identity
  CALL Run_Forward( 4 )
  CALL Run_TL( 4 )
  LHS = ZERO
  CALL CRTM_RTSolution_Zero( RTS_AD )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_AD(l,m)%Stokes(1:4) = RTS_TL(l,m)%Stokes(1:4)
      LHS = LHS + SUM( RTS_TL(l,m)%Stokes(1:4)**2 )
    END DO
  END DO
  CALL CRTM_Atmosphere_Zero( Atm_AD ) ; CALL CRTM_Surface_Zero( Sfc_AD )
  DO m = 1, N_PROFILES
    Options(m)%n_Stokes        = 4
    Options(m)%RT_Algorithm_Id = RT_ADA
  END DO
  Error_Status = CRTM_Adjoint( Atm, Sfc, RTS_AD, Geometry, ChannelInfo, &
                               Atm_AD, Sfc_AD, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL judge( 'CRTM_Adjoint call', .FALSE. )
  ELSE
    RHS = ZERO ; absum = ZERO
    DO m = 1, N_PROFILES
      DO k = 1, N_LAYERS
        RHS   = RHS   + o3_ref(k,m)*Atm_AD(m)%Absorber(k,IDX_O3)
        absum = absum + ABS(o3_ref(k,m)*Atm_AD(m)%Absorber(k,IDX_O3))
      END DO
    END DO
    ! Reported so a future failure can be attributed. If this is ~1 the
    ! residual is not being manufactured by cancellation in this summation,
    ! which is the first thing to rule out. Measured at 1.0.
    WRITE(*,'(7x,a,es12.4)') 'inner-product cancellation sum|t|/|sum t| = ', &
          absum/MAX(ABS(RHS),1.0e-300_fp)
    rel = ABS(LHS - RHS)/MAX(ABS(LHS),1.0e-300_fp)
    CALL judge( 'adjoint identity, full four components (U and V, sine tier)', &
                rel < TOL_ADJ_SIN .AND. LHS > ZERO )
    WRITE(*,'(7x,a,es22.15)') '<TL,y>  = ', LHS
    WRITE(*,'(7x,a,es22.15,a,es12.4)') '<x,ADy> = ', RHS, '   rel = ', rel
  END IF

  ! ...K against AD
  CALL CRTM_RTSolution_Zero( RTS_K )
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      RTS_K(l,m)%Stokes(1:4) = RTS_TL(l,m)%Stokes(1:4)
    END DO
    DO l = 1, n_Channels
      CALL CRTM_Atmosphere_Zero( Atm_K(l:l,m) )
    END DO
  END DO
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTS_K, Geometry, ChannelInfo, &
                                Atm_K, Sfc_K, RTS, Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    CALL judge( 'CRTM_K_Matrix call', .FALSE. )
  ELSE
    LHS = ZERO ; RHS = ZERO
    DO m = 1, N_PROFILES
      DO k = 1, N_LAYERS
        DO l = 1, n_Channels
          LHS = LHS + o3_ref(k,m)*Atm_K(l,m)%Absorber(k,IDX_O3)
        END DO
        RHS = RHS + o3_ref(k,m)*Atm_AD(m)%Absorber(k,IDX_O3)
      END DO
    END DO
    rel = ABS(LHS - RHS)/MAX(ABS(RHS),1.0e-300_fp)
    CALL judge( 'K-matrix summed over channels equals the adjoint', &
                rel < TOL_K .AND. ABS(RHS) > ZERO )
    WRITE(*,'(7x,a,es22.15,a,es12.4)') 'sum(K) = ', LHS, '   rel = ', rel
  END IF

  ! ------------------------------------------------------------------
  CALL CRTM_Atmosphere_Destroy(Atm)
  CALL CRTM_RTSolution_Destroy(RTS)
  Error_Status = CRTM_Destroy( ChannelInfo )

  IF ( all_ok ) THEN
    WRITE(*,'(/5x,a/)') 'ALL CHECKS PASSED'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'ONE OR MORE CHECKS FAILED'
    STOP 1
  END IF

CONTAINS

  SUBROUTINE Set_Geometry( dphi )
    REAL(fp), INTENT(IN) :: dphi
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Sfc(mm) = CRTM_Surface_type()
      Sfc(mm)%Water_Coverage    = ONE
      Sfc(mm)%Water_Type        = 1
      Sfc(mm)%Water_Temperature = 290.0_fp
      Sfc(mm)%Wind_Speed        = 6.25_fp
      Sfc(mm)%Salinity          = 33.0_fp
      CALL CRTM_Geometry_SetValue( Geometry(mm),                          &
                                   Sensor_Zenith_Angle  = ZENITH,         &
                                   Sensor_Azimuth_Angle = SENS_AZI,       &
                                   Source_Zenith_Angle  = SOLAR,          &
                                   Source_Azimuth_Angle = SENS_AZI + dphi )
    END DO
  END SUBROUTINE Set_Geometry

  SUBROUTINE Run_Forward( ns )
    INTEGER, INTENT(IN) :: ns
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Forward( Atm, Sfc, Geometry, ChannelInfo, RTS, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Forward failed', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_Forward

  SUBROUTINE Run_TL( ns )
    INTEGER, INTENT(IN) :: ns
    INTEGER :: mm
    DO mm = 1, N_PROFILES
      Options(mm)%n_Stokes        = ns
      Options(mm)%RT_Algorithm_Id = RT_ADA
    END DO
    Error_Status = CRTM_Tangent_Linear( Atm, Sfc, Atm_TL, Sfc_TL, Geometry, &
                                        ChannelInfo, RTS, RTS_TL, Options=Options )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_Tangent_Linear failed', FAILURE ); STOP 1
    END IF
  END SUBROUTINE Run_TL

  SUBROUTINE judge( label, ok )
    CHARACTER(*), INTENT(IN) :: label
    LOGICAL,      INTENT(IN) :: ok
    IF ( ok ) THEN
      WRITE(*,'(5x,"PASS  ",a)') label
    ELSE
      WRITE(*,'(5x,"FAIL  ",a)') label
      all_ok = .FALSE.
    END IF
  END SUBROUTINE judge

  INCLUDE 'Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_VectorRT_VIS_Rayleigh
