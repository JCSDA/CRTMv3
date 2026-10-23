!
! test_ADA_VectorDecoupling
!
! Ground-truth test for the ADA adding-doubling solver under polarimetric
! (n_Stokes > 1) operation, constructed so that it depends on no coefficient
! file, no cloud lookup table, and no external radiative transfer code.
!
! The property being tested
! -------------------------
! ADA is driven directly, with a synthetic scattering phase matrix whose
! polarized blocks are identically zero and whose intensity and Q blocks are
! identical:
!
!     Pff(I,I) = Pff(Q,Q) = P        Pff(I,Q) = Pff(Q,I) = 0
!
! Such a matrix cannot transfer energy between Stokes components. The Stokes
! vector therefore evolves as two independent copies of the same scalar problem,
! and a single n_Stokes = 2 solve must reproduce, exactly, two separate
! n_Stokes = 1 solves that differ only in their surface boundary condition:
!
!     I = ( Iv + Ih ) / 2
!     Q = ( Iv - Ih ) / 2
!
! with the vector run given the surface source ((eV+eH)/2, (eV-eH)/2) and the
! two scalar runs given eV and eH respectively.
!
! Crucially this holds with scattering fully active, which is what the earlier
! scalar-limit probe could not achieve: SCATTERING_ALBEDO_THRESHOLD is 1.0e-10,
! so there is no water content that keeps ADA engaged while disabling its
! scattering coupling. Driving ADA directly removes that obstacle.
!
! Why an unphysical phase matrix is acceptable
! --------------------------------------------
! This is a relative comparison. Whatever the synthetic phase matrix does to the
! intensity, it must do identically in the scalar and vector runs, so absolute
! physical fidelity is irrelevant. That is precisely the separation required
! when the available lookup tables are not yet trustworthy for full Stokes work:
! it tests whether the code is correct, independently of whether the data is.
!
! What a failure means
! --------------------
! Any leakage between Stokes components, any inconsistency in how the adding
! recursion indexes the Stokes slots, any asymmetry in the surface reflection
! matrix assembly, or any mis-striding of the thermal source will break the
! identity. A pass does not prove the polarized physics is right, since the
! polarized blocks are zero here by construction; it proves the vector machinery
! reduces correctly to the scalar case, which is the necessary foundation.
!

PROGRAM test_ADA_VectorDecoupling

  USE CRTM_Module
  USE RTV_Define , ONLY: RTV_type, RTV_Create, RTV_Destroy, RTV_Associated
  USE ADA_Module , ONLY: CRTM_ADA
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_ADA_VectorDecoupling'

  INTEGER , PARAMETER :: N_LAYERS  = 6
  INTEGER , PARAMETER :: N_STREAMS = 2
  INTEGER , PARAMETER :: N_ANGLES  = N_STREAMS + 1   ! CRTM appends the sensor angle
  INTEGER , PARAMETER :: N_LEG     = 8

  ! Surface: a strongly polarized ocean-like boundary so that any leakage
  ! between I and Q is numerically obvious.
  REAL(fp), PARAMETER :: EV = 0.60_fp
  REAL(fp), PARAMETER :: EH = 0.35_fp

  REAL(fp), PARAMETER :: ALBEDO   = 0.4_fp     ! genuine scattering, not a limit case
  REAL(fp), PARAMETER :: TAU_LAY  = 0.15_fp
  REAL(fp), PARAMETER :: PSCAT    = 1.0_fp     ! isotropic; normalisation irrelevant
  REAL(fp), PARAMETER :: B_SFC    = 290.0_fp
  REAL(fp), PARAMETER :: B_ATM    = 250.0_fp
  REAL(fp), PARAMETER :: CBR      = 2.7_fp

  ! The identity is exact algebra, so only round-off needs absorbing. It is
  ! stated relative to the intensity because the recursion accumulates over
  ! layers.
  REAL(fp), PARAMETER :: TOL = 1.0e-11_fp

  TYPE(RTV_type) :: RTV_v, RTV_s
  INTEGER  :: Error_Status, i, n1
  REAL(fp) :: Iv, Ih, Iexp, Qexp, Igot, Qgot, dI, dQ
  LOGICAL  :: ok_I, ok_Q

  WRITE(*,'(/5x,a)')  'ADA vector decoupling: n_Stokes=2 must reduce to two scalar solves'
  WRITE(*,'(5x,a/)')  'Synthetic phase matrix, no coefficient files, scattering ACTIVE'

  ! ------------------------------------------------------------------
  ! Two scalar reference solves, differing only in surface emissivity
  ! ------------------------------------------------------------------
  CALL solve_scalar( EV, Iv )
  CALL solve_scalar( EH, Ih )

  ! ------------------------------------------------------------------
  ! One vector solve with the Stokes-basis surface source
  ! ------------------------------------------------------------------
  CALL solve_vector( POINT_5*(EV+EH), POINT_5*(EV-EH), Igot, Qgot )

  Iexp = POINT_5*(Iv + Ih)
  Qexp = POINT_5*(Iv - Ih)
  dI   = ABS(Igot - Iexp)
  dQ   = ABS(Qgot - Qexp)

  ! Full-precision scalar references. These go through the same solver code as
  ! the vector run at n_Stokes=1, so they are the probe of choice for confirming
  ! bit-identity of the scalar path across any change to ADA: capture them
  ! before and after and compare all 17 digits.
  WRITE(*,'(5x,a,es24.16)') 'scalar-bit-reference Iv ', Iv
  WRITE(*,'(5x,a,es24.16)') 'scalar-bit-reference Ih ', Ih
  WRITE(*,'(5x,a,f14.8)') 'scalar solve, e = eV      : Iv = ', Iv
  WRITE(*,'(5x,a,f14.8)') 'scalar solve, e = eH      : Ih = ', Ih
  WRITE(*,'(5x,a,f14.8,a,f14.8)') 'expected  I = ', Iexp, '   Q = ', Qexp
  WRITE(*,'(5x,a,f14.8,a,f14.8)') 'ADA n_Stokes=2   I = ', Igot, '   Q = ', Qgot
  WRITE(*,'(/5x,a,es12.4)') '|I - (Iv+Ih)/2| = ', dI
  WRITE(*,'(5x,a,es12.4)')  '|Q - (Iv-Ih)/2| = ', dQ
  WRITE(*,'(5x,a,es12.4)')  'tolerance       = ', TOL*MAX(ONE,ABS(Iexp))

  ok_I = ( dI < TOL*MAX(ONE,ABS(Iexp)) )
  ok_Q = ( dQ < TOL*MAX(ONE,ABS(Iexp)) )

  IF ( ok_I .AND. ok_Q ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: ADA vector solve reduces exactly to the scalar solves'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: Stokes components are coupled or mis-indexed in ADA'
    STOP 1
  END IF

CONTAINS

  ! Populate the parts of RTV that ADA reads. n_Stokes MUST be set before
  ! RTV_Create, because RTV_Create sizes Pff and Pbb from RTV%n_Stokes
  ! (RTV_Define.f90 ~line 463).
  SUBROUTINE setup_rtv( RTV, ns )
    TYPE(RTV_type), INTENT(INOUT) :: RTV
    INTEGER,        INTENT(IN)    :: ns
    INTEGER :: ia, ja, i1, j1, k, nZ

    RTV%n_Stokes = ns
    CALL RTV_Create( RTV, N_ANGLES, N_LEG, N_LAYERS )
    IF ( .NOT. RTV_Associated(RTV) ) THEN
      CALL Display_Message( PROGRAM_NAME, 'RTV_Create failed', FAILURE ); STOP 1
    END IF

    RTV%n_Angles       = N_ANGLES
    RTV%n_Streams      = N_STREAMS
    RTV%n_Layers       = N_LAYERS
    RTV%mth_Azi        = 0
    RTV%Scattering_RT  = .TRUE.
    RTV%Solar_Flag_true= .FALSE.
    RTV%Diffuse_Surface= .FALSE.

    ! Quadrature: two hemispheric streams plus the sensor angle. The cosine is
    ! replicated across the Stokes slots of each angle, exactly as
    ! Common_RTSolution does (~line 360), so all components of one angle share a
    ! transmittance.
    RTV%COS_Angle(1)  = 0.80_fp ; RTV%COS_Weight(1) = 0.5_fp
    RTV%COS_Angle(2)  = 0.40_fp ; RTV%COS_Weight(2) = 0.5_fp
    RTV%COS_Angle(3)  = 0.60_fp ; RTV%COS_Weight(3) = 0.0_fp   ! sensor angle
    k = 0
    DO ia = 1, N_ANGLES
      DO ja = 1, ns
        k = k + 1
        RTV%COS_AngleS(k)  = RTV%COS_Angle(ia)
        RTV%COS_WeightS(k) = RTV%COS_Weight(ia)
      END DO
    END DO

    RTV%Planck_Surface          = B_SFC
    RTV%Planck_Atmosphere(0:N_LAYERS) = B_ATM
    RTV%Cosmic_Background_Radiance    = CBR

    ! Synthetic phase matrix: block diagonal in Stokes, identical I and Q
    ! blocks, zero cross terms. Cannot couple Stokes components by construction.
    nZ = N_ANGLES*ns
    RTV%Pff = ZERO
    RTV%Pbb = ZERO
    DO ia = 1, N_ANGLES
      i1 = (ia-1)*ns + 1
      DO ja = 1, N_ANGLES
        j1 = (ja-1)*ns + 1
        DO k = 1, N_LAYERS
          RTV%Pff(i1,j1,k) = PSCAT
          RTV%Pbb(i1,j1,k) = PSCAT
          IF ( ns > 1 ) THEN
            RTV%Pff(i1+1,j1+1,k) = PSCAT     ! Q block identical to I block
            RTV%Pbb(i1+1,j1+1,k) = PSCAT
          END IF
        END DO
      END DO
    END DO
  END SUBROUTINE setup_rtv

  ! Scalar solve with a single surface emissivity.
  SUBROUTINE solve_scalar( e, Iout )
    REAL(fp), INTENT(IN)  :: e
    REAL(fp), INTENT(OUT) :: Iout
    REAL(fp) :: w(N_LAYERS), tau(N_LAYERS)
    REAL(fp) :: emis(N_ANGLES), refl(N_ANGLES,N_ANGLES), dref(N_ANGLES)
    INTEGER  :: ia, ja

    CALL setup_rtv( RTV_s, 1 )
    w   = ALBEDO
    tau = TAU_LAY
    emis = e
    dref = ZERO
    refl = ZERO
    DO ia = 1, N_ANGLES
      refl(ia,ia) = ONE - e
    END DO

    CALL CRTM_ADA( N_LAYERS, w, tau, CBR, emis, refl, dref, RTV_s, Error_Status )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_ADA (scalar) failed', FAILURE ); STOP 1
    END IF
    Iout = RTV_s%s_Level_Rad_UP( N_ANGLES, 0 )   ! sensor angle, TOA
    CALL RTV_Destroy( RTV_s )
  END SUBROUTINE solve_scalar

  ! Vector solve with the Stokes-basis surface source (eI, eQ).
  SUBROUTINE solve_vector( eI, eQ, Iout, Qout )
    REAL(fp), INTENT(IN)  :: eI, eQ
    REAL(fp), INTENT(OUT) :: Iout, Qout
    REAL(fp) :: w(N_LAYERS), tau(N_LAYERS)
    REAL(fp) :: emis(2*N_ANGLES), refl(2*N_ANGLES,2*N_ANGLES), dref(2*N_ANGLES)
    REAL(fp) :: rI, rQ
    INTEGER  :: ia, i1

    CALL setup_rtv( RTV_v, 2 )
    w   = ALBEDO
    tau = TAU_LAY
    dref = ZERO

    ! Surface source and reflection in the Stokes basis, mirroring the
    ! conversion CRTM_SfcOptics applies on the coupled-polarization branch.
    rI = ONE - eI
    rQ =     - eQ
    emis = ZERO
    refl = ZERO
    DO ia = 1, N_ANGLES
      i1 = (ia-1)*2 + 1
      emis(i1)   = eI
      emis(i1+1) = eQ
      refl(i1  ,i1  ) = rI
      refl(i1+1,i1+1) = rI
      refl(i1  ,i1+1) = rQ
      refl(i1+1,i1  ) = rQ
    END DO

    CALL CRTM_ADA( N_LAYERS, w, tau, CBR, emis, refl, dref, RTV_v, Error_Status )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( PROGRAM_NAME, 'CRTM_ADA (vector) failed', FAILURE ); STOP 1
    END IF
    i1   = (N_ANGLES-1)*2 + 1
    Iout = RTV_v%s_Level_Rad_UP( i1  , 0 )
    Qout = RTV_v%s_Level_Rad_UP( i1+1, 0 )
    CALL RTV_Destroy( RTV_v )
  END SUBROUTINE solve_vector

END PROGRAM test_ADA_VectorDecoupling
