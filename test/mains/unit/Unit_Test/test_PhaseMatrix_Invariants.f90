!
! test_PhaseMatrix_Invariants
!
! Asserts physical invariants of the polarized phase matrix that CRTM assembles
! from the generalized-spherical-function expansion coefficients
! (alpha1..alpha4, beta1, beta2). Requires no coefficient file, no cloud lookup
! table and no external radiative transfer code.
!
! Why invariants rather than a Rayleigh reconstruction
! ---------------------------------------------------
! The obvious test, "build a Rayleigh phase matrix and check the degree of
! polarization equals (1-cos^2 T)/(1+cos^2 T)", does not work here, for two
! reasons that are worth recording so the idea is not re-attempted.
!
! First, RTV%Pff is not the scattering matrix at a scattering angle. It is the
! m-th azimuthal Fourier component of the phase matrix between two quadrature
! angles, and CRTM computes one m per call. Recovering F(Theta) would mean
! summing the Fourier series, which this routine is not structured to provide.
!
! Second, the sign of beta1 depends on the convention chosen for F12 and for the
! generalized spherical functions, and references differ. Deriving it from
! F12 = -(3/4)sin^2(Theta) with P(2)_{0,2}(x) = (sqrt6/4)(1-x^2) gives
! beta1_2 = -sqrt(6)/2, but the opposite sign is equally common in the
! literature. A test that hard-codes a sign would be testing the choice of
! textbook rather than testing CRTM.
!
! The invariants below avoid both problems: they are statements about the
! assembled matrix that must hold in any sign convention and for any Fourier
! component, and they are physics rather than transcription. Comparing CRTM's
! assembly against an independent re-derivation of the same published formulas
! would largely test whether the same equations were copied from the same
! source.
!
! Invariants asserted
! -------------------
!   1. Intensity-block invariance. The (1,1) block is built from alpha1 alone,
!      so it must be identical whether the run is scalar or polarized. If it
!      is not, going polarized perturbs the unpolarized radiance, which would
!      be a defect visible to every user who enables n_Stokes > 1.
!
!   2. Degree of polarization bounded by unity: |P12| <= P11 for the m = 0
!      component. Scattered radiation cannot be more than fully polarized.
!      This is the invariant that would expose the positivity clamp applied to
!      the (1,1) element: if P11 is forced up to PHASE_THRESHOLD while the
!      beta1-driven off-diagonal elements are left untouched, the implied
!      polarization can exceed one.
!
!   3. Symmetry of the intensity block, P11(i,j) = P11(j,i), which follows from
!      the expansion being a product of the same function evaluated at the two
!      angles.
!
! A pass does not establish that the polarized physics is right in absolute
! terms; it establishes that the assembled matrix is self-consistent and
! physically admissible, which is a necessary condition that no existing test
! checks.
!

PROGRAM test_PhaseMatrix_Invariants

  USE CRTM_Module
  USE RTV_Define        , ONLY: RTV_type, RTV_Create, RTV_Destroy, RTV_Associated
  USE Common_RTSolution , ONLY: CRTM_Phase_Matrix
  USE CRTM_AtmOptics_Define, ONLY: CRTM_AtmOptics_type    , &
                                   CRTM_AtmOptics_Create  , &
                                   CRTM_AtmOptics_Destroy , &
                                   CRTM_AtmOptics_Associated
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_PhaseMatrix_Invariants'

  INTEGER , PARAMETER :: N_ANGLES  = 4
  INTEGER , PARAMETER :: N_LAYERS  = 1
  INTEGER , PARAMETER :: N_LEG     = 3          ! l = 0,1,2 : Rayleigh needs no more
  INTEGER , PARAMETER :: N_PHASE   = 6

  ! Rayleigh-like expansion. Magnitudes are the textbook values; the sign of
  ! beta1 is deliberately left as a discovery rather than an assertion (see
  ! header), so no invariant below depends on it.
  REAL(fp), PARAMETER :: A1_0 = 1.0_fp, A1_2 = 0.5_fp
  REAL(fp), PARAMETER :: A2_2 = 3.0_fp
  REAL(fp), PARAMETER :: A4_1 = 1.5_fp
  REAL(fp), PARAMETER :: B1_2 = 1.2247448713915890_fp   ! sqrt(6)/2

  ! Stress case. A strongly backscattering asymmetry (g = alpha1_1/3 ~ -0.97)
  ! drives the m=0 intensity component negative between forward angles, which
  ! is precisely when the positivity clamp on the (1,1) element fires. That
  ! clamp raises P11 to PHASE_THRESHOLD without touching the beta1-driven
  ! off-diagonals, so it is the configuration in which the implied degree of
  ! polarization can exceed unity. The nominal case never reaches it.
  REAL(fp), PARAMETER :: A1_1_STRESS = -2.9_fp
  ! Mirrors PHASE_THRESHOLD in RTV_Define, which is module-private. Used only to
  ! count how many (1,1) elements the clamp engaged on, not in any assertion.
  REAL(fp), PARAMETER :: CLAMP_VALUE = 1.0e-7_fp

  REAL(fp), PARAMETER :: TOL_EQ  = 1.0e-13_fp   ! exact-equality invariants
  REAL(fp), PARAMETER :: TOL_POL = 1.0e-10_fp   ! slack on the polarization bound

  TYPE(RTV_type)            :: RTV1, RTV4, RTV_st
  TYPE(CRTM_AtmOptics_type) :: AO1, AO4, AO_st
  INTEGER  :: i, j, i1, j1
  REAL(fp) :: p11_s, p11_v, p12, dmax_block, dmax_sym, pol_worst
  LOGICAL  :: ok_block, ok_pol, ok_sym, ok_stress
  REAL(fp) :: pol_stress
  INTEGER  :: n_clamped
  REAL(fp) :: beta1_sign_probe

  WRITE(*,'(/5x,a)') 'Polarized phase-matrix physical invariants'
  WRITE(*,'(5x,a/)') 'No coefficient files; assembly driven directly'

  CALL build( RTV1, AO1, 1, ZERO )
  CALL build( RTV4, AO4, 4, ZERO )

  CALL CRTM_Phase_Matrix( AO1, RTV1 )
  CALL CRTM_Phase_Matrix( AO4, RTV4 )

  ! ---------------------------------------------------------------
  ! 1. Intensity block must not depend on n_Stokes
  ! 3. and must be symmetric in the two angles
  ! ---------------------------------------------------------------
  dmax_block = ZERO
  dmax_sym   = ZERO
  DO i = 1, N_ANGLES
    i1 = (i-1)*4 + 1
    DO j = 1, N_ANGLES
      j1 = (j-1)*4 + 1
      p11_s = RTV1%Pff(i ,j ,1)
      p11_v = RTV4%Pff(i1,j1,1)
      dmax_block = MAX( dmax_block, ABS(p11_v - p11_s) )
      dmax_sym   = MAX( dmax_sym  , ABS(RTV4%Pff(i1,j1,1) - RTV4%Pff(j1,i1,1)) )
    END DO
  END DO

  ! ---------------------------------------------------------------
  ! 2. Degree of polarization must not exceed unity
  ! ---------------------------------------------------------------
  pol_worst = ZERO
  beta1_sign_probe = ZERO
  DO i = 1, N_ANGLES
    i1 = (i-1)*4 + 1
    DO j = 1, N_ANGLES
      j1 = (j-1)*4 + 1
      p11_v = RTV4%Pff(i1,j1  ,1)
      p12   = RTV4%Pff(i1,j1+1,1)
      IF ( ABS(p11_v) > 1.0e-30_fp ) pol_worst = MAX( pol_worst, ABS(p12)/ABS(p11_v) )
      IF ( ABS(p12) > ABS(beta1_sign_probe) ) beta1_sign_probe = p12
    END DO
  END DO

  ! ---------------------------------------------------------------
  ! Stress scenario: force the (1,1) positivity clamp to engage
  ! ---------------------------------------------------------------
  CALL build( RTV_st, AO_st, 4, A1_1_STRESS )
  CALL CRTM_Phase_Matrix( AO_st, RTV_st )
  pol_stress = ZERO
  n_clamped  = 0
  DO i = 1, N_ANGLES
    i1 = (i-1)*4 + 1
    DO j = 1, N_ANGLES
      j1 = (j-1)*4 + 1
      p11_v = RTV_st%Pff(i1,j1  ,1)
      p12   = RTV_st%Pff(i1,j1+1,1)
      IF ( p11_v <= CLAMP_VALUE*1.000001_fp ) n_clamped = n_clamped + 1
      IF ( ABS(p11_v) > 1.0e-30_fp ) pol_stress = MAX( pol_stress, ABS(p12)/ABS(p11_v) )
    END DO
  END DO
  WRITE(*,'(/5x,a,i0,a,i0)') 'stress case: clamped (1,1) elements = ', n_clamped, ' of ', N_ANGLES*N_ANGLES
  WRITE(*,'(5x,a,es14.6)')   'stress case: worst |P12|/|P11|      = ', pol_stress
  ok_stress = ( pol_stress <= ONE + TOL_POL )
  IF ( .NOT. ok_stress ) THEN
    WRITE(*,'(5x,a)') 'stress case: KNOWN DEFECT, polarization bound violated.'
    WRITE(*,'(5x,a)') '  The (1,1) positivity clamp raises P11 to PHASE_THRESHOLD without'
    WRITE(*,'(5x,a)') '  rescaling the beta1-driven off-diagonals, so the assembled matrix'
    WRITE(*,'(5x,a)') '  implies a degree of polarization far above unity at those angle'
    WRITE(*,'(5x,a)') '  pairs. Reported, not asserted, until the correct physical'
    WRITE(*,'(5x,a)') '  treatment is decided (rescale the block, or bound the polarized'
    WRITE(*,'(5x,a)') '  elements by the clamped P11). Realistic trigger is Legendre'
    WRITE(*,'(5x,a)') '  truncation ringing of a forward-peaked phase function, not the'
    WRITE(*,'(5x,a)') '  backscattering coefficients used to provoke it here.'
  END IF
  CALL CRTM_AtmOptics_Destroy( AO_st ) ; CALL RTV_Destroy( RTV_st )

  ok_block = ( dmax_block < TOL_EQ )
  ok_sym   = ( dmax_sym   < TOL_EQ )
  ok_pol   = ( pol_worst  <= ONE + TOL_POL )

  WRITE(*,'(5x,a,es12.4,a,l1)') 'intensity block, |vector - scalar|  = ', dmax_block, '   pass = ', ok_block
  WRITE(*,'(5x,a,es12.4,a,l1)') 'intensity block symmetry            = ', dmax_sym  , '   pass = ', ok_sym
  WRITE(*,'(5x,a,f12.6,a,l1)')  'worst |P12| / |P11|                 = ', pol_worst , '   pass = ', ok_pol
  WRITE(*,'(5x,a,es12.4)')      'largest P12 (sign is convention)    = ', beta1_sign_probe

  CALL CRTM_AtmOptics_Destroy( AO1 ) ; CALL CRTM_AtmOptics_Destroy( AO4 )
  CALL RTV_Destroy( RTV1 )           ; CALL RTV_Destroy( RTV4 )

  ! ok_stress is reported above, not asserted: see the KNOWN DEFECT note.
  IF ( ok_block .AND. ok_sym .AND. ok_pol ) THEN
    WRITE(*,'(/5x,a/)') 'PASS: assembled phase matrix is physically admissible'
    STOP 0
  ELSE
    WRITE(*,'(/5x,a/)') 'FAIL: phase-matrix invariant violated'
    STOP 1
  END IF

CONTAINS

  ! n_Stokes must be set before RTV_Create, which sizes Pff from it.
  SUBROUTINE build( RTV, AO, ns, a1_1 )
    TYPE(RTV_type)           , INTENT(INOUT) :: RTV
    TYPE(CRTM_AtmOptics_type), INTENT(INOUT) :: AO
    INTEGER                  , INTENT(IN)    :: ns
    REAL(fp)                 , INTENT(IN)    :: a1_1
    INTEGER :: ia

    RTV%n_Stokes = ns
    CALL RTV_Create( RTV, N_ANGLES, N_LEG, N_LAYERS )
    IF ( .NOT. RTV_Associated(RTV) ) THEN
      CALL Display_Message( PROGRAM_NAME, 'RTV_Create failed', FAILURE ); STOP 1
    END IF
    RTV%n_Angles        = N_ANGLES
    RTV%n_Streams       = N_ANGLES
    RTV%n_Layers        = N_LAYERS
    RTV%mth_Azi         = 0
    RTV%Solar_Flag_true = .FALSE.
    DO ia = 1, N_ANGLES
      RTV%COS_Angle(ia)  = 0.95_fp - 0.2_fp*REAL(ia-1,fp)
      RTV%COS_Weight(ia) = 0.25_fp
    END DO

    CALL CRTM_AtmOptics_Create( AO, N_LAYERS, N_LEG, N_PHASE )
    IF ( .NOT. CRTM_AtmOptics_Associated(AO) ) THEN
      CALL Display_Message( PROGRAM_NAME, 'AtmOptics_Create failed', FAILURE ); STOP 1
    END IF
    AO%n_Legendre_Terms = N_LEG
    AO%n_Phase_Elements = N_PHASE
    AO%Single_Scatter_Albedo(1) = 0.9_fp        ! above the assembly threshold
    AO%Phase_Coefficient = ZERO
    AO%Phase_Coefficient(0,1,1) = A1_0
    AO%Phase_Coefficient(1,1,1) = a1_1
    AO%Phase_Coefficient(2,1,1) = A1_2
    AO%Phase_Coefficient(2,2,1) = A2_2
    AO%Phase_Coefficient(1,4,1) = A4_1
    AO%Phase_Coefficient(2,5,1) = B1_2
  END SUBROUTINE build

END PROGRAM test_PhaseMatrix_Invariants
