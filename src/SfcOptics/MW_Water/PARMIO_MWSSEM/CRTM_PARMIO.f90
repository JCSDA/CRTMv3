!
! CRTM_PARMIO
!
! Compute_PARMIO: forward microwave ocean surface emissivity / reflectivity
! from the PARMIOCoeff lookup table. Mirrors the Compute_FastemX positional
! signature so the dispatcher swap in CRTM_MW_Water_SfcOptics.f90 is a
! one-line edit.
!
! TL/AD live in CRTM_PARMIO_TL.f90 and CRTM_PARMIO_AD.f90 (separate compile
! units to keep this file scannable); they consume the iVar_type stashed
! by Compute_PARMIO here.

MODULE CRTM_PARMIO

  USE Type_Kinds, ONLY: fp
  USE PARMIOCoeff_Define, ONLY: PARMIOCoeff_type, N_PARMIO_HARMONIC_TERMS
  USE PARMIO_LUT_Interpolation, ONLY: &
        PARMIO_LUT_iVar_type,         &
        PARMIO_LUT_Interp_Forward
  USE PARMIO_Azimuth_Module, ONLY: PARMIO_Azimuth_Recombine
  USE PARMIO_RC_Interpolation, ONLY: &
        PARMIO_RC_iVar_type,         &
        PARMIO_RC_Interp_Forward
  USE Reflection_Correction_Module, ONLY: &
        RC_iVar_type => iVar_type,        &
        Reflection_Correction
  USE CRTM_MWwaterCoeff, ONLY: MWwaterC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: iVar_type
  PUBLIC :: Compute_PARMIO

  ! CRTM microwave surface-optics index convention. These mirror
  ! CRTM_FastemX.f90:121-122; slots 1/2 are decoupled V/H, not
  ! canonical Stokes I/Q.
  INTEGER, PARAMETER :: N_STOKES = 4
  INTEGER, PARAMETER :: Iv_IDX = 1
  INTEGER, PARAMETER :: Ih_IDX = 2

  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  REAL(fp), PARAMETER :: PI   = 3.141592653589793238462643383279_fp
  REAL(fp), PARAMETER :: DEGREES_TO_RADIANS = PI / 180.0_fp
  ! Generous physical band for the catastrophic-reflectivity guard (see the
  ! clamp in Compute_PARMIO). Reflectivity is physically in [0,1]; the band is
  ! wide enough to leave a small RT-tolerated grazing overshoot untouched but
  ! catches the gross blow-up. Matches CRTM_FastemX.
  REAL(fp), PARAMETER :: R_PHYS_LO = -0.5_fp
  REAL(fp), PARAMETER :: R_PHYS_HI =  1.5_fp

  ! ---------------------------------------------------------------
  ! Internal-state carrier: hands forward results into TL/AD.
  ! ---------------------------------------------------------------
  TYPE :: iVar_type
    ! LUT bracketing + interpolated coefficients
    TYPE(PARMIO_LUT_iVar_type) :: LUT_Var
    ! Coefficients evaluated at query point (foam-blended)
    REAL(fp) :: Coefficients(N_PARMIO_HARMONIC_TERMS) = ZERO
    REAL(fp) :: Foam_Fraction = ZERO
    ! Query inputs (saved for TL/AD chains)
    REAL(fp) :: Frequency      = ZERO
    REAL(fp) :: Zenith_Angle   = ZERO
    REAL(fp) :: cos_z          = ONE
    REAL(fp) :: Wind_Speed     = ZERO
    REAL(fp) :: Azimuth_Angle  = ZERO
    REAL(fp) :: Transmittance  = ONE
    LOGICAL  :: Has_Azimuth    = .FALSE.
    LOGICAL  :: Has_Transmittance = .FALSE.
    ! Emissivity (V-pol, H-pol, U, V_circ) and pre-RC reflectivity, plus the
    ! Reflection_Correction iVar so TL/AD can drive the FASTEM RC kernel.
    REAL(fp) :: Emissivity(N_STOKES)   = ZERO
    REAL(fp) :: Reflectivity(N_STOKES) = ZERO
    REAL(fp) :: Rv_Mod = ONE
    REAL(fp) :: Rh_Mod = ONE
    LOGICAL  :: Has_PARMIO_RC = .FALSE.
    ! Per-Stokes flag: the V/H reflection correction was clamped to the bare
    ! (1 - emissivity) because it left [0,1] (near-grazing quadrature angles).
    ! TL/AD read this to drop the (blown-up) correction derivative for that
    ! component and use the bare d(1 - emissivity) instead.
    LOGICAL  :: Reflectivity_Clamped(N_STOKES) = .FALSE.
    TYPE(RC_iVar_type) :: RC_Var
    TYPE(PARMIO_RC_iVar_type) :: PARMIO_RC_Var
  END TYPE iVar_type

CONTAINS

  !-----------------------------------------------------------------
  !  Compute_PARMIO
  !  Same positional signature as Compute_FastemX in
  !  CRTMv3/src/SfcOptics/MW_Water/FASTEM_MWSSEM/CRTM_FastemX.f90.
  !-----------------------------------------------------------------
  SUBROUTINE Compute_PARMIO( &
      PARMIOCoeff,    &  ! in  coefficient/LUT struct
      Frequency,      &  ! in  GHz
      n_Angles,       &  ! in  (kept for API symmetry; ignored — single-angle call)
      Zenith_Angle,   &  ! in  deg
      Temperature,    &  ! in  K  (SST)
      Salinity,       &  ! in  ppt
      Wind_Speed,     &  ! in  m/s
      iVar,           &  ! out internal state for TL/AD
      Emissivity,     &  ! out (4) V-pol, H-pol, U, V_circ
      Reflectivity,   &  ! out (4)
      Azimuth_Angle,  &  ! in,opt deg (relative to wind)
      Transmittance)     ! in,opt
    TYPE(PARMIOCoeff_type), INTENT(IN)  :: PARMIOCoeff
    REAL(fp),               INTENT(IN)  :: Frequency
    INTEGER,                INTENT(IN)  :: n_Angles
    REAL(fp),               INTENT(IN)  :: Zenith_Angle
    REAL(fp),               INTENT(IN)  :: Temperature
    REAL(fp),               INTENT(IN)  :: Salinity
    REAL(fp),               INTENT(IN)  :: Wind_Speed
    TYPE(iVar_type),        INTENT(OUT) :: iVar
    REAL(fp),               INTENT(OUT) :: Emissivity(:)
    REAL(fp),               INTENT(OUT) :: Reflectivity(:)
    REAL(fp), OPTIONAL,     INTENT(IN)  :: Azimuth_Angle
    REAL(fp), OPTIONAL,     INTENT(IN)  :: Transmittance
    REAL(fp) :: phi_deg, SST_C, SSS_psu
    REAL(fp) :: rdown(2)

    ! Save query inputs
    iVar%Frequency      = Frequency
    iVar%Zenith_Angle   = Zenith_Angle
    iVar%cos_z          = COS(Zenith_Angle * DEGREES_TO_RADIANS)
    iVar%Wind_Speed     = Wind_Speed
    iVar%Has_Azimuth       = PRESENT(Azimuth_Angle)
    iVar%Has_Transmittance = PRESENT(Transmittance)
    IF (iVar%Has_Azimuth)       iVar%Azimuth_Angle  = Azimuth_Angle
    IF (iVar%Has_Transmittance) iVar%Transmittance  = Transmittance
    phi_deg = ZERO
    IF (iVar%Has_Azimuth) phi_deg = Azimuth_Angle

    ! Convert SST K → C, salinity passthrough (PARMIO LUT uses C/psu)
    SST_C   = Temperature - 273.15_fp
    SSS_psu = Salinity

    ! 1) LUT interpolation → 14 dimensionless harmonic coefficients
    !    (already foam-blended) at (freq, theta, U10, sst, sss)
    CALL PARMIO_LUT_Interp_Forward( &
          PARMIOCoeff,                        &
          Frequency_GHz    = Frequency,       &
          Zenith_Angle_deg = Zenith_Angle,    &
          Wind_Speed_mps   = Wind_Speed,      &
          SST_C            = SST_C,           &
          SSS_psu          = SSS_psu,         &
          Coefficients     = iVar%Coefficients, &
          Foam_Fraction    = iVar%Foam_Fraction, &
          iVar             = iVar%LUT_Var)

    ! 2) Recombine the 14 coefficients into CRTM's microwave surface-
    !    optics basis at the requested azimuth (phi=0 if azimuth absent,
    !    matching FastemX).
    CALL PARMIO_Azimuth_Recombine( &
          Coefficients      = iVar%Coefficients, &
          Azimuth_Angle_deg = phi_deg,           &
          Emissivity        = iVar%Emissivity)

    ! 3) Bare reflectivity = 1 - emissivity (per polarization).
    !    The Reflection_Correction below scales V/H reflectivity.
    iVar%Reflectivity = ONE - iVar%Emissivity

    ! 4) Transmittance-dependent bistatic-scattering correction on V/H.
    !    Mirrors CRTM_FastemX.f90:432, but uses MWwaterC%RCCoeff for the
    !    FASTEM-fit polynomial. PARMIOCoeff currently does NOT carry its
    !    own RCCoeff; this is acceptable as a Phase-4 baseline (see plan
    !    §4.5). When Transmittance is absent or RCCoeff is unallocated
    !    we leave Reflectivity at the (1 - Emissivity) value.
    iVar%Rv_Mod = ONE
    iVar%Rh_Mod = ONE
    IF (iVar%Has_Transmittance) THEN
      CALL PARMIO_RC_Interp_Forward( &
            PARMIOCoeff,       &
            Frequency,         &
            Zenith_Angle,      &
            Wind_Speed,        &
            SST_C,             &
            SSS_psu,           &
            iVar%Foam_Fraction,&
            Transmittance,     &
            rdown,             &
            iVar%PARMIO_RC_Var,&
            iVar%Has_PARMIO_RC)
      IF (iVar%Has_PARMIO_RC) THEN
        iVar%Reflectivity(Iv_IDX) = rdown(1)
        iVar%Reflectivity(Ih_IDX) = rdown(2)
      ELSE
        CALL Reflection_Correction( &
              MWwaterC%RCCoeff, &
              Frequency,        &
              iVar%cos_z,       &
              Wind_Speed,       &
              Transmittance,    &
              iVar%Rv_Mod,      &
              iVar%Rh_Mod,      &
              iVar%RC_Var)
        iVar%Reflectivity(Iv_IDX) = iVar%Rv_Mod * (ONE - iVar%Emissivity(Iv_IDX))
        iVar%Reflectivity(Ih_IDX) = iVar%Rh_Mod * (ONE - iVar%Emissivity(Ih_IDX))
      END IF
    END IF

    ! Catastrophic-reflectivity guard. The V/H reflection correction above is a
    ! FASTEM-fit polynomial valid only for typical view angles; the scattering RT
    ! evaluates the surface optics at Gaussian quadrature angles up to ~86 deg
    ! (near grazing), where it can extrapolate to a wildly non-physical
    ! reflectivity (observed: V-pol ~ -1e35 at za=86 deg). Clamp only GROSSLY
    ! out-of-range values (a generous [R_PHYS_LO,R_PHYS_HI] band), falling back to
    ! the bare (1 - emissivity), which is physical by construction; a small
    ! RT-tolerated overshoot is left alone so validated results are unchanged.
    ! Without this the garbage reflectivity propagates into the adding-doubling
    ! surface boundary and blows the radiance up (only reachable on the
    ! cloudy/scattering path at
    ! >= 200 GHz, i.e. the PARMIO regime).
    iVar%Reflectivity_Clamped = .FALSE.
    IF ( iVar%Reflectivity(Iv_IDX) < R_PHYS_LO .OR. iVar%Reflectivity(Iv_IDX) > R_PHYS_HI ) THEN
      iVar%Reflectivity(Iv_IDX) = MIN( MAX( ONE - iVar%Emissivity(Iv_IDX), ZERO ), ONE )
      iVar%Reflectivity_Clamped(Iv_IDX) = .TRUE.
    END IF
    IF ( iVar%Reflectivity(Ih_IDX) < R_PHYS_LO .OR. iVar%Reflectivity(Ih_IDX) > R_PHYS_HI ) THEN
      iVar%Reflectivity(Ih_IDX) = MIN( MAX( ONE - iVar%Emissivity(Ih_IDX), ZERO ), ONE )
      iVar%Reflectivity_Clamped(Ih_IDX) = .TRUE.
    END IF

    ! Write outputs
    Emissivity   = iVar%Emissivity
    Reflectivity = iVar%Reflectivity

  END SUBROUTINE Compute_PARMIO

END MODULE CRTM_PARMIO
