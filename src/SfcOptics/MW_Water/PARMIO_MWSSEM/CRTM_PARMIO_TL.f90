!
! CRTM_PARMIO_TL
!
! Tangent-linear PARMIO microwave ocean surface emissivity / reflectivity.
!

MODULE CRTM_PARMIO_TL

  USE Type_Kinds, ONLY: fp
  USE PARMIOCoeff_Define, ONLY: PARMIOCoeff_type, N_PARMIO_HARMONIC_TERMS
  USE CRTM_PARMIO, ONLY: iVar_type
  USE PARMIO_LUT_Interpolation, ONLY: PARMIO_LUT_Interp_TL
  USE PARMIO_Azimuth_Module, ONLY: PARMIO_Azimuth_Recombine_TL, &
                                   PARMIO_AZ_HARMONIC_FIRST,    &
                                   PARMIO_AZ_HARMONIC_LAST
  USE PARMIO_RC_Interpolation, ONLY: PARMIO_RC_Interp_TL
  USE Reflection_Correction_Module, ONLY: Reflection_Correction_TL
  USE CRTM_MWwaterCoeff, ONLY: MWwaterC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Compute_PARMIO_TL

  INTEGER, PARAMETER :: N_STOKES = 4
  INTEGER, PARAMETER :: Iv_IDX = 1
  INTEGER, PARAMETER :: Ih_IDX = 2
  INTEGER, PARAMETER :: U_IDX  = 3
  INTEGER, PARAMETER :: V_IDX  = 4

  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp

CONTAINS

  SUBROUTINE Compute_PARMIO_TL( &
      PARMIOCoeff,     &  ! Input
      Temperature_TL,  &  ! TL input
      Salinity_TL,     &  ! TL input
      Wind_Speed_TL,   &  ! TL input
      iVar,            &  ! Internal variable input
      Emissivity_TL,   &  ! TL output
      Reflectivity_TL, &  ! TL output
      Azimuth_Angle_TL,&  ! Optional TL input
      Transmittance_TL )  ! Optional TL input
    TYPE(PARMIOCoeff_type), INTENT(IN)  :: PARMIOCoeff
    REAL(fp),               INTENT(IN)  :: Temperature_TL
    REAL(fp),               INTENT(IN)  :: Salinity_TL
    REAL(fp),               INTENT(IN)  :: Wind_Speed_TL
    TYPE(iVar_type),        INTENT(IN)  :: iVar
    REAL(fp),               INTENT(OUT) :: Emissivity_TL(:)
    REAL(fp),               INTENT(OUT) :: Reflectivity_TL(:)
    REAL(fp), OPTIONAL,     INTENT(IN)  :: Azimuth_Angle_TL
    REAL(fp), OPTIONAL,     INTENT(IN)  :: Transmittance_TL

    REAL(fp) :: coefficients_TL(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: foam_fraction_TL
    REAL(fp) :: azimuth_tl
    REAL(fp) :: transmittance_tl_local
    REAL(fp) :: Rv_Mod_TL, Rh_Mod_TL
    REAL(fp) :: rdown_TL(2)

    Emissivity_TL   = ZERO
    Reflectivity_TL = ZERO
    IF (iVar%LUT_Var%Group_ID < 1) RETURN

    CALL PARMIO_LUT_Interp_TL( &
          PARMIOCoeff,                 &
          Frequency_GHz_TL    = ZERO,  &
          Zenith_Angle_deg_TL = ZERO,  &
          Wind_Speed_mps_TL   = Wind_Speed_TL, &
          SST_C_TL            = Temperature_TL, &
          SSS_psu_TL          = Salinity_TL, &
          Coefficients_TL     = coefficients_TL, &
          Foam_Fraction_TL    = foam_fraction_TL, &
          iVar                = iVar%LUT_Var)

    azimuth_tl = ZERO
    IF (PRESENT(Azimuth_Angle_TL) .AND. iVar%Has_Azimuth) THEN
      azimuth_tl = Azimuth_Angle_TL
    END IF
    ! No valid azimuth: the forward dropped the harmonic slots, so their
    ! coefficient perturbations must not leak into the emissivity TL.
    IF (.NOT. iVar%Has_Azimuth) &
      coefficients_TL(PARMIO_AZ_HARMONIC_FIRST:PARMIO_AZ_HARMONIC_LAST) = ZERO

    CALL PARMIO_Azimuth_Recombine_TL( &
          Coefficients         = iVar%Coefficients, &
          Coefficients_TL      = coefficients_TL, &
          Azimuth_Angle_deg    = iVar%Azimuth_Angle, &
          Azimuth_Angle_deg_TL = azimuth_tl, &
          Emissivity_TL        = Emissivity_TL(1:N_STOKES))

    Rv_Mod_TL = ZERO
    Rh_Mod_TL = ZERO
    rdown_TL = ZERO
    transmittance_tl_local = ZERO
    IF (PRESENT(Transmittance_TL)) transmittance_tl_local = Transmittance_TL
    ! Gate the RC linearization on iVar%Has_Transmittance only, mirroring the
    ! forward exactly. transmittance_tl_local already defaults to ZERO when the
    ! optional Transmittance_TL is absent, so the wind-speed sensitivity of the
    ! reflection correction is still propagated when only the transmittance
    ! perturbation is missing (previously the whole kernel was skipped, silently
    ! dropping that term from the Jacobian).
    IF (iVar%Has_Transmittance) THEN
      IF (iVar%Has_PARMIO_RC) THEN
        CALL PARMIO_RC_Interp_TL( &
               PARMIOCoeff, &
               Frequency_GHz_TL     = ZERO, &
               Zenith_Angle_deg_TL  = ZERO, &
               Wind_Speed_mps_TL    = Wind_Speed_TL, &
               SST_C_TL             = Temperature_TL, &
               SSS_psu_TL           = Salinity_TL, &
               Foam_Fraction_TL     = foam_fraction_TL, &
               Transmittance_TL     = transmittance_tl_local, &
               Rdown_TL             = rdown_TL, &
               iVar                 = iVar%PARMIO_RC_Var)
      ELSE
        CALL Reflection_Correction_TL( &
               MWwaterC%RCCoeff, &
               Wind_Speed_TL, &
               transmittance_tl_local, &
               Rv_Mod_TL, &
               Rh_Mod_TL, &
               iVar%RC_Var)
      END IF
    END IF

    ! Base: d(1 - emissivity). The V/H correction overrides this UNLESS the
    ! forward clamped that component (the correction left [0,1] at a near-grazing
    ! quadrature angle), in which case the forward used the bare (1 - emissivity)
    ! and its derivative IS the base term -- so leave it.
    Reflectivity_TL(1:N_STOKES) = -Emissivity_TL(1:N_STOKES)
    ! 3rd/4th Stokes reflectivity is identically zero in the forward.
    Reflectivity_TL(U_IDX) = ZERO
    Reflectivity_TL(V_IDX) = ZERO
    IF (iVar%Has_PARMIO_RC) THEN
      IF (.NOT. iVar%Reflectivity_Clamped(Iv_IDX)) Reflectivity_TL(Iv_IDX) = rdown_TL(1)
      IF (.NOT. iVar%Reflectivity_Clamped(Ih_IDX)) Reflectivity_TL(Ih_IDX) = rdown_TL(2)
    ELSE
      IF (.NOT. iVar%Reflectivity_Clamped(Iv_IDX)) &
        Reflectivity_TL(Iv_IDX) = (ONE - iVar%Emissivity(Iv_IDX)) * Rv_Mod_TL &
                                - iVar%Rv_Mod * Emissivity_TL(Iv_IDX)
      IF (.NOT. iVar%Reflectivity_Clamped(Ih_IDX)) &
        Reflectivity_TL(Ih_IDX) = (ONE - iVar%Emissivity(Ih_IDX)) * Rh_Mod_TL &
                                - iVar%Rh_Mod * Emissivity_TL(Ih_IDX)
    END IF
  END SUBROUTINE Compute_PARMIO_TL

END MODULE CRTM_PARMIO_TL
