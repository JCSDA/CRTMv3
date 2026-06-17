!
! CRTM_PARMIO_AD
!
! Adjoint PARMIO microwave ocean surface emissivity / reflectivity.
!

MODULE CRTM_PARMIO_AD

  USE Type_Kinds, ONLY: fp
  USE PARMIOCoeff_Define, ONLY: PARMIOCoeff_type, N_PARMIO_HARMONIC_TERMS
  USE CRTM_PARMIO, ONLY: iVar_type
  USE PARMIO_LUT_Interpolation, ONLY: PARMIO_LUT_Interp_AD
  USE PARMIO_Azimuth_Module, ONLY: PARMIO_Azimuth_Recombine_AD
  USE PARMIO_RC_Interpolation, ONLY: PARMIO_RC_Interp_AD
  USE Reflection_Correction_Module, ONLY: Reflection_Correction_AD
  USE CRTM_MWwaterCoeff, ONLY: MWwaterC

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: Compute_PARMIO_AD

  INTEGER, PARAMETER :: N_STOKES = 4
  INTEGER, PARAMETER :: Iv_IDX = 1
  INTEGER, PARAMETER :: Ih_IDX = 2
  INTEGER, PARAMETER :: U_IDX  = 3
  INTEGER, PARAMETER :: V_IDX  = 4

  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp

CONTAINS

  SUBROUTINE Compute_PARMIO_AD( &
      PARMIOCoeff,     &  ! Input
      Emissivity_AD,   &  ! AD input
      Reflectivity_AD, &  ! AD input
      iVar,            &  ! Internal variable input
      Temperature_AD,  &  ! AD output
      Salinity_AD,     &  ! AD output
      Wind_Speed_AD,   &  ! AD output
      Azimuth_Angle_AD,&  ! Optional AD output
      Transmittance_AD )  ! Optional AD output
    TYPE(PARMIOCoeff_type), INTENT(IN)     :: PARMIOCoeff
    REAL(fp),               INTENT(IN OUT) :: Emissivity_AD(:)
    REAL(fp),               INTENT(IN OUT) :: Reflectivity_AD(:)
    TYPE(iVar_type),        INTENT(IN)     :: iVar
    REAL(fp),               INTENT(IN OUT) :: Temperature_AD
    REAL(fp),               INTENT(IN OUT) :: Salinity_AD
    REAL(fp),               INTENT(IN OUT) :: Wind_Speed_AD
    REAL(fp), OPTIONAL,     INTENT(IN OUT) :: Azimuth_Angle_AD
    REAL(fp), OPTIONAL,     INTENT(IN OUT) :: Transmittance_AD

    REAL(fp) :: coefficients_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: foam_fraction_AD
    REAL(fp) :: azimuth_AD
    REAL(fp) :: Rv_Mod_AD, Rh_Mod_AD
    REAL(fp) :: rdown_AD(2)
    REAL(fp) :: frequency_AD, theta_AD, sst_AD, sss_AD
    REAL(fp) :: transmittance_AD_local

    IF (iVar%LUT_Var%Group_ID < 1) THEN
      Emissivity_AD   = ZERO
      Reflectivity_AD = ZERO
      RETURN
    END IF

    Rv_Mod_AD = ZERO
    Rh_Mod_AD = ZERO
    rdown_AD = ZERO
    foam_fraction_AD = ZERO
    frequency_AD = ZERO
    theta_AD = ZERO
    sst_AD = ZERO
    sss_AD = ZERO
    transmittance_AD_local = ZERO
    IF (PRESENT(Transmittance_AD)) transmittance_AD_local = Transmittance_AD

    Emissivity_AD(U_IDX) = Emissivity_AD(U_IDX) - Reflectivity_AD(U_IDX)
    Reflectivity_AD(U_IDX) = ZERO
    Emissivity_AD(V_IDX) = Emissivity_AD(V_IDX) - Reflectivity_AD(V_IDX)
    Reflectivity_AD(V_IDX) = ZERO

    ! V/H components the forward clamped to the bare (1 - emissivity): apply that
    ! transpose and zero the seed, so the (blown-up) correction adjoint below is
    ! skipped for them. Transpose of the TL's "leave the base term" branch.
    IF (iVar%Reflectivity_Clamped(Iv_IDX)) THEN
      Emissivity_AD(Iv_IDX) = Emissivity_AD(Iv_IDX) - Reflectivity_AD(Iv_IDX)
      Reflectivity_AD(Iv_IDX) = ZERO
    END IF
    IF (iVar%Reflectivity_Clamped(Ih_IDX)) THEN
      Emissivity_AD(Ih_IDX) = Emissivity_AD(Ih_IDX) - Reflectivity_AD(Ih_IDX)
      Reflectivity_AD(Ih_IDX) = ZERO
    END IF

    IF (iVar%Has_PARMIO_RC) THEN
      rdown_AD(2) = Reflectivity_AD(Ih_IDX)
      Reflectivity_AD(Ih_IDX) = ZERO
      rdown_AD(1) = Reflectivity_AD(Iv_IDX)
      Reflectivity_AD(Iv_IDX) = ZERO
      CALL PARMIO_RC_Interp_AD( &
             PARMIOCoeff, &
             Rdown_AD          = rdown_AD, &
             iVar              = iVar%PARMIO_RC_Var, &
             Foam_Fraction_AD  = foam_fraction_AD, &
             Frequency_GHz_AD  = frequency_AD, &
             Zenith_Angle_deg_AD = theta_AD, &
             Wind_Speed_mps_AD = Wind_Speed_AD, &
             SST_C_AD          = sst_AD, &
             SSS_psu_AD        = sss_AD, &
             Transmittance_AD  = transmittance_AD_local)
      IF (PRESENT(Transmittance_AD)) Transmittance_AD = transmittance_AD_local
    ELSE
      Emissivity_AD(Ih_IDX) = Emissivity_AD(Ih_IDX) - &
                              iVar%Rh_Mod * Reflectivity_AD(Ih_IDX)
      Rh_Mod_AD = (ONE - iVar%Emissivity(Ih_IDX)) * Reflectivity_AD(Ih_IDX)
      Reflectivity_AD(Ih_IDX) = ZERO

      Emissivity_AD(Iv_IDX) = Emissivity_AD(Iv_IDX) - &
                              iVar%Rv_Mod * Reflectivity_AD(Iv_IDX)
      Rv_Mod_AD = (ONE - iVar%Emissivity(Iv_IDX)) * Reflectivity_AD(Iv_IDX)
      Reflectivity_AD(Iv_IDX) = ZERO

      ! Transpose of the TL: gate on iVar%Has_Transmittance only. Accumulate the
      ! transmittance adjoint through a local so an absent optional Transmittance_AD
      ! does not suppress the wind-speed adjoint of the reflection correction;
      ! write the local back to the caller's accumulator only when it is present.
      IF (iVar%Has_Transmittance) THEN
        CALL Reflection_Correction_AD( &
               MWwaterC%RCCoeff, &
               Rv_Mod_AD, &
               Rh_Mod_AD, &
               Wind_Speed_AD, &
               transmittance_AD_local, &
               iVar%RC_Var)
        IF (PRESENT(Transmittance_AD)) Transmittance_AD = transmittance_AD_local
      ELSE
        Rv_Mod_AD = ZERO
        Rh_Mod_AD = ZERO
      END IF
    END IF

    coefficients_AD = ZERO
    azimuth_AD = ZERO
    CALL PARMIO_Azimuth_Recombine_AD( &
          Coefficients          = iVar%Coefficients, &
          Emissivity_AD         = Emissivity_AD(1:N_STOKES), &
          Azimuth_Angle_deg     = iVar%Azimuth_Angle, &
          Coefficients_AD       = coefficients_AD, &
          Azimuth_Angle_deg_AD  = azimuth_AD)
    IF (PRESENT(Azimuth_Angle_AD) .AND. iVar%Has_Azimuth) THEN
      Azimuth_Angle_AD = Azimuth_Angle_AD + azimuth_AD
    END IF

    CALL PARMIO_LUT_Interp_AD( &
          PARMIOCoeff, &
          Coefficients_AD    = coefficients_AD, &
          Foam_Fraction_AD   = foam_fraction_AD, &
          iVar               = iVar%LUT_Var, &
          Frequency_GHz_AD   = frequency_AD, &
          Zenith_Angle_deg_AD = theta_AD, &
          Wind_Speed_mps_AD  = Wind_Speed_AD, &
          SST_C_AD           = sst_AD, &
          SSS_psu_AD         = sss_AD)
    Temperature_AD = Temperature_AD + sst_AD
    Salinity_AD    = Salinity_AD + sss_AD
  END SUBROUTINE Compute_PARMIO_AD

END MODULE CRTM_PARMIO_AD
