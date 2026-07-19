!
! PARMIO_LUT_Interpolation
!
! Multilinear interpolation kernel for the PARMIOCoeff lookup table.
! Selects the appropriate frequency group (sss_dependent / sss_nominal_m /
! sss_nominal_h), brackets each axis, and returns the 14 azimuthal-harmonic
! coefficients at the query point. Out-of-range inputs are clamped to the
! nearest grid edge and the clamp is recorded in iVar for diagnostics.
!
! Foam blending: PARMIO emits one set of coefficients with foam contribution
! applied (foam_on) and one without (foam_off). The runtime blends with
!   coeff(k) = (1 - F) * coeff_off(k) + F * coeff_on(k)
! where F = foam fraction (0..1) interpolated from the LUT's foam-on slot.
!
! Forward only for now; TL/AD will live in companion modules but reuse
! the bracketing weights stashed in iVar_type.

MODULE PARMIO_LUT_Interpolation

  USE Type_Kinds, ONLY: fp, Double
  USE PARMIOCoeff_Define, ONLY: &
    PARMIOCoeff_type, PARMIOCoeff_Group_type, &
    PARMIOCoeff_GroupName_For_Frequency,      &
    N_PARMIO_HARMONIC_TERMS,                  &
    PARMIO_FOAM_OFF, PARMIO_FOAM_ON,          &
    PARMIO_GROUP_SSS_DEPENDENT,               &
    PARMIO_N_GROUPS

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PARMIO_LUT_iVar_type
  PUBLIC :: PARMIO_LUT_Interp_Forward
  PUBLIC :: PARMIO_LUT_Interp_TL
  PUBLIC :: PARMIO_LUT_Interp_AD

  ! ---------------------------------------------------------------
  ! 1-D bracket descriptor: lo/hi indices and the linear weight
  !   value(query) = (1-w) * arr(lo) + w * arr(hi)
  ! ---------------------------------------------------------------
  TYPE :: Bracket_1D_type
    INTEGER  :: lo = 1
    INTEGER  :: hi = 1
    REAL(fp) :: w  = 0.0_fp
    LOGICAL  :: clamped_low  = .FALSE.
    LOGICAL  :: clamped_high = .FALSE.
  END TYPE Bracket_1D_type

  ! ---------------------------------------------------------------
  ! Internal-state record carried from forward into TL / AD.
  ! ---------------------------------------------------------------
  TYPE :: PARMIO_LUT_iVar_type
    INTEGER :: Group_ID  = 0
    LOGICAL :: SSS_Active = .FALSE.
    TYPE(Bracket_1D_type) :: B_Frequency
    TYPE(Bracket_1D_type) :: B_Theta
    TYPE(Bracket_1D_type) :: B_Wind
    TYPE(Bracket_1D_type) :: B_SST
    TYPE(Bracket_1D_type) :: B_SSS
    REAL(fp) :: Foam_Fraction = 0.0_fp   ! 0..1, interpolated from LUT
    REAL(fp) :: Coeff_Foam_Off(N_PARMIO_HARMONIC_TERMS) = 0.0_fp
    REAL(fp) :: Coeff_Foam_On (N_PARMIO_HARMONIC_TERMS) = 0.0_fp
  END TYPE PARMIO_LUT_iVar_type

CONTAINS


  !-----------------------------------------------------------------
  !  PARMIO_LUT_Interp_Forward
  !  Look up + interpolate the 14 harmonic coefficients at a query
  !  point. Out-of-range inputs are clamped silently; the clamp is
  !  recorded on iVar for downstream diagnostics.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_LUT_Interp_Forward( &
      LUT, Frequency_GHz, Zenith_Angle_deg, Wind_Speed_mps, &
      SST_C, SSS_psu, Coefficients, Foam_Fraction, iVar)
    TYPE(PARMIOCoeff_type),       INTENT(IN)  :: LUT
    REAL(fp),                     INTENT(IN)  :: Frequency_GHz
    REAL(fp),                     INTENT(IN)  :: Zenith_Angle_deg
    REAL(fp),                     INTENT(IN)  :: Wind_Speed_mps
    REAL(fp),                     INTENT(IN)  :: SST_C
    REAL(fp),                     INTENT(IN)  :: SSS_psu
    REAL(fp),                     INTENT(OUT) :: Coefficients(N_PARMIO_HARMONIC_TERMS)
    REAL(fp),                     INTENT(OUT) :: Foam_Fraction
    TYPE(PARMIO_LUT_iVar_type),   INTENT(OUT) :: iVar
    INTEGER :: g
    REAL(fp) :: cutoff_sss, cutoff_perm

    ! Select group from frequency
    cutoff_sss  = LUT%SSS_Cutoff_GHz
    cutoff_perm = LUT%Permittivity_Switch_GHz
    g = PARMIOCoeff_GroupName_For_Frequency(  &
          Frequency_GHz,                      &
          SSS_Cutoff_GHz_Override          = cutoff_sss, &
          Permittivity_Switch_GHz_Override = cutoff_perm)
    ! Unallocated group: return inert zeros with Group_ID=0, which the TL/AD
    ! already treat as "no contribution" (same guard they carry).
    IF (.NOT. LUT%Group(g)%Is_Allocated) THEN
      iVar%Group_ID = 0
      Coefficients  = 0.0_fp
      Foam_Fraction = 0.0_fp
      RETURN
    END IF
    iVar%Group_ID  = g
    iVar%SSS_Active = LUT%Group(g)%SSS_Axis_Active

    ! Bracket each axis (clamp to range)
    CALL Bracket(LUT%Group(g)%Frequency,  Frequency_GHz,    iVar%B_Frequency)
    CALL Bracket(LUT%Group(g)%Theta,      Zenith_Angle_deg, iVar%B_Theta)
    CALL Bracket(LUT%Group(g)%Wind_Speed, Wind_Speed_mps,   iVar%B_Wind)
    CALL Bracket(LUT%Group(g)%SST,        SST_C,            iVar%B_SST)
    IF (iVar%SSS_Active) THEN
      CALL Bracket(LUT%Group(g)%SSS, SSS_psu, iVar%B_SSS)
    ELSE
      iVar%B_SSS%lo = 1
      iVar%B_SSS%hi = 1
      iVar%B_SSS%w  = 0.0_fp
    END IF

    ! Interpolate coefficients on each foam slot
    CALL Interp_All_Harmonics(LUT%Group(g), PARMIO_FOAM_OFF, iVar, &
                              iVar%Coeff_Foam_Off)
    CALL Interp_All_Harmonics(LUT%Group(g), PARMIO_FOAM_ON,  iVar, &
                              iVar%Coeff_Foam_On)

    ! Foam fraction comes from the foam-on slot (foam-off has Foam=0).
    ! PARMIO emits Foam in percent; convert to fraction.
    iVar%Foam_Fraction = 0.01_fp * Interp_Foam( &
        LUT%Group(g), PARMIO_FOAM_ON, iVar)
    Foam_Fraction = iVar%Foam_Fraction

    ! Foam-fraction-weighted blend of the two foam slots
    Coefficients = (1.0_fp - iVar%Foam_Fraction) * iVar%Coeff_Foam_Off &
                 +           iVar%Foam_Fraction  * iVar%Coeff_Foam_On

  END SUBROUTINE PARMIO_LUT_Interp_Forward


  !-----------------------------------------------------------------
  !  PARMIO_LUT_Interp_TL
  !  Tangent-linear of PARMIO_LUT_Interp_Forward for the state
  !  variables exposed by CRTM_PARMIO_TL.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_LUT_Interp_TL( &
      LUT, Frequency_GHz_TL, Zenith_Angle_deg_TL, Wind_Speed_mps_TL, &
      SST_C_TL, SSS_psu_TL, Coefficients_TL, Foam_Fraction_TL, iVar)
    TYPE(PARMIOCoeff_type),     INTENT(IN)  :: LUT
    REAL(fp),                   INTENT(IN)  :: Frequency_GHz_TL
    REAL(fp),                   INTENT(IN)  :: Zenith_Angle_deg_TL
    REAL(fp),                   INTENT(IN)  :: Wind_Speed_mps_TL
    REAL(fp),                   INTENT(IN)  :: SST_C_TL
    REAL(fp),                   INTENT(IN)  :: SSS_psu_TL
    REAL(fp),                   INTENT(OUT) :: Coefficients_TL(N_PARMIO_HARMONIC_TERMS)
    REAL(fp),                   INTENT(OUT) :: Foam_Fraction_TL
    TYPE(PARMIO_LUT_iVar_type), INTENT(IN)  :: iVar
    INTEGER :: g
    REAL(fp) :: coef_off_TL(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: coef_on_TL(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: foam_raw_TL

    Coefficients_TL   = 0.0_fp
    Foam_Fraction_TL  = 0.0_fp
    g = iVar%Group_ID
    IF (g < 1 .OR. g > PARMIO_N_GROUPS) RETURN
    IF (.NOT. LUT%Group(g)%Is_Allocated) RETURN

    CALL Interp_All_Harmonics_TL( &
          LUT%Group(g), PARMIO_FOAM_OFF, iVar, &
          Frequency_GHz_TL, Zenith_Angle_deg_TL, Wind_Speed_mps_TL, &
          SST_C_TL, SSS_psu_TL, coef_off_TL)
    CALL Interp_All_Harmonics_TL( &
          LUT%Group(g), PARMIO_FOAM_ON, iVar, &
          Frequency_GHz_TL, Zenith_Angle_deg_TL, Wind_Speed_mps_TL, &
          SST_C_TL, SSS_psu_TL, coef_on_TL)
    CALL Interp_Foam_TL( &
          LUT%Group(g), PARMIO_FOAM_ON, iVar, &
          Frequency_GHz_TL, Zenith_Angle_deg_TL, Wind_Speed_mps_TL, &
          SST_C_TL, SSS_psu_TL, foam_raw_TL)

    Foam_Fraction_TL = 0.01_fp * foam_raw_TL
    Coefficients_TL = (1.0_fp - iVar%Foam_Fraction) * coef_off_TL &
                    +            iVar%Foam_Fraction  * coef_on_TL  &
                    + Foam_Fraction_TL * (iVar%Coeff_Foam_On - iVar%Coeff_Foam_Off)
  END SUBROUTINE PARMIO_LUT_Interp_TL


  !-----------------------------------------------------------------
  !  PARMIO_LUT_Interp_AD
  !  Adjoint of PARMIO_LUT_Interp_Forward for the state variables
  !  exposed by CRTM_PARMIO_AD. Adjoint inputs are zeroed on exit.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_LUT_Interp_AD( &
      LUT, Coefficients_AD, Foam_Fraction_AD, iVar, &
      Frequency_GHz_AD, Zenith_Angle_deg_AD, Wind_Speed_mps_AD, &
      SST_C_AD, SSS_psu_AD)
    TYPE(PARMIOCoeff_type),     INTENT(IN)     :: LUT
    REAL(fp),                   INTENT(IN OUT) :: Coefficients_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp),                   INTENT(IN OUT) :: Foam_Fraction_AD
    TYPE(PARMIO_LUT_iVar_type), INTENT(IN)     :: iVar
    REAL(fp),                   INTENT(IN OUT) :: Frequency_GHz_AD
    REAL(fp),                   INTENT(IN OUT) :: Zenith_Angle_deg_AD
    REAL(fp),                   INTENT(IN OUT) :: Wind_Speed_mps_AD
    REAL(fp),                   INTENT(IN OUT) :: SST_C_AD
    REAL(fp),                   INTENT(IN OUT) :: SSS_psu_AD
    INTEGER :: g
    REAL(fp) :: coef_off_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: coef_on_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp) :: foam_raw_AD
    REAL(fp) :: wfreq_AD(2), wth_AD(2), wu_AD(2), ws_AD(2), wq_AD(2)

    g = iVar%Group_ID
    IF (g < 1 .OR. g > PARMIO_N_GROUPS) THEN
      Coefficients_AD  = 0.0_fp
      Foam_Fraction_AD = 0.0_fp
      RETURN
    END IF
    IF (.NOT. LUT%Group(g)%Is_Allocated) THEN
      Coefficients_AD  = 0.0_fp
      Foam_Fraction_AD = 0.0_fp
      RETURN
    END IF

    coef_off_AD = (1.0_fp - iVar%Foam_Fraction) * Coefficients_AD
    coef_on_AD  =            iVar%Foam_Fraction  * Coefficients_AD
    Foam_Fraction_AD = Foam_Fraction_AD + &
      SUM((iVar%Coeff_Foam_On - iVar%Coeff_Foam_Off) * Coefficients_AD)
    Coefficients_AD = 0.0_fp

    foam_raw_AD = 0.01_fp * Foam_Fraction_AD
    Foam_Fraction_AD = 0.0_fp

    wfreq_AD = 0.0_fp
    wth_AD   = 0.0_fp
    wu_AD    = 0.0_fp
    ws_AD    = 0.0_fp
    wq_AD    = 0.0_fp

    CALL Interp_Foam_AD( &
          LUT%Group(g), PARMIO_FOAM_ON, iVar, foam_raw_AD, &
          wfreq_AD, wth_AD, wu_AD, ws_AD, wq_AD)
    CALL Interp_All_Harmonics_AD( &
          LUT%Group(g), PARMIO_FOAM_ON, iVar, coef_on_AD, &
          wfreq_AD, wth_AD, wu_AD, ws_AD, wq_AD)
    CALL Interp_All_Harmonics_AD( &
          LUT%Group(g), PARMIO_FOAM_OFF, iVar, coef_off_AD, &
          wfreq_AD, wth_AD, wu_AD, ws_AD, wq_AD)

    CALL Axis_Value_AD(LUT%Group(g)%Frequency,   iVar%B_Frequency, wfreq_AD, Frequency_GHz_AD)
    CALL Axis_Value_AD(LUT%Group(g)%Theta,       iVar%B_Theta,     wth_AD,   Zenith_Angle_deg_AD)
    CALL Axis_Value_AD(LUT%Group(g)%Wind_Speed,  iVar%B_Wind,      wu_AD,    Wind_Speed_mps_AD)
    CALL Axis_Value_AD(LUT%Group(g)%SST,         iVar%B_SST,       ws_AD,    SST_C_AD)
    IF (iVar%SSS_Active) THEN
      CALL Axis_Value_AD(LUT%Group(g)%SSS,       iVar%B_SSS,       wq_AD,    SSS_psu_AD)
    END IF
  END SUBROUTINE PARMIO_LUT_Interp_AD


  !-----------------------------------------------------------------
  !  Bracket: find lo/hi indices and linear weight for a 1-D axis.
  !  Clamps to range when the query falls outside.
  !-----------------------------------------------------------------
  SUBROUTINE Bracket(axis, query, b)
    REAL(Double),            INTENT(IN)  :: axis(:)
    REAL(fp),                INTENT(IN)  :: query
    TYPE(Bracket_1D_type),   INTENT(OUT) :: b
    INTEGER  :: n, k, lo, hi
    REAL(fp) :: ax_lo, ax_hi
    n = SIZE(axis)
    IF (n == 1) THEN
      b%lo = 1; b%hi = 1; b%w = 0.0_fp
      RETURN
    END IF
    IF (query <= REAL(axis(1), fp)) THEN
      b%lo = 1; b%hi = 1; b%w = 0.0_fp
      b%clamped_low = (query < REAL(axis(1), fp))
      RETURN
    END IF
    IF (query >= REAL(axis(n), fp)) THEN
      b%lo = n; b%hi = n; b%w = 0.0_fp
      b%clamped_high = (query > REAL(axis(n), fp))
      RETURN
    END IF
    ! Linear search (axes are small, < 100 nodes); avoids importing a
    ! search utility just for this.
    DO k = 1, n - 1
      ax_lo = REAL(axis(k),     fp)
      ax_hi = REAL(axis(k + 1), fp)
      IF (query >= ax_lo .AND. query <= ax_hi) THEN
        b%lo = k
        b%hi = k + 1
        IF (ax_hi > ax_lo) THEN
          b%w = (query - ax_lo) / (ax_hi - ax_lo)
        ELSE
          b%w = 0.0_fp
        END IF
        RETURN
      END IF
    END DO
    ! Unreachable (caught by the bounds checks above), but be safe.
    b%lo = n; b%hi = n; b%w = 0.0_fp
  END SUBROUTINE Bracket


  !-----------------------------------------------------------------
  !  Multilinear interpolation of all 14 harmonic terms at a single
  !  foam slot.  Visits 16 (with SSS) or 32 (without SSS) corners
  !  with the (1-w) / w weight expansion.
  !  Note: dim order is (k, foam, sss, sst, U10, theta, freq).
  !-----------------------------------------------------------------
  SUBROUTINE Interp_All_Harmonics(grp, foam_idx, iVar, out_coef)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN)  :: grp
    INTEGER,                      INTENT(IN)  :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN)  :: iVar
    REAL(fp),                     INTENT(OUT) :: out_coef(N_PARMIO_HARMONIC_TERMS)

    INTEGER :: ifreq, jfreq, ith, jth, iu, ju, is, js, iq, jq
    REAL(fp) :: wfreq0, wfreq1, wth0, wth1, wu0, wu1, ws0, ws1, wq0, wq1

    ifreq = iVar%B_Frequency%lo;  jfreq = iVar%B_Frequency%hi
    ith   = iVar%B_Theta%lo;      jth   = iVar%B_Theta%hi
    iu    = iVar%B_Wind%lo;       ju    = iVar%B_Wind%hi
    is    = iVar%B_SST%lo;        js    = iVar%B_SST%hi
    iq    = iVar%B_SSS%lo;        jq    = iVar%B_SSS%hi
    wfreq1 = iVar%B_Frequency%w; wfreq0 = 1.0_fp - wfreq1
    wth1   = iVar%B_Theta%w;     wth0   = 1.0_fp - wth1
    wu1    = iVar%B_Wind%w;      wu0    = 1.0_fp - wu1
    ws1    = iVar%B_SST%w;       ws0    = 1.0_fp - ws1
    wq1    = iVar%B_SSS%w;       wq0    = 1.0_fp - wq1

    ! Coefficients(k, foam, sss, sst, U10, theta, freq) — the explicit
    ! 32-term unfold is verbose but maps directly to the chain rule used
    ! by TL/AD; keep it as a single hot loop for cache friendliness.
    out_coef = wq0 * wfreq0 * ( &
        ws0 * wu0 * wth0 * grp%Coefficients(:, foam_idx, iq, is, iu, ith, ifreq) &
      + ws0 * wu0 * wth1 * grp%Coefficients(:, foam_idx, iq, is, iu, jth, ifreq) &
      + ws0 * wu1 * wth0 * grp%Coefficients(:, foam_idx, iq, is, ju, ith, ifreq) &
      + ws0 * wu1 * wth1 * grp%Coefficients(:, foam_idx, iq, is, ju, jth, ifreq) &
      + ws1 * wu0 * wth0 * grp%Coefficients(:, foam_idx, iq, js, iu, ith, ifreq) &
      + ws1 * wu0 * wth1 * grp%Coefficients(:, foam_idx, iq, js, iu, jth, ifreq) &
      + ws1 * wu1 * wth0 * grp%Coefficients(:, foam_idx, iq, js, ju, ith, ifreq) &
      + ws1 * wu1 * wth1 * grp%Coefficients(:, foam_idx, iq, js, ju, jth, ifreq))
    out_coef = out_coef + wq0 * wfreq1 * ( &
        ws0 * wu0 * wth0 * grp%Coefficients(:, foam_idx, iq, is, iu, ith, jfreq) &
      + ws0 * wu0 * wth1 * grp%Coefficients(:, foam_idx, iq, is, iu, jth, jfreq) &
      + ws0 * wu1 * wth0 * grp%Coefficients(:, foam_idx, iq, is, ju, ith, jfreq) &
      + ws0 * wu1 * wth1 * grp%Coefficients(:, foam_idx, iq, is, ju, jth, jfreq) &
      + ws1 * wu0 * wth0 * grp%Coefficients(:, foam_idx, iq, js, iu, ith, jfreq) &
      + ws1 * wu0 * wth1 * grp%Coefficients(:, foam_idx, iq, js, iu, jth, jfreq) &
      + ws1 * wu1 * wth0 * grp%Coefficients(:, foam_idx, iq, js, ju, ith, jfreq) &
      + ws1 * wu1 * wth1 * grp%Coefficients(:, foam_idx, iq, js, ju, jth, jfreq))
    IF (iVar%SSS_Active) THEN
      out_coef = out_coef + wq1 * wfreq0 * ( &
          ws0 * wu0 * wth0 * grp%Coefficients(:, foam_idx, jq, is, iu, ith, ifreq) &
        + ws0 * wu0 * wth1 * grp%Coefficients(:, foam_idx, jq, is, iu, jth, ifreq) &
        + ws0 * wu1 * wth0 * grp%Coefficients(:, foam_idx, jq, is, ju, ith, ifreq) &
        + ws0 * wu1 * wth1 * grp%Coefficients(:, foam_idx, jq, is, ju, jth, ifreq) &
        + ws1 * wu0 * wth0 * grp%Coefficients(:, foam_idx, jq, js, iu, ith, ifreq) &
        + ws1 * wu0 * wth1 * grp%Coefficients(:, foam_idx, jq, js, iu, jth, ifreq) &
        + ws1 * wu1 * wth0 * grp%Coefficients(:, foam_idx, jq, js, ju, ith, ifreq) &
        + ws1 * wu1 * wth1 * grp%Coefficients(:, foam_idx, jq, js, ju, jth, ifreq))
      out_coef = out_coef + wq1 * wfreq1 * ( &
          ws0 * wu0 * wth0 * grp%Coefficients(:, foam_idx, jq, is, iu, ith, jfreq) &
        + ws0 * wu0 * wth1 * grp%Coefficients(:, foam_idx, jq, is, iu, jth, jfreq) &
        + ws0 * wu1 * wth0 * grp%Coefficients(:, foam_idx, jq, is, ju, ith, jfreq) &
        + ws0 * wu1 * wth1 * grp%Coefficients(:, foam_idx, jq, is, ju, jth, jfreq) &
        + ws1 * wu0 * wth0 * grp%Coefficients(:, foam_idx, jq, js, iu, ith, jfreq) &
        + ws1 * wu0 * wth1 * grp%Coefficients(:, foam_idx, jq, js, iu, jth, jfreq) &
        + ws1 * wu1 * wth0 * grp%Coefficients(:, foam_idx, jq, js, ju, ith, jfreq) &
        + ws1 * wu1 * wth1 * grp%Coefficients(:, foam_idx, jq, js, ju, jth, jfreq))
    END IF
  END SUBROUTINE Interp_All_Harmonics


  !-----------------------------------------------------------------
  !  Multilinear interpolation of the foam-fraction scalar at a
  !  single foam slot. Same logic as Interp_All_Harmonics but
  !  without the leading harmonic-index dimension.
  !-----------------------------------------------------------------
  REAL(fp) FUNCTION Interp_Foam(grp, foam_idx, iVar) RESULT(F)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN) :: grp
    INTEGER,                      INTENT(IN) :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN) :: iVar
    INTEGER :: ifreq, jfreq, ith, jth, iu, ju, is, js, iq, jq
    REAL(fp) :: wfreq0, wfreq1, wth0, wth1, wu0, wu1, ws0, ws1, wq0, wq1
    ifreq = iVar%B_Frequency%lo;  jfreq = iVar%B_Frequency%hi
    ith   = iVar%B_Theta%lo;      jth   = iVar%B_Theta%hi
    iu    = iVar%B_Wind%lo;       ju    = iVar%B_Wind%hi
    is    = iVar%B_SST%lo;        js    = iVar%B_SST%hi
    iq    = iVar%B_SSS%lo;        jq    = iVar%B_SSS%hi
    wfreq1 = iVar%B_Frequency%w; wfreq0 = 1.0_fp - wfreq1
    wth1   = iVar%B_Theta%w;     wth0   = 1.0_fp - wth1
    wu1    = iVar%B_Wind%w;      wu0    = 1.0_fp - wu1
    ws1    = iVar%B_SST%w;       ws0    = 1.0_fp - ws1
    wq1    = iVar%B_SSS%w;       wq0    = 1.0_fp - wq1
    F = wq0 * wfreq0 * ( &
        ws0 * wu0 * wth0 * grp%Foam(foam_idx, iq, is, iu, ith, ifreq) &
      + ws0 * wu0 * wth1 * grp%Foam(foam_idx, iq, is, iu, jth, ifreq) &
      + ws0 * wu1 * wth0 * grp%Foam(foam_idx, iq, is, ju, ith, ifreq) &
      + ws0 * wu1 * wth1 * grp%Foam(foam_idx, iq, is, ju, jth, ifreq) &
      + ws1 * wu0 * wth0 * grp%Foam(foam_idx, iq, js, iu, ith, ifreq) &
      + ws1 * wu0 * wth1 * grp%Foam(foam_idx, iq, js, iu, jth, ifreq) &
      + ws1 * wu1 * wth0 * grp%Foam(foam_idx, iq, js, ju, ith, ifreq) &
      + ws1 * wu1 * wth1 * grp%Foam(foam_idx, iq, js, ju, jth, ifreq))
    F = F + wq0 * wfreq1 * ( &
        ws0 * wu0 * wth0 * grp%Foam(foam_idx, iq, is, iu, ith, jfreq) &
      + ws0 * wu0 * wth1 * grp%Foam(foam_idx, iq, is, iu, jth, jfreq) &
      + ws0 * wu1 * wth0 * grp%Foam(foam_idx, iq, is, ju, ith, jfreq) &
      + ws0 * wu1 * wth1 * grp%Foam(foam_idx, iq, is, ju, jth, jfreq) &
      + ws1 * wu0 * wth0 * grp%Foam(foam_idx, iq, js, iu, ith, jfreq) &
      + ws1 * wu0 * wth1 * grp%Foam(foam_idx, iq, js, iu, jth, jfreq) &
      + ws1 * wu1 * wth0 * grp%Foam(foam_idx, iq, js, ju, ith, jfreq) &
      + ws1 * wu1 * wth1 * grp%Foam(foam_idx, iq, js, ju, jth, jfreq))
    IF (iVar%SSS_Active) THEN
      F = F + wq1 * wfreq0 * ( &
          ws0 * wu0 * wth0 * grp%Foam(foam_idx, jq, is, iu, ith, ifreq) &
        + ws0 * wu0 * wth1 * grp%Foam(foam_idx, jq, is, iu, jth, ifreq) &
        + ws0 * wu1 * wth0 * grp%Foam(foam_idx, jq, is, ju, ith, ifreq) &
        + ws0 * wu1 * wth1 * grp%Foam(foam_idx, jq, is, ju, jth, ifreq) &
        + ws1 * wu0 * wth0 * grp%Foam(foam_idx, jq, js, iu, ith, ifreq) &
        + ws1 * wu0 * wth1 * grp%Foam(foam_idx, jq, js, iu, jth, ifreq) &
        + ws1 * wu1 * wth0 * grp%Foam(foam_idx, jq, js, ju, ith, ifreq) &
        + ws1 * wu1 * wth1 * grp%Foam(foam_idx, jq, js, ju, jth, ifreq))
      F = F + wq1 * wfreq1 * ( &
          ws0 * wu0 * wth0 * grp%Foam(foam_idx, jq, is, iu, ith, jfreq) &
        + ws0 * wu0 * wth1 * grp%Foam(foam_idx, jq, is, iu, jth, jfreq) &
        + ws0 * wu1 * wth0 * grp%Foam(foam_idx, jq, is, ju, ith, jfreq) &
        + ws0 * wu1 * wth1 * grp%Foam(foam_idx, jq, is, ju, jth, jfreq) &
        + ws1 * wu0 * wth0 * grp%Foam(foam_idx, jq, js, iu, ith, jfreq) &
        + ws1 * wu0 * wth1 * grp%Foam(foam_idx, jq, js, iu, jth, jfreq) &
        + ws1 * wu1 * wth0 * grp%Foam(foam_idx, jq, js, ju, ith, jfreq) &
        + ws1 * wu1 * wth1 * grp%Foam(foam_idx, jq, js, ju, jth, jfreq))
    END IF
  END FUNCTION Interp_Foam


  SUBROUTINE Interp_All_Harmonics_TL( &
      grp, foam_idx, iVar, Frequency_TL, Theta_TL, Wind_TL, SST_TL, SSS_TL, &
      out_coef_TL)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN)  :: grp
    INTEGER,                      INTENT(IN)  :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN)  :: iVar
    REAL(fp),                     INTENT(IN)  :: Frequency_TL
    REAL(fp),                     INTENT(IN)  :: Theta_TL
    REAL(fp),                     INTENT(IN)  :: Wind_TL
    REAL(fp),                     INTENT(IN)  :: SST_TL
    REAL(fp),                     INTENT(IN)  :: SSS_TL
    REAL(fp),                     INTENT(OUT) :: out_coef_TL(N_PARMIO_HARMONIC_TERMS)
    INTEGER :: ifreq(2), ith(2), iu(2), isst(2), iq(2)
    INTEGER :: af, at, au, asst, aq, nq
    REAL(fp) :: wf(2), wth(2), wu(2), ws(2), wq(2)
    REAL(fp) :: wf_TL(2), wth_TL(2), wu_TL(2), ws_TL(2), wq_TL(2)
    REAL(fp) :: cw_TL

    CALL Axis_Weights(iVar%B_Frequency, ifreq, wf)
    CALL Axis_Weights(iVar%B_Theta,     ith,   wth)
    CALL Axis_Weights(iVar%B_Wind,      iu,    wu)
    CALL Axis_Weights(iVar%B_SST,       isst,  ws)
    CALL Axis_Weights(iVar%B_SSS,       iq,    wq)
    CALL Axis_Weights_TL(grp%Frequency,  iVar%B_Frequency, Frequency_TL, wf_TL)
    CALL Axis_Weights_TL(grp%Theta,      iVar%B_Theta,     Theta_TL,     wth_TL)
    CALL Axis_Weights_TL(grp%Wind_Speed, iVar%B_Wind,      Wind_TL,      wu_TL)
    CALL Axis_Weights_TL(grp%SST,        iVar%B_SST,       SST_TL,       ws_TL)
    IF (iVar%SSS_Active) THEN
      CALL Axis_Weights_TL(grp%SSS,      iVar%B_SSS,       SSS_TL,       wq_TL)
      nq = 2
    ELSE
      wq_TL = 0.0_fp
      nq = 1
    END IF

    out_coef_TL = 0.0_fp
    DO aq = 1, nq
      DO af = 1, 2
        DO at = 1, 2
          DO au = 1, 2
            DO asst = 1, 2
              cw_TL = wq_TL(aq)*wf(af)*wth(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf_TL(af)*wth(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth_TL(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth(at)*wu_TL(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth(at)*wu(au)*ws_TL(asst)
              out_coef_TL = out_coef_TL + cw_TL * REAL( &
                grp%Coefficients(:, foam_idx, iq(aq), isst(asst), &
                                 iu(au), ith(at), ifreq(af)), fp)
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE Interp_All_Harmonics_TL


  SUBROUTINE Interp_Foam_TL( &
      grp, foam_idx, iVar, Frequency_TL, Theta_TL, Wind_TL, SST_TL, SSS_TL, &
      foam_TL)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN)  :: grp
    INTEGER,                      INTENT(IN)  :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN)  :: iVar
    REAL(fp),                     INTENT(IN)  :: Frequency_TL
    REAL(fp),                     INTENT(IN)  :: Theta_TL
    REAL(fp),                     INTENT(IN)  :: Wind_TL
    REAL(fp),                     INTENT(IN)  :: SST_TL
    REAL(fp),                     INTENT(IN)  :: SSS_TL
    REAL(fp),                     INTENT(OUT) :: foam_TL
    INTEGER :: ifreq(2), ith(2), iu(2), isst(2), iq(2)
    INTEGER :: af, at, au, asst, aq, nq
    REAL(fp) :: wf(2), wth(2), wu(2), ws(2), wq(2)
    REAL(fp) :: wf_TL(2), wth_TL(2), wu_TL(2), ws_TL(2), wq_TL(2)
    REAL(fp) :: cw_TL

    CALL Axis_Weights(iVar%B_Frequency, ifreq, wf)
    CALL Axis_Weights(iVar%B_Theta,     ith,   wth)
    CALL Axis_Weights(iVar%B_Wind,      iu,    wu)
    CALL Axis_Weights(iVar%B_SST,       isst,  ws)
    CALL Axis_Weights(iVar%B_SSS,       iq,    wq)
    CALL Axis_Weights_TL(grp%Frequency,  iVar%B_Frequency, Frequency_TL, wf_TL)
    CALL Axis_Weights_TL(grp%Theta,      iVar%B_Theta,     Theta_TL,     wth_TL)
    CALL Axis_Weights_TL(grp%Wind_Speed, iVar%B_Wind,      Wind_TL,      wu_TL)
    CALL Axis_Weights_TL(grp%SST,        iVar%B_SST,       SST_TL,       ws_TL)
    IF (iVar%SSS_Active) THEN
      CALL Axis_Weights_TL(grp%SSS,      iVar%B_SSS,       SSS_TL,       wq_TL)
      nq = 2
    ELSE
      wq_TL = 0.0_fp
      nq = 1
    END IF

    foam_TL = 0.0_fp
    DO aq = 1, nq
      DO af = 1, 2
        DO at = 1, 2
          DO au = 1, 2
            DO asst = 1, 2
              cw_TL = wq_TL(aq)*wf(af)*wth(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf_TL(af)*wth(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth_TL(at)*wu(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth(at)*wu_TL(au)*ws(asst) + &
                      wq(aq)*wf(af)*wth(at)*wu(au)*ws_TL(asst)
              foam_TL = foam_TL + cw_TL * REAL( &
                grp%Foam(foam_idx, iq(aq), isst(asst), &
                         iu(au), ith(at), ifreq(af)), fp)
            END DO
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE Interp_Foam_TL


  SUBROUTINE Interp_All_Harmonics_AD( &
      grp, foam_idx, iVar, out_coef_AD, &
      wfreq_AD, wth_AD, wu_AD, ws_AD, wq_AD)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN)     :: grp
    INTEGER,                      INTENT(IN)     :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN)     :: iVar
    REAL(fp),                     INTENT(IN OUT) :: out_coef_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp),                     INTENT(IN OUT) :: wfreq_AD(2), wth_AD(2)
    REAL(fp),                     INTENT(IN OUT) :: wu_AD(2), ws_AD(2), wq_AD(2)
    INTEGER :: ifreq(2), ith(2), iu(2), isst(2), iq(2)
    INTEGER :: af, at, au, asst, aq, nq
    REAL(fp) :: wf(2), wth(2), wu(2), ws(2), wq(2)
    REAL(fp) :: cw_AD

    CALL Axis_Weights(iVar%B_Frequency, ifreq, wf)
    CALL Axis_Weights(iVar%B_Theta,     ith,   wth)
    CALL Axis_Weights(iVar%B_Wind,      iu,    wu)
    CALL Axis_Weights(iVar%B_SST,       isst,  ws)
    CALL Axis_Weights(iVar%B_SSS,       iq,    wq)
    nq = 1
    IF (iVar%SSS_Active) nq = 2

    DO aq = 1, nq
      DO af = 1, 2
        DO at = 1, 2
          DO au = 1, 2
            DO asst = 1, 2
              cw_AD = SUM(REAL( &
                grp%Coefficients(:, foam_idx, iq(aq), isst(asst), &
                                 iu(au), ith(at), ifreq(af)), fp) * out_coef_AD)
              wq_AD(aq)    = wq_AD(aq)    + wf(af)*wth(at)*wu(au)*ws(asst) * cw_AD
              wfreq_AD(af) = wfreq_AD(af) + wq(aq)*wth(at)*wu(au)*ws(asst) * cw_AD
              wth_AD(at)   = wth_AD(at)   + wq(aq)*wf(af)*wu(au)*ws(asst)  * cw_AD
              wu_AD(au)    = wu_AD(au)    + wq(aq)*wf(af)*wth(at)*ws(asst) * cw_AD
              ws_AD(asst)  = ws_AD(asst)  + wq(aq)*wf(af)*wth(at)*wu(au)   * cw_AD
            END DO
          END DO
        END DO
      END DO
    END DO
    out_coef_AD = 0.0_fp
  END SUBROUTINE Interp_All_Harmonics_AD


  SUBROUTINE Interp_Foam_AD( &
      grp, foam_idx, iVar, foam_AD, &
      wfreq_AD, wth_AD, wu_AD, ws_AD, wq_AD)
    TYPE(PARMIOCoeff_Group_type), INTENT(IN)     :: grp
    INTEGER,                      INTENT(IN)     :: foam_idx
    TYPE(PARMIO_LUT_iVar_type),   INTENT(IN)     :: iVar
    REAL(fp),                     INTENT(IN OUT) :: foam_AD
    REAL(fp),                     INTENT(IN OUT) :: wfreq_AD(2), wth_AD(2)
    REAL(fp),                     INTENT(IN OUT) :: wu_AD(2), ws_AD(2), wq_AD(2)
    INTEGER :: ifreq(2), ith(2), iu(2), isst(2), iq(2)
    INTEGER :: af, at, au, asst, aq, nq
    REAL(fp) :: wf(2), wth(2), wu(2), ws(2), wq(2)
    REAL(fp) :: cw_AD

    CALL Axis_Weights(iVar%B_Frequency, ifreq, wf)
    CALL Axis_Weights(iVar%B_Theta,     ith,   wth)
    CALL Axis_Weights(iVar%B_Wind,      iu,    wu)
    CALL Axis_Weights(iVar%B_SST,       isst,  ws)
    CALL Axis_Weights(iVar%B_SSS,       iq,    wq)
    nq = 1
    IF (iVar%SSS_Active) nq = 2

    DO aq = 1, nq
      DO af = 1, 2
        DO at = 1, 2
          DO au = 1, 2
            DO asst = 1, 2
              cw_AD = REAL(grp%Foam(foam_idx, iq(aq), isst(asst), &
                                     iu(au), ith(at), ifreq(af)), fp) * foam_AD
              wq_AD(aq)    = wq_AD(aq)    + wf(af)*wth(at)*wu(au)*ws(asst) * cw_AD
              wfreq_AD(af) = wfreq_AD(af) + wq(aq)*wth(at)*wu(au)*ws(asst) * cw_AD
              wth_AD(at)   = wth_AD(at)   + wq(aq)*wf(af)*wu(au)*ws(asst)  * cw_AD
              wu_AD(au)    = wu_AD(au)    + wq(aq)*wf(af)*wth(at)*ws(asst) * cw_AD
              ws_AD(asst)  = ws_AD(asst)  + wq(aq)*wf(af)*wth(at)*wu(au)   * cw_AD
            END DO
          END DO
        END DO
      END DO
    END DO
    foam_AD = 0.0_fp
  END SUBROUTINE Interp_Foam_AD


  SUBROUTINE Axis_Weights(b, idx, w)
    TYPE(Bracket_1D_type), INTENT(IN)  :: b
    INTEGER,               INTENT(OUT) :: idx(2)
    REAL(fp),              INTENT(OUT) :: w(2)
    idx(1) = b%lo
    idx(2) = b%hi
    w(2) = b%w
    w(1) = 1.0_fp - w(2)
  END SUBROUTINE Axis_Weights


  SUBROUTINE Axis_Weights_TL(axis, b, query_TL, w_TL)
    REAL(Double),          INTENT(IN)  :: axis(:)
    TYPE(Bracket_1D_type), INTENT(IN)  :: b
    REAL(fp),              INTENT(IN)  :: query_TL
    REAL(fp),              INTENT(OUT) :: w_TL(2)
    REAL(fp) :: dw
    w_TL = 0.0_fp
    IF (b%hi > b%lo) THEN
      dw = query_TL / (REAL(axis(b%hi), fp) - REAL(axis(b%lo), fp))
      w_TL(1) = -dw
      w_TL(2) =  dw
    END IF
  END SUBROUTINE Axis_Weights_TL


  SUBROUTINE Axis_Value_AD(axis, b, w_AD, query_AD)
    REAL(Double),          INTENT(IN)     :: axis(:)
    TYPE(Bracket_1D_type), INTENT(IN)     :: b
    REAL(fp),              INTENT(IN)     :: w_AD(2)
    REAL(fp),              INTENT(IN OUT) :: query_AD
    IF (b%hi > b%lo) THEN
      query_AD = query_AD + (w_AD(2) - w_AD(1)) / &
        (REAL(axis(b%hi), fp) - REAL(axis(b%lo), fp))
    END IF
  END SUBROUTINE Axis_Value_AD

END MODULE PARMIO_LUT_Interpolation
