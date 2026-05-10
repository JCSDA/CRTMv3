!
! PARMIO_RC_Interpolation
!
! Interpolation kernel for optional PARMIO-native reflection correction.
! The table stores effective V/H reflectivity of downwelling atmospheric
! radiation, Rdown, on the same physical axes as the PARMIO emissivity LUT
! plus a transmittance axis.
!

MODULE PARMIO_RC_Interpolation

  USE Type_Kinds, ONLY: fp, Double
  USE PARMIOCoeff_Define, ONLY: &
    PARMIOCoeff_type, PARMIOCoeff_RC_Group_type, &
    PARMIOCoeff_GroupName_For_Frequency,         &
    PARMIO_RC_V_POL, PARMIO_RC_H_POL,            &
    PARMIO_N_GROUPS

  IMPLICIT NONE
  PRIVATE

  PUBLIC :: PARMIO_RC_iVar_type
  PUBLIC :: PARMIO_RC_Interp_Forward
  PUBLIC :: PARMIO_RC_Interp_TL
  PUBLIC :: PARMIO_RC_Interp_AD

  INTEGER, PARAMETER :: N_POL = 2

  TYPE :: Bracket_1D_type
    INTEGER  :: lo = 1
    INTEGER  :: hi = 1
    REAL(fp) :: w  = 0.0_fp
    LOGICAL  :: clamped_low  = .FALSE.
    LOGICAL  :: clamped_high = .FALSE.
  END TYPE Bracket_1D_type

  TYPE :: PARMIO_RC_iVar_type
    INTEGER :: Group_ID = 0
    LOGICAL :: SSS_Active = .FALSE.
    LOGICAL :: Is_Available = .FALSE.
    TYPE(Bracket_1D_type) :: B_Frequency
    TYPE(Bracket_1D_type) :: B_Theta
    TYPE(Bracket_1D_type) :: B_Wind
    TYPE(Bracket_1D_type) :: B_SST
    TYPE(Bracket_1D_type) :: B_SSS
    TYPE(Bracket_1D_type) :: B_Transmittance
    REAL(fp) :: Foam_Fraction = 0.0_fp
    REAL(fp) :: Rdown_Foam_Off(N_POL) = 0.0_fp
    REAL(fp) :: Rdown_Foam_On (N_POL) = 0.0_fp
  END TYPE PARMIO_RC_iVar_type

CONTAINS

  SUBROUTINE PARMIO_RC_Interp_Forward( &
      LUT, Frequency_GHz, Zenith_Angle_deg, Wind_Speed_mps, &
      SST_C, SSS_psu, Foam_Fraction, Transmittance, Rdown, iVar, Is_Available)
    TYPE(PARMIOCoeff_type),     INTENT(IN)  :: LUT
    REAL(fp),                   INTENT(IN)  :: Frequency_GHz
    REAL(fp),                   INTENT(IN)  :: Zenith_Angle_deg
    REAL(fp),                   INTENT(IN)  :: Wind_Speed_mps
    REAL(fp),                   INTENT(IN)  :: SST_C
    REAL(fp),                   INTENT(IN)  :: SSS_psu
    REAL(fp),                   INTENT(IN)  :: Foam_Fraction
    REAL(fp),                   INTENT(IN)  :: Transmittance
    REAL(fp),                   INTENT(OUT) :: Rdown(N_POL)
    TYPE(PARMIO_RC_iVar_type),  INTENT(OUT) :: iVar
    LOGICAL,                    INTENT(OUT) :: Is_Available
    INTEGER :: g

    Rdown = 0.0_fp
    Is_Available = .FALSE.
    g = PARMIOCoeff_GroupName_For_Frequency( &
          Frequency_GHz, LUT%SSS_Cutoff_GHz, LUT%Permittivity_Switch_GHz)
    iVar%Group_ID = g
    IF (g < 1 .OR. g > PARMIO_N_GROUPS) RETURN
    IF (.NOT. LUT%RC_Group(g)%Is_Allocated) RETURN

    iVar%Is_Available = .TRUE.
    iVar%SSS_Active = LUT%RC_Group(g)%SSS_Axis_Active
    iVar%Foam_Fraction = Foam_Fraction

    CALL Bracket(LUT%RC_Group(g)%Frequency,     Frequency_GHz,    iVar%B_Frequency)
    CALL Bracket(LUT%RC_Group(g)%Theta,         Zenith_Angle_deg, iVar%B_Theta)
    CALL Bracket(LUT%RC_Group(g)%Wind_Speed,    Wind_Speed_mps,   iVar%B_Wind)
    CALL Bracket(LUT%RC_Group(g)%SST,           SST_C,            iVar%B_SST)
    CALL Bracket(LUT%RC_Group(g)%Transmittance, Transmittance,    iVar%B_Transmittance)
    IF (iVar%SSS_Active) THEN
      CALL Bracket(LUT%RC_Group(g)%SSS, SSS_psu, iVar%B_SSS)
    ELSE
      iVar%B_SSS%lo = 1
      iVar%B_SSS%hi = 1
      iVar%B_SSS%w  = 0.0_fp
    END IF

    CALL Interp_Rdown(LUT%RC_Group(g), 1, iVar, iVar%Rdown_Foam_Off)
    CALL Interp_Rdown(LUT%RC_Group(g), 2, iVar, iVar%Rdown_Foam_On)
    Rdown = (1.0_fp - Foam_Fraction) * iVar%Rdown_Foam_Off + &
                       Foam_Fraction  * iVar%Rdown_Foam_On
    Is_Available = .TRUE.
  END SUBROUTINE PARMIO_RC_Interp_Forward


  SUBROUTINE PARMIO_RC_Interp_TL( &
      LUT, Frequency_GHz_TL, Zenith_Angle_deg_TL, Wind_Speed_mps_TL, &
      SST_C_TL, SSS_psu_TL, Foam_Fraction_TL, Transmittance_TL, &
      Rdown_TL, iVar)
    TYPE(PARMIOCoeff_type),    INTENT(IN)  :: LUT
    REAL(fp),                  INTENT(IN)  :: Frequency_GHz_TL
    REAL(fp),                  INTENT(IN)  :: Zenith_Angle_deg_TL
    REAL(fp),                  INTENT(IN)  :: Wind_Speed_mps_TL
    REAL(fp),                  INTENT(IN)  :: SST_C_TL
    REAL(fp),                  INTENT(IN)  :: SSS_psu_TL
    REAL(fp),                  INTENT(IN)  :: Foam_Fraction_TL
    REAL(fp),                  INTENT(IN)  :: Transmittance_TL
    REAL(fp),                  INTENT(OUT) :: Rdown_TL(N_POL)
    TYPE(PARMIO_RC_iVar_type), INTENT(IN)  :: iVar
    INTEGER :: g
    REAL(fp) :: off_TL(N_POL), on_TL(N_POL)

    Rdown_TL = 0.0_fp
    g = iVar%Group_ID
    IF (g < 1 .OR. g > PARMIO_N_GROUPS) RETURN
    IF (.NOT. iVar%Is_Available) RETURN
    IF (.NOT. LUT%RC_Group(g)%Is_Allocated) RETURN

    CALL Interp_Rdown_TL( &
      LUT%RC_Group(g), 1, iVar, Frequency_GHz_TL, Zenith_Angle_deg_TL, &
      Wind_Speed_mps_TL, SST_C_TL, SSS_psu_TL, Transmittance_TL, off_TL)
    CALL Interp_Rdown_TL( &
      LUT%RC_Group(g), 2, iVar, Frequency_GHz_TL, Zenith_Angle_deg_TL, &
      Wind_Speed_mps_TL, SST_C_TL, SSS_psu_TL, Transmittance_TL, on_TL)

    Rdown_TL = (1.0_fp - iVar%Foam_Fraction) * off_TL + &
                          iVar%Foam_Fraction  * on_TL  + &
               Foam_Fraction_TL * (iVar%Rdown_Foam_On - iVar%Rdown_Foam_Off)
  END SUBROUTINE PARMIO_RC_Interp_TL


  SUBROUTINE PARMIO_RC_Interp_AD( &
      LUT, Rdown_AD, iVar, Foam_Fraction_AD, &
      Frequency_GHz_AD, Zenith_Angle_deg_AD, Wind_Speed_mps_AD, &
      SST_C_AD, SSS_psu_AD, Transmittance_AD)
    TYPE(PARMIOCoeff_type),    INTENT(IN)     :: LUT
    REAL(fp),                  INTENT(IN OUT) :: Rdown_AD(N_POL)
    TYPE(PARMIO_RC_iVar_type), INTENT(IN)     :: iVar
    REAL(fp),                  INTENT(IN OUT) :: Foam_Fraction_AD
    REAL(fp),                  INTENT(IN OUT) :: Frequency_GHz_AD
    REAL(fp),                  INTENT(IN OUT) :: Zenith_Angle_deg_AD
    REAL(fp),                  INTENT(IN OUT) :: Wind_Speed_mps_AD
    REAL(fp),                  INTENT(IN OUT) :: SST_C_AD
    REAL(fp),                  INTENT(IN OUT) :: SSS_psu_AD
    REAL(fp),                  INTENT(IN OUT) :: Transmittance_AD
    INTEGER :: g
    REAL(fp) :: off_AD(N_POL), on_AD(N_POL)
    REAL(fp) :: wf_AD(2), wt_AD(2), wu_AD(2), ws_AD(2), wq_AD(2), wx_AD(2)

    g = iVar%Group_ID
    IF (g < 1 .OR. g > PARMIO_N_GROUPS) THEN
      Rdown_AD = 0.0_fp
      RETURN
    END IF
    IF (.NOT. iVar%Is_Available .OR. .NOT. LUT%RC_Group(g)%Is_Allocated) THEN
      Rdown_AD = 0.0_fp
      RETURN
    END IF

    off_AD = (1.0_fp - iVar%Foam_Fraction) * Rdown_AD
    on_AD  =            iVar%Foam_Fraction  * Rdown_AD
    Foam_Fraction_AD = Foam_Fraction_AD + &
      SUM((iVar%Rdown_Foam_On - iVar%Rdown_Foam_Off) * Rdown_AD)
    Rdown_AD = 0.0_fp

    wf_AD = 0.0_fp
    wt_AD = 0.0_fp
    wu_AD = 0.0_fp
    ws_AD = 0.0_fp
    wq_AD = 0.0_fp
    wx_AD = 0.0_fp

    CALL Interp_Rdown_AD(LUT%RC_Group(g), 2, iVar, on_AD,  wf_AD, wt_AD, wu_AD, ws_AD, wq_AD, wx_AD)
    CALL Interp_Rdown_AD(LUT%RC_Group(g), 1, iVar, off_AD, wf_AD, wt_AD, wu_AD, ws_AD, wq_AD, wx_AD)

    CALL Axis_Value_AD(LUT%RC_Group(g)%Frequency,     iVar%B_Frequency,     wf_AD, Frequency_GHz_AD)
    CALL Axis_Value_AD(LUT%RC_Group(g)%Theta,         iVar%B_Theta,         wt_AD, Zenith_Angle_deg_AD)
    CALL Axis_Value_AD(LUT%RC_Group(g)%Wind_Speed,    iVar%B_Wind,          wu_AD, Wind_Speed_mps_AD)
    CALL Axis_Value_AD(LUT%RC_Group(g)%SST,           iVar%B_SST,           ws_AD, SST_C_AD)
    CALL Axis_Value_AD(LUT%RC_Group(g)%Transmittance, iVar%B_Transmittance, wx_AD, Transmittance_AD)
    IF (iVar%SSS_Active) THEN
      CALL Axis_Value_AD(LUT%RC_Group(g)%SSS,         iVar%B_SSS,           wq_AD, SSS_psu_AD)
    END IF
  END SUBROUTINE PARMIO_RC_Interp_AD


  SUBROUTINE Bracket(axis, query, b)
    REAL(Double),          INTENT(IN)  :: axis(:)
    REAL(fp),              INTENT(IN)  :: query
    TYPE(Bracket_1D_type), INTENT(OUT) :: b
    INTEGER :: n, k
    REAL(fp) :: lo, hi
    n = SIZE(axis)
    IF (n == 1) THEN
      b%lo = 1; b%hi = 1; b%w = 0.0_fp
      RETURN
    END IF
    IF (query <= REAL(axis(1), fp)) THEN
      b%lo = 1; b%hi = 1; b%w = 0.0_fp
      b%clamped_low = query < REAL(axis(1), fp)
      RETURN
    END IF
    IF (query >= REAL(axis(n), fp)) THEN
      b%lo = n; b%hi = n; b%w = 0.0_fp
      b%clamped_high = query > REAL(axis(n), fp)
      RETURN
    END IF
    DO k = 1, n - 1
      lo = REAL(axis(k), fp)
      hi = REAL(axis(k + 1), fp)
      IF (query >= lo .AND. query <= hi) THEN
        b%lo = k
        b%hi = k + 1
        b%w = (query - lo) / (hi - lo)
        RETURN
      END IF
    END DO
    b%lo = n; b%hi = n; b%w = 0.0_fp
  END SUBROUTINE Bracket


  SUBROUTINE Axis_Weights(b, idx, w)
    TYPE(Bracket_1D_type), INTENT(IN)  :: b
    INTEGER,               INTENT(OUT) :: idx(2)
    REAL(fp),              INTENT(OUT) :: w(2)
    idx = (/b%lo, b%hi/)
    w = (/1.0_fp - b%w, b%w/)
    IF (b%lo == b%hi) w = (/1.0_fp, 0.0_fp/)
  END SUBROUTINE Axis_Weights


  SUBROUTINE Axis_Weights_TL(axis, b, value_TL, w_TL)
    REAL(Double),          INTENT(IN)  :: axis(:)
    TYPE(Bracket_1D_type), INTENT(IN)  :: b
    REAL(fp),              INTENT(IN)  :: value_TL
    REAL(fp),              INTENT(OUT) :: w_TL(2)
    REAL(fp) :: dw
    w_TL = 0.0_fp
    IF (b%lo == b%hi .OR. b%clamped_low .OR. b%clamped_high) RETURN
    dw = value_TL / REAL(axis(b%hi) - axis(b%lo), fp)
    w_TL = (/-dw, dw/)
  END SUBROUTINE Axis_Weights_TL


  SUBROUTINE Axis_Value_AD(axis, b, w_AD, value_AD)
    REAL(Double),          INTENT(IN)     :: axis(:)
    TYPE(Bracket_1D_type), INTENT(IN)     :: b
    REAL(fp),              INTENT(IN OUT) :: w_AD(2)
    REAL(fp),              INTENT(IN OUT) :: value_AD
    REAL(fp) :: dw_AD
    IF (b%lo /= b%hi .AND. .NOT. b%clamped_low .AND. .NOT. b%clamped_high) THEN
      dw_AD = w_AD(2) - w_AD(1)
      value_AD = value_AD + dw_AD / REAL(axis(b%hi) - axis(b%lo), fp)
    END IF
    w_AD = 0.0_fp
  END SUBROUTINE Axis_Value_AD


  SUBROUTINE Interp_Rdown(grp, foam_idx, iVar, out)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(IN)  :: grp
    INTEGER,                         INTENT(IN)  :: foam_idx
    TYPE(PARMIO_RC_iVar_type),       INTENT(IN)  :: iVar
    REAL(fp),                        INTENT(OUT) :: out(N_POL)
    INTEGER :: ipol, lf, lt, lu, ls, lq, lx
    INTEGER :: idf(2), idt(2), idu(2), ids(2), idq(2), idx(2)
    REAL(fp) :: wf(2), wt(2), wu(2), ws(2), wq(2), wx(2), w

    CALL Axis_Weights(iVar%B_Frequency,     idf, wf)
    CALL Axis_Weights(iVar%B_Theta,         idt, wt)
    CALL Axis_Weights(iVar%B_Wind,          idu, wu)
    CALL Axis_Weights(iVar%B_SST,           ids, ws)
    CALL Axis_Weights(iVar%B_SSS,           idq, wq)
    CALL Axis_Weights(iVar%B_Transmittance, idx, wx)

    out = 0.0_fp
    DO lx = 1, 2; DO lf = 1, 2; DO lt = 1, 2; DO lu = 1, 2; DO ls = 1, 2; DO lq = 1, 2
      w = wx(lx) * wf(lf) * wt(lt) * wu(lu) * ws(ls) * wq(lq)
      IF (w == 0.0_fp) CYCLE
      DO ipol = 1, N_POL
        IF (ipol == PARMIO_RC_V_POL) THEN
          out(ipol) = out(ipol) + w * REAL( &
            grp%Rdown_v(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        ELSE
          out(ipol) = out(ipol) + w * REAL( &
            grp%Rdown_h(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        END IF
      END DO
    END DO; END DO; END DO; END DO; END DO; END DO
  END SUBROUTINE Interp_Rdown


  SUBROUTINE Interp_Rdown_TL( &
      grp, foam_idx, iVar, Frequency_GHz_TL, Zenith_Angle_deg_TL, &
      Wind_Speed_mps_TL, SST_C_TL, SSS_psu_TL, Transmittance_TL, out_TL)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(IN)  :: grp
    INTEGER,                         INTENT(IN)  :: foam_idx
    TYPE(PARMIO_RC_iVar_type),       INTENT(IN)  :: iVar
    REAL(fp),                        INTENT(IN)  :: Frequency_GHz_TL
    REAL(fp),                        INTENT(IN)  :: Zenith_Angle_deg_TL
    REAL(fp),                        INTENT(IN)  :: Wind_Speed_mps_TL
    REAL(fp),                        INTENT(IN)  :: SST_C_TL
    REAL(fp),                        INTENT(IN)  :: SSS_psu_TL
    REAL(fp),                        INTENT(IN)  :: Transmittance_TL
    REAL(fp),                        INTENT(OUT) :: out_TL(N_POL)
    INTEGER :: ipol, lf, lt, lu, ls, lq, lx
    INTEGER :: idf(2), idt(2), idu(2), ids(2), idq(2), idx(2)
    REAL(fp) :: wf(2), wt(2), wu(2), ws(2), wq(2), wx(2)
    REAL(fp) :: wf_TL(2), wt_TL(2), wu_TL(2), ws_TL(2), wq_TL(2), wx_TL(2), w_TL
    REAL(fp) :: value

    CALL Axis_Weights(iVar%B_Frequency,     idf, wf)
    CALL Axis_Weights(iVar%B_Theta,         idt, wt)
    CALL Axis_Weights(iVar%B_Wind,          idu, wu)
    CALL Axis_Weights(iVar%B_SST,           ids, ws)
    CALL Axis_Weights(iVar%B_SSS,           idq, wq)
    CALL Axis_Weights(iVar%B_Transmittance, idx, wx)
    CALL Axis_Weights_TL(grp%Frequency,     iVar%B_Frequency,     Frequency_GHz_TL,    wf_TL)
    CALL Axis_Weights_TL(grp%Theta,         iVar%B_Theta,         Zenith_Angle_deg_TL, wt_TL)
    CALL Axis_Weights_TL(grp%Wind_Speed,    iVar%B_Wind,          Wind_Speed_mps_TL,   wu_TL)
    CALL Axis_Weights_TL(grp%SST,           iVar%B_SST,           SST_C_TL,            ws_TL)
    CALL Axis_Weights_TL(grp%Transmittance, iVar%B_Transmittance, Transmittance_TL,    wx_TL)
    IF (iVar%SSS_Active) THEN
      CALL Axis_Weights_TL(grp%SSS,         iVar%B_SSS,           SSS_psu_TL,          wq_TL)
    ELSE
      wq_TL = 0.0_fp
    END IF

    out_TL = 0.0_fp
    DO lx = 1, 2; DO lf = 1, 2; DO lt = 1, 2; DO lu = 1, 2; DO ls = 1, 2; DO lq = 1, 2
      w_TL = wx_TL(lx) * wf(lf)    * wt(lt)    * wu(lu)    * ws(ls)    * wq(lq) + &
             wx(lx)    * wf_TL(lf) * wt(lt)    * wu(lu)    * ws(ls)    * wq(lq) + &
             wx(lx)    * wf(lf)    * wt_TL(lt) * wu(lu)    * ws(ls)    * wq(lq) + &
             wx(lx)    * wf(lf)    * wt(lt)    * wu_TL(lu) * ws(ls)    * wq(lq) + &
             wx(lx)    * wf(lf)    * wt(lt)    * wu(lu)    * ws_TL(ls) * wq(lq) + &
             wx(lx)    * wf(lf)    * wt(lt)    * wu(lu)    * ws(ls)    * wq_TL(lq)
      IF (w_TL == 0.0_fp) CYCLE
      DO ipol = 1, N_POL
        IF (ipol == PARMIO_RC_V_POL) THEN
          value = REAL(grp%Rdown_v(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        ELSE
          value = REAL(grp%Rdown_h(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        END IF
        out_TL(ipol) = out_TL(ipol) + w_TL * value
      END DO
    END DO; END DO; END DO; END DO; END DO; END DO
  END SUBROUTINE Interp_Rdown_TL


  SUBROUTINE Interp_Rdown_AD(grp, foam_idx, iVar, out_AD, wf_AD, wt_AD, wu_AD, ws_AD, wq_AD, wx_AD)
    TYPE(PARMIOCoeff_RC_Group_type), INTENT(IN)     :: grp
    INTEGER,                         INTENT(IN)     :: foam_idx
    TYPE(PARMIO_RC_iVar_type),       INTENT(IN)     :: iVar
    REAL(fp),                        INTENT(IN OUT) :: out_AD(N_POL)
    REAL(fp),                        INTENT(IN OUT) :: wf_AD(2), wt_AD(2), wu_AD(2)
    REAL(fp),                        INTENT(IN OUT) :: ws_AD(2), wq_AD(2), wx_AD(2)
    INTEGER :: ipol, lf, lt, lu, ls, lq, lx
    INTEGER :: idf(2), idt(2), idu(2), ids(2), idq(2), idx(2)
    REAL(fp) :: wf(2), wt(2), wu(2), ws(2), wq(2), wx(2), value, adj

    CALL Axis_Weights(iVar%B_Frequency,     idf, wf)
    CALL Axis_Weights(iVar%B_Theta,         idt, wt)
    CALL Axis_Weights(iVar%B_Wind,          idu, wu)
    CALL Axis_Weights(iVar%B_SST,           ids, ws)
    CALL Axis_Weights(iVar%B_SSS,           idq, wq)
    CALL Axis_Weights(iVar%B_Transmittance, idx, wx)

    DO lx = 1, 2; DO lf = 1, 2; DO lt = 1, 2; DO lu = 1, 2; DO ls = 1, 2; DO lq = 1, 2
      DO ipol = 1, N_POL
        IF (ipol == PARMIO_RC_V_POL) THEN
          value = REAL(grp%Rdown_v(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        ELSE
          value = REAL(grp%Rdown_h(idx(lx), foam_idx, idq(lq), ids(ls), idu(lu), idt(lt), idf(lf)), fp)
        END IF
        adj = value * out_AD(ipol)
        wx_AD(lx) = wx_AD(lx) + wf(lf) * wt(lt) * wu(lu) * ws(ls) * wq(lq) * adj
        wf_AD(lf) = wf_AD(lf) + wx(lx) * wt(lt) * wu(lu) * ws(ls) * wq(lq) * adj
        wt_AD(lt) = wt_AD(lt) + wx(lx) * wf(lf) * wu(lu) * ws(ls) * wq(lq) * adj
        wu_AD(lu) = wu_AD(lu) + wx(lx) * wf(lf) * wt(lt) * ws(ls) * wq(lq) * adj
        ws_AD(ls) = ws_AD(ls) + wx(lx) * wf(lf) * wt(lt) * wu(lu) * wq(lq) * adj
        wq_AD(lq) = wq_AD(lq) + wx(lx) * wf(lf) * wt(lt) * wu(lu) * ws(ls) * adj
      END DO
    END DO; END DO; END DO; END DO; END DO; END DO
    out_AD = 0.0_fp
  END SUBROUTINE Interp_Rdown_AD

END MODULE PARMIO_RC_Interpolation
