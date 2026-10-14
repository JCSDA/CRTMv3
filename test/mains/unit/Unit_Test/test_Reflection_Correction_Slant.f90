!
! test_Reflection_Correction_Slant
!
! Validates the transmittance convention of the FASTEM reflection correction
! (Reflection_Correction_Module), used by FASTEM-4/5/6 and PARMIO.
!
! The regression is defined for the slant-path surface-to-space transmittance:
! its optical-depth predictor, -LOG(T)*cos_z, is the vertical optical depth
! only when T is the slant value, which is what RTTOV passes. CRTM's model
! modules pass the vertical transmittance EXP(-SUM(layer optical depth)), so
! the routine converts it to the slant value at the stream's angle. Before
! that conversion the correction was evaluated with the vertical value, which
! left CRTM 0.4 to 1.7 K warmer than RTTOV in off-nadir ocean window channels
! and agreed only at nadir.
!
! Checks, with synthetic regression coefficients (no coefficient files):
!   1. Forward: Reflection_Correction(T_vertical, cos_z) equals an independent
!      evaluation of the regression at T_slant = T_vertical**(1/cos_z); at
!      nadir this is the pre-fix result.
!   2. TL against a central finite difference of the forward, in transmittance
!      and wind speed.
!   3. AD against TL (dot-product test).
!   4. Below the transmittance floor the transmittance derivative is zero in
!      both TL and AD.
! Exit status: STOP 0 = success, STOP 1 = failure.
!
PROGRAM test_Reflection_Correction_Slant

  USE Type_Kinds,      ONLY: fp
  USE FitCoeff_Define, ONLY: FitCoeff_3D_type, FitCoeff_Create, FitCoeff_Destroy
  USE Slope_Variance,  ONLY: svVar_type => iVar_type, Compute_Slope_Variance
  USE CRTM_Parameters, ONLY: LIMIT_EXP
  USE Reflection_Correction_Module, ONLY: rcVar_type => iVar_type, &
                                          Reflection_Correction, &
                                          Reflection_Correction_TL, &
                                          Reflection_Correction_AD
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_Reflection_Correction_Slant'
  INTEGER,  PARAMETER :: N_PREDICTORS = 7, N_TERMS = 3
  REAL(fp), PARAMETER :: FREQ(3)  = [23.8_fp, 89.0_fp, 183.31_fp]
  REAL(fp), PARAMETER :: ZEN(4)   = [0.0_fp, 30.0_fp, 50.0_fp, 55.0_fp]
  REAL(fp), PARAMETER :: TV(3)    = [0.30_fp, 0.70_fp, 0.95_fp]
  REAL(fp), PARAMETER :: WIND(2)  = [2.0_fp, 10.0_fp]
  REAL(fp), PARAMETER :: DEG2RAD  = 3.141592653589793238462643_fp / 180.0_fp
  REAL(fp), PARAMETER :: TOL_FWD  = 1.0e-12_fp, TOL_FD = 1.0e-5_fp, TOL_DOT = 1.0e-10_fp

  TYPE(FitCoeff_3D_type) :: RC
  TYPE(rcVar_type) :: iVar
  REAL(fp) :: cos_z, Tv_in, Ts, w, f
  REAL(fp) :: Rv, Rh, Rv_ref, Rh_ref, Rv_p, Rh_p, Rv_m, Rh_m
  REAL(fp) :: Rv_TL, Rh_TL, Rv_AD, Rh_AD, T_AD, W_AD, h, fd, lhs, rhs
  INTEGER  :: i, k, p, jf, jz, jt, jw, n_fail, n_case

  ! Synthetic regression coefficients of FASTEM's magnitude and sign pattern
  CALL FitCoeff_Create(RC, [N_TERMS, N_PREDICTORS, 2])
  DO p = 1, 2
    DO i = 1, N_PREDICTORS
      DO k = 1, N_TERMS
        RC%C(k,i,p) = 0.02_fp * (-1.0_fp)**(i+k) * (1.0_fp + 0.1_fp*p) / REAL(i, fp)
      END DO
    END DO
  END DO

  n_fail = 0
  n_case = 0
  DO jf = 1, SIZE(FREQ)
    f = FREQ(jf)
    DO jz = 1, SIZE(ZEN)
      cos_z = COS(ZEN(jz)*DEG2RAD)
      DO jt = 1, SIZE(TV)
        Tv_in = TV(jt)
        DO jw = 1, SIZE(WIND)
          w = WIND(jw)
          n_case = n_case + 1

          ! 1. Forward against the regression evaluated at the slant value
          CALL Reflection_Correction(RC, f, cos_z, w, Tv_in, Rv, Rh, iVar)
          Ts = Tv_in**(1.0_fp/cos_z)
          CALL Reference(f, cos_z, w, Ts, Rv_ref, Rh_ref)
          CALL Check('forward Rv', Rv, Rv_ref, TOL_FWD)
          CALL Check('forward Rh', Rh, Rh_ref, TOL_FWD)

          ! 2a. TL in transmittance against a central difference
          h = 1.0e-6_fp * Tv_in
          CALL Reflection_Correction(RC, f, cos_z, w, Tv_in+h, Rv_p, Rh_p, iVar)
          CALL Reflection_Correction(RC, f, cos_z, w, Tv_in-h, Rv_m, Rh_m, iVar)
          CALL Reflection_Correction(RC, f, cos_z, w, Tv_in, Rv, Rh, iVar)
          CALL Reflection_Correction_TL(RC, 0.0_fp, 1.0_fp, Rv_TL, Rh_TL, iVar)
          fd = (Rv_p - Rv_m) / (2.0_fp*h);  CALL Check('TL dRv/dT', Rv_TL, fd, TOL_FD)
          fd = (Rh_p - Rh_m) / (2.0_fp*h);  CALL Check('TL dRh/dT', Rh_TL, fd, TOL_FD)

          ! 2b. TL in wind speed against a central difference
          h = 1.0e-4_fp
          CALL Reflection_Correction(RC, f, cos_z, w+h, Tv_in, Rv_p, Rh_p, iVar)
          CALL Reflection_Correction(RC, f, cos_z, w-h, Tv_in, Rv_m, Rh_m, iVar)
          CALL Reflection_Correction(RC, f, cos_z, w, Tv_in, Rv, Rh, iVar)
          CALL Reflection_Correction_TL(RC, 1.0_fp, 0.0_fp, Rv_TL, Rh_TL, iVar)
          fd = (Rv_p - Rv_m) / (2.0_fp*h);  CALL Check('TL dRv/dU', Rv_TL, fd, TOL_FD)
          fd = (Rh_p - Rh_m) / (2.0_fp*h);  CALL Check('TL dRh/dU', Rh_TL, fd, TOL_FD)

          ! 3. Dot-product test: <TL(dU,dT), (a,b)> = <(dU,dT), AD(a,b)>
          CALL Reflection_Correction_TL(RC, 0.3_fp, 0.7_fp, Rv_TL, Rh_TL, iVar)
          lhs = 1.3_fp*Rv_TL + 0.6_fp*Rh_TL
          Rv_AD = 1.3_fp;  Rh_AD = 0.6_fp;  T_AD = 0.0_fp;  W_AD = 0.0_fp
          CALL Reflection_Correction_AD(RC, Rv_AD, Rh_AD, W_AD, T_AD, iVar)
          rhs = 0.3_fp*W_AD + 0.7_fp*T_AD
          CALL Check('dot product', lhs, rhs, TOL_DOT)
        END DO
      END DO
    END DO
  END DO

  ! 4. Below the floor the transmittance derivative vanishes
  f     = 89.0_fp
  w     = 5.0_fp
  cos_z = COS(50.0_fp*DEG2RAD)
  Tv_in = EXP(-LIMIT_EXP)
  CALL Reflection_Correction(RC, f, cos_z, w, Tv_in, Rv, Rh, iVar)
  CALL Reflection_Correction_TL(RC, 0.0_fp, 1.0_fp, Rv_TL, Rh_TL, iVar)
  CALL Check('floor TL Rv', Rv_TL, 0.0_fp, TOL_FWD)
  CALL Check('floor TL Rh', Rh_TL, 0.0_fp, TOL_FWD)
  Rv_AD = 1.0_fp;  Rh_AD = 1.0_fp;  T_AD = 0.0_fp;  W_AD = 0.0_fp
  CALL Reflection_Correction_AD(RC, Rv_AD, Rh_AD, W_AD, T_AD, iVar)
  CALL Check('floor AD T', T_AD, 0.0_fp, TOL_FWD)

  CALL FitCoeff_Destroy(RC)
  WRITE(*,'(a,": ",i0," cases, ",i0," failures")') PROGRAM_NAME, n_case, n_fail
  IF ( n_fail > 0 ) STOP 1
  STOP 0

CONTAINS

  ! The FASTEM reflection-correction regression evaluated directly at a
  ! slant-path transmittance (the pre-fix code path with slant input).
  SUBROUTINE Reference(f, cos_z, w, Ts, Rv, Rh)
    REAL(fp), INTENT(IN)  :: f, cos_z, w, Ts
    REAL(fp), INTENT(OUT) :: Rv, Rh
    TYPE(svVar_type) :: sv
    REAL(fp) :: variance, od, odx(2), zx(N_PREDICTORS), rough(2)
    INTEGER  :: i, p
    CALL Compute_Slope_Variance(f, w, sv, variance)
    od     = -LOG(Ts) * cos_z
    odx(1) = LOG(od)
    odx(2) = odx(1)**2
    zx(1) = 1.0_fp
    zx(2) = variance
    zx(4) = 1.0_fp / cos_z
    zx(3) = zx(2) * zx(4)
    zx(5) = zx(3) * zx(3)
    zx(6) = zx(4) * zx(4)
    zx(7) = zx(2) * zx(2)
    DO p = 1, 2
      rough(p) = 1.0_fp
      DO i = 1, N_PREDICTORS
        rough(p) = rough(p) + zx(i) * (RC%C(1,i,p) + odx(1)*RC%C(2,i,p) + odx(2)*RC%C(3,i,p))
      END DO
    END DO
    Rv = (1.0_fp - Ts**rough(1)) / (1.0_fp - Ts)
    Rh = (1.0_fp - Ts**rough(2)) / (1.0_fp - Ts)
  END SUBROUTINE Reference

  SUBROUTINE Check(what, got, want, tol)
    CHARACTER(*), INTENT(IN) :: what
    REAL(fp),     INTENT(IN) :: got, want, tol
    IF ( ABS(got - want) > tol * MAX(1.0_fp, ABS(want)) ) THEN
      n_fail = n_fail + 1
      WRITE(*,'(a,": FAIL ",a," f=",f7.2," cos_z=",f6.4," T=",es10.3," U=",f5.1," got ",es22.14," want ",es22.14)') &
        PROGRAM_NAME, what, f, cos_z, Tv_in, w, got, want
    END IF
  END SUBROUTINE Check

END PROGRAM test_Reflection_Correction_Slant
