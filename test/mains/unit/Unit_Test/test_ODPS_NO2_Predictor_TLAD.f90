!
! test_ODPS_NO2_Predictor_TLAD
!
! Machine-precision TL/AD transpose check for the GROUP_UV_NO2 (Group_Index=8)
! ODPS predictor mapping, directly at the ODPS_Compute_Predictor level with no
! radiative transfer in the loop.
!
! Rationale: the RT-level adjoint dot-product in test_UV_NO2_TLAD is bounded
! by float64 accumulation roundoff of the UV solar scattering solver (~1e-11
! relative), which cannot distinguish a tiny transpose slip in the new NO2
! predictor block from RT summation noise. Here the dot-product
!     <X_TL, X_AD>  ==  <(dT, dAbsorber), (T_AD, Absorber_AD)>
! covers the complete group-8 predictor mapping (all 6 components, including
! the new 3-predictor NO2 block and its DT/DT2 temperature coupling) and must
! hold to machine epsilon: any real transcription error fails by many orders.
!
! Inputs are synthetic but physical (positive absorbers, realistic T range);
! duals are deterministic pseudo-random. No coefficient files are required.
!
! Exit: STOP 0 on pass, STOP 1 otherwise.
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson, 25-Jul-2026
!
PROGRAM test_ODPS_NO2_Predictor_TLAD

  USE Type_Kinds           , ONLY: fp
  USE ODPS_Predictor_Define, ONLY: ODPS_Predictor_type, &
                                   ODPS_Predictor_Create, &
                                   ODPS_Predictor_Destroy, &
                                   ODPS_Predictor_Associated, &
                                   PAFV_Create, &
                                   PAFV_Associated
  USE ODPS_Predictor       , ONLY: GROUP_UV_NO2, &
                                   ODPS_Compute_Predictor, &
                                   ODPS_Compute_Predictor_TL, &
                                   ODPS_Compute_Predictor_AD
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_ODPS_NO2_Predictor_TLAD'
  INTEGER,  PARAMETER :: N_LAYERS     = 100
  INTEGER,  PARAMETER :: N_COMPONENTS = 6    ! group-8 [Dry,WLO,WCO,OZO,CO2,NO2]
  INTEGER,  PARAMETER :: N_ABSORBERS  = 4    ! group-8 [H2O,O3,CO2,NO2]
  INTEGER,  PARAMETER :: MAX_N_PRED   = 15
  ! Canonical group-8 rosters [Dry,WLO,WCO,OZO,CO2,NO2] over [H2O,O3,CO2,NO2]
  INTEGER,  PARAMETER :: COMPONT_IDS(N_COMPONENTS) = (/ 20, 101, 15, 114, 121, 122 /)
  INTEGER,  PARAMETER :: ABSORBR_IDS(N_ABSORBERS)  = (/ 1, 3, 2, 10 /)
  REAL(fp), PARAMETER :: ZERO = 0.0_fp, ONE = 1.0_fp
  REAL(fp), PARAMETER :: TOL  = 1.0e-13_fp   ! ~1e3 * eps: pure-arithmetic bound

  TYPE(ODPS_Predictor_type) :: prd, prd_TL, prd_AD
  REAL(fp) :: ref_level_p(N_LAYERS+1)
  REAL(fp) :: ref_t(N_LAYERS), t(N_LAYERS), t_tl(N_LAYERS), t_ad(N_LAYERS)
  REAL(fp) :: ref_abs(N_LAYERS,N_ABSORBERS)
  REAL(fp) :: abs_prof(N_LAYERS,N_ABSORBERS)
  REAL(fp) :: abs_tl(N_LAYERS,N_ABSORBERS), abs_ad(N_LAYERS,N_ABSORBERS)
  REAL(fp) :: secang(N_LAYERS)
  REAL(fp) :: x_ad_save(N_LAYERS,MAX_N_PRED,N_COMPONENTS)
  REAL(fp) :: lhs, rhs, rel
  INTEGER  :: k, i, j

  ! Synthetic column: 101-level log grid 0.005..1100 hPa, realistic T,
  ! positive absorber profiles with layer-scale structure.
  DO k = 1, N_LAYERS+1
    ref_level_p(k) = 0.005_fp * (1100.0_fp/0.005_fp)**(REAL(k-1,fp)/REAL(N_LAYERS,fp))
  END DO
  DO k = 1, N_LAYERS
    ref_t(k)      = 230.0_fp + 60.0_fp*SIN( 0.06_fp*REAL(k,fp) )
    t(k)          = ref_t(k) + 12.0_fp*SIN( 0.21_fp*REAL(k,fp) + 0.5_fp )
    ref_abs(k,1)  = 1.0e-3_fp + 8.0_fp*EXP( -REAL(N_LAYERS-k,fp)/12.0_fp )      ! H2O
    ref_abs(k,2)  = 0.03_fp + 7.0_fp*EXP( -((REAL(k,fp)-35.0_fp)/12.0_fp)**2 )  ! O3
    ref_abs(k,3)  = 380.0_fp + 3.0_fp*SIN( 0.1_fp*REAL(k,fp) )                  ! CO2
    ref_abs(k,4)  = 3.0e-4_fp + 6.7e-3_fp*EXP( -((REAL(k,fp)-15.0_fp)/8.0_fp)**2 ) ! NO2
    DO j = 1, N_ABSORBERS
      abs_prof(k,j) = ref_abs(k,j) * ( ONE + 0.35_fp*SIN( 0.17_fp*REAL(k,fp) + 0.9_fp*REAL(j,fp) ) )
    END DO
    secang(k) = 1.5_fp + 0.3_fp*SIN( 0.05_fp*REAL(k,fp) )
  END DO

  CALL ODPS_Predictor_Create( prd,    N_LAYERS, N_LAYERS, N_COMPONENTS, MAX_N_PRED, No_OPTRAN=.TRUE. )
  CALL ODPS_Predictor_Create( prd_TL, N_LAYERS, N_LAYERS, N_COMPONENTS, MAX_N_PRED, No_OPTRAN=.TRUE. )
  CALL ODPS_Predictor_Create( prd_AD, N_LAYERS, N_LAYERS, N_COMPONENTS, MAX_N_PRED, No_OPTRAN=.TRUE. )
  IF ( .NOT. ( ODPS_Predictor_Associated(prd) .AND. &
               ODPS_Predictor_Associated(prd_TL) .AND. &
               ODPS_Predictor_Associated(prd_AD) ) ) THEN
    WRITE(*,*) 'Predictor allocation failed'
    STOP 1
  END IF
  ! The TL/AD predictor routines read the forward-saved integrated variables
  ! (Tz_ref, GAzp_ref, PDP, ...) from Predictor%PAFV, so the forward pass must
  ! run with PAFV allocated (mirrors CRTM_Predictor_Define with SaveFWV).
  CALL PAFV_Create( prd%PAFV, N_LAYERS, N_LAYERS, N_ABSORBERS, No_OPTRAN=.TRUE. )
  IF ( .NOT. PAFV_Associated(prd%PAFV) ) THEN
    WRITE(*,*) 'PAFV allocation failed'
    STOP 1
  END IF

  ! Forward (fills prd%n_CP and the predictor values the AD recomputation uses)
  CALL ODPS_Compute_Predictor( GROUP_UV_NO2, COMPONT_IDS, ABSORBR_IDS, &
                               t, abs_prof, ref_level_p, &
                               ref_t, ref_abs, secang, prd )

  ! TL input: relative structure on every absorber, additive on T
  DO k = 1, N_LAYERS
    t_tl(k) = SIN( 0.31_fp*REAL(k,fp) + 0.2_fp )
    DO j = 1, N_ABSORBERS
      abs_tl(k,j) = 0.1_fp * abs_prof(k,j) * COS( 0.23_fp*REAL(k,fp) + 1.1_fp*REAL(j,fp) )
    END DO
  END DO
  CALL ODPS_Compute_Predictor_TL( GROUP_UV_NO2, COMPONT_IDS, ABSORBR_IDS, &
                                  t, abs_prof, ref_t, ref_abs, &
                                  secang, prd, t_tl, abs_tl, prd_TL )

  ! AD dual: deterministic pseudo-random weights on every active predictor slot
  prd_AD%X = ZERO
  x_ad_save = ZERO
  DO j = 1, N_COMPONENTS
    DO i = 1, prd%n_CP(j)
      DO k = 1, N_LAYERS
        x_ad_save(k,i,j) = SIN( 0.13_fp*REAL(k,fp) + 0.7_fp*REAL(i,fp) + 1.7_fp*REAL(j,fp) )
        prd_AD%X(k,i,j)  = x_ad_save(k,i,j)
      END DO
    END DO
  END DO

  t_ad   = ZERO
  abs_ad = ZERO
  CALL ODPS_Compute_Predictor_AD( GROUP_UV_NO2, COMPONT_IDS, ABSORBR_IDS, &
                                  t, abs_prof, ref_t, ref_abs, &
                                  secang, prd, prd_AD, t_ad, abs_ad )

  ! <X_TL, X_AD>  vs  <(dT,dAbs), (T_AD, Abs_AD)>
  lhs = ZERO
  DO j = 1, N_COMPONENTS
    DO i = 1, prd%n_CP(j)
      DO k = 1, N_LAYERS
        lhs = lhs + prd_TL%X(k,i,j) * x_ad_save(k,i,j)
      END DO
    END DO
  END DO
  rhs = DOT_PRODUCT( t_tl, t_ad )
  DO j = 1, N_ABSORBERS
    rhs = rhs + DOT_PRODUCT( abs_tl(:,j), abs_ad(:,j) )
  END DO

  rel = ABS(lhs-rhs) / MAX( ABS(lhs), TINY(ONE) )
  WRITE(*,'(/5x,a)') 'GROUP_UV_NO2 predictor-level TL/AD transpose check'
  WRITE(*,'(7x,"n_CP = ",6(i0,:,", "))') prd%n_CP
  WRITE(*,'(7x,"<X_TL,X_AD> = ",es20.13)') lhs
  WRITE(*,'(7x,"<x,   AD  > = ",es20.13)') rhs
  WRITE(*,'(7x,"relative difference = ",es11.4,"   (tol ",es8.1,")")') rel, TOL

  CALL ODPS_Predictor_Destroy( prd )
  CALL ODPS_Predictor_Destroy( prd_TL )
  CALL ODPS_Predictor_Destroy( prd_AD )

  IF ( rel < TOL ) THEN
    WRITE(*,'(7x,a/)') 'PASS'
    STOP 0
  ELSE
    WRITE(*,'(7x,a/)') 'FAIL'
    STOP 1
  END IF

END PROGRAM test_ODPS_NO2_Predictor_TLAD
