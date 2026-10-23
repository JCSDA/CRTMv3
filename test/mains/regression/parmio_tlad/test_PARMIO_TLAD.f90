!
! test_PARMIO_TLAD
!
! Focused finite-difference and adjoint-consistency test for the PARMIO
! microwave ocean surface emissivity model.
!

PROGRAM test_PARMIO_TLAD

  USE Type_Kinds, ONLY: fp
  USE Message_Handler, ONLY: SUCCESS
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy, PARMIOC
  USE CRTM_PARMIO, ONLY: PARMIO_iVar_type => iVar_type, Compute_PARMIO
  USE CRTM_PARMIO_TL, ONLY: Compute_PARMIO_TL
  USE CRTM_PARMIO_AD, ONLY: Compute_PARMIO_AD

  IMPLICIT NONE

  INTEGER, PARAMETER :: N_STOKES = 4
  REAL(fp), PARAMETER :: FD_STEP = 1.0e-5_fp
  REAL(fp), PARAMETER :: FWD_TOL = 5.0e-8_fp
  REAL(fp), PARAMETER :: TL_TOL  = 2.0e-6_fp
  REAL(fp), PARAMETER :: AD_TOL  = 2.0e-8_fp

  CHARACTER(512) :: lut_file
  INTEGER :: err_stat, nfail

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, lut_file)
  ELSE
    lut_file = './testinput/PARMIO.MWwater.EmisCoeff.nc'
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load FASTEM6 MWwater coefficients'

  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load PARMIO coefficient LUT'

  nfail = 0
  CALL Check_Forward_Grid_Point( &
       label       = 'sss_active_grid', &
       frequency   = 6.9_fp, &
       theta       = 40.0_fp, &
       temperature = 288.15_fp, &
       salinity    = 35.0_fp, &
       wind_speed  = 7.0_fp, &
       azimuth     = 37.0_fp, &
       expected_e  = (/ &
         4.5444981127404382e-01_fp, &
         3.0886369657373791e-01_fp, &
        -2.4314563947960824e-03_fp, &
         4.9242013158089660e-04_fp /), &
       nfail       = nfail)
  ! Label kept as 'sss_nominal_h_grid' for historical continuity. 89 GHz now
  ! routes to sss_nominal_m (the Meissner-permittivity threshold moved to
  ! 200 GHz); expected_e was refreshed when the LUT's m-group rows were
  ! regenerated with the correct Meissner permittivity (previously the m-group
  ! held data computed with the high-frequency tabulated dielectric, mislabelled).
  CALL Check_Forward_Grid_Point( &
       label       = 'sss_nominal_h_grid', &
       frequency   = 89.0_fp, &
       theta       = 45.0_fp, &
       temperature = 288.15_fp, &
       salinity    = 35.0_fp, &
       wind_speed  = 10.0_fp, &
       azimuth     = -22.0_fp, &
       expected_e  = (/ &
         7.1240464322802577e-01_fp, &
         4.9766145521667254e-01_fp, &
         5.0850443213180914e-03_fp, &
        -2.7572192979976429e-04_fp /), &
       nfail       = nfail)

  CALL Run_Case( &
       label       = 'sss_active', &
       frequency   = 6.9_fp, &
       theta       = 41.2_fp, &
       temperature = 289.15_fp, &
       salinity    = 34.2_fp, &
       wind_speed  = 8.4_fp, &
       azimuth     = 37.0_fp, &
       trans       = 0.72_fp, &
       nfail       = nfail)
  CALL Run_Case( &
       label       = 'sss_nominal_h', &
       frequency   = 89.0_fp, &
       theta       = 43.7_fp, &
       temperature = 291.15_fp, &
       salinity    = 35.0_fp, &
       wind_speed  = 11.2_fp, &
       azimuth     = -22.0_fp, &
       trans       = 0.80_fp, &
       nfail       = nfail)
  ! >= 200 GHz: the only LUT group the integrated CRTM dispatcher actually
  ! routes to (CRTM_MW_Water_SfcOptics sends MW-water channels >= 200 GHz to
  ! PARMIO). AWS 325 GHz sideband conditions; previously this group had no
  ! direct kernel FD/adjoint coverage.
  CALL Run_Case( &
       label       = 'sss_nominal_m_325', &
       frequency   = 325.15_fp, &
       theta       = 53.0_fp, &
       temperature = 290.15_fp, &
       salinity    = 34.5_fp, &
       wind_speed  = 9.0_fp, &
       azimuth     = 15.0_fp, &
       trans       = 0.55_fp, &
       nfail       = nfail)

  ! CRTM marks "no sensor azimuth" with an out-of-range sentinel (Geometry
  ! default 999.9). An invalid relative azimuth must behave exactly like an
  ! absent one: azimuthal-mean emissivity (harmonic slots dropped), zero
  ! 3rd/4th-Stokes emissivity and reflectivity.
  CALL Check_Sentinel_Azimuth( &
       label       = 'sentinel_azimuth', &
       frequency   = 325.15_fp, &
       theta       = 53.0_fp, &
       temperature = 290.15_fp, &
       salinity    = 34.5_fp, &
       wind_speed  = 9.0_fp, &
       nfail       = nfail)

  CALL CRTM_PARMIOCoeff_Destroy()

  IF (nfail > 0) THEN
    WRITE(*,'("PARMIO TL/AD consistency failures: ",i0)') nfail
    ERROR STOP 'PARMIO TL/AD consistency test failed'
  END IF

  WRITE(*,'("PARMIO TL/AD consistency test passed")')

CONTAINS

  SUBROUTINE Check_Forward_Grid_Point(label, frequency, theta, temperature, &
                                      salinity, wind_speed, azimuth, &
                                      expected_e, nfail)
    CHARACTER(*), INTENT(IN)     :: label
    REAL(fp),     INTENT(IN)     :: frequency
    REAL(fp),     INTENT(IN)     :: theta
    REAL(fp),     INTENT(IN)     :: temperature
    REAL(fp),     INTENT(IN)     :: salinity
    REAL(fp),     INTENT(IN)     :: wind_speed
    REAL(fp),     INTENT(IN)     :: azimuth
    REAL(fp),     INTENT(IN)     :: expected_e(N_STOKES)
    INTEGER,      INTENT(IN OUT) :: nfail

    TYPE(PARMIO_iVar_type) :: ivar
    REAL(fp) :: e(N_STOKES), r(N_STOKES)
    REAL(fp) :: expected_r(N_STOKES)
    REAL(fp) :: err

    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, temperature, salinity, wind_speed, &
         ivar, e, r, Azimuth_Angle=azimuth)

    expected_r = 1.0_fp - expected_e
    ! 3rd/4th Stokes reflectivity is identically zero (FastemX convention):
    ! the U/circular emissivities are azimuthal harmonics, not (1 - r) pairs.
    expected_r(3:4) = 0.0_fp
    err = MAX(MAXVAL(ABS(e - expected_e)), MAXVAL(ABS(r - expected_r)))
    IF (err > FWD_TOL) THEN
      WRITE(*,'(a,": forward grid-point mismatch: ",es13.5)') TRIM(label), err
      nfail = nfail + 1
    END IF
  END SUBROUTINE Check_Forward_Grid_Point

  SUBROUTINE Check_Sentinel_Azimuth(label, frequency, theta, temperature, &
                                    salinity, wind_speed, nfail)
    CHARACTER(*), INTENT(IN)     :: label
    REAL(fp),     INTENT(IN)     :: frequency
    REAL(fp),     INTENT(IN)     :: theta
    REAL(fp),     INTENT(IN)     :: temperature
    REAL(fp),     INTENT(IN)     :: salinity
    REAL(fp),     INTENT(IN)     :: wind_speed
    INTEGER,      INTENT(IN OUT) :: nfail

    TYPE(PARMIO_iVar_type) :: ivar_inv, ivar_abs
    REAL(fp) :: e_inv(N_STOKES), r_inv(N_STOKES)
    REAL(fp) :: e_abs(N_STOKES), r_abs(N_STOKES)
    REAL(fp) :: err

    ! Wind_Direction(0) - Sensor_Azimuth_Angle(999.9) as the dispatcher forms it
    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, temperature, salinity, wind_speed, &
         ivar_inv, e_inv, r_inv, Azimuth_Angle=-999.9_fp)
    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, temperature, salinity, wind_speed, &
         ivar_abs, e_abs, r_abs)

    err = MAX(MAXVAL(ABS(e_inv - e_abs)), MAXVAL(ABS(r_inv - r_abs)))
    IF (err > FWD_TOL) THEN
      WRITE(*,'(a,": sentinel azimuth differs from absent azimuth: ",es13.5)') &
        TRIM(label), err
      nfail = nfail + 1
    END IF
    err = MAX(MAXVAL(ABS(e_inv(3:4))), MAXVAL(ABS(r_inv(3:4))))
    IF (err > FWD_TOL) THEN
      WRITE(*,'(a,": nonzero 3rd/4th Stokes for sentinel azimuth: ",es13.5)') &
        TRIM(label), err
      nfail = nfail + 1
    END IF
  END SUBROUTINE Check_Sentinel_Azimuth


  SUBROUTINE Run_Case(label, frequency, theta, temperature, salinity, &
                      wind_speed, azimuth, trans, nfail)
    CHARACTER(*), INTENT(IN)     :: label
    REAL(fp),     INTENT(IN)     :: frequency
    REAL(fp),     INTENT(IN)     :: theta
    REAL(fp),     INTENT(IN)     :: temperature
    REAL(fp),     INTENT(IN)     :: salinity
    REAL(fp),     INTENT(IN)     :: wind_speed
    REAL(fp),     INTENT(IN)     :: azimuth
    REAL(fp),     INTENT(IN)     :: trans
    INTEGER,      INTENT(IN OUT) :: nfail

    TYPE(PARMIO_iVar_type) :: ivar, ivar_plus, ivar_minus
    REAL(fp) :: e(N_STOKES), r(N_STOKES)
    REAL(fp) :: e_plus(N_STOKES), r_plus(N_STOKES)
    REAL(fp) :: e_minus(N_STOKES), r_minus(N_STOKES)
    REAL(fp) :: e_tl(N_STOKES), r_tl(N_STOKES)
    REAL(fp) :: e_fd(N_STOKES), r_fd(N_STOKES)
    REAL(fp) :: e_ad(N_STOKES), r_ad(N_STOKES)
    REAL(fp) :: e_seed(N_STOKES), r_seed(N_STOKES)
    REAL(fp) :: temperature_tl, salinity_tl, wind_speed_tl
    REAL(fp) :: azimuth_tl, trans_tl
    REAL(fp) :: temperature_ad, salinity_ad, wind_speed_ad
    REAL(fp) :: azimuth_ad, trans_ad
    REAL(fp) :: output_inner, input_inner, tl_err, ad_err

    temperature_tl = 0.70_fp
    salinity_tl    = -0.40_fp
    wind_speed_tl  = 0.30_fp
    azimuth_tl     = 1.20_fp
    trans_tl       = -0.02_fp

    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, temperature, salinity, wind_speed, &
         ivar, e, r, Azimuth_Angle=azimuth, Transmittance=trans)
    CALL Compute_PARMIO_TL( &
         PARMIOC, temperature_tl, salinity_tl, wind_speed_tl, ivar, &
         e_tl, r_tl, Azimuth_Angle_TL=azimuth_tl, Transmittance_TL=trans_tl)

    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, &
         temperature + FD_STEP*temperature_tl, &
         salinity    + FD_STEP*salinity_tl, &
         wind_speed  + FD_STEP*wind_speed_tl, &
         ivar_plus, e_plus, r_plus, &
         Azimuth_Angle=azimuth + FD_STEP*azimuth_tl, &
         Transmittance=trans + FD_STEP*trans_tl)
    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, theta, &
         temperature - FD_STEP*temperature_tl, &
         salinity    - FD_STEP*salinity_tl, &
         wind_speed  - FD_STEP*wind_speed_tl, &
         ivar_minus, e_minus, r_minus, &
         Azimuth_Angle=azimuth - FD_STEP*azimuth_tl, &
         Transmittance=trans - FD_STEP*trans_tl)

    e_fd = (e_plus - e_minus) / (2.0_fp * FD_STEP)
    r_fd = (r_plus - r_minus) / (2.0_fp * FD_STEP)
    tl_err = MAX(MAXVAL(ABS(e_tl - e_fd)), MAXVAL(ABS(r_tl - r_fd)))
    IF (tl_err > TL_TOL) THEN
      WRITE(*,'(a,": TL finite-difference mismatch: ",es13.5)') TRIM(label), tl_err
      nfail = nfail + 1
    END IF

    e_seed = (/ 0.70_fp, -1.10_fp, 0.30_fp, -0.20_fp /)
    r_seed = (/ 0.40_fp, -0.50_fp, 0.25_fp, -0.15_fp /)
    e_ad = e_seed
    r_ad = r_seed
    temperature_ad = 0.0_fp
    salinity_ad    = 0.0_fp
    wind_speed_ad  = 0.0_fp
    azimuth_ad     = 0.0_fp
    trans_ad       = 0.0_fp

    output_inner = SUM(e_seed * e_tl) + SUM(r_seed * r_tl)
    CALL Compute_PARMIO_AD( &
         PARMIOC, e_ad, r_ad, ivar, &
         temperature_ad, salinity_ad, wind_speed_ad, &
         Azimuth_Angle_AD=azimuth_ad, Transmittance_AD=trans_ad)
    input_inner = temperature_ad * temperature_tl &
                + salinity_ad    * salinity_tl &
                + wind_speed_ad  * wind_speed_tl &
                + azimuth_ad     * azimuth_tl &
                + trans_ad       * trans_tl
    ad_err = ABS(output_inner - input_inner)
    IF (ad_err > AD_TOL) THEN
      WRITE(*,'(a,": AD inner-product mismatch: ",es13.5)') TRIM(label), ad_err
      nfail = nfail + 1
    END IF
  END SUBROUTINE Run_Case

END PROGRAM test_PARMIO_TLAD
