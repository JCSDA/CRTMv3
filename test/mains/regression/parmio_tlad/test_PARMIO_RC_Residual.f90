!
! test_PARMIO_RC_Residual
!
! Spot-check the Phase-4 PARMIO + FASTEM-RCCoeff reflection correction
! against PARMIO's own atmosphere-on brightness-temperature output.
!

PROGRAM test_PARMIO_RC_Residual

  USE Type_Kinds, ONLY: fp
  USE Message_Handler, ONLY: SUCCESS
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy, PARMIOC
  USE CRTM_PARMIO, ONLY: PARMIO_iVar_type => iVar_type, Compute_PARMIO

  IMPLICIT NONE

  INTEGER, PARAMETER :: N_STOKES = 4
  INTEGER, PARAMETER :: N_CASES = 8
  REAL(fp), PARAMETER :: GATE_TOL_K = 0.5_fp

  CHARACTER(512) :: lut_file
  INTEGER :: err_stat, i
  INTEGER :: n_over_v, n_over_h
  REAL(fp) :: emissivity(N_STOKES), reflectivity(N_STOKES)
  REAL(fp) :: pred_v, pred_h, trans
  REAL(fp) :: residual_v(N_CASES), residual_h(N_CASES)
  REAL(fp) :: max_residual
  TYPE(PARMIO_iVar_type) :: ivar

  CHARACTER(16), PARAMETER :: label(N_CASES) = (/ &
       'c06_th40_u07    ', &
       'c06_th55_u20    ', &
       'c18_th40_u07    ', &
       'c18_th55_u20    ', &
       'c36_th40_u07    ', &
       'c36_th55_u20    ', &
       'c89_th40_u07    ', &
       'c89_th55_u20    ' /)
  REAL(fp), PARAMETER :: frequency(N_CASES) = (/ &
       6.900000000000000e+00_fp, &
       6.900000000000000e+00_fp, &
       1.870000000000000e+01_fp, &
       1.870000000000000e+01_fp, &
       3.650000000000000e+01_fp, &
       3.650000000000000e+01_fp, &
       8.900000000000000e+01_fp, &
       8.900000000000000e+01_fp /)
  REAL(fp), PARAMETER :: theta(N_CASES) = (/ &
       4.000000000000000e+01_fp, &
       5.500000000000000e+01_fp, &
       4.000000000000000e+01_fp, &
       5.500000000000000e+01_fp, &
       4.000000000000000e+01_fp, &
       5.500000000000000e+01_fp, &
       4.000000000000000e+01_fp, &
       5.500000000000000e+01_fp /)
  REAL(fp), PARAMETER :: wind_speed(N_CASES) = (/ &
       7.000000000000000e+00_fp, &
       2.000000000000000e+01_fp, &
       7.000000000000000e+00_fp, &
       2.000000000000000e+01_fp, &
       7.000000000000000e+00_fp, &
       2.000000000000000e+01_fp, &
       7.000000000000000e+00_fp, &
       2.000000000000000e+01_fp /)
  REAL(fp), PARAMETER :: sst_k(N_CASES) = (/ &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp, &
       2.881500000000000e+02_fp /)
  REAL(fp), PARAMETER :: sss(N_CASES) = (/ &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp, &
       3.500000000000000e+01_fp /)
  REAL(fp), PARAMETER :: phi(N_CASES) = 0.0_fp
  REAL(fp), PARAMETER :: tb_up(N_CASES) = (/ &
       3.240000000000000e+00_fp, &
       4.310000000000000e+00_fp, &
       8.199999999999999e+00_fp, &
       1.089000000000000e+01_fp, &
       1.881000000000000e+01_fp, &
       2.477000000000000e+01_fp, &
       3.313000000000000e+01_fp, &
       4.322000000000000e+01_fp /)
  REAL(fp), PARAMETER :: tb_down(N_CASES) = (/ &
       5.900000000000000e+00_fp, &
       6.960000000000000e+00_fp, &
       1.082000000000000e+01_fp, &
       1.348000000000000e+01_fp, &
       2.134000000000000e+01_fp, &
       2.728000000000000e+01_fp, &
       3.556000000000000e+01_fp, &
       4.561000000000000e+01_fp /)
  REAL(fp), PARAMETER :: tau(N_CASES) = (/ &
       1.307000000000000e-02_fp, &
       1.745000000000000e-02_fp, &
       3.287000000000000e-02_fp, &
       4.387000000000000e-02_fp, &
       7.783000000000000e-02_fp, &
       1.039000000000000e-01_fp, &
       1.387000000000000e-01_fp, &
       1.852000000000000e-01_fp /)
  REAL(fp), PARAMETER :: target_v(N_CASES) = (/ &
       1.364552450000000e+02_fp, &
       1.670052420000000e+02_fp, &
       1.620781040000000e+02_fp, &
       1.916080430000000e+02_fp, &
       1.903238910000000e+02_fp, &
       2.155549790000000e+02_fp, &
       2.341035940000000e+02_fp, &
       2.495075990000000e+02_fp /)
  REAL(fp), PARAMETER :: target_h(N_CASES) = (/ &
       9.567352400000001e+01_fp, &
       9.125832299999999e+01_fp, &
       1.235214280000000e+02_fp, &
       1.304386340000000e+02_fp, &
       1.530302360000000e+02_fp, &
       1.643161350000000e+02_fp, &
       2.103029000000000e+02_fp, &
       2.250442120000000e+02_fp /)

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, lut_file)
  ELSE
    lut_file = '../parmio/Outputs/sweep/lut_production/PARMIO.production.MWwater.EmisCoeff.nc'
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load FASTEM6 MWwater coefficients'

  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load PARMIO coefficient LUT'

  n_over_v = 0
  n_over_h = 0
  WRITE(*,'("PARMIO atmosphere-on residual gate")')
  WRITE(*,'("case              freq theta wind  res_v(K)  res_h(K)")')
  DO i = 1, N_CASES
    trans = EXP(-tau(i))
    CALL Compute_PARMIO( &
         PARMIOC, frequency(i), 1, theta(i), sst_k(i), sss(i), wind_speed(i), &
         ivar, emissivity, reflectivity, &
         Azimuth_Angle = phi(i), Transmittance = trans)
    pred_v = tb_up(i) + trans * (sst_k(i) * emissivity(1) + tb_down(i) * reflectivity(1))
    pred_h = tb_up(i) + trans * (sst_k(i) * emissivity(2) + tb_down(i) * reflectivity(2))
    residual_v(i) = pred_v - target_v(i)
    residual_h(i) = pred_h - target_h(i)
    IF (ABS(residual_v(i)) > GATE_TOL_K) n_over_v = n_over_v + 1
    IF (ABS(residual_h(i)) > GATE_TOL_K) n_over_h = n_over_h + 1
    WRITE(*,'(a16,1x,f5.1,1x,f5.1,1x,f5.1,2(1x,f9.3))') &
      label(i), frequency(i), theta(i), wind_speed(i), residual_v(i), residual_h(i)
  END DO

  max_residual = MAX(MAXVAL(ABS(residual_v)), MAXVAL(ABS(residual_h)))
  WRITE(*,'("max_abs_residual(K): ",f9.3)') max_residual
  WRITE(*,'("cases_over_0.5K: V=",i0,"/",i0," H=",i0,"/",i0)') n_over_v, N_CASES, n_over_h, N_CASES

  CALL CRTM_PARMIOCoeff_Destroy()

  IF (max_residual > GATE_TOL_K) THEN
    ERROR STOP 'PARMIO+FASTEM-RCCoeff atmosphere-on residual exceeds 0.5 K gate'
  END IF
END PROGRAM test_PARMIO_RC_Residual
