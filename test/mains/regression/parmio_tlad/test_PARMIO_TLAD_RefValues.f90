!
! test_PARMIO_TLAD_RefValues
!
! One-off helper that prints Compute_PARMIO's emissivity at the two
! grid points used by test_PARMIO_TLAD's Check_Forward_Grid_Point cases.
! Use to refresh the hardcoded expected_e values when the LUT changes.
!
PROGRAM test_PARMIO_TLAD_RefValues
  USE Type_Kinds, ONLY: fp
  USE Message_Handler, ONLY: SUCCESS
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy, PARMIOC
  USE CRTM_PARMIO, ONLY: Compute_PARMIO, PARMIO_iVar_type => iVar_type
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  IMPLICIT NONE
  CHARACTER(512) :: lut
  INTEGER :: err
  REAL(fp) :: e(4), r(4)
  TYPE(PARMIO_iVar_type) :: iv

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, lut)
  ELSE
    lut = '../parmio/Outputs/sweep/lut_production/PARMIO.production.MWwater.EmisCoeff.nc'
  END IF

  err = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err /= SUCCESS) ERROR STOP 'FASTEM load failed'
  err = CRTM_PARMIOCoeff_Load(TRIM(lut), Quiet=.TRUE.)
  IF (err /= SUCCESS) ERROR STOP 'PARMIO LUT load failed'

  CALL Compute_PARMIO(PARMIOC, 6.9_fp, 1, 40.0_fp, 288.15_fp, 35.0_fp, 7.0_fp, &
                      iv, e, r, Azimuth_Angle=37.0_fp)
  WRITE(*,'("sss_active_grid expected_e:")')
  WRITE(*,'(4(es25.17,",",/))') e

  CALL Compute_PARMIO(PARMIOC, 89.0_fp, 1, 45.0_fp, 288.15_fp, 35.0_fp, 10.0_fp, &
                      iv, e, r, Azimuth_Angle=-22.0_fp)
  WRITE(*,'("sss_nominal_h_grid expected_e:")')
  WRITE(*,'(4(es25.17,",",/))') e

  CALL CRTM_PARMIOCoeff_Destroy()
END PROGRAM test_PARMIO_TLAD_RefValues
