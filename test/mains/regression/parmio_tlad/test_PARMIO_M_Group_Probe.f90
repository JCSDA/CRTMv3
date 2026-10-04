!
! test_PARMIO_M_Group_Probe
!
! Diagnostic probe for the PARMIO sss_nominal_m group. Calls Compute_PARMIO
! and Compute_FastemX at the same query point and writes a CSV that can be
! plotted alongside the standalone PARMIO reference for visual comparison.
!
! Usage: test_PARMIO_M_Group_Probe <lut_file> [csv_file]
!

PROGRAM test_PARMIO_M_Group_Probe

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS
  USE PARMIOCoeff_Define,     ONLY: PARMIOCoeff_Inspect, &
                                    PARMIO_GROUP_SSS_NOMINAL_M, &
                                    N_PARMIO_HARMONIC_TERMS
  USE CRTM_PARMIOCoeff,       ONLY: CRTM_PARMIOCoeff_Load, &
                                    CRTM_PARMIOCoeff_Destroy, &
                                    PARMIOC
  USE CRTM_PARMIO,            ONLY: Compute_PARMIO, &
                                    PARMIO_iVar_type => iVar_type
  USE CRTM_FastemX,           ONLY: Compute_FastemX, &
                                    FastemX_iVar_type => iVar_type
  USE CRTM_MWwaterCoeff,      ONLY: CRTM_MWwaterCoeff_Load_FASTEM, MWwaterC

  IMPLICIT NONE

  CHARACTER(512) :: lut_file, csv_file
  INTEGER :: err_stat
  INTEGER :: ifreq, csv_unit
  REAL(fp) :: freq, theta, sst, u10, sss
  REAL(fp) :: emis_p(4), refl_p(4), emis_f(4), refl_f(4)
  TYPE(PARMIO_iVar_type)   :: cp_ivar
  TYPE(FastemX_iVar_type)  :: cf_ivar

  ! Dense frequency grid spanning the m-group plus the ATMS channel set
  ! (drawn from atms_npp.SpcCoeff). Designed so smooth interpolation curves
  ! plus discrete ATMS markers can be drawn from one CSV.
  REAL(fp), PARAMETER :: PROBE_FREQS(*) = (/ &
       15.0_fp,  18.7_fp,  21.0_fp,  23.8_fp,  25.0_fp,  29.0_fp,  31.4_fp, &
       35.0_fp,  40.0_fp,  45.0_fp,  50.3_fp,  51.76_fp, 52.8_fp,  53.596_fp, &
       54.4_fp,  54.94_fp, 55.5_fp,  57.29_fp, 60.0_fp,  70.0_fp,  80.0_fp, &
       88.2_fp, 100.0_fp, 120.0_fp, 140.0_fp, 165.5_fp, 175.0_fp, 183.31_fp /)

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, lut_file)
  ELSE
    lut_file = './testinput/PARMIO.MWwater.EmisCoeff.nc'
  END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, csv_file)
  ELSE
    csv_file = 'parmio_m_group_probe.csv'
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load FASTEM6'

  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load PARMIO LUT'

  WRITE(*,'(/,"=== PARMIOCoeff_Inspect dump ===")')
  CALL PARMIOCoeff_Inspect(PARMIOC)

  ! Match the conditions in the original parmio_comparison.png plot:
  !   SST = -0.4 C, U10 = 3.7 m/s, theta = 20.8 deg
  u10   = 3.7_fp
  sss   = 35.0_fp
  theta = 20.8_fp
  sst   = 272.75_fp

  OPEN(NEWUNIT=csv_unit, FILE=TRIM(csv_file), STATUS='REPLACE', ACTION='WRITE')
  WRITE(csv_unit,'(a)') 'freq_ghz,sst_k,u10_mps,theta_deg,sss_psu,' // &
                        'emis_v_parmio,emis_h_parmio,' // &
                        'emis_v_fastem,emis_h_fastem'

  WRITE(*,'(/,"=== Probe scan at SST=",f5.1," K, U10=",f4.1," m/s, theta=",f4.1," deg ===")') &
    sst, u10, theta
  WRITE(*,'(a8,1x,a8,1x,a8,1x,a8,1x,a8)') &
    'freq', 'V_PARMIO', 'H_PARMIO', 'V_FASTEM', 'H_FASTEM'

  DO ifreq = 1, SIZE(PROBE_FREQS)
    freq = PROBE_FREQS(ifreq)

    CALL Compute_PARMIO( &
      PARMIOCoeff   = PARMIOC,             &
      Frequency     = freq,                &
      n_Angles      = 1,                   &
      Zenith_Angle  = theta,               &
      Temperature   = sst,                 &
      Salinity      = sss,                 &
      Wind_Speed    = u10,                 &
      iVar          = cp_ivar,             &
      Emissivity    = emis_p,              &
      Reflectivity  = refl_p)

    CALL Compute_FastemX( &
      MWwaterCoeff  = MWwaterC,            &
      Frequency     = freq,                &
      n_Angles      = 1,                   &
      Zenith_Angle  = theta,               &
      Temperature   = sst,                 &
      Salinity      = sss,                 &
      Wind_Speed    = u10,                 &
      iVar          = cf_ivar,             &
      Emissivity    = emis_f,              &
      Reflectivity  = refl_f)

    WRITE(*,'(f8.3,1x,f8.5,1x,f8.5,1x,f8.5,1x,f8.5)') &
      freq, emis_p(1), emis_p(2), emis_f(1), emis_f(2)
    WRITE(csv_unit,'(f0.4,",",f0.4,",",f0.4,",",f0.4,",",f0.4,",",f0.6,",",f0.6,",",f0.6,",",f0.6)') &
      freq, sst, u10, theta, sss, emis_p(1), emis_p(2), emis_f(1), emis_f(2)
  END DO

  CLOSE(csv_unit)
  WRITE(*,'(/,"Wrote ",a)') TRIM(csv_file)

  CALL CRTM_PARMIOCoeff_Destroy()

END PROGRAM test_PARMIO_M_Group_Probe
