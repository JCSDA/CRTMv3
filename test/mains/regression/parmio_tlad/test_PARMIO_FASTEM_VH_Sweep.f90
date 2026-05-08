!
! test_PARMIO_FASTEM_VH_Sweep
!
! Sensor-agnostic V-pol vs H-pol surface-emissivity comparison between
! CRTM PARMIO (LUT) and CRTM FASTEM-6 at four frequencies that matter
! for current and near-future polarimetric MW imagers/sounders:
! 89, 166, 183, and 325 GHz. Sweeps a fine SST x U10 x theta grid and
! writes a CSV the plotter consumes.
!
! Build-only (not registered with CTest); driven externally.
!

PROGRAM test_PARMIO_FASTEM_VH_Sweep

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS
  USE CRTM_PARMIOCoeff,       ONLY: CRTM_PARMIOCoeff_Load,    &
                                    CRTM_PARMIOCoeff_Destroy, &
                                    PARMIOC
  USE CRTM_PARMIO,            ONLY: Compute_PARMIO,           &
                                    PARMIO_iVar_type => iVar_type
  USE CRTM_FastemX,           ONLY: Compute_FastemX,          &
                                    FastemX_iVar_type => iVar_type
  USE CRTM_MWwaterCoeff,      ONLY: CRTM_MWwaterCoeff_Load_FASTEM, MWwaterC

  IMPLICIT NONE

  CHARACTER(512) :: lut_file, csv_file
  INTEGER :: err_stat
  INTEGER :: ifreq, isst, iu, ith, csv_unit
  REAL(fp) :: freq, sst, u10, theta, sss
  REAL(fp) :: emis_p(4), refl_p(4), emis_f(4), refl_f(4)
  TYPE(PARMIO_iVar_type)   :: cp_ivar
  TYPE(FastemX_iVar_type)  :: cf_ivar

  REAL(fp), PARAMETER :: PROBE_FREQS(4) = (/ 89.0_fp, 166.0_fp, 183.31_fp, 325.0_fp /)
  REAL(fp), PARAMETER :: PROBE_SSTS(6)  = &
    (/ 273.15_fp, 278.15_fp, 283.15_fp, 288.15_fp, 293.15_fp, 298.15_fp /)
  REAL(fp), PARAMETER :: PROBE_U10S(7)  = &
    (/ 1.0_fp, 3.0_fp, 5.0_fp, 7.0_fp, 10.0_fp, 15.0_fp, 20.0_fp /)
  REAL(fp), PARAMETER :: PROBE_THETAS(3) = (/ 30.0_fp, 45.0_fp, 55.0_fp /)
  REAL(fp), PARAMETER :: FIXED_SSS = 35.0_fp

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, lut_file)
  ELSE
    lut_file = '../parmio/Outputs/sweep/lut_production/PARMIO.production.MWwater.EmisCoeff.nc'
  END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, csv_file)
  ELSE
    csv_file = 'parmio_fastem_vh_sweep.csv'
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load FASTEM6'
  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load PARMIO LUT'

  sss = FIXED_SSS

  OPEN(NEWUNIT=csv_unit, FILE=TRIM(csv_file), STATUS='REPLACE', ACTION='WRITE')
  WRITE(csv_unit,'(a)') 'freq_ghz,sst_k,u10_mps,theta_deg,sss_psu,' // &
                        'emis_v_parmio,emis_h_parmio,' // &
                        'emis_v_fastem,emis_h_fastem'

  WRITE(*,'("PARMIO/FASTEM V/H sweep: ",i0," frequencies x ",i0," SSTs x ",i0, &
            " U10s x ",i0," thetas = ",i0," points")') &
    SIZE(PROBE_FREQS), SIZE(PROBE_SSTS), SIZE(PROBE_U10S), SIZE(PROBE_THETAS), &
    SIZE(PROBE_FREQS) * SIZE(PROBE_SSTS) * SIZE(PROBE_U10S) * SIZE(PROBE_THETAS)

  DO ifreq = 1, SIZE(PROBE_FREQS)
    freq = PROBE_FREQS(ifreq)
    DO isst = 1, SIZE(PROBE_SSTS)
      sst = PROBE_SSTS(isst)
      DO iu = 1, SIZE(PROBE_U10S)
        u10 = PROBE_U10S(iu)
        DO ith = 1, SIZE(PROBE_THETAS)
          theta = PROBE_THETAS(ith)

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

          WRITE(csv_unit,'(f0.4,",",f0.4,",",f0.4,",",f0.4,",",f0.4,",", &
                           f0.6,",",f0.6,",",f0.6,",",f0.6)') &
            freq, sst, u10, theta, sss, &
            emis_p(1), emis_p(2), emis_f(1), emis_f(2)
        END DO
      END DO
    END DO
  END DO

  CLOSE(csv_unit)
  WRITE(*,'("Wrote ",a)') TRIM(csv_file)

  CALL CRTM_PARMIOCoeff_Destroy()

END PROGRAM test_PARMIO_FASTEM_VH_Sweep
