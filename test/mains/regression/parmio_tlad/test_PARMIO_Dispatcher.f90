!
! test_PARMIO_Dispatcher
!
! Exercise the CRTM_MW_Water_SfcOptics PARMIO branch through forward,
! tangent-linear, and adjoint dispatch and compare against direct PARMIO calls.
!

PROGRAM test_PARMIO_Dispatcher

  USE Type_Kinds, ONLY: fp
  USE Message_Handler, ONLY: SUCCESS
  USE CRTM_Surface_Define, ONLY: CRTM_Surface_type, CRTM_Surface_Zero
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                      CRTM_GeometryInfo_SetValue
  USE CRTM_SfcOptics_Define, ONLY: CRTM_SfcOptics_type, &
                                   CRTM_SfcOptics_Create, &
                                   CRTM_SfcOptics_Destroy, &
                                   CRTM_SfcOptics_Zero
  USE CRTM_SpcCoeff, ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy, PARMIOC
  USE CRTM_PARMIO, ONLY: PARMIO_iVar_type => iVar_type, Compute_PARMIO
  USE CRTM_PARMIO_TL, ONLY: Compute_PARMIO_TL
  USE CRTM_PARMIO_AD, ONLY: Compute_PARMIO_AD
  USE CRTM_MW_Water_SfcOptics, ONLY: MW_iVar_type => iVar_type, &
                                     Compute_MW_Water_SfcOptics, &
                                     Compute_MW_Water_SfcOptics_TL, &
                                     Compute_MW_Water_SfcOptics_AD

  IMPLICIT NONE

  INTEGER, PARAMETER :: N_STOKES = 4
  REAL(fp), PARAMETER :: TOL = 1.0e-12_fp

  CHARACTER(512) :: spc_path, lut_file
  CHARACTER(32)  :: sensor_id(1)
  INTEGER :: err_stat, nfail

  IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
    CALL GET_COMMAND_ARGUMENT(1, spc_path)
  ELSE
    spc_path = './testinput/'
  END IF
  IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
    CALL GET_COMMAND_ARGUMENT(2, lut_file)
  ELSE
    lut_file = './testinput/PARMIO.production.MWwater.EmisCoeff.nc'
  END IF

  sensor_id = (/ 'atms_npp' /)
  err_stat = CRTM_SpcCoeff_Load( &
       Sensor_ID = sensor_id, File_Path = TRIM(spc_path), &
       netCDF = .FALSE., Quiet = .TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load SpcCoeff'

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load FASTEM6 MWwater coefficients'

  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) ERROR STOP 'Failed to load PARMIO coefficient LUT'

  nfail = 0
  CALL Check_Dispatcher(nfail)

  CALL CRTM_PARMIOCoeff_Destroy()
  err_stat = CRTM_SpcCoeff_Destroy()

  IF (nfail > 0) THEN
    WRITE(*,'("PARMIO dispatcher failures: ",i0)') nfail
    ERROR STOP 'PARMIO dispatcher test failed'
  END IF
  WRITE(*,'("PARMIO dispatcher test passed")')

CONTAINS

  SUBROUTINE Check_Dispatcher(nfail)
    INTEGER, INTENT(IN OUT) :: nfail

    TYPE(CRTM_Surface_type)      :: surface, surface_tl, surface_ad, surface_ad_ref
    TYPE(CRTM_GeometryInfo_type) :: geometry
    TYPE(CRTM_SfcOptics_type)    :: sfcoptics, sfcoptics_tl, sfcoptics_ad
    TYPE(MW_iVar_type)           :: mw_ivar
    TYPE(PARMIO_iVar_type)       :: parmio_ivar
    REAL(fp) :: e_ref(N_STOKES), r_ref(N_STOKES)
    REAL(fp) :: e_tl_ref(N_STOKES), r_tl_ref(N_STOKES)
    REAL(fp) :: e_ad_ref(N_STOKES), r_ad_ref(N_STOKES)
    REAL(fp) :: azimuth_ad_ref, trans_ad_ref
    REAL(fp) :: frequency
    INTEGER :: channel_index

    channel_index = 1
    frequency = SC(1)%Frequency(channel_index)

    CALL CRTM_Surface_Zero(surface)
    surface%Water_Coverage    = 1.0_fp
    surface%Water_Temperature = 289.15_fp
    surface%Salinity          = 34.5_fp
    surface%Wind_Speed        = 8.0_fp
    surface%Wind_Direction    = 47.0_fp

    CALL CRTM_GeometryInfo_SetValue( &
         geometry, Source_Azimuth_Angle = 0.0_fp, &
         Sensor_Azimuth_Angle = 10.0_fp)

    CALL CRTM_SfcOptics_Create(sfcoptics, 1, N_STOKES)
    CALL CRTM_SfcOptics_Create(sfcoptics_tl, 1, N_STOKES)
    CALL CRTM_SfcOptics_Create(sfcoptics_ad, 1, N_STOKES)
    sfcoptics%Use_PARMIO_Model = .TRUE.
    sfcoptics%Angle(1) = 40.0_fp
    sfcoptics%Weight(1) = 1.0_fp
    sfcoptics%Transmittance = 0.74_fp

    err_stat = Compute_MW_Water_SfcOptics( &
         surface, geometry, 1, channel_index, sfcoptics, mw_ivar)
    IF (err_stat /= SUCCESS) THEN
      WRITE(*,'("Dispatcher forward returned failure")')
      nfail = nfail + 1
    END IF

    CALL Compute_PARMIO( &
         PARMIOC, frequency, 1, sfcoptics%Angle(1), &
         surface%Water_Temperature, surface%Salinity, surface%Wind_Speed, &
         parmio_ivar, e_ref, r_ref, &
         Azimuth_Angle = surface%Wind_Direction - 10.0_fp, &
         Transmittance = sfcoptics%Transmittance)
    IF (MAXVAL(ABS(sfcoptics%Emissivity(1,:) - e_ref)) > TOL .OR. &
        MAXVAL(ABS((/ sfcoptics%Reflectivity(1,1,1,1), &
                     sfcoptics%Reflectivity(1,2,1,2), &
                     sfcoptics%Reflectivity(1,3,1,3), &
                     sfcoptics%Reflectivity(1,4,1,4) /) - r_ref)) > TOL) THEN
      WRITE(*,'("Dispatcher forward differs from direct PARMIO")')
      nfail = nfail + 1
    END IF

    CALL CRTM_Surface_Zero(surface_tl)
    surface_tl%Water_Temperature = 0.4_fp
    surface_tl%Salinity          = -0.2_fp
    surface_tl%Wind_Speed        = 0.3_fp
    surface_tl%Wind_Direction    = 0.7_fp
    CALL CRTM_SfcOptics_Zero(sfcoptics_tl)
    sfcoptics_tl%Transmittance = -0.01_fp

    err_stat = Compute_MW_Water_SfcOptics_TL( &
         sfcoptics, surface_tl, geometry, 1, channel_index, sfcoptics_tl, mw_ivar)
    IF (err_stat /= SUCCESS) THEN
      WRITE(*,'("Dispatcher TL returned failure")')
      nfail = nfail + 1
    END IF
    CALL Compute_PARMIO_TL( &
         PARMIOC, surface_tl%Water_Temperature, surface_tl%Salinity, &
         surface_tl%Wind_Speed, parmio_ivar, e_tl_ref, r_tl_ref, &
         Azimuth_Angle_TL = surface_tl%Wind_Direction, &
         Transmittance_TL = sfcoptics_tl%Transmittance)
    IF (MAXVAL(ABS(sfcoptics_tl%Emissivity(1,:) - e_tl_ref)) > TOL .OR. &
        MAXVAL(ABS((/ sfcoptics_tl%Reflectivity(1,1,1,1), &
                     sfcoptics_tl%Reflectivity(1,2,1,2), &
                     sfcoptics_tl%Reflectivity(1,3,1,3), &
                     sfcoptics_tl%Reflectivity(1,4,1,4) /) - r_tl_ref)) > TOL) THEN
      WRITE(*,'("Dispatcher TL differs from direct PARMIO")')
      nfail = nfail + 1
    END IF

    CALL CRTM_Surface_Zero(surface_ad)
    CALL CRTM_Surface_Zero(surface_ad_ref)
    CALL CRTM_SfcOptics_Zero(sfcoptics_ad)
    sfcoptics_ad%Emissivity(1,:) = (/ 0.6_fp, -0.5_fp, 0.2_fp, -0.1_fp /)
    sfcoptics_ad%Reflectivity(1,1,1,1) = 0.3_fp
    sfcoptics_ad%Reflectivity(1,2,1,2) = -0.4_fp
    sfcoptics_ad%Reflectivity(1,3,1,3) = 0.1_fp
    sfcoptics_ad%Reflectivity(1,4,1,4) = -0.2_fp
    sfcoptics_ad%Transmittance = 0.05_fp

    e_ad_ref = sfcoptics_ad%Emissivity(1,:)
    r_ad_ref = (/ sfcoptics_ad%Reflectivity(1,1,1,1), &
                  sfcoptics_ad%Reflectivity(1,2,1,2), &
                  sfcoptics_ad%Reflectivity(1,3,1,3), &
                  sfcoptics_ad%Reflectivity(1,4,1,4) /)
    azimuth_ad_ref = 0.0_fp
    trans_ad_ref = sfcoptics_ad%Transmittance
    CALL Compute_PARMIO_AD( &
         PARMIOC, e_ad_ref, r_ad_ref, parmio_ivar, &
         surface_ad_ref%Water_Temperature, surface_ad_ref%Salinity, &
         surface_ad_ref%Wind_Speed, &
         Azimuth_Angle_AD = azimuth_ad_ref, Transmittance_AD = trans_ad_ref)
    surface_ad_ref%Wind_Direction = surface_ad_ref%Wind_Direction + azimuth_ad_ref

    err_stat = Compute_MW_Water_SfcOptics_AD( &
         sfcoptics, sfcoptics_ad, geometry, 1, channel_index, surface_ad, mw_ivar)
    IF (err_stat /= SUCCESS) THEN
      WRITE(*,'("Dispatcher AD returned failure")')
      nfail = nfail + 1
    END IF
    IF (ABS(surface_ad%Water_Temperature - surface_ad_ref%Water_Temperature) > TOL .OR. &
        ABS(surface_ad%Salinity          - surface_ad_ref%Salinity         ) > TOL .OR. &
        ABS(surface_ad%Wind_Speed        - surface_ad_ref%Wind_Speed       ) > TOL .OR. &
        ABS(surface_ad%Wind_Direction    - surface_ad_ref%Wind_Direction   ) > TOL .OR. &
        ABS(sfcoptics_ad%Transmittance   - trans_ad_ref                    ) > TOL) THEN
      WRITE(*,'("Dispatcher AD differs from direct PARMIO")')
      nfail = nfail + 1
    END IF

    CALL CRTM_SfcOptics_Destroy(sfcoptics)
    CALL CRTM_SfcOptics_Destroy(sfcoptics_tl)
    CALL CRTM_SfcOptics_Destroy(sfcoptics_ad)
  END SUBROUTINE Check_Dispatcher

END PROGRAM test_PARMIO_Dispatcher
