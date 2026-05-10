!
! test_PARMIO_Lifecycle
!
! Exercise the Phase-4 lifecycle integration:
!   - CRTM_Init(..., PARMIOCoeff_File=...) loads the LUT.
!   - Options%Use_PARMIO_Model toggles dispatch in MW_Water surface optics.
!   - With LUT loaded + flag on, CRTM_Forward routes through PARMIO.
!   - With LUT loaded + flag off, behaviour is byte-identical to FASTEM.
!   - With NO LUT loaded + flag on, dispatcher silently falls back to FASTEM.
!   - CRTM_Destroy frees the LUT (no leak; re-init works).
!
! Single ATMS-NPP scene; the broader grid sweep lives in
! test_PARMIO_FASTEM_DeltaSweep.
!

PROGRAM test_PARMIO_Lifecycle

  USE CRTM_Module
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_IsLoaded

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME       = 'test_PARMIO_Lifecycle'
  CHARACTER(*), PARAMETER :: DEFAULT_COEFF_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: DEFAULT_LUT_FILE   = 'PARMIO.MWwater.EmisCoeff.nc'
  CHARACTER(*), PARAMETER :: SENSOR_ID          = 'atms_npp'

  INTEGER, PARAMETER  :: N_ATM_PROFILES = 2  ! Load_ECMWF84_Atm_Data.inc fills two
  INTEGER, PARAMETER  :: N_PROFILES   = 1
  INTEGER, PARAMETER  :: N_LAYERS     = 100
  INTEGER, PARAMETER  :: N_ABSORBERS  = 6
  INTEGER, PARAMETER  :: N_CLOUDS     = 0
  INTEGER, PARAMETER  :: N_AEROSOLS   = 0
  INTEGER, PARAMETER  :: N_SENSORS    = 1
  REAL(fp), PARAMETER :: SST          = 288.15_fp
  REAL(fp), PARAMETER :: U10          = 10.0_fp
  REAL(fp), PARAMETER :: ZENITH       = 45.0_fp
  REAL(fp), PARAMETER :: SALINITY_PSU = 35.0_fp
  REAL(fp), PARAMETER :: MIN_DELTA_TB = 0.05_fp  ! V or H ATMS channel must differ at least this much

  CHARACTER(512) :: coeff_path
  CHARACTER(512) :: lut_file
  INTEGER :: err_stat
  INTEGER :: allocate_status
  INTEGER :: n_channels
  INTEGER :: l
  INTEGER :: nfail
  REAL(fp) :: max_abs_delta_loaded
  REAL(fp) :: max_abs_delta_fallback

  TYPE(CRTM_ChannelInfo_type)             :: channel_info(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: atm(N_ATM_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: opt(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rt_a(:,:), rt_b(:,:), rt_c(:,:)

  CALL Parse_Arguments(coeff_path, lut_file)
  nfail = 0

  ! ---- Phase A: CRTM_Init WITH PARMIOCoeff_File; LUT should be loaded.
  err_stat = CRTM_Init( (/ SENSOR_ID /), channel_info, &
                       File_Path        = TRIM(coeff_path), &
                       PARMIOCoeff_File = TRIM(lut_file), &
                       Quiet            = .TRUE. )
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Init with PARMIOCoeff_File failed', FAILURE)
    STOP 1
  END IF
  IF (.NOT. CRTM_PARMIOCoeff_IsLoaded()) THEN
    CALL Display_Message(PROGRAM_NAME, &
         'PARMIOCoeff_File supplied but CRTM_PARMIOCoeff_IsLoaded() is .FALSE.', FAILURE)
    STOP 1
  END IF

  n_channels = SUM(CRTM_ChannelInfo_n_Channels(channel_info))
  ALLOCATE(rt_a(n_channels, N_PROFILES), &
           rt_b(n_channels, N_PROFILES), &
           rt_c(n_channels, N_PROFILES), STAT=allocate_status)
  IF (allocate_status /= 0) STOP 'RTSolution allocation failed'
  CALL CRTM_RTSolution_Create(rt_a, N_LAYERS)
  CALL CRTM_RTSolution_Create(rt_b, N_LAYERS)
  CALL CRTM_RTSolution_Create(rt_c, N_LAYERS)

  CALL CRTM_Atmosphere_Create(atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  CALL Load_ECMWF84_Atm_Data()
  CALL Configure_Surface_And_Geometry(sfc, geometry)

  opt(:)%Use_PARMIO_Model = .TRUE.
  err_stat = CRTM_Forward(atm(1:N_PROFILES), sfc, geometry, channel_info, rt_a, Options=opt)
  IF (err_stat /= SUCCESS) STOP 'CRTM_Forward with PARMIO failed'

  opt(:)%Use_PARMIO_Model = .FALSE.
  err_stat = CRTM_Forward(atm(1:N_PROFILES), sfc, geometry, channel_info, rt_b, Options=opt)
  IF (err_stat /= SUCCESS) STOP 'CRTM_Forward with FASTEM (LUT loaded) failed'

  max_abs_delta_loaded = 0.0_fp
  DO l = 1, n_channels
    max_abs_delta_loaded = MAX(max_abs_delta_loaded, &
      ABS(rt_a(l,1)%Brightness_Temperature - rt_b(l,1)%Brightness_Temperature))
  END DO
  IF (max_abs_delta_loaded < MIN_DELTA_TB) THEN
    WRITE(*,'("FAIL: dispatch did not route through PARMIO. max |dT| = ",f10.6, &
            " K; expected >= ",f10.6)') max_abs_delta_loaded, MIN_DELTA_TB
    nfail = nfail + 1
  END IF

  err_stat = CRTM_Destroy(channel_info)
  IF (err_stat /= SUCCESS) STOP 'CRTM_Destroy (Phase A) failed'
  IF (CRTM_PARMIOCoeff_IsLoaded()) THEN
    CALL Display_Message(PROGRAM_NAME, &
         'PARMIO LUT still loaded after CRTM_Destroy', FAILURE)
    nfail = nfail + 1
  END IF

  ! ---- Phase B: CRTM_Init WITHOUT PARMIOCoeff_File; flag-on should fall back.
  err_stat = CRTM_Init( (/ SENSOR_ID /), channel_info, &
                       File_Path = TRIM(coeff_path), &
                       Quiet     = .TRUE. )
  IF (err_stat /= SUCCESS) STOP 'CRTM_Init without PARMIOCoeff_File failed'
  IF (CRTM_PARMIOCoeff_IsLoaded()) THEN
    CALL Display_Message(PROGRAM_NAME, &
         'PARMIO LUT loaded despite PARMIOCoeff_File not being passed', FAILURE)
    STOP 1
  END IF

  opt(:)%Use_PARMIO_Model = .TRUE.
  err_stat = CRTM_Forward(atm(1:N_PROFILES), sfc, geometry, channel_info, rt_c, Options=opt)
  IF (err_stat /= SUCCESS) STOP 'CRTM_Forward fallback path failed'

  max_abs_delta_fallback = 0.0_fp
  DO l = 1, n_channels
    max_abs_delta_fallback = MAX(max_abs_delta_fallback, &
      ABS(rt_c(l,1)%Brightness_Temperature - rt_b(l,1)%Brightness_Temperature))
  END DO
  IF (max_abs_delta_fallback /= 0.0_fp) THEN
    WRITE(*,'("FAIL: fallback path not byte-identical to FASTEM. max |dT| = ",es12.4)') &
      max_abs_delta_fallback
    nfail = nfail + 1
  END IF

  err_stat = CRTM_Destroy(channel_info)
  IF (err_stat /= SUCCESS) STOP 'CRTM_Destroy (Phase B) failed'

  CALL CRTM_Atmosphere_Destroy(atm)
  DEALLOCATE(rt_a, rt_b, rt_c)

  IF (nfail > 0) THEN
    WRITE(*,'("test_PARMIO_Lifecycle: ",i0," failure(s)")') nfail
    ERROR STOP 'PARMIO lifecycle integration test failed'
  END IF

  WRITE(*,'("test_PARMIO_Lifecycle passed: max PARMIO/FASTEM delta = ",f8.4, &
          " K; fallback delta = ",es10.2)') max_abs_delta_loaded, max_abs_delta_fallback

CONTAINS

  SUBROUTINE Parse_Arguments(coeff_path, lut_file)
    CHARACTER(*), INTENT(OUT) :: coeff_path
    CHARACTER(*), INTENT(OUT) :: lut_file
    IF (COMMAND_ARGUMENT_COUNT() >= 1) THEN
      CALL GET_COMMAND_ARGUMENT(1, coeff_path)
    ELSE
      coeff_path = DEFAULT_COEFF_PATH
    END IF
    IF (COMMAND_ARGUMENT_COUNT() >= 2) THEN
      CALL GET_COMMAND_ARGUMENT(2, lut_file)
    ELSE
      lut_file = DEFAULT_LUT_FILE
    END IF
  END SUBROUTINE Parse_Arguments

  SUBROUTINE Configure_Surface_And_Geometry(sfc, geometry)
    TYPE(CRTM_Surface_type),  INTENT(IN OUT) :: sfc(:)
    TYPE(CRTM_Geometry_type), INTENT(IN OUT) :: geometry(:)
    CALL CRTM_Surface_Zero(sfc)
    sfc(1)%Water_Coverage    = 1.0_fp
    sfc(1)%Water_Temperature = SST
    sfc(1)%Wind_Speed        = U10
    sfc(1)%Wind_Direction    = 0.0_fp
    sfc(1)%Salinity          = SALINITY_PSU
    CALL CRTM_Geometry_SetValue( &
         geometry, &
         Sensor_Zenith_Angle  = ZENITH, &
         Sensor_Scan_Angle    = ZENITH, &
         Sensor_Azimuth_Angle = 0.0_fp, &
         Source_Zenith_Angle  = 100.0_fp, &
         Source_Azimuth_Angle = 0.0_fp)
  END SUBROUTINE Configure_Surface_And_Geometry

  INCLUDE '../../unit/Unit_Test/Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_PARMIO_Lifecycle
