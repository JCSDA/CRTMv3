!
! test_PARMIO_FASTEM_DeltaSweep
!
! Compare FASTEM6 and PARMIO top-of-atmosphere microwave brightness
! temperatures through the full CRTM forward path for a small ATMS-NPP
! ocean-state grid.
!

PROGRAM test_PARMIO_FASTEM_DeltaSweep

  USE CRTM_Module
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy
  USE CRTM_SpcCoeff, ONLY: SC

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_PARMIO_FASTEM_DeltaSweep'
  CHARACTER(*), PARAMETER :: DEFAULT_COEFF_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: DEFAULT_LUT_FILE = &
       './testinput/PARMIO.MWwater.EmisCoeff.nc'
  CHARACTER(*), PARAMETER :: DEFAULT_SENSOR_ID = 'atms_npp'

  INTEGER, PARAMETER :: N_ATM_PROFILES = 2
  INTEGER, PARAMETER :: N_PROFILES = 1
  INTEGER, PARAMETER :: N_LAYERS = 100
  INTEGER, PARAMETER :: N_ABSORBERS = 6
  INTEGER, PARAMETER :: N_CLOUDS = 0
  INTEGER, PARAMETER :: N_AEROSOLS = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  INTEGER, PARAMETER :: N_SST = 4
  INTEGER, PARAMETER :: N_U10 = 5
  INTEGER, PARAMETER :: N_ZENITH = 2

  REAL(fp), PARAMETER :: SST_GRID(N_SST) = (/ 275.0_fp, 285.0_fp, 295.0_fp, 300.0_fp /)
  REAL(fp), PARAMETER :: U10_GRID(N_U10) = (/ 2.0_fp, 5.0_fp, 10.0_fp, 15.0_fp, 20.0_fp /)
  REAL(fp), PARAMETER :: ZENITH_GRID(N_ZENITH) = (/ 30.0_fp, 55.0_fp /)
  REAL(fp), PARAMETER :: SALINITY = 35.0_fp
  REAL(fp), PARAMETER :: WIND_DIRECTION = 0.0_fp
  REAL(fp), PARAMETER :: MIN_VALID_TB = 2.7_fp
  ! Loosened above ATMS's 30 K to accommodate AWS-class 325 GHz channels in
  ! cold/windy regimes; PARMIO and FASTEM physics genuinely diverge more there.
  REAL(fp), PARAMETER :: MAX_DELTA_TB = 50.0_fp

  CHARACTER(512) :: coeff_path
  CHARACTER(512) :: lut_file
  CHARACTER(256) :: message
  CHARACTER(256) :: version
  CHARACTER(32)  :: sensor_id_arg
  CHARACTER(64)  :: csv_file
  INTEGER :: err_stat
  INTEGER :: allocate_status
  INTEGER :: n_channels
  INTEGER :: csv_unit
  INTEGER :: case_id
  INTEGER :: i_sst, i_u10, i_zenith, l
  INTEGER :: n_cases
  REAL(fp) :: tb_fastem
  REAL(fp) :: tb_parmio
  REAL(fp) :: delta_tb
  REAL(fp) :: frequency

  TYPE(CRTM_ChannelInfo_type)             :: channel_info(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: atm(N_ATM_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: opt(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rt(:,:)
  ! Per-(channel, case) brightness-temperature buffers for the two phases.
  ! FASTEM-vs-PARMIO is now determined inside CRTM by frequency (>=200 GHz
  ! routes to PARMIO when the LUT is loaded), so we run two passes over the
  ! same grid: phase 1 with no LUT loaded (pure FASTEM), phase 2 with the
  ! LUT loaded (PARMIO at high frequency channels, FASTEM elsewhere).
  REAL(fp), ALLOCATABLE :: tb_fastem_grid(:,:)
  REAL(fp), ALLOCATABLE :: tb_parmio_grid(:,:)

  CALL Parse_Arguments(coeff_path, lut_file, sensor_id_arg)
  csv_file = 'delta_sweep_'//TRIM(sensor_id_arg)//'.csv'

  CALL CRTM_Version(version)
  CALL Program_Message( &
       PROGRAM_NAME, &
       'FASTEM6-vs-PARMIO '//TRIM(sensor_id_arg)//' brightness-temperature delta sweep.', &
       'CRTM Version: '//TRIM(version))

  err_stat = CRTM_Init((/ sensor_id_arg /), channel_info, File_Path=TRIM(coeff_path))
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Init failed', FAILURE)
    STOP 1
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Failed to load FASTEM6 MWwater coefficients', FAILURE)
    STOP 1
  END IF

  n_channels = SUM(CRTM_ChannelInfo_n_Channels(channel_info))
  n_cases    = N_SST * N_U10 * N_ZENITH
  ALLOCATE(rt(n_channels, N_PROFILES), &
           tb_fastem_grid(n_channels, n_cases), &
           tb_parmio_grid(n_channels, n_cases), STAT=allocate_status)
  IF (allocate_status /= 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'RTSolution / TB-grid allocation failed', FAILURE)
    STOP 1
  END IF
  CALL CRTM_RTSolution_Create(rt, N_LAYERS)
  IF (ANY(.NOT. CRTM_RTSolution_Associated(rt))) THEN
    CALL Display_Message(PROGRAM_NAME, 'RTSolution create failed', FAILURE)
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Create(atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  IF (ANY(.NOT. CRTM_Atmosphere_Associated(atm))) THEN
    CALL Display_Message(PROGRAM_NAME, 'Atmosphere allocation failed', FAILURE)
    STOP 1
  END IF
  CALL Load_ECMWF84_Atm_Data()

  ! ---- Phase 1: FASTEM-only sweep (PARMIO LUT not loaded) ----
  CALL Run_Grid_Sweep(tb_fastem_grid)

  ! ---- Phase 2: PARMIO sweep (LUT loaded, dispatcher routes >=200 GHz channels) ----
  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Failed to load PARMIO coefficient LUT', FAILURE)
    STOP 1
  END IF
  CALL Run_Grid_Sweep(tb_parmio_grid)
  CALL CRTM_PARMIOCoeff_Destroy()

  ! ---- Compare phases, write CSV, gate on sanity bounds ----
  OPEN(NEWUNIT=csv_unit, FILE=TRIM(csv_file), STATUS='REPLACE', ACTION='WRITE')
  WRITE(csv_unit,'(a)') &
       'case_id,channel,freq_GHz,sst,u10,zenith,tb_fastem,tb_parmio,delta_tb'

  case_id = 0
  DO i_zenith = 1, N_ZENITH
    DO i_sst = 1, N_SST
      DO i_u10 = 1, N_U10
        case_id = case_id + 1
        DO l = 1, n_channels
          tb_fastem = tb_fastem_grid(l, case_id)
          tb_parmio = tb_parmio_grid(l, case_id)
          delta_tb  = tb_parmio - tb_fastem
          frequency = SC(channel_info(1)%Sensor_Index)%Frequency(channel_info(1)%Channel_Index(l))

          CALL Check_Row( &
               case_id, channel_info(1)%Sensor_Channel(l), SST_GRID(i_sst), U10_GRID(i_u10), &
               ZENITH_GRID(i_zenith), tb_fastem, tb_parmio, delta_tb)

          WRITE(csv_unit,'(i0,",",i0,",",f10.4,",",f8.2,",",f8.2,",",f8.2,",", &
                          f12.6,",",f12.6,",",f12.6)') &
               case_id, channel_info(1)%Sensor_Channel(l), frequency, SST_GRID(i_sst), &
               U10_GRID(i_u10), ZENITH_GRID(i_zenith), tb_fastem, tb_parmio, delta_tb
        END DO
      END DO
    END DO
  END DO

  CLOSE(csv_unit)

  WRITE(*,'("PARMIO FASTEM delta sweep passed: sensor=",a,", ",i0," cases, ",i0, &
            " channels, CSV=",a)') TRIM(sensor_id_arg), case_id, n_channels, TRIM(csv_file)

  err_stat = CRTM_Destroy(channel_info)
  CALL CRTM_Atmosphere_Destroy(atm)
  DEALLOCATE(rt, tb_fastem_grid, tb_parmio_grid)

CONTAINS

  ! Walk the SST x U10 x ZENITH grid; store brightness temperatures into
  ! tb_out(:,:) keyed by (channel, case_id). Whether the call resolves to
  ! FASTEM or PARMIO is determined solely by whether the PARMIO LUT was
  ! loaded prior to invocation (and per-channel frequency >= 200 GHz).
  SUBROUTINE Run_Grid_Sweep(tb_out)
    REAL(fp), INTENT(OUT) :: tb_out(:,:)
    INTEGER :: i_sst_l, i_u10_l, i_zenith_l, ll, case_local

    case_local = 0
    DO i_zenith_l = 1, N_ZENITH
      DO i_sst_l = 1, N_SST
        DO i_u10_l = 1, N_U10
          case_local = case_local + 1
          CALL Configure_Case( &
               sfc, geometry, &
               SST_GRID(i_sst_l), U10_GRID(i_u10_l), ZENITH_GRID(i_zenith_l))
          err_stat = CRTM_Forward( &
               atm(1:N_PROFILES), sfc, geometry, channel_info, rt, Options=opt)
          IF (err_stat /= SUCCESS) THEN
            CALL Display_Message(PROGRAM_NAME, 'CRTM_Forward failed in Run_Grid_Sweep', FAILURE)
            STOP 1
          END IF
          DO ll = 1, SIZE(tb_out, 1)
            tb_out(ll, case_local) = rt(ll, 1)%Brightness_Temperature
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE Run_Grid_Sweep

  SUBROUTINE Parse_Arguments(coeff_path, lut_file, sensor_id)
    CHARACTER(*), INTENT(OUT) :: coeff_path
    CHARACTER(*), INTENT(OUT) :: lut_file
    CHARACTER(*), INTENT(OUT) :: sensor_id

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

    IF (COMMAND_ARGUMENT_COUNT() >= 3) THEN
      CALL GET_COMMAND_ARGUMENT(3, sensor_id)
    ELSE
      sensor_id = DEFAULT_SENSOR_ID
    END IF
  END SUBROUTINE Parse_Arguments

  SUBROUTINE Configure_Case(sfc, geometry, sst, u10, zenith)
    TYPE(CRTM_Surface_type),  INTENT(IN OUT) :: sfc(:)
    TYPE(CRTM_Geometry_type), INTENT(IN OUT) :: geometry(:)
    REAL(fp),                 INTENT(IN)     :: sst
    REAL(fp),                 INTENT(IN)     :: u10
    REAL(fp),                 INTENT(IN)     :: zenith

    CALL CRTM_Surface_Zero(sfc)
    sfc(1)%Water_Coverage    = 1.0_fp
    sfc(1)%Water_Temperature = sst
    sfc(1)%Wind_Speed        = u10
    sfc(1)%Wind_Direction    = WIND_DIRECTION
    sfc(1)%Salinity          = SALINITY

    CALL CRTM_Geometry_SetValue( &
         geometry, &
         Sensor_Zenith_Angle  = zenith, &
         Sensor_Scan_Angle    = zenith, &
         Sensor_Azimuth_Angle = 0.0_fp, &
         Source_Zenith_Angle  = 100.0_fp, &
         Source_Azimuth_Angle = 0.0_fp)
  END SUBROUTINE Configure_Case

  SUBROUTINE Check_Row(case_id, channel, sst, u10, zenith, tb_fastem, tb_parmio, delta_tb)
    INTEGER,  INTENT(IN) :: case_id
    INTEGER,  INTENT(IN) :: channel
    REAL(fp), INTENT(IN) :: sst
    REAL(fp), INTENT(IN) :: u10
    REAL(fp), INTENT(IN) :: zenith
    REAL(fp), INTENT(IN) :: tb_fastem
    REAL(fp), INTENT(IN) :: tb_parmio
    REAL(fp), INTENT(IN) :: delta_tb

    IF (tb_fastem <= MIN_VALID_TB .OR. tb_parmio <= MIN_VALID_TB .OR. &
        tb_fastem >= sst + 5.0_fp .OR. tb_parmio >= sst + 5.0_fp .OR. &
        ABS(delta_tb) >= MAX_DELTA_TB) THEN
      WRITE(message,'("Sanity gate failed: case=",i0,", channel=",i0, &
                     ", sst=",f7.2,", u10=",f6.2,", zenith=",f6.2, &
                     ", tb_fastem=",f10.4,", tb_parmio=",f10.4, &
                     ", delta_tb=",f10.4)') &
           case_id, channel, sst, u10, zenith, tb_fastem, tb_parmio, delta_tb
      CALL Display_Message(PROGRAM_NAME, TRIM(message), FAILURE)
      ERROR STOP 'PARMIO FASTEM delta sweep sanity gate failed'
    END IF
  END SUBROUTINE Check_Row

  INCLUDE '../../unit/Unit_Test/Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_PARMIO_FASTEM_DeltaSweep
