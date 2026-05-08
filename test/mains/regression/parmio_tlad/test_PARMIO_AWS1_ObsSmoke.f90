!
! test_PARMIO_AWS1_ObsSmoke
!
! Compare CRTM FASTEM6 and PARMIO brightness-temperature residuals against
! AWS-1 observed TBs for a scene CSV.  The CSV must provide observation
! geometry and background surface state; this program uses the ECMWF84
! atmosphere as a smoke-test placeholder until a full collocation pipeline is
! supplied.
!

PROGRAM test_PARMIO_AWS1_ObsSmoke

  USE CRTM_Module
  USE CRTM_MWwaterCoeff, ONLY: CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff, ONLY: CRTM_PARMIOCoeff_Load, CRTM_PARMIOCoeff_Destroy
  USE CRTM_SpcCoeff, ONLY: SC

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_PARMIO_AWS1_ObsSmoke'
  CHARACTER(*), PARAMETER :: DEFAULT_SENSOR_ID = 'mwr_aws'
  CHARACTER(*), PARAMETER :: DEFAULT_COEFF_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: DEFAULT_LUT_FILE = &
       '../parmio/Outputs/sweep/lut_production/PARMIO.production.MWwater.EmisCoeff.nc'
  CHARACTER(*), PARAMETER :: DEFAULT_SCENE_FILE = 'aws1_scenes.csv'
  CHARACTER(*), PARAMETER :: DEFAULT_RESIDUAL_FILE = 'aws1_residuals.csv'
  CHARACTER(*), PARAMETER :: DEFAULT_SUMMARY_FILE = 'aws1_summary.csv'

  INTEGER, PARAMETER :: N_ATM_PROFILES = 2
  INTEGER, PARAMETER :: N_PROFILES = 1
  INTEGER, PARAMETER :: N_LAYERS = 100
  INTEGER, PARAMETER :: N_ABSORBERS = 6
  INTEGER, PARAMETER :: N_CLOUDS = 0
  INTEGER, PARAMETER :: N_AEROSOLS = 0
  INTEGER, PARAMETER :: N_SENSORS = 1
  INTEGER, PARAMETER :: EXPECTED_AWS_CHANNELS = 19
  REAL(fp), PARAMETER :: WIND_DIRECTION = 0.0_fp

  TYPE :: Scene_type
    INTEGER :: scene_id
    INTEGER :: scan
    INTEGER :: fov
    REAL(fp) :: lat
    REAL(fp) :: lon
    REAL(fp) :: scan_angle
    REAL(fp) :: zenith
    REAL(fp) :: azimuth
    REAL(fp) :: sst
    REAL(fp) :: u10
    REAL(fp) :: sss
    REAL(fp) :: obs_tb(EXPECTED_AWS_CHANNELS)
  END TYPE Scene_type

  CHARACTER(512) :: coeff_path
  CHARACTER(512) :: lut_file
  CHARACTER(512) :: scene_file
  CHARACTER(512) :: residual_file
  CHARACTER(512) :: summary_file
  CHARACTER(512) :: message
  CHARACTER(256) :: version
  CHARACTER(32)  :: sensor_id(N_SENSORS)
  INTEGER :: err_stat
  INTEGER :: allocate_status
  INTEGER :: n_channels
  INTEGER :: scene_unit
  INTEGER :: residual_unit
  INTEGER :: summary_unit
  INTEGER :: processed
  INTEGER :: l
  REAL(fp) :: obs_tb
  REAL(fp) :: tb_fastem
  REAL(fp) :: tb_parmio
  REAL(fp) :: omf_fastem
  REAL(fp) :: omf_parmio
  REAL(fp) :: frequency
  TYPE(Scene_type) :: scene

  TYPE(CRTM_ChannelInfo_type)             :: channel_info(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: atm(N_ATM_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: opt(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rt_fastem(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rt_parmio(:,:)

  INTEGER :: n(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: sum_omf_fastem(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: sum_omf_parmio(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: sum_sq_fastem(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: sum_sq_parmio(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: sum_abs_delta(EXPECTED_AWS_CHANNELS)
  REAL(fp) :: max_abs_delta

  CALL Parse_Arguments(coeff_path, lut_file, scene_file, residual_file, summary_file)

  CALL CRTM_Version(version)
  CALL Program_Message( &
       PROGRAM_NAME, &
       'FASTEM6-vs-PARMIO AWS-1 residual smoke comparison.', &
       'CRTM Version: '//TRIM(version))

  sensor_id = (/ DEFAULT_SENSOR_ID /)
  err_stat = CRTM_Init( &
       sensor_id, channel_info, File_Path=TRIM(coeff_path), &
       SpcCoeff_Format='netCDF', TauCoeff_Format='netCDF')
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Init failed for '//DEFAULT_SENSOR_ID, FAILURE)
    STOP 1
  END IF

  err_stat = CRTM_MWwaterCoeff_Load_FASTEM('FASTEM6', Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Failed to load FASTEM6 MWwater coefficients', FAILURE)
    STOP 1
  END IF

  err_stat = CRTM_PARMIOCoeff_Load(TRIM(lut_file), Quiet=.TRUE.)
  IF (err_stat /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Failed to load PARMIO coefficient LUT', FAILURE)
    STOP 1
  END IF

  n_channels = SUM(CRTM_ChannelInfo_n_Channels(channel_info))
  IF (n_channels /= EXPECTED_AWS_CHANNELS) THEN
    WRITE(message,'("Expected ",i0," AWS channels but CRTM initialized ",i0)') &
         EXPECTED_AWS_CHANNELS, n_channels
    CALL Display_Message(PROGRAM_NAME, TRIM(message), FAILURE)
    STOP 1
  END IF

  ALLOCATE(rt_fastem(n_channels, N_PROFILES), &
           rt_parmio(n_channels, N_PROFILES), STAT=allocate_status)
  IF (allocate_status /= 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'RTSolution allocation failed', FAILURE)
    STOP 1
  END IF
  CALL CRTM_RTSolution_Create(rt_fastem, N_LAYERS)
  CALL CRTM_RTSolution_Create(rt_parmio, N_LAYERS)
  IF (ANY(.NOT. CRTM_RTSolution_Associated(rt_fastem)) .OR. &
      ANY(.NOT. CRTM_RTSolution_Associated(rt_parmio))) THEN
    CALL Display_Message(PROGRAM_NAME, 'RTSolution create failed', FAILURE)
    STOP 1
  END IF

  CALL CRTM_Atmosphere_Create(atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  IF (ANY(.NOT. CRTM_Atmosphere_Associated(atm))) THEN
    CALL Display_Message(PROGRAM_NAME, 'Atmosphere allocation failed', FAILURE)
    STOP 1
  END IF
  CALL Load_ECMWF84_Atm_Data()

  n = 0
  sum_omf_fastem = 0.0_fp
  sum_omf_parmio = 0.0_fp
  sum_sq_fastem = 0.0_fp
  sum_sq_parmio = 0.0_fp
  sum_abs_delta = 0.0_fp
  max_abs_delta = 0.0_fp
  processed = 0

  OPEN(NEWUNIT=scene_unit, FILE=TRIM(scene_file), STATUS='OLD', ACTION='READ', IOSTAT=err_stat)
  IF (err_stat /= 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'Unable to open scene CSV: '//TRIM(scene_file), FAILURE)
    STOP 1
  END IF
  CALL Skip_Header(scene_unit)

  OPEN(NEWUNIT=residual_unit, FILE=TRIM(residual_file), STATUS='REPLACE', ACTION='WRITE', IOSTAT=err_stat)
  IF (err_stat /= 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'Unable to open residual CSV: '//TRIM(residual_file), FAILURE)
    STOP 1
  END IF
  WRITE(residual_unit,'(a)') &
       'scene_id,scan,fov,channel,freq_GHz,obs_tb,tb_fastem,tb_parmio,omf_fastem,omf_parmio,abs_omf_delta'

  DO
    IF (.NOT. Read_Scene(scene_unit, scene)) EXIT
    processed = processed + 1
    CALL Configure_Scene(sfc, geometry, scene)

    opt(:)%Use_PARMIO_Model = .FALSE.
    err_stat = CRTM_Forward( &
         atm(1:N_PROFILES), sfc, geometry, channel_info, rt_fastem, Options=opt)
    IF (err_stat /= SUCCESS) THEN
      CALL Display_Message(PROGRAM_NAME, 'FASTEM CRTM_Forward failed', FAILURE)
      STOP 1
    END IF

    opt(:)%Use_PARMIO_Model = .TRUE.
    err_stat = CRTM_Forward( &
         atm(1:N_PROFILES), sfc, geometry, channel_info, rt_parmio, Options=opt)
    IF (err_stat /= SUCCESS) THEN
      CALL Display_Message(PROGRAM_NAME, 'PARMIO CRTM_Forward failed', FAILURE)
      STOP 1
    END IF

    DO l = 1, n_channels
      obs_tb = scene%obs_tb(l)
      tb_fastem = rt_fastem(l,1)%Brightness_Temperature
      tb_parmio = rt_parmio(l,1)%Brightness_Temperature
      omf_fastem = obs_tb - tb_fastem
      omf_parmio = obs_tb - tb_parmio
      frequency = SC(channel_info(1)%Sensor_Index)%Frequency(channel_info(1)%Channel_Index(l))

      n(l) = n(l) + 1
      sum_omf_fastem(l) = sum_omf_fastem(l) + omf_fastem
      sum_omf_parmio(l) = sum_omf_parmio(l) + omf_parmio
      sum_sq_fastem(l) = sum_sq_fastem(l) + omf_fastem**2
      sum_sq_parmio(l) = sum_sq_parmio(l) + omf_parmio**2
      sum_abs_delta(l) = sum_abs_delta(l) + (ABS(omf_parmio) - ABS(omf_fastem))
      max_abs_delta = MAX(max_abs_delta, ABS(tb_parmio - tb_fastem))

      WRITE(residual_unit,'(i0,",",i0,",",i0,",",i0,",",f10.4,",", &
                           f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
           scene%scene_id, scene%scan, scene%fov, rt_fastem(l,1)%Sensor_Channel, frequency, &
           obs_tb, tb_fastem, tb_parmio, omf_fastem, omf_parmio, &
           ABS(omf_parmio) - ABS(omf_fastem)
    END DO
  END DO

  CLOSE(scene_unit)
  CLOSE(residual_unit)

  IF (processed == 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'No scenes were read from '//TRIM(scene_file), FAILURE)
    STOP 1
  END IF

  OPEN(NEWUNIT=summary_unit, FILE=TRIM(summary_file), STATUS='REPLACE', ACTION='WRITE', IOSTAT=err_stat)
  IF (err_stat /= 0) THEN
    CALL Display_Message(PROGRAM_NAME, 'Unable to open summary CSV: '//TRIM(summary_file), FAILURE)
    STOP 1
  END IF
  WRITE(summary_unit,'(a)') &
       'channel,freq_GHz,n,bias_fastem,bias_parmio,rmse_fastem,rmse_parmio,rmse_delta,mean_abs_omf_delta'
  DO l = 1, n_channels
    frequency = SC(channel_info(1)%Sensor_Index)%Frequency(channel_info(1)%Channel_Index(l))
    WRITE(summary_unit,'(i0,",",f10.4,",",i0,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6,",",f12.6)') &
         rt_fastem(l,1)%Sensor_Channel, frequency, n(l), &
         sum_omf_fastem(l)/REAL(n(l),fp), &
         sum_omf_parmio(l)/REAL(n(l),fp), &
         SQRT(sum_sq_fastem(l)/REAL(n(l),fp)), &
         SQRT(sum_sq_parmio(l)/REAL(n(l),fp)), &
         SQRT(sum_sq_parmio(l)/REAL(n(l),fp)) - SQRT(sum_sq_fastem(l)/REAL(n(l),fp)), &
         sum_abs_delta(l)/REAL(n(l),fp)
  END DO
  CLOSE(summary_unit)

  WRITE(*,'("PARMIO AWS-1 obs smoke comparison complete: scenes=",i0, &
            ", channels=",i0,", max_abs_model_delta=",f10.4," K")') &
       processed, n_channels, max_abs_delta
  WRITE(*,'("Residual CSV: ",a)') TRIM(residual_file)
  WRITE(*,'("Summary CSV: ",a)') TRIM(summary_file)

  CALL CRTM_PARMIOCoeff_Destroy()
  err_stat = CRTM_Destroy(channel_info)
  CALL CRTM_Atmosphere_Destroy(atm)
  DEALLOCATE(rt_fastem, rt_parmio)

CONTAINS

  SUBROUTINE Parse_Arguments(coeff_path, lut_file, scene_file, residual_file, summary_file)
    CHARACTER(*), INTENT(OUT) :: coeff_path
    CHARACTER(*), INTENT(OUT) :: lut_file
    CHARACTER(*), INTENT(OUT) :: scene_file
    CHARACTER(*), INTENT(OUT) :: residual_file
    CHARACTER(*), INTENT(OUT) :: summary_file

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
      CALL GET_COMMAND_ARGUMENT(3, scene_file)
    ELSE
      scene_file = DEFAULT_SCENE_FILE
    END IF
    IF (COMMAND_ARGUMENT_COUNT() >= 4) THEN
      CALL GET_COMMAND_ARGUMENT(4, residual_file)
    ELSE
      residual_file = DEFAULT_RESIDUAL_FILE
    END IF
    IF (COMMAND_ARGUMENT_COUNT() >= 5) THEN
      CALL GET_COMMAND_ARGUMENT(5, summary_file)
    ELSE
      summary_file = DEFAULT_SUMMARY_FILE
    END IF
  END SUBROUTINE Parse_Arguments

  SUBROUTINE Skip_Header(unit)
    INTEGER, INTENT(IN) :: unit
    CHARACTER(4096) :: line
    READ(unit,'(a)',IOSTAT=err_stat) line
    IF (err_stat /= 0) THEN
      CALL Display_Message(PROGRAM_NAME, 'Scene CSV is empty', FAILURE)
      STOP 1
    END IF
  END SUBROUTINE Skip_Header

  LOGICAL FUNCTION Read_Scene(unit, scene)
    INTEGER, INTENT(IN) :: unit
    TYPE(Scene_type), INTENT(OUT) :: scene
    CHARACTER(4096) :: line

    DO
      READ(unit,'(a)',IOSTAT=err_stat) line
      IF (err_stat /= 0) THEN
        Read_Scene = .FALSE.
        RETURN
      END IF
      IF (LEN_TRIM(line) /= 0) EXIT
    END DO

    READ(line,*,IOSTAT=err_stat) &
         scene%scene_id, scene%scan, scene%fov, scene%lat, scene%lon, &
         scene%scan_angle, scene%zenith, scene%azimuth, scene%sst, scene%u10, &
         scene%sss, scene%obs_tb
    IF (err_stat /= 0) THEN
      CALL Display_Message(PROGRAM_NAME, 'Malformed scene CSV row: '//TRIM(line), FAILURE)
      STOP 1
    END IF
    Read_Scene = .TRUE.
  END FUNCTION Read_Scene

  SUBROUTINE Configure_Scene(sfc, geometry, scene)
    TYPE(CRTM_Surface_type),  INTENT(IN OUT) :: sfc(:)
    TYPE(CRTM_Geometry_type), INTENT(IN OUT) :: geometry(:)
    TYPE(Scene_type),         INTENT(IN)     :: scene

    CALL CRTM_Surface_Zero(sfc)
    sfc(1)%Water_Coverage    = 1.0_fp
    sfc(1)%Water_Temperature = scene%sst
    sfc(1)%Wind_Speed        = scene%u10
    sfc(1)%Wind_Direction    = WIND_DIRECTION
    sfc(1)%Salinity          = scene%sss

    CALL CRTM_Geometry_SetValue( &
         geometry, &
         iFOV                 = scene%fov + 1, &
         Longitude            = scene%lon, &
         Latitude             = scene%lat, &
         Sensor_Scan_Angle    = scene%scan_angle, &
         Sensor_Zenith_Angle  = scene%zenith, &
         Sensor_Azimuth_Angle = scene%azimuth, &
         Source_Zenith_Angle  = 100.0_fp, &
         Source_Azimuth_Angle = 0.0_fp, &
         Year                 = 2026, &
         Month                = 4, &
         Day                  = 9)
  END SUBROUTINE Configure_Scene

  INCLUDE '../../unit/Unit_Test/Load_ECMWF84_Atm_Data.inc'

END PROGRAM test_PARMIO_AWS1_ObsSmoke
