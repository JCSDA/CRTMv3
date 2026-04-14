!
! test_taucoeff_emulator
!
! Compare an LBL-emulator-generated TauCoeff netCDF file against the
! reference binary TauCoeff for the same sensor:
!   (1) Structural check via ODPS netCDF/binary inquire (dimensions, IDs,
!       absorbers, reference profile).
!   (2) Forward-model brightness-temperature comparison over a standard
!       atmosphere set; assert max |dBT| and RMSE under tolerances.
!
! Both loads use the same testinput/ directory:
!   - Reference: <sensor>.TauCoeff.bin (Binary, default format)
!   - Emulator:  <sensor>.TauCoeff.nc  (netCDF, via TauCoeff_Format arg)
!
! Invocation: test_TauCoeff_Emulator <sensor_id>
!

PROGRAM test_taucoeff_emulator

  USE CRTM_Module
  USE ODPS_Define,     ONLY: ODPS_type, Destroy_ODPS
  USE ODPS_netCDF_IO,  ONLY: Read_ODPS_netCDF
  USE ODPS_Binary_IO,  ONLY: Read_ODPS_Binary
  USE UnitTest_Define, ONLY: UnitTest_type, UnitTest_Init, UnitTest_Setup, &
                             UnitTest_Assert, UnitTest_Passed, UnitTest_Report

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_taucoeff_emulator'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'

  ! Profile setup mirrors test_Simple (UMBC 92-layer, 6-gas, 2 profiles)
  INTEGER,  PARAMETER :: N_PROFILES  = 2
  INTEGER,  PARAMETER :: N_LAYERS    = 92
  INTEGER,  PARAMETER :: N_ABSORBERS = 2
  INTEGER,  PARAMETER :: N_CLOUDS    = 1
  INTEGER,  PARAMETER :: N_AEROSOLS  = 1
  INTEGER,  PARAMETER :: N_SENSORS   = 1
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! Tolerances (K). Tune after baseline run.
  REAL(fp), PARAMETER :: BT_MAX_TOL  = 1.0_fp
  REAL(fp), PARAMETER :: BT_RMSE_TOL = 0.3_fp

  CHARACTER(256) :: Message, Sensor_Id, ref_file, emu_file
  INTEGER :: Error_Status, n_Channels, l, m, ac, n_valid
  REAL(fp), ALLOCATABLE :: BT_ref(:,:), BT_emu(:,:), dBT(:,:)
  REAL(fp) :: max_abs, rmse, sse

  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  TYPE(ODPS_type) :: ODPS_ref, ODPS_emu
  TYPE(UnitTest_type) :: ut
  LOGICAL :: ok

  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
    WRITE(*,*) PROGRAM_NAME//': ERROR — expected one argument: <sensor_id>'
    STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1, Sensor_Id)
  Sensor_Id = ADJUSTL(Sensor_Id)

  CALL UnitTest_Init(ut, .TRUE.)
  CALL UnitTest_Setup(ut, 'TauCoeff_Emulator_'//TRIM(Sensor_Id), PROGRAM_NAME, .TRUE.)

  WRITE(*,'(/5x,"Comparing emulator vs. reference TauCoeff for ",a)') TRIM(Sensor_Id)

  ! =========================================================================
  ! 1. STRUCTURAL CHECK
  ! =========================================================================
  ref_file = COEFFICIENTS_PATH//TRIM(Sensor_Id)//'.TauCoeff.bin'
  emu_file = COEFFICIENTS_PATH//TRIM(Sensor_Id)//'.TauCoeff.nc'

  Error_Status = Read_ODPS_Binary( TRIM(ref_file), ODPS_ref, Quiet=1 )
  CALL UnitTest_Assert(ut, (Error_Status == SUCCESS))
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Error reading reference ODPS: '//TRIM(ref_file), FAILURE)
    STOP 1
  END IF

  Error_Status = Read_ODPS_netCDF( TRIM(emu_file), ODPS_emu, Quiet=1 )
  CALL UnitTest_Assert(ut, (Error_Status == SUCCESS))
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Error reading emulator ODPS: '//TRIM(emu_file), FAILURE)
    STOP 1
  END IF

  WRITE(*,'(5x,"  ref: n_Layers=",i0,"  n_Channels=",i0,"  n_Absorbers=",i0)') &
        ODPS_ref%n_Layers, ODPS_ref%n_Channels, ODPS_ref%n_Absorbers
  WRITE(*,'(5x,"  emu: n_Layers=",i0,"  n_Channels=",i0,"  n_Absorbers=",i0)') &
        ODPS_emu%n_Layers, ODPS_emu%n_Channels, ODPS_emu%n_Absorbers

  CALL UnitTest_Assert(ut, ODPS_emu%n_Layers    == ODPS_ref%n_Layers)
  CALL UnitTest_Assert(ut, ODPS_emu%n_Channels  == ODPS_ref%n_Channels)
  CALL UnitTest_Assert(ut, ODPS_emu%n_Absorbers == ODPS_ref%n_Absorbers)
  CALL UnitTest_Assert(ut, ODPS_emu%Sensor_Type == ODPS_ref%Sensor_Type)
  CALL UnitTest_Assert(ut, TRIM(ODPS_emu%Sensor_Id) == TRIM(ODPS_ref%Sensor_Id))
  CALL UnitTest_Assert(ut, ODPS_emu%WMO_Sensor_Id     == ODPS_ref%WMO_Sensor_Id)
  CALL UnitTest_Assert(ut, ODPS_emu%WMO_Satellite_Id  == ODPS_ref%WMO_Satellite_Id)

  IF (ODPS_emu%n_Absorbers == ODPS_ref%n_Absorbers) THEN
    DO ac = 1, ODPS_ref%n_Absorbers
      CALL UnitTest_Assert(ut, ODPS_emu%Absorber_ID(ac) == ODPS_ref%Absorber_ID(ac))
    END DO
  END IF

  Error_Status = Destroy_ODPS(ODPS_ref)
  Error_Status = Destroy_ODPS(ODPS_emu)

  ! =========================================================================
  ! 2. FORWARD-MODEL COMPARISON
  ! =========================================================================
  !
  ! Pass A: reference (binary TauCoeff)
  !
  WRITE(*,'(/5x,"Initializing CRTM with REFERENCE TauCoeff...")')
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Error initializing CRTM (reference)', FAILURE)
    STOP 1
  END IF
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  WRITE(*,'(5x,"  n_Channels = ",i0)') n_Channels

  ALLOCATE(RTSolution(n_Channels, N_PROFILES))
  CALL CRTM_RTSolution_Create(RTSolution, N_LAYERS)
  CALL CRTM_Atmosphere_Create(Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS)
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  CALL CRTM_Geometry_SetValue(Geometry, Sensor_Zenith_Angle=ZENITH_ANGLE, &
                              Sensor_Scan_Angle=SCAN_ANGLE)

  Error_Status = CRTM_Forward(Atm, Sfc, Geometry, ChannelInfo, RTSolution)
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Forward failed (reference)', FAILURE)
    STOP 1
  END IF

  ALLOCATE(BT_ref(n_Channels, N_PROFILES), BT_emu(n_Channels, N_PROFILES), &
           dBT(n_Channels, N_PROFILES))
  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      BT_ref(l,m) = RTSolution(l,m)%Brightness_Temperature
    END DO
  END DO

  CALL CRTM_RTSolution_Destroy(RTSolution)
  DEALLOCATE(RTSolution)
  Error_Status = CRTM_Destroy(ChannelInfo)
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Destroy failed (reference)', FAILURE)
    STOP 1
  END IF

  !
  ! Pass B: emulator (netCDF TauCoeff)
  !
  WRITE(*,'(/5x,"Initializing CRTM with EMULATOR TauCoeff...")')
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, &
                            File_Path       = COEFFICIENTS_PATH, &
                            TauCoeff_Format = 'netCDF' )
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'Error initializing CRTM (emulator)', FAILURE)
    STOP 1
  END IF

  ALLOCATE(RTSolution(n_Channels, N_PROFILES))
  CALL CRTM_RTSolution_Create(RTSolution, N_LAYERS)
  Error_Status = CRTM_Forward(Atm, Sfc, Geometry, ChannelInfo, RTSolution)
  IF (Error_Status /= SUCCESS) THEN
    CALL Display_Message(PROGRAM_NAME, 'CRTM_Forward failed (emulator)', FAILURE)
    STOP 1
  END IF

  DO m = 1, N_PROFILES
    DO l = 1, n_Channels
      BT_emu(l,m) = RTSolution(l,m)%Brightness_Temperature
    END DO
  END DO

  ! -------------------------------------------------------------------------
  ! 3. COMPARISON METRICS
  ! -------------------------------------------------------------------------
  dBT = BT_emu - BT_ref
  max_abs = MAXVAL(ABS(dBT))
  n_valid = n_Channels * N_PROFILES
  sse     = SUM(dBT*dBT)
  rmse    = SQRT(sse / REAL(n_valid, fp))

  WRITE(*,'(/5x,"Brightness temperature comparison:")')
  WRITE(*,'(5x,"  max|dBT| = ",f10.4," K   (tol ",f6.3,")")') max_abs, BT_MAX_TOL
  WRITE(*,'(5x,"  RMSE     = ",f10.4," K   (tol ",f6.3,")")') rmse,    BT_RMSE_TOL
  WRITE(*,'(/5x,"Per-channel summary (profile 1):")')
  DO l = 1, n_Channels
    WRITE(*,'(7x,"ch ",i3,":  ref=",f8.3," K  emu=",f8.3," K  d=",f7.3," K")') &
          l, BT_ref(l,1), BT_emu(l,1), dBT(l,1)
  END DO

  CALL UnitTest_Assert(ut, max_abs < BT_MAX_TOL)
  CALL UnitTest_Assert(ut, rmse    < BT_RMSE_TOL)

  ! Cleanup
  CALL CRTM_RTSolution_Destroy(RTSolution)
  DEALLOCATE(RTSolution, BT_ref, BT_emu, dBT)
  CALL CRTM_Atmosphere_Destroy(Atm)
  Error_Status = CRTM_Destroy(ChannelInfo)

  CALL UnitTest_Report(ut)
  ok = UnitTest_Passed(ut)
  IF (.NOT. ok) STOP 1
  STOP 0

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_taucoeff_emulator
