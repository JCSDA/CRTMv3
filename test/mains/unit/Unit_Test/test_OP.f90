!
! test_OP
!
! Test program for the CRTM Forward function with user defined input optical profiles



PROGRAM test_OP

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ==========================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_OP'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/unit/'
  CHARACTER(*), PARAMETER :: TOP_FILE = 'TOP_SingleProfile.nc'

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 1
  INTEGER, PARAMETER :: N_AEROSOLS  = 1
  ! ...but only ONE Sensor at a time
  INTEGER, PARAMETER :: N_SENSORS = 1

  ! Test GeometryInfo angles. The test scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! ============================================================================

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message
  CHARACTER(256) :: Version
  CHARACTER(256) :: Sensor_Id
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels, n_Stokes
  INTEGER :: l, m
  ! Declarations for RTSolution comparison
  INTEGER :: n_l, n_m, n_k, n_s
  CHARACTER(256) :: rts_File
  CHARACTER(256) :: op_File
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  TYPE(OP_Input_type), ALLOCATABLE :: OP

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_OP(:,:)

  ! Define Options
  TYPE(CRTM_Options_type)                 :: Options(N_PROFILES)
  ! ============================================================================


  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Test program for the CRTM Forward function with user defined aerosol optical profiles.', &
    'CRTM Version: '//TRIM(Version) )

  ! Sensor_Id
  Sensor_Id = 'v.abi_gr'
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)
  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. Initialise the requested sensor
  ! -----------------------------------
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  ! ============================================================================

  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
  ! 3a. Allocate the ARRAYS
  ! -----------------------
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ALLOCATE( RTSolution_OP( n_Channels, N_PROFILES ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3a-2. Allocate N_Layers for layered outputs
  CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  CALL CRTM_RTSolution_Create( RTSolution_OP, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3b. Allocate the Structures
  ! ---------------------------
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  CALL CRTM_Options_Create( Options, n_Channels )
  IF ( ANY(.NOT. CRTM_Options_Associated(Options)) ) THEN
    Message = 'Error allocating CRTM Options structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! ============================================================================


  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! 4a. Atmosphere and Surface input
  ! --------------------------------
  CALL Load_Atm_Data_SingleProfile()
  CALL Load_Sfc_Data_SingleProfile()


  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )

  ! 4c. Optional varibles
  CALL CRTM_Options_SetValue( Options, &
                              n_Streams     = 6       , &
                              Use_Total_OP  = .TRUE.    )
  op_File = COEFFICIENTS_PATH//TOP_FILE

  !op_File = '/Users/dangch/Documents/CRTM/CRTM_dev/crtm_code_review/CRTMv3_ITF/teST/mains/unit/Unit_Test/top_test.nc'
  print *, 'start reading optical profile from file : ', op_File
  Error_Status = OP_Input_ReadFile(op_File, OP)
  IF ( Error_Status /= SUCCESS ) THEN
     Message = 'Error reading OP file'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
  END IF

  ! ...set total op values
  Options(1)%TOP%tau    = OP%tau
  Options(1)%TOP%bs     = OP%bs
  Options(1)%TOP%kb     = OP%kb
  Options(1)%TOP%pcoeff = OP%pcoeff

  !PRINT *, 'RT_Algorithm_Id', Options(1)%RT_Algorithm_Id
  !print*, 'TOP%n_Legendre_Terms', OP%n_Legendre_Terms

  !PRINT *,   Options(1)%TOP%tau
  !PRINT *,   Options(1)%n_Stokes


  ! ============================================================================

  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  ! CRTM Default Interface
    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! n_Stokes = RTSolution(1,1)%n_Stokes+1
    PRINT *, 'FINISH CRTM Default Interface'

    ! CRTM TOP Interface
    Error_Status = CRTM_Forward( Atm        , &
                                 Sfc        , &
                                 Geometry   , &
                                 ChannelInfo, &
                                 RTSolution_OP, &
                                 Options )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error in CRTM Forward Model with Options'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! n_Stokes = RTSolution_OP(1,1)%n_Stokes+1
    PRINT *, 'CRTM TOP Interface'


  ! ============================================================================

  ! ============================================================================
  ! 8. **** COMPARE RTSolution RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! ! 8a. Create the output file if it does not exist
  ! ! -----------------------------------------------
  ! ! ...Generate a filename
  ! rts_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.nc'
  ! ! ...Check if the file exists
  ! IF ( .NOT. File_Exists(rts_File) ) THEN
  !   Message = 'RTSolution save file does not exist. Creating...'
  !   CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  !   ! ...File not found, so write RTSolution structure to file
  !   Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, NetCDF=.TRUE., Quiet=.TRUE. )
  !   IF ( Error_Status /= SUCCESS ) THEN
  !     Message = 'Error creating RTSolution save file'
  !     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !     STOP 1
  !   END IF
  ! END IF


  ! 8b. Inquire the saved file
  ! ! --------------------------
  ! Error_Status = CRTM_RTSolution_InquireFile( rts_File, &
  !                                             NetCDF=.TRUE.,    &
  !                                             n_Profiles = n_m, &
  !                                             n_Layers   = n_k, &
  !                                             n_Channels = n_l, &
  !                                             n_Stokes   = n_s )
  ! IF ( Error_Status /= SUCCESS ) THEN
  !   Message = 'Error inquiring RTSolution save file'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF

  ! 8c. Compare the dimensions
  ! --------------------------
  ! IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES .OR. n_k /= n_Layers .OR. n_s /= n_Stokes ) THEN
  !   Message = 'Dimensions of saved data different from that calculated!'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF


  ! 8d. Allocate the structure to read in saved data
  ! ------------------------------------------------
  ! ALLOCATE( rts( n_l, n_m ), STAT=Allocate_Status )
  ! IF ( Allocate_Status /= 0 ) THEN
  !   Message = 'Error allocating RTSolution saved data array'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF
  ! CALL CRTM_RTSolution_Create( rts, n_k )
  ! IF ( ANY(.NOT. CRTM_RTSolution_Associated(rts)) ) THEN
  !   Message = 'Error allocating CRTM RTSolution structures'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF
  !
  !
  ! ! 8e. Read the saved data
  ! ! -----------------------
  ! Error_Status = CRTM_RTSolution_ReadFile( rts_File, rts, NetCDF=.TRUE., Quiet=.TRUE. )
  ! IF ( Error_Status /= SUCCESS ) THEN
  !   Message = 'Error reading RTSolution save file'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF


  ! 8f. Compare the structures
  ! --------------------------
  IF ( ALL(CRTM_RTSolution_Compare(RTSolution, RTSolution_OP)) ) THEN
    Message = 'RTSolution results are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'RTSolution results are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current RTSolution results to file
    rts_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.nc'
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, NetCDF=.TRUE., Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary RTSolution save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    STOP 1
  END IF


  ! ============================================================================

  ! ============================================================================
  ! 7. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  Error_Status = CRTM_Destroy( ChannelInfo )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 9. **** CLEAN UP ****
  !
  ! 9a. Deallocate the structures
  ! -----------------------------
  CALL CRTM_Atmosphere_Destroy(Atm)

  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(RTSolution, RTSolution_OP, STAT=Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data_SingleProfile.inc'
  INCLUDE 'Load_Sfc_Data_SingleProfile.inc'

END PROGRAM test_OP
