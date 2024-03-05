!
! test_function CRTM_Atmosphere_Aerosol_Bypass
!
! Test program for the CRTM Forward function including clouds and aerosols.
!
!

PROGRAM test_Aerosol_Bypass

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_Aerosol_Bypass'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/unit/'

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 1

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
  INTEGER :: N_AEROSOLS, ICASE
  ! Declarations for RTSolution comparison
  INTEGER :: n_l, n_m, n_k, n_s
  CHARACTER(256) :: rts_File
  ! CHARACTER(*), PARAMETER :: Outcase(2) = (/'case1', 'case2'/)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)
  ! ============================================================================

  !First, make sure the right number of inputs have been provided
  Sensor_Id = 'atms_npp'

  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Test program for the CRTM Forward function including clouds and aerosols.', &
    'CRTM Version: '//TRIM(Version) )


  ! Get sensor id from user
  ! -----------------------
  Sensor_Id = ADJUSTL(Sensor_Id)
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


   DO ICASE = 1, 2

     ! --------------------------------
     IF (ICASE == 1) THEN
       N_AEROSOLS = 2
     END IF
     IF (ICASE == 2) THEN
       N_AEROSOLS = 3
     END IF

    CALL CRTM_Atmosphere_Destroy(Atm)
    DEALLOCATE(RTSolution, STAT=Allocate_Status)

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

   ! 3a-2. Allocate N_Layers for layered outputs
   CALL CRTM_RTSolution_Create( RTSolution, N_LAYERS )
   IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution)) ) THEN
     Message = 'Error allocating CRTM RTSolution structures'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
   END IF

   ! 3b. Allocate the STRUCTURES
   ! ---------------------------
   CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
   IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
     Message = 'Error allocating CRTM Atmosphere structures'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
   END IF
  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! 4a. Atmosphere and Surface input
  ! --------------------------------
  IF (ICASE == 1) THEN
    CALL Load_Atm_Data()
  END IF
  IF (ICASE == 2) THEN
    CALL Load_Atm_Data_Bypass_Aerosol()
  END IF

  CALL Load_Sfc_Data()


  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value
  !  The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )
  ! ============================================================================




  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  Error_Status = CRTM_Forward( Atm        , &
                               Sfc        , &
                               Geometry   , &
                               ChannelInfo, &
                               RTSolution  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  n_Stokes = RTSolution(1,1)%n_Stokes+1
  ! ============================================================================




  ! ============================================================================
  ! 6. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  ! DO m = 1, N_PROFILES
  !   WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id)
  !   DO l = 1, n_Channels
  !     WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
  !     CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
  !   END DO
  ! END DO
  ! ============================================================================

  ! ============================================================================
  ! 8. **** COMPARE RTSolution RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 8a. Create the output file if it does not exist
  ! -----------------------------------------------
  ! ...Generate a filename
  rts_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution.nc'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(rts_File) ) THEN
    Message = 'RTSolution save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write RTSolution structure to file
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution, NetCDF=.TRUE., Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating RTSolution save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  END DO !ICASE

  ! 8b. Inquire the saved file
  ! --------------------------
  Error_Status = CRTM_RTSolution_InquireFile( rts_File, &
                                              NetCDF=.TRUE.,    &
                                              n_Profiles = n_m, &
                                              n_Layers   = n_k, &
                                              n_Channels = n_l, &
                                              n_Stokes   = n_s )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF


  ! 8c. Compare the dimensions
  ! --------------------------
  IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES .OR. n_k /= n_Layers .OR. n_s /= n_Stokes ) THEN
    Message = 'Dimensions of saved data different from that calculated!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF


  ! 8d. Allocate the structure to read in saved data
  ! ------------------------------------------------
  ALLOCATE( rts( n_l, n_m ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating RTSolution saved data array'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  CALL CRTM_RTSolution_Create( rts, n_k )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(rts)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF


  ! 8e. Read the saved data
  ! -----------------------
  Error_Status = CRTM_RTSolution_ReadFile( rts_File, rts, NetCDF=.TRUE., Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading RTSolution save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF


  ! 8f. Compare the structures
  ! --------------------------
  IF ( ALL(CRTM_RTSolution_Compare(RTSolution, rts)) ) THEN
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
  DEALLOCATE(RTSolution, rts, STAT=Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'
  INCLUDE 'Load_Atm_Data_Bypass_Aerosol.inc'

END PROGRAM test_Aerosol_Bypass
