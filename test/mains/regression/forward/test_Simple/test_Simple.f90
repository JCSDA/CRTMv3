!
! test_Simple
!
! Test program for the CRTM Forward function including clouds and aerosols.
!
!

PROGRAM test_Simple

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_Simple'

  ! Directory location of coefficients
  CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='./testinput/'
  CHARACTER(*), PARAMETER :: NC_COEFFICIENT_PATH='./testinput/'

  ! Aerosol/Cloud coefficient format
  CHARACTER(*), PARAMETER :: Coeff_Format = 'Binary'
  !CHARACTER(*), PARAMETER :: Coeff_Format = 'netCDF'
  
  ! Aerosol/Cloud coefficient scheme
  CHARACTER(*), PARAMETER :: Aerosol_Model = 'CRTM'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'CMAQ'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'GOCART-GEOS5'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'NAAPS'
  CHARACTER(*), PARAMETER :: Cloud_Model   = 'CRTM'
  
  ! Directory location of results for comparison [NOT USED YET]
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/forward/'
  

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 1
  INTEGER, PARAMETER :: N_AEROSOLS  = 1


  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS =1 

  ! Some pretend geometry angles. The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp
  ! ============================================================================

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: Message, Version
  CHARACTER(256) :: AerosolCoeff_File
  CHARACTER(256) :: AerosolCoeff_Format
  CHARACTER(256) :: CloudCoeff_File
  CHARACTER(256) :: CloudCoeff_Format
  CHARACTER(256) :: Aerosol_Scheme
  CHARACTER(256) :: Cloud_Scheme
  INTEGER        :: Error_Status, Allocate_Status
  INTEGER        :: n_Channels
  INTEGER        :: l, m, n, nc, n_l, n_m
  CHARACTER(256) :: rts_File, Sensor_Id
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 1a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)


  ! 1b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)



  !First, make sure the right number of inputs have been provided
  IF(COMMAND_ARGUMENT_COUNT().NE.1)THEN
     WRITE(*,*)'test_Simple.f90: ERROR, only one command-line argument required, returning'
     STOP 1
  ENDIF
  CALL GET_COMMAND_ARGUMENT(1,Sensor_Id)   !read in the value

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
  ! STEP 2. **** INITIALIZE THE CRTM ****
  !
  ! 2a. Initialise all the sensors at once
  ! --------------------------------------
  ! ... Cloud coefficient information
  IF ( Cloud_Model /= 'CRTM' ) THEN
     Cloud_Scheme = Cloud_Model//'.'
  ELSE
     Cloud_Scheme = ' '
  END IF
  ! ... Aerosol coefficient information
  IF ( Aerosol_Model /= 'CRTM' ) THEN
     Aerosol_Scheme = Aerosol_Model//'.'
  ELSE
     Aerosol_Scheme = ' '
  END IF
  ! ... Coefficient table format
  IF ( Coeff_Format == 'Binary' ) THEN
     AerosolCoeff_Format = 'Binary'
     AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'bin'
     CloudCoeff_Format   = 'Binary'
     CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'bin'
  ELSE IF ( Coeff_Format == 'netCDF' ) THEN
     AerosolCoeff_Format = 'netCDF'
     AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'nc4'
     CloudCoeff_Format   = 'netCDF'
     CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'nc4'
  END IF
 
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                        ChannelInfo, &
                        Aerosol_Model, &
                        AerosolCoeff_Format, &
                        AerosolCoeff_File, &
                        Cloud_Model, &
                        CloudCoeff_Format, &
                        CloudCoeff_File, &
                        File_Path=COEFFICIENT_PATH, &
                        NC_File_Path=NC_COEFFICIENT_PATH, &
                        Quiet=.TRUE.)

  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 2b. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_Channels
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
  CALL CRTM_RTSolution_Create(RTSolution,N_LAYERS)

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




  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! 4a. Atmosphere and Surface input
  ! --------------------------------
  CALL Load_Atm_Data()
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
  CALL CRTM_RTSolution_Destroy(RTSolution)
  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(rts, STAT=Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_Simple
