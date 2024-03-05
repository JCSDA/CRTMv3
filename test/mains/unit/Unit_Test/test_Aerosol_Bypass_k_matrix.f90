!
! test_Simple
!
! Test program for the CRTM K-Matrix function including clouds and aerosols.
!
!

PROGRAM test_Aerosol_Bypass_k_matrix

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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_Aerosol_Bypass_k_matrix'
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
  INTEGER :: n_Channels
  INTEGER :: l, m
  INTEGER :: N_AEROSOLS, ICASE

  ! Declarations for Jacobian comparisons
  INTEGER :: n_la, n_ma
  INTEGER :: n_ls, n_ms
  CHARACTER(256) :: atmk_File, sfck_File
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)



  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)

  ! Define the FORWARD variables
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)

  ! Define the K-MATRIX variables
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atmosphere_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: Surface_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_K(:,:)
  ! ============================================================================

  !First, make sure the right number of inputs have been provided
  Sensor_Id = 'atms_npp'

  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Test program for the CRTM K-Matrix function including clouds and aerosols.', &
    'CRTM Version: '//TRIM(Version) )


  ! Get sensor id from user
  ! -----------------------
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)



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
  ALLOCATE( RTSolution( n_Channels, N_PROFILES ), &
            atm_K( n_Channels, N_PROFILES ), &
            Atmosphere_K( n_Channels, N_PROFILES ), &
            Surface_K( n_Channels, N_PROFILES ), &
            RTSolution_K( n_Channels, N_PROFILES ), &
            STAT = Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  DO ICASE = 1, 2

    ! --------------------------------
    IF (ICASE == 1) THEN
      N_AEROSOLS = 2
    END IF
    IF (ICASE == 2) THEN
      N_AEROSOLS = 3
    END IF

  ! 3b. Allocate the STRUCTURES
  ! ---------------------------
  ! The input FORWARD structure
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! The output K-MATRIX structure
  CALL CRTM_Atmosphere_Create( Atmosphere_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atmosphere_K)) ) THEN
    Message = 'Error allocating CRTM Atmosphere_K structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! Deleted in V2.4.1
  ! The comparative K-MATRIX structure inside the results file
  !CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  !IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
  !  Message = 'Error allocating CRTM atm_K structure'
  !  CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !  STOP 1
  !END IF

  ! ============================================================================




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
  ! 5. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
  !
  ! 5a. Zero the K-matrix OUTPUT structures
  ! ---------------------------------------
  CALL CRTM_Atmosphere_Zero( Atmosphere_K )
  !CALL CRTM_Atmosphere_Zero( atm_K ) ! deleted in V2.4.1
  CALL CRTM_Surface_Zero( Surface_K )

  ! 5b. Inintialize the K-matrix INPUT so
  !     that the IR/MW results are dTb/dx
  !     and the visible results are dR/dx
  ! -------------------------------------
  DO l = 1, n_Channels
    IF ( ChannelInfo(1)%Sensor_Type == INFRARED_SENSOR .OR. &
         ChannelInfo(1)%Sensor_Type == MICROWAVE_SENSOR ) THEN
      RTSolution_K(l,:)%Radiance               = ZERO
      RTSolution_K(l,:)%Brightness_Temperature = ONE
    ELSE
      RTSolution_K(l,:)%Radiance               = ONE
      RTSolution_K(l,:)%Brightness_Temperature = ZERO
    END IF
  END DO
  ! ============================================================================




  ! ============================================================================
  ! 6. **** CALL THE CRTM K-MATRIX MODEL ****
  !
  Error_Status = CRTM_K_Matrix( Atm         , &
                                Sfc         , &
                                RTSolution_K, &
                                Geometry    , &
                                ChannelInfo , &
                                Atmosphere_K, &
                                Surface_K   , &
                                RTSolution  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM K_Matrix Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================




  ! ============================================================================
  ! 7. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  DO m = 1, N_PROFILES
    WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)
    DO l = 1, n_Channels
      WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
      ! FWD output
      WRITE( *, '(/3x,"FORWARD OUTPUT")')
      CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
      ! K-MATRIX output
      WRITE( *, '(/3x,"K-MATRIX OUTPUT")')
      CALL CRTM_Surface_Inspect(Surface_K(l,m))
      CALL CRTM_Atmosphere_Inspect(Atmosphere_K(l,m))
    END DO
  END DO
  ! ============================================================================







  ! ============================================================================
  ! 9. **** COMPARE Atmosphere_K and Surface_K RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 9a. Create the output files if they do not exist
  ! ------------------------------------------------
  ! 9a.1 Atmosphere file
  ! ...Generate filename
  atmk_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Atmosphere.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(atmk_File) ) THEN
    Message = 'Atmosphere_K save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write Atmosphere_K structure to file
    Error_Status = CRTM_Atmosphere_WriteFile( atmk_file, Atmosphere_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating Atmosphere_K save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF
  ! 9a.2 Surface file
  ! ...Generate filename
  sfck_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Surface.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(sfck_File) ) THEN
    Message = 'Surface_K save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write Surface_K structure to file
    Error_Status = CRTM_Surface_WriteFile( sfck_file, Surface_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating Surface_K save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

   ! IF (ICASE == 1) THEN
   !   DEALLOCATE(RTSolution, RTSolution_K, &
   !              STAT = Allocate_Status)
   ! END IF

  END DO !ICASE

  ! 9b. Inquire the saved files
  ! ---------------------------
  ! 9b.1 Atmosphere file
  Error_Status = CRTM_Atmosphere_InquireFile( atmk_File, &
                                              n_Channels = n_la, &
                                              n_Profiles = n_ma )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring Atmosphere_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! 9b.2 Surface file
  Error_Status = CRTM_Surface_InquireFile( sfck_File, &
                                           n_Channels = n_ls, &
                                           n_Profiles = n_ms )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring Surface_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9c. Compare the dimensions
  ! --------------------------
  IF ( n_la /= n_Channels .OR. n_ma /= N_PROFILES .OR. &
       n_ls /= n_Channels .OR. n_ms /= N_PROFILES      ) THEN
    Message = 'Dimensions of saved data different from that calculated!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9d. Read the saved data
  ! -----------------------
  ! 9d.1 Atmosphere file
  Error_Status = CRTM_Atmosphere_ReadFile( atmk_File, atm_k, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading Atmosphere_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! 9d.2 Surface file
  Error_Status = CRTM_Surface_ReadFile( sfck_File, sfc_k, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading Surface_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9e. Compare some Jacobians
  ! --------------------------
  ! 9e.1 Atmosphere
  ! IF ( ALL(CRTM_Atmosphere_Compare(Atmosphere_K, atm_k)) ) THEN
  !   Message = 'Atmosphere_K Jacobians are the same!'
  !   CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ! ELSE
  !   Message = 'Atmosphere_K Jacobians are different!'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   ! Write the current Atmosphere_K results to file
  !   atmk_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Atmosphere.bin'
  !   Error_Status = CRTM_Atmosphere_WriteFile( atmk_file, Atmosphere_K, Quiet=.TRUE. )
  !   IF ( Error_Status /= SUCCESS ) THEN
  !     Message = 'Error creating temporary Atmosphere_K save file for failed comparison'
  !     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   END IF
  ! END IF
  ! 9e.2 Surface
  IF ( ALL(CRTM_Surface_Compare(Surface_K, sfc_k, n_SigFig=5)) ) THEN
    Message = 'Surface_K Jacobians are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'Surface_K Jacobians are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current Surface_K results to file
    sfck_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Surface.bin'
    Error_Status = CRTM_Surface_WriteFile( sfck_file, Surface_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary Surface_K save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
  END IF
  ! ============================================================================

  ! ============================================================================
  ! 8. **** DESTROY THE CRTM ****
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
  ! 10. **** CLEAN UP ****
  !
  ! 10a. Deallocate the structures
  ! ------------------------------
  CALL CRTM_Atmosphere_Destroy(atm_K)
  CALL CRTM_Atmosphere_Destroy(Atmosphere_K)
  CALL CRTM_Atmosphere_Destroy(Atm)




  ! 10b. Deallocate the arrays
  ! --------------------------
  DEALLOCATE(RTSolution, RTSolution_K, &
             Surface_K, Atmosphere_K, &
             sfc_K, atm_K, &
             STAT = Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'
  INCLUDE 'Load_Atm_Data_Bypass_Aerosol.inc'

END PROGRAM test_Aerosol_Bypass_k_matrix
