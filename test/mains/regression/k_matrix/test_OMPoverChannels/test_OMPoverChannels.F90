!
! test_OMPoverChannels
!
! Test program for the CRTM K-Matrix function including clouds and aerosols.
!
!

PROGRAM test_OMPoverChannels

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
#ifdef _OPENMP
  USE OMP_LIB
#endif
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'test_OMPoverChannels'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/k_matrix/'
  INTEGER, PARAMETER :: TARGET_CHANNEL_COUNT = 256


  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions...
  INTEGER, PARAMETER :: N_PROFILES  = 2
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
  CHARACTER(32) :: threads_arg
  INTEGER :: Error_Status
  INTEGER :: env_status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels
  INTEGER :: n_Channels_Total
  INTEGER :: subset_stride
  INTEGER :: subset_count
  INTEGER :: l, m
  INTEGER :: obs_level
  INTEGER :: arg_count
  INTEGER :: n_threads
  INTEGER :: ios
  INTEGER :: clock_start
  INTEGER :: clock_end
  INTEGER :: clock_rate
  ! Declarations for Jacobian comparisons
  INTEGER :: n_la, n_ma
  INTEGER :: n_ls, n_ms
  INTEGER :: n_l, n_m, i
  CHARACTER(256) :: atmk_File, sfck_File, rtsk_File
  INTEGER, ALLOCATABLE :: channel_subset(:)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_k(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_k(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_k(:,:)



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
  TYPE(CRTM_Options_type)                :: Options(N_PROFILES)
  ! ============================================================================

  !First, make sure the right number of inputs have been provided
  arg_count = COMMAND_ARGUMENT_COUNT()
  IF ( arg_count < 1 .OR. arg_count > 2 ) THEN
     WRITE(*,*) TRIM(PROGRAM_NAME)//': ERROR, 1 or 2 command-line arguments required, returning'
     STOP 1
  END IF
  CALL GET_COMMAND_ARGUMENT(1,Sensor_Id)   !read in the value
  IF ( arg_count == 2 ) THEN
    CALL GET_COMMAND_ARGUMENT(2,threads_arg)
    READ(threads_arg,*,IOSTAT=ios) n_threads
    IF ( ios /= 0 .OR. n_threads < 1 ) THEN
      WRITE(*,*) TRIM(PROGRAM_NAME)//': ERROR, invalid OpenMP thread count, returning'
      STOP 1
    END IF
  END IF


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

  ! Set number of OpenMP threads, if present:
  ! -----------------------------------------
#ifdef _OPENMP
  CALL GET_ENVIRONMENT_VARIABLE("OMP_NUM_THREADS", STATUS=env_status)
  IF (env_status /= 0) THEN
     CALL OMP_SET_NUM_THREADS(1)
  END IF
  IF ( arg_count == 2 ) THEN
    ! CALL OMP_SET_NUM_THREADS(n_threads)
  END IF
!$OMP PARALLEL
!$OMP SINGLE
  WRITE(*,'(5x,"OpenMP threads: ",i0)') OMP_GET_NUM_THREADS()
!$OMP END SINGLE
!$OMP END PARALLEL
#endif



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

  ! 2b. Specify an instrument channel subset for processing
  ! -------------------------------------------------------
  ! n_Channels_Total = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  ! subset_stride = MAX(1, n_Channels_Total / TARGET_CHANNEL_COUNT)
  ! subset_count = ((n_Channels_Total - 1) / subset_stride) + 1
  ! ALLOCATE( channel_subset(subset_count), STAT=Allocate_Status )
  ! IF ( Allocate_Status /= 0 ) THEN
  !   Message = 'Error allocating channel subset list'
  !   CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
  !   STOP 1
  ! END IF
  ! channel_subset = [(l, l=1, n_Channels_Total, subset_stride)]
  Error_Status = CRTM_ChannelInfo_Subset( ChannelInfo(1), &
                                          Channel_Subset = (/ (i, i = 1, 16920, 2) /) )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Selecting channel subset unsuccessful!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 2c. Determine the total number of channels
  !     for which the CRTM was initialized
  ! ------------------------------------------
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  WRITE( *,'(/5x,"Processing ",i0," of ",i0," channels (stride ",i0,")")' ) &
    n_Channels, n_Channels_Total, subset_stride
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
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  obs_level = N_LAYERS / 2
  DO m = 1, N_PROFILES
    Options(m)%Obs_4_downward_P = Atm(m)%Level_Pressure(obs_level)
  END DO


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
    IF ( Options(1)%Obs_4_downward_P > ZERO ) THEN
      RTSolution_K(l,:)%Radiance               = ZERO
      RTSolution_K(l,:)%Down_Radiance          = ONE
      RTSolution_K(l,:)%Brightness_Temperature = ZERO
    ELSE IF ( ChannelInfo(1)%Sensor_Type == INFRARED_SENSOR .OR. &
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
  CALL SYSTEM_CLOCK( clock_start, clock_rate )
  IF ( clock_rate == 0 ) clock_rate = 1
  Error_Status = CRTM_K_Matrix( Atm         , &
                                Sfc         , &
                                RTSolution_K, &
                                Geometry    , &
                                ChannelInfo , &
                                Atmosphere_K, &
                                Surface_K   , &
                                RTSolution  , &
                                Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM K_Matrix Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  CALL SYSTEM_CLOCK( clock_end, clock_rate )
  IF ( clock_rate == 0 ) clock_rate = 1
  WRITE( *,'(/5x,"CRTM_K_Matrix wall time (s): ",f10.3)' ) &
    REAL(clock_end - clock_start, fp) / REAL(clock_rate, fp)
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
      CALL CRTM_RTSolution_Inspect(RTSolution_K(l,m))
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
  ! 9a.3 RTSolution_K file
  ! ...Generate filename
  rtsk_file = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution_K.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(rtsk_file) ) THEN
    Message = 'RTSolution_K save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write RTSolution_K structure to file
    Error_Status = CRTM_RTSolution_WriteFile( rtsk_file, RTSolution_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating RTSolution_K save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

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
  ! 9b.3 RTSolution_K file
  Error_Status = CRTM_RTSolution_InquireFile(rtsk_file, &
                                             n_Channels = n_l, &
                                             n_Profiles = n_m )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring RTSolution_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9c. Compare the dimensions
  ! --------------------------
  IF ( n_la /= n_Channels .OR. n_ma /= N_PROFILES .OR. &
       n_ls /= n_Channels .OR. n_ms /= N_PROFILES .OR. &
       n_l  /= n_Channels .OR. n_m  /= N_PROFILES ) THEN
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
  ! 9d.3 RTSolution_K file
  ALLOCATE( rts_k( n_l, n_m ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating RTSolution_K saved data array'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  Error_Status = CRTM_RTSolution_ReadFile( rtsk_file, rts_k, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading RTSolution_K save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9e. Compare some Jacobians
  ! --------------------------
  ! 9e.1 Atmosphere
  IF ( ALL(CRTM_Atmosphere_Compare(Atmosphere_K, atm_k)) ) THEN
    Message = 'Atmosphere_K Jacobians are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'Atmosphere_K Jacobians are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current Atmosphere_K results to file
    atmk_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.Atmosphere.bin'
    Error_Status = CRTM_Atmosphere_WriteFile( atmk_file, Atmosphere_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary Atmosphere_K save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    STOP 1
  END IF
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
    STOP 1
  END IF
  ! 9e.3 RTSolution_K
  IF ( ALL(CRTM_RTSolution_Compare(RTSolution_K, rts_k, n_SigFig=5)) ) THEN
    Message = 'RTSolution_K results are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'RTSolution_K results are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    rtsk_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution_K.bin'
    Error_Status = CRTM_RTSolution_WriteFile( rtsk_File, RTSolution_K, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary RTSolution_K save file for failed comparison'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    END IF
    STOP 1
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
             sfc_k, atm_k, rts_k, channel_subset, &
             STAT = Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_OMPoverChannels
