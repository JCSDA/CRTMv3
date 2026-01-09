!
! test_OMPoverChannels
!
! Test program for the CRTM Tangent-Linear function including clouds and aerosols.
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
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/tangent_linear/'
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
  ! Declarations for RTSolution comparison
  INTEGER :: n_l, n_m, i
  CHARACTER(256) :: rts_File
  INTEGER, ALLOCATABLE :: channel_subset(:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_TL(:,:)


  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES), Atm_TL(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES), Sfc_TL(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_TL(:,:)
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
    'Test program for the CRTM Tangent-Linear function including clouds and aerosols.', &
    'CRTM Version: '//TRIM(Version) )


  ! Get sensor id from user
  ! -----------------------
  Sensor_Id = ADJUSTL(Sensor_Id)
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)

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
            RTSolution_TL( n_Channels, N_PROFILES ), &
            STAT=Allocate_Status )
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
    Message = 'Error allocating CRTM Atmosphere FWD structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! The input TANGENT-LINEAR structure
  CALL CRTM_Atmosphere_Create( Atm_TL, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm_TL)) ) THEN
    Message = 'Error allocating CRTM Atmosphere TL structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  ! ============================================================================




  ! ============================================================================
  ! 4. **** ASSIGN INPUT DATA ****
  !
  ! 4a. Atmosphere and Surface FWD input
  ! ------------------------------------
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
  ! 5. **** INITIALIZE THE TANGENT-LINEAR ARGUMENTS ****
  !
  ! 5a. Zero the tangent-liner INPUT structures
  ! ---------------------------------------
  ! Copy...
  Atm_TL = Atm
  ! ...zero...
  CALL CRTM_Atmosphere_Zero(Atm_TL)
  ! ...and perturb temperature by 0.5K
  DO m = 1, N_PROFILES
    Atm_TL(m)%Temperature = 0.5_fp
  END DO

  ! Copy...
  Sfc_TL = Sfc
  ! ...and zero.
  CALL CRTM_Surface_Zero(Sfc_TL)
  ! ============================================================================




  ! ============================================================================
  ! 6. **** CALL THE CRTM TANGENT-LINEAR MODEL ****
  !
  CALL SYSTEM_CLOCK( clock_start, clock_rate )
  IF ( clock_rate == 0 ) clock_rate = 1
  Error_Status = CRTM_Tangent_Linear( Atm          , &
                                      Sfc          , &
                                      Atm_TL       , &
                                      Sfc_TL       , &
                                      Geometry     , &
                                      ChannelInfo  , &
                                      RTSolution   , &
                                      RTSolution_TL , &
                                      Options=Options )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Tangent-Linear Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  CALL SYSTEM_CLOCK( clock_end, clock_rate )
  IF ( clock_rate == 0 ) clock_rate = 1
  WRITE( *,'(/5x,"CRTM_Tangent_Linear wall time (s): ",f10.3)' ) &
    REAL(clock_end - clock_start, fp) / REAL(clock_rate, fp)
  ! ============================================================================




  ! ============================================================================
  ! 7. **** OUTPUT THE RESULTS TO SCREEN ****
  !
  DO m = 1, N_PROFILES
    WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id)
    DO l = 1, n_Channels
      WRITE( *, '(/5x,"Channel ",i0," results")') RTSolution(l,m)%Sensor_Channel
      ! FWD output
      WRITE( *, '(/3x,"FORWARD OUTPUT")')
      CALL CRTM_RTSolution_Inspect(RTSolution(l,m))
      ! TL output
      WRITE( *, '(/3x,"TANGENT-LINEAR OUTPUT")')
      CALL CRTM_RTSolution_Inspect(RTSolution_TL(l,m))
    END DO
  END DO
  ! ============================================================================








  ! ============================================================================
  ! 9. **** COMPARE RTSolution_TL RESULTS TO SAVED VALUES ****
  !
  WRITE( *, '( /5x, "Comparing calculated results with saved ones..." )' )

  ! 9a. Create the output file if it does not exist
  ! -----------------------------------------------
  ! ...Generate a filename
  rts_File = RESULTS_PATH//TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.RTSolution_TL.bin'
  ! ...Check if the file exists
  IF ( .NOT. File_Exists(rts_File) ) THEN
    Message = 'RTSolution_TL save file does not exist. Creating...'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
    ! ...File not found, so write RTSolution_TL structure to file
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution_TL, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating RTSolution_TL save file'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
  END IF

  ! 9b. Inquire the saved file
  ! --------------------------
  Error_Status = CRTM_RTSolution_InquireFile( rts_File, &
                                              n_Channels = n_l, &
                                              n_Profiles = n_m )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error inquiring RTSolution_TL save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9c. Compare the dimensions
  ! --------------------------
  IF ( n_l /= n_Channels .OR. n_m /= N_PROFILES ) THEN
    Message = 'Dimensions of saved data different from that calculated!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9d. Allocate the structure to read in saved data
  ! ------------------------------------------------
  ALLOCATE( rts_TL( n_l, n_m ), STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating RTSolution_TL saved data array'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9e. Read the saved data
  ! -----------------------
  Error_Status = CRTM_RTSolution_ReadFile( rts_File, rts_TL, Quiet=.TRUE. )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error reading RTSolution_TL save file'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 9f. Compare the structures
  ! --------------------------
  IF ( ALL(CRTM_RTSolution_Compare(RTSolution_TL, rts_TL)) ) THEN
    Message = 'RTSolution_TL results are the same!'
    CALL Display_Message( PROGRAM_NAME, Message, INFORMATION )
  ELSE
    Message = 'RTSolution_TL results are different!'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    ! Write the current RTSolution results to file
    rts_File = TRIM(Sensor_Id)//'.RTSolution_TL.bin'
    Error_Status = CRTM_RTSolution_WriteFile( rts_File, RTSolution_TL, Quiet=.TRUE. )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error creating temporary RTSolution_TL save file for failed comparison'
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
  CALL CRTM_Atmosphere_Destroy(Atm)
  CALL CRTM_Atmosphere_Destroy(Atm_TL)

  ! 9b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(RTSolution, RTSolution_TL, rts_TL, channel_subset, STAT=Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data.inc'
  INCLUDE 'Load_Sfc_Data.inc'

END PROGRAM test_OMPoverChannels
