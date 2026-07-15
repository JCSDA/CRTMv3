!
! test_OP
!
! Test program for the CRTM Forward function with user defined input optical
! profiles. It checks that supplying externally-computed aerosol-only (AOP),
! cloud-only (COP), and combined (TOP) optical profiles via CRTM_Options
! produces RTSolution results identical to the default CRTM_Forward
! calculation, which computes the cloud/aerosol scattering internally.
! The AOP/COP/TOP netCDF input files are produced by Generate_OP.f90.



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
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'test_OP'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: RESULTS_PATH      = './results/unit/'
  CHARACTER(*), PARAMETER :: AOP_FILE = 'AOP_SingleProfile.nc'
  CHARACTER(*), PARAMETER :: COP_FILE = 'COP_SingleProfile.nc'
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
  CHARACTER(256) :: op_File
  INTEGER :: Error_Status
  INTEGER :: Allocate_Status
  INTEGER :: n_Channels

  ! Optical profiles read back from the netCDF files produced by Generate_OP.f90
  TYPE(OP_Input_type), ALLOCATABLE :: OP_AOP
  TYPE(OP_Input_type), ALLOCATABLE :: OP_COP
  TYPE(OP_Input_type), ALLOCATABLE :: OP_TOP

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)             :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)              :: Atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: Sfc(N_PROFILES)

  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:)      ! Default (no Options)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_AOP(:,:)  ! Options%AOP
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_COP(:,:)  ! Options%COP
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution_TOP(:,:)  ! Options%TOP

  ! Define Options: one set per optical-profile source, so each carries its
  ! own Use_*_OP flag independently of the others. Options_Default carries
  ! only the fixed n_Streams=16 setting (no Use_*_OP flags), so the "default"
  ! run below uses the SAME stream count as the AOP/COP/TOP runs and as
  ! Generate_OP.f90 - all four calculations must agree on n_Streams for the
  ! comparison in section 6 to be meaningful.
  TYPE(CRTM_Options_type)                 :: Options_Default(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: Options_AOP(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: Options_COP(N_PROFILES)
  TYPE(CRTM_Options_type)                 :: Options_TOP(N_PROFILES)
  ! ============================================================================


  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Test program for the CRTM Forward function with user defined aerosol/cloud/total optical profiles.', &
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
  ALLOCATE( RTSolution(     n_Channels, N_PROFILES ), &
            RTSolution_AOP( n_Channels, N_PROFILES ), &
            RTSolution_COP( n_Channels, N_PROFILES ), &
            RTSolution_TOP( n_Channels, N_PROFILES ), &
            STAT=Allocate_Status )
  IF ( Allocate_Status /= 0 ) THEN
    Message = 'Error allocating structure arrays'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3a-2. Allocate N_Layers for layered outputs
  CALL CRTM_RTSolution_Create( RTSolution,     N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_AOP, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_COP, N_LAYERS )
  CALL CRTM_RTSolution_Create( RTSolution_TOP, N_LAYERS )
  IF ( ANY(.NOT. CRTM_RTSolution_Associated(RTSolution))     .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_AOP)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_COP)) .OR. &
       ANY(.NOT. CRTM_RTSolution_Associated(RTSolution_TOP)) ) THEN
    Message = 'Error allocating CRTM RTSolution structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3b. Allocate the Atmosphere structure
  ! --------------------------------------
  CALL CRTM_Atmosphere_Create( Atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(Atm)) ) THEN
    Message = 'Error allocating CRTM Atmosphere structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 3c. Allocate the Options structures
  ! -------------------------------------
  CALL CRTM_Options_Create( Options_Default, n_Channels )
  CALL CRTM_Options_Create( Options_AOP,     n_Channels )
  CALL CRTM_Options_Create( Options_COP,     n_Channels )
  CALL CRTM_Options_Create( Options_TOP,     n_Channels )
  IF ( ANY(.NOT. CRTM_Options_Associated(Options_Default)) .OR. &
       ANY(.NOT. CRTM_Options_Associated(Options_AOP))     .OR. &
       ANY(.NOT. CRTM_Options_Associated(Options_COP))     .OR. &
       ANY(.NOT. CRTM_Options_Associated(Options_TOP))     ) THEN
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

  ! Fix the "default" run's stream count at 16 (the maximum CRTM supports)
  ! instead of letting CRTM_Compute_nStreams pick 4/6 per channel, so it can
  ! be compared against the AOP/COP/TOP runs below, which were generated by
  ! Generate_OP.f90 using that same fixed 16-stream setting.
  CALL CRTM_Options_SetValue( Options_Default, &
                               n_Streams = 16 )


  ! 4c. Load the aerosol-only optical profile (AOP)
  ! -------------------------------------------------
  op_File = COEFFICIENTS_PATH//AOP_FILE
  PRINT *, 'Reading optical profile from file : ', op_File
  Error_Status = OP_Input_ReadFile(op_File, OP_AOP)
  IF ( Error_Status /= SUCCESS ) THEN
     Message = 'Error reading AOP file'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
  END IF
  Options_AOP(1)%AOP%n_Channels       = OP_AOP%n_Channels
  Options_AOP(1)%AOP%n_Layers         = OP_AOP%n_Layers
  Options_AOP(1)%AOP%n_Phase_Elements = OP_AOP%n_Phase_Elements
  Options_AOP(1)%AOP%n_Legendre_Terms = OP_AOP%n_Legendre_Terms
  Options_AOP(1)%AOP%Channel_Index    = OP_AOP%Channel_Index
  Options_AOP(1)%AOP%tau    = OP_AOP%tau
  Options_AOP(1)%AOP%bs     = OP_AOP%bs
  Options_AOP(1)%AOP%kb     = OP_AOP%kb
  Options_AOP(1)%AOP%pcoeff = OP_AOP%pcoeff
  ! n_Streams=16 matches Options_Default above and Generate_OP.f90's fixed
  ! stream count - all runs being compared must agree on this value.
  CALL CRTM_Options_SetValue( Options_AOP, &
                               n_Streams      = 16     , &
                               Use_Aerosol_OP = .TRUE.   )

  ! 4d. Load the cloud-only optical profile (COP)
  ! -------------------------------------------------
  op_File = COEFFICIENTS_PATH//COP_FILE
  PRINT *, 'Reading optical profile from file : ', op_File
  Error_Status = OP_Input_ReadFile(op_File, OP_COP)
  IF ( Error_Status /= SUCCESS ) THEN
     Message = 'Error reading COP file'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
  END IF
  Options_COP(1)%COP%n_Channels       = OP_COP%n_Channels
  Options_COP(1)%COP%n_Layers         = OP_COP%n_Layers
  Options_COP(1)%COP%n_Phase_Elements = OP_COP%n_Phase_Elements
  Options_COP(1)%COP%n_Legendre_Terms = OP_COP%n_Legendre_Terms
  Options_COP(1)%COP%Channel_Index = OP_COP%Channel_Index
  Options_COP(1)%COP%tau    = OP_COP%tau
  Options_COP(1)%COP%bs     = OP_COP%bs
  Options_COP(1)%COP%kb     = OP_COP%kb
  Options_COP(1)%COP%pcoeff = OP_COP%pcoeff
  CALL CRTM_Options_SetValue( Options_COP, &
                               n_Streams    = 16    , &
                               Use_Cloud_OP = .TRUE.  )

  ! 4e. Load the combined total optical profile (TOP)
  ! -------------------------------------------------
  op_File = COEFFICIENTS_PATH//TOP_FILE
  PRINT *, 'Reading optical profile from file : ', op_File
  Error_Status = OP_Input_ReadFile(op_File, OP_TOP)
  IF ( Error_Status /= SUCCESS ) THEN
     Message = 'Error reading TOP file'
     CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
     STOP 1
  END IF
  Options_TOP(1)%TOP%n_Channels       = OP_TOP%n_Channels
  Options_TOP(1)%TOP%n_Layers         = OP_TOP%n_Layers
  Options_TOP(1)%TOP%n_Phase_Elements = OP_TOP%n_Phase_Elements
  Options_TOP(1)%TOP%n_Legendre_Terms = OP_TOP%n_Legendre_Terms
  Options_TOP(1)%TOP%Channel_Index    = OP_TOP%Channel_Index
  Options_TOP(1)%TOP%tau    = OP_TOP%tau
  Options_TOP(1)%TOP%bs     = OP_TOP%bs
  Options_TOP(1)%TOP%kb     = OP_TOP%kb
  Options_TOP(1)%TOP%pcoeff = OP_TOP%pcoeff
  CALL CRTM_Options_SetValue( Options_TOP, &
                               n_Streams    = 16    , &
                               Use_Total_OP = .TRUE.  )
  ! ============================================================================


  ! ============================================================================
  ! 5. **** CALL THE CRTM FORWARD MODEL ****
  !
  ! 5a. Default CRTM interface: cloud/aerosol scattering computed internally,
  !     but with the stream count fixed at 16 via Options_Default so it can
  !     be compared against the AOP/COP/TOP runs below.
  ! ----------------------------------------------------------------------------
  Error_Status = CRTM_Forward( Atm            , &
                               Sfc            , &
                               Geometry       , &
                               ChannelInfo    , &
                               RTSolution     , &
                               Options_Default  )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  PRINT *, 'FINISH CRTM Default Interface'

  ! 5b. User-defined aerosol-only optical profile (AOP)
  !     Cloud scattering is still computed internally.
  ! ----------------------------------------------------------------------------
  Error_Status = CRTM_Forward( Atm           , &
                               Sfc           , &
                               Geometry      , &
                               ChannelInfo   , &
                               RTSolution_AOP, &
                               Options_AOP     )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model with Options%AOP'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  PRINT *, 'FINISH CRTM AOP Interface'

  ! 5c. User-defined cloud-only optical profile (COP)
  !     Aerosol scattering is still computed internally.
  ! ----------------------------------------------------------------------------
  Error_Status = CRTM_Forward( Atm           , &
                               Sfc           , &
                               Geometry      , &
                               ChannelInfo   , &
                               RTSolution_COP, &
                               Options_COP     )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model with Options%COP'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  PRINT *, 'FINISH CRTM COP Interface'

  ! 5d. User-defined combined total optical profile (TOP)
  !     Neither cloud nor aerosol scattering is computed internally.
  ! ----------------------------------------------------------------------------
  Error_Status = CRTM_Forward( Atm           , &
                               Sfc           , &
                               Geometry      , &
                               ChannelInfo   , &
                               RTSolution_TOP, &
                               Options_TOP     )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error in CRTM Forward Model with Options%TOP'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF
  PRINT *, 'FINISH CRTM TOP Interface'
  ! ============================================================================


  ! ============================================================================
  ! 6. **** COMPARE RTSolution RESULTS TO THE DEFAULT FORWARD CALCULATION ****
  !
  WRITE( *, '( /5x, "Comparing AOP/COP/TOP results with the default forward calculation..." )' )

  CALL Compare_And_Report( RTSolution, RTSolution_AOP, 'AOP' )
  CALL Compare_And_Report( RTSolution, RTSolution_COP, 'COP' )
  CALL Compare_And_Report( RTSolution, RTSolution_TOP, 'TOP' )
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
  ! 8. **** CLEAN UP ****
  !
  ! 8a. Deallocate the structures
  ! -----------------------------
  CALL CRTM_Atmosphere_Destroy(Atm)

  ! 8b. Deallocate the arrays
  ! -------------------------
  DEALLOCATE(RTSolution, RTSolution_AOP, RTSolution_COP, RTSolution_TOP, STAT=Allocate_Status)
  ! ============================================================================

  ! Signal the completion of the program. It is not a necessary step for running CRTM.

CONTAINS

  INCLUDE 'Load_Atm_Data_SingleProfile.inc'
  INCLUDE 'Load_Sfc_Data_SingleProfile.inc'

  ! Compare a test RTSolution array against the reference (default) one and
  ! report the outcome. On mismatch, the test RTSolution is written to a
  ! netCDF file (named after Label) for offline inspection.
  SUBROUTINE Compare_And_Report( rts_Ref, rts_Test, Label )
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: rts_Ref(:,:), rts_Test(:,:)
    CHARACTER(*),                INTENT(IN) :: Label
    CHARACTER(256) :: local_Message
    CHARACTER(256) :: local_File

    IF ( ALL(CRTM_RTSolution_Compare(rts_Ref, rts_Test)) ) THEN
      local_Message = Label//': RTSolution results are the same as the default calculation!'
      CALL Display_Message( PROGRAM_NAME, local_Message, INFORMATION )
    ELSE
      local_Message = Label//': RTSolution results are different from the default calculation!'
      CALL Display_Message( PROGRAM_NAME, local_Message, FAILURE )
      ! Write the mismatching RTSolution results to file
      local_File = TRIM(PROGRAM_NAME)//'_'//TRIM(Sensor_Id)//'.'//Label//'.RTSolution.nc'
      Error_Status = CRTM_RTSolution_WriteFile( local_File, rts_Test, NetCDF=.TRUE., Quiet=.TRUE. )
      IF ( Error_Status /= SUCCESS ) THEN
        local_Message = 'Error creating temporary RTSolution save file for failed '//Label//' comparison'
        CALL Display_Message( PROGRAM_NAME, local_Message, FAILURE )
      END IF
      STOP 1
    END IF
  END SUBROUTINE Compare_And_Report

END PROGRAM test_OP
