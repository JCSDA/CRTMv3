!
! Generate_OP
!
! Utility program to generate user-defined aerosol/cloud/total optical
! profile (AOP/COP/TOP) netCDF files for a single test profile, using the
! same atmosphere/surface data as test_OP.f90. The resulting TOP file is
! the one read back by test_OP.f90 via OP_Input_ReadFile.



PROGRAM Generate_OP

  ! ============================================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  ! ...Additional modules needed to compute AtmOptics directly (not
  !    re-exported by CRTM_Module)
  USE CRTM_AtmOptics_Define,    ONLY: CRTM_AtmOptics_type   , &
                                       CRTM_AtmOptics_Create , &
                                       CRTM_AtmOptics_Destroy, &
                                       CRTM_AtmOptics_Zero
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type, &
                                       CRTM_GeometryInfo_SetValue
  USE CRTM_GeometryInfo,        ONLY: CRTM_GeometryInfo_Compute
  USE CRTM_CloudScatter,        ONLY: CRTM_Compute_CloudScatter
  USE CRTM_AerosolScatter,      ONLY: CRTM_Compute_AerosolScatter
  USE CSvar_Define,             ONLY: CSvar_type, CSvar_Create, CSvar_Destroy
  USE ASvar_Define,             ONLY: ASvar_type, ASvar_Create, ASvar_Destroy
  USE CRTM_CloudCoeff,          ONLY: CloudC
  USE CRTM_RTSolution,          ONLY: CRTM_Compute_nStreams
  USE CRTM_Atmosphere,          ONLY: CRTM_Atmosphere_AddLayers
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ==========================================================================

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME      = 'Generate_OP'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  CHARACTER(*), PARAMETER :: AOP_FILE = 'AOP_SingleProfile.nc'
  CHARACTER(*), PARAMETER :: COP_FILE = 'COP_SingleProfile.nc'
  CHARACTER(*), PARAMETER :: TOP_FILE = 'TOP_SingleProfile.nc'

  ! ============================================================================
  ! 0. **** SOME SET UP PARAMETERS FOR THIS TEST ****
  !
  ! Profile dimensions... must match test_OP.f90
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
  INTEGER :: n_Channels
  INTEGER :: n_Phase_Elements
  INTEGER :: n_Legendre_Terms_OP
  INTEGER :: n_Full_Streams
  INTEGER :: SensorIndex, ChannelIndex
  INTEGER :: l
  TYPE(CRTM_RTSolution_type) :: RTS_Dummy

  ! ============================================================================
  ! 1. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  TYPE(CRTM_ChannelInfo_type)   :: ChannelInfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)      :: Geometry(N_PROFILES)
  TYPE(CRTM_Atmosphere_type)    :: Atm(N_PROFILES)
  ! CRTM_Forward_Module internally extends the user's profile above its top
  ! level with a climatological reference atmosphere (CRTM_Atmosphere_AddLayers)
  ! before computing cloud/aerosol scattering, so Atm_x (not Atm) is what must
  ! be used for the scattering calls below - otherwise Atm_x%n_Layers ends up
  ! LARGER than N_LAYERS, and CRTM_Forward_Module reads opt%AOP/COP/TOP out of
  ! bounds for the added layers when it later consumes these OP_Input files.
  TYPE(CRTM_Atmosphere_type)    :: Atm_x(N_PROFILES)
  TYPE(CRTM_Surface_type)       :: Sfc(N_PROFILES)
  TYPE(CRTM_GeometryInfo_type)  :: GeometryInfo(N_PROFILES)

  ! Per-channel scattering optics and their internal (Jacobian) variables
  TYPE(CRTM_AtmOptics_type) :: AtmOptics_Cloud
  TYPE(CRTM_AtmOptics_type) :: AtmOptics_Aerosol
  TYPE(CRTM_AtmOptics_type) :: AtmOptics_Total
  TYPE(CSvar_type) :: CSvar
  TYPE(ASvar_type) :: ASvar

  ! Output optical profile objects (aerosol-only, cloud-only, combined)
  TYPE(OP_Input_type) :: AOP
  TYPE(OP_Input_type) :: COP
  TYPE(OP_Input_type) :: TOP
  ! ============================================================================


  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Program to generate user-defined aerosol/cloud/total optical profile (AOP/COP/TOP) files.', &
    'CRTM Version: '//TRIM(Version) )

  ! Sensor_Id
  Sensor_Id = 'v.abi_gr'
  WRITE( *,'(//5x,"Running CRTM for ",a," sensor...")' ) TRIM(Sensor_Id)
  ! ============================================================================
  ! 2. **** INITIALIZE THE CRTM ****
  !
  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  Error_Status = CRTM_Init( (/Sensor_Id/), &
                            ChannelInfo, &
                            File_Path=COEFFICIENTS_PATH)
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! Determine the total number of channels
  n_Channels = SUM(CRTM_ChannelInfo_n_Channels(ChannelInfo))
  ! ============================================================================

  ! ============================================================================
  ! 3. **** ALLOCATE STRUCTURE ARRAYS ****
  !
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
  CALL Load_Atm_Data_SingleProfile()
  CALL Load_Sfc_Data_SingleProfile()

  ! ...Extend the profile above its top level, exactly as CRTM_Forward does
  ! internally, so the scattering computed below lines up with the number
  ! of layers CRTM_Forward_Module will actually use.
  Error_Status = CRTM_Atmosphere_AddLayers( Atm(1), Atm_x(1) )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error adding extra layers to the CRTM Atmosphere structure'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  ! 4b. GeometryInfo input
  ! ----------------------
  ! All profiles are given the same value. The Sensor_Scan_Angle is optional.
  CALL CRTM_Geometry_SetValue( Geometry, &
                               Sensor_Zenith_Angle = ZENITH_ANGLE, &
                               Sensor_Scan_Angle   = SCAN_ANGLE )
  CALL CRTM_GeometryInfo_SetValue( GeometryInfo, Geometry=Geometry )
  CALL CRTM_GeometryInfo_Compute( GeometryInfo )
  ! ============================================================================


  ! ============================================================================
  ! 5. **** COMPUTE THE AEROSOL/CLOUD/TOTAL OPTICAL PROFILES ****
  !
  ! Dimensions of the scattering coefficient tables. The same
  ! n_Phase_Elements is used for cloud, aerosol, and total optics so that
  ! their AtmOptics objects share the same array shapes and can be added
  ! together directly (this mirrors CRTM_Forward_Module, which accumulates
  ! cloud and aerosol contributions into a single shared AtmOptics object).
  n_Phase_Elements = CloudC%N_PHASE_ELEMENTS
  ! OP_Input%pcoeff stores the Legendre index 0:n_Full_Streams shifted by
  ! one (see Pack_AtmOptics below), so its Legendre dimension is one larger
  ! than the largest per-channel AtmOptics%n_Legendre_Terms can be.
  n_Legendre_Terms_OP = MAX_N_STREAMS + 1

  ! Allocate the OP_Input objects that will hold all channels' worth of
  ! aerosol-only (AOP), cloud-only (COP), and combined (TOP) optics. The
  ! layer dimension must match the EXTENDED atmosphere (Atm_x), not the
  ! user-supplied N_LAYERS, since that is what CRTM_Forward_Module will
  ! actually loop over when it reads these files back via Options.
  CALL OP_Input_Create( AOP, n_Channels, Atm_x(1)%n_Layers, n_Phase_Elements, n_Legendre_Terms_OP )
  CALL OP_Input_Create( COP, n_Channels, Atm_x(1)%n_Layers, n_Phase_Elements, n_Legendre_Terms_OP )
  CALL OP_Input_Create( TOP, n_Channels, Atm_x(1)%n_Layers, n_Phase_Elements, n_Legendre_Terms_OP )
  IF ( .NOT. ( OP_Input_Associated(AOP) .AND. OP_Input_Associated(COP) .AND. OP_Input_Associated(TOP) ) ) THEN
    Message = 'Error allocating OP_Input structures'
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  SensorIndex = 1
  Channel_Loop: DO l = 1, n_Channels
    ChannelIndex = l

    ! Number of RT streams to use for this channel's scattering computation.
    ! Fixed at 16 (the maximum CRTM supports) to match Options%n_Streams=16
    ! forced in test_OP.f90 for the default/AOP/COP/TOP calls alike - the
    ! stream count selects AtmOptics%n_Legendre_Terms and hence the lOffset
    ! slice of the Legendre-coefficient lookup table (see CRTM_Compute_CloudScatter/
    ! AerosolScatter), so it must be identical on the generating side (here)
    ! and the consuming side (test_OP.f90) or the phase coefficients come
    ! from the wrong table slice.
    !
    ! CRTM_Compute_nStreams would instead pick 4 or 6 per channel here (the
    ! same choice CRTM_Forward_Module makes when Options%Use_n_Streams is
    ! left off) - kept below, commented out, in case a dynamic-stream-count
    ! version of this program is needed again.
    ! n_Full_Streams = CRTM_Compute_nStreams( Atm_x(1), SensorIndex, ChannelIndex, RTS_Dummy )
    n_Full_Streams = MAX_N_STREAMS

    ! 2. Cloud optical properties (COP)
    ! ----------------------------------
    ! CRTM_AtmOptics_Create's internal Zero-before-allocation-flag-is-set is a
    ! no-op the first time an object is (re)created after CRTM_AtmOptics_Destroy
    ! (which never deallocates, only flips Is_Allocated to .FALSE.), so without
    ! this explicit Zero, Phase_Coefficient etc. would silently retain values
    ! left over from the previous channel's computation.
    CALL CRTM_AtmOptics_Create( AtmOptics_Cloud, Atm_x(1)%n_Layers, n_Full_Streams, n_Phase_Elements )
    CALL CRTM_AtmOptics_Zero( AtmOptics_Cloud )
    AtmOptics_Cloud%Include_Scattering = .TRUE.
    IF ( Atm_x(1)%n_Clouds > 0 ) THEN
      CALL CSvar_Create( CSvar, n_Full_Streams, n_Phase_Elements, Atm_x(1)%n_Layers, Atm_x(1)%n_Clouds )
      Error_Status = CRTM_Compute_CloudScatter(  Atm_x(1)       , &  ! Input
                                                 GeometryInfo(1), &  ! Input
                                                 SensorIndex    , &  ! Input
                                                 ChannelIndex   , &  ! Input
                                                 AtmOptics_Cloud, &  ! Output
                                                 CSvar            )  ! Internal variable output
      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message,'("Error computing CloudScatter for channel ",i0)' ) l
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP 1
      END IF
      CALL CSvar_Destroy( CSvar )
    END IF

    ! 1. Aerosol optical properties (AOP)
    ! ------------------------------------
    CALL CRTM_AtmOptics_Create( AtmOptics_Aerosol, Atm_x(1)%n_Layers, n_Full_Streams, n_Phase_Elements )
    CALL CRTM_AtmOptics_Zero( AtmOptics_Aerosol )
    AtmOptics_Aerosol%Include_Scattering = .TRUE.
    IF ( Atm_x(1)%n_Aerosols > 0 ) THEN
      CALL ASvar_Create( ASvar, n_Full_Streams, n_Phase_Elements, Atm_x(1)%n_Layers, Atm_x(1)%n_Aerosols )
      Error_Status = CRTM_Compute_AerosolScatter(  Atm_x(1)         , &  ! Input
                                                   SensorIndex      , &  ! Input
                                                   ChannelIndex     , &  ! Input
                                                   AtmOptics_Aerosol, &  ! Output
                                                   ASvar              )  ! Internal variable output
      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message,'("Error computing AerosolScatter for channel ",i0)' ) l
        CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
        STOP 1
      END IF
      CALL ASvar_Destroy( ASvar )
    END IF

    ! 3. Combine AOP + COP into TOP
    ! ------------------------------
    CALL CRTM_AtmOptics_Create( AtmOptics_Total, Atm_x(1)%n_Layers, n_Full_Streams, n_Phase_Elements )
    CALL CRTM_AtmOptics_Zero( AtmOptics_Total )
    AtmOptics_Total%Include_Scattering    = .TRUE.
    AtmOptics_Total%Optical_Depth         = AtmOptics_Cloud%Optical_Depth         + AtmOptics_Aerosol%Optical_Depth
    AtmOptics_Total%Single_Scatter_Albedo = AtmOptics_Cloud%Single_Scatter_Albedo + AtmOptics_Aerosol%Single_Scatter_Albedo
    AtmOptics_Total%Backscat_Coefficient  = AtmOptics_Cloud%Backscat_Coefficient  + AtmOptics_Aerosol%Backscat_Coefficient
    AtmOptics_Total%Phase_Coefficient     = AtmOptics_Cloud%Phase_Coefficient     + AtmOptics_Aerosol%Phase_Coefficient

    ! Pack this channel's optics into the three OP_Input objects
    CALL Pack_AtmOptics( AtmOptics_Cloud,   ChannelIndex, COP )
    CALL Pack_AtmOptics( AtmOptics_Aerosol, ChannelIndex, AOP )
    CALL Pack_AtmOptics( AtmOptics_Total,   ChannelIndex, TOP )

    CALL CRTM_AtmOptics_Destroy( AtmOptics_Cloud )
    CALL CRTM_AtmOptics_Destroy( AtmOptics_Aerosol )
    CALL CRTM_AtmOptics_Destroy( AtmOptics_Total )

  END DO Channel_Loop
  ! ============================================================================


  ! ============================================================================
  ! 6. **** WRITE THE AOP/COP/TOP DATA TO SEPARATE NETCDF FILES ****
  !
  WRITE( *, '( /5x, "Writing optical profile files..." )' )

  Error_Status = OP_Input_WriteFile( AOP_FILE, AOP )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error writing AOP file '//AOP_FILE
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Error_Status = OP_Input_WriteFile( COP_FILE, COP )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error writing COP file '//COP_FILE
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
    STOP 1
  END IF

  Error_Status = OP_Input_WriteFile( TOP_FILE, TOP )
  IF ( Error_Status /= SUCCESS ) THEN
    Message = 'Error writing TOP file '//TOP_FILE
    CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
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
  ! 8. **** CLEAN UP ****
  !
  CALL CRTM_Atmosphere_Destroy(Atm)
  CALL CRTM_Atmosphere_Destroy(Atm_x)
  ! ============================================================================

CONTAINS

  INCLUDE 'Load_Atm_Data_SingleProfile.inc'
  INCLUDE 'Load_Sfc_Data_SingleProfile.inc'

  ! Copy one channel's computed AtmOptics into the corresponding channel
  ! slice of an OP_Input object. The Legendre index in Phase_Coefficient is
  ! 0-based (0:n_Legendre_Terms) while OP_Input%pcoeff's Legendre dimension
  ! is 1-based, hence the ileg+1 shift - this matches the convention used
  ! when CRTM_Forward_Module reads Options%TOP/COP/AOP back into AtmOptics.
  SUBROUTINE Pack_AtmOptics( AO, ChannelIndex, OP )
    TYPE(CRTM_AtmOptics_type), INTENT(IN)    :: AO
    INTEGER,                   INTENT(IN)    :: ChannelIndex
    TYPE(OP_Input_type),       INTENT(INOUT) :: OP
    INTEGER :: ilay, iphas, ileg, ileg1

    ! Record which SpcCoeff ChannelIndex this column corresponds to, so
    ! CRTM_Forward_Module's OP_Input_Channel_Position lookup can find it
    OP%Channel_Index(ChannelIndex) = ChannelIndex

    DO ilay = 1, AO%n_Layers
      OP%tau(ChannelIndex,ilay) = AO%Optical_Depth(ilay)
      OP%bs(ChannelIndex,ilay)  = AO%Single_Scatter_Albedo(ilay)
      OP%kb(ChannelIndex,ilay)  = AO%Backscat_Coefficient(ilay)
      DO iphas = 1, AO%n_Phase_Elements
        DO ileg = 0, AO%n_Legendre_Terms
          ileg1 = ileg + 1
          OP%pcoeff(ChannelIndex,ilay,iphas,ileg1) = AO%Phase_Coefficient(ileg,iphas,ilay)
        END DO
      END DO
    END DO
  END SUBROUTINE Pack_AtmOptics

END PROGRAM Generate_OP
