!
! check_crtm
!
! Check/example program for the CRTM Forward and K-Matrix functions.
!
! NOTE: No results are output or compared in this program.
!       At this stage, it is included as an example only,
!       and to determine that the CRTM library built
!       correctly such that it can be linked to create
!       an executable.
!

! Aerosol-Index mapping:
!  --- CRTM ---
!  1 = Dust
!  2 = Sea salt-SSAM
!  3 = Sea salt-SSCM1
!  4 = Sea salt-SSCM2
!  5 = Sea salt-SSCM3
!  6 = Organic carbon
!  7 = Black carbon
!  8 = Sulfate
!  --- CMAQ ---
!  1 = Dust
!  2 = Soot
!  3 = Water soluble
!  4 = Sulfate
!  5 = Sea salt
!  6 = Water
!  7 = Insoluble
!  8 = dust-like
!  --- GOCART-GEOS5 ---
!  1, 2, 3, 4, 5  = Dust 1, 2, 3, 4, 5
!  6, 7, 8, 9, 10 = Sea salt 1, 2, 3, 4, 5
!  11, 12 = Organic carbon 1, 2
!  13, 14 = Black carbon 1, 2
!  15, 16 = Sulfate 1, 2
!  17, 18, 19 = Nitrate 1, 2, 3
!  --- NAAPS ---
!  1 = Dust
!  2 = Smoke
!  3 = Sea Salt
!  4 = Anthropogenic and Biogenic Fine Particles

PROGRAM check_crtm

  ! ============================================================================
  ! STEP 1. **** ENVIRONMENT SETUP FOR CRTM USAGE ****
  !
  ! Module usage
  USE CRTM_Module
  USE CRTM_RTSolution_Define, ONLY: CRTM_RTSolution_WriteFile_netCDF
  USE netcdf
  ! Disable all implicit typing
  IMPLICIT NONE
  ! ============================================================================



  ! --------------------------
  ! Some non-CRTM-y Parameters
  ! --------------------------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'check_crtm'

  ! ============================================================================
  ! STEP 2. **** SET UP SOME PARAMETERS FOR THE CRTM RUN ****
  !
  ! Directory location of coefficients
#ifdef LITTLE_ENDIAN
  CHARACTER(*), PARAMETER :: ENDIAN_TYPE='little_endian'
#else
  CHARACTER(*), PARAMETER :: ENDIAN_TYPE='big_endian'
#endif
  CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='coefficients/'//ENDIAN_TYPE//'/'
  CHARACTER(*), PARAMETER :: NC_COEFFICIENT_PATH='coefficients/netcdf/'

  ! Aerosol/Cloud coefficient format
  CHARACTER(*), PARAMETER :: Coeff_Format = 'Binary'
  !CHARACTER(*), PARAMETER :: Coeff_Format = 'netCDF'

  ! Aerosol/Cloud coefficient scheme
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'CRTM'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'CMAQ'
  CHARACTER(*), PARAMETER :: Aerosol_Model = 'GOCART-GEOS5'
  !CHARACTER(*), PARAMETER :: Aerosol_Model = 'NAAPS'
  CHARACTER(*), PARAMETER :: Cloud_Model   = 'CRTM'

  ! Directory location of results for comparison [NOT USED YET]
  CHARACTER(*), PARAMETER :: RESULTS_PATH = './results/'

  ! Profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_LAYERS    = 30
  INTEGER, PARAMETER :: N_ABSORBERS = 2
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 20  !tbd
  INTEGER, PARAMETER :: N_TIME = 1
  INTEGER, PARAMETER :: N_LAT = 1
  INTEGER, PARAMETER :: N_LON = 1
  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'v.modis_aqua'/)

  ! Some pretend geometry angles. The scan angle is based
  ! on the default Re (earth radius) and h (satellite height)
  REAL(fp), PARAMETER :: ZENITH_ANGLE = 30.0_fp
  REAL(fp), PARAMETER :: SCAN_ANGLE   = 26.37293341421_fp

  ! WRF-Chem input files
  INTEGER :: FileId
  INTEGER :: NF90_Status
  INTEGER :: VarId
  CHARACTER(*), PARAMETER :: Input_File = '/Users/dangch/Documents/WRF/WRF_output/WRFChem_map2_CRTM.nc4'
  REAL(Double), DIMENSION(N_LAYERS) :: rh, Reff_i, Reff_j, Reff_k, &
                                       PresLayer, Temperature
  REAL(Double), DIMENSION(0:N_LAYERS) :: PresLevel
  REAL(Double), DIMENSION(N_LON, N_LAT, N_LAYERS, N_TIME) :: DUk, SSi, SSj, SSk, &
                                       OCiHP, OCjHP, OCkHP, OCiHB, OCjHB, OCkHB, &
                                       ECiHP, ECjHP, ECkHP, ECiHB, ECjHB, ECkHB, &
                                       SFi, SFj, NIi, NIj

  ! ============================================================================

  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: message, version
  CHARACTER(256) :: AerosolCoeff_File
  CHARACTER(256) :: AerosolCoeff_Format
  CHARACTER(256) :: CloudCoeff_File
  CHARACTER(256) :: CloudCoeff_Format
  CHARACTER(256) :: Aerosol_Scheme
  CHARACTER(256) :: Cloud_Scheme
  CHARACTER(256) :: output_nc_file
  INTEGER :: err_stat, alloc_stat
  INTEGER :: n_channels
  INTEGER :: l, m, n, nc
  ! ============================================================================
  ! STEP 3. **** DEFINE THE CRTM INTERFACE STRUCTURES ****
  !
  ! 3a. Define the "non-demoninational" arguments
  ! ---------------------------------------------
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)
  TYPE(CRTM_Geometry_type)                :: geo(N_PROFILES)


  ! 3b. Define the FORWARD variables
  ! --------------------------------
  TYPE(CRTM_Atmosphere_type)              :: atm(N_PROFILES)
  TYPE(CRTM_Surface_type)                 :: sfc(N_PROFILES)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)

  ! 3c. Define the K-MATRIX variables
  ! ---------------------------------
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: atm_K(:,:)
  TYPE(CRTM_Surface_type)   , ALLOCATABLE :: sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts_K(:,:)
  ! ============================================================================


  ! Program header
  ! --------------
  CALL CRTM_Version( Version )
  CALL Program_Message( PROGRAM_NAME, &
    'Check/example program for the CRTM Forward and K-Matrix functions using '//&
    ENDIAN_TYPE//' coefficient datafiles', &
    'CRTM Version: '//TRIM(Version) )



  ! ============================================================================
  ! STEP 4. **** INITIALIZE THE CRTM ****
  !
  ! 4a. Initialise all the sensors at once
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
  err_stat = CRTM_Init( SENSOR_ID, &
                        chinfo, &
                        Aerosol_Model, &
                        AerosolCoeff_Format, &
                        AerosolCoeff_File, &
                        Cloud_Model, &
                        CloudCoeff_Format, &
                        CloudCoeff_File, &
                        File_Path=COEFFICIENT_PATH, &
                        NC_File_Path=NC_COEFFICIENT_PATH, &
                        Quiet=.TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  ! 4b. Output some channel information
  ! -----------------------------------
  n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  DO n = 1, N_SENSORS
    WRITE( *,'(7x,i0," from ",a)' ) &
      CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(SENSOR_ID(n))
  END DO
  ! ============================================================================



  ! Begin loop over sensors
  ! ----------------------
  Sensor_Loop: DO n = 1, N_SENSORS

    ! ==========================================================================
    ! STEP 5. **** ALLOCATE STRUCTURE ARRAYS ****
    !
    ! 5a. Determine the number of channels
    !     for the current sensor
    ! ------------------------------------
    n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

    ! 5b. Allocate the ARRAYS
    ! -----------------------
    ALLOCATE( rts( n_channels, N_PROFILES ), &
              atm_K( n_channels, N_PROFILES ), &
              sfc_K( n_channels, N_PROFILES ), &
              rts_K( n_channels, N_PROFILES ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) THEN
      message = 'Error allocating structure arrays'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF


    ! 5c. Allocate the STRUCTURE INTERNALS
    !     NOTE: Only the Atmosphere structures
    !           are allocated in this example
    ! ----------------------------------------
    ! The input FORWARD structure
    CALL CRTM_Atmosphere_Create( atm, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm)) ) THEN
      message = 'Error allocating CRTM Forward Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! The output K-MATRIX structure
    CALL CRTM_Atmosphere_Create( atm_K, N_LAYERS, N_ABSORBERS, N_CLOUDS, N_AEROSOLS )
    IF ( ANY(.NOT. CRTM_Atmosphere_Associated(atm_K)) ) THEN
      message = 'Error allocating CRTM K-matrix Atmosphere structure'
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! The output structure
    CALL CRTM_RTSolution_Create( rts, N_LAYERS )
    IF ( ANY(.NOT. CRTM_RTSolution_Associated(rts)) ) THEN
      Message = 'Error allocating CRTM RTSolution structures'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! ==========================================================================

    ! ==========================================================================
    ! STEP 6. **** ASSIGN INPUT DATA ****
    !
    ! 6a. Atmosphere and Surface input
    !     NOTE: that this is the hard part (in my opinion :o). The mechanism by
    !     by which the atmosphere and surface data are loaded in to their
    !     respective structures below was done purely to keep the step-by-step
    !     instructions in this program relatively "clean".
    ! ------------------------------------------------------------------------
    PRINT *, 'Start loading ATM and SFC data'
    CALL Load_Atm_Data()
    CALL Load_Sfc_Data()
    PRINT *, 'ATM and SFC data loaded'

    ! 6b. Geometry input
    ! ------------------
    ! All profiles are given the same value
    !  The Sensor_Scan_Angle is optional.
    CALL CRTM_Geometry_SetValue( geo, &
                                 Sensor_Zenith_Angle = ZENITH_ANGLE, &
                                 Sensor_Scan_Angle   = SCAN_ANGLE )
    ! ==========================================================================




    ! ==========================================================================
    ! STEP 7. **** INITIALIZE THE K-MATRIX ARGUMENTS ****
    !
    ! 7a. Zero the K-matrix OUTPUT structures
    ! ---------------------------------------
    CALL CRTM_Atmosphere_Zero( atm_K )
    CALL CRTM_Surface_Zero( sfc_K )


    ! 7b. Inintialize the K-matrix INPUT so
    !     that the results are dTb/dx
    ! -------------------------------------
    rts_K%Radiance               = ZERO
    rts_K%Brightness_Temperature = ONE
    ! ==========================================================================

    ! ==========================================================================
    ! STEP 8. **** CALL THE CRTM FUNCTIONS FOR THE CURRENT SENSOR ****
    !
    WRITE( *, '( /5x, "Calling the CRTM functions for ",a,"..." )' ) TRIM(SENSOR_ID(n))

    ! 8a. The forward model
    ! ---------------------
    err_stat = CRTM_Forward( atm        , &  ! Input
                             sfc        , &  ! Input
                             geo        , &  ! Input
                             chinfo(n:n), &  ! Input
                             rts          )  ! Output
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM Forward Model for '//TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF


    ! 8b. The K-matrix model
    ! ----------------------
    err_stat = CRTM_K_Matrix( atm        , &  ! FORWARD  Input
                              sfc        , &  ! FORWARD  Input
                              rts_K      , &  ! K-MATRIX Input
                              geo        , &  ! Input
                              chinfo(n:n), &  ! Input
                              atm_K      , &  ! K-MATRIX Output
                              sfc_K      , &  ! K-MATRIX Output
                              rts          )  ! FORWARD  Output
    IF ( err_stat /= SUCCESS ) THEN
      message = 'Error calling CRTM K-Matrix Model for '//TRIM(SENSOR_ID(n))
      CALL Display_Message( PROGRAM_NAME, message, FAILURE )
      STOP
    END IF

    ! 8b. The AOD model
    ! ----------------------
    err_stat = CRTM_AOD( atm        , &
                         chinfo(n:n), &
                         rts )
    IF ( err_stat /= SUCCESS ) THEN
      Message = 'Error in CRTM AOD Model'
      CALL Display_Message( PROGRAM_NAME, Message, FAILURE )
      STOP 1
    END IF
    ! ==========================================================================

   ! ============================================================================
   ! 8c. **** OUTPUT THE RESULTS TO SCREEN ****
   !
   ! User should read the user guide or the source code of the routine
   ! CRTM_RTSolution_Inspect in the file CRTM_RTSolution_Define.f90 to
   ! select the needed variables for outputs.  These variables are contained
   ! in the structure RTSolution.
   ! DO m = 1, N_PROFILES
   !   WRITE( *,'(//7x,"Profile ",i0," output for ",a )') m, TRIM(Sensor_Id(n))
   !   DO l = 1, n_Channels
   !     WRITE( *, '(/5x,"Channel ",i0," results")') chinfo(n)%Sensor_Channel(l)
   !     CALL CRTM_RTSolution_Inspect(rts(l,m))
   !     CALL CRTM_RTSolution_Inspect(rts_K(l,m))
   !     CALL CRTM_Atmosphere_Inspect(atm_K(l,m))
   !     CALL CRTM_Surface_Inspect(sfc_K(l,m))
   !
   !   END DO
   !   CALL CRTM_Atmosphere_Inspect(atm(m))
   !   CALL CRTM_Surface_Inspect(sfc(m))
   ! END DO

    ! 8d. **** OUTPUT THE RESULTS TO NetCDF files ****
    output_nc_file = TRIM(SENSOR_ID(n))
    err_stat = CRTM_RTSolution_WriteFile_netCDF( &
                       output_nc_file     , &  ! Input
                       rts                  )  ! Input
     IF ( err_stat /= SUCCESS ) THEN
       message = 'Error calling CRTM_RTSolution_WriteFile_netCDF for '//TRIM(SENSOR_ID(n))
       CALL Display_Message( PROGRAM_NAME, message, FAILURE )
       STOP
     END IF

    ! ==========================================================================
    ! STEP 9. **** CLEAN UP FOR NEXT SENSOR ****
    !
    ! 9a. Deallocate the structures
    ! -----------------------------
    CALL CRTM_Atmosphere_Destroy(atm_K)
    CALL CRTM_Atmosphere_Destroy(atm)


    ! 9b. Deallocate the arrays
    ! -------------------------
    DEALLOCATE(rts, rts_K, sfc_k, atm_k, STAT = alloc_stat)
    ! ==========================================================================

  END DO Sensor_Loop


  ! ==========================================================================
  ! 10. **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF
  ! ==========================================================================




  ! ==========================================================================
  ! 11. **** CREATE A SIGNAL FILE FOR TESTING SUCCESS ****
  !
  ! This step is just to allow the CRTM library build process
  ! to detect success or failure at the shell level
  CALL SignalFile_Create()
  ! ==========================================================================


CONTAINS


  ! ==========================================================================
  !                Below are some internal procedures that load the
  !                necessary input structures with some pretend data
  ! ==========================================================================

  !
  ! Internal subprogam to load some test profile data
  !
  SUBROUTINE Load_Atm_Data()
    ! Local variables
    INTEGER :: nc
    INTEGER :: k1, k2

    ! Read WRF-Chem output
    ! ... Read RH
    print *, 'Check input file: ', Input_File
    IF ( .NOT. File_Exists(Input_File) ) THEN
       WRITE(*,*) 'File_Exists(Input_File)'
    END IF
    NF90_Status = NF90_OPEN( Input_File, NF90_NOWRITE, FileId )
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: read input file'

    NF90_Status = NF90_INQ_VARID( FileId, 'RH', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, Relative_Humidity'
    NF90_Status = NF90_GET_VAR( FileId, VarId, RH)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, Relative_Humidity'

    ! ... Reff_i
    NF90_Status = NF90_INQ_VARID( FileId, 'Reff_i', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, Reff_i'
    NF90_Status = NF90_GET_VAR( FileId, VarId, Reff_i)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, Reff_i'

    ! ... Reff_j
    NF90_Status = NF90_INQ_VARID( FileId, 'Reff_j', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, Reff_j'
    NF90_Status = NF90_GET_VAR( FileId, VarId, Reff_j)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, Reff_j'

    ! ... Reff_k
    NF90_Status = NF90_INQ_VARID( FileId, 'Reff_k', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, Reff_k'
    NF90_Status = NF90_GET_VAR( FileId, VarId, Reff_k)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, Reff_k'

    ! ... PressureLev
    NF90_Status = NF90_INQ_VARID( FileId, 'PresLevel', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, PresLevel'
    NF90_Status = NF90_GET_VAR( FileId, VarId, PresLevel)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, PresLevel'

    ! ... Pressure layer
    NF90_Status = NF90_INQ_VARID( FileId, 'PresLayer', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, PresLayer'
    NF90_Status = NF90_GET_VAR( FileId, VarId, PresLayer)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, PresLayer'

    ! ... Temperature
    NF90_Status = NF90_INQ_VARID( FileId, 'Temperature', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, Temperature'
    NF90_Status = NF90_GET_VAR( FileId, VarId, Temperature)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, Temperature'

    ! ... DUk
    NF90_Status = NF90_INQ_VARID( FileId, 'DUk', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, DUk'
    NF90_Status = NF90_GET_VAR( FileId, VarId, DUk)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, DUk'

    ! ... SSi
    NF90_Status = NF90_INQ_VARID( FileId, 'SSi', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, SSi'
    NF90_Status = NF90_GET_VAR( FileId, VarId, SSi)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, SSi'

    ! ... SSj
    NF90_Status = NF90_INQ_VARID( FileId, 'SSj', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, SSj'
    NF90_Status = NF90_GET_VAR( FileId, VarId, SSj)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, SSj'

    ! ... SSk
    NF90_Status = NF90_INQ_VARID( FileId, 'SSk', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, SSk'
    NF90_Status = NF90_GET_VAR( FileId, VarId, SSk)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, SSk'

    ! ... OCiHP
    NF90_Status = NF90_INQ_VARID( FileId, 'OCiHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCiHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCiHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCiHP'

    ! ... OCjHP
    NF90_Status = NF90_INQ_VARID( FileId, 'OCjHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCjHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCjHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCjHP'

    ! ... OCkHP
    NF90_Status = NF90_INQ_VARID( FileId, 'OCkHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCkHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCkHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCkHP'

    ! ... OCiHB
    NF90_Status = NF90_INQ_VARID( FileId, 'OCiHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCiHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCiHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCiHB'

    ! ... OCjHB
    NF90_Status = NF90_INQ_VARID( FileId, 'OCjHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCjHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCjHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCjHB'

    ! ... OCkHB
    NF90_Status = NF90_INQ_VARID( FileId, 'OCkHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, OCkHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, OCkHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, OCkHB'

    ! ... ECiHP
    NF90_Status = NF90_INQ_VARID( FileId, 'ECiHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECiHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECiHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECiHP'

    ! ... ECjHP
    NF90_Status = NF90_INQ_VARID( FileId, 'ECjHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECjHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECjHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECjHP'

    ! ... ECkHP
    NF90_Status = NF90_INQ_VARID( FileId, 'ECkHP', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECkHP'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECkHP)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECkHP'

    ! ... ECiHB
    NF90_Status = NF90_INQ_VARID( FileId, 'ECiHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECiHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECiHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECiHB'

    ! ... ECjHB
    NF90_Status = NF90_INQ_VARID( FileId, 'ECjHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECjHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECjHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECjHB'

    ! ... ECkHB
    NF90_Status = NF90_INQ_VARID( FileId, 'ECkHB', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, ECkHB'
    NF90_Status = NF90_GET_VAR( FileId, VarId, ECkHB)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, ECkHB'

    ! ... SFi
    NF90_Status = NF90_INQ_VARID( FileId, 'SFi', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, SFi'
    NF90_Status = NF90_GET_VAR( FileId, VarId, SFi)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, SFi'

    ! ... SFj
    NF90_Status = NF90_INQ_VARID( FileId, 'SFj', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, SFj'
    NF90_Status = NF90_GET_VAR( FileId, VarId, SFj)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, SFj'

    ! ... NIi
    NF90_Status = NF90_INQ_VARID( FileId, 'NIi', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, NIi'
    NF90_Status = NF90_GET_VAR( FileId, VarId, NIi)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, NIi'

    ! ... NIj
    NF90_Status = NF90_INQ_VARID( FileId, 'NIj', VarId )
    IF ( NF90_Status /= NF90_NOERR) WRITE(*,*) 'ERROR: NF90_INQ_VARID, NIj'
    NF90_Status = NF90_GET_VAR( FileId, VarId, NIj)
    IF ( NF90_Status /= NF90_NOERR ) WRITE(*,*) 'ERROR: NF90_GET_VAR, NIj'

    Print *, 'Finish reading wrf-chem data'

    ! 4a.1 Profile #1
    ! ---------------
    ! ...Profile and absorber definitions
    atm(1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID                 , O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    ! ...Profile data
    atm(1)%Level_Pressure = PresLevel
    atm(1)%Pressure = PresLayer
    atm(1)%Temperature = Temperature

    ! print *, atm(1)%Level_Pressure


    ! atm(1)%Level_Pressure = &
    ! (/7.14000000e-01, 3.51721379e+01, 6.96302759e+01, 1.04088414e+02, &
    !    1.38546552e+02, 1.73004690e+02, 2.07462828e+02, 2.41920966e+02, &
    !    2.76379103e+02, 3.10837241e+02, 3.45295379e+02, 3.79753517e+02, &
    !    4.14211655e+02, 4.48669793e+02, 4.83127931e+02, 5.17586069e+02, &
    !    5.52044207e+02, 5.86502345e+02, 6.20960483e+02, 6.55418621e+02, &
    !    6.89876759e+02, 7.24334897e+02, 7.58793034e+02, 7.93251172e+02, &
    !    8.27709310e+02, 8.62167448e+02, 8.96625586e+02, 9.31083724e+02, &
    !    9.65541862e+02, 1.00000000e+03/)

    atm(1)%Relative_Humidity = RH  ! convert to fraction
    atm(1)%Absorber(:,1) = &
    (/4.187E-03_fp,4.401E-03_fp,4.250E-03_fp,3.688E-03_fp,3.516E-03_fp,3.739E-03_fp,3.694E-03_fp,3.449E-03_fp, &
      3.228E-03_fp,3.212E-03_fp,3.245E-03_fp,3.067E-03_fp,2.886E-03_fp,2.796E-03_fp,2.704E-03_fp,2.617E-03_fp, &
      2.568E-03_fp,2.536E-03_fp,2.506E-03_fp,2.468E-03_fp,2.427E-03_fp,2.438E-03_fp,2.493E-03_fp,2.543E-03_fp, &
      2.586E-03_fp,2.632E-03_fp,2.681E-03_fp,2.703E-03_fp,2.636E-03_fp,2.512E-03_fp/)

    atm(1)%Absorber(:,2) = &
    (/3.035E+00_fp,3.943E+00_fp,4.889E+00_fp,5.812E+00_fp,6.654E+00_fp,7.308E+00_fp,7.660E+00_fp,7.745E+00_fp, &
      7.696E+00_fp,7.573E+00_fp,7.413E+00_fp,7.246E+00_fp,7.097E+00_fp,6.959E+00_fp,6.797E+00_fp,6.593E+00_fp, &
      6.359E+00_fp,6.110E+00_fp,5.860E+00_fp,5.573E+00_fp,5.253E+00_fp,4.937E+00_fp,4.625E+00_fp,4.308E+00_fp, &
      3.986E+00_fp,3.642E+00_fp,3.261E+00_fp,2.874E+00_fp,2.486E+00_fp,2.102E+00_fp/)


    ! Load CO2 absorber data if there are three absorrbers
    ! IF ( atm(1)%n_Absorbers > 2 ) THEN
    !   atm(1)%Absorber_Id(3)    = CO2_ID
    !   atm(1)%Absorber_Units(3) = VOLUME_MIXING_RATIO_UNITS
    !   atm(1)%Absorber(:,3)     = 380.0_fp
    ! END IF


    ! Cloud data
    ! IF ( atm(1)%n_Clouds > 0 ) THEN
    !   k1 = 75
    !   k2 = 79
    !   DO nc = 1, atm(1)%n_Clouds
    !     atm(1)%Cloud(nc)%Type = WATER_CLOUD
    !     atm(1)%Cloud(nc)%Effective_Radius(k1:k2) = 20.0_fp ! microns
    !     atm(1)%Cloud(nc)%Water_Content(k1:k2)    = 5.0_fp  ! kg/m^2
    !   END DO
    ! END IF

    !Assign aerosols
    atm(1)%Aerosol(1)%Type = 3 !or 4, test both
    atm(1)%Aerosol(1)%Effective_Radius = Reff_k
    atm(1)%Aerosol(1)%Concentration = DUk(1,1,:,1)

    ! sea salt
    atm(1)%Aerosol(2)%Type = 6
    atm(1)%Aerosol(2)%Effective_Radius = Reff_i
    atm(1)%Aerosol(2)%Concentration = SSi(1,1,:,1)
    atm(1)%Aerosol(3)%Type = 7
    atm(1)%Aerosol(3)%Effective_Radius = Reff_j
    atm(1)%Aerosol(3)%Concentration = SSj(1,1,:,1)
    atm(1)%Aerosol(4)%Type = 9
    atm(1)%Aerosol(4)%Effective_Radius = Reff_k
    atm(1)%Aerosol(4)%Concentration = SSk(1,1,:,1)

    ! Organic carbon
    atm(1)%Aerosol(5)%Type = 12
    atm(1)%Aerosol(5)%Effective_Radius = Reff_i
    atm(1)%Aerosol(5)%Concentration = OCiHP(1,1,:,1)
    atm(1)%Aerosol(6)%Type = 12
    atm(1)%Aerosol(6)%Effective_Radius = Reff_j
    atm(1)%Aerosol(6)%Concentration = OCjHP(1,1,:,1)
    atm(1)%Aerosol(7)%Type = 12
    atm(1)%Aerosol(7)%Effective_Radius = Reff_k
    atm(1)%Aerosol(7)%Concentration = OCkHP(1,1,:,1)
    atm(1)%Aerosol(8)%Type = 11
    atm(1)%Aerosol(8)%Effective_Radius = Reff_i
    atm(1)%Aerosol(8)%Concentration = OCiHB(1,1,:,1)
    atm(1)%Aerosol(9)%Type = 11
    atm(1)%Aerosol(9)%Effective_Radius = Reff_j
    atm(1)%Aerosol(9)%Concentration = OCjHB(1,1,:,1)
    atm(1)%Aerosol(10)%Type = 11
    atm(1)%Aerosol(10)%Effective_Radius = Reff_k
    atm(1)%Aerosol(10)%Concentration = OCkHB(1,1,:,1)

    ! Black carbon
    atm(1)%Aerosol(11)%Type = 14
    atm(1)%Aerosol(11)%Effective_Radius = Reff_i
    atm(1)%Aerosol(11)%Concentration = ECiHP(1,1,:,1)
    atm(1)%Aerosol(12)%Type = 14
    atm(1)%Aerosol(12)%Effective_Radius = Reff_j
    atm(1)%Aerosol(12)%Concentration = ECjHP(1,1,:,1)
    atm(1)%Aerosol(13)%Type = 14
    atm(1)%Aerosol(13)%Effective_Radius = Reff_k
    atm(1)%Aerosol(13)%Concentration = ECkHP(1,1,:,1)
    atm(1)%Aerosol(14)%Type = 13
    atm(1)%Aerosol(14)%Effective_Radius = Reff_i
    atm(1)%Aerosol(14)%Concentration = ECiHB(1,1,:,1)
    atm(1)%Aerosol(15)%Type = 13
    atm(1)%Aerosol(15)%Effective_Radius = Reff_j
    atm(1)%Aerosol(15)%Concentration = ECjHB(1,1,:,1)
    atm(1)%Aerosol(16)%Type = 13
    atm(1)%Aerosol(16)%Effective_Radius = Reff_k
    atm(1)%Aerosol(16)%Concentration = ECkHB(1,1,:,1)

    ! Sulfate
    atm(1)%Aerosol(17)%Type = 15
    atm(1)%Aerosol(17)%Effective_Radius = Reff_i
    atm(1)%Aerosol(17)%Concentration = SFi(1,1,:,1)
    atm(1)%Aerosol(18)%Type = 16
    atm(1)%Aerosol(18)%Effective_Radius = Reff_j
    atm(1)%Aerosol(18)%Concentration = SFj(1,1,:,1)

    ! Nitrate
    atm(1)%Aerosol(19)%Type = 17
    atm(1)%Aerosol(19)%Effective_Radius = Reff_i
    atm(1)%Aerosol(19)%Concentration = NIi(1,1,:,1)
    atm(1)%Aerosol(20)%Type = 18
    atm(1)%Aerosol(20)%Effective_Radius = Reff_j
    atm(1)%Aerosol(20)%Concentration = NIi(1,1,:,1)

  END SUBROUTINE Load_Atm_Data



  !
  ! Internal subprogam to load some test surface data
  !
  SUBROUTINE Load_Sfc_Data()


    ! 4a.0 Surface type definitions for default SfcOptics definitions
    !      For IR and VIS, this is the NPOESS reflectivities.
    ! ---------------------------------------------------------------
    INTEGER, PARAMETER :: TUNDRA_SURFACE_TYPE         = 10  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: SCRUB_SURFACE_TYPE          =  7  ! NPOESS Land surface type for IR/VIS Land SfcOptics
    INTEGER, PARAMETER :: COARSE_SOIL_TYPE            =  1  ! Soil type                for MW land SfcOptics
    INTEGER, PARAMETER :: GROUNDCOVER_VEGETATION_TYPE =  7  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: BARE_SOIL_VEGETATION_TYPE   = 11  ! Vegetation type          for MW Land SfcOptics
    INTEGER, PARAMETER :: SEA_WATER_TYPE              =  1  ! Water type               for all SfcOptics
    INTEGER, PARAMETER :: FRESH_SNOW_TYPE             =  2  ! NPOESS Snow type         for IR/VIS SfcOptics
    INTEGER, PARAMETER :: FRESH_ICE_TYPE              =  1  ! NPOESS Ice type          for IR/VIS SfcOptics



    ! 4a.1 Profile #1
    ! ---------------
    ! ...Land surface characteristics
    sfc(1)%Land_Coverage     = 0.1_fp
    sfc(1)%Land_Type         = TUNDRA_SURFACE_TYPE
    sfc(1)%Land_Temperature  = 272.0_fp
    sfc(1)%Lai               = 0.17_fp
    sfc(1)%Soil_Type         = COARSE_SOIL_TYPE
    sfc(1)%Vegetation_Type   = GROUNDCOVER_VEGETATION_TYPE
    ! ...Water surface characteristics
    sfc(1)%Water_Coverage    = 0.5_fp
    sfc(1)%Water_Type        = SEA_WATER_TYPE
    sfc(1)%Water_Temperature = 275.0_fp
    ! ...Snow coverage characteristics
    sfc(1)%Snow_Coverage    = 0.25_fp
    sfc(1)%Snow_Type        = FRESH_SNOW_TYPE
    sfc(1)%Snow_Temperature = 265.0_fp
    ! ...Ice surface characteristics
    sfc(1)%Ice_Coverage    = 0.15_fp
    sfc(1)%Ice_Type        = FRESH_ICE_TYPE
    sfc(1)%Ice_Temperature = 269.0_fp

  END SUBROUTINE Load_Sfc_Data


  !
  ! Internal subprogam to create a signal file
  !
  SUBROUTINE SignalFile_Create()
    CHARACTER(256) :: Filename
    INTEGER :: fid
    Filename = '.signal'
    fid = Get_Lun()
    OPEN( fid, FILE = Filename )
    WRITE( fid,* ) TRIM(Filename)
    CLOSE( fid )
  END SUBROUTINE SignalFile_Create

END PROGRAM check_crtm
