!
! TB_to_Radiance
!
! Example program to convert input TBs to Radiance for a given sensor
!

PROGRAM TB_to_Radiance

  ! Module usage
  USE CRTM_Module
  USE CRTM_Planck_Functions,      ONLY: CRTM_Planck_Radiance, CRTM_Planck_Temperature
  USE CRTM_SpcCoeff  , ONLY: SC
  ! Disable all implicit typing
  IMPLICIT NONE

  ! ============================================================================
  CHARACTER(*), PARAMETER :: PROGRAM_NAME   = 'TB_to_Radiance'

  ! Directory location of coefficients
  CHARACTER(*), PARAMETER :: COEFFICIENT_PATH='../test/testinput/'

  ! Aerosol/Cloud coefficient format
  CHARACTER(LEN=128) :: Coeff_Format, SpcCoeff_Format, TauCoeff_Format
  !CHARACTER(*), PARAMETER :: Coeff_Format = 'netCDF'

  ! Aerosol/Cloud coefficient scheme
  CHARACTER(*), PARAMETER :: Aerosol_Model = 'CRTM'
  CHARACTER(*), PARAMETER :: Cloud_Model   = 'CRTM'

  ! Profile dimensions
  INTEGER, PARAMETER :: N_PROFILES  = 1
  INTEGER, PARAMETER :: N_LAYERS    = 92
  INTEGER, PARAMETER :: N_ABSORBERS = 1
  INTEGER, PARAMETER :: N_CLOUDS    = 0
  INTEGER, PARAMETER :: N_AEROSOLS  = 0

  ! Sensor information
  INTEGER     , PARAMETER :: N_SENSORS = 1
  CHARACTER(*), PARAMETER :: SENSOR_ID(N_SENSORS) = (/'gmi_gpm'/)  !** change to your sensor here
  REAL(fp)                :: BT(13)  !** this is your input TB array size (n_channels), 13 in the GMI example.
  

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
  INTEGER :: err_stat
  INTEGER :: n_channels
  INTEGER :: n, i, SensorIndex, ChannelIndex
  REAL(fp) :: Radiance, tmp
  TYPE(CRTM_ChannelInfo_type)             :: chinfo(N_SENSORS)

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

  Coeff_Format = 'Binary'
     
  ! ... Coefficient table format
  IF ( Coeff_Format == 'Binary' ) THEN
    AerosolCoeff_Format = 'Binary'
    AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'bin'
    CloudCoeff_Format   = 'Binary'
    CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'bin'
    SpcCoeff_Format     = 'Binary'
    TauCoeff_Format     = 'Binary'
  ELSE IF ( Coeff_Format == 'netCDF' ) THEN
    AerosolCoeff_Format = 'netCDF'
    AerosolCoeff_File   = 'AerosolCoeff.'//TRIM(Aerosol_Scheme)//'nc4'
    CloudCoeff_Format   = 'netCDF'
    CloudCoeff_File     = 'CloudCoeff.'//TRIM(Cloud_Scheme)//'nc4'
    SpcCoeff_Format     = 'netCDF'
    TauCoeff_Format     = 'netCDF'
  END IF

  WRITE( *,'(/5x,"Initializing the CRTM...")' )
  err_stat = CRTM_Init( SENSOR_ID          , &
                        chinfo             , &
                        Aerosol_Model      , &
                        AerosolCoeff_Format, &
                        AerosolCoeff_File  , &
                        Cloud_Model        , &
                        CloudCoeff_Format  , &
                        CloudCoeff_File    , &
                        SpcCoeff_Format    , &
                        TauCoeff_Format    , &
                        File_Path=COEFFICIENT_PATH      , &
                        Quiet=.FALSE.)
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF

  ! Output some channel information
  ! -----------------------------------
  n_channels = SUM(CRTM_ChannelInfo_n_Channels(chinfo))
  WRITE( *,'(/5x,"Processing a total of ",i0," channels...")' ) n_channels
  DO n = 1, N_SENSORS
    WRITE( *,'(7x,i0," from ",a)' ) &
      CRTM_ChannelInfo_n_Channels(chinfo(n)), TRIM(SENSOR_ID(n))
  END DO
  ! ============================================================================

  Sensor_Loop: DO n = 1, N_SENSORS
     
     n_channels = CRTM_ChannelInfo_n_Channels(chinfo(n))

     ! example TB data
     !         1     2     3     4     5     6     7     8     9     10    11    12     13
     BT = (/243.0,253.0,263.2,223.2,253.3,201.9,210.8,215.0,211.5,217.2,211.0,200.3, 186.3/) !** 13 channels for GMI in this example
     
     SensorIndex = chinfo(n)%Sensor_Index
     
     DO i = 1,n_channels
        ChannelIndex = chinfo(n)%Channel_Index(i)
        
        CALL CRTM_Planck_Radiance( &
             SensorIndex  , &
             ChannelIndex , &
             BT(i)        , &  !* input
             Radiance     )    !* output
             
        PRINT *, 'A',  i, BT(i), Radiance

        !** check to ensure it works both ways
        CALL CRTM_Planck_Temperature( &
             SensorIndex  , &
             ChannelIndex , &
             Radiance     , &    !* input
             tmp           )  !* output TB

        PRINT *, 'B', i, Radiance, tmp
        
     END DO
  END DO Sensor_Loop

  ! ==========================================================================
  ! **** DESTROY THE CRTM ****
  !
  WRITE( *, '( /5x, "Destroying the CRTM..." )' )
  err_stat = CRTM_Destroy( chinfo )
  IF ( err_stat /= SUCCESS ) THEN
    message = 'Error destroying CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  END IF
  ! ==========================================================================

END PROGRAM TB_to_Radiance
