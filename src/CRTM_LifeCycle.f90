!
! CRTM_LifeCycle
!
! Module containing CRTM life cycle functions to initialize and destroy
! the CRTM space.
!
! Written by:     Paul van Delst, 21-May-2004
!                 paul.vandelst@noaa.gov
!
!
! Record of Revisions:
! ====================
!
!
! Date:           Name:                    Description:
! =====           =====                    ============
!
! 2020-08-01      Cheng Dang               Add optional arguments for
!                                          aerosol/cloud scheme/file/format.
!
! 2021-07-26      Patrick Stegmann         Add optional format input for
!                                          TauCoeff files.
!
! 2022-03-09      Cheng Dang               Add optional format input for
!                                          EmisCoeff files.
!
! 2022-05-27      Cheng Dang               Add optional input file for
!                                          Snow Emissivity files.

MODULE CRTM_LifeCycle


  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Message_Handler
  USE File_Utility           , ONLY: File_Exists
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type, &
                                     CRTM_ChannelInfo_Associated, &
                                     CRTM_ChannelInfo_Destroy, &
                                     CRTM_ChannelInfo_Create
  ! ...Spectral coefficients
  USE CRTM_SpcCoeff          , ONLY: SC, &
                                     CRTM_SpcCoeff_Load   , &
                                     CRTM_SpcCoeff_Destroy, &
                                     SpcCoeff_IsMicrowaveSensor  , &
                                     SpcCoeff_IsInfraredSensor   , &
                                     SpcCoeff_IsVisibleSensor    , &
                                     SpcCoeff_IsUltravioletSensor
  ! ...Transmittance model coefficients
  USE CRTM_TauCoeff
  ! ...Aerosol optical properties
  USE CRTM_AerosolCoeff      , ONLY: CRTM_AerosolCoeff_Load, &
                                     CRTM_AerosolCoeff_Destroy
  ! ...Cloud optical properties
  USE CRTM_CloudCoeff        , ONLY: CRTM_CloudCoeff_Load, &
                                     CRTM_CloudCoeff_Destroy
  ! ...Infrared surface emissivities
  USE CRTM_IRwaterCoeff      , ONLY: CRTM_IRwaterCoeff_Load, &
                                     CRTM_IRwaterCoeff_Destroy
  USE CRTM_IRlandCoeff       , ONLY: CRTM_IRlandCoeff_Load, &
                                     CRTM_IRlandCoeff_Destroy
  USE CRTM_IRsnowCoeff       , ONLY: CRTM_IRsnowCoeff_Load, &
                                     CRTM_IRsnowCoeff_Destroy
  USE CRTM_IRiceCoeff        , ONLY: CRTM_IRiceCoeff_Load, &
                                     CRTM_IRiceCoeff_Destroy
  ! ...Visible surface emissivities
  USE CRTM_VISwaterCoeff     , ONLY: CRTM_VISwaterCoeff_Load, &
                                     CRTM_VISwaterCoeff_Destroy
  USE CRTM_VISlandCoeff      , ONLY: CRTM_VISlandCoeff_Load, &
                                     CRTM_VISlandCoeff_Destroy
  USE CRTM_VISsnowCoeff      , ONLY: CRTM_VISsnowCoeff_Load, &
                                     CRTM_VISsnowCoeff_Destroy
  USE CRTM_VISiceCoeff       , ONLY: CRTM_VISiceCoeff_Load, &
                                     CRTM_VISiceCoeff_Destroy
  ! ...Microwave surface emissivities
  USE CRTM_MWwaterCoeff      , ONLY: CRTM_MWwaterCoeff_Load, &
                                     CRTM_MWwaterCoeff_Destroy, &
                                     CRTM_MWwaterCoeff_Load_FASTEM
  USE CRTM_PARMIOCoeff       , ONLY: CRTM_PARMIOCoeff_Load, &
                                     CRTM_PARMIOCoeff_Destroy, &
                                     CRTM_PARMIOCoeff_IsLoaded
  ! ...OpenMP API
#ifdef _OPENMP
  USE OMP_LIB
#endif
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Public procedures
  PUBLIC :: CRTM_Init
  PUBLIC :: CRTM_Destroy
  PUBLIC :: CRTM_IsInitialized

  ! -----------------
  ! Module parameters
  ! -----------------
  ! String lengths
  INTEGER, PARAMETER :: ML = 256   ! Error message length
  INTEGER, PARAMETER :: SL = 5000  ! Maximum length for path+filenames


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Init
!
! PURPOSE:
!       Function to initialise the CRTM.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Init( Sensor_ID  , &
!                                 ChannelInfo, &
!                                 Aerosol_Model       = Aerosol_Model       , &
!                                 AerosolCoeff_Format = AerosolCoeff_Format , &
!                                 AerosolCoeff_File   = AerosolCoeff_File   , &
!                                 Cloud_Model         = Cloud_Model         , &
!                                 CloudCoeff_Format   = CloudCoeff_Format   , &
!                                 CloudCoeff_File     = CloudCoeff_File     , &
!                                 SpcCoeff_Format     = SpcCoeff_Format     , &
!                                 TauCoeff_Format     = TauCoeff_Format     , &
!                                 Load_CloudCoeff     = Load_CloudCoeff     , &
!                                 Load_AerosolCoeff   = Load_AerosolCoeff   , &
!                                 IRwaterCoeff_File   = IRwaterCoeff_File   , &
!                                 IRlandCoeff_File    = IRlandCoeff_File    , &
!                                 IRsnowCoeff_File    = IRsnowCoeff_File    , &
!                                 IRiceCoeff_File     = IRiceCoeff_File     , &
!                                 VISwaterCoeff_File  = VISwaterCoeff_File  , &
!                                 VISlandCoeff_File   = VISlandCoeff_File   , &
!                                 VISsnowCoeff_File   = VISsnowCoeff_File   , &
!                                 VISiceCoeff_File    = VISiceCoeff_File    , &
!                                 MWwaterCoeff_File   = MWwaterCoeff_File   , &
!                                 IRwaterCoeff_Format = IRwaterCoeff_Format , &
!                                 IRlandCoeff_Format  = IRlandCoeff_Format  , &
!                                 IRiceCoeff_Format   = IRiceCoeff_Format   , &
!                                 VISwaterCoeff_Format= VISwaterCoeff_Format, &
!                                 VISlandCoeff_Format = VISlandCoeff_Format , &
!                                 VISsnowCoeff_Format = VISsnowCoeff_Format , &
!                                 VISiceCoeff_Format  = VISiceCoeff_Format  , &
!                                 File_Path           = File_Path           , &
!                                 NC_File_Path        = NC_File_Path        , &
!                                 Quiet               = Quiet               , &
!                                 Process_ID          = Process_ID          , &
!                                 Output_Process_ID   = Output_Process_ID   )
!
! INPUTS:
!       Sensor_ID:          List of the sensor IDs (e.g. hirs3_n17, amsua_n18,
!                           ssmis_f16, etc) with which the CRTM is to be
!                           initialised. These sensor ids are used to construct
!                           the sensor specific SpcCoeff and TauCoeff filenames
!                           containing the necessary coefficient data, i.e.
!                             <Sensor_ID>.SpcCoeff.bin
!                           and
!                             <Sensor_ID>.TauCoeff.bin
!                           for each sensor Id in the list.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Rank-1 (n_Sensors)
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OUTPUTS:
!       ChannelInfo:        ChannelInfo structure array populated based on
!                           the contents of the coefficient files and the
!                           user inputs.
!                           UNITS:      N/A
!                           TYPE:       CRTM_ChannelInfo_type
!                           DIMENSION:  Same as input Sensor_Id argument
!                           ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUTS:
!       Aerosol_Model:     Name of the aerosol scheme for scattering calculation
!                          Available aerosol scheme:
!                          - CRTM  [DEFAULT]
!                          - CMAQ
!                          - GOCART-GEOS5
!                          - NAAPS
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AerosolCoeff_Format:    Format of the aerosol optical properties data
!                               Available options:
!                               - Binary  [DEFAULT]
!                               - netCDF
!                               UNITS:      N/A
!                               TYPE:       CHARACTER(*)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AerosolCoeff_File:  Name of the data file containing the aerosol optical
!                           properties data for scattering calculations.
!                           Available datafiles:
!                           CRTM:
!                           - AerosolCoeff.bin      [DEFAULT, Binary]
!                           - AerosolCoeff.nc/nc4   [netCDF-Classic/4]
!                           CMAQ:
!                           - AerosolCoeff.CMAQ.bin      [Binary]
!                           - AerosolCoeff.CMAQ.nc/nc4   [netCDF-Classic/4]
!                           GOCART-GEOS5:
!                           - AerosolCoeff.GOCART-GEOS5.bin      [Binary]
!                           - AerosolCoeff.GOCART-GEOS5.nc/nc4   [netCDF-Classic/4]
!                           NAAPS:
!                           - AerosolCoeff.NAAPS.bin      [Binary]
!                           - AerosolCoeff.NAAPS.nc/nc4   [netCDF-Classic/4]
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Cloud_Model:       Name of the cloud scheme for scattering calculation
!                          Available cloud scheme:
!                          - CRTM  [DEFAULT]
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       CloudCoeff_Format:     Format of the cloud optical properties data
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       CloudCoeff_File:    Name of the data file containing the cloud optical
!                           properties data for scattering calculations.
!                           Available datafiles:
!                           - CloudCoeff.bin  [DEFAULT, Binary]
!                           - CloudCoeff.nc     [netCDF-Classic/4]
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       SpcCoeff_Format:     Format of the CRTM spectral coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       TauCoeff_Format:     Format of the CRTM transmittance coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Load_CloudCoeff:    Set this logical argument for not loading the CloudCoeff data
!                           to save memory space under the clear conditions
!                           If == .FALSE., the CloudCoeff data will not be loaded;
!                              == .TRUE.,  the CloudCoeff data will be loaded.
!                           If not specified, default is .TRUE. (will be loaded)
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Load_AerosolCoeff:  Set this logical argument for not loading the AerosolCoeff data
!                           to save memory space under the clear conditions
!                           If == .FALSE., the AerosolCoeff data will not be loaded;
!                              == .TRUE.,  the AerosolCoeff data will be loaded.
!                           If not specified, default is .TRUE. (will be loaded)
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       MWwaterCoeff_File:  Name of the data file containing the coefficient
!                           data for the microwave water emissivity model.
!                           Available datafiles:
!                           - FASTEM6.MWwater.EmisCoeff.bin  [DEFAULT]
!                           - FASTEM5.MWwater.EmisCoeff.bin
!                           - FASTEM4.MWwater.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       IRwaterCoeff_File:  Name of the data file containing the coefficient
!                           data for the infrared water emissivity model.
!                           Available datafiles:
!                           - Nalli.IRwater.EmisCoeff.bin  [DEFAULT]
!                           - WuSmith.IRwater.EmisCoeff.bin
!                           If not specified the Nalli datafile is read.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       IRlandCoeff_File:   Name of the data file containing the coefficient
!                           data for the infrared land emissivity model.
!                           Available datafiles:
!                           - NPOESS.IRland.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.IRland.EmisCoeff.bin
!                           - USGS.IRland.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       IRsnowCoeff_File:   Name of the data file containing the coefficient
!                           data for the infrared snow emissivity model.
!                           Available datafiles:
!                           - NPOESS.IRsnow.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.IRsnow.EmisCoeff.bin
!                           - USGS.IRsnow.EmisCoeff.bin
!                           - Nalli.IRsnow.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       IRiceCoeff_File:    Name of the data file containing the coefficient
!                           data for the infrared ice emissivity model.
!                           Available datafiles:
!                           - NPOESS.IRice.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.IRice.EmisCoeff.bin
!                           - USGS.IRice.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       VISwaterCoeff_File: Name of the data file containing the coefficient
!                           data for the visible water emissivity model.
!                           Available datafiles:
!                           - NPOESS.VISwater.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.VISwater.EmisCoeff.bin
!                           - USGS.VISwater.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       VISlandCoeff_File:  Name of the data file containing the coefficient
!                           data for the visible land emissivity model.
!                           Available datafiles:
!                           - NPOESS.VISland.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.VISland.EmisCoeff.bin
!                           - USGS.VISland.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       VISsnowCoeff_File:  Name of the data file containing the coefficient
!                           data for the visible snow emissivity model.
!                           Available datafiles:
!                           - NPOESS.VISsnow.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.VISsnow.EmisCoeff.bin
!                           - USGS.VISsnow.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       VISiceCoeff_File:   Name of the data file containing the coefficient
!                           data for the visible ice emissivity model.
!                           Available datafiles:
!                           - NPOESS.VISice.EmisCoeff.bin  [DEFAULT]
!                           - IGBP.VISice.EmisCoeff.bin
!                           - USGS.VISice.EmisCoeff.bin
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!     IRwaterCoeff_Format:   Format of the CRTM IRwater coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!      IRlandCoeff_Format:   Format of the CRTM IRland coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!      IRsnowCoeff_Format:   Format of the CRTM IRsnow coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       IRiceCoeff_Format:   Format of the CRTM IRice coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!    VISwaterCoeff_Format:   Format of the CRTM VISwater coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!     VISlandCoeff_Format:   Format of the CRTM VISland coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!     VISsnowCoeff_Format:   Format of the CRTM VISsnow coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!     VISiceCoeff_Format:   Format of the CRTM VISice coefficients
!                              Available options
!                              - Binary  [DEFAULT]
!                              - netCDF
!                              UNITS:      N/A
!                              TYPE:       CHARACTER(*)
!                              DIMENSION:  Scalar
!                              ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       File_Path:          Character string specifying a file path for the
!                           input data files in Binary format. If not specified,
!                           the current directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       NC_File_Path:       Character string specifying a file path for the
!                           input data files in netCDF format. If not specified,
!                           the current directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATION message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION message!       File_Path:          Character string specifying a file path for the
!                           input data files. If not specified, the current
!                           directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONALs are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output.
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the CRTM initialisation was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       All public data arrays accessed by this module and its dependencies
!       are overwritten.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Init( &
    Sensor_ID           , &  ! Input
    ChannelInfo         , &  ! Output
    Aerosol_Model       , &  ! Optional input
    AerosolCoeff_Format , &  ! Optional input
    AerosolCoeff_File   , &  ! Optional input
    Cloud_Model         , &  ! Optional input
    CloudCoeff_Format   , &  ! Optional input
    CloudCoeff_File     , &  ! Optional input
    SpcCoeff_Format     , &  ! Optional input
    TauCoeff_Format     , &  ! Optional input
    EmisCoeff_File      , &  ! Optional input  ! *** DEPRECATED. Replaced by IRwaterCoeff_File
    IRwaterCoeff_File   , &  ! Optional input
    IRlandCoeff_File    , &  ! Optional input
    IRsnowCoeff_File    , &  ! Optional input
    IRiceCoeff_File     , &  ! Optional input
    VISwaterCoeff_File  , &  ! Optional input
    VISlandCoeff_File   , &  ! Optional input
    VISsnowCoeff_File   , &  ! Optional input
    VISiceCoeff_File    , &  ! Optional input
    MWwaterCoeff_File   , &  ! Optional input
    MWwaterCoeff_Scheme , &  ! Optional input
    PARMIOCoeff_File    , &  ! Optional input
    IRwaterCoeff_Format , &  ! Optional input
    IRlandCoeff_Format  , &  ! Optional input
    IRsnowCoeff_Format  , &  ! Optional input
    IRiceCoeff_Format   , &  ! Optional input
    VISwaterCoeff_Format, &  ! Optional input
    VISlandCoeff_Format , &  ! Optional input
    VISsnowCoeff_Format , &  ! Optional input
    VISiceCoeff_Format  , &  ! Optional input
    IRsnow_Model        , &  ! Optional input
    File_Path           , &  ! Optional input
    NC_File_Path        , &  ! Optional input
    Load_CloudCoeff     , &  ! Optional input
    Load_AerosolCoeff   , &  ! Optional input
    Quiet               , &  ! Optional input
    Process_ID          , &  ! Optional input
    Output_Process_ID )   &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*)               , INTENT(IN)  :: Sensor_ID(:)
    TYPE(CRTM_ChannelInfo_type), INTENT(OUT) :: ChannelInfo(:)
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: Aerosol_Model
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: AerosolCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: AerosolCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: Cloud_Model
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: CloudCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: CloudCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: SpcCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: TauCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: EmisCoeff_File  ! *** DEPRECATED. Replaced by IRwaterCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRwaterCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRlandCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRsnowCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRiceCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISwaterCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISlandCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISsnowCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISiceCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: MWwaterCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: MWwaterCoeff_Scheme
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: PARMIOCoeff_File
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRwaterCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRlandCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRsnowCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRiceCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISwaterCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISlandCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISsnowCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: VISiceCoeff_Format
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: IRsnow_Model
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: File_Path
    CHARACTER(*),      OPTIONAL, INTENT(IN)  :: NC_File_Path
    LOGICAL     ,      OPTIONAL, INTENT(IN)  :: Load_CloudCoeff
    LOGICAL     ,      OPTIONAL, INTENT(IN)  :: Load_AerosolCoeff
    LOGICAL     ,      OPTIONAL, INTENT(IN)  :: Quiet
    INTEGER     ,      OPTIONAL, INTENT(IN)  :: Process_ID
    INTEGER     ,      OPTIONAL, INTENT(IN)  :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Init'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(SL) :: Default_Aerosol_Model
    CHARACTER(SL) :: Default_AerosolCoeff_Format
    CHARACTER(SL) :: Default_AerosolCoeff_File
    CHARACTER(SL) :: Default_Cloud_Model
    CHARACTER(SL) :: Default_CloudCoeff_Format
    CHARACTER(SL) :: Default_CloudCoeff_File
    CHARACTER(SL) :: Default_SpcCoeff_Format
    CHARACTER(SL) :: Default_TauCoeff_Format
    CHARACTER(SL) :: Default_IRwaterCoeff_File
    CHARACTER(SL) :: Default_IRlandCoeff_File
    CHARACTER(SL) :: Default_IRsnow_Model
    CHARACTER(SL) :: Default_IRsnowCoeff_File
    CHARACTER(SL) :: Default_IRiceCoeff_File
    CHARACTER(SL) :: Default_VISwaterCoeff_File
    CHARACTER(SL) :: Default_VISlandCoeff_File
    CHARACTER(SL) :: Default_VISsnowCoeff_File
    CHARACTER(SL) :: Default_VISiceCoeff_File
    CHARACTER(SL) :: Default_MWwaterCoeff_File
    CHARACTER(SL) :: Default_MWwaterCoeff_Scheme
    CHARACTER(SL) :: Resolved_PARMIOCoeff_File
    CHARACTER(SL) :: Default_IRwaterCoeff_Format
    CHARACTER(SL) :: Default_IRlandCoeff_Format
    CHARACTER(SL) :: Default_IRsnowCoeff_Format
    CHARACTER(SL) :: Default_IRsnowCoeff_Model
    CHARACTER(SL) :: Default_IRiceCoeff_Format
    CHARACTER(SL) :: Default_VISwaterCoeff_Format
    CHARACTER(SL) :: Default_VISlandCoeff_Format
    CHARACTER(SL) :: Default_VISsnowCoeff_Format
    CHARACTER(SL) :: Default_VISiceCoeff_Format
    CHARACTER(SL) :: Default_File_Path
    CHARACTER(SL) :: Effective_Bin_Path
    CHARACTER(SL) :: Effective_NC_Path
    CHARACTER(SL) :: Effective_Coeff_Path


    INTEGER :: l, n, n_Sensors
    LOGICAL :: Local_Load_CloudCoeff
    LOGICAL :: Local_Load_AerosolCoeff
    LOGICAL :: netCDF, isSEcategory
    LOGICAL :: Quiet_
    INTEGER :: iQuiet ! TODO: iQuiet should be removed once load routine interfaces have been modified
    Quiet_ = .TRUE.
    IF ( PRESENT(Quiet) ) Quiet_ = Quiet
    iQuiet = 0
    IF ( Quiet_ ) THEN
      iQuiet = 1
    END IF

    ! Honor OMP_NUM_THREADS from the runtime environment.
    ! If the variable is unset OR set but empty, default to a single thread.
    ! libgomp can corrupt the heap when OMP_NUM_THREADS is the empty string,
    ! so an empty value must be coerced to 1 the same as unset.
#ifdef _OPENMP
    BLOCK
      CHARACTER(32) :: omp_env_value
      INTEGER       :: omp_env_status
      CALL GET_ENVIRONMENT_VARIABLE( 'OMP_NUM_THREADS', &
                                     VALUE  = omp_env_value, &
                                     STATUS = omp_env_status )
      IF ( omp_env_status /= 0 .OR. LEN_TRIM(omp_env_value) == 0 ) &
        CALL OMP_SET_NUM_THREADS( 1 )
    END BLOCK
#endif

    ! Set up
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF
    ! ...Check coefficient loading flags
    Local_Load_CloudCoeff = .TRUE.
    IF( PRESENT(Load_CloudCoeff) ) Local_Load_CloudCoeff = Load_CloudCoeff
    Local_Load_AerosolCoeff = .TRUE.
    IF( PRESENT(Load_AerosolCoeff) ) Local_Load_AerosolCoeff = Load_AerosolCoeff
    ! ...Check dimensionality
    n_Sensors = SIZE(Sensor_ID)
    IF ( SIZE(ChannelInfo) /= n_Sensors ) THEN
      err_stat = FAILURE
      msg = 'Inconsistent Sensor_ID and ChannelInfo dimensions'
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
      RETURN
    END IF
    ! ...Check for deprecated EmisCoeff_File argument
    IF ( PRESENT(EmisCoeff_File) ) THEN
      err_stat = FAILURE
      msg = 'The EmisCoeff_File optional argument is deprecated. Use IRwaterCoeff_File instead'
      CALL Display_Message( ROUTINE_NAME, TRIM(msg)//TRIM(pid_msg), err_stat )
      RETURN
    END IF


    ! Specify sensor-independent coefficient filenames
    ! ...Default File_Path
    Default_File_Path = ''
    ! ...Default filenames
    Default_Aerosol_Model       = 'CRTM'
    Default_AerosolCoeff_File   = 'AerosolCoeff.nc'
    Default_Cloud_Model         = 'CRTM'
    Default_CloudCoeff_File     = 'CloudCoeff.nc'
    Default_IRwaterCoeff_File   = 'Nalli.IRwater.EmisCoeff.nc'
    Default_IRlandCoeff_File    = 'NPOESS.IRland.EmisCoeff.nc'
    Default_IRsnow_Model        = 'SEcategory'
    Default_IRsnowCoeff_File    = 'NPOESS.IRsnow.EmisCoeff.nc'
    Default_IRiceCoeff_File     = 'NPOESS.IRice.EmisCoeff.nc'
    Default_VISwaterCoeff_File  = 'NPOESS.VISwater.EmisCoeff.nc'
    Default_VISlandCoeff_File   = 'NPOESS.VISland.EmisCoeff.nc'
    Default_VISsnowCoeff_File   = 'NPOESS.VISsnow.EmisCoeff.nc'
    Default_VISiceCoeff_File    = 'NPOESS.VISice.EmisCoeff.nc'
    Default_MWwaterCoeff_File   = 'FASTEM6.MWwater.EmisCoeff.nc'
    Default_MWwaterCoeff_Scheme = 'FASTEM6'
    ! ... Default file formats (NetCDF is the canonical format from REL-3.2.0;
    !     Resolve_Coeff_Format below will fall back to Binary if a NetCDF file
    !     is missing on disk but the .bin equivalent is present.)
    Default_AerosolCoeff_Format = 'netCDF'
    Default_CloudCoeff_Format   = 'netCDF'
    Default_SpcCoeff_Format     = 'netCDF'
    Default_TauCoeff_Format     = 'netCDF'
    Default_IRwaterCoeff_Format = 'netCDF'
    Default_IRlandCoeff_Format  = 'netCDF'
    Default_IRsnowCoeff_Format  = 'netCDF'
    Default_IRiceCoeff_Format   = 'netCDF'
    Default_VISwaterCoeff_Format= 'netCDF'
    Default_VISlandCoeff_Format = 'netCDF'
    Default_VISsnowCoeff_Format = 'netCDF'
    Default_VISiceCoeff_Format  = 'netCDF'
    ! ...Were coefficient models specified?
    IF ( PRESENT(Aerosol_Model       ) ) Default_Aerosol_Model       = TRIM(ADJUSTL(Aerosol_Model))
    IF ( PRESENT(Cloud_Model         ) ) Default_Cloud_Model         = TRIM(ADJUSTL(Cloud_Model))
    IF ( PRESENT(IRsnow_Model        ) ) Default_IRsnow_Model        = TRIM(ADJUSTL(IRsnow_Model))
    ! ...Were other filenames specified?
    IF ( PRESENT(AerosolCoeff_File   ) ) Default_AerosolCoeff_File   = TRIM(ADJUSTL(AerosolCoeff_File))
    IF ( PRESENT(CloudCoeff_File     ) ) Default_CloudCoeff_File     = TRIM(ADJUSTL(CloudCoeff_File))
    IF ( PRESENT(IRwaterCoeff_File   ) ) Default_IRwaterCoeff_File   = TRIM(ADJUSTL(IRwaterCoeff_File))
    IF ( PRESENT(IRlandCoeff_File    ) ) Default_IRlandCoeff_File    = TRIM(ADJUSTL(IRlandCoeff_File))
    IF ( PRESENT(IRsnowCoeff_File    ) ) Default_IRsnowCoeff_File    = TRIM(ADJUSTL(IRsnowCoeff_File))
    IF ( PRESENT(IRiceCoeff_File     ) ) Default_IRiceCoeff_File     = TRIM(ADJUSTL(IRiceCoeff_File))
    IF ( PRESENT(VISwaterCoeff_File  ) ) Default_VISwaterCoeff_File  = TRIM(ADJUSTL(VISwaterCoeff_File))
    IF ( PRESENT(VISlandCoeff_File   ) ) Default_VISlandCoeff_File   = TRIM(ADJUSTL(VISlandCoeff_File))
    IF ( PRESENT(VISsnowCoeff_File   ) ) Default_VISsnowCoeff_File   = TRIM(ADJUSTL(VISsnowCoeff_File))
    IF ( PRESENT(VISiceCoeff_File    ) ) Default_VISiceCoeff_File    = TRIM(ADJUSTL(VISiceCoeff_File))
    IF ( PRESENT(MWwaterCoeff_File   ) ) Default_MWwaterCoeff_File   = TRIM(ADJUSTL(MWwaterCoeff_File))
    ! ...Were data formats specificed?
    IF ( PRESENT(AerosolCoeff_Format ) ) Default_AerosolCoeff_Format  = TRIM(ADJUSTL(AerosolCoeff_Format))
    IF ( PRESENT(CloudCoeff_Format   ) ) Default_CloudCoeff_Format    = TRIM(ADJUSTL(CloudCoeff_Format))
    IF ( PRESENT(SpcCoeff_Format     ) ) Default_SpcCoeff_Format      = TRIM(ADJUSTL(SpcCoeff_Format))
    IF ( PRESENT(TauCoeff_Format     ) ) Default_TauCoeff_Format      = TRIM(ADJUSTL(TauCoeff_Format))
    IF ( PRESENT(IRwaterCoeff_Format ) ) Default_IRwaterCoeff_Format  = TRIM(ADJUSTL(IRwaterCoeff_Format))
    IF ( PRESENT(IRlandCoeff_Format  ) ) Default_IRlandCoeff_Format   = TRIM(ADJUSTL(IRlandCoeff_Format))
    IF ( PRESENT(IRsnowCoeff_Format  ) ) Default_IRsnowCoeff_Format   = TRIM(ADJUSTL(IRsnowCoeff_Format))
    IF ( PRESENT(IRiceCoeff_Format   ) ) Default_IRiceCoeff_Format    = TRIM(ADJUSTL(IRiceCoeff_Format))
    IF ( PRESENT(VISwaterCoeff_Format) ) Default_VISwaterCoeff_Format = TRIM(ADJUSTL(VISwaterCoeff_Format))
    IF ( PRESENT(VISlandCoeff_Format ) ) Default_VISlandCoeff_Format  = TRIM(ADJUSTL(VISlandCoeff_Format))
    IF ( PRESENT(VISsnowCoeff_Format ) ) Default_VISsnowCoeff_Format  = TRIM(ADJUSTL(VISsnowCoeff_Format))
    IF ( PRESENT(VISiceCoeff_Format  ) ) Default_VISiceCoeff_Format   = TRIM(ADJUSTL(VISiceCoeff_Format))
    ! ...MW water emissivity scheme
    IF ( PRESENT(MWwaterCoeff_Scheme ) ) Default_MWwaterCoeff_Scheme  = TRIM(ADJUSTL(MWwaterCoeff_Scheme))

    ! ...Was a path specified?
    IF ( PRESENT(File_Path) ) THEN
      Default_MWwaterCoeff_File  = TRIM(ADJUSTL(File_Path)) // TRIM(Default_MWwaterCoeff_File)
    END IF

    ! ...Effective search paths for the format-resolver. NC_File_Path wins for
    !    NetCDF when it is supplied separately; otherwise both formats are
    !    expected to live under File_Path (the typical layout for the CRTM
    !    test-data tree and most external callers).
    Effective_Bin_Path = ''
    IF ( PRESENT(File_Path)    ) Effective_Bin_Path = File_Path
    IF ( PRESENT(NC_File_Path) ) THEN
      Effective_NC_Path = NC_File_Path
    ELSE
      Effective_NC_Path = Effective_Bin_Path
    END IF

    ! Load the spectral coefficients
    netCDF = .FALSE.
    IF (Default_SpcCoeff_Format == 'netCDF' ) THEN
        netCDF = .TRUE.
    END IF
    IF ( .NOT. Quiet_ ) THEN
      WRITE(*,*) "Loading "//TRIM(Default_SpcCoeff_Format)//" spectral coefficients."
    END IF
    err_stat = CRTM_SpcCoeff_Load( &
                 Sensor_ID                                 , &
                 File_Path         = Effective_NC_Path     , &
                 netCDF            = netCDF                , &
                 Quiet             = Quiet                 , &
                 Process_ID        = Process_ID            , &
                 Output_Process_ID = Output_Process_ID       )
    IF ( err_stat /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME,'Error loading SpcCoeff data'//TRIM(pid_msg),err_stat )
      RETURN
    END IF


    ! Load the transmittance model coefficients
    netCDF = .FALSE.
    IF (Default_TauCoeff_Format == 'netCDF' ) THEN
        netCDF = .TRUE.
    END IF
    IF ( .NOT. Quiet_ ) THEN
      WRITE(*,*) "Loading "//TRIM(Default_TauCoeff_Format)//" transmittance coefficients."
    END IF
    err_stat = CRTM_Load_TauCoeff( &
                 Sensor_ID         = Sensor_ID             , &
                 File_Path         = Effective_NC_Path     , &
                 Quiet             = iQuiet                , &  ! *** Use of iQuiet temporary
                 netCDF            = netCDF                , &
                 Process_ID        = Process_ID            , &
                 Output_Process_ID = Output_Process_ID       )
    IF ( err_stat /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME,'Error loading TauCoeff data'//TRIM(pid_msg),err_stat )
      RETURN
    END IF


    ! Load the cloud coefficients
    IF ( Local_Load_CloudCoeff ) THEN
      CALL Resolve_Coeff_Format( Default_CloudCoeff_File, Default_CloudCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_CloudCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading cloud coefficients: ", a) ') TRIM(Default_CloudCoeff_File)
      END IF
      err_stat = CRTM_CloudCoeff_Load( &
                   Default_Cloud_Model                  , &
                   Default_CloudCoeff_File              , &
                   File_Path         = Default_File_Path, &
                   netCDF            = netCDF           , &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading CloudCoeff data from '//TRIM(Default_CloudCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF


    ! Load the aerosol coefficients
    IF ( Local_Load_AerosolCoeff ) THEN
      CALL Resolve_Coeff_Format( Default_AerosolCoeff_File, Default_AerosolCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_AerosolCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading aerosol coefficients: ", a) ') TRIM(Default_AerosolCoeff_File)
      END IF
      err_stat = CRTM_AerosolCoeff_Load( &
                   Default_Aerosol_Model                , &
                   Default_AerosolCoeff_File            , &
                   File_Path         = Default_File_Path, &
                   netCDF            = netCDF           , &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading AerosolCoeff data from '//TRIM(Default_AerosolCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF


    ! Load the emissivity model coefficients
    ! ...Infrared
    Infrared_Sensor: IF ( ANY(SpcCoeff_IsInfraredSensor(SC)) ) THEN
      ! ...IR land
      CALL Resolve_Coeff_Format( Default_IRlandCoeff_File, Default_IRlandCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_IRlandCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading IR land emissivity coefficients: ", a) ') TRIM(Default_IRlandCoeff_File)
      END IF
      err_stat = CRTM_IRlandCoeff_Load( &
                   Default_IRlandCoeff_File, &
                   File_Path         = Default_File_Path, &
                   netCDF            = netCDF           , &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading IRlandCoeff data from '//TRIM(Default_IRlandCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...IR Water
      CALL Resolve_Coeff_Format( Default_IRwaterCoeff_File, Default_IRwaterCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_IRwaterCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading IR water emissivity coefficients: ", a) ') TRIM(Default_IRwaterCoeff_File)
      END IF
      err_stat = CRTM_IRwaterCoeff_Load( &
                   Default_IRwaterCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading IRwaterCoeff data from '//TRIM(Default_IRwaterCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...IR snow
      CALL Resolve_Coeff_Format( Default_IRsnowCoeff_File, Default_IRsnowCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_IRsnowCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF (Default_IRsnow_Model == 'SEcategory') THEN
        isSEcategory = .TRUE.
      ELSE
        isSEcategory = .FALSE.
      END IF
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading IR snow emissivity coefficients: ", a) ')  TRIM(Default_IRsnowCoeff_File)
      END IF
      err_stat = CRTM_IRsnowCoeff_Load( &
                   Default_IRsnowCoeff_File, &
                   netCDF            = netCDF           , &
                   isSEcategory      = isSEcategory     , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading IRsnowCoeff data from '//TRIM(Default_IRsnowCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...IR ice
      CALL Resolve_Coeff_Format( Default_IRiceCoeff_File, Default_IRiceCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_IRiceCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading IR ice emissivity coefficients: ", a) ') TRIM(Default_IRiceCoeff_File)
      END IF
      err_stat = CRTM_IRiceCoeff_Load( &
                   Default_IRiceCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading IRiceCoeff data from '//TRIM(Default_IRiceCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF Infrared_Sensor

    ! ...Visible
    Visible_Sensor: IF ( ANY(SpcCoeff_IsVisibleSensor(SC)) ) THEN
      ! ...VIS land
      CALL Resolve_Coeff_Format( Default_VISlandCoeff_File, Default_VISlandCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_VISlandCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading VIS land emissivity coefficients: ", a) ') TRIM(Default_VISlandCoeff_File)
      END IF
      err_stat = CRTM_VISlandCoeff_Load( &
                   Default_VISlandCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading VISlandCoeff data from '//TRIM(Default_VISlandCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...VIS water
      CALL Resolve_Coeff_Format( Default_VISwaterCoeff_File, Default_VISwaterCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_VISwaterCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading VIS water emissivity coefficients: ", a) ') TRIM(Default_VISwaterCoeff_File)
      END IF
      err_stat = CRTM_VISwaterCoeff_Load( &
                   Default_VISwaterCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading VISwaterCoeff data from '//TRIM(Default_VISwaterCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...VIS snow
      CALL Resolve_Coeff_Format( Default_VISsnowCoeff_File, Default_VISsnowCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_VISsnowCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading VIS snow emissivity coefficients: ", a) ') TRIM(Default_VISsnowCoeff_File)
      END IF
        err_stat = CRTM_VISsnowCoeff_Load( &
                   Default_VISsnowCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading VISsnowCoeff data from '//TRIM(Default_VISsnowCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...VIS ice
      CALL Resolve_Coeff_Format( Default_VISiceCoeff_File, Default_VISiceCoeff_Format, &
                                 Effective_NC_Path, Effective_Bin_Path, &
                                 Effective_Coeff_Path, Quiet=Quiet_ )
      netCDF = ( TRIM(Default_VISiceCoeff_Format) == 'netCDF' )
      Default_File_Path = Effective_Coeff_Path
      IF ( .NOT. Quiet_ ) THEN
        WRITE(*, '("Loading VIS ice emissivity coefficients: ", a) ') TRIM(Default_VISiceCoeff_File)
      END IF
      err_stat = CRTM_VISiceCoeff_Load( &
                   Default_VISiceCoeff_File, &
                   netCDF            = netCDF           , &
                   File_Path         = Default_File_Path, &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading VISiceCoeff data from '//TRIM(Default_VISiceCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
    END IF Visible_Sensor

    ! ...Microwave
    IF ( .NOT. Quiet_ ) THEN
      WRITE(*, '("Loading MW water emissivity coefficients: ", a) ') TRIM(Default_MWwaterCoeff_File)
    END IF
    Microwave_Sensor: IF ( ANY(SpcCoeff_IsMicrowaveSensor(SC)) ) THEN

      ! Note: this module is only needed for FASTEM5 scheme and the corresponding
      ! Binary lookup tables.
      ! ! ...MW water
      ! err_stat = CRTM_MWwaterCoeff_Load( &
      !              Default_MWwaterCoeff_File, &
      !              Quiet             = Quiet            , &
      !              Process_ID        = Process_ID       , &
      !              Output_Process_ID = Output_Process_ID  )
      ! IF ( err_stat /= SUCCESS ) THEN
      !   msg = 'Error loading MWwaterCoeff data from '//TRIM(Default_MWwaterCoeff_File)
      !   CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
      !   RETURN
      ! END IF

      ! ...MW water without LUT, support FASTEM4 and FASTEM6
      err_stat = CRTM_MWwaterCoeff_Load_FASTEM( &
                   Default_MWwaterCoeff_Scheme          , &
                   Quiet             = Quiet            , &
                   Process_ID        = Process_ID       , &
                   Output_Process_ID = Output_Process_ID  )
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error loading MWwaterCoeff module '//TRIM(Default_MWwaterCoeff_Scheme)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF

      ! ...Optional PARMIO LUT. Loaded only when the caller supplies a filename;
      !    Options%Use_PARMIO_Model toggles dispatch at Forward/TL/AD time. When
      !    the LUT is not loaded the dispatcher silently falls back to FASTEM,
      !    so omitting PARMIOCoeff_File is byte-identical to a pre-PARMIO build.
      IF ( PRESENT(PARMIOCoeff_File) ) THEN
        Resolved_PARMIOCoeff_File = TRIM(ADJUSTL(PARMIOCoeff_File))
        IF ( LEN_TRIM(Resolved_PARMIOCoeff_File) > 0 ) THEN
          IF ( PRESENT(File_Path) ) THEN
            Resolved_PARMIOCoeff_File = TRIM(ADJUSTL(File_Path)) // TRIM(Resolved_PARMIOCoeff_File)
          END IF
          IF ( .NOT. Quiet_ ) THEN
            WRITE(*, '("Loading PARMIO MW water emissivity LUT: ", a)') TRIM(Resolved_PARMIOCoeff_File)
          END IF
          err_stat = CRTM_PARMIOCoeff_Load( TRIM(Resolved_PARMIOCoeff_File), Quiet=Quiet )
          IF ( err_stat /= SUCCESS ) THEN
            msg = 'Error loading PARMIOCoeff data from '//TRIM(Resolved_PARMIOCoeff_File)
            CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
            RETURN
          END IF
        END IF
      END IF
    END IF Microwave_Sensor


    ! Load the ChannelInfo structure
    DO n = 1, n_Sensors
      ! ...Allocate the ChannelInfo structure
      CALL CRTM_ChannelInfo_Create( ChannelInfo(n), SC(n)%n_Channels )
      IF ( .NOT. CRTM_ChannelInfo_Associated(ChannelInfo(n)) ) THEN
        err_stat = FAILURE
        msg = 'ChannelInfo allocation failed for '//TRIM(Sensor_Id(n))//' sensor'
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      ! ...Fill the structure for the current sensor
      ChannelInfo(n)%Sensor_Index     = n
      ChannelInfo(n)%Channel_Index    = (/(l, l=1,SC(n)%n_Channels)/)
      ChannelInfo(n)%Sensor_ID        = SC(n)%Sensor_Id
      ChannelInfo(n)%Sensor_Type      = SC(n)%Sensor_Type
      ChannelInfo(n)%WMO_Satellite_ID = SC(n)%WMO_Satellite_ID
      ChannelInfo(n)%WMO_Sensor_ID    = SC(n)%WMO_Sensor_ID
      ChannelInfo(n)%Sensor_Channel   = SC(n)%Sensor_Channel
    END DO

  END FUNCTION CRTM_Init


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Destroy
!
! PURPOSE:
!       Function to deallocate all the shared data arrays allocated and
!       populated during the CRTM initialization.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Destroy( ChannelInfo            , &
!                                    Process_ID = Process_ID  )
!
! OUTPUTS:
!       ChannelInfo:  Reinitialized ChannelInfo structure.
!                     UNITS:      N/A
!                     TYPE:       CRTM_ChannelInfo_type
!                     DIMENSION:  Rank-1
!                     ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       Process_ID:   Set this argument to the MPI process ID that this
!                     function call is running under. This value is used
!                     solely for controlling message output. If MPI is not
!                     being used, ignore this argument.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status: The return value is an integer defining the error
!                     status. The error codes are defined in the
!                     Message_Handler module.
!                     If == SUCCESS the CRTM deallocations were successful
!                        == FAILURE an unrecoverable error occurred.
!                     UNITS:      N/A
!                     TYPE:       INTEGER
!                     DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       All CRTM shared data arrays and structures are deallocated.
!
! COMMENTS:
!       Note the INTENT on the output ChannelInfo argument is IN OUT rather than
!       just OUT. This is necessary because the argument may be defined upon
!       input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_Destroy( &
    ChannelInfo, &  ! Output
    Process_ID ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_ChannelInfo_type), INTENT(IN OUT) :: ChannelInfo(:)
    INTEGER,           OPTIONAL, INTENT(IN)     :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Destroy'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    INTEGER :: Destroy_Status

    ! Set up
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF


    ! Destroy all the ChannelInfo structures
    CALL CRTM_ChannelInfo_Destroy( ChannelInfo )
    IF ( ANY(CRTM_ChannelInfo_Associated(ChannelInfo)) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating ChannelInfo structure(s)'
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF


    ! Destroy the shared data structure
    Destroy_Status = CRTM_VISiceCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared VISiceCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_VISsnowCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared VISsnowCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_VISwaterCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared VISwaterCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_VISlandCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared VISlandCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_IRiceCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared IRiceCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_IRsnowCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared IRsnowCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_IRwaterCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared IRwaterCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_IRlandCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared IRlandCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_AerosolCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared AerosolCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_CloudCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared CloudCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_Destroy_TauCoeff( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared TauCoeff data structure(s)'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_SpcCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared SpcCoeff data structure(s)'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    Destroy_Status = CRTM_MWwaterCoeff_Destroy( Process_ID = Process_ID )
    IF ( Destroy_Status /= SUCCESS ) THEN
      err_stat = Destroy_Status
      msg = 'Error deallocating shared MWwaterCoeff data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
    END IF

    IF ( CRTM_PARMIOCoeff_IsLoaded() ) CALL CRTM_PARMIOCoeff_Destroy()

  END FUNCTION CRTM_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_IsInitialized
!
! PURPOSE:
!       Logical function to test if the CRTM has been correctly initialized.
!
! CALLING SEQUENCE:
!       status = CRTM_IsInitialized( ChannelInfo )
!
! INPUTS:
!       ChannelInfo:  ChannelInfo structure array.
!                     UNITS:      N/A
!                     TYPE:       CRTM_ChannelInfo_type
!                     DIMENSION:  Rank-1
!                     ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:       The return value is a logical result indicating if the
!                     CRTM has been correctly initialised.
!                     If == .TRUE., all the ChannelInfo entries are valid.
!                        == .FALSE., any of the ChannelInfo entries are invalid.
!                     UNITS:      N/A
!                     TYPE:       LOGICAL
!                     DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION CRTM_IsInitialized( ChannelInfo ) RESULT( Status )
    TYPE(CRTM_ChannelInfo_type), INTENT(IN) :: ChannelInfo(:)
    LOGICAL :: Status
    Status = ALL(CRTM_ChannelInfo_Associated(ChannelInfo))
  END FUNCTION CRTM_IsInitialized


!------------------------------------------------------------------------------
! Resolve_Coeff_Format
!
! Decide which format/file/directory a coefficient loader should actually use
! given a requested format and the user-supplied paths. If the requested file
! exists, return it unchanged. Otherwise, probe the alternate-format
! equivalent (NetCDF <-> Binary) and rewrite Filename, Format, and
! Resolved_Path in place so the caller transparently loads the available
! format. Supports the REL-3.2.0 Binary -> NetCDF transition.
!
! NC_Path and Bin_Path are the user's preferred directories for each format
! (typically NC_File_Path and File_Path arguments to CRTM_Init). When only
! one path is meaningful, callers can pass the same string for both.
!
! When falling back from .bin to NetCDF the helper prefers .nc4 then .nc,
! matching the on-disk extensions used by different coefficient categories.
!
! When neither file exists, the in/out arguments are left unchanged so the
! downstream loader's existing "file not found" error path triggers normally.
!------------------------------------------------------------------------------
  SUBROUTINE Resolve_Coeff_Format( Filename, Format, NC_Path, Bin_Path, &
                                   Resolved_Path, Quiet )
    CHARACTER(*),           INTENT(IN OUT) :: Filename
    CHARACTER(*),           INTENT(IN OUT) :: Format
    CHARACTER(*),           INTENT(IN)     :: NC_Path
    CHARACTER(*),           INTENT(IN)     :: Bin_Path
    CHARACTER(*),           INTENT(OUT)    :: Resolved_Path
    LOGICAL,      OPTIONAL, INTENT(IN)     :: Quiet
    ! Locals
    CHARACTER(LEN(Filename)) :: trial_name, requested_name
    CHARACTER(LEN(Format))   :: requested_format
    LOGICAL :: noisy
    INTEGER :: dot

    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet

    requested_name   = Filename
    requested_format = Format

    ! Default the resolved path to whichever directory matches the request.
    IF ( TRIM(Format) == 'netCDF' ) THEN
      Resolved_Path = NC_Path
    ELSE
      Resolved_Path = Bin_Path
    END IF

    ! Already on disk: nothing to do.
    IF ( File_Exists(TRIM(Resolved_Path)//TRIM(Filename)) ) RETURN

    dot = INDEX(Filename, '.', BACK=.TRUE.)
    IF ( dot < 1 ) RETURN  ! No extension; can't infer alternate.

    SELECT CASE ( TRIM(Filename(dot:)) )
    CASE ( '.nc', '.nc4' )
      ! Requested NetCDF; try the binary alternate in Bin_Path.
      trial_name = Filename(:dot-1) // '.bin'
      IF ( File_Exists(TRIM(Bin_Path)//TRIM(trial_name)) ) THEN
        Filename      = trial_name
        Format        = 'Binary'
        Resolved_Path = Bin_Path
      END IF

    CASE ( '.bin' )
      ! Requested Binary; try .nc4 then .nc in NC_Path.
      trial_name = Filename(:dot-1) // '.nc4'
      IF ( File_Exists(TRIM(NC_Path)//TRIM(trial_name)) ) THEN
        Filename      = trial_name
        Format        = 'netCDF'
        Resolved_Path = NC_Path
      ELSE
        trial_name = Filename(:dot-1) // '.nc'
        IF ( File_Exists(TRIM(NC_Path)//TRIM(trial_name)) ) THEN
          Filename      = trial_name
          Format        = 'netCDF'
          Resolved_Path = NC_Path
        END IF
      END IF

    END SELECT

    IF ( noisy .AND. Filename /= requested_name ) THEN
      CALL Display_Message( 'CRTM_Init', &
        'Requested '//TRIM(requested_format)//' file '//TRIM(requested_name)// &
        ' not found; falling back to '//TRIM(Format)//' file '//TRIM(Filename), &
        INFORMATION )
    END IF
  END SUBROUTINE Resolve_Coeff_Format

END MODULE CRTM_LifeCycle
