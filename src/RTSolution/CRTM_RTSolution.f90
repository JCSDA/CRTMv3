!
! CRTM_RTSolution
!
! Module containing the CRTM radiative transfer solution routines.
!
!
! CREATION HISTORY:
!       Written by:     Quanhua Liu,    QSS at JCSDA;    Quanhua.Liu@noaa.gov
!                       Yong Han,       NOAA/NESDIS;     Yong.Han@noaa.gov
!                       Paul van Delst, CIMSS/SSEC;      paul.vandelst@ssec.wisc.edu
!                       08-Jun-2004

MODULE CRTM_RTSolution

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds              , ONLY: fp
  USE Message_Handler         , ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters         , ONLY: ZERO, ONE, TWO, PI, &
                                      MAX_N_ANGLES, MAX_N_STOKES, &
                                      RT_ADA, RT_SOI, RT_AMOM, RT_VMOM
  USE Common_RTSolution       , ONLY: Assign_Common_Input, &
                                      Assign_Common_Output, &
                                      Assign_Common_Input_TL, &
                                      Assign_Common_Output_TL, &
                                      Assign_Common_Input_AD, &
                                      Assign_Common_Output_AD
  USE CRTM_SpcCoeff           , ONLY: SC
  USE CRTM_Atmosphere_Define  , ONLY: CRTM_Atmosphere_type
  USE CRTM_Surface_Define     , ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE CRTM_AtmOptics_Define   , ONLY: CRTM_AtmOptics_type
  USE CRTM_SfcOptics_Define   , ONLY: CRTM_SfcOptics_type
  USE CRTM_RTSolution_Define  , ONLY: CRTM_RTSolution_type
  USE CRTM_Utility
  USE RTV_Define
  ! RT modules
  USE SOI_Module
  USE ADA_Module
  USE Emission_Module
  ! Disable all implicit typing
  IMPLICIT NONE

  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  ! RTSolution structure entities
  ! ...Datatypes
  PUBLIC :: CRTM_RTSolution_type
  ! RTV structure entities
  ! ...Datatypes
  PUBLIC :: RTV_type
  ! Module procedures
  PUBLIC :: CRTM_Compute_RTSolution
  PUBLIC :: CRTM_Compute_RTSolution_TL
  PUBLIC :: CRTM_Compute_RTSolution_AD
  PUBLIC :: CRTM_Compute_nStreams
  PUBLIC :: CRTM_RTSolution_Version

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*),  PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_RTSolution.f90 60152 2015-08-13 19:19:13Z paul.vandelst@noaa.gov $'

CONTAINS

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_RTSolution
!
! PURPOSE:
!       Function to solve the radiative transfer equation.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_RTSolution( Atmosphere  , &  ! Input
!                                               Surface     , &  ! Input
!                                               AtmOptics   , &  ! Input
!                                               SfcOptics   , &  ! Input
!                                               GeometryInfo, &  ! Input
!                                               SensorIndex , &  ! Input
!                                               ChannelIndex, &  ! Input
!                                               RTSolution  , &  ! Output
!                                               RTV           )  ! Internal variable output
!
! INPUT ARGUMENTS:
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Surface:        Structure containing the surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:      Structure containing the surface optical properties
!                       data. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface
!                       optics properties before this routine is called.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!                       See the SensorIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution:     Structure containing the soluition to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       RTV:            Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output RTSolution argument is IN OUT rather than
!       just OUT. This is necessary because the argument is defined upon
!       input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_RTSolution( &
    Atmosphere  , &  ! Input
    Surface     , &  ! Input
    AtmOptics   , &  ! Input
    SfcOptics   , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    RTSolution  , &  ! Output
    RTV         ) &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN OUT) :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution
    TYPE(RTV_type),               INTENT(IN OUT) :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_RTSolution'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: i, nZ

    Error_Status = SUCCESS

    ! Populate the RTV structure
    Error_Status = Assign_Common_Input( Atmosphere  , &
                                        Surface     , &
                                        AtmOptics   , &
                                        SfcOptics   , &
                                        GeometryInfo, &
                                        SensorIndex , &
                                        ChannelIndex, &
                                        RTSolution  , &
                                        nz          , &
                                        RTV           )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning input for RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF


    ! Direct reflectivity normalisation fix for visible sensors
    IF( RTV%Visible_Flag_true ) THEN
      DO i = 1, nZ
! incorrect        SfcOptics%Direct_Reflectivity(i,1) = SfcOptics%Direct_Reflectivity(i,1) * PI
        ! ...Apply the UW limiter
        IF (SfcOptics%Direct_Reflectivity(i,1) > ONE) THEN
          SfcOptics%Direct_Reflectivity(i,1) = ONE
        END IF
      END DO
    END IF

    IF( RTV%n_Stokes > 1 ) THEN
      CALL Reshape_Surf_Opt(RTV%n_Angles, RTV%n_Stokes, SfcOptics%Emissivity, SfcOptics%Direct_Reflectivity, &
        SfcOptics%Reflectivity, SfcOptics%S_Emissivity, SfcOptics%S_Direct_Ref, SfcOptics%S_Reflectivity)
    ! ------------------------------
    ! Perform the radiative transfer
    ! ------------------------------
    ! Select the RT model
        IF( RTV%Scattering_RT ) THEN
!!          RTSolution%RT_Algorithm_Name = 'VMOM'
          ! NESDIS advanced adding-doubling method   

          CALL CRTM_ADA( &
               Atmosphere%n_Layers                       , & ! Input, number of atmospheric layers
               AtmOptics%Single_Scatter_Albedo           , & ! Input, layer single scattering albedo
               AtmOptics%Optical_Depth                   , & ! Input, layer optical depth
               RTV%Cosmic_Background_Radiance            , & ! Input, cosmic background radiation
               SfcOptics%S_Emissivity(:),   & ! Input, surface emissivity
               SfcOptics%S_Reflectivity(:,:), & ! Input, surface reflectivity for a point source
               SfcOptics%S_Direct_Ref(:),   & ! Input, surface reflectivity
               RTV, Error_Status                 ) ! Output, Internal variables
          IF( Error_Status /= SUCCESS  ) THEN
            WRITE( Message,'("Error in  CALL CRTM_ADA 2")' ) 
            CALL Display_Message( ROUTINE_NAME,  &                                                    
                            TRIM(Message), &                                                   
                            Error_Status   )                                          
            RETURN                                                                                    
          END IF
        ELSE
        RTSolution%RT_Algorithm_Name = 'Emission'
        CALL CRTM_Emission( &
             Atmosphere%n_Layers,                   & ! Input, number of atmospheric layers
             RTV%n_Angles,                          & ! Input, number of discrete zenith angles
             RTV%Diffuse_Surface,                   & ! Input, surface behavior
             GeometryInfo%Cosine_Sensor_Zenith,     & ! Input, cosine of sensor zenith angle
             AtmOptics%Optical_Depth,               & ! Input, layer optical depth
             RTV%Planck_Atmosphere,                 & ! Input, layer radiances
             RTV%Planck_Surface,                    & ! Input, surface radiance
             SfcOptics%S_Emissivity(:),          & ! Input, surface emissivity
             SfcOptics%S_Reflectivity(:,:), & ! Input, surface reflectivity
             SfcOptics%S_Direct_Ref(:), & ! Input, surface reflectivity for a point source
             RTV%Cosmic_Background_Radiance,        & ! Input, cosmic background radiation
             RTV%Solar_Irradiance,                  & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                  & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,     & ! Input, Source zenith angle
             RTV                                    ) ! Output, Internal variables     
        END IF
    
    ! ------------------------------
    ! Perform the radiative transfer
    ! ------------------------------
    ! Select the RT model
    ELSE IF( RTV%Scattering_RT ) THEN

      ! Select the scattering RT model
      SELECT CASE(RTV%RT_Algorithm_Id)
        CASE (RT_ADA, RT_VMOM) ;
          RTSolution%RT_Algorithm_Name = 'ADA'
          ! NESDIS advanced adding-doubling method
          CALL CRTM_ADA( &
               Atmosphere%n_Layers                       , & ! Input, number of atmospheric layers
               AtmOptics%Single_Scatter_Albedo           , & ! Input, layer single scattering albedo
               AtmOptics%Optical_Depth                   , & ! Input, layer optical depth
               RTV%Cosmic_Background_Radiance            , & ! Input, cosmic background radiation
               SfcOptics%Emissivity( 1:nZ, 1 )           , & ! Input, surface emissivity
               SfcOptics%Reflectivity( 1:nZ, 1, 1:nZ, 1 ), & ! Input, surface reflectivity
               SfcOptics%Direct_Reflectivity(1:nZ,1)     , & ! Input, surface reflectivity for a point source
               RTV, Error_Status        ) ! Output, Internal variables
        CASE (RT_SOI) ;
          RTSolution%RT_Algorithm_Name = 'SOI'
          ! UW SOI RT solver
          CALL CRTM_SOI( &
                 Atmosphere%n_Layers                       , & ! Input, number of atmospheric layers
                 AtmOptics%Single_Scatter_Albedo           , & ! Input, layer single scattering albedo
                 AtmOptics%Optical_Depth                   , & ! Input, layer optical depth
                 RTV%Cosmic_Background_Radiance            , & ! Input, cosmic background radiation
                 SfcOptics%Emissivity( 1:nZ, 1 )           , & ! Input, surface emissivity
                 SfcOptics%Reflectivity( 1:nZ, 1, 1:nZ, 1 ), & ! Input, surface reflectivity
                 SfcOptics%Index_Sat_Ang                   , & ! Input, Satellite angle index
                 RTV                                         ) ! Output, Internal variables
      CASE DEFAULT
      Error_Status = FAILURE
      WRITE(Message,'("Incorrect RT_Algorithm_ID, ",i0,", do not fit model")') &
                    RTV%RT_Algorithm_Id
      CALL Display_Message( ROUTINE_NAME,Message,Error_Status )
      RETURN      
      END SELECT

    ELSE

      ! -----------------
      ! Emission model RT
      ! -----------------
      RTSolution%RT_Algorithm_Name = 'Emission'
      CALL CRTM_Emission( &
             Atmosphere%n_Layers,                   & ! Input, number of atmospheric layers
             RTV%n_Angles,                          & ! Input, number of discrete zenith angles
             RTV%Diffuse_Surface,                   & ! Input, surface behavior
             GeometryInfo%Cosine_Sensor_Zenith,     & ! Input, cosine of sensor zenith angle
             AtmOptics%Optical_Depth,               & ! Input, layer optical depth
             RTV%Planck_Atmosphere,                 & ! Input, layer radiances
             RTV%Planck_Surface,                    & ! Input, surface radiance
             SfcOptics%Emissivity(1:nZ,1),          & ! Input, surface emissivity
             SfcOptics%Reflectivity(1:nZ,1,1:nZ,1), & ! Input, surface reflectivity
             SfcOptics%Direct_Reflectivity(1:nZ,1), & ! Input, surface reflectivity for a point source
             RTV%Cosmic_Background_Radiance,        & ! Input, cosmic background radiation
             RTV%Solar_Irradiance,                  & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                  & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,     & ! Input, Source zenith angle
             RTV                                    ) ! Output, Internal variables

    END IF

    Error_Status = Assign_Common_Output( Atmosphere,   &
                                         SfcOptics,    &
                                         GeometryInfo, &
                                         SensorIndex,  &
                                         ChannelIndex, &
                                         RTV,          &
                                         RTSolution    )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning output for RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF

  END FUNCTION CRTM_Compute_RTSolution

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_RTSolution_TL
!
! PURPOSE:
!       Function to solve the tangent-linear radiative transfer equation.
!
! CALLING SEQUENCE:
!      Error_Status = CRTM_Compute_RTSolution_TL( Atmosphere   , &  ! FWD Input
!                                                 Surface      , &  ! FWD Input
!                                                 AtmOptics    , &  ! FWD Input
!                                                 SfcOptics    , &  ! FWD Input
!                                                 RTSolution   , &  ! FWD Input
!                                                 Atmosphere_TL, &  ! TL Input
!                                                 Surface_TL   , &  ! TL Input
!                                                 AtmOptics_TL , &  ! TL Input
!                                                 SfcOptics_TL , &  ! TL Input
!                                                 GeometryInfo , &  ! Input
!                                                 SensorIndex  , &  ! Input
!                                                 ChannelIndex , &  ! Input
!                                                 RTSolution_TL, &  ! TL Output
!                                                 RTV            )  ! Internal variable input
!
! INPUT ARGUMENTS:
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Surface:        Structure containing the surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:      Structure containing the surface optical properties
!                       data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution:     Structure containing the solution to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:  Structure containing the tangent-linear atmospheric
!                       state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:     Structure containing the tangent-linear surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics_TL:   Structure containing the tangent-linear atmospheric
!                       optical properties.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_TL:   Structure containing the tangent-linear surface optical
!                       properties. Argument is defined as INTENT (IN OUT ) as
!                       different RT algorithms may compute the surface optics
!                       properties before this routine is called.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!                       See the SensorIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTV:            Structure containing internal forward model variables
!                       required for subsequent tangent-linear or adjoint model
!                       calls. The contents of this structure are NOT accessible
!                       outside of the CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
! OUTPUT ARGUMENTS:
!       RTSolution_TL:  Structure containing the solution to the tangent-linear
!                       RT equation for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output RTSolution_TL argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!--------------------------------------------------------------------------------




  FUNCTION CRTM_Compute_RTSolution_TL( &
    Atmosphere   , &  ! FWD Input
    Surface      , &  ! FWD Input
    AtmOptics    , &  ! FWD Input
    SfcOptics    , &  ! FWD Input
    RTSolution   , &  ! FWD Input
    Atmosphere_TL, &  ! TL Input
    Surface_TL   , &  ! TL Input
    AtmOptics_TL , &  ! TL Input
    SfcOptics_TL , &  ! TL Input
    GeometryInfo , &  ! Input
    SensorIndex  , &  ! Input
    ChannelIndex , &  ! Input
    RTSolution_TL, &  ! TL Output
    RTV          ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_RTSolution_type),   INTENT(IN)     :: RTSolution
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere_TL
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface_TL
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics_TL
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_TL
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution_TL
    TYPE(RTV_type),               INTENT(IN)     :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_RTSolution_TL'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: nZ
    REAL(fp) :: User_Emissivity_TL, Direct_Reflectivity_TL
    REAL(fp)                                     :: Planck_Surface_TL    ! Surface TL radiance
    REAL(fp), DIMENSION( 0:Atmosphere%n_Layers ) :: Planck_Atmosphere_TL ! *LAYER* TL radiances

    ! The following variables are RT model specific
    REAL(fp), DIMENSION( RTV%n_Angles * RTV%n_Stokes, &
                         (RTV%n_Angles+1) * RTV%n_Stokes, &
                         Atmosphere%n_Layers ) :: Pff_TL ! Forward scattering TL phase matrix
    REAL(fp), DIMENSION( RTV%n_Angles * RTV%n_Stokes, &
                         (RTV%n_Angles+1) * RTV%n_Stokes, &
                         Atmosphere%n_Layers ) :: Pbb_TL ! Backward scattering TL phase matrix
    REAL(fp), DIMENSION( RTV%n_Angles * RTV%n_Stokes ) :: Scattering_Radiance_TL
    REAL(fp) :: Radiance_TL

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS

    RTSolution_TL%RT_Algorithm_Name = RTSolution%RT_Algorithm_Name

    Error_Status = Assign_Common_Input_TL( Atmosphere             , &  ! FWD Input
                                           Surface                , &  ! FWD Input
                                           AtmOptics              , &  ! FWD Input
                                           SfcOptics              , &  ! FWD Input
                                           Atmosphere_TL          , &  ! TL Input
                                           Surface_TL             , &  ! TL Input
                                           AtmOptics_TL           , &  ! TL Input
                                           SfcOptics_TL           , &  ! TL Input Output
                                           GeometryInfo           , &  ! Input
                                           SensorIndex            , &  ! Input
                                           ChannelIndex           , &  ! Input
                                           RTSolution_TL          , &  ! TL Output
                                           nz                     , &  ! Output
                                           User_Emissivity_TL     , &  ! Output
                                           Direct_Reflectivity_TL , &  ! Output
                                           Planck_Surface_TL      , &  ! Output
                                           Planck_Atmosphere_TL   , &  ! Output
                                           Pff_TL                 , &  ! Output
                                           Pbb_TL                 , &  ! Output
                                           RTV          )   ! Internal variable input

    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning input for TL RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF

    nZ = RTV%n_Angles * RTV%n_Stokes
    IF( RTV%n_Stokes > 1 ) THEN
       CALL Reshape_Surf_Opt(RTV%n_Angles, RTV%n_Stokes, SfcOptics_TL%Emissivity, SfcOptics_TL%Direct_Reflectivity, &
        SfcOptics_TL%Reflectivity, SfcOptics_TL%S_Emissivity, SfcOptics_TL%S_Direct_Ref, SfcOptics_TL%S_Reflectivity)
    ! ---------------------------------------------
    ! Perform the tangent-linear radiative transfer
    ! ---------------------------------------------
    ! Select the RT model
      IF( RTV%Scattering_RT ) THEN
          CALL CRTM_ADA_TL( &
                 Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
                 AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
                 AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
                 RTV%Cosmic_Background_Radiance,           & ! cosmic background radiation
                 SfcOptics%S_Emissivity(1:nZ),             & ! Input, FWD surface emissivity
                 SfcOptics%S_Direct_Ref(1:nZ),    & ! Input, surface direct reflectivity
                 RTV,                                      & ! Input, structure containing forward results
                 Planck_Atmosphere_TL,                     & ! Input, TL layer radiances
                 Planck_Surface_TL,                        & ! Input, TL surface radiance
                 AtmOptics_TL%Single_Scatter_Albedo,       & ! Input, TL layer single scattering albedo
                 AtmOptics_TL%Optical_Depth,               & ! Input, TL layer optical depth
                 SfcOptics_TL%S_Emissivity(1:nZ),          & ! Input, TL surface emissivity
                 SfcOptics_TL%S_Reflectivity(1:nZ,1:nZ), & ! Input, TL surface reflectivity
                 SfcOptics_TL%S_Direct_Ref(1:nZ), & ! Input, TL surface direct reflectivity
                 Pff_TL(1:nZ,1:(nZ+1),:),                  & ! Input, TL layer forward phase matrix
                 Pbb_TL(1:nZ,1:(nZ+1),:),                  & ! Input, TL layer backward phase matrix
                 Scattering_Radiance_TL(1:nZ)              ) ! Output, TL radiances
      ELSE
        CALL CRTM_Emission_TL( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             RTV%n_Angles,                             & ! Input, number of discrete zenith angles
             GeometryInfo%Cosine_Sensor_Zenith,        & ! Input, cosine of sensor zenith angle
             RTV%Planck_Atmosphere,                    & ! Input, FWD layer radiances
             RTV%Planck_Surface,                       & ! Input, FWD surface radiance
             SfcOptics%S_Emissivity(1:nZ),             & ! Input, FWD surface emissivity
             SfcOptics%S_Reflectivity(1:nZ,1:nZ),    & ! Input, FWD surface reflectivity
             SfcOptics%S_Direct_Ref(1:nZ),    & ! Input, FWD surface reflectivity for a point source
             RTV%Solar_Irradiance,                     & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                     & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,        & ! Input, Source zenith angle
             RTV,                                      & ! Input, internal variables
             AtmOptics_TL%Optical_Depth,               & ! Input, TL layer optical depth
             Planck_Atmosphere_TL,                     & ! Input, TL layer radiances
             Planck_Surface_TL,                        & ! Input, TL surface radiance
             SfcOptics_TL%S_Emissivity(1:nZ),          & ! Input, TL surface emissivity
             SfcOptics_TL%S_Reflectivity(1:nZ,1:nZ), & ! Input, TL surface reflectivity
             SfcOptics_TL%S_Direct_Ref(1:nZ), & ! Input, TL surface reflectivity for a point source
             Radiance_TL                               ) ! Output, TL radiances    
      END IF

   ELSE IF( RTV%Scattering_RT ) THEN
      ! Select the scattering RT model
      SELECT CASE(RTV%RT_Algorithm_Id)
        CASE (RT_ADA, RT_VMOM)
          ! NESDIS advanced adding-doubling method    
          CALL CRTM_ADA_TL( &
                 Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
                 AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
                 AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
                 RTV%Cosmic_Background_Radiance,           & ! cosmic background radiation
                 SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
                 SfcOptics%Direct_Reflectivity(1:nZ,1),    & ! Input, surface direct reflectivity
                 RTV,                                      & ! Input, structure containing forward results
                 Planck_Atmosphere_TL,                     & ! Input, TL layer radiances
                 Planck_Surface_TL,                        & ! Input, TL surface radiance
                 AtmOptics_TL%Single_Scatter_Albedo,       & ! Input, TL layer single scattering albedo
                 AtmOptics_TL%Optical_Depth,               & ! Input, TL layer optical depth
                 SfcOptics_TL%Emissivity(1:nZ,1),          & ! Input, TL surface emissivity
                 SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1), & ! Input, TL surface reflectivity
                 SfcOptics_TL%Direct_Reflectivity(1:nZ,1), & ! Input, TL surface direct reflectivity
                 Pff_TL(1:nZ,1:(nZ+1),:),                  & ! Input, TL layer forward phase matrix
                 Pbb_TL(1:nZ,1:(nZ+1),:),                  & ! Input, TL layer backward phase matrix
                 Scattering_Radiance_TL(1:nZ)              ) ! Output, TL radiances

        CASE (RT_SOI)
          ! UW SOI RT solver
          CALL CRTM_SOI_TL( &
                 Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
                 AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
                 AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
                 SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
                 SfcOptics%Reflectivity(1:nZ,1,1:nZ,1),    & ! Input, surface reflectivity
                 SfcOptics%Index_Sat_Ang,                  & ! Input, Satellite angle index
                 RTV,                                      & ! Input, structure containing forward results
                 Planck_Atmosphere_TL,                     & ! Input, TL layer radiances
                 Planck_Surface_TL,                        & ! Input, TL surface radiance
                 AtmOptics_TL%Single_Scatter_Albedo,       & ! Input, TL layer single scattering albedo
                 AtmOptics_TL%Optical_Depth,               & ! Input, TL layer optical depth
                 SfcOptics_TL%Emissivity(1:nZ,1),          & ! Input, TL surface emissivity
                 SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1), & ! Input, TL surface reflectivity
                 Pff_TL(1:nZ,1:nZ,:),                      & ! Input, TL layer forward phase matrix
                 Pbb_TL(1:nZ,1:nZ,:),                      & ! Input, TL layer backward phase matrix
                 Scattering_Radiance_TL(1:nZ)              ) ! Output, TL radiances
        CASE DEFAULT
          Error_Status = FAILURE
          WRITE(Message,'("Incorrect TL RT_Algorithm_ID, ",i0,", do not fit model")') &
                        RTV%RT_Algorithm_Id
          CALL Display_Message( ROUTINE_NAME,Message,Error_Status )
        RETURN      
      END SELECT
    ELSE
      ! -----------------
      ! Emission model RT
      ! -----------------
      CALL CRTM_Emission_TL( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             RTV%n_Angles,                             & ! Input, number of discrete zenith angles
             GeometryInfo%Cosine_Sensor_Zenith,        & ! Input, cosine of sensor zenith angle
             RTV%Planck_Atmosphere,                    & ! Input, FWD layer radiances
             RTV%Planck_Surface,                       & ! Input, FWD surface radiance
             SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
             SfcOptics%Reflectivity(1:nZ,1,1:nZ,1),    & ! Input, FWD surface reflectivity
             SfcOptics%Direct_Reflectivity(1:nZ,1),    & ! Input, FWD surface reflectivity for a point source
             RTV%Solar_Irradiance,                     & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                     & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,        & ! Input, Source zenith angle
             RTV,                                      & ! Input, internal variables
             AtmOptics_TL%Optical_Depth,               & ! Input, TL layer optical depth
             Planck_Atmosphere_TL,                     & ! Input, TL layer radiances
             Planck_Surface_TL,                        & ! Input, TL surface radiance
             SfcOptics_TL%Emissivity(1:nZ,1),          & ! Input, TL surface emissivity
             SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1), & ! Input, TL surface reflectivity
             SfcOptics_TL%Direct_Reflectivity(1:nZ,1), & ! Input, TL surface reflectivity for a point source
             Radiance_TL                               ) ! Output, TL radiances
    END IF

    Error_Status = Assign_Common_Output_TL( SfcOptics             , &
                                            RTSolution            , &
                                            GeometryInfo          , &
                                            Radiance_TL           , &
                                            Scattering_Radiance_TL, &
                                            SensorIndex           , &
                                            ChannelIndex          , &
                                            RTV                   , &
                                            RTSolution_TL           )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning output for TL RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF

  END FUNCTION CRTM_Compute_RTSolution_TL
!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_RTSolution_AD
!
! PURPOSE:
!       Function to solve the adjoint radiative transfer equation.
!
! CALLING SEQUENCE:
!      Error_Status = CRTM_Compute_RTSolution_AD( Atmosphere   , &  ! FWD Input
!                                                 Surface      , &  ! FWD Input
!                                                 AtmOptics    , &  ! FWD Input
!                                                 SfcOptics    , &  ! FWD Input
!                                                 RTSolution   , &  ! FWD Input
!                                                 RTSolution_AD, &  ! AD Input
!                                                 GeometryInfo , &  ! Input
!                                                 SensorIndex  , &  ! Input
!                                                 ChannelIndex , &  ! Input
!                                                 Atmosphere_AD, &  ! AD Output
!                                                 Surface_AD   , &  ! AD Output
!                                                 AtmOptics_AD , &  ! AD Output
!                                                 SfcOptics_AD , &  ! AD Output
!                                                 RTV            )  ! Internal variable input
!
! INPUT ARGUMENTS:
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Surface:        Structure containing the surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:      Structure containing the surface optical properties
!                       data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution:     Structure containing the solution to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution_AD:  Structure containing the RT solution adjoint inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!                       See the SensorIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTV:            Structure containing internal forward model variables
!                       required for subsequent tangent-linear or adjoint model
!                       calls. The contents of this structure are NOT accessible
!                       outside of the CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Atmosphere_AD:  Structure containing the adjoint atmospheric
!                       state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       Surface_AD:     Structure containing the adjoint surface state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Surface_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       AtmOptics_AD:   Structure containing the adjoint combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       SfcOptics_AD:   Structure containing the adjoint surface optical
!                       properties data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_SfcOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT( IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the computation was sucessful
!                          == FAILURE an unrecoverable error occurred
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on all of the adjoint arguments (whether input or output)
!       is IN OUT rather than just OUT. This is necessary because the Input
!       adjoint arguments are modified, and the Output adjoint arguments must
!       be defined prior to entry to this routine. So, anytime a structure is
!       to be output, to prevent memory leaks the IN OUT INTENT is a must.
!
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_RTSolution_AD( &
    Atmosphere   , &  ! FWD Input
    Surface      , &  ! FWD Input
    AtmOptics    , &  ! FWD Input
    SfcOptics    , &  ! FWD Input
    RTSolution   , &  ! FWD Input
    RTSolution_AD, &  ! AD Input
    GeometryInfo , &  ! Input
    SensorIndex  , &  ! Input
    ChannelIndex , &  ! Input
    Atmosphere_AD, &  ! AD Output
    Surface_AD   , &  ! AD Output
    AtmOptics_AD , &  ! AD Output
    SfcOptics_AD , &  ! AD Output
    RTV          ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_RTSolution_type),   INTENT(IN)     :: RTSolution
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Atmosphere_type),   INTENT(IN OUT) :: Atmosphere_AD
    TYPE(CRTM_Surface_type),      INTENT(IN OUT) :: Surface_AD
    TYPE(CRTM_AtmOptics_type),    INTENT(IN OUT) :: AtmOptics_AD
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    TYPE(RTV_type),               INTENT(IN)     :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_RTSolution_AD'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: nz, n1
    REAL(fp)                                     :: Planck_Surface_AD    ! Surface AD radiance
    REAL(fp), DIMENSION( 0:Atmosphere%n_Layers ) :: Planck_Atmosphere_AD ! *LAYER* AD radiances
    REAL(fp) :: User_Emissivity_AD    ! Temporary adjoint variable for SfcOptics calcs.
    ! The following variables are RT model specific
    REAL(fp), DIMENSION( RTV%n_Angles * RTV%n_Stokes, &
                         (RTV%n_Angles+1) * RTV%n_Stokes, &
                         Atmosphere%n_Layers ) :: Pff_AD ! Forward scattering AD phase matrix
    REAL(fp), DIMENSION( RTV%n_Angles * RTV%n_Stokes, &
                         (RTV%n_Angles+1) * RTV%n_Stokes, &
                         Atmosphere%n_Layers ) :: Pbb_AD ! Backward scattering AD phase matrix
    REAL (fp),DIMENSION( RTV%n_Angles * RTV%n_Stokes ) :: Scattering_Radiance_AD
    REAL (fp) :: Radiance_AD(MAX_N_STOKES)


    ! -----
    ! Setup
    ! -----
    Error_Status = SUCCESS

    Pff_AD = ZERO
    Pbb_AD = ZERO


    RTSolution_AD%RT_Algorithm_Name = RTSolution%RT_Algorithm_Name

    Error_Status = Assign_Common_Input_AD( SfcOptics            , &  ! FWD Input
                                           RTSolution           , &  ! FWD Input
                                           GeometryInfo         , &  ! Input
                                           SensorIndex          , &  ! Input
                                           ChannelIndex         , &  ! Input
                                           RTSolution_AD        , &  ! AD Output/Input
                                           SfcOptics_AD         , &  ! AD Output
                                           Planck_Surface_AD    , &  ! AD Output
                                           Planck_Atmosphere_AD , &  ! AD Output
                                           Radiance_AD          , &  ! AD Output
                                           nz                   , &  ! Output
                                           RTV                    )  ! Internal variable input
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning input for AD RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF

    ! --------------------------------------
    ! Perform the adjoint radiative transfer
    ! --------------------------------------
    ! Select the RT model
    nZ = RTV%n_Angles * RTV%n_Stokes
    n1 = (SfcOptics%Index_Sat_Ang-1)*RTV%n_Stokes + 1
        
    IF( RTV%Scattering_RT ) THEN

      ! Initialise the input adjoint radiance
      Scattering_Radiance_AD = ZERO
      Scattering_Radiance_AD(n1:n1-1+RTV%n_Stokes) = Radiance_AD(1:RTV%n_Stokes)    ! qliu 1/11
    END IF

    IF( RTV%n_Stokes > 1 ) THEN
        ! Select the RT model
      IF( RTV%Scattering_RT ) THEN
        CALL CRTM_ADA_AD( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
             AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
             RTV%Cosmic_Background_Radiance,           & ! Input, cosmic background radiation
             SfcOptics%S_Emissivity(1:nZ),                   & ! Input, FWD surface emissivity
             SfcOptics%S_Direct_Ref(1:nZ),    & ! Input, FWD surface reflectivity for a point source
             RTV,                                      & ! In/Output, internal variables
             Scattering_Radiance_AD(1:nZ),             & ! Input, AD radiances
             Planck_Atmosphere_AD,                     & ! Output, AD layer radiances
             Planck_Surface_AD,                        & ! Output, AD surface radiance
             AtmOptics_AD%Single_Scatter_Albedo,       & ! Output, AD layer single scattering albedo
             AtmOptics_AD%Optical_Depth,               & ! Output, AD layer optical depth
             SfcOptics_AD%S_Emissivity(1:nZ),          & ! Output, AD surface emissivity
             SfcOptics_AD%S_Reflectivity(1:nZ,1:nZ), & ! Output, AD surface reflectivity
             SfcOptics_AD%S_Direct_Ref(1:nZ), & ! Output, AD surface reflectivity for a point source
             Pff_AD(1:nZ,1:(nZ+1),:),                  & ! Output, AD layer forward phase matrix
             Pbb_AD(1:nZ,1:(nZ+1),:)                   ) ! Output, AD layer backward phase matrix      
      ELSE
        CALL CRTM_Emission_AD( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             RTV%n_Angles,                             & ! Input, number of discrete zenith angles
             GeometryInfo%Cosine_Sensor_Zenith,        & ! Input, cosine of sensor zenith angle
             RTV%Planck_Atmosphere,                    & ! Input, FWD layer radiances
             RTV%Planck_Surface,                       & ! Input, FWD surface radiance
             SfcOptics%S_Emissivity(1:nZ),             & ! Input, FWD surface emissivity
             SfcOptics%S_Reflectivity(1:nZ,1:nZ),    & ! Input, FWD surface reflectivity
             SfcOptics%S_Direct_Ref(1:nZ),    & ! Input, FWD surface reflectivity for a point source
             RTV%Solar_Irradiance,                     & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                     & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,        & ! Input, Source zenith angle
             RTV,                                      & ! Input, internal variables
             Radiance_AD(1),                         & ! Input, AD radiance
             AtmOptics_AD%Optical_Depth,               & ! Output, AD layer optical depth
             Planck_Atmosphere_AD,                     & ! Output, AD layer radiances
             Planck_Surface_AD,                        & ! Output, AD surface radiance
             SfcOptics_AD%S_Emissivity(1:nZ),          & ! Output, AD surface emissivity
             SfcOptics_AD%S_Reflectivity(1:nZ,1:nZ), & ! Output, AD surface reflectivity
             SfcOptics_AD%S_Direct_Ref(1:nZ)  ) ! Output, AD surface reflectivity for a point source      
      END IF
      CALL Reshape_Surf_Opt_AD(RTV%n_Angles, RTV%n_Stokes, SfcOptics_AD%Emissivity, SfcOptics_AD%Direct_Reflectivity, &
        SfcOptics_AD%Reflectivity, SfcOptics_AD%S_Emissivity, SfcOptics_AD%S_Direct_Ref, SfcOptics_AD%S_Reflectivity)
!!      Scattering_Radiance_AD( SfcOptics%Index_Sat_Ang ) = Radiance_AD(1)    ! qliu 1/11
!!      RTSolution_AD%Radiance = ZERO
   ELSE IF( RTV%Scattering_RT ) THEN
      ! Select the scattering RT model
      SELECT CASE(RTV%RT_Algorithm_Id)

        CASE (RT_ADA, RT_VMOM)
        ! NESDIS advanced adding-doubling method
        CALL CRTM_ADA_AD( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
             AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
             RTV%Cosmic_Background_Radiance,           & ! Input, cosmic background radiation
             SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
             SfcOptics%Direct_Reflectivity(1:nZ,1),    & ! Input, FWD surface reflectivity for a point source
             RTV,                                      & ! In/Output, internal variables
             Scattering_Radiance_AD(1:nZ),             & ! Input, AD radiances
             Planck_Atmosphere_AD,                     & ! Output, AD layer radiances
             Planck_Surface_AD,                        & ! Output, AD surface radiance
             AtmOptics_AD%Single_Scatter_Albedo,       & ! Output, AD layer single scattering albedo
             AtmOptics_AD%Optical_Depth,               & ! Output, AD layer optical depth
             SfcOptics_AD%Emissivity(1:nZ,1),          & ! Output, AD surface emissivity
             SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1), & ! Output, AD surface reflectivity
             SfcOptics_AD%Direct_Reflectivity(1:nZ,1), & ! Output, AD surface reflectivity for a point source
             Pff_AD(1:nZ,1:(nZ+1),:),                  & ! Output, AD layer forward phase matrix
             Pbb_AD(1:nZ,1:(nZ+1),:)                   ) ! Output, AD layer backward phase matrix

        CASE (RT_SOI)
        ! UW SOI RT solver
        CALL CRTM_SOI_AD( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             AtmOptics%Single_Scatter_Albedo,          & ! Input, FWD layer single scattering albedo
             AtmOptics%Optical_Depth,                  & ! Input, FWD layer optical depth
             SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
             SfcOptics%Reflectivity(1:nZ,1,1:nZ,1),    & ! Input, FWD surface reflectivity
             SfcOptics%Index_Sat_Ang,                  & ! Input, Satellite angle index
             RTV,                                      & ! In/Output, internal variables
             Scattering_Radiance_AD(1:nZ),             & ! Input, AD radiances
             Planck_Atmosphere_AD,                     & ! Output AD atmospheric layer Planck radiance
             Planck_Surface_AD,                        & ! Output AD surface Planck radiance
             AtmOptics_AD%Single_Scatter_Albedo,       & ! Output, AD layer single scattering albedo
             AtmOptics_AD%Optical_Depth,               & ! Output, AD layer optical depth
             SfcOptics_AD%Emissivity(1:nZ,1),          & ! Output, AD surface emissivity
             SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1), & ! Output, AD surface reflectivity
             Pff_AD(1:nZ,1:(nZ+1),:),                  & ! Output, AD layer forward phase matrix
             Pbb_AD(1:nZ,1:(nZ+1),:)                   ) ! Output, AD layer backward phase matrix
      CASE DEFAULT
      Error_Status = FAILURE
      WRITE(Message,'("Incorrect AD RT_Algorithm_ID, ",i0,", do not fit model")') &
                    RTV%RT_Algorithm_Id
      CALL Display_Message( ROUTINE_NAME,Message,Error_Status )
      RETURN      
      END SELECT

    ELSE

      ! -----------------
      ! Emission model RT
      ! -----------------
      CALL CRTM_Emission_AD( &
             Atmosphere%n_Layers,                      & ! Input, number of atmospheric layers
             RTV%n_Angles,                             & ! Input, number of discrete zenith angles
             GeometryInfo%Cosine_Sensor_Zenith,        & ! Input, cosine of sensor zenith angle
             RTV%Planck_Atmosphere,                    & ! Input, FWD layer radiances
             RTV%Planck_Surface,                       & ! Input, FWD surface radiance
             SfcOptics%Emissivity(1:nZ,1),             & ! Input, FWD surface emissivity
             SfcOptics%Reflectivity(1:nZ,1,1:nZ,1),    & ! Input, FWD surface reflectivity
             SfcOptics%Direct_Reflectivity(1:nZ,1),    & ! Input, FWD surface reflectivity for a point source
             RTV%Solar_Irradiance,                     & ! Input, Source irradiance at TOA
             RTV%Is_Solar_Channel,                     & ! Input, Source sensitive channel info.
             GeometryInfo%Source_Zenith_Radian,        & ! Input, Source zenith angle
             RTV,                                      & ! Input, internal variables
             Radiance_AD(1),                         & ! Input, AD radiance
             AtmOptics_AD%Optical_Depth,               & ! Output, AD layer optical depth
             Planck_Atmosphere_AD,                     & ! Output, AD layer radiances
             Planck_Surface_AD,                        & ! Output, AD surface radiance
             SfcOptics_AD%Emissivity(1:nZ,1),          & ! Output, AD surface emissivity
             SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1), & ! Output, AD surface reflectivity
             SfcOptics_AD%Direct_Reflectivity(1:nZ,1)  ) ! Output, AD surface reflectivity for a point source
    END IF

    Error_Status = Assign_Common_Output_AD( Atmosphere           , & ! Input
                                            Surface              , & ! Input
                                            AtmOptics            , & ! Input
                                            SfcOptics            , & ! Input
                                            Pff_AD               , & ! Input
                                            Pbb_AD               , & ! Input
                                            GeometryInfo         , & ! Input
                                            SensorIndex          , & ! Input
                                            ChannelIndex         , & ! Input
                                            nZ                   , & ! Input
                                            AtmOptics_AD         , & ! Output
                                            SfcOptics_AD         , & ! Output
                                            Planck_Surface_AD    , & ! Output
                                            Planck_Atmosphere_AD , & ! Output
                                            User_Emissivity_AD   , & ! Output
                                            Atmosphere_AD        , & ! Output
                                            Surface_AD           , & ! Output
                                            RTSolution_AD        , & ! Output
                                            RTV                  )   ! Input
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Error assigning output for AD RTSolution algorithms'
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
      RETURN
    END IF


  END FUNCTION CRTM_Compute_RTSolution_AD

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_n_Streams
!
! PURPOSE:
!       Function to compute the number of streams required for subsequent
!       radiative transfer calculations
!
! CALLING SEQUENCE:
!       nStreams = CRTM_Compute_n_Streams( Atmosphere,   &  ! Input
!                                          SensorIndex,  &  ! Input
!                                          ChannelIndex, &  ! Input
!                                          RTSolution    )  ! Output
!
! INPUT ARGUMENTS:
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id. This is a unique index associated
!                       with a (supported) sensor used to access the
!                       shared coefficient data for a particular sensor.
!                       See the ChannelIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id. This is a unique index associated
!                       with a (supported) sensor channel used to access the
!                       shared coefficient data for a particular sensor's
!                       channel.
!                       See the SensorIndex argument.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution:     Structure containing the scattering flag to be set
!                       for the RT calcs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       nStreams:       The number of RT streams required to perform radiative
!                       transfer in a scattering atmosphere.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!--------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_nStreams( &
    Atmosphere  , &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    RTSolution  ) &  ! Output
  RESULT( nStreams )
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atmosphere
    INTEGER,                    INTENT(IN)     :: SensorIndex
    INTEGER,                    INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTSolution
    ! Function result
    INTEGER :: nStreams
    ! Local variables
    REAL(fp) :: maxReff, Reff, MieParameter
    INTEGER :: n

    ! Set up
    nStreams = 0
    RTSolution%n_full_Streams = nStreams
    RTSolution%Scattering_FLAG = .FALSE.

    ! If no clouds and no aerosols, no scattering, so return
    IF ( Atmosphere%n_Clouds   == 0 .AND. &
         Atmosphere%n_Aerosols == 0       ) RETURN

    ! Determine the maximum cloud particle size
    maxReff = ZERO
    DO n = 1, Atmosphere%n_Clouds
      Reff = MAXVAL(Atmosphere%Cloud(n)%Effective_Radius)
      IF( Reff > maxReff) maxReff = Reff
    END DO
    DO n = 1, Atmosphere%n_Aerosols
      Reff = MAXVAL(Atmosphere%Aerosol(n)%Effective_Radius)
      IF( Reff > maxReff) maxReff = Reff
    END DO

    ! Compute the Mie parameter, 2.pi.Reff/lambda
    MieParameter = TWO * PI * maxReff * SC(SensorIndex)%Wavenumber(ChannelIndex)/10000.0_fp

    ! Determine the number of streams based on Mie parameter
    IF ( MieParameter < 0.01_fp ) THEN
      nStreams = 2
    ELSE IF( MieParameter < ONE ) THEN
      nStreams = 4
    ELSE
      nStreams = 6
    END IF

! Hardcode number of streams for testing purposes
!    nStreams = 6

    ! Set RTSolution scattering info
    RTSolution%Scattering_Flag = .TRUE.
    RTSolution%n_full_Streams  = nStreams + 2

  END FUNCTION CRTM_Compute_nStreams
  
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_RTSolution_Version
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CRTM_RTSolution_Version( Id )
!
! OUTPUT ARGUMENTS:
!       Id:            Character string containing the version Id information
!                      for the module.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_RTSolution_Version( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_RTSolution_Version
  
END MODULE CRTM_RTSolution
