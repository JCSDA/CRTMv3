!
! CRTM_SfcOptics
!
! Module to compute the surface optical properties required for
! determining the surface contribution to the radiative transfer.
!
!
! CREATION HISTORY:
!       Written by:     Yong Han,       NOAA/NESDIS;     Yong.Han@noaa.gov
!                       Quanhua Liu,    QSS Group, Inc;  Quanhua.Liu@noaa.gov
!                       Paul van Delst, CIMSS/SSEC;      paul.vandelst@ssec.wisc.edu
!                       02-Apr-2004
!
! It has been modified replace the CRTM surface modules with the
! stand-along, expanded and integrated surface RT modling system
! CSEM (Community Surface Emissivity Models)
!
! Since the CRTM surface data structures haven't been modified for
! the new surface modeling system,"adptor" subroutines are used to
! convert the surface I/O streams between CRTM and CSEM.
!
!
! CREATION HISTORY:
!       Written by:     Yong Han,       NOAA/NESDIS;     Yong.Han@noaa.gov
!                       Quanhua Liu,    QSS Group, Inc;  Quanhua.Liu@noaa.gov
!                       Paul van Delst, CIMSS/SSEC;      paul.vandelst@ssec.wisc.edu
!                       02-Apr-2004
!
!       Modifed by:     Ming Chen,       UMD-CICS;     ming.chen@noaa.gov
!                       08-12-2018
! initial commit to jcsda/crtm repository  01/26/2021
! Patrick Stegmann 2021-01-22     Added CONST_MIXED_POLARIZATION scheme.
!
! Patrick Stegmann 2021-08-31     Added PRA_POLARIZATION scheme for GEMS-1.
!


MODULE CRTM_SfcOptics

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE CSEM_Type_Kinds,               ONLY: fp => CSEM_fp
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,          ONLY: ZERO, POINT_5, ONE, DEGREES_TO_RADIANS, MAX_N_STOKES
  USE CRTM_SpcCoeff,            ONLY: SC, &
                                      SpcCoeff_IsSolar,                   &
                                      SpcCoeff_IsMicrowaveSensor,         &
                                      SpcCoeff_IsInfraredSensor,          &
                                      SpcCoeff_IsVisibleSensor,           &
                                      SpcCoeff_IsUltravioletSensor,       &
                                      UNPOLARIZED,                        &
                                      INTENSITY,                          &
                                      FIRST_STOKES_COMPONENT,             &
                                      SECOND_STOKES_COMPONENT,            &
                                      THIRD_STOKES_COMPONENT,             &
                                      FOURTH_STOKES_COMPONENT,            &
                                      VL_POLARIZATION,                    &
                                      HL_POLARIZATION,                    &
                                      plus45L_POLARIZATION,               &
                                      minus45L_POLARIZATION,              &
                                      VL_MIXED_POLARIZATION,              &
                                      HL_MIXED_POLARIZATION,              &
                                      RC_POLARIZATION,                    &
                                      LC_POLARIZATION,                    &
                                      CONST_MIXED_POLARIZATION,           &
                                      PRA_POLARIZATION
 
  USE CRTM_Surface_Define,      ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE CRTM_SfcOptics_Define,    ONLY: CRTM_SfcOptics_type,                &
                                      OPERATOR(==),                       &
                                      CRTM_SfcOptics_Associated,          &
                                      CRTM_SfcOptics_Destroy   ,          &
                                      CRTM_SfcOptics_Create
      
  USE CSEM_LandMW_SfcOptics,    ONLY: CSEM_MWLSOVar_type => iVar_type,    &
                                      CSEM_Compute_LandMW_SfcOptics,      &
                                      CSEM_Compute_LandMW_SfcOptics_TL,   &
                                      CSEM_Compute_LandMW_SfcOptics_AD
  USE CSEM_WaterMW_SfcOptics,   ONLY: CSEM_MWWSOVar_type => iVar_type,    &
                                      CSEM_Compute_WaterMW_SfcOptics,     &
                                      CSEM_Compute_WaterMW_SfcOptics_TL,  &
                                      CSEM_Compute_WaterMW_SfcOptics_AD
  USE CSEM_SnowMW_SfcOptics,    ONLY: CSEM_MWSSOVar_type => iVar_type,    &
                                      CSEM_Compute_SnowMW_SfcOptics,      &
                                      CSEM_Compute_SnowMW_SfcOptics_TL,   &
                                      CSEM_Compute_SnowMW_SfcOptics_AD
  USE CSEM_IceMW_SfcOptics,     ONLY: CSEM_MWISOVar_type => iVar_type,    &
                                      CSEM_Compute_IceMW_SfcOptics,       &
                                      CSEM_Compute_IceMW_SfcOptics_TL,    &
                                      CSEM_Compute_IceMW_SfcOptics_AD

  USE CSEM_LandIR_SfcOptics,    ONLY: CSEM_IRLSOVar_type => iVar_type,    &
                                      CSEM_Compute_LandIR_SfcOptics,      &
                                      CSEM_Compute_LandIR_SfcOptics_TL,   &
                                      CSEM_Compute_LandIR_SfcOptics_AD
  USE CSEM_WaterIR_SfcOptics,   ONLY: CSEM_IRWSOVar_type => iVar_type,    &
                                      CSEM_Compute_WaterIR_SfcOptics,     &
                                      CSEM_Compute_WaterIR_SfcOptics_TL,  &
                                      CSEM_Compute_WaterIR_SfcOptics_AD
  USE CSEM_SnowIR_SfcOptics,    ONLY: CSEM_IRSSOVar_type => iVar_type,    &
                                      CSEM_Compute_SnowIR_SfcOptics,      &
                                      CSEM_Compute_SnowIR_SfcOptics_TL,   &
                                      CSEM_Compute_SnowIR_SfcOptics_AD
  USE CSEM_IceIR_SfcOptics,     ONLY: CSEM_IRISOVar_type => iVar_type,    &
                                      CSEM_Compute_IceIR_SfcOptics,       &
                                      CSEM_Compute_IceIR_SfcOptics_TL,    &
                                      CSEM_Compute_IceIR_SfcOptics_AD
  USE CSEM_LandVIS_SfcOptics,   ONLY: CSEM_VISLSOVar_type => iVar_type,   &
                                      CSEM_Compute_LandVIS_SfcOptics,     &
                                      CSEM_Compute_LandVIS_SfcOptics_TL,  &
                                      CSEM_Compute_LandVIS_SfcOptics_AD
  USE CSEM_WaterVIS_SfcOptics,  ONLY: CSEM_VISWSOVar_type => iVar_type,   &
                                      CSEM_Compute_WaterVIS_SfcOptics,    &
                                      CSEM_Compute_WaterVIS_SfcOptics_TL, &
                                      CSEM_Compute_WaterVIS_SfcOptics_AD
  USE CSEM_SnowVIS_SfcOptics,   ONLY: CSEM_VISSSOVar_type => iVar_type,   &
                                      CSEM_Compute_SnowVIS_SfcOptics,     &
                                      CSEM_Compute_SnowVIS_SfcOptics_TL,  &
                                      CSEM_Compute_SnowVIS_SfcOptics_AD
  USE CSEM_IceVIS_SfcOptics,    ONLY: CSEM_VISISOVar_type => iVar_type,   &
                                      CSEM_Compute_IceVIS_SfcOptics,      &
                                      CSEM_Compute_IceVIS_SfcOptics_TL,   &
                                      CSEM_Compute_IceVIS_SfcOptics_AD
     
 
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type,             &
                                      CRTM_GeometryInfo_GetValue
      
  USE CSEM_Define
  USE CSEM_Model_Manager
     
     
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  PUBLIC :: iVar_type
  ! Procedures
  PUBLIC :: CRTM_Compute_SurfaceT
  PUBLIC :: CRTM_Compute_SurfaceT_TL
  PUBLIC :: CRTM_Compute_SurfaceT_AD
  PUBLIC :: CRTM_Compute_SfcOptics
  PUBLIC :: CRTM_Compute_SfcOptics_TL
  PUBLIC :: CRTM_Compute_SfcOptics_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_SfcOptics.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov &'
  ! Message length
  INTEGER, PARAMETER :: ML = 256

  LOGICAL :: COM_CRTM_CSEM = .FALSE.

  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Microwave
    TYPE(CSEM_MWLSOVar_type)  ::  CSEM_MWLSOV  ! Land
    TYPE(CSEM_MWWSOVar_type)  ::  CSEM_MWWSOV  ! Water
    TYPE(CSEM_MWSSOVar_type)  ::  CSEM_MWSSOV  ! Snow
    TYPE(CSEM_MWISOVar_type)  ::  CSEM_MWISOV  ! Ice
    ! Infrared
    TYPE(CSEM_IRLSOVar_type)  ::  CSEM_IRLSOV  ! Land
    TYPE(CSEM_IRWSOVar_type)  ::  CSEM_IRWSOV  ! Water
    TYPE(CSEM_IRSSOVar_type)  ::  CSEM_IRSSOV  ! Snow
    TYPE(CSEM_IRISOVar_type)  ::  CSEM_IRISOV  ! Ice
    ! Visible
    TYPE(CSEM_VISLSOVar_type) ::  CSEM_VISLSOV ! Land
    TYPE(CSEM_VISWSOVar_type) ::  CSEM_VISWSOV ! Water
    TYPE(CSEM_VISSSOVar_type) ::  CSEM_VISSSOV ! Snow
    TYPE(CSEM_VISISOVar_type) ::  CSEM_VISISOV ! Ice
  END TYPE iVar_type

  TYPE(CSEM_Model_ID)  :: CSEM_Model
 
CONTAINS


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_SurfaceT
!
! PURPOSE:
!       Subroutine to compute the average of the various surface type
!       temperatures weighted by their coverage fraction.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_SurfaceT( Surface,  &  ! Input
!                                   SfcOptics )  ! Output
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        temperature required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_SurfaceT( Surface,  &  ! Input
                                    SfcOptics )  ! Output
    ! Arguments
    TYPE(CRTM_Surface_type),   INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type), INTENT(IN OUT) :: SfcOptics

    ! The weighted average surface temperature
    SfcOptics%Surface_Temperature = &
      ( Surface%Land_Coverage  * Surface%Land_Temperature  ) + &
      ( Surface%Water_Coverage * Surface%Water_Temperature ) + &
      ( Surface%Snow_Coverage  * Surface%Snow_Temperature  ) + &
      ( Surface%Ice_Coverage   * Surface%Ice_Temperature   )

  END SUBROUTINE CRTM_Compute_SurfaceT


!----------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_SurfaceT_TL
!
! PURPOSE:
!       Subroutine to compute the tangent-linear average of the various
!       surface type temperatures weighted by their coverage fraction.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_SurfaceT_TL( Surface,     &  ! Input
!                                      Surface_TL,  &  ! Input
!                                      SfcOptics_TL )  ! In/Output
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:      CRTM_Surface structure containing the tangent-linerar
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface temperature required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!----------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_SurfaceT_TL( Surface,     &  ! Input
                                       Surface_TL,  &  ! Input
                                       SfcOptics_TL )  ! Output
    ! Arguments
    TYPE(CRTM_Surface_type),   INTENT(IN)     :: Surface
    TYPE(CRTM_Surface_type),   INTENT(IN)     :: Surface_TL
    TYPE(CRTM_SfcOptics_type), INTENT(IN OUT) :: SfcOptics_TL

    ! The weighted average tangent-linear surface temperature
    SfcOptics_TL%Surface_Temperature = &
      ( Surface%Land_Coverage  * Surface_TL%Land_Temperature  ) + &
      ( Surface%Water_Coverage * Surface_TL%Water_Temperature ) + &
      ( Surface%Snow_Coverage  * Surface_TL%Snow_Temperature  ) + &
      ( Surface%Ice_Coverage   * Surface_TL%Ice_Temperature   )

  END SUBROUTINE CRTM_Compute_SurfaceT_TL


!----------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_SurfaceT_AD
!
! PURPOSE:
!       Subroutine to compute the adjoint of the average of the various
!       surface type temperatures weighted by their coverage fraction.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_SurfaceT_AD( Surface,      &  ! Input
!                                      SfcOptics_AD, &  ! Input
!                                      Surface_AD    )  ! Output
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_AD:    CRTM_SfcOptics structure containing the adjoint
!                        surface temperature required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUTS:
!       Surface_AD:      CRTM_Surface structure containing the adjoint surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! SIDE EFFECTS:
!       Even though the SfcOptics_AD argument is listed as an INPUT, its
!       INTENT is ( IN OUT ) as it is modified on output since the
!       Surface_Temperature component is set to zero after the adjoint
!       calculation.
!
!       Even though the Surface_AD argument is listed as an OUTPUT, its
!       INTENT is ( IN OUT ) as the components of the adjoint calculation
!       in this routine may already have a value from a previous adjoint
!       calculation performed on the structure.
!
! COMMENTS:
!       In addition to the input/output requirements described in the SIDE
!       EFFECTS section, the SfcOptics_AD and Surface_AD arguments require
!       an INTENT of IN OUT to prevent memory leaks.
!
!----------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_SurfaceT_AD( Surface,      &  ! Input
                                       SfcOptics_AD, &  ! Input
                                       Surface_AD    )  ! Output
    ! Arguments
    TYPE(CRTM_Surface_type),   INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type), INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_Surface_type),   INTENT(IN OUT) :: Surface_AD

    ! The adjoint of the weighted average surface temperature
    Surface_AD%Land_Temperature  = Surface_AD%Land_Temperature + &
                                   (Surface%Land_Coverage *SfcOptics_AD%Surface_Temperature)
    Surface_AD%Water_Temperature = Surface_AD%Water_Temperature + &
                                   (Surface%Water_Coverage*SfcOptics_AD%Surface_Temperature)
    Surface_AD%Snow_Temperature  = Surface_AD%Snow_Temperature  + &
                                   (Surface%Snow_Coverage *SfcOptics_AD%Surface_Temperature)
    Surface_AD%Ice_Temperature   = Surface_AD%Ice_Temperature   + &
                                   (Surface%Ice_Coverage  *SfcOptics_AD%Surface_Temperature)
    SfcOptics_AD%Surface_Temperature = ZERO

  END SUBROUTINE CRTM_Compute_SurfaceT_AD


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_SfcOptics
!
! PURPOSE:
!       Function to compute the surface optical properties and populate
!       the output SfcOptics structure for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_SfcOptics( &
!                        Surface     , &  ! Input
!                        GeometryInfo, &  ! Input
!                        SensorIndex , &  ! Input
!                        ChannelIndex, &  ! Input
!                        SfcOptics   , &  ! Output
!                        iVar          )  ! Internal variable output
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       CRTM_GeometryInfo_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        On Input:  The Secant_Angle component is assumed to
!                                   contain data.
!                        On Output: The Emissivity and Reflectivity components
!                                   will contain the required data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT. This is necessary because the argument should be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_SfcOptics( &
    Surface     , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   , &  ! Output
    iVar        ) &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type)     , INTENT(IN)     :: Surface
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN OUT) :: SfcOptics
    TYPE(iVar_type)             , INTENT(OUT)    :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_SfcOptics'
    ! Local variables
    CHARACTER(ML) :: Message
    INTEGER :: i
    INTEGER :: nL, nZ
    REAL(fp) :: SIN2_Angle
    REAL(fp) :: pv, ph
    REAL(fp) :: phi, theta_f
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Emissivity
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES, &
                        SfcOptics%n_Angles,MAX_N_STOKES) :: Reflectivity
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Direct_Reflectivity
    INTEGER :: Polarization

    TYPE(CSEM_Land_Surface)  :: CSEM_Land
    TYPE(CSEM_Snow_Surface)  :: CSEM_Snow
    TYPE(CSEM_Water_Surface) :: CSEM_Water
    TYPE(CSEM_Ice_Surface)   :: CSEM_Ice
    TYPE(CSEM_Options_Type)  :: CSEM_Options

    TYPE(CSEM_SfcOptics_Type):: CSEM_SfcOptics
   
    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    nL = SfcOptics%n_Stokes
    nZ = SfcOptics%n_Angles
    Polarization = SC(SensorIndex)%Polarization(ChannelIndex)
    ! Initialise the local emissivity and reflectivities
    Emissivity   = ZERO
    Reflectivity = ZERO
    Direct_Reflectivity = ZERO

    Error_Status = CRTM_CSEM_Input_Adaptor( &
                   Surface,     &  ! Input
                   CSEM_Land,   &  ! Output
                   CSEM_Water,  &  ! Output
                   CSEM_Snow,   &  ! Output
                   CSEM_Ice)    

    ! Data stream from CRTM to CSEM
    CALL CSEM_Options%SensorObs%init(size(Surface%SensorData%Tb))
    CALL CSEM_SfcOptics%init(N_Angles=nZ)

    CSEM_SfcOptics%Frequency     =    SC(SensorIndex)%Frequency(ChannelIndex)
    CSEM_SfcOptics%Wavenumber    =    SC(SensorIndex)%Wavenumber(ChannelIndex)
    CSEM_SfcOptics%Is_Solar      =    SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex)
    CSEM_SfcOptics%Angle(1:nZ)   =    SfcOptics%Angle(1:nZ)
    CSEM_SfcOptics%Weight(1:nZ)  =    SfcOptics%Weight(1:nZ)

    CSEM_Options%SensorObs%Sensor_ID =  TRIM(SC(SensorIndex)%Sensor_ID)

    CALL CRTM_GeometryInfo_GetValue ( &
         GeometryInfo,                                                      &
         Longitude               =    CSEM_Options%GeoInfo%Longitude,       &
         Latitude                =    CSEM_Options%GeoInfo%Latitude,        &
         Month                   =    CSEM_Options%GeoInfo%Month,           &
         Sensor_Scan_Angle       =    CSEM_SfcOptics%Sensor_Scan_Angle,     &
         Sensor_Zenith_Angle     =    CSEM_SfcOptics%Sensor_Zenith_Angle,   &
         Sensor_Azimuth_Angle    =    CSEM_SfcOptics%Sensor_Azimuth_Angle,  &
         Source_Zenith_Angle     =    CSEM_SfcOptics%Source_Zenith_Angle,   &
         Source_Azimuth_Angle    =    CSEM_SfcOptics%Source_Azimuth_Angle)

    IF(CSEM_Options%SensorObs%Is_Allocated) THEN
         CSEM_Options%SensorObs%Tb = Surface%SensorData%Tb
         CSEM_Options%SensorObs%Channel_Frequency=SC(SensorIndex)%Frequency
         CSEM_Options%SensorObs%Channel_Polarization=SC(SensorIndex)%Polarization
         !CALL Set_Ref_Channel(CSEM_Options%SensorObs)
    END IF
    
    CSEM_Options%Atmos%Transmittance = SfcOptics%Transmittance
    !CSEM_Options%Atmos%Transmittance = 0.6_fp
    
    
    
      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                     ## MICROWAVE CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      Sensor_Select: IF ( SpcCoeff_IsMicrowaveSensor( SC(SensorIndex) ) ) THEN

        ! --------------------------------------
        ! Microwave LAND emissivity/reflectivity
        ! --------------------------------------
        Microwave_Land: IF( Surface%Land_Coverage > ZERO) THEN

          ! Compute the surface optics
          CSEM_Model = Inq_Model_Option("MW_LAND") 
          !CSEM_Model%NAME = "TELSEM_ATLAS"
          !CALL SET_CSEM_Model(CSEM_Model)
          Error_Status = CSEM_Compute_LandMW_SfcOptics(  &
                         CSEM_Land,              &  ! Input
                         CSEM_SfcOptics,         &  ! output
                         CSEM_Options,           &  ! input
                         iVar%CSEM_MWLSOV   )       ! Output

          IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Land',CSEM_SfcOptics%Emissivity(1,1:2)
          END IF
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW land SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          ! Accumulate the surface optics properties
          ! based on land coverage fraction
          DO i = 1, nZ
             Emissivity(i,1:2)           =  CSEM_SfcOptics%Emissivity(i,1:2)         *  Surface%Land_Coverage
             Reflectivity(i,1:2,i,1:2)   =  CSEM_SfcOptics%Reflectivity(i,1:2,i,1:2) *  Surface%Land_Coverage
          END DO

        END IF Microwave_Land


        ! ---------------------------------------
        ! Microwave WATER emissivity/reflectivity
        ! ---------------------------------------
        Microwave_Water: IF( Surface%Water_Coverage > ZERO ) THEN
          CSEM_Model = Inq_Model_Option("MW_WATER") 
          !CSEM_Model%NAME = "NESDIS_FASTEM_V6"
          !CALL SET_CSEM_Model(CSEM_Model)

          ! Compute the surface optics
          Error_Status =  CSEM_Compute_WaterMW_SfcOptics(  &
                          CSEM_Water,                      &  ! Input
                          CSEM_SfcOptics,                  &  ! output
                          CSEM_Options,ivar%CSEM_MWWSOV   )    ! input
          IF(COM_CRTM_CSEM) THEN
            print*, TRIM(CSEM_Model%NAME)
            WRITE(*,*)'CSEM_Water',CSEM_SfcOptics%Emissivity(1,1:4)
          END IF
         IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW water SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          ! Accumulate the surface optics properties
          ! based on water coverage fraction
          DO i = 1, nZ
             Emissivity(i,1:4)           =  Emissivity(i,1:4)         +  &
                 CSEM_SfcOptics%Emissivity(i,1:4)         *  Surface%Water_Coverage
             Reflectivity(i,1:4,i,1:4)   =  Reflectivity(i,1:4,i,1:4) +  &
                 CSEM_SfcOptics%Reflectivity(i,1:4,i,1:4) *  Surface%Water_Coverage
          END DO


         END IF Microwave_Water


        ! --------------------------------------
        ! Microwave SNOW emissivity/reflectivity
        ! --------------------------------------
        Microwave_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN
          ! Compute the surface optics
        Error_Status = CSEM_Compute_SnowMW_SfcOptics( &
                       CSEM_Snow,             &  ! Input
                       CSEM_SfcOptics,        &  ! output
                       CSEM_Options   )          ! intput
          IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Snow',CSEM_SfcOptics%Emissivity(1,1:4)
          END IF

          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW snow SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity(i,1:2)           =  Emissivity(i,1:2)         +  &
                 CSEM_SfcOptics%Emissivity(i,1:2)         *  Surface%Snow_Coverage
             Reflectivity(i,1:2,i,1:2)   =  Reflectivity(i,1:2,i,1:2) +  &
                 CSEM_SfcOptics%Reflectivity(i,1:2,i,1:2) *  Surface%Snow_Coverage
          END DO
  
        END IF Microwave_Snow


        ! -------------------------------------
        ! Microwave ICE emissivity/reflectivity
        ! -------------------------------------
        Microwave_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN
         ! Compute the surface optics

         Error_Status = CSEM_Compute_IceMW_SfcOptics(  &
                        CSEM_Ice,              &  ! Input
                        CSEM_SfcOptics,        &  ! Output
                        CSEM_Options   )          ! input
         IF(COM_CRTM_CSEM) THEN
             !WRITE(*,*)'CSEM_Ice',CSEM_Channel%frequency,CSEM_Geometry(1)%zenith_angle
             WRITE(*,*)'CSEM_Ice',CSEM_SfcOptics%Emissivity(i,1:2)
          END IF

          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW ice SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity(i,1:2)           =  Emissivity(i,1:2)         +  &
                 CSEM_SfcOptics%Emissivity(i,1:2)         *  Surface%Ice_Coverage
             Reflectivity(i,1:2,i,1:2)   =  Reflectivity(i,1:2,i,1:2) +  &
                 CSEM_SfcOptics%Reflectivity(i,1:2,i,1:2) *  Surface%Ice_Coverage
          END DO
   
        END IF Microwave_Ice



        !#----------------------------------------------------------------------#
        !#                 -- HANDLE THE DECOUPLED POLARISATION --              #
        !#                                                                      #
        !# The SfcOptics n_Stokes dimension determines whether the surface      #
        !# optics takes into account the second order effect of cross           #
        !# polarisation, e.g. if the surface optics for a purely vertically     #
        !# polarised channel has a horizontal (or other) component due to       #
        !# scattering at the surface.                                           #
        !#                                                                      #
        !# If the SfcOptics n_Stokes dimension == 1, the polarisations are      #
        !# decoupled.                                                           #
        !#----------------------------------------------------------------------#

        Decoupled_Polarization: IF( SfcOptics%n_Stokes == 1 ) THEN


          ! ------------------------------------------------------
          ! Decoupled polarisation. Branch on channel polarisation
          ! ------------------------------------------------------
          Polarization_Type: SELECT CASE( Polarization )

            ! The unpolarised case, I
            ! e = (eV + eH)/2
            ! r = (rV + rH)/2
            ! Note: INTENSITY == UNPOLARIZED == FIRST_STOKES_COMPONENT
            CASE( INTENSITY )
              SfcOptics%Emissivity(1:nZ,1) = &
                POINT_5 * ( Emissivity(1:nZ,1) + Emissivity(1:nZ,2) )
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = &
                POINT_5 * ( Reflectivity(1:nZ,1,1:nZ,1) + Reflectivity(1:nZ,2,1:nZ,2) )

            ! The second Stokes component, Q, the polarisation difference.
            ! e = (eV - eH)/2
            ! r = (rV - rH)/2
            CASE( SECOND_STOKES_COMPONENT )
              SfcOptics%Emissivity(1:nZ,1) = &
                POINT_5 * ( Emissivity(1:nZ,1) - Emissivity(1:nZ,2) )
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = &
                POINT_5 * ( Reflectivity(1:nZ,1,1:nZ,1) - Reflectivity(1:nZ,2,1:nZ,2) )

            ! The third Stokes component, U.
            CASE ( THIRD_STOKES_COMPONENT )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,3)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,3,1:nZ,3)

            ! The fourth Stokes component, V.
            CASE ( FOURTH_STOKES_COMPONENT )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,4)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,4,1:nZ,4)

            ! Vertical linear polarisation
            CASE ( VL_POLARIZATION )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)

            ! Horizontal linear polarisation
            CASE ( HL_POLARIZATION )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,2)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,2,1:nZ,2)

            ! +45deg. linear polarisation
            CASE ( plus45L_POLARIZATION )

              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)

            ! -45deg. linear polarisation
            CASE ( minus45L_POLARIZATION )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)

            ! Vertical, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! vertical, but the instrument scans cross-track.
            ! e = eV * (1-SIN^2(z))  +  eH * SIN^2(z)
            ! r = rV * (1-SIN^2(z))  +  rH * SIN^2(z)
            CASE ( VL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                SfcOptics%Emissivity(i,1) = (Emissivity(i,1)*(ONE-SIN2_Angle)) + &
                                            (Emissivity(i,2)*SIN2_Angle)
                SfcOptics%Reflectivity(i,1,i,1) = (Reflectivity(i,1,i,1)*(ONE-SIN2_Angle)) + &
                                                  (Reflectivity(i,2,i,2)*SIN2_Angle)
              END DO

            ! Horizontal, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! horizontal, but the instrument scans cross-track.
            ! e = eV * SIN^2(z)  +  eH * (1-SIN^2(z))
            ! r = rV * SIN^2(z)  +  rH * (1-SIN^2(z))
            CASE ( HL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                SfcOptics%Emissivity(i,1) = (Emissivity(i,1)*SIN2_Angle) + &
                                            (Emissivity(i,2)*(ONE-SIN2_Angle))
                SfcOptics%Reflectivity(i,1,i,1) = (Reflectivity(i,1,i,1)*SIN2_Angle) + &
                                                  (Reflectivity(i,2,i,2)*(ONE-SIN2_Angle))
              END DO

            ! Right circular polarisation
            CASE ( RC_POLARIZATION )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)

            ! Left circular polarisation
            CASE ( LC_POLARIZATION )
              SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
              SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)
            !
            ! Description:
            ! ============
            ! Polarization mixing with constant offset angle for TROPICS
            !
            ! Reference:
            ! ==========
            ! Leslie, V. (2020): TROPICS Polarization Description, 20 November 2020.
            ! (Personal Communication)
            ! 
            CASE ( CONST_MIXED_POLARIZATION )
              SIN2_Angle = (GeometryInfo%Distance_Ratio * &
                           SIN(DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)))**2
              DO i = 1, nZ
                SfcOptics%Emissivity(i,1) = (Emissivity(i,1)*(SIN2_Angle)) + &
                                              (Emissivity(i,2)*(ONE-SIN2_Angle))
                SfcOptics%Reflectivity(i,1,i,1) = (Reflectivity(i,1,i,1)*SIN2_Angle) + &
                                                  (Reflectivity(i,2,i,2)*(ONE-SIN2_Angle))
              END DO
            
            !
            ! Description:
            ! ============
            ! Polarization changing with a defined polarization rotation angle 
            ! as instrument zenith angle changes. Implemented for GEMS-1 SmallSat.
            !
            CASE ( PRA_POLARIZATION )
              DO i = 1, nZ
                ! Alias for the sensor scan angle:
                phi = GeometryInfo%Sensor_Scan_Radian
                ! Instrument offset angle: 
                theta_f = DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)
                ph = SIN(phi) * ( COS(phi) + SIN(theta_f)*(1.0_fp - COS(phi))  ) &
                   ! --------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                pv = - ( SIN(phi)**2 - SIN(theta_f)*(1.0_fp - COS(phi))*COS(phi) ) &
                   ! ---------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                ! Sine square of Polarization Rotation Angle (PRA)
                SIN2_Angle = SIN(ATAN( -pv/ph ))**2
                SfcOptics%Emissivity(i,1) = (Emissivity(i,1)*(SIN2_Angle)) + &
                                               (Emissivity(i,2)*(ONE-SIN2_Angle))
                SfcOptics%Reflectivity(i,1,i,1) = (Reflectivity(i,1,i,1)*SIN2_Angle) + &
                                                  (Reflectivity(i,2,i,2)*(ONE-SIN2_Angle))
              END DO
 
            ! Serious problem if we got to this points
            CASE DEFAULT
               Error_Status = FAILURE
               WRITE( Message,'("Unrecognised polarization flag for microwave ",&
                               &"channel index ",i0)' ) ChannelIndex
               CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
               RETURN

           END SELECT Polarization_Type

        ELSE


          ! ------------------------------------
          ! Coupled polarization from atmosphere
          ! considered. Simply copy the data
          ! ------------------------------------
          SfcOptics%Emissivity(1:nZ,1:nL)             = Emissivity(1:nZ,1:nL)
          SfcOptics%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) = Reflectivity(1:nZ,1:nL,1:nZ,1:nL)

        END IF Decoupled_Polarization



      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                      ## INFRARED CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE IF ( SpcCoeff_IsInfraredSensor( SC(SensorIndex) ) ) THEN

        ! -------------------------------------
        ! Infrared LAND emissivity/reflectivity
        ! -------------------------------------
        Infrared_Land: IF( Surface%Land_Coverage > ZERO ) THEN
          ! Compute the surface optics

          Error_Status = CSEM_Compute_LandIR_SfcOptics ( &
                         CSEM_Land,              &  ! Input
                         CSEM_SfcOptics,         &  ! Input
                         CSEM_Options  )            ! Output
          IF(COM_CRTM_CSEM) THEN
             !WRITE(*,*)'CSEM_Land',CSEM_Channel%wavenumber,CSEM_Geometry(1)%zenith_angle
             WRITE(*,*)'CSEM_Land',CSEM_SfcOptics%Emissivity(1,1:2)
          END IF

         IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR land SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
  
          ! Accumulate the surface optics properties
          ! based on land coverage fraction
          DO i = 1, nZ
             Emissivity(i,1)            =  CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Land_Coverage
             Reflectivity(i,1,i,1)      =  CSEM_SfcOptics%Reflectivity(i,1,i,1)     *  Surface%Land_Coverage
             Direct_Reflectivity(i,1)   =  CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Land_Coverage
          END DO
 
         END IF Infrared_Land


        ! --------------------------------------
        ! Infrared WATER emissivity/reflectivity
        ! --------------------------------------
        Infrared_Water: IF( Surface%Water_Coverage > ZERO ) THEN
          ! Compute the surface optics
          !CALL SET_AlgID("IR_WATER", NESDIS_IRW_Nalli)
          Error_Status = CSEM_Compute_WaterIR_SfcOptics ( &
                         CSEM_Water,                      &  ! Input
                         CSEM_SfcOptics,                  &  ! Input
                         CSEM_Options,ivar%CSEM_IRWSOV   )   ! Output
          IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Water,wavenumber, angle',CSEM_SfcOptics%wavenumber,CSEM_SfcOptics%Angle
              WRITE(*,*)'CSEM_Water,Emissivity',CSEM_SfcOptics%Emissivity(1,1:1)
         END IF

          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR water SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
 
          ! Accumulate the surface optics properties
          ! based on water coverage fraction
          DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Water_Coverage
             Reflectivity(i,1,i,1)      =  Reflectivity(i,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(i,1,i,1)     *  Surface%Water_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Water_Coverage
          END DO
        END IF Infrared_Water


        ! -------------------------------------
        ! Infrared SNOW emissivity/reflectivity
        ! -------------------------------------
        Infrared_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN
        
          ! Compute the surface optics
          Error_Status = CSEM_Compute_SnowIR_SfcOptics ( &
                         CSEM_Snow,               &  ! Input
                         CSEM_SfcOptics,          &  ! Input
                         CSEM_Options  )             ! input
          IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Snow',CSEM_SfcOptics%Emissivity(1,1:2)
          END IF

          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR snow SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Snow_Coverage
             Reflectivity(i,1,i,1)      =  Reflectivity(i,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(i,1,i,1)     *  Surface%Snow_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Snow_Coverage
          END DO

        ENDIF Infrared_Snow


        ! ------------------------------------
        ! Infrared ICE emissivity/reflectivity
        ! ------------------------------------
        Infrared_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN
         ! Compute the surface optics

         Error_Status = CSEM_Compute_IceIR_SfcOptics ( &
                        CSEM_Ice,              &  ! Input
                        CSEM_SfcOptics,        &  ! output
                        CSEM_Options  )           ! input
          IF(COM_CRTM_CSEM) THEN
            WRITE(*,*)'CSEM_Ice',CSEM_SfcOptics%Emissivity(1,1:2)
          END IF

          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR ice SfcOptics at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on Ice coverage fraction
          DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Ice_Coverage
             Reflectivity(i,1,i,1)      =  Reflectivity(i,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(i,1,i,1)     *  Surface%Ice_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Ice_Coverage
          END DO
   
        END IF Infrared_Ice


        ! -----------------------
        ! Assign the final result
        ! -----------------------
        SfcOptics%Emissivity(1:nZ,1)          = Emissivity(1:nZ,1)
        SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity(1:nZ,1,1:nZ,1)
        SfcOptics%Direct_Reflectivity(1:nZ,1) = Direct_Reflectivity(1:nZ,1)


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                       ## VISIBLE CALCULATIONS ##                     ##
      !## Visible part shares using the IR code, in which visible              ##
      !## lambertian emissivity/reflectivity can be computed for visible       ##
      !## wavenumber.                                                          ##
      !##########################################################################
      !##########################################################################

!old      ELSE IF ( SpcCoeff_IsVisibleSensor( SC(SensorIndex) ) ) THEN
      ELSE IF ( SpcCoeff_IsVisibleSensor( SC(SensorIndex)).or. &
          SpcCoeff_IsUltravioletSensor( SC(SensorIndex)) ) THEN

        mth_Azi_Test: IF( SfcOptics%mth_Azi == 0 ) THEN

          !  ==================
          !  Lambertian surface
          !  ==================

          ! -------------------------------------
          ! Visible LAND emissivity/reflectivity
          ! -------------------------------------
          Visible_Land: IF( Surface%Land_Coverage > ZERO ) THEN
            ! Compute the surface optics
            Error_Status = CSEM_Compute_LandVIS_SfcOptics(    &
                           CSEM_Land,                 &  ! Input
                           CSEM_SfcOptics,            &  ! output
                           CSEM_Options  )               ! input
            IF(COM_CRTM_CSEM) THEN
              WRITE(*,*)'CSEM_Land',CSEM_SfcOptics%wavenumber,CSEM_SfcOptics%Angle
              WRITE(*,*)'CSEM_Land',CSEM_SfcOptics%Direct_Reflectivity(1,1)
              END IF
 
            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing VIS land SfcOptics at ", &
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF

            ! Accumulate the surface optics properties
            ! based on land coverage fraction
            DO i = 1, nZ
               Emissivity(i,1)            =  CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Land_Coverage
               Reflectivity(1:nZ,1,i,1)   =  CSEM_SfcOptics%Reflectivity(1:nZ,1,i,1)  *  Surface%Land_Coverage
               Direct_Reflectivity(i,1)   =  CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Land_Coverage
            END DO

          END IF Visible_Land


          ! -------------------------------------
          ! Visible WATER emissivity/reflectivity
          ! -------------------------------------
          Visible_Water: IF( Surface%Water_Coverage > ZERO ) THEN

            ! Compute the surface optics
            Error_Status = CSEM_Compute_WaterVIS_SfcOptics (   &
                           CSEM_Water,                         &  ! Input
                           CSEM_SfcOptics,                     &  ! output
                           CSEM_Options ,ivar%CSEM_VISWSOV )    ! input
            IF(COM_CRTM_CSEM) THEN
              WRITE(*,*)'CSEM_Water',CSEM_SfcOptics%wavenumber,CSEM_SfcOptics%Angle
              WRITE(*,*)'CSEM_Water',CSEM_SfcOptics%Direct_Reflectivity(1,1) 
            END IF

            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing VIS water SfcOptics at ",&
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF

            ! Accumulate the surface optics properties
            ! based on water coverage fraction
            DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Water_Coverage
             Reflectivity(1:nZ,1,i,1)   =  Reflectivity(1:nZ,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(1:nZ,1,i,1)  *  Surface%Water_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Water_Coverage
           END DO

          END IF Visible_Water


          ! ------------------------------------
          ! Visible SNOW emissivity/reflectivity
          ! ------------------------------------
          Visible_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN

            ! Compute the surface optics
            Error_Status = CSEM_Compute_SnowVIS_SfcOptics (   &
                           CSEM_Snow,                 &  ! Input
                           CSEM_SfcOptics,            &  ! output
                           CSEM_Options  )               ! input
            IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Snow',CSEM_SfcOptics%wavenumber,CSEM_SfcOptics%Angle
             WRITE(*,*)'CSEM_Snow',CSEM_SfcOptics%Emissivity(1,1:2)
            END IF

            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing VIS snow SfcOptics at ",&
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF

            ! Accumulate the surface optics properties
            ! based on snow coverage fraction
            DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Snow_Coverage
             Reflectivity(1:nZ,1,i,1)   =  Reflectivity(1:nZ,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(1:nZ,1,i,1)  *  Surface%Snow_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Snow_Coverage
            END DO
          ENDIF Visible_Snow


          ! -----------------------------------
          ! Visible ICE emissivity/reflectivity
          ! -----------------------------------
          Visible_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN

           ! Compute the surface optics
            Error_Status = CSEM_Compute_IceVIS_SfcOptics (   &
                           CSEM_Ice,                 &  ! Input
                           CSEM_SfcOptics,           &  ! output
                           CSEM_Options  )              ! input
          IF(COM_CRTM_CSEM) THEN
             WRITE(*,*)'CSEM_Ice',CSEM_SfcOptics%wavenumber,CSEM_SfcOptics%Angle
             WRITE(*,*)'CSEM_Ice',CSEM_SfcOptics%Emissivity(1,1:2)
          END IF

            IF ( Error_Status /= SUCCESS ) THEN
              WRITE( Message,'("Error computing VIS ice SfcOptics at ",&
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN
            END IF

            ! Accumulate the surface optics properties
            ! based on Ice coverage fraction
      
            DO i = 1, nZ
             Emissivity(i,1)            =  Emissivity(i,1)           +   &
                 CSEM_SfcOptics%Emissivity(i,1)           *  Surface%Ice_Coverage
             Reflectivity(1:nZ,1,i,1)   =  Reflectivity(1:nZ,1,i,1)     +   &
                 CSEM_SfcOptics%Reflectivity(1:nZ,1,i,1)  *  Surface%Ice_Coverage
             Direct_Reflectivity(i,1)   =  Direct_Reflectivity(i,1)  +   &
                 CSEM_SfcOptics%Direct_Reflectivity(i,1)  *  Surface%Ice_Coverage
            END DO

          END IF Visible_Ice


          ! -----------------------
          ! Assign the final result
          ! -----------------------
          SfcOptics%Emissivity(1:nZ,1)          =  Emissivity(1:nZ,1)
          SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) =  Reflectivity(1:nZ,1,1:nZ,1)
          SfcOptics%Direct_Reflectivity(1:nZ,1) =  Direct_Reflectivity(1:nZ,1)

        ELSE

          SfcOptics%Emissivity(1:nZ,1)          =  ZERO
          SfcOptics%Reflectivity(1:nZ,1,1:nZ,1) =  ZERO
          SfcOptics%Direct_Reflectivity         =  ZERO

        END IF mth_Azi_Test



      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                        ## INVALID SENSOR TYPE ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE Sensor_Select

        Error_Status = FAILURE
        WRITE( Message,'("Unrecognised sensor type for channel index ",i0)' ) &
                       ChannelIndex
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN

      END IF Sensor_Select

  END FUNCTION CRTM_Compute_SfcOptics


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface optical properties
!       and populate the output SfcOptics_TL structure for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_SfcOptics_TL( &
!                       Surface     , &  ! Input
!                       SfcOptics   , &  ! Input
!                       Surface_TL  , &  ! Input
!                       GeometryInfo, &  ! Input
!                       SensorIndex , &  ! Input
!                       ChannelIndex, &  ! Input
!                       SfcOptics_TL, &  ! In/Output
!                       iVar          )  ! Internal variable input
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:      CRTM_Surface structure containing the tangent-linear
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       CRTM_GeometryInfo_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface optical properties required for the radiative
!                        transfer calculation.
!                        On Input:  The Secant_Angle component is assumed to
!                                   contain data.
!                        On Output: The Emissivity and Reflectivity components
!                                   will contain the required data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics_TL argument is IN OUT rather
!       than just OUT. This is necessary because the argument should be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_SfcOptics_TL( &
    Surface     , &  ! Input
    SfcOptics   , &  ! Input
    Surface_TL  , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics_TL, &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type)     , INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN)     :: SfcOptics
    TYPE(CRTM_Surface_type)     , INTENT(IN)     :: Surface_TL
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN OUT) :: SfcOptics_TL
    TYPE(iVar_type)             , INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_SfcOptics_TL'
    ! Local variables
    CHARACTER(ML) :: Message
    INTEGER :: i
    INTEGER :: nL, nZ
    INTEGER :: Polarization
    REAL(fp) :: SIN2_Angle
    REAL(fp) :: pv, ph
    REAL(fp) :: phi, theta_f
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Emissivity_TL
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES, &
                        SfcOptics%n_Angles,MAX_N_STOKES) :: Reflectivity_TL
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Direct_Reflectivity_TL

    TYPE(CSEM_Land_Surface)  :: CSEM_Land_TL
    TYPE(CSEM_Snow_Surface)  :: CSEM_Snow_TL
    TYPE(CSEM_Water_Surface) :: CSEM_Water_TL
    TYPE(CSEM_Ice_Surface)   :: CSEM_Ice_TL
    TYPE(CSEM_SfcOptics_Type):: CSEM_SfcOptics_TL
    TYPE(CSEM_Atmosphere_Parameters) :: Atmos_TL 


    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    nL = SfcOptics%n_Stokes
    nZ = SfcOptics%n_Angles
    Polarization = SC(SensorIndex)%Polarization( ChannelIndex )
    ! Initialise the local emissivity and reflectivities
    Emissivity_TL   = ZERO
    Reflectivity_TL = ZERO
    Direct_Reflectivity_TL = ZERO

    CALL CSEM_SfcOptics_TL%Init(N_Angles=nZ)
    CSEM_SfcOptics_TL%Frequency   =  SC(SensorIndex)%Frequency(ChannelIndex)
    CSEM_SfcOptics_TL%Wavenumber  =  SC(SensorIndex)%Wavenumber(ChannelIndex)
    CSEM_SfcOptics_TL%Is_Solar    =  SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex)
    CSEM_SfcOptics_TL%Angle(1:nZ) =  SfcOptics%Angle(1:nZ)
    CSEM_SfcOptics_TL%Weight(1:nZ)=  SfcOptics%Weight(1:nZ)

    Error_Status = CRTM_CSEM_Input_Adaptor( &
                   Surface_TL,      &  ! Input
                   CSEM_Land_TL,    &  ! Output
                   CSEM_Water_TL ,  &  ! Output
                   CSEM_Snow_TL,    &  ! Output
                   CSEM_Ice_TL)    



      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                     ## MICROWAVE CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      Sensor_Select: IF ( SpcCoeff_IsMicrowaveSensor( SC(SensorIndex) ) ) THEN

        ! --------------------------------------
        ! Microwave LAND emissivity/reflectivity
        ! --------------------------------------
        Microwave_Land: IF( Surface%Land_Coverage > ZERO) THEN

          ! Compute the surface optics
          Error_Status = CSEM_Compute_LandMW_SfcOptics_TL (     &
                         CSEM_Land_TL,                  &
                         CSEM_SfcOptics_TL,                     &
                         iVar%CSEM_MWLSOV   )          
          IF(COM_CRTM_CSEM) THEN
            WRITE(*,*)'CSEM_Land_TL',CSEM_SfcOptics_TL%Emissivity(1,1:4)
          END IF
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW land SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on land coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1:2)         =  CSEM_SfcOptics_TL%Emissivity(i,1:2)         *  Surface%Land_Coverage
             Reflectivity_TL(i,1:2,i,1:2) =  CSEM_SfcOptics_TL%Reflectivity(i,1:2,i,1:2) *  Surface%Land_Coverage
          END DO
          
        END IF Microwave_Land


        ! ---------------------------------------
        ! Microwave WATER emissivity/reflectivity
        ! ---------------------------------------
        Microwave_Water: IF( Surface%Water_Coverage > ZERO ) THEN
      
          Atmos_TL%Transmittance = SfcOptics_TL%Transmittance

          ! Compute the surface optics
          Error_Status =  CSEM_Compute_WaterMW_SfcOptics_TL( &
                          CSEM_Water_TL                     ,&  ! Input
                          Atmos_TL                          ,&
                          CSEM_SfcOptics_TL                 ,&
                          iVar%CSEM_MWWSOV)                          

          IF(COM_CRTM_CSEM) THEN
            print*,'CSEM_Water_Surface_TL',CSEM_Water_TL
            WRITE(*,*)'CSEM_water_TL',CSEM_SfcOptics_TL%Emissivity(1,1:4)
          END IF

     
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW water SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on water coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1:4)         =   Emissivity_TL(i,1:4)         +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1:4)         *  Surface%Water_Coverage
             Reflectivity_TL(i,1:4,i,1:4) =   Reflectivity_TL(i,1:4,i,1:4) +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1:4,i,1:4) *  Surface%Water_Coverage
          END DO

        END IF Microwave_Water


        ! --------------------------------------
        ! Microwave SNOW emissivity/reflectivity
        ! --------------------------------------
        Microwave_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN

          ! Compute the surface optics
          Error_Status = CSEM_Compute_SnowMW_SfcOptics_TL( CSEM_SfcOptics_TL )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW snow SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1:2)         =   Emissivity_TL(i,1:2)         +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1:2)         *  Surface%Snow_Coverage
             Reflectivity_TL(i,1:2,i,1:2) =   Reflectivity_TL(i,1:2,i,1:2) +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1:2,i,1:2) *  Surface%Snow_Coverage
          END DO

        ENDIF Microwave_Snow


        ! -------------------------------------
        ! Microwave ICE emissivity/reflectivity
        ! -------------------------------------

        Microwave_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN

          ! Compute the surface optics
          Error_Status = CSEM_Compute_IceMW_SfcOptics_TL( CSEM_SfcOptics_TL )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW ice SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1:2)         =   Emissivity_TL(i,1:2)         +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1:2)         *  Surface%Ice_Coverage
             Reflectivity_TL(i,1:2,i,1:2) =   Reflectivity_TL(i,1:2,i,1:2) +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1:2,i,1:2) *  Surface%Ice_Coverage
          END DO

        ENDIF Microwave_Ice



        !#----------------------------------------------------------------------#
        !#                 -- HANDLE THE DECOUPLED POLARISATION --              #
        !#                                                                      #
        !# The SfcOptics n_Stokes dimension determines whether the surface      #
        !# optics takes into account the second order effect of cross           #
        !# polarisation, e.g. if the surface optics for a purely vertically     #
        !# polarised channel has a horizontal (or other) component due to       #
        !# scattering at the surface.                                           #
        !#                                                                      #
        !# If the SfcOptics n_Stokes dimension == 1, the polarisations are      #
        !# decoupled.                                                           #
        !#----------------------------------------------------------------------#

        Decoupled_Polarization: IF( SfcOptics%n_Stokes == 1 ) THEN


          ! ------------------------------------------------------
          ! Decoupled polarisation. Branch on channel polarisation
          ! ------------------------------------------------------
          Polarization_Type: SELECT CASE( Polarization )

            ! The unpolarised case, I
            ! e = (eV + eH)/2
            ! r = (rV + rH)/2
            ! Note: INTENSITY == UNPOLARIZED == FIRST_STOKES_COMPONENT
            CASE( INTENSITY )
              SfcOptics_TL%Emissivity(1:nZ,1) = &
                POINT_5 * ( Emissivity_TL(1:nZ,1) + Emissivity_TL(1:nZ,2) )
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = &
                POINT_5 * ( Reflectivity_TL(1:nZ,1,1:nZ,1) + Reflectivity_TL(1:nZ,2,1:nZ,2) )

            ! The second Stokes component, Q, the polarisation difference.
            ! e = (eV - eH)/2
            ! r = (rV - rH)/2
            CASE( SECOND_STOKES_COMPONENT )
              SfcOptics_TL%Emissivity(1:nZ,1) = &
                POINT_5 * ( Emissivity_TL(1:nZ,1) - Emissivity_TL(1:nZ,2) )
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = &
                POINT_5 * ( Reflectivity_TL(1:nZ,1,1:nZ,1) - Reflectivity_TL(1:nZ,2,1:nZ,2) )

            ! The third Stokes component, U.
            CASE ( THIRD_STOKES_COMPONENT )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,3)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,3,1:nZ,3)

            ! The fourth Stokes component, V.
            CASE ( FOURTH_STOKES_COMPONENT )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,4)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,4,1:nZ,4)

            ! Vertical linear polarisation
            CASE ( VL_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)

            ! Horizontal linear polarisation
            CASE ( HL_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,2)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(:,2,:,2)

            ! +45deg. linear polarisation
            CASE ( plus45L_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)

            ! -45deg. linear polarisation
            CASE ( minus45L_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)

            ! Vertical, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! vertical, but the instrument scans cross-track.
            ! e = eV * (1-SIN^2(z))  +  eH * SIN^2(z)
            ! r = rV * (1-SIN^2(z))  +  rH * SIN^2(z)
            CASE ( VL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                SfcOptics_TL%Emissivity(i,1) = (Emissivity_TL(i,1)*(ONE-SIN2_Angle)) + &
                                               (Emissivity_TL(i,2)*SIN2_Angle)
                SfcOptics_TL%Reflectivity(i,1,i,1) = (Reflectivity_TL(i,1,i,1)*(ONE-SIN2_Angle)) + &
                                                     (Reflectivity_TL(i,2,i,2)*SIN2_Angle)
              END DO

            ! Horizontal, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! horizontal, but the instrument scans cross-track.
            ! e = eV * SIN^2(z)  +  eH * (1-SIN^2(z))
            ! r = rV * SIN^2(z)  +  rH * (1-SIN^2(z))
            CASE ( HL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                SfcOptics_TL%Emissivity(i,1) = (Emissivity_TL(i,1)*SIN2_Angle) + &
                                               (Emissivity_TL(i,2)*(ONE-SIN2_Angle))
                SfcOptics_TL%Reflectivity(i,1,i,1) = (Reflectivity_TL(i,1,i,1)*SIN2_Angle) + &
                                                     (Reflectivity_TL(i,2,i,2)*(ONE-SIN2_Angle))
              END DO

            ! Right circular polarisation
            CASE ( RC_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)

            ! Left circular polarisation
            CASE ( LC_POLARIZATION )
              SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
              SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)
 
            ! Polarization mixing with constant offset angle for TROPICS
            CASE ( CONST_MIXED_POLARIZATION )

              SIN2_Angle = (GeometryInfo%Distance_Ratio * &
                           SIN(DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)))**2
              DO i = 1, nZ
                SfcOptics_TL%Emissivity(1:nZ,1) = (Emissivity_TL(i,1)*(SIN2_Angle)) + &
                                                  (Emissivity_TL(i,2)*(ONE-SIN2_Angle))
                SfcOptics_TL%Reflectivity(i,1,i,1) =(Reflectivity_TL(i,1,i,1)*SIN2_Angle) + &
                                                     (Reflectivity_TL(i,2,i,2)*(ONE-SIN2_Angle))
              END DO

            !
            ! Description:
            ! ============
            ! Polarization changing with a defined polarization rotation angle 
            ! as instrument zenith angle changes. Implemented for GEMS-1 SmallSat.
            !
            CASE ( PRA_POLARIZATION )
              DO i = 1, nZ
                ! Alias for the sensor scan angle:
                phi = GeometryInfo%Sensor_Scan_Radian
                ! Instrument offset angle: 
                theta_f = DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)
                ph = SIN(phi) * ( COS(phi) + SIN(theta_f)*(1.0_fp - COS(phi))  ) &
                   ! --------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                pv = - ( SIN(phi)**2 - SIN(theta_f)*(1.0_fp - COS(phi))*COS(phi) ) &
                   ! ---------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                ! Sine square of Polarization Rotation Angle (PRA)
                SIN2_Angle = SIN(ATAN( -pv/ph ))**2
                SfcOptics_TL%Emissivity(i,1) = (Emissivity_TL(i,1)*(SIN2_Angle)) + &
                                               (Emissivity_TL(i,2)*(ONE-SIN2_Angle))
                SfcOptics_TL%Reflectivity(i,1,i,1) = (Reflectivity_TL(i,1,i,1)*SIN2_Angle) + &
                                                  (Reflectivity_TL(i,2,i,2)*(ONE-SIN2_Angle))
              END DO

            ! Serious problem if we got to this point
            CASE DEFAULT
              Error_Status = FAILURE
              WRITE( Message,'("Unrecognised polarization flag for microwave ",&
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN

           END SELECT Polarization_Type


        ELSE


          ! ------------------------------------
          ! Coupled polarization from atmosphere
          ! considered. Simply copy the data
          ! ------------------------------------
          SfcOptics_TL%Emissivity   = Emissivity_TL(1:nZ,1:nL)
          SfcOptics_TL%Reflectivity = Reflectivity_TL(1:nZ,1:nL,1:nZ,1:nL)

        END IF Decoupled_Polarization



      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                      ## INFRARED CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE IF ( SpcCoeff_IsInfraredSensor( SC(SensorIndex) ) ) THEN


        ! -------------------------------------
        ! Infrared LAND emissivity/reflectivity
        ! -------------------------------------
        Infrared_Land: IF( Surface%Land_Coverage > ZERO ) THEN

          ! Compute the surface optics
          ! **STUB PROCEDURE**
          Error_Status = CSEM_Compute_LandIR_SfcOptics_TL( CSEM_SfcOptics_TL )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR land SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on land coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1)          =  CSEM_SfcOptics_TL%Emissivity(i,1)           *   Surface%Land_Coverage
             Reflectivity_TL(i,1,i,1)    =  CSEM_SfcOptics_TL%Reflectivity(i,1,i,1)     *   Surface%Land_Coverage
             Direct_Reflectivity_TL(i,1) =  CSEM_SfcOptics_TL%Direct_Reflectivity(i,1)  *   Surface%Land_Coverage
          END DO

        END IF Infrared_Land


        ! --------------------------------------
        ! Infrared WATER emissivity/reflectivity
        ! --------------------------------------
        Infrared_Water: IF( Surface%Water_Coverage > ZERO ) THEN

          ! Compute the surface optics
           Error_Status =  CSEM_Compute_WaterIR_SfcOptics_TL(   &
             CSEM_Water_TL                             ,&  
             CSEM_SfcOptics_TL, iVar%CSEM_IRWSOV)                         
          IF(COM_CRTM_CSEM) THEN
            WRITE(*,*)'CSEM_water_TL',CSEM_SfcOptics_TL%Emissivity(1,1)
            WRITE(*,*)'CSEM_water_TL',CSEM_SfcOptics_TL%Direct_Reflectivity(1,1)
          END IF
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR water SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on water coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1)          =  Emissivity_TL(i,1)           +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1)          *   Surface%Water_Coverage
             Reflectivity_TL(i,1,i,1)    =  Reflectivity_TL(i,1,i,1)     +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1,i,1)    *   Surface%Water_Coverage
             Direct_Reflectivity_TL(i,1) =  Direct_Reflectivity_TL(i,1)  +   &
                 CSEM_SfcOptics_TL%Direct_Reflectivity(i,1) *   Surface%Water_Coverage
          END DO
  
        END IF Infrared_Water


        ! -------------------------------------
        ! Infrared SNOW emissivity/reflectivity
        ! -------------------------------------
        Infrared_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN

          ! Compute the surface optics
          Error_Status = CSEM_Compute_SnowIR_SfcOptics_TL( CSEM_SfcOptics_TL )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR snow SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on snow coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1)          =  Emissivity_TL(i,1)           +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1)          *   Surface%Snow_Coverage
             Reflectivity_TL(i,1,i,1)    =  Reflectivity_TL(i,1,i,1)     +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1,i,1)    *   Surface%Snow_Coverage
             Direct_Reflectivity_TL(i,1) =  Direct_Reflectivity_TL(i,1)  +   &
                 CSEM_SfcOptics_TL%Direct_Reflectivity(i,1) *   Surface%Snow_Coverage
          END DO

        END IF Infrared_Snow


        ! ------------------------------------
        ! Infrared ICE emissivity/reflectivity
        ! ------------------------------------
        Infrared_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN

          ! Compute the surface optics
          Error_Status = CSEM_Compute_IceIR_SfcOptics_TL( CSEM_SfcOptics_TL )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR ice SfcOptics_TL at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

          ! Accumulate the surface optics properties
          ! based on Ice coverage fraction
          DO i = 1, nZ
             Emissivity_TL(i,1)          =  Emissivity_TL(i,1)           +   &
                 CSEM_SfcOptics_TL%Emissivity(i,1)          *   Surface%Ice_Coverage
             Reflectivity_TL(i,1,i,1)    =  Reflectivity_TL(i,1,i,1)     +   &
                 CSEM_SfcOptics_TL%Reflectivity(i,1,i,1)    *   Surface%Ice_Coverage
             Direct_Reflectivity_TL(i,1) =  Direct_Reflectivity_TL(i,1)  +   &
                 CSEM_SfcOptics_TL%Direct_Reflectivity(i,1) *   Surface%Ice_Coverage
          END DO

        END IF Infrared_Ice


        ! -----------------------
        ! Assign the final result
        ! -----------------------
        SfcOptics_TL%Emissivity(1:nZ,1)          = Emissivity_TL(1:nZ,1)
        SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = Reflectivity_TL(1:nZ,1,1:nZ,1)
        SfcOptics_TL%Direct_Reflectivity(1:nZ,1) = Direct_Reflectivity_TL(1:nZ,1)


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                       ## VISIBLE CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      !ELSE IF ( SpcCoeff_IsVisibleSensor( SC(SensorIndex) ) ) THEN
      ELSE IF (SpcCoeff_IsVisibleSensor(SC(SensorIndex)).or. &
           SpcCoeff_IsUltravioletSensor(SC(SensorIndex))) THEN


        ! -------------------
        ! Default values only
        ! -------------------
        SfcOptics_TL%Emissivity(1:nZ,1)          = ZERO
        SfcOptics_TL%Reflectivity(1:nZ,1,1:nZ,1) = ZERO
        SfcOptics_TL%Direct_Reflectivity = ZERO


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                        ## INVALID SENSOR TYPE ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE Sensor_Select

        Error_Status = FAILURE
        WRITE( Message,'("Unrecognised sensor type for channel index ",i0)' ) &
                       ChannelIndex
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN

      END IF Sensor_Select

  END FUNCTION CRTM_Compute_SfcOptics_TL


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface optical properties
!       for a single channel.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_Compute_SfcOptics_AD( &
!                        Surface     , &  ! Input
!                        SfcOptics   , &  ! Input
!                        SfcOptics_AD, &  ! Input
!                        GeometryInfo, &  ! Input
!                        SensorIndex , &  ! Input
!                        ChannelIndex, &  ! Input
!                        Surface_AD  , &  ! Output
!                        iVar          )  ! Internal variable input
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_AD:    CRTM_SfcOptics structure containing the adjoint
!                        surface optical properties.
!                        **NOTE: On EXIT from this function, the contents of
!                                this structure may be modified (e.g. set to
!                                zero.)
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:    CRTM_GeometryInfo structure containing the
!                        view geometry information.
!                        UNITS:      N/A
!                        TYPE:       CRTM_GeometryInfo_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Surface_AD:      CRTM_Surface structure containing the adjoint
!                        surface state data.
!                        **NOTE: On ENTRY to this function, the contents of
!                                this structure should be defined (e.g.
!                                initialized to some value based on the
!                                position of this function in the call chain.)
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the ERROR_HANDLER module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on all of the adjoint arguments (whether input or output)
!       is IN OUT rather than just OUT. This is necessary because the INPUT
!       adjoint arguments are modified, and the OUTPUT adjoint arguments must
!       be defined prior to entry to this routine. So, anytime a structure is
!       to be output, to prevent memory leaks the IN OUT INTENT is a must.
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_Compute_SfcOptics_AD( &
    Surface     , &  ! Input
    SfcOptics   , &  ! Input
    SfcOptics_AD, &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Surface_AD  , &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type)     , INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN)     :: SfcOptics
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Surface_type)     , INTENT(IN OUT) :: Surface_AD
    TYPE(iVar_type)             , INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Compute_SfcOptics_AD'
    ! Local variables
    CHARACTER(256)  :: Message
    INTEGER :: i
    INTEGER :: nL, nZ
    INTEGER :: Polarization
    REAL(fp) :: SIN2_Angle
    REAL(fp) :: pv, ph
    REAL(fp) :: phi, theta_f
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Emissivity_AD
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES, &
                        SfcOptics%n_Angles,MAX_N_STOKES) :: Reflectivity_AD
    REAL(fp), DIMENSION(SfcOptics%n_Angles,MAX_N_STOKES) :: Direct_Reflectivity_AD
    
    TYPE(CSEM_Land_Surface)  :: CSEM_Land_AD
    TYPE(CSEM_Snow_Surface)  :: CSEM_Snow_AD
    TYPE(CSEM_Water_Surface) :: CSEM_Water_AD
    TYPE(CSEM_Ice_Surface)   :: CSEM_Ice_AD
    TYPE(CSEM_SfcOptics_Type):: CSEM_SfcOptics_AD

    TYPE(CSEM_Atmosphere_Parameters)  :: Atmos_AD 
 
   

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
    nL = SfcOptics%n_Stokes
    nZ = SfcOptics%n_Angles
    Polarization = SC(SensorIndex)%Polarization( ChannelIndex )
    ! Initialise the local emissivity and reflectivity adjoints
    Emissivity_AD = ZERO
    Reflectivity_AD = ZERO
    Direct_Reflectivity_AD = ZERO

    CALL CSEM_SfcOptics_AD%Init(N_Angles=nZ)
    CSEM_SfcOptics_AD%Frequency    =  SC(SensorIndex)%Frequency(ChannelIndex)
    CSEM_SfcOptics_AD%Wavenumber   =  SC(SensorIndex)%Wavenumber(ChannelIndex)
    CSEM_SfcOptics_AD%Is_Solar     =  SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex)
    CSEM_SfcOptics_AD%Angle(1:nZ)  =  SfcOptics%Angle(1:nZ)
    CSEM_SfcOptics_AD%Weight(1:nZ) = SfcOptics%Weight(1:nZ)
 
    Error_Status = CRTM_CSEM_Input_Adaptor(   &
                   Surface_AD,                &  ! Input
                   CSEM_Land_AD,      &  ! Output
                   CSEM_Water_AD,     &  ! Output
                   CSEM_Snow_AD,      &  ! Output
                   CSEM_Ice_AD)    


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                     ## MICROWAVE CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      Sensor_Select: IF ( SpcCoeff_IsMicrowaveSensor( SC(SensorIndex) ) ) THEN

        !#----------------------------------------------------------------------#
        !#                 -- HANDLE THE DECOUPLED POLARISATION --              #
        !#                                                                      #
        !# The SfcOptics n_Stokes dimension determines whether the surface      #
        !# optics takes into account the second order effect of cross           #
        !# polarisation, e.g. if the surface optics for a purely vertically     #
        !# polarised channel has a horizontal (or other) component due to       #
        !# scattering at the surface.                                           #
        !#                                                                      #
        !# If the SfcOptics n_Stokes dimension == 1, the polarisations are      #
        !# decoupled.                                                           #
        !#----------------------------------------------------------------------#
        Decoupled_Polarization: IF( SfcOptics%n_Stokes == 1 ) THEN


          ! ------------------------------------------------------
          ! Decoupled polarisation. Branch on channel polarisation
          ! ------------------------------------------------------
          Polarization_Type: SELECT CASE( Polarization )

            ! The unpolarised case, I
            ! e = (eV + eH)/2
            ! r = (rV + rH)/2
            ! Note: INTENSITY == UNPOLARIZED == FIRST_STOKES_COMPONENT
            CASE( INTENSITY )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              Emissivity_AD(1:nZ,2) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              Reflectivity_AD(1:nZ,2,1:nZ,2) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! The second Stokes component, Q, the polarisation difference.
            ! e = (eV - eH)/2
            ! r = (rV - rH)/2
            CASE( SECOND_STOKES_COMPONENT )
              Emissivity_AD(1:nZ,1) =  SfcOptics_AD%Emissivity(1:nZ,1)
              Emissivity_AD(1:nZ,2) = -SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) =  SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              Reflectivity_AD(1:nZ,2,1:nZ,2) = -SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! The third Stokes component, U.
            CASE ( THIRD_STOKES_COMPONENT )
              Emissivity_AD(1:nZ,3) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,3,1:nZ,3) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! The fourth Stokes component, V.
            CASE ( FOURTH_STOKES_COMPONENT )
              Emissivity_AD(1:nZ,4) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,4,1:nZ,4) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! Vertical linear polarisation
            CASE ( VL_POLARIZATION )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! Horizontal linear polarisation
            CASE ( HL_POLARIZATION )
              Emissivity_AD(1:nZ,2) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,2,1:nZ,2) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! +45deg. linear polarisation
            CASE ( plus45L_POLARIZATION )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! -45deg. linear polarisation
            CASE ( minus45L_POLARIZATION )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! Vertical, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! vertical, but the instrument scans cross-track.
            ! e = eV * (1-SIN^2(z))  +  eH * SIN^2(z)
            ! r = rV * (1-SIN^2(z))  +  rH * SIN^2(z)
            CASE ( VL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                Emissivity_AD(i,1) = SfcOptics_AD%Emissivity(i,1)*(ONE-SIN2_Angle)
                Emissivity_AD(i,2) = SfcOptics_AD%Emissivity(i,1)*SIN2_Angle
                Reflectivity_AD(i,1,i,1) = SfcOptics_AD%Reflectivity(i,1,i,1)*(ONE-SIN2_Angle)
                Reflectivity_AD(i,2,i,2) = SfcOptics_AD%Reflectivity(i,1,i,1)*SIN2_Angle
              END DO
              SfcOptics_AD%Emissivity   = ZERO
              SfcOptics_AD%Reflectivity = ZERO

            ! Horizontal, mixed polarisation. This category of polarisation is
            ! for those microwave channels where the nadir polarisation is
            ! horizontal, but the instrument scans cross-track.
            ! e = eV * SIN^2(z)  +  eH * (1-SIN^2(z))
            ! r = rV * SIN^2(z)  +  rH * (1-SIN^2(z))
            CASE ( HL_MIXED_POLARIZATION )
              DO i = 1, nZ
                SIN2_Angle = (GeometryInfo%Distance_Ratio*SIN(DEGREES_TO_RADIANS*SfcOptics%Angle(i)))**2
                Emissivity_AD(i,1) = SfcOptics_AD%Emissivity(i,1)*SIN2_Angle
                Emissivity_AD(i,2) = SfcOptics_AD%Emissivity(i,1)*(ONE-SIN2_Angle)
                Reflectivity_AD(i,1,i,1) = SfcOptics_AD%Reflectivity(i,1,i,1)*SIN2_Angle
                Reflectivity_AD(i,2,i,2) = SfcOptics_AD%Reflectivity(i,1,i,1)*(ONE-SIN2_Angle)
              END DO
              SfcOptics_AD%Emissivity = ZERO
              SfcOptics_AD%Reflectivity = ZERO

            ! Right circular polarisation
            CASE ( RC_POLARIZATION )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! Left circular polarisation
            CASE ( LC_POLARIZATION )
              Emissivity_AD(1:nZ,1) = SfcOptics_AD%Emissivity(1:nZ,1)
              SfcOptics_AD%Emissivity = ZERO
              Reflectivity_AD(1:nZ,1,1:nZ,1) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1)
              SfcOptics_AD%Reflectivity = ZERO

            ! Polarization mixing with constant offset angle for TROPICS
            CASE ( CONST_MIXED_POLARIZATION )
              SIN2_Angle = (GeometryInfo%Distance_Ratio * &
                           SIN(DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)))**2
              DO i = 1, nZ
                ! PS: The adjoint is the transpose of the TL relationship:
                ! eV_AD = e_AD * SIN^2(theta)
                ! eH_AD = e_AD * COS^2(theta)
                Emissivity_AD(i,1) = SfcOptics_AD%Emissivity(i,1)*SIN2_Angle
                Emissivity_AD(i,2) = SfcOptics_AD%Emissivity(i,1)*(ONE-SIN2_Angle)
                Reflectivity_AD(i,1,i,1) = SfcOptics_AD%Reflectivity(i,1,i,1)*SIN2_Angle
                Reflectivity_AD(i,2,i,2) = SfcOptics_AD%Reflectivity(i,1,i,1)*(ONE-SIN2_Angle)
              END DO
              SfcOptics_AD%Emissivity   = ZERO
              SfcOptics_AD%Reflectivity = ZERO

            !
            ! Description:
            ! ============
            ! Polarization changing with a defined polarization rotation angle 
            ! as instrument zenith angle changes. Implemented for GEMS-1 SmallSat.
            !
            CASE ( PRA_POLARIZATION )
              DO i = 1, nZ
                ! Alias for the sensor scan angle:
                phi = GeometryInfo%Sensor_Scan_Radian
                ! Instrument offset angle: 
                theta_f = DEGREES_TO_RADIANS*SC(SensorIndex)%PolAngle(ChannelIndex)
                ph = SIN(phi) * ( COS(phi) + SIN(theta_f)*(1.0_fp - COS(phi))  ) &
                   ! --------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                pv = - ( SIN(phi)**2 - SIN(theta_f)*(1.0_fp - COS(phi))*COS(phi) ) &
                   ! ---------------------------------------------------------------
                     / SQRT( SIN(phi)**2 + SIN(theta_f)**2*(1.0_fp - COS(phi)**2) )
                ! Sine square of Polarization Rotation Angle (PRA)
                SIN2_Angle = SIN(ATAN( -pv/ph ))**2
                ! PS: The adjoint is the transpose of the TL relationship:
                ! eV_AD = e_AD * SIN^2(theta)
                ! eH_AD = e_AD * COS^2(theta)
                Emissivity_AD(i,1) = SfcOptics_AD%Emissivity(i,1)*SIN2_Angle
                Emissivity_AD(i,2) = SfcOptics_AD%Emissivity(i,1)*(ONE-SIN2_Angle)
                Reflectivity_AD(i,1,i,1) = SfcOptics_AD%Reflectivity(i,1,i,1)*SIN2_Angle
                Reflectivity_AD(i,2,i,2) = SfcOptics_AD%Reflectivity(i,1,i,1)*(ONE-SIN2_Angle)
              END DO
              SfcOptics_AD%Emissivity   = ZERO
              SfcOptics_AD%Reflectivity = ZERO

  
            ! Serious problem if we got to this point
            CASE DEFAULT
              Error_Status = FAILURE
              WRITE( Message,'("Unrecognised polarization flag for microwave ",&
                              &"channel index ",i0)' ) ChannelIndex
              CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
              RETURN

          END SELECT Polarization_Type


        ELSE


          ! ------------------------------------
          ! Coupled polarization from atmosphere
          ! considered. Simply copy the data
          ! ------------------------------------
          Emissivity_AD(1:nZ,1:nL) = SfcOptics_AD%Emissivity(1:nZ,1:nL)
          SfcOptics_AD%Emissivity = ZERO
          Reflectivity_AD(1:nZ,1:nL,1:nZ,1:nL) = SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL)
          SfcOptics_AD%Reflectivity = ZERO

        END IF Decoupled_Polarization

 

        ! -------------------------------------
        ! Microwave ICE emissivity/reflectivity
        ! -------------------------------------
        Microwave_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN

          ! The surface optics properties based on ice coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:2) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:2) + &
            (Emissivity_AD(1:nZ,1:2)*Surface%Ice_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) + &
            (Reflectivity_AD(1:nZ,1:2,1:nZ,1:2)*Surface%Ice_Coverage)

          ! Compute the surface optics adjoints
          Error_Status = CSEM_Compute_IceMW_SfcOptics_AD( CSEM_Ice_AD )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW ice SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          SfcOptics_AD%Emissivity = ZERO
          SfcOptics_AD%Reflectivity = ZERO
        END IF Microwave_Ice


        ! --------------------------------------
        ! Microwave SNOW emissivity/reflectivity
        ! --------------------------------------

        Microwave_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN

          ! The surface optics properties based on snow coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:2) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:2) + &
            (Emissivity_AD(1:nZ,1:2)*Surface%Snow_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) + &
            (Reflectivity_AD(1:nZ,1:2,1:nZ,1:2)*Surface%Snow_Coverage)

          ! Compute the surface optics adjoints
          Error_Status = CSEM_Compute_SnowMW_SfcOptics_AD( CSEM_Snow_AD )
         IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW snow SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          SfcOptics_AD%Emissivity = ZERO
          SfcOptics_AD%Reflectivity = ZERO

        END IF Microwave_Snow


        ! ---------------------------------------
        ! Microwave WATER emissivity/reflectivity
        ! ---------------------------------------
        Microwave_Water: IF( Surface%Water_Coverage > ZERO ) THEN

          !Surface_AD%Wind_Speed=0.0_fp;Surface_AD%Water_Temperature=0.0_fp
          !Surface_AD%Salinity =0.0_fp      

          ! The surface optics properties based on water coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:4) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:2) + &
            (Emissivity_AD(1:nZ,1:4)*Surface%Water_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:4,1:nZ,1:4) = &
           ! SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) + &
            (Reflectivity_AD(1:nZ,1:4,1:nZ,1:4)*Surface%Water_Coverage)
    
          ! Compute the surface optics adjoints
          DO i = 1, nZ
            CSEM_SfcOptics_AD%Emissivity(i,1:4)         = SfcOptics_AD%Emissivity(i,1:4)
            CSEM_SfcOptics_AD%Reflectivity(i,1:4,i,1:4) = SfcOptics_AD%Reflectivity(i,1:4,i,1:4)
          END DO 
         
          Atmos_AD%Transmittance = SfcOptics_AD%Transmittance

          Error_Status =  CSEM_Compute_WaterMW_SfcOptics_AD( &
                          CSEM_SfcOptics_AD,   &  
                          CSEM_Water_AD,       &  
                          Atmos_AD,            &  
                         iVar%CSEM_MWWSOV)                          
         IF(COM_CRTM_CSEM) THEN   
           WRITE(*,*) 'CSEM_Water_AD',CSEM_Water_AD%Wind_Speed, &
             CSEM_Water_AD%Water_Temperature,&
             CSEM_Water_AD%Salinity
         ENDIF
         IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW water SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
  
          Error_Status =  CSEM_CRTM_Input_Adaptor(Surface_AD, CSEM_Water=CSEM_Water_AD) 
          SfcOptics_AD%Transmittance = Atmos_AD%Transmittance
          SfcOptics_AD%Emissivity = ZERO
          SfcOptics_AD%Reflectivity = ZERO
       END IF Microwave_Water


        ! --------------------------------------
        ! Microwave LAND emissivity/reflectivity
        ! --------------------------------------
        Microwave_Land: IF( Surface%Land_Coverage > ZERO ) THEN

          ! The surface optics properties based on land coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:2) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:2) + &
            (Emissivity_AD(1:nZ,1:2)*Surface%Land_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:2,1:nZ,1:2) + &
            (Reflectivity_AD(1:nZ,1:2,1:nZ,1:2)*Surface%Land_Coverage)
          ! Compute the surface optics adjoints
          DO i = 1, nZ
            CSEM_SfcOptics_AD%Emissivity(i,1:2)         = SfcOptics_AD%Emissivity(i,1:2)
            CSEM_SfcOptics_AD%Reflectivity(i,1:2,i,1:2) = SfcOptics_AD%Reflectivity(i,1:2,i,1:2)
          END DO 

          ! Compute the surface optics adjoints
          Error_Status = CSEM_Compute_LandMW_SfcOptics_AD (     &
                         CSEM_SfcOptics_AD,                     &
                         CSEM_Land_AD,                  &
                         iVar%CSEM_MWLSOV   )      
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing MW land SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
          Error_Status =  CSEM_CRTM_Input_Adaptor(Surface_AD,CSEM_Land=CSEM_Land_AD) 
          SfcOptics_AD%Emissivity = ZERO
          SfcOptics_AD%Reflectivity = ZERO

        END IF Microwave_Land


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                      ## INFRARED CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE IF ( SpcCoeff_IsInfraredSensor( SC(SensorIndex) ) ) THEN
        Reflectivity_AD(1:nZ,1,1:nZ,1:nL) = SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1:nL)
        SfcOptics_AD%Reflectivity = ZERO        
        Emissivity_AD(1:nZ,1:nL) = SfcOptics_AD%Emissivity(1:nZ,1:nL)
        SfcOptics_AD%Emissivity = ZERO
        Direct_Reflectivity_AD(1:nZ,1) = SfcOptics_AD%Direct_Reflectivity(1:nZ,1)
        SfcOptics_AD%Direct_Reflectivity(1:nZ,1) = ZERO

        ! ------------------------------------
        ! Infrared ICE emissivity/reflectivity
        ! ------------------------------------
        Infrared_Ice: IF( Surface%Ice_Coverage > ZERO ) THEN

          ! The surface optics properties based on ice coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:nL) + &
            (Emissivity_AD(1:nZ,1:nL)*Surface%Ice_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) + &
            (Reflectivity_AD(1:nZ,1:nL,1:nZ,1:nL)*Surface%Ice_Coverage)
          SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) + &
            (Direct_Reflectivity_AD(1:nZ,1:nL)*Surface%Ice_Coverage) 

          ! Compute the surface optics adjoints
          Error_Status = CSEM_Compute_IceIR_SfcOptics_AD( CSEM_Ice_AD )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR ice SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

        END IF Infrared_Ice

        ! -------------------------------------
        ! Infrared SNOW emissivity/reflectivity
        ! -------------------------------------
        Infrared_Snow: IF( Surface%Snow_Coverage > ZERO ) THEN

          ! The surface optics properties based on snow coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:nL) + &
            (Emissivity_AD(1:nZ,1:nL)*Surface%Snow_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) + &
            (Reflectivity_AD(1:nZ,1:nL,1:nZ,1:nL)*Surface%Snow_Coverage)
          SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) + &
            (Direct_Reflectivity_AD(1:nZ,1:nL)*Surface%Snow_Coverage) 

          ! Compute the surface optics adjoints
          Error_Status = CSEM_Compute_SnowIR_SfcOptics_AD( CSEM_Snow_AD )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR snow SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

        END IF Infrared_Snow

        ! --------------------------------------
        ! Infrared WATER emissivity/reflectivity
        ! --------------------------------------
        Infrared_Water: IF ( Surface%Water_Coverage > ZERO ) THEN

          ! The surface optics properties based on water coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:nL) + &
            (Emissivity_AD(1:nZ,1:nL)*Surface%Water_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) + &
            (Reflectivity_AD(1:nZ,1:nL,1:nZ,1:nL)*Surface%Water_Coverage)
          SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) + &
            (Direct_Reflectivity_AD(1:nZ,1:nL)*Surface%Water_Coverage) 

          DO i = 1, nZ
             CSEM_SfcOptics_AD%Emissivity(i,1:2)           =  SfcOptics_AD%Emissivity(i,1:2)
             CSEM_SfcOptics_AD%Reflectivity(i,1:2,i,1:2)   =  SfcOptics_AD%Reflectivity(i,1:2,i,1:2)
             CSEM_SfcOptics_AD%Direct_Reflectivity(i,1:2)  =  SfcOptics_AD%Direct_Reflectivity(i,1:2)
          END DO
          !Surface_AD%Wind_Speed=0.0_fp;Surface_AD%Water_Temperature=0.0_fp
          !Surface_AD%Salinity =0.0_fp      
  
          Error_Status = CSEM_Compute_WaterIR_SfcOptics_AD(       &
                         CSEM_SfcOptics_AD,                       &  
                         CSEM_Water_AD,                   &
                         iVar%CSEM_IRWSOV)  

          IF(COM_CRTM_CSEM) THEN   
             WRITE(*,*) 'CSEM_Water_AD',CSEM_Water_AD%Wind_Speed
          END IF
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR water SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF
  
          Error_Status =  CSEM_CRTM_Input_Adaptor(Surface_AD, CSEM_Water=CSEM_Water_AD) 

        END IF Infrared_Water


        ! --------------------------------------
        ! Infrared LAND emissivity/reflectivity
        ! --------------------------------------
        Infrared_Land: IF( Surface%Land_Coverage > ZERO ) THEN

          ! The surface optics properties based on land coverage fraction
          ! Note that the Emissivity_AD and Reflectivity_AD local adjoints
          ! are NOT zeroed here.
          SfcOptics_AD%Emissivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Emissivity(1:nZ,1:nL) + &
            (Emissivity_AD(1:nZ,1:nL)*Surface%Land_Coverage)
          SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) = &
            !SfcOptics_AD%Reflectivity(1:nZ,1:nL,1:nZ,1:nL) + &
            (Reflectivity_AD(1:nZ,1:nL,1:nZ,1:nL)*Surface%Land_Coverage)
          SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) = &
            !SfcOptics_AD%Direct_Reflectivity(1:nZ,1:nL) + &
            (Direct_Reflectivity_AD(1:nZ,1:nL)*Surface%Land_Coverage) 

          ! Compute the surface optics adjoints
          ! **STUB PROCEDURE**
          Error_Status = CSEM_Compute_LandIR_SfcOptics_AD( CSEM_Land_AD )
          IF ( Error_Status /= SUCCESS ) THEN
            WRITE( Message,'("Error computing IR land SfcOptics_AD at ",&
                            &"channel index ",i0)' ) ChannelIndex
            CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
            RETURN
          END IF

        END IF Infrared_Land



      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                       ## VISIBLE CALCULATIONS ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      !ELSE IF ( SpcCoeff_IsVisibleSensor( SC(SensorIndex) ) ) THEN
      ELSE IF (SpcCoeff_IsVisibleSensor(SC(SensorIndex)).or. &
           SpcCoeff_IsUltravioletSensor(SC(SensorIndex)) ) THEN


        ! -------------------
        ! Default values only
        ! -------------------
        SfcOptics_AD%Emissivity(1:nZ,1)          = ZERO
        SfcOptics_AD%Reflectivity(1:nZ,1,1:nZ,1) = ZERO
        SfcOptics_AD%Direct_Reflectivity         = ZERO


      !##########################################################################
      !##########################################################################
      !##                                                                      ##
      !##                        ## INVALID SENSOR TYPE ##                     ##
      !##                                                                      ##
      !##########################################################################
      !##########################################################################

      ELSE Sensor_Select
        Error_Status = FAILURE
        WRITE( Message,'("Unrecognised sensor type for channel index ",i0)' ) &
                       ChannelIndex
        CALL Display_Message( ROUTINE_NAME, Message, Error_Status )
        RETURN

      END IF Sensor_Select

  END FUNCTION CRTM_Compute_SfcOptics_AD
  
  
  
!===============================================================================
!
!===============================================================================
 FUNCTION CRTM_CSEM_Input_Adaptor( &
    Surface,     &  ! Input
    CSEM_Land,   &  ! Output
    CSEM_Water,  &  ! Output
    CSEM_Snow,   &  ! Output
    CSEM_Ice)    &  ! Output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type)          , INTENT(IN)   :: Surface
    TYPE(CSEM_Land_Surface) ,OPTIONAL, INTENT(OUT)  :: CSEM_Land
    TYPE(CSEM_Water_Surface),OPTIONAL, INTENT(OUT)  :: CSEM_Water
    TYPE(CSEM_Snow_Surface) ,OPTIONAL, INTENT(OUT)  :: CSEM_Snow
    TYPE(CSEM_Ice_Surface)  ,OPTIONAL, INTENT(OUT)  :: CSEM_Ice
 
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_CSEM_Input_Adaptor'
    ! Local variables
    CHARACTER(ML) :: Message


    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS

 
    IF( PRESENT(CSEM_Land) ) THEN
        CSEM_Land%Top_Soil_Moisture     = Surface%Soil_Moisture_Content
        CSEM_Land%vegetation_Fraction   = Surface%Vegetation_Fraction
        CSEM_Land%Top_Soil_Temperature  = Surface%Soil_Temperature
        CSEM_Land%Land_Skin_Temperature = Surface%Land_Temperature
        CSEM_Land%LAI                   = Surface%Lai
        CSEM_Land%soil_Type             = Surface%Soil_Type
        CSEM_Land%Vegetation_Type       = Surface%Vegetation_Type
        CSEM_Land%Land_Cover_Type       = Surface%Land_Type

    END IF
     
    IF( PRESENT(CSEM_Snow) ) THEN
   
        CSEM_Snow%Snow_Type             = Surface%Snow_Type
        CSEM_Snow%Snow_Temperature      = Surface%Snow_Temperature
        CSEM_Snow%Snow_Depth            = Surface%Snow_Depth
        CSEM_Snow%Snow_Density          = Surface%Snow_Density 
        CSEM_Snow%Snow_Grain_Size       = Surface%Snow_Grain_Size
        CSEM_Snow%soil_Type             = Surface%Soil_Type
        CSEM_Snow%Top_Soil_Temperature  = Surface%Soil_Temperature
        CSEM_Snow%Top_Soil_Moisture_Content = Surface%Soil_Moisture_Content
        CSEM_Land%LAI                   = Surface%Lai
        CSEM_Land%Vegetation_Type       = Surface%Vegetation_Type
    END IF


    IF( PRESENT(CSEM_Ice) ) THEN
        CSEM_Ice%Ice_Type               = Surface%Ice_Type
        CSEM_Ice%Ice_Temperature        = Surface%Ice_Temperature
        CSEM_Ice%Ice_Thickness          = Surface%Ice_Thickness
        CSEM_Ice%Ice_Density            = Surface%Ice_Density 
        CSEM_Ice%Ice_Roughness          = Surface%Ice_Roughness
        CSEM_Ice%Salinity               = Surface%Salinity

    END IF


    IF( PRESENT(CSEM_Water) ) THEN
        CSEM_Water%Water_Type           = Surface%Water_Type 
        CSEM_Water%Water_Temperature    = Surface%Water_Temperature
        CSEM_Water%Wind_Speed           = Surface%Wind_Speed
        CSEM_Water%Wind_Direction       = Surface%Wind_Direction
        CSEM_Water%Salinity             = Surface%Salinity
    END IF


   END FUNCTION CRTM_CSEM_Input_Adaptor

!===============================================================================
!
!===============================================================================
 FUNCTION CSEM_CRTM_Input_Adaptor( &
    Surface,     &  ! output
    CSEM_Land,   &  ! input
    CSEM_Water,  &  ! input
    CSEM_Snow,   &  ! input
    CSEM_Ice)    &  ! input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Surface_type)     , INTENT(INOUT)     :: Surface
   
    TYPE(CSEM_Land_Surface) ,OPTIONAL, INTENT(IN)  :: CSEM_Land
    TYPE(CSEM_Water_Surface),OPTIONAL, INTENT(IN)  :: CSEM_Water
    TYPE(CSEM_Snow_Surface) ,OPTIONAL, INTENT(IN)  :: CSEM_Snow
    TYPE(CSEM_Ice_Surface)  ,OPTIONAL, INTENT(IN)  :: CSEM_Ice
 
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_CRTM_Input_Adaptor'
    ! Local variables
    CHARACTER(ML) :: Message

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS
 
    IF( PRESENT(CSEM_Land) ) THEN
        Surface%Soil_Moisture_Content  =   CSEM_Land%Top_Soil_Moisture     
        Surface%Vegetation_Fraction    =   CSEM_Land%vegetation_Fraction  
        Surface%Soil_Temperature       =   CSEM_Land%Top_Soil_Temperature  
        Surface%Land_Temperature       =   CSEM_Land%Land_Skin_Temperature 
        Surface%Lai                    =   CSEM_Land%LAI                  
        Surface%Soil_Type              =   CSEM_Land%soil_Type            
        Surface%Vegetation_Type        =   CSEM_Land%Vegetation_Type  
        Surface%Land_Type              =   CSEM_Land%Land_Cover_Type 

    END IF
     
    IF( PRESENT(CSEM_Snow) ) THEN
   
        Surface%Snow_Type              =   CSEM_Snow%Snow_Type       
        Surface%Snow_Temperature       =   CSEM_Snow%Snow_Temperature     
        Surface%Snow_Depth             =   CSEM_Snow%Snow_Depth          
        Surface%Snow_Density           =   CSEM_Snow%Snow_Density          
        Surface%Snow_Grain_Size        =   CSEM_Snow%Snow_Grain_Size      
        Surface%Soil_Type              =   CSEM_Snow%soil_Type            
        Surface%Soil_Temperature       =   CSEM_Snow%Top_Soil_Temperature  
        Surface%Soil_Moisture_Content  =   CSEM_Snow%Top_Soil_Moisture_Content 
        Surface%Lai                    =   CSEM_Land%LAI                  
        Surface%Vegetation_Type        =   CSEM_Land%Vegetation_Type  
    END IF


    IF( PRESENT(CSEM_Ice) ) THEN
        Surface%Ice_Type              =    CSEM_Ice%Ice_Type  
        Surface%Ice_Temperature       =    CSEM_Ice%Ice_Temperature  
        Surface%Ice_Thickness         =    CSEM_Ice%Ice_Thickness  
        Surface%Ice_Density           =    CSEM_Ice%Ice_Density  
        Surface%Ice_Roughness         =    CSEM_Ice%Ice_Roughness  
        Surface%Salinity              =    CSEM_Ice%Salinity 

    END IF


    IF( PRESENT(CSEM_Water) ) THEN
        Surface%Water_Type            =    CSEM_Water%Water_Type    
        Surface%Water_Temperature     =    CSEM_Water%Water_Temperature   
        Surface%Wind_Speed            =    CSEM_Water%Wind_Speed  
        Surface%Wind_Direction        =    CSEM_Water%Wind_Direction 
        Surface%Salinity              =    CSEM_Water%Salinity            
    END IF


   END FUNCTION CSEM_CRTM_Input_Adaptor

END MODULE CRTM_SfcOptics
