!
! Common_RTSolution
!
! Module containing common ancillary functions for
! available radiative transfer algorithms.
!
! CREATION HISTORY:
!       Written by:     David Neil Groff,    IMSG at EMC;    david.groff@noaa.gov
!                       09-Sep-2010


MODULE Common_RTSolution

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE Type_Kinds,                ONLY: fp
  USE CRTM_Parameters,           ONLY: ONE, ZERO, TWO, PI, &
                                       DEGREES_TO_RADIANS, &
                                       SECANT_DIFFUSIVITY, &
                                       SCATTERING_ALBEDO_THRESHOLD
  USE Message_Handler,           ONLY: SUCCESS, Display_Message
  USE CRTM_Atmosphere_Define,    ONLY: CRTM_Atmosphere_type
  USE CRTM_Surface_Define,       ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define,  ONLY: CRTM_GeometryInfo_type
  USE CRTM_Planck_Functions,     ONLY: CRTM_Planck_Radiance , &
                                       CRTM_Planck_Temperature, &
                                       CRTM_Planck_Radiance_TL, &
                                       CRTM_Planck_Temperature_TL, &
                                       CRTM_Planck_Temperature_AD, &
                                       CRTM_Planck_Radiance_AD
  USE CRTM_SpcCoeff
  USE CRTM_AtmOptics_Define,     ONLY: CRTM_AtmOptics_type
  USE CRTM_AtmOptics,            ONLY: CRTM_Include_Scattering
  USE RTV_Define
  USE CRTM_SfcOptics
  USE CRTM_SfcOptics_Define
  USE CRTM_Utility
  USE CRTM_RTSolution_Define

  ! Disable all implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE

  PUBLIC :: Assign_Common_Input
  PUBLIC :: Assign_Common_Output
  PUBLIC :: Assign_Common_Input_TL
  PUBLIC :: Assign_Common_Output_TL
  PUBLIC :: Assign_Common_Input_AD
  PUBLIC :: Assign_Common_Output_AD

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: $'

CONTAINS

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!-------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Input
!
! PURPOSE:
!       Function to assign input that is used when calling both the ADA and
!       SOI RT algorithms.
!
! CALLING SEQUENCE:
!       Error_Status = Assign_Common_Input( Atmosphere  , &  ! Input
!                                           Surface     , &  ! Input
!                                           AtmOptics   , &  ! Input
!                                           SfcOptics   , &  ! Input
!                                           GeometryInfo, &  ! Input
!                                           SensorIndex , &  ! Input
!                                           ChannelIndex, &  ! Input
!                                           RTSolution  , &  ! Output
!                                           RTV           )  ! Internal variable output
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
!                       ATTRIBUTES: INTENT(IN OUT)
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
!       nz              Assignment of RTV%n_Angles for
!                       existing RTSolution interfaces
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
!
!       RTV:            Structure containing internal variables required for
!                       tangent-linear or adjoint model calls in CRTM_RTSolution.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
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

  FUNCTION Assign_Common_Input( &
    Atmosphere  , &  ! Input
    Surface     , &  ! Input
    AtmOptics   , &  ! Input
    SfcOptics   , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    RTSolution  , &  ! Output
    nz          , &  ! Input
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
    INTEGER,                      INTENT(OUT)    :: nz
    TYPE(RTV_type),               INTENT(IN OUT) :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Input'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: no, na, nt
    INTEGER :: i, j, k
    REAL(fp) :: Factor  ! SfcOptics quadrature weights normalisation factor
    REAL(fp) :: User_Emissivity, Direct_Reflectivity

    ! Set Up
    Error_Status = SUCCESS
    GeometryInfo%Cosine_Sensor_Zenith = ONE / GeometryInfo%Secant_Sensor_Zenith

    RTV%n_Added_Layers = Atmosphere%n_Added_Layers
    RTV%n_Layers       = Atmosphere%n_Layers
    RTV%n_Streams      = RTSolution%n_Full_Streams/2

    ! Save the optical depth if required
    IF ( CRTM_RTSolution_Associated( RTSolution ) ) THEN
      ! Shorter names for indexing
      no = RTSolution%n_Layers  ! Original no. of layers
      na = RTV%n_Added_Layers   ! No. of added layers
      nt = RTV%n_Layers         ! Current total no. of layers
      ! Assign only the optical depth profile
      ! defined by the user input layering
      RTSolution%Layer_Optical_Depth(1:no) = AtmOptics%Optical_Depth(na+1:nt)
    END IF
    ! Required SpcCoeff components
    RTV%Cosmic_Background_Radiance = SC(SensorIndex)%Cosmic_Background_Radiance(ChannelIndex)
    RTV%Sensor_Type                = SC(SensorIndex)%Sensor_Type
    RTV%Is_Solar_Channel           = SpcCoeff_IsSolar( SC(SensorIndex), ChannelIndex=ChannelIndex )
    RTV%COS_SUN                = cos(GeometryInfo%Source_Zenith_Radian)
    RTV%Solar_Irradiance       = SC(SensorIndex)%Solar_Irradiance(ChannelIndex) * GeometryInfo%AU_ratio2
    !
    RTSolution%Solar_Irradiance = RTV%Solar_Irradiance

    ! Determine the surface emission behavior
    !   By default, surface is SPECULAR.
    RTV%Diffuse_Surface = .FALSE.
    IF ( SpcCoeff_IsVisibleSensor( SC(SensorIndex)).or.SpcCoeff_IsUltravioletSensor( SC(SensorIndex)) ) THEN
      RTV%Diffuse_Surface = .TRUE.
    END IF
    ! -------------------------------------------
    ! Determine the quadrature angles and weights
    ! -------------------------------------------
    ! Default is to assume scattering RT
    RTV%Scattering_RT = .FALSE.

    ! Check for scattering conditions
    Determine_Stream_Angles: IF( RTSolution%Scattering_FLAG .AND. CRTM_Include_Scattering(AtmOptics) ) THEN

      ! --------------------------
      ! Set the scattering RT flag
      ! --------------------------
      RTV%Scattering_RT = .TRUE.

      ! ---------------------------------------
      ! Compute the Gaussian angles and weights
      ! ---------------------------------------
      CALL Double_Gauss_Quadrature( RTV%n_Streams, RTV%COS_Angle, RTV%COS_Weight )

      ! ----------------------------------------------------------------
      ! Is an additional stream required for the satellite zenith angle?
      ! That is, is the satellite zenith angle sufficiently different
      ! from the stream angles?
      ! ----------------------------------------------------------------
      ! Search for a matching angle
      DO i = 1, RTV%n_Streams
        IF( ABS( RTV%COS_Angle(i) - GeometryInfo%Cosine_Sensor_Zenith ) < ANGLE_THRESHOLD ) THEN
          SfcOptics%Index_Sat_Ang = i         ! Identify the sensor zenith angle
          RTV%n_Angles = RTV%n_Streams
          ! *****FLAW*****
          ! CAN WE REPLACE THIS CONSTRUCT WITH ANOTHER
          GO TO 100
          ! *****FLAW*****
        END IF
      END DO

      ! No match found, so flag an additional
      ! spot for the satellite zenith angle
      RTV%n_Angles = RTV%n_Streams + 1
      SfcOptics%Index_Sat_Ang = RTV%n_Angles      ! Identify the sensor zenith angle
      RTV%COS_Angle( RTV%n_Angles )  = GeometryInfo%Cosine_Sensor_Zenith
      RTV%COS_Weight( RTV%n_Angles ) = ZERO

      100 CONTINUE

      ! ------------------------
      ! Compute the phase matrix
      ! ------------------------
      CALL CRTM_Phase_Matrix( AtmOptics, & ! Input
                              RTV        ) ! Internal variable

    ELSE IF( RTV%Diffuse_Surface) THEN

      ! ------------------
      ! Lambertian surface
      ! ------------------
      ! The default diffusivity angle
      RTV%n_Streams       = 1
      RTV%COS_Angle( 1 )  = ONE/SECANT_DIFFUSIVITY
      RTV%COS_Weight( 1 ) = ONE
      ! The satellite zenith angle
      RTV%n_Angles                   = RTV%n_Streams + 1
      SfcOptics%Index_Sat_Ang        = RTV%n_Angles      ! Identify the sensor zenith angle
      RTV%COS_Angle( RTV%n_Angles )  = GeometryInfo%Cosine_Sensor_Zenith
      RTV%COS_Weight( RTV%n_Angles ) = ZERO

    ELSE

      ! ----------------
      ! Specular surface
      ! ----------------
      RTV%n_Streams = 0
      RTV%n_Angles  = 1
      RTV%COS_Angle( RTV%n_Angles )  = GeometryInfo%Cosine_Sensor_Zenith
      RTV%COS_Weight( RTV%n_Angles ) = ONE
      SfcOptics%Index_Sat_Ang        = RTV%n_Angles

    END IF Determine_Stream_Angles
    ! -------------------------------------------------
    ! Copy the RTV angles to a short name for use below
    ! -------------------------------------------------
    nZ = RTV%n_Angles

    ! ---------------------------------------------------------------------
    ! Assign the number of Stokes parameters. Currently, this is set to 1
    ! for decoupled polarization between surface and atmosphere.
    ! Remember for polarised microwave instruments, each frequency's Stokes
    ! parameter is treated as a separate channel.
    ! ---------------------------------------------------------------------
    ! ----------------------------------------------------
    ! Assign the angles and weights to SfcOptics structure
    ! ----------------------------------------------------
    SfcOptics%n_Angles = RTV%n_Angles
    Factor = ZERO
    k = 0
    DO i = 1, nZ
      SfcOptics%Angle(i)  = ACOS( RTV%COS_Angle(i) ) / DEGREES_TO_RADIANS
      SfcOptics%Weight(i) = RTV%COS_Angle(i) * RTV%COS_Weight(i)
      ! Compute the weight normalisation factor
      Factor = Factor + SfcOptics%Weight(i)
      DO j = 1, RTV%n_Stokes
        k = k + 1
        RTV%COS_WeightS(k) = RTV%COS_Weight(i)
        RTV%COS_AngleS(k) = RTV%COS_Angle(i)
      END DO
    END DO
    ! Normalise the quadrature weights
    IF( Factor > ZERO) SfcOptics%Weight(1:nZ) = SfcOptics%Weight(1:nZ) / Factor
    ! --------------------------------
    ! Populate the SfcOptics structure
    ! --------------------------------
    IF ( SfcOptics%Compute ) THEN
!    print *,Surface%Water_Coverage,Surface%Land_Coverage,Surface%Snow_Coverage,Surface%Ice_Coverage
      Error_Status = CRTM_Compute_SfcOptics( &
                       Surface     , & ! Input
                       GeometryInfo, & ! Input
                       SensorIndex , & ! Input
                       ChannelIndex, & ! Input
                       SfcOptics   , & ! In/Output
                       RTV%SOV       ) ! Internal variable output

      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message,'("Error computing SfcOptics for ",a," channel ",i0)' ) &
               TRIM(SC(SensorIndex)%Sensor_Id), &
               SC(SensorIndex)%Sensor_Channel(ChannelIndex)
        CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
        RETURN
      END IF

    ELSE
      IF( RTV%Scattering_RT ) THEN
        ! Replicate the user emissivity for all angles
        SfcOptics%Reflectivity = ZERO
        User_Emissivity = SfcOptics%Emissivity(1,1)
        SfcOptics%Emissivity(1,1) = ZERO
        Direct_Reflectivity = SfcOptics%Direct_Reflectivity(1,1)
        SfcOptics%Emissivity(1:nZ,1) = User_Emissivity
        ! Replicate the user reflectivities for all angles

        SfcOptics%Direct_Reflectivity(1:nZ,1) = Direct_Reflectivity
        IF( RTV%Diffuse_Surface) THEN
          IF( RTV%mth_Azi == 0 ) THEN
          DO i = 1, nZ
            SfcOptics%Reflectivity(1:nZ, 1, i, 1) = (ONE-SfcOptics%Emissivity(i,1))*SfcOptics%Weight(i)
          END DO
          ELSE
            SfcOptics%Reflectivity(1:nZ, 1, 1:nZ, 1) = ZERO
            SfcOptics%Direct_Reflectivity(1:nZ,1) = ZERO
            SfcOptics%Emissivity(1:nZ,1) = ZERO
          END IF


        ELSE ! Specular surface
          DO i = 1, nZ
            SfcOptics%Reflectivity(i, 1, i, 1) = (ONE-SfcOptics%Emissivity(i,1))
          END DO
        END IF
      ELSE
          IF( RTV%mth_Azi == 0 ) THEN
            User_Emissivity = SfcOptics%Emissivity(1,1)
            SfcOptics%Emissivity( SfcOptics%Index_Sat_Ang,1 ) = User_Emissivity
            SfcOptics%Reflectivity(1,1,1,1) = ONE - User_Emissivity
            Direct_Reflectivity = SfcOptics%Direct_Reflectivity(1,1)
            SfcOptics%Direct_Reflectivity(1,1) = Direct_Reflectivity
          ELSE
            User_Emissivity = ZERO
            SfcOptics%Emissivity( SfcOptics%Index_Sat_Ang,1 ) = ZERO
            SfcOptics%Reflectivity(1,1,1,1) = ZERO
            SfcOptics%Direct_Reflectivity(1,1) = ZERO
          END IF
      END IF

    END IF

    ! ---------------------------
    ! Save the surface properties
    ! ---------------------------
    RTSolution%Surface_Emissivity   = SfcOptics%Emissivity( SfcOptics%Index_Sat_Ang, 1 )
    RTSolution%Surface_Reflectivity = SfcOptics%Reflectivity( SfcOptics%Index_Sat_Ang, 1, SfcOptics%Index_Sat_Ang, 1 )

    ! ------------------------
    ! Compute Planck radiances
    ! ------------------------
    IF( RTV%mth_Azi == 0 ) THEN
      DO k = 1, Atmosphere%n_Layers
        CALL CRTM_Planck_Radiance( &
               SensorIndex              , & ! Input
               ChannelIndex             , & ! Input
               Atmosphere%Temperature(k), & ! Input
               RTV%Planck_Atmosphere(k)   ) ! Output
      END DO
      ! Surface radiance
      CALL CRTM_Planck_Radiance( &
             SensorIndex                  , & ! Input
             ChannelIndex                 , & ! Input
             SfcOptics%Surface_Temperature, & ! Input
             RTV%Planck_Surface             ) ! Output
    ELSE
      RTV%Planck_Atmosphere = ZERO
      RTV%Planck_Surface = ZERO
    END IF

  END FUNCTION Assign_Common_Input

!--------------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Input_TL
!
! PURPOSE:
!       Function to assign tangent-linear input common to both the ADA and
!       SOI RT models.
!
! CALLING SEQUENCE:
!      Error_Status = Assign_Common_Input_TL( Atmosphere   , &  ! FWD Input
!                                             Surface      , &  ! FWD Input
!                                             AtmOptics    , &  ! FWD Input
!                                             SfcOptics    , &  ! FWD Input
!                                             RTSolution   , &  ! FWD Input
!                                             Atmosphere_TL, &  ! TL Input
!                                             Surface_TL   , &  ! TL Input
!                                             AtmOptics_TL , &  ! TL Input
!                                             SfcOptics_TL , &  ! TL Input
!                                             GeometryInfo , &  ! Input
!                                             SensorIndex  , &  ! Input
!                                             ChannelIndex , &  ! Input
!                                             RTSolution_TL, &  ! TL Output
!                                             RTV            )  ! Internal variable input
!
! INPUT ARGUMENTS:
!       Atmosphere:             Structure containing the atmospheric state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Atmosphere_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       Surface:                Structure containing the surface state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Surface_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:              Structure containing the combined atmospheric
!                               optical properties for gaseous absorption, clouds,
!                               and aerosols.
!                               UNITS:      N/A
!                               TYPE:       CRTM_AtmOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       SfcOptics:              Structure containing the surface optical properties
!                               data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_SfcOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       Atmosphere_TL:          Structure containing the tangent-linear atmospheric
!                               state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Atmosphere_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       Surface_TL:             Structure containing the tangent-linear surface state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Surface_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       AtmOptics_TL:           Structure containing the tangent-linear atmospheric
!                               optical properties.
!                               UNITS:      N/A
!                               TYPE:       CRTM_AtmOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       SfcOptics_TL:           Structure containing the tangent-linear surface optical
!                               properties. Argument is defined as INTENT (IN OUT ) as
!                               different RT algorithms may compute the surface optics
!                               properties before this routine is called.
!                               UNITS:      N/A
!                               TYPE:       CRTM_SfcOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!       GeometryInfo:           Structure containing the view geometry data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_GeometryInfo_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:            Sensor index id. This is a unique index associated
!                               with a (supported) sensor used to access the
!                               shared coefficient data for a particular sensor.
!                               See the ChannelIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:           Channel index id. This is a unique index associated
!                               with a (supported) sensor channel used to access the
!                               shared coefficient data for a particular sensor's
!                               channel.
!                               See the SensorIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       RTV:                    Structure containing internal forward model variables
!                               required for subsequent tangent-linear or adjoint model
!                               calls. The contents of this structure are NOT accessible
!                               outside of the CRTM_RTSolution module.
!                               UNITS:      N/A
!                               TYPE:       RTV_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!
!
!       RTSolution_TL:          Structure containing the solution to the tangent-linear
!                               RT equation for the given inputs.
!                               UNITS:      N/A
!                               TYPE:       CRTM_RTSolution_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!       nz:                     Integer assigned from RTV%n_Angles
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT)
!
!       User_Emissivity_TL:     Tangent linear of emissivity
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT)
!
!       Direct_Reflectivity_TL: Tangent linear of direct reflectivity
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT)
!
!       Planck_Surface_TL:      Tangent linear of Planck Surface
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT)
!
!       Planck_Atmosphere_TL:   Tangent linear of Planck Atmosphere
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-1
!                               ATTRIBUTES: INTENT(OUT)
!
!       Pff_TL:                 Tangent linear of forward scattering
!                               phase function
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-3
!                               ATTRIBUTES: INTENT(OUT)
!
!       Pbb_TL:                 Tangent linear of backward scattering
!                               phase function
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-3
!                               ATTRIBUTES: INTENT(OUT)
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
  FUNCTION Assign_Common_Input_TL( &
    Atmosphere            , &  ! FWD Input
    Surface               , &  ! FWD Input
    AtmOptics             , &  ! FWD Input
    SfcOptics             , &  ! FWD Input
    Atmosphere_TL         , &  ! TL Input
    Surface_TL            , &  ! TL Input
    AtmOptics_TL          , &  ! TL Input
    SfcOptics_TL          , &  ! TL Input/Output
    GeometryInfo          , &  ! Input
    SensorIndex           , &  ! Input
    ChannelIndex          , &  ! Input
    RTSolution_TL         , &  ! TL Output
    nz                    , &  ! Output
    User_Emissivity_TL    , &  ! Output
    Direct_Reflectivity_TL, &  ! Output
    Planck_Surface_TL     , &  ! Output
    Planck_Atmosphere_TL  , &  ! Output
    Pff_TL                , &  ! Output
    Pbb_TL                , &  ! Output
    RTV                   ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_Atmosphere_type),   INTENT(IN)     :: Atmosphere_TL
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface_TL
    TYPE(CRTM_AtmOptics_type),    INTENT(IN)     :: AtmOptics_TL
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_TL
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution_TL
    INTEGER,                      INTENT(OUT)    :: nz
    REAL(fp),                     INTENT(OUT)    :: User_Emissivity_TL
    REAL(fp),                     INTENT(OUT)    :: Direct_Reflectivity_TL
    REAL(fp),                     INTENT(OUT)    :: Planck_Surface_TL
    REAL(fp),                     INTENT(OUT)    :: Planck_Atmosphere_TL(0:)
    REAL(fp),                     INTENT(OUT)    :: Pff_TL(:,:,:)
    REAL(fp),                     INTENT(OUT)    :: Pbb_TL(:,:,:)
    TYPE(RTV_type),               INTENT(IN)     :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Input_TL'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: i, k
    INTEGER :: no, na, nt

    ! ------
    ! Set up
    ! ------
    Error_Status = SUCCESS

    nz = RTV%n_Angles
    ! Set the sensor zenith angle index
    SfcOptics_TL%Index_Sat_Ang = SfcOptics%Index_Sat_Ang
    ! Save the TL optical depth if required
    IF ( CRTM_RTSolution_Associated( RTSolution_TL ) ) THEN
      ! Shorter names for indexing
      no = RTSolution_TL%n_Layers  ! Original no. of layers
      na = RTV%n_Added_Layers      ! No. of added layers
      nt = RTV%n_Layers            ! Current total no. of layers
      ! Assign only the optical depth profile
      ! defined by the user input layering
      RTSolution_TL%Layer_Optical_Depth(1:no) = AtmOptics_TL%Optical_Depth(na+1:nt)
    END IF

    ! ---------------------------------------------------
    ! Compute the tangent-linear phase matrix if required
    ! ---------------------------------------------------
    IF ( RTV%Scattering_RT ) THEN
      CALL CRTM_Phase_Matrix_TL( &
             AtmOptics,    & ! Input
             AtmOptics_TL, & ! Input
             Pff_TL,       & ! Output - TL forward scattering
             Pbb_TL,       & ! Output - TL backward scattering
             RTV           ) ! Internal variable
    END IF


    ! ---------------------------------------------------------------------
    ! Assign the number of Stokes parameters. Currently, this is set to 1
    ! for decoupled polarization between surface and atmosphere.
    ! Remember for polarised microwave instruments, each frequency's Stokes
    ! parameter is treated as a separate channel.
    ! ---------------------------------------------------------------------
!    SfcOptics_TL%n_Stokes = 1

    ! ----------------------------------------------------
    ! Assign the angles and weights to SfcOptics structure
    ! ----------------------------------------------------
    SfcOptics_TL%n_Angles = SfcOptics%n_Angles
    SfcOptics_TL%Angle    = SfcOptics%Angle
    SfcOptics_TL%Weight   = SfcOptics%Weight


    ! -----------------------------------------
    ! Populate the SfcOptics_TL structure based
    ! on FORWARD model SfcOptics Compute_Switch
    ! -----------------------------------------
    IF ( SfcOptics%Compute ) THEN
      Error_Status = CRTM_Compute_SfcOptics_TL( &
                       Surface     , & ! Input
                       SfcOptics   , & ! Input
                       Surface_TL  , & ! Input
                       GeometryInfo, & ! Input
                       SensorIndex , & ! Input
                       ChannelIndex, & ! Input
                       SfcOptics_TL, & ! In/Output
                       RTV%SOV       ) ! Internal variable input
      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message,'("Error computing SfcOptics_TL for ",a," channel ",i0)' ) &
               TRIM(SC(SensorIndex)%Sensor_Id), &
               SC(SensorIndex)%Sensor_Channel(ChannelIndex)
        CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
        RETURN
      END IF
    ELSE
      IF( RTV%Scattering_RT ) THEN
        ! Replicate the user emissivity for all angles
        SfcOptics_TL%Reflectivity       = ZERO
        User_Emissivity_TL              = SfcOptics_TL%Emissivity(1,1)
        SfcOptics_TL%Emissivity(1,1)    = ZERO
        Direct_Reflectivity_TL          = SfcOptics_TL%Direct_Reflectivity(1,1)/PI
        SfcOptics_TL%Emissivity(1:nZ,1) = User_Emissivity_TL
        ! Replicate the user reflectivities for all angles
        SfcOptics_TL%Direct_Reflectivity(1:nZ,1) = Direct_Reflectivity_TL
        IF( RTV%Diffuse_Surface) THEN

          IF( RTV%mth_Azi == 0 ) THEN
          DO i = 1, nZ
            SfcOptics_TL%Reflectivity(1:nZ, 1, i, 1) = -SfcOptics_TL%Emissivity(i,1)*SfcOptics%Weight(i)
          END DO
          ELSE
            SfcOptics_TL%Reflectivity(1:nZ, 1, 1:nZ, 1) = ZERO
            SfcOptics_TL%Direct_Reflectivity(1:nZ,1) = ZERO
            SfcOptics_TL%Emissivity(1:nZ,1) = ZERO
          END IF

        ELSE ! Specular surface
          DO i = 1, nZ
            SfcOptics_TL%Reflectivity(i, 1, i, 1) = -SfcOptics_TL%Emissivity(i,1)
          END DO
        END IF
      ELSE
        IF( RTV%mth_Azi == 0 ) THEN
          User_Emissivity_TL = SfcOptics_TL%Emissivity(1,1)
          SfcOptics_TL%Emissivity( SfcOptics%Index_Sat_Ang,1 ) = User_Emissivity_TL
          SfcOptics_TL%Reflectivity(1,1,1,1) = - User_Emissivity_TL
        ELSE
          User_Emissivity_TL = ZERO
          SfcOptics_TL%Emissivity( SfcOptics%Index_Sat_Ang,1 ) = ZERO
          SfcOptics_TL%Reflectivity(1,1,1,1) = ZERO
        END IF

      END IF
    END IF

    ! ------------------------------
    ! Save the TL surface properties
    ! ------------------------------
    RTSolution_TL%Surface_Emissivity   = SfcOptics_TL%Emissivity( SfcOptics_TL%Index_Sat_Ang, 1 )
    RTSolution_TL%Surface_Reflectivity = SfcOptics_TL%Reflectivity( SfcOptics_TL%Index_Sat_Ang, 1, &
                                                                    SfcOptics_TL%Index_Sat_Ang, 1  )

    ! -------------------------------------------
    ! Compute the tangent-linear planck radiances
    ! -------------------------------------------
    IF( RTV%mth_Azi == 0 ) THEN
      ! Atmospheric layer TL radiances
      DO k = 1, Atmosphere%n_Layers
        CALL CRTM_Planck_Radiance_TL( &
               SensorIndex                 , & ! Input
               ChannelIndex                , & ! Input
               Atmosphere%Temperature(k)   , & ! Input
               Atmosphere_TL%Temperature(k), & ! Input
               Planck_Atmosphere_TL(k)       ) ! Output
      END DO
      ! Surface TL radiance
      CALL CRTM_Planck_Radiance_TL( &
             SensorIndex                     , & ! Input
             ChannelIndex                    , & ! Input
             SfcOptics%Surface_Temperature   , & ! Input
             SfcOptics_TL%Surface_Temperature, & ! Input
             Planck_Surface_TL                 ) ! Output
    ELSE
      Planck_Atmosphere_TL = ZERO
      Planck_Surface_TL = ZERO
    END IF

  END FUNCTION Assign_Common_Input_TL



!--------------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Input_AD
!
! PURPOSE:
!       Function to assign common input before calling the adjoint RT algorithms.
!
! CALLING SEQUENCE:
!      Error_Status = Assign_Common_Input_AD( SfcOptics    , &  ! FWD Input
!                                             RTSolution   , &  ! FWD Input
!                                             GeometryInfo , &  ! Input
!                                             SensorIndex  , &  ! Input
!                                             ChannelIndex , &  ! Input
!                                             RTSolution_AD, &  ! AD Input/Output
!                                             Atmosphere_AD, &  ! AD Output
!                                             Surface_AD   , &  ! AD Output
!                                             AtmOptics_AD , &  ! AD Output
!                                             SfcOptics_AD , &  ! AD Output
!                                             RTV            )  ! Internal variable input
!
! INPUT ARGUMENTS:
!
!       SfcOptics:            Structure containing the surface optical properties
!                             data.
!                             UNITS:      N/A
!                             TYPE:       CRTM_SfcOptics_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       RTSolution:           Structure containing the solution to the RT equation
!                             for the given inputs.
!                             UNITS:      N/A
!                             TYPE:       CRTM_RTSolution_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:         Structure containing the view geometry data.
!                             UNITS:      N/A
!                             TYPE:       CRTM_GeometryInfo_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:          Sensor index id. This is a unique index associated
!                             with a (supported) sensor used to access the
!                             shared coefficient data for a particular sensor.
!                             See the ChannelIndex argument.
!                             UNITS:      N/A
!                             TYPE:       INTEGER
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:         Channel index id. This is a unique index associated
!                             with a (supported) sensor channel used to access the
!                             shared coefficient data for a particular sensor's
!                             channel.
!                             See the SensorIndex argument.
!                             UNITS:      N/A
!                             TYPE:       INTEGER
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
!       RTV:                  Structure containing internal forward model variables
!                             required for subsequent tangent-linear or adjoint model
!                             calls. The contents of this structure are NOT accessible
!                             outside of the CRTM_RTSolution module.
!                             UNITS:      N/A
!                             TYPE:       RTV_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!
!       RTSolution_AD:        Structure containing the RT solution adjoint inputs.
!                             UNITS:      N/A
!                             TYPE:       CRTM_RTSolution_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
!       SfcOptics_AD:         Structure containing the adjoint surface optical
!                             properties data.
!                             UNITS:      N/A
!                             TYPE:       CRTM_SfcOptics_type
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(IN OUT)
!
!       Planck_Surface_AD:    Adjoint of Planck surface.
!                             UNITS:      N/A
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(OUT)
!
!       Planck_Atmosphere_AD: Adjoint of Planck atmosphere
!                             UNITS:      N/A
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Rank-1 (0:n_Layers)
!                             ATTRIBUTES: INTENT(OUT)
!
!       Radiance_AD:          Adjoint of Radiance
!                             UNITS:      N/A
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(OUT)
!
!       nz:                   Adjoint of Planck atmosphere
!                             UNITS:      N/A
!                             TYPE:       INTEGER
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: INTENT(OUT)
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

  FUNCTION Assign_Common_Input_AD( &
    SfcOptics            , &  ! FWD Input
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
    RTV                  ) &  ! Internal variable input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_RTSolution_type),   INTENT(IN)     :: RTSolution
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_RTSolution_type),   INTENT(IN OUT) :: RTSolution_AD
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    REAL(fp),                     INTENT(OUT)    :: Planck_Surface_AD
    REAL(fp),                     INTENT(OUT)    :: Planck_Atmosphere_AD(0:)
    REAL(fp), DIMENSION(:), INTENT(OUT)    :: Radiance_AD
    INTEGER,                      INTENT(OUT)    :: nz
    TYPE(RTV_type),               INTENT(IN)     :: RTV
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Input_AD'

    ! -----
    ! Setup
    ! -----
    Error_Status = SUCCESS
    ! Short named local copies
    nz = RTV%n_Angles
    ! Set the sensor zenith angle index
    SfcOptics_AD%Index_Sat_Ang = SfcOptics%Index_Sat_Ang

    ! Initialise local adjoint variables
    Planck_Surface_AD    = ZERO
    Planck_Atmosphere_AD = ZERO
    Radiance_AD = ZERO

    ! ------------------------------------------
    ! Compute the brightness temperature adjoint
    ! ------------------------------------------
    IF ( SpcCoeff_IsInfraredSensor( SC(SensorIndex) ) .OR. &
         SpcCoeff_IsMicrowaveSensor( SC(SensorIndex) ) ) THEN
      IF( RTV%mth_Azi == 0 ) THEN
        CALL CRTM_Planck_Temperature_AD( &
               SensorIndex                         , & ! Input
               ChannelIndex                        , & ! Input
               RTSolution%Radiance                 , & ! Input
               RTSolution_AD%Brightness_Temperature, & ! Input
               Radiance_AD(1)                        ) ! Output
        RTSolution_AD%Brightness_Temperature = ZERO
      END IF
    END IF

    ! accumulate Fourier component
    IF( RTV%n_Stokes == 1 ) THEN
      Radiance_AD(1) = Radiance_AD(1) + RTSolution_AD%Radiance * &
         COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
    ELSE
      Radiance_AD(1:2) = Radiance_AD(1:2) + RTSolution_AD%Stokes(1:2) * &
         COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )

      IF(  RTV%n_Stokes > 2 ) THEN
       Radiance_AD(3:RTV%n_Stokes) = Radiance_AD(3:RTV%n_Stokes) + RTSolution_AD%Stokes(3:RTV%n_Stokes)* &
       SIN( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
      END IF
    END IF


    ! ---------------------------------------------------------------
    ! Compute the corresponding tangent-linear brightness temperature
    ! ---------------------------------------------------------------

  END FUNCTION Assign_Common_Input_AD

!-------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Output
!
! PURPOSE:
!       Function to assign output from CRTM_ADA, CRTM_SOI and CRTM_Emission
!
! CALLING SEQUENCE:
!       Error_Status = Assign_Common_Output( Atmosphere   , & ! Input
!                                            SfcOptics    , & ! Input
!                                            GeometryInfo , & ! Input
!                                            SensorIndex  , & ! Input
!                                            ChannelIndex , & ! Input
!                                            RTV          , & ! Input
!                                            RTSolution     ) ! Output
!
! INPUT ARGUMENTS:
!
!       Atmosphere:     Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
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
!                       ATTRIBUTES: INTENT(IN)
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
!       RTV:            Structure containing internal variables required for
!                       subsequent tangent-linear or adjoint model calls.
!                       The contents of this structure are NOT accessible
!                       outside of the CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
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

  FUNCTION Assign_Common_Output( &
    Atmosphere   , & ! Input
    SfcOptics    , & ! Input
    GeometryInfo , & ! Input
    SensorIndex  , & ! Input
    ChannelIndex , & ! Input
    RTV          , & ! Input
    RTSolution   ) & ! Output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type)  , INTENT(IN)     :: Atmosphere
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN)     :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    TYPE(RTV_type)              , INTENT(IN)     :: RTV
    TYPE(CRTM_RTSolution_type)  , INTENT(IN OUT) :: RTSolution

    ! Function Result
    INTEGER :: Error_Status
    ! Local Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Output'
    ! Local variables
    INTEGER :: no, na, nt, n1
    REAL(fp) :: Radiance(RTV%n_Stokes)

    Error_Status = SUCCESS
    n1 = (SfcOptics%Index_Sat_Ang-1)*RTV%n_Stokes + 1
    ! ADA and SOI specific assignments
    IF( RTV%Scattering_RT ) THEN

      ! Assign radiance output from ADA or SOI
      IF ( RTV%aircraft%rt ) THEN
        Radiance(:) = RTV%s_Level_Rad_UP(n1:n1-1+RTV%n_Stokes, RTV%aircraft%idx)
      ELSE
        Radiance(:) = RTV%s_Level_Rad_UP(n1:n1-1+RTV%n_Stokes, 0)
      END IF

      ! Output downwelling radiance
      IF ( RTV%obs_4_downward%rt ) THEN
        Radiance = RTV%s_Level_Rad_DOWN(n1:n1-1+RTV%n_Stokes, RTV%obs_4_downward%idx)
      END IF

    ! Emission specific assignments
    ELSE

      ! The output radiance at TOA or at Aircraft height
      IF ( RTV%aircraft%rt ) THEN
        Radiance(1) = RTV%e_Level_Rad_UP(RTV%aircraft%idx)
      ELSE
        Radiance(1) = RTV%e_Level_Rad_UP(0)
      END IF

      ! Output downwelling radiance
      IF ( RTV%obs_4_downward%rt ) THEN
        Radiance = RTV%e_Level_Rad_DOWN(RTV%obs_4_downward%idx)
      END IF

      ! Other emission-only output
      RTSolution%Up_Radiance             = RTV%Up_Radiance
      RTSolution%Down_Radiance           = RTV%e_Level_Rad_DOWN(Atmosphere%n_Layers)
      RTSolution%Down_Solar_Radiance     = RTV%Down_Solar_Radiance
      RTSolution%Surface_Planck_Radiance = RTV%Planck_Surface
      IF ( CRTM_RTSolution_Associated( RTSolution ) ) THEN
        ! Shorter names for indexing
        no = RTSolution%n_Layers  ! Original no. of layers
        na = RTV%n_Added_Layers   ! No. of added layers
        nt = RTV%n_Layers         ! Current total no. of layers
        ! Assign only the upwelling radiance profile
        ! defined by the user input layering
        RTSolution%Upwelling_Radiance(1:no) = RTV%e_Level_Rad_UP(na+1:nt)
        RTSolution%Upwelling_Overcast_Radiance(1:no) = RTV%e_Cloud_Radiance_UP(na+1:nt)
      END IF
    END IF

    ! accumulate Fourier component
    IF( RTV%n_Stokes == 1 ) THEN
     RTSolution%Radiance = RTSolution%Radiance + Radiance(1)*  &
       COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
      RTSolution%Stokes(1) = RTSolution%Radiance
    ELSE
     RTSolution%Stokes(1:2) = RTSolution%Stokes(1:2) + Radiance(1:2)*  &
       COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
       IF(  RTV%n_Stokes > 2 ) THEN
     RTSolution%Stokes(3:RTV%n_Stokes) = RTSolution%Stokes(3:RTV%n_Stokes) + Radiance(3:RTV%n_Stokes)*  &
       SIN( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
       END IF
       RTSolution%Radiance = RTSolution%Stokes(1)
    END IF

    ! ------------------------------------------------
    ! Compute the corresponding brightness temperature
    ! ------------------------------------------------
    IF( RTV%mth_Azi == 0 ) &
      CALL CRTM_Planck_Temperature( &
             SensorIndex,                      & ! Input
             ChannelIndex,                     & ! Input
             RTSolution%Radiance,              & ! Input
             RTSolution%Brightness_Temperature ) ! Output

    ! -------------------------------------------------------
    ! Compute the corresponding solar reflectance (TOA only)
    ! -------------------------------------------------------
    RTSolution%Reflectance = ZERO
    IF( RTV%Solar_Flag_true .AND. (.NOT.(RTV%aircraft%rt))) THEN
         RTSolution%Reflectance = RTSolution%Radiance*PI/RTV%Solar_irradiance
    END IF

    ! -------------------------------------------------------
    ! Same thing but assuming clear sky
    ! -------------------------------------------------------
    RTSolution%Reflectance_clear = ZERO
    IF( RTV%Solar_Flag_true .AND. (.NOT.(RTV%aircraft%rt))) THEN
         RTSolution%Reflectance_clear = RTSolution%R_clear*PI/RTV%Solar_irradiance
    END IF

  END FUNCTION Assign_Common_Output
!-------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Output_TL
!
! PURPOSE:
!       Function to assign output from CRTM_ADA, CRTM_SOI and CRTM_Emission
!
! CALLING SEQUENCE:
!       Error_Status = Assign_Common_Output_TL( SfcOptics              , & ! Input
!                                               RTSolution             , & ! Input
!                                               GeometryInfo           , & ! Input
!                                               Radiance_TL            , & ! Input
!                                               Scattering_Radiance_TL , & ! Input
!                                               SensorIndex            , & ! Input
!                                               ChannelIndex           , & ! Input
!                                               RTV                    , & ! Input
!                                               RTSolution_TL            ) ! Output
!
! INPUT ARGUMENTS:
!
!       SfcOptics:              Structure containing the surface optical properties
!                               data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_SfcOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       RTSolution:             Structure containing the soluition to the RT equation
!                               UNITS:      N/A
!                               TYPE:       CRTM_RTSolution_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:           Structure containing the view geometry data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_GeometryInfo_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       Radiance_TL:            TL radiance
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       Scattering_Radiance_TL: Scattering TL radiance
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-1
!                               ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:            Sensor index id. This is a unique index associated
!                               with a (supported) sensor used to access the
!                               shared coefficient data for a particular sensor.
!                               See the ChannelIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:           Channel index id. This is a unique index associated
!                               with a (supported) sensor channel used to access the
!                               shared coefficient data for a particular sensor's
!                               channel.
!                               See the SensorIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!       RTV:                    Structure containing internal variables required for
!                               subsequent tangent-linear or adjoint model calls.
!                               The contents of this structure are NOT accessible
!                               outside of the CRTM_RTSolution module.
!                               UNITS:      N/A
!                               TYPE:       RTV_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution_TL:          Structure containing the TL soluition to the RT equation.
!                               UNITS:      N/A
!                               TYPE:       CRTM_RTSolution_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:           The return value is an integer defining the error status.
!                               The error codes are defined in the Message_Handler module.
!                               If == SUCCESS the computation was sucessful
!                               == FAILURE an unrecoverable error occurred
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output RTSolution argument is IN OUT rather than
!       just OUT. This is necessary because the argument is defined upon
!       input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!--------------------------------------------------------------------------------

  FUNCTION Assign_Common_Output_TL( &
    SfcOptics              , & ! Input
    RTSolution             , & ! Input
    GeometryInfo           , & ! Input
    Radiance_TL            , & ! Input
    Scattering_Radiance_TL , & ! Input
    SensorIndex            , & ! Input
    ChannelIndex           , & ! Input
    RTV                    , & ! Input
    RTSolution_TL          ) & ! Output
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN)     :: SfcOptics
    TYPE(CRTM_RTSolution_type)  , INTENT(IN)     :: RTSolution
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    REAL(fp)                    , INTENT(IN)     :: Radiance_TL
    REAL(fp)                    , INTENT(IN)     :: Scattering_Radiance_TL(:)
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    TYPE(RTV_type)              , INTENT(IN)     :: RTV
    TYPE(CRTM_RTSolution_type)  , INTENT(IN OUT) :: RTSolution_TL

    ! Function Result
    INTEGER :: Error_Status,n1
    REAL(fp) :: SRadiance_TL(RTV%n_Stokes)
    ! Local Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Output_TL'

    Error_Status = SUCCESS

    n1 = (SfcOptics%Index_Sat_Ang-1)*RTV%n_Stokes + 1
    ! No TL/AD for the simulations at aircraft
    IF( RTV%Scattering_RT ) THEN
        SRadiance_TL(:) = Scattering_Radiance_TL(n1:n1-1+RTV%n_Stokes)
    ! Emission specific assignments
    ELSE
        SRadiance_TL(1) = Radiance_TL
    END IF

    ! accumulate Fourier component
    IF( RTV%n_Stokes == 1 ) THEN
      RTSolution_TL%Radiance = RTSolution_TL%Radiance + SRadiance_TL(1)*  &
         COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
      RTSolution_TL%Stokes(1) = RTSolution_TL%Radiance
    ELSE
      RTSolution_TL%Stokes(1:2) = RTSolution_TL%Stokes(1:2) + SRadiance_TL(1:2)*  &
         COS( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
      RTSolution_TL%Radiance = RTSolution_TL%Stokes(1)
      IF(  RTV%n_Stokes > 2 ) THEN
      RTSolution_TL%Stokes(3:RTV%n_Stokes) = RTSolution_TL%Stokes(3:RTV%n_Stokes) + SRadiance_TL(3:RTV%n_Stokes)* &
       SIN( RTV%mth_Azi*(GeometryInfo%Sensor_Azimuth_Radian-GeometryInfo%Source_Azimuth_Radian) )
     END IF
    END IF


    ! ---------------------------------------------------------------
    ! Compute the corresponding tangent-linear brightness temperature
    ! ---------------------------------------------------------------
    IF( RTV%mth_Azi == 0 ) &
      CALL CRTM_Planck_Temperature_TL( &
             SensorIndex                        , & ! Input
             ChannelIndex                       , & ! Input
             RTSolution%Radiance                , & ! Input
             RTSolution_TL%Radiance             , & ! Input
             RTSolution_TL%Brightness_Temperature ) ! Output

  END FUNCTION Assign_Common_Output_TL

!------------------------------------------------------------------------------
!
! NAME:
!       Assign_Common_Output_AD
!
! PURPOSE:
!       Function to assign output from CRTM_ADA, CRTM_SOI and CRTM_Emission
!       adjoint calls
!
! CALLING SEQUENCE:
!       Error_Status = Assign_Common_Output_AD( Atmosphere           , & ! Input
!                                               Surface              , & ! Input
!                                               AtmOptics            , & ! Input
!                                               SfcOptics            , & ! Input
!                                               Pff_AD               , & ! Input
!                                               Pbb_AD               , & ! Input
!                                               RTSolution           , & ! Input
!                                               GeometryInfo         , & ! Input
!                                               SensorIndex          , & ! Input
!                                               ChannelIndex         , & ! Input
!                                               nZ                   , & ! Input
!                                               AtmOptics_AD         , & ! Output
!                                               SfcOptics_AD         , & ! Output
!                                               Planck_Surface_AD    , & ! Output
!                                               Planck_Atmosphere_AD , & ! Output
!                                               User_Emissivity_AD   , & ! Output
!                                               Atmosphere_AD        , & ! Output
!                                               Surface_AD           , & ! Output
!                                               RTSolution_AD        , & ! Output
!                                               RTV                  ) & ! Input
!
! INPUT ARGUMENTS:
!
!      Atmosphere:              Structure containing the atmospheric state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Atmosphere_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      Surface:                 Structure containing the Surface state data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_Surface_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      AtmOptics:               Structure containing the combined atmospheric
!                               optical properties for gaseous absorption, clouds,
!                               and aerosols.
!                               UNITS:      N/A
!                               TYPE:       CRTM_AtmOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      SfcOptics:               Structure containing the surface
!                               optical properties data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_SfcOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      GeometryInfo:            Structure containing the view geometry data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_GeometryInfo_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      SensorIndex:             Sensor index id. This is a unique index associated
!                               with a (supported) sensor used to access the
!                               shared coefficient data for a particular sensor.
!                               See the ChannelIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      ChannelIndex:            Channel index id. This is a unique index associated
!                               with a (supported) sensor channel used to access the
!                               shared coefficient data for a particular sensor's
!                               channel.
!                               See the SensorIndex argument.
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      nz:                      Integer assigned from RTV%n_Angles
!                               UNITS:      N/A
!                               TYPE:       INTEGER
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!      RTV:                     Structure containing internal variables required for
!                               subsequent tangent-linear or adjoint model calls.
!                               The contents of this structure are NOT accessible
!                               outside of the CRTM_RTSolution module.
!                               UNITS:      N/A
!                               TYPE:       RTV_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN)
!
!  OUTPUT ARGUMENTS
!
!      Pff_AD:                  Forward scattering AD phase matrix
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      Pbb_AD:                  Backward scattering AD phase matrix
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      AtmOptics_AD:            Structure containing the adjoint combined atmospheric
!                               optical properties for gaseous absorption, clouds,
!                               and aerosols.
!                               UNITS:      N/A
!                               TYPE:       CRTM_AtmOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      SfcOptics_AD:            Structure containing the adjoint surface optical
!                               properties data.
!                               UNITS:      N/A
!                               TYPE:       CRTM_SfcOptics_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      Planck_Surface_AD:       Adjoint of Planck surface.
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      Planck_Atmosphere_AD:    Adjoint of Planck atmosphere
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Rank-1 (0:n_Layers)
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      User_Emissivity_AD:      Adjoint of user emissivity
!                               UNITS:      N/A
!                               TYPE:       REAL(fp)
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(OUT)
!
!      Atmosphere_AD:           Adjoint of atmosphere
!                               UNITS:      N/A
!                               TYPE:       CRTM_Atmosphere_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      Surface_AD:              Adjoint of surface
!                               UNITS:      N/A
!                               TYPE:       CRTM_Surface_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
!
!      RTSolution_AD:           Structure containing the RT solution adjoint inputs.
!                               UNITS:      N/A
!                               TYPE:       CRTM_RTSolution_type
!                               DIMENSION:  Scalar
!                               ATTRIBUTES: INTENT(IN OUT)
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

  FUNCTION Assign_Common_Output_AD( &
    Atmosphere             , & ! Input
    Surface                , & ! Input
    AtmOptics              , & ! Input
    SfcOptics              , & ! Input
    Pff_AD                 , & ! Input
    Pbb_AD                 , & ! Input
    GeometryInfo           , & ! Input
    SensorIndex            , & ! Input
    ChannelIndex           , & ! Input
    nZ                     , & ! Input
    AtmOptics_AD           , & ! Output
    SfcOptics_AD           , & ! Output
    Planck_Surface_AD      , & ! Output
    Planck_Atmosphere_AD   , & ! Output
    User_Emissivity_AD     , & ! Output
    Atmosphere_AD          , & ! Output
    Surface_AD             , & ! Output
    RTSolution_AD          , & ! Output
    RTV                    ) & ! Input
  RESULT( Error_Status )
    ! Arguments
    TYPE(CRTM_Atmosphere_type)  , INTENT(IN)     :: Atmosphere
    TYPE(CRTM_Surface_type)     , INTENT(IN)     :: Surface
    TYPE(CRTM_AtmOptics_type)   , INTENT(IN)     :: AtmOptics
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN)     :: SfcOptics
    REAL(fp)                    , INTENT(IN OUT) :: Pff_AD(:,:,:)
    REAL(fp)                    , INTENT(IN OUT) :: Pbb_AD(:,:,:)
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                     , INTENT(IN)     :: SensorIndex
    INTEGER                     , INTENT(IN)     :: ChannelIndex
    INTEGER                     , INTENT(IN)     :: nZ
    TYPE(CRTM_AtmOptics_type)   , INTENT(IN OUT) :: AtmOptics_AD
    TYPE(CRTM_SfcOptics_type)   , INTENT(IN OUT) :: SfcOptics_AD
    REAL(fp)                    , INTENT(IN OUT) :: Planck_Surface_AD
    REAL(fp)                    , INTENT(IN OUT) :: Planck_Atmosphere_AD(0:)
    REAL(fp)                    , INTENT(OUT)    :: User_Emissivity_AD
    TYPE(CRTM_Atmosphere_type)  , INTENT(IN OUT) :: Atmosphere_AD
    TYPE(CRTM_Surface_type)     , INTENT(IN OUT) :: Surface_AD
    TYPE(CRTM_RTSolution_type)  , INTENT(IN OUT) :: RTSolution_AD
    TYPE(RTV_type)              , INTENT(IN)     :: RTV

    ! Function Result
    INTEGER :: Error_Status
    ! Local Parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Assign_Common_Output_AD'
    ! Local variables
    CHARACTER(256) :: Message
    INTEGER :: no, na, nt
    INTEGER :: k, i

    Error_Status = SUCCESS

    ! ------------------------------------
    ! Compute the adjoint planck radiances
    ! ------------------------------------

    IF( RTV%mth_Azi == 0 ) THEN
      ! Surface AD radiance
      CALL CRTM_Planck_Radiance_AD( &
             SensorIndex                    , & ! Input
             ChannelIndex                   , & ! Input
             SfcOptics%Surface_Temperature  , & ! Input
             Planck_Surface_AD              , & ! Input
             SfcOptics_AD%Surface_Temperature ) ! In/Output
      Planck_Surface_AD = ZERO
      ! Atmospheric layer AD radiances
      DO k = 1, Atmosphere%n_Layers
        CALL CRTM_Planck_Radiance_AD( &
               SensorIndex                , & ! Input
               ChannelIndex               , & ! Input
               Atmosphere%Temperature(k)  , & ! Input
               Planck_Atmosphere_AD(k)    , & ! Input
               Atmosphere_AD%Temperature(k) ) ! In/Output
        Planck_Atmosphere_AD(k) = ZERO
      END DO
    ELSE
      Planck_Surface_AD = ZERO
      Planck_Atmosphere_AD = ZERO
    END IF

    ! ---------------------------------------------------------------------
    ! Assign the number of Stokes parameters. Currently, this is set to 1
    ! for decoupled polarization between surface and atmosphere.
    ! Remember for polarised microwave instruments, each frequency's Stokes
    ! parameter is treated as a separate channel.
    ! ---------------------------------------------------------------------
    SfcOptics_AD%n_Stokes = RTV%n_Stokes


    ! ----------------------------------------------------
    ! Assign the angles and weights to SfcOptics structure
    ! ----------------------------------------------------
    SfcOptics_AD%n_Angles = SfcOptics%n_Angles
    SfcOptics_AD%Angle    = SfcOptics%Angle
    SfcOptics_AD%Weight   = SfcOptics%Weight


    ! --------------------------------------------------------------------------------
    ! This part is designed for specific requirement for the total sensitivity
    ! of emissivity. The requirement is regardless whether having user input or not.
    ! Therefore, this part is common part with/without optional input.
    ! The outputs of CRTM_Compute_RTSolution are the partial derivative
    ! of emissivity and reflectivity. Since SfcOptics_AD is used as input only in
    ! CRTM_Compute_SfcOptics_AD, the total sensitivity of the emissivity is taken into
    ! account for here.
    ! --------------------------------------------------------------------------------

    ! -----------------------------------------
    ! Populate the SfcOptics_AD structure based
    ! on FORWARD model SfcOptics Compute_Switch
    ! -----------------------------------------
    IF ( SfcOptics%Compute ) THEN
      Error_Status = CRTM_Compute_SfcOptics_AD( &
                       Surface     , & ! Input
                       SfcOptics   , & ! Input
                       SfcOptics_AD, & ! Input
                       GeometryInfo, & ! Input
                       SensorIndex , & ! Input
                       ChannelIndex, & ! Input
                       Surface_AD  , & ! In/Output
                       RTV%SOV       ) ! Internal variable input

      IF ( Error_Status /= SUCCESS ) THEN
        WRITE( Message,'("Error computing SfcOptics_AD for ",a," channel ",i0)' ) &
               TRIM(SC(SensorIndex)%Sensor_Id), &
               SC(SensorIndex)%Sensor_Channel(ChannelIndex)
        CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status )
        RETURN
      END IF

    ELSE

      IF( RTV%Scattering_RT ) THEN
        User_Emissivity_AD = SUM(SfcOptics_AD%Emissivity(1:nZ,1))
        IF( RTV%Diffuse_Surface) THEN
          IF( RTV%mth_Azi == 0 ) THEN
        ! Assuming Lambertian surface isn't polarized!
          DO i = SfcOptics%n_Angles, 1, -1
            User_Emissivity_AD = User_Emissivity_AD - &
            (SUM(SfcOptics_AD%Reflectivity(1:SfcOptics%n_Angles,1,i,1))*SfcOptics%Weight(i))
          END DO
          ELSE
            SfcOptics_AD%Reflectivity(1:SfcOptics%n_Angles, 1, 1:SfcOptics%n_Angles, 1) = ZERO
            SfcOptics_AD%Direct_Reflectivity(1:SfcOptics%n_Angles,1) = ZERO
            SfcOptics_AD%Emissivity(1:SfcOptics%n_Angles,1) = ZERO
          END IF

        ELSE ! Specular surface
        DO i = nZ, 1, -1
          User_Emissivity_AD = User_Emissivity_AD - SfcOptics_AD%Reflectivity(i,1,i,1)
        END DO
        END IF
!      Direct_Reflectivity_AD = SUM(SfcOptics_AD%Direct_Reflectivity(1:nZ,1))
!      SfcOptics_AD%Direct_Reflectivity(1,1) = SfcOptics_AD%Direct_Reflectivity(1,1) +
!                                              (Direct_Reflectivity_AD/PI)
        RTSolution_AD%Surface_Emissivity = User_Emissivity_AD
      ELSE
!! need to check !!
        RTSolution_AD%Surface_Emissivity = SfcOptics_AD%Emissivity(SfcOptics_AD%Index_Sat_Ang,1) - &
                                         SfcOptics_AD%Reflectivity(1,1,1,1)
      END IF

    END IF

    ! --------------------------------------------
    ! Compute the adjoint phase matrix if required
    ! --------------------------------------------
    IF ( RTV%Scattering_RT ) THEN
      CALL CRTM_Phase_Matrix_AD( &
             AtmOptics,    &  ! Input - FWD
             Pff_AD,       &  ! Input - AD forward scattering
             Pbb_AD,       &  ! Input - AD backward scattering
             AtmOptics_AD, &  ! Output - AD
             RTV           )  ! Internal variable
    END IF

    ! ------------------------------------------
    ! Save the adjoint optical depth if required
    ! ------------------------------------------
    IF ( CRTM_RTSolution_Associated( RTSolution_AD ) ) THEN
      ! Shorter names for indexing
      no = RTSolution_AD%n_Layers  ! Original no. of layers
      na = RTV%n_Added_Layers      ! No. of added layers
      nt = RTV%n_Layers            ! Current total no. of layers
      ! Assign only the optical depth profile
      ! defined by the user input layering
      RTSolution_AD%Layer_Optical_Depth(1:no) = AtmOptics_AD%Optical_Depth(na+1:nt)
    END IF

    ! ----------------------------------------------------------------
    ! Compute the corresponding adjoint solar reflectance (TOA only)
    ! ----------------------------------------------------------------
    ! RTSolution_AD%Reflectance = ZERO
    ! IF( RTV%Solar_Flag_true .AND. (.NOT.(RTV%aircraft%rt))) THEN
    !  RTSolution_AD%Reflectance = RTSolution_AD%Radiance*PI/RTV%Solar_irradiance
    ! END IF

  END FUNCTION Assign_Common_Output_AD

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PRIVATE MODULE ROUTINES ##                      ##
!##                                                                            ##
!################################################################################
!################################################################################

! -------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Phase_Matrix
!
! PURPOSE:
!       Subroutine to calculate the phase function for the scattering model.
!
! CALLING SEQUENCE:
!       CALL CRTM_Phase_Matrix( AtmOptics,  &  ! Input
!                               RTV         )  ! Internal variable
!
! INPUT ARGUMENTS:
!       AtmOptics:      Structure containing the atmospheric optical
!                       parameters
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       RTV:            Structure containing internal forward model variables
!                       required for tangent-linear or adjoint model calls in
!                       CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUT ARGUMENTS:
!       RTV:            Structure containing internal forward model variables
!                       required for tangent-linear or adjoint model
!                       calls in CRTM_RTSolution module.
!                       UNITS:      N/A
!                       TYPE:       RTV_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!       Note that the INTENT on the RTV         argument is IN OUT as it
!       contains data prior to this call and is filled with data within this
!       routine.
!
!S-
!--------------------------------------------------------------------------------
  SUBROUTINE CRTM_Phase_Matrix( &
    AtmOptics, &  ! Input
    RTV        )  ! Internal variable
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: AtmOptics
    TYPE(RTV_type),            INTENT(IN OUT) :: RTV
    ! Local variables
    INTEGER :: i, j, k, l, ifac, jn, i1, j1

    !#--------------------------------------------------------------------------#
    !#                   -- COMPUTING LEGENDRE FUNCTION --                      #
    !#                   Pplus  for positive cosine of angles                   #
    !#                   Pminus for negative cosine of angles                   #
    !#--------------------------------------------------------------------------#

    IF( RTV%n_Stokes == 1 ) THEN

        !#--------------------------------------------------------------------------#
    !#                   -- COMPUTING LEGENDRE FUNCTION --                      #
    !#                   Pplus  for positive cosine of angles                   #
    !#                   Pminus for negative cosine of angles                   #
    !#--------------------------------------------------------------------------#
    DO i = 1, RTV%n_Angles
       CALL Legendre_M( RTV%mth_Azi               , &
                        AtmOptics%n_Legendre_Terms, &
                        RTV%COS_Angle(i)          , &
                        RTV%Pleg(0:,i)              )
    END DO

    IF(RTV%Solar_Flag_true) &
      CALL Legendre_M( RTV%mth_Azi                , &
                       AtmOptics%n_Legendre_Terms , &
                       RTV%COS_SUN                , &
                       RTV%Pleg(0:,RTV%n_Angles+1)  )

    IF( RTV%Solar_Flag_true ) THEN
      jn = RTV%n_Angles + 1
    ELSE
      jn = RTV%n_Angles
    END IF

    !#--------------------------------------------------------------------------#
    !#                    -- COMPUTE THE PHASE MATRICES --                      #
    !#--------------------------------------------------------------------------#

    Layer_Loop1: DO  k = 1, RTV%n_Layers


      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering1: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN

        DO j = 1, jn
          ! add solar angle
          DO i = 1, RTV%n_Angles

            RTV%Off(i,j,k)=ZERO
            RTV%Obb(i,j,k)=ZERO

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              RTV%Off(i,j,k) = RTV%Off(i,j,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              RTV%Obb(i,j,k) = RTV%Obb(i,j,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac )
            END DO

            RTV%Pff(i,j,k) = RTV%Off(i,j,k)
            RTV%Pbb(i,j,k) = RTV%Obb(i,j,k)

            ! For intensity, the phase matrix element must >= ZERO
            IF ( RTV%mth_Azi == 0 ) THEN
              IF(RTV%Pff(i,j,k) < ZERO) RTV%Pff(i,j,k) = PHASE_THRESHOLD
              IF(RTV%Pbb(i,j,k) < ZERO) RTV%Pbb(i,j,k) = PHASE_THRESHOLD
            END IF

          END DO
        END DO

        ! Normalization to ensure energy conservation
        IF (  RTV%mth_Azi == 0 ) CALL Normalize_Phase( k, RTV )
      END IF Significant_Scattering1

    END DO Layer_Loop1

      RETURN
    END IF

    L = AtmOptics%n_Legendre_Terms-1
    DO i = 1, RTV%n_Angles
       CALL Legendre_M( RTV%mth_Azi               , &
                        L, &
                        RTV%COS_Angle(i)          , &
                        RTV%Pleg(0:,i)              )
       IF(RTV%n_Stokes > 1 ) THEN
       !  Rlm
         RTV%Pplus(0:L,i)=-(Gl2n(RTV%mth_Azi,L,-2,RTV%COS_Angle(i))+Gl2n(RTV%mth_Azi,L,2,RTV%COS_Angle(i)))/TWO
       ! Tlm
         RTV%Pminus(0:L,i)=(Gl2n(RTV%mth_Azi,L,-2,RTV%COS_Angle(i))-Gl2n(RTV%mth_Azi,L,2,RTV%COS_Angle(i)))/TWO
       END IF
    END DO

    IF(RTV%Solar_Flag_true) THEN
      CALL Legendre_M( RTV%mth_Azi                , &
                       L, &
                       RTV%COS_SUN                , &
                       RTV%Pleg(0:,RTV%n_Angles+1)  )

       IF(RTV%n_Stokes > 1 ) THEN
       !  Rlm
         RTV%Pplus(0:L,RTV%n_Angles+1)=-(Gl2n(RTV%mth_Azi,L,-2,RTV%COS_SUN)+Gl2n(RTV%mth_Azi,L,2,RTV%COS_SUN))/TWO
       ! Tlm
         RTV%Pminus(0:L,RTV%n_Angles+1)=(Gl2n(RTV%mth_Azi,L,-2,RTV%COS_SUN)-Gl2n(RTV%mth_Azi,L,2,RTV%COS_SUN))/TWO
       END IF

    END IF

    IF( RTV%Solar_Flag_true ) THEN
      jn = RTV%n_Angles + 1
    ELSE
      jn = RTV%n_Angles
    END IF

    !#--------------------------------------------------------------------------#
    !#                    -- COMPUTE THE PHASE MATRICES --                      #
    !#--------------------------------------------------------------------------#

    Layer_Loop: DO  k = 1, RTV%n_Layers

      RTV%Pff(:,:,k) = ZERO
      RTV%Pbb(:,:,k) = ZERO
      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN

        DO j = 1, RTV%n_Angles
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1
            RTV%Off(i,j,k)=ZERO
            RTV%Obb(i,j,k)=ZERO

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              RTV%Off(i,j,k) = RTV%Off(i,j,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              RTV%Obb(i,j,k) = RTV%Obb(i,j,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac )

            IF( RTV%n_Stokes > 1 ) THEN
              !  (1,2)
              RTV%Pff(i1,j1+1,k) = RTV%Pff(i1,j1+1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j) )
              RTV%Pbb(i1,j1+1,k) = RTV%Pbb(i1,j1+1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j)*ifac )
              ! beta1  (2,1)
              RTV%Pff(i1+1,j1,k) = RTV%Pff(i1+1,j1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i) )
              RTV%Pbb(i1+1,j1,k) = RTV%Pbb(i1+1,j1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i)*ifac )
              ! alpha2, alpha3  (2,2)
              RTV%Pff(i1+1,j1+1,k) = RTV%Pff(i1+1,j1+1,k) + ( AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j) ) &
                                   + ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j) )
              RTV%Pbb(i1+1,j1+1,k) = RTV%Pbb(i1+1,j1+1,k) +AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac &
                                   - ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac )
            END IF
            IF( RTV%n_Stokes > 2 ) THEN
              ! alpha3, alpha2  (3,3)
              RTV%Pff(i1+2,j1+2,k) = RTV%Pff(i1+2,j1+2,k) + ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j) ) &
                                   + ( AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j) )
              RTV%Pbb(i1+2,j1+2,k) = RTV%Pbb(i1+2,j1+2,k)+AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac &
                                   - ( AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac )
              ! beta1  (1,3)
              RTV%Pff(i1,j1+2,k) = RTV%Pff(i1,j1+2,k) - ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j) )
              RTV%Pbb(i1,j1+2,k) = RTV%Pbb(i1,j1+2,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j)*ifac )

              ! beta1  (3,1)
              RTV%Pff(i1+2,j1,k) = RTV%Pff(i1+2,j1,k) - ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i) )
              RTV%Pbb(i1+2,j1,k) = RTV%Pbb(i1+2,j1,k) - ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i)*ifac )
              ! alpha2, alpha3  (2,3)
              RTV%Pff(i1+1,j1+2,k) = RTV%Pff(i1+1,j1+2,k) -AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j) &
                                 - ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j) )
              RTV%Pbb(i1+1,j1+2,k) = RTV%Pbb(i1+1,j1+2,k) + AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j)*ifac &
                                 - ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j)*ifac )

              ! alpha2, alpha3  (3,2)
              RTV%Pff(i1+2,j1+1,k) = RTV%Pff(i1+2,j1+1,k) - ( AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i) ) &
                                 - ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i) )
              RTV%Pbb(i1+2,j1+1,k) = RTV%Pbb(i1+2,j1+1,k) -AtmOptics%Phase_Coefficient(l,2,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i)*ifac &
                                   + ( AtmOptics%Phase_Coefficient(l,3,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i)*ifac )
            END IF
            IF( RTV%n_Stokes == 4 ) THEN
              ! beta2  (3,4)
              RTV%Pff(i1+2,j1+3,k) = RTV%Pff(i1+2,j1+3,k) - ( AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j) )
              RTV%Pbb(i1+2,j1+3,k) = RTV%Pbb(i1+2,j1+3,k) -AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j)*ifac

              ! beta2  (4,3)
              RTV%Pff(i1+3,j1+2,k) = RTV%Pff(i1+3,j1+2,k) + ( AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i) )
              RTV%Pbb(i1+3,j1+2,k) = RTV%Pbb(i1+3,j1+2,k) +AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i)*ifac
              ! alpha4 (4,4)
              RTV%Pff(i1+3,j1+3,k) = RTV%Pff(i1+3,j1+3,k) + ( AtmOptics%Phase_Coefficient(l,4,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              RTV%Pbb(i1+3,j1+3,k) = RTV%Pbb(i1+3,j1+3,k) +AtmOptics%Phase_Coefficient(l,4,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac

              ! beta2  (2,4)
              RTV%Pff(i1+1,j1+3,k) = RTV%Pff(i1+1,j1+3,k) + ( AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j) )
              RTV%Pbb(i1+1,j1+3,k) = RTV%Pbb(i1+1,j1+3,k) +AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j)*ifac
             ! beta2  (4,2)
             RTV%Pff(i1+3,j1+1,k) = RTV%Pff(i1+3,j1+1,k) - ( AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i) )
             RTV%Pbb(i1+3,j1+1,k) = RTV%Pbb(i1+3,j1+1,k) +AtmOptics%Phase_Coefficient(l,6,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i)*ifac
           END IF

            END DO

            RTV%Pff(i1,j1,k) = RTV%Off(i,j,k)
            RTV%Pbb(i1,j1,k) = RTV%Obb(i,j,k)

            ! For intensity, the phase matrix element must >= ZERO
            IF ( RTV%mth_Azi == 0 ) THEN
              IF(RTV%Pff(i1,j1,k) < ZERO) RTV%Pff(i1,j1,k) = PHASE_THRESHOLD
              IF(RTV%Pbb(i1,j1,k) < ZERO) RTV%Pbb(i1,j1,k) = PHASE_THRESHOLD
            END IF

            ! qliu   set P' = P D
            IF( RTV%n_Stokes > 2 ) THEN
              RTV%Pbb(i1,j1+2,k) = -RTV%Pbb(i1,j1+2,k)
              RTV%Pbb(i1+1,j1+2,k) = -RTV%Pbb(i1+1,j1+2,k)
              RTV%Pbb(i1+2,j1+2,k) = -RTV%Pbb(i1+2,j1+2,k)
            END IF
            IF( RTV%n_Stokes == 4 ) THEN
              RTV%Pbb(i1,j1+3,k) = -RTV%Pbb(i1,j1+3,k)
              RTV%Pbb(i1+1,j1+3,k) = -RTV%Pbb(i1+1,j1+3,k)
              RTV%Pbb(i1+2,j1+3,k) = -RTV%Pbb(i1+2,j1+3,k)

              RTV%Pbb(i1+3,j1+2,k) = -RTV%Pbb(i1+3,j1+2,k)
              RTV%Pbb(i1+3,j1+3,k) = -RTV%Pbb(i1+3,j1+3,k)
            END IF

          END DO
       END DO


       IF( jn > RTV%n_Angles ) THEN
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = RTV%n_Angles*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1
            RTV%Off(i,jn,k)=ZERO
            RTV%Obb(i,jn,k)=ZERO

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              RTV%Off(i,jn,k) = RTV%Off(i,jn,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn) )
              RTV%Obb(i,jn,k) = RTV%Obb(i,jn,k) + ( AtmOptics%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn)*ifac )
            IF( RTV%n_Stokes > 1 ) THEN
              ! beta1  (2,1)
              RTV%Pff(i1+1,j1,k) = RTV%Pff(i1+1,j1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i) )
              RTV%Pbb(i1+1,j1,k) = RTV%Pbb(i1+1,j1,k) + ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i)*ifac )
            END IF
            IF( RTV%n_Stokes > 2 ) THEN
              ! beta1  (3,1)
              RTV%Pff(i1+2,j1,k) = RTV%Pff(i1+2,j1,k) - ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i) )
              RTV%Pbb(i1+2,j1,k) = RTV%Pbb(i1+2,j1,k) - ( AtmOptics%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i)*ifac)
            END IF
            END DO
            RTV%Pff(i1,j1,k) = RTV%Off(i,jn,k)
            RTV%Pbb(i1,j1,k) = RTV%Obb(i,jn,k)

            ! For intensity, the phase matrix element must >= ZERO
            IF ( RTV%mth_Azi == 0 ) THEN
              IF(RTV%Pff(i1,j1,k) < ZERO) RTV%Pff(i1,j1,k) = PHASE_THRESHOLD
              IF(RTV%Pbb(i1,j1,k) < ZERO) RTV%Pbb(i1,j1,k) = PHASE_THRESHOLD
            END IF

          END DO

       END IF

       ! Normalization to ensure energy conservation
       IF (  RTV%mth_Azi == 0 ) CALL Normalize_Phase( k, RTV )  ! normalize RTV%Pff and RTV%Pbb
                                  ! (n_Stokes*n_Angles, n_Stokes*n_Angles)
      END IF Significant_Scattering

    END DO Layer_Loop

  END SUBROUTINE CRTM_Phase_Matrix

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Phase_Matrix_TL
!
! PURPOSE:
!       Subroutine to calculate the tangent-linear phase function for the
!       scattering model.
!
! CALLING SEQUENCE:
!       CALL CRTM_Phase_Matrix_TL( AtmOptics,    &  ! FWD Input
!                                  AtmOptics_TL, &  ! TL  Input
!                                  Pff_TL,       &  ! TL Output
!                                  Pff_TL,       &  ! TL Output
!                                  RTV           )  ! Internal variable
!
! INPUT ARGUMENTS:
!       AtmOptics:      Structure containing the atmospheric optical
!                       parameters
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics_TL:   Structure containing the tangent-linear atmospheric
!                       optical parameters
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
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
!       Pff_TL:         Array containing the tangent-linear of the
!                       forward phase matrix.
!                       UNITS:      N/A
!                       TYPE:       REAL
!                       DIMENSION:  Rank-3, n_Angles x n_Angles x n_Layers
!                       ATTRIBUTES: INTENT(OUT)
!
!       Pbb_TL:         Array containing the tangent-linear of the
!                       backward phase matrix.
!                       UNITS:      N/A
!                       TYPE:       REAL
!                       DIMENSION:  Rank-3, n_Angles x n_Angles x n_Layers
!                       ATTRIBUTES: INTENT(OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Phase_Matrix_TL( &
    AtmOptics,    &  ! FWD Input
    AtmOptics_TL, &  ! TL  Input
    Pff_TL,       &  ! TL Output
    Pbb_TL,       &  ! TL Output
    RTV           )  ! Internal variable
    ! Arguments
    TYPE(CRTM_AtmOptics_type)  , INTENT(IN)  :: AtmOptics
    TYPE(CRTM_AtmOptics_type)  , INTENT(IN)  :: AtmOptics_TL
    REAL(fp)                   , INTENT(OUT) :: Pff_TL(:,:,:) ! n_Angles x n_Angles x n_Layers
    REAL(fp)                   , INTENT(OUT) :: Pbb_TL(:,:,:) ! n_Angles x n_Angles x n_Layers
    TYPE(RTV_type)             , INTENT(IN)  :: RTV
    ! Local variables
    INTEGER :: i, j, k, l, nZ, ifac, jn, i1, j1
    REAL(fp), DIMENSION(RTV%n_Angles,RTV%n_Angles+1) :: Lff, Lbb, Lff_TL, Lbb_TL



    !#--------------------------------------------------------------------------#
    !#             -- COMPUTE THE TANGENT-LINEAR PHASE MATRICES --              #
    !#--------------------------------------------------------------------------#

    Pff_TL = ZERO
    Pbb_TL = ZERO

    nZ = RTV%n_Angles
    IF( RTV%Solar_Flag_true ) THEN
      jn = RTV%n_Angles + 1
    ELSE
      jn = RTV%n_Angles
    END IF

    IF(RTV%n_Stokes == 1 ) THEN

    Layer_Loop1: DO  k = 1, RTV%n_Layers

      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering1: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN

        Lff(1:nZ,1:jn) = RTV%Off(1:nZ,1:jn,k)
        Lbb(1:nZ,1:jn) = RTV%Obb(1:nZ,1:jn,k)


        DO j = 1, jn
          DO i = 1, RTV%n_Angles

            Pff_TL(i,j,k) = ZERO
            Pbb_TL(i,j,k) = ZERO

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              Pff_TL(i,j,k) = Pff_TL(i,j,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              Pbb_TL(i,j,k) = Pbb_TL(i,j,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac )
            END DO

            ! For intensity, the FWD phase matrix element must >= ZERO
            ! so the TL form is always == zero.
            IF ( RTV%mth_Azi == 0 ) THEN
              IF ( RTV%Off(i,j,k) < ZERO ) THEN
                Pff_TL(i,j,k) = ZERO
                Lff(i,j)      = PHASE_THRESHOLD
              END IF

              IF ( RTV%Obb(i,j,k) < ZERO ) THEN
                Pbb_TL(i,j,k) = ZERO
                Lbb(i,j) = PHASE_THRESHOLD
              END IF
            END IF

          END DO
        END DO

        IF ( RTV%mth_Azi == 0 ) THEN
          ! Normalisation for energy conservation
          CALL Normalize_Phase_TL( &
                 k, RTV, &
                 Lff,           & ! FWD Input
                 Lbb,           & ! FWD Input
                 Pff_TL(:,:,k), & ! TL  Output
                 Pbb_TL(:,:,k)  ) ! TL  Output
        END IF
      END IF Significant_Scattering1

    END DO Layer_Loop1
    RETURN

    END IF





    Layer_Loop: DO  k = 1, RTV%n_Layers

      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN

        DO j = 1, RTV%n_Angles
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1

          n_Legendre_Terms_loop:  DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)

              Pff_TL(i1,j1,k) = Pff_TL(i1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              Pbb_TL(i1,j1,k) = Pbb_TL(i1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac )

            IF( RTV%n_Stokes > 1) THEN
              !  (1,2)
              Pff_TL(i1,j1+1,k) = Pff_TL(i1,j1+1,k) + ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j) )
              Pbb_TL(i1,j1+1,k) = Pbb_TL(i1,j1+1,k) + ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j)*ifac )

              ! beta1  (2,1)
              Pff_TL(i1+1,j1,k) = Pff_TL(i1+1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i) )
              Pbb_TL(i1+1,j1,k) = Pbb_TL(i1+1,j1,k) +AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i)*ifac

              ! alpha2, alpha3  (2,2)
              Pff_TL(i1+1,j1+1,k) = Pff_TL(i1+1,j1+1,k) +AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j) &
                                   + ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j) )
              Pbb_TL(i1+1,j1+1,k) = Pbb_TL(i1+1,j1+1,k)+AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac &
                                   - ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac )
             END IF

             IF( RTV%n_Stokes > 2 ) THEN
              ! alpha3, alpha2  (3,3)
              Pff_TL(i1+2,j1+2,k) = Pff_TL(i1+2,j1+2,k) +AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j) &
                                   + ( AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j) )
              Pbb_TL(i1+2,j1+2,k) = Pbb_TL(i1+2,j1+2,k)+AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac &
                                   - ( AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac )
              ! beta1  (1,3)
              Pff_TL(i1,j1+2,k) = Pff_TL(i1,j1+2,k) -AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j)
              Pbb_TL(i1,j1+2,k) = Pbb_TL(i1,j1+2,k) +AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j)*ifac

              ! beta1  (3,1)
              Pff_TL(i1+2,j1,k) = Pff_TL(i1+2,j1,k) - AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i)
              Pbb_TL(i1+2,j1,k) = Pbb_TL(i1+2,j1,k) - AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i)*ifac

              ! alpha2, alpha3  (2,3)
              Pff_TL(i1+1,j1+2,k) = Pff_TL(i1+1,j1+2,k) -AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j) &
                                 - ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j) )
              Pbb_TL(i1+1,j1+2,k)=Pbb_TL(i1+1,j1+2,k)+AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j)*ifac &
                                 - ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j)*ifac )

              ! alpha2, alpha3  (3,2)
              Pff_TL(i1+2,j1+1,k) = Pff_TL(i1+2,j1+1,k) -AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i) &
                                 - ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i) )
              Pbb_TL(i1+2,j1+1,k) = Pbb_TL(i1+2,j1+1,k)-AtmOptics_TL%Phase_Coefficient(l,2,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i)*ifac &
                                   + ( AtmOptics_TL%Phase_Coefficient(l,3,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i)*ifac )
             END IF
             IF( RTV%n_Stokes == 4 ) THEN
              ! beta2  (3,4)
              Pff_TL(i1+2,j1+3,k) = Pff_TL(i1+2,j1+3,k) -AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j)
              Pbb_TL(i1+2,j1+3,k) = Pbb_TL(i1+2,j1+3,k) -AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j)*ifac

              ! beta2  (4,3)
              Pff_TL(i1+3,j1+2,k) = Pff_TL(i1+3,j1+2,k) +AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i)
              Pbb_TL(i1+3,j1+2,k) = Pbb_TL(i1+3,j1+2,k) +AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i)*ifac
              ! alpha4 (4,4)
              Pff_TL(i1+3,j1+3,k) = Pff_TL(i1+3,j1+3,k) +AtmOptics_TL%Phase_Coefficient(l,4,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)
              Pbb_TL(i1+3,j1+3,k) = Pbb_TL(i1+3,j1+3,k) +AtmOptics_TL%Phase_Coefficient(l,4,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac

              ! beta2  (2,4)
              Pff_TL(i1+1,j1+3,k) = Pff_TL(i1+1,j1+3,k) +AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j)
              Pbb_TL(i1+1,j1+3,k) = Pbb_TL(i1+1,j1+3,k) +AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j)*ifac

              ! beta2  (4,2)
              Pff_TL(i1+3,j1+1,k) = Pff_TL(i1+3,j1+1,k) -AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i)
              Pbb_TL(i1+3,j1+1,k) = Pbb_TL(i1+3,j1+1,k) +AtmOptics_TL%Phase_Coefficient(l,6,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i)*ifac
             END IF

            END DO n_Legendre_Terms_loop

        ! qliu   set P' = P D
           IF( RTV%n_Stokes > 2 ) THEN
             Pbb_TL(i1,j1+2,k)   = -Pbb_TL(i1,j1+2,k)
             Pbb_TL(i1+1,j1+2,k) = -Pbb_TL(i1+1,j1+2,k)
             Pbb_TL(i1+2,j1+2,k) = -Pbb_TL(i1+2,j1+2,k)
           END IF
           IF( RTV%n_Stokes == 4 ) THEN
             Pbb_TL(i1,j1+3,k)   = -Pbb_TL(i1,j1+3,k)
             Pbb_TL(i1+1,j1+3,k) = -Pbb_TL(i1+1,j1+3,k)
             Pbb_TL(i1+2,j1+3,k) = -Pbb_TL(i1+2,j1+3,k)
             Pbb_TL(i1+3,j1+2,k) = -Pbb_TL(i1+3,j1+2,k)
             Pbb_TL(i1+3,j1+3,k) = -Pbb_TL(i1+3,j1+3,k)
           END IF

           IF( RTV%mth_Azi == 0 ) THEN
              Lff_TL(i,j) = Pff_TL(i1,j1,k)
              Lbb_TL(i,j) = Pbb_TL(i1,j1,k)
           END IF

          END DO
        END DO

       IF( jn > RTV%n_Angles ) THEN
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = RTV%n_Angles*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              Pff_TL(i1,j1,k) = Pff_TL(i1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn) )
              Pbb_TL(i1,j1,k) = Pbb_TL(i1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn)*ifac )

            IF( RTV%n_Stokes > 1 ) THEN
              ! beta1  (2,1)
              Pff_TL(i1+1,j1,k)  = Pff_TL(i1+1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i) )
              Pbb_TL(i1+1,j1,k)  = Pbb_TL(i1+1,j1,k) + ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i)*ifac )
            END IF
            IF( RTV%n_Stokes > 2 ) THEN
              ! beta1  (3,1)
              Pff_TL(i1+2,j1,k) = Pff_TL(i1+2,j1,k) - ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i) )
              Pbb_TL(i1+2,j1,k) = Pbb_TL(i1+2,j1,k) - ( AtmOptics_TL%Phase_Coefficient(l,5,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i)*ifac)
            END IF
            END DO
          END DO
       END IF

      IF ( RTV%mth_Azi == 0 ) THEN
        DO j = 1, jn
          DO i = 1, RTV%n_Angles
            Lff(i,j) = RTV%Off(i,j,k)
            Lbb(i,j) = RTV%Obb(i,j,k)
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1

            Lff_TL(i,j) = Pff_TL(i1,j1,k)
            Lbb_TL(i,j) = Pbb_TL(i1,j1,k)
            Lff(i,j) =  RTV%Off(i,j,k)
            Lbb(i,j) =  RTV%Obb(i,j,k)
            ! For intensity, the FWD phase matrix element must >= ZERO
            ! so the TL form is always == zero.
            IF ( RTV%Off(i,j,k) < ZERO ) THEN
              Lff_TL(i,j) = ZERO
              Lff(i,j)      = PHASE_THRESHOLD
            END IF

            IF ( RTV%Obb(i,j,k) < ZERO ) THEN
              Lbb_TL(i,j) = ZERO
              Lbb(i,j) = PHASE_THRESHOLD
            END IF
          END DO
        END DO

          ! Normalisation for energy conservation
          ! Using FWD results Lff, Lbb and normalize Pff_TL, Pbb_TL (n_Angles, n_Angles)
          CALL Normalize_Phase_TL( &
                 k, RTV, &
                 Lff,           & ! FWD Input
                 Lbb,           & ! FWD Input
                 Lff_TL(:,:), & ! TL  Output
                 Lbb_TL(:,:)  ) ! TL  Output

        DO j = 1, jn
          DO i = 1, RTV%n_Angles
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1
            Pff_TL(i1,j1,k) = Lff_TL(i,j)
            Pbb_TL(i1,j1,k) = Lbb_TL(i,j)
          END DO
        END DO
      END IF


      END IF Significant_Scattering

    END DO Layer_Loop

  END SUBROUTINE CRTM_Phase_Matrix_TL

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Phase_Matrix_AD
!
! PURPOSE:
!       Subroutine to calculate the adjoint of the phase function for the
!       scattering model.
!
! CALLING SEQUENCE:
!       CALL CRTM_Phase_Matrix_AD( AtmOptics,    &  ! FWD Input
!                                  Pff_AD,       &  ! AD Input
!                                  Pbb_AD,       &  ! AD Input
!                                  AtmOptics_AD, &  ! AD Output
!                                  RTV           )  ! Internal variable
!
! INPUT ARGUMENTS:
!       AtmOptics:      Structure containing the atmospheric optical
!                       parameters
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Pff_AD:         Array containing the adjoint of the
!                       forward phase matrix.
!                       ** NOTE: This argument will be zeroed upon exit
!                       **       from this routine.
!                       UNITS:      N/A
!                       TYPE:       REAL
!                       DIMENSION:  Rank-3, n_Angles x n_Angles x n_Layers
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       Pbb_AD:         Array containing the adjoint of the
!                       backward phase matrix.
!                       ** NOTE: This argument will be zeroed upon exit
!                       **       from this routine.
!                       UNITS:      N/A
!                       TYPE:       REAL
!                       DIMENSION:  Rank-3, n_Angles x n_Angles x n_Layers
!                       ATTRIBUTES: INTENT(IN OUT)
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
!       AtmOptics_AD:   Structure containing the adjoint atmospheric optical
!                       parameters
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
! COMMENTS:
!       Note the INTENT on the output AtmOptics argument is IN OUT rather than
!       just OUT. This is necessary because the argument MUST be defined upon
!       input. To prevent memory leaks, and in this case errors in accessing
!       unallocated memory, the IN OUT INTENT is a must.
!
!S-
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Phase_Matrix_AD( &
    AtmOptics,    &  ! FWD Input
    Pff_AD,       &  ! AD Input
    Pbb_AD,       &  ! AD Input
    AtmOptics_AD, &  ! AD Output
    RTV           )  ! Internal variable
    ! Arguments
    TYPE(CRTM_AtmOptics_type), INTENT(IN)     :: AtmOptics
    REAL(fp)                 , INTENT(IN OUT) :: Pff_AD(:,:,:) ! n_Angles x n_Angles x n_Layers
    REAL(fp)                 , INTENT(IN OUT) :: Pbb_AD(:,:,:) ! n_Angles x n_Angles x n_Layers
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_AD
    TYPE(RTV_type)           , INTENT(IN)     :: RTV
    ! Local variables
    INTEGER :: i, j, k, l, nZ, ifac, jn, i1, j1
    REAL(fp), DIMENSION(RTV%n_Angles,RTV%n_Angles+1) :: Lff, Lbb,  Lff_AD, Lbb_AD

    !#--------------------------------------------------------------------------#
    !#             -- COMPUTE THE ADJOINT OF THE PHASE MATRICES --              #
    !#--------------------------------------------------------------------------#

     Lff_AD = ZERO
     Lbb_AD = ZERO

    nZ = RTV%n_Angles
    IF( RTV%Solar_Flag_true ) THEN
      jn = RTV%n_Angles + 1
    ELSE
      jn = RTV%n_Angles
    END IF


    IF( RTV%n_Stokes == 1 ) THEN

    Layer_Loop1: DO  k = 1, RTV%n_Layers


      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering1: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN


        ! AD normalization to ensure energy conservation
     !!   AtmOptics_AD%Phase_Coefficient(0:,1,k) = ZERO

        Lff(1:nZ,1:jn) = RTV%Off(1:nZ,1:jn,k)
        Lbb(1:nZ,1:jn) = RTV%Obb(1:nZ,1:jn,k)

        IF ( RTV%mth_Azi == 0 ) THEN
          DO j = 1, jn
            DO i = 1, RTV%n_Angles

              ! For intensity, the FWD phase matrix element must >= ZERO
              ! so the TL, and thus the AD, for is always == zero.

                IF ( RTV%Off(i,j,k) < ZERO) Lff(i,j) = PHASE_THRESHOLD
                IF ( RTV%Obb(i,j,k) < ZERO) Lbb(i,j) = PHASE_THRESHOLD

            END DO
          END DO

          CALL Normalize_Phase_AD( &
                 k, RTV, &
                 Lff, Lbb,      & ! FWD Input
                 Pff_AD(:,:,k), & ! AD  Output
                 Pbb_AD(:,:,k)  ) ! AD  Output
        END IF

        DO j = 1, jn
          DO i = 1, RTV%n_Angles

            ! For intensity, the FWD phase matrix element must >= ZERO
            ! so the TL, and thus the AD, for is always == zero.
            IF ( RTV%mth_Azi == 0 ) THEN
              IF ( RTV%Off(i,j,k) < ZERO) Pff_AD(i,j,k) = ZERO
              IF ( RTV%Obb(i,j,k) < ZERO) Pbb_AD(i,j,k) = ZERO
            END IF

            DO l = RTV%mth_Azi, AtmOptics%n_Legendre_Terms - 1
              ifac = (-1) ** (l - RTV%mth_Azi)
              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                                                      ( Pff_AD(i,j,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j) )
              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                                                      ( Pbb_AD(i,j,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac )
            END DO

            Pff_AD(i,j,k) = ZERO
            Pbb_AD(i,j,k) = ZERO

          END DO
        END DO

      END IF Significant_Scattering1

    END DO Layer_Loop1

    RETURN
    END IF
!    AtmOptics_AD%Phase_Coefficient(0:,:,:) = ZERO

    Layer_Loop: DO  k = 1, RTV%n_Layers

      ! ------------------------------
      ! Only proceed if the scattering
      ! coefficient is significant
      ! ------------------------------

      Significant_Scattering: IF( AtmOptics%Single_Scatter_Albedo(k) > SCATTERING_ALBEDO_THRESHOLD) THEN

      IF ( RTV%mth_Azi == 0 ) THEN
        DO j = 1, jn
          DO i = 1, RTV%n_Angles
            Lff(i,j) = RTV%Off(i,j,k)
            Lbb(i,j) = RTV%Obb(i,j,k)
            IF ( RTV%Off(i,j,k) < ZERO) Lff(i,j) = PHASE_THRESHOLD
            IF ( RTV%Obb(i,j,k) < ZERO) Lbb(i,j) = PHASE_THRESHOLD
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1
            Lff_AD(i,j) = Pff_AD(i1,j1,k)
            Lbb_AD(i,j) = Pbb_AD(i1,j1,k)
          END DO
        END DO

        CALL Normalize_Phase_AD( &
              k, RTV, &
              Lff, Lbb,      & ! FWD Input
              Lff_AD, & ! AD  Output
              Lbb_AD  ) ! AD  Output

        DO j = 1, RTV%n_Angles
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1
            IF ( RTV%Off(i,j,k) < ZERO) Lff_AD(i,j) = ZERO
            IF ( RTV%Obb(i,j,k) < ZERO) Lbb_AD(i,j) = ZERO
            Pff_AD(i1,j1,k) = Lff_AD(i,j)
            Pbb_AD(i1,j1,k) = Lbb_AD(i,j)
          END DO
        END DO

    END IF

       IF( jn > RTV%n_Angles ) THEN
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = RTV%n_Angles*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1

            DO l = AtmOptics%n_Legendre_Terms - 1, RTV%mth_Azi, -1
              ifac = (-1) ** (l - RTV%mth_Azi)

            IF( RTV%n_Stokes > 2 ) THEN
            ! beta1  (3,1)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) - &
                 ( Pbb_AD(i1+2,j1,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i)*ifac)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) - &
                 ( Pff_AD(i1+2,j1,k)*RTV%Pleg(l,jn)*RTV%Pminus(l,i))
            END IF

            IF( RTV%n_Stokes > 1 ) THEN
              ! beta1  (2,1)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                 ( Pbb_AD(i1+1,j1,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i)*ifac)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                 ( Pff_AD(i1+1,j1,k)*RTV%Pleg(l,jn)*RTV%Pplus(l,i) )
            END IF

              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                ( Pbb_AD(i1,j1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn)*ifac )
              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                  ( Pff_AD(i1,j1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,jn) )
            END DO
          END DO
       END IF


        DO j = 1, RTV%n_Angles
          ! add solar angle
          DO i = 1, RTV%n_Angles
            j1 = (j-1)*RTV%n_Stokes + 1
            i1 = (i-1)*RTV%n_Stokes + 1

!           IF( RTV%mth_Azi == 0 ) THEN
!              Lff_TL(i,j) = Pff_TL(i1,j1,k)
!              Lbb_TL(i,j) = Pbb_TL(i1,j1,k)
!           END IF

        ! qliu   set P' = P D
           IF( RTV%n_Stokes == 4 ) THEN
             Pbb_AD(i1,j1+3,k)   = -Pbb_AD(i1,j1+3,k)
             Pbb_AD(i1+1,j1+3,k) = -Pbb_AD(i1+1,j1+3,k)
             Pbb_AD(i1+2,j1+3,k) = -Pbb_AD(i1+2,j1+3,k)
             Pbb_AD(i1+3,j1+2,k) = -Pbb_AD(i1+3,j1+2,k)
             Pbb_AD(i1+3,j1+3,k) = -Pbb_AD(i1+3,j1+3,k)
           END IF

           IF( RTV%n_Stokes > 2 ) THEN
             Pbb_AD(i1+2,j1+2,k) = -Pbb_AD(i1+2,j1+2,k)
             Pbb_AD(i1+1,j1+2,k) = -Pbb_AD(i1+1,j1+2,k)
             Pbb_AD(i1,j1+2,k)   = -Pbb_AD(i1,j1+2,k)
           END IF

         n_Legendre_Terms_loop: DO l = AtmOptics%n_Legendre_Terms - 1, RTV%mth_Azi, -1
              ifac = (-1) ** (l - RTV%mth_Azi)

            IF( RTV%n_Stokes == 4 ) THEN
             ! beta2  (4,2)
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) + &
               Pbb_AD(i1+3,j1+1,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i)*ifac

              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) - &
               Pff_AD(i1+3,j1+1,k)*RTV%Pminus(l,j)*RTV%Pleg(l,i)

              ! beta2  (2,4)
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) + &
              Pbb_AD(i1+1,j1+3,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) + &
              Pff_AD(i1+1,j1+3,k)*RTV%Pminus(l,i)*RTV%Pleg(l,j)

              ! alpha4 (4,4)
              AtmOptics_AD%Phase_Coefficient(l,4,k) = AtmOptics_AD%Phase_Coefficient(l,4,k) + &
               Pbb_AD(i1+3,j1+3,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,4,k) = AtmOptics_AD%Phase_Coefficient(l,4,k) + &
               Pff_AD(i1+3,j1+3,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)

              ! beta2  (4,3)
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) + &
               Pbb_AD(i1+3,j1+2,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i)*ifac
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) + &
               Pff_AD(i1+3,j1+2,k)*RTV%Pplus(l,j)*RTV%Pleg(l,i)

              ! beta2  (3,4)
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) - &
                Pbb_AD(i1+2,j1+3,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,6,k) = AtmOptics_AD%Phase_Coefficient(l,6,k) - &
                Pff_AD(i1+2,j1+3,k)*RTV%Pplus(l,i)*RTV%Pleg(l,j)
            END IF


            IF( RTV%n_Stokes > 2 ) THEN
              ! alpha2, alpha3  (3,2)
              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) + &
                Pbb_AD(i1+2,j1+1,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i)*ifac

              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) - &
                Pbb_AD(i1+2,j1+1,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i)*ifac

              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) - &
                Pff_AD(i1+2,j1+1,k)*RTV%Pminus(l,j)*RTV%Pplus(l,i)
              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) - &
                Pff_AD(i1+2,j1+1,k)*RTV%Pplus(l,j)*RTV%Pminus(l,i)

              ! alpha2, alpha3  (2,3)
              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) - &
                Pbb_AD(i1+1,j1+2,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j)*ifac

              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) + &
                Pbb_AD(i1+1,j1+2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j)*ifac

              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) - &
                Pff_AD(i1+1,j1+2,k)*RTV%Pminus(l,i)*RTV%Pplus(l,j)
              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) - &
                Pff_AD(i1+1,j1+2,k)*RTV%Pplus(l,i)*RTV%Pminus(l,j)

              ! beta1  (3,1)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) - &
                Pbb_AD(i1+2,j1,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i)*ifac
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) - &
                Pff_AD(i1+2,j1,k)*RTV%Pleg(l,j)*RTV%Pminus(l,i)

              ! beta1  (1,3)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                 Pbb_AD(i1,j1+2,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) - &
                 Pff_AD(i1,j1+2,k)*RTV%Pleg(l,i)*RTV%Pminus(l,j)

              ! alpha3, alpha2  (3,3)
              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) - &
                Pbb_AD(i1+2,j1+2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac

             AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) + &
                Pbb_AD(i1+2,j1+2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac

              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) + &
                Pff_AD(i1+2,j1+2,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)
              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) + &
                Pff_AD(i1+2,j1+2,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)
            END IF

            ! alpha2, alpha3  (2,2)
            IF( RTV%n_Stokes > 1 ) THEN
              AtmOptics_AD%Phase_Coefficient(l,3,k) = AtmOptics_AD%Phase_Coefficient(l,3,k) - &
                Pbb_AD(i1+1,j1+1,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) + &
                Pbb_AD(i1+1,j1+1,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)*ifac

              AtmOptics_AD%Phase_Coefficient(l,3,k) =  AtmOptics_AD%Phase_Coefficient(l,3,k) + &
                Pff_AD(i1+1,j1+1,k)*RTV%Pminus(l,i)*RTV%Pminus(l,j)
              AtmOptics_AD%Phase_Coefficient(l,2,k) = AtmOptics_AD%Phase_Coefficient(l,2,k) + &
                Pff_AD(i1+1,j1+1,k)*RTV%Pplus(l,i)*RTV%Pplus(l,j)

              ! beta1  (2,1)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                Pbb_AD(i1+1,j1,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i)*ifac
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                Pff_AD(i1+1,j1,k)*RTV%Pleg(l,j)*RTV%Pplus(l,i)

              !  (1,2)
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                Pbb_AD(i1,j1+1,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,5,k) = AtmOptics_AD%Phase_Coefficient(l,5,k) + &
                Pff_AD(i1,j1+1,k)*RTV%Pleg(l,i)*RTV%Pplus(l,j)
            END IF

              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                Pbb_AD(i1,j1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)*ifac
              AtmOptics_AD%Phase_Coefficient(l,1,k) = AtmOptics_AD%Phase_Coefficient(l,1,k) + &
                Pff_AD(i1,j1,k)*RTV%Pleg(l,i)*RTV%Pleg(l,j)

            END DO n_Legendre_Terms_loop

          END DO
        END DO

      END IF Significant_Scattering

    END DO Layer_Loop

  END SUBROUTINE CRTM_Phase_Matrix_AD

!**************************************************************
!                 Normalize Phase Subroutines
!**************************************************************

  SUBROUTINE Normalize_Phase( k, RTV )
    ! Arguments
    INTEGER,        INTENT(IN)     :: k
    TYPE(RTV_type), INTENT(IN OUT) :: RTV
    ! Local variables
    INTEGER :: i, j, nZ, i1, j1

    nZ = RTV%n_Angles

    IF( RTV%n_Stokes == 1 ) THEN
        ! Normalisation for stream angles
      RTV%Sum_Fac(0,k)=ZERO
      DO i = 1, RTV%n_Streams
        RTV%n_Factor(i,k)=ZERO
        DO j = i,nZ
          RTV%n_Factor(i,k)=RTV%n_Factor(i,k)+(RTV%Pff(i,j,k)+RTV%Pbb(i,j,k))*RTV%COS_Weight(j)
        END DO
        DO j=i,nZ
          RTV%Pff(i,j,k)=RTV%Pff(i,j,k)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))
          RTV%Pbb(i,j,k)=RTV%Pbb(i,j,k)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))
        END DO
        RTV%Sum_Fac(i,k)=ZERO
        IF( i < nZ ) THEN
          DO j=i+1,nZ
            RTV%Pff(j,i,k)=RTV%Pff(i,j,k)
            RTV%Pbb(j,i,k)=RTV%Pbb(i,j,k)
          END DO
          DO j = 1, i
            RTV%Sum_Fac(i,k)=RTV%Sum_Fac(i,k) + (RTV%Pff(i+1,j,k)+RTV%Pbb(i+1,j,k))*RTV%COS_Weight(j)
          END DO
        END IF
      END DO

      IF( RTV%n_Streams < nZ ) THEN
        ! Sensor viewing angle differs from the Gaussian angles
        RTV%n_Factor(nZ,k) =  RTV%Sum_Fac(nZ-1,k)
        DO j = 1, nZ
          RTV%Pff(j,nZ,k) = RTV%Pff(j,nZ,k)/RTV%n_Factor(nZ,k)
          RTV%Pbb(j,nZ,k) = RTV%Pbb(j,nZ,k)/RTV%n_Factor(nZ,k)
        ! Symmetric condition
          IF( j < nZ ) THEN
            RTV%Pff(nZ,j,k) = RTV%Pff(j,nZ,k)
            RTV%Pbb(nZ,j,k) = RTV%Pbb(j,nZ,k)
          END IF
        END DO
      END IF

      RETURN
    END IF


    ! Normalisation for stream angles
    RTV%Sum_Fac(0,k)=ZERO
    DO i = 1, RTV%n_Streams
      RTV%n_Factor(i,k)=ZERO
      DO j = i,nZ
        j1 = (j-1)*RTV%n_Stokes + 1
        i1 = (i-1)*RTV%n_Stokes + 1
        RTV%n_Factor(i,k)=RTV%n_Factor(i,k)+(RTV%Pff(i1,j1,k)+RTV%Pbb(i1,j1,k))*RTV%COS_Weight(j)
      END DO
      DO j=i,nZ
        j1 = (j-1)*RTV%n_Stokes + 1
        RTV%Pff(i1,j1,k)=RTV%Pff(i1,j1,k)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))
        RTV%Pbb(i1,j1,k)=RTV%Pbb(i1,j1,k)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))
      END DO
      RTV%Sum_Fac(i,k)=ZERO
      IF( i < nZ ) THEN
        DO j=i+1,nZ
          j1 = (j-1)*RTV%n_Stokes + 1
          RTV%Pff(j1,i1,k)=RTV%Pff(i1,j1,k)
          RTV%Pbb(j1,i1,k)=RTV%Pbb(i1,j1,k)
        END DO
        DO j = 1, i
          j1 = (j-1)*RTV%n_Stokes + 1
          RTV%Sum_Fac(i,k)=RTV%Sum_Fac(i,k) + (RTV%Pff(i1+RTV%n_Stokes,j1,k)+RTV%Pbb(i1+RTV%n_Stokes,j1,k))*RTV%COS_Weight(j)
        END DO
      END IF
    END DO

    IF( RTV%n_Streams < nZ ) THEN
      ! Sensor viewing angle differs from the Gaussian angles
      RTV%n_Factor(nZ,k) =  RTV%Sum_Fac(nZ-1,k)
      DO j = 1, nZ
        j1 = (j-1)*RTV%n_Stokes + 1
        RTV%Pff(j1,nZ*RTV%n_Stokes,k) = RTV%Pff(j1,nZ*RTV%n_Stokes,k)/RTV%n_Factor(nZ,k)
        RTV%Pbb(j1,nZ*RTV%n_Stokes,k) = RTV%Pbb(j1,nZ*RTV%n_Stokes,k)/RTV%n_Factor(nZ,k)
        ! Symmetric condition
        IF( j < nZ ) THEN
          RTV%Pff(nZ*RTV%n_Stokes,j1,k) = RTV%Pff(j1,nZ*RTV%n_Stokes,k)
          RTV%Pbb(nZ*RTV%n_Stokes,j1,k) = RTV%Pbb(j1,nZ*RTV%n_Stokes,k)
        END IF
      END DO
    END IF

  END SUBROUTINE Normalize_Phase

  SUBROUTINE Normalize_Phase_TL( k, RTV, Pff, Pbb, Pff_TL, Pbb_TL )
    ! Arguments
    INTEGER       , INTENT(IN)     :: k
    TYPE(RTV_type), INTENT(IN)     :: RTV
    REAL(fp)      , INTENT(IN)     :: Pff(:,:)
    REAL(fp)      , INTENT(IN)     :: Pbb(:,:)
    REAL(fp)      , INTENT(IN OUT) :: Pff_TL(:,:)
    REAL(fp)      , INTENT(IN OUT) :: Pbb_TL(:,:)
    ! Local variables
    REAL(fp) :: n_Factor_TL
    REAL(fp) :: Sum_Fac_TL(0:RTV%n_Angles)
    INTEGER :: i, j, nZ

    nZ = RTV%n_Angles

    ! Normalisation for stream angles
    Sum_Fac_TL(0) = ZERO
    DO i = 1, RTV%n_Streams
      n_Factor_TL = ZERO
      DO j = i,nZ
        n_Factor_TL=n_Factor_TL+(Pff_TL(i,j)+Pbb_TL(i,j))*RTV%COS_Weight(j)
      END DO
      DO j=i,nZ
        Pff_TL(i,j) = Pff_TL(i,j)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k)) - &
                      Pff(i,j)/RTV%n_Factor(i,k)/RTV%n_Factor(i,k)*n_Factor_TL*(ONE-RTV%Sum_Fac(i-1,k)) - &
                      Pff(i,j)/RTV%n_Factor(i,k)*Sum_Fac_TL(i-1)

        Pbb_TL(i,j) = Pbb_TL(i,j)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k)) - &
                      Pbb(i,j)/RTV%n_Factor(i,k)/RTV%n_Factor(i,k)*n_Factor_TL*(ONE-RTV%Sum_Fac(i-1,k)) - &
                      Pbb(i,j)/RTV%n_Factor(i,k)*Sum_Fac_TL(i-1)
      END DO
      Sum_Fac_TL(i)=ZERO
      ! Symmetric condition
      IF( i < nZ ) THEN
        DO j = i+1, nZ
          Pff_TL(j,i) = Pff_TL(i,j)
          Pbb_TL(j,i) = Pbb_TL(i,j)
        END DO
        DO j = 1, i
          Sum_Fac_TL(i)=Sum_Fac_TL(i) + (Pff_TL(j,i+1)+Pbb_TL(j,i+1))*RTV%COS_Weight(j)
        END DO
     END IF
    END DO

    IF( RTV%n_Streams < nZ ) THEN
      ! Sensor viewing angle differs from the Gaussian angles
      n_Factor_TL = Sum_Fac_TL(nZ-1)
      DO j = 1, nZ
        Pff_TL(j,nZ) = Pff_TL(j,nZ)/RTV%n_Factor(nZ,k) - &
                       Pff(j,nZ)/RTV%n_Factor(nZ,k)/RTV%n_Factor(nZ,k)*n_Factor_TL

        Pbb_TL(j,nZ) = Pbb_TL(j,nZ)/RTV%n_Factor(nZ,k) - &
                       Pbb(j,nZ)/RTV%n_Factor(nZ,k)/RTV%n_Factor(nZ,k)*n_Factor_TL
        ! Symmetric condition
        IF( j < nZ ) THEN
          Pff_TL(nZ,j) = Pff_TL(j,nZ)
          Pbb_TL(nZ,j) = Pbb_TL(j,nZ)
        END IF
      END DO
    END IF

  END SUBROUTINE Normalize_Phase_TL

  SUBROUTINE Normalize_Phase_AD( k, RTV, Pff, Pbb, Pff_AD, Pbb_AD )
    ! Arguments
    INTEGER       , INTENT(IN)     :: k
    TYPE(RTV_type), INTENT(IN)     :: RTV
    REAL(fp)      , INTENT(IN)     :: Pff(:,:)
    REAL(fp)      , INTENT(IN)     :: Pbb(:,:)
    REAL(fp)      , INTENT(IN OUT) :: Pff_AD(:,:)
    REAL(fp)      , INTENT(IN OUT) :: Pbb_AD(:,:)
    ! Local variables
    INTEGER :: i, j, nZ
    REAL(fp) :: n_Factor_AD
    REAL(fp) :: Sum_Fac_AD(0:RTV%n_Angles)


    nZ = RTV%n_Angles
    Sum_Fac_AD = ZERO

    n_Factor_AD = ZERO
    IF( RTV%n_Streams < nZ ) THEN
      ! Sensor viewing angle diffs from the Gaussian angles
      DO j = nZ, 1, -1
        ! Symmetric condition
        IF( j < nZ ) THEN
          Pff_AD(j,nZ) = Pff_AD(j,nZ) + Pff_AD(nZ,j)
          Pff_AD(nZ,j) = ZERO
          Pbb_AD(j,nZ) = Pbb_AD(j,nZ) + Pbb_AD(nZ,j)
          Pbb_AD(nZ,j) = ZERO
        END IF

        n_Factor_AD = n_Factor_AD - Pff(j,nZ)/RTV%n_Factor(nZ,k)/RTV%n_Factor(nZ,k)*Pff_AD(j,nZ)
        Pff_AD(j,nZ) = Pff_AD(j,nZ)/RTV%n_Factor(nZ,k)

        n_Factor_AD = n_Factor_AD - Pbb(j,nZ)/RTV%n_Factor(nZ,k)/RTV%n_Factor(nZ,k)*Pbb_AD(j,nZ)
        Pbb_AD(j,nZ) = Pbb_AD(j,nZ)/RTV%n_Factor(nZ,k)
      END DO
      Sum_Fac_AD(nZ-1) = n_Factor_AD
      n_Factor_AD = ZERO
    END IF

    DO i = RTV%n_Streams, 1, -1
      ! Symmetric condition
      IF( i < nZ ) THEN
        DO j = i, 1, -1
          Pbb_AD(j,i+1) = Pbb_AD(j,i+1) + Sum_Fac_AD(i)*RTV%COS_Weight(j)
          Pff_AD(j,i+1) = Pff_AD(j,i+1) + Sum_Fac_AD(i)*RTV%COS_Weight(j)
        END DO
        DO j = nZ,i+1,-1
          Pff_AD(i,j) = Pff_AD(i,j) + Pff_AD(j,i)
          Pff_AD(j,i) = ZERO
          Pbb_AD(i,j) = Pbb_AD(i,j) + Pbb_AD(j,i)
          Pbb_AD(j,i) = ZERO
        END DO
      END IF
      Sum_Fac_AD(i) = ZERO
      DO j = nZ, i, -1
        Sum_Fac_AD(i-1) = Sum_Fac_AD(i-1) - Pbb(i,j)/RTV%n_Factor(i,k)*Pbb_AD(i,j)
        n_Factor_AD = n_Factor_AD -Pbb(i,j)/RTV%n_Factor(i,k)/RTV%n_Factor(i,k) * &
                      Pbb_AD(i,j)*(ONE-RTV%Sum_Fac(i-1,k))
        Pbb_AD(i,j) = Pbb_AD(i,j)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))

        Sum_Fac_AD(i-1) = Sum_Fac_AD(i-1) - Pff(i,j)/RTV%n_Factor(i,k)*Pff_AD(i,j)
        n_Factor_AD = n_Factor_AD -Pff(i,j)/RTV%n_Factor(i,k)/RTV%n_Factor(i,k) * &
                      Pff_AD(i,j)*(ONE-RTV%Sum_Fac(i-1,k))
        Pff_AD(i,j) = Pff_AD(i,j)/RTV%n_Factor(i,k)*(ONE-RTV%Sum_Fac(i-1,k))
      END DO
      DO j = nZ, i, -1
        Pbb_AD(i,j) = Pbb_AD(i,j) + n_Factor_AD*RTV%COS_Weight(j)
        Pff_AD(i,j) = Pff_AD(i,j) + n_Factor_AD*RTV%COS_Weight(j)
      END DO
      n_Factor_AD = ZERO
    END DO
    Sum_Fac_AD(0) = ZERO
  END SUBROUTINE Normalize_Phase_AD

END MODULE Common_RTSolution
