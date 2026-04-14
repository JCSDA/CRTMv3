!
! CRTM_AtmAbsorption
!
! Module containing routines to compute the optical depth profile
! due to gaseous absorption.
!
!
! CREATION HISTORY:
!       Modifed by:     Yong Han, NESDIS/STAR 25-June-2008
!                       yong.han@noaa.gov


MODULE CRTM_AtmAbsorption

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds,                 ONLY: fp
  USE Message_Handler,            ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,            ONLY: ZERO, &
                                        ODAS_ALGORITHM,  ODPS_ALGORITHM, ODSSU_ALGORITHM
  USE CRTM_Atmosphere_Define,     ONLY: CRTM_Atmosphere_type
  USE CRTM_TauCoeff,              ONLY: TC
  USE CRTM_AncillaryInput_Define, ONLY: CRTM_AncillaryInput_type
  USE CRTM_GeometryInfo_Define,   ONLY: CRTM_GeometryInfo_type, &
                                        CRTM_GeometryInfo_GetValue
  USE CRTM_AtmOptics_Define,      ONLY: CRTM_AtmOptics_type
  USE CRTM_Predictor_Define,      ONLY: CRTM_Predictor_type
  USE CRTM_ChannelInfo_Define,    ONLY: CRTM_ChannelInfo_type
  USE iso_c_binding,               ONLY: c_float, c_size_t, c_int
  USE crtm_onnx_interface,         ONLY: crtm_onnx_predict
  ! ODAS modules
  USE ODAS_AtmAbsorption,         ONLY: ODAS_AAVar_type => iVar_type , &
                                        ODAS_Compute_AtmAbsorption   , &
                                        ODAS_Compute_AtmAbsorption_TL, &
                                        ODAS_Compute_AtmAbsorption_AD
  ! ODPS modules
  USE ODPS_AtmAbsorption,         ONLY: ODPS_AAVar_type => iVar_type , &
                                        ODPS_Compute_AtmAbsorption   , &
                                        ODPS_Compute_AtmAbsorption_TL, &
                                        ODPS_Compute_AtmAbsorption_AD
  ! ODSSU modules
  USE ODSSU_AtmAbsorption,        ONLY: ODSSU_AAVar_type => iVar_type , &
                                        ODSSU_Compute_Weights         , &
                                        ODSSU_Compute_AtmAbsorption   , &
                                        ODSSU_Compute_AtmAbsorption_TL, &
                                        ODSSU_Compute_AtmAbsorption_AD
  ! ODZeeman modules
  USE ODZeeman_AtmAbsorption,     ONLY: Zeeman_Compute_AtmAbsorption,    &
                                        Zeeman_Compute_AtmAbsorption_TL, &
                                        Zeeman_Compute_AtmAbsorption_AD, &
                                        Is_Zeeman_Channel

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: iVar_type
  ! Procedures
  PUBLIC :: CRTM_Compute_AtmAbsorption
  PUBLIC :: CRTM_Compute_AtmAbsorption_TL
  PUBLIC :: CRTM_Compute_AtmAbsorption_AD
  PUBLIC :: CRTM_Compute_AtmAbsorption_ONNX_Batch


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


  ! ---------------------------------------------
  ! Structure to hold AtmAbsorption forward model
  ! variables across FWD, TL, and AD calls
  ! ---------------------------------------------
  !:tdoc+:
  TYPE :: iVar_type
    PRIVATE
    TYPE(ODAS_AAVar_type)   :: ODAS
    TYPE(ODPS_AAVar_type)   :: ODPS
    TYPE(ODSSU_AAVar_type)  :: ODSSU
  END TYPE iVar_type
  !:tdoc-:


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_AtmAbsorption
!
! PURPOSE:
!       Subroutine to calculate the layer optical depths due to gaseous
!       absorption for a given sensor and channel and atmospheric profile.
!       It is a wrapper which calls the algorithm-specific routine.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_AtmAbsorption( &
!              SensorIndex   , &  ! Input
!              ChannelIndex  , &  ! Input
!              AncillaryInput, &  ! Input
!              Predictor     , &  ! Input
!              AtmOptics     , &  ! Output
!              iVar            )  ! Internal variable output
!
! INPUTS:
!       SensorIndex:
!         Sensor index id. This is a unique index associated
!         with a (supported) sensor used to access the
!         shared coefficient data for a particular sensor.
!         See the ChannelIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:
!         Channel index id. This is a unique index associated
!         with a (supported) sensor channel used to access the
!         shared coefficient data for a particular sensor's
!         channel.
!         See the SensorIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       AncillaryInput:
!         Structure holding ancillary inputs
!         UNITS:      N/A
!         TYPE:       AncillaryInput_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       Predictor:
!         Structure containing the integrated absorber and
!         predictor profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_Predictor_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       AtmOptics:
!         Structure containing computed optical depth profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_AtmOptics_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:
!         Structure containing internal variables required for
!         subsequent tangent-linear or adjoint model calls.
!         The contents of this structure are NOT accessible
!         outside of this module.
!         UNITS:      N/A
!         TYPE:       iVar_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_AtmAbsorption( &
    SensorIndex   , &  ! Input
    ChannelIndex  , &  ! Input
    AncillaryInput, &  ! Input
    Predictor     , &  ! Input
    AtmOptics     , &  ! Output
    iVar            )  ! Internal variable output
    ! Arguments
    INTEGER                       , INTENT(IN)     :: SensorIndex
    INTEGER                       , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AncillaryInput_type), INTENT(IN)     :: AncillaryInput
    TYPE(CRTM_Predictor_type)     , INTENT(IN OUT) :: Predictor
    TYPE(CRTM_AtmOptics_type)     , INTENT(IN OUT) :: AtmOptics
    TYPE(iVar_type)               , INTENT(IN OUT) :: iVar
    ! Local variables
    INTEGER :: idx

    ! Is this a Zeeman channel?
    idx = TC%ZSensor_LoIndex(SensorIndex)
    IF( idx > 0 )THEN
      IF( Is_Zeeman_Channel(TC%ODZeeman(idx), ChannelIndex) )THEN
        CALL Zeeman_Compute_AtmAbsorption( &
               TC%ODZeeman(idx)  , &  ! Input
               ChannelIndex      , &  ! Input
               Predictor%ODZeeman, &  ! Input
               AtmOptics           )  ! Output
        RETURN
      END IF
    END IF


    ! Call required model
    idx = TC%Sensor_LoIndex(SensorIndex)
    SELECT CASE( TC%Algorithm_ID(SensorIndex) )

      ! ODAS transmittance model
      CASE( ODAS_ALGORITHM )
        CALL ODAS_Compute_AtmAbsorption( &
               TC%ODAS(idx)  , &  ! Input
               ChannelIndex  , &  ! Input
               Predictor%ODAS, &  ! Input
               AtmOptics     , &  ! Output
               iVar%ODAS       )  ! Internal variable output

      ! ODPS transmittance model
      CASE( ODPS_ALGORITHM )
        CALL ODPS_Compute_AtmAbsorption( &
               TC%ODPS(idx)  , &  ! Input
               ChannelIndex  , &  ! Input
               Predictor%ODPS, &  ! Input
               AtmOptics       )  ! Output

      ! SSU instrument specific
      CASE( ODSSU_ALGORITHM )
        CALL ODSSU_Compute_Weights( &
               AncillaryInput%SSU, &  ! Input
               SensorIndex       , &  ! Input
               ChannelIndex      , &  ! Input
               iVar%ODSSU          )  ! Internal variable output

        ! ...Select particular transmittance algorithm for this instrument
        SELECT CASE( TC%ODSSU(idx)%subAlgorithm )
          CASE( ODAS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODAS                , &  ! Input
                   AtmOptics                     , &  ! Output
                   iVar%ODSSU                      )  ! Internal variable output
          CASE( ODPS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODPS                , &  ! Input
                   AtmOptics                     , &  ! Output
                   iVar%ODSSU                      )  ! Internal variable output
        END SELECT
    END SELECT

  END SUBROUTINE CRTM_Compute_AtmAbsorption


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_Compute_AtmAbsorption_TL
!
! PURPOSE:
!       Subroutine to calculate the tangent-linear layer optical depths due
!       to gaseous absorption for a given sensor and channel and atmospheric
!       profile. It is a wrapper which calls the algorithm specific routine.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_AtmAbsorption_TL( &
!              SensorIndex   , &  ! Input
!              ChannelIndex  , &  ! Input
!              Predictor     , &  ! FWD Input
!              Predictor_TL  , &  ! TL Input
!              AtmOptics_TL  , &  ! TL Output
!              iVar            )  ! Internal variable input
!
! INPUTS:
!       SensorIndex:
!         Sensor index id. This is a unique index associated
!         with a (supported) sensor used to access the
!         shared coefficient data for a particular sensor.
!         See the ChannelIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:
!         Channel index id. This is a unique index associated
!         with a (supported) sensor channel used to access the
!         shared coefficient data for a particular sensor's
!         channel.
!         See the SensorIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       Predictor:
!         Structure containing the integrated absorber and
!         predictor profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_Predictor_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       Predictor_TL:
!         Structure containing the tangent-linearintegrated absorber and
!         predictor profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_Predictor_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       iVar:
!         Structure containing internal variables required for
!         subsequent tangent-linear or adjoint model calls.
!         The contents of this structure are NOT accessible
!         outside of this module.
!         UNITS:      N/A
!         TYPE:       iVar_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       AtmOptics_TL:
!         Structure containing the computed tangent-linear
!         optical depth profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_AtmOptics_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN OUT)
!
! :sdoc-:
!------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_AtmAbsorption_TL( &
    SensorIndex   , &  ! Input
    ChannelIndex  , &  ! Input
    Predictor     , &  ! Input
    Predictor_TL  , &  ! Input
    AtmOptics_TL  , &  ! Output
    iVar            )  ! Internal variable input
    ! Arguments
    INTEGER                  , INTENT(IN)     :: SensorIndex
    INTEGER                  , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Predictor_type), INTENT(IN)     :: Predictor
    TYPE(CRTM_Predictor_type), INTENT(IN OUT) :: Predictor_TL
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_TL
    TYPE(iVar_type)          , INTENT(IN)     :: iVar
    ! Local variables
    INTEGER :: idx

    ! Is this a Zeeman channel?
    idx = TC%ZSensor_LoIndex(SensorIndex)
    IF( idx > 0 )THEN
      IF( Is_Zeeman_Channel(TC%ODZeeman(idx), ChannelIndex) )THEN
        CALL Zeeman_Compute_AtmAbsorption_TL( &
               TC%ODZeeman(idx)     , &  ! Input
               ChannelIndex         , &  ! Input
               Predictor%ODZeeman   , &  ! Input
               Predictor_TL%ODZeeman, &  ! Input
               AtmOptics_TL           )  ! Output
        RETURN
      END IF
    END IF


    ! Call required model
    idx = TC%Sensor_LoIndex(SensorIndex)
    SELECT CASE( TC%Algorithm_ID(SensorIndex) )

      ! ODAS transmittance model
      CASE( ODAS_ALGORITHM )
        CALL ODAS_Compute_AtmAbsorption_TL( &
               TC%ODAS(idx)     , &  ! Input
               ChannelIndex     , &  ! Input
               Predictor%ODAS   , &  ! Input
               Predictor_TL%ODAS, &  ! Input
               AtmOptics_TL     , &  ! Output
               iVar%ODAS          )  ! Internal variable input

      ! ODPS transmittance model
      CASE( ODPS_ALGORITHM )
        CALL ODPS_Compute_AtmAbsorption_TL( &
               TC%ODPS(idx)     , &  ! Input
               ChannelIndex     , &  ! Input
               Predictor%ODPS   , &  ! Input
               Predictor_TL%ODPS, &  ! Input
               AtmOptics_TL       )  ! Output

      ! SSU instrument specific
      CASE( ODSSU_ALGORITHM )

        ! ...Select particular transmittance algorithm for this instrument
        SELECT CASE( TC%ODSSU(idx)%subAlgorithm )
          CASE( ODAS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption_TL( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODAS                , &  ! Input
                   Predictor_TL%ODAS             , &  ! Input
                   AtmOptics_TL                  , &  ! Output
                   iVar%ODSSU                      )  ! Internal variable input
          CASE( ODPS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption_TL( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODPS                , &  ! Input
                   Predictor_TL%ODPS             , &  ! Input
                   AtmOptics_TL                  , &  ! Output
                   iVar%ODSSU                      )  ! Internal variable input
        END SELECT
    END SELECT

  END SUBROUTINE CRTM_Compute_AtmAbsorption_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_AtmAbsorption_AD
!
! PURPOSE:
!       Subroutine to calculate the layer optical depth adjoints due to
!       gaseous absorption for a given sensor and channel and atmospheric
!       profile. It is a wrapper which calls the algorithm specific routine.
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_AtmAbsorption_AD( &
!              SensorIndex , &  ! Input
!              ChannelIndex, &  ! Input
!              Predictor   , &  ! FWD Input
!              AtmOptics_AD, &  ! AD  Input
!              Predictor_AD, &  ! AD  Output
!              iVar          )  ! Internal variable input
!
! INPUT ARGUMENTS:
!       SensorIndex:
!         Sensor index id. This is a unique index associated
!         with a (supported) sensor used to access the
!         shared coefficient data for a particular sensor.
!         See the ChannelIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:
!         Channel index id. This is a unique index associated
!         with a (supported) sensor channel used to access the
!         shared coefficient data for a particular sensor's
!         channel.
!         See the SensorIndex argument.
!         UNITS:      N/A
!         TYPE:       INTEGER
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       Predictor:
!         Structure containing the integrated absorber and
!         predictor profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_Predictor_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
!       AtmOptics_AD:
!         Structure containing the computed adjoint optical depth profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_AtmOptics_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:
!         Structure containing internal variables required for
!         subsequent tangent-linear or adjoint model calls.
!         The contents of this structure are NOT accessible
!         outside of this module.
!         UNITS:      N/A
!         TYPE:       iVar_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       Predictor_AD:
!         Structure containing the adjoint integrated absorber and
!         predictor profile data.
!         UNITS:      N/A
!         TYPE:       CRTM_Predictor_type
!         DIMENSION:  Scalar
!         ATTRIBUTES: INTENT(IN OUT)
!
! SIDE EFFECTS:
!       Components of the AtmOptics_AD structure argument are modified
!       in this function.
!
!:sdoc-:
!------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_AtmAbsorption_AD( &
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Predictor   , &  ! FWD Input
    AtmOptics_AD, &  ! AD  Input
    Predictor_AD, &  ! AD  Output
    iVar          )  ! Internal variable input
    ! Arguments
    INTEGER                  , INTENT(IN)     :: SensorIndex
    INTEGER                  , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Predictor_type), INTENT(IN)     :: Predictor
    TYPE(CRTM_AtmOptics_type), INTENT(IN OUT) :: AtmOptics_AD
    TYPE(CRTM_Predictor_type), INTENT(IN OUT) :: Predictor_AD
    TYPE(iVar_type)          , INTENT(IN)     :: iVar
    ! Local variables
    INTEGER :: idx

    ! Is this a Zeeman channel?
    idx = TC%ZSensor_LoIndex(SensorIndex)
    IF( idx > 0 )THEN
      IF( Is_Zeeman_Channel(TC%ODZeeman(idx), ChannelIndex) )THEN
        CALL Zeeman_Compute_AtmAbsorption_AD( &
               TC%ODZeeman(idx)     , &  ! Input
               ChannelIndex         , &  ! Input
               Predictor%ODZeeman   , &  ! Input
               AtmOptics_AD         , &  ! AD Input
               Predictor_AD%ODZeeman  )  ! AD Output
        RETURN
      END IF
    END IF

    ! Call required model
    idx = TC%Sensor_LoIndex(SensorIndex)
    SELECT CASE( TC%Algorithm_ID(SensorIndex) )

      ! ODAS transmittance model
      CASE( ODAS_ALGORITHM )
        CALL ODAS_Compute_AtmAbsorption_AD( &
               TC%ODAS(idx)     , &  ! Input
               ChannelIndex     , &  ! Input
               Predictor%ODAS   , &  ! FWD Input
               AtmOptics_AD     , &  ! AD Input
               Predictor_AD%ODAS, &  ! AD Output
               iVar%ODAS          )  ! Internal variable input

      ! ODPS transmittance model
      CASE( ODPS_ALGORITHM )
        CALL ODPS_Compute_AtmAbsorption_AD( &
               TC%ODPS(idx)     , &  ! Input
               ChannelIndex     , &  ! Input
               Predictor%ODPS   , &  ! FWD Input
               AtmOptics_AD     , &  ! AD Input
               Predictor_AD%ODPS  )  ! AD Output

      ! SSU instrument specific
      CASE( ODSSU_ALGORITHM )

        ! Select particular transmittance algorithm for this instrument
        SELECT CASE( TC%ODSSU(idx)%subAlgorithm )
          CASE( ODAS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption_AD( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODAS                , &  ! FWD Input
                   AtmOptics_AD                  , &  ! AD Input
                   Predictor_AD%ODAS             , &  ! AD Output
                   iVar%ODSSU                      )  ! Internal variable input
          CASE( ODPS_ALGORITHM )
            CALL ODSSU_Compute_AtmAbsorption_AD( &
                   TC%Sensor_LoIndex(SensorIndex), &  ! Input
                   ChannelIndex                  , &  ! Input
                   Predictor%ODPS                , &  ! FWD Input
                   AtmOptics_AD                  , &  ! AD Input
                   Predictor_AD%ODPS             , &  ! AD Output
                   iVar%ODSSU                      )  ! Internal variable input
        END SELECT
    END SELECT

  END SUBROUTINE CRTM_Compute_AtmAbsorption_AD



  SUBROUTINE CRTM_Compute_AtmAbsorption_ONNX_Batch(SensorIndex, ChannelInfo, Atmosphere, GeometryInfo, Optical_Depth)
    INTEGER, INTENT(IN) :: SensorIndex
    TYPE(CRTM_ChannelInfo_type), INTENT(IN) :: ChannelInfo
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: Atmosphere
    TYPE(CRTM_GeometryInfo_type), INTENT(IN) :: GeometryInfo
    REAL(fp), INTENT(OUT) :: Optical_Depth(:,:) ! [n_channels, n_layers]
    
    INTEGER(c_int) :: status
    REAL(c_float), ALLOCATABLE :: batch_features(:,:)
    REAL(c_float), ALLOCATABLE :: batch_outputs(:,:)
    INTEGER :: n_layers, n_channels, j, k, H2O_idx, CO2_idx, O3_idx
    REAL(fp) :: Sensor_Zenith
    REAL(c_float) :: secant

    n_layers = Atmosphere%n_Layers
    n_channels = ChannelInfo%n_Channels
    ALLOCATE(batch_features(6, n_layers))
    ALLOCATE(batch_outputs(n_channels, n_layers))

    CALL CRTM_GeometryInfo_GetValue(GeometryInfo, Sensor_Zenith_Angle=Sensor_Zenith)
    secant = REAL(1.0_fp / COS(Sensor_Zenith * 0.017453292519943_fp), c_float)
    
    H2O_idx = 0; CO2_idx = 0; O3_idx = 0
    DO j = 1, Atmosphere%n_Absorbers
       IF (Atmosphere%Absorber_ID(j) == 1) H2O_idx = j
       IF (Atmosphere%Absorber_ID(j) == 2) CO2_idx = j
       IF (Atmosphere%Absorber_ID(j) == 3) O3_idx = j
    END DO

    DO k = 1, n_layers
       batch_features(1, k) = REAL(LOG10(MAX(Atmosphere%Pressure(k), 1.0e-4_fp)), c_float)
       batch_features(2, k) = REAL(Atmosphere%Temperature(k), c_float)
       batch_features(3, k) = -12.0_c_float; batch_features(4, k) = -12.0_c_float; batch_features(5, k) = -12.0_c_float
       IF (H2O_idx > 0) batch_features(3, k) = REAL(LOG10(MAX(Atmosphere%Absorber(k, H2O_idx), 1.0e-12_fp)), c_float)
       IF (CO2_idx > 0) batch_features(4, k) = REAL(LOG10(MAX(Atmosphere%Absorber(k, CO2_idx), 1.0e-12_fp)), c_float)
       IF (O3_idx > 0) batch_features(5, k) = REAL(LOG10(MAX(Atmosphere%Absorber(k, O3_idx), 1.0e-12_fp)), c_float)
       batch_features(6, k) = secant
       DO j = 1, 6
          batch_features(j, k) = (batch_features(j, k) - REAL(ChannelInfo%ONNX_Mean(j), c_float)) / REAL(ChannelInfo%ONNX_Std(j), c_float)
       END DO
    END DO

    status = crtm_onnx_predict(batch_features, INT(n_layers, c_size_t), 6_c_size_t, &
                               batch_outputs, INT(n_channels, c_size_t))
    
    IF (status == 0) THEN
       DO k = 1, n_layers
          DO j = 1, n_channels
             Optical_Depth(j, k) = REAL(-LOG(MAX(batch_outputs(j, k), 1.0e-20_c_float)), fp)
          END DO
       END DO
    END IF
    DEALLOCATE(batch_features, batch_outputs)
  END SUBROUTINE CRTM_Compute_AtmAbsorption_ONNX_Batch

END MODULE CRTM_AtmAbsorption
