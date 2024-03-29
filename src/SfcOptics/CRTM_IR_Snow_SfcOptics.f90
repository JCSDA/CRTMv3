!
! CRTM_IR_Snow_SfcOptics
!
! Module to compute the surface optical properties for snow surfaces at
! infrared frequencies required for determining the snow surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to "wrap" their existing
! codes inside the provided functions to simplify integration into
! the main CRTM_SfcOptics module.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 23-Jun-2005
!                       paul.vandelst@noaa.gov
!
!       Modified by:   Cheng Dang, 31-May-2022
!                      dangch@ucar.edu
!

MODULE CRTM_IR_Snow_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds               , ONLY: fp
  USE Message_Handler          , ONLY: SUCCESS, Display_Message
  USE Spectral_Units_Conversion, ONLY: Inverse_cm_to_Micron
  USE CRTM_Parameters          , ONLY: ZERO, ONE, MAX_N_ANGLES
  USE CRTM_SpcCoeff            , ONLY: SC, SpcCoeff_IsSolar
  USE CRTM_Surface_Define      , ONLY: CRTM_Surface_type
  USE CRTM_GeometryInfo_Define , ONLY: CRTM_GeometryInfo_type
  USE CRTM_SfcOptics_Define    , ONLY: CRTM_SfcOptics_type
  USE CRTM_SEcategory          , ONLY: SEVar_type => iVar_type, &
                                       SEcategory_Emissivity
  USE CRTM_IRsnowCoeff         , ONLY: CRTM_IRsnowCoeff_IsLoaded, &
                                       CRTM_IRsnowCoeff_SE_IsLoaded, &
                                       IRsnowC, &
                                       IRsnowC_SE
  USE CRTM_IRSnowEM            , ONLY: IRsnowVar_type => iVar_type, &
                                       CRTM_Compute_IRSnowEM, &
                                       CRTM_Compute_IRSnowEM_TL, &
                                       CRTM_Compute_IRSnowEM_AD
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  PUBLIC :: iVar_SE_type, iVar_type
  ! Science routines
  PUBLIC :: Compute_IR_Snow_SfcOptics
  PUBLIC :: Compute_IR_Snow_SfcOptics_TL
  PUBLIC :: Compute_IR_Snow_SfcOptics_AD


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_SE_type
    PRIVATE
    TYPE(SEVar_type) :: sevar
  END TYPE iVar_SE_type

  TYPE :: iVar_type
    PRIVATE
    TYPE(IRsnowVar_type) :: irsnowvar
  END TYPE iVar_type



CONTAINS


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_IR_Snow_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at infrared
!       frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_IR_Snow_SfcOptics( &
!                        Surface     , &
!                        SensorIndex , &
!                        ChannelIndex, &
!                        SfcOptics   , &
!                        iVar          )
!
! INPUTS:
!       Surface:         Structure containing the surface state data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
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
!                        transfer calculation. On input the Angle component
!                        is assumed to contain data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the module containing this procedure.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT as it is assumed to contain some data upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_IR_Snow_SfcOptics( &
    Surface     , &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics   , &  ! Output
    iVar_SE     , &  ! Internal variable output
    iVar        ) &  ! Internal variable output
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics
    TYPE(iVar_SE_type),           INTENT(IN OUT) :: iVar_SE
    TYPE(iVar_type),              INTENT(IN OUT) :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_IR_Snow_SfcOptics'
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER :: j, nZ
    REAL(fp) :: frequency, emissivity
    LOGICAL  :: isSEcategory, isIRsnowC

    ! Set up
    err_stat     = SUCCESS
    ! ...Retrieve frequency data from structures
    frequency    = SC(SensorIndex)%Wavenumber(ChannelIndex)
    ! ...Short name for angle dimensions
    nZ = SfcOptics%n_Angles
    isSEcategory = CRTM_IRsnowCoeff_SE_IsLoaded()
    isIRsnowC    = CRTM_IRsnowCoeff_IsLoaded()


    ! Compute Lambertian surface emissivity
    IF ( isSEcategory ) THEN
      err_stat = SEcategory_Emissivity( &
                   IRsnowC_SE       , &  ! Input
                   frequency        , &  ! Input
                   Surface%Snow_Type, &  ! Input
                   emissivity       , &  ! Output
                   iVar_SE%sevar      )  ! Internal variable output
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error occurred in SEcategory_Emissivity()'
        CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
      END IF

      ! Solar direct component
      IF ( SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex) ) THEN
        SfcOptics%Direct_Reflectivity(:,1) = ONE - emissivity
      END IF

      ! Fill the return emissivity arrays
      SfcOptics%Emissivity(1:SfcOptics%n_Angles,1) = emissivity

    ELSE IF ( isIRsnowC ) THEN
      err_stat = CRTM_Compute_IRSnowEM(&
                   IRsnowC                      , &  ! Input
                   Surface%Snow_Temperature     , &  ! Input
                   Surface%Snow_Grain_Size      , &  ! Input
                   frequency                    , &  ! Input
                   SfcOptics%Angle(1:nZ)        , &  ! Input
                   iVar%irsnowvar               , &  ! Internal variable output
                   SfcOptics%Emissivity(1:nZ,1)   )  ! Output
      IF ( err_stat /= SUCCESS ) THEN
        msg = 'Error occurred in CRTM_Compute_IRSnowEM()'
        CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
      END IF

      ! Compute the solar direct BRDF
      IF ( SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex) ) THEN
        ! Cheng: placeholder for BRDF
        SfcOptics%Direct_Reflectivity(1:nZ,1) = ZERO
      END IF

    END IF

    ! Fill the return reflectivity arrays
    DO j = 1, nZ
      SfcOptics%Reflectivity(j,1,j,1) = ONE - SfcOptics%Emissivity(j,1)
    END DO

  END FUNCTION Compute_IR_Snow_SfcOptics


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_IR_Snow_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at infrared frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       err_stat = Compute_IR_Snow_SfcOptics_TL( &
!                           Surface     , &  ! Input
!                           SfcOptics   , &  ! Input
!                           Surface_TL  , &  ! Input
!                           GeometryInfo, &  ! Input
!                           SensorIndex , &  ! Input
!                           ChannelIndex, &  ! Input
!                           SfcOptics_TL, &  ! Output
!                           iVar        ) &  ! Internal variable input
!
! OUTPUTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! INPUTS:
!       Surface:         CRTM_Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
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
!       SfcOptics:       CRTM_SfcOptics structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
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
!                        outside of the CRTM_IR_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       SfcOptics_TL:    CRTM_SfcOptics structure containing the tangent-linear
!                        surface optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CRTM_SfcOptics_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       err_stat:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the input SfcOptics_TL argument is IN OUT rather
!       than just OUT as it may be defined upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_IR_Snow_SfcOptics_TL( &
    Surface     , &  ! Input
    SfcOptics   , &  ! Input
    Surface_TL  , &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    SfcOptics_TL, &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface_TL
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_TL
    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat, j, nZ
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_IR_Snow_SfcOptics_TL'
    ! Local variables
    LOGICAL  :: isIRsnowC

    ! Set up
    err_stat  = SUCCESS
    ! ...Short name for angle dimensions
    nZ = SfcOptics%n_Angles
    ! ...Snow coefficient model
    isIRsnowC = CRTM_IRsnowCoeff_IsLoaded()

    ! Compute the tangent-linear surface optical parameters
    IF ( isIRsnowC ) THEN

      ! Compute tangent-linear IR snow surface emissivity
      err_stat = CRTM_Compute_IRSnowEM_TL( &
                   IRsnowC                        , &  ! Input model coefficients
                   Surface_TL%Snow_Temperature    , &  ! Input
                   Surface_TL%Snow_Grain_Size     , &  ! Input
                   iVar%irsnowvar                 , &  ! Internal variable input
                   SfcOptics_TL%Emissivity(1:nZ,1)  )  ! Output
      IF ( err_stat /= SUCCESS ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Error computing Tangent_linear IR snow surface emissivity', &
                              err_stat  )
        RETURN
      END IF

      ! Compute the tangent-linear solar direct BRDF
      IF ( SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex) ) THEN
        ! Cheng: placeholder for BRDF
          SfcOptics_TL%Direct_Reflectivity(1:nZ,1) = ZERO
      END IF

      ! Surface reflectance (currently assumed to be specular ALWAYS)
      DO j = 1, nZ
        SfcOptics_TL%Reflectivity(j,1,j,1) = -SfcOptics_TL%Emissivity(j,1)
      END DO

    ELSE

      ! No TL component for SEcategory files
      SfcOptics_TL%Reflectivity        = ZERO
      SfcOptics_TL%Direct_Reflectivity = ZERO
      SfcOptics_TL%Emissivity          = ZERO

    END IF

  END FUNCTION Compute_IR_Snow_SfcOptics_TL


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_IR_Snow_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at infrared frequencies over a snow surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       Error_Status = Compute_IR_Snow_SfcOptics_AD( &
!                        Surface     , &
!                        SfcOptics   , &
!                        SfcOptics_AD, &
!                        GeometryInfo, &
!                        SensorIndex , &
!                        ChannelIndex, &
!                        Surface_AD  , &
!                        iVar          )
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
!                        surface optical properties required for the adjoint
!                        radiative transfer calculation.
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
!                        outside of the CRTM_IR_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Surface_AD:      CRTM_Surface structure containing the adjoint
!                        surface state data.
!                        UNITS:      N/A
!                        TYPE:       CRTM_Surface_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       Error_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
!
! COMMENTS:
!       Note the INTENT on the input adjoint arguments are IN OUT regardless
!       of their specification as "input" or "output". This is because these
!       arguments may contain information on input, or need to be zeroed on
!       output (or both).
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_IR_Snow_SfcOptics_AD( &
    Surface     , &  ! Input
    SfcOptics   , &  ! Input
    SfcOptics_AD, &  ! Input
    GeometryInfo, &  ! Input
    SensorIndex , &  ! Input
    ChannelIndex, &  ! Input
    Surface_AD  , &  ! Output
    iVar        ) &  ! Internal variable input
  RESULT( err_stat )
    ! Arguments
    TYPE(CRTM_Surface_type),      INTENT(IN)     :: Surface
    TYPE(CRTM_SfcOptics_type),    INTENT(IN)     :: SfcOptics
    TYPE(CRTM_SfcOptics_type),    INTENT(IN OUT) :: SfcOptics_AD
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER,                      INTENT(IN)     :: SensorIndex
    INTEGER,                      INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_Surface_type),      INTENT(IN OUT) :: Surface_AD
    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat, j, nZ
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_IR_Snow_SfcOptics_AD'
    ! Local variables
    LOGICAL  :: isIRsnowC

    ! Set up
    err_stat = SUCCESS
    ! ...Short name for angle dimensions
    nZ = SfcOptics%n_Angles
    ! ...Snow coefficient model
    isIRsnowC = CRTM_IRsnowCoeff_IsLoaded()

    ! Compute the adjoint surface optical parameters
    IF ( isIRsnowC ) THEN
      ! Surface reflectance (currently assumed to be specular ALWAYS)
      DO j = nZ, 1, -1
        SfcOptics_AD%Emissivity(j,1) = SfcOptics_AD%Emissivity(j,1) - &
                                       SfcOptics_AD%Reflectivity(j,1,j,1)
        SfcOptics_AD%Reflectivity(j,1,j,1) = ZERO
      END DO

      ! Compute the adjoint solar direct BRDF
      IF ( SpcCoeff_IsSolar(SC(SensorIndex), ChannelIndex=ChannelIndex) ) THEN
        ! Cheng: placeholder for BRDF
        SfcOptics_AD%Direct_Reflectivity(1:nZ,1) = ZERO
      END IF

      ! Compute sdjoint IRSSEM sea surface emissivity
      err_stat = CRTM_Compute_IRSnowEM_AD( &
                   IRsnowC                       , &  ! Input model coefficients
                   SfcOptics_AD%Emissivity(1:nZ,1), &  ! Input
                   iVar%irsnowvar                 , &  ! Internal Variable Input
                   Surface_AD%Snow_Grain_Size     , &  ! Output
                   Surface_AD%Snow_Temperature      )  ! Output
      IF ( err_stat /= SUCCESS ) THEN
        CALL Display_Message( ROUTINE_NAME, &
                              'Error computing Adjoint IR sea surface emissivity', &
                              err_stat )
        RETURN
      END IF

    ELSE

      ! No AD component for SEcategory files
      SfcOptics_AD%Reflectivity        = ZERO
      SfcOptics_AD%Direct_Reflectivity = ZERO
      SfcOptics_AD%Emissivity          = ZERO

    END IF

  END FUNCTION Compute_IR_Snow_SfcOptics_AD

END MODULE CRTM_IR_Snow_SfcOptics
