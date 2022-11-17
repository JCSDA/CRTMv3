!----------------------------------------------------------------------------------
!:sdoc+:
!
!
! CRTM_FASTEM_MODULE
!
! This module is provided to "wrap" all the existing CRTM FASTEM
! versions and provide a general interface to simplify integration into
! the main CRTM_SfcOptics module.
!
!
! CREATION HISTORY:
!      Written by:      Paul van Delst, 25-Jun-2005
!                       paul.vandelst@noaa.gov
!      Modified by:     Ming Chen, 08-Feb-2015
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!
!

MODULE CRTM_FASTEM_MODULE

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds,          ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler,   ONLY: SUCCESS
  USE CRTM_FASTEM_Parameters,         ONLY: SET, NOT_SET, &
                                      ZERO, ONE, &
                                      MAX_N_ANGLES, &
                                      N_STOKES => MAX_N_STOKES
  USE CRTM_LowFrequency_MWSSEM, ONLY: LF_MWSSEM_type => iVar_type, &
                                      LowFrequency_MWSSEM, &
                                      LowFrequency_MWSSEM_TL, &
                                      LowFrequency_MWSSEM_AD
  USE CRTM_Fastem1,             ONLY: Fastem1
  USE CRTM_Fastem3,             ONLY: Fastem3
  USE CRTM_FastemXX,            ONLY: FastemX_type => iVar_type, &
                                     Compute_FastemX    => Compute_FastemXX,   &
                                     Compute_FastemX_TL => Compute_FastemXX_TL,&
                                     Compute_FastemX_AD => Compute_FastemXX_AD
  USE FASTEM_Coeff_Reader,   MWwaterC => CSEM_MWwaterC

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  PUBLIC :: iVar_type
  ! Science routines
  PUBLIC :: CRTM_FASTEM_Emiss
  PUBLIC :: Compute_FASTEM_SfcOptics
  PUBLIC :: Compute_FASTEM_SfcOptics_TL
  PUBLIC :: Compute_FASTEM_SfcOptics_AD
  PUBLIC :: CRTM_FASTEM_Init 
  PUBLIC :: CRTM_FASTEM_Destroy
  PUBLIC :: CSEM_MWwaterCoeff_INIT 
  
  ! -----------------
  ! Module parameters
  ! -----------------
  ! RCS Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_MW_Water_SfcOptics.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Low frequency model threshold
  REAL(fp), PARAMETER :: LOW_F_THRESHOLD = 20.0_fp ! GHz

  LOGICAL,SAVE :: CSEM_MWwaterCoeff_INIT = .FALSE.
  
  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! FastemX model internal variable structure
    TYPE(FastemX_type), DIMENSION(MAX_N_ANGLES) :: FastemX_Var
    ! Low frequency model internal variable structure
    TYPE(LF_MWSSEM_type), DIMENSION(MAX_N_ANGLES) :: LF_MWSSEM_Var
    LOGICAL :: Is_LOW_Freq = .FALSE. 
    ! Fastem outputs
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dTs        = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEH_dWindSpeed = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dTs        = ZERO
    REAL(fp), DIMENSION(MAX_N_ANGLES) :: dEV_dWindSpeed = ZERO
  END TYPE iVar_type

CONTAINS
!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_FASTEM_Emiss
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a water surface and at SINGLE frequency channel and SINGLE
!       receiving angle
!
!       This function is a wrapper of different FASTEM versions
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_FASTEM_Emiss(    
!                              Frequency                              , &  ! Input
!                              Angle                                  , &  ! Input
!                              Water_Temperature                      , &  ! Input
!                              Salinity                               , &  ! Input
!                              Wind_Speed                             , &  ! Input
!                              Wind_Direction                         , &  ! Input
!                              Emissivity                             , &  ! Output
!                              Reflectivity                           , &  ! Output
!                              FASTEM_Version                         , &  ! Input
!                              Sensor_Azimuth_Angle                   , &  ! Input
!                              Transmittance)                           &  ! Input
! INPUTS:
!      Frequency:        Frequency
!                        UNITS:      GHz
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Angle:            Receiving Angle
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Water_Temperature:Water surface temperature
!                        UNITS:      K
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Salinity:         Water salinity
!                        UNITS:      Per thousand
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Speed:       Ocean surface wind speed
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Direction:   Wind direction
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      FASTEM_Version:   FASTEM version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Sensor_Azimuth_Angle: Sensor azimuth angle
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: OPTIONAL INTENT(IN)
!          
!      Transmittance :   Atmospheric transmittance
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: OPTIONAL INTENT(IN)

! OUTPUTS:
!      Emissivity:       Full Stokes surface emissivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Array(4)
!                        ATTRIBUTES: INTENT(OUT)
!
!      Reflectivity:     Full Stokes surface reflectivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  ARRAY(4)
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
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_FASTEM_Emiss( &
           Frequency                              , &  ! Input
           Angle                                  , &  ! Input
           Water_Temperature                      , &  ! Input
           Salinity                               , &  ! Input
           Wind_Speed                             , &  ! Input
           Wind_Direction                         , &  ! Input
           Emissivity                             , &  ! Output
           Reflectivity                           , &  ! Output
           FASTEM_Version                         , &  ! Input
           Sensor_Azimuth_Angle                   , &  ! Input
           Transmittance)                           &  ! Input
  RESULT( err_stat )
    
    ! Arguments
    REAL(fp), INTENT(IN)   :: Frequency 
    REAL(fp), INTENT(IN)   :: Angle
    REAL(fp), INTENT(IN)   :: Water_Temperature
    REAL(fp), INTENT(IN)   :: Salinity
    REAL(fp), INTENT(IN)   :: Wind_Speed
    REAL(fp), INTENT(IN)   :: Wind_Direction
    REAL(fp), INTENT(OUT)  :: Reflectivity(N_STOKES)
    REAL(fp), INTENT(OUT)  :: Emissivity(N_STOKES)
    INTEGER,  INTENT(IN)   :: FASTEM_Version 
    REAL(fp), INTENT(IN)   :: Sensor_Azimuth_Angle
    REAL(fp), INTENT(IN)   :: Transmittance

    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics'
    ! Local variables
    REAL(fp) :: Azimuth_Angle 
    TYPE(FastemX_type )  :: iVar_FastemX
    TYPE(LF_MWSSEM_type) :: iVar_LF_MWSSEM
    REAL(fp) :: dEH_dWindSpeed , dEV_dWindSpeed              
    INTEGER :: nZ = 2

    ! Set up
    err_stat = SUCCESS
    Reflectivity = ZERO
    
    ! Compute the surface optical parameters
    IF( FASTEM_Version == 6   .OR. FASTEM_Version == 5 ) THEN
       ! FastemX model
       Azimuth_Angle = Wind_Direction - Sensor_Azimuth_Angle

       CALL Compute_FastemX(                    &
               MWwaterC                       , & ! Input model coefficients
               Frequency                      , & ! Input
               nZ                             , &  ! Input
               Angle                          , & ! Input
               Water_Temperature              , & ! Input
               Salinity                       , & ! Input
               Wind_Speed                     , & ! Input
               iVar_FastemX                   , & ! Internal variable output
               Emissivity(:)                  , & ! Output
               Reflectivity(:)                , & ! Output
               Azimuth_Angle                  , & ! Optional input
               Transmittance                    )  ! Optional input
    ELSE IF( FASTEM_Version == 3 ) THEN
        CALL Fastem3(                            &
                     Frequency,                 & ! INPUT
                     Angle,                     & ! INPUT
                     Sensor_Azimuth_Angle,      & ! INPUT in degree
                     Water_Temperature,         & ! INPUT
                     Wind_Speed,                & ! INPUT
                     Wind_Direction,            & ! INPUT in degree
                     Transmittance,             & ! INPUT
                     Fastem_Version,            & ! INPUT
                     Emissivity(:),             & ! OUTPUT
                     Reflectivity(:))             ! OUTPUT

    ELSE
        ! Low frequency model coupled with Fastem1
        IF( Frequency < LOW_F_THRESHOLD ) THEN
           ! Call the low frequency model
           CALL LowFrequency_MWSSEM( &
                 Frequency                    , & ! Input
                 Angle                        , & ! Input
                 Water_Temperature            , & ! Input
                 Salinity                     , & ! Input
                 Wind_Speed                   , & ! Input
                 Emissivity                   , & ! Output
                 iVar_LF_MWSSEM                 ) ! Internal variable output	 
            Reflectivity = ONE-Emissivity
        ELSE
          ! Call Fastem1
           CALL Fastem1( Frequency            , & ! Input
                  Angle                       , & ! Input
                  Water_Temperature           , & ! Input
                  Wind_Speed                  , & ! Input
                  Emissivity                  , & ! Output
                  dEH_dWindSpeed              , & ! Output
                  dEV_dWindSpeed                ) ! Output
            Reflectivity = ONE-Emissivity 

        END IF

    END IF

  END FUNCTION CRTM_FASTEM_Emiss

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      Compute_FASTEM_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a water surface and at SINGLE frequency channel and MULTIPLE
!       receiving angle
!
!       This function is a wrapper of different FASTEM versions
!
! CALLING SEQUENCE:
!       Error_Status = Compute_FASTEM_SfcOptics(    
!                                 Frequency                              , &  ! Input
!                                 Angle                                  , &  ! Input
!                                 Water_Temperature                      , &  ! Input
!                                 Salinity                               , &  ! Input
!                                 Wind_Speed                             , &  ! Input
!                                 Wind_Direction                         , &  ! Input
!                                 iVar                                   , &  ! Input
!                                 Emissivity                             , &  ! Output
!                                 Reflectivity                           , &  ! Output
!                                 FASTEM_Version                         , &  ! Input
!                                 Sensor_Azimuth_Angle                   , &  ! Input
!                                 Transmittance)                           &  ! Input
! INPUTS:
!      Frequency:        Frequency
!                        UNITS:      GHz
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Angle:            Receiving Angle
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(IN)
!
!      Water_Temperature:Water surface temperature
!                        UNITS:      K
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Salinity:         Water salinity
!                        UNITS:      Per thousand
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Speed:       Ocean surface wind speed
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Direction:   Wind direction
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      FASTEM_Version:   FASTEM version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Sensor_Azimuth_Angle: Sensor azimuth angle
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: OPTIONAL INTENT(IN) 
!          
!      Transmittance :   Atmospheric transmittance
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: OPTIONAL INTENT(IN)
!
!
! OUTPUTS:
!      Emissivity:       Full Stokes surface emissivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(OUT)
!
!      Reflectivity:     Full Stokes surface reflectivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(OUT)
!
!      iVar:             Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!      iVar:             Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
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
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_FASTEM_SfcOptics(  &
           Frequency                              , &  ! Input
           Angles                                 , &  ! Input
           Water_Temperature                      , &  ! Input
           Salinity                               , &  ! Input
           Wind_Speed                             , &  ! Input
           Wind_Direction                         , &  ! Input
           iVar                                   , &  ! InOUT
           Emissivity                             , &  ! Output
           Reflectivity                           , &  ! Output
           FASTEM_Version                         , &  ! Input
           Sensor_Azimuth_Angle                   , &  ! Input
           Transmittance)                           &  ! Input
  RESULT( err_stat )
    
    ! Arguments
    REAL(fp), INTENT(IN)   :: Frequency 
    REAL(fp), INTENT(IN)   :: Angles(:)
    REAL(fp), INTENT(IN)   :: Water_Temperature
    REAL(fp), INTENT(IN)   :: Salinity
    REAL(fp), INTENT(IN)   :: Wind_Speed
    REAL(fp), INTENT(IN)   :: Wind_Direction
    REAL(fp), INTENT(OUT)  :: Reflectivity(:,:)
    REAL(fp), INTENT(OUT)  :: Emissivity(:,:)
    INTEGER,  INTENT(IN)   :: FASTEM_Version 
    REAL(fp), INTENT(IN)   :: Sensor_Azimuth_Angle
    REAL(fp), INTENT(IN)   :: Transmittance
    TYPE(iVar_type),  INTENT(OUT) :: iVar

    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics'
    ! Local variables
    INTEGER  :: i, nZ
    REAL(fp) :: Azimuth_Angle 

    ! Set up
    err_stat = SUCCESS
    Reflectivity = ZERO
    nZ = size(Angles)
    ! Compute the surface optical parameters
    IF( FASTEM_Version == 6 .OR. FASTEM_Version == 5) THEN
      ! FastemX model
      Azimuth_Angle = Wind_Direction - Sensor_Azimuth_Angle

      DO i = 1, nZ
        CALL Compute_FastemX( &
               MWwaterC                               , &  ! Input model coefficients
               Frequency                              , &  ! Input
               nZ                                     , &  ! Input
               Angles(i)                              , &  ! Input
               Water_Temperature                      , &  ! Input
               Salinity                               , &  ! Input
               Wind_Speed                             , &  ! Input
               iVar%FastemX_Var(i)                    , &  ! Internal variable output
               Emissivity(i,:)                        , &  ! Output
               Reflectivity(i,:)                      , &  ! Output
               Azimuth_Angle                          , &  ! Optional input
               Transmittance                            )  ! Optional input
      END DO
    ELSE IF( FASTEM_Version == 3 ) THEN

       DO i = 1, nZ
          CALL Fastem3(                         &
                     Frequency,                 & ! INPUT
                     Angles(i),                 & ! INPUT
                     Sensor_Azimuth_Angle,      & ! INPUT in degree
                     Water_Temperature,         & ! INPUT
                     Wind_Speed,                & ! INPUT
                     Wind_Direction,            & ! INPUT in degree
                     Transmittance,             & ! INPUT
                     Fastem_Version,            & ! INPUT
                     Emissivity(i,:),           & ! OUTPUT
                     Reflectivity(i,:))           ! OUTPUT

       END DO
   ELSE

      ! Low frequency model coupled with Fastem1
      IF( Frequency < LOW_F_THRESHOLD ) THEN
        iVar%Is_LOW_Freq = .TRUE.
        ! Call the low frequency model
        DO i = 1, nZ
          CALL LowFrequency_MWSSEM( &
                 Frequency                            , &  ! Input
                 Angles(i)                            , &  ! Input
                 Water_Temperature                    , &  ! Input
                 Salinity                             , &  ! Input
                 Wind_Speed                           , &  ! Input
                 Emissivity(i,:)                      , &  ! Output
                 iVar%LF_MWSSEM_Var(i)                  )  ! Internal variable output
          Reflectivity(i,1:2) = ONE - Emissivity(i,1:2)
        END DO
      ELSE
        ! Call Fastem1
        DO i = 1, nZ
          CALL Fastem1( Frequency                     , &  ! Input
                        Angles(i)                     , &  ! Input
                        Water_Temperature             , &  ! Input
                        Wind_Speed                    , &  ! Input
                        Emissivity(i,:)               , &  ! Output
                        iVar%dEH_dWindSpeed(i)        , &  ! Output
                        iVar%dEV_dWindSpeed(i)          )  ! Output
          Reflectivity(i,1:2) = ONE- Emissivity(i,1:2)
        END DO
      END IF

    END IF



  END FUNCTION Compute_FASTEM_SfcOptics

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      Compute_FASTEM_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a water surface.
!
!       This function is a wrapper of different FASTEM versions
!
! CALLING SEQUENCE:
!       Error_Status = Compute_FASTEM_SfcOptics_TL(    
!                            Water_Temperature_TL                   , &  ! Input
!                            Salinity_TL                            , &  ! Input
!                            Wind_Speed_TL                          , &  ! Input
!                            Wind_Direction_TL                      , &  ! Input
!                            iVar                                   , &  ! Input
!                            EmissivityTL                           , &  ! Output
!                            Reflectivity_TL                        , &  ! Output
!                            FASTEM_Version)                          &  ! Input
! INPUTS:
!
!      Water_Temperature_TL:Water surface temperature tangent-linear 
!                        UNITS:      K
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Salinity_TL:      Water salinity tangent-linear 
!                        UNITS:      Per thousand
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Speed_TL:    Ocean surface wind speed tangent-linear
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      Wind_Direction_TL:Wind direction tangent-linear
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      FASTEM_Version:   FASTEM version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      iVar:             Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)

! OUTPUTS:
!      Emissivity_TL:    Full Stokes surface emissivity tangent-linear
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(OUT)
!
!      Reflectivity_TL:  Full Stokes surface reflectivity tangent-linear
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
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
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_FASTEM_SfcOptics_TL(   &
           Water_Temperature_TL                 , &  ! Input
           Salinity_TL                          , &  ! Input
           Wind_Speed_TL                        , &  ! Input
           Wind_Direction_TL                    , &  ! Input
           Transmittance_TL                     , &  ! Input
           iVar                                 , &  ! InOUT
           Emissivity_TL                        , &  ! Output
           Reflectivity_TL                      , &  ! Output
           FASTEM_Version )                       &  ! Input
  RESULT( err_stat )
    ! Arguments
    REAL(fp), INTENT(IN)   :: Water_Temperature_TL
    REAL(fp), INTENT(IN)   :: Salinity_TL
    REAL(fp), INTENT(IN)   :: Wind_Speed_TL
    REAL(fp), INTENT(IN)   :: Wind_Direction_TL
    REAL(fp), INTENT(IN)   :: Transmittance_TL
    REAL(fp), INTENT(OUT)  :: Reflectivity_TL(:,:)
    REAL(fp), INTENT(OUT)  :: Emissivity_TL(:,:)
    INTEGER,  INTENT(IN)   :: FASTEM_Version 

    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics_TL'
    ! Local variables
    INTEGER :: i, nZ
   
    ! Set up
    err_stat = SUCCESS
    nZ = size(Emissivity_TL,1) 
    Reflectivity_TL = ZERO
    ! ...Retrieve data from structures

    ! Compute the tangent-linear surface optical parameters
    IF( FASTEM_Version == 6 .OR. FASTEM_Version == 5) THEN

      ! FastemX model
      DO i = 1, nZ
         CALL Compute_FastemX_TL( &
               MWwaterC                            , &  ! Input model coefficients
               Water_Temperature_TL                , &  ! TL Input
               Salinity_TL                         , &  ! TL Input
               Wind_Speed_TL                       , &  ! TL Input
               iVar%FastemX_Var(i)                 , &  ! Internal variable input
               Emissivity_TL(i,:)                  , &  ! TL Output
               Reflectivity_TL(i,:)                , &  ! TL Output
               Azimuth_Angle_TL = Wind_Direction_TL, &  ! Optional TL input
               Transmittance_TL = Transmittance_TL   )  ! Optional TL input

      END DO
    ELSE IF( FASTEM_Version == 3 ) THEN
        Emissivity_TL =  ZERO                 
        Reflectivity_TL = ZERO 

    ELSE

      ! Low frequency model coupled with Fastem1
      IF( iVar%Is_LOW_Freq) THEN
        ! Call the low frequency model
        DO i = 1, nZ
            CALL LowFrequency_MWSSEM_TL(              &
                 Water_Temperature_TL             , &  ! TL  Input
                 Salinity_TL                      , &  ! TL  Input
                 Wind_Speed_TL                    , &  ! TL  Input
                 Emissivity_TL(i,:)               , &  ! TL  Output
                 iVar%LF_MWSSEM_Var(i)              )  ! Internal variable input
          Reflectivity_TL(i,1) = -Emissivity_TL(i,1)
          Reflectivity_TL(i,2) = -Emissivity_TL(i,2)
        END DO
      ELSE
        ! Call Fastem1
        DO i = 1, nZ
          Emissivity_TL(i,2) = (iVar%dEH_dTs(i)*Water_Temperature_TL) + &
                                        (iVar%dEH_dWindSpeed(i)*Wind_Speed_TL)
          Emissivity_TL(i,1) = (iVar%dEV_dTs(i)*Water_Temperature_TL) + &
                                         (iVar%dEV_dWindSpeed(i)*Wind_Speed_TL)
          Reflectivity_TL(i,1) = -Emissivity_TL(i,1)
          Reflectivity_TL(i,2) = -Emissivity_TL(i,2)
        END DO
      END IF
    END IF

  END FUNCTION Compute_FASTEM_SfcOptics_TL

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      Compute_FASTEM_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a water surface.
!
!
!       This function is a wrapper of different FASTEM versions
!
! CALLING SEQUENCE:
!       Error_Status = Compute_FASTEM_SfcOptics_AD(    
!                              Emissivity_AD                          , &  ! Output
!                              Reflectivity_AD                        , &  ! Output
!                              Water_Temperature_AD                   , &  ! Input
!                              Salinity_AD                            , &  ! Input
!                              Wind_Speed_AD                          , &  ! Input
!                              Wind_Direction_AD                      , &  ! Input
!                              iVar                                   , &  ! Input
!                              FASTEM_Version)                          &  ! Input
! INPUTS:
!
!      Emissivity_AD:    Full Stokes surface emissivity adjoint 
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(INOUT)
!
!      Reflectivity_AD:  Full Stokes surface reflectivity adjoint 
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape ARRAY
!                        ATTRIBUTES: INTENT(INOUT)

!      FASTEM_Version:   FASTEM version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!      iVar:             Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_MW_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!
!      Water_Temperature_AD:Water surface temperature adjoint 
!                        UNITS:      K
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
!
!      Salinity_AD:      Water salinity adjoint 
!                        UNITS:      Per thousand
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
!
!      Wind_Speed_AD:    Ocean surface wind speed djoint 
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
!
!      Wind_Direction_AD:Wind direction adjoint 
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
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
!!      Note the INTENT on the input SfcOptics_AD argument is IN OUT rather
!       than just OUT. This is necessary because components of this argument
!       may need to be zeroed out upon output.
!
!       Note the INTENT on the output Surface_AD argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input. To prevent memory leaks, the IN OUT INTENT is a must.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION Compute_FASTEM_SfcOptics_AD(  &
           Emissivity_AD                        , &  ! Input
           Reflectivity_AD                      , &  ! Input
           Water_Temperature_AD                 , &  ! Output
           Salinity_AD                          , &  ! Output
           Wind_Speed_AD                        , &  ! Output
           Wind_Direction_AD                    , &  ! Output
           Transmittance_AD                     , &  ! Output
           iVar                                 , &  ! InOUT
           FASTEM_Version)                        &  ! Input
  RESULT( err_stat )
  
    ! Arguments
    REAL(fp), INTENT(INOUT) :: Reflectivity_AD(:,:)
    REAL(fp), INTENT(INOUT) :: Emissivity_AD(:,:)
    REAL(fp), INTENT(INOUT) :: Transmittance_AD
    REAL(fp), INTENT(OUT)   :: Water_Temperature_AD
    REAL(fp), INTENT(OUT)   :: Salinity_AD
    REAL(fp), INTENT(OUT)   :: Wind_Speed_AD
    REAL(fp), INTENT(OUT)   :: Wind_Direction_AD
    INTEGER , INTENT(IN)    :: FASTEM_Version 


    TYPE(iVar_type),              INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_MW_Water_SfcOptics_AD'
    ! Local variables
    INTEGER :: i, j, nZ
    REAL(fp) :: Azimuth_Angle_AD

    ! Set up
    err_stat = SUCCESS
    nZ = size(Emissivity_AD,1) 


    ! Compute the adjoint surface optical parameters
    IF(FASTEM_Version == 6 .OR. FASTEM_Version == 5) THEN

      ! FastemX model
      Azimuth_Angle_AD = ZERO
      DO i = 1, nZ
        CALL Compute_FastemX_AD( &
               MWwaterC                             , &  ! Input model coefficients
               Emissivity_AD(i,:)                   , &  ! AD Input
               Reflectivity_AD(i,:)                 , &  ! AD Input
               iVar%FastemX_Var(i)                        , &  ! Internal variable input
               Water_Temperature_AD                 , &  ! AD Output
               Salinity_AD                          , &  ! AD Output
               Wind_Speed_AD                        , &  ! AD Output
               Azimuth_Angle_AD = Azimuth_Angle_AD  , &  ! Optional AD Output
               Transmittance_AD = Transmittance_AD    )  ! Optional AD Output
      END DO
      Wind_Direction_AD = Wind_Direction_AD + Azimuth_Angle_AD

    ELSE IF( FASTEM_Version == 3 ) THEN
        Emissivity_AD = ZERO               
        Reflectivity_AD = ZERO
 
    ELSE

      ! Low frequency model coupled with Fastem1
      IF( iVar%Is_LOW_Freq ) THEN
        ! Call the low frequency model
        DO i = 1, nZ
          Emissivity_AD(i,1) = Emissivity_AD(i,1)-Reflectivity_AD(i,1)
          Emissivity_AD(i,2) = Emissivity_AD(i,2)-Reflectivity_AD(i,2)
          CALL LowFrequency_MWSSEM_AD(                &
                 Emissivity_AD(i,:)                 , &  ! AD  Input
                 Water_Temperature_AD               , &  ! AD  Output
                 Salinity_AD                        , &  ! AD  Output
                 Wind_Speed_AD                      , &  ! AD  Output
                 iVar%LF_MWSSEM_Var(i)                ) ! Internal variable input
        END DO
      ELSE
        ! Call Fastem1
        DO i = nZ, 1, -1
          DO j = 1, 2
            Emissivity_AD(i,1:2) = Emissivity_AD(i,1:2) - Reflectivity_AD(i,1:2)
            Reflectivity_AD(i,1:2) = ZERO
          END DO
          ! Vertical polarisation component
          Water_Temperature_AD  = Water_Temperature_AD + &
                                         (iVar%dEV_dTs(i)*Emissivity_AD(i,1))
          Wind_Speed_AD         = Wind_Speed_AD + &
                                         (iVar%dEV_dWindSpeed(i)*Emissivity_AD(i,1))
          Emissivity_AD(i,1)  = ZERO
          ! Horizontal polarization component
          Water_Temperature_AD  = Water_Temperature_AD + &
                                         (iVar%dEH_dTs(i)*Emissivity_AD(i,2))
          Wind_Speed_AD         = Wind_Speed_AD + &
                                         (iVar%dEH_dWindSpeed(i)*Emissivity_AD(i,2))
          Emissivity_AD(i,2)  = ZERO
        END DO
      END IF
    END IF

    Reflectivity_AD = ZERO

  END FUNCTION Compute_FASTEM_SfcOptics_AD
!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      CRTM_FASTEM_Init
!
! PURPOSE:
!       Function to load FASTEM coefficient NETDCF files
!
!       This function must be called before calling other FASTEM functions
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_FASTEM_Init(    
!                              MWwaterCoeff_File                      , &  ! Output
!                              Version)                                 &  ! Input
! INPUTS:
!
!      MWwaterCoeff_File:FASTEM coefficient file (full-path)
!                        UNITS:      N/A
!                        TYPE:       CHARACTER
!                        DIMENSION:  SCALAR
!                        ATTRIBUTES: INTENT(INOUT)
!
!
!      Version:          version index
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!
!      CSEM_MWwaterCoeff_INIT: Global-scope variable
!                        UNITS:      N/A
!                        TYPE:       LOGICAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: Global 
!
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
!
!:sdoc-:
!----------------------------------------------------------------------------------
  
  FUNCTION CRTM_FASTEM_Init( &
    MWwaterCoeff_File,       &
    Version)                 &
  RESULT( Error_Status  )
    CHARACTER(LEN=*)  :: MWwaterCoeff_File
    INTEGER,INTENT(IN)  :: Version
    INTEGER  :: Error_Status
    LOGICAL  :: Quiet=.TRUE.
    !PRINT*, 'loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File), Version
    Error_Status = CSEM_MWwaterCoeff_Load( &
                 TRIM(MWwaterCoeff_File), &
                 Version=Version, Quiet = Quiet )

    IF ( Error_Status /= SUCCESS ) THEN
        PRINT*, 'Error loading MWwaterCoeff data from '//TRIM(MWwaterCoeff_File)
        STOP
    END IF
    MWwaterC%Version = Version 
    CSEM_MWwaterCoeff_INIT = .TRUE.
  END FUNCTION CRTM_FASTEM_Init 
  
  FUNCTION CRTM_FASTEM_Destroy () RESULT( Error_Status  )
    INTEGER  :: Error_Status
    Error_Status = SUCCESS
    !PRINT*, ' Closeing MWwaterCoeff data ....'
    IF(.NOT. CSEM_MWwaterCoeff_INIT) THEN 
       RETURN
    ENDIF
    Error_Status  = CSEM_MWwaterCoeff_CleanUp( )

    IF ( Error_Status /= SUCCESS ) THEN
        PRINT*, 'Error Closeing MWwaterCoeff data ....'
        STOP
    END IF
    CSEM_MWwaterCoeff_INIT = .FALSE.
  END FUNCTION CRTM_FASTEM_Destroy
 
END MODULE CRTM_FASTEM_MODULE
