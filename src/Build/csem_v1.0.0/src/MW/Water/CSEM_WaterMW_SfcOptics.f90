!----------------------------------------------------------------------------------
!:sdoc+:
!
! CSEM_WaterMW_SfcOptics
!
! Module to compute the surface optical properties of Ocean surfaces at microwave
! frequencies, which is required for determining the surface contribution
! to the overall radiative transfer. 
!
! This module provide a general container for developers to integrate individual 
! models into the model repository, and a generic interface for 
! the upper-level applications to access the model options in the repository.
! 
! This module implementes the FWD, TL and AD functions for variational data 
! assimilation application and retrieval applications. It replacess the similar
! functions as those in the original CRTM surface module "CRTM_WaterMW_SfcOptics", 
! which was Written by Paul van Delst, 23-Jun-2005
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 28-Feb-2016
!                       ming.chen@noaa.gov
!
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!
!

MODULE CSEM_WaterMW_SfcOptics
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler

  USE CSEM_Define
  USE CSEM_Model_Manager
  USE CRTM_FASTEM_Parameters,ONLY: MAX_n_Angles, MAX_N_STOKES
  
  USE CRTM_FASTEM_MODULE,ONLY: FASTEM_iVar_Type => iVar_type, &
                         Compute_FASTEM_SfcOptics, &
                         Compute_FASTEM_SfcOptics_TL, &
                         Compute_FASTEM_SfcOptics_AD, &
                         CRTM_FASTEM_Init, & 
                         CRTM_FASTEM_Destroy, &
                         CSEM_MWwaterCoeff_INIT
  USE RTTOV_FASTEM_MODULE, ONLY :  RTTOV_iVar_Type => iVar_type, &
                         Compute_RTTOV_Fastem, &
                         Compute_RTTOV_Fastem_TL, &
                         Compute_RTTOV_Fastem_AD 
  
 !Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  REAL(fp), PARAMETER ::  ONE = 1.0_fp, ZERO = 0.0_fp
  ! Data types
  ! Science routines
  PUBLIC :: CSEM_Compute_WaterMW_SfcOptics
  PUBLIC :: CSEM_Compute_WaterMW_SfcOptics_TL
  PUBLIC :: CSEM_Compute_WaterMW_SfcOptics_AD
  PUBLIC :: iVar_Type
    
  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    TYPE(FASTEM_iVar_Type) :: FASTEM
    TYPE(RTTOV_iVar_Type)  :: RTTOV(MAX_n_Angles)
    REAL(fp) :: Transmittance_TL
    REAL(fp) :: Transmittance_AD
  END TYPE iVar_type

  TYPE(CSEM_Model_ID), SAVE  :: MODEL

CONTAINS

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_WaterMW_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at Infrared
!       frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
! PURPOSE:
!       Function to compute the surface emissivity  at microwave
!       frequencies over a Water surface.
!
!       This function encapsulates all available  CSEM WaterMW emissivity models,and
!       provides a genereic interface for the upper-level user application.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterMW_SfcOptics( 
!                                        Surface                  ,&  ! Input
!                                        SfcOptics                ,&  ! Output
!                                        Options                     ,&  ! Input extension 
!                                        iVar)                     &  ! output
!
! INPUTS:
!
!       Surface:         CSEM Water Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Water_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Options:         CSEM input extension  structure containing additional inputs,
!                        e.g., Gelocation/Time Metadata,  Sensor observations, et al,
!                        which may be used by emprirical and semi-empirical models.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_WaterMW_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!
! OUTPUTS:
!       SfcOptics:       CSEM output structure containing the surface optical
!                        properties required for the radiative transfer calculation,
!                        including wavelength and geometry structure as the inputs
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(IN OUT)
! FUNCTION RESULT:
!       IO_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the CSEM_Exception_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
!
! COMMENTS:
!                       CSEM infrastructure data structures are defined in CSEM_DEFINE.f90
!       
!
!:sdoc-:
!-----------------------------------------------------------------------------------------------


  FUNCTION CSEM_Compute_WaterMW_SfcOptics( &
      Surface,                  &  ! Input
      SfcOptics,                &  ! Output
      Options,                  &  ! AdditionalInput 
      iVar)                     &  ! Output
    RESULT (IO_Status)
    
    TYPE(CSEM_Water_Surface),  INTENT(IN)     :: Surface
    TYPE(CSEM_Options_Type),   INTENT(IN)     :: Options
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT)  :: SfcOptics
    TYPE(iVar_type),           INTENT(OUT)    :: ivar
  
    ! Function result
    INTEGER ::  IO_Status 

    ! Local parameters
    CHARACTER(LEN=256) :: DB_PATH
    CHARACTER(LEN=256) :: MWwaterCoeff_File
    INTEGER  :: FASTEM_Version 
    REAL(fp) :: Emiss(4), Refl(4)
    REAL(fp) :: Emissivity(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Reflectivity(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Rel_Azimuth
    INTEGER  ::  i, nZ
 
    IO_Status  = 0
    nZ = SfcOptics%n_Angles
 
  
    Model= Inq_Model_Option("MW_WATER") 
    ! use path in CSEM alg registor file
    DB_PATH = TRIM(MODEL%DATA_PATH)
    
    !Rel_Azimuth = Surface%Wind_Direction - Options%SensorObs%Azimuth_Angle + 180.0
    Rel_Azimuth = Surface%Wind_Direction - SfcOptics%Sensor_Azimuth_Angle 
  
    IF(INDEX(TRIM(MODEL%NAME),'RTTOV_FASTEM') > 0) THEN
         FASTEM_Version = 6 
         IF(TRIM(MODEL%NAME) .EQ. "RTTOV_FASTEM_V5") THEN
           FASTEM_Version = 5
         ENDIF
         IF(TRIM(MODEL%NAME) .EQ. "RTTOV_FASTEM_V4") THEN
           FASTEM_Version = 4
         ENDIF
         !Rel_Azimuth = Surface%Wind_Direction - Options%SensorObs%Azimuth_Angle

         DO i = 1, nZ
            CALL Compute_RTTOV_Fastem(fastem_version,       &         ! Input
                            SfcOptics%Frequency,            &         ! Input
                            SfcOptics%Angle(i),             &         ! Input
                            Surface%Water_Temperature,      &         ! Input
                            Surface%Salinity,               &         ! Input
                            Surface%Wind_Speed,             &         ! Input
                            Emiss,                          &         ! Output
                            Refl,                           &         ! Output
                            Options%Atmos%Transmittance,       &         ! Input, may not be used
                            Rel_Azimuth, iVar=iVar%RTTOV(i))          ! Input, may not be used
                            !Supply_Foam_Fraction, & ! Optional input
                            !Foam_Fraction)          ! Optional input
 
            SfcOptics%Emissivity(i, 1:4)     =  Emiss(:)      
            SfcOptics%Reflectivity(i,1,i,1)  =  Refl(1)          
            SfcOptics%Reflectivity(i,2,i,2)  =  Refl(2)          
            SfcOptics%Reflectivity(i,3,i,3)  =  Refl(3)          
            SfcOptics%Reflectivity(i,4,i,4)  =  Refl(4)          
         END DO

       ELSE IF(INDEX(TRIM(MODEL%NAME),'NESDIS_FASTEM') > 0) THEN
 
         IF(TRIM(MODEL%NAME) .EQ. "NESDIS_FASTEM_V3") THEN
           FASTEM_Version = 3
         ENDIF

         IF(TRIM(MODEL%NAME) .EQ. "NESDIS_FASTEM_V5") THEN
           FASTEM_Version = 5
           MWwaterCoeff_File = 'FASTEM5.MWwater.EmisCoeff.nc'
         ENDIF
         IF(TRIM(MODEL%NAME) .EQ. "NESDIS_FASTEM_V6") THEN
           FASTEM_Version = 6
           MWwaterCoeff_File = 'FASTEM6.MWwater.EmisCoeff.nc'
         ENDIF
         IF((.NOT. CSEM_MWwaterCoeff_INIT) .AND. FASTEM_Version /= 3) THEN
           IO_Status = CRTM_FASTEM_Init( &
                TRIM(ADJUSTL(DB_PATH))//TRIM(ADJUSTL(MWwaterCoeff_File)), FASTEM_Version )

           IF ( IO_Status /= SUCCESS ) THEN
             PRINT*, 'Error loading MWwaterCoeff data from '// &
             TRIM(ADJUSTL(DB_PATH))//TRIM(ADJUSTL(MWwaterCoeff_File))
             STOP
           END IF
         ENDIF
 
         IO_Status = Compute_FASTEM_SfcOptics(                  &
                       SfcOptics%Frequency                    , &  ! Input
                       SfcOptics%Angle                        , &  ! Input
                       Surface%Water_Temperature              , &  ! Input
                       Surface%Salinity                       , &  ! Input
                       Surface%Wind_Speed                     , &  ! Input
                       Surface%Wind_Direction                 , &  ! Input
                       iVar%FASTEM                            , &  ! InOUT
                       Emissivity                             , &  ! Output
                       Reflectivity                           , &  ! Output
                       FASTEM_Version                         , &  ! Input
                       SfcOptics%Sensor_Azimuth_Angle         , &  ! Input
                       Options%Atmos%Transmittance                 )  ! Input
         DO i = 1, nZ
           SfcOptics%Emissivity(i, 1:4)     =  Emissivity(i,:)     
           SfcOptics%Reflectivity(i,1,i,1)  =  Reflectivity(i,1)   
           SfcOptics%Reflectivity(i,2,i,2)  =  Reflectivity(i,2)      
           SfcOptics%Reflectivity(i,3,i,3)  =  Reflectivity(i,3)        
           SfcOptics%Reflectivity(i,4,i,4)  =  Reflectivity(i,4)       
         END DO
      ENDIF
 

  END FUNCTION CSEM_Compute_WaterMW_SfcOptics

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_WaterMW_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterMW_SfcOptics_TL( 
!                                    Surface_TL,      &
!                                    SfcOptics_TL,    &
!                                    iVar)
!
! INPUTS:
!      Surface_TL:   Structure containing the tangent-linear surface
!                        state inputs required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Water_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_WaterMW_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!      SfcOptics_TL: Structure containing the tangent-linear surface
!                        optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
! FUNCTION RESULT:
!       IO_Status:    The return value is an integer defining the error status.
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

  FUNCTION CSEM_Compute_WaterMW_SfcOptics_TL( &
      Surface_TL,               &  ! Input
      Atmos_TL,                 &  ! Input
      SfcOptics_TL,             &  ! Output
      iVar        )             &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Water_Surface),         INTENT(IN)     :: Surface_TL
    TYPE(CSEM_Atmosphere_Parameters), INTENT(IN)     :: Atmos_TL
    TYPE(CSEM_SfcOptics_Type),        INTENT(INOUT)  :: SfcOptics_TL
    TYPE(iVar_type),                  INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: IO_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_WaterMW_SfcOptics_TL'
    ! Local variables
    INTEGER ::  FASTEM_Version = 5

    REAL(fp) :: Emissivity_TL(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Reflectivity_TL(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Emiss_TL(4), Refl_TL(4)
    INTEGER ::  i,nZ


    ! Set up
    IO_Status = SUCCESS
    nZ = SfcOptics_TL%n_Angles

    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity = ZERO

    Model= Inq_Model_Option("MW_WATER") 

    IF(INDEX(TRIM(MODEL%NAME),'RTTOV_FASTEM') > 0) THEN
         DO i = 1, nZ
            CALL Compute_RTTOV_Fastem_TL(         &  ! Input
                 Surface_TL%Water_Temperature,    &  ! Input
                 Surface_TL%Salinity,             &  ! Input
                 Surface_TL%Wind_Speed,           &  ! Input
                 Emiss_TL,                        &  ! Output
                 Refl_TL,                         &  ! Output
                 Atmos_TL%Transmittance ,         &  ! Input, may not be used
                 Surface_TL%Wind_Direction,       &  ! Input, may not be used
                 Surface_TL%Foam_Fraction,        &  ! Input, may not be used
                 ivar=iVar%RTTOV(i))                 ! Optional input

 
           SfcOptics_TL%Emissivity(i,1:4)      =  Emiss_TL(1:4)      
           SfcOptics_TL%Reflectivity(i,1,i,1)  =  Refl_TL(1)          
           SfcOptics_TL%Reflectivity(i,2,i,2)  =  Refl_TL(2)          
           SfcOptics_TL%Reflectivity(i,3,i,3)  =  Refl_TL(3)          
           SfcOptics_TL%Reflectivity(i,4,i,4)  =  Refl_TL(4)          

          END DO
 
  
    ELSE IF (INDEX(TRIM(MODEL%NAME),'NESDIS_FASTEM') > 0) THEN
         IF(TRIM(MODEL%NAME) .EQ. "NESDIS_FASTEM_V6") THEN
            FASTEM_Version = 6 
          END IF
          IO_Status = Compute_FASTEM_SfcOptics_TL( &
               Surface_TL%Water_Temperature                 , &  ! Input
               Surface_TL%Salinity                          , &  ! Input
               Surface_TL%Wind_Speed                        , &  ! Input
               Surface_TL%Wind_Direction                    , &  ! Input
               Atmos_TL%Transmittance                       , &  ! Input
               iVar%FASTEM                                  , &  ! InOUT
               Emissivity_TL                                , &  ! Output
               Reflectivity_TL                              , &  ! Output
               FASTEM_Version )                                  ! Input
         DO i = 1, nZ
           SfcOptics_TL%Emissivity(i, 1:4)     =  Emissivity_TL(i,1:4)     
           SfcOptics_TL%Reflectivity(i,1,i,1)  =  Reflectivity_TL(i,1)   
           SfcOptics_TL%Reflectivity(i,2,i,2)  =  Reflectivity_TL(i,2)      
           SfcOptics_TL%Reflectivity(i,3,i,3)  =  Reflectivity_TL(i,3)        
           SfcOptics_TL%Reflectivity(i,4,i,4)  =  Reflectivity_TL(i,4)       
           END DO
    END IF
 
 END FUNCTION CSEM_Compute_WaterMW_SfcOptics_TL

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      CSEM_Compute_WaterMW_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
!       NB: CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO AD
!           COMPONENTS IN THE MW Water SFCOPTICS COMPUTATIONS.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterMW_SfcOptics_AD(
!                                  SfcOptics_AD,         &
!                                  Surface_AD,           &
!                                  Atmos_AD,             &  
!                                  iVar)
!
! INPUTS:
!       SfcOptics_AD:    Structure containing the adjoint surface optical
!                        properties required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_WaterMW_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
! OUTPUTS:
!
!       Surface_AD:     Structure containing the adjoint surface state
!                        properties required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_Water_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
!       Atmos_AD:        Structure containing the adjoint atmosphere otpical
!                        properties required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_Atmosphere_Parameters
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       IO_Status:       The return value is an integer defining the error status.
!                        The error codes are defined in the Message_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the input adjoint arguments are IN OUT regardless
!       of their specification as "input" or "output". This is because these
!       arguments may contain information on input, or need to be zeroed on
!       output (or both).
!
!:sdoc-:
!----------------------------------------------------------------------------------


  FUNCTION CSEM_Compute_WaterMW_SfcOptics_AD( &
      SfcOptics_AD,             &  ! Input
      Surface_AD,               &  ! Output
      Atmos_AD,                 &  ! Output
      iVar        )             &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Water_Surface),         INTENT(INOUT) :: Surface_AD
    TYPE(CSEM_Atmosphere_Parameters), INTENT(INOUT) :: Atmos_AD
    TYPE(CSEM_SfcOptics_Type),        INTENT(INOUT) :: SfcOptics_AD
    TYPE(iVar_type),                  INTENT(IN)    :: iVar
 
    ! Function result
    INTEGER :: IO_Status 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_WaterMW_SfcOptics_AD'
    ! Local variables
    INTEGER ::  FASTEM_Version = 5
    REAL(fp) :: Emissivity_AD(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Reflectivity_AD(MAX_n_Angles, MAX_N_STOKES)
    REAL(fp) :: Emiss_AD(4), Refl_AD(4)
    REAL(fp) :: Temperature_ad, Salinity_ad, Wind_Speed_ad
    REAL(fp) :: Transmittance_ad, Rel_Azimuth_ad
    REAL(fp) :: Foam_Fraction_ad
    INTEGER  :: i, nZ

    ! Set up
    IO_Status = SUCCESS
    nZ = SfcOptics_AD%n_Angles
 
    Model = Inq_Model_Option("MW_WATER") 

    IF (INDEX(TRIM(MODEL%NAME),'RTTOV_FASTEM') > 0) THEN
         DO i = 1, nZ
           Emiss_AD   = SfcOptics_AD%Emissivity(i,1:4) 
           Refl_AD(1) = SfcOptics_AD%Reflectivity(i,1,i,1)
           Refl_AD(2) = SfcOptics_AD%Reflectivity(i,2,i,2)
           Refl_AD(3) = SfcOptics_AD%Reflectivity(i,3,i,3)
           Refl_AD(4) = SfcOptics_AD%Reflectivity(i,4,i,4)
   
           temperature_ad   = ZERO; salinity_ad    = ZERO; wind_Speed_ad    = ZERO
           transmittance_ad = ZERO; Rel_Azimuth_ad = ZERO; foam_Fraction_ad = ZERO
   
            CALL Compute_RTTOV_Fastem_AD(         &  ! Input
                 Emiss_AD,                        &  ! Inout
                 Refl_AD,                         &  ! Inout
                 temperature_ad,                  &  ! Inout
                 salinity_ad,                     &  ! Inout
                 wind_Speed_ad,                   &  ! Inout
                 transmittance_ad,                &  ! inout, may not be used
                 Rel_Azimuth_ad,                  &  ! Inout, may not be used
                 foam_Fraction_ad,                &  ! Inout, may not be used
                 ivar=iVar%RTTOV(i))                 ! Optional input
 
                 Surface_AD%Water_Temperature = Surface_AD%Water_Temperature + Temperature_ad
                 Surface_AD%Salinity          = Surface_AD%Salinity + Salinity_ad
                 Surface_AD%Wind_Speed        = Surface_AD%Wind_Speed +  Wind_Speed_ad
                 Atmos_AD%Transmittance       = Atmos_AD%Transmittance + Transmittance_ad
                 Surface_AD%Wind_Direction    = Surface_AD%Wind_Direction  + Rel_Azimuth_ad
                 Surface_AD%Foam_Fraction     = Surface_AD%Foam_Fraction + Foam_Fraction_ad
          END DO
             
    ELSE IF (INDEX(TRIM(MODEL%NAME),'NESDIS_FASTEM') > 0) THEN
       IF(TRIM(MODEL%NAME) .EQ. "NESDIS_FASTEM_V6") THEN
           FASTEM_Version = 6 
        END IF
       
        DO i = 1, nZ
           Emissivity_AD(i,1:4) = SfcOptics_AD%Emissivity(i,1:4)       
           Reflectivity_AD(i,1) = SfcOptics_AD%Reflectivity(i,1,i,1)        
           Reflectivity_AD(i,2) = SfcOptics_AD%Reflectivity(i,2,i,2)      
           Reflectivity_AD(i,3) = SfcOptics_AD%Reflectivity(i,3,i,3)    
           Reflectivity_AD(i,4) = SfcOptics_AD%Reflectivity(i,4,i,4)
        END DO
     
        IO_Status = Compute_FASTEM_SfcOptics_AD(       &
               Emissivity_AD(1:nZ, 1:4)                 , &  ! INOUT
               Reflectivity_AD(1:nZ, 1:4)               , &  ! INOUT
               Surface_AD%Water_Temperature             , &  ! INOUT
               Surface_AD%Salinity                      , &  ! INOUT
               Surface_AD%Wind_Speed                    , &  ! INOUT
               Surface_AD%Wind_Direction                , &  ! INOUT
               Atmos_AD%Transmittance                   , &  ! INOUT
               iVar%FASTEM                              , &  ! InOUT
               FASTEM_Version)                               ! Input
    END IF
 
    SfcOptics_AD%Emissivity      = ZERO      
    SfcOptics_AD%Reflectivity    = ZERO        
     
 
  END FUNCTION CSEM_Compute_WaterMW_SfcOptics_AD

END MODULE CSEM_WaterMW_SfcOptics
