!----------------------------------------------------------------------------------
!:sdoc+:
!
! CSEM_WaterIR_SfcOptics
!
! Module to compute the surface optical properties of Water surfaces at infrared 
! frequencies, which is required for determining the Water surface contribution
! to the overall radiative transfer. 
!
! This module provide a general container for developers to integrate individual 
! models into the model repository, and a generic interface for 
! the upper-level applications to access the model options in the repository.
! 
! This module implementes the FWD, TL and AD functions for variational data 
! assimilation application and retrieval applications. It replacess the similar
! functions as those in the original CRTM surface module "CRTM_WaterIR_SfcOptics", 
! which was Written by Paul van Delst, 23-Jun-2005
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 28-Feb-2016
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!
!
!

MODULE CSEM_WaterIR_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Define
  USE CSEM_Model_Manager
  USE CSEM_Exception_Handler
  
  USE NPOESS_LUT_Module
  USE NESDIS_WaterIR_PhyModel, ONLY: IRSSEM_BRDF_Type =>iVar_type, &
                                NESDIS_IRSSEM_BRDF,&
                                NESDIS_IRSSEM_BRDF_TL,&
                                NESDIS_IRSSEM_BRDF_AD,&
                                NESDIS_IRSSEM_Setup, &
                                NESDIS_IRSSEM_Close, &
                                NESDIS_IRSSEM_Initialized
  USE NESDIS_WaterIR_PhyModel_V2, ONLY: IRSSEM_BRDF_V2_Type =>iVar_type, &
                                NESDIS_IRSSEM_BRDF_V2,&
                                NESDIS_IRSSEM_BRDF_V2_TL,&
                                NESDIS_IRSSEM_BRDF_V2_AD,&
                                NESDIS_IRSSEM_V2_Setup, &
                                NESDIS_IRSSEM_V2_Close, &
                                NESDIS_IRSSEM_V2_Initialized
  !Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  ! Science routines
  PUBLIC :: CSEM_Compute_WaterIR_SfcOptics
  PUBLIC :: CSEM_Compute_WaterIR_SfcOptics_TL
  PUBLIC :: CSEM_Compute_WaterIR_SfcOptics_AD
  PUBLIC :: iVar_type

 
  REAL(fp), PARAMETER ::  ONE = 1.0_fp, ZERO = 0.0_fp
  REAL(fp), PARAMETER ::  PI  = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER ::  D2R  = PI/180.0_fp
  INTEGER,  PARAMETER ::  MAX_n_Angles = 15
  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    INTEGER :: n_Angles
    ! Variables in routines rough sea BRDF
    TYPE(IRSSEM_BRDF_V2_type)   :: IRSSEM_BRDF_V2
    TYPE(IRSSEM_BRDF_type)      :: IRSSEM_BRDF
    CHARACTER(LEN=256)          :: AlgID
    LOGICAL                     :: VALID_RTTOV_ALG = .FALSE.
  END TYPE iVar_type

  TYPE(CSEM_MODEL_ID),SAVE :: MODEL

CONTAINS

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_WaterIR_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at Infrared
!       frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
! PURPOSE:
!       Function to compute the surface emissivity  at Ifrared
!       frequencies over a Water surface.
!
!       This function encapsulates all available  CSEM WaterIR emissivity models,and
!       provides a genereic interface for the upper-level user application.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterIR_SfcOptics( 
!                                        Surface                  ,&  ! Input
!                                        SfcOptics                ,&  ! Output
!                                        Options )                       ! OptionalInput 
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
!       Options:            CSEM input extension structure containing additional inputs,
!                        e.g., Gelocation/Time Metadata,  Sensor observations, et al,
!                        which may be used by emprirical and semi-empirical models.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
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

  FUNCTION CSEM_Compute_WaterIR_SfcOptics(  &
      Surface,                              &  ! Input
      SfcOptics,                            &  ! Output
      Options,                              &  ! Optional Input
      iVar)                                 &  ! Input
    RESULT (IO_Status)

    TYPE(CSEM_Water_Surface),  INTENT(IN)    ::  Surface
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT) ::  SfcOptics
    TYPE(CSEM_Options_Type),   INTENT(IN)    ::  Options
    TYPE(iVar_type),           INTENT(OUT)   ::  iVar
    ! Function result
    INTEGER ::  IO_Status
 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_WaterIR_SfcOptics'
    CHARACTER(LEN=256) :: Message 
    ! Local variables
    CHARACTER(LEN=256) :: DB_PATH, IRSSEM_File
    REAL(fp) ::  Emiss,  Direct_Refl
    REAL(fp) ::  Sensor_Zenith_Radian, Sensor_Azimuth_Radian
    REAL(fP) ::  Solar_Zenith_Radian,  Solar_Azimuth_Radian

    INTEGER  ::  i, nZ, stype   
    CHARACTER(LEN=100)    ::   CRTM_Sensor_Id  
    
    REAL(fP) ::  Emissivity(MAX_n_Angles), Reflectivity(MAX_n_Angles)
    REAL(fP) ::  Direct_Reflectivity(MAX_n_Angles)
    
    IO_Status  = SUCCESS
    
    Sensor_Zenith_Radian   =  D2R * SfcOptics%Sensor_Zenith_Angle
    Sensor_Azimuth_Radian  =  D2R * SfcOptics%Sensor_Azimuth_Angle
    Solar_Zenith_Radian    =  D2R * SfcOptics%Source_Zenith_Angle
    Solar_Azimuth_Radian   =  D2R * SfcOptics%Source_Azimuth_Angle
    Direct_Refl    =  ZERO

    nZ = SfcOptics%n_Angles
   
    MODEL = Inq_Model_Option("IR_WATER") 
    ! use path in CSEM alg registor file
    DB_PATH = TRIM(MODEL%DATA_PATH)
    iVar%AlgID = TRIM(MODEL%NAME)

    CRTM_Sensor_Id  = TRIM(Options%SensorObs%Sensor_Id)
 
    IF(TRIM(MODEL%NAME) .EQ. 'NPOESS_LUT' ) THEN
        IF(Surface%Water_Type /= 1) THEN
          PRINT*,'NPOESS_LUT Water Type 1'
          IO_Status  = FAILURE 
          RETURN
        ENDIF
        Stype =  Surface%Water_Type + 20
        IO_Status = NPOESS_LUT_Emiss(    &
            SfcOptics%Wavenumber,        & ! INPUT, wavelength in micrometer
            Emiss ,                      & ! OUTPUT, surface emissivity (0 - 1)
            Stype)                         ! INPUT, surface type (1)

        SfcOptics%Emissivity           =   Emiss
        SfcOptics%Reflectivity         =   ONE-Emiss
        SfcOptics%Direct_Reflectivity  =   Direct_Refl

    ELSE IF((TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2CM") .OR. &
          (TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2EK") ) THEN

       IF(.NOT. NESDIS_IRSSEM_V2_Initialized()) THEN
          IF(TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2CM") THEN
            IRSSEM_File = 'Nalli.IRwater.EmisCoeff.V2.CM.nc'
          ELSE
          !"NESDIS_IRW_WuSmith" 
            IRSSEM_File = 'Nalli.IRwater.EmisCoeff.V2.EK.nc'
          END IF
          IO_Status = NESDIS_IRSSEM_V2_Setup(TRIM(DB_PATH)//'/'//TRIM(IRSSEM_File))
       END IF
       
       IF(SfcOptics%Is_Solar) THEN
       
          IO_Status = NESDIS_IRSSEM_BRDF_V2(               &
                  SfcOptics%Wavenumber                   , &  ! Input
                  Surface%Wind_Speed                     , &  ! Input
                  Surface%Water_Temperature              , &  ! Input
                  SfcOptics%Angle                        , &  ! Input
                  Sensor_Zenith_Radian                   , &  ! Input
                  Sensor_Azimuth_Radian                  , &  ! Input
                  Solar_Zenith_Radian                    , &  ! Input
                  Solar_Azimuth_Radian                   , &  ! Input
                  Emissivity(1:nZ)                       , &  ! Output
                  Reflectivity(1:nZ)                     , &  ! Output
                  Direct_Reflectivity (1:nZ)             , &  ! Output
                  iVar%IRSSEM_BRDF_V2 )                    
        ELSE
          IO_Status = NESDIS_IRSSEM_BRDF_V2( &
                  SfcOptics%Wavenumber                   , &  ! Input
                  Surface%Wind_Speed                     , &  ! Input
                  Surface%Water_Temperature              , &  ! Input
                  SfcOptics%Angle                        , &  ! Input
                  Emissivity(1:nZ)                       , &  ! Output
                  Reflectivity(1:nZ)                     , &  ! Output
                  iVar%IRSSEM_BRDF_V2)
           
       END IF
       
       
       DO I = 1, NZ
        SfcOptics%Emissivity(i,1)   = Emissivity(i)
        SfcOptics%Reflectivity(i,1,i,1) = Reflectivity(i)
        SfcOptics%Direct_Reflectivity(i,1) = Direct_Reflectivity(i)
       ENDDO
       
         
       IF ( IO_Status /= SUCCESS ) THEN
            Message = 'Error occurred in forward'
            CALL Display_Message( ROUTINE_NAME, Message, IO_Status );
       END IF
        
    ELSE
       IF((TRIM(MODEL%NAME) .NE. "NESDIS_IRW_WuSmith") .AND. &
          (TRIM(MODEL%NAME) .NE. "NESDIS_IRW_Nalli") ) THEN
         iVar%AlgID = "NESDIS_IRW_WuSmith"
         DB_PATH = TRIM(GET_DATA_PATH('IR_WATER','NESDIS_IRW_WuSmith'))
        ENDIF

       IF(.NOT. NESDIS_IRSSEM_Initialized()) THEN
          IF(TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_Nalli") THEN
            IRSSEM_File = 'Nalli.IRwater.EmisCoeff.nc'
          ELSE
          !"NESDIS_IRW_WuSmith" 
            IRSSEM_File = 'WuSmith.IRwater.EmisCoeff.nc'
          END IF
          IO_Status = NESDIS_IRSSEM_Setup(TRIM(DB_PATH)//'/'//TRIM(IRSSEM_File))
       END IF
        
       IF(SfcOptics%Is_Solar) THEN
          IO_Status = NESDIS_IRSSEM_BRDF(                  &
                  SfcOptics%Wavenumber                   , &  ! Input
                  Surface%Wind_Speed                     , &  ! Input
                  SfcOptics%Angle                        , &  ! Input
                  Sensor_Zenith_Radian                   , &  ! Input
                  Sensor_Azimuth_Radian                  , &  ! Input
                  Solar_Zenith_Radian                    , &  ! Input
                  Solar_Azimuth_Radian                   , &  ! Input
                  Emissivity(1:nZ)                       , &  ! Output
                  Reflectivity(1:nZ)                     , &  ! Output
                  Direct_Reflectivity (1:nZ)             , &  ! Output
                  iVar%IRSSEM_BRDF )                    
        ELSE
          IO_Status = NESDIS_IRSSEM_BRDF( &
                  SfcOptics%Wavenumber                   , &  ! Input
                  Surface%Wind_Speed                     , &  ! Input
                  SfcOptics%Angle                        , &  ! Input
                  Emissivity(1:nZ)                       , &  ! Output
                  Reflectivity(1:nZ)                     , &  ! Output
                  iVar%IRSSEM_BRDF)
           
       END IF
       
       DO I = 1, NZ
        SfcOptics%Emissivity(i,1)   = Emissivity(i)
        SfcOptics%Reflectivity(i,1,i,1) = Reflectivity(i)
        SfcOptics%Direct_Reflectivity(i,1) = Direct_Reflectivity(i)
       ENDDO
       
         
       IF ( IO_Status /= SUCCESS ) THEN
            Message = 'Error occurred in forward'
            CALL Display_Message( ROUTINE_NAME, Message, IO_Status );
       END IF

    END IF

 
  END FUNCTION CSEM_Compute_WaterIR_SfcOptics

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_WaterIR_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at infrared frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterIR_SfcOptics_TL( SfcOptics_TL )
!
! INPUTS:
!
!         Surface_TL:    Structure containing the tagent-linear surface state
!                        variables required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
! OUTPUTS:
!
!        SfcOptics_TL:   Structure containing the tangent-linear surface
!                        optical properties required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Assumed-shape Array
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
!
!:sdoc-:
!----------------------------------------------------------------------------------

 FUNCTION CSEM_Compute_WaterIR_SfcOptics_TL(   &
       Surface_TL,                            &  ! Input
       SfcOptics_TL,                          &  ! Output
       iVar)                                  &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Water_Surface),  INTENT(IN)     :: Surface_TL
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT)  :: SfcOptics_TL
    TYPE(iVar_type),           INTENT(IN)     :: iVar
 
   ! Function result
    INTEGER :: IO_Status 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_WaterIR_SfcOptics_TL'
    ! Local variables
    INTEGER :: i,nZ

    REAL(fp) :: Emissivity_TL(MAX_n_Angles), Reflectivity_TL(MAX_n_Angles)
    REAL(fp) :: Direct_Reflectivity_TL(MAX_n_Angles)
    
    ! Set up
    IO_Status  = SUCCESS

    nZ = SfcOptics_TL%n_Angles
 
    Emissivity_TL          =  ZERO
    Reflectivity_TL        =  ZERO
    Direct_Reflectivity_TL =  ZERO

    IF(TRIM(ivar%AlgID) .EQ. 'NPOESS_LUT' ) THEN


    ELSE IF((TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2CM") .OR. &
          (TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2EK") ) THEN
        IO_Status = NESDIS_IRSSEM_BRDF_V2_TL(            &
                        Surface_TL%Wind_Speed          , &  ! Input
                        Surface_TL%Water_Temperature   , &  ! Input
                        Emissivity_TL(1:nZ)            , &  ! Output
                        Reflectivity_TL(1:nZ)          , &  ! Output
                        Direct_Reflectivity_TL(1:nZ)   , &  ! output
                        iVar%IRSSEM_BRDF_V2        )  ! Internal variable input
        IF ( IO_Status /= SUCCESS ) THEN
           CALL Display_Message( ROUTINE_NAME, &
                           'Error computing Tangent_linear IR sea surface emissivity', &
                            IO_Status  )
           RETURN
        END IF
 

    ELSE
        !Compute tangent-linear IR sea surface emissivity
        IO_Status = NESDIS_IRSSEM_BRDF_TL(               &
                        Surface_TL%Wind_Speed          , &  ! Input
                        Emissivity_TL(1:nZ)            , &  ! Output
                        Reflectivity_TL(1:nZ)          , &  ! Output
                        Direct_Reflectivity_TL(1:nZ)   , &  ! output
                        iVar%IRSSEM_BRDF        )  ! Internal variable input
        IF ( IO_Status /= SUCCESS ) THEN
           CALL Display_Message( ROUTINE_NAME, &
                           'Error computing Tangent_linear IR sea surface emissivity', &
                            IO_Status  )
           RETURN
        END IF
 
    END IF
    DO i=1,nZ
      SfcOptics_TL%Emissivity(i,1)          =  Emissivity_TL(i)
      SfcOptics_TL%Reflectivity(i,1,i,1)    =  Reflectivity_TL(i)
      SfcOptics_TL%Direct_Reflectivity(i,1) =  Direct_Reflectivity_TL(i)
    END DO

 
  END FUNCTION CSEM_Compute_WaterIR_SfcOptics_TL


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_WaterIR_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at infrared frequencies over a Water surface.
!
!       This function is a wrapper for third party code.
!
!       NB: CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO AD
!           COMPONENTS IN THE IR Water SFCOPTICS COMPUTATIONS.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_WaterIR_SfcOptics_AD( SfcOptics_AD )
!
! INPUTS:
!       SfcOptics_AD:  Structure containing the adjoint surface optical
!                      properties required for the adjoint radiative
!                      transfer calculation.
!                      *** COMPONENTS MODIFIED UPON OUTPUT ***
!                      UNITS:      N/A
!                      TYPE:       CSEM_SfcOptics_Type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
! OUTPUTS:
!
!         Surface_AD:  Structure containing the adjoint surface state
!                      variables required for the adjoint radiative
!                      transfer calculation.
!                      *** COMPONENTS MODIFIED UPON OUTPUT ***
!                      UNITS:      N/A
!                      TYPE:       CSEM_SfcOptics_Type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       IO_Status:  The return value is an integer defining the error status.
!                      The error codes are defined in the Message_Handler module.
!                      If == SUCCESS the computation was sucessful
!                      == FAILURE an unrecoverable error occurred
!                      UNITS:      N/A
!                      TYPE:       INTEGER
!                      DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the input adjoint arguments are IN OUT regardless
!       of their specification as "input" or "output". This is because these
!       arguments may contain information on input, or need to be zeroed on
!       output (or both).
!
!:sdoc-:
!----------------------------------------------------------------------------------


  FUNCTION CSEM_Compute_WaterIR_SfcOptics_AD(  &
       SfcOptics_AD,                           &  ! Input
       Surface_AD,                             &  ! Output
       iVar)                                   &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Water_Surface),  INTENT(INOUT) :: Surface_AD
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT) :: SfcOptics_AD
    TYPE(iVar_type) :: iVar
 
   ! Function result
    INTEGER :: IO_Status 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_WaterIR_SfcOptics_AD'
    ! Local variables
    INTEGER :: i, nZ

    REAL(fp) :: Emissivity_AD(MAX_n_Angles), Reflectivity_AD(MAX_n_Angles)
    REAL(fp) :: Direct_Reflectivity_AD(MAX_n_Angles)


    ! Set up
    IO_Status = SUCCESS

    nZ = SfcOptics_AD%n_Angles
 
    !ALLOCATE(Emissivity_AD(nZ), Reflectivity_AD(nZ),Direct_Reflectivity_AD(nZ))
    DO i=1,nZ
      Emissivity_AD(i)          = SfcOptics_AD%Emissivity(i,1)
      Reflectivity_AD(i)        = SfcOptics_AD%Reflectivity(i,1,i,1)
      Direct_Reflectivity_AD(i) = SfcOptics_AD%Direct_Reflectivity(i,1)
    ENDDO

    IF(TRIM(iVar%AlgID) .EQ. 'NPOESS_LUT' ) THEN

       Surface_AD%Wind_Speed         =  Surface_AD%Wind_Speed + ZERO
       Surface_AD%Water_Temperature  = Surface_AD%Water_Temperature + ZERO

    ELSE IF((TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2CM") .OR. &
          (TRIM(iVar%AlgID) .EQ. "NESDIS_IRW_NalliV2EK") ) THEN
        IO_Status = NESDIS_IRSSEM_BRDF_V2_AD(            &
                       Emissivity_AD(1:nz)             , &  ! Input
                       Reflectivity_AD (1:nz)          , &  ! Input
                       Direct_Reflectivity_AD (1:nz)   , &  ! Input
                       Surface_AD%Wind_Speed           , &  ! Output
                       Surface_AD%Water_Temperature    , &  ! Outout
                       iVar%IRSSEM_BRDF_V2)                   
        IF ( IO_Status /= SUCCESS ) THEN
           CALL Display_Message( ROUTINE_NAME, &
                           'Error computing Tangent_linear IR sea surface emissivity', &
                            IO_Status  )
            Surface_AD%Wind_Speed = ZERO
            RETURN
        END IF

    ELSE
        !Compute tangent-linear IR sea surface emissivity
        IO_Status = NESDIS_IRSSEM_BRDF_AD(               &
                       Emissivity_AD(1:nz)             , &  ! Input
                       Reflectivity_AD (1:nz)          , &  ! Input
                       Direct_Reflectivity_AD (1:nz)   , &  ! Input
                       Surface_AD%Wind_Speed           , &  ! Output
                       iVar%IRSSEM_BRDF)                   
        IF ( IO_Status /= SUCCESS ) THEN
           CALL Display_Message( ROUTINE_NAME, &
                           'Error computing Tangent_linear IR sea surface emissivity', &
                            IO_Status  )
            Surface_AD%Wind_Speed = ZERO
            RETURN
        END IF

    END IF
    SfcOptics_AD%Emissivity    =  ZERO
    SfcOptics_AD%Reflectivity  =  ZERO
    SfcOptics_AD%Direct_Reflectivity =  ZERO


  END FUNCTION CSEM_Compute_WaterIR_SfcOptics_AD

END MODULE CSEM_WaterIR_SfcOptics
