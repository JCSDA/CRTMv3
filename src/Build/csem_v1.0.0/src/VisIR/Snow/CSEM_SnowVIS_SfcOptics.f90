!----------------------------------------------------------------------------------
!:sdoc+:
!
! CSEM_SnowVIS_SfcOptics
!
! Module to compute the surface optical properties of Snow surfaces at VISIBLE 
! frequencies, which is required for determining the land surface contribution
! to the overall radiative transfer. 
!
! This module provide a general container for developers to integrate individual 
! models into the model repository, and a generic interface for 
! the upper-level applications to access the model options in the repository.
! 
! This module implementes the FWD, TL and AD functions for variational data 
! assimilation application and retrieval applications. It replacess the similar
! functions as those in the original CRTM surface module "CRTM_SnowVIS_SfcOptics", 
! which was Written by Paul van Delst, 23-Jun-2005
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 28-Feb-2016
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!

MODULE CSEM_SnowVIS_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler
  USE CSEM_Define
  USE CSEM_Model_Manager

  USE NESDIS_SnowVIS_PhyModel
  USE NPOESS_LUT_Module

  !Disable implicit typing
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  ! Science routines
  PUBLIC :: CSEM_Compute_SnowVIS_SfcOptics
  PUBLIC :: CSEM_Compute_SnowVIS_SfcOptics_TL
  PUBLIC :: CSEM_Compute_SnowVIS_SfcOptics_AD
  PUBLIC :: iVar_type
  
  REAL(fp), PARAMETER ::  ONE = 1.0_fp, ZERO = 0.0_fp
  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    INTEGER :: Dummy = 0
  END TYPE iVar_type

  TYPE(CSEM_MODEL_ID),SAVE  :: MODEL

CONTAINS

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_SnowVIS_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at visible
!       frequencies over a Snow surface.
!
!       This function is a wrapper for third party code.
!
! PURPOSE:
!       Function to compute the surface emissivity  at visible
!       frequencies over a Snow surface.
!
!       This function encapsulates all available  CSEM SnowVIS emissivity models,and
!       provides a genereic interface for the upper-level user application.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_SnowVIS_SfcOptics( 
!                                        Surface                  ,&  ! Input
!                                        SfcOptics                ,&  ! Output
!                                        Options )                       ! OptionalInput 
!
! INPUTS:
!
!       Surface:         CSEM Snow Surface structure containing the surface state
!                        inputs.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Snow_Surface
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
!       IO_Status:       The return value is an integer defining the error status.
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

  FUNCTION CSEM_Compute_SnowVIS_SfcOptics(  &
      Surface,                              &  ! Input
      SfcOptics,                            &  ! Output
      Options)                              &  ! Optional Input
    RESULT (IO_Status)

    TYPE(CSEM_Snow_Surface),   INTENT(IN)    :: Surface
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT) :: SfcOptics
    ! Optional input    
    TYPE(CSEM_Options_Type), OPTIONAL, INTENT(IN)  :: Options
    ! Function result
    INTEGER ::  IO_Status

    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_SnowVIS_SfcOptics'
    ! Local variables
    CHARACTER(LEN=256) :: Message
    INTEGER :: Stype 
    REAL(fp):: Emissivity
    INTEGER :: i, n_Angles
    REAL(fp) :: rdummy

    IF(.FALSE.)rdummy =  SfcOptics%Sensor_Zenith_Angle
   
    IF(PRESENT(Options)) THEN
    
    ENDIF 
 
    IO_Status = SUCCESS
    n_Angles = SfcOptics%n_Angles

    MODEL = Inq_Model_Option("VIS_SNOW")        
    SELECT CASE (TRIM(MODEL%NAME))
     
      CASE ("ATLAS")
      
      CASE DEFAULT
        IF(Surface%Snow_Type > 2) THEN
           Message = 'NPOESS_LUT  Snow Type 1-2'
           CALL Display_Message( ROUTINE_NAME, Message, IO_Status )
        ENDIF
        Stype =  Surface%Snow_Type + 21
        IO_Status = NPOESS_LUT_Emiss(        &
            SfcOptics%Wavenumber,            & ! INPUT, wavelength in micrometer
            Emissivity,                      & ! OUTPUT, surface emissivity (0 - 1)
            Stype)                             ! INPUT, surface type (1 - 2)

        DO i = 1, n_Angles
          SfcOptics%Emissivity(i,1)           =  Emissivity
          SfcOptics%Direct_Reflectivity(i,1)  =  ONE-Emissivity
          SfcOptics%Reflectivity(1:n_Angles,1,i,1)     =  (ONE-Emissivity)*SfcOptics%Weight(i)
        END DO

    END SELECT


  END FUNCTION CSEM_Compute_SnowVIS_SfcOptics

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_SnowVIS_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at visible frequencies over a Snow surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_SnowVIS_SfcOptics_TL( SfcOptics_TL )
!
! OUTPUTS:
!       SfcOptics_TL:    Structure containing the tangent-linear surface
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
!       CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO TL
!        COMPONENTS IN THE VIS Snow SFCOPTICS COMPUTATIONS.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CSEM_Compute_SnowVIS_SfcOptics_TL(  &
       SfcOptics_TL)                           &  ! TL  Output
    RESULT ( IO_Status )
    
    ! Arguments
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT)   :: SfcOptics_TL
    ! Function result
    INTEGER :: IO_Status
    
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_SnowVIS_SfcOptics_TL'


    ! Set up
    IO_Status = SUCCESS

    ! Compute the tangent-linear surface optical parameters
    ! ***No TL models yet, so default TL output is zero***
    SfcOptics_TL%Emissivity          =  ZERO
    SfcOptics_TL%Reflectivity        =  ZERO
    SfcOptics_TL%Direct_Reflectivity =  ZERO


  END FUNCTION CSEM_Compute_SnowVIS_SfcOptics_TL



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_SnowVIS_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at visible frequencies over a Snow surface.
!
!       This function is a wrapper for third party code.
!
!       NB: CURRENTLY THIS IS A STUB FUNCTION AS THERE ARE NO AD
!           COMPONENTS IN THE VIS Snow SFCOPTICS COMPUTATIONS.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_SnowVIS_SfcOptics_AD( SfcOptics_AD )
!
! INPUTS:
!
!      SfcOptics_AD:    Structure containing the adjoint surface optical
!                        properties required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! OUTPUTS:
!      Surface_AD:       Structure containing the adjoint surface state
!                        variables required for the adjoint radiative
!                        transfer calculation.
!                        *** COMPONENTS MODIFIED UPON OUTPUT ***
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
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

  FUNCTION CSEM_Compute_SnowVIS_SfcOptics_AD(  &
         Surface_AD)                           &  ! AD  Input
    RESULT( IO_Status )
    
    ! Arguments
    TYPE(CSEM_Snow_Surface), INTENT(INOUT)   :: Surface_AD
    ! Function result
    INTEGER :: IO_Status

    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_SnowVIS_SfcOptics_AD'

    ! Set up
    IO_Status = SUCCESS

    ! Compute the adjoint surface optical parameters
    ! ***No AD models yet, so there is no impact on AD result***
    Surface_AD%Snow_Temperature           = Surface_AD%Snow_Temperature     + ZERO
    Surface_AD%Snow_Depth                 = Surface_AD%Snow_Depth           + ZERO
    Surface_AD%Snow_Density               = Surface_AD%Snow_Density         + ZERO
    Surface_AD%Snow_Grain_Size            = Surface_AD%Snow_Grain_Size      + ZERO
    Surface_AD%Top_Soil_Temperature       = Surface_AD%Top_Soil_Temperature + ZERO
    Surface_AD%Top_Soil_Moisture_Content  = Surface_AD%Top_Soil_Moisture_Content + ZERO

  END FUNCTION CSEM_Compute_SnowVIS_SfcOptics_AD

END MODULE CSEM_SnowVIS_SfcOptics
