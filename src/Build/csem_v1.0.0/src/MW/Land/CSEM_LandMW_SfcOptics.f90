!----------------------------------------------------------------------------------
!:sdoc+:
!
! CSEM_LandMW_SfcOptics
!
! Module to compute the surface optical properties of Land surfaces at microwave
! frequencies, which is required for determining the surface contribution
! to the overall radiative transfer. 
!
! This module provide a general template for developers to integrate individual 
! models into the model repository, and a generic interface for 
! the upper-level applications to access the model options in the model repository.
! 
! This module implementes the FWD, TL and AD functions for variational data 
! assimilation application and retrieval applications. It replacess the similar
! functions as those in the original CRTM surface module "CRTM_LandMW_SfcOptics", 
! which was Written by Paul van Delst, 23-Jun-2005
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 28-Feb-2015
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!

MODULE CSEM_LandMW_SfcOptics

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Define
  USE CSEM_Model_Manager
  USE CSEM_Exception_Handler
 
  USE NESDIS_LandMW_PhyModel, ONLY: phy_iVar =>  iVar_type, &
                                    NESDIS_LandMW_Emiss,    &
                                    NESDIS_LandMW_Emiss_TL, &
                                    NESDIS_LandMW_Emiss_AD
  USE NESDIS_LandEM_Module,   ONLY: NESDIS_LandEM_213 =>  NESDIS_LandEM_213, &
                                    NESDIS_LandEM_205 =>  NESDIS_LandEM_OLD
  USE TELSEM_Atlas_Module
  USE TELSEM2_Atlas_Module
  USE CNRM_Atlas_Module

  !Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Data types
  ! Science routines
  PUBLIC :: CSEM_Compute_LandMW_SfcOptics
  PUBLIC :: CSEM_Compute_LandMW_SfcOptics_TL
  PUBLIC :: CSEM_Compute_LandMW_SfcOptics_AD
  PUBLIC :: iVar_type
  PUBLIC :: Get_Ref_Index
  
  REAL(fp), PARAMETER  ::  ONE = 1.0_fp, ZERO = 0.0_fp
  INTEGER,  PARAMETER  ::  MAX_N_ANGLES = 16
  ! Valid type indices for the microwave land emissivity model
  ! ...The soil types
  INTEGER, PARAMETER :: N_VALID_SOIL_TYPES = 8
  INTEGER, PARAMETER :: INVALID_SOIL    =  0
  INTEGER, PARAMETER :: COARSE          =  1
  INTEGER, PARAMETER :: MEDIUM          =  2
  INTEGER, PARAMETER :: FINE            =  3
  INTEGER, PARAMETER :: COARSE_MEDIUM   =  4
  INTEGER, PARAMETER :: COARSE_FINE     =  5
  INTEGER, PARAMETER :: MEDIUM_FINE     =  6
  INTEGER, PARAMETER :: COARSE_MED_FINE =  7
  INTEGER, PARAMETER :: ORGANIC         =  8
  ! ...The vegetation types
  INTEGER, PARAMETER :: N_VALID_VEGETATION_TYPES       = 12
  INTEGER, PARAMETER :: INVALID_VEGETATION             =  0
  INTEGER, PARAMETER :: BROADLEAF_EVERGREEN_TREES      =  1
  INTEGER, PARAMETER :: BROADLEAF_DECIDUOUS_TREES      =  2
  INTEGER, PARAMETER :: BROADLEAF_NEEDLELEAF_TREES     =  3
  INTEGER, PARAMETER :: NEEDLELEAF_EVERGREEN_TREES     =  4
  INTEGER, PARAMETER :: NEEDLELEAF_DECIDUOUS_TREES     =  5
  INTEGER, PARAMETER :: BROADLEAF_TREES_GROUNDCOVER    =  6
  INTEGER, PARAMETER :: GROUNDCOVER                    =  7
  INTEGER, PARAMETER :: GROADLEAF_SHRUBS_GROUNDCOVER   =  8
  INTEGER, PARAMETER :: BROADLEAF_SHRUBS_BARE_SOIL     =  9
  INTEGER, PARAMETER :: DWARF_TREES_SHRUBS_GROUNDCOVER = 10
  INTEGER, PARAMETER :: BARE_SOIL                      = 11
  INTEGER, PARAMETER :: CULTIVATIONS                   = 12
  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    TYPE(phy_iVar) :: phyVar(MAX_N_ANGLES)
    LOGICAL :: FREQ_CUTOFF = .FALSE.
  END TYPE iVar_type

  TYPE(CSEM_Model_ID), SAVE  :: MODEL

CONTAINS

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_LandMW_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a Land surface.
!
!       This function is a wrapper for third party code.
!
! PURPOSE:
!       Function to compute the surface emissivity  at microwave
!       frequencies over a Land surface.
!
!       This function encapsulates all available  CSEM LandMW emissivity models,and
!       provides a genereic interface for the upper-level user application.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_LandMW_SfcOptics( 
!                                        Surface                  ,&  ! Input
!                                        SfcOptics                ,&  ! Output
!                                        Options                  ,&  ! OptionalInput 
!                                        iVar)                     &  ! output
!
! INPUTS:
!
!       Surface:         CSEM Land Surface structure containing the surface state
!                        inputs.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Land_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Options:         CSEM input extension structure containing additional inputs,
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
!                        outside of the CRTM_LandMW_SfcOptics module.
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
!                       CSEM infrastructure data structures are defined in CSEM_Define.f90
!       
!
!:sdoc-:
!-----------------------------------------------------------------------------------------------


  FUNCTION CSEM_Compute_LandMW_SfcOptics( &
      Surface,                  &  ! Input
      SfcOptics,                &  ! Output
      Options,                  &  ! OptionalInput 
      iVar)                     &  ! Output
    RESULT (IO_Status)
    
    TYPE(CSEM_Land_Surface),   INTENT(IN)     :: Surface
    TYPE(CSEM_Options_Type),   INTENT(IN)     :: Options
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT)  :: SfcOptics
    TYPE(iVar_type),           INTENT(OUT)    :: iVar

    ! Function result
    INTEGER :: IO_Status 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_LandMW_SfcOptics'
    REAL(fp),     PARAMETER :: FREQUENCY_CUTOFF   = 200.0_fp  ! GHz
    REAL(fp),     PARAMETER :: DEFAULT_EMISSIVITY = 0.95_fp
    CHARACTER(LEN=256):: DB_PATH
    ! Local variables
    INTEGER :: i
    INTEGER :: imonth
    INTEGER :: N_Angles
    CHARACTER(256) :: msg

    IO_Status = SUCCESS
    N_Angles = SfcOptics%N_Angles
    SfcOptics%Reflectivity = ZERO
    SfcOptics%Emissivity   = ZERO
    
    MODEL = Inq_Model_Option("MW_LAND") 
    ! use telsem_atlas_path in CSEM alg registor file
    DB_PATH = TRIM(MODEL%DATA_PATH)
    IF(TRIM(MODEL%NAME) .EQ. 'TELSEM_ATLAS' ) THEN

        imonth = Options%GeoInfo%Month
        IF (.NOT. (TELSEM_Atlas_Initialized(imonth))) THEN
          IO_Status = TELSEM_Atlas_Setup(imonth, path=TRIM(DB_PATH))
        ENDIF
        DO i = 1, N_Angles
          IO_Status = TELSEM_Atlas_Emiss(        &
                SfcOptics%Frequency,                &  ! Input
                SfcOptics%Angle(i),                 &  ! Input
                Options%GeoInfo%Latitude,           &  ! Input
                Options%GeoInfo%Longitude,          &  ! Input
                Options%GeoInfo%Month,              &  ! Input
                SfcOptics%Emissivity(i,2),          &  ! Output
                SfcOptics%Emissivity(i,1))             ! Output
                SfcOptics%Emissivity(i,3:4) = ZERO
        END DO

    ELSE IF(TRIM(MODEL%NAME) .EQ. 'TELSEM2_ATLAS' ) THEN

        imonth = Options%GeoInfo%Month
        IF (.NOT. (TELSEM_Atlas_Initialized(imonth))) THEN
          IO_Status = TELSEM2_Atlas_Setup(imonth, path=TRIM(DB_PATH))
        ENDIF
        DO i = 1, N_Angles
          IO_Status = TELSEM2_Atlas_Emiss(       &
                SfcOptics%Frequency,                &  ! Input
                SfcOptics%Angle(i),                 &  ! Input
                Options%GeoInfo%Latitude,           &  ! Input
                Options%GeoInfo%Longitude,          &  ! Input
                Options%GeoInfo%Month,              &  ! Input
                SfcOptics%Emissivity(i,2),          &  ! Output
                SfcOptics%Emissivity(i,1))             ! Output
                SfcOptics%Emissivity(i,3:4) = ZERO
        END DO
    ELSE IF(TRIM(MODEL%NAME) .EQ. 'CNRM_AMSUA_ATLAS' .AND. &
       INDEX(TRIM(Options%SensorObs%Sensor_Id),'amsu') > 0) THEN

        imonth = Options%GeoInfo%Month
        IF (.NOT. (CNRM_Atlas_Initialized(imonth))) THEN
          IO_Status = CNRM_Atlas_Setup(imonth, path=TRIM(DB_PATH))
        ENDIF

        DO i = 1, N_Angles
          IO_Status = CNRM_Atlas_Emiss(         &
                SfcOptics%Frequency,               &  ! Input
                SfcOptics%Angle(i),                &  ! Input
                Options%GeoInfo%Latitude,          &  ! Input
                Options%GeoInfo%Longitude,         &  ! Input
                Options%GeoInfo%Month,             &  ! Input
                SfcOptics%Emissivity(i,2),         &  ! Output
                SfcOptics%Emissivity(i,1))            ! Output
          SfcOptics%Emissivity(i,3:4) = ZERO
       END DO

    ELSE IF(TRIM(MODEL%NAME) .EQ. 'NESDIS_Land_MW213' ) THEN
       IF ( Surface%Soil_Type < 1 .OR. &
         Surface%Soil_Type > N_VALID_SOIL_TYPES ) THEN
         SfcOptics%Emissivity   = ZERO
         SfcOptics%Reflectivity = ZERO
         IO_Status = FAILURE
         msg = 'Invalid soil type index specified'
         WRITE(*,*) TRIM(msg);  RETURN
       END IF
       ! ...and the vegetation type
       IF ( Surface%Vegetation_Type < 1 .OR. &
         Surface%Vegetation_Type > N_VALID_VEGETATION_TYPES ) THEN
         SfcOptics%Emissivity   = ZERO
         SfcOptics%Reflectivity = ZERO
         IO_Status = FAILURE
         msg = 'Invalid vegetation type index specified'
         WRITE(*,*) TRIM(msg);  RETURN
       END IF

       IF(SfcOptics%Frequency < 80.0_fp) THEN
 
         DO i = 1, N_Angles
           CALL  NESDIS_LandEM_213(                &  ! Input
                 SfcOptics%Angle(i),               &  ! Input
                 SfcOptics%Frequency,              &  ! Input
                 Surface%Top_Soil_Moisture,        &  ! Input
                 Surface%Vegetation_Fraction,      &  ! Input
                 Surface%Top_Soil_Temperature,     &  ! Input
                 Surface%Land_Skin_Temperature,    &  ! Input
                 Surface%LAI,                      &  ! Input
                 Surface%Soil_Type,                &  ! Input
                 Surface%Vegetation_Type,          &  ! Input
                 SfcOptics%Emissivity(i,2),        &  ! Output
                 SfcOptics%Emissivity(i,1))            ! Output

           ! Assume specular surface
           SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
           SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
         END DO
       ELSE
         ! Frequency is too high for model. Use default.
         DO i = 1, SfcOptics%n_Angles
            SfcOptics%Emissivity(i,1:2)         = DEFAULT_EMISSIVITY
            SfcOptics%Reflectivity(i,1:2,i,1:2) = ONE-DEFAULT_EMISSIVITY
         END DO
       END IF

    ELSE IF(TRIM(MODEL%NAME) .EQ. 'NESDIS_Land_MW' ) THEN
       IF(SfcOptics%Frequency < 80.0_fp) THEN
 
         DO i = 1, N_Angles
           IO_Status = NESDIS_LandMW_Emiss(    &
                SfcOptics%Frequency,              &  ! Input
                SfcOptics%Angle(i),               &  ! Input
                Surface%Land_Skin_Temperature,    &  ! Input
                Surface%Top_Soil_Temperature,     &  ! Input
                Surface%Top_Soil_Moisture,        &  ! Input
                Surface%Vegetation_Fraction,      &  ! Input
                Surface%LAI,                      &  ! Input
                Surface%Vegetation_Type,          &  ! Input
                Surface%Soil_Type,                &  ! Input
                SfcOptics%Emissivity(i,2),        &  ! Output
                SfcOptics%Emissivity(i,1),        &  ! Output
                iVar%phyVar(i))                      ! Output
            ! Assume specular surface
            SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
            SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
        END DO
       ELSE
         ! Frequency is too high for model. Use default.
         ivar%FREQ_CUTOFF = .TRUE.
         DO i = 1, SfcOptics%n_Angles
            SfcOptics%Emissivity(i,1:2)         = DEFAULT_EMISSIVITY
            SfcOptics%Reflectivity(i,1:2,i,1:2) = ONE-DEFAULT_EMISSIVITY
         END DO

       ENDIF
    ELSE 
       IO_Status = FAILURE
        msg = 'Invalid MW_LAND model name'
       WRITE(*,*) TRIM(msg);  RETURN

    END IF

    IO_Status = SUCCESS


  END FUNCTION CSEM_Compute_LandMW_SfcOptics

  SUBROUTINE Get_Ref_Index(Frequency, Polarization, i_ref_h,i_ref_v)
 
    REAL(fp):: Frequency(:)
    INTEGER :: Polarization(:)
    INTEGER :: i, i_ref_h,i_ref_v
    !INTEGER, PARAMETER :: N_POLARIZATION_TYPES    = 12
    !INTEGER, PARAMETER :: INVALID_POLARIZATION    = 0
    !INTEGER, PARAMETER :: UNPOLARIZED             = 1
    !INTEGER, PARAMETER :: INTENSITY               = UNPOLARIZED
    !INTEGER, PARAMETER :: FIRST_STOKES_COMPONENT  = UNPOLARIZED
    !INTEGER, PARAMETER :: SECOND_STOKES_COMPONENT = 2
    !INTEGER, PARAMETER :: THIRD_STOKES_COMPONENT  = 3
    !INTEGER, PARAMETER :: FOURTH_STOKES_COMPONENT = 4
    !INTEGER, PARAMETER :: VL_POLARIZATION         = 5
    !INTEGER, PARAMETER :: HL_POLARIZATION         = 6
    !INTEGER, PARAMETER :: plus45L_POLARIZATION    = 7
    !INTEGER, PARAMETER :: minus45L_POLARIZATION   = 8
    !INTEGER, PARAMETER :: VL_MIXED_POLARIZATION   = 9
    !INTEGER, PARAMETER :: HL_MIXED_POLARIZATION   = 10
    !INTEGER, PARAMETER :: RC_POLARIZATION         = 11
    !INTEGER, PARAMETER :: LC_POLARIZATION         = 12

    i_ref_h = -1 ; i_ref_v = -1
    DO i = 1, size(Frequency)
       IF(Frequency(i) > 15.0_fp .AND. Frequency(i) <= 20.0_fp ) THEN
           IF(Polarization(i) == 5 .OR. Polarization(i) == 9)  i_ref_v = i
           IF(Polarization(i) == 6 .OR. Polarization(i) == 10) i_ref_h = i
       ENDIF 
    ENDDO
      
    IF(i_ref_v > 0 .AND. i_ref_h > 0)  RETURN
    i_ref_h = -1 ; i_ref_v = -1
    DO i = 1, size(Frequency)
       IF(Frequency(i) > 20.0_fp .AND. Frequency(i) <= 25.0_fp ) THEN
          i_ref_v = i
          i_ref_h = i 
          EXIT
       ENDIF
    ENDDO
      
  END SUBROUTINE Get_Ref_Index

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_LandMW_SfcOptics_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at microwave frequencies over a Land surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_LandMW_SfcOptics_TL( 
!                                    Surface_TL,      &
!                                    SfcOptics_TL,    &
!                                    iVar)
!
! INPUTS:
!      Surface_TL:   Structure containing the tangent-linear surface
!                        state inputs required for the tangent-
!                        linear radiative transfer calculation.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Land_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_LandMW_SfcOptics module.
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

  FUNCTION CSEM_Compute_LandMW_SfcOptics_TL( &
      Surface_TL,               &  ! Input
      SfcOptics_TL,             &  ! Output
      iVar        )             &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Land_Surface),   INTENT(IN)     :: Surface_TL
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT)  :: SfcOptics_TL
    TYPE(iVar_type),           INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: IO_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_LandMW_SfcOptics_TL'
    ! Local variables
    REAL(fp) :: Emissivity_H_TL, Emissivity_V_TL
    INTEGER :: i, N_Angles

    IO_Status = SUCCESS
    N_Angles = SfcOptics_TL%N_Angles
   
    SfcOptics_TL%Direct_Reflectivity = ZERO
    SfcOptics_TL%Reflectivity = ZERO
    SfcOptics_TL%Emissivity   = ZERO   

    MODEL = Inq_Model_Option("MW_LAND") 

    IF(TRIM(MODEL%NAME) .EQ. 'NESDIS_Land_MW' .AND. &
       .NOT. ivar%FREQ_CUTOFF) THEN
        DO i = 1, N_Angles
          IO_Status = NESDIS_LandMW_Emiss_TL(        &
                Surface_TL%Land_Skin_Temperature,       &  ! Input
                Surface_TL%Top_Soil_Temperature,        &  ! Input
                Surface_TL%Top_Soil_Moisture,           &  ! Input
                Surface_TL%Vegetation_Fraction,         &  ! Input
                Emissivity_H_TL,                        &  ! Output
                Emissivity_V_TL,                        &  ! Output
                iVar%phyVar(i))                            ! Internal variable output
          SfcOptics_TL%Emissivity(i,1)     = Emissivity_V_TL
          SfcOptics_TL%Emissivity(i,2)     = Emissivity_H_TL
          SfcOptics_TL%Reflectivity(i,1,i,1) = -SfcOptics_TL%Emissivity(i,1)
          SfcOptics_TL%Reflectivity(i,2,i,2) = -SfcOptics_TL%Emissivity(i,2)

        END DO

      ELSE
          SfcOptics_TL%Reflectivity = ZERO
          SfcOptics_TL%Emissivity   = ZERO
    END IF

  END FUNCTION CSEM_Compute_LandMW_SfcOptics_TL



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      CSEM_Compute_LandMW_SfcOptics_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at microwave frequencies over a land surface.
!
!       This function is a wrapper for third party code.
!
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_LandMW_SfcOptics_AD(
!                                  SfcOptics_AD,         &
!                                  Surface_AD,           &
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
!                        outside of the CRTM_LandMW_SfcOptics module.
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
!                        TYPE:       CSEM_Land_Surface
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

  FUNCTION CSEM_Compute_LandMW_SfcOptics_AD( &
      SfcOptics_AD,             &  ! Input
      Surface_AD,               &  ! Output
      iVar        )             &  ! Internal variable input
   RESULT ( IO_Status  )
    ! Arguments
    TYPE(CSEM_Land_Surface),   INTENT(INOUT) :: Surface_AD
    TYPE(CSEM_SfcOptics_Type), INTENT(INOUT) :: SfcOptics_AD
    TYPE(iVar_type),           INTENT(IN)    :: iVar
 
    ! Function result
    INTEGER :: IO_Status 
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_Compute_LandMW_SfcOptics_AD'
    ! Local variables

    ! Set up
    REAL(fp) :: Emissivity_H_AD, Emissivity_V_AD
    INTEGER :: i, N_Angles

    IO_Status = SUCCESS
    N_Angles = SfcOptics_AD%N_Angles

    MODEL = Inq_Model_Option("MW_LAND") 
    IF(TRIM(MODEL%NAME) .EQ. 'NESDIS_Land_MW' .AND. &
        .NOT. ivar%FREQ_CUTOFF) THEN

        DO i = 1, N_Angles
           Emissivity_V_AD = SfcOptics_AD%Emissivity(i,1)     
           Emissivity_H_AD = SfcOptics_AD%Emissivity(i,2)  
           Emissivity_V_AD = Emissivity_V_AD - SfcOptics_AD%Reflectivity(i,1,i,1)   
           Emissivity_H_AD = Emissivity_H_AD - SfcOptics_AD%Reflectivity(i,2,i,2)
           IO_Status = NESDIS_LandMW_Emiss_AD(       &
                Surface_AD%Land_Skin_Temperature,       &  ! Input
                Surface_AD%Top_Soil_Temperature,        &  ! Input
                Surface_AD%Top_Soil_Moisture,           &  ! Input
                Surface_AD%Vegetation_Fraction,         &  ! Input
                Emissivity_H_AD,                        &  ! Output
                Emissivity_V_AD,                        &  ! Output
                iVar%phyVar(i))                            ! Internal variable output

        END DO

      ELSE
        Surface_AD%Land_Skin_Temperature     =  Surface_AD%Land_Skin_Temperature     + ZERO
        Surface_AD%Top_Soil_Temperature      =  Surface_AD%Top_Soil_Temperature      + ZERO
        Surface_AD%Top_Soil_Moisture         =  Surface_AD%Top_Soil_Moisture         + ZERO
        Surface_AD%Vegetation_Fraction       =  Surface_AD%Vegetation_Fraction       + ZERO

    END IF

    SfcOptics_AD%Direct_Reflectivity = ZERO
    SfcOptics_AD%Reflectivity = ZERO
    SfcOptics_AD%Emissivity   = ZERO
 
  END FUNCTION CSEM_Compute_LandMW_SfcOptics_AD

END MODULE CSEM_LandMW_SfcOptics
