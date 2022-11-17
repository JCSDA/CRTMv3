!
! NESDIS_WaterIR_PhyModel_V2
!
! Module containing the NESDIS water emissivity model of infrared channels.
! The phsyicall model computes the surface optical properties for WATER surfaces at
! infrared frequencies required for determining the WATER surface
! contribution to the radiative transfer.
!
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, ESSIC_UMD, 25-May-2020
!                       ming.chen@noaa.gov
!

MODULE NESDIS_WaterIR_PhyModel_V2

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler,   ONLY: SUCCESS, Display_Message
  USE NESDIS_WaterIR_Emiss_V2_Module,   ONLY: IRSSEM_V2_type=>Einterp_V2_type, &
                                      IRSSEM_V2_Setup,                &
                                      IRSSEM_V2_CleanUp,              &
                                      IRSSEM_V2_Initialized,          &
                                      NESDIS_WaterIR_Emiss_V2,        &
                                      NESDIS_WaterIR_Emiss_V2_TL,     &
                                      NESDIS_WaterIR_Emiss_V2_AD
  USE NESDIS_WaterIR_BRDF_Module,     ONLY: BRDF_type=>iVar_type,     &
                                      NESDIS_IRWater_BRDF,            &
                                      NESDIS_IRWater_BRDF_TL,         &
                                      NESDIS_IRWater_BRDF_AD

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
  PUBLIC :: NESDIS_IRSSEM_BRDF_V2
  PUBLIC :: NESDIS_IRSSEM_BRDF_V2_TL
  PUBLIC :: NESDIS_IRSSEM_BRDF_V2_AD
  PUBLIC :: NESDIS_IRSSEM_V2_Setup
  PUBLIC :: NESDIS_IRSSEM_V2_Close
  PUBLIC :: NESDIS_IRSSEM_V2_Initialized
  
  INTERFACE NESDIS_IRSSEM_BRDF_V2
    MODULE PROCEDURE Compute_WaterIR_Emiss_V2
    MODULE PROCEDURE Compute_WaterIR_Emiss_BRDF_V2
  END INTERFACE NESDIS_IRSSEM_BRDF_V2

  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: CRTM_IR_Water_SfcOptics.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Coefficients for Sigma**2 in the Cox & Munk slope probability density function
  REAL(fp), PARAMETER :: CM_1 = 0.003_fp, CM_2 = 5.12e-3_fp
  REAL(fp), PARAMETER :: ZERO = 0.0_fp, ONE=1.0_fp, TWO = 2.0_fp
  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643_fp


  ! --------------------------------------
  ! Structure definition to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Variables in routines rough sea BRDF
    TYPE(BRDF_type)   :: BRDF
    ! IRSSEM data structure
    TYPE(IRSSEM_V2_type) :: IRSSEM
  END TYPE iVar_type


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_WaterIR_Emiss_V2
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity of one infrared
!       channel at several view angles over a water surface without solar contribution.
!
!
! CALLING SEQUENCE:
!       Error_Status = Compute_WaterIR_Emiss_V2( &
!                        Wavenumber           , &  ! Input
!                        Wind_Speed           , &  ! Input
!                        Temperature          , &  ! Input
!                        Angle                , &  ! Input
!                        Emissivity           , &  ! Output
!                        Reflectivity         , &  ! Output
!                        iVar        )          &  ! Internal variable output
!
! INPUTS:
!       Wavenumber:      wavenumber
!                        data.
!                        UNITS:      cm-1
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Wind_speed:      Ocean surface speed
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Temperature:     Ocean surface temperature
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!       Angle:           Receving angle
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(IN)
!
!
! OUTPUTS:
!       Emissivity:      Surface non-polarized emissivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!       Reflectivity:    Surface non-polarized reflectivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_IR_Water_SfcOptics module.
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
  FUNCTION Compute_WaterIR_Emiss_V2( &
    Wavenumber                    , &  ! Input
    Wind_Speed                    , &  ! Input
    Temperature                   , &  ! Input
    Angle                         , &  ! Input
    Emissivity                    , &  ! Output
    Reflectivity                  , &  ! Output
    iVar)                           &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    REAL(fp), INTENT(IN)     :: Wavenumber
    REAL(fp), INTENT(IN)     :: Wind_Speed
    REAL(fp), INTENT(IN)     :: Temperature
    REAL(fp), INTENT(IN)     :: Angle(:)
    REAL(fp), INTENT(OUT)    :: Emissivity(:)
    REAL(fp), INTENT(OUT)    :: Reflectivity(:)
    TYPE(iVar_type), INTENT(INOUT)  :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_WaterIR_Emiss'
    ! Local variables
    INTEGER  :: nZ
 

    ! Set up
    Error_Status = SUCCESS
    ! ...Short name for angle dimensions
    nZ = size(Angle)
    ! Compute IR sea surface emissivity
    Error_Status = NESDIS_WaterIR_Emiss_V2(Wind_Speed,       &
                                        Temperature,      &
                                        Wavenumber,       &
                                        Angle(1:nZ),      &
                                        Emissivity(1:nZ), &
                                        iVar%IRSSEM )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error computing IR sea surface emissivity', &
                            Error_Status )
      RETURN
    END IF
    
    iVar%BRDF%isSolar  = .FALSE.
    ! Surface reflectance (currently assumed to be specular ALWAYS)
    Reflectivity(1:nZ) = ONE-Emissivity(1:nZ)
    
  
  END FUNCTION Compute_WaterIR_Emiss_V2
  
  
!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       Compute_WaterIR_Emiss_BRDF_V2
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity of one infrared
!       channel at several view angles over a water surface WITH solar contribution.
!
!
! CALLING SEQUENCE:
!       Error_Status = Compute_WaterIR_Emiss_BRD_V2F( &
!                        Wavenumber             , &  ! Input
!                        Wind_Speed             , &  ! Input
!                        Angle                  , &  ! Input
!                        Sensor_Zenith_Radian   , &  ! Input
!                        Sensor_Azimuth_Radian  , &  ! Input
!                        Source_Zenith_Radian   , &  ! Input
!                        Source_Azimuth_Radian  , &  ! Input
!                        Emissivity             , &  ! Output
!                        Reflectivity           , &  ! Output
!                        Direct_Reflectivity    , &  ! Output
!                        iVar)                    &  ! Internal variable output
!
! INPUTS:
!          Wavenumber:   wavenumber
!                        data.
!                        UNITS:      cm-1
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!          Wind_speed:   Ocean surface speed
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!       Temperature:     Ocean surface temperature
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!               Angle:   Receving angle
!                        UNITS:      Degree
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(IN)
!
!  Sensor_Zenith_Radian: Sensor zenith angle 
!                        UNITS:      Radiian
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! Sensor_Azimuth_Radian: Sensor azimuth angle 
!                        UNITS:      Radian
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!  Source_Zenith_Radian: Solar zenith angle 
!                        UNITS:      Radiian
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! Source_Azimuth_Radian: Solar azimuth angle 
!                        UNITS:      Radian
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!        Emissivity:     Surface non-polarized emissivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!       Reflectivity:    Surface non-polarized reflectivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!  Direct_Reflectivity:  Surface solar reflectivity
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_IR_Water_SfcOptics module.
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

  FUNCTION Compute_WaterIR_Emiss_BRDF_V2( &
      Wavenumber                    , &  ! Input
      Wind_Speed                    , &  ! Input
      Temperature                   , &  ! Input
      Angle                         , &  ! Input
      Sensor_Zenith_Radian          , &  ! Input
      Sensor_Azimuth_Radian         , &  ! Input
      Source_Zenith_Radian          , &  ! Input
      Source_Azimuth_Radian         , &  ! Input
      Emissivity                    , &  ! Output
      Reflectivity                  , &  ! Output
      Direct_Reflectivity           , &  ! Output
      iVar)                           &  ! Internal variable output
  RESULT( Error_Status )
    ! Arguments
    REAL(fp),  INTENT(IN)    ::   Wavenumber
    REAL(fp),  INTENT(IN)    ::   Wind_Speed
    REAL(fp),  INTENT(IN)    ::   Temperature
    REAL(fp),  INTENT(IN)    ::   Angle(:)
    REAL(fp),  INTENT(OUT)   ::   Emissivity(:)
    REAL(fp),  INTENT(OUT)   ::   Reflectivity(:)
    REAL(fp),  INTENT(IN)    ::   Sensor_Zenith_Radian
    REAL(fp),  INTENT(IN)    ::   Sensor_Azimuth_Radian
    REAL(fp),  INTENT(IN)    ::   Source_Zenith_Radian
    REAL(fp),  INTENT(IN)    ::   Source_Azimuth_Radian
    REAL(fp),  INTENT(OUT)   ::   Direct_Reflectivity(:)
    TYPE(iVar_type), INTENT(INOUT)  :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Compute_WaterIR_Emiss_BRDF'
    ! Local variables
    INTEGER  :: nZ
    REAL(fp) :: Relative_Azimuth_Radian
    REAL(fp) :: BRDF
    TYPE(BRDF_type)   :: iVar_BRDF


    ! Set up
    Error_Status = SUCCESS
    ! ...Short name for angle dimensions
    nZ = size(Angle)
    ! Compute IR sea surface emissivity
    Error_Status = NESDIS_WaterIR_Emiss_V2(Wind_Speed,       &
                                        Temperature,      &
                                        Wavenumber,       &
                                        Angle(1:nZ),      &
                                        Emissivity(1:nZ), &
                                        iVar%IRSSEM )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error computing IR sea surface emissivity', &
                            Error_Status )
      RETURN
    END IF
    
    
    ! Compute the solar direct BRDF
    Direct_Reflectivity(1:nZ) = ZERO
    IF(Source_Zenith_Radian < PI/TWO ) THEN
       Relative_Azimuth_Radian = Sensor_Azimuth_Radian - &
                            Source_Azimuth_Radian
       Error_Status =  NESDIS_IRWater_BRDF(        &
                            Wavenumber           , &  ! Input
                            Wind_Speed           , &  ! Input
                            Sensor_Zenith_Radian , &  ! Input
                            Sensor_Azimuth_Radian, &  ! Input
                            Source_Zenith_Radian , &  ! Input
                            Source_Azimuth_Radian, &  ! Input
                            BRDF                 , &  ! Output
                            iVar_BRDF            ) ! Internal variable output
       Direct_Reflectivity(1:nZ) = BRDF
       iVar%BRDF = iVar_BRDF 
    END IF
    iVar%BRDF%isSolar  = .TRUE.
     
    ! Surface reflectance (currently assumed to be specular ALWAYS)
    Reflectivity(1:nZ) = ONE-Emissivity(1:nZ)
 
  END FUNCTION Compute_WaterIR_Emiss_BRDF_V2
  
  
!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NESDIS_IRSSEM_BRDF_V2_TL
!
! PURPOSE:
!       Function to compute the tangent-linear surface emissivity and
!       reflectivity at infrared frequencies over a water surface.
!
!
! CALLING SEQUENCE:
!       Error_Status = NESDIS_IRSSEM_BRDF_V2_TL(      &
!                        Wind_Speed_TL           , &  ! Input
!                        Temperature_TL          , &  ! Input
!                        Emissivity_TL           , &  ! Output
!                        Reflectivity_TL         , &  ! Output
!                        Direct_Reflectivity_TL  , &  ! Output
!                        iVar        )             &  ! Internal variable output
!
! INPUTS:
!
!       Wind_speed_TL:   Ocean surface speed tangent-linear
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!       Temperature_TL:  Ocean surface temperature tangent-linear
!                        UNITS:      m/s
!                        TYPE:       REAL
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

! OUTPUTS:
!      Emissivity_TL:    Surface non-polarized emissivity tangent-linear
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!    Reflectivity_TL:    Surface non-polarized reflectivity tangent-linear
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
! Direct_Reflectivity_TL:Surface solar reflectivity tangent-linear
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
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


  FUNCTION NESDIS_IRSSEM_BRDF_V2_TL( &
    Wind_Speed_TL               , &  ! Input
    Temperature_TL              , &  ! Input
    Emissivity_TL               , &  ! Output
    Reflectivity_TL             , &  ! Output
    Direct_Reflectivity_TL      , &  ! output
    iVar                        ) &  ! Internal variable input
  RESULT ( Error_Status )
    REAL(fp),              INTENT(IN)     :: Wind_Speed_TL
    REAL(fp),              INTENT(IN)     :: Temperature_TL
    REAL(fp),              INTENT(OUT)    :: Emissivity_TL(:)
    REAL(fp),              INTENT(OUT)    :: Reflectivity_TL(:)
    REAL(fp),              INTENT(OUT)    :: Direct_Reflectivity_TL(:)
    TYPE(iVar_type),       INTENT(IN)     :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = ' NESDIS_IRSSEM_BRDF_TL'
    ! Local variables
    INTEGER  :: nZ
  
    TYPE(BRDF_type)   :: iVar_BRDF

    ! Set up
    Error_Status = SUCCESS
    ! ...Short name for angle dimensions
    nZ = size(Emissivity_TL)
 
    ! Compute tangent-linear IR sea surface emissivity
    Error_Status = NESDIS_WaterIR_Emiss_V2_TL(Wind_Speed_TL, Temperature_TL,&
                                           Emissivity_TL(1:nZ), &
                                           iVar%IRSSEM )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error computing Tangent_linear IR sea surface emissivity', &
                            Error_Status  )
      RETURN
    END IF

    Direct_Reflectivity_TL(1:nZ) = ZERO
    ! Compute the tangent-linear solar direct BRDF
    IF ( iVar%BRDF%isSolar ) THEN
        iVar_BRDF = iVar%BRDF
        Error_Status = NESDIS_IRWater_BRDF_TL(                &
                               Wind_Speed_TL,                 &
                               Direct_Reflectivity_TL(1:nZ),  &
                               iVar_BRDF)
    END IF

    ! Surface reflectance (currently assumed to be specular ALWAYS)
    Reflectivity_TL(1:nZ) = -Emissivity_TL(1:nZ)

  END FUNCTION NESDIS_IRSSEM_BRDF_V2_TL

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NESDIS_IRSSEM_BRDF_V2_AD
!
! PURPOSE:
!       Function to compute the adjoint surface emissivity and
!       reflectivity at infrared frequencies over a water surface.
!
!
!
! CALLING SEQUENCE:
!       Error_Status = NESDIS_IRSSEM_BRDF_V2_AD(      &
!                        Emissivity_AD           , &  ! Output
!                        Reflectivity_AD         , &  ! Output
!                        Direct_Reflectivity_AD  , &  ! Output
!                        Wind_Speed_AD           , &  ! Input
!                        iVar        )             &  ! Internal variable output
!
! INPUTS:
!
!      Emissivity_AD:    Surface non-polarized emissivity  adjoint
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!    Reflectivity_AD:    Surface non-polarized reflectivity adjoint
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
! Direct_Reflectivity_AD:Surface solar reflectivity adjoint
!                        UNITS:      N/A
!                        TYPE:       REAL
!                        DIMENSION:  Assumed-shape Array
!                        ATTRIBUTES: INTENT(OUT)
!
!       iVar:            Structure containing internal variables required for
!                        subsequent tangent-linear or adjoint model calls.
!                        The contents of this structure are NOT accessible
!                        outside of the CRTM_IR_Water_SfcOptics module.
!                        UNITS:      N/A
!                        TYPE:       iVar_type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! OUTPUTS:
!       Wind_speed_AD:   Ocean surface speed adjoint
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Temperature_AD:  Ocean surface temperature adjoint
!                        UNITS:      m/s
!                        TYPE:       REAL
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
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
!       Note the INTENT on the input SfcOptics_AD argument is IN OUT rather
!       than just OUT. This is necessary because components of this argument
!       may need to be zeroed out upon output.
!
!       Note the INTENT on the output Surface_AD argument is IN OUT rather
!       than just OUT. This is necessary because the argument may be defined
!       upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------


  FUNCTION NESDIS_IRSSEM_BRDF_V2_AD( &
    Emissivity_AD               , &  ! Input
    Reflectivity_AD             , &  ! Input
    Direct_Reflectivity_AD      , &  ! Input
    Wind_Speed_AD               , &  ! Output
    Temperature_AD              , &  ! Output
    iVar                        ) &  ! Internal variable input
  RESULT ( Error_Status )
    ! Arguments
    REAL(fp),                  INTENT(INOUT)  :: Emissivity_AD(:)
    REAL(fp),                  INTENT(INOUT)  :: Reflectivity_AD(:)
    REAL(fp),                  INTENT(INOUT)  :: Direct_Reflectivity_AD(:)
    REAL(fp),                  INTENT(INOUT)  :: Wind_Speed_AD
    REAL(fp),                  INTENT(INOUT)  :: Temperature_AD
    TYPE(iVar_type),           INTENT(INOUT)  :: iVar
    ! Function result
    INTEGER :: Error_Status
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NESDIS_IRSSEM_BRDF_AD'
    ! Local variables
    INTEGER  :: nZ
    TYPE(BRDF_type)   :: iVar_BRDF

    ! Set up
    Error_Status = SUCCESS
    ! ...Short name for angle dimensions
    nZ = size(Emissivity_AD)

    ! Surface reflectance (currently assumed to be specular ALWAYS)
    Emissivity_AD = Emissivity_AD - Reflectivity_AD
    Reflectivity_AD = ZERO

    ! Solar direct BRDF
    IF (iVar%BRDF%isSolar ) THEN
        iVar_BRDF =  iVar%BRDF
        Error_Status = NESDIS_IRWater_BRDF_AD( &
                               Direct_Reflectivity_AD, &  ! Input
                               Wind_Speed_AD,          &  ! Output
                               iVar_BRDF  )               ! Input  
 
    END IF

    ! Compute sdjoint IRSSEM sea surface emissivity
    Error_Status = NESDIS_WaterIR_Emiss_V2_AD( Emissivity_AD(1:nZ), &
                                           Wind_Speed_AD,        &
                                           Temperature_AD,        &
                                           iVar%IRSSEM )
    IF ( Error_Status /= SUCCESS ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error computing Adjoint IR sea surface emissivity', &
                            Error_Status  )
      RETURN
    END IF

  END FUNCTION NESDIS_IRSSEM_BRDF_V2_AD

  FUNCTION NESDIS_IRSSEM_V2_Setup( File_Name)  RESULT( Error_Status )

    CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: File_Name     ! IRSSEM netcdf coeff file name
    
    ! Function result
    INTEGER :: Error_Status
    
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NESDIS_IRSSEM_Setup'
    ! Local variables

  
    ! set path
    IF (PRESENT(File_Name) ) THEN
        Error_Status=IRSSEM_V2_Setup(TRIM(File_Name))
    ELSE     
        Error_Status = IRSSEM_V2_Setup()
    END IF    


  END FUNCTION NESDIS_IRSSEM_V2_Setup
  
  SUBROUTINE NESDIS_IRSSEM_V2_Close( )
    INTEGER :: Error_Status
    IF (IRSSEM_V2_Initialized()) THEN
      Error_Status = IRSSEM_V2_CleanUP()
    ENDIF
   
  END SUBROUTINE NESDIS_IRSSEM_V2_Close
  
  FUNCTION NESDIS_IRSSEM_V2_Initialized() RESULT( LUT_Status )
    LOGICAL :: LUT_Status
    LUT_Status=IRSSEM_V2_Initialized() 
  END FUNCTION NESDIS_IRSSEM_V2_Initialized

END MODULE NESDIS_WaterIR_PhyModel_V2
