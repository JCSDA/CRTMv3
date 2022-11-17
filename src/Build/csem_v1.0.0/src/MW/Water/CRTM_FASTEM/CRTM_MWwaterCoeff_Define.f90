!
! MWwaterCoeff_Define
!
! Module defining the MWwaterCoeff object.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 10-Nov-2011
!                       paul.vandelst@noaa.gov
!       Modified by:    Ming Chen, 10-Nov-2015
!                       ming.chen@noaa.gov

MODULE CRTM_MWwaterCoeff_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds         , ONLY: Long   => CSEM_Long, fp => CSEM_fp, &
                                      Double => CSEM_Double
  USE CSEM_Exception_Handler  , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message

  USE CSEM_FitCoeff_Define
  USE CRTM_MWwaterLUT_Define
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: CRTM_MWwaterCoeff_type

  ! Procedures
  PUBLIC :: CRTM_MWwaterCoeff_Associated
  PUBLIC :: CRTM_MWwaterCoeff_Destroy
  PUBLIC :: CRTM_MWwaterCoeff_Create
  PUBLIC :: CRTM_MWwaterCoeff_Inspect
  PUBLIC :: CRTM_MWwaterCoeff_ValidRelease
  PUBLIC :: CRTM_MWwaterCoeff_Info
  PUBLIC :: CRTM_MWwaterCoeff_DefineVersion


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
    '$Id: MWwaterCoeff_Define.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Release and version
  INTEGER, PARAMETER :: MWWATERCOEFF_RELEASE = 1  ! This determines structure and file formats.
  INTEGER, PARAMETER :: MWWATERCOEFF_VERSION = 1  ! This is just the default data version.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Data indicators
  INTEGER, PARAMETER :: DATA_MISSING = 0
  INTEGER, PARAMETER :: DATA_PRESENT = 1
  ! String lengths
  INTEGER,  PARAMETER :: ML = 256 ! Message length
  INTEGER,  PARAMETER :: SL =  80 ! String length


  ! ---------------------------------
  ! MWwaterCoeff data type definition
  ! ---------------------------------
  !:tdoc+:
  TYPE :: CRTM_MWwaterCoeff_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = MWWATERCOEFF_RELEASE
    INTEGER(Long) :: Version = MWWATERCOEFF_VERSION
    ! Derived type components
    TYPE(CSEM_FitCoeff_1D_type) :: FCCoeff     ! Foam coverage          fitting coefficients
    TYPE(CSEM_FitCoeff_1D_type) :: FRCoeff     ! Foam reflectivity      fitting coefficients
    TYPE(CSEM_FitCoeff_3D_type) :: RCCoeff     ! Reflection correction  fitting coefficients
    TYPE(CSEM_FitCoeff_3D_type) :: AZCoeff     ! Azimuth emissivity     fitting coefficients
    TYPE(CSEM_FitCoeff_1D_type) :: SSCCoeff    ! Small-scale correction fitting coefficients
    TYPE(CSEM_FitCoeff_3D_type) :: LSCCoeff    ! Large-scale correction fitting coefficients
    TYPE(MWwaterLUT_type)  :: LUT         ! Emissivity look-up table
  END TYPE CRTM_MWwaterCoeff_type
  !:tdoc-:


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                           ## PUBLIC PROCEDURES ##                          ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_Associated
!
! PURPOSE:
!       Pure function to test the status of the allocatable components
!       of the MWwaterCoeff structure.
!
! CALLING SEQUENCE:
!       Status = MWwaterCoeff_Associated( MWwaterCoeff )
!
! OBJECTS:
!       MWwaterCoeff:  Structure which is to have its member's
!                      status tested.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:        The return value is a logical value indicating the
!                      status of the components.
!                       .TRUE.  - if ALL of the MWwaterCoeff allocatable members
!                                 are in use.
!                       .FALSE. - if ANY of the MWwaterCoeff allocatable members
!                                 are not in use.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE FUNCTION CRTM_MWwaterCoeff_Associated( self ) RESULT( Status )
    TYPE(CRTM_MWwaterCoeff_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION CRTM_MWwaterCoeff_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_Destroy
!
! PURPOSE:
!       Pure subroutine to re-initialize MWwaterCoeff objects.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_Destroy( MWwaterCoeff )
!
! OBJECTS:
!       MWwaterCoeff: Re-initialized MWwaterCoeff structure.
!                     UNITS:      N/A
!                     TYPE:       MWwaterCoeff_type
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE CRTM_MWwaterCoeff_Destroy( self )
    TYPE(CRTM_MWwaterCoeff_type), INTENT(OUT) :: self
    CALL CSEM_FitCoeff_Destroy(self%FCCoeff)
    CALL CSEM_FitCoeff_Destroy(self%FRCoeff)
    CALL CSEM_FitCoeff_Destroy(self%RCCoeff)
    CALL CSEM_FitCoeff_Destroy(self%AZCoeff)
    CALL CSEM_FitCoeff_Destroy(self%SSCCoeff)
    CALL CSEM_FitCoeff_Destroy(self%LSCCoeff)
    self%Is_Allocated = .FALSE.
  END SUBROUTINE CRTM_MWwaterCoeff_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_Create
!
! PURPOSE:
!       Pure subroutine to create a valid instance of an MWwaterCoeff object.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_Create( MWwaterCoeff )
!
! OBJECTS:
!       MWwaterCoeff:       MWwaterCoeff object structure.
!                           UNITS:      N/A
!                           TYPE:       MWwaterCoeff_type
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE CRTM_MWwaterCoeff_Create( &
     self,           &
     ndim_subgrp,    &
     dims_subgrp)

    ! Arguments

    INTEGER,   INTENT(IN)  :: ndim_subgrp(:)
    INTEGER,   INTENT(IN)  :: dims_subgrp(:,:)
    TYPE(CRTM_MWwaterCoeff_type), INTENT(IN OUT) :: self
  
    CALL CSEM_FitCoeff_Create(self%FCCoeff, dims_subgrp(1:ndim_subgrp(1),1))
    CALL CSEM_FitCoeff_Create(self%FRCoeff, dims_subgrp(1:ndim_subgrp(2),2))
    CALL CSEM_FitCoeff_Create(self%RCCoeff, dims_subgrp(1:ndim_subgrp(3),3))
    CALL CSEM_FitCoeff_Create(self%AZCoeff, dims_subgrp(1:ndim_subgrp(4),4))
    CALL CSEM_FitCoeff_Create(self%SSCCoeff,dims_subgrp(1:ndim_subgrp(5),5))
    CALL CSEM_FitCoeff_Create(self%LSCCoeff,dims_subgrp(1:ndim_subgrp(6),6))

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.
  END SUBROUTINE CRTM_MWwaterCoeff_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a MWwaterCoeff object to stdout.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_Inspect( MWwaterCoeff )
!
! OBJECTS:
!       MWwaterCoeff:  MWwaterCoeff object to display.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_MWwaterCoeff_Inspect( self, pause )
    TYPE(CRTM_MWwaterCoeff_type), INTENT(IN) :: self
    LOGICAL,       OPTIONAL, INTENT(IN) :: pause
    LOGICAL :: wait

    wait = .FALSE.
    IF ( PRESENT(pause) ) wait = pause

    WRITE(*,'(1x,"MWwaterCoeff OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version : ",i0,".",i0)') self%Release, self%Version
    IF ( .NOT. CRTM_MWwaterCoeff_Associated(self) ) RETURN
    ! Derived types
    IF ( CSEM_FitCoeff_Associated(self%FCCoeff ) ) THEN
      WRITE(*,'(1x,"Foam coverage cofficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%FCCoeff )
    END IF
    IF ( CSEM_FitCoeff_Associated(self%FRCoeff ) ) THEN
      WRITE(*,'(1x,"Foam reflectivity cofficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%FRCoeff )
    END IF
    IF ( CSEM_FitCoeff_Associated(self%RCCoeff) ) THEN
      WRITE(*,'(1x,"Reflection correction cofficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%RCCoeff)
    END IF
    IF ( CSEM_FitCoeff_Associated(self%AZCoeff) ) THEN
      WRITE(*,'(1x,"Azimuth emissivity coefficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%AZCoeff)
    END IF
    IF ( CSEM_FitCoeff_Associated(self%SSCCoeff) ) THEN
      WRITE(*,'(1x,"Small-scale correction coefficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%SSCCoeff)
    END IF
    IF ( CSEM_FitCoeff_Associated(self%LSCCoeff) ) THEN
      WRITE(*,'(1x,"Large-scale correction coefficients, ")',ADVANCE='NO')
      CALL CSEM_FitCoeff_Inspect(self%LSCCoeff)
    END IF
    IF ( MWwaterLUT_Associated(self%LUT) ) THEN
      WRITE(*,'(1x,"Emissivity look-up table, ")',ADVANCE='NO')
      CALL MWwaterLUT_Inspect(self%LUT,pause=pause)
    END IF
  END SUBROUTINE CRTM_MWwaterCoeff_Inspect



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_ValidRelease
!
! PURPOSE:
!       Function to check the MWwaterCoeff Release value.
!
! CALLING SEQUENCE:
!       IsValid = MWwaterCoeff_ValidRelease( MWwaterCoeff )
!
! INPUTS:
!       MWwaterCoeff:  MWwaterCoeff object for which the Release component
!                      is to be checked.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       IsValid:       Logical value defining the release validity.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Scalar
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_MWwaterCoeff_ValidRelease( self ) RESULT( IsValid )
    ! Arguments
    TYPE(CRTM_MWwaterCoeff_type), INTENT(IN) :: self
    ! Function result
    LOGICAL :: IsValid
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'MWwaterCoeff_ValidRelease'
    ! Local variables
    CHARACTER(ML) :: msg

    ! Set up
    IsValid = .TRUE.


    ! Check release is not too old
    IF ( self%Release < MWWATERCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An MWwaterCoeff data update is needed. ", &
                  &"MWwaterCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, MWWATERCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF


    ! Check release is not too new
    IF ( self%Release > MWWATERCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An MWwaterCoeff software update is needed. ", &
                  &"MWwaterCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, MWWATERCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF

  END FUNCTION CRTM_MWwaterCoeff_ValidRelease


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_Info
!
! PURPOSE:
!       Subroutine to return a string containing version and dimension
!       information about a MWwaterCoeff object.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_Info( MWwaterCoeff, Info )
!
! OBJECTS:
!       MWwaterCoeff:  MWwaterCoeff object about which info is required.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Info:          String containing version and dimension information
!                      about the MWwaterCoeff object.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_MWwaterCoeff_Info( self, Info )
    ! Arguments
    TYPE(CRTM_MWwaterCoeff_type), INTENT(IN)  :: self
    CHARACTER(*),            INTENT(OUT) :: Info
    ! Parameters
    INTEGER, PARAMETER :: CARRIAGE_RETURN = 13
    INTEGER, PARAMETER :: LINEFEED = 10
    ! Local variables
    CHARACTER(2000) :: Long_String

    ! Write the required data to the local string
    WRITE( Long_String, &
           '(a,1x,"MWwaterCoeff RELEASE.VERSION: ",i0,".",i0 )' ) &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%Release, self%Version

    ! Trim the output based on the
    ! dummy argument string length
    Info = Long_String(1:MIN(LEN(Info), LEN_TRIM(Long_String)))

  END SUBROUTINE CRTM_MWwaterCoeff_Info


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_DefineVersion
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_DefineVersion( Id )
!
! OUTPUTS:
!       Id:    Character string containing the version Id information
!              for the module.
!              UNITS:      N/A
!              TYPE:       CHARACTER(*)
!              DIMENSION:  Scalar
!              ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_MWwaterCoeff_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CRTM_MWwaterCoeff_DefineVersion


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_SetValue
!
! PURPOSE:
!       Subroutine to set the contents of a valid MWwaterCoeff object.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_SetValue( MWwaterCoeff, &
!                                   FCCoeff    = FCCoeff   , &
!                                   FRCoeff    = FRCoeff   , &
!                                   RCCoeff    = RCCoeff   , &
!                                   AZCoeff    = AZCoeff   , &
!                                   SSCCoeff   = SSCCoeff  , &
!                                   LSCCoeff   = LSCCoeff  , &
!                                   MWwaterLUT = MWwaterLUT  )
! OBJECTS:
!       MWwaterCoeff:  Valid, allocated MWwaterCoeff object for which
!                      values are to be set.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL INPUTS:
!       FCCoeff:       Object containing the foam coverage fitting coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       FRCoeff:       Object containing the foam reflectivity fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       RCCoeff:       Object containing the reflection correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       AZCoeff:       Object containing the azimuth emissivity fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       SSCCoeff:      Object containing the small-scale correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       LSCCoeff:      Object containing the large-scale correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       MWwaterLUT:    Object containing the emissivity look-up table.
!                      UNITS:      N/A
!                      TYPE:       MWwaterLUT_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN), OPTIONAL
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_MWwaterCoeff_SetValue( &
    self      , &  ! In/Output
    FCCoeff   , &  ! Optional input
    FRCoeff   , &  ! Optional input
    RCCoeff   , &  ! Optional input
    AZCoeff   , &  ! Optional input
    SSCCoeff  , &  ! Optional input
    LSCCoeff  , &  ! Optional input
    MWwaterLUT  )  ! Optional input
    ! Arguments
    TYPE(CRTM_MWwaterCoeff_type)         , INTENT(IN OUT) :: self
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(IN)     :: FCCoeff
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(IN)     :: FRCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(IN)     :: RCCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(IN)     :: AZCoeff
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(IN)     :: SSCCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(IN)     :: LSCCoeff
    TYPE(MWwaterLUT_type) , OPTIONAL, INTENT(IN)     :: MWwaterLUT

    IF ( PRESENT(FCCoeff   ) ) self%FCCoeff  = FCCoeff
    IF ( PRESENT(FRCoeff   ) ) self%FRCoeff  = FRCoeff
    IF ( PRESENT(RCCoeff   ) ) self%RCCoeff  = RCCoeff
    IF ( PRESENT(AZCoeff   ) ) self%AZCoeff  = AZCoeff
    IF ( PRESENT(SSCCoeff  ) ) self%SSCCoeff = SSCCoeff
    IF ( PRESENT(LSCCoeff  ) ) self%LSCCoeff = LSCCoeff
    IF ( PRESENT(MWwaterLUT) ) self%LUT      = MWwaterLUT

  END SUBROUTINE CRTM_MWwaterCoeff_SetValue


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       MWwaterCoeff_GetValue
!
! PURPOSE:
!       Subroutine to get the contents of a valid MWwaterCoeff object.
!
! CALLING SEQUENCE:
!       CALL MWwaterCoeff_GetValue( MWwaterCoeff, &
!                                   FCCoeff    = FCCoeff   , &
!                                   FRCoeff    = FRCoeff   , &
!                                   RCCoeff    = RCCoeff   , &
!                                   AZCoeff    = AZCoeff   , &
!                                   SSCCoeff   = SSCCoeff  , &
!                                   LSCCoeff   = LSCCoeff  , &
!                                   MWwaterLUT = MWwaterLUT  )
!
! OBJECTS:
!       MWwaterCoeff:  Valid, allocated MWwaterCoeff object from which
!                      values are to be retrieved.
!                      UNITS:      N/A
!                      TYPE:       MWwaterCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN OUT)
!
! OPTIONAL OUTPUTS:
!       FCCoeff:       Object containing the foam coverage fitting coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       FRCoeff:       Object containing the foam reflectivity fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       RCCoeff:       Object containing the reflection correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       AZCoeff:       Object containing the azimuth emissivity fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       SSCCoeff:      Object containing the small-scale correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_1D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       LSCCoeff:      Object containing the large-scale correction fitting
!                      coefficients.
!                      UNITS:      N/A
!                      TYPE:       CSEM_FitCoeff_3D_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       MWwaterLUT:    Object containing the emissivity look-up table.
!                      UNITS:      N/A
!                      TYPE:       MWwaterLUT_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_MWwaterCoeff_GetValue( &
    self      , &  ! Input
    FCCoeff   , &  ! Optional output
    FRCoeff   , &  ! Optional output
    RCCoeff   , &  ! Optional output
    AZCoeff   , &  ! Optional output
    SSCCoeff  , &  ! Optional output
    LSCCoeff  , &  ! Optional output
    MWwaterLUT  )  ! Optional output
    ! Arguments
    TYPE(CRTM_MWwaterCoeff_type)         , INTENT(IN)  :: self
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(OUT) :: FCCoeff
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(OUT) :: FRCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(OUT) :: RCCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(OUT) :: AZCoeff
    TYPE(CSEM_FitCoeff_1D_type), OPTIONAL, INTENT(OUT) :: SSCCoeff
    TYPE(CSEM_FitCoeff_3D_type), OPTIONAL, INTENT(OUT) :: LSCCoeff
    TYPE(MWwaterLUT_type) , OPTIONAL, INTENT(OUT) :: MWwaterLUT

    IF ( PRESENT(FCCoeff   ) ) FCCoeff    = self%FCCoeff
    IF ( PRESENT(FRCoeff   ) ) FRCoeff    = self%FRCoeff
    IF ( PRESENT(RCCoeff   ) ) RCCoeff    = self%RCCoeff
    IF ( PRESENT(AZCoeff   ) ) AZCoeff    = self%AZCoeff
    IF ( PRESENT(SSCCoeff  ) ) SSCCoeff   = self%SSCCoeff
    IF ( PRESENT(LSCCoeff  ) ) LSCCoeff   = self%LSCCoeff
    IF ( PRESENT(MWwaterLUT) ) MWwaterLUT = self%LUT

  END SUBROUTINE CRTM_MWwaterCoeff_GetValue




END MODULE CRTM_MWwaterCoeff_Define
