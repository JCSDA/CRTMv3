!
! CSEM_FitCoeff_Define
!
! Module defining the FitCoeff objects.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 17-Nov-2011
!                       paul.vandelst@noaa.gov
 
MODULE CSEM_FitCoeff_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds,        ONLY: Long   => CSEM_Long, fp => CSEM_fp, &
                                    Double => CSEM_Double
  USE CSEM_Exception_Handler, ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message

  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Parameters
  PUBLIC :: FITCOEFF_MAX_N_DIMENSIONS
  ! Datatypes
  PUBLIC :: CSEM_FitCoeff_1D_type
  PUBLIC :: CSEM_FitCoeff_2D_type
  PUBLIC :: CSEM_FitCoeff_3D_type

  ! Procedures
  PUBLIC :: CSEM_FitCoeff_Associated
  PUBLIC :: CSEM_FitCoeff_Destroy
  PUBLIC :: CSEM_FitCoeff_Create
  PUBLIC :: CSEM_FitCoeff_Inspect
  PUBLIC :: CSEM_FitCoeff_ValidRelease
  PUBLIC :: CSEM_FitCoeff_DefineVersion


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE CSEM_FitCoeff_Associated
    MODULE PROCEDURE FitCoeff_1D_Associated
    MODULE PROCEDURE FitCoeff_2D_Associated
    MODULE PROCEDURE FitCoeff_3D_Associated
  END INTERFACE CSEM_FitCoeff_Associated
  
  INTERFACE CSEM_FitCoeff_Destroy
    MODULE PROCEDURE FitCoeff_1D_Destroy
    MODULE PROCEDURE FitCoeff_2D_Destroy
    MODULE PROCEDURE FitCoeff_3D_Destroy
  END INTERFACE CSEM_FitCoeff_Destroy
  
  INTERFACE CSEM_FitCoeff_Create
    MODULE PROCEDURE FitCoeff_1D_Create
    MODULE PROCEDURE FitCoeff_2D_Create
    MODULE PROCEDURE FitCoeff_3D_Create
  END INTERFACE CSEM_FitCoeff_Create
  
  INTERFACE CSEM_FitCoeff_Inspect
    MODULE PROCEDURE FitCoeff_1D_Inspect
    MODULE PROCEDURE FitCoeff_2D_Inspect
    MODULE PROCEDURE FitCoeff_3D_Inspect
  END INTERFACE CSEM_FitCoeff_Inspect
  
  INTERFACE CSEM_FitCoeff_ValidRelease
    MODULE PROCEDURE FitCoeff_1D_ValidRelease
    MODULE PROCEDURE FitCoeff_2D_ValidRelease
    MODULE PROCEDURE FitCoeff_3D_ValidRelease
  END INTERFACE CSEM_FitCoeff_ValidRelease
  


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
    '$Id: FitCoeff_Define.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'
  ! Release and version
  INTEGER, PARAMETER :: FITCOEFF_RELEASE = 1  ! This determines structure and file formats.
  INTEGER, PARAMETER :: FITCOEFF_VERSION = 1  ! This is just the default data version.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! String lengths
  INTEGER, PARAMETER :: ML = 256 ! Message length
  INTEGER, PARAMETER :: SL =  80 ! String length
  ! Maximum number of dimensions
  INTEGER, PARAMETER :: FITCOEFF_MAX_N_DIMENSIONS = 3  ! Only implemented up to 3-D arrays so far


  ! ----------------------------------
  ! FitCoeff data type definitions
  ! ----------------------------------
  !:tdoc+:
  TYPE :: CSEM_FitCoeff_1D_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = FITCOEFF_RELEASE
    INTEGER(Long) :: Version = FITCOEFF_VERSION
    ! Dimensions
    INTEGER(Long) :: Dimensions(1) = 0
    ! Data
    REAL(Double), ALLOCATABLE :: C(:)
  END TYPE CSEM_FitCoeff_1D_type
  !:tdoc-:

  !:tdoc+:
  TYPE :: CSEM_FitCoeff_2D_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = FITCOEFF_RELEASE
    INTEGER(Long) :: Version = FITCOEFF_VERSION
    ! Dimensions
    INTEGER(Long) :: Dimensions(2) = 0
    ! Data
    REAL(Double), ALLOCATABLE :: C(:,:)
  END TYPE CSEM_FitCoeff_2D_type
  !:tdoc-:

  !:tdoc+:
  TYPE :: CSEM_FitCoeff_3D_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = FITCOEFF_RELEASE
    INTEGER(Long) :: Version = FITCOEFF_VERSION
    ! Dimensions
    INTEGER(Long) :: Dimensions(3) = 0
    ! Data
    REAL(Double), ALLOCATABLE :: C(:,:,:)
  END TYPE CSEM_FitCoeff_3D_type
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
!       CSEM_FitCoeff_Associated
!
! PURPOSE:
!       Pure function to test the status of the allocatable components
!       of the FitCoeff structure.
!
! CALLING SEQUENCE:
!       Status = CSEM_FitCoeff_Associated( FitCoeff )
!
! OBJECTS:
!       FitCoeff:      Structure which is to have its member's
!                      status tested.
!                      UNITS:      N/A
!                      TYPE:       Any FitCoeff type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:        The return value is a logical value indicating the
!                      status of the components.
!                       .TRUE.  - if ANY of the FitCoeff allocatable members
!                                 are in use.
!                       .FALSE. - if ALL of the FitCoeff allocatable members
!                                 are not in use.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE FUNCTION FitCoeff_1D_Associated( self ) RESULT( Status )
    TYPE(CSEM_FitCoeff_1D_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION FitCoeff_1D_Associated

  PURE FUNCTION FitCoeff_2D_Associated( self ) RESULT( Status )
    TYPE(CSEM_FitCoeff_2D_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION FitCoeff_2D_Associated

  PURE FUNCTION FitCoeff_3D_Associated( self ) RESULT( Status )
    TYPE(CSEM_FitCoeff_3D_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION FitCoeff_3D_Associated

  
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_FitCoeff_Destroy
!
! PURPOSE:
!       Pure subroutine to re-initialize FitCoeff objects.
!
! CALLING SEQUENCE:
!       CALL CSEM_FitCoeff_Destroy( FitCoeff )
!
! OBJECTS:
!       FitCoeff:     Re-initialized FitCoeff structure.
!                     UNITS:      N/A
!                     TYPE:       Any FitCoeff type
!                     DIMENSION:  Scalar
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE FitCoeff_1D_Destroy( self )
    TYPE(CSEM_FitCoeff_1D_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%Dimensions = 0

  END SUBROUTINE FitCoeff_1D_Destroy

  PURE SUBROUTINE FitCoeff_2D_Destroy( self )
    TYPE(CSEM_FitCoeff_2D_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%Dimensions = 0
  END SUBROUTINE FitCoeff_2D_Destroy

  PURE SUBROUTINE FitCoeff_3D_Destroy( self )
    TYPE(CSEM_FitCoeff_3D_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%Dimensions = 0
  END SUBROUTINE FitCoeff_3D_Destroy



!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_FitCoeff_Create
!
! PURPOSE:
!       Pure subroutine to create an instance of a FitCoeff object.
!
! CALLING SEQUENCE:
!       CALL CSEM_FitCoeff_Create( FitCoeff, Dimensions )
!
! OBJECTS:
!       FitCoeff:   FitCoeff object structure.
!                   UNITS:      N/A
!                   TYPE:       Any FitCoeff type
!                   DIMENSION:  Scalar
!                   ATTRIBUTES: INTENT(OUT)
!
! INPUTS:
!       Dimensions: Dimension vector for the fitting coefficient array.
!                   The number of elements of this array must agree with
!                   the rank of the FitCoeff datatype specified, e.g. 2D
!                   type requires 2 dimensions specified.
!                   Values must be > 0.
!                   UNITS:      N/A
!                   TYPE:       INTEGER
!                   DIMENSION:  Rank
!                   ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  PURE SUBROUTINE FitCoeff_1D_Create( &
    self      , &  ! Output
    dimensions  )  ! Input
    ! Arguments
    TYPE(CSEM_FitCoeff_1D_type), INTENT(OUT) :: self
    INTEGER               , INTENT(IN)  :: dimensions(1)
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( ANY(dimensions < 1) ) RETURN

    ! Perform the allocation
    ALLOCATE( self%C(dimensions(1)), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    ! ...Dimensions
    self%Dimensions = dimensions
    ! ...Arrays
    self%C = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE FitCoeff_1D_Create


  PURE SUBROUTINE FitCoeff_2D_Create( &
    self      , &  ! Output
    dimensions  )  ! Input
    ! Arguments
    TYPE(CSEM_FitCoeff_2D_type), INTENT(OUT) :: self
    INTEGER               , INTENT(IN)  :: dimensions(2)
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( ANY(dimensions < 1) ) RETURN

    ! Perform the allocation
    ALLOCATE( self%C(dimensions(1), dimensions(2)), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    ! ...Dimensions
    self%Dimensions = dimensions
    ! ...Arrays
    self%C = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE FitCoeff_2D_Create


  PURE SUBROUTINE FitCoeff_3D_Create( &
    self      , &  ! Output
    dimensions  )  ! Input
    ! Arguments
    TYPE(CSEM_FitCoeff_3D_type), INTENT(OUT) :: self
    INTEGER               , INTENT(IN)  :: dimensions(3)
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( ANY(dimensions < 1) ) RETURN

    ! Perform the allocation
    ALLOCATE( self%C(dimensions(1), dimensions(2), dimensions(3)), STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    ! ...Dimensions
    self%Dimensions = dimensions
    ! ...Arrays
    self%C = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE FitCoeff_3D_Create




!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_FitCoeff_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a FitCoeff object to stdout.
!
! CALLING SEQUENCE:
!       CALL CSEM_FitCoeff_Inspect( FitCoeff )
!
! OBJECTS:
!       FitCoeff:      FitCoeff object to display.
!                      UNITS:      N/A
!                      TYPE:       Any FitCoeff type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE FitCoeff_1D_Inspect( self )
    TYPE(CSEM_FitCoeff_1D_type), INTENT(IN) :: self
    WRITE(*,'(1x,"FitCoeff 1D OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version : ",i0,".",i0)') self%Release, self%Version
    ! Dimensions
    WRITE(*,'(3x,"Dimensions : ",10(i5,:))') self%Dimensions
    IF ( .NOT. FitCoeff_1D_Associated(self) ) RETURN
    ! Coefficient data
    WRITE(*,'(3x,"Coefficients:")')
    WRITE(*,'(5(1x,es13.6,:))') self%C
  END SUBROUTINE FitCoeff_1D_Inspect


  SUBROUTINE FitCoeff_2D_Inspect( self )
    TYPE(CSEM_FitCoeff_2D_type), INTENT(IN) :: self
    INTEGER :: i
    WRITE(*,'(1x,"FitCoeff 2D OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version : ",i0,".",i0)') self%Release, self%Version
    ! Dimensions
    WRITE(*,'(3x,"Dimensions : ",10(i5,:))') self%Dimensions
    IF ( .NOT. FitCoeff_2D_Associated(self) ) RETURN
    ! Coefficient data
    WRITE(*,'(3x,"Coefficients:")')
    DO i = 1, self%Dimensions(2)
      WRITE(*,'(5x,"Outer dimension = ",i0," of ",i0)') i, self%Dimensions(2)
      WRITE(*,'(5(1x,es13.6,:))') self%C(:,i)
    END DO
  END SUBROUTINE FitCoeff_2D_Inspect


  SUBROUTINE FitCoeff_3D_Inspect( self )
    TYPE(CSEM_FitCoeff_3D_type), INTENT(IN) :: self
    INTEGER :: i, j
    WRITE(*,'(1x,"FitCoeff 3D OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version : ",i0,".",i0)') self%Release, self%Version
    ! Dimensions
    WRITE(*,'(3x,"Dimensions : ",10(i5,:))') self%Dimensions
    IF ( .NOT. FitCoeff_3D_Associated(self) ) RETURN
    ! Coefficient data
    WRITE(*,'(3x,"Coefficients:")')
    DO j = 1, self%Dimensions(3)
      WRITE(*,'(5x,"Outer dimension = ",i0," of ",i0)') j, self%Dimensions(3)
      DO i = 1, self%Dimensions(2)
        WRITE(*,'(7x,"Middle dimension = ",i0," of ",i0)') i, self%Dimensions(2)
        WRITE(*,'(5(1x,es13.6,:))') self%C(:,i,j)
      END DO
    END DO
  END SUBROUTINE FitCoeff_3D_Inspect


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_FitCoeff_ValidRelease
!
! PURPOSE:
!       Function to check the FitCoeff Release value.
!
! CALLING SEQUENCE:
!       IsValid = CSEM_FitCoeff_ValidRelease( FitCoeff )
!
! INPUTS:
!       FitCoeff:      FitCoeff object for which the Release component
!                      is to be checked.
!                      UNITS:      N/A
!                      TYPE:       Any FitCoeff type
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

  FUNCTION ValidRelease( Release ) RESULT( IsValid )
    ! Arguments
    INTEGER, INTENT(IN) :: Release
    ! Function result
    LOGICAL :: IsValid
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'FitCoeff_ValidRelease'
    ! Local variables
    CHARACTER(ML) :: msg

    ! Set up
    IsValid = .TRUE.

    ! Check release is not too old
    IF ( Release < FITCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("A FitCoeff data update is needed. ", &
                  &"FitCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  Release, FITCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF

    ! Check release is not too new
    IF ( Release > FITCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("A FitCoeff software update is needed. ", &
                  &"FitCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  Release, FITCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF
  END FUNCTION ValidRelease


  FUNCTION FitCoeff_1D_ValidRelease( self ) RESULT( IsValid )
    TYPE(CSEM_FitCoeff_1D_type), INTENT(IN) :: self
    LOGICAL :: IsValid
    IsValid = ValidRelease( self%Release )
  END FUNCTION FitCoeff_1D_ValidRelease


  FUNCTION FitCoeff_2D_ValidRelease( self ) RESULT( IsValid )
    TYPE(CSEM_FitCoeff_2D_type), INTENT(IN) :: self
    LOGICAL :: IsValid
    IsValid = ValidRelease( self%Release )
  END FUNCTION FitCoeff_2D_ValidRelease


  FUNCTION FitCoeff_3D_ValidRelease( self ) RESULT( IsValid )
    TYPE(CSEM_FitCoeff_3D_type), INTENT(IN) :: self
    LOGICAL :: IsValid
    IsValid = ValidRelease( self%Release )
  END FUNCTION FitCoeff_3D_ValidRelease


 
!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_FitCoeff_DefineVersion
!
! PURPOSE:
!       Subroutine to return the module version information.
!
! CALLING SEQUENCE:
!       CALL CSEM_FitCoeff_DefineVersion( Id )
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

  SUBROUTINE CSEM_FitCoeff_DefineVersion( Id )
    CHARACTER(*), INTENT(OUT) :: Id
    Id = MODULE_VERSION_ID
  END SUBROUTINE CSEM_FitCoeff_DefineVersion


  

END MODULE CSEM_FitCoeff_Define
