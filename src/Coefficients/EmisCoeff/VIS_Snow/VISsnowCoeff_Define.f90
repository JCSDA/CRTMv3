!
! VISsnowCoeff_Define
!
! Module defining the VISsnowCoeff object to hold coefficient
! data for the visible and near-IR snow surface reflectivity models.
!
!
! CREATION HISTORY:
!       Written by:     Cheng Dang, May-2026
!                       dangch@ucar.edu

MODULE VISsnowCoeff_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds           , ONLY: fp, Long, Double
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)
  USE File_Utility         , ONLY: File_Open, File_Exists
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Datatypes
  PUBLIC :: VISsnowCoeff_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Procedures
  PUBLIC :: VISsnowCoeff_Associated
  PUBLIC :: VISsnowCoeff_Destroy
  PUBLIC :: VISsnowCoeff_Create
  PUBLIC :: VISsnowCoeff_Inspect
  PUBLIC :: VISsnowCoeff_ValidRelease
  PUBLIC :: VISsnowCoeff_Info


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE VISsnowCoeff_Equal
  END INTERFACE OPERATOR(==)


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Current valid release and version
  INTEGER, PARAMETER :: VISsnowCOEFF_RELEASE = 1  ! This determines structure and file formats.
  INTEGER, PARAMETER :: VISsnowCOEFF_VERSION = 1  ! This is just the default data version.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Literal constants
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  ! String lengths
  INTEGER,  PARAMETER :: ML = 256 ! Message length


  ! ----------------------------------
  ! VISsnowCoeff_type data type definitions
  ! ----------------------------------
  !:tdoc+:
  TYPE :: VISsnowCoeff_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = VISsnowCOEFF_RELEASE
    INTEGER(Long) :: Version = VISsnowCOEFF_VERSION
    ! Surface classification name
    CHARACTER(ML) :: Classification_Name = ''
    ! Dimensions
    INTEGER(Long) :: n_Angles      = 0   ! I dimension
    INTEGER(Long) :: n_Frequencies = 0   ! L dimension
    INTEGER(Long) :: n_Grain_Sizes = 0   ! G dimension
    INTEGER(Long) :: n_Depths      = 0   ! T dimension
    INTEGER(Long) :: n_Densities   = 0   ! J dimension
    ! Dimensional vectors
    REAL(Double), ALLOCATABLE :: Angle(:)        ! I
    REAL(Double), ALLOCATABLE :: Frequency(:)    ! L
    REAL(Double), ALLOCATABLE :: Grain_Size(:)   ! G
    REAL(Double), ALLOCATABLE :: Depth(:)        ! T
    REAL(Double), ALLOCATABLE :: Density(:)      ! J
    ! Reflectance LUT data
    REAL(Double), ALLOCATABLE :: Reflectance(:,:,:,:,:)  ! I x L x G x T x J
  END TYPE VISsnowCoeff_type
  !:tdoc-:


CONTAINS


!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_Associated
!
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of the VISsnowCoeff structure.
!
! CALLING SEQUENCE:
!       Status = VISsnowCoeff_Associated( VISsnowCoeff )
!
! OBJECTS:
!       VISsnowCoeff:  Structure which is to have its member's
!                      status tested.
!                      UNITS:      N/A
!                      TYPE:       VISsnowCoeff_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       Status:        The return value is a logical value indicating the
!                      status of the NLTE members.
!                       .TRUE.  - if ANY of the VISsnowCoeff allocatable members
!                                 are in use.
!                       .FALSE. - if ALL of the VISsnowCoeff allocatable members
!                                 are not in use.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as input
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION VISsnowCoeff_Associated( self ) RESULT( Status )
    TYPE(VISsnowCoeff_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION VISsnowCoeff_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_Destroy
!
! PURPOSE:
!       Elemental subroutine to re-initialize VISsnowCoeff objects.
!
! CALLING SEQUENCE:
!       CALL VISsnowCoeff_Destroy( VISsnowCoeff )
!
! OBJECTS:
!       VISsnowCoeff: Re-initialized VISsnowCoeff structure.
!                     UNITS:      N/A
!                     TYPE:       VISsnowCoeff_type
!                     DIMENSION:  Scalar or any rank
!                     ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE VISsnowCoeff_Destroy( self )
    TYPE(VISsnowCoeff_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
  END SUBROUTINE VISsnowCoeff_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_Create
!
! PURPOSE:
!       Elemental subroutine to create an instance of an VISsnowCoeff object.
!
! CALLING SEQUENCE:
!       CALL VISsnowCoeff_Create( VISsnowCoeff   , &
!                                 n_Angles     , &
!                                 n_Frequencies, &
!                                 n_Grain_Sizes, &
!                                 n_Temperature  )
!
! OBJECTS:
!       VISsnowCoeff:   VISsnowCoeff object structure.
!                       UNITS:      N/A
!                       TYPE:       VISsnowCoeff_type
!                       DIMENSION:  Scalar or any rank
!                       ATTRIBUTES: INTENT(OUT)
!
! INPUTS:
!       n_Angles:       Number of angles dimension.
!                       Must be > 0.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Conformable with the VISsnowCoeff object
!                       ATTRIBUTES: INTENT(IN)
!
!       n_Frequencies:  Number of frequencies dimension.
!                       Must be > 0.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Conformable with the VISsnowCoeff object
!                       ATTRIBUTES: INTENT(IN)
!
!       n_Grain_Sizes:  Number of Grain Sizes dimension.
!                       Must be > 0.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Conformable with the VISsnowCoeff object
!                       ATTRIBUTES: INTENT(IN)
!
!       n_Temperature:  Number oftemperature dimension.
!                       Must be > 0.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Conformable with the VISsnowCoeff object
!                       ATTRIBUTES: INTENT(IN)
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE VISsnowCoeff_Create( &
    self         , &  ! Output
    n_Angles     , &  ! Input
    n_Frequencies, &  ! Input
    n_Grain_Sizes, &  ! Input
    n_Depths      , &  ! Input
    n_Densities   )  ! Input
    ! Arguments
    TYPE(VISsnowCoeff_type) , INTENT(OUT) :: self
    INTEGER                 , INTENT(IN)  :: n_Angles
    INTEGER                 , INTENT(IN)  :: n_Frequencies
    INTEGER                 , INTENT(IN)  :: n_Grain_Sizes
    INTEGER                 , INTENT(IN)  :: n_Depths
    INTEGER                 , INTENT(IN)  :: n_Densities
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( self%Is_Allocated .OR. &
         n_Angles      < 1 .OR. &
         n_Frequencies < 1 .OR. &
         n_Grain_Sizes < 1 .OR. &
         n_Depths      < 1 .OR. &
         n_Densities   < 1) RETURN

    ! Perform the allocation
    ALLOCATE( self%Angle( n_Angles ), &
              self%Frequency( n_Frequencies ), &
              self%Grain_Size( n_Grain_Sizes ), &
              self%Depth( n_Depths ), &
              self%Density( n_Densities ), &
              self%Reflectance( n_Angles, n_Frequencies, n_Grain_Sizes, n_Depths, n_Densities ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN


    ! Initialise
    ! ...Dimensions
    self%n_Angles      = n_Angles
    self%n_Frequencies = n_Frequencies
    self%n_Grain_Sizes = n_Grain_Sizes
    self%n_Depths      = n_Depths
    self%n_Densities   = n_Densities
    ! ...Arrays
    self%Angle        = ZERO
    self%Frequency    = ZERO
    self%Grain_Size   = ZERO
    self%Depth        = ZERO
    self%Density      = ZERO
    self%Reflectance  = ZERO

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.

  END SUBROUTINE VISsnowCoeff_Create


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of a VISsnowCoeff object to stdout.
!
! CALLING SEQUENCE:
!       CALL VISsnowCoeff_Inspect( VISsnowCoeff )
!
! OBJECTS:
!       VISsnowCoeff:  VISsnowCoeff object to display.
!                      UNITS:      N/A
!                      TYPE:       VISsnowCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE VISsnowCoeff_Inspect( self )
    TYPE(VISsnowCoeff_type), INTENT(IN) :: self
    INTEGER :: i2, i3, i4, i5
    WRITE(*,'(1x,"VISsnowCoeff OBJECT")')
    ! Release/version info
    WRITE(*,'(3x,"Release.Version :",1x,i0,".",i0)') self%Release, self%Version
    ! Surface classification name
    WRITE(*,'(3x,"Classification_Name :",1x,a)') TRIM(self%Classification_Name)
    ! Dimensions
    WRITE(*,'(3x,"n_Angles        :",1x,i0)') self%n_Angles
    WRITE(*,'(3x,"n_Frequencies   :",1x,i0)') self%n_Frequencies
    WRITE(*,'(3x,"n_Grain_Sizes   :",1x,i0)') self%n_Grain_Sizes
    WRITE(*,'(3x,"n_Depths        :",1x,i0)') self%n_Depths
    WRITE(*,'(3x,"n_Densities     :",1x,i0)') self%n_Densities
    IF ( .NOT. VISsnowCoeff_Associated(self) ) RETURN
    ! Dimension arrays
    WRITE(*,'(3x,"Angle      :")')
    WRITE(*,'(5(1x,es22.15,:))') self%Angle
    WRITE(*,'(3x,"Frequency  :")')
    WRITE(*,'(5(1x,es22.15,:))') self%Frequency
    WRITE(*,'(3x,"Grain_Size :")')
    WRITE(*,'(5(1x,es22.15,:))') self%Grain_Size
    WRITE(*,'(3x,"Depth :")')
    WRITE(*,'(5(1x,es22.15,:))') self%Depth
    WRITE(*,'(3x,"Density :")')
    WRITE(*,'(5(1x,es22.15,:))') self%Density
    ! Reflectance array
    WRITE(*,'(3x,"Reflectance :")')
    DO i5 = 1, self%n_Densities
      WRITE(*,'(5x,"DENSITY :",es22.15)') self%Density(i5)
      DO i4 = 1, self%n_Depths
        WRITE(*,'(5x,"DEPTH :",es22.15)') self%Depth(i4)
        DO i3 = 1, self%n_Grain_Sizes
          WRITE(*,'(5x,"GRAIN_SIZE :",es22.15)') self%Grain_Size(i3)
          DO i2 = 1, self%n_Frequencies
            WRITE(*,'(5x,"FREQUENCY  :",es22.15)') self%Frequency(i2)
            WRITE(*,'(5(1x,es22.15,:))') self%Reflectance(:,i2,i3,i4,i5)
          END DO
        END DO
      END DO
    END DO
  END SUBROUTINE VISsnowCoeff_Inspect



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_ValidRelease
!
! PURPOSE:
!       Function to check the VISsnowCoeff Release value.
!
! CALLING SEQUENCE:
!       IsValid = VISsnowCoeff_ValidRelease( VISsnowCoeff )
!
! INPUTS:
!       VISsnowCoeff:  VISsnowCoeff object for which the Release component
!                      is to be checked.
!                      UNITS:      N/A
!                      TYPE:       VISsnowCoeff_type
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

  FUNCTION VISsnowCoeff_ValidRelease( self ) RESULT( IsValid )
    ! Arguments
    TYPE(VISsnowCoeff_type), INTENT(IN) :: self
    ! Function result
    LOGICAL :: IsValid
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'VISsnowCoeff_ValidRelease'
    ! Local variables
    CHARACTER(ML) :: msg

    ! Set up
    IsValid = .TRUE.


    ! Check release is not too old
    IF ( self%Release < VISsnowCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An VISsnowCoeff data update is needed. ", &
                  &"VISsnowCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, VISsnowCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF


    ! Check release is not too new
    IF ( self%Release > VISsnowCOEFF_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("An VISsnowCoeff software update is needed. ", &
                  &"VISsnowCoeff release is ",i0,". Valid release is ",i0,"." )' ) &
                  self%Release, VISsnowCOEFF_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF

  END FUNCTION VISsnowCoeff_ValidRelease


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_Info
!
! PURPOSE:
!       Subroutine to return a string containing version and dimension
!       information about a VISsnowCoeff object.
!
! CALLING SEQUENCE:
!       CALL VISsnowCoeff_Info( VISsnowCoeff, Info )
!
! OBJECTS:
!       VISsnowCoeff:  VISsnowCoeff object about which info is required.
!                      UNITS:      N/A
!                      TYPE:       VISsnowCoeff_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Info:          String containing version and dimension information
!                      about the VISsnowCoeff object.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE VISsnowCoeff_Info( self, Info )
    ! Arguments
    TYPE(VISsnowCoeff_type),  INTENT(IN)  :: self
    CHARACTER(*),             INTENT(OUT) :: Info
    ! Parameters
    INTEGER, PARAMETER :: CARRIAGE_RETURN = 13
    INTEGER, PARAMETER :: LINEFEED = 10
    ! Local variables
    CHARACTER(2000) :: Long_String

    ! Write the required data to the local string
    WRITE( Long_String, &
           '( a,1x,"VISsnowCoeff RELEASE.VERSION: ", i2, ".", i2.2,a,3x, &
              &"CLASSIFICATION: ",a,",",2x,&
              &"N_ANGLES=",i3,2x,&
              &"N_FREQUENCIES=",i5,2x,&
              &"N_GRAIN_SIZES=",i3,2x,&
              &"N_DEPTHS=",i3,2x,&
              &"N_DENSITIES=",i3 )' ) &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           self%Release, self%Version, &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           TRIM(self%Classification_Name), &
           self%n_Angles, &
           self%n_Frequencies, &
           self%n_Grain_Sizes, &
           self%n_Depths, &
           self%n_Densities

    ! Trim the output based on the
    ! dummy argument string length
    Info = Long_String(1:MIN(LEN(Info), LEN_TRIM(Long_String)))

  END SUBROUTINE VISsnowCoeff_Info


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

!------------------------------------------------------------------------------
!
! NAME:
!       VISsnowCoeff_Equal
!
! PURPOSE:
!       Elemental function to test the equality of two VISsnowCoeff objects.
!       Used in OPERATOR(==) interface block.
!
! CALLING SEQUENCE:
!       is_equal = VISsnowCoeff_Equal( x, y )
!
!         or
!
!       IF ( x == y ) THEN
!         ...
!       END IF
!
! OBJECTS:
!       x, y:          Two VISsnowCoeff objects to be compared.
!                      UNITS:      N/A
!                      TYPE:       VISsnowCoeff_type
!                      DIMENSION:  Scalar or any rank
!                      ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       is_equal:      Logical value indicating whether the inputs are equal.
!                      UNITS:      N/A
!                      TYPE:       LOGICAL
!                      DIMENSION:  Same as inputs.
!
!------------------------------------------------------------------------------

  ELEMENTAL FUNCTION VISsnowCoeff_Equal( x, y ) RESULT( is_equal )
    TYPE(VISsnowCoeff_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal

    ! Set up
    is_equal = .FALSE.

    ! Check the object association status
    IF ( (.NOT. VISsnowCoeff_Associated(x)) .OR. &
         (.NOT. VISsnowCoeff_Associated(y))      ) RETURN

    ! Check contents
    ! ...Release/version info
    IF ( (x%Release /= y%Release) .OR. &
         (x%Version /= y%Version) ) RETURN
    ! ...Classification name
    IF ( (x%Classification_Name /= y%Classification_Name) ) RETURN
    ! ...Dimensions
    IF ( (x%n_Angles      /= y%n_Angles      ) .OR. &
         (x%n_Frequencies /= y%n_Frequencies ) .OR. &
         (x%n_Grain_Sizes /= y%n_Grain_Sizes ) .OR. &
         (x%n_Depths      /= y%n_Depths      ) .OR. &
         (x%n_Densities   /= y%n_Densities ) ) RETURN
    ! ...Arrays
    IF ( ALL(x%Angle       .EqualTo. y%Angle       )  .AND. &
         ALL(x%Frequency   .EqualTo. y%Frequency   )  .AND. &
         ALL(x%Grain_Size  .EqualTo. y%Grain_Size  )  .AND. &
         ALL(x%Depth       .EqualTo. y%Depth       )  .AND. &
         ALL(x%Density     .EqualTo. y%Density     )  .AND. &
         ALL(x%Reflectance .EqualTo. y%Reflectance ) ) &
      is_equal = .TRUE.

  END FUNCTION VISsnowCoeff_Equal

END MODULE VISsnowCoeff_Define