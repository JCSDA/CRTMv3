!
! TELSEM2Atlas_Define
!
! Module defining the TELSEM2Atlas object: the TELSEM2 (Tool to Estimate Land
! Surface Emissivities at Microwave to millimetre frequencies) monthly
! climatology atlas of microwave land surface emissivity.
!
! The atlas is held on the native TELSEM2 0.25-degree equal-area grid. For each
! populated grid cell and month it stores the mean emissivity at the seven SSM/I
! channels (19V, 19H, 22V, 37V, 37H, 85V, 85H) together with two surface-class
! indices: class1 drives the angular-correction regression and class2 drives the
! frequency interpolation above 85 GHz. Only the data required to reconstruct the
! emissivity are stored; the atlas uncertainty (std/covariance) used by RTTOV's
! optional outputs is not retained because the CRTM surface optics do not use it.
!
! All twelve monthly atlases are stored in a single object (stacked along n_Data
! with per-month offsets) because CRTM loads coefficients once at initialisation
! while the month is supplied per profile through the Geometry structure.
!
! Reference:
!   Aires, F., C. Prigent, F. Bernardo, C. Jimenez, R. Saunders, P. Brunel, 2011.
!   A Tool to Estimate Land-Surface Emissivities at Microwave frequencies
!   (TELSEM) for use in numerical weather prediction. Q.J.R. Meteorol. Soc.,
!   137, 690-699. doi:10.1002/qj.803
!

MODULE TELSEM2Atlas_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds           , ONLY: fp, Long, Double
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Parameters
  PUBLIC :: TELSEM2ATLAS_DATATYPE
  ! Datatypes
  PUBLIC :: TELSEM2Atlas_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Procedures
  PUBLIC :: TELSEM2Atlas_Associated
  PUBLIC :: TELSEM2Atlas_Destroy
  PUBLIC :: TELSEM2Atlas_Create
  PUBLIC :: TELSEM2Atlas_Inspect
  PUBLIC :: TELSEM2Atlas_ValidRelease
  PUBLIC :: TELSEM2Atlas_Info
  PUBLIC :: TELSEM2Atlas_Name


  ! ---------------------
  ! Procedure overloading
  ! ---------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE TELSEM2Atlas_Equal
  END INTERFACE OPERATOR(==)


  ! -----------------
  ! Module parameters
  ! -----------------
  CHARACTER(*), PARAMETER :: TELSEM2ATLAS_DATATYPE = 'TELSEM2Atlas'
  ! Current valid release and version
  INTEGER(Long), PARAMETER :: TELSEM2ATLAS_RELEASE = 1  ! Structure/file format
  INTEGER(Long), PARAMETER :: TELSEM2ATLAS_VERSION = 1  ! Default data version
  ! Literal constants
  REAL(Double), PARAMETER :: ZERO = 0.0_Double
  ! String lengths
  INTEGER, PARAMETER :: ML = 256 ! Message length
  INTEGER, PARAMETER :: SL =  80 ! String length


  ! ------------------------------
  ! TELSEM2Atlas data type
  ! Components are public: the atlas is a large read-only lookup table that the
  ! interpolation module accesses field-by-field; routing the multi-hundred-MB
  ! arrays through accessors would be pointlessly expensive.
  ! ------------------------------
  !:tdoc+:
  TYPE :: TELSEM2Atlas_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Datatype information
    CHARACTER(SL) :: Datatype_Name = TELSEM2ATLAS_DATATYPE
    ! Release and version information
    INTEGER(Long) :: Release = TELSEM2ATLAS_RELEASE
    INTEGER(Long) :: Version = TELSEM2ATLAS_VERSION
    ! Atlas description / dimensions
    REAL(Double)  :: Resolution       = ZERO  ! Equal-area grid resolution (deg)
    INTEGER(Long) :: n_Channels       = 0     ! Number of atlas channels (7 SSM/I)
    INTEGER(Long) :: n_Latitude_Bands = 0     ! Equal-area latitude bands (180/Resolution)
    INTEGER(Long) :: n_Cells          = 0     ! Total cells over the globe (max cell number)
    INTEGER(Long) :: n_Months         = 0     ! Number of monthly atlases (12)
    INTEGER(Long) :: n_Data           = 0     ! Total populated cells stacked over all months
    ! Per-month indexing into the stacked data arrays
    INTEGER(Long), ALLOCATABLE :: Month_Data_Count(:) ! (n_Months) populated cells in each month
    INTEGER(Long), ALLOCATABLE :: Month_Offset(:)     ! (n_Months) data elements preceding each month
    ! Stacked sparse atlas data (1:n_Data)
    INTEGER(Long), ALLOCATABLE :: Cell_Number(:)  ! (n_Data) equal-area cell number
    INTEGER(Long), ALLOCATABLE :: Class1(:)       ! (n_Data) angular-correction class
    INTEGER(Long), ALLOCATABLE :: Class2(:)       ! (n_Data) frequency-interpolation class
    REAL(Double),  ALLOCATABLE :: Emissivity(:,:) ! (n_Data x n_Channels) mean emissivity
    ! Derived equal-area grid geometry and reverse lookup (rebuilt on read)
    INTEGER(Long), ALLOCATABLE :: Cells_Per_Band(:) ! (n_Latitude_Bands)
    INTEGER(Long), ALLOCATABLE :: First_Cell(:)     ! (n_Latitude_Bands) first cell number of band
    INTEGER(Long), ALLOCATABLE :: Correspondence(:,:) ! (n_Cells x n_Months) cell number -> data index (0 if empty)
  END TYPE TELSEM2Atlas_type
  !:tdoc-:


CONTAINS


!################################################################################
!##                          ## PUBLIC PROCEDURES ##                           ##
!################################################################################

!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Associated
! PURPOSE:
!       Elemental function to test the status of the allocatable components
!       of the TELSEM2Atlas structure.
!:sdoc-:
!--------------------------------------------------------------------------------
  ELEMENTAL FUNCTION TELSEM2Atlas_Associated( self ) RESULT( Status )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: self
    LOGICAL :: Status
    Status = self%Is_Allocated
  END FUNCTION TELSEM2Atlas_Associated


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Destroy
! PURPOSE:
!       Elemental subroutine to re-initialize TELSEM2Atlas objects.
!:sdoc-:
!--------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE TELSEM2Atlas_Destroy( self )
    TYPE(TELSEM2Atlas_type), INTENT(OUT) :: self
    self%Is_Allocated = .FALSE.
    self%n_Channels       = 0
    self%n_Latitude_Bands = 0
    self%n_Cells          = 0
    self%n_Months         = 0
    self%n_Data           = 0
  END SUBROUTINE TELSEM2Atlas_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Create
! PURPOSE:
!       Elemental subroutine to create an instance of a TELSEM2Atlas object.
!
! CALLING SEQUENCE:
!       CALL TELSEM2Atlas_Create( self            , &
!                                 n_Channels      , &
!                                 n_Latitude_Bands, &
!                                 n_Cells         , &
!                                 n_Months        , &
!                                 n_Data            )
!:sdoc-:
!--------------------------------------------------------------------------------
  ELEMENTAL SUBROUTINE TELSEM2Atlas_Create( &
    self            , &
    n_Channels      , &
    n_Latitude_Bands, &
    n_Cells         , &
    n_Months        , &
    n_Data            )
    ! Arguments
    TYPE(TELSEM2Atlas_type), INTENT(OUT) :: self
    INTEGER(Long),           INTENT(IN)  :: n_Channels
    INTEGER(Long),           INTENT(IN)  :: n_Latitude_Bands
    INTEGER(Long),           INTENT(IN)  :: n_Cells
    INTEGER(Long),           INTENT(IN)  :: n_Months
    INTEGER(Long),           INTENT(IN)  :: n_Data
    ! Local variables
    INTEGER :: alloc_stat

    ! Check input
    IF ( n_Channels       < 1 .OR. &
         n_Latitude_Bands < 1 .OR. &
         n_Cells          < 1 .OR. &
         n_Months         < 1 .OR. &
         n_Data           < 1 ) RETURN

    ! Perform the allocation
    ALLOCATE( self%Month_Data_Count( n_Months ), &
              self%Month_Offset( n_Months ), &
              self%Cell_Number( n_Data ), &
              self%Class1( n_Data ), &
              self%Class2( n_Data ), &
              self%Emissivity( n_Data, n_Channels ), &
              self%Cells_Per_Band( n_Latitude_Bands ), &
              self%First_Cell( n_Latitude_Bands ), &
              self%Correspondence( n_Cells, n_Months ), &
              STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Initialise
    self%n_Channels       = n_Channels
    self%n_Latitude_Bands = n_Latitude_Bands
    self%n_Cells          = n_Cells
    self%n_Months         = n_Months
    self%n_Data           = n_Data
    self%Month_Data_Count = 0
    self%Month_Offset     = 0
    self%Cell_Number      = 0
    self%Class1           = 0
    self%Class2           = 0
    self%Emissivity       = ZERO
    self%Cells_Per_Band   = 0
    self%First_Cell       = 0
    self%Correspondence   = 0

    ! Set allocation indicator
    self%Is_Allocated = .TRUE.
  END SUBROUTINE TELSEM2Atlas_Create


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Inspect
! PURPOSE:
!       Subroutine to print the contents of a TELSEM2Atlas object to stdout.
!:sdoc-:
!--------------------------------------------------------------------------------
  SUBROUTINE TELSEM2Atlas_Inspect( self )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: self
    INTEGER :: m
    WRITE(*,'(1x,"TELSEM2Atlas OBJECT")')
    ! Release/version
    WRITE(*,'(3x,"Release           :",1x,i0)') self%Release
    WRITE(*,'(3x,"Version           :",1x,i0)') self%Version
    ! Dimensions
    WRITE(*,'(3x,"Resolution (deg)  :",1x,es13.6)') self%Resolution
    WRITE(*,'(3x,"n_Channels        :",1x,i0)') self%n_Channels
    WRITE(*,'(3x,"n_Latitude_Bands  :",1x,i0)') self%n_Latitude_Bands
    WRITE(*,'(3x,"n_Cells           :",1x,i0)') self%n_Cells
    WRITE(*,'(3x,"n_Months          :",1x,i0)') self%n_Months
    WRITE(*,'(3x,"n_Data            :",1x,i0)') self%n_Data
    IF ( .NOT. TELSEM2Atlas_Associated(self) ) RETURN
    WRITE(*,'(3x,"Populated cells per month:")')
    DO m = 1, self%n_Months
      WRITE(*,'(5x,"Month ",i2.2," :",1x,i0)') m, self%Month_Data_Count(m)
    END DO
  END SUBROUTINE TELSEM2Atlas_Inspect


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_ValidRelease
! PURPOSE:
!       Function to check the TELSEM2Atlas Release value is valid.
!:sdoc-:
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_ValidRelease( self ) RESULT( IsValid )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: self
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'TELSEM2Atlas_ValidRelease'
    CHARACTER(ML) :: msg

    IsValid = .TRUE.
    ! Release in the future?
    IF ( self%Release > TELSEM2ATLAS_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("A newer release is needed to read this data. ", &
                  &"Data Release=",i0,", valid release=",i0,"." )' ) &
                  self%Release, TELSEM2ATLAS_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF
    ! Release too old?
    IF ( self%Release < TELSEM2ATLAS_RELEASE ) THEN
      IsValid = .FALSE.
      WRITE( msg,'("This data is for an old release. ", &
                  &"Data Release=",i0,", valid release=",i0,"." )' ) &
                  self%Release, TELSEM2ATLAS_RELEASE
      CALL Display_Message( ROUTINE_NAME, msg, INFORMATION ); RETURN
    END IF
  END FUNCTION TELSEM2Atlas_ValidRelease


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Info
! PURPOSE:
!       Subroutine to return a string containing version and dimension
!       information about a TELSEM2Atlas object.
!:sdoc-:
!--------------------------------------------------------------------------------
  SUBROUTINE TELSEM2Atlas_Info( self, Info )
    TYPE(TELSEM2Atlas_type), INTENT(IN)  :: self
    CHARACTER(*),            INTENT(OUT) :: Info
    CHARACTER(2000) :: long_string
    WRITE( long_string, &
           '(a,1x,"TELSEM2Atlas RELEASE.VERSION: ",i2,".",i2.2,2x, &
           &"N_CHANNELS=",i0,2x,"N_MONTHS=",i0,2x,"N_DATA=",i0 )' ) &
           ACHAR(13)//ACHAR(10), &
           self%Release, self%Version, &
           self%n_Channels, self%n_Months, self%n_Data
    Info = long_string(1:MIN(LEN(Info), LEN_TRIM(long_string)))
  END SUBROUTINE TELSEM2Atlas_Info


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       TELSEM2Atlas_Name
! PURPOSE:
!       Function to return the data type name string.
!:sdoc-:
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_Name( self ) RESULT( Name )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: self
    CHARACTER(LEN(self%Datatype_Name)) :: Name
    Name = self%Datatype_Name
  END FUNCTION TELSEM2Atlas_Name


!################################################################################
!##                          ## PRIVATE PROCEDURES ##                          ##
!################################################################################

  ELEMENTAL FUNCTION TELSEM2Atlas_Equal( x, y ) RESULT( is_equal )
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal

    is_equal = .FALSE.
    ! Check allocation status
    IF ( (x%Is_Allocated .NEQV. y%Is_Allocated) ) RETURN
    ! Check scalars
    IF ( (x%Release /= y%Release) .OR. &
         (x%Version /= y%Version) ) RETURN
    IF ( .NOT. (x%Resolution .EqualTo. y%Resolution) ) RETURN
    ! Check dimensions
    IF ( (x%n_Channels       /= y%n_Channels      ) .OR. &
         (x%n_Latitude_Bands /= y%n_Latitude_Bands) .OR. &
         (x%n_Cells          /= y%n_Cells         ) .OR. &
         (x%n_Months         /= y%n_Months        ) .OR. &
         (x%n_Data           /= y%n_Data          ) ) RETURN
    ! Check the data (both unallocated counts as equal here)
    IF ( TELSEM2Atlas_Associated(x) .AND. TELSEM2Atlas_Associated(y) ) THEN
      IF ( ALL(x%Cell_Number == y%Cell_Number) .AND. &
           ALL(x%Class1      == y%Class1     ) .AND. &
           ALL(x%Class2      == y%Class2     ) .AND. &
           ALL(x%Emissivity .EqualTo. y%Emissivity) ) is_equal = .TRUE.
    ELSE
      is_equal = .TRUE.
    END IF
  END FUNCTION TELSEM2Atlas_Equal

END MODULE TELSEM2Atlas_Define
