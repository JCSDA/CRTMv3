!
! OP_Input_Define
!
! Module containing the structure definition and associated routines
! for CRTM optional inputs specific to user-defined aerosol/cloud/total optical profiles
!
!
! CREATION HISTORY:
!       Written by:     Cheng Dang, Jan, 2024
!                       dangch@ucar.edu
!

MODULE OP_Input_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds           , ONLY: fp, Long, Double
  USE Message_Handler      , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE Compare_Float_Numbers, ONLY: OPERATOR(.EqualTo.)
  USE File_Utility         , ONLY: File_Open, File_Exists
  USE netcdf
  ! Disable all implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Datatypes
  PUBLIC :: OP_Input_type
  ! Operators
  PUBLIC :: OPERATOR(==)
  ! Procedures
  PUBLIC :: OP_Input_IsValid
  PUBLIC :: OP_Input_Inspect
  PUBLIC :: OP_Input_Create
  PUBLIC :: OP_Input_Associated
  PUBLIC :: OP_Input_InquireFile
  PUBLIC :: OP_Input_ReadFile
  PUBLIC :: OP_Input_WriteFile

  ! -------------------
  ! Procedure overloads
  ! -------------------
  INTERFACE OPERATOR(==)
    MODULE PROCEDURE OP_Input_Equal
  END INTERFACE OPERATOR(==)

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Release and version
  INTEGER, PARAMETER :: OP_INPUT_RELEASE = 1  ! This determines structure and file formats.
  ! Close status for write errors
  CHARACTER(*), PARAMETER :: WRITE_ERROR_STATUS = 'DELETE'
  ! Literal constants
  REAL(Double), PARAMETER :: ZERO = 0.0_Double
  ! Message length
  INTEGER,  PARAMETER :: ML = 256

  ! NetCDF attributes
  ! Global attribute names. Case sensitive
  CHARACTER(*), PARAMETER :: RELEASE_GATTNAME     = 'Release'

  ! Dimension names
  CHARACTER(*), PARAMETER :: CHANNEL_DIMNAME    = 'n_Channels'
  CHARACTER(*), PARAMETER :: LAYER_DIMNAME      = 'n_Layers'
  CHARACTER(*), PARAMETER :: LEGENDRE_DIMNAME   = 'n_Legendre_Terms'
  CHARACTER(*), PARAMETER :: PHASE_DIMNAME      = 'n_Phase_Elements'

  ! Variable names
  CHARACTER(*), PARAMETER :: TAU_VARNAME        = 'tau'
  CHARACTER(*), PARAMETER :: BS_VARNAME         = 'bs'
  CHARACTER(*), PARAMETER :: PCOEFF_VARNAME     = 'pcoeff'
  CHARACTER(*), PARAMETER :: KB_VARNAME         = 'kb'

  ! Variable description attribute.
  CHARACTER(*), PARAMETER :: DESCRIPTION_ATTNAME = 'description'
  CHARACTER(*), PARAMETER :: TAU_DESCRIPTION     = 'Layer optical depth'
  CHARACTER(*), PARAMETER :: BS_DESCRIPTION      = 'Layer volume scattering coefficient'
  CHARACTER(*), PARAMETER :: PCOEFF_DESCRIPTION  = 'Layer phase function coefficients for scatters'
  CHARACTER(*), PARAMETER :: KB_DESCRIPTION      = 'Layer backward scattering coefficient'


  ! Variable units attribute.
  CHARACTER(*), PARAMETER :: UNITS_ATTNAME = 'units'
  CHARACTER(*), PARAMETER :: TAU_UNITS     = 'unit'
  CHARACTER(*), PARAMETER :: BS_UNITS      = 'unit'
  CHARACTER(*), PARAMETER :: PCOEFF_UNITS  = 'unit'
  CHARACTER(*), PARAMETER :: KB_UNITS      = 'Metres squared per kilogram (m^2.kg^-1)'

  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'
  REAL(Double), PARAMETER :: FILL_FLOAT = -999.0_fp

  ! Variable types
  INTEGER, PARAMETER :: FLOAT_TYPE = NF90_DOUBLE

  !--------------------
  ! Structure defintion
  !--------------------
  !:tdoc+:
  TYPE :: OP_Input_type
    ! Allocation indicator
    LOGICAL :: Is_Allocated = .FALSE.
    ! Release and version information
    INTEGER(Long) :: Release = OP_INPUT_RELEASE

    ! Dimensions
    INTEGER :: n_Channels         = 0  ! K dimension
    INTEGER :: n_Layers           = 0  ! L dimension
    INTEGER :: n_Phase_Elements   = 0  ! Ip dimension
    INTEGER :: n_Legendre_Terms   = 0  ! Il dimension

    ! Scalar components
    ! LOGICAL  :: Include_Scattering = .TRUE.
    ! INTEGER  :: lOffset = 0    ! Start position in array for Legendre coefficients
    ! REAL(fp) :: Scattering_Optical_Depth = ZERO
    ! REAL(fp) :: depolarization = 0.0279_fp

    REAL(fp), ALLOCATABLE :: tau(:,:)         ! K * L
    REAL(fp), ALLOCATABLE :: bs(:,:)          ! K * L
    REAL(fp), ALLOCATABLE :: kb(:,:)          ! K * L
    REAL(fp), ALLOCATABLE :: pcoeff(:,:,:,:)  ! K * L * Il * Ip

  END TYPE OP_Input_type
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
!   OP_Input_Associated
!
! PURPOSE:
!   Elemental function to test the status of the allocatable components
!   of a CRTM RTSolution object.
!
! CALLING SEQUENCE:
!   Status = OP_Input_Associated( OP )
!
! OBJECTS:
!   RTSolution:   OP structure which is to have its member's
!                 status tested.
!                 UNITS:      N/A
!                 TYPE:       OP
!                 DIMENSION:  Scalar or any rank
!                 ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!   Status:       The return value is a logical value indicating the
!                 status of the RTSolution members.
!                   .TRUE.  - if the array components are allocated.
!                   .FALSE. - if the array components are not allocated.
!                 UNITS:      N/A
!                 TYPE:       LOGICAL
!                 DIMENSION:  Same as input RTSolution argument
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL FUNCTION OP_Input_Associated( op ) RESULT( Status )
    TYPE(OP_Input_type), INTENT(IN) :: op
    LOGICAL :: Status
    Status = op%Is_Allocated
  END FUNCTION OP_Input_Associated

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       OP_Input_IsValid
!
! PURPOSE:
!       Non-pure function to perform some simple validity checks on a
!       OP_Input object.
!
!       If invalid data is found, a message is printed to stdout.
!
! CALLING SEQUENCE:
!       result = OP_Input_IsValid( OP )
!
!         or
!
!       IF ( OP_Input_IsValid( OP ) ) THEN....
!
! OBJECTS:
!       op:        OP_Input object which is to have its
!                  contents checked.
!                  UNITS:      N/A
!                  TYPE:       OP_Input_type
!                  DIMENSION:  Scalar
!                  ATTRIBUTES: INTENT(IN)
!
! FUNCTION RESULT:
!       result:    Logical variable indicating whether or not the input
!                  passed the check.
!                  If == .FALSE., object is unused or contains
!                                 invalid data.
!                     == .TRUE.,  object can be used.
!                  UNITS:      N/A
!                  TYPE:       LOGICAL
!                  DIMENSION:  Scalar
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION OP_Input_IsValid( op ) RESULT( IsValid )
    TYPE(OP_Input_type), INTENT(IN) :: op
    LOGICAL :: IsValid
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_IsValid'
    CHARACTER(ML) :: msg

    ! Setup
    IsValid = .TRUE.

    ! CD: placeholder for now

  END FUNCTION OP_Input_IsValid


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       OP_Input_Inspect
!
! PURPOSE:
!       Subroutine to print the contents of an Zeeman_Input object to stdout.
!
! CALLING SEQUENCE:
!       CALL OP_Input_Inspect( op )
!
! INPUTS:
!       z:             Zeeman_Input object to display.
!                      UNITS:      N/A
!                      TYPE:       OP_Input_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE OP_Input_Inspect(op)
    ! Arguments
    TYPE(OP_Input_type), INTENT(IN) :: op
    !INTEGER,    OPTIONAL, INTENT(IN) :: Unit
    ! Local variables
    INTEGER :: fid
    CHARACTER(len=*), PARAMETER :: fmt64 = '(3x,a,es22.15)'  ! print in 64-bit precision
    CHARACTER(len=*), PARAMETER :: fmt32 = '(3x,a,es13.6)'   ! print in 32-bit precision
    CHARACTER(len=*), PARAMETER :: fmt = fmt64               ! choose 64-bit precision

    WRITE(*,'(1x,"OP_Input OBJECT")')
    WRITE(*,'(3x,"n_Layers             :",i0)') op%n_Layers
    WRITE(*,'(3x,"n_Phase_Elements     :",i0)') op%n_Phase_Elements
    WRITE(*,'(3x,"n_Legendre_Terms     :",i0)') op%n_Legendre_Terms
    WRITE(fid,fmt) "Optical_Depth          : ", op%tau
    WRITE(fid,fmt) "Backscat_Coefficient   : ", op%kb
    WRITE(fid,fmt) "Scattering coefficient : ", op%bs
    WRITE(fid,fmt) "Phase_Coefficient      : ", op%pcoeff

  END SUBROUTINE OP_Input_Inspect


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       OP_Input_WriteFile
!
! PURPOSE:
!       Function to write OP_Input object files in netCDF format.
!
! CALLING SEQUENCE:
!       Error_Status = OP_Input_WriteFile( &
!                        Filename, &
!                        OP, &
!                        Quiet   = Quiet )
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the data write was successful
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
! INPUTS:

!       Filename:       Character string specifying the name of the
!                       OP_Input data file to write.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       OP:            Object containing the OP_Input data.
!                       UNITS:      N/A
!                       TYPE:       TYPE(OP_Input_type)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION OP_Input_WriteFile( &
    Filename     , &  ! Input
    OP           , &  ! Input
    Quiet        ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN) :: Filename
    TYPE(OP_Input_type)    , INTENT(IN) :: OP
    LOGICAL,       OPTIONAL, INTENT(IN) :: Quiet
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_WriteFile(netCDF)'
    ! Local variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    LOGICAL :: Noisy
    INTEGER :: NF90_Status
    INTEGER :: FileId
    INTEGER :: VarId

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    ! ...Check structure pointer association status
    IF ( .NOT. OP_Input_Associated( OP ) ) THEN
      msg = 'OP_Input structure is empty. Nothing to do!'
      CALL Write_CleanUp(); RETURN
    END IF
    ! ...Check if release is valid
    ! IF ( .NOT. OP_Input_ValidRelease( OP ) ) THEN
    !   msg = 'OP_Input Release check failed.'
    !   CALL Write_Cleanup(); RETURN
    ! END IF
    ! ...Check Quiet argument
    Noisy = .TRUE.
    IF ( PRESENT(Quiet) ) Noisy = .NOT. Quiet

    ! Create the output file
    err_stat = CreateFile( &
                 Filename                , &  ! Input
                 OP%n_Channels           , &  ! Input
                 OP%n_Layers             , &  ! Input
                 OP%n_Phase_Elements     , &  ! Input
                 OP%n_Legendre_Terms     , &  ! Input
                 OP%Release              , &  ! Input
                 FileId                   )    ! Output
    IF ( err_stat /= SUCCESS ) THEN
       msg = 'Error creating output file '//TRIM(Filename)
       CALL Write_Cleanup(); RETURN
    END IF

    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Write the data items
    ! ...tau variable
    NF90_Status = NF90_INQ_VARID( FileId,TAU_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//TAU_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarID,OP%tau )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//TAU_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...bs variable
    NF90_Status = NF90_INQ_VARID( FileId,BS_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//BS_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarID,OP%bs )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//BS_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...kb variable
    NF90_Status = NF90_INQ_VARID( FileId,KB_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//KB_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarID,OP%kb )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//KB_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...pcoeff variable
    NF90_Status = NF90_INQ_VARID( FileId,PCOEFF_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//PCOEFF_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarID,OP%pcoeff )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//PCOEFF_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId )
    Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing output file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF

    ! Output an info message
    IF ( Noisy ) THEN
      CALL OP_Input_Info( OP, msg )
      CALL Display_Message( ROUTINE_NAME, 'FILE: '//TRIM(Filename)//'; '//TRIM(msg), INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Write_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing output file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Write_CleanUp

  END FUNCTION OP_Input_WriteFile


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!   OP_Input_ReadFile
!
! PURPOSE:
!   Function to read OP_Input object files.
!
! CALLING SEQUENCE:
!   Error_Status = OP_Input_ReadFile( Filename      , &
!                                      OP           , &
!                                      noisy     )
!
! INPUTS:
!   Filename:     Character string specifying the name of an
!                 RTSolution format data file to read.
!                 UNITS:      N/A
!                 TYPE:       CHARACTER(*)
!                 DIMENSION:  Scalar
!                 ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!   OP:          OP_Input object array containing the OP_Input
!                 data.
!                 UNITS:      N/A
!                 TYPE:       OP_Input
!                 ATTRIBUTES: INTENT(OUT), ALLOCATABLE
!
! OPTIONAL INPUTS:
!   noisy:        Set this logical argument to suppress INFORMATION
!                 messages being printed to stdout
!                 If == .TRUE., INFORMATION messages are OUTPUT [DEFAULT].
!                    == .FALSE.,INFORMATION messages are SUPPRESSED.
!                 If not specified, default is .TRUE.
!                 UNITS:      N/A
!                 TYPE:       LOGICAL
!                 DIMENSION:  Scalar
!                 ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!   Error_Status: The return value is an integer defining the error status.
!                 The error codes are defined in the Message_Handler module.
!                 If == SUCCESS, the file read was successful
!                    == FAILURE, an unrecoverable error occurred.
!                 UNITS:      N/A
!                 TYPE:       INTEGER
!                 DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION OP_Input_ReadFile( &
    Filename   , &  ! Input
    OP        , &  ! Output
    Quiet      ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),                     INTENT(IN)  :: Filename
    TYPE(OP_Input_type), ALLOCATABLE, INTENT(OUT) :: OP
    LOGICAL,             OPTIONAL,    INTENT(IN)  :: Quiet
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: io_msg
    CHARACTER(ML) :: alloc_msg
    INTEGER :: io_stat
    INTEGER :: alloc_stat
    INTEGER :: fid
    INTEGER :: l, m, s, c
    LOGICAL :: Close_File
    INTEGER :: NF90_Status, FileId, VarId, Allocate_Status
    INTEGER :: n_Channels
    INTEGER :: n_Layers
    INTEGER :: n_Phase_Elements
    INTEGER :: n_Legendre_Terms
    INTEGER :: Release
    LOGICAL :: noisy
    REAL(fp), ALLOCATABLE :: tau(:,:), bs(:,:), kb(:,:),pcoeff(:,:,:,:)


    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    ! ...Check that the file exists
    IF ( .NOT. File_Exists(Filename) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet

    ! Inquire the file to get the dimensions
    err_stat = OP_Input_InquireFile( &
                 Filename                             , &
                 n_Channels        = n_Channels       , &
                 n_Layers          = n_Layers         , &
                 n_Phase_Elements  = n_Phase_Elements , &
                 n_Legendre_Terms  = n_Legendre_Terms , &
                 Release           = Release            )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error obtaining OP_Input dimensions from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF

    ! Perform the allocations for local variables
    ALLOCATE(tau( n_Channels, n_Layers ), &
             bs(  n_Channels, n_Layers ), &
             kb(  n_Channels, n_Layers ), &
             pcoeff( n_Channels         , &
                     n_Layers           , &
                     n_Phase_Elements   , &
                     n_Legendre_Terms   ), &
             STAT = alloc_stat )
    IF ( alloc_stat /= 0 ) RETURN

    ! Open the file for reading
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read access - '//&
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! ...Close the file if any error from here on
    Close_File = .TRUE.


    ! Read the OP Input data
    ! ...tau variable
    NF90_Status = NF90_INQ_VARID( FileId,TAU_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//TAU_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarID,tau )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//TAU_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...bs variable
    NF90_Status = NF90_INQ_VARID( FileId,BS_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//BS_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarID,bs )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//BS_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...kb variable
    NF90_Status = NF90_INQ_VARID( FileId,KB_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//KB_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarID,kb )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//KB_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...pcoeff variable
    NF90_Status = NF90_INQ_VARID( FileId,PCOEFF_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//PCOEFF_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarID,pcoeff )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//PCOEFF_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Assign variables
    ALLOCATE(OP)
    OP%n_Channels = n_Channels
    OP%n_Layers   = n_Layers
    OP%n_Phase_Elements = n_Phase_Elements
    OP%n_Legendre_Terms = n_Legendre_Terms
    OP%tau = tau
    OP%bs  = bs
    OP%kb  = kb
    OP%pcoeff = pcoeff

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing output file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Output an info message
    IF ( noisy ) THEN
      CALL OP_Input_Info( OP, msg )
      CALL Display_Message( ROUTINE_NAME, 'FILE: '//TRIM(Filename)//'; '//TRIM(msg), INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Read_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup- '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      CALL OP_Input_Destroy( OP )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION OP_Input_ReadFile


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       OP_Input_InquireFile
!
! PURPOSE:
!       Subroutine to print the contents of an SSU_Input object to stdout.
!
! CALLING SEQUENCE:
!       CALL   OP_Input_InquireFile( &
!     Filename         , &  ! Input
!     n_Channels       , &  ! Output
!     n_Layers         , &  ! Output
!     n_Phase_Elements , &  ! Output
!     n_Legendre_Terms , &  ! Output
!     Release          ) &  ! Output
!
! INPUTS:
!       ssu:           SSU_Input object to display.
!                      UNITS:      N/A
!                      TYPE:       SSU_Input_type
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  FUNCTION OP_Input_InquireFile( &
    Filename         , &  ! Input
    n_Channels       , &  ! Output
    n_Layers         , &  ! Output
    n_Phase_Elements , &  ! Output
    n_Legendre_Terms , &  ! Output
    Release          ) &  ! Output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),  INTENT(IN)  :: Filename
    INTEGER     ,  INTENT(OUT) :: n_Channels
    INTEGER     ,  INTENT(OUT) :: n_Layers
    INTEGER     ,  INTENT(OUT) :: n_Phase_Elements
    INTEGER     ,  INTENT(OUT) :: n_Legendre_Terms
    INTEGER     ,  INTENT(OUT) :: Release
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: GAttName
    INTEGER :: io_stat
    INTEGER :: fid
    LOGICAL :: Close_File
    INTEGER :: NF90_Status, FileId, VarId, DimId

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.

    ! Open the file
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read access - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF

    ! ...Close the file if any error from here on
    Close_File = .TRUE.
    ! Get the dimensions
    ! ...n_Channels dimension
    NF90_Status = NF90_INQ_DIMID( FileId,CHANNEL_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//CHANNEL_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=n_Channels )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//CHANNEL_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    ! ...n_Layers dimension
    NF90_Status = NF90_INQ_DIMID( FileId,LAYER_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//LAYER_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=n_Layers )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//LAYER_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    ! ...n_Legendre_Terms dimension
    NF90_Status = NF90_INQ_DIMID( FileId, LEGENDRE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//LEGENDRE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=n_Legendre_Terms )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//LEGENDRE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    ! ...n_Phase_Elements dimension
    NF90_Status = NF90_INQ_DIMID( FileId,PHASE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//PHASE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=n_Phase_Elements )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//PHASE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_CleanUp(); RETURN
    END IF

    ! Get the global attributes
    GAttName = RELEASE_GATTNAME
    NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),Release )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      CALL ReadGAtts_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId )
    Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing input file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE ReadGAtts_CleanUp()
      err_stat = FAILURE
      msg = 'Error reading '//TRIM(GAttName)//' attribute from '//TRIM(Filename)//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ) )
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE ReadGAtts_CleanUp

    SUBROUTINE Inquire_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup.'
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION  OP_Input_InquireFile

!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!   OP_Input_Destroy
!
! PURPOSE:
!   Elemental subroutine to re-initialize CRTM RTSolution objects.
!
! CALLING SEQUENCE:
!   CALL OP_Input_Destroy( OP )
!
! OBJECTS:
!   RTSolution:   Re-initialized RTSolution structure.
!                 UNITS:      N/A
!                 TYPE:       CRTM_RTSolution_type
!                 DIMENSION:  Scalar OR any rank
!                 ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE OP_Input_Destroy( op )
    TYPE(OP_Input_type), INTENT(OUT) :: op
    op%Is_Allocated = .FALSE.
    op%n_Channels = 0
  END SUBROUTINE OP_Input_Destroy


!--------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       OP_Input_Info
!
! PURPOSE:
!       Subroutine to return a string containing version and dimension
!       information about a AerosolCoeff object.
!
! CALLING SEQUENCE:
!       CALL OP_Input_Info( OP, Info )
!
! INPUTS:
!       OP:            OP_Input object about which info is required.
!                      UNITS:      N/A
!                      TYPE:       TYPE(OP_Input)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Info:          String containing version and dimension information
!                      about the passed AerosolCoeff object.
!                      UNITS:      N/A
!                      TYPE:       CHARACTER(*)
!                      DIMENSION:  Scalar
!                      ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!--------------------------------------------------------------------------------

  SUBROUTINE OP_Input_Info( OP, Info )
    ! Arguments
    TYPE(OP_Input_type), INTENT(IN)  :: OP
    CHARACTER(*),        INTENT(OUT) :: Info
    ! Parameters
    INTEGER, PARAMETER :: CARRIAGE_RETURN = 13
    INTEGER, PARAMETER :: LINEFEED = 10
    ! Local variables
    CHARACTER(2000) :: Long_String

    ! Write the required data to the local string
    WRITE( Long_String, &
           '(a,1x,"OP_Input RELEASE : ",i3, &
           &"N_CHANNELS=",i4,2x,&
           &"N_LAYERS=",i3,2x,&
           &"N_PHASE_ELEMENTS=",i2,2x,&
           &"N_LEGENDRE_TERMS=",i2 )' ) &
           ACHAR(CARRIAGE_RETURN)//ACHAR(LINEFEED), &
           OP%Release          , &
           OP%n_Channels       , &
           OP%n_Layers         , &
           OP%n_Phase_Elements , &
           OP%n_Legendre_Terms

    ! Trim the output based on the
    ! dummy argument string length
    Info = Long_String(1:MIN(LEN(Info), LEN_TRIM(Long_String)))

  END SUBROUTINE OP_Input_Info

!--------------------------------------------------------------------------------
!
! NAME:
!       OP_Input_Create
! PURPOSE:
!       Elemental subroutine to create an instance of a OP_Input object.
!
! CALLING SEQUENCE:
!       CALL OP_Input_Create( OP             , &
!                              n_Channels      , &
!                              n_Layers        , &
!                              n_Phase_Elements, &
!                              n_Legendre_Terms  )
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE OP_Input_Create( &
      OP             , &
      n_Channels      , &
      n_Layers        , &
      n_Phase_Elements, &
      n_Legendre_Terms  )
      ! Arguments
      TYPE(OP_Input_type), INTENT(OUT) :: OP
      INTEGER,              INTENT(IN)  :: n_Channels
      INTEGER,              INTENT(IN)  :: n_Layers
      INTEGER,              INTENT(IN)  :: n_Phase_Elements
      INTEGER,              INTENT(IN)  :: n_Legendre_Terms
      ! Local parameters
      CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_Create'
      ! Local variables
      INTEGER :: alloc_stat

      ! Check input
      IF ( n_Channels       < 1 .OR. &
           n_Layers         < 1 .OR. &
           n_Legendre_Terms < 0 .OR. &
           n_Phase_Elements < 1      ) RETURN

      ! Perform the allocations.
      ALLOCATE(OP%tau( n_Channels, n_Layers ), &
               OP%bs(  n_Channels, n_Layers ), &
               OP%kb(  n_Channels, n_Layers ), &
               OP%pcoeff( n_Channels         , &
                           n_Layers           , &
                           n_Phase_Elements   , &
                           n_Legendre_Terms   ), &
               STAT = alloc_stat )
      IF ( alloc_stat /= 0 ) RETURN

      ! Initialise
      ! ...Dimensions
      OP%n_Channels          = n_Channels
      OP%n_Layers            = n_Layers
      OP%n_Phase_Elements    = n_Phase_Elements
      OP%n_Legendre_Terms    = n_Legendre_Terms

      ! ...Arrays
      OP%tau    = ZERO
      OP%bs     = ZERO
      OP%kb     = ZERO
      OP%pcoeff = ZERO

      ! Set allocationindicator
      OP%Is_Allocated = .TRUE.

  END SUBROUTINE OP_Input_Create

!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################
  FUNCTION CreateFile( &
    Filename        , &  ! Input
    n_Channels      , &  ! Input
    n_Layers        , &  ! Input
    n_Phase_Elements, &  ! Input
    n_Legendre_Terms, &  ! Input
    Release         , &  ! Input
    FileId          ) &  ! Output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     ,           INTENT(IN)  :: n_Channels
    INTEGER     ,           INTENT(IN)  :: n_Layers
    INTEGER     ,           INTENT(IN)  :: n_Phase_Elements
    INTEGER     ,           INTENT(IN)  :: n_Legendre_Terms
    INTEGER     ,           INTENT(IN)  :: Release
    INTEGER     ,           INTENT(OUT) :: FileId
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'OP_Input_CreateFile(netCDF)'
    ! Local variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status, VarID
    INTEGER :: n_Channels_DimID
    INTEGER :: n_Layers_DimID
    INTEGER :: n_Phase_Elements_DimID
    INTEGER :: n_Legendre_Terms_DimID
    INTEGER :: Put_Status(3)

    ! Setup
    err_stat = SUCCESS
    Close_File = .FALSE.

    ! Create the data file
    NF90_Status = NF90_CREATE( Filename,NF90_CLOBBER,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error creating '//TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Define the dimensions
    ! ...Number of Channels
    NF90_Status = NF90_DEF_DIM( FileID,CHANNEL_DIMNAME,n_Channels,n_Channels_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//CHANNEL_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...Number of Layers
    NF90_Status = NF90_DEF_DIM( FileID,LAYER_DIMNAME,n_Layers,n_Layers_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//LAYER_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...Number of Phase_Elements
    NF90_Status = NF90_DEF_DIM( FileID,PHASE_DIMNAME,n_Phase_Elements,n_Phase_Elements_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//PHASE_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...Number of Legendre_Terms
    NF90_Status = NF90_DEF_DIM( FileID,LEGENDRE_DIMNAME,n_Legendre_Terms,n_Legendre_Terms_DimID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//LEGENDRE_DIMNAME//' dimension in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

    ! Write the global attributes
    NF90_Status = NF90_PUT_ATT( FileId, NF90_GLOBAL,TRIM(RELEASE_GATTNAME),Release )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error setting '//RELEASE_GATTNAME//' global attribute in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF


    ! Write the variables
    ! ...tau variable
    NF90_Status = NF90_DEF_VAR( FileID, &
      TAU_VARNAME, &
      FLOAT_TYPE, &
      dimIDs=(/n_Channels_DimID, n_Layers_DimID/), &
      varID=VarID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//TAU_VARNAME//' variable in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    Put_Status(1) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME ,TAU_DESCRIPTION)
    Put_Status(2) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME       ,TAU_UNITS)
    Put_Status(3) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME   ,FILL_FLOAT)
    IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
      msg = 'Error writing '//TAU_VARNAME//' variable attributes to '//TRIM(Filename)
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...bs variable
    NF90_Status = NF90_DEF_VAR( FileID, &
      BS_VARNAME, &
      FLOAT_TYPE, &
      dimIDs=(/n_Channels_DimID, n_Layers_DimID/), &
      varID=VarID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//BS_VARNAME//' variable in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    Put_Status(1) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME ,BS_DESCRIPTION)
    Put_Status(2) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME       ,BS_UNITS)
    Put_Status(3) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME   ,FILL_FLOAT)
    IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
      msg = 'Error writing '//BS_VARNAME//' variable attributes to '//TRIM(Filename)
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...kb variable
    NF90_Status = NF90_DEF_VAR( FileID, &
      KB_VARNAME, &
      FLOAT_TYPE, &
      dimIDs=(/n_Channels_DimID, n_Layers_DimID/), &
      varID=VarID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//KB_VARNAME//' variable in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    Put_Status(1) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME ,KB_DESCRIPTION)
    Put_Status(2) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME       ,KB_UNITS)
    Put_Status(3) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME   ,FILL_FLOAT)
    IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
      msg = 'Error writing '//KB_VARNAME//' variable attributes to '//TRIM(Filename)
      CALL Create_Cleanup(); RETURN
    END IF
    ! ...pcoeff variable
    NF90_Status = NF90_DEF_VAR( FileID, &
      PCOEFF_VARNAME, &
      FLOAT_TYPE, &
      dimIDs=(/n_Channels_DimID, n_Layers_DimID, n_Phase_Elements_DimID, n_Legendre_Terms_DimID/), &
      varID=VarID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error defining '//PCOEFF_VARNAME//' variable in '//&
            TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF
    Put_Status(1) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME ,PCOEFF_DESCRIPTION)
    Put_Status(2) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME       ,PCOEFF_UNITS)
    Put_Status(3) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME   ,FILL_FLOAT)
    IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
      msg = 'Error writing '//PCOEFF_VARNAME//' variable attributes to '//TRIM(Filename)
      CALL Create_Cleanup(); RETURN
    END IF

    ! Take netCDF file out of define mode
    NF90_Status = NF90_ENDDEF( FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error taking file '//TRIM(Filename)// &
            ' out of define mode - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Create_Cleanup(); RETURN
    END IF

  CONTAINS

    SUBROUTINE Create_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileID )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup - '//&
                TRIM(NF90_STRERROR( NF90_Status ))
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Create_CleanUp

  END FUNCTION CreateFile


  ! ---------------------------------------------------------------------------

  ELEMENTAL FUNCTION OP_Input_Equal(x, y) RESULT(is_equal)
    TYPE(OP_Input_type), INTENT(IN) :: x, y
    LOGICAL :: is_equal

    ! Setup
    is_equal = .FALSE.

    is_equal = (x%n_Layers              == y%n_Layers             ) .AND. &
               (x%n_Phase_Elements      == y%n_Phase_Elements     ) .AND. &
               (x%n_Legendre_Terms      == y%n_Legendre_Terms     ) !.AND. &
            ! ALL(x%Optical_Depth         .EqualTo. y%Optical_Depth        ) .AND. &
            ! ALL(x%Single_Scatter_Albedo .EqualTo. y%Single_Scatter_Albedo) .AND. &
            ! ALL(x%Asymmetry_Factor      .EqualTo. y%Asymmetry_Factor     ) .AND. &
            ! ALL(x%Backscat_Coefficient  .EqualTo. y%Backscat_Coefficient ) .AND. &
            ! ALL(x%Delta_Truncation      .EqualTo. y%Delta_Truncation     ) .AND. &
            ! ALL(x%Phase_Coefficient     .EqualTo. y%Phase_Coefficient    )
  END FUNCTION OP_Input_Equal

END MODULE OP_Input_Define
