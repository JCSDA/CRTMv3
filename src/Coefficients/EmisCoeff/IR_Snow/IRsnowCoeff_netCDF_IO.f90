!
! IRsnowCoeff_netCDF_IO
!
! Module containing routines to read and write IRsnowCoeff netCDF
! format files.
!
!
! CREATION HISTORY:
!
!       Written by:    Cheng Dang, 23-May-2022
!                      dangch@ucar.edu

MODULE IRsnowCoeff_netCDF_IO

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds         , ONLY: fp, Double, Long
  USE Message_Handler    , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility       , ONLY: File_Exists
  USE String_Utility     , ONLY: StrClean
  USE IRsnowCoeff_Define , ONLY: IRsnowCoeff_type, &
                                 IRsnowCoeff_Associated, &
                                 IRsnowCoeff_Create, &
                                 IRsnowCoeff_Inspect, &
                                 IRsnowCoeff_Destroy, &
                                 IRsnowCoeff_ValidRelease, &
                                 IRsnowCoeff_Info
  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: IRsnowCoeff_netCDF_InquireFile
  PUBLIC :: IRsnowCoeff_netCDF_ReadFile
  PUBLIC :: IRsnowCoeff_netCDF_WriteFile

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Default msg string length
  INTEGER, PARAMETER :: ML = 1024
  ! Literal constants
  REAL(fp), PARAMETER :: FILL_FLOAT = -999.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp

  ! Global attribute names. Case sensitive
  CHARACTER(*), PARAMETER :: RELEASE_GATTNAME     = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_GATTNAME     = 'Version'
  CHARACTER(*), PARAMETER :: DATA_SOURCE_GATTNAME = 'Data_Source'
  CHARACTER(*), PARAMETER :: TITLE_GATTNAME       = 'Title'
  CHARACTER(*), PARAMETER :: HISTORY_GATTNAME     = 'History'
  CHARACTER(*), PARAMETER :: COMMENT_GATTNAME     = 'Comment'
  CHARACTER(*), PARAMETER :: CLASSIFICATION_NAME_GATTNAME = 'Classification_Name'

  ! Dimension names
  CHARACTER(*), PARAMETER :: TNSL_DIMNAME         = 'String_Length'
  CHARACTER(*), PARAMETER :: FREQUENCY_DIMNAME    = 'n_Frequencies'
  CHARACTER(*), PARAMETER :: ANGLE_DIMNAME        = 'n_Angles'
  CHARACTER(*), PARAMETER :: GRAINSIZE_DIMNAME    = 'n_Grain_Sizes'
  CHARACTER(*), PARAMETER :: TEMPERATURE_DIMNAME  = 'n_Temperature'

  ! Variable names
  CHARACTER(*), PARAMETER :: ANGLE_VARNAME        = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_VARNAME    = 'Frequency'
  CHARACTER(*), PARAMETER :: GRAINSIZE_VARNAME    = 'Grain_Size'
  CHARACTER(*), PARAMETER :: TEMPERATURE_VARNAME  = 'Temperature'
  CHARACTER(*), PARAMETER :: EMISSIVITY_VARNAME   = 'Emissivity'

  ! Variable long name attribute.
  CHARACTER(*), PARAMETER :: LONGNAME_ATTNAME = 'long_name'
  CHARACTER(*), PARAMETER :: ANGLE_LONGNAME        = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_LONGNAME    = 'Frequency'
  CHARACTER(*), PARAMETER :: GRAINSIZE_LONGNAME    = 'Grain Size'
  CHARACTER(*), PARAMETER :: TEMPERATURE_LONGNAME  = 'Temperature'
  CHARACTER(*), PARAMETER :: EMISSIVITY_LONGNAME   = 'Emissivity'

  ! Variable description attribute.
  CHARACTER(*), PARAMETER :: DESCRIPTION_ATTNAME = 'description'
  CHARACTER(*), PARAMETER :: ANGLE_DESCRIPTION        = 'Angle dimension values for emissivity data'
  CHARACTER(*), PARAMETER :: FREQUENCY_DESCRIPTION    = 'Frequency dimension values for emissivity data'
  CHARACTER(*), PARAMETER :: GRAINSIZE_DESCRIPTION    = 'Grain Size dimension values for emissivity data'
  CHARACTER(*), PARAMETER :: TEMPERATURE_DESCRIPTION  = 'Temperature dimension values for emissivity data'
  CHARACTER(*), PARAMETER :: EMISSIVITY_DESCRIPTION   = 'Spectral water surface emissivity data'

  ! Variable units attribute.
  CHARACTER(*), PARAMETER :: UNITS_ATTNAME = 'units'
  CHARACTER(*), PARAMETER :: ANGLE_UNITS        = 'degrees from vertical'
  CHARACTER(*), PARAMETER :: FREQUENCY_UNITS    = 'inverse centimeters (cm^-1)'
  CHARACTER(*), PARAMETER :: GRAINSIZE_UNITS    = 'effective radius in microns (um)'
  CHARACTER(*), PARAMETER :: TEMPERATURE_UNITS  = 'Kelvins'
  CHARACTER(*), PARAMETER :: EMISSIVITY_UNITS   = 'N/A'

  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'
  REAL(Double), PARAMETER :: ANGLE_FILLVALUE       = FILL_FLOAT
  REAL(Double), PARAMETER :: FREQUENCY_FILLVALUE   = FILL_FLOAT
  REAL(Double), PARAMETER :: GRAINSIZE_FILLVALUE   = FILL_FLOAT
  REAL(Double), PARAMETER :: TEMPERATURE_FILLVALUE = FILL_FLOAT
  REAL(Double), PARAMETER :: EMISSIVITY_FILLVALUE  = FILL_FLOAT

  ! Variable types
  INTEGER, PARAMETER :: ANGLE_TYPE       = NF90_DOUBLE
  INTEGER, PARAMETER :: FREQUENCY_TYPE   = NF90_DOUBLE
  INTEGER, PARAMETER :: GRAINSIZE_TYPE   = NF90_DOUBLE
  INTEGER, PARAMETER :: TEMPERATURE_TYPE = NF90_DOUBLE
  INTEGER, PARAMETER :: EMISSIVITY_TYPE  = NF90_DOUBLE


CONTAINS

!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################
!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       IRsnowCoeff_netCDF_InquireFile
!
! PURPOSE:
!       Function to inquire IRsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = IRsnowCoeff_netCDF_InquireFile( &
!                        Filename, &
!                        n_Angles         = n_Angles        , &
!                        n_Frequencies    = n_Frequencies   , &
!                        n_Grain_Sizes    = n_Grain_Sizes   , &
!                        n_Temperature    = n_Temperature   , &
!                        Release          = Release         , &
!                        Version          = Version         , &
!                        Title            = Title           , &
!                        History          = History         , &
!                        Comment          = Comment           )
!
! INPUTS:
!       Filename:          Character string specifying the name of a
!                          IRsnowCoeff data file to read.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
! OPTIONAL OUTPUTS:
!       n_Angles:          The number of angles in the look-up
!                          table (LUT). Must be > 0.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
!       n_Frequencies:     The number of frequencies in the LUT.
!                          Must be > 0.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
!       n_Grain_Sizes:     The number of grain size in
!                          the LUT. Must be > 0.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
!       n_Grain_Sizes:     The number of temperature in
!                          the LUT. Must be > 0.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
!       Release:           The release number of the IRsnowCoeff file.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Version:           The version number of the IRsnowCoeff file.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Title:             Character string written into the TITLE global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       History:           Character string written into the HISTORY global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Comment:           Character string written into the COMMENT global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
! FUNCTION RESULT:
!       Error_Status:      The return value is an integer defining the error status.
!                          The error codes are defined in the Message_Handler module.
!                          If == SUCCESS, the file inquire was successful
!                             == FAILURE, an unrecoverable error occurred.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION IRsnowCoeff_netCDF_InquireFile( &
    Filename     , &  ! Input
    n_Angles     , &  ! Optional output
    n_Frequencies, &  ! Optional output
    n_Grain_Sizes, &  ! Optional output
    n_Temperature, &  ! Optional output
    Release      , &  ! Optional output
    Version      , &  ! Optional output
    Title        , &  ! Optional output
    History      , &  ! Optional output
    Comment      ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Angles
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Frequencies
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Grain_Sizes
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Temperature
    INTEGER     , OPTIONAL, INTENT(OUT) :: Release
    INTEGER     , OPTIONAL, INTENT(OUT) :: Version
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: History
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_netCDF_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status
    INTEGER :: FileId
    INTEGER :: DimId
    TYPE(IRsnowCoeff_type) :: IRsnowCoeff

    ! Setup
    err_stat = SUCCESS
    Close_File = .FALSE.

    ! Open the file
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read access - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Get the dimensions
    ! ...n_Angles dimension
    NF90_Status = NF90_INQ_DIMID( FileId,ANGLE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//ANGLE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=IRsnowCoeff%n_Angles )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//ANGLE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    ! ...n_Frequencies dimension
    NF90_Status = NF90_INQ_DIMID( FileId,FREQUENCY_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//FREQUENCY_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=IRsnowCoeff%n_Frequencies )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//FREQUENCY_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    ! ...n_Grain_Sizes dimension
    NF90_Status = NF90_INQ_DIMID( FileId,GRAINSIZE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//GRAINSIZE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=IRsnowCoeff%n_Grain_Sizes )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//GRAINSIZE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    ! ...n_Temperature dimension
    NF90_Status = NF90_INQ_DIMID( FileId,TEMPERATURE_DIMNAME,DimId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring dimension ID for '//TEMPERATURE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_INQUIRE_DIMENSION( FileId,DimId,Len=IRsnowCoeff%n_Temperature )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading dimension value for '//TEMPERATURE_DIMNAME//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Get the global attributes
    err_stat = ReadGAtts( Filename, &
                          FileId  , &
                          Release = IRsnowCoeff%Release, &
                          Version = IRsnowCoeff%Version, &
                          Title   = Title  , &
                          History = History, &
                          Comment = Comment  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading global attributes from '//TRIM(Filename)
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId )
    Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing input file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Inquire_Cleanup(); RETURN
    END IF

    ! Set the return values
    IF ( PRESENT(n_Angles     ) ) n_Angles       = IRsnowCoeff%n_Angles
    IF ( PRESENT(n_Frequencies) ) n_Frequencies  = IRsnowCoeff%n_Frequencies
    IF ( PRESENT(n_Grain_Sizes) ) n_Grain_Sizes  = IRsnowCoeff%n_Grain_Sizes
    IF ( PRESENT(n_Temperature) ) n_Temperature  = IRsnowCoeff%n_Temperature
    IF ( PRESENT(Release      ) ) Release        = IRsnowCoeff%Release
    IF ( PRESENT(Version      ) ) Version        = IRsnowCoeff%Version

  CONTAINS

    SUBROUTINE Inquire_CleanUp()
      IF ( Close_File ) THEN
        NF90_Status = NF90_CLOSE( FileId )
        IF ( NF90_Status /= NF90_NOERR ) &
          msg = TRIM(msg)//'; Error closing input file during error cleanup.'
      END IF
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Inquire_CleanUp

  END FUNCTION IRsnowCoeff_netCDF_InquireFile

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       IRsnowCoeff_netCDF_WriteFile
!
! PURPOSE:
!       Function to write IRsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = IRsnowCoeff_netCDF_WriteFile( &
!                        IRsnowCoeff, &
!                        Filename, &
!                        Quiet   = Quiet  , &
!                        Title   = Title  , &
!                        History = History, &
!                        Comment = Comment, &
!                        Debug  )
!
! INPUTS:
!       IRsnowCoeff:   Object containing the IRsnow coefficient data.
!                       UNITS:      N/A
!                       TYPE:       TYPE(IRsnowCoeff_type)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       Filename:       Character string specifying the name of the
!                       IRsnowCoeff data file to write.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
! OPTIONAL INPUTS:
!       Quiet:          Set this logical argument to suppress INFORMATION
!                       messages being printed to stdout
!                       If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                          == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                       If not specified, default is .FALSE.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Title:          Character string written into the TITLE global
!                       attribute field of the IRsnowCoeff file.
!                       This argument is ignored if the netCDF argument
!                       is not supplied or set.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       History:        Character string written into the HISTORY global
!                       attribute field of the IRsnowCoeff file.
!                       This argument is ignored if the netCDF argument
!                       is not supplied or set.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Comment:        Character string written into the COMMENT global
!                       attribute field of the IRsnowCoeff file.
!                       This argument is ignored if the netCDF argument
!                       is not supplied or set.
!                       UNITS:      N/A
!                       TYPE:       CHARACTER(*)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:   The return value is an integer defining the error status.
!                       The error codes are defined in the Message_Handler module.
!                       If == SUCCESS the data write was successful
!                          == FAILURE an unrecoverable error occurred.
!                       UNITS:      N/A
!                       TYPE:       INTEGER
!                       DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION IRsnowCoeff_netCDF_WriteFile( &
    IRsnowCoeff , &  ! Input
    Filename    , &  ! Input
    Quiet       , &  ! Optional input
    Title       , &  ! Optional input
    History     , &  ! Optional input
    Comment     , &  ! Optional input
    Debug       ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(IRsnowCoeff_type), INTENT(IN) :: IRsnowCoeff
    CHARACTER(*),            INTENT(IN) :: Filename
    LOGICAL     ,  OPTIONAL, INTENT(IN) :: Quiet
    CHARACTER(*),  OPTIONAL, INTENT(IN) :: Title
    CHARACTER(*),  OPTIONAL, INTENT(IN) :: History
    CHARACTER(*),  OPTIONAL, INTENT(IN) :: Comment
    LOGICAL     ,  OPTIONAL, INTENT(IN) :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_netCDF_WriteFile'
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
    IF ( .NOT. IRsnowCoeff_Associated( IRsnowCoeff ) ) THEN
      msg = 'IRsnowCoeff structure is empty. Nothing to do!'
      CALL Write_CleanUp(); RETURN
    END IF
    ! ...Check if release is valid
    IF ( .NOT. IRsnowCoeff_ValidRelease( IRsnowCoeff ) ) THEN
      msg = 'IRsnowCoeff Release check failed.'
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Check Quiet argument
    Noisy = .TRUE.
    IF ( PRESENT(Quiet) ) Noisy = .NOT. Quiet

    ! Create the output file
    err_stat = CreateFile( &
                 Filename                                              , &  ! Input
                 IRsnowCoeff%n_Angles                                  , &  ! Input
                 IRsnowCoeff%n_Frequencies                             , &  ! Input
                 IRsnowCoeff%n_Grain_Sizes                             , &  ! Input
                 IRsnowCoeff%n_Temperature                             , &  ! Input
                 FileId                                                , &  ! Input
                 Release             = IRsnowCoeff%Release             , &  ! Optional input
                 Version             = IRsnowCoeff%Version             , &  ! Optional input
                 Classification_Name = IRsnowCoeff%Classification_Name , &  ! Optional input
                 Title               = Title                           , &  ! Optional input
                 History             = History                         , &  ! Optional input
                 Comment             = Comment                           )  ! Optional input
    IF ( err_stat /= SUCCESS ) THEN
       msg = 'Error creating output file '//TRIM(Filename)
       CALL Write_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Write the data items
    ! ...Angle variable
    NF90_Status = NF90_INQ_VARID( FileId,ANGLE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//ANGLE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,IRsnowCoeff%Angle )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//ANGLE_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Frequency variable
    NF90_Status = NF90_INQ_VARID( FileId,FREQUENCY_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//FREQUENCY_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,IRsnowCoeff%Frequency )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//FREQUENCY_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Grain Size variable
    NF90_Status = NF90_INQ_VARID( FileId,GRAINSIZE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//GRAINSIZE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,IRsnowCoeff%Grain_Size )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//GRAINSIZE_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Temperature variable
    NF90_Status = NF90_INQ_VARID( FileId,TEMPERATURE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//TEMPERATURE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,IRsnowCoeff%Temperature )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//TEMPERATURE_VARNAME//' to '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    ! ...Emissivity variable
    NF90_Status = NF90_INQ_VARID( FileId,EMISSIVITY_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//EMISSIVITY_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Write_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_PUT_VAR( FileId,VarId,IRsnowCoeff%Emissivity )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error writing '//EMISSIVITY_VARNAME//' to '//TRIM(Filename)//&
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
      CALL IRsnowCoeff_Info( IRsnowCoeff, msg )
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

  END FUNCTION IRsnowCoeff_netCDF_WriteFile

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       IRsnowCoeff_netCDF_ReadFile
!
! PURPOSE:
!       Function to read IRsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = IRsnowCoeff_netCDF_ReadFile( &
!                        IRsnowCoeff, &
!                        Filename, &
!                        Quiet   = Quiet  , &
!                        Title   = Title  , &
!                        History = History, &
!                        Comment = Comment  )
!
! INPUTS:
!       Filename:          Character string specifying the name of a
!                          IRsnowCoeff data file to read.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       IRsnowCoeff:   Object containing the IRsnow coefficient data.
!                       UNITS:      N/A
!                       TYPE:       TYPE(IRsnowCoeff_type)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUTS:
!       Quiet:          Set this logical argument to suppress INFORMATION
!                       messages being printed to stdout
!                       If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                          == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                       If not specified, default is .FALSE.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
! OPTIONAL OUTPUTS:
!       Title:             Character string written into the TITLE global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       History:           Character string written into the HISTORY global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Comment:           Character string written into the COMMENT global
!                          attribute field of the IRsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
! FUNCTION RESULT:
!       Error_Status:      The return value is an integer defining the error status.
!                          The error codes are defined in the Message_Handler module.
!                          If == SUCCESS, the file inquire was successful
!                             == FAILURE, an unrecoverable error occurred.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION IRsnowCoeff_netCDF_ReadFile( &
    IRsnowCoeff , &  ! Output
    Filename    , &  ! Input
    Quiet       , &  ! Optional input
    Title       , &  ! Optional output
    History     , &  ! Optional output
    Comment     , &  ! Optional output
    Debug       ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(IRsnowCoeff_type) , INTENT(OUT) :: IRsnowCoeff
    CHARACTER(*),            INTENT(IN)  :: Filename
    LOGICAL     ,  OPTIONAL, INTENT(IN)  :: Quiet
    CHARACTER(*),  OPTIONAL, INTENT(OUT) :: Title
    CHARACTER(*),  OPTIONAL, INTENT(OUT) :: History
    CHARACTER(*),  OPTIONAL, INTENT(OUT) :: Comment
    LOGICAL     ,  OPTIONAL, INTENT(IN)  :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_netCDF_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    LOGICAL :: Noisy
    INTEGER :: NF90_Status
    INTEGER :: FileId
    INTEGER :: n_Angles
    INTEGER :: n_Frequencies
    INTEGER :: n_Grain_Sizes
    INTEGER :: n_Temperature
    INTEGER :: VarId

    ! Set up
    err_stat = SUCCESS
    Close_File = .FALSE.
    ! ...Check that the file exists
    IF ( .NOT. File_Exists(Filename) ) THEN
      msg = 'File '//TRIM(Filename)//' not found.'
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check Quiet argument
    Noisy = .TRUE.
    IF ( PRESENT(Quiet) ) Noisy = .NOT. Quiet

    ! Inquire the file to get the dimensions
    err_stat = IRsnowCoeff_netCDF_InquireFile( &
                 Filename                       , &
                 n_Angles      = n_Angles       , &
                 n_Frequencies = n_Frequencies  , &
                 n_Grain_Sizes = n_Grain_Sizes  , &
                 n_Temperature = n_Temperature   )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error obtaining IRsnowCoeff dimensions from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF

    ! Allocate the output structure
    CALL IRsnowCoeff_Create( &
           IRsnowCoeff   , &
           n_Angles       , &
           n_Frequencies  , &
           n_Grain_Sizes  , &
           n_Temperature   )
    IF ( .NOT. IRsnowCoeff_Associated( IRsnowCoeff ) ) THEN
      msg = 'IRsnowCoeff object allocation failed.'
      CALL Read_Cleanup(); RETURN
    END IF

    ! Open the file for reading
    NF90_Status = NF90_OPEN( Filename,NF90_NOWRITE,FileId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error opening '//TRIM(Filename)//' for read access - '//&
            TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Close the file if any error from here on
    Close_File = .TRUE.

    ! Read the global attributes
    err_stat = ReadGAtts( Filename, &
                          FileID  , &
                          Release                 = IRsnowCoeff%Release             , &
                          Version                 = IRsnowCoeff%Version             , &
                          Classification_Name     = IRsnowCoeff%Classification_Name , &
                          Title                   = Title                           , &
                          History                 = History                         , &
                          Comment                 = Comment                           )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading global attribute from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check if release is valid
    IF ( .NOT. IRsnowCoeff_ValidRelease( IRsnowCoeff ) ) THEN
      msg = 'IRsnowCoeff Release check failed.'
      CALL Read_Cleanup(); RETURN
    END IF

    ! Read the IRsnowCoeff data
    ! ...Angle variable
    NF90_Status = NF90_INQ_VARID( FileId,ANGLE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//ANGLE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,IRsnowCoeff%Angle )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//ANGLE_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Frequency variable
    NF90_Status = NF90_INQ_VARID( FileId,FREQUENCY_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//FREQUENCY_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,IRsnowCoeff%Frequency )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//FREQUENCY_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Grain Size variable
    NF90_Status = NF90_INQ_VARID( FileId,GRAINSIZE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//GRAINSIZE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,IRsnowCoeff%Grain_Size )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//GRAINSIZE_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Temperature variable
    NF90_Status = NF90_INQ_VARID( FileId,TEMPERATURE_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//TEMPERATURE_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,IRsnowCoeff%Temperature )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//TEMPERATURE_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Emissivity variable
    NF90_Status = NF90_INQ_VARID( FileId,EMISSIVITY_VARNAME,VarId )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error inquiring '//TRIM(Filename)//' for '//EMISSIVITY_VARNAME//&
            ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF
    NF90_Status = NF90_GET_VAR( FileId,VarId,IRsnowCoeff%Emissivity )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error reading '//EMISSIVITY_VARNAME//' from '//TRIM(Filename)//&
            ' - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing output file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Output an info message
    IF ( Noisy ) THEN
      CALL IRsnowCoeff_Info( IRsnowCoeff, msg )
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
      CALL IRsnowCoeff_Destroy( IRsnowCoeff )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Read_CleanUp

  END FUNCTION IRsnowCoeff_netCDF_ReadFile


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

  ! Function to write the global attributes to a IRsnowCoeff data file.

  FUNCTION WriteGAtts( &
    Filename            , &  ! Input
    FileId              , &  ! Input
    Release             , &  ! Optional input
    Version             , &  ! Optional input
    Classification_Name , &  ! Optional input
    Title               , &  ! Optional input
    History             , &  ! Optional input
    Comment             ) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    INTEGER     ,           INTENT(IN) :: FileId
    INTEGER     , OPTIONAL, INTENT(IN) :: Release
    INTEGER     , OPTIONAL, INTENT(IN) :: Version
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Classification_Name
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Title
    CHARACTER(*), OPTIONAL, INTENT(IN) :: History
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Comment
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_WriteGAtts(netCDF)'
    CHARACTER(*), PARAMETER :: WRITE_MODULE_HISTORY_GATTNAME   = 'write_module_history'
    CHARACTER(*), PARAMETER :: CREATION_DATE_AND_TIME_GATTNAME = 'creation_date_and_time'
    ! Local variables
    CHARACTER(ML) :: msg
    CHARACTER(ML) :: GAttName
    CHARACTER(8)  :: cdate
    CHARACTER(10) :: ctime
    CHARACTER(5)  :: czone
    CHARACTER(ML) :: CLSname
    INTEGER :: Ver
    INTEGER :: NF90_Status
    TYPE(IRsnowCoeff_type) :: IRsnowCoeff

    ! Set up
    err_stat = SUCCESS
    msg = ' '

    ! Mandatory global attributes
    ! ...Software ID
    !GAttName = WRITE_MODULE_HISTORY_GATTNAME
    !NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),MODULE_VERSION_ID )
    !IF ( NF90_Status /= NF90_NOERR ) THEN
    !  CALL WriteGAtts_Cleanup(); RETURN
    !END IF
    ! ...Creation date
    CALL DATE_AND_TIME( cdate, ctime, czone )
    GAttName = CREATION_DATE_AND_TIME_GATTNAME
    NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName), &
                                cdate(1:4)//'/'//cdate(5:6)//'/'//cdate(7:8)//', '// &
                                ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//' '// &
                                czone//'UTC' )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      CALL WriteGAtts_Cleanup(); RETURN
    END IF

    ! Optional global attributes
    ! ...The Release
    IF ( PRESENT(Release) ) THEN
      GAttName = RELEASE_GATTNAME
      NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),IRsnowCoeff%Release )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL WriteGAtts_Cleanup(); RETURN
      END IF
    END IF

    ! ...The Version
    IF ( PRESENT(Version) ) THEN
      Ver = Version
    ELSE
      Ver = IRsnowCoeff%Version
    END IF
    GAttName = VERSION_GATTNAME
    NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),Ver )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      CALL WriteGAtts_Cleanup(); RETURN
    END IF
    ! ...The Classification_Name
    IF ( PRESENT(Classification_Name) ) THEN
      CLSname = Classification_Name
    ELSE
      CLSname = IRsnowCoeff%Classification_Name
    END IF
    GAttName = CLASSIFICATION_NAME_GATTNAME
    NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),CLSname )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      CALL WriteGAtts_Cleanup(); RETURN
    END IF
    ! ...The title
    IF ( PRESENT(title) ) THEN
      GAttName = TITLE_GATTNAME
      NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),title )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL WriteGAtts_Cleanup(); RETURN
      END IF
    END IF
    ! ...The history
    IF ( PRESENT(history) ) THEN
      GAttName = HISTORY_GATTNAME
      NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),history )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL WriteGAtts_Cleanup(); RETURN
      END IF
    END IF
    ! ...The comment
    IF ( PRESENT(comment) ) THEN
      GAttName = COMMENT_GATTNAME
      NF90_Status = NF90_PUT_ATT( FileId,NF90_GLOBAL,TRIM(GAttName),comment )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL WriteGAtts_Cleanup(); RETURN
      END IF
    END IF

  CONTAINS

     SUBROUTINE WriteGAtts_CleanUp()
       NF90_Status = NF90_CLOSE( FileId )
       IF ( NF90_Status /= NF90_NOERR ) &
         msg = '; Error closing input file during error cleanup - '//&
               TRIM(NF90_STRERROR( NF90_Status ) )
       err_stat = FAILURE
       CALL Display_Message( ROUTINE_NAME, &
                             'Error writing '//TRIM(GAttName)//' attribute to '//&
                             TRIM(Filename)//' - '// &
                             TRIM(NF90_STRERROR( NF90_Status ) )//TRIM(msg), &
                             err_stat )
     END SUBROUTINE WriteGAtts_CleanUp

   END FUNCTION WriteGAtts

   ! Function to read the global attributes from a IRsnowCoeff data file.

   FUNCTION ReadGAtts( &
     Filename            , &  ! Input
     FileId              , &  ! Input
     Release             , &  ! Optional output
     Version             , &  ! Optional output
     Classification_Name , &  ! Optional output
     Title               , &  ! Optional output
     History             , &  ! Optional output
     Comment             ) &  ! Optional output
   RESULT( err_stat )
     ! Arguments
     CHARACTER(*),           INTENT(IN)  :: Filename
     INTEGER     ,           INTENT(IN)  :: FileId
     INTEGER     , OPTIONAL, INTENT(OUT) :: Release
     INTEGER     , OPTIONAL, INTENT(OUT) :: Version
     CHARACTER(*), OPTIONAL, INTENT(OUT) :: Classification_Name
     CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title
     CHARACTER(*), OPTIONAL, INTENT(OUT) :: History
     CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment
     ! Function result
     INTEGER :: err_stat
     ! Local parameters
     CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_ReadGAtts(netCDF)'
     ! Local variables
     CHARACTER(ML)   :: msg
     CHARACTER(256)  :: GAttName
     CHARACTER(5000) :: GAttString
     INTEGER :: NF90_Status

     ! Set up
     err_stat = SUCCESS

     ! The global attributes
     ! ...The Release
     IF ( PRESENT(Release) ) THEN
       GAttName = RELEASE_GATTNAME
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),Release )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
     END IF
     ! ...The Version
     IF ( PRESENT(Version) ) THEN
       GAttName = VERSION_GATTNAME
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),Version )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
     END IF
     ! ...The Classification Name
     IF ( PRESENT(Classification_Name) ) THEN
       GAttName = CLASSIFICATION_NAME_GATTNAME; GAttString = ''
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),GAttString )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
       CALL StrClean( GAttString )
       Classification_Name = GAttString(1:MIN(LEN(Classification_Name), LEN_TRIM(GAttString)))
     END IF
     ! ...The title
     IF ( PRESENT(title) ) THEN
       GAttName = TITLE_GATTNAME; GAttString = ''
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),GAttString )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
       CALL StrClean( GAttString )
       title = GAttString(1:MIN(LEN(title), LEN_TRIM(GAttString)))
     END IF
     ! ...The history
     IF ( PRESENT(history) ) THEN
       GAttName = HISTORY_GATTNAME; GAttString = ''
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),GAttString )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
       CALL StrClean( GAttString )
       history = GAttString(1:MIN(LEN(history), LEN_TRIM(GAttString)))
     END IF
     ! ...The comment
     IF ( PRESENT(comment) ) THEN
       GAttName = COMMENT_GATTNAME; GAttString = ''
       NF90_Status = NF90_GET_ATT( FileID,NF90_GLOBAL,TRIM(GAttName),GAttString )
       IF ( NF90_Status /= NF90_NOERR ) THEN
         CALL ReadGAtts_Cleanup(); RETURN
       END IF
       CALL StrClean( GAttString )
       comment = GAttString(1:MIN(LEN(comment), LEN_TRIM(GAttString)))
     END IF

   CONTAINS

     SUBROUTINE ReadGAtts_CleanUp()
       err_stat = FAILURE
       msg = 'Error reading '//TRIM(GAttName)//' attribute from '//TRIM(Filename)//' - '// &
             TRIM(NF90_STRERROR( NF90_Status ) )
       CALL Display_Message( ROUTINE_NAME, msg, err_stat )
     END SUBROUTINE ReadGAtts_CleanUp

   END FUNCTION ReadGAtts

   FUNCTION CreateFile( &
     Filename            , &  ! Input
     n_Angles            , &  ! Input
     n_Frequencies       , &  ! Input
     n_Grain_Sizes       , &  ! Input
     n_Temperature       , &  ! Input
     FileId              , &  ! Output
     Release             , &  ! Optional input
     Version             , &  ! Optional input
     Classification_Name , &  ! Optional input
     Title               , &  ! Optional input
     History             , &  ! Optional input
     Comment             ) &  ! Optional input
   RESULT( err_stat )
     ! Arguments
     CHARACTER(*),           INTENT(IN)  :: Filename
     INTEGER     ,           INTENT(IN)  :: n_Angles
     INTEGER     ,           INTENT(IN)  :: n_Frequencies
     INTEGER     ,           INTENT(IN)  :: n_Grain_Sizes
     INTEGER     ,           INTENT(IN)  :: n_Temperature
     INTEGER     ,           INTENT(OUT) :: FileId
     INTEGER     , OPTIONAL, INTENT(IN)  :: Release
     INTEGER     , OPTIONAL, INTENT(IN)  :: Version
     CHARACTER(*), OPTIONAL, INTENT(IN)  :: Classification_Name
     CHARACTER(*), OPTIONAL, INTENT(IN)  :: Title
     CHARACTER(*), OPTIONAL, INTENT(IN)  :: History
     CHARACTER(*), OPTIONAL, INTENT(IN)  :: Comment
     ! Function result
     INTEGER :: err_stat
     ! Local parameters
     CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'IRsnowCoeff_CreateFile(netCDF)'
     ! Local variables
     CHARACTER(ML) :: msg
     LOGICAL :: Close_File
     INTEGER :: NF90_Status
     INTEGER :: n_Angles_DimID
     INTEGER :: n_Frequencies_DimID
     INTEGER :: n_Grain_Sizes_DimID
     INTEGER :: n_Temperature_DimID
     INTEGER :: varID
     INTEGER :: Put_Status(4)
     TYPE(IRsnowCoeff_type) :: dummy

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
     ! ...Number of angles
     NF90_Status = NF90_DEF_DIM( FileID,ANGLE_DIMNAME,n_Angles,n_Angles_DimID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//ANGLE_DIMNAME//' dimension in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Number of frequencies
     NF90_Status = NF90_DEF_DIM( FileID,FREQUENCY_DIMNAME,n_Frequencies,n_Frequencies_DimID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//FREQUENCY_DIMNAME//' dimension in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Number of Grain Size
     NF90_Status = NF90_DEF_DIM( FileID,GRAINSIZE_DIMNAME,n_Grain_Sizes,n_Grain_Sizes_DimID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//GRAINSIZE_DIMNAME//' dimension in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Number of temperature
     NF90_Status = NF90_DEF_DIM( FileID,TEMPERATURE_DIMNAME,n_Temperature,n_Temperature_DimID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//TEMPERATURE_DIMNAME//' dimension in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     CALL IRsnowCoeff_Destroy(dummy)

     ! Write the global attributes
     err_stat = WriteGAtts( Filename, &
                            FileId  , &
                            Release             = Release             , &
                            Version             = Version             , &
                            Classification_Name = Classification_Name , &
                            Title               = Title               , &
                            History             = History             , &
                            Comment             = Comment               )
     IF ( err_stat /= SUCCESS ) THEN
       msg = 'Error writing global attribute to '//TRIM(Filename)
       CALL Create_Cleanup(); RETURN
     END IF

     ! Define the variables
     ! ...Angle variable
     NF90_Status = NF90_DEF_VAR( FileID, &
       ANGLE_VARNAME, &
       ANGLE_TYPE, &
       dimIDs=(/n_Angles_DimID/), &
       varID=VarID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//ANGLE_VARNAME//' variable in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     Put_Status(1) = NF90_PUT_ATT( FileID,VarID,LONGNAME_ATTNAME   ,ANGLE_LONGNAME    )
     Put_Status(2) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME,ANGLE_DESCRIPTION )
     Put_Status(3) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME      ,ANGLE_UNITS       )
     Put_Status(4) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME  ,ANGLE_FILLVALUE   )
     IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
       msg = 'Error writing '//ANGLE_VARNAME//' variable attributes to '//TRIM(Filename)
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Frequency variable
     NF90_Status = NF90_DEF_VAR( FileID, &
       FREQUENCY_VARNAME, &
       FREQUENCY_TYPE, &
       dimIDs=(/n_Frequencies_DimID/), &
       varID=VarID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//FREQUENCY_VARNAME//' variable in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     Put_Status(1) = NF90_PUT_ATT( FileID,VarID,LONGNAME_ATTNAME   ,FREQUENCY_LONGNAME    )
     Put_Status(2) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME,FREQUENCY_DESCRIPTION )
     Put_Status(3) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME      ,FREQUENCY_UNITS       )
     Put_Status(4) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME  ,FREQUENCY_FILLVALUE   )
     IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
       msg = 'Error writing '//FREQUENCY_VARNAME//' variable attributes to '//TRIM(Filename)
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Grain Size variable
     NF90_Status = NF90_DEF_VAR( FileID, &
       GRAINSIZE_VARNAME, &
       GRAINSIZE_TYPE, &
       dimIDs=(/n_Grain_Sizes_DimID/), &
       varID=VarID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//GRAINSIZE_VARNAME//' variable in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     Put_Status(1) = NF90_PUT_ATT( FileID,VarID,LONGNAME_ATTNAME   ,GRAINSIZE_LONGNAME    )
     Put_Status(2) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME,GRAINSIZE_DESCRIPTION )
     Put_Status(3) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME      ,GRAINSIZE_UNITS       )
     Put_Status(4) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME  ,GRAINSIZE_FILLVALUE   )
     IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
       msg = 'Error writing '//GRAINSIZE_VARNAME//' variable attributes to '//TRIM(Filename)
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Temperature variable
     NF90_Status = NF90_DEF_VAR( FileID, &
       TEMPERATURE_VARNAME, &
       TEMPERATURE_TYPE, &
       dimIDs=(/n_Temperature_DimID/), &
       varID=VarID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//TEMPERATURE_VARNAME//' variable in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     Put_Status(1) = NF90_PUT_ATT( FileID,VarID,LONGNAME_ATTNAME   ,TEMPERATURE_LONGNAME    )
     Put_Status(2) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME,TEMPERATURE_DESCRIPTION )
     Put_Status(3) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME      ,TEMPERATURE_UNITS       )
     Put_Status(4) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME  ,TEMPERATURE_FILLVALUE   )
     IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
       msg = 'Error writing '//TEMPERATURE_VARNAME//' variable attributes to '//TRIM(Filename)
       CALL Create_Cleanup(); RETURN
     END IF
     ! ...Emissivity variable
     NF90_Status = NF90_DEF_VAR( FileID, &
       EMISSIVITY_VARNAME, &
       EMISSIVITY_TYPE, &
       dimIDs=(/n_Angles_DimID, n_Frequencies_DimID, n_Grain_Sizes_DimID, n_Temperature_DimID/), &
       varID=VarID )
     IF ( NF90_Status /= NF90_NOERR ) THEN
       msg = 'Error defining '//EMISSIVITY_VARNAME//' variable in '//&
             TRIM(Filename)//' - '//TRIM(NF90_STRERROR( NF90_Status ))
       CALL Create_Cleanup(); RETURN
     END IF
     Put_Status(1) = NF90_PUT_ATT( FileID,VarID,LONGNAME_ATTNAME   ,EMISSIVITY_LONGNAME    )
     Put_Status(2) = NF90_PUT_ATT( FileID,VarID,DESCRIPTION_ATTNAME,EMISSIVITY_DESCRIPTION )
     Put_Status(3) = NF90_PUT_ATT( FileID,VarID,UNITS_ATTNAME      ,EMISSIVITY_UNITS       )
     Put_Status(4) = NF90_PUT_ATT( FileID,VarID,FILLVALUE_ATTNAME  ,EMISSIVITY_FILLVALUE   )
     IF ( ANY(Put_Status /= NF90_NOERR) ) THEN
       msg = 'Error writing '//EMISSIVITY_VARNAME//' variable attributes to '//TRIM(Filename)
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

 END MODULE IRsnowCoeff_netCDF_IO
