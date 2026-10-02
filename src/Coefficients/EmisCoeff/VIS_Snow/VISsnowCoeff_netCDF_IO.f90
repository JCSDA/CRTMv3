!
! VISsnowCoeff_netCDF_IO
!
! Module containing routines to read and write VISsnowCoeff netCDF
! format files.
!
!
! CREATION HISTORY:
!
!       Written by:    Cheng Dang, May 2026
!                      dangch@ucar.edu

MODULE VISsnowCoeff_netCDF_IO

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds          , ONLY: fp, Double, Long
  USE Message_Handler     , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility        , ONLY: File_Exists
  USE String_Utility      , ONLY: StrClean
  USE VISsnowCoeff_Define , ONLY: VISsnowCoeff_type, &
                                  VISsnowCoeff_Associated, &
                                  VISsnowCoeff_Create, &
                                  VISsnowCoeff_Inspect, &
                                  VISsnowCoeff_Destroy, &
                                  VISsnowCoeff_ValidRelease, &
                                  VISsnowCoeff_Info
  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Procedures
  PUBLIC :: VISsnowCoeff_netCDF_InquireFile
  PUBLIC :: VISsnowCoeff_netCDF_ReadFile

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
  CHARACTER(*), PARAMETER :: TITLE_GATTNAME       = 'Title'
  CHARACTER(*), PARAMETER :: HISTORY_GATTNAME     = 'History'
  CHARACTER(*), PARAMETER :: COMMENT_GATTNAME     = 'Comment'
  CHARACTER(*), PARAMETER :: CLASSIFICATION_NAME_GATTNAME = 'Classification_Name'

  ! Dimension names
  CHARACTER(*), PARAMETER :: TNSL_DIMNAME         = 'String_Length'
  CHARACTER(*), PARAMETER :: FREQUENCY_DIMNAME    = 'n_Frequencies'
  CHARACTER(*), PARAMETER :: ANGLE_DIMNAME        = 'n_Angles'
  CHARACTER(*), PARAMETER :: GRAINSIZE_DIMNAME    = 'n_Grain_Sizes'
  CHARACTER(*), PARAMETER :: DEPTH_DIMNAME        = 'n_Depths'
  CHARACTER(*), PARAMETER :: DENSITY_DIMNAME      = 'n_Densities'

  ! Variable names
  CHARACTER(*), PARAMETER :: ANGLE_VARNAME        = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_VARNAME    = 'Frequency'
  CHARACTER(*), PARAMETER :: GRAINSIZE_VARNAME    = 'Grain_Size'
  CHARACTER(*), PARAMETER :: DEPTH_VARNAME        = 'Depth'
  CHARACTER(*), PARAMETER :: DENSITY_VARNAME      = 'Density'
  CHARACTER(*), PARAMETER :: REFLECTANCE_VARNAME  = 'Reflectance'

  ! Variable long name attribute.
  CHARACTER(*), PARAMETER :: LONGNAME_ATTNAME = 'long_name'
  CHARACTER(*), PARAMETER :: ANGLE_LONGNAME        = 'Angle'
  CHARACTER(*), PARAMETER :: FREQUENCY_LONGNAME    = 'Frequency'
  CHARACTER(*), PARAMETER :: GRAINSIZE_LONGNAME    = 'Grain Size'
  CHARACTER(*), PARAMETER :: DEPTH_LONGNAME        = 'Depth'
  CHARACTER(*), PARAMETER :: DENSITY_LONGNAME      = 'Density'
  CHARACTER(*), PARAMETER :: REFLECTANCE_LONGNAME  = 'Reflectance'

  ! Variable description attribute.
  CHARACTER(*), PARAMETER :: DESCRIPTION_ATTNAME = 'description'
  CHARACTER(*), PARAMETER :: ANGLE_DESCRIPTION        = 'Angle dimension values for reflectance data'
  CHARACTER(*), PARAMETER :: FREQUENCY_DESCRIPTION    = 'Frequency dimension values for reflectance data'
  CHARACTER(*), PARAMETER :: GRAINSIZE_DESCRIPTION    = 'Grain Size dimension values for reflectance data'
  CHARACTER(*), PARAMETER :: DEPTH_DESCRIPTION        = 'Depth dimension values for reflectance data'
  CHARACTER(*), PARAMETER :: DENSITY_DESCRIPTION      = 'Density   dimension values for reflectance data'
  CHARACTER(*), PARAMETER :: REFLECTANCE_DESCRIPTION  = 'Spectral snow surface reflectance data'

  ! Variable units attribute.
  CHARACTER(*), PARAMETER :: UNITS_ATTNAME = 'units'
  CHARACTER(*), PARAMETER :: ANGLE_UNITS        = 'degrees from vertical'
  CHARACTER(*), PARAMETER :: FREQUENCY_UNITS    = 'inverse centimeters (cm^-1)'
  CHARACTER(*), PARAMETER :: GRAINSIZE_UNITS    = 'effective radius in microns (um)'
  CHARACTER(*), PARAMETER :: DEPTH_UNITS        = 'meters'
  CHARACTER(*), PARAMETER :: DENSITY_UNITS      = 'kg/m^3'
  CHARACTER(*), PARAMETER :: REFLECTANCE_UNITS  = 'N/A'

  ! Variable _FillValue attribute.
  CHARACTER(*), PARAMETER :: FILLVALUE_ATTNAME = '_FillValue'
  REAL(Double), PARAMETER :: ANGLE_FILLVALUE        = FILL_FLOAT
  REAL(Double), PARAMETER :: FREQUENCY_FILLVALUE    = FILL_FLOAT
  REAL(Double), PARAMETER :: GRAINSIZE_FILLVALUE    = FILL_FLOAT
  REAL(Double), PARAMETER :: DEPTH_FILLVALUE        = FILL_FLOAT
  REAL(Double), PARAMETER :: DENSITY_FILLVALUE      = FILL_FLOAT
  REAL(Double), PARAMETER :: REFLECTANCE_FILLVALUE  = FILL_FLOAT

  ! Variable types
  INTEGER, PARAMETER :: ANGLE_TYPE        = NF90_DOUBLE
  INTEGER, PARAMETER :: FREQUENCY_TYPE    = NF90_DOUBLE
  INTEGER, PARAMETER :: GRAINSIZE_TYPE    = NF90_DOUBLE
  INTEGER, PARAMETER :: DEPTH_TYPE        = NF90_DOUBLE
  INTEGER, PARAMETER :: DENSITY_TYPE      = NF90_DOUBLE
  INTEGER, PARAMETER :: REFLECTANCE_TYPE  = NF90_DOUBLE


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
!       VISsnowCoeff_netCDF_InquireFile
!
! PURPOSE:
!       Function to inquire VISsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = VISsnowCoeff_netCDF_InquireFile( &
!                        Filename, &
!                        n_Angles         = n_Angles        , &
!                        n_Frequencies    = n_Frequencies   , &
!                        n_Grain_Sizes    = n_Grain_Sizes   , &
!                        n_Depths         = n_Depths        , &
!                        n_Densities      = n_Densities     , &
!                        Release          = Release         , &
!                        Version          = Version         , &
!                        Title            = Title           , &
!                        History          = History         , &
!                        Comment          = Comment           )
!
! INPUTS:
!       Filename:          Character string specifying the name of a
!                          VISsnowCoeff data file to read.
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
!       Release:           The release number of the VISsnowCoeff file.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Version:           The version number of the VISsnowCoeff file.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Title:             Character string written into the TITLE global
!                          attribute field of the VISsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       History:           Character string written into the HISTORY global
!                          attribute field of the VISsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Comment:           Character string written into the COMMENT global
!                          attribute field of the VISsnowCoeff file.
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

  FUNCTION VISsnowCoeff_netCDF_InquireFile( &
    Filename            , &  ! Input
    n_Angles            , &  ! Optional output
    n_Frequencies       , &  ! Optional output
    n_Grain_Sizes       , &  ! Optional output
    n_Depths            , &  ! Optional output
    n_Densities         , &  ! Optional output
    Release             , &  ! Optional output
    Version             , &  ! Optional output
    Classification_Name , &  ! Optional output
    Title               , &  ! Optional output
    History             , &  ! Optional output
    Comment             ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Angles
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Frequencies
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Grain_Sizes
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Depths
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Densities
    INTEGER     , OPTIONAL, INTENT(OUT) :: Release
    INTEGER     , OPTIONAL, INTENT(OUT) :: Version
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Classification_Name
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: History
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'VISsnowCoeff_netCDF_InquireFile'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    INTEGER :: NF90_Status
    INTEGER :: FileId
    INTEGER :: DimId
    TYPE(VISsnowCoeff_type) :: VISsnowCoeff

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
    CALL Read_Dim( ANGLE_DIMNAME,     VISsnowCoeff%n_Angles,      err_stat ); IF (err_stat/=SUCCESS) RETURN
    CALL Read_Dim( FREQUENCY_DIMNAME, VISsnowCoeff%n_Frequencies, err_stat ); IF (err_stat/=SUCCESS) RETURN
    CALL Read_Dim( GRAINSIZE_DIMNAME, VISsnowCoeff%n_Grain_Sizes, err_stat ); IF (err_stat/=SUCCESS) RETURN
    CALL Read_Dim( DEPTH_DIMNAME,     VISsnowCoeff%n_Depths,      err_stat ); IF (err_stat/=SUCCESS) RETURN
    CALL Read_Dim( DENSITY_DIMNAME,   VISsnowCoeff%n_Densities,   err_stat ); IF (err_stat/=SUCCESS) RETURN

    ! Get the global attributes
    err_stat = ReadGAtts( Filename, &
                          FileId  , &
                          Release = VISsnowCoeff%Release, &
                          Version = VISsnowCoeff%Version, &
                          Classification_Name = VISsnowCoeff%Classification_Name )
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
    IF ( PRESENT(n_Angles     ) ) n_Angles       = VISsnowCoeff%n_Angles
    IF ( PRESENT(n_Frequencies) ) n_Frequencies  = VISsnowCoeff%n_Frequencies
    IF ( PRESENT(n_Grain_Sizes) ) n_Grain_Sizes  = VISsnowCoeff%n_Grain_Sizes
    IF ( PRESENT(n_Depths     ) ) n_Depths       = VISsnowCoeff%n_Depths
    IF ( PRESENT(n_Densities  ) ) n_Densities    = VISsnowCoeff%n_Densities
    IF ( PRESENT(Release      ) ) Release        = VISsnowCoeff%Release
    IF ( PRESENT(Version      ) ) Version        = VISsnowCoeff%Version
    IF ( PRESENT(Classification_Name) ) Classification_Name = VISsnowCoeff%Classification_Name

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

    SUBROUTINE Read_Dim( DimName, DimValue, Error_Status )
      CHARACTER(*), INTENT(IN)  :: DimName
      INTEGER,      INTENT(OUT) :: DimValue
      INTEGER,      INTENT(OUT) :: Error_Status
      Error_Status = SUCCESS
      NF90_Status = NF90_INQ_DIMID( FileId, DimName, DimId )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error inquiring dimension ID for '//DimName//' - '// &
              TRIM(NF90_STRERROR( NF90_Status ))
        CALL Inquire_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
      NF90_Status = NF90_INQUIRE_DIMENSION( FileId, DimId, Len=DimValue )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error reading dimension value for '//DimName//' - '// &
              TRIM(NF90_STRERROR( NF90_Status ))
        CALL Inquire_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
    END SUBROUTINE Read_Dim

  END FUNCTION VISsnowCoeff_netCDF_InquireFile



!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_netCDF_ReadFile
!
! PURPOSE:
!       Function to read VISsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = VISsnowCoeff_netCDF_ReadFile( &
!                        VISsnowCoeff, &
!                        Filename, &
!                        Quiet   = Quiet  , &
!                        Title   = Title  , &
!                        History = History, &
!                        Comment = Comment  )
!
! INPUTS:
!       Filename:          Character string specifying the name of a
!                          VISsnowCoeff data file to read.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       VISsnowCoeff:   Object containing the VISsnow coefficient data.
!                       UNITS:      N/A
!                       TYPE:       TYPE(VISsnowCoeff_type)
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
!                          attribute field of the VISsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       History:           Character string written into the HISTORY global
!                          attribute field of the VISsnowCoeff file.
!                          This argument is ignored if the netCDF argument
!                          is not supplied or set.
!                          UNITS:      N/A
!                          TYPE:       CHARACTER(*)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(OUT), OPTIONAL
!
!       Comment:           Character string written into the COMMENT global
!                          attribute field of the VISsnowCoeff file.
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
  FUNCTION VISsnowCoeff_netCDF_ReadFile( &
    VISsnowCoeff , &  ! Output
    Filename    , &  ! Input
    Quiet       , &  ! Optional input
    Title       , &  ! Optional output
    History     , &  ! Optional output
    Comment     , &  ! Optional output
    Debug       ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(VISsnowCoeff_type) , INTENT(OUT) :: VISsnowCoeff
    CHARACTER(*),            INTENT(IN)   :: Filename
    LOGICAL     ,  OPTIONAL, INTENT(IN)   :: Quiet
    CHARACTER(*),  OPTIONAL, INTENT(OUT)  :: Title
    CHARACTER(*),  OPTIONAL, INTENT(OUT)  :: History
    CHARACTER(*),  OPTIONAL, INTENT(OUT)  :: Comment
    LOGICAL     ,  OPTIONAL, INTENT(IN)   :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'VISsnowCoeff_netCDF_ReadFile'
    ! Function variables
    CHARACTER(ML) :: msg
    LOGICAL :: Close_File
    LOGICAL :: Noisy
    INTEGER :: NF90_Status
    INTEGER :: FileId
    INTEGER :: n_Angles
    INTEGER :: n_Frequencies
    INTEGER :: n_Grain_Sizes
    INTEGER :: n_Depths
    INTEGER :: n_Densities
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
    err_stat = VISsnowCoeff_netCDF_InquireFile( &
                 Filename                       , &
                 n_Angles      = n_Angles       , &
                 n_Frequencies = n_Frequencies  , &
                 n_Grain_Sizes = n_Grain_Sizes  , &
                 n_Depths      = n_Depths       , &
                 n_Densities   = n_Densities      )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error obtaining VISsnowCoeff dimensions from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF

    ! Allocate the output structure
    CALL VISsnowCoeff_Create( &
           VISsnowCoeff   , &
           n_Angles       , &
           n_Frequencies  , &
           n_Grain_Sizes  , &
           n_Depths       , &
           n_Densities      )
    IF ( .NOT. VISsnowCoeff_Associated( VISsnowCoeff ) ) THEN
      msg = 'VISsnowCoeff object allocation failed.'
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
                          Release             = VISsnowCoeff%Release            , &
                          Version             = VISsnowCoeff%Version            , &
                          Classification_Name = VISsnowCoeff%Classification_Name  )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading global attribute from '//TRIM(Filename)
      CALL Read_Cleanup(); RETURN
    END IF
    ! ...Check if release is valid
    IF ( .NOT. VISsnowCoeff_ValidRelease( VISsnowCoeff ) ) THEN
      msg = 'VISsnowCoeff Release check failed.'
      CALL Read_Cleanup(); RETURN
    END IF

    ! Read the VISsnowCoeff data
    CALL Read_Var_1D_Real( ANGLE_VARNAME,       VISsnowCoeff%Angle,       err_stat ); IF ( err_stat /= SUCCESS ) RETURN
    CALL Read_Var_1D_Real( FREQUENCY_VARNAME,   VISsnowCoeff%Frequency,   err_stat ); IF ( err_stat /= SUCCESS ) RETURN
    CALL Read_Var_1D_Real( GRAINSIZE_VARNAME,   VISsnowCoeff%Grain_Size,  err_stat ); IF ( err_stat /= SUCCESS ) RETURN
    CALL Read_Var_1D_Real( DEPTH_VARNAME,       VISsnowCoeff%Depth,       err_stat ); IF ( err_stat /= SUCCESS ) RETURN
    CALL Read_Var_1D_Real( DENSITY_VARNAME,     VISsnowCoeff%Density,     err_stat ); IF ( err_stat /= SUCCESS ) RETURN    
    CALL Read_Var_5D_Real( REFLECTANCE_VARNAME, VISsnowCoeff%Reflectance, err_stat ); IF ( err_stat /= SUCCESS ) RETURN

    ! Close the file
    NF90_Status = NF90_CLOSE( FileId ); Close_File = .FALSE.
    IF ( NF90_Status /= NF90_NOERR ) THEN
      msg = 'Error closing output file - '//TRIM(NF90_STRERROR( NF90_Status ))
      CALL Read_Cleanup(); RETURN
    END IF

    ! Output an info message
    IF ( Noisy ) THEN
      CALL VISsnowCoeff_Info( VISsnowCoeff, msg )
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
      CALL VISsnowCoeff_Destroy( VISsnowCoeff )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
    END SUBROUTINE Read_CleanUp

    SUBROUTINE Read_Var_1D_Real( VarName, VarData, Error_Status )
      CHARACTER(*), INTENT(IN)  :: VarName
      REAL(fp),     INTENT(OUT) :: VarData(:)
      INTEGER,      INTENT(OUT) :: Error_Status
      Error_Status = SUCCESS
      NF90_Status = NF90_INQ_VARID( FileId, VarName, VarId )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error inquiring '//TRIM(Filename)//' for '//VarName// &
              ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
        CALL Read_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
      NF90_Status = NF90_GET_VAR( FileId, VarId, VarData )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error reading '//VarName//' from '//TRIM(Filename)// &
              ' - '//TRIM(NF90_STRERROR( NF90_Status ))
        CALL Read_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
    END SUBROUTINE Read_Var_1D_Real

    SUBROUTINE Read_Var_5D_Real( VarName, VarData, Error_Status )
      CHARACTER(*), INTENT(IN)  :: VarName
      REAL(fp),     INTENT(OUT) :: VarData(:,:,:,:,:)
      INTEGER,      INTENT(OUT) :: Error_Status
      Error_Status = SUCCESS
      NF90_Status = NF90_INQ_VARID( FileId, VarName, VarId )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error inquiring '//TRIM(Filename)//' for '//VarName// &
              ' variable ID - '//TRIM(NF90_STRERROR( NF90_Status ))
        CALL Read_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
      NF90_Status = NF90_GET_VAR( FileId, VarId, VarData )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        msg = 'Error reading '//VarName//' from '//TRIM(Filename)// &
              ' - '//TRIM(NF90_STRERROR( NF90_Status ))
        CALL Read_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
    END SUBROUTINE Read_Var_5D_Real

  END FUNCTION VISsnowCoeff_netCDF_ReadFile


!##################################################################################
!##################################################################################
!##                                                                              ##
!##                          ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                              ##
!##################################################################################
!##################################################################################

  ! Function to read the global attributes from a VISsnowCoeff data file.

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
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'VISsnowCoeff_ReadGAtts(netCDF)'
    ! Local variables
    CHARACTER(ML)   :: msg
    CHARACTER(256)  :: GAttName
    CHARACTER(5000) :: GAttString
    INTEGER :: NF90_Status

    ! Set up
    err_stat = SUCCESS

    ! The global attributes
    IF ( PRESENT(Release) ) THEN
      CALL Read_GAtt_Int( RELEASE_GATTNAME, Release, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF
    IF ( PRESENT(Version) ) THEN
      CALL Read_GAtt_Int( VERSION_GATTNAME, Version, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF
    IF ( PRESENT(Classification_Name) ) THEN
      CALL Read_GAtt_Str( CLASSIFICATION_NAME_GATTNAME, Classification_Name, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF
    IF ( PRESENT(title) ) THEN
      CALL Read_GAtt_Str( TITLE_GATTNAME, title, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF
    IF  ( PRESENT(history) ) THEN
      CALL Read_GAtt_Str( HISTORY_GATTNAME, history, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF
    IF  ( PRESENT(comment) ) THEN
      CALL Read_GAtt_Str( COMMENT_GATTNAME, comment, err_stat )
      IF ( err_stat/=SUCCESS ) RETURN
    END IF

  CONTAINS

    SUBROUTINE ReadGAtts_CleanUp()
      err_stat = FAILURE
      msg = 'Error reading '//TRIM(GAttName)//' attribute from '//TRIM(Filename)//' - '// &
            TRIM(NF90_STRERROR( NF90_Status ) )
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE ReadGAtts_CleanUp

    SUBROUTINE Read_GAtt_Int( AttName, AttValue, Error_Status )
      CHARACTER(*), INTENT(IN)  :: AttName
      INTEGER,      INTENT(OUT) :: AttValue
      INTEGER,      INTENT(OUT) :: Error_Status
      Error_Status = SUCCESS
      GAttName    = AttName
      NF90_Status = NF90_GET_ATT( FileID, NF90_GLOBAL, TRIM(GAttName), AttValue )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL ReadGAtts_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
    END SUBROUTINE Read_GAtt_Int

    SUBROUTINE Read_GAtt_Str( AttName, AttValue, Error_Status )
      CHARACTER(*), INTENT(IN)  :: AttName
      CHARACTER(*), INTENT(OUT) :: AttValue
      INTEGER,      INTENT(OUT) :: Error_Status
      Error_Status = SUCCESS
      GAttName    = AttName
      GAttString  = ''
      NF90_Status = NF90_GET_ATT( FileID, NF90_GLOBAL, TRIM(GAttName), GAttString )
      IF ( NF90_Status /= NF90_NOERR ) THEN
        CALL ReadGAtts_Cleanup(); Error_Status = FAILURE; RETURN
      END IF
      CALL StrClean( GAttString )
      AttValue = GAttString(1:MIN(LEN(AttValue), LEN_TRIM(GAttString)))
    END SUBROUTINE Read_GAtt_Str

  END FUNCTION ReadGAtts

 END MODULE VISsnowCoeff_netCDF_IO
