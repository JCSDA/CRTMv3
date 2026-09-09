!
! VISsnowCoeff_IO
!
! Container module for Binary and netCDF VISsnowCoeff I/O modules.
! All Binary related modules are placeholder for now.
!
! CREATION HISTORY:
!
!       Written by:    Cheng Dang, May 2026
!                      dangch@ucar.edu

MODULE VISsnowCoeff_IO

    ! -----------------
    ! Environment setup
    ! -----------------
    ! Module use
    USE Type_Kinds             , ONLY: fp
    USE Message_Handler        , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
    USE Compare_Float_Numbers  , ONLY: OPERATOR(.EqualTo.)
    USE File_Utility           , ONLY: File_Exists
    USE VISsnowCoeff_Define    , ONLY: VISsnowCoeff_type, &
                                       OPERATOR(==), &
                                       VISsnowCoeff_Associated
    USE VISsnowCoeff_netCDF_IO  , ONLY: VISsnowCoeff_netCDF_InquireFile , &
                                        VISsnowCoeff_netCDF_ReadFile
                                        
    ! Disable implicit typing
    IMPLICIT NONE

    ! ------------
    ! Visibilities
    ! ------------
    PRIVATE
    PUBLIC :: VISsnowCoeff_InquireFile_IO
    PUBLIC :: VISsnowCoeff_ReadFile_IO


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
!       VISsnowCoeff_InquireFile
!
! PURPOSE:
!       Function to inquire VISsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = VISsnowCoeff_InquireFile( &
!                        Filename, &
!                        netCDF           = netCDF          , &
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
! OPTIONAL INPUTS:
!       netCDF:            Set this logical argument to access netCDF format
!                          VISsnowCoeff datafiles.
!                          If == .FALSE., file format is BINARY [DEFAULT].
!                             == .TRUE.,  file format is NETCDF.
!                          If not specified, default is .FALSE.
!                          UNITS:      N/A
!                          TYPE:       LOGICAL
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN), OPTIONAL
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
!       n_Depths:          The number of depths in
!                          the LUT. Must be > 0.
!                          UNITS:      N/A
!                          TYPE:       INTEGER
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: INTENT(IN)

!       n_Densities:       The number of densities in
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

  FUNCTION VISsnowCoeff_InquireFile_IO( &
    Filename       , &  ! Input
    netCDF         , &  ! Optional input
    n_Angles       , &  ! Optional output
    n_Frequencies  , &  ! Optional output
    n_Grain_Sizes  , &  ! Optional output
    n_Depths         , &  ! Optional output
    n_Densities      , &  ! Optional output
    Release        , &  ! Optional output
    Version        , &  ! Optional output
    Title          , &  ! Optional output
    History        , &  ! Optional output
    Comment        ) &  ! Optional output
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN)  :: Filename
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Angles
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Frequencies
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Grain_Sizes
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Depths
    INTEGER     , OPTIONAL, INTENT(OUT) :: n_Densities
    LOGICAL     , OPTIONAL, INTENT(IN)  :: netCDF
    INTEGER     , OPTIONAL, INTENT(OUT) :: Release
    INTEGER     , OPTIONAL, INTENT(OUT) :: Version
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Title
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: History
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Comment
    ! Function result
    INTEGER :: err_stat
    ! Function variables
    LOGICAL :: Binary

    ! Set up
    err_stat = SUCCESS
    ! ...Check netCDF argument
    Binary = .FALSE.
    IF ( PRESENT(netCDF) ) Binary = .NOT. netCDF


    ! Call the appropriate function
    err_stat = VISsnowCoeff_netCDF_InquireFile( &
                Filename                         , &
                n_Angles        = n_Angles       , &
                n_Frequencies   = n_Frequencies  , &
                n_Grain_Sizes   = n_Grain_Sizes  , &
                n_Depths        = n_Depths       , &
                n_Densities     = n_Densities    , &
                Release         = Release        , &
                Version         = Version          )
    

  END FUNCTION VISsnowCoeff_InquireFile_IO

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       VISsnowCoeff_ReadFile
!
! PURPOSE:
!       Function to read VISsnowCoeff object files.
!
! CALLING SEQUENCE:
!       Error_Status = VISsnowCoeff_ReadFile( &
!                        VISsnowCoeff, &
!                        Filename, &
!                        netCDF  = netCDF , &
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
!       netCDF:         Set this logical argument to access netCDF format
!                       VISsnowCoeff datafiles.
!                       If == .FALSE., file format is BINARY [DEFAULT].
!                          == .TRUE.,  file format is NETCDF.
!                       If not specified, default is .FALSE.
!                       UNITS:      N/A
!                       TYPE:       LOGICAL
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN), OPTIONAL
!
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
  FUNCTION VISsnowCoeff_ReadFile_IO( &
    VISsnowCoeff , &  ! Output
    Filename    , &  ! Input
    netCDF      , &  ! Optional input
    No_Close    , &  ! Optional input
    Quiet       , &  ! Optional input
    Title       , &  ! Optional output
    History     , &  ! Optional output
    Comment     , &  ! Optional output
    Debug       ) &  ! Optional input (Debug output control)
  RESULT( err_stat )
    ! Arguments
    TYPE(VISsnowCoeff_type), INTENT(OUT) :: VISsnowCoeff
    CHARACTER(*),           INTENT(IN)   :: Filename
    LOGICAL,      OPTIONAL, INTENT(IN)   :: netCDF
    LOGICAL,      OPTIONAL, INTENT(IN)   :: No_Close
    LOGICAL,      OPTIONAL, INTENT(IN)   :: Quiet
    CHARACTER(*), OPTIONAL, INTENT(OUT)  :: Title
    CHARACTER(*), OPTIONAL, INTENT(OUT)  :: History
    CHARACTER(*), OPTIONAL, INTENT(OUT)  :: Comment
    LOGICAL,      OPTIONAL, INTENT(IN)   :: Debug
    ! Function result
    INTEGER :: err_stat
    ! Function variables
    LOGICAL :: Binary

    ! Set up
    err_stat = SUCCESS
    ! ...Check netCDF argument
    Binary = .FALSE.
    IF ( PRESENT(netCDF) ) Binary = .NOT. netCDF

    ! Call the appropriate function
    err_stat = VISsnowCoeff_netCDF_ReadFile( &
                  VISsnowCoeff , &
                  Filename    , &
                  Quiet       , &
                  Title       , &
                  History     , &
                  Comment     , &
                  Debug       )
    
  END FUNCTION VISsnowCoeff_ReadFile_IO


END MODULE VISsnowCoeff_IO
