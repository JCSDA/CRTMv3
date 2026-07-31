!
! CRTM_MWwaterCoeff
!
! Module containing the shared CRTM microwave water surface emissivity
! data and their load/destruction routines.
!
! PUBLIC DATA:
!   MWwaterC:  Data structure containing the microwave water surface
!              emissivity data.
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure MWwaterC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CRTM initialisation.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 12-Nov-2011
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_MWwaterCoeff

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Message_Handler    , ONLY: SUCCESS, FAILURE, Display_Message
  USE File_Utility          , ONLY: Join_Path
  USE MWwaterCoeff_Define, ONLY: MWwaterCoeff_type      , &
                                 MWwaterCoeff_Associated, &
                                 MWwaterCoeff_Destroy   , &
                                 MWwaterCoeff_Readfile
  USE MWwaterCoeff_FASTEM6, ONLY: Load_FASTEM6 => MWwaterCoeff_LoadCoeffs
  USE MWwaterCoeff_FASTEM4, ONLY: Load_FASTEM4 => MWwaterCoeff_LoadCoeffs

  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared data
  PUBLIC :: MWwaterC
  ! Procedures
  PUBLIC :: CRTM_MWwaterCoeff_Load_FASTEM
  PUBLIC :: CRTM_MWwaterCoeff_Load
  PUBLIC :: CRTM_MWwaterCoeff_Destroy
  PUBLIC :: CRTM_MWwaterCoeff_IsLoaded


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 512


  ! --------------------------------------------------
  ! The shared microwave water surface emissivity data
  ! --------------------------------------------------
  TYPE(MWwaterCoeff_type), TARGET, SAVE :: MWwaterC


CONTAINS

!---------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWwaterCoeff_Load_FASTEM
!
! PURPOSE:
!       Function to load the microwave water surface emissivity data into
!       the public data structure MWwaterC
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWwaterCoeff_Load_FASTEM( &
!                        FASTEM_Scheme,                         &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!   FASTEM_Scheme:          Name of the FASTEM Scheme
!                           OPTIONS: FASTEM4 or FASTEM6
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATIOn message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION messages are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output.
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_MWwaterCoeff_Load_FASTEM( &
      FASTEM_Scheme    , &  ! Input
      Quiet            , &  ! Optional input
      Process_ID       , &  ! Optional input
      Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: FASTEM_Scheme
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWwaterCoeff_Load_FASTEM'
    CHARACTER(ML) :: msg, pid_msg
    LOGICAL :: noisy

    ! Setup
    err_stat = SUCCESS
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Check the MPI Process Ids
    IF ( noisy .AND. PRESENT(Process_ID) .AND. PRESENT(Output_Process_ID) ) THEN
      IF ( Process_Id /= Output_Process_Id ) noisy = .FALSE.
    END IF
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF

    ! Discard anything already loaded before loading a different scheme.
    ! FitCoeff_SetValue only allocates when the target is unassociated, and
    ! otherwise rejects a shape mismatch by DESTROYING the structure and
    ! returning no status. The FASTEM4 and FASTEM6 azimuth coefficients have
    ! different shapes, so switching scheme without this left the shared
    ! MWwaterC deallocated while this function still reported SUCCESS. That
    ! matters for polarimetric work specifically: FASTEM6 is the default and
    ! has no third or fourth Stokes azimuth model, so anyone wanting a
    ! polarimetric surface has to switch to FASTEM4 or FASTEM5.
    IF ( MWwaterCoeff_Associated( MWwaterC ) ) CALL MWwaterCoeff_Destroy( MWwaterC )

    ! Load MWwaterCoeff data
    SELECT CASE ( FASTEM_Scheme )
      CASE ( 'FASTEM6' )
        CALL Load_FASTEM6( MWwaterC )
      CASE ( 'FASTEM4' )
        CALL Load_FASTEM4( MWwaterC )
      CASE DEFAULT
        err_stat = FAILURE
        msg = 'This option not yet implemented: '//TRIM(FASTEM_Scheme)//TRIM(pid_msg)
        CALL Display_Message( ROUTINE_NAME, msg, err_stat )
        RETURN
    END SELECT

    ! The loaders return no status of their own, and FitCoeff_SetValue signals
    ! failure by leaving the structure unassociated. Check rather than assume.
    IF ( .NOT. MWwaterCoeff_Associated( MWwaterC ) ) THEN
      err_stat = FAILURE
      msg = 'MWwaterCoeff structure is unassociated after loading '// &
            TRIM(FASTEM_Scheme)//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF

  END FUNCTION CRTM_MWwaterCoeff_Load_FASTEM

!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWwaterCoeff_Load
!
! PURPOSE:
!       Function to load the microwave water surface emissivity data into
!       the public data structure MWwaterC
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWwaterCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the MWwaterCoeff file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          Character string specifying a file path for the
!                           input data file. If not specified, the current
!                           directory is the default.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Quiet:              Set this logical argument to suppress INFORMATION
!                           messages being printed to stdout
!                           If == .FALSE., INFORMATION messages are OUTPUT [DEFAULT].
!                              == .TRUE.,  INFORMATION messages are SUPPRESSED.
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Process_ID:         Set this argument to the MPI process ID that this
!                           function call is running under. This value is used
!                           solely for controlling INFORMATIOn message output.
!                           If MPI is not being used, ignore this argument.
!                           This argument is ignored if the Quiet argument is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Output_Process_ID:  Set this argument to the MPI process ID in which
!                           all INFORMATION messages are to be output. If
!                           the passed Process_ID value agrees with this value
!                           the INFORMATION messages are output.
!                           This argument is ignored if the Quiet argument
!                           is set.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:       The return value is an integer defining the error
!                           status. The error codes are defined in the
!                           Message_Handler module.
!                           If == SUCCESS the data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_MWwaterCoeff_Load( &
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWwaterCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(:), ALLOCATABLE :: MWwaterCoeff_File
    LOGICAL :: noisy

    ! Setup
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    MWwaterCoeff_File = TRIM(ADJUSTL(Filename))
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) MWwaterCoeff_File = Join_Path(File_Path, MWwaterCoeff_File)
    ! ...Check Quiet argument
    noisy = .TRUE.
    IF ( PRESENT(Quiet) ) noisy = .NOT. Quiet
    ! ...Check the MPI Process Ids
    IF ( noisy .AND. PRESENT(Process_ID) .AND. PRESENT(Output_Process_ID) ) THEN
      IF ( Process_Id /= Output_Process_Id ) noisy = .FALSE.
    END IF
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF


    ! Read the MWwaterCoeff data file
    err_stat = MWwaterCoeff_ReadFile( &
                 MWwaterC, &
                 MWwaterCoeff_File, &
                 Quiet = .NOT. noisy )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading MWwaterCoeff file '//TRIM(MWwaterCoeff_File)//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
      RETURN
    END IF

  END FUNCTION CRTM_MWwaterCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWwaterCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure MWwaterC containing
!       the CRTM microwave water surface emissivity data.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWwaterCoeff_Destroy( Process_ID = Process_ID )
!
! OPTIONAL INPUTS:
!       Process_ID:       Set this argument to the MPI process ID that this
!                         function call is running under. This value is used
!                         solely for controlling message output. If MPI is not
!                         being used, ignore this argument.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: INTENT(IN), OPTIONAL
!
! FUNCTION RESULT:
!       Error_Status:     The return value is an integer defining the error
!                         status. The error codes are defined in the
!                         Message_Handler module.
!                         If == SUCCESS the deallocation of the public data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data
!       structure MWwaterC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_MWwaterCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWwaterCoeff_Destroy'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg

    ! Setup
    err_stat = SUCCESS
    ! ...Create a process ID message tag for error messages
    IF ( PRESENT(Process_Id) ) THEN
      WRITE( pid_msg,'("; Process ID: ",i0)' ) Process_ID
    ELSE
      pid_msg = ''
    END IF

    ! Destroy the structure
    CALL MWwaterCoeff_Destroy( MWwaterC )
    IF ( MWwaterCoeff_Associated( MWwaterC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating MWwaterCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF

  END FUNCTION CRTM_MWwaterCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWwaterCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if microwave water surface emissivity data has
!       been loaded into the public data structure MWwaterC.
!
! CALLING SEQUENCE:
!       status = CRTM_MWwaterCoeff_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_MWwaterCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = MWwaterCoeff_Associated( MWwaterC )
  END FUNCTION CRTM_MWwaterCoeff_IsLoaded

END MODULE CRTM_MWwaterCoeff
