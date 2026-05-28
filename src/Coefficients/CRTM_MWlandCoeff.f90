!
! CRTM_MWlandCoeff
!
! Module containing the shared CRTM microwave land surface emissivity atlas
! (TELSEM2) data and its load/destruction routines.
!
! PUBLIC DATA:
!   MWlandC:  TELSEM2Atlas structure containing the microwave land surface
!             emissivity climatology atlas.
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public data
!       structure MWlandC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the CRTM
!       initialisation.
!

MODULE CRTM_MWlandCoeff

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE TELSEM2Atlas_Define   , ONLY: TELSEM2Atlas_type     , &
                                    TELSEM2Atlas_Associated, &
                                    TELSEM2Atlas_Destroy
  USE TELSEM2Atlas_netCDF_IO, ONLY: TELSEM2Atlas_netCDF_ReadFile
  USE TELSEM2_Atlas_Module  , ONLY: TELSEM2_Setup_Grid
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! The shared data
  PUBLIC :: MWlandC
  ! Procedures
  PUBLIC :: CRTM_MWlandCoeff_Load
  PUBLIC :: CRTM_MWlandCoeff_Destroy
  PUBLIC :: CRTM_MWlandCoeff_IsLoaded


  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PARAMETER :: ML = 512


  ! ----------------------------------------------------
  ! The shared microwave land surface emissivity atlas
  ! ----------------------------------------------------
  TYPE(TELSEM2Atlas_type), SAVE :: MWlandC


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_Load
!
! PURPOSE:
!       Function to load the TELSEM2 microwave land surface emissivity atlas
!       into the public data structure MWlandC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWlandCoeff_Load( &
!                        Filename,                              &
!                        File_Path         = File_Path        , &
!                        netCDF            = netCDF           , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Filename:           Name of the TELSEM2 atlas coefficient file.
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       File_Path:          File path for the input data file. Default is the
!                           current directory.
!       netCDF:             Present for interface consistency with the other
!                           coefficient loaders. The TELSEM2 atlas is only
!                           distributed in netCDF, so this is effectively always
!                           treated as netCDF.
!       Quiet:              Suppress INFORMATION messages if .TRUE.
!       Process_ID:         MPI process ID (message control only).
!       Output_Process_ID:  MPI process ID that emits messages.
!
! FUNCTION RESULT:
!       Error_Status:       SUCCESS or FAILURE.
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_Load( &
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    netCDF           , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL,      OPTIONAL, INTENT(IN) :: netCDF
    LOGICAL,      OPTIONAL, INTENT(IN) :: Quiet
    INTEGER,      OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER,      OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWlandCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(ML) :: MWlandCoeff_File
    LOGICAL :: noisy

    ! Setup
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    MWlandCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) MWlandCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(MWlandCoeff_File)
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

    ! Read the TELSEM2 atlas file (netCDF)
    err_stat = TELSEM2Atlas_netCDF_ReadFile( TRIM(MWlandCoeff_File), MWlandC )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading MWlandCoeff TELSEM2 atlas file '//TRIM(MWlandCoeff_File)//TRIM(pid_msg)
      CALL Load_Cleanup(); RETURN
    END IF

    ! Build the equal-area grid geometry and reverse lookup
    CALL TELSEM2_Setup_Grid( MWlandC )

    IF ( noisy ) THEN
      WRITE( msg,'("TELSEM2 MW land emissivity atlas loaded: ",i0," cells, ",i0," months")' ) &
             MWlandC%n_Data, MWlandC%n_Months
      CALL Display_Message( ROUTINE_NAME, TRIM(msg)//TRIM(pid_msg), INFORMATION )
    END IF

  CONTAINS

    SUBROUTINE Load_CleanUp()
      CALL TELSEM2Atlas_Destroy( MWlandC )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Load_CleanUp

  END FUNCTION CRTM_MWlandCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure MWlandC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_MWlandCoeff_Destroy( Process_ID = Process_ID )
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_MWlandCoeff_Destroy'
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
    CALL TELSEM2Atlas_Destroy( MWlandC )
    IF ( TELSEM2Atlas_Associated( MWlandC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating MWlandCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME, msg, err_stat ); RETURN
    END IF

  END FUNCTION CRTM_MWlandCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_MWlandCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if the TELSEM2 MW land emissivity atlas has been
!       loaded into the public data structure MWlandC.
!
!:sdoc-:
!------------------------------------------------------------------------------
  FUNCTION CRTM_MWlandCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = TELSEM2Atlas_Associated( MWlandC )
  END FUNCTION CRTM_MWlandCoeff_IsLoaded

END MODULE CRTM_MWlandCoeff
