!
! CRTM_CloudCoeff
!
! Module containing the shared CRTM scattering coefficient data
! (CloudCoeff) and their load/destruction routines.
!
! PUBLIC DATA:
!       CloudC:  Data structure containing the cloud bulk optical
!                properties data
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure CloudC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CRTM initialisation.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 24-Jun-2004
!                       paul.vandelst@noaa.gov
!

MODULE CRTM_CloudCoeff

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use
  USE Message_Handler,      ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE CloudCoeff_Define,    ONLY: CloudCoeff_type, &
                                  CloudCoeff_Associated, &
                                  CloudCoeff_Destroy, &
                                  INVALID_CLOUDCOEFF, &
                                  MIE_TAMU_CLOUDCOEFF, &
                                  DDA_ARTS_CLOUDCOEFF
  USE CloudCoeff_IO       , ONLY: CloudCoeff_ReadFile
  ! Experimental ('CRTM-Exp') opt-in scheme
  USE CloudCoeff_Exp_Define   , ONLY: CloudCoeff_Exp_type, &
                                      CloudCoeff_Exp_Destroy, &
                                      CloudCoeff_Exp_Associated, &
                                      CRTM_EXP_CLOUDCOEFF
  USE CloudCoeff_Exp_netCDF_IO, ONLY: CloudCoeff_Exp_netCDF_ReadFile
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared variables
  PUBLIC :: INVALID_CLOUDCOEFF
  PUBLIC :: MIE_TAMU_CLOUDCOEFF
  PUBLIC :: DDA_ARTS_CLOUDCOEFF

  ! The shared data
  PUBLIC :: CloudC
  PUBLIC :: CloudC_Exp
  PUBLIC :: Active_Cloud_Scheme
  PUBLIC :: CRTM_EXP_CLOUDCOEFF
  PUBLIC :: SCHEME_LEGACY
  ! Procedures
  PUBLIC :: CRTM_CloudCoeff_Load
  PUBLIC :: CRTM_CloudCoeff_Destroy
  PUBLIC :: CRTM_CloudCoeff_IsLoaded


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


  ! Active cloud-optics scheme selector (set by CRTM_CloudCoeff_Load).
  ! Legacy (MIE_TAMU / DDA_ARTS, auto-detected from CloudC) is the default.
  INTEGER, PARAMETER :: SCHEME_LEGACY = 0


  ! ---------------------------------
  ! The shared cloud coefficient data
  ! ---------------------------------
  TYPE(CloudCoeff_type), TARGET, SAVE :: CloudC
  ! Experimental 'CRTM-Exp' shared data + which scheme is active
  TYPE(CloudCoeff_Exp_type), TARGET, SAVE :: CloudC_Exp
  INTEGER, SAVE :: Active_Cloud_Scheme = SCHEME_LEGACY


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_Load
!
! PURPOSE:
!       Function to load the CloudCoeff scattering coefficient data into
!       the public data structure CloudC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_CloudCoeff_Load( &
!                        Cloud_Model , &
!                        Filename , &
!                        File_Path         = File_Path        , &
!                        netCDF            = netCDF           , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Cloud_Model:        Name of the cloud scheme for scattering calculation
!                           Available cloud scheme:
!                           - CRTM  [DEFAULT]
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       CloudCoeff_IO:      Format of the cloud optical properties data
!                           Available options:
!                           - Binary  [DEFAULT]
!                           - netCDF
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Filename:           Name of the data file containing the cloud optical
!                           properties data for scattering calculations.
!                           Available datafiles:
!                           - CloudCoeff.bin  [DEFAULT, Binary]
!                           - CloudCoeff.nc   [netCDF-Classic/4]
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
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
!       netCDF:             Set this logical argument to specify file format.
!                           If == .FALSE., Binary [DEFAULT].
!                              == .TRUE.,  netCDF
!                           If not specified, default is .FALSE.
!                           UNITS:      N/A
!                           TYPE:       LOGICAL
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
!                           If == SUCCESS the CloudCoeff data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure CloudC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_Load( &
    Cloud_Model      , &  ! Input
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    netCDF           , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Cloud_Model
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL,      OPTIONAL, INTENT(IN) :: netCDF
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_CloudCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(ML) :: CloudCoeff_File
    LOGICAL :: noisy
    ! Function variables
    LOGICAL :: Binary

    ! Setup
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    CloudCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) CloudCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(CloudCoeff_File)
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
    ! ...Check netCDF argument
    Binary = .TRUE.
    IF ( PRESENT(netCDF) ) Binary = .NOT. netCDF


    ! Read the CloudCoeff data file.
    ! The experimental scheme is selected EXPLICITLY via Cloud_Model=='CRTM-Exp'
    ! (never auto-detected from file contents) so a black-box user cannot trip
    ! into it by swapping a coefficient file. Default => legacy, unchanged.
    IF ( TRIM(ADJUSTL(Cloud_Model)) == 'CRTM-Exp' ) THEN
      err_stat = CloudCoeff_Exp_netCDF_ReadFile( &
                   CloudCoeff_File, &
                   CloudC_Exp, &
                   Quiet = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading experimental CloudCoeff file ",a)') TRIM(CloudCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      Active_Cloud_Scheme = CRTM_EXP_CLOUDCOEFF
      ! Mirror the phase-element count onto the (otherwise empty) legacy CloudC
      ! scalar so the existing AtmOptics/CSvar allocations (sized from
      ! CloudC%n_Phase_Elements) are correct for the experimental scheme. The
      ! legacy CloudC arrays remain unallocated and unused.
      CloudC%n_Phase_Elements = CloudC_Exp%n_Phase_Elements
      IF ( noisy ) CALL Display_Message( ROUTINE_NAME, &
        'Active cloud-optics scheme: CRTM-Exp (experimental)'//TRIM(pid_msg), INFORMATION )
    ELSE
      err_stat = CloudCoeff_ReadFile( &
                   CloudCoeff_File, &
                   CloudC, &
                   netCDF = .NOT. Binary, &
                   Quiet  = .NOT. noisy )
      IF ( err_stat /= SUCCESS ) THEN
        WRITE( msg,'("Error reading CloudCoeff file ",a)') TRIM(CloudCoeff_File)
        CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
        RETURN
      END IF
      Active_Cloud_Scheme = SCHEME_LEGACY
    END IF

     ! Derive the cloud-optics scheme from the loaded data, once, so the MW cloud-scatter dispatch is
     ! explicit. The Mie-TAMU tables carry a positive MW effective-radius axis (Reff_MW); the DDA-ARTS
     ! database has none (it interpolates on water content), so Reff_MW is left zero on read.
     IF ( ALL(CloudC%Reff_MW > 0.0) ) THEN
       CloudC%Data_Type = MIE_TAMU_CLOUDCOEFF
     ELSE
       CloudC%Data_Type = DDA_ARTS_CLOUDCOEFF
     END IF

  CONTAINS

    SUBROUTINE Load_CleanUp()
      CALL CloudCoeff_Destroy( CloudC )
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, msg, err_stat )
    END SUBROUTINE Load_CleanUp

  END FUNCTION CRTM_CloudCoeff_Load


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure CloudC containing
!       the CRTM CloudCoeff scattering coefficient data.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_CloudCoeff_Destroy( Process_ID =Process_ID )
!
! OPTIONAL INPUT ARGUMENTS:
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
!                         If == SUCCESS the deallocation of the public CloudC data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure CloudC.
!
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN) :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_CloudCoeff_Destroy'
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

    ! Destroy the structures (both schemes) and reset the active-scheme flag
    CALL CloudCoeff_Exp_Destroy( CloudC_Exp )
    Active_Cloud_Scheme = SCHEME_LEGACY
    CALL CloudCoeff_Destroy( CloudC )
    IF ( CloudCoeff_Associated( CloudC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating CloudCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
      RETURN
    END IF

  END FUNCTION CRTM_CloudCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CloudCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if the CloudCoeff scattering coefficient data has
!       loaded into the public data structure CloudC.
!
! CALLING SEQUENCE:
!       status = CRTM_CloudCoeff_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_CloudCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    ! Loaded if EITHER the legacy or the experimental scheme data is present
    IsLoaded = CloudCoeff_Associated( CloudC ) .OR. &
               CloudCoeff_Exp_Associated( CloudC_Exp )
  END FUNCTION CRTM_CloudCoeff_IsLoaded

END MODULE CRTM_CloudCoeff
