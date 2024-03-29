!
! CRTM_AerosolCoeff
!
! Module containing the shared CRTM aerosol coefficients (AerosolCoeff)
! and their load/destruction routines.
!
! PUBLIC DATA:
!       AeroC:  Data structure containing the aerosol bulk optical
!               properties data
!
! SIDE EFFECTS:
!       Routines in this module modify the contents of the public
!       data structure AeroC.
!
! RESTRICTIONS:
!       Routines in this module should only be called during the
!       CRTM initialisation.
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, CIMSS/SSEC 24-Jun-2004
!                       paul.vandelst@noaa.gov
!       Modified by     Yingtao Ma, 2020/6/11
!                       yingtao.ma@noaa.gov
!                       Implemented CMAQ aerosol
!

MODULE CRTM_AerosolCoeff

  ! ----------------
  ! Enviroment setup
  ! ----------------
  ! Module use
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, Display_Message
  USE AerosolCoeff_Define   , ONLY: AerosolCoeff_type, &
                                    AerosolCoeff_Associated, &
                                    AerosolCoeff_Destroy
  USE AerosolCoeff_IO       , ONLY: AerosolCoeff_ReadFile
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE

  ! The shared data
  PUBLIC :: AeroC

  ! Public routines defined in this module
  PUBLIC :: CRTM_AerosolCoeff_Load
  PUBLIC :: CRTM_AerosolCoeff_Destroy
  PUBLIC :: CRTM_AerosolCoeff_IsLoaded

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Message string length
  INTEGER, PARAMETER :: ML = 256


  ! -----------------------------------
  ! The shared aerosol coefficient data
  ! -----------------------------------
  TYPE(AerosolCoeff_type), TARGET, SAVE :: AeroC


CONTAINS


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AerosolCoeff_Load
!
! PURPOSE:
!       Function to load the AerosolCoeff scattering coefficient data into
!       the public data structure AerosolC.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AerosolCoeff_Load( &
!                        Aerosol_Model , &
!                        Filename , &
!                        File_Path         = File_Path        , &
!                        netCDF            = netCDF           , &
!                        Quiet             = Quiet            , &
!                        Process_ID        = Process_ID       , &
!                        Output_Process_ID = Output_Process_ID  )
!
! INPUT ARGUMENTS:
!       Aerosol_Model:      Name of the aerosol scheme for scattering calculation
!                           Available aerosol scheme:
!                           - CRTM  [DEFAULT]
!                           - CMAQ
!                           - GOCART-GEOS5
!                           - NAAPS
!                           UNITS:      N/A
!                           TYPE:       CHARACTER(*)
!                           DIMENSION:  Scalar
!                           ATTRIBUTES: INTENT(IN), OPTIONAL
!
!       Filename:           Name of the data file containing the aerosol optical
!                           properties data for scattering calculations.
!                           Available datafiles:
!                           CRTM:
!                           - AerosolCoeff.bin  [DEFAULT, Binary]
!                           - AerosolCoeff.nc   [netCDF-Classic/4]
!                           CMAQ:
!                           - AerosolCoeff.CMAQ.bin  [Binary]
!                           - AerosolCoeff.CMAQ.nc   [netCDF-Classic/4]
!                           GOCART-GEOS5:
!                           - AerosolCoeff.GOCART-GEOS5.bin  [Binary]
!                           - AerosolCoeff.GOCART-GEOS5.nc   [netCDF-Classic/4]
!                           NAAPS:
!                           - AerosolCoeff.NAAPS.bin  [Binary]
!                           - AerosolCoeff.NAAPS.nc   [netCDF-Classic/4]
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
!                           If == SUCCESS the AerosolCoeff data load was successful
!                              == FAILURE an unrecoverable error occurred.
!                           UNITS:      N/A
!                           TYPE:       INTEGER
!                           DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure AerosolC.
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_AerosolCoeff_Load( &
    Aerosol_Model    , &  ! Input
    Filename         , &  ! Input
    File_Path        , &  ! Optional input
    netCDF           , &  ! Optional input
    Quiet            , &  ! Optional input
    Process_ID       , &  ! Optional input
    Output_Process_ID) &  ! Optional input
  RESULT( err_stat )
    ! Arguments
    CHARACTER(*),           INTENT(IN) :: Aerosol_Model
    CHARACTER(*),           INTENT(IN) :: Filename
    CHARACTER(*), OPTIONAL, INTENT(IN) :: File_Path
    LOGICAL,      OPTIONAL, INTENT(IN) :: netCDF
    LOGICAL     , OPTIONAL, INTENT(IN) :: Quiet
    INTEGER     , OPTIONAL, INTENT(IN) :: Process_ID
    INTEGER     , OPTIONAL, INTENT(IN) :: Output_Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AerosolCoeff_Load'
    ! Local variables
    CHARACTER(ML) :: msg, pid_msg
    CHARACTER(ML) :: AerosolCoeff_File
    LOGICAL :: noisy
    ! Function variables
    LOGICAL :: Binary

    ! Setup
    err_stat = SUCCESS
    ! ...Assign the filename to local variable
    AerosolCoeff_File = ADJUSTL(Filename)
    ! ...Add the file path
    IF ( PRESENT(File_Path) ) AerosolCoeff_File = TRIM(ADJUSTL(File_Path))//TRIM(AerosolCoeff_File)
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


    ! Read the AerosolCoeff file
    err_stat = AerosolCoeff_ReadFile( &
                 Aerosol_Model, &
                 AerosolCoeff_File, &
                 AeroC, &
                 netCDF = .NOT. Binary, &
                 Quiet  = .NOT. noisy )
     IF ( err_stat /= SUCCESS ) THEN
       WRITE( msg,'("Error reading AerosolCoeff file ",a)') TRIM(AerosolCoeff_File)
       CALL Display_Message( ROUTINE_NAME,TRIM(msg)//TRIM(pid_msg),err_stat )
       RETURN
     END IF

   CONTAINS

     SUBROUTINE Load_CleanUp()
       CALL AerosolCoeff_Destroy( AeroC )
       err_stat = FAILURE
       CALL Display_Message( ROUTINE_NAME, msg, err_stat )
     END SUBROUTINE Load_CleanUp

  END FUNCTION CRTM_AerosolCoeff_Load


!------------------------------------------------------------------------------
!
! NAME:
!       CRTM_AerosolCoeff_Destroy
!
! PURPOSE:
!       Function to deallocate the public data structure AeroC containing
!       the CRTM AerosolCoeff aerosol coefficient data.
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_AerosolCoeff_Destroy( Process_ID = Process_ID )
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
!                         If == SUCCESS the deallocation of the public AeroC data
!                                       structure was successful
!                            == FAILURE an unrecoverable error occurred.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!
! SIDE EFFECTS:
!       This function modifies the contents of the public data structure AeroC.
!
!------------------------------------------------------------------------------

  FUNCTION CRTM_AerosolCoeff_Destroy( Process_ID ) RESULT( err_stat )
    ! Arguments
    INTEGER, OPTIONAL, INTENT(IN)  :: Process_ID
    ! Function result
    INTEGER :: err_stat
    ! Local parameters
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_AerosolCoeff_Destroy'
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
    CALL AerosolCoeff_Destroy( AeroC )
    IF ( AerosolCoeff_Associated( AeroC ) ) THEN
      err_stat = FAILURE
      msg = 'Error deallocating AerosolCoeff shared data structure'//TRIM(pid_msg)
      CALL Display_Message( ROUTINE_NAME,msg,err_stat )
      RETURN
    END IF

  END FUNCTION CRTM_AerosolCoeff_Destroy


!------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_AerosolCoeff_IsLoaded
!
! PURPOSE:
!       Function to test if the AerosolCoeff scattering coefficient data has
!       loaded into the public data structure AerosolC.
!
! CALLING SEQUENCE:
!       status = CRTM_AerosolCoeff_IsLoaded()
!
!:sdoc-:
!------------------------------------------------------------------------------

  FUNCTION CRTM_AerosolCoeff_IsLoaded() RESULT( IsLoaded )
    LOGICAL :: IsLoaded
    IsLoaded = AerosolCoeff_Associated( AeroC )
  END FUNCTION CRTM_AerosolCoeff_IsLoaded



END MODULE CRTM_AerosolCoeff
