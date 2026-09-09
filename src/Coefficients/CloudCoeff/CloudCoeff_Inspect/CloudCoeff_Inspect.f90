!
! CloudCoeff_Inspect
!
! Program to inspect the contents of a CRTM CloudCoeff file (binary or netCDF).
!
! The file format is selected automatically from the filename: any file whose
! name contains ".nc" is read as netCDF, otherwise it is read as Binary.
!
! For netCDF files the schema is also auto-detected: a file carrying the
! 'CRTM-Exp' scheme (or an n_Habit dimension) is read with the experimental
! CloudCoeff reader and displayed via CloudCoeff_Exp_Inspect; all other files
! use the standard CloudCoeff reader/inspector.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 20-Jun-2006
!                       paul.vandelst@noaa.gov
!
!       Modified by:    Benjamin Johnson, 05-Jun-2026
!                       Added netCDF support via CloudCoeff_IO, plus auto-detect
!                       and display of the experimental ('CRTM-Exp') schema.
!

PROGRAM CloudCoeff_Inspect

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE File_Utility,             ONLY: File_Exists
  USE Message_Handler,          ONLY: SUCCESS, FAILURE, Program_Message, Display_Message
  USE CloudCoeff_Define,        ONLY: CloudCoeff_type, CloudCoeff_Destroy, &
                                      Inspect => CloudCoeff_Inspect
  USE CloudCoeff_IO,            ONLY: CloudCoeff_ReadFile
  ! ...Experimental ('CRTM-Exp') schema support
  USE CloudCoeff_Exp_Define,    ONLY: CloudCoeff_Exp_type, CloudCoeff_Exp_Destroy, &
                                      Exp_Inspect => CloudCoeff_Exp_Inspect
  USE CloudCoeff_Exp_netCDF_IO, ONLY: CloudCoeff_Exp_netCDF_ReadFile, &
                                      CloudCoeff_Exp_netCDF_InquireFile
  ! Disable implicit typing
  IMPLICIT NONE


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'CloudCoeff_Inspect'
  CHARACTER(*), PARAMETER :: PROGRAM_RCS_ID = ''


  ! ---------
  ! Variables
  ! ---------
  CHARACTER(256) :: msg
  CHARACTER(256) :: filename, arg_string
  CHARACTER(32)  :: scheme
  INTEGER :: err_stat
  INTEGER :: n, n_cmd_args
  INTEGER :: n_habit
  LOGICAL :: is_nc, is_exp
  LOGICAL :: pause
  TYPE(CloudCoeff_type)     :: coeffs
  TYPE(CloudCoeff_Exp_type) :: exp_coeffs


  ! Output program header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to display the contents of a CRTM '//&
                        'Binary/netCDF format CloudCoeff file (standard or '//&
                        'experimental CRTM-Exp schema) to stdout.', &
                        '$Revision$' )


  ! Parse the command line arguments.
  ! ...The first non-"pause" argument is taken to be the filename; the keyword
  !    "pause" (case sensitive) enables paging between sections.
  filename = ''
  pause    = .FALSE.
  n_cmd_args = COMMAND_ARGUMENT_COUNT()
  DO n = 1, n_cmd_args
    CALL GET_COMMAND_ARGUMENT(n, arg_string)
    arg_string = ADJUSTL(arg_string)
    IF ( TRIM(arg_string) == 'pause' ) THEN
      pause = .TRUE.
    ELSE IF ( LEN_TRIM(filename) == 0 ) THEN
      filename = arg_string
    END IF
  END DO


  ! Prompt for the filename if it was not supplied on the command line
  IF ( LEN_TRIM(filename) == 0 ) THEN
    WRITE( *,FMT='(/5x,"Enter the CloudCoeff filename: ")',ADVANCE='NO' )
    READ( *,'(a)' ) filename
  END IF
  filename = ADJUSTL(filename)
  IF ( .NOT. File_Exists( TRIM(filename) ) ) THEN
    msg = 'File '//TRIM(filename)//' not found.'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF


  ! Select the file format from the filename: ".nc" (or ".nc4") => netCDF
  is_nc = ( INDEX(TRIM(filename), '.nc', BACK=.TRUE.) > 0 )


  ! For netCDF files, probe for the experimental ('CRTM-Exp') schema. The probe
  ! is best-effort: a standard CloudCoeff file has neither the 'CRTM-Exp' scheme
  ! attribute nor an n_Habit dimension, so it falls through to the standard path.
  is_exp = .FALSE.
  IF ( is_nc ) THEN
    scheme  = ''
    n_habit = 0
    err_stat = CloudCoeff_Exp_netCDF_InquireFile( TRIM(filename), &
                                                  n_Habit = n_habit, &
                                                  Scheme  = scheme )
    IF ( err_stat == SUCCESS ) &
      is_exp = ( TRIM(scheme) == 'CRTM-Exp' ) .OR. ( n_habit > 0 )
  END IF


  ! Read and display the contents using the appropriate reader/inspector
  IF ( is_exp ) THEN
    ! ...Experimental 'CRTM-Exp' schema
    err_stat = CloudCoeff_Exp_netCDF_ReadFile( TRIM(filename), exp_coeffs )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading experimental CloudCoeff file '//TRIM(filename)
      CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
    END IF
    CALL Exp_Inspect( exp_coeffs )
    CALL CloudCoeff_Exp_Destroy( exp_coeffs )
  ELSE
    ! ...Standard schema (binary or netCDF)
    err_stat = CloudCoeff_ReadFile( filename, coeffs, netCDF=is_nc )
    IF ( err_stat /= SUCCESS ) THEN
      msg = 'Error reading CloudCoeff file '//TRIM(filename)
      CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
    END IF
    CALL Inspect( coeffs, Pause=pause )
    CALL CloudCoeff_Destroy( coeffs )
  END IF

END PROGRAM CloudCoeff_Inspect
