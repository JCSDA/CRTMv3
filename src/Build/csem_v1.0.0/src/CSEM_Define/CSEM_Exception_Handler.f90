!
! CSEM_Exception_Handler
!
! This module is currently used to define simple error/exit codes and output messages, but will be
! expanded to include more functions to handle the CSEM exceptions. It is made up of several 
! subrotuines from the NOAA/NESDIS CRTM package
!
!
! CREATION HISTORY:
!       Modified by:    Ming Chen, 28-Feb-2015
!                       ming.chen@noaa.gov
!
!
MODULE CSEM_Exception_Handler

  ! Module use statements
  ! Disable all implicit typing
  IMPLICIT NONE

  ! Visibilities
  PRIVATE
  ! Module parameters
  PUBLIC :: SUCCESS    
  PUBLIC :: INFORMATION
  PUBLIC :: WARNING    
  PUBLIC :: FAILURE    
  PUBLIC :: EOF        
  PUBLIC :: UNDEFINED
  ! Module procedures  
  PUBLIC :: Display_Message
  PUBLIC :: Open_Message_Log

  ! Integer values that define the error or exit state.
  ! Note: These values are totally arbitrary.
  INTEGER, PARAMETER :: SUCCESS     = 0
  INTEGER, PARAMETER :: INFORMATION = 1
  INTEGER, PARAMETER :: WARNING     = 2
  INTEGER, PARAMETER :: FAILURE     = 3
  INTEGER, PARAMETER :: EOF         = 4
  INTEGER, PARAMETER :: UNDEFINED   = 5

  ! Character descriptors of the error states
  INTEGER,      PARAMETER :: MAX_N_STATES = 5
  CHARACTER(*), PARAMETER, DIMENSION( 0:MAX_N_STATES ) :: &
    STATE_DESCRIPTOR = (/ 'SUCCESS    ', &
                          'INFORMATION', &
                          'WARNING    ', &
                          'FAILURE    ', &
                          'END-OF-FILE', &
                          'UNDEFINED  ' /)


CONTAINS



  ! Subroutine to display messages.
  !
  ! This routine calls itself if the optional argument Message_Log 
  ! is passed and an error occurs opening the output log file.
  !
  RECURSIVE SUBROUTINE Display_Message(Routine_Name, &
                                       Message,      &
                                       Error_State,  &
                                       Message_Log   )
    ! Arguments
    CHARACTER(*), INTENT(IN)           :: Routine_Name
    CHARACTER(*), INTENT(IN)           :: Message
    INTEGER,      INTENT(IN)           :: Error_State
    CHARACTER(*), INTENT(IN), OPTIONAL :: Message_Log
    ! Local parameters
    CHARACTER(*), PARAMETER :: THIS_ROUTINE_NAME = 'Display_Message'
    CHARACTER(*), PARAMETER :: FMT_STRING = '( 1x, a, "(", a, ") : ", a )'
    ! Local variables
    INTEGER :: Error_State_To_Use
    LOGICAL :: Log_To_StdOut
    INTEGER :: File_ID
    INTEGER :: Error_Status

    ! Check the input error state
    Error_State_To_Use = Error_State
    IF ( Error_State < 0 .OR. Error_State > MAX_N_STATES ) THEN
      Error_State_To_Use = UNDEFINED
    END IF

    ! Set the message log. Default is output to stdout
    Log_To_StdOut = .TRUE.
    IF ( PRESENT( Message_Log ) ) THEN
      Log_To_StdOut = .FALSE.
      Error_Status = Open_Message_Log( TRIM( Message_Log ), File_ID )
      IF ( Error_Status /= 0 ) THEN
        CALL Display_Message( THIS_ROUTINE_NAME, &
                              'Error opening message log file', &
                              FAILURE )
        Log_To_StdOut = .TRUE.
      END IF
    END IF

    ! Output the message
    IF ( Log_To_StdOut ) THEN
      WRITE( *, FMT = FMT_STRING ) &
                TRIM( Routine_Name ), &
                TRIM( STATE_DESCRIPTOR( Error_State_To_Use ) ), &
                TRIM( Message )
    ELSE
      WRITE( File_ID, FMT = FMT_STRING ) &
                      TRIM( Routine_Name ), &
                      TRIM( STATE_DESCRIPTOR( Error_State_To_Use ) ), &
                      TRIM( Message )
      CLOSE( File_ID )
    END IF

  END SUBROUTINE Display_Message

  ! Function to open the message log file.
  !
  ! SIDE EFFECTS:
  !   The file is opened for SEQUENTIAL, FORMATTED access with
  !   UNKNOWN status, position of APPEND, and action of READWRITE.
  !
  !   Hopefully all of these options will not cause an existing file
  !   to be inadvertantly overwritten.
  !
  FUNCTION Open_Message_Log(Message_Log, File_ID) RESULT(Error_Status)
    ! Arguments
    CHARACTER(*), INTENT(IN)  :: Message_Log
    INTEGER,      INTENT(OUT) :: File_ID
    ! Function result
    INTEGER :: Error_Status
    ! Local variables
    INTEGER :: Lun
    INTEGER :: IO_Status

    ! Set successful return status
    Error_Status = SUCCESS

    ! Get a file unit number
    Lun = Get_Lun()
    IF ( Lun < 0 ) THEN
      Error_Status = FAILURE
      RETURN
    END IF

    ! Open the file
    OPEN( Lun, FILE     = TRIM( Message_Log ), &
               ACCESS   = 'SEQUENTIAL', &
               FORM     = 'FORMATTED', &
               STATUS   = 'UNKNOWN', &
               POSITION = 'APPEND', &
               ACTION   = 'READWRITE', &
               IOSTAT   = IO_Status )
    IF ( IO_Status /= 0 ) THEN
      Error_Status = FAILURE
      RETURN
    END IF

    ! Return the file ID
    File_ID = Lun

  END FUNCTION Open_Message_Log

!
! Get_Lun
!
! Function to obtain a free logical unit number for file access
!
! CALLING SEQUENCE:
!       Lun = Get_Lun()
!
! FUNCTION RESULT:
!       Lun:   Logical unit number that may be used for file access.
!              If Lun > 0 it can be used as a logical unit number to open
!                         and access a file.
!                 Lun < 0 a non-existant logical unit number was reached
!                         during the search.
!              UNITS:      N/A
!              TYPE:       INTEGER
!              DIMENSION:  Scalar
!

  FUNCTION Get_Lun() RESULT( Lun )
    INTEGER :: Lun
    LOGICAL :: Exist_Open_Unit

    ! Initialise logical unit number
    Lun = 9

    ! Start open loop for Lun Search
    Lun_Search: DO
      Lun = Lun + 1
      INQUIRE( UNIT = Lun, EXIST = Exist_Open_Unit )
      IF ( .NOT. Exist_Open_Unit ) THEN
        Lun = -1
        EXIT Lun_Search
      END IF
      INQUIRE( UNIT = Lun, OPENED = Exist_Open_Unit )
      IF ( .NOT. Exist_Open_Unit) EXIT Lun_Search
    END DO Lun_Search

  END FUNCTION Get_Lun

END MODULE CSEM_Exception_Handler
