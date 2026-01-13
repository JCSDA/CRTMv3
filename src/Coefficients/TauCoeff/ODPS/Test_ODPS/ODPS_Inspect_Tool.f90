PROGRAM ODPS_Inspect_Tool
  USE ODPS_Define, ONLY: ODPS_type, Info_ODPS, Destroy_ODPS
  USE ODPS_Binary_IO, ONLY: Read_ODPS_Binary
  USE Message_Handler, ONLY: SUCCESS, FAILURE, Display_Message
  IMPLICIT NONE

  CHARACTER(256) :: filename
  CHARACTER(2048) :: info_string
  TYPE(ODPS_type) :: ODPS
  INTEGER :: err_stat

  IF (COMMAND_ARGUMENT_COUNT() /= 1) THEN
     WRITE(*,*) "Usage: ODPS_Inspect_Tool <filename>"
     STOP
  END IF

  CALL GET_COMMAND_ARGUMENT(1, filename)

  err_stat = Read_ODPS_Binary( TRIM(filename), ODPS )
  IF (err_stat /= SUCCESS) THEN
     WRITE(*,*) "Error reading file: ", TRIM(filename)
     STOP
  END IF

  CALL Info_ODPS(ODPS, info_string)
  WRITE(*,*) TRIM(info_string)

  err_stat = Destroy_ODPS(ODPS)
  IF (err_stat /= SUCCESS) THEN
     WRITE(*,*) "Error destroying ODPS structure."
  END IF

END PROGRAM ODPS_Inspect_Tool
