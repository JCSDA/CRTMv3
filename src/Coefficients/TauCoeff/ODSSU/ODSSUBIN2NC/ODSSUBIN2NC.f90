!
! ODSSUBIN2NC
!
! Program to convert a CRTM ODSSU TauCoeff file from Binary to netCDF format.
!
! Usage:
!   ODSSUBIN2NC <input.TauCoeff.bin> [output.TauCoeff.nc]
!
! If the output filename is omitted, it is derived from the input by
! replacing the trailing ".bin" with ".nc".
!

PROGRAM ODSSUBIN2NC

  USE File_Utility    , ONLY: File_Exists
  USE Message_Handler , ONLY: SUCCESS, FAILURE, INFORMATION, &
                              Program_Message, Display_Message
  USE ODSSU_Define    , ONLY: ODSSU_type, Destroy_ODSSU
  USE ODSSU_Binary_IO , ONLY: Read_ODSSU_Binary
  USE ODSSU_netCDF_IO , ONLY: Write_ODSSU_netCDF

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'ODSSUBIN2NC'

  INTEGER :: err_stat, n_args
  CHARACTER(512) :: bin_filename
  CHARACTER(512) :: nc_filename
  CHARACTER(256) :: msg
  TYPE(ODSSU_type) :: ODSSU

  CALL Program_Message( PROGRAM_NAME, &
                        'Convert a CRTM ODSSU TauCoeff file from Binary to netCDF.', &
                        'CRTM v3 REL-3.2.0' )

  n_args = COMMAND_ARGUMENT_COUNT()
  IF ( n_args < 1 ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Usage: ODSSUBIN2NC <input.TauCoeff.bin> [output.TauCoeff.nc]', FAILURE )
    STOP 1
  END IF

  CALL GET_COMMAND_ARGUMENT(1, bin_filename)
  bin_filename = ADJUSTL(bin_filename)

  IF ( n_args >= 2 ) THEN
    CALL GET_COMMAND_ARGUMENT(2, nc_filename)
    nc_filename = ADJUSTL(nc_filename)
  ELSE
    ! Derive: strip trailing ".bin" (4 chars) and append ".nc"
    IF ( LEN_TRIM(bin_filename) > 4 .AND. &
         bin_filename(LEN_TRIM(bin_filename)-3:LEN_TRIM(bin_filename)) == '.bin' ) THEN
      nc_filename = bin_filename(1:LEN_TRIM(bin_filename)-4) // '.nc'
    ELSE
      nc_filename = TRIM(bin_filename) // '.nc'
    END IF
  END IF

  IF ( TRIM(bin_filename) == TRIM(nc_filename) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Input and output filenames are the same.', FAILURE )
    STOP 1
  END IF

  IF ( .NOT. File_Exists( TRIM(bin_filename) ) ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Input file '//TRIM(bin_filename)//' not found.', FAILURE )
    STOP 1
  END IF

  WRITE(*,'(/5x,"Reading Binary ODSSU file ",a)') TRIM(bin_filename)
  err_stat = Read_ODSSU_Binary( TRIM(bin_filename), ODSSU )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Read_ODSSU_Binary failed for '//TRIM(bin_filename), FAILURE )
    STOP 1
  END IF

  WRITE(*,'(/5x,"Writing netCDF ODSSU file ",a)') TRIM(nc_filename)
  err_stat = Write_ODSSU_netCDF( TRIM(nc_filename), ODSSU )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Write_ODSSU_netCDF failed for '//TRIM(nc_filename), FAILURE )
    STOP 1
  END IF

  err_stat = Destroy_ODSSU( ODSSU )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
      'Destroy_ODSSU failed (non-fatal).', INFORMATION )
  END IF

  msg = 'ODSSU Binary -> netCDF conversion successful: '//TRIM(nc_filename)
  CALL Display_Message( PROGRAM_NAME, TRIM(msg), INFORMATION )

END PROGRAM ODSSUBIN2NC
