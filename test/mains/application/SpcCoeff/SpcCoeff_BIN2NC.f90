!
! SpcCoeff_BIN2NC
!
! Program to convert netCDF format SpcCoeff files to the Binary format.
!
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 27-Jul-2002
!                       paul.vandelst@noaa.gov
!       Updated:        Benjamin Johnson 1-22-2024 (bjohns@ucar.edu)
!                       Converted from SpcCoeff_NC2BIN.f90, includes ACCoeff and NLTE support
!

PROGRAM SpcCoeff_BIN2NC

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE Message_Handler   , ONLY: SUCCESS, FAILURE, Program_Message, Display_Message
  USE String_Utility    , ONLY: StrLowCase
  USE SignalFile_Utility, ONLY: Create_SignalFile
  USE SpcCoeff_Define   , ONLY: SpcCoeff_type
  USE SpcCoeff_IO       , ONLY: SpcCoeff_Binary_to_netCDF
  ! Disable implicit typing
  IMPLICIT NONE


  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'SpcCoeff_BIN2NC'

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(256) :: msg
  CHARACTER(256) :: NC_filename
  CHARACTER(256) :: BIN_filename
  CHARACTER(256) :: answer
  INTEGER :: version

  ! Program header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to convert a CRTM SpcCoeff data file '//&
                        'from netCDF to Binary format.')


  ! Get the filenames
  WRITE(*,FMT='(/5x,"Enter the INPUT Binary SpcCoeff filename : ")', ADVANCE='NO')
  READ(*,'(a)') BIN_filename
  BIN_filename = ADJUSTL(BIN_filename)
  WRITE(*,FMT='(/5x,"Enter the OUTPUT netCDF SpcCoeff filename: ")', ADVANCE='NO')
  READ(*,'(a)') NC_filename
  NC_filename = ADJUSTL(NC_filename)

  ! ...Sanity check that they're not the same
  IF ( bin_filename == nc_filename ) THEN
    msg = 'SpcCoeff netCDF and Binary filenames are the same!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF

  ! Ask if version increment required
  WRITE(*,FMT='(/5x,"Increment the OUTPUT version number? [y/n]: ")', ADVANCE='NO')
  READ(*,'(a)') answer
  answer = StrLowCase(ADJUSTL(answer))
  SELECT CASE( TRIM(answer) )
    CASE('y','yes')
      version = -1
    CASE DEFAULT
      version = 0
  END SELECT
  

  ! Perform the conversion
  err_stat = SpcCoeff_Binary_to_netCDF( BIN_filename, NC_filename, Version = version )
  IF ( err_stat /= SUCCESS ) THEN
    msg = 'SpcCoeff Binary -> netCDF conversion failed!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  ELSE
    msg = 'SpcCoeff Binary -> netCDF conversion successful!'
    CALL Display_Message( PROGRAM_NAME, msg, err_stat )
  END IF
  
  
  ! Create a signal file indicating success
  err_stat = Create_SignalFile( bin_filename )
  IF ( err_stat /= SUCCESS ) THEN
    msg = 'Error creating signal file for '//TRIM(bin_filename)
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE )
  END IF

END PROGRAM SpcCoeff_NC2BIN
