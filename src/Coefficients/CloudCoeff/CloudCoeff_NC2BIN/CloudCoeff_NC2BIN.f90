!
! CloudCoeff_NC2BIN
!
! Program to convert a CRTM CloudCoeff data file
! from netCDF to Binary format
!
!
! CREATION HISTORY:
!       Written by: Paul van Delst, 27-Apr-2007
!                   paul.vandelst@noaa.gov
!

PROGRAM CloudCoeff_NC2BIN

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module usage
  USE Message_Handler   , ONLY: SUCCESS, FAILURE, Program_Message, Display_Message
!  USE SignalFile_Utility, ONLY: Create_SignalFile
  USE CloudCoeff_Define , ONLY: CloudCoeff_type
  USE CloudCoeff_IO     , ONLY: CloudCoeff_netCDF_to_Binary
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'CloudCoeff_NC2BIN'
  CHARACTER(*), PARAMETER :: PROGRAM_VERSION_ID = &
  '$Id: CloudCoeff_NC2BIN.f90 7080 2010-03-15 17:19:38Z david.groff@noaa.gov $'
  
  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(256) :: nc_filename, bin_filename, tmp_filename, msg
  INTEGER :: n_args
 
  ! Program header
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to convert a CRTM CloudCoeff data file '//&
                        'from netCDF to Binary format.', &
                        '$Revision: 7080 $')

  ! Get the filename
  n_args = COMMAND_ARGUMENT_COUNT()
  IF ( n_args > 0 ) THEN
     CALL GET_COMMAND_ARGUMENT(1, nc_filename)
     !** automatically generate the binary filename based on the command line netcdf filename
     tmp_filename = nc_filename(1:LEN_TRIM(nc_filename) - 3)//".bin"
     bin_filename = TRIM(ADJUSTL(tmp_filename))
     PRINT *, "Output filename:", bin_filename
  ELSE
     ! Get the filenames
     WRITE(*,FMT='(/5x,"Enter the INPUT netCDF CloudCoeff filename : ")', ADVANCE='NO')
     READ(*,'(a)') nc_filename
     nc_filename = ADJUSTL(nc_filename)
     WRITE(*,FMT='(/5x,"Enter the OUTPUT Binary CloudCoeff filename: ")', ADVANCE='NO')
     READ(*,'(a)') bin_filename
     bin_filename = ADJUSTL(bin_filename)
     ! ...Sanity check that they're not the same
     IF ( bin_filename == nc_filename ) THEN
        msg = 'CloudCoeff netCDF and Binary filenames are the same!'
        CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
     END IF
  END IF

  
  
  ! Perform the conversion
  err_stat = CloudCoeff_netCDF_to_Binary( NC_Filename, BIN_Filename )
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( PROGRAM_NAME, &
                          'CloudCoeff netCDF -> Binary conversion failed!', &
                          FAILURE )
    STOP
  END IF
  
  
!!$  ! Create a signal file indicating success
!!$  err_stat = Create_SignalFile( BIN_Filename )
!!$  IF ( err_stat /= SUCCESS ) THEN
!!$    CALL Display_Message( PROGRAM_NAME, &
!!$                          'Error creating signal file for '//TRIM(BIN_Filename), &
!!$                          FAILURE )
!!$  END IF
!!$
END PROGRAM CloudCoeff_NC2BIN

