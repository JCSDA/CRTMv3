!
! CloudCoeff_BIN2NC
!
! Program to convert Binary format CloudCoeff files to the netCDF format.
!
!
! CREATION HISTORY:
!       Written by:     Benjamin Johnson
!                       benjamin.t.johnson@noaa.gov / bjohns@ucar.edu
!       Updated:        Benjamin Johnson 8-9-2024 (bjohns@ucar.edu)
!                       Converted from CloudCoeff_BIN2NC.f90
!

PROGRAM CloudCoeff_BIN2NC

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage
  USE CRTM_Module
  USE CloudCoeff_Define   , ONLY: CloudCoeff_type
  USE CloudCoeff_IO       , ONLY: CloudCoeff_Binary_to_netCDF
  USE Endian_Utility    , ONLY: Big_Endian
  ! Disable implicit typing
  IMPLICIT NONE

  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'CloudCoeff_BIN2NC'

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: err_stat
  CHARACTER(256) :: msg
  CHARACTER(256) :: bin_filename
  CHARACTER(256) :: nc_filename, tmp_filename
  CHARACTER(256) :: answer
  CHARACTER(256) :: version_str
  INTEGER :: n_args
  
  ! Program header
  CALL CRTM_Version(version_str)
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to convert a CRTM CloudCoeff data file from Binary to netCDF format.', &
                        'CRTM Version: '//TRIM(version_str) )

  ! Get the filename                                                                                                                                                                                   ! Check for command line argument
  n_args = COMMAND_ARGUMENT_COUNT()
  IF ( n_args > 0 ) THEN
     CALL GET_COMMAND_ARGUMENT(1, BIN_filename)
     tmp_filename = BIN_filename(1:LEN_TRIM(bin_filename) - 4)//".nc"
     nc_filename = TRIM(ADJUSTL(tmp_filename))
     PRINT *, "nc_filename:", nc_filename
  ELSE
     ! Get the filenames
     IF ( Big_Endian() ) THEN
        WRITE(*,FMT='(/5x,"Enter the INPUT Binary [Big Endian] CloudCoeff filename : ")', ADVANCE='NO')
     ELSE
        WRITE(*,FMT='(/5x,"Enter the INPUT Binary [Little Endian] CloudCoeff filename : ")', ADVANCE='NO')
     END IF
     READ(*,'(a)') bin_filename
     bin_filename = ADJUSTL(bin_filename)
     WRITE(*,FMT='(/5x,"Enter the OUTPUT netCDF CloudCoeff filename: ")', ADVANCE='NO')
     READ(*,'(a)') nc_filename
     nc_filename = ADJUSTL(nc_filename)
  END IF



  
  ! ...Sanity check that they're not the same
  IF ( bin_filename == nc_filename ) THEN
    msg = 'CloudCoeff netCDF and Binary filenames are the same!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  END IF  

  ! Perform the conversion
  err_stat = CloudCoeff_Binary_to_netCDF( BIN_filename, NC_filename, Quiet=.FALSE.)
  IF ( err_stat /= SUCCESS ) THEN
    msg = 'CloudCoeff Binary -> netCDF conversion failed!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  ELSE
    msg = 'CloudCoeff Binary -> netCDF conversion successful!'
    CALL Display_Message( PROGRAM_NAME, msg, err_stat )
  END IF
  
END PROGRAM CloudCoeff_BIN2NC
