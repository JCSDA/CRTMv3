!
! SpcCoeff_BIN2NC
!
! Program to convert Binary format SpcCoeff files to the netCDF format.
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
  USE CRTM_Module
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
  CHARACTER(256) :: BIN_filename
  CHARACTER(256) :: NC_filename
  CHARACTER(256) :: answer
  CHARACTER(256) :: version_str
  
  ! Program header
  CALL CRTM_Version(version_str)
  CALL Program_Message( PROGRAM_NAME, &
                        'Program to convert a CRTM SpcCoeff data file from Binary to netCDF format.', &
                        'CRTM Version: '//TRIM(version_str) )

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

  ! Perform the conversion
  err_stat = SpcCoeff_Binary_to_netCDF( BIN_filename, NC_filename, Quiet=.FALSE.)
  IF ( err_stat /= SUCCESS ) THEN
    msg = 'SpcCoeff Binary -> netCDF conversion failed!'
    CALL Display_Message( PROGRAM_NAME, msg, FAILURE ); STOP
  ELSE
    msg = 'SpcCoeff Binary -> netCDF conversion successful!'
    CALL Display_Message( PROGRAM_NAME, msg, err_stat )
  END IF
  
END PROGRAM SpcCoeff_BIN2NC
