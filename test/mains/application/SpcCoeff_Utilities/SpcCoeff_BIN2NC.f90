!
! SpcCoeff_BIN2NC
!
! Program to convert binary format SpcCoeff files to the netCDF 
! format.
!
!
! Record of Revisions:
! ====================
!
! Date:          Author:                  Description:
! =====          =======                  ============
! 2021-08-31     Ben Johnson (JCSDA)      Initial Implementation, based on ODPS_BIN2NC
!

PROGRAM SpcCoeff_BIN2NC

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module usage

  USE SpcCoeff_Define    , ONLY: SpcCoeff_type, SpcCoeff_Destroy, SpcCoeff_Equal
  USE SpcCoeff_Binary_IO , ONLY: SpcCoeff_Binary_ReadFile, SpcCoeff_Binary_WriteFile
  USE SpcCoeff_netCDF_IO , ONLY: SpcCoeff_netCDF_ReadFile, SpcCoeff_netCDF_WriteFile

  USE CRTM_Module


  ! Disable implicit typing
  IMPLICIT NONE
  
  ! ----------
  ! Parameters
  ! ----------
  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'SpcCoeff_BIN2NC'

  ! ---------
  ! Variables
  ! ---------
  INTEGER :: Error_Status
  CHARACTER(256)      :: Version
  CHARACTER(256)      :: NC_Filename
  CHARACTER(256)      :: BIN_Filename
  TYPE(SpcCoeff_type) :: SpcCoeff
  TYPE(SpcCoeff_type) :: SpcCoeff_Test
  

  ! Output prgram header
  ! Program header                                                                                                                                                                                    
  ! --------------
  call CRTM_Version(Version)
  CALL Program_Message( PROGRAM_NAME, 'Program to convert binary format SpcCoeff files to netCDF format.', 'CRTM Version:'//trim(Version))

  ! Get the input and output filenames
  WRITE( *, FMT     = '( /5x, "Enter the INPUT Binary SpcCoeff file: " )', ADVANCE = 'NO' )
  READ( *, '( a )' ) BIN_Filename
  
  BIN_Filename = ADJUSTL( BIN_FileNAME )
  IF ( .NOT. File_Exists( TRIM( BIN_Filename ) ) ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'File '//TRIM( BIN_Filename )//' not found.', &
          FAILURE )
     STOP
  END IF
  
  WRITE( *, FMT     = '( /5x, "Enter the OUTPUT netCDF SpcCoeff file: " )', &
       ADVANCE = 'NO' )
  READ( *, '( a )' ) NC_Filename
  NC_Filename = ADJUSTL( NC_Filename )
  
  ! Check that the BIN file isn't accidentally overwritten
  IF ( TRIM( BIN_Filename ) == TRIM( NC_Filename ) ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'Output filename is the same as the input filename!', &
          FAILURE )
     STOP
  END IF
  
  ! Read the input binary file
  WRITE( *, '( /5x, "Reading Binary SpcCoeff data ..." )' )
  Error_Status = SpcCoeff_Binary_ReadFile( BIN_Filename, SpcCoeff )
  IF ( Error_Status /= SUCCESS ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'Error reading Binary SpcCoeff file '//&
          TRIM( BIN_Filename ), &
          Error_Status )
     STOP
  END IF
  
  ! Write the netCDF file
  WRITE( *, '( /5x, "Writeing netCDF SpcCoeff data ..." )' )
  Error_Status = SpcCoeff_netCDF_WriteFile( NC_Filename, SpcCoeff )
  IF ( Error_Status /= SUCCESS ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'Error writing netCDF SpcCoeff file '//&
          TRIM( NC_Filename ), &
          Error_Status )
     STOP
  END IF
  
  
  ! Test read the netCDF data file
  WRITE( *, '( /5x, "Test reading the netCDF SpcCoeff data file ..." )' )
  Error_Status = SpcCoeff_netCDF_ReadFile( NC_Filename, SpcCoeff_Test )
  IF ( Error_Status /= SUCCESS ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'Error reading netCDF SpcCoeff file '//&
          TRIM( NC_Filename ), &
          Error_Status )
     STOP
  END IF
  
  ! Compare the two structures
  WRITE( *, '( /5x, "Comparing the Binary and netCDF SpcCoeff structures ..." )' )
  Error_Status = SpcCoeff_Equal( SpcCoeff_Test, SpcCoeff)
  IF ( Error_Status /= SUCCESS ) THEN
     CALL Display_Message( PROGRAM_NAME, &
          'Differences found in Binary and netCDF'//&
          'file SpcCoeff structure comparison.', &
          Error_Status )
  ELSE
     CALL Display_Message( PROGRAM_NAME, &
          'Binary and netCDF file SpcCoeff structures are equal.', &
          INFORMATION )
  END IF
  
  ! Destroy the structures
  CALL SpcCoeff_Destroy( SpcCoeff )
  CALL SpcCoeff_Destroy( SpcCoeff_Test )
  
END PROGRAM SpcCoeff_BIN2NC
