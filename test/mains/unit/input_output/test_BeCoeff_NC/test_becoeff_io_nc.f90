!-------------------------------------------------------
!
! Description:
!	Simple test program to inspect the CRTM BeCoeff
!	files.
!
!-------------------------------------------------------

PROGRAM test_becoeff_io_nc

  ! ====================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !

  ! Module usage
  USE UnitTest_Define, ONLY: UnitTest_type,   &
                             UnitTest_Init,   &
                             UnitTest_Setup,  &
                             UnitTest_Assert, &
                             UnitTest_Passed
  USE CRTM_BeCoeff
  USE Message_Handler, ONLY: SUCCESS, Display_Message

  ! Disable all implicit typing
  IMPLICIT NONE

  ! Data dictionary:
  CHARACTER(2000)         :: info
  CHARACTER(*), PARAMETER :: BeCoeff_File = 'BeCoeff.nc'
  CHARACTER(*), PARAMETER :: File_Path   = './testinput/'
  LOGICAL,      PARAMETER :: netCDF = .TRUE.
  LOGICAL,      PARAMETER :: Quiet  = .TRUE.
  INTEGER                 :: err_stat
  TYPE(UnitTest_type)     :: ioTest
  LOGICAL                 :: testPassed
  CHARACTER(*), PARAMETER :: Program_Name = 'Test_BeCoeff_IO_NetCDF'

  ! Initialize Unit test:
  CALL UnitTest_Init(ioTest, .TRUE.)
  CALL UnitTest_Setup(ioTest, 'BeCoeff_IO_Test_NetCDF', Program_Name, .TRUE.)

  ! Greeting:
  WRITE(*,*) 'HELLO, THIS IS A TEST CODE TO INSPECT BeCoeff files.'
  WRITE(*,*) 'test_becoeff_io_nc'

  ! Load the BeCoeff look-up table:
  err_stat = 3
  err_stat = CRTM_BeCoeff_Load( &
                BeCoeff_File         , &
                File_Path    = File_Path, &
                netCDF       = netCDF, &
                Quiet        = Quiet )
  CALL UnitTest_Assert(ioTest, (err_stat==SUCCESS) )
  testPassed = UnitTest_Passed(ioTest)

  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( 'CRTM_BeCoeff_Load' ,'Error loading BeCoeff data', err_stat )
    STOP 1
  END IF
  STOP 0

END PROGRAM test_becoeff_io_nc
