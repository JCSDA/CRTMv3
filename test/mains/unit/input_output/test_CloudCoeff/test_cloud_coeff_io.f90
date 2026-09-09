!-------------------------------------------------------
!
! Description:
!	Simple test program to inspect the CRTM CloudCoeff
!	files.
!
!	Date: 2018-08-14	Author: P. Stegmann
!
! MODIFICATION HISTORY:
! =====================
!
! Author:           Date:          Description:
! =======           =====          ============
! Patrick Stegmann  2021-02-05     Refactored as a CRTM
!                                  unit test.
! Cheng Dang        2021-07-28     Modified for Cloud
!                                  Coeff look-up table
!-------------------------------------------------------

PROGRAM test_cloud_coeff_io

  ! ====================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !

  ! Module usage
  USE UnitTest_Define, ONLY: UnitTest_type,   &
                             UnitTest_Init,   &
                             UnitTest_Setup,  &
                             UnitTest_Assert, &
                             UnitTest_Passed
  !USE CloudCoeff_Define, ONLY: CloudCoeff_type
  USE CRTM_CloudCoeff
  USE Message_Handler, ONLY: SUCCESS, Display_Message

  ! Disable all implicit typing
  IMPLICIT NONE

  ! Data dictionary:
  CHARACTER(2000)         :: info
  CHARACTER(*), PARAMETER :: Cloud_Model = 'CRTM'
  CHARACTER(*), PARAMETER :: CloudCoeff_File = 'CloudCoeff.nc'
  CHARACTER(*), PARAMETER :: File_Path = './testinput/'
  LOGICAL,      PARAMETER :: Quiet  = .TRUE.
  LOGICAL,      PARAMETER :: netCDF = .TRUE.
  INTEGER                 :: err_stat
  TYPE(UnitTest_type)     :: ioTest
  LOGICAL                 :: testPassed
  CHARACTER(*), PARAMETER :: Program_Name = 'Test_Cloud_Coeff_IO_Binary'

  ! Initialize Unit test:
  CALL ioTest%Init(.TRUE.)
  CALL ioTest%Setup('Test_Cloud_Coeff_IO_Binary', Program_Name, .TRUE.)

  ! Greeting:
  WRITE(*,*) 'HELLO, THIS IS A TEST CODE TO INSPECT CloudCoefff files.'
  WRITE(*,*) 'test_cloud_coeff_io', 'The following Cloud scheme is investigated: ', Cloud_Model
  ! Load the Cloud coefficient look-up table:
  err_stat = 3
  err_stat = CRTM_CloudCoeff_Load( &
                Cloud_Model              , &
                CloudCoeff_File          , &
                File_Path                , &
                netCDF          = netCDF , &
                Quiet           = Quiet    )
  CALL ioTest%Assert((err_stat==SUCCESS) )
  testPassed = ioTest%Passed()

  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( 'CRTM_Load_Cloud_Coeff' ,'Error loading CloudCoeff data', err_stat )
    STOP 1
  END IF
  STOP 0

END PROGRAM test_cloud_coeff_io
