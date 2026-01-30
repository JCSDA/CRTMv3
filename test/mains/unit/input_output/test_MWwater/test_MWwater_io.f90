!-------------------------------------------------------
!
! Description:
!       Test replacement module for MW water emissivity FASTEM
!
!       Date: 2024-10-07        Author: Cheng Dang

! MODIFICATION HISTORY:
! =====================
!
!-------------------------------------------------------

PROGRAM test_MWwater_io

  ! ====================================================
  ! **** ENVIRONMENT SETUP FOR RTM USAGE ****
  !

  ! Module usage
  USE UnitTest_Define, ONLY: UnitTest_type,   &
                             UnitTest_Init,   &
                             UnitTest_Setup,  &
                             UnitTest_Assert, &
                             UnitTest_Passed
  ! ...Infrared surface emissivities
  USE CRTM_MWwaterCoeff  , ONLY: CRTM_MWwaterCoeff_Load, &
                                 CRTM_MWwaterCoeff_Load_FASTEM, &
                                 MWwaterC
  USE MWwaterCoeff_Define, ONLY: MWwaterCoeff_type, &
                                 MWwaterCoeff_Create, &
                                 MWwaterCoeff_Destroy, &
                                 MWwaterCoeff_Equal
  USE Message_Handler    , ONLY: SUCCESS, Display_Message

  ! Disable all implicit typing
  IMPLICIT NONE

  ! Data dictionary:
  LOGICAL,      PARAMETER :: Quiet = .TRUE.
  CHARACTER(*), PARAMETER :: Program_Name = 'test_MWwater_io'
  CHARACTER(*), PARAMETER :: File_Path = './testinput/'
  CHARACTER(*), PARAMETER :: File_FASTEM6   = 'FASTEM6.MWwater.EmisCoeff.bin'
  CHARACTER(*), PARAMETER :: File_FASTEM5   = 'FASTEM5.MWwater.EmisCoeff.bin'
  CHARACTER(*), PARAMETER :: File_FASTEM4   = 'FASTEM4.MWwater.EmisCoeff.bin'
  CHARACTER(*), PARAMETER :: Scheme_FASTEM6 = 'FASTEM6'
  CHARACTER(*), PARAMETER :: Scheme_FASTEM5 = 'FASTEM5'
  CHARACTER(*), PARAMETER :: Scheme_FASTEM4 = 'FASTEM4'

  TYPE(MWwaterCoeff_type) :: MWwaterC_LUT, MWwaterC_NEW
  INTEGER                 :: err_stat
  TYPE(UnitTest_type)     :: ioTest
  LOGICAL                 :: testPassed, is_equal

  ! Initialize Unit test:
  CALL UnitTest_Init(ioTest, .TRUE.)
  CALL UnitTest_Setup(ioTest, 'MWwater_Coeff_IO_Test', Program_Name, .TRUE.)

  ! Load the default emissivity coefficient look-up table:

  ! FASTEM 6
  WRITE(*,*) 'CRTM_MWwaterCoeff_Load ...LOADING: ', File_FASTEM6
  CALL MWwaterCoeff_Create (MWwaterC)
  CALL MWwaterCoeff_Create (MWwaterC_LUT)
  err_stat = CRTM_MWwaterCoeff_Load(&
                  TRIM(TRIM(File_Path)//File_FASTEM6), &
                  Quiet = .TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( 'CRTM_MWwaterCoeff_Load' ,'Error loading MWwaterCoeff data', err_stat )
    STOP 1
  END IF
  MWwaterC_LUT = MWwaterC
  CALL MWwaterCoeff_Destroy (MWwaterC)

  WRITE(*,*) 'CRTM_MWwaterCoeff_Load_FASTEM ...LOADING: ', Scheme_FASTEM6
  CALL MWwaterCoeff_Create (MWwaterC)
  CALL MWwaterCoeff_Create (MWwaterC_NEW)
  
  err_stat = CRTM_MWwaterCoeff_Load_FASTEM( &
                     Scheme_FASTEM6, &
                     Quiet = .TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
   CALL Display_Message( 'CRTM_MWwaterCoeff_Load_FASTEM' ,'Error loading MWwaterCoeff data', err_stat )
   STOP 1
  END IF
  MWwaterC_NEW = MWwaterC
  CALL MWwaterCoeff_Destroy (MWwaterC)

  is_equal = MWwaterCoeff_Equal(MWwaterC_LUT, MWwaterC_NEW)
  IF ( .NOT. is_equal ) THEN
    CALL Display_Message( 'MWwaterCoeff_Equal' ,'MWwaterCoeff are different', err_stat )
    STOP 1
  END IF
  CALL MWwaterCoeff_Destroy (MWwaterC_LUT)
  CALL MWwaterCoeff_Destroy (MWwaterC_NEW)

  ! FASTEM 5 - not supported by CRTM_MWwaterCoeff_Load_FASTEM

  ! FASTEM 4
  WRITE(*,*) 'CRTM_MWwaterCoeff_Load ...LOADING: ', File_FASTEM4
  CALL MWwaterCoeff_Create (MWwaterC)
  err_stat = CRTM_MWwaterCoeff_Load(&
                  TRIM(TRIM(File_Path)//File_FASTEM4), &
                  Quiet = .TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( 'CRTM_MWwaterCoeff_Load' ,'Error loading MWwaterCoeff data', err_stat )
    STOP 1
  END IF
  MWwaterC_LUT = MWwaterC
  CALL MWwaterCoeff_Destroy (MWwaterC)

  WRITE(*,*) 'CRTM_MWwaterCoeff_Load_FASTEM ...LOADING: ', Scheme_FASTEM4
  CALL MWwaterCoeff_Create (MWwaterC)
  err_stat = CRTM_MWwaterCoeff_Load_FASTEM( &
                     Scheme_FASTEM4, &
                     Quiet = .TRUE.)
  IF ( err_stat /= SUCCESS ) THEN
    CALL Display_Message( 'CRTM_MWwaterCoeff_Load_FASTEM' ,'Error loading MWwaterCoeff data', err_stat )
    STOP 1
  END IF
  MWwaterC_NEW = MWwaterC
  CALL MWwaterCoeff_Destroy (MWwaterC)

  is_equal = MWwaterCoeff_Equal(MWwaterC_LUT, MWwaterC_NEW)
  IF ( .NOT. is_equal ) THEN
    CALL Display_Message( 'MWwaterCoeff_Equal' ,'MWwaterCoeff are different', err_stat )
    STOP 1
  END IF
  CALL MWwaterCoeff_Destroy (MWwaterC_LUT)
  CALL MWwaterCoeff_Destroy (MWwaterC_NEW)

  STOP 0

END PROGRAM test_MWwater_io
