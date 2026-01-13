
!
! test_AerosolScatter
!
! Unit test for CRTM_Compute_AerosolScatter
!

PROGRAM test_AerosolScatter

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE
  USE CRTM_SpcCoeff,          ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_AerosolCoeff,      ONLY: AeroC, CRTM_AerosolCoeff_Load, CRTM_AerosolCoeff_Destroy
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, &
                                    CRTM_Atmosphere_Create, &
                                    CRTM_Atmosphere_Destroy
  USE CRTM_AtmOptics_Define,  ONLY: CRTM_AtmOptics_type, &
                                    CRTM_AtmOptics_Create, &
                                    CRTM_AtmOptics_Destroy
  USE CRTM_AerosolScatter,    ONLY: CRTM_Compute_AerosolScatter
  USE ASvar_Define,           ONLY: ASvar_type, &
                                    ASvar_Create, ASvar_Destroy
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter'
  CHARACTER(256) :: Message
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id
  CHARACTER(256) :: File_Path

  TYPE(CRTM_Atmosphere_type) :: Atm
  TYPE(CRTM_AtmOptics_type)  :: AScat
  TYPE(ASvar_type)           :: ASvar
  INTEGER :: SensorIndex = 1
  INTEGER :: ChannelIndex = 1
  INTEGER :: n_Layers = 1
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Legendre = 4
  INTEGER :: n_Phase = 2

  ! Setup
  Sensor_Id = 'atms_npp'
  ! Assuming test is run from build/test/mains/unit/Unit_Test
  ! and data is in build/test_data/... which is symlinked to ./testinput
  File_Path = './testinput/' 

  PRINT *, 'Starting test_AerosolScatter...'

  ! Load Coefficients
  PRINT *, 'Loading SpcCoeff...'
  Error_Status = CRTM_SpcCoeff_Load( (/Sensor_Id/), File_Path=File_Path )
  IF ( Error_Status /= SUCCESS ) THEN
     PRINT *, 'CRTM_SpcCoeff_Load FAILED'
     STOP 1
  END IF
  
  PRINT *, 'Loading AerosolCoeff...'
  Error_Status = CRTM_AerosolCoeff_Load('CRTM', 'AerosolCoeff.bin', File_Path=File_Path)
  IF ( Error_Status /= SUCCESS ) THEN
     PRINT *, 'CRTM_AerosolCoeff_Load FAILED'
     STOP 1
  END IF

  ! ... (rest of code) ...

  ! Clean up
  PRINT *, 'Cleaning up...'
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_AtmOptics_Destroy( AScat )
  CALL ASvar_Destroy( ASvar )
  Error_Status = CRTM_SpcCoeff_Destroy()
  Error_Status = CRTM_AerosolCoeff_Destroy()

  PRINT *, 'Test Complete.'

END PROGRAM test_AerosolScatter
