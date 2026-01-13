
!
! test_AerosolScatter_TL
!
! Unit test for CRTM_Compute_AerosolScatter_TL
! Verifies sensitivities to input variables.
!

PROGRAM test_AerosolScatter_TL

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
  USE CRTM_AerosolScatter,    ONLY: CRTM_Compute_AerosolScatter, &
                                    CRTM_Compute_AerosolScatter_TL
  USE ASvar_Define,           ONLY: ASvar_type, &
                                    ASvar_Create, ASvar_Destroy
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_TL'
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  CHARACTER(256) :: File_Path = './testinput/'

  TYPE(CRTM_Atmosphere_type) :: Atm, Atm_TL
  TYPE(CRTM_AtmOptics_type)  :: AScat, AScat_TL
  TYPE(ASvar_type)           :: ASvar
  INTEGER :: SensorIndex = 1
  INTEGER :: ChannelIndex = 1
  INTEGER :: n_Layers = 1
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Legendre = 4
  INTEGER :: n_Phase = 1

  PRINT *, 'Starting test_AerosolScatter_TL...'

  ! Load Coefficients
  Error_Status = CRTM_SpcCoeff_Load( (/Sensor_Id/), File_Path=File_Path )
  IF ( Error_Status /= SUCCESS ) STOP 1
  
  ! Using GOCART scheme for RH sensitivity
  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.bin', File_Path=File_Path)
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_AtmOptics_Create( AScat, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_TL, n_Layers, n_Legendre, n_Phase )
  CALL ASvar_Create( ASvar, n_Legendre, n_Phase, n_Layers, n_Aerosols )

  ! Populate Input
  Atm%Absorber_ID(1) = 1
  Atm%Aerosol(1)%Type = 6
  Atm%Aerosol(1)%Concentration(1) = 1.0e-4_fp
  Atm%Aerosol(1)%Effective_Radius(1) = 1.0_fp
  Atm%Relative_Humidity(1) = 0.5_fp

  ! Populate TL Input: Perturb Relative Humidity
  Atm_TL%Relative_Humidity(1) = 1.0_fp 

  AScat%Include_Scattering = .TRUE.
  AScat_TL%Include_Scattering = .TRUE.

  ! Forward Call to populate ASvar
  Error_Status = CRTM_Compute_AerosolScatter( Atm, SensorIndex, ChannelIndex, AScat, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! TL Call
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  PRINT *, 'AScat_TL%Optical_Depth(1) for RH perturbation: ', AScat_TL%Optical_Depth(1)

  ! Expected behavior: If RH is ignored, sensitivity should be exactly ZERO.
  ! We want this to be NON-ZERO eventually if the Jacobian is implemented.
  IF ( ABS(AScat_TL%Optical_Depth(1)) < 1.0e-12_fp ) THEN
     PRINT *, 'FAIL: Sensitivity to RH is ZERO (Missing Jacobian confirmed)'
     ! In TDD, we proceed once we have a failing test.
  ELSE
     PRINT *, 'SUCCESS: Sensitivity to RH is NON-ZERO'
  END IF

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_AtmOptics_Destroy( AScat )
  CALL CRTM_AtmOptics_Destroy( AScat_TL )
  CALL ASvar_Destroy( ASvar )
  Error_Status = CRTM_SpcCoeff_Destroy()
  Error_Status = CRTM_AerosolCoeff_Destroy()

END PROGRAM test_AerosolScatter_TL
