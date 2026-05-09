
!
! test_AerosolScatter_AD
!
! Unit test for CRTM_Compute_AerosolScatter_AD
! Verifies Adjoint model consistency with Tangent Linear model.
!

PROGRAM test_AerosolScatter_AD

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE
  USE CRTM_SpcCoeff,          ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_AerosolCoeff,      ONLY: AeroC, CRTM_AerosolCoeff_Load, CRTM_AerosolCoeff_Destroy
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, &
                                    CRTM_Atmosphere_Create, &
                                    CRTM_Atmosphere_Destroy, &
                                    CRTM_Atmosphere_Zero
  USE CRTM_AtmOptics_Define,  ONLY: CRTM_AtmOptics_type, &
                                    CRTM_AtmOptics_Create, &
                                    CRTM_AtmOptics_Destroy, &
                                    CRTM_AtmOptics_Zero
  USE CRTM_AerosolScatter,    ONLY: CRTM_Compute_AerosolScatter, &
                                    CRTM_Compute_AerosolScatter_TL, &
                                    CRTM_Compute_AerosolScatter_AD
  USE ASvar_Define,           ONLY: ASvar_type, &
                                    ASvar_Create, ASvar_Destroy
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_AD'
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  CHARACTER(256) :: File_Path = './testinput/'

  TYPE(CRTM_Atmosphere_type) :: Atm, Atm_TL, Atm_AD
  TYPE(CRTM_AtmOptics_type)  :: AScat, AScat_TL, AScat_AD
  TYPE(ASvar_type)           :: ASvar
  INTEGER :: SensorIndex = 1
  INTEGER :: ChannelIndex = 1
  INTEGER :: n_Layers = 1
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Legendre = 4
  INTEGER :: n_Phase = 1
  
  REAL(fp) :: sum_tl, sum_ad
  INTEGER :: k, n, l, m

  PRINT *, 'Starting test_AerosolScatter_AD...'

  ! Load Coefficients
  Error_Status = CRTM_SpcCoeff_Load( (/Sensor_Id/), File_Path=File_Path )
  IF ( Error_Status /= SUCCESS ) STOP 1
  
  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.nc4', File_Path=File_Path, netCDF=.TRUE.)
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_TL, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_AD, n_Layers, 1, 0, n_Aerosols )
  CALL CRTM_AtmOptics_Create( AScat, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_TL, n_Layers, n_Legendre, n_Phase )
  CALL CRTM_AtmOptics_Create( AScat_AD, n_Layers, n_Legendre, n_Phase )
  CALL ASvar_Create( ASvar, n_Legendre, n_Phase, n_Layers, n_Aerosols )

  ! Populate Input
  Atm%Absorber_ID(1) = 1
  Atm%Aerosol(1)%Type = 6 ! Sea Salt
  Atm%Aerosol(1)%Concentration(1) = 1.0e-4_fp
  Atm%Aerosol(1)%Effective_Radius(1) = 1.0_fp
  Atm%Relative_Humidity(1) = 0.5_fp

  ! Populate TL Input: Perturb Relative Humidity and Concentration
  CALL CRTM_Atmosphere_Zero( Atm_TL )
  Atm_TL%Relative_Humidity(1) = 0.1_fp 
  Atm_TL%Aerosol(1)%Concentration(1) = 1.0e-5_fp

  AScat%Include_Scattering = .TRUE.
  AScat_TL%Include_Scattering = .TRUE.
  AScat_AD%Include_Scattering = .TRUE.

  ! Forward Call
  Error_Status = CRTM_Compute_AerosolScatter( Atm, SensorIndex, ChannelIndex, AScat, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! TL Call
  CALL CRTM_AtmOptics_Zero( AScat_TL )
  Error_Status = CRTM_Compute_AerosolScatter_TL( Atm, AScat, Atm_TL, SensorIndex, ChannelIndex, AScat_TL, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Adjoint Input: Set AScat_AD = AScat_TL
  CALL CRTM_AtmOptics_Zero( AScat_AD )
  AScat_AD%Optical_Depth = AScat_TL%Optical_Depth
  AScat_AD%Single_Scatter_Albedo = AScat_TL%Single_Scatter_Albedo
  AScat_AD%Asymmetry_Factor = AScat_TL%Asymmetry_Factor
  AScat_AD%Phase_Coefficient = AScat_TL%Phase_Coefficient

  ! AD Call
  CALL CRTM_Atmosphere_Zero( Atm_AD )
  Error_Status = CRTM_Compute_AerosolScatter_AD( Atm, AScat, AScat_AD, SensorIndex, ChannelIndex, Atm_AD, ASvar )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Adjoint Test: sum(AScat_TL * AScat_AD_input) == sum(Atm_TL * Atm_AD_output)
  ! sum_tl = sum(AScat_TL**2)
  sum_tl = SUM(AScat_TL%Optical_Depth**2) + &
           SUM(AScat_TL%Single_Scatter_Albedo**2) + &
           SUM(AScat_TL%Asymmetry_Factor**2) + &
           SUM(AScat_TL%Phase_Coefficient**2)

  ! sum_ad = sum(Atm_TL * Atm_AD)
  sum_ad = SUM(Atm_TL%Relative_Humidity * Atm_AD%Relative_Humidity) + &
           SUM(Atm_TL%Aerosol(1)%Concentration * Atm_AD%Aerosol(1)%Concentration) + &
           SUM(Atm_TL%Aerosol(1)%Effective_Radius * Atm_AD%Aerosol(1)%Effective_Radius) + &
           SUM(Atm_TL%Aerosol(1)%Effective_Variance * Atm_AD%Aerosol(1)%Effective_Variance)

  PRINT *, 'LHS (TL sum): ', sum_tl
  PRINT *, 'RHS (AD sum): ', sum_ad
  PRINT *, 'Diff:         ', ABS(sum_tl - sum_ad)
  
  IF ( sum_tl > ZERO ) THEN
    PRINT *, 'Rel. Diff:    ', ABS(sum_tl - sum_ad) / sum_tl
    IF ( ABS(sum_tl - sum_ad) / sum_tl < 1.0e-12_fp ) THEN
       PRINT *, 'SUCCESS: Adjoint test passed.'
    ELSE
       PRINT *, 'FAIL: Adjoint test failed.'
       STOP 1
    END IF
  ELSE
    PRINT *, 'FAIL: LHS is ZERO, test is invalid.'
    STOP 1
  END IF

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_TL )
  CALL CRTM_Atmosphere_Destroy( Atm_AD )
  CALL CRTM_AtmOptics_Destroy( AScat )
  CALL CRTM_AtmOptics_Destroy( AScat_TL )
  CALL CRTM_AtmOptics_Destroy( AScat_AD )
  CALL ASvar_Destroy( ASvar )
  Error_Status = CRTM_SpcCoeff_Destroy()
  Error_Status = CRTM_AerosolCoeff_Destroy()

END PROGRAM test_AerosolScatter_AD
