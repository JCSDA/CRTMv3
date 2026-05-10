
!
! test_AerosolScatter_K
!
! Unit test for CRTM_K_Matrix focusing on Aerosol sensitivities (Jacobians).
!

PROGRAM test_AerosolScatter_K

  USE Type_Kinds,             ONLY: fp
  USE Message_Handler,        ONLY: SUCCESS, FAILURE, Display_Message
  USE CRTM_Parameters,        ONLY: ZERO, ONE
  USE CRTM_SpcCoeff,          ONLY: SC, CRTM_SpcCoeff_Load, CRTM_SpcCoeff_Destroy
  USE CRTM_AerosolCoeff,      ONLY: AeroC, CRTM_AerosolCoeff_Load, CRTM_AerosolCoeff_Destroy
  USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, &
                                    CRTM_Atmosphere_Create, &
                                    CRTM_Atmosphere_Destroy, &
                                    CRTM_Atmosphere_Zero, &
                                    H2O_ID, O3_ID, CO2_ID, N2O_ID, CH4_ID, CO_ID, &
                                    MASS_MIXING_RATIO_UNITS
  USE CRTM_AtmOptics_Define,  ONLY: CRTM_AtmOptics_type, &
                                    CRTM_AtmOptics_Create, &
                                    CRTM_AtmOptics_Destroy, &
                                    CRTM_AtmOptics_Zero
  USE CRTM_Surface_Define,    ONLY: CRTM_Surface_type, &
                                    CRTM_Surface_Create, &
                                    CRTM_Surface_Destroy, &
                                    CRTM_Surface_Zero
  USE CRTM_ChannelInfo_Define, ONLY: CRTM_ChannelInfo_type
  USE CRTM_Geometry_Define,    ONLY: CRTM_Geometry_type
  USE CRTM_RTSolution_Define,  ONLY: CRTM_RTSolution_type
  USE CRTM_LifeCycle,          ONLY: CRTM_Init, CRTM_Destroy
  USE CRTM_K_Matrix_Module,    ONLY: CRTM_K_Matrix

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: PROGRAM_NAME = 'test_AerosolScatter_K'
  CHARACTER(*), PARAMETER :: COEFFICIENTS_PATH = './testinput/'
  
  INTEGER :: Error_Status
  CHARACTER(256) :: Sensor_Id = 'modis_aqua'
  
  TYPE(CRTM_ChannelInfo_type) :: ChannelInfo(1)
  TYPE(CRTM_Geometry_type)    :: Geometry(1)
  TYPE(CRTM_Atmosphere_type)  :: Atm(1)
  TYPE(CRTM_Atmosphere_type), ALLOCATABLE :: Atm_K(:,:)
  TYPE(CRTM_Surface_type)     :: Sfc(1)
  TYPE(CRTM_Surface_type), ALLOCATABLE :: Sfc_K(:,:)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: RTSolution(:,:), RTSolution_K(:,:)
  
  INTEGER :: n_Channels
  INTEGER :: n_Layers = 10
  INTEGER :: n_Absorbers = 6
  INTEGER :: n_Aerosols = 1
  INTEGER :: n_Clouds = 0
  
  INTEGER :: l, k

  PRINT *, 'Starting test_AerosolScatter_K...'

  ! Initialize CRTM
  Error_Status = CRTM_Init( (/Sensor_Id/), ChannelInfo, File_Path=COEFFICIENTS_PATH )
  IF ( Error_Status /= SUCCESS ) STOP 1
  
  ! Load GOCART-GEOS5
  Error_Status = CRTM_AerosolCoeff_Load('GOCART-GEOS5', 'AerosolCoeff.GOCART-GEOS5.nc', File_Path=COEFFICIENTS_PATH, netCDF=.TRUE.)
  IF ( Error_Status /= SUCCESS ) STOP 1

  n_Channels = ChannelInfo(1)%n_Channels
  ALLOCATE( RTSolution( n_Channels, 1 ), &
            RTSolution_K( n_Channels, 1 ), &
            Atm_K( n_Channels, 1 ), &
            Sfc_K( n_Channels, 1 ) )

  ! Create structures
  CALL CRTM_Atmosphere_Create( Atm, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Atmosphere_Create( Atm_K, n_Layers, n_Absorbers, n_Clouds, n_Aerosols )
  CALL CRTM_Surface_Create( Sfc_K, 0 ) ! Allocate Surface_K too!
  
  ! Populate Input
  Atm(1)%Absorber_ID(1) = H2O_ID
  Atm(1)%Absorber_ID(2) = O3_ID
  Atm(1)%Absorber_ID(3) = CO2_ID
  Atm(1)%Absorber_ID(4) = N2O_ID
  Atm(1)%Absorber_ID(5) = CH4_ID
  Atm(1)%Absorber_ID(6) = CO_ID
  Atm(1)%Aerosol(1)%Type = 6 ! Sea Salt (hygroscopic)
  
  ! Standard Atmosphere (simplified)
  Atm(1)%Pressure(1:n_Layers) = (/ (1000.0_fp - REAL(k-1,fp)*100.0_fp, k=1,n_Layers) /)
  Atm(1)%Temperature(1:n_Layers) = 280.0_fp
  Atm(1)%Absorber_Units = MASS_MIXING_RATIO_UNITS
  Atm(1)%Absorber(1:n_Layers, 1) = 1.0e-3_fp ! H2O
  Atm(1)%Absorber(1:n_Layers, 2) = 1.0e-6_fp ! O3
  Atm(1)%Absorber(1:n_Layers, 3:6) = 1.0e-7_fp
  Atm(1)%Aerosol(1)%Concentration(1:n_Layers) = 1.0e-5_fp
  Atm(1)%Aerosol(1)%Effective_Radius(1:n_Layers) = 1.0_fp
  Atm(1)%Relative_Humidity(1:n_Layers) = 0.5_fp ! 50% as fraction

  ! Geometry
  Geometry(1)%Sensor_Zenith_Angle = 30.0_fp
  
  ! Surface
  CALL CRTM_Surface_Create( Sfc, 0 ) ! Land
  Sfc(1)%Land_Coverage = 1.0_fp
  Sfc(1)%Land_Temperature = 290.0_fp

  ! Initialize K-Matrix inputs: Sensitivity to Brightness Temperature
  CALL CRTM_Atmosphere_Zero( Atm_K )
  CALL CRTM_Surface_Zero( Sfc_K )
  DO l = 1, n_Channels
    RTSolution_K(l,1)%Radiance = ZERO
    RTSolution_K(l,1)%Brightness_Temperature = ONE
  END DO

  ! Call K-Matrix
  Error_Status = CRTM_K_Matrix( Atm, Sfc, RTSolution_K, Geometry, ChannelInfo, Atm_K, Sfc_K, RTSolution )
  IF ( Error_Status /= SUCCESS ) STOP 1

  ! Inspect Aerosol Jacobians for Channel 1 (VIS)
  PRINT *, 'Aerosol Jacobians for Channel 1 (MODIS VIS 0.65um):'
  PRINT *, 'Layer, RH_Jacobian, Concentration_Jacobian'
  DO k = 1, n_Layers
    WRITE(*, '(i3, 2(1x,es12.5))') k, Atm_K(1,1)%Relative_Humidity(k), Atm_K(1,1)%Aerosol(1)%Concentration(k)
  END DO

  ! Check if RH Jacobian is non-zero
  IF ( ANY(ABS(Atm_K(1,1)%Relative_Humidity) > 1.0e-12_fp) ) THEN
    PRINT *, 'SUCCESS: Aerosol RH Jacobians are NON-ZERO.'
  ELSE
    PRINT *, 'FAIL: Aerosol RH Jacobians are ZERO.'
    STOP 1
  END IF

  ! Clean up
  CALL CRTM_Atmosphere_Destroy( Atm )
  CALL CRTM_Atmosphere_Destroy( Atm_K )
  CALL CRTM_Surface_Destroy( Sfc )
  CALL CRTM_Surface_Destroy( Sfc_K )
  DEALLOCATE( RTSolution, RTSolution_K )
  Error_Status = CRTM_Destroy( ChannelInfo )
  Error_Status = CRTM_AerosolCoeff_Destroy()

END PROGRAM test_AerosolScatter_K
