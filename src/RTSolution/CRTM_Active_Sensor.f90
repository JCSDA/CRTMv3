!
! CRTM_Active_Sensor
!
! Module containing the active sensor routines.
!
!
! CREATION HISTORY:
!       Written by:     Isaac Moradi - NASA GSFC March 10, 2022
!                       Isaac.Moradi@NASA.GOV
!

MODULE CRTM_Active_Sensor

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use statements
  USE Type_Kinds,      ONLY: FP
  USE CRTM_Parameters, ONLY: ZERO, ONE, TWO, TEN, ONE_THOUSAND, MISSING_REFL, EPSILON_FP
  USE CRTM_SpcCoeff,            ONLY: SC, &
                                      SpcCoeff_IsMicrowaveSensor , &
                                      SpcCoeff_IsInfraredSensor  , &
                                      SpcCoeff_IsVisibleSensor   , &
                                      SpcCoeff_IsUltravioletSensor
  USE CRTM_Atmosphere_Define,   ONLY: CRTM_Atmosphere_type, &
                                      H2O_ID
  USE CRTM_AtmOptics_Define,     ONLY: CRTM_AtmOptics_type
  USE CRTM_RTSolution_Define,    ONLY: CRTM_RTSolution_type
  USE Spectral_Units_Conversion, ONLY: GHz_to_inverse_cm
  USE Fundamental_Constants,     ONLY: PI
  USE ODPS_CoordinateMapping,     ONLY: Geopotential_Height
  USE CRTM_GeometryInfo_Define, ONLY: CRTM_GeometryInfo_type
  USE Message_Handler         , ONLY: Display_Message
  ! Disable all implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC  :: CRTM_Compute_Reflectivity
  PUBLIC  :: CRTM_Compute_Reflectivity_TL
  PUBLIC  :: CRTM_Compute_Reflectivity_AD
  PUBLIC  :: Calculate_Height
  PUBLIC  :: Calculate_Cloud_Water_Density

  ! ----------
  ! Parameters
  ! ----------
  REAL(fp), PARAMETER :: POINT_01 = 0.01_fp
  REAL(fp), PARAMETER :: Kw_2_fixed = 0.93_fp
  ! 1.0d818 converts from m^3 to mm^6/m^3 (standard radar reflectivity units)
  REAL(fp), PARAMETER :: M6_MM6 = 1.0d18
  REAL(fp), PARAMETER :: REFLECTIVITY_THRESHOLD = TINY(REAL(fp))
  REAL(fp), PARAMETER :: MM6_TO_CM6 = 1.0e-6_fp
CONTAINS

!--------------------------------------------------------------------------------
!
! NAME:
!       Calculate_Height
!
! PURPOSE:
!       Subroutine to calculate the height usig atmospheric input profiles
!
! CALLING SEQUENCE:
!       Height = Calculate_Height(Atm)
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
! OUTPUT ARGUMENTS:
!       Height:         Array containing height
!                       UNITS:      km
!                       TYPE:       REAL
!                       DIMENSION:  Array
!
!--------------------------------------------------------------------------------

Function Calculate_Height(Atm) RESULT (Height)

   TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
   REAL(fp) :: Height(0:Atm%n_Layers)

   INTEGER  :: j, H2O_idx

   !** locate the WV absorber ID
   H2O_idx = 1
   DO j = 1,Atm%n_Absorbers
        IF (Atm%Absorber_ID(j) == H2O_ID) H2O_idx = j
   END DO
   ! Calculate the geometric heights of the pressure levels in km
   CALL Geopotential_Height(Atm%Level_Pressure      , & ! Input
                            Atm%Temperature         , & ! Input
                            Atm%Absorber(:, H2O_idx), & ! Input
                            ZERO                    , & ! Input - surface height
                            Height                ) ! Output in km
   
   Height = Height * ONE_THOUSAND ! convert to meters

END FUNCTION Calculate_Height

!--------------------------------------------------------------------------------
!
! NAME:
!       Calculate_Cloud_Water_Density
!
! PURPOSE:
!       Subroutine to calculate cloud water density (mass/volume) and also height if not set yet
!
! CALLING SEQUENCE:
!       CALL Calculate_Cloud_Water_Density(Atm, GeometryInfo)
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN/OUT)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!--------------------------------------------------------------------------------

SUBROUTINE Calculate_Cloud_Water_Density(Atm, &
                                         GeometryInfo)
   TYPE(CRTM_Atmosphere_type), INTENT(IN OUT)     :: Atm
   TYPE(CRTM_GeometryInfo_type), OPTIONAL, INTENT(IN)     :: GeometryInfo
   INTEGER :: n_Layers
   REAL(fp) :: Height(0:Atm%n_Layers), dZ_m(Atm%n_Layers)
   Integer :: n

   n_Layers = Atm%n_Layers

   ! Calculate heights if hasn't been set already
   IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
       Atm%Height = Calculate_Height(Atm)
   END IF

   dZ_m = (Atm%Height(0:n_Layers-1) - Atm%Height(1:n_Layers)) 

   DO n = 1, Atm%n_Clouds
      Atm%Cloud(n)%Water_Density = Atm%Cloud(n)%Water_Content / dZ_m
   END DO

END SUBROUTINE Calculate_Cloud_Water_Density

!--------------------------------------------------------------------------------
!
! NAME:
!       Water_Permittivity_Turner_2016
!
! PURPOSE:
!       Function to calculate water permittivity using Turnet et al. 2016 Method
!
! CALLING SEQUENCE:
!       CALL Water_Permittivity_Turner_2016(freq_Hz, & ! Input
!                                           temp_K)  ! Input
!                                           RESULT(perm) ! Output
!
! INPUT ARGUMENTS:
!
!       freq_Hz:        Scalar representing frequency in Hz.
!                       UNITS:      Hz
!                       TYPE:       REAL(fp)
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       temp_K:         Array representing temperature in Kelvin.
!                       UNITS:      Kelvin
!                       TYPE:       REAL(fp)
!                       DIMENSION:  1-D Array
!                       ATTRIBUTES: INTENT(IN)
!
! OOUTUT ARGUMENTS:
!
!       perm:        Array representing water permittivity in a complex form of (real, imaginary).
!                       UNITS:      N/A
!                       TYPE:       REAL(fp)
!                       DIMENSION:  1-D Complex Array
!
!--------------------------------------------------------------------------------

FUNCTION Water_Permittivity_Turner_2016(freq_Hz, temp_K) RESULT(perm)

  REAL(fp), INTENT(IN) :: freq_Hz, temp_K(:)
  REAL(fp), DIMENSION(SIZE(temp_K)) :: perm_re, perm_im
  COMPLEX :: perm(SIZE(temp_K))

  REAL(fp) :: a1, a2, b1, b2, c1, c2, d1, d2
  REAL(fp) :: e0, e1, e2, e3, tc
  REAL(fp), DIMENSION(SIZE(temp_K)) :: delta1, delta2, tau1, tau2
  REAL(fp), DIMENSION(SIZE(temp_K)) :: term11, term12, denom1, denom2
  REAL(fp), DIMENSION(SIZE(temp_K)) :: term21, term22, temp_c, eps_s

  temp_c = temp_k - 273.15_fp

  IF (freq_Hz .GT. 500.d9) THEN
     STOP
  END IF

  a1 = 8.111E+1
  a2 = 2.025E+0

  b1 = 4.434E-3
  b2 = 1.073E-2

  c1 = 1.302E-13
  c2 = 1.012E-14

  d1 = 6.627E+2
  d2 = 6.089E+2

  tc = 133.1383_fp

  e0 = 87.9144_fp
  e1 = 0.404399_fp
  e2 = 9.58726d-4
  e3 = -1.32802d-6

  ! static dielectric constant
  eps_s = e0 + e1 * temp_c + e2 * (temp_c**TWO) + e3 * (temp_c**3.0_fp)

  delta1 = a1 * EXP(-b1 * temp_c)
  delta2 = a2 * EXP(-b2 * temp_c)

  tau1   = c1 * EXP(d1 / (temp_c + tc))
  tau2   = c2 * EXP(d2 / (temp_c + tc))

  denom1 = ONE + (TWO * PI * freq_Hz * tau1)**TWO
  denom2 = ONE + (TWO * PI * freq_Hz * tau2)**TWO

  term11 = ((tau1**TWO) * delta1) / denom1
  term12 = ((tau2**TWO) * delta2) / denom2

  term21 = (tau1 * delta1) / denom1
  term22 = (tau2 * delta2) / denom2

  perm_re = eps_s - ((TWO * PI * freq_Hz)**TWO) * (term11 + term12)
  perm_im = (TWO * PI * freq_Hz) * (term21 + term22)

  perm = cmplx(perm_re, perm_im, FP)

END FUNCTION Water_Permittivity_Turner_2016

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_Reflectivity
!
! PURPOSE:
!       Subroutine to calculate the reflectivity for active sensors
!
! CALLING SEQUENCE:
!       CALL CRTM_Compute_Reflectivity(Atm         , &  ! Input
!                                      AtmOptics   , &  ! Input
!                                      GeometryInfo , &  ! Input
!                                      SensorIndex , &  ! Input
!                                      ChannelIndex, &  ! Input
!                                      RTSolution )     ! Input/Output
!
! INPUT ARGUMENTS:
!
!       Atm:            Structure containing the atmospheric state data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_Atmosphere_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       UNITS:      N/A
!                       TYPE:       CRTM_AtmOptics_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       UNITS:      N/A
!                       TYPE:       CRTM_GeometryInfo_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:     Sensor index id. This is a unique index associated
!                        with a (supported) sensor used to access the
!                        shared coefficient data for a particular sensor.
!                        See the ChannelIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:    Channel index id. This is a unique index associated
!                        with a (supported) sensor channel used to access the
!                        shared coefficient data for a particular sensor's
!                        channel.
!                        See the SensorIndex argument.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution:     Structure containing the soluition to the RT equation
!                       for the given inputs.
!                       UNITS:      N/A
!                       TYPE:       CRTM_RTSolution_type
!                       DIMENSION:  Scalar
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Compute_Reflectivity(Atm         , &  ! Input
                                       AtmOptics   , &  ! Input
                                       GeometryInfo , &  ! Input
                                       SensorIndex , &  ! Input
                                       ChannelIndex, &  ! Input
                                       RTSolution )     ! Input/Output
    ! Arguments
    TYPE(CRTM_Atmosphere_type), INTENT(IN)     :: Atm
    TYPE(CRTM_GeometryInfo_type), INTENT(IN)     :: GeometryInfo
    INTEGER                   , INTENT(IN)     :: SensorIndex
    INTEGER                   , INTENT(IN)     :: ChannelIndex
    TYPE(CRTM_AtmOptics_type) , INTENT(IN OUT)     :: AtmOptics
    TYPE(CRTM_RTSolution_type), INTENT(IN OUT) :: RTSolution

    REAL(fp) :: Frequency, Wavenumber, Wavelength_m
    REAL(fp) :: Reflectivity(AtmOptics%n_Layers)
    REAL(fp) :: Reflectivity_Attenuated(AtmOptics%n_Layers)
    REAL(fp) :: Transmittance(AtmOptics%n_Layers)
    REAL(fp) :: Optical_Depth(AtmOptics%n_Layers)
    REAL(fp) :: P1(AtmOptics%n_Layers)
    REAL(fp) :: Height(0:AtmOptics%n_Layers)
    REAL(fp) :: dZ_m(AtmOptics%n_Layers)
    REAL(fp) :: dS_m(AtmOptics%n_Layers)
    REAL(fp) :: Temp_K(AtmOptics%n_Layers)
    REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Kw_2, perm_re, perm_im
    COMPLEX :: perm(AtmOptics%n_Layers), kw(AtmOptics%n_Layers)
    INTEGER :: k, n_Layers
    REAL(fp) :: temp_sum
    REAL(fp) :: ext_coef(AtmOptics%n_Layers) ! extinction coefficient  
    character(len=16) :: ms_env_val

    n_Layers = Atm%n_Layers

    ! Calculate heights if hasn't been set already
    IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
        Height = Calculate_Height(Atm)
    ELSE
        Height = Atm%Height
    ENDIF
    dZ_m = (Height(0:n_Layers-1) - Height(1:n_Layers))
    dS_m = dZ_m / GeometryInfo%Cosine_Sensor_Zenith

    IF ( SpcCoeff_IsMicrowaveSensor(SC(SensorIndex)) ) THEN
       Frequency = SC(SensorIndex)%Frequency(ChannelIndex) ! GHz
       Wavelength_m = POINT_01 / GHz_to_Inverse_cm( Frequency )
    ELSE IF( SpcCoeff_IsInfraredSensor(SC(SensorIndex)) ) THEN
       Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex) ! 1/cm
       Wavelength_m = POINT_01 / Wavenumber
    END IF

    Temp_K = 273.15_fp
    perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, & ! Input
                                          Temp_K) ! Input

    perm_re = REAL(REAL(perm))  ! Double REAL is required to avoud issues with GNU Fortran
    perm_im = REAL(AIMAG(perm))
    ! perm = cmplx(perm_re, perm_im, FP)
    !kw = (perm**2 - ONE )/(perm**2 + TWO)
    kw = (perm - ONE )/(perm + ONE)
    Kw_2 = ABS(kw)**TWO      
    P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)

    ! This avoids some outliers when transmittance is low
    ! don't multiple OD by dS
    Optical_Depth = AtmOptics%optical_depth / dS_m
    DO k = 1, n_layers
       temp_sum = SUM(Optical_Depth(1:k))
       Transmittance(k) = EXP(-TWO * temp_sum)
    END DO
    
    ! Divide by dZ_m as water conntent calculations are based on
    ! vertical thickness of layers
    Reflectivity =  P1 * (AtmOptics%Backscat_Coefficient / dZ_m) ! mm^6 m^-3
    Reflectivity_Attenuated = Transmittance * Reflectivity ! mm^6 m^-3
    ! Convert the unit to dBz
    WHERE (Reflectivity .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity       = TEN * LOG10(Reflectivity) ! [dBZ]  
        RTSolution%ReflectivityLinear = MM6_TO_CM6 * Reflectivity  ! cm^6m^-3  
    ELSE WHERE
        RTSolution%Reflectivity       = MISSING_REFL
        RTSolution%ReflectivityLinear = MISSING_REFL
    END WHERE

    ! Convert the unit to dBz
    ! Note that Reflectivity can be greater than zero but Reflectivity_Attenuated
    ! can be still zero if Transmittance is zero
    WHERE (Reflectivity_Attenuated .GT.  REFLECTIVITY_THRESHOLD)
        RTSolution%Reflectivity_Attenuated = TEN * LOG10(Reflectivity_Attenuated)      ! dBZ
        RTSolution%Reflectivity_AttenuatedLinear = MM6_TO_CM6 * Reflectivity_attenuated ! cm^6m^-3   
    ELSE WHERE
        RTSolution%Reflectivity_Attenuated       = MISSING_REFL
        RTSolution%Reflectivity_AttenuatedLinear = MISSING_REFL
    END WHERE
    
  END SUBROUTINE CRTM_Compute_Reflectivity

SUBROUTINE CRTM_Compute_Reflectivity_TL(Atm, &
                                       AtmOptics, &
                                       AtmOptics_TL, &
                                       GeometryInfo, &
                                       SensorIndex, &
                                       ChannelIndex, &
                                       RTSolution_TL)
  TYPE(CRTM_Atmosphere_type),     INTENT(IN)     :: Atm
  TYPE(CRTM_AtmOptics_type),      INTENT(IN)     :: AtmOptics, AtmOptics_TL
  TYPE(CRTM_GeometryInfo_type),   INTENT(IN)     :: GeometryInfo
  INTEGER,                        INTENT(IN)     :: SensorIndex, ChannelIndex
  TYPE(CRTM_RTSolution_type),     INTENT(IN OUT) :: RTSolution_TL

  ! Locals
  REAL(fp) :: Frequency, Wavenumber, Wavelength_m
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity, Reflectivity_Attenuated
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity_TL, Reflectivity_Attenuated_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear, Reflectivity_AttenuatedLinear
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear_TL, Reflectivity_AttenuatedLinear_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance, Transmittance_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Optical_Depth, Optical_Depth_TL
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: P1, Kw_2, perm_re, perm_im
  REAL(fp), DIMENSION(0:AtmOptics%n_Layers) :: Height
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: dZ_m, dS_m
  COMPLEX :: perm(AtmOptics%n_Layers)
  INTEGER :: k, n_Layers
  REAL(fp) :: temp_sum

  n_Layers = Atm%n_Layers

  ! Compute layer interface heights if missing
  IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
     Height = Calculate_Height(Atm)
  ELSE
     Height = Atm%Height
  END IF

  dZ_m = Height(0:n_Layers-1) - Height(1:n_Layers)
  dS_m = dZ_m / GeometryInfo%Cosine_Sensor_Zenith

  ! Compute wavelength
  IF (SpcCoeff_IsMicrowaveSensor(SC(SensorIndex))) THEN
     Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
     Wavelength_m = POINT_01 / GHz_to_Inverse_cm(Frequency)
  ELSE IF (SpcCoeff_IsInfraredSensor(SC(SensorIndex))) THEN
     Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
     Wavelength_m = POINT_01 / Wavenumber
  END IF

  ! Compute permittivity and related factor
  perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, Atm%Temperature)
  perm_re = REAL(perm)
  perm_im = AIMAG(perm)
  Kw_2 = ABS((perm_re - ONE) / (perm_re + ONE)) ** TWO
  P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)

  ! Forward (base state) calculations
  Optical_Depth = AtmOptics%optical_depth / dS_m
  DO k = 1, n_Layers
     temp_sum = SUM(Optical_Depth(1:k))
     Transmittance(k) = EXP(-TWO * temp_sum)
  END DO

  Reflectivity = P1 * (AtmOptics%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated = Transmittance * Reflectivity

  ! Tangent-linear calculations
  Optical_Depth_TL = AtmOptics_TL%optical_depth / dS_m
  DO k = 1, n_Layers
     temp_sum = SUM(Optical_Depth_TL(1:k))
     Transmittance_TL(k) = -TWO * Transmittance(k) * temp_sum
  END DO

  Reflectivity_TL = P1 * (AtmOptics_TL%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated_TL = Transmittance * Reflectivity_TL + Reflectivity * Transmittance_TL

  ! Linear reflectivity TL
  ReflectivityLinear = MM6_TO_CM6 * Reflectivity
  ReflectivityLinear_TL = MM6_TO_CM6 * Reflectivity_TL
  Reflectivity_AttenuatedLinear = MM6_TO_CM6 * Reflectivity_Attenuated
  Reflectivity_AttenuatedLinear_TL = MM6_TO_CM6 * Reflectivity_Attenuated_TL

  ! Convert to dBZ tangent linear if above threshold
  WHERE (Reflectivity > REFLECTIVITY_THRESHOLD)
     Reflectivity_TL = TEN * Reflectivity_TL / (Reflectivity * LOG(TEN))
     RTSolution_TL%Reflectivity = Reflectivity_TL
     RTSolution_TL%ReflectivityLinear = ReflectivityLinear_TL
  ELSEWHERE
     RTSolution_TL%Reflectivity = ZERO
     RTSolution_TL%ReflectivityLinear = ZERO
  END WHERE

  WHERE (Reflectivity_Attenuated > REFLECTIVITY_THRESHOLD)
     Reflectivity_Attenuated_TL = TEN * Reflectivity_Attenuated_TL / (Reflectivity_Attenuated * LOG(TEN))
     RTSolution_TL%Reflectivity_Attenuated = Reflectivity_Attenuated_TL
     RTSolution_TL%Reflectivity_AttenuatedLinear = Reflectivity_AttenuatedLinear_TL
  ELSEWHERE
     RTSolution_TL%Reflectivity_Attenuated = ZERO
     RTSolution_TL%Reflectivity_AttenuatedLinear = ZERO
  END WHERE

END SUBROUTINE CRTM_Compute_Reflectivity_TL

!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Compute_Reflectivity_AD
!
! PURPOSE:
!       Subroutine to calculate the adjoint reflectivity for active instruments.
!       Includes both logarithmic (dBZ) and linear reflectivity adjoints:
!       ReflectivityLinear_AD and Reflectivity_AttenuatedLinear_AD.
!
! CALLING SEQUENCE:
!        CALL CRTM_Compute_Reflectivity_AD(Atm, &
!                                          AtmOptics     , &  ! Input
!                                          RTSolution    , &  ! Input
!                                          GeometryInfo  , &  ! Input
!                                          SensorIndex   , &  ! Input
!                                          ChannelIndex  , &  ! Input
!                                          AtmOptics_AD  , &  ! Input/Output
!                                          RTSolution_AD )    ! Input/Output
!
! INPUT ARGUMENTS:
!       Atm:            Structure containing the atmospheric state data.
!                       TYPE:       CRTM_Atmosphere_type
!                       ATTRIBUTES: INTENT(IN)
!
!       AtmOptics:      Structure containing the combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols.
!                       TYPE:       CRTM_AtmOptics_type
!                       ATTRIBUTES: INTENT(IN)
!
!       RTSolution:     Structure containing the solution to the RT equation
!                       for the given inputs.
!                       TYPE:       CRTM_RTSolution_type
!                       ATTRIBUTES: INTENT(IN)
!
!       GeometryInfo:   Structure containing the view geometry data.
!                       TYPE:       CRTM_GeometryInfo_type
!                       ATTRIBUTES: INTENT(IN)
!
!       SensorIndex:    Sensor index id used to access shared coefficient data.
!                       TYPE:       INTEGER
!                       ATTRIBUTES: INTENT(IN)
!
!       ChannelIndex:   Channel index id used to access shared coefficient data.
!                       TYPE:       INTEGER
!                       ATTRIBUTES: INTENT(IN)
!
! OUTPUT ARGUMENTS:
!       RTSolution_AD:  Structure containing the RT solution adjoint inputs.
!                       TYPE:       CRTM_RTSolution_type
!                       ATTRIBUTES: INTENT(IN OUT)
!
!       AtmOptics_AD:   Structure containing the adjoint combined atmospheric
!                       optical properties for gaseous absorption, clouds,
!                       and aerosols, including:
!                       - ReflectivityLinear_AD
!                       - Reflectivity_AttenuatedLinear_AD
!                       TYPE:       CRTM_AtmOptics_type
!                       ATTRIBUTES: INTENT(IN OUT)
!
!--------------------------------------------------------------------------------
SUBROUTINE CRTM_Compute_Reflectivity_AD(Atm, &
                                       AtmOptics, &
                                       RTSolution, &
                                       GeometryInfo, &
                                       SensorIndex, &
                                       ChannelIndex, &
                                       AtmOptics_AD, &
                                       RTSolution_AD)
  !-----------------------------------------------------------------------
  !  Compute adjoint sensitivities of reflectivity (linear or dBZ)
  !  depending on which adjoint variable is set to one in RTSolution_AD.
  !
  !  Supported adjoint targets:
  !    - RTSolution_AD%Reflectivity_Attenuated        (logarithmic dBZ)
  !    - RTSolution_AD%Reflectivity_AttenuatedLinear  (linear cm^6/m^3)
  !    - RTSolution_AD%Reflectivity                   (logarithmic dBZ)
  !    - RTSolution_AD%ReflectivityLinear             (linear cm^6/m^3)
  !
  !-----------------------------------------------------------------------

  ! Arguments
  TYPE(CRTM_Atmosphere_type),      INTENT(IN)     :: Atm
  TYPE(CRTM_GeometryInfo_type),    INTENT(IN)     :: GeometryInfo
  INTEGER,                         INTENT(IN)     :: SensorIndex
  INTEGER,                         INTENT(IN)     :: ChannelIndex
  TYPE(CRTM_AtmOptics_type),       INTENT(IN)     :: AtmOptics
  TYPE(CRTM_RTSolution_type),      TARGET, INTENT(IN) :: RTSolution
  TYPE(CRTM_AtmOptics_type),       INTENT(IN OUT) :: AtmOptics_AD
  TYPE(CRTM_RTSolution_type),      INTENT(IN OUT) :: RTSolution_AD

  ! Locals
  INTEGER :: k, j, n_Layers
  REAL(fp) :: Frequency, Wavenumber, Wavelength_m
  REAL(fp) :: temp_sum

  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Reflectivity, Reflectivity_Attenuated
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear, Reflectivity_AttenuatedLinear
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance, Optical_Depth
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: P1, Kw_2, perm_re, perm_im
  REAL(fp), DIMENSION(0:AtmOptics%n_Layers) :: Height
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: dZ_m, dS_m
  COMPLEX :: perm(AtmOptics%n_Layers)

  ! Adjoint variables
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: Transmittance_AD, Reflectivity_AD, Reflectivity_Attenuated_AD, Optical_Depth_AD
  REAL(fp), DIMENSION(AtmOptics%n_Layers) :: ReflectivityLinear_AD, Reflectivity_AttenuatedLinear_AD

  !-----------------------------------------------------------------------
  ! Initialization
  !-----------------------------------------------------------------------
  n_Layers = Atm%n_Layers

  Transmittance_AD           = ZERO
  Reflectivity_AD            = ZERO
  Reflectivity_Attenuated_AD = ZERO
  ReflectivityLinear_AD            = ZERO
  Reflectivity_AttenuatedLinear_AD = ZERO
  Optical_Depth_AD           = ZERO
  AtmOptics_AD%Backscat_Coefficient = ZERO
  AtmOptics_AD%optical_depth       = ZERO

  !-----------------------------------------------------------------------
  ! Heights and geometry
  !-----------------------------------------------------------------------
  IF (ALL(Atm%Height .LT. EPSILON_FP)) THEN
      Height = Calculate_Height(Atm)
  ELSE
      Height = Atm%Height
  END IF

  dZ_m = Height(0:n_Layers-1) - Height(1:n_Layers)
  dS_m = dZ_m / GeometryInfo%Cosine_Sensor_Zenith

  !-----------------------------------------------------------------------
  ! Frequency / wavelength
  !-----------------------------------------------------------------------
  IF (SpcCoeff_IsMicrowaveSensor(SC(SensorIndex))) THEN
      Frequency = SC(SensorIndex)%Frequency(ChannelIndex)
      Wavelength_m = POINT_01 / GHz_to_Inverse_cm(Frequency)
  ELSE IF (SpcCoeff_IsInfraredSensor(SC(SensorIndex))) THEN
      Wavenumber = SC(SensorIndex)%Wavenumber(ChannelIndex)
      Wavelength_m = POINT_01 / Wavenumber
  END IF

  !-----------------------------------------------------------------------
  ! Permittivity and forward reflectivity computation
  !-----------------------------------------------------------------------
  perm = Water_Permittivity_Turner_2016(Frequency * 1.0d9, Atm%Temperature)
  perm_re = REAL(perm)
  perm_im = AIMAG(perm)
  Kw_2 = ABS((perm_re - ONE) / (perm_re + ONE))**TWO

  Optical_Depth = AtmOptics%optical_depth / dS_m

  DO k = 1, n_Layers
      temp_sum = SUM(Optical_Depth(1:k))
      Transmittance(k) = EXP(-TWO * temp_sum)
  END DO

  P1 = (M6_MM6 * Wavelength_m**4.0_fp) / (PI**5.0_fp * Kw_2)
  Reflectivity = P1 * (AtmOptics%Backscat_Coefficient / dZ_m)
  Reflectivity_Attenuated = Transmittance * Reflectivity

  ! Linear reflectivities (cm^6/m^3)
  ReflectivityLinear = MM6_TO_CM6 * Reflectivity
  Reflectivity_AttenuatedLinear = MM6_TO_CM6 * Reflectivity_Attenuated

  !-----------------------------------------------------------------------
  ! ADJOINT PROPAGATION: determine which output adjoint is active
  !-----------------------------------------------------------------------

  ! === Case 1: Logarithmic attenuated reflectivity (dBZ)
  IF (ANY(RTSolution_AD%Reflectivity_Attenuated /= ZERO)) THEN
      WHERE (Reflectivity_Attenuated > REFLECTIVITY_THRESHOLD)
        Reflectivity_Attenuated_AD = Reflectivity_Attenuated_AD + &
             TEN * RTSolution_AD%Reflectivity_Attenuated / (Reflectivity_Attenuated * LOG(TEN))
      END WHERE

  ! === Case 2: Linear attenuated reflectivity (cm^6/m^3)
  ELSEIF (ANY(RTSolution_AD%Reflectivity_AttenuatedLinear /= ZERO)) THEN
      Reflectivity_AttenuatedLinear_AD = Reflectivity_AttenuatedLinear_AD + RTSolution_AD%Reflectivity_AttenuatedLinear
      Reflectivity_Attenuated_AD = Reflectivity_Attenuated_AD + MM6_TO_CM6 * Reflectivity_AttenuatedLinear_AD
  END IF


  ! === Case 3: Logarithmic reflectivity (dBZ)
  IF (ANY(RTSolution_AD%Reflectivity /= ZERO)) THEN
      WHERE (Reflectivity > REFLECTIVITY_THRESHOLD)
        Reflectivity_AD = Reflectivity_AD + &
             TEN * RTSolution_AD%Reflectivity / (Reflectivity * LOG(TEN))
      END WHERE

  ! === Case 4: Linear reflectivity (cm^6/m^3)
  ELSEIF (ANY(RTSolution_AD%ReflectivityLinear /= ZERO)) THEN
      ReflectivityLinear_AD = ReflectivityLinear_AD + RTSolution_AD%ReflectivityLinear
      Reflectivity_AD = Reflectivity_AD + MM6_TO_CM6 * ReflectivityLinear_AD
  END IF

  !-----------------------------------------------------------------------
  ! Common adjoint propagation for reflectivity relationships
  !-----------------------------------------------------------------------
  Transmittance_AD = Transmittance_AD + Reflectivity * Reflectivity_Attenuated_AD
  Reflectivity_AD  = Reflectivity_AD  + Transmittance * Reflectivity_Attenuated_AD

  ! Backscatter coefficient adjoint
  AtmOptics_AD%Backscat_Coefficient = AtmOptics_AD%Backscat_Coefficient + Reflectivity_AD * P1 / dZ_m  

  !-----------------------------------------------------------------------
  ! Adjoint of transmittance → optical depth
  !-----------------------------------------------------------------------
  DO k = n_Layers, 1, -1
      DO j = 1, k
          Optical_Depth_AD(j) = Optical_Depth_AD(j) - TWO * Transmittance(k) * Transmittance_AD(k)
      END DO
  END DO

  !-----------------------------------------------------------------------
  ! Adjoint of optical depth normalization
  !-----------------------------------------------------------------------
  AtmOptics_AD%optical_depth = AtmOptics_AD%optical_depth + Optical_Depth_AD / dS_m

END SUBROUTINE CRTM_Compute_Reflectivity_AD

END MODULE CRTM_Active_Sensor
