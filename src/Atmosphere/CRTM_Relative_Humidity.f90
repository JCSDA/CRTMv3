!
! CRTM_RH_Calculation
!
! Compute relative humidity based on CRTM Atmosphere sturctures
!
! CREATION HISTORY:
!       Written by:     Cheng Dang, 17-Aug-2022
!                       dangch@ucar.edu
!
! Contains submodules from Profile_Utility package:
!    MR_PPMV.f90 :  PPMV_to_MR
!    RH_MR.f90   :  MR_to_RH [RH unit in one not %]
!    Atmospheric_Properties.f90:  Saturation_Mixing_Ratio, SVP_Water, SVP_Ice
!

MODULE CRTM_Relative_Humidity
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE Type_Kinds            , ONLY: fp
  USE Message_Handler       , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, Display_Message
  USE CRTM_Parameters       , ONLY: ZERO


  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! Module procedures
  PUBLIC :: PPMV_to_MR, MR_to_RH

  ! -----------------
  ! Module parameters
  ! -----------------
  ! Parameters for unit conversion
  REAL(fp), PARAMETER :: CELSIUS_TO_KELVIN = 273.15_fp
  REAL(fp), PARAMETER :: KG_TO_G           = 1.0e+03_fp
  REAL(fp), PARAMETER :: TO_PERCENT        = 1.0e+02_fp
  REAL(fp), PARAMETER :: PPMV_TO_PPV       = 1.0e-06_fp
  REAL(fp), PARAMETER :: PPMV_TO_MR_SCALE_FACTOR = KG_TO_G * PPMV_TO_PPV
  !
  ! Maximum number of molecular species (for MW_H2O only in this module)
  INTEGER,  PARAMETER :: MAX_N_MOLECULAR_SPECIES = 32
  ! Molecular weights of first seven HITRAN molecular species
  REAL(fp), PUBLIC, PARAMETER :: MW_H2O = 18.01528_fp
  REAL(fp), PARAMETER :: MW_CO2 = 44.00950_fp
  REAL(fp), PARAMETER :: MW_O3  = 47.99820_fp
  REAL(fp), PARAMETER :: MW_N2O = 44.01288_fp
  REAL(fp), PARAMETER :: MW_CO  = 28.01010_fp
  REAL(fp), PARAMETER :: MW_CH4 = 16.04246_fp
  REAL(fp), PARAMETER :: MW_O2  = 31.99880_fp
  REAL(fp), PARAMETER :: MW_N2  = 28.01348_fp
  ! Weights of all 32 HITRAN molecular species
  REAL(fp), PARAMETER :: MOLECULAR_WEIGHT(MAX_N_MOLECULAR_SPECIES) = &
    (/       MW_H2O,        MW_CO2,         MW_O3,        MW_N2O, &
             MW_CO ,        MW_CH4,         MW_O2,   30.00614_fp, &
        64.06480_fp,   46.00554_fp,   17.03056_fp,   63.01288_fp, &
        17.00734_fp,   20.00634_fp,   36.46064_fp,   80.91194_fp, &
       127.91241_fp,   51.45210_fp,   60.07610_fp,   30.02598_fp, &
        52.46004_fp,        MW_N2 ,   27.02538_fp,   50.48722_fp, &
        34.01468_fp,   26.03728_fp,   30.06904_fp,   33.99758_fp, &
        66.00690_fp,  146.05643_fp,   34.08188_fp,   46.02538_fp /)

  ! Average molecular weight of dry air
  REAL(fp), PARAMETER :: MW_DRYAIR = 28.9648_fp
  ! Ratio of water vapor and dry air weights for conversion routines
  REAL(fp), PARAMETER :: EPS       = MW_H2O / MW_DRYAIR
  ! Minimum pressure for saturation mixing ratio calculation
  REAL(fp), PARAMETER :: MIN_SMR_PRESSURE = 50.0_fp

  ! Coefficients for saturation vapor pressure over water
  INTEGER,  PARAMETER :: N_SVPW_COEFFICIENTS = 8
  REAL(fp), PARAMETER :: SVPW_COEFFICIENTS(0:N_SVPW_COEFFICIENTS) = &
    (/-3.21582393e-16_fp, 3.79534310e-14_fp, 7.02620698e-11_fp, &
       2.03154182e-08_fp, 2.99291081e-06_fp, 2.64224321e-04_fp, &
       1.43177157e-02_fp, 4.44606896e-01_fp, 6.11583699e+00_fp /)
  ! ...Valid temperature range
  REAL(fp), PARAMETER :: MIN_SVPW_TEMPERATURE = 188.15_fp
  REAL(fp), PARAMETER :: MAX_SVPW_TEMPERATURE = 343.15_fp

  ! Coefficients for saturation vapor pressure over ice
  INTEGER,  PARAMETER :: N_SVPI_COEFFICIENTS = 8
  REAL(fp), PARAMETER :: SVPI_COEFFICIENTS(0:N_SVPI_COEFFICIENTS) = &
    (/ 1.61444444e-15_fp, 1.05785160e-12_fp, 3.07839583e-10_fp, &
       5.21693933e-08_fp, 5.65392987e-06_fp, 4.02737184e-04_fp, &
       1.84672631e-02_fp, 4.99320233e-01_fp, 6.09868993e+00_fp /)
  ! ...Valid temperature range
  REAL(fp), PARAMETER :: MIN_SVPI_TEMPERATURE = 183.15_fp
  REAL(fp), PARAMETER :: MAX_SVPI_TEMPERATURE = 273.15_fp


CONTAINS



!################################################################################
!################################################################################
!##                                                                            ##
!##                         ## PUBLIC MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################
!------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       PPMV_to_MR
!
! PURPOSE:
!       Subroutine to convert gas concentrations from volume mixing
!       ratio in ppmv to mass mixing ratio in g/kg.
!
! CALLING SEQUENCE:
!       CALL PPMV_to_MR( ppmv                     , &  ! Input
!                        Mixing_Ratio             , &  ! Output
!                        Molecule_ID = Molecule_ID  )  ! Optional input
!
! INPUTS:
!       ppmv:             Volume mixing ratio of gas.
!                         Must be > or = 0.0
!                         UNITS:      ppmv
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or any rank
!                         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Mixing_Ratio:     Mass mixing ratio of gas.
!                         Set to 0.0 if input < 0.0
!                         UNITS:      g/kg
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as ppmv argument
!                         ATTRIBUTES: INTENT(IN)
!
! OPTIONAL INPUT ARGUMENTS:
!       Molecule_ID:      HITRAN molecular designation identifying the
!                         molecule for which the concentration units
!                         conversion is required. If not specified, the
!                         default value is that for water vapor.
!                         Valid values are:
!                           1: H2O       9: SO2      17: HI       25: H2O2
!                           2: CO2      10: NO2      18: ClO      26: C2H2
!                           3: O3       11: NH3      19: OCS      27: C2H6
!                           4: N2O      12: HNO3     20: H2CO     28: PH3
!                           5: CO       13: OH       21: HOCl     29: COF2
!                           6: CH4      14: HF       22: N2       30: SF6
!                           7: O2       15: HCl      23: HCN      31: H2S
!                           8: NO       16: HBr      24: CH3Cl    32: HCOOH
!                         Output is set to zero if an invalid Molecule_Id
!                         is supplied.
!                         UNITS:      N/A
!                         TYPE:       INTEGER
!                         DIMENSION:  Scalar
!                         ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! PROCEDURE:
!       To convert ppmv to mixing ratio, the following is used:
!
!                                          MW(MOL)
!         mr(MOL) = 0.001 . ppmv(MOL) . -------------
!                                        MW(Dry Air)
!
!       where MW(Dry Air) = Average molecular weight of dry air
!             MW(MOL)     = Molecular weight of the gas in question.
!
!       The factor of 0.001 derives from the product of the g/g to g/kg
!       scale factor (1000) and the "parts-per-million" to "parts-per"
!       scale factor (1.0e-06)
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 03-May-2000
!                       paul.vandelst@noaa.gov
!:sdoc-:
!------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE PPMV_to_MR( &
    ppmv        , &  ! Input
    Mixing_Ratio, &  ! Output
    Molecule_ID   )  ! Optional input
    ! Arguments
    REAL(fp),          INTENT(IN)  :: ppmv
    REAL(fp),          INTENT(OUT) :: Mixing_Ratio
    INTEGER, OPTIONAL, INTENT(IN)  :: Molecule_ID
    ! Local variables
    INTEGER :: Id

    ! Error checks
    ! ...Zero output for -ve input
    IF ( ppmv < ZERO ) THEN
      Mixing_Ratio = ZERO
      RETURN
    ENDIF
    ! ...Zero output for invalid id
    IF ( PRESENT(Molecule_ID) ) THEN
      IF ( Molecule_ID < 1 .OR. Molecule_ID > MAX_N_MOLECULAR_SPECIES ) THEN
        Mixing_Ratio = ZERO
        RETURN
      END IF
      Id = Molecule_ID
    ELSE
      Id = 1  ! Default value is for water vapor
    END IF

    ! Convert volume to mass mixing ratio
    Mixing_Ratio = ppmv * PPMV_TO_MR_SCALE_FACTOR * MOLECULAR_WEIGHT(Id) / MW_DRYAIR

  END SUBROUTINE PPMV_to_MR

!------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       MR_to_RH
!
! PURPOSE:
!       Subroutine to convert water vapor mass mixing ratio
!       to relative humidity
!
! CALLING SEQUENCE:
!       CALL MR_to_RH( Pressure                         , &  ! Input
!                      Temperature                      , &  ! Input
!                      Mixing_Ratio                     , &  ! Input
!                      Relative_Humidity                , &  ! Output
!                      Ice_Temperature = Ice_Temperature, &  ! Optional input
!                      Min_Pressure    = Min_Pressure     )  ! Optional input
!
! INPUTS:
!       Pressure:             Total atmospheric pressure.
!                             Must be > 0.
!                             UNITS:      hectoPascals, hPa
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar or Rank-1 (K x 1)
!                             ATTRIBUTES: INTENT(IN)
!
!       Temperature:          Atmospheric temperature.
!                             Must be > 0.
!                             UNITS:      Kelvin, K
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Same as Pressure
!                             ATTRIBUTES: INTENT(IN)
!
!       Mixing_Ratio:         Water vapor mixing ratio.
!                             Must be > 0.
!                             UNITS:      g/kg
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Same as Pressure
!                             ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Relative_Humidity:    Relative humidity.
!                             Set to zero for invalid input.
!                             UNITS:      %
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Same as Pressure
!                             ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUT ARGUMENTS:
!       Ice_Temperature:      Temperature below which the saturation vapor
!                             pressure over ice is used in the conversion.
!                             By default, only the saturation vapor pressure
!                             over water is used.
!                             UNITS:      Kelvin, K
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Min_Pressure:         Pressure value below which the saturation
!                             mixing ratio is not calculated. The default
!                             is 50mb. Saturation mixing ratios below the
!                             minimum pressure are set to zero. This is
!                             because at pressures less than 50mb, the
!                             saturation vapour Pressure, which is based
!                             only on temperature, can exceed the total
!                             air Pressure.
!                             UNITS:      hectoPascals, hPa
!                             TYPE:       REAL(fp)
!                             DIMENSION:  Scalar
!                             ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! PROCEDURE:
!       Once the saturation mixing ratio is calculated the relative humidity
!       corresponding to the input mixing ratio is determined using:
!
!                                         Mixing_Ratio
!       Relative_Humidity = 100.0 * -------------------------
!                                    Saturation_Mixing_Ratio
!
! CREATION HISTORY:
!       Written by:     Paul van Delst, 02-Mar-1999
!                       paul.vandelst@noaa.gov
!:sdoc-:
!------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE MR_to_RH( &
    P              , &  ! Input
    T              , &  ! Input
    mr             , &  ! Input
    rh             , &  ! Output
    Ice_Temperature, &  ! Optional Input
    Min_Pressure     )  ! Optional Input
    ! Arguments
    REAL(fp),           INTENT(IN)  :: P
    REAL(fp),           INTENT(IN)  :: T
    REAL(fp),           INTENT(IN)  :: mr
    REAL(fp),           INTENT(OUT) :: rh
    REAL(fp), OPTIONAL, INTENT(IN)  :: Ice_Temperature
    REAL(fp), OPTIONAL, INTENT(IN)  :: Min_Pressure
    ! Local variables
    REAL(fp) :: smr

    ! Setup
    IF ( P < ZERO .OR. T < ZERO .OR. mr < ZERO ) THEN
      rh = ZERO
      RETURN
    ENDIF


    ! Calculate saturation mixing ratio in g/kg
    CALL Saturation_Mixing_Ratio( P, &
                                  T, &
                                  smr, &
                                  Ice_Temperature=Ice_Temperature, &
                                  Min_Pressure   =Min_Pressure )

    ! Calculate relative humidity from 0 to 1
    IF ( smr > ZERO ) THEN
      rh =  mr / smr
    ELSE
      rh = ZERO
    END IF

  END SUBROUTINE MR_to_RH

!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       Saturation_Mixing_Ratio
!
! PURPOSE:
!       Subroutine to calculate the saturation mixing ratio for
!       a given pressure and temperature
!
! CALLING SEQUENCE:
!       CALL Saturation_Mixing_Ratio( Pressure                         , &  ! Input
!                                     Temperature                      , &  ! Input
!                                     Mixing_Ratio                     , &  ! Output
!                                     Ice_Temperature = Ice_Temperature, &  ! Optional input
!                                     Min_Pressure    = Min_Pressure     )  ! Optional input
!
! INPUTS:
!       Pressure:          Total atmospheric pressure.
!                          Valid pressures are > 50hPa.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar or any rank.
!                          ATTRIBUTES: INTENT(IN)
!
!       Temperature:       Atmospheric Temperature.
!                          Valid temperature ranges for saturation vapor
!                          pressure calculation are:
!                            Over ice:   183K - 273K (-90C - 0C).
!                            Over water: 188K - 343K (-85C - +70C).
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Mixing_Ratio:      The saturation mixing ratio for the supplied
!                          pressure and temperature.
!                          Value is set to zero for invalid input.
!                          UNITS:      g/kg
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Same as Pressure
!                          ATTRIBUTES: INTENT(OUT)
!
! OPTIONAL INPUTS:
!       Ice_Temperature:   Temperature below which the saturation vapor
!                          pressure over ice is used in the conversion.
!                          By default, only the saturation vapor pressure
!                          over water is used.
!                          UNITS:      Kelvin, K
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
!       Min_Pressure:      Pressure value below which the saturation
!                          mixing ratio is not calculated. The default, and
!                          absolute, minimum value used in this routine is
!                          50hPa. Saturation mixing ratios at pressures
!                          less than the minimum pressure are set to zero.
!                          This is because at pressures less than 50mb, the
!                          saturation vapour pressure, which is based only on
!                          temperature, can exceed the total air pressure.
!                          UNITS:      hectoPascals, hPa
!                          TYPE:       REAL(fp)
!                          DIMENSION:  Scalar
!                          ATTRIBUTES: OPTIONAL, INTENT(IN)
!
! PROCEDURE:
!       The saturation mixing ratio can be defined as:
!
!                rho_ws
!          ws = --------     .....(1)
!                rho_d
!
!       where rho_ws = the partial density of water vapour required to
!                      saturate air with respect to water at a Temperature, T
!             rho_d  = the partial density of dry air.
!
!       Equation (1) can be rewritten as:
!
!                   es
!               ---------
!                R_w . T
!         ws = ------------
!                p - es
!               ---------
!                R_d . T
!
!               R_d       es
!            = ----- . --------
!               R_w     p - es
!
!               M_w       es
!            = ----- . --------     .....(2)
!               M_d     p - es
!
!       where M_w = molecular weight of water
!             M_d = molecular weight of dry air
!             es  = water vapor partial pressure
!             p   = total air pressure
!             R_d = gas constant for dry air
!             R_w = gas constant for water vapor
!
!       The units of equation (2) are:
!
!               g     hPa
!         ws = --- . -----
!               g     hPa
!
!                      g
!            = 1000.0 ----
!                      kg
!
!       A factor of 1000 is used to return values in units of g/kg.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE Saturation_Mixing_Ratio( &
    Pressure       , &  ! Input
    Temperature    , &  ! Input
    Mixing_Ratio   , &  ! Output
    Ice_Temperature, &  ! Optional Input
    Min_Pressure     )  ! Optional Input
    ! Arguments
    REAL(fp),           INTENT(IN)  :: Pressure
    REAL(fp),           INTENT(IN)  :: Temperature
    REAL(fp),           INTENT(OUT) :: Mixing_Ratio
    REAL(fp), OPTIONAL, INTENT(IN)  :: Ice_Temperature
    REAL(fp), OPTIONAL, INTENT(IN)  :: Min_Pressure
    ! Local variables
    REAL(fp) :: Pmin, Tmin, Tice
    REAL(fp) :: svp
    REAL(fp) :: dp

    ! Setup
    ! ...Check optional arguments
    IF ( PRESENT(Min_Pressure) ) THEN
      Pmin = MAX(Min_Pressure, MIN_SMR_PRESSURE)
    ELSE
      Pmin = MIN_SMR_PRESSURE
    END IF
    IF ( PRESENT(Ice_Temperature) ) THEN
      Tice = Ice_Temperature
      Tmin = MIN_SVPI_TEMPERATURE
    ELSE
      Tice = ZERO
      Tmin = MIN_SVPW_TEMPERATURE
    END IF
    ! ...Check input
    IF ( Pressure < Pmin .OR. Temperature < Tmin ) THEN
      Mixing_Ratio = ZERO
      RETURN
    ENDIF


    ! Calculate saturation vapor pressure
    IF ( Temperature > Tice ) THEN
      CALL SVP_Water( Temperature, svp )
    ELSE
      CALL SVP_Ice( Temperature, svp )
    END IF


    ! Calculate saturation mixing ratio only if the
    ! total pressure is greater than the saturation
    ! vapor pressure.
    dp = Pressure - svp
    IF ( dp > ZERO ) THEN
      Mixing_Ratio = KG_TO_G * EPS * svp / dp
    ELSE
      Mixing_Ratio = ZERO
    END IF

  END SUBROUTINE Saturation_Mixing_Ratio


!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       SVP_Water
!
! PURPOSE:
!       Subroutine to calculate the saturation vapor pressure
!       over water.
!
! CALLING SEQUENCE:
!       CALL SVP_Water( Temperature   , &  ! Input
!                       Vapor_Pressure  )  ! Output
!
! INPUTS:
!       Temperature:      Temperatures for which the saturation vapor
!                         pressure is required.
!                         Valid temperature range is 188K - 343K (-85C - +70C).
!                         UNITS:      Kelvin, K
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or any rank.
!                         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Vapor_Pressure:   The saturation vapor pressure over water.
!                         Value is set to zero if input temperatures are
!                         outside the valid range.
!                         UNITS:      hectoPascals, hPa
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as input Temperature
!                         ATTRIBUTES: INTENT(OUT)
!
! PROCEDURE:
!       Flatau,P.J., R.L.Walko, and W.R.Cotton, 1992: "Polynomial fits to
!         saturation vapor pressure", J.Appl.Met., v31, pp1507-1513
!
!                           __ N
!                          \            i
!         SVP_Water = c0 +  >   c(i) . T
!                          /__
!                             i=1
!
!       where the c(i) are the relative error norm coefficients obtained
!       from the reference above.
!
!       Horner's method is used to evaluate the above polynomial.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE SVP_Water( &
    Temperature,   &  ! Input
    Vapor_Pressure )  ! Output
    ! Arguments
    REAL(fp), INTENT(IN)  :: Temperature
    REAL(fp), INTENT(OUT) :: Vapor_Pressure
    ! Local variables
    INTEGER :: i
    REAL(fp) :: T

    ! Setup
    IF ( Temperature < MIN_SVPW_TEMPERATURE .OR. Temperature > MAX_SVPW_TEMPERATURE ) THEN
      Vapor_Pressure = ZERO
      RETURN
    END IF

    ! Calculate saturation vapor pressure
    T = Temperature - CELSIUS_TO_KELVIN
    Vapor_Pressure = SVPW_COEFFICIENTS(0)
    DO i = 1, N_SVPW_COEFFICIENTS
      Vapor_Pressure = (Vapor_Pressure * T) + SVPW_COEFFICIENTS(i)
    END DO

  END SUBROUTINE SVP_Water

!--------------------------------------------------------------------------------
!:sdoc+:
! NAME:
!       SVP_Ice
!
! PURPOSE:
!       Subroutine to calculate the saturation vapor pressure
!       over ice.
!
! CALLING SEQUENCE:
!       CALL SVP_Ice( Temperature   , &  ! Input
!                     Vapor_Pressure  )  ! Output
!
! INPUTS:
!       Temperature:      Temperatures for which the saturation vapor
!                         pressure is required.
!                          Valid temperature range is 183K - 273K (-90C - 0C).
!                         UNITS:      Kelvin, K
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Scalar or any rank.
!                         ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       Vapor_Pressure:   The saturation vapor pressure over ice.
!                         Value is set to zero if input temperatures are
!                         outside the valid range.
!                         UNITS:      hectoPascals, hPa
!                         TYPE:       REAL(fp)
!                         DIMENSION:  Same as input Temperature
!                         ATTRIBUTES: INTENT(OUT)
!
! PROCEDURE:
!       Flatau,P.J., R.L.Walko, and W.R.Cotton, 1992: "Polynomial fits to
!         saturation vapor pressure", J.Appl.Met., v31, pp1507-1513
!
!                         __ N
!                        \            i
!         SVP_Ice = c0 +  >   c(i) . T
!                        /__
!                           i=1
!
!       where the c(i) are the relative error norm coefficients obtained
!       from the reference above.
!
!       Horner's method is used to evaluate the above polynomial.
!
!:sdoc-:
!--------------------------------------------------------------------------------

  ELEMENTAL SUBROUTINE SVP_Ice( &
    Temperature,   &  ! Input
    Vapor_Pressure )  ! Output
    ! Arguments
    REAL(fp), INTENT(IN)  :: Temperature
    REAL(fp), INTENT(OUT) :: Vapor_Pressure
    ! Local variables
    INTEGER  :: i
    REAL(fp) :: T

    ! Setup
    IF ( Temperature < MIN_SVPI_TEMPERATURE .OR. Temperature > MAX_SVPI_TEMPERATURE ) THEN
      Vapor_Pressure = ZERO
      RETURN
    END IF

    ! Calculate saturation vapor pressure
    T = Temperature - CELSIUS_TO_KELVIN
    Vapor_Pressure = SVPI_COEFFICIENTS(0)
    DO i = 1, N_SVPI_COEFFICIENTS
      Vapor_Pressure = (Vapor_Pressure * T) + SVPI_COEFFICIENTS(i)
    END DO

  END SUBROUTINE SVP_Ice

END MODULE CRTM_Relative_Humidity
