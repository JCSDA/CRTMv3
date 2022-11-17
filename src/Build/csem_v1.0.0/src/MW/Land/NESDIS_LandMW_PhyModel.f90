!
! NESDIS_LandMW_PhyModel
!
! Module containing the NESDIS physical NON-SNOW (bare soil, desert, vegetation-covered)
! land emissivity model of microwave 
! channels. It wraps all the available model versions and provide a general interface
! to couple with upper-level applications.
!
! This module encloses the similar forward function of the original NESDIS_LandEM_Module 
! in the earlier CRTM releases, but with several additional features in model physics and
! model implementation techniques. 
!
! Unlike the module in the ealier CRTM  releases, functions of Snow and Sea ice emissivity
! and reflectivity are not enclosed in the this module dedicated for non-snow land. Instead,
! separate modules are designed and created for the Snow and sea ice model implementation.
!
! Soil and canopy models are also implemented in their individual modules, which provides
! soil and canopy optical parameters required by this module. Since different models
! are implemented as options in Soil module (MW_Soil_Optics) and canopy module(MW_Canopy_Optics)
! users need to choose and specify the soil model and the canopy model to be used in this module.
! In other words, users may cretae their own model by different combination of model options. 
!
!
! Non-isothermal two-stream model is enclosed in this module to account for the temperature
! difference between the canopy and the underlying soil. 
!
! Tangent-linear and adjoint functions are developed for the applications in data assimilation
! and surface property retrieval systems.
!
! REFERENCES:
!       Weng, F., B. Yan, and N. Grody, 2001: "A microwave land emissivity model",
!         J. Geophys. Res., 106, 20, 115-20, 123
 
! 
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       Ming.Chen@noaa.gov
!                   



MODULE NESDIS_LandMW_PhyModel
 
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds,  ONLY: fp => CSEM_fp
  USE MW_Leaf_Optics,   ONLY: Leaf_iVar_Type => iVar_type,   &
                              CSEM_LeafMW_Optics 
  USE MW_Canopy_Optics, ONLY: Canopy_iVar_Type => iVar_type, &
                              CSEM_CanopyMW_Optics,          &
                              CSEM_CanopyMW_Optics_TL,       &
                              CSEM_CanopyMW_Optics_AD      
  USE MW_Soil_Optics,   ONLY: Soil_iVar_Type => iVar_type,   &
                              Max_Soil_Layers,               &
                              CSEM_SoilMW_Optics,            &
                              CSEM_SoilMW_Optics_TL,         &
                              CSEM_SoilMW_Optics_AD      
  USE SoilMW_roughness_Correction, ONLY:  Rsoil_iVar_Type =>iVar_type,  &      
                              CSEM_SoilMW_roughness_Correction,   &
                              CSEM_SoilMW_roughness_Correction_TL,&
                              CSEM_SoilMW_roughness_Correction_AD
  USE SoilMW_DSM_Module
  
  
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC  :: NESDIS_LandMW_Emiss
  PUBLIC  :: NESDIS_LandMW_Emiss_TL
  PUBLIC  :: NESDIS_LandMW_Emiss_AD
  
  PUBLIC  :: iVar_type
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: ZERO = 0.0_fp
  REAL(fp), PARAMETER :: ONE  = 1.0_fp
  REAL(fp), PARAMETER :: TWO  = 2.0_fp
  REAL(fp), PARAMETER :: PI   = 3.141592653589793238462643_fp

  ! First and Second Planck function constants  
  REAL(fp), PARAMETER :: C_1 = 1.191042722543248E-016_fp   
  REAL(fp), PARAMETER :: C_2 = 1.438775246065196E-002_fp   
   
  REAL(fp), PARAMETER :: EMISSH_DEFAULT = 0.25_fp
  REAL(fp), PARAMETER :: EMISSV_DEFAULT = 0.30_fp

 
  TYPE Soil_Parameters
    INTEGER  :: soil_type
    REAL(fp) :: rhos     ! density of the soil solids (2.65 g.cm^3 for solid soil material)
    REAL(fp) :: rhob     ! soil bulk volume density of the soil (1.18-1.12)
    REAL(fp) :: sand     ! sand fraction (sand + clay = 1.0)
    REAL(fp) :: clay     ! caly fraction 
  END TYPE Soil_Parameters

  TYPE Veg_Parameters
    INTEGER  :: veg_type
    REAL(fp) :: lthick   ! leaf thickness (mm)
    REAL(fp) :: langl    ! canopy dominant leaf inclination angle
    REAL(fp) :: vrho     ! Bulk density of the dry vegetation material 
    REAL(fp) :: vmge     ! leaf gravimetric water content
    REAL(fp) :: maxlai   ! maximum leaf area index
    REAL(fp) :: minlai   ! minimum  leaf area index
    REAL(fp) :: sigma    ! surface roughness formed between medium 1 and 2,
  END TYPE Veg_Parameters
 
  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE iVar_2stream
    PRIVATE
    ! Forward model input values  
    REAL(fp) :: alfa(2)            = ZERO
    REAL(fp) :: beta(2)            = ZERO 
    REAL(fp) :: gamma(2)           = ZERO
    REAL(fp) :: kk(2)              = ZERO
    REAL(fp) :: fact0              = ZERO
    REAL(fp) :: fact1              = ZERO
    REAL(fp) :: fact(2,2)          = ZERO
  END TYPE iVar_2stream

  TYPE iVar_type
    PRIVATE
    ! Forward model input values
    REAL(fp) :: Frequency         =  ZERO
    REAL(fp) :: Theta             =  ZERO
    REAL(fp) :: Tskin             =  ZERO
    REAL(fp) :: Teff              =  ZERO
    REAL(fp) :: Tsoil             =  ZERO
    REAL(fp) :: smc               =  ZERO
    REAL(fp) :: LAI               =  ZERO
    REAL(fp) :: fveg              =  ZERO
    REAL(fp) :: r21(2)            =  ZERO
    REAL(fp) :: t21(2)            =  ZERO
    REAL(fp) :: r23(2)            =  ZERO
    REAL(fp) :: ssalb(2)          =  ZERO
    REAL(fp) :: tau(2)            =  ZERO
    REAL(fp) :: g(2)              =  ZERO
    REAL(fp) :: smc_frac          =  1
    TYPE(Leaf_iVar_Type)     ::  iVar_leaf
    TYPE(Soil_iVar_Type)     ::  iVar_soil
    TYPE(Rsoil_iVar_Type)    ::  iVar_rough
    TYPE(Canopy_iVar_Type)   ::  iVar_canopy
    TYPE(iVar_2stream)       ::  iVar2
  END TYPE iVar_type
 
 
CONTAINS

!
!--------------------------------------------------------------------------------
!
! NAME:
!       NESDIS_LandMW_Emiss
!
! PURPOSE:
!       Physical simulation of microwave emissivity over non-snow land conditions with
!       non-isothermal two-stream radiative transfer model. This is the version in 
!       CRTM-REL2.1.3
!
!
! CALLING SEQUENCE:
!              IO_Status = NESDIS_LandMW_Emiss(                 &
!                                  Frequency,                   &  ! Input
!                                  Angle,                       &  ! Input
!                                  Land_Skin_Temperature,       &  ! Input
!                                  Soil_Temperature,            &  ! Input
!                                  Soil_Moisture_Content,       &  ! Input
!                                  Vegetation_Fraction,         &  ! Input
!                                  LAI,                         &  ! Input
!                                  Vegetation_Type,             &  ! Input
!                                  Soil_Type,                   &  ! Input
!                                  Emissivity_H,                &  ! Output
!                                  Emissivity_V)                   ! Output
! 	    
!
! INPUT ARGUMENTS:
!
!         Frequency                Frequency User defines
!                                  This is the "I" dimension
!                                  UNITS:      GHz
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         Angle                    The angle values in degree.
!                                  UNITS:      Degrees
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Rank-1, (I)
!                                  ATTRIBUTES: INTENT(IN)
!
!         Land_Skin_Temperature:   The land surface temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         Soil_Temperature:        The soil temperature.
!                                  UNITS:      Kelvin, K
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!
!         Soil_Moisture_Content:   The volumetric water content of the soil (0:1).
!                                  UNITS:      cm-3/cm-3
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         Vegetation_Fraction:     The vegetation fraction of the surface (0:1).
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         LAI:                     Leaf area index.
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         Vegetation_Type:         Land surface vegetation cover type (1-13)
!                                  UNITS:      N/A
!                                  TYPE:       INTEGER
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!         Soil_Type:               Soil type(1-9)
!                                  UNITS:      N/A
!                                  TYPE:       INTEGER
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(IN)
!
!
! OUTPUT ARGUMENTS:
!         Emissivity_H:            The surface emissivity at a horizontal
!                                  polarization.
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(OUT)
!
!         Emissivity_V:            The surface emissivity at a vertical polarization.
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(OUT)
!
!         iVar:                    Structure containing internal variables required for
!                                  subsequent tangent-linear or adjoint model calls.
!                                  The contents of this structure are NOT accessible
!                                  outside of this module.
!                                  UNITS:      N/A
!                                  TYPE:       iVar_type
!                                  DIMENSION:  Scalar
!                                  ATTRIBUTES: INTENT(OUT)
!
!
! CREATION HISTORY:
!
!       Writeen by:     Ming Chen, Jul 2014
!                       Ming.Chen@noaa.gov
!
!------------------------------------------------------------------------------------------------------------

  FUNCTION NESDIS_LandMW_Emiss(  &
      & Frequency,                   &  ! Input
      & Angle,                       &  ! Input
      & Land_Skin_Temperature,       &  ! Input
      & Soil_Temperature,            &  ! Input
      & Soil_Moisture_Content,       &  ! Input
      & Vegetation_Fraction,         &  ! Input
      & LAI,                         &  ! Input
      & Vegetation_Type,             &  ! Input
      & Soil_Type,                   &  ! Input
      & Emissivity_H,                &  ! Output
      & Emissivity_V,                &  ! Output
      & iVar)                        &  ! Output
    RESULT (IO_Status)
     
    ! Arguments
    REAL(fp), INTENT(IN)  :: Frequency
    REAL(fp), INTENT(IN)  :: Angle
    REAL(fp), INTENT(IN)  :: Soil_Moisture_Content
    REAL(fp), INTENT(IN)  :: Vegetation_Fraction
    REAL(fp), INTENT(IN)  :: Soil_Temperature
    REAL(fp), INTENT(IN)  :: Land_Skin_Temperature
    REAL(fp), INTENT(IN)  :: LAI
    INTEGER,  INTENT(IN)  :: Soil_Type
    INTEGER,  INTENT(IN)  :: Vegetation_Type
    REAL(fp), INTENT(OUT) :: Emissivity_V, Emissivity_H
    TYPE(iVar_type)       :: iVar
    INTEGER               :: IO_Status
    
    ! Local variables
    REAL(fp) :: theta
    REAL(fp) :: roughness
    REAL(fp) :: fveg, vlai
    REAL(fp) :: Tskin, Tsoil(1), smc(1)
    REAL(fp) :: Leaf_Refls(2), Leaf_Trans(2) 
    REAL(fp) :: r23_smooth(2), emiss(2)

    COMPLEX(fp) ::  eveg
    TYPE(Soil_Parameters) :: SoilPar 
    TYPE(Veg_Parameters)  :: VegPar 

    IO_Status = 0
    SoilPar%soil_type = Soil_Type
    CALL Load_Soil_Params("CRTM_REL213", SoilPar)
    VegPar%veg_type   = Vegetation_Type
    CALL Load_Vege_Params("CRTM_REL213", VegPar)

    iVar%Frequency       = Frequency
    iVar%theta           = Angle*PI/180.0_fp
    iVar%Tskin           = Land_Skin_Temperature
    iVar%Tsoil           = Soil_Temperature
    iVar%smc             = Soil_Moisture_Content
    iVar%LAI             = LAI
    iVar%fveg            = Vegetation_Fraction
    
    iVar%iVar_leaf%Model_Option   = 2
    iVar%iVar_Soil%Model_Option   = 2
    iVar%iVar_canopy%Model_Option = 2
    iVar%iVar_rough%Model_Option  = 2
 
    iVar%iVar_Soil%rho_b   = SoilPar%rhob
    iVar%iVar_Leaf%lthick  = VegPar%lthick
    iVar%iVar_Leaf%rhoveg  = VegPar%vrho
 
    Tskin    = iVar%tskin
    Tsoil(1) = iVar%tsoil
    smc(1)   = iVar%smc
    fveg     = iVar%fveg
   
    ! Check soil moisture content range
    smc(1) = MAX(MIN(smc(1),ONE), ZERO)
    ! Limit for vegetation fraction
    fveg   = MAX(MIN(fveg,ONE), ZERO)
   
    vlai  = iVar%LAI*fveg 
    theta = iVar%theta 
    roughness = 0.5_fp

    iVar%r21 = ZERO ; iVar%t21 = ONE

    CALL CSEM_SoilMW_Optics(frequency, theta, Tskin, Tsoil, smc, &
         SoilPar%sand, SoilPar%clay, r23_smooth, iVar%Teff, iVar%iVar_Soil)
    CALL CSEM_SoilMW_roughness_Correction(frequency, roughness, &
         iVar%theta, r23_smooth, iVar%r23, iVar%iVar_rough)
    CALL CSEM_LeafMW_Optics(frequency, angle, VegPar%vmge,  &
         Leaf_Refls, Leaf_Trans, eveg, iVar%iVar_Leaf)
    CALL CSEM_CanopyMW_Optics(vlai, Leaf_Refls, Leaf_Trans, &
         iVar%g, iVar%ssalb, iVar%tau, iVar%iVar_canopy)
    CALL Two_Stream_Solution(emiss,iVar)
 
    Emissivity_H = emiss(1) ; Emissivity_V = emiss(2) 
  
  END FUNCTION NESDIS_LandMW_Emiss


 FUNCTION NESDIS_LandMW_Emiss_TL(       &
      & Land_Skin_Temperature_TL,       &  ! Input
      & Soil_Temperature_TL,            &  ! Input
      & Soil_Moisture_Content_TL,       &  ! Input
      & Vegetation_Fraction_TL,         &  ! Input
      & Emissivity_H_TL,                &  ! Output
      & Emissivity_V_TL, iVar)          &  ! Output
    RESULT (IO_Status)
     
    ! Arguments
    REAL(fp), INTENT(IN) :: Soil_Moisture_Content_TL
    REAL(fp), INTENT(IN) :: Vegetation_Fraction_TL
    REAL(fp), INTENT(IN) :: Soil_Temperature_TL
    REAL(fp), INTENT(IN) :: Land_Skin_Temperature_TL
    REAL(fp), INTENT(OUT):: Emissivity_V_TL, Emissivity_H_TL
    TYPE(iVar_type)  :: iVar
    INTEGER :: IO_Status
    

    ! Local variables
    REAL(fp) :: Tsoil_TL, Tskin_TL, smc_TL
    REAL(fp) :: vlai_TL
    REAL(fp) :: Teff_TL, rh23_TL, rv23_TL 
    REAL(fp) :: emiss_TL(2), ssalb_TL(2),tau_TL(2), r23_TL(2)

    IO_Status = 0
    
    vlai_TL  = iVar%LAI*Vegetation_Fraction_TL
    Tsoil_TL = Soil_Temperature_TL
    Tskin_TL = Land_Skin_Temperature_TL
    smc_TL   = Soil_Moisture_Content_TL 

    CALL CSEM_SoilMW_Optics_TL(TSkin_TL,Tsoil_TL,smc_TL, &
         rh23_TL, rv23_TL, teff_TL, iVar%iVar_soil)
    CALL CSEM_SoilMW_roughness_Correction_TL( rh23_TL, rv23_TL, &
         r23_TL(1), r23_TL(2), iVar%iVar_rough)
    CALL CSEM_CanopyMW_Optics_TL(vLAI_TL,ssalb_TL,tau_TL,iVar%iVar_canopy)
    CALL Two_Stream_Solution_TL(ssalb_TL, tau_TL, r23_TL, Tskin_TL,Teff_TL, &
         emiss_TL, iVar)
    Emissivity_H_TL = emiss_TL(1) ; Emissivity_V_TL = emiss_TL(2) 

  END FUNCTION NESDIS_LandMW_Emiss_TL

 FUNCTION NESDIS_LandMW_Emiss_AD(       &
      & Land_Skin_Temperature_AD,       &  ! Input
      & Soil_Temperature_AD,            &  ! Input
      & Soil_Moisture_Content_AD,       &  ! Input
      & Vegetation_Fraction_AD,         &  ! Input
      & Emissivity_H_AD,                &  ! Output
      & Emissivity_V_AD, iVar)                &  ! Output
    RESULT (IO_Status)
     
    ! Arguments
    REAL(fp), INTENT(INOUT) :: Soil_Moisture_Content_AD
    REAL(fp), INTENT(INOUT) :: Vegetation_Fraction_AD
    REAL(fp), INTENT(INOUT) :: Soil_Temperature_AD
    REAL(fp), INTENT(INOUT) :: Land_Skin_Temperature_AD
    REAL(fp), INTENT(INOUT) :: Emissivity_V_AD, Emissivity_H_AD
    TYPE(iVar_type)  :: iVar

    INTEGER :: IO_Status
    
    ! Local variables
    REAL(fp) :: Tsoil_AD, Tskin_AD, smc_AD
    REAL(fp) :: vlai_AD
    REAL(fp) :: Teff_AD, rh23_AD, rv23_AD
    REAL(fp) :: ssalb_AD(2), tau_AD(2), r23_AD(2)
    REAL(fp) :: emiss_AD(2)
  

    IO_Status = 0

    emiss_AD(1) = Emissivity_H_AD ; emiss_AD(2) = Emissivity_V_AD   

    Tskin_AD = ZERO ; Tsoil_AD = ZERO ; smc_AD  = ZERO
    rh23_AD  = ZERO ; rv23_AD  = ZERO ; r23_AD  = ZERO
    ssalb_AD = ZERO ; tau_AD   = ZERO ; teff_AD = ZERO
    vLAI_AD = ZERO

    CALL Two_Stream_Solution_AD(ssalb_AD, tau_AD, r23_AD,      &
         Tskin_AD, Teff_AD, emiss_AD, iVar)
    CALL CSEM_CanopyMW_Optics_AD(vLAI_AD, ssalb_AD, tau_AD,    &
         iVar%iVar_canopy)
    CALL CSEM_SoilMW_roughness_Correction_AD(rh23_AD, rv23_AD, &
         r23_AD(1), r23_AD(2), iVar%iVar_rough)
    CALL CSEM_SoilMW_Optics_AD(TSkin_AD, Tsoil_AD, smc_AD,     &
         rh23_AD, rv23_AD, teff_AD, iVar%iVar_soil)

    Soil_Moisture_Content_AD = Soil_Moisture_Content_AD + smc_AD 
    Land_Skin_Temperature_AD = Land_Skin_Temperature_AD + Tskin_AD
    Soil_Temperature_AD      = Soil_Temperature_AD      + Tsoil_AD
    Vegetation_Fraction_AD   = Vegetation_Fraction_AD   + iVar%LAI * vlai_AD 

    Emissivity_H_AD= ZERO
    Emissivity_V_AD= ZERO
 
  END FUNCTION NESDIS_LandMW_Emiss_AD


 SUBROUTINE Two_Stream_Solution(emiss, iVar)

    REAL(fp), INTENT(OUT) :: emiss(2)
    TYPE(iVar_type) :: iVar
    ! local
    REAL(fp) :: f1, f2, f3, f4, f5
    REAL(fp) :: mu, beta
   
    INTEGER :: i
    
    mu = cos(iVar%Theta)
    iVar%iVar2%fact0   =   EXP(C_2*iVar%frequency/iVar%Tskin) - ONE
    iVar%iVar2%fact1   =   EXP(C_2*iVar%frequency/iVar%Teff)  - ONE

    DO i = 1, 2
      iVar%iVar2%alfa(i)    = SQRT((ONE - iVar%ssalb(i))/(ONE -  iVar%g(i)*iVar%ssalb(i)))
      iVar%iVar2%kk(i)      = SQRT((ONE - iVar%ssalb(i))*(ONE -  iVar%g(i)*iVar%ssalb(i)))/mu
      iVar%iVar2%beta(i)    = (ONE - iVar%iVar2%alfa(i) )/(ONE + iVar%iVar2%alfa(i) )
      iVar%iVar2%gamma(i)   = (iVar%iVar2%beta(i)  -iVar%r23(i))/(ONE-iVar%iVar2%beta(i) *iVar%r23(i))
      iVar%iVar2%fact(1,i)  = EXP(-iVar%iVar2%kk(i) *iVar%tau(i))
      iVar%iVar2%fact(2,i)  = ((ONE-(iVar%iVar2%beta(i))**2 )/(ONE-iVar%iVar2%beta(i) *iVar%r23(i)))

      f1   = iVar%iVar2%gamma(i) *  iVar%iVar2%fact(1,i)**2
      f2   = (ONE-iVar%r23(i))*(iVar%iVar2%fact0/iVar%iVar2%fact1-ONE)
      f3   =  iVar%iVar2%fact(1,i)*iVar%iVar2%fact(2,i)

      beta = iVar%iVar2%beta(i)
      f4 = (ONE - beta)*(ONE + f1)+f2*f3
      f5 =  ONE-beta*iVar%r21(i)-(beta-iVar%r21(i))*f1
     
      emiss(i)  = iVar%t21(i)* f4 / f5

      if (emiss(i) < EMISSH_DEFAULT) emiss(i) = EMISSH_DEFAULT
      if (emiss(i) > ONE) emiss(i) = ONE

    END DO


  END SUBROUTINE Two_Stream_Solution

 SUBROUTINE Two_Stream_Solution_TL(ssalb_TL, tau_TL, r23_TL, Tskin_TL,Tsoil_TL, &
                                   emiss_TL, iVar)

    REAL(fp), INTENT(IN)  :: ssalb_TL(2), tau_TL(2), r23_TL(2)
    REAL(fp), INTENT(IN)  :: Tskin_TL,Tsoil_TL
    REAL(fp), INTENT(OUT) :: emiss_TL(2)
    !local
    REAL(fp) :: mu, beta
    REAL(fp) :: Tskin, Tsoil
    REAL(fp) :: x1, x2, x3
    REAL(fp) :: f1, f2, f3, f4, f5
    REAL(fp) :: f1_TL, f2_TL, f3_TL, f4_TL, f5_TL
    REAL(fp) :: alfa_TL, beta_TL, gamma_TL, kk_TL
    REAL(fp) :: fact0_TL, fact1_TL, fact_TL(2)

    TYPE(iVar_type)      :: iVar
    TYPE(iVar_2stream)   :: iVar2
    INTEGER :: i

    iVar2 = iVar%iVar2
    
    mu = cos(iVar%Theta)
    Tskin = iVar%Tskin ; Tsoil = iVar%Teff 
    fact0_TL    = -C_2*iVar%frequency*EXP(C_2*iVar%frequency/Tskin)/Tskin**2 * Tskin_TL
    fact1_TL    = -C_2*iVar%frequency*EXP(C_2*iVar%frequency/Tsoil)/Tsoil**2 * Tsoil_TL

    DO i = 1, 2
      alfa_TL  = (ONE/(2.0_fp*iVar2%alfa(i))) * (iVar%g(i)-ONE) / &
                 (ONE - iVar%g(i)*iVar%ssalb(i))**2 * ssalb_TL(i)
      kk_TL    = (2.0_fp*iVar%ssalb(i) - iVar%g(i) - ONE) / &
                 (2.0_fp*iVar2%kk(i)*mu**2) * ssalb_TL(i)
      beta     = ivar2%beta(i)
      beta_TL  = -2.0_fp / ( ONE + iVar2%alfa(i) )**2 * alfa_TL
    
      gamma_TL = (ONE - iVar%r23(i)**2)/(ONE - beta*iVar%r23(i))**2 * beta_TL  + &
                 (beta**2 - ONE)/(ONE-beta*iVar%r23(i))**2 * r23_TL(i)
      fact_TL(1) = -EXP(-iVar2%kk(i) * iVar%tau(i)) * &
                   (iVar%tau(i) * kk_TL + iVar2%kk(i) * tau_TL(i))
 
      x1 = iVar%r23(i) * (ONE - beta**2) / (ONE - beta *iVar%r23(i))**2 -  &
           2.0_fp*beta / (ONE - beta *iVar%r23(i))
      x2 = beta * (ONE - beta**2)/(ONE - beta *iVar%r23(i))**2
      fact_TL(2)    = x1 * beta_TL +  x2 * r23_TL(i)
 
      f1 = iVar%iVar2%gamma(i) *  iVar%iVar2%fact(1,i)**2
      f2 = (ONE - iVar%r23(i)) * (iVar2%fact0/iVar2%fact1 - ONE)
      f3 = iVar2%fact(1,i)*iVar2%fact(2,i)

      f1_TL = iVar2%fact(1,i)**2 * gamma_TL + &
         TWO*iVar2%gamma(i)*iVar2%fact(1,i) * fact_TL(1) 

      x1 = -(iVar2%fact0/iVar2%fact1 - ONE)
      x2 =  (ONE - iVar%r23(i)) * ONE/iVar2%fact1
      x3 = -(ONE - iVar%r23(i)) * iVar2%fact0/(iVar2%fact1**2)
      f2_TL = x1 * r23_TL(i) + x2 * fact0_TL +x3 * fact1_TL 
      f3_TL = iVar2%fact(2,i) * fact_TL(1) + iVar2%fact(1,i) * fact_TL(2)

      f4 = (ONE - beta)*(ONE + f1) + f2 * f3
      f5 =  ONE - beta*iVar%r21(i) - (beta-iVar%r21(i))*f1
      f4_TL = (ONE - beta)*f1_TL +f3*f2_TL + f2*f3_TL - (ONE + f1) * beta_TL
      f5_TL = -(iVar%r21(i)+f1) * beta_TL - (beta-iVar%r21(i)) * f1_TL
      emiss_TL(i) = iVar%t21(i) * ONE/f5 * f4_TL - iVar%t21(i) * f4/f5**2 * f5_TL
    END DO

  END SUBROUTINE Two_Stream_Solution_TL

 SUBROUTINE Two_Stream_Solution_AD(ssalb_AD, tau_AD, r23_AD, Tskin_AD,Tsoil_AD, &
                                   emiss_AD, iVar)

    REAL(fp), INTENT(INOUT)  :: ssalb_AD(2), tau_AD(2), r23_AD(2)
    REAL(fp), INTENT(INOUT)  :: Tskin_AD, Tsoil_AD
    REAL(fp), INTENT(INOUT)  :: emiss_AD(2)
    !local
    REAL(fp) :: mu, beta
    REAL(fp) :: Tskin, Tsoil
    REAL(fp) :: x1, x2, x3
    REAL(fp) :: f1, f2, f3, f4, f5
    REAL(fp) :: f1_AD, f2_AD, f3_AD, f4_AD, f5_AD
    REAL(fp) :: alfa_AD, beta_AD, gamma_AD, kk_AD
    REAL(fp) :: fact0_AD, fact1_AD, fact_AD(2)
 
    TYPE(iVar_type)      :: iVar
    TYPE(iVar_2stream)   :: iVar2
    
    INTEGER :: i       
       
    iVar2 = iVar%iVar2
    mu = cos(iVar%Theta)
    
    Tskin = iVar%Tskin ; Tsoil = iVar%Teff

    fact0_AD = ZERO ; fact1_AD = ZERO
    DO i = 1, 2
      beta   = ivar2%beta(i)
      f1   = iVar%iVar2%gamma(i) *  iVar%iVar2%fact(1,i)**2
      f2   = (ONE-iVar%r23(i)) * (iVar2%fact0/iVar2%fact1 - ONE)
      f3   = iVar2%fact(1,i)*iVar2%fact(2,i)
      f4 = (ONE - beta)*(ONE + f1)+f2*f3
      f5 =  ONE-beta*iVar%r21(i)-(beta-iVar%r21(i))*f1

      f5_AD = - iVar%t21(i) * f4/f5**2  * emiss_AD(i)
      f4_AD =   iVar%t21(i) * ONE/f5    * emiss_AD(i)
      beta_AD = -(iVar%r21(i)+f1) * f5_AD
      f1_AD = - (beta-iVar%r21(i))* f5_AD

      beta_AD = beta_AD - (ONE + f1)  * f4_AD
      f1_AD = f1_AD + (ONE - beta)* f4_AD
      f2_AD = f3 * f4_AD
      f3_AD = f2 * f4_AD
      fact_AD(1) = iVar2%fact(2,i)*f3_AD
      fact_AD(2) = iVar2%fact(1,i)*f3_AD

      x1 = -(iVar2%fact0/iVar2%fact1-ONE)
      x2 = (ONE-iVar%r23(i))*ONE/iVar2%fact1
      x3 = - (ONE-iVar%r23(i)) * iVar2%fact0/(iVar2%fact1**2)
      fact0_AD  = fact0_AD + x2 *f2_AD
      fact1_AD  = fact1_AD + x3 *f2_AD
      r23_AD(i) = r23_AD(i)+ x1 *f2_AD
      fact_AD(1) = fact_AD(1)  + TWO*iVar2%gamma(i)*iVar2%fact(1,i) * f1_AD
      gamma_AD = iVar2%fact(1,i)**2 * f1_AD

      x1 = iVar%r23(i)*(ONE-beta**2)/(ONE-beta *iVar%r23(i))**2 - &
           2.0_fp*beta /(ONE-beta *iVar%r23(i))
      x2 = beta*(ONE-beta**2)/(ONE-beta *iVar%r23(i))**2
      beta_AD   = beta_AD   + x1 * fact_AD(2)
      r23_AD(i) = r23_AD(i) + x2 * fact_AD(2)
      
      kk_AD     =  -EXP(-iVar2%kk(i) *iVar%tau(i))*iVar%tau(i) * fact_AD(1) 
      tau_AD(i) = tau_AD(i) - EXP(-iVar2%kk(i) *iVar%tau(i))*iVar2%kk(i) * fact_AD(1) 
      
       
      beta_AD   = beta_AD   + (ONE-iVar%r23(i)**2)/(ONE-beta*iVar%r23(i))**2 * gamma_AD
      r23_AD(i) = r23_AD(i) + (beta**2-ONE)/(ONE-beta*iVar%r23(i))**2 * gamma_AD

      alfa_AD   = -2.0_fp/(ONE + iVar2%alfa(i))**2 * beta_AD

      ssalb_AD(i) = ssalb_AD(i) + &
            (2.0_fp*iVar%ssalb(i)-iVar%g(i)-ONE)/(2.0_fp*iVar2%kk(i)*mu**2) * kk_AD 
      ssalb_AD(i) = ssalb_AD(i) +&
            ONE/(2.0_fp*iVar2%alfa(i))*(iVar%g(i)-ONE) & 
                   /(ONE-iVar%g(i)*iVar%ssalb(i))**2 * alfa_AD

      emiss_AD(i) = ZERO
    END DO
    Tskin_AD  = Tskin_AD - &
        C_2*iVar%frequency*EXP(C_2*iVar%frequency/Tskin)/Tskin**2 * fact0_AD
    Tsoil_AD  = Tsoil_AD - &
        C_2*iVar%frequency*EXP(C_2*iVar%frequency/Tsoil)/Tsoil**2 * fact1_AD
    
  END SUBROUTINE Two_Stream_Solution_AD


  SUBROUTINE Load_Vege_Params(tbl, VegPar)
    CHARACTER(*)              :: tbl ! define your own lut
    TYPE(Veg_Parameters)      :: VegPar 
    REAL(fp), DIMENSION(0:13) :: veg_rho         ! Veg Specific Density
    REAL(fp), DIMENSION(0:13) :: veg_mge         ! Veg water content
    REAL(fp), DIMENSION(0:13) :: lai_min         ! Minimum LAI	  
    REAL(fp), DIMENSION(0:13) :: lai_max         ! Maximum LAI
    
    REAL(fp), DIMENSION(0:13) :: leaf_thick      ! Leaf_thickness
    REAL(fp), DIMENSION(0:13) :: leaf_angls      ! Dominant leaf 
                                                 ! inclination angle 
    REAL(fp), DIMENSION(0:13) :: sigmas          ! surface roughness length
    INTEGER :: iveg = 0

    SELECT CASE (TRIM(tbl))
      CASE ("CRTM_REL213")
         lai_min(0:13)    = (/ &
         0.52_fp, 3.08_fp, 1.85_fp, 2.80_fp, 5.00_fp, 1.00_fp, 0.50_fp, &
         0.52_fp, 0.60_fp, 0.50_fp, 0.60_fp, 0.10_fp, 1.56_fp, 0.01_fp /)
         lai_max(0:13)    = (/ &
         2.90_fp, 6.48_fp, 3.31_fp, 5.50_fp, 6.40_fp, 5.16_fp, 3.66_fp, &
         2.90_fp, 2.60_fp, 3.66_fp, 2.60_fp, 0.75_fp, 5.68_fp, 0.01_fp /)
         leaf_thick(0:13) = (/ &
         0.07_fp, 0.18_fp, 0.18_fp, 0.18_fp, 0.18_fp, 0.18_fp, 0.12_fp, &
         0.12_fp, 0.12_fp, 0.12_fp, 0.12_fp, 0.12_fp, 0.15_fp, 0.12_fp /)
         leaf_angls(0:13) = (/ &
         30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, &
         45.0_fp, 50.0_fp, 0.00_fp, 0.00_fp, 0.00_fp, 0.00_fp, 10.0_fp /)
         veg_rho(0:13)    = (/ &
         0.33_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.25_fp, &
         0.25_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.33_fp, 0.33_fp /)
         veg_mge(0:13)    = (/ &
         0.50_fp, 0.45_fp, 0.45_fp, 0.45_fp, 0.40_fp, 0.40_fp, 0.30_fp, &
         0.35_fp, 0.30_fp, 0.30_fp, 0.40_fp, 0.30_fp, 0.50_fp, 0.40_fp /)

      CASE DEFAULT
         lai_min(0:13)    = (/ &
         0.52_fp, 3.08_fp, 1.85_fp, 2.80_fp, 5.00_fp, 1.00_fp, 0.50_fp, &
         0.52_fp, 0.60_fp, 0.50_fp, 0.60_fp, 0.10_fp, 1.56_fp, 0.01_fp /)
         lai_max(0:13)    = (/ &
         2.90_fp, 6.48_fp, 3.31_fp, 5.50_fp, 6.40_fp, 5.16_fp, 3.66_fp, &
         2.90_fp, 2.60_fp, 3.66_fp, 2.60_fp, 0.75_fp, 5.68_fp, 0.01_fp /)
         leaf_thick(0:13) = (/ &
         0.07_fp, 0.20_fp, 0.18_fp, 0.15_fp, 0.12_fp, 0.12_fp, 0.18_fp, &
         0.18_fp, 0.18_fp, 0.18_fp, 0.18_fp, 0.15_fp, 0.20_fp, 0.12_fp /)
         leaf_angls(0:13) = (/ &
         30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, 30.0_fp, &
         45.0_fp, 50.0_fp, 0.00_fp, 0.00_fp, 0.00_fp, 0.00_fp, 10.0_fp /)
         veg_rho(0:13)    = (/ &
         0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.45_fp, 0.45_fp, 0.40_fp, &
         0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.33_fp /)
         veg_mge(0:13)    = (/ &
         0.50_fp, 0.50_fp, 0.50_fp, 0.50_fp, 0.50_fp, 0.50_fp, 0.50_fp, &
         0.60_fp, 0.50_fp, 0.40_fp, 0.30_fp, 0.50_fp, 0.50_fp, 0.40_fp /)
         sigmas(0:13)     = (/ &
         1.00_fp, 1.00_fp, 1.00_fp, 1.00_fp, 2.00_fp, 2.50_fp, 0.80_fp, &
         0.60_fp, 0.40_fp, 0.25_fp, 0.25_fp, 0.25_fp, 0.80_fp, 0.30_fp/) 

    END SELECT
      
    iveg = VegPar%veg_type
    
    VegPar%maxlai = lai_max(iveg)    ;  VegPar%minlai = lai_min(iveg)
    VegPar%lthick = leaf_thick(iveg) ;  VegPar%langl  = leaf_angls(iveg)
    VegPar%vrho   = veg_rho(iveg)    ;  VegPar%vmge   = veg_mge(iveg)
    VegPar%sigma  = sigmas(iveg)
     
  END SUBROUTINE Load_Vege_Params
  
  SUBROUTINE Load_Soil_Params(tbl, SoilPar)
  
    CHARACTER(*)              ::  tbl ! define your own lut
    TYPE(Soil_Parameters)     ::  SoilPar 
    REAL(fp), DIMENSION(0:9)  ::  frac_sand      ! Sand fraction
    REAL(fp), DIMENSION(0:9)  ::  frac_clay      ! clay fraction
    REAL(fp), DIMENSION(0:9)  ::  rhob_soil      ! soil bulk desnsity
    INTEGER :: isoil = 0
   
    SELECT CASE (TRIM(tbl))
      CASE ("CRTM_REL213")
         frac_sand(0:9) = (/ 0.80_fp, 0.92_fp, 0.10_fp, &
         0.20_fp, 0.51_fp, 0.50_fp, 0.35_fp, 0.60_fp, 0.42_fp, 0.92_fp /)
         frac_clay(0:9) = (/ 0.20_fp, 0.06_fp, 0.34_fp, &
         0.63_fp, 0.14_fp, 0.43_fp, 0.34_fp, 0.28_fp, 0.085_fp,0.06_fp /)
         rhob_soil(0:9) = (/ 1.48_fp, 1.68_fp, 1.27_fp, &
         1.21_fp, 1.48_fp, 1.31_fp, 1.32_fp, 1.40_fp, 1.54_fp, 1.68_fp /)

      CASE DEFAULT
         frac_sand(0:9) = (/ 0.80_fp, 0.92_fp, 0.10_fp, &
         0.20_fp, 0.51_fp, 0.50_fp, 0.35_fp, 0.60_fp, 0.42_fp, 0.92_fp /)  
         frac_clay(0:9) = (/ 0.20_fp, 0.06_fp, 0.34_fp, &
         0.63_fp, 0.14_fp, 0.43_fp, 0.34_fp, 0.28_fp, 0.085_fp,0.06_fp /)
         rhob_soil(0:9) = (/ 1.48_fp, 1.90_fp, 1.27_fp, &
         1.21_fp, 1.48_fp, 1.31_fp, 1.32_fp, 1.40_fp, 1.54_fp, 1.68_fp /)  
  
    END SELECT
   
    isoil = SoilPar%soil_type

    SoilPar%rhos = 2.65_fp
    SoilPar%rhob = rhob_soil(isoil)
    SoilPar%sand = frac_sand(isoil)
    SoilPar%clay = frac_clay(isoil)
    
  END SUBROUTINE Load_Soil_Params


  
END MODULE NESDIS_LandMW_PhyModel
