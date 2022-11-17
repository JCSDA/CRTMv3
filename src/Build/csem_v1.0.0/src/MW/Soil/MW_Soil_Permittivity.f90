!
! MW_Soil_Permittivity
!
! Module to compute the soil dielectric properties for LAND surfaces at
! microwave frequencies required for determining the LAND surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to add their soil
! optical model codes and to simplify integration into
! the CSEM_LandMW_Emiss module.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!



MODULE MW_Soil_Permittivity

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE MW_SoilWater_Permittivity

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CSEM_SoilMW_Permittivity
  PUBLIC :: CSEM_SoilMW_Permittivity_TL
  PUBLIC :: CSEM_SoilMW_Permittivity_AD
  PUBLIC :: iVar_type
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: PI     = 3.141592653589793238462643_fp 
  REAL(fp), PARAMETER :: TWOPI  = 2.0_fp*Pi
  REAL(fp), PARAMETER :: ZERO   = 0.0_fp, ONE = 1.0_fp
  ! set soil model common constants
  
  REAL(fp),    PARAMETER ::  rho_s = 2.65_fp        !   rho_s : soil specific density (g/cm3), ATBD 2.664
  REAL(fp),    PARAMETER ::  eps_winf  = 4.9_fp     !   eps_winf : dielectric constant at infinite frequency (Stogryn 1971)
                                                    !   also called: high frequency dielectric constant (4.9 Lane and Saxton )
  REAL(fp),    PARAMETER ::  eps_0 = 8.854e-12_fp   !   permittivity of free space (Klein and Swift 1977) [Farads/meter]
  COMPLEX(fp), PARAMETER ::  ea = (1.0_fp, 0.0_fp)  !   dielectric constant of air
  COMPLEX(fp), PARAMETER ::  er = (5.5_fp, 0.2_fp)  !   dielectric constant of rock
  COMPLEX(fp), PARAMETER ::  ef = (5.0_fp, 0.5_fp)  !   dielectric constant of frozen soil (Hallikainen et al. 1984, 1985)
  COMPLEX(fp), PARAMETER ::  ei = (3.2_fp, 0.1_fp)  !   dielectric constant of ice

  INTEGER, PARAMETER  :: N_Diel_Models = 5
  CHARACTER(LEN=20), PARAMETER :: Soil_Diel_Models(N_Diel_Models) = &
         [ CHARACTER(LEN=20) :: 'CRTM3.0.0', 'CRTM2.1.3', 'Dobson', 'Wang', 'Mironov']

  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Forward model input values
    REAL(fp)    :: freq            = ZERO
    REAL(fp)    :: tsoil           = ZERO
    REAL(fp)    :: smc             = ZERO
    INTEGER     :: Soil_Type       = 1
    REAL(fp)    :: alpha, beta, ess, sigma_eff
    REAL(fp)    :: eswo, tauw, esof, eswi
    COMPLEX(fp) :: esw, esm, esm1, esoil

    INTEGER,   PUBLIC :: Model_Option    = 1
    REAL(fp),  PUBLIC :: rhob, rhos=2.65_fp, salinity=0.65_fp
  END TYPE iVar_type

CONTAINS


  SUBROUTINE CSEM_SoilMW_Permittivity(frequency, tsoil, smc, sand, clay, eps, ivar)
    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: tsoil, smc
    REAL(fp), INTENT(IN)  :: sand, clay

    COMPLEX(fp),INTENT(OUT) :: eps
    TYPE(iVar_type) :: iVar
    
    iVar%rhos = rho_s 

    SELECT CASE (TRIM(Soil_Diel_Models(iVar%Model_Option)))
 
      CASE('CRTM3.0.0')
        CALL CRTM_SoilMW_Diel_V3(frequency, tsoil, smc, sand, clay, eps, iVar)
   
      CASE('CRTM2.1.3')
        CALL CRTM_SoilMW_Diel(frequency, tsoil, smc, sand, clay, iVar%rhob, eps, iVar)

      CASE('Bobson')
        CALL Dobson_SoilMW_Diel(frequency, tsoil, smc, sand, clay, iVar%salinity, eps)

      CASE('Wang')
        CALL Wang_SoilMW_Diel(frequency, tsoil, smc, sand, clay,  eps)

      CASE('Mironov')
        CALL Mironov_SoilMW_Diel(frequency, smc, clay, eps)

      CASE Default
        WRITE(*,*)"MW_Soil_Permittivity: wrong model option ...."
        RETURN
     
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Permittivity
  

  SUBROUTINE CSEM_SoilMW_Permittivity_TL(tsoil_TL,smc_TL, eps_TL, iVar)
    REAL(fp),   INTENT(IN)   :: tsoil_TL, smc_TL
    COMPLEX(fp),INTENT(OUT)  :: eps_TL
    TYPE(iVar_type) :: iVar
    
    SELECT CASE (TRIM(Soil_Diel_Models(iVar%Model_Option)))
    
      CASE('CRTM3.0.0')
        CALL CRTM_SoilMW_Diel_V3_TL(tsoil_TL,smc_TL,eps_TL, iVar)
      CASE('CRTM2.1.3')
        CALL CRTM_SoilMW_Diel_TL(tsoil_TL,smc_TL,eps_TL, iVar)
      CASE Default
        eps_TL = 0.0
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Permittivity_TL
  
  SUBROUTINE CSEM_SoilMW_Permittivity_AD( tsoil_AD, smc_AD, eps_AD, iVar)
    COMPLEX(fp), INTENT(INOUT)  :: eps_AD
    REAL(fp),    INTENT(INOUT)  :: tsoil_AD, smc_AD
    TYPE(iVar_type) :: iVar

    
    SELECT CASE (TRIM(Soil_Diel_Models(iVar%Model_Option)))
    
      CASE('CRTM3.0.0')
        CALL CRTM_SoilMW_Diel_V3_AD(tsoil_AD,smc_AD, eps_AD, iVar)
      CASE('CRTM2.1.3')
        CALL CRTM_SoilMW_Diel_AD(tsoil_AD,smc_AD, eps_AD, iVar)
      CASE Default
        RETURN
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Permittivity_AD
  
  SUBROUTINE CRTM_SoilMW_Diel(freq, t_soil, vmc, sand, clay, rhob, esm, iVar)

  !----------------------------------------------------------------------------------
  !$$$  subprogram documentation block
  !                .      .    .                                       .
  ! subprogram:    Soil_Diel   calculate the dielectric properties of soil
  !
  !   prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
  !
  ! abstract: compute the dilectric constant of the bare soil
  !
  ! program history log:
  !
  ! input argument list:
  !
  !      theta        -  local zenith angle (degree)
  !      frequency    -  frequency (ghz)
  !      t_soil       -  soil temperature
  !      vmc          -  volumetric moisture content (demensionless)
  !      rhob         -  bulk volume density of the soil (1.18-1.12)
  !      rhos         -  density of the solids (2.65 g.cm^3 for
  !                       solid soil material)
  !      sand         -  sand fraction (sand + clay = 1.0)
  !      clay         -  clay fraction
  !
  ! output argument list:
  !
  !      esm          -  dielectric constant for bare soil
  !
  ! important internal variables:
  !
  !      esof         -  the permittivity of free space
  !      eswo         -  static dieletric constant
  !      tauw         -  relaxation time of water
  !      s            -  salinity
  !
  ! remarks:
  !
  ! attributes:
  !   language: f90
  !   machine:  ibm rs/6000 sp
  !
  !----------------------------------------------------------------------------------

    REAL(fp)    :: freq, t_soil, vmc
    REAL(fp)    :: rhob, sand, clay
    COMPLEX(fp) :: esm
    COMPLEX(fp) :: es1, es2
    REAL(fp)    :: t, es1_img 

    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)
    TYPE(iVar_type) :: iVar
    
    iVar%freq  = freq*1.0e9_fp  
    iVar%tsoil = t_soil - 273.0_fp
    iVar%smc   = vmc
    iVar%alpha = 0.65_fp
    iVar%beta  = 1.09_fp - 0.11_fp*sand + 0.18_fp*clay
    iVar%ess   = (1.01_fp + 0.44_fp*iVar%rhos)**2 - 0.062_fp

    iVar%sigma_eff = -1.645_fp + 1.939_fp*rhob - 0.020213_fp*sand + 0.01594_fp*clay
    ! the permittivity at the high frequency limit
    iVar%eswi = 5.5_fp
    ! the permittivity of free space (esof)
    iVar%esof = 8.854e-12_fp
    
    ! static dieletric constant (eswo)
    t = iVar%tsoil
    iVar%eswo = a(1)+a(2)*t+a(3)*t**2+a(4)*t**3
    iVar%tauw = b(1)+b(2)*t+b(3)*t**2+b(4)*t**3
   
    IF (iVar%smc > ZERO) THEN
      es1_img = -iVar%sigma_eff*(iVar%rhos-iVar%rhob)/ &
         (twopi*iVar%freq*iVar%esof*iVar%rhos*iVar%smc)
      es1 = CMPLX(iVar%eswi, es1_img, fp)
    ELSE
      es1 = CMPLX(iVar%eswi, ZERO, fp)
    ENDIF

    es2 = CMPLX(iVar%eswo-iVar%eswi, ZERO, fp)/CMPLX(ONE, iVar%freq*iVar%tauw, fp)
    iVar%esw = es1 + es2
    iVar%esm1 = ONE + (iVar%ess**(iVar%alpha) - ONE)*iVar%rhob/iVar%rhos + &
        iVar%smc**(iVar%beta)*iVar%esw**(iVar%alpha) -iVar%smc
   
    esm = iVar%esm1**(ONE/iVar%alpha)
    iVar%esm  = esm
    IF (AIMAG(iVar%esm) >= ZERO) esm = CMPLX(REAL(iVar%esm,fp),-0.0001_fp, fp)
    iVar%esoil = esm
  END SUBROUTINE CRTM_SoilMW_Diel

  SUBROUTINE CRTM_SoilMW_Diel_TL(t_soil_tl, vmc_tl, esm_tl, iVar)

    REAL(fp),    INTENT(IN)    :: t_soil_tl, vmc_tl
    COMPLEX(fp), INTENT(OUT)   :: esm_tl
    TYPE(iVar_type),INTENT(IN) :: iVar

    COMPLEX(fp) :: ac, bc
    COMPLEX(fp) :: ac_tl, bc_tl
    COMPLEX(fp) :: es1_tl, es2_tl
    COMPLEX(fp) :: esw_tl, esm1_tl
    REAL(fp)    :: eswo_tl, tauw_tl, es1_img_tl
    REAL(fp)    :: t, t_tl

    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)
 
    t = iVar%tsoil
    t_tl = t_soil_tl
    ! static dieletric constant (eswo)
    eswo_tl = (a(2)+2.0_fp*a(3)*t+3.0_fp*a(4)*t**2)*t_tl
    tauw_tl = (b(2)+2.0_fp*b(3)*t+3.0_fp*b(4)*t**2)*t_tl

    IF (iVar%smc > ZERO) THEN
      es1_img_tl =  iVar%sigma_eff*(iVar%rhos-iVar%rhob)/ &
        (twopi*iVar%freq*iVar%esof*iVar%rhos*iVar%smc**2)*vmc_tl
      es1_tl = CMPLX(ZERO, es1_img_tl, fp)
    ELSE
      es1_tl = CMPLX(ZERO, ZERO, fp)
   ENDIF

    ac = CMPLX(iVar%eswo-iVar%eswi, ZERO, fp)
    bc = CMPLX(ONE, iVar%freq*iVar%tauw, fp)
    ac_tl = CMPLX(eswo_tl, ZERO, fp)
    bc_tl = CMPLX(ZERO,  iVar%freq*tauw_tl, fp)
    es2_tl = ONE/bc*ac_tl - ac/bc**2*bc_tl
    esw_tl = es1_tl + es2_tl

    esm1_tl = (iVar%beta* iVar%smc**(iVar%beta-ONE)*iVar%esw**iVar%alpha - ONE)* vmc_tl+ &
             (iVar%smc**iVar%beta*iVar%alpha*iVar%esw**(iVar%alpha-ONE))* esw_tl  

    esm_tl = (ONE/iVar%alpha)*iVar%esm1**(ONE/iVar%alpha-ONE)*esm1_tl

    IF (AIMAG(iVar%esm) >= ZERO) esm_tl = CMPLX(REAL(esm_tl,fp), ZERO, fp)

  END SUBROUTINE CRTM_SoilMW_Diel_TL

  SUBROUTINE CRTM_SoilMW_Diel_AD(t_soil_ad, vmc_ad, esm_ad, iVar)

    REAL(fp),    INTENT(INOUT)  :: t_soil_ad, vmc_ad
    COMPLEX(fp), INTENT(INOUT)  :: esm_ad
    TYPE(iVar_type), INTENT(INOUT) :: iVar

    COMPLEX(fp) :: ac, bc, ac_ad,  bc_ad
    COMPLEX(fp) :: es1_ad, es2_ad, esw_ad
    COMPLEX(fp) :: esm1_ad
    REAL(fp)    :: eswo_ad, tauw_ad 
    REAL(fp)    :: es1_img_ad
    REAL(fp)    :: t, t_ad

    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)
  
    IF (AIMAG(iVar%esm) >= ZERO) esm_ad = CMPLX(REAL(esm_ad,fp),ZERO, fp)
   
    esm1_ad = CONJG(ONE/iVar%alpha*iVar%esm1**(ONE/iVar%alpha-ONE))*esm_ad
    esw_ad = CONJG(iVar%smc**iVar%beta*iVar%alpha*iVar%esw**(iVar%alpha-ONE))*esm1_ad
    vmc_ad = vmc_ad + REAL(&
             CONJG(iVar%beta* iVar%smc**(iVar%beta-ONE)*iVar%esw**iVar%alpha - ONE)*esm1_ad,fp)
    esm_ad = CMPLX(ZERO, ZERO,fp)
    
    es1_ad = esw_ad
    es2_ad = esw_ad
    esw_ad = CMPLX(ZERO, ZERO,fp)
    
    ac = CMPLX(iVar%eswo-iVar%eswi, ZERO, fp)
    bc = CMPLX(ONE, iVar%freq*iVar%tauw, fp)
    ac_ad =  CONJG(ONE/bc)*es2_ad
    bc_ad =  - CONJG(ac/bc**2)*es2_ad
    eswo_ad = REAL(ac_ad)
    tauw_ad =  iVar%freq*AIMAG(bc_ad)
    
    IF (iVar%smc > ZERO) THEN
      es1_img_ad = AIMAG( es1_ad)
      vmc_ad = vmc_ad + iVar%sigma_eff*(iVar%rhos-iVar%rhob) / &
        (twopi*iVar%freq*iVar%esof*iVar%rhos*iVar%smc**2)*es1_img_ad
    ENDIF
    es1_img_ad = ZERO
    t = iVar%tsoil

    t_ad = ZERO
    t_ad= t_ad + (b(2)+2.0_fp*b(3)*t+3.0_fp*b(4)*t**2)*tauw_ad
    t_ad= t_ad + (a(2)+2.0_fp*a(3)*t+3.0_fp*a(4)*t**2)*eswo_ad
    t_soil_ad = t_soil_ad + t_ad
    esm_ad = ZERO
  END SUBROUTINE CRTM_SoilMW_Diel_AD

  SUBROUTINE CRTM_SoilMW_Diel_V3(freq, t_soil, vmc, sand, clay, esm, iVar)

    REAL(fp),    INTENT(IN)    :: freq, t_soil, vmc, sand, clay 
    COMPLEX(fp), INTENT(OUT)   :: esm
    TYPE(iVar_type)  :: iVar
    
    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)
   
    
    REAL(fp) :: t
    REAL(fp) :: omega,ftau2, feps
    REAL(fp) :: esw_real, esw_img
    COMPLEX(fp) :: epsa

    iVar%freq  = freq*1.0e9_fp  
    iVar%tsoil = t_soil - 273.0_fp
    iVar%smc   = vmc
    iVar%alpha = 0.60_fp
    iVar%beta  = 1.09_fp - 0.11_fp*sand + 0.18_fp*clay
    
    iVar%rhos = rho_s
    iVar%rhob = sand*1.6_fp +  clay*1.1_fp + (1._fp-sand-clay)*1.2_fp
    iVar%ess = (1.01_fp + 0.44_fp*iVar%rhos)**2 - 0.062_fp
    iVar%sigma_eff = -1.645_fp + 1.939_fp*iVar%rhob - 2.0213_fp*sand + 1.594_fp*clay

    IF(iVar%smc <= 1.0e-5) THEN
      iVar%esm = (ONE + (iVar%ess**iVar%alpha-ONE)*iVar%rhob/iVar%rhos)**(ONE/iVar%alpha)
      RETURN
    ENDIF
    
    ! the permittivity at the high frequency limit
    iVar%eswi = 4.9_fp
    ! the permittivity of free space (esof)
    iVar%esof = 8.854e-12_fp

    ! static dieletric constant (eswo)
    t = iVar%tsoil
    iVar%eswo = a(1)+a(2)*t+a(3)*t**2+a(4)*t**3
    iVar%tauw = b(1)+b(2)*t+b(3)*t**2+b(4)*t**3

   
    omega =  1.2*iVar%freq
    ftau2 = one+(omega*iVar%tauw)**2
    feps  =  omega*iVar%tauw*(iVar%eswo- iVar%eswi)
    esw_real = iVar%eswi + (iVar%eswo-iVar%eswi)/ftau2
    esw_img  = omega*iVar%tauw*(iVar%eswo- iVar%eswi)/ftau2 + &
          iVar%sigma_eff/(TWOPI*iVar%freq*iVar%esof)*(iVar%rhos-iVar%rhob)/(iVar%rhos*iVar%smc)
    
    iVar%esw   = CMPLX(esw_real, -esw_img,fp)
    epsa  = ONE + (iVar%ess**iVar%alpha-ONE)*iVar%rhob/iVar%rhos + &
        (iVar%smc**iVar%beta)*iVar%esw**iVar%alpha-iVar%smc
    esm   = epsa**(ONE/iVar%alpha) 
    iVar%esm   = esm
    
    IF (AIMAG(iVar%esm) >= ZERO) esm = CMPLX(REAL(iVar%esm,fp),-0.0001_fp, fp)
    iVar%esoil = esm

  END SUBROUTINE CRTM_SoilMW_Diel_V3 

  SUBROUTINE CRTM_SoilMW_Diel_V3_TL(t_soil_TL,vmc_TL,esm_TL, iVar)

    REAL(fp),    INTENT(IN)   :: t_soil_TL,vmc_TL
    COMPLEX(fp), INTENT(OUT)  :: esm_TL
    TYPE(iVar_type)  :: iVar
      
    REAL(fp)   :: t
    REAL(fp)   :: omega, ftau2, feps
    COMPLEX(fp):: epsa
    REAL(fp)   :: t_TL, tauw_TL, ftau2_TL, feps_TL
    REAL(fp)   :: eswo_TL, esw_real_TL, esw_img_TL
    COMPLEX(fp):: esw_TL, epsa_TL

    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)

    IF(iVar%smc <= 1.0e-5) THEN
      esm_TL = 0.0_fp
      RETURN
    ENDIF
    
    omega = 1.2_fp*iVar%freq
    ftau2 = 1.0_fp+(omega*iVar%tauw)**2
    feps  = omega*iVar%tauw*(iVar%eswo - iVar%eswi)
    epsa  = 1.0_fp+iVar%rhob/iVar%rhos*(iVar%ess**iVar%alpha-1.0_fp)+ &
            iVar%smc**iVar%beta*iVar%esw**iVar%alpha-iVar%smc

    t = iVar%tsoil ; t_TL = t_soil_TL
    eswo_TL  = (a(2)+2.0_fp*a(3)*t+3.0_fp*a(4)*t**2)*t_TL
    tauw_TL  = (b(2)+2.0_fp*b(3)*t+3.0_fp*b(4)*t**2)*t_TL
    ftau2_TL = 2.0_fp*omega**2*iVar%tauw*tauw_TL
    feps_TL  = omega*(iVar%eswo- iVar%eswi)*tauw_TL + omega*iVar%tauw*eswo_TL
        
    esw_real_TL = 1.0_fp/ftau2*eswo_TL - (iVar%eswo-iVar%eswi)/ftau2**2*ftau2_TL
    esw_img_TL  = 1.0_fp/ftau2 * feps_TL - feps/ftau2**2 * ftau2_TL - &
                  iVar%sigma_eff/(TWOPI*iVar%freq*iVar%esof)* &
                  (iVar%rhos-iVar%rhob)/(iVar%rhos*iVar%smc**2)*vmc_TL

    esw_TL  = CMPLX(esw_real_TL, -esw_img_TL, fp)
    epsa_TL = (iVar%beta*iVar%smc**(iVar%beta-1.0_fp)*iVar%esw**iVar%alpha-1.0_fp)*vmc_TL + &
               iVar%alpha*iVar%smc**iVar%beta*iVar%esw**(iVar%alpha-1)*esw_TL
    
    esm_TL  = ONE/iVar%alpha*epsa**(ONE/iVar%alpha-ONE)*epsa_TL
    
    IF (AIMAG(iVar%esm) >= ZERO) esm_TL = CMPLX(REAL(esm_TL,fp), ZERO, fp)

  END SUBROUTINE CRTM_SoilMW_Diel_V3_TL
  
 
  SUBROUTINE CRTM_SoilMW_Diel_V3_AD(t_soil_AD, vmc_AD, esm_AD, iVar)

    COMPLEX(fp), INTENT(INOUT):: esm_AD
    REAL(fp),    INTENT(INOUT):: t_soil_AD, vmc_AD
    TYPE(iVar_type)  :: iVar
       
    REAL(fp)   :: t
    REAL(fp)   :: omega,  ftau2, feps
    COMPLEX(fp):: epsa
    REAL(fp)   :: t_AD, tauw_AD, ftau2_AD, feps_AD
    REAL(fp)   :: eswo_AD, esw_real_AD, esw_img_AD
    COMPLEX(fp):: esw_AD, epsa_AD

    REAL(fp) :: a(4) = (/ 87.134_fp, -1.949e-1_fp,-1.276e-2_fp, 2.491e-4_fp/)
    REAL(fp) :: b(4) = (/ 1.1109e-10_fp, -3.824e-12_fp, 6.938e-14_fp, -5.096e-16_fp/)

 
    IF(iVar%smc <= 1.0e-5) THEN
      t_soil_AD = 0.0
      vmc_AD = 0.0
      esm_AD = 0.0
      RETURN
    ENDIF

    omega = 1.2_fp*iVar%freq
    ftau2 = 1.0_fp+(omega*iVar%tauw)**2
    feps  = omega*iVar%tauw*(iVar%eswo - iVar%eswi)
    epsa  = 1.0_fp+iVar%rhob/iVar%rhos*(iVar%ess**iVar%alpha-1.0_fp)+ &
            iVar%smc**iVar%beta*iVar%esw**iVar%alpha-iVar%smc

    esw_AD  = zero; esw_img_AD = zero; esw_real_AD  = zero; eswo_AD = zero
    feps_AD = zero; epsa_AD = zero; ftau2_AD = zero; t_AD  = zero; tauw_AD = zero
    
    IF (AIMAG(iVar%esm) >= ZERO) esm_AD = CMPLX(REAL(esm_AD,fp), ZERO, fp)

    epsa_AD = epsa_AD + CONJG(ONE/iVar%alpha*epsa**(ONE/iVar%alpha-ONE))*esm_AD
    esm_AD  = zero
    
    vmc_AD  = vmc_AD + REAL(CONJG(iVar%beta*iVar%smc**(iVar%beta-ONE)* &
                                 iVar%esw**iVar%alpha-ONE)* epsa_AD, fp)
    esw_AD  = esw_AD + CONJG(iVar%alpha*iVar%smc**iVar%beta* &
                            iVar%esw**(iVar%alpha-ONE)) * epsa_AD
    esw_real_AD  = esw_real_AD + REAL(esw_AD)
    esw_img_AD   = esw_img_AD - AIMAG(esw_AD)
        
    vmc_AD   = vmc_AD - iVar%sigma_eff/(TWOPI*iVar%freq*iVar%esof)* &
               (iVar%rhos-iVar%rhob)/(iVar%rhos*iVar%smc**2)*esw_img_AD
    ftau2_AD = ftau2_AD - feps/ftau2**2*esw_img_AD
    feps_AD  = feps_AD + ONE/ftau2*esw_img_AD
    
    ftau2_AD = ftau2_AD - (iVar%eswo-iVar%eswi)/ftau2**2*esw_real_AD 
    eswo_AD  = eswo_AD + ONE/ftau2*esw_real_AD 
    
    eswo_AD  = eswo_AD + omega*iVar%tauw*feps_AD 
    tauw_AD  = tauw_AD + omega*(iVar%eswo- iVar%eswi)*feps_AD 
    tauw_AD  = tauw_AD + 2.0_fp*omega**2*iVar%tauw*ftau2_AD
    
    t = iVar%tsoil
    t_AD = t_AD +(b(2)+2.0_fp*b(3)*t+3.0_fp*b(4)*t**2)*tauw_AD
    t_AD = t_AD +(a(2)+2.0_fp*a(3)*t+3.0_fp*a(4)*t**2)*eswo_AD
    t_soil_AD = t_soil_AD + t_AD
    esm_ad = ZERO
  END SUBROUTINE CRTM_SoilMW_Diel_V3_AD

!===========================================================================

  SUBROUTINE Dobson_SoilMW_Diel(freq, T, wc, sand, clay, sal, eps )

! Purpose :
!   Calculate the dielectric constant of a wet soil
!   Developed and validated for 1.4 and 18 GHz.

! Reference:
!  Dobson et al., 1985: Microwave Dielectric behavior of
!    wet soil - part II: Dielectric mixing models,
!    IEEE Trans. Geosc. Rem. Sens., GE-23, No. 1, 35-46.
!    sal_soil : soil water salinity (psu = ppt(weight) ) , ATBD 0, (~0.65)
!    sal_sea : sea water salinity (psu = ppt(weight) ) , lmeb 32.5, range 6-60psu
! alphas : constant shape factor
!---------------------------------------------------------------------------
    REAL(fp),    INTENT(IN)  :: freq
    REAL(fp),    INTENT(IN)  :: T, wc, sand, clay
    COMPLEX(fp), INTENT(OUT) :: eps

    REAL(fp), PARAMETER :: alphas = 0.65_fp
    REAL(fp) :: eps_s, beta, epsr, epsi, eaa
    REAL(fp) :: wcm, rho_b, sal
    COMPLEX(fp) :: ew

    wcm = MAX(wc,0.001_fp) ! to avoid dividing by zero
    ! bulk density frho_b= 1.4
    rho_b= sand*1.6_fp +  clay*1.1_fp + (1._fp-sand-clay)*1.2_fp
 
    eps_s = (1.01_fp + 0.44_fp * rho_s)**2._fp - 0.062_fp  ! dry soil dielec
    CALL MW_Soil_Water_Diel(freq, T, wc, sal, sand, clay, ew)

    beta = 1.2748_fp - 0.519_fp * sand - 0.152_fp * clay ! Eq.(30)
    eaa = 1.0_fp + (rho_b / rho_s) * (eps_s ** alphas - 1.0_fp)   &
         + (wcm ** beta) * (real(ew)) ** alphas - wcm    ! Eq. 28
    epsr = eaa ** (1._fp/alphas)

    beta = 1.33797_fp - 0.603_fp * sand - 0.166_fp * clay ! Eq.(31)
    eaa= (wcm ** beta) * (abs(aimag(ew)) ** alphas)        ! Eq. 28
    epsi = eaa ** (1._fp/alphas)

    eps = cmplx(epsr,epsi,fp)


  END SUBROUTINE Dobson_SoilMW_Diel


  SUBROUTINE Wang_SoilMW_Diel(freq, T,  wc,  sand, clay,  eps )

! Purpose :
!   Calculate the dielectric constant of a wet soil
!   Developed and validated for 1.4 and 5 GHz.

! Reference:
!  Wang and Schmugge, 1980: An empirical model for the
!    complex dielectric permittivity of soils as a function of water
!    content. IEEE Trans. Geosci. Rem. Sens., GE-18, No. 4, 288-295.

! External :
!   dielsal_sub.F90

! Input/Output Variables:
!   wc :  volumetric soil water content of layer (cm3/cm3)
!   wp : wilting point (cm3/cm3)
!   p  : porosity of dry soil (cm3/cm3)
!   wt : transition moisture point (cm3/cm3)
!   gamma : fitting parameter
!   ecl : conductivity loss
!   ex : dielectric constant of the initially absorbed water
!   eps : Dielectric constant of the soil layer
!---------------------------------------------------------------------------

    REAL(fp),    INTENT(IN)  :: freq
    REAL(fp),    INTENT(IN)  :: wc, T, sand, clay
    COMPLEX(fp), INTENT(OUT) :: eps

    REAL(fp) :: rho_b
    REAL(fp) :: p, wp, wt
    REAL(fp) :: alpha,gamma
    REAL(fp) :: ecl
    COMPLEX(fp)  :: ex,ew
    COMPLEX(fp)  :: j = (0._fp,1._fp)
     
    ! bulk density frho_b= 1.4
    rho_b= sand*1.6_fp +  clay*1.1_fp + (1._fp-sand-clay)*1.2_fp

    ! porosity refers to the volume of soil voids that can be filled by water and/or air
    ! rho_b/rho_s = soild fraction   porosity = 1-soild fraction CEME was wrong
    p = 1.0_fp-rho_b/rho_s  !(Eq.7)

    ! Wilting point
    wp = 0.06774_fp - 0.064_fp * sand + 0.478_fp * clay !Eq.1

    ! alpha, for conductivity loss (CMEM)
    ! this parameter may be improved with 60.*wavelength*sigma
    ! sigma_eff of Dobson or ION_CONDUCT (T,sal,sigma)
    alpha = 0.0_fp
    IF ( freq <= 2.5_fp ) THEN
      alpha=100._fp * wp
      alpha=min(alpha,26.0_fp)
    ENDIF
    ! Wang's paper didn't mention free water type, but
    ! ion correction was added, so the water may be pure water
    CALL MW_Pure_Water_Diel(freq, T, ew)
    !CALL MW_Soil_Water_Diel(freq, T, wc, sal, sand, clay, ew)

    ! Calculate dielectric constant of soil-water mixture
    gamma = -0.57_fp * wp + 0.481_fp !(Eq.8)
    wt = 0.49_fp * wp + 0.165_fp     !(Eq.9)
    IF (wc <= wt) THEN
      ex = ei + (ew-ei) * (wc/wt) * gamma
      eps = wc*ex + (p-wc) * ea + (1.-p) * er                  !Eq 3
    ELSE
      ex = ei + (ew-ei) * gamma
      eps = wt*ex + (wc-wt) * ew + (p-wc) * ea + (1.-p) * er  !Eq 4
    ENDIF

    ! add conductivity loss

    ecl = alpha * wc**2._fp      !Eq.6
    eps = eps + j * ecl

  END SUBROUTINE Wang_SoilMW_Diel


 
  SUBROUTINE Mironov_SoilMW_Diel(freq, wc, clay, eps)

! Purpose :
!   Calculate the dielectric constant of a wet soil
!   Developed and validated from 1 to 10 GHz.
!   adapted for a large range of soil moisture

! Reference:
!   Mironov et al: Generalized Refractive Mixing Dielectric Model for moist soil
!    IEEE Trans. Geosc. Rem. Sens., vol 42 (4), 773-785. 2004.
!
! Adapted from the Matlab version of JP Wigneron
!
! Patricia de Rosnay, October 9 2007

    REAL(fp),    INTENT(IN)  :: freq
    REAL(fp),    INTENT(IN)  :: wc, clay
    COMPLEX(fp), INTENT(OUT) :: eps

    REAL(fp) :: f
    REAL(fp) :: mv,mvt
    REAL(fp) :: nd,kd,nb,kb,nu,ku,nm,km
    REAL(fp) :: eps0b,taub,sigmab
    REAL(fp) :: eps0u,tauu,sigmau
    REAL(fp) :: zb,epsb1,epsb2
    REAL(fp) :: zu,epsu1,epsu2
    REAL(FP) :: epsm1,epsm2
    REAL(FP) :: zflag


    f = freq *1.0E9_fp

    ! RI & NAC of dry soils
    nd = 1.634_fp - 0.539_fp * clay + 0.2748_fp * clay**2._fp
    kd = 0.03952_fp - 0.04038_fp * clay

    ! Maximum bound water fraction
    mvt = 0.02863_fp + 0.30673_fp * clay

    ! Bound water parameters
    eps0b = 79.8_fp - 85.4_fp*clay  + 32.7_fp *clay**2._fp
    taub = 1.062e-11_fp + 3.450e-12_fp * clay
    sigmab = 0.3112_fp + 0.467_fp * clay 

    ! Unbound (free) water parameters
    eps0u = 100._fp
    tauu = 8.5e-12_fp
    sigmau = 0.3631_fp + 1.217_fp * clay


    ! Computation of epsilon water (bound & unbound)
    zb = (eps0b - eps_winf) / (1._fp + (2._fp*pi*f*taub)**2._fp)
    epsb1 = eps_winf + zb
    epsb2 =  zb * (2._fp*pi*f*taub) + sigmab / (2._fp*pi*eps_0*f)

    zu = (eps0u - eps_winf) / (1._fp + (2._fp*pi*f*tauu)**2._fp)
    epsu1 = eps_winf + zu
    epsu2 =  zu * (2._fp*pi*f*tauu) + sigmau/(2._fp*pi*eps_0*f)

    ! Computation of refractive index of water (bound & unbound)
    nb= sqrt( sqrt( epsb1**2._fp + epsb2**2._fp) + epsb1 ) / sqrt(2._fp)
    kb= sqrt( sqrt( epsb1**2._fp + epsb2**2._fp) - epsb1 ) / sqrt(2._fp)

    nu= sqrt( sqrt( epsu1**2._fp + epsu2**2._fp) + epsu1 ) / sqrt(2._fp)
    ku= sqrt( sqrt( epsu1**2._fp + epsu2**2._fp) - epsu1 ) / sqrt(2._fp)

    ! Computation of soil refractive index (nm & km):
    ! xmv can be a vector
    mv= min (wc, mvt)
    zflag = 0._fp
    IF ( wc >= mvt ) zflag = 1._fp

    nm = nd + (nb - 1._fp) * mv + (nu - 1._fp) * (wc-mvt) * zflag
    km = kd + kb * mv + ku * (wc-mvt) * zflag

    ! computation of soil dielectric constant:
    epsm1= nm ** 2._fp - km ** 2._fp
    epsm2= nm * km * 2._fp

    eps = cmplx(epsm1,epsm2,fp)
    
  END SUBROUTINE Mironov_SoilMW_Diel

 
END MODULE MW_Soil_Permittivity



