!
! SoilMW_Roughness_Correction
!
! Module to compute the soil surface roughness correction for LAND surfaces at
! microwave frequencies required for determining the LAND surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to add their 
! model codes and to simplify integration into
! the CSEM_LandMW_Emiss module.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!



MODULE SoilMW_Roughness_Correction

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CSEM_SoilMW_Roughness_Correction
  PUBLIC :: CSEM_SoilMW_Roughness_Correction_TL
  PUBLIC :: CSEM_SoilMW_Roughness_Correction_AD
  PUBLIC :: iVar_type
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER  :: PI   = 3.141592653589793238462643_fp 
  REAL(fp), PARAMETER  :: ONE  = 1.0_fp
  REAL(fp), PARAMETER  :: TWO  = 2.0_fp
  REAL(fp), PARAMETER  :: ZERO = 0.0_fp
   ! set soil model common constants
  
  REAL(fp),    PARAMETER ::  rho_s = 2.66_fp        !   rho_s : soil specific density (g/cm3), ATBD 2.664
  REAL(fp),    PARAMETER ::  eps_winf  = 4.9_fp     !   eps_winf : dielectric constant at infinite frequency (Stogryn 1971)
                                                    !   also called: high frequency dielectric constant (4.9 Lane and Saxton )
  REAL(fp),    PARAMETER ::  eps_0 = 8.854e-12_fp   !   permittivity of free space (Klein and Swift 1977) [Farads/meter]
  COMPLEX(fp), PARAMETER ::  ea = (1.0_fp, 0.0_fp)  !   dielectric constant of air
  COMPLEX(fp), PARAMETER ::  er = (5.5_fp, 0.2_fp)  !   dielectric constant of rock
  COMPLEX(fp), PARAMETER ::  ef = (5.0_fp, 0.5_fp)  !   dielectric constant of frozen soil (Hallikainen et al. 1984, 1985)
  COMPLEX(fp), PARAMETER ::  ei = (3.2_fp, 0.1_fp)  !   dielectric constant of ice


  INTEGER, PARAMETER  :: N_Roughness_Models = 7
  CHARACTER(LEN=20), PARAMETER :: Soil_Roughness_Models(N_Roughness_Models) = &
         [CHARACTER(LEN=20) :: 'Chen', 'CRTM', 'Wang','Coppo','Wigneron','Wigneron_2011','Wegmuller']

  TYPE iVar_type
    PRIVATE
    ! Forward model input values
    REAL(fp) :: ph, pv, q

    INTEGER, PUBLIC  :: Model_Option = 1
   END TYPE iVar_type

CONTAINS

  SUBROUTINE CSEM_SoilMW_Roughness_Correction(frequency, sigma, theta, &
             refl_smooth, refl_rough, iVar)
    REAL(fp), INTENT(IN)    :: frequency
    REAL(fp), INTENT(IN)    :: sigma,theta
    REAL(fp), INTENT(IN)    :: refl_smooth(2)
    REAL(fp), INTENT(OUT)   :: refl_rough(2)
    TYPE(iVar_type)         :: iVar
    
      
    SELECT CASE (TRIM(Soil_Roughness_Models(iVar%Model_Option)))
      CASE('CRTM')
        CALL Roughness_PQ_CRTM(frequency,sigma,theta,iVar%ph,iVar%pv,iVar%q )
        CALL Roughness_Correction(iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
      CASE('Wang')
        CALL Roughness_PQ_Wang(frequency,sigma,theta,iVar%ph,iVar%pv,iVar%q)
        CALL Roughness_Correction(iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
      CASE('Coppo')
        CALL Roughness_PQ_Coppo(frequency,sigma, iVar%ph, iVar%pv, iVar%q)
        CALL Roughness_Correction( iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
      CASE('Wigneron')
        CALL Roughness_PQ_Wigneron(frequency,sigma, iVar%ph,iVar% pv, iVar%q)
        CALL Roughness_Correction(iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
      CASE('Wigneron_2011')
        CALL Roughness_PQ_Wigneron_2011(frequency,sigma,theta,iVar%ph,iVar%pv,iVar%q)
        CALL Roughness_Correction(iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
      CASE('Wegmuller')
        CALL Roughness_Wegmuller(frequency,sigma,theta, &
                                 refl_smooth(1), refl_rough(1), refl_rough(2))
      CASE Default
        CALL Roughness_PQ_Chen(frequency,sigma,theta,iVar%ph,iVar%pv,iVar%q )
        CALL Roughness_Correction(iVar%ph, iVar%pv, iVar%q, refl_smooth, refl_rough)
  
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Roughness_Correction
  
  SUBROUTINE CSEM_SoilMW_Roughness_Correction_TL(rh_s_TL, rv_s_TL, rh_TL, rv_TL, iVar)
    REAL(fp), INTENT(IN)    :: rh_s_TL,rv_s_TL
    REAL(fp), INTENT(OUT)   :: rh_TL, rv_TL
    TYPE(iVar_type)         :: iVar
    
    SELECT CASE (TRIM(Soil_Roughness_Models(iVar%Model_Option)))
      CASE('CRTM')
        CALL Roughness_Correction_TL(rv_s_TL, rh_s_TL, rv_TL, rh_TL, iVar)

      CASE Default
        CALL Roughness_Correction_TL(rv_s_TL, rh_s_TL, rv_TL,rh_TL, iVar)
  
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Roughness_Correction_TL
  
  SUBROUTINE CSEM_SoilMW_Roughness_Correction_AD(rh_s_AD,rv_s_AD,rh_AD,rv_AD, iVar )

    REAL(fp), INTENT(INOUT) :: rh_s_AD,rv_s_AD
    REAL(fp), INTENT(INOUT) :: rh_AD, rv_AD

    TYPE(iVar_type)         :: iVar
     
    SELECT CASE (TRIM(Soil_Roughness_Models(iVar%Model_Option)))
      CASE('CRTM')
        CALL Roughness_Correction_AD(rv_s_AD, rh_s_AD, rv_AD,rh_AD, iVar)

      CASE Default
        CALL Roughness_Correction_AD(rv_s_AD, rh_s_AD, rv_AD,rh_AD, iVar)
  
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Roughness_Correction_AD
  
  SUBROUTINE  Roughness_Correction( ph, pv, q, refl_smooth, refl_rough)
    REAL(fp), INTENT(IN)  :: ph, pv, q
    REAL(fp), INTENT(IN)  :: refl_smooth(2)
    REAL(fp), INTENT(OUT) :: refl_rough(2)
    REAL(fp) :: rh_p, rv_p
    
    rh_p = ph*refl_smooth(1)
    rv_p = pv*refl_smooth(2)
    refl_rough(1) = rh_p + q*(rv_p-rh_p)
    refl_rough(2) = rv_p + q*(rh_p-rv_p)
    
  END SUBROUTINE Roughness_Correction

  SUBROUTINE  Roughness_Correction_TL( rv_s_TL, rh_s_TL, rv_TL, rh_TL, iVar)
    REAL(fp), INTENT(IN)  :: rh_s_TL, rv_s_TL
    REAL(fp), INTENT(OUT) :: rh_TL, rv_TL
    TYPE(iVar_type)       :: iVar
    
    REAL(fp) :: rh_p, rv_p
   
    rh_p  = iVar%ph*rh_s_TL; rv_p = iVar%pv*rv_s_TL
    rh_TL = rh_p + iVar%q*(rv_p-rh_p)
    rv_TL = rv_p + iVar%q*(rh_p-rv_p)
    
  END SUBROUTINE Roughness_Correction_TL

  SUBROUTINE  Roughness_Correction_AD( rv_s_AD, rh_s_AD, rv_AD, rh_AD, iVar)
    REAL(fp), INTENT(INOUT) :: rh_AD, rv_AD
    REAL(fp), INTENT(INOUT) :: rh_s_AD, rv_s_AD
    TYPE(iVar_type)         :: iVar
    
    REAL(fp) :: rh_p, rv_p
   
    rh_p = zero; rv_p = zero
    
    rh_p = rh_p + iVar%q * rv_AD
    rv_p = rv_p + (one-iVar%q) * rv_AD
    rh_p = rh_p + (one-iVar%q) * rh_AD
    rv_p = rv_p + iVar%q * rh_AD
    rv_AD = 0.0_fp; rh_AD = 0.0_fp
    rh_s_AD = rh_s_AD + iVar%ph * rh_p
    rv_s_AD = rv_s_AD + iVar%pv * rv_p
    rh_p = zero; rv_p = zero
   
  END SUBROUTINE Roughness_Correction_AD

  SUBROUTINE Roughness_PQ_Chen(frequency,sigma,theta,ph,pv,q )

    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma,theta
    REAL(fp), INTENT(OUT) :: q,pv,ph
    REAL(fp) :: a, b,c,ac
    REAL(fp) :: wv, ks,cs,hpp,kkk
    REAL(fp) :: h, hf, lf,xhf,xlf

    REAL(fp) :: smcp 

    wv = 30.0_fp/frequency
    ks = 2.0_fp*PI/wv
 
    cs = cos(theta)
    hpp = 1.0_fp
    h = (Sigma/wv)**hpp
    a = 1.4_fp;  b=0.15_fp; ac=0.43_fp
    smcp = 0.3_fp
 
    lf = 0.0_fp  !lf=(sigma/60.0)**hpp   ! 60, 90
    hf = (sigma/3.0_fp)**hpp     !1.0 in the old test
    xlf = (h-lf)*cs**1.0_fp/smcp ! 0.1 0.3
    xhf = (h-hf)*cs**1.0_fp/0.5_fp
     
    kkk = 2.0_fp*Pi*6.0_fp/30.0_fp*sigma  !f/30=1/w
    ! use the asmpototic value of Wegmuller0.1
    c = exp(-1.0_fp * kkk**(sqrt(0.30_fp * cs)))
   
    ph = (a+c+(b-a)*tanh(xlf)+(c-b)*tanh(xhf))*0.5_fp
    pv = ph
    q  = ac*(1.0_fp - EXP(-15.00_fp*(h)**1.0_fp*cs**2.0_fp))!q=ac=0.35


  END SUBROUTINE Roughness_PQ_Chen
  
  SUBROUTINE  Roughness_PQ_CRTM(frequency,sigma,theta, ph, pv, q)

    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma, theta
    REAL(fp), INTENT(OUT) :: ph, pv, q
    REAL(fp) :: wv, ks, cs, h
 
    wv = 30.0_fp/frequency
    ks = 2.0_fp*PI/wv
 
    h  = (2.0_fp*Ks*sigma)**2
    cs = cos(PI*theta/180.0_fp)
    ph = 0.3_fp
    pv = ph
    
    q  = 0.35_fp*(1.0_fp - EXP(-0.60_fp*frequency*sigma**2))
     
    
  END SUBROUTINE Roughness_PQ_CRTM

 
  SUBROUTINE  Roughness_PQ_Wang(frequency,sigma,theta, ph, pv, q)

    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma, theta
    REAL(fp), INTENT(OUT) :: ph, pv, q
    REAL(fp) :: wv, ks, cs, h
 
    wv = 30.0_fp/frequency
    ks = 2.0_fp*PI/wv
 
    h  = (2.0_fp*Ks*sigma)**2
    cs = cos(PI*theta/180.0_fp)

    ph = 1.0_fp/exp(h*cs**2)
    pv = ph
    q  = 0.0_fp 
    
    
  END SUBROUTINE Roughness_PQ_Wang


  SUBROUTINE  Roughness_PQ_Wigneron(frequency,sigma, ph, pv, q)
    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma
    REAL(fp), INTENT(OUT) :: ph, pv, q
    REAL(fp) :: wv, m, h
 
    wv = 30.0_fp/frequency  
    m  = 0.02_fp*sigma**0.85_fp
    h  = 1.3972_fp*(m)**0.5879_fp

    ph = 1.0_fp/exp(h)
    pv = ph
    q  = 0.0_fp
  
    
  END SUBROUTINE  Roughness_PQ_Wigneron
  
  SUBROUTINE  Roughness_PQ_Wigneron_2011(frequency,sigma,theta,ph, pv, q)
    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma, theta
    REAL(fp), INTENT(OUT) :: pv, ph, q
    REAL(fp) :: wv, ks, cs, h
    REAL(fp) :: a1, a2, a3, a4, b1, b2
     
    wv = 30.0_fp/frequency
    ks = 2.0_fp*PI/wv
    cs = cos(theta)
  
    a1 = 0.880_fp;  a2 = 0.775_fp
    a3 = 3.748_fp;  a4 =-0.100_fp
    b1 = 1.619_fp;  b2 = 0.134_fp
    h  = (a1*sigma*10.0_fp/(a2*sigma*10.0_fp+a3))**6.0_fp

    ph = 1.0_fp/exp(h*cs**b1)
    pv = 1.0_fp/exp(h*cs**b2)
    q  = a4  
    
    
  END SUBROUTINE Roughness_PQ_Wigneron_2011


  SUBROUTINE Roughness_PQ_Coppo(frequency,sigma,ph, pv, q)
    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma
    REAL(fp), INTENT(OUT) :: ph, pv, q
    REAL(fp) :: wv, h
     
    wv = 30.0_fp/frequency
    h  = 3.0_fp*SQRT(sigma/wv)
    ph = 1.0_fp/exp(h)
    pv = ph
    q  = 0.0_fp
  END SUBROUTINE Roughness_PQ_Coppo
  
  SUBROUTINE Roughness_Wegmuller(frequency,sigma,theta,rh_s, rh,rv)
    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: sigma, theta
    REAL(fp), INTENT(IN)  :: rh_s
    REAL(fp), INTENT(OUT) :: rh, rv
    REAL(fp) :: wv, ks, cs, rhm, rvm
    REAL(fp) :: angle
      
    wv = 30.0_fp/frequency
    ks = 2.0_fp*PI/wv
    cs = cos(theta)
    angle = theta/PI*180.0_fp
    rhm = rh_s * exp(-1.0_fp * (ks*sigma)**(sqrt(0.10_fp * cs)))
    IF (angle <= 60.0_fp) THEN 
       rvm = rhm * (cs**0.655_fp)
    ELSE
       rvm = rhm * (0.635_fp - 0.0014_fp*(angle-60.0_fp))
    ENDIF
    rh=rhm ; rv=rvm
  END SUBROUTINE Roughness_Wegmuller


END MODULE SoilMW_Roughness_Correction



