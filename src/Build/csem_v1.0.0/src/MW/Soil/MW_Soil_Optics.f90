!
! MW_Soil_Optics
!
! Module to compute the soil optical properties for LAND surfaces at
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



MODULE MW_Soil_Optics

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE MW_Soil_Permittivity, ONLY : iVar_soil_diel=>iVar_type, &
                                   CSEM_SoilMW_Permittivity, &
                                   CSEM_SoilMW_Permittivity_TL, &
                                   CSEM_SoilMW_Permittivity_AD
  USE CSEM_Fresnel

  IMPLICIT NONE
  PRIVATE
  PUBLIC :: CSEM_SoilMW_Optics
  PUBLIC :: CSEM_SoilMW_Optics_TL
  PUBLIC :: CSEM_SoilMW_Optics_AD
  !PUBLIC :: CSEM_SoilMW_Reflectance_TL
  !PUBLIC :: CSEM_SoilMW_Reflectance_AD
  !PUBLIC :: CSEM_SoilMW_Reflectance
  PUBLIC :: iVar_type
  PUBLIC :: Max_Soil_Layers
 
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: ZERO = 0.0, ONE = 1.0
  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643_fp 
  ! set soil model common constants
  

  INTEGER, PARAMETER  :: Max_Soil_Layers = 1

  INTEGER, PARAMETER  :: N_Soil_Models = 3
  CHARACTER(LEN=20), PARAMETER :: Soil_Models(N_Soil_Models) =   &
         [CHARACTER(LEN=20) :: 'CRTM3.0.0', 'CRTM2.1.3', 'Burke']
 
  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! Forward model input values
    REAL(fp) :: freq  = ZERO
    REAL(fp) :: tskin, teff
    REAL(fp) :: tsoil(Max_Soil_Layers)
    REAL(fp) :: smc(Max_Soil_Layers)
    REAL(fp) :: sand, clay
    REAL(fp) :: theta_t, theta_i, arcs
    COMPLEX(fp) :: eair, esoil
    INTEGER  ::  n_layers = 1
    TYPE(iVar_soil_diel) :: diel

    INTEGER,  PUBLIC :: Model_Option    = 1
    REAL(fp), PUBLIC :: rho_b
   END TYPE iVar_type

CONTAINS

  SUBROUTINE CSEM_SoilMW_Optics(frequency, theta, Tskin, Tsoil, smc, sand, clay, &
                                refl_smooth, teff,  iVar)

    REAL(fp), INTENT(IN)  :: frequency
    REAL(fp), INTENT(IN)  :: theta
    REAL(fp), INTENT(IN)  :: tskin
    REAL(fp), INTENT(IN)  :: tsoil(:), smc(:)
    REAL(fp), INTENT(IN)  :: sand, clay
    REAL(fp), INTENT(OUT) :: refl_smooth(2), teff

    TYPE(iVar_type)  ::  iVar

    INTEGER  :: nLayer_Burke = 120

    iVar%n_layers  = size(tsoil)
    IF(iVar%n_layers > Max_Soil_Layers) THEN
      WRITE(*,*) "Soil layers larger than maximum allowed ...."
      RETURN
    ENDIF
    iVar%freq      =  frequency
    iVar%theta_i   =  theta
    iVar%tskin     =  Tskin
    iVar%sand      =  sand
    iVar%clay      =  clay
    iVar%tsoil(1:iVar%n_layers) =  Tsoil(:)
    iVar%smc(1:iVar%n_layers)   =  smc(:)
    ivar%teff =  Tsoil(1)
    
    SELECT CASE (TRIM(Soil_Models(iVar%Model_Option)))
      CASE('Burke')
        iVar%diel%Model_Option = 3
        CALL ML_Burke_SoilMW_Reflectance(nLayer_Burke, &
            refl_smooth(1), refl_smooth(2),  iVar)
      CASE('CRTM2.1.3')
        iVar%diel%Model_Option = 2
        iVar%diel%rhob = iVar%rho_b
        teff = Tsoil(1)
        iVar%teff = teff

        CALL CRTM_SoilMW_Reflectance_V2( refl_smooth(1), refl_smooth(2), iVar)
      CASE('CRTM3.0.0')
        iVar%diel%Model_Option = 1
        teff=iVar%Tsoil(1)+(Tskin-iVar%Tsoil(1))*(iVar%smc(1)/0.40)**0.5
        iVar%teff = teff
        CALL CRTM_SoilMW_Reflectance_V3( refl_smooth(1), refl_smooth(2), iVar)

      CASE Default
         WRITE(*,*)"CSEM_SoilMW_Optics: wrong model option ...."
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Optics
  
  SUBROUTINE CSEM_SoilMW_Optics_TL(Tskin_TL, Tsoil_TL, smc_TL, &
         refl_h_TL, refl_v_TL, teff_TL, iVar)
    REAL(fp),INTENT(IN)  :: tskin_TL, tsoil_TL, smc_TL
    REAL(fp),INTENT(OUT) :: refl_h_TL, refl_v_TL, teff_TL
    TYPE(iVar_type) ::  iVar
    

    SELECT CASE (TRIM(Soil_Models(iVar%Model_Option)))
      CASE('CRTM3.0.0')
        teff_tl  = (1.0-(iVar%smc(1)/0.40)**0.5)*Tsoil_tl + &
               (iVar%smc(1)/0.40)**0.5*Tskin_tl + &
               (iVar%Tskin-iVar%Tsoil(1))*(1.0/(0.4*iVar%smc(1)))**0.5*smc_tl    
        CALL CRTM_SoilMW_Reflectance_V3_TL(smc_TL, teff_TL, refl_h_TL, refl_v_TL, iVar)
      CASE('CRTM2.1.3')
        teff_tl  = Tsoil_tl 
        CALL CRTM_SoilMW_Reflectance_v2_TL(smc_TL, teff_TL, refl_h_TL, refl_v_TL, iVar)
         
      CASE Default
        teff_TL = 0.0 ; refl_h_TL = 0.0 ; refl_v_TL= 0.0
    
    END SELECT
    
 
  END SUBROUTINE CSEM_SoilMW_Optics_TL
  
   SUBROUTINE CSEM_SoilMW_Optics_AD(Tskin_AD, Tsoil_AD, smc_AD,&
         refl_h_AD, refl_v_AD, teff_AD, iVar)
    REAL(fp),INTENT(INOUT) :: tskin_AD, tsoil_AD, smc_AD
    REAL(fp),INTENT(INOUT) :: refl_h_AD, refl_v_AD, teff_AD
    TYPE(iVar_type) ::  iVar
    
    SELECT CASE (TRIM(Soil_Models(iVar%Model_Option)))
      CASE('CRTM3.0.0')
        CALL CRTM_SoilMW_Reflectance_V3_AD(smc_AD, teff_AD, refl_h_AD, refl_v_AD,iVar)
        Tsoil_ad = Tsoil_ad + (1.0-(iVar%smc(1)/0.40)**0.5)* teff_AD
        Tskin_ad = Tskin_ad + (iVar%smc(1)/0.40)**0.5 * teff_AD
        smc_ad   = smc_ad + (iVar%Tskin-iVar%Tsoil(1))*(1.0/(0.4*iVar%smc(1)))**0.5* teff_AD
      CASE('CRTM2.1.3')
        CALL CRTM_SoilMW_Reflectance_V2_AD(smc_AD, teff_AD, refl_h_AD, refl_v_AD,iVar)
        Tsoil_ad = Tsoil_ad +  teff_AD
        Tskin_ad = Tskin_ad + ZERO
        smc_ad   = smc_ad   + ZERO
      CASE Default
        Tsoil_ad = Tsoil_ad + ZERO
        Tskin_ad = Tskin_ad + ZERO
        smc_ad   = smc_ad   + ZERO
    END SELECT
    
    teff_ad = ZERO
    refl_h_AD = ZERO ; refl_v_AD = ZERO
   
 
  END SUBROUTINE CSEM_SoilMW_Optics_AD
 
  SUBROUTINE CRTM_SoilMW_Reflectance_V2(refl_h, refl_v, iVar)
    REAL(fp),INTENT(OUT) :: refl_h,refl_v
    TYPE(iVar_type) ::  iVar
    iVar%eair = CMPLX(1.0,0.0,fp)
    CALL CSEM_SoilMW_Permittivity(iVar%freq, iVar%teff, iVar%smc(1), &
       iVar%sand, iVar%clay, iVar%esoil, iVar%diel)
    iVar%theta_t = ASIN(REAL(SIN(iVar%theta_i)*SQRT(iVar%eair)/SQRT(iVar%esoil),fp))
    CALL Fresnel_Reflectance(iVar%eair, iVar%esoil, iVar%theta_i, iVar%theta_t, refl_v, refl_h)
  END SUBROUTINE CRTM_SoilMW_Reflectance_V2

  SUBROUTINE CRTM_SoilMW_Reflectance_V2_TL(smc_TL,tsoil_TL,refl_h_TL,refl_v_TL, iVar)
    REAL(fp),INTENT(IN)  :: smc_TL, tsoil_TL
    REAL(fp),INTENT(OUT) :: refl_h_TL, refl_v_TL
    COMPLEX(fp) :: eair_TL, esoil_TL
    REAL(fp) :: theta_t_tl
    TYPE(iVar_Type) :: iVar
   
    eair_TL = 0.0
    CALL CSEM_SoilMW_Permittivity_TL(tsoil_TL, smc_TL,  esoil_tl, iVar%diel)
    iVar%arcs = ONE/SQRT(ONE-(REAL(SIN(iVar%theta_i)*SQRT(iVar%eair)/SQRT(iVar%esoil)))**2)
    theta_t_tl = -2.0/3.0*iVar%arcs*REAL(SIN(iVar%theta_i)*SQRT(iVar%eair)/iVar%esoil*esoil_tl)
    
    CALL Fresnel_Reflectance_TL(iVar%eair, iVar%esoil, iVar%theta_i, eair_TL, esoil_TL, refl_v_TL, refl_h_TL)

  END SUBROUTINE CRTM_SoilMW_Reflectance_V2_TL

  SUBROUTINE CRTM_SoilMW_Reflectance_V2_AD(smc_AD,tsoil_AD,refl_h_AD, refl_v_AD, iVar)
    REAL(fp),INTENT(INOUT)  :: smc_AD, tsoil_AD
    REAL(fp),INTENT(INOUT)  :: refl_h_AD,refl_v_AD
    COMPLEX(fp) :: eair_AD, esoil_AD
    REAL(fp)    :: theta_t_ad
    TYPE(iVar_Type) :: iVar
      
    eair_AD = 0.0 ; esoil_AD = 0.0
    CALL Fresnel_Reflectance_AD(iVar%eair, iVar%esoil, iVar%theta_i,&
       eair_AD, esoil_AD, refl_v_AD, refl_h_AD)
    theta_t_ad = -2.0/3.0*iVar%arcs*REAL(SIN(iVar%theta_i)*SQRT(iVar%eair)/iVar%esoil*esoil_ad)
    CALL CSEM_SoilMW_Permittivity_AD(tsoil_AD,smc_AD,  esoil_AD, iVar%diel)
    refl_h_AD = 0.0; refl_v_AD = 0.0
 
  END SUBROUTINE CRTM_SoilMW_Reflectance_V2_AD
  
  SUBROUTINE CRTM_SoilMW_Reflectance_V3(refl_h, refl_v, iVar)
    REAL(fp),INTENT(OUT) :: refl_h,refl_v
    TYPE(iVar_type) ::  iVar
   
    iVar%eair = CMPLX(1.0,0.0,fp)
    CALL CSEM_SoilMW_Permittivity(iVar%freq, iVar%teff, iVar%smc(1), &
       iVar%sand, iVar%clay, iVar%esoil, iVar%diel)
    CALL Fresnel_Reflectance(iVar%eair, iVar%esoil, iVar%theta_i,  refl_v, refl_h)

  END SUBROUTINE CRTM_SoilMW_Reflectance_V3
  
  SUBROUTINE CRTM_SoilMW_Reflectance_V3_TL(smc_TL,tsoil_TL,refl_h_TL,refl_v_TL,iVar)

    REAL(fp),INTENT(IN)  :: smc_TL, tsoil_TL
    REAL(fp),INTENT(OUT) :: refl_h_TL, refl_v_TL
    COMPLEX(fp) :: eair_TL, esoil_TL
    TYPE(iVar_type) ::  iVar
   
    eair_TL = 0.0
    CALL CSEM_SoilMW_Permittivity_TL(tsoil_TL,smc_TL,  esoil_tl, iVar%diel)
    CALL Fresnel_Reflectance_TL(iVar%eair, iVar%esoil, iVar%theta_i, &
           eair_TL, esoil_TL, refl_v_TL, refl_h_TL)

  END SUBROUTINE CRTM_SoilMW_Reflectance_V3_TL
  
  SUBROUTINE CRTM_SoilMW_Reflectance_V3_AD(smc_AD,tsoil_AD,refl_h_AD,refl_v_AD,iVar)
    REAL(fp),INTENT(INOUT)  :: smc_AD, tsoil_AD
    REAL(fp),INTENT(INOUT)  :: refl_h_AD,refl_v_AD
    COMPLEX(fp) :: eair_AD, esoil_AD
    TYPE(iVar_type) ::  iVar
   
    eair_AD = 0.0 ; esoil_AD = 0.0
    CALL Fresnel_Reflectance_AD(iVar%eair, iVar%esoil, iVar%theta_i, eair_AD, esoil_AD, refl_v_AD, refl_h_AD)
    CALL CSEM_SoilMW_Permittivity_AD(tsoil_AD,smc_AD,  esoil_AD, iVar%diel)
    refl_h_AD = 0.0; refl_v_AD = 0.0
 
  END SUBROUTINE CRTM_SoilMW_Reflectance_V3_AD
  
  SUBROUTINE ML_Burke_SoilMW_Reflectance( Nlayer, refl_h, refl_v, iVar)
      
    REAL(fp), INTENT(OUT) :: refl_h, refl_v
    INTEGER , INTENT(IN)  :: Nlayer
    TYPE(iVar_type)       :: iVar
    
    REAL(fp) :: tl, vmc
    REAL(fp) :: tlsm(2),wlsm(2)
    REAL(fp) :: k_real, k_img, k0
    REAL(fp) :: rv, rh, tb_h, tb_v

    REAL(fp),   DIMENSION(Nlayer) :: Q, T
    REAL(fp),   DIMENSION(Nlayer) :: Rp_H, Rp_V
    REAL(fp),   DIMENSION(Nlayer) :: gama, deltaz, zsoil

    COMPLEX(fp):: esm1,esm2
    INTEGER  :: I
    
    tlsm=(/iVar%tskin,iVar%tsoil(1)/)
    wlsm=(/iVar%smc(1)*0.15,iVar%smc(1)/) !0.62
    IF(nlayer > 2 ) THEN
      CALL SOIL_PROFILE(nlayer,T,Q,deltaz,zsoil,tlsm,wlsm) 
      deltaz=deltaz*10.0
    ELSE 
      T = tlsm
      Q = MAX(MIN(iVar%smc(1)*1.5_fp,1.0_fp),0.0_fp)
      deltaz=(/10.0,100.0/)! mm
    ENDIF

    k0 =2.0*Pi*iVar%freq/300.0
    esm1 = CMPLX(1.0, 0.0, fp)
    
    DO I=1,Nlayer
      tl=T(I)
      vmc = MAX(MIN(Q(I),1.0_fp),0.0_fp)
      !CALL CSEM_SoilMW_Diel(freq,tl,vmc,sand,clay,esm2)
      CALL CSEM_SoilMW_Permittivity(iVar%freq, tl, vmc, &
       iVar%sand, iVar%clay,  esm2, iVar%diel)
      CALL Dispersion(iVar%theta_i,esm2,k_real,k_img)
      gama(I)=abs(2.0*K0*k_img)

      CALL Fresnel_Reflectance_liou(iVar%theta_i,esm1,esm2,rv,rh)
      Rp_H(I) = Rh ; Rp_V(I) = Rv
      
      esm1 = esm2
    ENDDO
      
    !deltaz(1)=1.0/gama(1)*10.0 !mm
    CALL Burke_Soil_Model(T,deltaz,gama,Rp_H,Tb_H,Nlayer)
    CALL Burke_Soil_Model(T,deltaz,gama,Rp_V,Tb_V,Nlayer)
    
    iVar%teff=max(iVar%tskin,iVar%tsoil(1))
    refl_h=1.0-Tb_H/iVar%teff
    refl_v=1.0-Tb_V/iVar%teff
  
  END SUBROUTINE ML_Burke_SoilMW_Reflectance



  SUBROUTINE Burke_Soil_Model(T,deltaz,gama,Rp,Tb,N)
    REAL(fp) :: TB,Ti
    REAL(fp) :: RR,Rpj,Rpi,Epj,Epi
    INTEGER  :: N,I,J
    REAL(fp) :: gama(:),Rp(:),T(:),deltaz(:)

    TB=0.0_fp
    DO I=1,N 
      RR=1.0
      DO J=1,I 
        IF (J == 1) THEN
          Epj=1.0 ; Rpj=0.0
        ELSE
          Epj=EXP(-1.0*gama(j-1)*deltaz(j-1))
          Rpj=Rp(j-1) 
        ENDIF
        RR=RR*(1.0-Rpj)*Epj
      ENDDO 
       
      IF (I==N) THEN
        Epi=0.0
        Rpi=0.0
      ELSE
        Epi=exp(-1.0*gama(i)*deltaz(i))
        Rpi=Rp(i+1)
      ENDIF
      Ti=T(i)*(1.0-Epi)*(1.0-Rp(i))*(1.0+Rpi*Epi)*RR   
      Tb=Tb+Ti
     
    END DO

  END SUBROUTINE Burke_Soil_Model



  SUBROUTINE SOIL_PROFILE(NL,tsoil,wsoil,dzsoil,zsoil, tlsm, wlsm)

    IMPLICIT NONE
  
    INTEGER, PARAMETER:: nsoil=120
    INTEGER  :: NL
    REAL(fp) :: dzsoil(NL), zsoil(NL),tsoil(NL), wsoil(NL)
    REAL(fp) :: tlsm(2), wlsm(2)
    REAL(fp) :: thicks(nsoil)
    REAL(fp) :: sumdz,e1,e2
    REAL(fp) :: ta,tb,tc
    REAL(fp) :: wa,wb,wc
    INTEGER  :: i, iz
      
    thicks(1:10)=(/(0.1,i=1,10)/)
    thicks(11:30)=(/(0.5,i=1,20)/)
    thicks(31:nsoil)=(/(1.0,i=1,90)/)
    
    IF (NL/= nsoil) THEN
      WRITE(*,*)'Number of Soil Layers should be ',nsoil
      STOP
    ENDIF
      
    ! Declare soil moisture layer parameters
    tc=30.0; wc=30.0
    e1=1.0-exp(-3.0*0.1/tc)
    e2=1.0-exp(-3.0*10.0/tc)
    tb=(tlsm(1)-tlsm(2))/(e1-e2)
    ta=tlsm(1)-tb*e2
    wb=(wlsm(1)-wlsm(2))/(e1-e2)
    wa=wlsm(1)-wb*e2

    sumdz  = 0.0                                                      
    DO iz = 1, nsoil
      sumdz   = sumdz  + thicks(iz)
      zsoil(iz) = sumdz  - thicks(iz)/2. ! depth to the center of layer
      dzsoil(iz) = thicks(iz) 
      tsoil(iz) = ta + tb*(1.0-exp(-3.0*zsoil(iz)/tc))
      wsoil(iz) = wa + wb*(1.0-exp(-3.0*zsoil(iz)/wc)) 
    ENDDO  
    WHERE (wsoil <0.0) wsoil=0.001
 
  END SUBROUTINE SOIL_PROFILE


END MODULE MW_Soil_Optics



