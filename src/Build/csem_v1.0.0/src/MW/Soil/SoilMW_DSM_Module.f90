!
! SoilMW_DSM_Module
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



MODULE SoilMW_DSM_Module

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Fresnel

  IMPLICIT NONE
  PUBLIC
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER, PRIVATE :: PI  = 3.141592653589793238462643_fp 
  REAL(fp), PARAMETER, PRIVATE :: ONE = 1.0_fp  
  REAL(fp), PARAMETER, PRIVATE :: TWO = 2.0_fp  
  REAL(fp), PARAMETER, PRIVATE :: POINT5=0.5_fp  
  REAL(fp), PARAMETER, PRIVATE :: ZERO = 0.0_fp 
  ! set soil model common constants
  
  REAL(fp),    PARAMETER ::  rho_s = 2.66_fp        !   rho_s : soil specific density (g/cm3), ATBD 2.664
  REAL(fp),    PARAMETER ::  eps_winf  = 4.9_fp     !   eps_winf : dielectric constant at infinite frequency (Stogryn 1971)
                                                    !   also called: high frequency dielectric constant (4.9 Lane and Saxton )
  REAL(fp),    PARAMETER ::  eps_0 = 8.854e-12_fp   !   permittivity of free space (Klein and Swift 1977) [Farads/meter]
  COMPLEX(fp), PARAMETER ::  ea = (1.0_fp, 0.0_fp)  !   dielectric constant of air
  COMPLEX(fp), PARAMETER ::  er = (5.5_fp, 0.2_fp)  !   dielectric constant of rock
  COMPLEX(fp), PARAMETER ::  ef = (5.0_fp, 0.5_fp)  !   dielectric constant of frozen soil (Hallikainen et al. 1984, 1985)
  COMPLEX(fp), PARAMETER ::  ei = (3.2_fp, 0.1_fp)  !   dielectric constant of ice

CONTAINS


  SUBROUTINE Soil_DSM(frequency,angle,vmc,Tsoil,TSkin,rh,rv,rhos,rhob,sand,clay,stype)
     REAL(fp)  :: frequency,angle, TSkin,Tsoil,vmc,rhob,rhos,sand,clay
     REAL(fp)  :: theta,mu,salb,gh,gv,theta_1
     COMPLEX(fp) :: esoil_1,eair
     REAL(fp)  ::   ssalb_h,ssalb_v, tau_h,tau_v
     REAL(fp)  ::   r_h,r_v
     REAL(fp)  ::   afactor
     REAL(fp)  ::   eh,ev,rh,rv
     INTEGER   ::   stype
     REAL(fp) ::    XX,YY1,YY2,RR
  
       
     theta=Pi/180.0_fp*angle
     mu = COS(Pi/180.0_fp*theta)
     eair=CMPLX(1.0_fp,-0.0_fp,fp)
     
     IF(Tsoil < TSkin) THEN
        IF(stype==9) vmc=1.0_fp*vmc !as amsr2 0.8
        IF(stype==11)vmc=1.0_fp*vmc
      vmc=0.05_fp !0.06 0.05
     ENDIF
     IF(Tsoil >= TSkin) THEN
        IF(stype==9) vmc=0.8_fp*vmc  !0.5
        IF(stype==11)vmc=0.8_fp*vmc
      vmc=0.03_fp !0.05 0.03
     ENDIF
   
      CALL CRTM_SandMW_Diel_N(Frequency, tsoil, vmc, rhob, rhos, sand, clay, esoil_1)
     CALL Soil_Optic_DSM(frequency,esoil_1,salb,gh,afactor)
  
     gv=gh
     ssalb_h=salb ; ssalb_v=ssalb_h ;  tau_h=0.5_fp ;  tau_v=tau_h

     theta_1 = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(esoil_1),fp))
     CALL Fresnel_Reflectance(eair, esoil_1, theta,  theta_1, r_v, r_h)
     
     XX=(frequency-12.0_fp)/5.0_fp
     YY1=(1.1_fp-1.0_fp)*tanh(XX)
     XX=(frequency-50.0_fp)/30.0_fp
     YY2=(1.5_fp-1.1_fp)*tanh(XX)
     RR=0.5_fp*(1.0_fp+1.5_fp+YY1+YY2)
   
     ev = ( one - r_v ) * ( two * afactor*RR / ( ( one+ afactor*RR) - ( one -afactor*RR )*r_v )  )
     
     XX=(frequency-12.0_fp)/5.0_fp
     YY1=(0.90_fp-1.2_fp)*tanh(XX)
     XX=(frequency-50.0_fp)/30.0_fp
     YY2=(0.35_fp-0.90_fp)*tanh(XX)
     RR=0.5_fp*(1.2_fp+0.35_fp+YY1+YY2)
     
     eh= ( one - r_h ) * ( two * afactor*RR / ( ( one + afactor*RR) - ( one -afactor*RR )*r_h )  )
     
     rh=1.0_fp-eh  ; rv=1.0_fp-ev
 

  END SUBROUTINE Soil_DSM


  SUBROUTINE Soil_Optic_DSM(frequency,esoil,salb,gh,afactor)

    REAL(fp) :: frequency,radius,f,salb,gh,pdepth
    REAL(fp) :: k0,Ka,Ks,r,Kr,YR,YI,ff,Qa,Qs,n0,wv
    REAL(fp) :: ftmp,omega,g,afactor
  
    COMPLEX(fp) :: esoil

    !radius=0.3 ; f=0.8 !f=0.6
    radius=0.3_fp ; f=0.6_fp !f=0.6,0.8 for sand, 0.1 for normal soil
    k0 = 2.0_fp *Pi*frequency/300.0_fp 
    r=radius
    Kr=K0*r
    n0=f*(3.0_fp /4.0_fp )/(Pi*r**3)
    YR=(REAL(esoil)-1.0_fp )/(REAL(esoil)+2.0_fp )
    YI=3.0_fp *abs(AIMAG(esoil))/(REAL(esoil)+2.0)**2

    ff=1.0_fp /(((1.0_fp -f*YR)**1.5_fp )*(one+two*f*YR)**0.5_fp )
    Qa=4.0_fp*(k0*r)*YI*ff
    Qs=8.0_fp/3.0_fp*(k0*r)**4*YR**2*(one-f)**4.0_fp/(one+two*f)**2*ff
 
    Ka=n0*Pi*r**2*Qa
    Ks=n0*Pi*r**2*Qs
 
    salb=Ks/(Ka+Ks)
    gh=0.23_fp*(k0*r)**2
    wv = 30.0_fp/frequency
  
    pdepth=  wv*SQRT(REAL(esoil))/(2.0_fp*Pi*ABS(AIMAG(esoil))*f**1.7_fp)
  
  !
    ftmp = ( 1.0_fp - f )**4
    ftmp = ftmp * kr**3 * yR * yR 
    omega = ftmp / (ftmp  + 1.5_fp * ( one + two * f )**2 * yI )   ! Eq. A16
    g = 0.23_fp * kr * kr
    afactor = Sqrt ( ( one - omega) / (one - omega * g) )    ! Eq. 3b

  
   END SUBROUTINE Soil_Optic_DSM

   SUBROUTINE CRTM_SandMW_Diel(freq,esm)

    REAL(fp) :: f,freq
    REAL(fp) :: eswi,esof,ap, f0
    COMPLEX(fp) :: esm,esf,esa
    REAL(fp) :: ONE=1.0_fp,ZERO=0.0_fp
    
 
    ap=0.002_fp
    f0=0.27_fp
    esa=CMPLX(ap, ZERO, fp)
    f = freq

    ! the permittivity at the high frequency limit
    eswi = 2.53_fp

    ! the permittivity of free space (esof)
    esof = 2.79_fp


    esf = CMPLX(esof-eswi, ZERO, fp)/CMPLX(ONE, -1.0_fp *f/f0, fp)
    esm = CMPLX(eswi, ZERO, fp)+esf+esa

    IF (AIMAG(esm) >= ZERO) esm = CMPLX(REAL(esm,fp),-0.0001_fp, fp)

   END SUBROUTINE CRTM_SandMW_Diel

   SUBROUTINE CRTM_SandMW_Diel_N(freq,t_soil,vmc,rhob,rhos,sand,clay,esm)

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

    REAL(fp) :: f,tauw,freq,t_soil,vmc,rhob,rhos,sand,clay
    REAL(fp) :: alpha,beta,ess,rhoef,t,eswi,eswo
    REAL(fp) :: esof,vmu,vmi
    COMPLEX(fp) :: esm,esw,es1,es2,esi
    REAL(fp) :: ONE=1.0_fp,ZERO=0.0_fp,TWOPI=2.0_fp*Pi
    REAL(fp) :: A=0.0_fp,B=1.0_fp,XX,YY,RR
 
      
    esi=CMPLX(3.15_fp , ZERO, fp)
    
    XX=(t_soil-275.0_fp )/1.2_fp 
    YY=(B-A)*tanh(XX)
    RR=0.5_fp *(A+B+YY)
    
    vmu=vmc*RR
    vmi=vmc-vmu
    !rhob = sand*1.6_fp +  clay*1.1_fp + (1._fp-sand-clay)*1.2_fp
     
    alpha = 0.65_fp
    alpha = 0.60_fp
    beta  = 1.09_fp - 0.11_fp*sand + 0.18_fp*clay
    ess = (1.01_fp + 0.44_fp*rhos)**2 - 0.062_fp
    rhoef = -1.645_fp + 1.939_fp*rhob - 0.020213_fp*sand + 0.01594_fp*clay
    !rhoef = -1.645_fp + 1.939_fp*rhob - 2.0213_fp*sand + 1.594_fp*clay
    t = t_soil - 273.0_fp
    f = freq*1.0e9_fp

    ! the permittivity at the high frequency limit
    eswi = 5.5_fp
    !eswi = 4.9_fp

    ! the permittivity of free space (esof)
    esof = 8.854e-12_fp

    ! static dieletric constant (eswo)
    eswo = 87.134_fp+(-1.949e-1_fp+(-1.276e-2_fp+2.491e-4_fp*t)*t)*t
    tauw = 1.1109e-10_fp+(-3.824e-12_fp+(6.938e-14_fp-5.096e-16_fp*t)*t)*t
    eswo = 87.134_fp+(-1.949e-1_fp+(-1.276e-2_fp+2.491e-4_fp*t)*t)*t
    tauw = 1.1109e-10_fp+(-3.824e-12_fp+(6.938e-14_fp-5.096e-16_fp*t)*t)*t


    IF (vmu > ZERO) THEN
      es1 = CMPLX(eswi, -rhoef*(rhos-rhob)/(TWOPI*f*esof*rhos*vmu), fp)
    ELSE
      es1 = CMPLX(eswi, ZERO, fp)
    ENDIF

    es2 = CMPLX(eswo-eswi, ZERO, fp)/CMPLX(ONE, f*tauw, fp)
    esw = es1 + es2
    esm = ONE + (ess**alpha - ONE)*rhob/rhos + vmu**beta*esw**alpha - vmu+vmi*esi**alpha
    esm = esm**(ONE/alpha)

    IF (AIMAG(esm) >= ZERO) esm = CMPLX(REAL(esm,fp),-0.0001_fp, fp)

  END SUBROUTINE CRTM_SandMW_Diel_N


END MODULE SoilMW_DSM_Module



