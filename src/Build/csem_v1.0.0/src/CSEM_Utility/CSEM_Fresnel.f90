!
! CSEM_Fresnel
!
! Module containing several algorithms for the calculation of 
! Fresnel Reflectance and transmittance. The associated Tangent linear (TL)
! and adjoint(AD) modes are provided to support the CSEM applications 
! in variational data assimlation and retrieval calculations
!
!
! CREATION HISTORY:
!
!       Written by:    Ming Chen 08-May-2016
!                       Ming.Chen@noaa.gov
!

MODULE CSEM_Fresnel

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp

  IMPLICIT NONE
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER, PRIVATE :: PI     = 3.141592653589793238462643_fp
  COMPLEX(fp), PARAMETER, PRIVATE :: CZERO = CMPLX(0.0, 0.0,fp)
  REAL(fp),    PARAMETER, PRIVATE :: ZERO  = 0.0_fp, ONE  = 1.0_fp
  
  INTERFACE Fresnel_Reflectance
      MODULE PROCEDURE  Fresnel_Reflectance_1
      MODULE PROCEDURE Fresnel_Reflectance_2
  END INTERFACE Fresnel_Reflectance
  
  INTERFACE Fresnel_Reflectance_TL
      MODULE PROCEDURE Fresnel_Reflectance_TL_1
      MODULE PROCEDURE Fresnel_Reflectance_TL_2
  END INTERFACE Fresnel_Reflectance_TL
  
  INTERFACE Fresnel_Reflectance_AD
      MODULE PROCEDURE Fresnel_Reflectance_AD_1
      MODULE PROCEDURE Fresnel_Reflectance_AD_2
  END INTERFACE Fresnel_Reflectance_AD
  
CONTAINS

  SUBROUTINE Fresnel_Reflectance_1(em1, em2, theta_i, theta_t, rv, rh)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Reflectance compute the surface reflectivity
!
!
! abstract: compute the surface reflectivety using fresnel equations
!    for a rough surface having a standard deviation of height of sigma
!
! program history log:
!
! input argument list:
!      theta_i      -  incident angle (degree)
!      theta_t      -  transmitted angle (degree)
!      em1          -  dielectric constant of the medium 1
!      em2          -  dielectric constant of the medium 2
!
! output argument list:
!
!      rv           -  reflectivity at vertical polarization
!      rh           -  reflectivity at horizontal polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------
    REAL(fp),    INTENT(IN)  :: theta_i, theta_t
    REAL(fp),    INTENT(OUT) :: rh, rv
    COMPLEX(fp), INTENT(IN)  :: em1, em2

    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t

    cos_i = COS(theta_i); cos_t = COS(theta_t)

    m1 = SQRT(em1); m2 = SQRT(em2)

    rv = (ABS((m1*cos_t-m2*cos_i)/(m1*cos_t+m2*cos_i)))**2
    rh = (ABS((m1*cos_i-m2*cos_t)/(m1*cos_i+m2*cos_t)))**2

  END SUBROUTINE Fresnel_Reflectance_1


  SUBROUTINE Fresnel_Reflectance_TL_1(em1, em2, theta_i, theta_t, &
             em1_TL, em2_TL, rv_TL, rh_TL)

    REAL(fp),    INTENT(IN)  :: theta_i, theta_t
    COMPLEX(fp), INTENT(IN)  :: em1, em2, em1_TL, em2_TL
    REAL(fp),    INTENT(OUT) :: rh_TL, rv_TL

    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t

    COMPLEX(fp) :: A, A_TL, B, B_TL, C, C_TL, D, D_TL
    COMPLEX(fp) :: m1_TL, m2_TL, frv_TL, frh_TL

    cos_i = COS(theta_i); cos_t = COS(theta_t)

    m1 = SQRT(em1) ; m2 = SQRT(em2)
    m1_TL = 0.5_fp/SQRT(em1)*em1_TL
    m2_TL = 0.5_fp/SQRT(em2)*em2_TL
    
    A = m1*cos_t-m2*cos_i
    B = m1*cos_t+m2*cos_i
    A_TL = m1_TL*cos_t-m2_TL*cos_i
    B_TL = m1_TL*cos_t+m2_TL*cos_i
    frv_TL = A_TL/B-A*B_TL/B**2
    Rv_TL = 2.0_fp * (REAL(A/B)*REAL(frv_TL)+ AIMAG(A/B)*AIMAG(frv_TL))
        
    C = m1*cos_i-m2*cos_t
    D = m1*cos_i+m2*cos_t
    C_TL = m1_TL*cos_i-m2_TL*cos_t
    D_TL = m1_TL*cos_i+m2_TL*cos_t
    frh_TL = C_TL/D-C*D_TL/D**2
    Rh_TL = 2.0_fp * (REAL(C/D)*REAL(frh_TL)+ AIMAG(C/D)*AIMAG(frh_TL))
    
  END SUBROUTINE Fresnel_Reflectance_TL_1
  
 
  SUBROUTINE Fresnel_Reflectance_AD_1(em1, em2, theta_i, theta_t, &
             em1_AD, em2_AD, rv_AD, rh_AD)

    REAL(fp),    INTENT(IN)    :: theta_i, theta_t
    COMPLEX(fp), INTENT(IN)    :: em1, em2
    REAL(fp),    INTENT(INOUT) :: rh_AD, rv_AD
    COMPLEX(fp), INTENT(INOUT) :: em1_AD, em2_AD

    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t

    COMPLEX(fp) :: A, A_AD, B, B_AD, C, C_AD, D, D_AD
    COMPLEX(fp) :: m1_AD, m2_AD, frv_AD, frh_AD

    cos_i = COS(theta_i) ; cos_t = COS(theta_t)

    m1 = SQRT(em1) ; m2 = SQRT(em2)
    m1_AD = CZERO; m2_AD = CZERO
     
    A = m1*cos_t-m2*cos_i
    B = m1*cos_t+m2*cos_i
    frv_AD = 2.0_fp *CMPLX(REAL(A/B)*Rv_AD,AIMAG(A/B)*Rv_AD,fp)
    Rv_AD =  ZERO
    A_AD  =  CONJG(1.0_fp/B)  * frv_AD
    B_AD  = -CONJG(A/B**2) * frv_AD
    m1_AD = m1_AD + cos_t * B_AD
    m2_AD = m2_AD + cos_i * B_AD
    m1_AD = m1_AD + cos_t * A_AD
    m2_AD = m2_AD - cos_i * A_AD
    A_AD  = CZERO; B_AD = CZERO
    
    C = m1*cos_i-m2*cos_t
    D = m1*cos_i+m2*cos_t
    frh_AD =  2.0_fp *CMPLX(REAL(C/D)*Rh_AD,AIMAG(C/D)*Rh_AD,fp)
    Rh_AD =  ZERO
    
    C_AD  =  CONJG(1.0_fp/D)  * frh_AD
    D_AD  = -CONJG(C/D**2) * frh_AD
    m1_AD = m1_AD + cos_i * D_AD
    m2_AD = m2_AD + cos_t * D_AD
    m1_AD = m1_AD + cos_i * C_AD
    m2_AD = m2_AD - cos_t * C_AD
    C_AD  = CZERO; D_AD = CZERO
  
    em1_AD = em1_AD +  CONJG(0.5_fp/SQRT(em1)) * m1_AD
    em2_AD = em2_AD +  CONJG(0.5_fp/SQRT(em2)) * m2_AD
    m1_AD = CZERO; m2_AD = CZERO
  END SUBROUTINE Fresnel_Reflectance_AD_1

  SUBROUTINE Fresnel_Reflectance_2(em1, em2, theta_i, rv, rh)

    REAL(fp),    INTENT(IN)  :: theta_i
    REAL(fp),    INTENT(OUT) :: rh, rv
    COMPLEX(fp), INTENT(IN)  :: em1, em2

    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t
    REAL(fp)    :: theta_t

    theta_t = ASIN(REAL(SIN(theta_i)*SQRT(em1)/SQRT(em2),fp))
    cos_i = COS(theta_i); cos_t = COS(theta_t)

    m1 = SQRT(em1); m2 = SQRT(em2)

    rv = (ABS((m1*cos_t-m2*cos_i)/(m1*cos_t+m2*cos_i)))**2
    rh = (ABS((m1*cos_i-m2*cos_t)/(m1*cos_i+m2*cos_t)))**2

  END SUBROUTINE Fresnel_Reflectance_2
  
  SUBROUTINE Fresnel_Reflectance_TL_2(em1, em2, theta_i, em1_TL, &
             em2_TL, rv_TL, rh_TL)

    REAL(fp),    INTENT(IN)  :: theta_i
    COMPLEX(fp), INTENT(IN)  :: em1, em2, em1_TL, em2_TL
    REAL(fp),    INTENT(OUT) :: rh_TL, rv_TL
   
    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t, theta_t
    REAL(fp)    :: theta_t_TL, cos_t_TL

    COMPLEX(fp) :: A, A_TL, B, B_TL, C, C_TL, D, D_TL
    COMPLEX(fp) :: m1_TL, m2_TL, frv_TL, frh_TL
    COMPLEX(fp) :: ctheta, ctheta_TL 
    
 
    ctheta  = SIN(theta_i)*SQRT(em1)/SQRT(em2)
    theta_t = ASIN(REAL(ctheta,fp))
    cos_i = COS(theta_i); cos_t = COS(theta_t)
    
    ctheta_TL  =  SIN(theta_i)/2.0_fp*SQRT(em2/em1)*(1.0_fp / &
                  em2*em1_TL-em1/em2**2*em2_TL)
    theta_t_TL =  1.0_fp/SQRT(1.0_fp-(REAL(ctheta,fp))**2) * REAL(ctheta_TL,fp)
    cos_t_TL = -REAL(ctheta,fp)*theta_t_TL
    
    m1 = SQRT(em1) ; m2 = SQRT(em2)
    m1_TL = 0.5_fp/SQRT(em1)*em1_TL
    m2_TL = 0.5_fp/SQRT(em2)*em2_TL
    
    A = m1*cos_t-m2*cos_i
    B = m1*cos_t+m2*cos_i
    A_TL = m1_TL*cos_t-m2_TL*cos_i+m1*cos_t_TL
    B_TL = m1_TL*cos_t+m2_TL*cos_i+m1*cos_t_TL
    frv_TL = A_TL/B-A*B_TL/B**2
    Rv_TL = 2.0_fp * (REAL(A/B)*REAL(frv_TL)+ AIMAG(A/B)*AIMAG(frv_TL))
        
    C = m1*cos_i-m2*cos_t
    D = m1*cos_i+m2*cos_t
    C_TL = m1_TL*cos_i-m2_TL*cos_t-m2*cos_t_TL
    D_TL = m1_TL*cos_i+m2_TL*cos_t+m2*cos_t_TL
    frh_TL = C_TL/D-C*D_TL/D**2
    Rh_TL = 2.0_fp * (REAL(C/D)*REAL(frh_TL)+ AIMAG(C/D)*AIMAG(frh_TL))
    
  END SUBROUTINE Fresnel_Reflectance_TL_2
  
  SUBROUTINE Fresnel_Reflectance_AD_2(em1, em2, theta_i,  &
             em1_AD, em2_AD, rv_AD, rh_AD)

    REAL(fp),    INTENT(IN)    :: theta_i
    COMPLEX(fp), INTENT(IN)    :: em1, em2
    REAL(fp),    INTENT(INOUT) :: rh_AD, rv_AD
    COMPLEX(fp), INTENT(INOUT) :: em1_AD, em2_AD

    COMPLEX(fp) :: m1, m2 
    REAL(fp)    :: cos_i, cos_t, theta_t
    REAL(fp)    :: cos_t_AD, theta_t_AD

    COMPLEX(fp) :: A, A_AD, B, B_AD, C, C_AD, D, D_AD
    COMPLEX(fp) :: m1_AD, m2_AD, frv_AD, frh_AD
    COMPLEX(fp) :: ctheta, ctheta_AD 

    ctheta  = SIN(theta_i)*SQRT(em1)/SQRT(em2)
    theta_t = ASIN(REAL(ctheta,fp))
    cos_i = COS(theta_i) ; cos_t = COS(theta_t)

    cos_t_AD = 0.0_fp; theta_t_AD = 0.0_fp
    ctheta_AD = CZERO
    
    m1 = SQRT(em1) ; m2 = SQRT(em2)
    m1_AD = CZERO; m2_AD = CZERO
    
    A = m1*cos_t-m2*cos_i
    B = m1*cos_t+m2*cos_i
    frv_AD = 2.0_fp *CMPLX(REAL(A/B)*Rv_AD,AIMAG(A/B)*Rv_AD,fp)
    Rv_AD =  ZERO
    A_AD  =  CONJG(1.0/B)  * frv_AD
    B_AD  = -CONJG(A/B**2) * frv_AD
    m1_AD = m1_AD + cos_t * B_AD
    m2_AD = m2_AD + cos_i * B_AD
    m1_AD = m1_AD + cos_t * A_AD
    m2_AD = m2_AD - cos_i * A_AD
    cos_t_AD = cos_t_AD + REAL(CONJG(m1) * A_AD,fp)
    cos_t_AD = cos_t_AD + REAL(CONJG(m1) * B_AD,fp)
   
    A_AD  = CZERO; B_AD = CZERO
    
    C = m1*cos_i-m2*cos_t
    D = m1*cos_i+m2*cos_t
    frh_AD =  2.0_fp *CMPLX(REAL(C/D)*Rh_AD,AIMAG(C/D)*Rh_AD,fp)
    Rh_AD =  ZERO
    
    C_AD  =  CONJG(1.0/D)  * frh_AD
    D_AD  = -CONJG(C/D**2) * frh_AD
    m1_AD = m1_AD + cos_i * D_AD
    m2_AD = m2_AD + cos_t * D_AD
    m1_AD = m1_AD + cos_i * C_AD
    m2_AD = m2_AD - cos_t * C_AD
    cos_t_AD = cos_t_AD - REAL(CONJG(m2) * C_AD,fp)
    cos_t_AD = cos_t_AD + REAL(CONJG(m2) * D_AD,fp)
    C_AD  = CZERO; D_AD = CZERO
  
    theta_t_AD = theta_t_AD - REAL(ctheta,fp)*cos_t_AD
    ctheta_AD  = ctheta_AD + CMPLX(1.0_fp/SQRT(1.0_fp-(REAL(ctheta,fp))**2)*&
                 theta_t_AD,0.0_fp,fp)
    
    em1_AD = em1_AD +  CONJG(0.5_fp/SQRT(em1)) * m1_AD
    em2_AD = em2_AD +  CONJG(0.5_fp/SQRT(em2)) * m2_AD
    em1_AD = em1_AD +  CONJG(SIN(theta_i)/2.0_fp*SQRT(em2/em1)*(1.0_fp/em2)) * ctheta_AD
    em2_AD = em2_AD -  CONJG(SIN(theta_i)/2.0_fp*SQRT(em2/em1)*(em1/em2**2)) * ctheta_AD
    
    m1_AD = CZERO; m2_AD = CZERO

  END SUBROUTINE Fresnel_Reflectance_AD_2

  SUBROUTINE Fresnel_Transmittance(em1,em2,theta_i,theta_t,tv,th)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Transmittance    calculate Transmittance
!
!
! abstract: compute Transmittance
!
! program history log:
!
! input argument list:
!
!      theta        -  local zenith angle (degree)
!      theta_i      -  incident angle (degree)
!      theta_t      -  transmitted angle (degree)
!      em1          -  dielectric constant of the medium 1
!      em2          -  dielectric constant of the medium 2
!
! output argument list:
!
!      tv           -  transmisivity at vertical polarization
!      th           -  transmisivity at horizontal polarization
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

    REAL(fp),    INTENT(IN)  :: theta_i, theta_t
    COMPLEX(fp), INTENT(IN)  :: em1, em2 
    REAL(fp),    INTENT(OUT) :: th, tv 

    COMPLEX(fp) ::  m1, m2 
    REAL(fp)    ::  rr, cos_i,cos_t

    cos_i = COS(theta_i); cos_t = COS(theta_t)

    m1 = SQRT(em1); m2 = SQRT(em2)

    rr = ABS(m2/m1)*cos_t/cos_i
    tv = rr*(ABS(2.0_fp*m1*cos_i/(m1*cos_t + m2*cos_i)))**2
    th = rr*(ABS(2.0_fp*m1*cos_i/(m1*cos_i + m2*cos_t)))**2

  END SUBROUTINE Fresnel_Transmittance


  SUBROUTINE Fresnel_Reflectance_Liou(theta_i,em1,em2,rv,rh)
    ! K,N,liou pp221
    REAL(fp),    INTENT(IN)    :: theta_i
    COMPLEX(fp), INTENT(IN)    :: em1,em2
    REAL(fp),    INTENT(OUT)   :: rv,rh
   
    REAL(fp)   :: k_real,k_img 
    COMPLEX(fp):: kz1,kz2,rvc,rhc
 
    CALL Dispersion(theta_i,em1,k_real,k_img)
    kz1=CMPLX(k_real,k_img,fp)
    CALL Dispersion(theta_i,em2,k_real,k_img)
    kz2=CMPLX(k_real,k_img,fp)
      
    ! wave field
    rhc = (kz1 - kz2)/(kz1 + kz2)
    rvc = (em2*kz1 - em1*kz2)/(em2*kz1 +em1*kz2)
    ! wave power
    rv = (ABS(rvc))**2
    rh = (ABS(rhc))**2
 
  END SUBROUTINE Fresnel_Reflectance_Liou
  
  SUBROUTINE Dispersion(theta_i,emc,k_real,k_img)

    REAL(fp),    INTENT(IN)    :: theta_i 
    COMPLEX(fp), INTENT(IN)    :: emc
    REAL(fp),    INTENT(OUT)   :: k_real,k_img

    REAL(fp) :: mr,mi,Nr,Ni,mn
     
    mr = DBLE(emc) ; mi = AIMAG(emc)
    mn = mr-(sin(theta_i))**2
    Nr = SQRT(0.5_fp*mn*(1.0_fp+SQRT(1.0_fp+mi**2/mn**2)))
    Ni = mi/(2.0_fp*Nr)
    k_real = Nr ; k_img = Ni
      
  END SUBROUTINE Dispersion

END MODULE CSEM_Fresnel



