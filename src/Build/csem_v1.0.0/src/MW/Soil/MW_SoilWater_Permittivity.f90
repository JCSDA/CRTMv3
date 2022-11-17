MODULE MW_SoilWater_Permittivity
! Purpose :
!   Calculate dielectric constant of water in soil media :

! Reference:
!  Dielectric constant of pure water
!   Ulaby p 2020
!  Dielectric constant of saline water
!   1) Stogryn, A. (1971): Equations for calculating the dielectric constant of
!    saline water, IEEE Transactions on Microwave Theory and Techniques,
!    Vol. MTT-19, 733-736.
!   2) Klein, L. A. and C. T. Swift (1977): An improved model
!    for the dielectric constant of saline water at microwave
!    frequencies, IEEE Transactions on  Antennas and Propagation,
!    Vol. AP-25, No. 1, 104-111.
!  Dielectric constant of soil water
!   1) Dobson '85. Modified Debye expression
!         Stern_Gouy double layer theory
!   2) Ulaby p 2024

! Interface :
!  medium = pure water(0) sea water(1) soil water(2)
!  isal = Stogryn (1) Klein and Swift (2)

! local variables :
!  N : normality from salinity (Stogryn, modified by Klein and Swift 1977)
!  T : temperature of water (C)
!  ew  : dielectric constant of water
!  sal : water salinity (psu = ppt(weight) )
!  eps_w0 : static dielectric constant of pure water (Klein and Swift 1977)
!  eps_sw0 : Static dielectric constant of soil water
!  tau_w : relaxation time of pure water (stogryn, 1970)
!  tau_sw : relaxation time of saline water
!  sigma : ionic conductivity
!  sigma_eff : effective conductivity of water (S/m)
!---------------------------------------------------------------------------

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  !USE Soil_Dielectric_Module
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MW_Pure_Water_Diel
  PUBLIC :: MW_Soil_Water_Diel
  
  ! set common constants
  REAL(fp),    PARAMETER ::  pi = 3.1415926_fp
  REAL(fp),    PARAMETER ::  eps_0 = 8.854e-12_fp   !   permittivity of free space (Klein and Swift 1977) [Farads/meter]
  REAL(fp),    PARAMETER ::  eps_winf  = 4.9_fp     !   eps_winf : dielectric constant at infinite frequency (Stogryn 1971)
  REAL(fp),    PARAMETER ::  rho_s = 2.66_fp        !   rho_s : soil specific density (g/cm3), ATBD 2.664
                                                    !   also called: high frequency dielectric constant (4.9 Lane and Saxton )
  COMPLEX(fp), PARAMETER :: j = (0._fp,1._fp)
CONTAINS


  SUBROUTINE MW_Pure_Water_Diel(freq, Tkelvin, ew)

    REAL(fp),    INTENT(IN)  :: freq
    REAL(fp),    INTENT(IN)  :: Tkelvin
    COMPLEX(fp), INTENT(OUT) :: ew

    REAL(fp) :: tau_w, eps_w0,sigma
    REAL(fp) :: f,omega
    REAL(fp) :: T
    COMPLEX(fp)  :: j = (0._fp,1._fp)

    f = freq*1.0E9_fp
    omega = 2._fp * pi * f
    T = Tkelvin- 273.15_fp 
    CALL Debye_Water_Parameters(T,tau_w,eps_w0,sigma)
    ew = eps_winf + (eps_w0 - eps_winf) / (1. - j * omega * tau_w) &
        + j * sigma / (omega * eps_0)

  END SUBROUTINE MW_Pure_Water_Diel


  SUBROUTINE MW_Soil_Water_Diel(freq, Tkelvin, wc, sal, sand, clay, ew)
    REAL(fp),    INTENT(IN)  :: freq
    REAL(fp),    INTENT(IN)  :: Tkelvin, wc, sal, sand, clay
    COMPLEX(fp), INTENT(OUT) :: ew
    ! 1 stogryn 2 Klein_Swift (just for 1.4Ghz)
    CHARACTER(LEN=20) :: DebyeAlg = 'Stogryn'

    REAL(fp) :: sigma_eff
    REAL(fp) :: tau_w, eps_w0
    REAL(fp) :: T, f, omega
  
    REAL(fp) :: rho_b
    
    ! Dobson model 1985
    f = freq*1.0E9_fp
    omega = 2._fp * pi * f
    T = Tkelvin- 273.15_fp
    ! bulk density frho_b= 1.4
    rho_b = sand*1.6_fp +  clay*1.1_fp + (1._fp-sand-clay)*1.2_fp
    
    CALL Debye_Water_Parameters(T,tau_w,eps_w0,sigma_eff,sal,DebyeAlg)
    ! Dobson Eq.32 2.013*sand CMEM 2.256_fp
    sigma_eff = -1.645_fp  + 1.939_fp * rho_b  - 2.256_fp * sand  + 1.594_fp * clay
    !  sigma_eff gets negative for very sandy soils
    !   with low bulk densities. If this happens set sigma_eff to zero.
    IF (sigma_eff < 0.) sigma_eff = 0._fp
  
    ! Modified Debye expression, Dobson '85
    ew = eps_winf + (eps_w0 - eps_winf) / (1._fp - j * omega * tau_w)  + &
         j * sigma_eff / (omega * eps_0) * (rho_s - rho_b) / (rho_s * MAX(0.001_fp, wc))

  END SUBROUTINE  MW_Soil_Water_Diel


!===========================================================================
! Private subroutines
!===========================================================================

  SUBROUTINE Debye_Water_Parameters(T,tau_w,eps_w0,sigma,sal,AlgTag)
    REAL(fp), INTENT(IN) :: T    
    REAL(fp), INTENT(IN), OPTIONAL  :: sal
    CHARACTER(LEN=*),  INTENT(IN),OPTIONAL  :: AlgTag
    REAL(fp), INTENT(OUT) :: eps_w0, tau_w, sigma
    CHARACTER(LEN=20)     :: AlgOpt = 'Stogryn'
   
    REAL(fp) :: a, bb, N
  
    tau_w = 1.768e-11_fp  - 6.068e-13_fp * T  + 1.104e-14_fp * T**2_fp  - 8.111e-17_fp * T**3_fp

    IF (.NOT. PRESENT(sal)) THEN
      eps_w0 = 88.045_fp- 0.4147_fp * T + 6.295e-4_fp* T**2_fp + 1.075e-5_fp * T**3_fp
      sigma = 0.0_fp
      RETURN
    ENDIF

    CALL ION_CONDUCT (T,sal,sigma)
    
    IF (PRESENT(AlgTag)) AlgOpt = AlgTag
    SELECT CASE (TRIM(AlgOpt))

      CASE ( 'Klein_Swift' ) !Klein, L. A. and C. T. Swift (1977)
        eps_w0 = 87.134_fp  - 1.949e-1_fp * T  - 1.276e-2_fp * T**2_fp  + 2.491e-4_fp * T**3_fp
        a = 1.000_fp  + 1.613e-5_fp * sal * T  - 3.656e-3_fp * sal  + 3.210e-5_fp * sal**2_fp  &
                       -  4.232e-7_fp * sal**3_fp
        eps_w0 = eps_w0 * a
        bb = 1.000_fp  + 2.282e-5_fp * sal * T  - 7.638e-4_fp * sal  - 7.760e-6_fp * sal**2_fp  &
                       + 1.105e-8_fp * sal**3_fp
        tau_w  = tau_w * bb

      CASE DEFAULT  !Stogryn
        N = 0.9141_fp * sal * (1.707e-2_fp + 1.205e-5_fp * sal + 4.058e-9_fp * sal**2_fp)
        eps_w0 = 87.74_fp  - 0.4008_fp * T  + 9.398e-4_fp * T**2_fp  + 1.410e-6_fp * T**3_fp
        a = 1.0_fp  - 0.2551_fp * N  + 5.151e-2_fp * N**2_fp  - 6.889e-3_fp * N**3_fp
        eps_w0 = eps_w0 * a
        bb = 1.0_fp - 0.04896_fp * N - 0.02967_fp * N**2_fp + 5.644e-3_fp &
             * N**3_fp + 0.1463e-2_fp * N * T
        tau_w  = tau_w * bb

    END SELECT


  END SUBROUTINE Debye_Water_Parameters


  SUBROUTINE ION_CONDUCT (T,sal,sigma)

! Purpose :
!  Calculate ionic conductivity for saline water
! Reference:
!  Stogryn, 1971

! local variables :
!  T : temperature of water (C)
!  sigma25 : ionic conductivity for sea water at 25 degree C (S/m)
!  sigma : ionic conductivity for saline water (S/m)
!---------------------------------------------------------------------------

    REAL(fp), INTENT(IN)  :: T
    REAL(fp), INTENT(IN)  :: sal
    REAL(fp), INTENT(OUT) :: sigma
    REAL(fp) :: alphac, delta, sigma25

    delta   =  25._fp - T
    alphac  =  2.033e-2_fp  +  1.266e-4_fp * delta  +  2.464e-6_fp * delta ** 2.0_fp -  &
               sal * (1.849e-5_fp  -  2.551e-7_fp * delta  +  2.551e-8_fp * delta ** 2.0_fp)
    sigma25 =  sal * (0.182521_fp  -  1.46192e-3_fp * sal  +  &
               2.09324e-5_fp * sal ** 2.0_fp  -  1.28205e-7_fp * sal ** 3.0_fp)

    sigma   =  sigma25 * exp(-1._fp * delta * alphac)

  END SUBROUTINE ION_CONDUCT
  
END MODULE MW_SoilWater_Permittivity
