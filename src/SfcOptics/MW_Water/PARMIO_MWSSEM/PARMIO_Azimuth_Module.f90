!
! PARMIO_Azimuth_Module
!
! Recombine PARMIO's azimuthal-harmonic emissivity coefficients into the
! CRTM microwave surface-optics emissivity basis:
!   (V-pol, H-pol, U, circular/Stokes-V).
! The first two entries are the decoupled V/H slots expected by
! CRTM_SfcOptics and FASTEM, not canonical Stokes I/Q.
!
! PARMIO stores per (freq, theta, U10, sst, sss) cell a 14-vector of
! dimensionless coefficients = Tb / SST_K:
!     1: evN     (V-pol specular)            }-- azimuth-independent
!     2: ehN     (H-pol specular)            /
!     3: ev0     (V-pol roughness 0th)       }-- azimuth-independent
!     4: eh0     (H-pol roughness 0th)       /
!     5: ev1     (V-pol cos(phi))            }-- 1st harmonic (cos in V/H,
!     6: eh1     (H-pol cos(phi))            /                  sin in U/V_S)
!     7: eU1     (3rd Stokes sin(phi))
!     8: eV1     (4th Stokes sin(phi))
!     9: ev2     (V-pol cos(2 phi))          }-- 2nd harmonic (same convention)
!    10: eh2     (H-pol cos(2 phi))          /
!    11: eU2     (3rd Stokes sin(2 phi))
!    12: eV2     (4th Stokes sin(2 phi))
!    13: edv_MR  (V-pol multi-reflection)    }-- azimuth-independent
!    14: edh_MR  (H-pol multi-reflection)    /
!
! At runtime:
!   e_V(phi) = evN + ev0 + ev1 cos(phi)   + ev2 cos(2 phi) + edv_MR
!   e_H(phi) = ehN + eh0 + eh1 cos(phi)   + eh2 cos(2 phi) + edh_MR
!   e_U(phi) =       eU1 sin(phi)         + eU2 sin(2 phi)
!   e_V_S(phi) =     eV1 sin(phi)         + eV2 sin(2 phi)
!
! Phi is the relative azimuth between the wind direction and the sensor
! viewing direction.

MODULE PARMIO_Azimuth_Module

  USE Type_Kinds, ONLY: fp
  USE PARMIOCoeff_Define, ONLY: N_PARMIO_HARMONIC_TERMS
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: PARMIO_Azimuth_Recombine
  PUBLIC :: PARMIO_Azimuth_Recombine_TL
  PUBLIC :: PARMIO_Azimuth_Recombine_AD

  REAL(fp), PARAMETER :: PI = 3.141592653589793238462643383279_fp
  REAL(fp), PARAMETER :: DEGREES_TO_RADIANS = PI / 180.0_fp

CONTAINS

  !-----------------------------------------------------------------
  !  PARMIO_Azimuth_Recombine
  !  Combine the 14 harmonic coefficients into CRTM's microwave
  !  surface-optics basis (V-pol, H-pol, U, circular/Stokes-V) at the
  !  requested relative azimuth.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_Azimuth_Recombine( &
      Coefficients, Azimuth_Angle_deg, Emissivity)
    REAL(fp), INTENT(IN)  :: Coefficients(N_PARMIO_HARMONIC_TERMS)
    REAL(fp), INTENT(IN)  :: Azimuth_Angle_deg
    REAL(fp), INTENT(OUT) :: Emissivity(4)
    REAL(fp) :: phi, c1, s1, c2, s2

    phi = Azimuth_Angle_deg * DEGREES_TO_RADIANS
    c1 = COS(phi)
    s1 = SIN(phi)
    c2 = COS(2.0_fp * phi)
    s2 = SIN(2.0_fp * phi)

    ! V-pol = evN + ev0 + ev1 c1 + ev2 c2 + edv_MR
    Emissivity(1) = Coefficients(1) + Coefficients(3)              &
                  + Coefficients(5)  * c1                          &
                  + Coefficients(9)  * c2                          &
                  + Coefficients(13)
    ! H-pol = ehN + eh0 + eh1 c1 + eh2 c2 + edh_MR
    Emissivity(2) = Coefficients(2) + Coefficients(4)              &
                  + Coefficients(6)  * c1                          &
                  + Coefficients(10) * c2                          &
                  + Coefficients(14)
    ! 3rd Stokes (U) = eU1 s1 + eU2 s2
    Emissivity(3) = Coefficients(7)  * s1 + Coefficients(11) * s2
    ! 4th Stokes (V_Stokes) = eV1 s1 + eV2 s2
    Emissivity(4) = Coefficients(8)  * s1 + Coefficients(12) * s2
  END SUBROUTINE PARMIO_Azimuth_Recombine


  !-----------------------------------------------------------------
  !  PARMIO_Azimuth_Recombine_TL
  !  Tangent-linear recombination of the 14 harmonic coefficients.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_Azimuth_Recombine_TL( &
      Coefficients, Coefficients_TL, Azimuth_Angle_deg, &
      Azimuth_Angle_deg_TL, Emissivity_TL)
    REAL(fp), INTENT(IN)  :: Coefficients(N_PARMIO_HARMONIC_TERMS)
    REAL(fp), INTENT(IN)  :: Coefficients_TL(N_PARMIO_HARMONIC_TERMS)
    REAL(fp), INTENT(IN)  :: Azimuth_Angle_deg
    REAL(fp), INTENT(IN)  :: Azimuth_Angle_deg_TL
    REAL(fp), INTENT(OUT) :: Emissivity_TL(4)
    REAL(fp) :: phi, phi_TL
    REAL(fp) :: c1, s1, c2, s2
    REAL(fp) :: c1_TL, s1_TL, c2_TL, s2_TL

    phi    = Azimuth_Angle_deg    * DEGREES_TO_RADIANS
    phi_TL = Azimuth_Angle_deg_TL * DEGREES_TO_RADIANS
    c1 = COS(phi)
    s1 = SIN(phi)
    c2 = COS(2.0_fp * phi)
    s2 = SIN(2.0_fp * phi)
    c1_TL = -s1 * phi_TL
    s1_TL =  c1 * phi_TL
    c2_TL = -2.0_fp * s2 * phi_TL
    s2_TL =  2.0_fp * c2 * phi_TL

    Emissivity_TL(1) = Coefficients_TL(1) + Coefficients_TL(3)     &
                     + Coefficients_TL(5)  * c1                    &
                     + Coefficients(5)     * c1_TL                 &
                     + Coefficients_TL(9)  * c2                    &
                     + Coefficients(9)     * c2_TL                 &
                     + Coefficients_TL(13)
    Emissivity_TL(2) = Coefficients_TL(2) + Coefficients_TL(4)     &
                     + Coefficients_TL(6)  * c1                    &
                     + Coefficients(6)     * c1_TL                 &
                     + Coefficients_TL(10) * c2                    &
                     + Coefficients(10)    * c2_TL                 &
                     + Coefficients_TL(14)
    Emissivity_TL(3) = Coefficients_TL(7)  * s1                    &
                     + Coefficients(7)     * s1_TL                 &
                     + Coefficients_TL(11) * s2                    &
                     + Coefficients(11)    * s2_TL
    Emissivity_TL(4) = Coefficients_TL(8)  * s1                    &
                     + Coefficients(8)     * s1_TL                 &
                     + Coefficients_TL(12) * s2                    &
                     + Coefficients(12)    * s2_TL
  END SUBROUTINE PARMIO_Azimuth_Recombine_TL


  !-----------------------------------------------------------------
  !  PARMIO_Azimuth_Recombine_AD
  !  Adjoint recombination of the 14 harmonic coefficients.
  !-----------------------------------------------------------------
  SUBROUTINE PARMIO_Azimuth_Recombine_AD( &
      Coefficients, Emissivity_AD, Azimuth_Angle_deg, &
      Coefficients_AD, Azimuth_Angle_deg_AD)
    REAL(fp), INTENT(IN)     :: Coefficients(N_PARMIO_HARMONIC_TERMS)
    REAL(fp), INTENT(IN OUT) :: Emissivity_AD(4)
    REAL(fp), INTENT(IN)     :: Azimuth_Angle_deg
    REAL(fp), INTENT(IN OUT) :: Coefficients_AD(N_PARMIO_HARMONIC_TERMS)
    REAL(fp), INTENT(IN OUT) :: Azimuth_Angle_deg_AD
    REAL(fp) :: phi
    REAL(fp) :: c1, s1, c2, s2
    REAL(fp) :: c1_AD, s1_AD, c2_AD, s2_AD, phi_AD

    phi = Azimuth_Angle_deg * DEGREES_TO_RADIANS
    c1 = COS(phi)
    s1 = SIN(phi)
    c2 = COS(2.0_fp * phi)
    s2 = SIN(2.0_fp * phi)

    c1_AD  = 0.0_fp
    s1_AD  = 0.0_fp
    c2_AD  = 0.0_fp
    s2_AD  = 0.0_fp
    phi_AD = 0.0_fp

    Coefficients_AD(8)  = Coefficients_AD(8)  + s1 * Emissivity_AD(4)
    s1_AD              = s1_AD               + Coefficients(8)  * Emissivity_AD(4)
    Coefficients_AD(12) = Coefficients_AD(12) + s2 * Emissivity_AD(4)
    s2_AD              = s2_AD               + Coefficients(12) * Emissivity_AD(4)
    Emissivity_AD(4) = 0.0_fp

    Coefficients_AD(7)  = Coefficients_AD(7)  + s1 * Emissivity_AD(3)
    s1_AD              = s1_AD               + Coefficients(7)  * Emissivity_AD(3)
    Coefficients_AD(11) = Coefficients_AD(11) + s2 * Emissivity_AD(3)
    s2_AD              = s2_AD               + Coefficients(11) * Emissivity_AD(3)
    Emissivity_AD(3) = 0.0_fp

    Coefficients_AD(2)  = Coefficients_AD(2)  + Emissivity_AD(2)
    Coefficients_AD(4)  = Coefficients_AD(4)  + Emissivity_AD(2)
    Coefficients_AD(6)  = Coefficients_AD(6)  + c1 * Emissivity_AD(2)
    c1_AD              = c1_AD               + Coefficients(6)  * Emissivity_AD(2)
    Coefficients_AD(10) = Coefficients_AD(10) + c2 * Emissivity_AD(2)
    c2_AD              = c2_AD               + Coefficients(10) * Emissivity_AD(2)
    Coefficients_AD(14) = Coefficients_AD(14) + Emissivity_AD(2)
    Emissivity_AD(2) = 0.0_fp

    Coefficients_AD(1)  = Coefficients_AD(1)  + Emissivity_AD(1)
    Coefficients_AD(3)  = Coefficients_AD(3)  + Emissivity_AD(1)
    Coefficients_AD(5)  = Coefficients_AD(5)  + c1 * Emissivity_AD(1)
    c1_AD              = c1_AD               + Coefficients(5)  * Emissivity_AD(1)
    Coefficients_AD(9)  = Coefficients_AD(9)  + c2 * Emissivity_AD(1)
    c2_AD              = c2_AD               + Coefficients(9)  * Emissivity_AD(1)
    Coefficients_AD(13) = Coefficients_AD(13) + Emissivity_AD(1)
    Emissivity_AD(1) = 0.0_fp

    phi_AD = phi_AD - s1 * c1_AD
    phi_AD = phi_AD + c1 * s1_AD
    phi_AD = phi_AD - 2.0_fp * s2 * c2_AD
    phi_AD = phi_AD + 2.0_fp * c2 * s2_AD
    Azimuth_Angle_deg_AD = Azimuth_Angle_deg_AD + phi_AD * DEGREES_TO_RADIANS
  END SUBROUTINE PARMIO_Azimuth_Recombine_AD

END MODULE PARMIO_Azimuth_Module
