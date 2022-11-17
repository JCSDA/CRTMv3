
!
! SnowMW_Optical_Model
!
! Module containing functions to simuate snow optics properties
! Currently, this model is not well established,
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, Fuzhong Weng, 02-28-2015
!                       Ming.Chen@noaa.gov
!                       Fuzhong.Weng@noaa.gov
!
MODULE SnowMW_Optical_Model
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Fresnel
  
  IMPLICIT NONE
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER, PRIVATE :: ZERO   = 0.0_fp
  REAL(fp), PARAMETER, PRIVATE :: ONE    = 1.0_fp
  REAL(fp), PARAMETER, PRIVATE :: TWO    = 2.0_fp
  REAL(fp), PARAMETER, PRIVATE :: PI     = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER, PRIVATE :: TWOPI  = TWO*PI
 
CONTAINS

!----------------------------------------------------------------------------------
!$$$  Snow_Diel
!                .      .    .                                       .
! subprogram:    Snow_Diel   compute dielectric constant of snow
!
! prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
!
! abstract: compute dielectric constant of snow
!
!
! program history log:
!
! input argument list:
!
!       frequency   -  frequency (ghz)
!       ep_real     -  real part of dielectric constant of particle
!       ep_imag     -  imaginary part of dielectric constant of particle
!       rad         -  particle radiu (mm)
!       frac        -  fraction volume of snow (0.0 - 1.0)
!
! output argument list:
!
!       ep_eff      -  dielectric constant of the dense medium
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!  Copyright (C) 2005 Fuzhong Weng and Banghua Yan
!
!  This program is free software; you can redistribute it and/or modify it under the terms of the GNU
!  General Public License as published by the Free Software Foundation; either version 2 of the License,
!  or (at your option) any later version.
!
!  This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even
!  the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public
!  License for more details.
!
!  You should have received a copy of the GNU General Public License along with this program; if not, write
!  to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
!----------------------------------------------------------------------------------

  SUBROUTINE Snow_Diel(frequency,ep_real,ep_imag,rad,frac,ep_eff)

    REAL(fp) :: ep_imag,ep_real
    REAL(fp) :: frequency,rad,frac,k0,yr,yi
    COMPLEX(fp) :: y,ep_r,ep_i,ep_eff,fracy

    k0 = TWOPI/(300.0_fp/frequency)

    yr = (ep_real - ONE)/(ep_real + TWO)
    yi = ep_imag/(ep_real + TWO)

    y = CMPLX(yr, yi, fp)
    fracy=frac*y

    ep_r = (ONE + TWO*fracy)/(ONE - fracy)
    ep_i = TWO*fracy*y*(k0*rad)**3*(ONE-frac)**4/((ONE-fracy)**2*(ONE+TWO*frac)**2)
    ep_eff = ep_r - CMPLX(ZERO,ONE,fp)*ep_i

    IF (AIMAG(ep_eff) >= ZERO) ep_eff = CMPLX(REAL(ep_eff), -0.0001_fp, fp)

  END SUBROUTINE Snow_Diel

!-------------------------------------------------------------------------------------------------------------
!$$$  Snow_Optic
!                .      .    .                                       .
! subprogram:    landem      comput optic parameters for snow
!
! prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
!
! abstract: compute optic parameters for snow
!
! program history log:
!
! input argument list:
!
!      theta        -  local zenith angle (degree)
!      frequency    -  frequency (ghz)
!      ep_real      -  real part of dielectric constant of particles
!      ep_imag      -  imaginary part of dielectric constant of particles
!      a            -  particle radiu (mm)
!      h            -  snow depth(mm)
!      f            -  fraction volume of snow (0.0 - 1.0)
!
! output argument list:
!
!       ssalb       -  single scattering albedo
!       tau         -  optical depth
!       g           -  asymmetry factor
!
!   important internal variables:
!
!       ks          -  scattering coeffcient (/mm)
!       ka          -  absorption coeffient (/mm)
!       kp          -  eigenvalue of two-stream approximation
!       y           -  = yr+iyi
!
! remarks:
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  SUBROUTINE Snow_Optic(frequency,a,h,f,ep_real,ep_imag,gv,gh, ssalb_v,ssalb_h,tau_v,tau_h)

    REAL(fp) :: yr,yi,ep_real,ep_imag
    REAL(fp) :: frequency,a,h,f,ssalb_v,ssalb_h,tau_v,tau_h,gv,gh,k
    REAL(fp) :: ks1,ks2,ks3,ks,kr1,kr2,kr3,kr,ki1,ki2,ki3,ki
    REAL(fp) :: fact1,fact2,fact3,fact4,fact5

    k = TWOPI/(300._fp/frequency)

    yr = (ep_real - ONE)/(ep_real + TWO)
    yi = -ep_imag/(ep_real + TWO)

    fact1 = (ONE+TWO*f)**2
    fact2 = ONE-f*yr
    fact3 = (ONE-f)**4
    fact4 = f*(k*a)**3
    fact5 = ONE+TWO*f*yr

    ks1 = k*SQRT(fact2/fact5)
    ks2 = fact4*fact3/fact1
    ks3 = (yr/fact2)**2
    ks = ks1*ks2*ks3

    kr1 = fact5/fact2
    kr2 = TWO*ks2
    kr3 = TWO*yi*yr/(fact2**3)
    kr = k*SQRT(kr1+kr2*kr3)

    ki1 = 3.0_fp*f*yi/fact2**2
    ki2 = kr2
    ki3 = ks3
    ki  = k**2/(TWO*kr)*(ki1+ki2*ki3)

    gv = 0.5_fp
    gh = 0.5_fp

    ssalb_v = MIN(ks/ki, 0.999_fp)
    ssalb_h = ssalb_v
    tau_v = TWO*ki*h
    tau_h = tau_v

  END SUBROUTINE Snow_Optic


END MODULE SnowMW_Optical_Model



