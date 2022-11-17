MODULE RTTOV_FASTEM5R1_TL_MODULE

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  ! Disable implicit typing
  IMPLICIT NONE

CONTAINS

! Description:
!> @file
!!   TL of FASTEM-4,5,6 emissivity and reflectance calculation
!
!> @brief
!!   TL of FASTEM-4,5,6 emissivity and reflectance calculation
!!
!!
!! @param[in]     fastem_version        FASTEM version to compute (4, 5 or 6)
!! @param[in]     frequency             channel frequency (GHz)
!! @param[in]     zenith_angle          profile zenith angle (degrees)
!! @param[in]     temperature           profile skin temperature (K)
!! @param[in]     salinity              profile salinity (practical salinity units)
!! @param[in]     wind_speed            profile wind speed (m/s)
!! @param[in,out] emissivity_tl         emissivity perturbation (4 Stokes components)
!! @param[in,out] reflectivity_tl       reflectivity perturbation (4 Stokes components)
!! @param[in]     temperature_tl        profile skin temperature perturbation
!! @param[in]     salinity_tl           profile salinity perturbation
!! @param[in]     wind_speed_tl         profile wind speed perturbation
!! @param[out]    emissivity            calculated emissivity (4 Stokes components)
!! @param[out]    reflectivity          calculated reflectivity (4 Stokes components)
!! @param[in]     transmittance         surface-to-space transmittance
!! @param[in]     rel_azimuth           relative azimuth angle
!! @param[in]     transmittance_tl      surface-to-space transmittance perturbation
!! @param[in]     rel_azimuth_tl        relative azimuth angle perturbation
!! @param[in]     supply_foam_fraction  flag to indicate user is supplying foam fraction, optional
!! @param[in]     foam_fraction         user supplied foam fraction, optional
!! @param[in]     foam_fraction_tl      user foam fraction perturbation, optional
!
! Copyright:
!    This software was developed within the context of
!    the EUMETSAT Satellite Application Facility on
!    Numerical Weather Prediction (NWP SAF), under the
!    Cooperation Agreement dated 25 November 1998, between
!    EUMETSAT and the Met Office, UK, by one or more partners
!    within the NWP SAF. The partners in the NWP SAF are
!    the Met Office, ECMWF, KNMI and MeteoFrance.
!
!    Copyright 2015, EUMETSAT, All Rights Reserved.
!

  SUBROUTINE rttov_fastem5r1_tl(fastem_version,       &  ! Input
                              Frequency   ,         &  ! Input
                              Zenith_Angle,         &  ! Input
                              Temperature ,         &  ! Input
                              Salinity    ,         &  ! Input
                              Wind_Speed  ,         &  ! Input
                              Temperature_tl,       &  ! Input
                              Salinity_tl ,         &  ! Input
                              Wind_Speed_tl,        &  ! Input
                              Emissivity  ,         &  ! Output
                              Reflectivity,         &  ! Output
                              Emissivity_tl,        &  ! Output
                              Reflectivity_tl,      &  ! Output
                              Transmittance,        &  ! Input, may not be used
                              Rel_Azimuth,          &  ! Input, may not be used
                              Transmittance_tl,     &  ! Input, may not be used
                              Rel_Azimuth_tl,       &  ! Input, may not be used
                              Supply_Foam_Fraction, &  ! Optional input
                              Foam_Fraction,        &  ! Optional input
                              Foam_Fraction_tl)        ! Optional input

!INTF_OFF
    USE mod_rttov_fastem5r1_coef, ONLY : FresnelVariables_type, PermittivityVariables_type,&
        ZERO, POINT_5, ONE, TWO, THREE,PI,DEGREES_TO_RADIANS,transmittance_limit_lower,&
        transmittance_limit_upper, e0_4, e0_5, min_f, max_f, min_wind, max_wind, A_COEF, Lcoef4, Lcoef5,&
        Scoef, t_c4, t_c5, b_coef, FR_COEFF, x, y, coef_mk_azi
!INTF_ON
!    USE mod_rttov_fastem5_coef, ONLY: fp
!    USE parkind1, ONLY : jpim, jplm
    IMPLICIT NONE
    ! Arguments
!   INTEGER(jpim),   INTENT(IN)     :: fastem_version
    INTEGER,         INTENT(IN)     :: fastem_version
    REAL(fp),        INTENT(IN)     :: Frequency
    REAL(fp),        INTENT(IN)     :: Zenith_Angle
    REAL(fp),        INTENT(IN)     :: Temperature
    REAL(fp),        INTENT(IN)     :: Salinity
    REAL(fp),        INTENT(IN)     :: Wind_Speed
    REAL(fp),        INTENT(OUT)    :: Emissivity(4), Reflectivity(4)
    REAL(fp), INTENT(IN)  :: Transmittance
    REAL(fp), INTENT(IN)  :: Rel_Azimuth
  ! TL part
    REAL(fp),                    INTENT(IN)     :: Temperature_tl
    REAL(fp),                    INTENT(IN)     :: Salinity_tl
    REAL(fp),                    INTENT(IN)     :: Wind_Speed_tl
    REAL(fp),        INTENT(INOUT)    :: Emissivity_tl(4), Reflectivity_tl(4)
    REAL(fp), OPTIONAL, INTENT(IN)  :: Transmittance_tl
    REAL(fp), OPTIONAL, INTENT(IN)  :: Rel_Azimuth_tl
!   LOGICAL(jplm),  OPTIONAL, INTENT(IN)  :: Supply_Foam_Fraction
    LOGICAL,        OPTIONAL, INTENT(IN)  :: Supply_Foam_Fraction
    REAL(fp),       OPTIONAL, INTENT(IN)  :: Foam_Fraction
    REAL(fp),       OPTIONAL, INTENT(IN)  :: Foam_Fraction_tl

!INTF_END

  !local variable

    REAL(fp) :: e0
    REAL(fp) :: cos_z, Foam_Cover
    REAL(fp) :: scor, small_corr, Azimuth_Emi(4),RV_Fresnel,RH_Fresnel
    REAL(fp) :: Ev,Eh,RvL,RhL,RvS,RhS
    REAL(fp) :: zreflmod_v,zreflmod_h,zrough_v,zrough_h
    INTEGER :: i, j, L, m
    LOGICAL  :: lcalc_foam_fraction = .TRUE.
    REAL(fp) :: Foam_Rv,Foam_Rh

    ! Local variables for the permittivity model
    REAL(fp) :: einf, sigma25
    REAL(fp) :: tau1, tau2, es, e1
    REAL(fp) :: perm_Real, perm_imag
!    REAL(fp) ::
    TYPE(PermittivityVariables_type) :: iVar
    COMPLEX( fp ) :: Permittivity

    ! Local variables for Fresnel reflectance
    COMPLEX(fp) :: zRv ! Vertical
    COMPLEX(fp) :: zRh ! Horizontal
    TYPE(FresnelVariables_type) :: frVar

    ! Local variables for small-scale
    REAL(fp) :: windspeed, freq_S
    ! Local variables for large-scale
    REAL(fp) :: seczen, zc(12)
    ! Local variables for including transmittance
    REAL(fp) :: variance,varm,opdpsfc,zx(9)
    ! Foam reflectance
    REAL(fp) :: Fh, Foam_ref
    ! Local variables for azimuth angle
    REAL(fp) :: ac, sc, fre_c, phi, wind10
    ! Local arrays to hold FASTEM-4/5 coefs
    REAL(fp) :: Lcoef(size(Lcoef5)), t_c(size(t_c5))

! TL internal variables
    REAL(fp) :: sigma_TL
    REAL(fp) :: einf_TL, e1_TL, es_TL, ce1_TL, ces_TL
    REAL(fp) :: t_TL, t_sq_TL, t_cu_TL, S_TL, beta_TL, delta_TL, sigma25_TL
    REAL(fp) :: tau1_TL, tau2_TL, ctau1_TL, ctau2_TL, f1_TL, f2_TL, del1_TL, del2_TL
    REAL(fp) :: perm_Real_TL, perm_imag_TL
    COMPLEX(fp) :: Permittivity_TL

    COMPLEX(fp) :: z1_TL, z2_TL
    COMPLEX(fp) :: zRv_TL           ! Vertical
    COMPLEX(fp) :: zRh_TL           ! Horizontal
    REAL(fp)    :: rzRv_TL,izRv_TL  ! Vertical
    REAL(fp)    :: rzRh_TL,izRh_TL  ! Horizontal
    REAL(fp) :: zreflmod_v_tl,zreflmod_h_tl, windspeed_tl
    REAL(fp) :: scor_tl, small_corr_tl, Azimuth_Emi_tl(4)
    REAL(fp) :: Foam_Cover_tl,RV_Fresnel_tl,RH_Fresnel_tl
    REAL(fp) :: Ev_tl,Eh_tl,RvL_tl,RhL_tl,RvS_tl,RhS_tl
    REAL(fp) :: variance_tl,varm_tl,opdpsfc_tl,zx_tl(9),zrough_v_tl,zrough_h_tl
    REAL(fp) :: ac_tl, sc_tl, wind10_tl, phi_tl

    INTEGER  :: ifreq,ipol
    REAL(fp) :: azimuth_component(2,6) ! pol,  freq
    REAL(fp) :: azimuth_component_tl(2,6)
    REAL(fp),parameter :: xs11=2
    REAL(fp),parameter :: xs12=2
    REAL(fp),parameter :: xs21=1
    REAL(fp),parameter :: xs22=4
    REAL(fp),parameter :: theta_ref=55.2d0
    REAL(fp) :: theta

    REAL(fp),dimension(6) :: A1v,A1h,A2v,A2h
    REAL(fp),dimension(6) :: A1s1,A1s2,A2s1,A2s2
    REAL(fp),dimension(6) :: A2s2_theta0, A1s1_theta,A2s1_theta,A1s2_theta,A2s2_theta
    REAL(fp),dimension(6) :: A1v_theta, A1h_theta,A2v_theta,A2h_theta
    REAL(fp),dimension(6) :: A1v_tl,A1h_tl,A2v_tl,A2h_tl
    REAL(fp),dimension(6) :: A1s1_tl,A1s2_tl,A2s1_tl,A2s2_tl
    REAL(fp),dimension(6) :: A2s2_theta0_tl, A1s1_theta_tl,A2s1_theta_tl,A1s2_theta_tl,A2s2_theta_tl
    REAL(fp),dimension(6) :: A1v_theta_tl, A1h_theta_tl,A2v_theta_tl,A2h_theta_tl
    REAL(fp) :: fratio

    IF (fastem_version == 4) THEN
      e0 = e0_4
      Lcoef = Lcoef4
      t_c = t_c4
    ELSE
      e0 = e0_5
      Lcoef = Lcoef5
      t_c = t_c5
    ENDIF
    cos_z = cos( Zenith_Angle*DEGREES_TO_RADIANS )

  ! Permittivity Calculation
  ! ------------------------
    !1.2 calculate permittivity using double-debye formula
    !-----------------------------------------------------
    !Set values for temperature polynomials (convert from kelvin to celsius)
    iVar%t = Temperature - 273.15_fp
    iVar%t_sq = iVar%t * iVar%t     !quadratic
    iVar%t_cu = iVar%t_sq * iVar%t  !cubic
    iVar%S = Salinity
    !-----------------------------------------------------
    !1.2 Pure or fresh water
    !-----------------------------------------------------
    einf = A_COEF(0) + A_COEF(1)*iVar%t
    es   = A_COEF(2) + A_COEF(3)*iVar%t  + A_COEF(4)*iVar%t_sq + A_COEF(5)*iVar%t_cu
    e1   = A_COEF(9) + A_COEF(10)*iVar%t + A_COEF(11)*iVar%t_sq
    tau1 = A_COEF(15) + A_COEF(16)*iVar%t + A_COEF(17)*iVar%t_sq + A_COEF(18)*iVar%t_cu
    tau2 = A_COEF(22) + A_COEF(23)*iVar%t + A_COEF(24)*iVar%t_sq + A_COEF(25)*iVar%t_cu

    iVar%es_k = es
    iVar%e1_k = e1
    iVar%tau1_k = tau1
    iVar%tau2_k = tau2
    perm_imag = ZERO

    IF( iVar%S > ZERO ) THEN
      iVar%delta = 25.0_fp - iVar%t
      iVar%beta  = A_COEF(29) +A_COEF(30)*iVar%delta +A_COEF(31)*iVar%delta**2  &
            + iVar%S*(A_COEF(32) +A_COEF(33)*iVar%delta +A_COEF(34)*iVar%delta**2)
      sigma25 = iVar%S*(A_COEF(35) +A_COEF(36)*iVar%S +A_COEF(37)*iVar%S**2  &
              +A_COEF(38)*iVar%S**3)
      iVar%sigma = sigma25*exp(-iVar%delta*iVar%beta)

      iVar%ces = ONE + iVar%S*(A_COEF(6) + A_COEF(7)*iVar%S + A_COEF(8)*iVar%t )
      iVar%ce1 = ONE + iVar%S*(A_COEF(12) + A_COEF(13)*iVar%S +A_COEF(14)*iVar%t )
      iVar%ctau1 = ONE + iVar%S*(A_COEF(19) +A_COEF(20)*iVar%t + A_COEF(21)*iVar%t_sq)
      iVar%ctau2 = ONE + iVar%S*(A_COEF(26) + A_COEF(27)*iVar%t + A_COEF(28)*iVar%S**2 )
      es = iVar%es_k * iVar%ces
      e1 = iVar%e1_k * iVar%ce1
      tau1 = iVar%tau1_k * iVar%ctau1
      tau2 = iVar%tau2_k * iVar%ctau2
      perm_imag = -iVar%sigma/(TWO*PI*e0*Frequency)
    END IF
    !Define two relaxation frequencies, f1 and f2
    iVar%f1 = Frequency*tau1
    iVar%f2 = Frequency*tau2
    iVar%del1 = es - e1
    iVar%del2 = e1 - einf

  ! TL
    t_TL = Temperature_TL
    S_TL = Salinity_TL
    t_sq_TL = TWO * iVar%t * t_TL     !quadratic
    t_cu_TL = THREE * iVar%t_sq * t_TL  !cubic
    !-----------------------------------------------------
    !1.2 Pure or fresh water
    !-----------------------------------------------------
    einf_TL = A_COEF(1)*t_TL
    es_TL   = A_COEF(3)*t_TL  + A_COEF(4)*t_sq_TL + A_COEF(5)*t_cu_TL
    e1_TL   = A_COEF(10)*t_TL + A_COEF(11)*t_sq_TL
    tau1_TL = A_COEF(16)*t_TL + A_COEF(17)*t_sq_TL + A_COEF(18)*t_cu_TL
    tau2_TL = A_COEF(23)*t_TL + A_COEF(24)*t_sq_TL + A_COEF(25)*t_cu_TL
    perm_imag_TL = ZERO

    IF( iVar%S > ZERO ) THEN
      delta_TL = - t_TL
      beta_TL  = A_COEF(30)*delta_TL +TWO*A_COEF(31)*iVar%delta*delta_TL &
            + S_TL*(A_COEF(32) +A_COEF(33)*iVar%delta +A_COEF(34)*iVar%delta**2)  &
            + iVar%S*(A_COEF(33)*delta_TL +TWO*A_COEF(34)*iVar%delta*delta_TL)
      sigma25_TL = S_TL*(A_COEF(35) +A_COEF(36)*iVar%S +A_COEF(37)*iVar%S**2   &
                 +A_COEF(38)*iVar%S**3) + iVar%S*(A_COEF(36)*S_TL &
                 +TWO*A_COEF(37)*iVar%S*S_TL +THREE*A_COEF(38)*iVar%S**2*S_TL)
      sigma_TL = sigma25_TL*exp(-iVar%delta*iVar%beta)  &
               -(delta_TL*iVar%beta+iVar%delta*beta_TL)*iVar%sigma
      ces_TL = A_COEF(6)*S_TL + TWO*A_COEF(7)*iVar%S*S_TL  &
             + A_COEF(8)*S_TL*iVar%t +A_COEF(8)*iVar%S*t_TL
      ce1_TL = A_COEF(12)*S_TL + TWO*A_COEF(13)*iVar%S*S_TL  &
             +A_COEF(14)*S_TL*iVar%t +A_COEF(14)*iVar%S*t_TL
      ctau1_TL = S_TL*(A_COEF(19) +A_COEF(20)*iVar%t + A_COEF(21)*iVar%t_sq)  &
               + iVar%S*(A_COEF(20)*t_TL + A_COEF(21)*t_sq_TL)
      ctau2_TL = S_TL*(A_COEF(26) + A_COEF(27)*iVar%t + A_COEF(28)*iVar%S**2) &
               + iVar%S*(A_COEF(27)*t_TL + TWO*A_COEF(28)*iVar%S*S_TL)
      es_TL = es_TL * iVar%ces + iVar%es_k * ces_TL
      e1_TL = e1_TL * iVar%ce1 + iVar%e1_k * ce1_TL
      tau1_TL = tau1_TL * iVar%ctau1 + iVar%tau1_k * ctau1_TL
      tau2_TL = tau2_TL * iVar%ctau2 + iVar%tau2_k * ctau2_TL
      perm_imag_TL = -sigma_TL/(TWO*PI*e0*Frequency)
    END IF
    !Define two relaxation frequencies, f1 and f2
    f1_TL = Frequency*tau1_TL
    f2_TL = Frequency*tau2_TL
    del1_TL = es_TL - e1_TL
    del2_TL = e1_TL - einf_TL

    perm_Real = einf + iVar%del1/(ONE + iVar%f1**2) + iVar%del2/(ONE + iVar%f2**2)
    perm_imag = -perm_imag + iVar%del1*iVar%f1/(ONE + iVar%f1**2)  &
              + iVar%del2*iVar%f2/(ONE + iVar%f2**2)
    Permittivity = Cmplx(perm_Real,-perm_imag,fp)

    perm_Real_TL = einf_TL + del1_TL/(ONE + iVar%f1**2)  &
                 -TWO*iVar%del1*iVar%f1*f1_TL/(ONE + iVar%f1**2)**2 &
    + del2_TL/(ONE + iVar%f2**2) -TWO*iVar%del2*iVar%f2*f2_TL/(ONE + iVar%f2**2)**2

    perm_imag_TL = -perm_imag_TL + (del1_TL*iVar%f1+iVar%del1*f1_TL)/(ONE + iVar%f1**2) &
      - TWO*iVar%del1*iVar%f1**2*f1_TL/(ONE + iVar%f1**2)**2  &
      + (del2_TL*iVar%f2+iVar%del2*f2_TL)/(ONE + iVar%f2**2)  &
      - TWO*iVar%del2*iVar%f2**2*f2_TL/(ONE + iVar%f2**2)**2
    Permittivity_TL = Cmplx(perm_Real_TL,-perm_imag_TL,fp)
  
  ! Compute Fresnel reflectance code, adopted from Masahiro Kazumori, JMA
  ! -------------------------------------------
    ! Compute the complex reflectivity components
    frVar%z1 = SQRT(permittivity - ONE + (cos_z*cos_z))
    frVar%z2 = permittivity * cos_z
    zRh = (cos_z  -frVar%z1) / (cos_z  +frVar%z1)
    zRv = (frVar%z2-frVar%z1) / (frVar%z2+frVar%z1)

    ! The square of the vertical abs value
    frVar%rzRv = REAL(zRv,fp)
    frVar%izRv = AIMAG(zRv)
    Rv_Fresnel = frVar%rzRv**2 + frVar%izRv**2

    ! The square of the horizontal abs value
    frVar%rzRh = REAL(zRh,fp)
    frVar%izRh = AIMAG(zRh)
    Rh_Fresnel = frVar%rzRh**2 + frVar%izRh**2

    ! Compute the tangent-linear complex reflectivity components
    z1_TL = POINT_5 * permittivity_TL / frVar%z1
    z2_TL = cos_z * permittivity_TL
    zRh_TL = -TWO * cos_z * z1_TL / (cos_z+frVar%z1)**2
    zRv_TL =  TWO * (frVar%z1*z2_TL - frVar%z2*z1_TL) / (frVar%z2+frVar%z1)**2

    ! The square of the tangent-linear vertical abs value
    rzRv_TL = REAL(zRv_TL,fp)
    izRv_TL = AIMAG(zRv_TL)
    Rv_Fresnel_tl = TWO * (frVar%rzRv*rzRv_TL + frVar%izRv*izRv_TL)

    ! The square of the tangent-linear horizontal abs value
    rzRh_TL = REAL(zRh_TL,fp)
    izRh_TL = AIMAG(zRh_TL)
    Rh_Fresnel_tl = TWO * (frVar%rzRh*rzRh_TL + frVar%izRh*izRh_TL)


  ! Apply small-scale correction
  ! --------------------------------
    windspeed = Wind_Speed
    IF( windspeed < min_wind ) windspeed = min_wind
    IF( windspeed > max_wind ) windspeed = max_wind

    freq_S = Frequency
    IF( freq_S < min_f ) freq_S = min_f
    IF( freq_S > max_f ) freq_S = max_f

    scor = Scoef(1) *windspeed*freq_S +Scoef(2) *windspeed*freq_S**2 &
           + Scoef(3) *windspeed**2* freq_S +Scoef(4) *windspeed**2* freq_S**2 &
           + Scoef(5) *windspeed**2 /freq_S +Scoef(6) *windspeed**2 /freq_S**2 &
           + Scoef(7) *windspeed + Scoef(8) *windspeed**2

    IF( Wind_Speed < min_wind .OR. Wind_Speed > max_wind ) THEN
      windspeed_tl = 0._fp
    ELSE
      windspeed_tl = wind_speed_tl
    ENDIF
    scor_tl = ( Scoef(1)*freq_S +Scoef(2)*freq_S**2 &
            +TWO*Scoef(3) *windspeed*freq_S +TWO*Scoef(4)*windspeed* freq_S**2 &
            +TWO*Scoef(5) *windspeed /freq_S +TWO*Scoef(6) *windspeed/freq_S**2 &
            +Scoef(7) + TWO*Scoef(8) *windspeed  )*windspeed_tl

    small_corr = exp(-scor*cos_z*cos_z )
    small_corr_tl = -scor_tl*cos_z*cos_z *small_corr
    RvS = Rv_Fresnel * small_corr
    RvS_tl = (Rv_Fresnel_tl * small_corr + Rv_Fresnel * small_corr_tl)

    RhS = Rh_Fresnel * small_corr
    RhS_tl = (Rh_Fresnel_tl * small_corr + Rh_Fresnel * small_corr_tl)


  ! Large Scale Correction Calculation
  ! ----------------------------------
    seczen = ONE/cos_z
    ! compute fitting coefficients for a given frequency
    DO j = 1, 12
      zc(j) = Lcoef(j*3-2) + Lcoef(j*3-1)*frequency + Lcoef(j*3)*frequency**2
    END DO

    RvL = zc(1) + zc(2)*seczen + zc(3)*seczen**2 + zc(4)*Wind_Speed &
      + zc(5)*Wind_Speed**2 + zc(6)*Wind_Speed*seczen
    RhL = zc(7) + zc(8)*seczen + zc(9)*seczen**2 + zc(10)*Wind_Speed &
      + zc(11)*Wind_Speed**2 + zc(12)*Wind_Speed*seczen

    RvL_tl = zc(4)*Wind_Speed_tl &
      + TWO*zc(5)*Wind_Speed*Wind_Speed_tl + zc(6)*Wind_Speed_tl*seczen

    RhL_tl = zc(10)*Wind_Speed_tl &
      + TWO*zc(11)*Wind_Speed*Wind_Speed_tl + zc(12)*Wind_Speed_tl*seczen     
   
    ! Check wether the foam fraction computation is needed or not
    IF (PRESENT(Supply_Foam_Fraction)) THEN
      lcalc_foam_fraction = .NOT.(Supply_Foam_Fraction .AND. &
                                & PRESENT(Foam_Fraction) .AND. &
                                & PRESENT(Foam_Fraction_tl))
    ENDIF

    ! change foam coverage back to FASTEM1,2,3 for FASTEM5
    IF (lcalc_foam_fraction) THEN
      IF (fastem_version == 4) THEN
        ! Compute foam coverage after Tang, 1974
        Foam_Cover = 7.75E-06_fp * Wind_Speed ** 3.231_fp
        Foam_Cover_tl = 3.231_fp*7.75E-06_fp *Wind_Speed**(3.231_fp-ONE)*Wind_Speed_tl
      ELSE
        ! Monahan et al., 1986 without surface stability term
        Foam_Cover = 1.95E-05_fp * Wind_Speed ** 2.55_fp
        Foam_Cover_tl = 2.55_fp*1.95E-05_fp *Wind_Speed**(2.55_fp-ONE)*Wind_Speed_tl
      ENDIF
    ELSE
      Foam_Cover = Foam_Fraction
      Foam_Cover_tl = Foam_Fraction_tl
    ENDIF

  ! The foam vertical and horizontal reflectanc codes, adopted from Masahiro Kazumori, JMA
  ! ----------------------------------
    Foam_Rv = FR_COEFF(1)
    Fh = ONE + Zenith_Angle*(FR_COEFF(2) +  Zenith_Angle*(FR_COEFF(3)  &
       + Zenith_Angle*FR_COEFF(4)))
    Foam_Rh = ONE + FR_COEFF(5)*Fh

    ! Added frequency dependence derived from Stogry model, 1971
    Foam_ref = 0.4_fp * exp(-0.05_fp*Frequency )
    Foam_Rv = Foam_Rv * Foam_ref
    Foam_Rh = Foam_Rh * Foam_ref

    Ev = (ONE-Foam_Cover)*(ONE - RvS + RvL) + Foam_Cover*(ONE-Foam_Rv)
    Ev_tl = Foam_Cover_tl*(RvS - RvL-Foam_Rv) + (ONE-Foam_Cover)*(-RvS_tl + RvL_tl)

    Eh = (ONE-Foam_Cover)*(ONE - RhS + RhL) + Foam_Cover*(ONE-Foam_Rh)
    Eh_tl = Foam_Cover_tl*(RhS - RhL-Foam_Rh) + (ONE-Foam_Cover)*(-RhS_tl + RhL_tl)

    Emissivity(1) = Ev
    Emissivity(2) = Eh
    Emissivity_tl(1) = Ev_tl
    Emissivity_tl(2) = Eh_tl

    zreflmod_v = ONE
    zreflmod_h = ONE
    zreflmod_v_tl = ZERO
    zreflmod_h_tl = ZERO


  ! correction for anisotropic downward radiation, adopted from the FASTEM3
  ! ----------------------------------
    IF( Transmittance > transmittance_limit_lower .and. Transmittance < transmittance_limit_upper) THEN
        !Using the Cox and Munk model to compute slope variance
        variance = 0.00512_fp * Wind_Speed + 0.0030_fp
        varm     = variance * t_c(43)
        variance = varm * ( t_c(44) * Frequency + t_c(45) )

        variance_tl = 0.00512_fp * wind_speed_tl
        varm_tl     = variance_tl * t_c(43)
        variance_tl = varm_tl * ( t_c(44) * Frequency + t_c(45) )

        IF ( variance >= varm ) THEN
          variance = varm
          variance_tl = varm_tl
        END IF
        IF ( variance <= ZERO  ) THEN
          variance = ZERO
          variance_tl = ZERO
        END IF
        !Compute surface to space optical depth
        opdpsfc = -log(Transmittance ) * cos_z
        opdpsfc_tl = -Transmittance_TL * cos_z/ ( Transmittance )
        !Define nine predictors for the effective angle calculation
        zx(1) = ONE
        zx(2) = variance
        zx(4) = ONE / cos_z
        zx(3) = zx(2) * zx(4)
        zx(5) = zx(3) * zx(3)
        zx(6) = zx(4) * zx(4)
        zx(7) = zx(2) * zx(2)
        zx(8) = log(opdpsfc)
        zx(9) = zx(8) * zx(8)

        zx_tl = ZERO
        zx_tl(1) = ZERO
        zx_tl(2) = variance_tl
        zx_tl(4) = ZERO
        zx_tl(3) = zx_tl(2) * zx(4)
        zx_tl(5) = 2 * zx_tl(3) * zx(3)
        zx_tl(7) = 2 * zx_tl(2) * zx(2)
        zx_tl(8) = opdpsfc_tl / opdpsfc
        zx_tl(9) = 2 * zx_tl(8) * zx(8)
        zrough_v = ONE
        zrough_h = ONE
        zrough_v_tl = ZERO
        zrough_h_tl = ZERO

        DO i = 1, 7
           j = i-1
          !Switched h to v Deblonde SSMIS june 7, 2001
          zrough_h = zrough_h + zx(i) *(t_c(1+j*3) + zx(8)*t_c(2+j*3) + zx(9)*t_c(3+j*3) )
          zrough_v = zrough_v + zx(i) *(t_c(22+j*3)+ zx(8)*t_c(23+j*3)+ zx(9)*t_c(24+j*3))
          zrough_h_tl = zrough_h_tl + zx_tl(i) *(t_c(1+j*3) + zx(8)*t_c(2+j*3) &
            + zx(9)*t_c(3+j*3) ) + zx(i) *( zx_tl(8)*t_c(2+j*3) + zx_tl(9)*t_c(3+j*3) )
          zrough_v_tl = zrough_v_tl + zx_tl(i) *(t_c(22+j*3)+ zx(8)*t_c(23+j*3)  &
            + zx(9)*t_c(24+j*3) ) + zx(i) *( zx_tl(8)*t_c(23+j*3)+ zx_tl(9)*t_c(24+j*3) )

        END DO
        zreflmod_v = (ONE-Transmittance ** zrough_v)/(ONE-Transmittance )
        zreflmod_v_tl = Transmittance_tl/(ONE-Transmittance)*(-zrough_v  &
                      *Transmittance**(zrough_v-ONE) + zreflmod_v )
        zreflmod_v_tl = zreflmod_v_tl -Transmittance**zrough_v *log(Transmittance)  &
            * zrough_v_tl/(ONE-Transmittance)

        zreflmod_h = (ONE-Transmittance ** zrough_h)/(ONE-Transmittance )
        zreflmod_h_tl = Transmittance_tl/(ONE-Transmittance)*(-zrough_h  &
                      *Transmittance**(zrough_h-ONE) + zreflmod_h )
        zreflmod_h_tl = zreflmod_h_tl -Transmittance**zrough_h *log(Transmittance)  &
            *  zrough_h_tl/(ONE-Transmittance)
    END IF


  ! azimuthal component
  ! --------------------------------
    Azimuth_Emi = ZERO
    Azimuth_Emi_tl = ZERO

    IF( abs(Rel_Azimuth) <= 360.0_fp ) THEN
      if(fastem_version == 6)then ! M.Kazumori azimuth model function

        azimuth_component = ZERO
        azimuth_component_tl = ZERO

        phi = Rel_Azimuth * DEGREES_TO_RADIANS
        phi_tl = Rel_Azimuth_tl * DEGREES_TO_RADIANS
        wind10 = Wind_Speed
        wind10_tl = Wind_Speed_tl
        theta = Zenith_Angle

        ! freq.
        DO ifreq = 1,6
          IF(wind10>18.0_fp)THEN
            ipol=1
            A1v(ifreq) = coef_mk_azi(1,ifreq,ipol) * ( exp(-coef_mk_azi(5,ifreq,ipol) * 18.0_fp * 18.0_fp ) - ONE ) * &
                      ( coef_mk_azi(2,ifreq,ipol) * 18.0_fp + coef_mk_azi(3,ifreq,ipol) * 18.0_fp * 18.0_fp + &
                        coef_mk_azi(4,ifreq,ipol) * 18.0_fp * 18.0_fp * 18.0_fp )
            A1v_tl(ifreq) = ZERO

            A2v(ifreq) = coef_mk_azi(6,ifreq,ipol) * 18.0_fp
            A2v_tl(ifreq) =  ZERO

            ipol=2
            A1h(ifreq) = coef_mk_azi(1,ifreq,ipol) * 18.0_fp
            A1h_tl(ifreq) = ZERO
            A2h(ifreq) = coef_mk_azi(2,ifreq,ipol) * ( exp(-coef_mk_azi(6,ifreq,ipol) * 18.0_fp * 18.0_fp ) - ONE ) * &
                      ( coef_mk_azi(3,ifreq,ipol) * 18.0_fp + coef_mk_azi(4,ifreq,ipol) * 18.0_fp * 18.0_fp + &
                        coef_mk_azi(5,ifreq,ipol) * 18.0_fp * 18.0_fp * 18.0_fp )
            A2h_tl(ifreq) = ZERO
          ELSE
            ipol=1
            A1v(ifreq) = coef_mk_azi(1,ifreq,ipol) * ( exp(-coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                      ( coef_mk_azi(2,ifreq,ipol) * wind10 + coef_mk_azi(3,ifreq,ipol) * wind10 * wind10 + &
                        coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 * wind10 )
            A1v_tl(ifreq) = coef_mk_azi(1,ifreq,ipol) *  &
  (&
                    -TWO*wind10*coef_mk_azi(5,ifreq,ipol)*( exp(-coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 ) )*wind10_tl * &
                    ( coef_mk_azi(2,ifreq,ipol) * wind10 + coef_mk_azi(3,ifreq,ipol) * wind10 * wind10 + &
                      coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 * wind10 ) &
  + &
                    ( exp(-coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                    ( coef_mk_azi(2,ifreq,ipol) * wind10_tl + TWO*coef_mk_azi(3,ifreq,ipol) * wind10 * wind10_tl + &
                      THREE*coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 * wind10_tl ) &
  )
            A2v(ifreq) = coef_mk_azi(6,ifreq,ipol) * wind10
            A2v_tl(ifreq) = coef_mk_azi(6,ifreq,ipol) * wind10_tl

            ipol=2
            A1h(ifreq) = coef_mk_azi(1,ifreq,ipol) * wind10
            A1h_tl(ifreq) = coef_mk_azi(1,ifreq,ipol) * wind10_tl

            A2h(ifreq) = coef_mk_azi(2,ifreq,ipol) * ( exp(-coef_mk_azi(6,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                    ( coef_mk_azi(3,ifreq,ipol) * wind10 + coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 + &
                      coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 * wind10 )

            A2h_tl(ifreq) = coef_mk_azi(2,ifreq,ipol) * &
  ( &
                    ( -TWO*coef_mk_azi(6,ifreq,ipol) * wind10*exp(-coef_mk_azi(6,ifreq,ipol) * wind10 * wind10 ) *wind10_tl ) * &
                    ( coef_mk_azi(3,ifreq,ipol) * wind10 + coef_mk_azi(4,ifreq,ipol) * wind10 * wind10 + &
                      coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 * wind10 ) &
  + &
                    ( exp(-coef_mk_azi(6,ifreq,ipol) * wind10 * wind10 ) - ONE ) * &
                    ( coef_mk_azi(3,ifreq,ipol) * wind10_tl + TWO*coef_mk_azi(4,ifreq,ipol) * wind10 * wind10_tl + &
                      THREE*coef_mk_azi(5,ifreq,ipol) * wind10 * wind10 * wind10_tl ) &
  )

          END IF
          A1s1(ifreq) = (A1v(ifreq) + A1h(ifreq))/TWO
          A1s1_tl(ifreq) = (A1v_tl(ifreq) + A1h_tl(ifreq))/TWO
          A1s2(ifreq) =  A1v(ifreq) - A1h(ifreq)
          A1s2_tl(ifreq) =  A1v_tl(ifreq) - A1h_tl(ifreq)
          A2s1(ifreq) = (A2v(ifreq) + A2h(ifreq))/TWO
          A2s1_tl(ifreq) = (A2v_tl(ifreq) + A2h_tl(ifreq))/TWO
          A2s2(ifreq) =  A2v(ifreq) - A2h(ifreq)
          A2s2_tl(ifreq) =  A2v_tl(ifreq) - A2h_tl(ifreq)

          IF(Frequency>37.0_fp)THEN
            IF(wind10>15.0_fp)THEN
              A2s2_theta0(ifreq) = (15.0_fp*15.0_fp - 15.0_fp*15.0_fp*15.0_fp/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/37.0_fp) )
              A2s2_theta0_tl(ifreq) = 0.0_fp
            ELSE
              A2s2_theta0(ifreq) = (wind10*wind10 - wind10*wind10*wind10/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/37.0_fp) )
              A2s2_theta0_tl(ifreq) = (TWO*wind10*wind10_tl - THREE*wind10*wind10*wind10_tl/22.5d0)/55.5556d0 * &
                                      (2.d0/290.d0)*(1.0d0 - log10(30.0d0/37.0_fp) )
            END IF
          ELSE
            IF(wind10>15.0_fp)THEN
              A2s2_theta0(ifreq) = (15.0_fp*15.0_fp - 15.0_fp*15.0_fp*15.0_fp/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/Frequency) )
              A2s2_theta0_tl(ifreq) = 0.0_fp
            ELSE
              A2s2_theta0(ifreq) = (wind10*wind10 - wind10*wind10*wind10/22.5d0)/55.5556d0 * &
                                   (2.d0/290.d0)*(1.0d0 - log10(30.0d0/Frequency) )
              A2s2_theta0_tl(ifreq) = (TWO*wind10*wind10_tl - THREE*wind10*wind10*wind10_tl/22.5d0)/55.5556d0 * &
                                      (2.d0/290.d0)*(1.0d0 - log10(30.0d0/Frequency) )
            END IF
          END IF

          A1s1_theta(ifreq)= A1s1(ifreq)*((theta/theta_ref)**xs11)
          A1s1_theta_tl(ifreq)= A1s1_tl(ifreq)*((theta/theta_ref)**xs11)
          A2s1_theta(ifreq)= A2s1(ifreq)*((theta/theta_ref)**xs12)
          A2s1_theta_tl(ifreq)= A2s1_tl(ifreq)*((theta/theta_ref)**xs12)
          A1s2_theta(ifreq)= A1s2(ifreq)*((theta/theta_ref)**xs21)
          A1s2_theta_tl(ifreq)= A1s2_tl(ifreq)*((theta/theta_ref)**xs21)
          A2s2_theta(ifreq)= A2s2_theta0(ifreq) + (A2s2(ifreq) - A2s2_theta0(ifreq))*((theta/theta_ref)**xs22)
          A2s2_theta_tl(ifreq)= A2s2_theta0_tl(ifreq) + (A2s2_tl(ifreq) - A2s2_theta0_tl(ifreq))*((theta/theta_ref)**xs22)

          A1v_theta(ifreq) = 0.5d0*(2.d0*A1s1_theta(ifreq) + A1s2_theta(ifreq))
          A1v_theta_tl(ifreq) = 0.5d0*(2.d0*A1s1_theta_tl(ifreq) + A1s2_theta_tl(ifreq))
          A1h_theta(ifreq) = 0.5d0*(2.d0*A1s1_theta(ifreq) - A1s2_theta(ifreq))
          A1h_theta_tl(ifreq) = 0.5d0*(2.d0*A1s1_theta_tl(ifreq) - A1s2_theta_tl(ifreq))
          A2v_theta(ifreq) = 0.5d0*(2.d0*A2s1_theta(ifreq) + A2s2_theta(ifreq))
          A2v_theta_tl(ifreq) = 0.5d0*(2.d0*A2s1_theta_tl(ifreq) + A2s2_theta_tl(ifreq))
          A2h_theta(ifreq) = 0.5d0*(2.d0*A2s1_theta(ifreq) - A2s2_theta(ifreq))
          A2h_theta_tl(ifreq) = 0.5d0*(2.d0*A2s1_theta_tl(ifreq) - A2s2_theta_tl(ifreq))

          azimuth_component(1,ifreq) = A1v_theta(ifreq) * cos (real(1,8)*phi) + A2v_theta(ifreq) * cos(real(2,8)*phi)
          azimuth_component_tl(1,ifreq) = A1v_theta_tl(ifreq) * cos (real(1,8)*phi) + &
                                          A1v_theta(ifreq) * (-real(1,8))*sin (real(1,8)*phi)*phi_tl + &
                                          A2v_theta_tl(ifreq) * cos (real(2,8)*phi) + &
                                          A2v_theta(ifreq) * (-real(2,8))*sin(real(2,8)*phi) *phi_tl

          azimuth_component(2,ifreq) = A1h_theta(ifreq) * cos (real(1,8)*phi) + A2h_theta(ifreq) * cos(real(2,8)*phi)
          azimuth_component_tl(2,ifreq) = A1h_theta_tl(ifreq) * cos (real(1,8)*phi) + &
                                          A1h_theta(ifreq) * (-real(1,8))*sin (real(1,8)*phi)*phi_tl + &
                                          A2h_theta_tl(ifreq) * cos (real(2,8)*phi) + &
                                          A2h_theta(ifreq) * (-real(2,8))*sin(real(2,8)*phi) *phi_tl
        END DO

        IF( Frequency >= 1.4_fp .and. Frequency < 6.925_fp ) THEN
          Azimuth_Emi(1) = azimuth_component(1,1)
          Azimuth_Emi(2) = azimuth_component(2,1)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,1)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,1)
        ELSE IF( Frequency >= 6.925_fp .and. Frequency < 10.65_fp ) THEN
          fratio = ONE-(Frequency - 6.925_fp)/(10.65_fp - 6.925_fp)
          Azimuth_Emi(1) = azimuth_component(1,1)*fratio + (ONE-fratio)*azimuth_component(1,2)
          Azimuth_Emi(2) = azimuth_component(2,1)*fratio + (ONE-fratio)*azimuth_component(2,2)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,1)*fratio + (ONE-fratio)*azimuth_component_tl(1,2)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,1)*fratio + (ONE-fratio)*azimuth_component_tl(2,2)
        ELSE IF( Frequency > 10.65_fp .and. Frequency <= 18.7_fp ) THEN
          fratio = ONE-(Frequency - 10.65_fp)/(18.7_fp - 10.65_fp)
          Azimuth_Emi(1) = azimuth_component(1,2)*fratio + (ONE-fratio)*azimuth_component(1,3)
          Azimuth_Emi(2) = azimuth_component(2,2)*fratio + (ONE-fratio)*azimuth_component(2,3)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,2)*fratio + (ONE-fratio)*azimuth_component_tl(1,3)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,2)*fratio + (ONE-fratio)*azimuth_component_tl(2,3)
        ELSE IF( Frequency > 18.7_fp .and. Frequency <= 23.8_fp ) THEN
          fratio = ONE-(Frequency - 18.7_fp)/(23.8_fp - 18.7_fp)
          Azimuth_Emi(1) = azimuth_component(1,3)*fratio + (ONE-fratio)*azimuth_component(1,4)
          Azimuth_Emi(2) = azimuth_component(2,3)*fratio + (ONE-fratio)*azimuth_component(2,4)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,3)*fratio + (ONE-fratio)*azimuth_component_tl(1,4)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,3)*fratio + (ONE-fratio)*azimuth_component_tl(2,4)
        ELSE IF( Frequency > 23.8_fp .and. Frequency <= 36.5_fp ) THEN
          fratio = ONE-(Frequency - 23.8_fp)/(36.5_fp - 23.8_fp)
          Azimuth_Emi(1) = azimuth_component(1,4)*fratio + (ONE-fratio)*azimuth_component(1,5)
          Azimuth_Emi(2) = azimuth_component(2,4)*fratio + (ONE-fratio)*azimuth_component(2,5)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,4)*fratio + (ONE-fratio)*azimuth_component_tl(1,5)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,4)*fratio + (ONE-fratio)*azimuth_component_tl(2,5)
        ELSE IF( Frequency > 36.5_fp .and. Frequency <= 89.0_fp ) THEN
          fratio = ONE-(Frequency - 36.5_fp)/(89.0_fp - 36.5_fp)
          Azimuth_Emi(1) = azimuth_component(1,5)*fratio + (ONE-fratio)*azimuth_component(1,6)
          Azimuth_Emi(2) = azimuth_component(2,5)*fratio + (ONE-fratio)*azimuth_component(2,6)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,5)*fratio + (ONE-fratio)*azimuth_component_tl(1,6)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,5)*fratio + (ONE-fratio)*azimuth_component_tl(2,6)
        ELSE IF( Frequency > 89.0_fp .and. Frequency <= 200.0_fp ) THEN
          Azimuth_Emi(1) = azimuth_component(1,6)
          Azimuth_Emi(2) = azimuth_component(2,6)
          Azimuth_Emi_tl(1) = azimuth_component_tl(1,6)
          Azimuth_Emi_tl(2) = azimuth_component_tl(2,6)
        END IF

      else                     ! M.Liu      azimuth model function 

        Fre_C = ZERO
        IF( Frequency >= min_f .or. Frequency <= max_f ) THEN
          DO i = 1, 8
            IF( Frequency >= x(i) .and. Frequency < x(i+1) ) THEN
              Fre_C = y(i) + (y(i+1)-y(i))/(x(i+1)-x(i))*(Frequency-x(i))
            END IF
          END DO
        END IF

        phi = Rel_Azimuth * DEGREES_TO_RADIANS
        phi_tl = Rel_Azimuth_tl * DEGREES_TO_RADIANS
        wind10 = Wind_Speed
        wind10_tl = Wind_Speed_tl

        DO m = 1, 3
          L = 10*(m-1)
          ac = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(1) = Azimuth_Emi(1) + ac*cos(m*phi)
          ac_tl = ( b_coef(L+5) +b_coef(L+6)*Frequency ) *wind10_tl +( TWO*(b_coef(L+7) &
              +b_coef(L+8)*Frequency)*wind10 +b_coef(L+9)*seczen   &
              +b_coef(L+10)*seczen*Frequency )*wind10_tl
          Azimuth_Emi_tl(1) = Azimuth_Emi_tl(1) + ac_tl*cos(m*phi) -m*phi_tl*ac*sin(m*phi)

          L = 10*(m-1) + 30
          ac = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(2) = Azimuth_Emi(2) + ac*cos(m*phi)
          ac_tl = ( b_coef(L+5) +b_coef(L+6)*Frequency ) *wind10_tl +( TWO*(b_coef(L+7) &
              +b_coef(L+8)*Frequency)*wind10 +b_coef(L+9)*seczen   &
              +b_coef(L+10)*seczen*Frequency )*wind10_tl
          Azimuth_Emi_tl(2) = Azimuth_Emi_tl(2) + ac_tl*cos(m*phi) -m*phi_tl*ac*sin(m*phi)

          L = 10*(m-1) + 60
          sc = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(3) = Azimuth_Emi(3) + sc*sin(m*phi)
          sc_tl = ( b_coef(L+5) +b_coef(L+6)*Frequency ) *wind10_tl +( TWO*(b_coef(L+7) &
              +b_coef(L+8)*Frequency)*wind10 +b_coef(L+9)*seczen   &
              +b_coef(L+10)*seczen*Frequency )*wind10_tl
          Azimuth_Emi_tl(3) = Azimuth_Emi_tl(3) + sc_tl*sin(m*phi) +m*phi_tl*sc*cos(m*phi)

          L = 10*(m-1) + 90
          sc = b_coef(L+1) +b_coef(L+2)*Frequency +b_coef(L+3)*seczen   &
            +b_coef(L+4)*seczen*Frequency &
            +b_coef(L+5)*wind10 +b_coef(L+6)*wind10*Frequency +b_coef(L+7)*wind10**2  &
            +b_coef(L+8)*Frequency*wind10**2 +b_coef(L+9)*wind10*seczen   &
            +b_coef(L+10)*wind10*seczen*Frequency
          Azimuth_Emi(4) = Azimuth_Emi(4) + sc*sin(m*phi)
          sc_tl = ( b_coef(L+5) +b_coef(L+6)*Frequency ) *wind10_tl +( TWO*(b_coef(L+7) &
              +b_coef(L+8)*Frequency)*wind10 +b_coef(L+9)*seczen   &
              +b_coef(L+10)*seczen*Frequency )*wind10_tl
          Azimuth_Emi_tl(4) = Azimuth_Emi_tl(4) + sc_tl*sin(m*phi) +m*phi_tl*sc*cos(m*phi)
        END DO
        Azimuth_Emi = Azimuth_Emi * Fre_C
        Azimuth_Emi_tl = Azimuth_Emi_tl * Fre_C

      endif
    END IF

    Emissivity(1) = Emissivity(1) + Azimuth_Emi(1)
    Emissivity_tl(1) = Emissivity_tl(1) + Azimuth_Emi_tl(1)
    Emissivity(2) = Emissivity(2) + Azimuth_Emi(2)
    Emissivity_tl(2) = Emissivity_tl(2) + Azimuth_Emi_tl(2)
    Emissivity(3) = Azimuth_Emi(3)
    Emissivity_tl(3) = Azimuth_Emi_tl(3)
    Emissivity(4) = Azimuth_Emi(4)
    Emissivity_tl(4) = Azimuth_Emi_tl(4)
    Reflectivity(1)  = zreflmod_v * (ONE-Emissivity(1))
    Reflectivity_tl(1)  = zreflmod_v_tl*(ONE-Emissivity(1)) -zreflmod_v*Emissivity_tl(1)
    Reflectivity(2)  = zreflmod_h * (ONE-Emissivity(2))
    Reflectivity_tl(2)  = zreflmod_h_tl * (ONE-Emissivity(2)) -zreflmod_h*Emissivity_tl(2)
    ! Reflectivities not computed for 3rd or 4th elements of Stokes vector, 
    ! as never used subsequently, as atmospheric source term = zero.
    Reflectivity(3:4)  = ZERO
    Reflectivity_tl(3:4)  = ZERO

   RETURN

  END SUBROUTINE rttov_fastem5r1_tl

END MODULE RTTOV_FASTEM5R1_TL_MODULE

