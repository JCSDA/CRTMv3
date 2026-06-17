
MODULE NESDIS_LandEM_Module
 
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE Type_Kinds, ONLY: fp
  USE Fundamental_Constants, ONLY: C_2
  USE NESDIS_SnowEM_Parameters
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Procedures
  PUBLIC  :: NESDIS_LandEM
  ! Parameters
  PUBLIC :: ZERO
  PUBLIC :: POINT1
  PUBLIC :: POINT5
  PUBLIC :: ONE
  PUBLIC :: TWO
  PUBLIC :: THREE
  PUBLIC :: FOUR
  PUBLIC :: PI
  PUBLIC :: EMISSH_DEFAULT
  PUBLIC :: EMISSV_DEFAULT
  
  PUBLIC :: ONE_TENTH
  PUBLIC :: HALF
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: ZERO   = 0.0_fp
  REAL(fp), PARAMETER :: POINT1 = 0.1_fp
  REAL(fp), PARAMETER :: POINT5 = 0.5_fp
  REAL(fp), PARAMETER :: ONE    = 1.0_fp
  REAL(fp), PARAMETER :: TWO    = 2.0_fp
  REAL(fp), PARAMETER :: THREE  = 3.0_fp
  REAL(fp), PARAMETER :: FOUR   = 4.0_fp
  REAL(fp), PARAMETER :: PI     = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER :: TWOPI  = TWO*PI
  REAL(fp), PARAMETER :: EMISSH_DEFAULT = 0.25_fp
  REAL(fp), PARAMETER :: EMISSV_DEFAULT = 0.30_fp

  REAL(fp), PARAMETER :: ONE_TENTH = POINT1
  REAL(fp), PARAMETER :: HALF      = POINT5
  
  
CONTAINS




  SUBROUTINE NESDIS_LandEM(Angle,                 &   ! Input
                           Frequency,             &   ! Input
                           Soil_Moisture_Content, &   ! Input
                           Vegetation_Fraction,   &   ! Input
                           Soil_Temperature,      &   ! Input
                           t_skin,                &   ! Input
                           Lai,                   &   ! Input
                           Soil_Type,             &   ! Input
                           Vegetation_Type,       &   ! Input
                           Snow_Depth,            &   ! Input
                           Emissivity_H,          &   ! Output
                           Emissivity_V,          &   ! Output
                           dEV_dvlai,             &   ! Optional Output
                           dEH_dvlai,             &   ! Optional Output
                           dEV_dmv,               &   ! Optional Output
                           dEH_dmv)                   ! Optional Output
    ! Arguments
    REAL(fp), intent(in) :: Angle
    REAL(fp), intent(in) :: Frequency
    REAL(fp), intent(in) :: Soil_Moisture_Content
    REAL(fp), intent(in) :: Vegetation_Fraction
    REAL(fp), intent(in) :: Soil_Temperature
    REAL(fp), intent(in) :: t_skin
    REAL(fp), intent(in) :: Lai
    INTEGER,  intent(in) :: Soil_Type
    INTEGER,  intent(in) :: Vegetation_Type
    REAL(fp), intent(in) :: Snow_Depth
    REAL(fp), intent(out):: Emissivity_V,Emissivity_H
    ! Optional tangent-linear/adjoint support, returned only for the (no-snow)
    ! canopy path and set to zero otherwise (snow/ice callers omit them):
    !   dEV_dvlai/dEH_dvlai : d(emissivity)/d(vlai), vlai = Lai*Vegetation_Fraction
    !   dEV_dmv/dEH_dmv     : d(emissivity)/d(Soil_Moisture_Content)
    REAL(fp), OPTIONAL, intent(out) :: dEV_dvlai, dEH_dvlai
    REAL(fp), OPTIONAL, intent(out) :: dEV_dmv, dEH_dmv
    ! Local parameters
    REAL(fp), PARAMETER :: snow_depth_c     = 10.0_fp
    REAL(fp), PARAMETER :: tsoilc_undersnow = 280.0_fp
    REAL(fp), PARAMETER :: rhos = 2.65_fp
    REAL(fp)            :: sand, clay, rhob
    REAL(fp), PARAMETER, dimension(0:9) :: frac_sand = (/ 0.80_fp,     &
                          0.92_fp, 0.10_fp, 0.20_fp, 0.51_fp, 0.50_fp, &
                          0.35_fp, 0.60_fp, 0.42_fp,  0.92_fp  /)
    REAL(fp), PARAMETER, dimension(0:9) :: frac_clay = (/ 0.20_fp,     &
                          0.06_fp, 0.34_fp, 0.63_fp, 0.14_fp, 0.43_fp, &
                          0.34_fp, 0.28_fp, 0.085_fp, 0.06_fp /)
    REAL(fp), PARAMETER, dimension(0:9) :: rhob_soil = (/ 1.48_fp,     &
                          1.68_fp, 1.27_fp, 1.21_fp, 1.48_fp, 1.31_fp, &
                          1.32_fp, 1.40_fp, 1.54_fp, 1.68_fp /)
    REAL(fp), PARAMETER, dimension(0:13) :: veg_rho  = (/ 0.33_fp,     &
                          0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
                          0.25_fp, 0.25_fp, 0.40_fp, 0.40_fp, 0.40_fp, &
                          0.40_fp, 0.33_fp, 0.33_fp            /)
    REAL(fp), PARAMETER, dimension(0:13) :: veg_mge  = (/ 0.50_fp,     &
                          0.45_fp, 0.45_fp, 0.45_fp, 0.40_fp, 0.40_fp, &
                          0.30_fp, 0.35_fp, 0.30_fp, 0.30_fp, 0.40_fp, &
                          0.30_fp, 0.50_fp, 0.40_fp            /)
    REAL(fp), PARAMETER, dimension(0:13) :: lai_min  = (/ 0.52_fp,     &
                          3.08_fp, 1.85_fp, 2.80_fp, 5.00_fp, 1.00_fp, &
                          0.50_fp, 0.52_fp, 0.60_fp, 0.50_fp, 0.60_fp, &
                          0.10_fp, 1.56_fp, 0.01_fp            /)
    REAL(fp), PARAMETER, dimension(0:13) :: lai_max  = (/ 2.90_fp,     &
                          6.48_fp, 3.31_fp, 5.50_fp, 6.40_fp, 5.16_fp, &
                          3.66_fp, 2.90_fp, 2.60_fp, 3.66_fp, 2.60_fp, &
                          0.75_fp, 5.68_fp, 0.01_fp            /)
    REAL(fp), PARAMETER, dimension(0:13) :: leaf_th  = (/ 0.07_fp,     &
                          0.18_fp, 0.18_fp, 0.18_fp, 0.18_fp, 0.18_fp, &
                          0.12_fp, 0.12_fp, 0.12_fp, 0.12_fp, 0.12_fp, &
                          0.12_fp, 0.15_fp, 0.12_fp            /)
    ! Local variables
    REAL(fp) :: mv,veg_frac,theta,theta_i,theta_t,mu,r21_h,r21_v,r23_h,r23_v,  &
                t21_v,t21_h,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v,mge, &
                    leaf_thick,rad,sigma,va,ep_real,ep_imag
    REAL(fp) :: t_soil
    REAL(fp) :: rhoveg, vlai
    REAL(fp) :: local_snow_depth
    REAL(fp) :: dtau_dvlai_l, desv_dtauv_l, desh_dtauh_l
    REAL(fp) :: desv_dr23v_l, desh_dr23h_l
    REAL(fp) :: dtheta_t_dmv, ds_dmv, cos_i_l, cos_tt, sin_tt, qr_l
    REAL(fp) :: drv0_dmv, drh0_dmv, dr23v_dmv, dr23h_dmv
    COMPLEX(fp) :: esoil, eveg, esnow, eair
    COMPLEX(fp) :: desoil_dmv, dqc_dmv, m1c, m2c, dm2_dmv, angle_i_c, angle_t_c, dangle_t_dmv
    COMPLEX(fp) :: nv_c, dv_c, dnv_c, ddv_c, zv_c, dzv_c
    COMPLEX(fp) :: nh_c, dh_c, dnh_c, ddh_c, zh_c, dzh_c
    LOGICAL :: SnowEM_Physical_Model

    eair = CMPLX(ONE,-ZERO,fp)
    theta = Angle*PI/180.0_fp

    ! Default the optional derivative outputs (filled only on the canopy path)
    IF (PRESENT(dEV_dvlai)) dEV_dvlai = ZERO
    IF (PRESENT(dEH_dvlai)) dEH_dvlai = ZERO
    IF (PRESENT(dEV_dmv))   dEV_dmv   = ZERO
    IF (PRESENT(dEH_dmv))   dEH_dmv   = ZERO

    ! By default use the 
    ! Assign local variable
    mv               = Soil_Moisture_Content
    veg_frac         = Vegetation_Fraction
    t_soil           = Soil_Temperature
    sand = frac_sand(Soil_Type)
    clay = frac_clay(Soil_Type )
    rhob = rhob_soil(Soil_Type )
    local_snow_depth = Snow_Depth

    ! Check soil/skin temperature
    if ( (t_soil <= 100.0_fp .OR.  t_soil >= 350.0_fp) .AND. &
         (t_skin >= 100.0_fp .AND. t_skin <= 350.0_fp) ) t_soil = t_skin

    ! Check soil moisture content range
    mv = MAX(MIN(mv,ONE),ZERO)

    ! Surface type based on snow depth
    IF (local_snow_depth > POINT1) THEN

      ! O.k.; we're going to compute snow emissivities....
      
      ! By default use the physical model for snow
      SnowEM_Physical_Model = .TRUE.
      if (local_snow_depth > snow_depth_c) SnowEM_Physical_Model = .FALSE.

      ! Compute the snow emissivity
      IF ( SnowEM_Physical_Model ) THEN

        ep_real = 3.2_fp
        ep_imag = -0.0005_fp
        sigma = ONE

        ! For deep snow, the performance of the model is poor
        local_snow_depth = MIN(local_snow_depth,1000.0_fp)

        ! The fraction volume of dense medium
        ! scatterers must be in the range (0-1)
        va = 0.4_fp + 0.0004_fp*local_snow_depth
        va = MAX(MIN(va,ONE),ZERO)

        ! Limit for snow grain size
        rad = MIN((POINT5 + 0.005_fp*local_snow_depth),ONE)

        ! Limit for soil temperature
        t_soil = MIN(t_soil,tsoilc_undersnow)

        CALL Snow_Diel(Frequency, ep_real, ep_imag, rad, va, esnow)
        CALL Soil_Diel(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil)

        theta_i = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(esnow),fp))

        CALL Reflectance(esnow, eair, theta_i,  theta, r21_v, r21_h)
        CALL Transmittance(esnow, eair, theta_i, theta, t21_v, t21_h)

        mu      = COS(theta_i)
        theta_t = ASIN(REAL(SIN(theta_i)*SQRT(esnow)/SQRT(esoil),fp))

        CALL Reflectance(esnow, esoil, theta_i, theta_t, r23_v, r23_h)
        CALL Roughness_Reflectance(Frequency, sigma, r23_v, r23_h)
        CALL Snow_Optic(Frequency,rad,local_snow_depth,va,ep_real, ep_imag,gv,gh,&
                        ssalb_v,ssalb_h,tau_v,tau_h)
        CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                                 r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                               frequency, t_soil, t_skin)
      ELSE
        ! Use the empirical method 
        CALL SnowEM_Default(Frequency,t_skin, local_snow_depth,Emissivity_V,Emissivity_H)
      END IF

    ELSE

      ! No snow, so we're going to compute canopy emissivities....
      
      ! Limit for vegetation fraction
      veg_frac = MAX(MIN(veg_frac,ONE),ZERO)

      mu  = COS(theta)
      sigma = POINT5
    
      vlai = Lai*veg_frac
      mge = veg_mge(Vegetation_Type)
      rhoveg = veg_rho(Vegetation_Type)
      leaf_thick = leaf_th(Vegetation_Type)

      r21_h    = ZERO
      r21_v    = ZERO
      t21_h    = ONE
      t21_v    = ONE

      CALL Soil_Diel(Frequency, t_soil, mv, rhob, rhos, sand, clay, esoil, &
                     desm_dvmc=desoil_dmv)
      theta_t = ASIN(REAL(SIN(theta)*SQRT(eair)/SQRT(esoil),fp))
      CALL Reflectance(eair, esoil, theta, theta_t, r23_v, r23_h)
      CALL Roughness_Reflectance(Frequency, sigma, r23_v, r23_h)
      CALL Canopy_Diel(Frequency, mge, eveg, rhoveg)
      CALL Canopy_Optic(vlai,Frequency,theta,eveg,leaf_thick,gv,gh,ssalb_v,ssalb_h,tau_v,tau_h, &
                        dtau_dvlai=dtau_dvlai_l)
      CALL Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
                               r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,Emissivity_V,Emissivity_H, &
                               frequency, t_soil, t_skin, &
                               desv_dtauv=desv_dtauv_l, desh_dtauh=desh_dtauh_l, &
                               desv_dr23v=desv_dr23v_l, desh_dr23h=desh_dr23h_l)

      ! Chain rule: d(emissivity)/d(vlai) = d(emissivity)/d(tau) * d(tau)/d(vlai)
      IF (PRESENT(dEV_dvlai)) dEV_dvlai = desv_dtauv_l*dtau_dvlai_l
      IF (PRESENT(dEH_dvlai)) dEH_dvlai = desh_dtauh_l*dtau_dvlai_l

      ! Chain rule for soil moisture: mv -> esoil (Soil_Diel) -> {theta_t, Fresnel
      ! reflectances rv0,rh0} -> roughened r23 -> emissivity (two-stream).
      IF ( PRESENT(dEV_dmv) .OR. PRESENT(dEH_dmv) ) THEN
        ! theta_t = ASIN(arg), arg = SIN(theta)*SQRT(eair)/SQRT(esoil); sin(theta_t)=arg
        ! d/dmv of arg = SIN(theta)*SQRT(eair)*esoil^(-1/2): -1/2*esoil^(-3/2)*desoil_dmv
        dqc_dmv      = SIN(theta)*SQRT(eair)*(-POINT5)*esoil**(-1.5_fp)*desoil_dmv
        ds_dmv       = REAL(dqc_dmv, fp)
        cos_tt       = COS(theta_t)
        sin_tt       = SIN(theta_t)
        dtheta_t_dmv = ds_dmv/cos_tt
        ! Fresnel reflectance derivatives (medium 1 = air, medium 2 = soil)
        cos_i_l      = COS(theta)
        m1c          = SQRT(eair)
        m2c          = SQRT(esoil)
        dm2_dmv      = POINT5*esoil**(-POINT5)*desoil_dmv
        angle_i_c    = CMPLX(cos_i_l, ZERO, fp)
        angle_t_c    = CMPLX(cos_tt,  ZERO, fp)
        dangle_t_dmv = CMPLX(-sin_tt*dtheta_t_dmv, ZERO, fp)
        ! Vertical polarisation: rv0 = |nv/dv|^2
        nv_c  = m1c*angle_t_c - m2c*angle_i_c
        dv_c  = m1c*angle_t_c + m2c*angle_i_c
        dnv_c = m1c*dangle_t_dmv - dm2_dmv*angle_i_c
        ddv_c = m1c*dangle_t_dmv + dm2_dmv*angle_i_c
        zv_c  = nv_c/dv_c
        dzv_c = (dnv_c*dv_c - nv_c*ddv_c)/(dv_c*dv_c)
        drv0_dmv = TWO*REAL(CONJG(zv_c)*dzv_c, fp)
        ! Horizontal polarisation: rh0 = |nh/dh|^2
        nh_c  = m1c*angle_i_c - m2c*angle_t_c
        dh_c  = m1c*angle_i_c + m2c*angle_t_c
        dnh_c = -(dm2_dmv*angle_t_c + m2c*dangle_t_dmv)
        ddh_c =  (dm2_dmv*angle_t_c + m2c*dangle_t_dmv)
        zh_c  = nh_c/dh_c
        dzh_c = (dnh_c*dh_c - nh_c*ddh_c)/(dh_c*dh_c)
        drh0_dmv = TWO*REAL(CONJG(zh_c)*dzh_c, fp)
        ! Roughness mixing (linear), matching Roughness_Reflectance
        qr_l = 0.35_fp*(ONE - EXP(-0.60_fp*Frequency*sigma**TWO))
        dr23h_dmv = 0.3_fp*drh0_dmv + qr_l*(0.3_fp*drv0_dmv - 0.3_fp*drh0_dmv)
        dr23v_dmv = 0.3_fp*drv0_dmv + qr_l*(0.3_fp*drh0_dmv - 0.3_fp*drv0_dmv)
        IF (PRESENT(dEV_dmv)) dEV_dmv = desv_dr23v_l*dr23v_dmv
        IF (PRESENT(dEH_dmv)) dEH_dmv = desh_dr23h_l*dr23h_dmv
      END IF
    END IF

  END SUBROUTINE NESDIS_LandEM






subroutine SnowEM_Default(frequency,ts, depth,Emissivity_V,Emissivity_H)


  ! Arguments
  REAL(fp) :: frequency,ts, depth,Emissivity_V,Emissivity_H
  ! Local parameters  
  INTEGER , PARAMETER :: new = 7
  INTEGER , PARAMETER :: NFRESH_SHALLOW_SNOW = 1
  INTEGER , PARAMETER :: NPOWDER_SNOW        = 2
  INTEGER , PARAMETER :: NWET_SNOW           = 3
  INTEGER , PARAMETER :: NDEEP_SNOW          = 4
  REAL(fp), PARAMETER :: twet    = 270.0_fp
  REAL(fp), PARAMETER :: tcrust  = 235.0_fp
  REAL(fp), PARAMETER :: depth_s =  50.0_fp
  REAL(fp), PARAMETER :: depth_c = 100.0_fp
  ! Local variables
  INTEGER :: ich,basic_snow_type
  REAL(fp), DIMENSION(new) :: ev, eh, freq
  REAL(fp) :: df, df0


  freq = FREQUENCY_AMSRE(1:new)
  
  ! Determine the snow type based on temperatures
  basic_snow_type = NFRESH_SHALLOW_SNOW
  if (ts >= twet .and. depth <= depth_s) then
    basic_snow_type = NWET_SNOW
  else
    if (depth <= depth_s) then
      basic_snow_type = NFRESH_SHALLOW_SNOW
    else
      basic_snow_type = NPOWDER_SNOW
    endif
  endif
  if (ts <= tcrust .and. depth >= depth_c) basic_snow_type = NDEEP_SNOW

  ! Assign the emissivity spectrum
  SELECT CASE (basic_snow_type)
    CASE (NFRESH_SHALLOW_SNOW)
      ev = GRASS_AFTER_SNOW_EV_AMSRE(1:new)
      eh = GRASS_AFTER_SNOW_EH_AMSRE(1:new)
    CASE (NPOWDER_SNOW)
      ev = POWDER_SNOW_EV_AMSRE(1:new)
      eh = POWDER_SNOW_EH_AMSRE(1:new)
    CASE (NWET_SNOW)
      ev = WET_SNOW_EV_AMSRE(1:new)
      eh = WET_SNOW_EH_AMSRE(1:new)
    CASE (NDEEP_SNOW)
      ev = DEEP_SNOW_EV_AMSRE(1:new)
      eh = DEEP_SNOW_EH_AMSRE(1:new)
  END SELECT

  ! Handle possible extrapolation
  if (frequency <= freq(1)) then
    Emissivity_H = eh(1)
    Emissivity_V = ev(1)
    return
  end if
  if (frequency >= freq(new)) then
    Emissivity_H = eh(new)
    Emissivity_V = ev(new)
    return
  end if

  ! Interpolate emissivity at a certain frequency
  Channel_loop: do ich=2,new
    if (frequency <= freq(ich)) then
      df  = frequency-freq(ich-1)
      df0 = freq(ich)-freq(ich-1)
      Emissivity_H = eh(ich-1) + (df*(eh(ich)-eh(ich-1))/df0)
      Emissivity_V = ev(ich-1) + (df*(ev(ich)-ev(ich-1))/df0)
      exit Channel_loop
    end if
  end do Channel_loop

end subroutine SnowEM_Default


subroutine Canopy_Optic(vlai,frequency,theta,esv,d,gv,gh,&
                        ssalb_v,ssalb_h,tau_v, tau_h, dtau_dvlai)


  REAL(fp) :: frequency,theta,d,vlai,ssalb_v,ssalb_h,tau_v,tau_h,gv, gh, mu
  COMPLEX(fp) :: ix,k0,kz0,kz1,rhc,rvc,esv,expval1,factt,factrvc,factrhc
  REAL(fp) :: rh,rv,th,tv
  REAL(fp), PARAMETER :: threshold = 0.999_fp
  ! Optional tangent-linear/adjoint support: d(tau)/d(vlai). The single-scatter
  ! albedo and asymmetry factor do not depend on vlai, so only the optical
  ! depth carries the LAI/vegetation sensitivity.
  REAL(fp), OPTIONAL, INTENT(OUT) :: dtau_dvlai

  mu = COS(theta)
  ix = CMPLX(ZERO, ONE, fp)

  k0  = CMPLX(TWOPI*frequency/300.0_fp, ZERO, fp)   ! 1/mm
  kz0 = k0*mu
  kz1 = k0*SQRT(esv - SIN(theta)**2)

  rhc = (kz0 - kz1)/(kz0 + kz1)
  rvc = (esv*kz0 - kz1)/(esv*kz0 + kz1)

  expval1 = EXP(-TWO*ix*kz1*d)
  factrvc = ONE-rvc**2*expval1
  factrhc = ONE-rhc**2*expval1
  factt   = FOUR*kz0*kz1*EXP(ix*(kz0-kz1)*d)

  rv = ABS(rvc*(ONE - expval1)/factrvc)**2
  rh = ABS(rhc*(ONE - expval1)/factrhc)**2

  th = ABS(factt/((kz1+kz0)**2*factrhc))**2
  tv = ABS(esv*factt/((kz1+esv*kz0)**2*factrvc))**2

  gv = POINT5
  gh = POINT5

  tau_v = POINT5*vlai*(TWO-tv-th)
  tau_h = tau_v

  ssalb_v = MIN((rv+rh)/(TWO-tv-th),threshold)
  ssalb_h = ssalb_v

  ! tau_v = tau_h = 0.5*vlai*(2-tv-th) is linear in vlai
  IF (PRESENT(dtau_dvlai)) dtau_dvlai = POINT5*(TWO-tv-th)

end subroutine Canopy_Optic


subroutine Snow_Optic(frequency,a,h,f,ep_real,ep_imag,gv,gh, ssalb_v,ssalb_h,tau_v,tau_h)


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

  ki1 = THREE*f*yi/fact2**2
  ki2 = kr2
  ki3 = ks3
  ki  = k**2/(TWO*kr)*(ki1+ki2*ki3)

  gv = POINT5
  gh = POINT5

  ssalb_v = MIN(ks/ki, 0.999_fp)
  ssalb_h = ssalb_v
  tau_v = TWO*ki*h
  tau_h = tau_v

end subroutine Snow_Optic


subroutine Soil_Diel(freq,t_soil,vmc,rhob,rhos,sand,clay,esm,desm_dvmc)


  REAL(fp) :: f,tauw,freq,t_soil,vmc,rhob,rhos,sand,clay
  REAL(fp) :: alpha,beta,ess,rhoef,t,eswi,eswo
  REAL(fp) :: esof
  COMPLEX(fp) :: esm,esw,es1,es2
  ! Optional tangent-linear/adjoint support: d(esm)/d(vmc) (complex)
  COMPLEX(fp), OPTIONAL, INTENT(OUT) :: desm_dvmc
  COMPLEX(fp) :: des1_dvmc, dinner_dvmc

  alpha = 0.65_fp
  beta  = 1.09_fp - 0.11_fp*sand + 0.18_fp*clay
  ess = (1.01_fp + 0.44_fp*rhos)**2 - 0.062_fp
  rhoef = -1.645_fp + 1.939_fp*rhob - 0.020213_fp*sand + 0.01594_fp*clay
  t = t_soil - 273.0_fp
  f = freq*1.0e9_fp

  ! the permittivity at the high frequency limit
  eswi = 5.5_fp

  ! the permittivity of free space (esof)
  esof = 8.854e-12_fp

  ! static dieletric constant (eswo)
  eswo = 87.134_fp+(-1.949e-1_fp+(-1.276e-2_fp+2.491e-4_fp*t)*t)*t
  tauw = 1.1109e-10_fp+(-3.824e-12_fp+(6.938e-14_fp-5.096e-16_fp*t)*t)*t

  if (vmc > ZERO) then
     es1 = CMPLX(eswi, -rhoef*(rhos-rhob)/(TWOPI*f*esof*rhos*vmc), fp)
  else
     es1 = CMPLX(eswi, ZERO, fp)
  endif

  es2 = CMPLX(eswo-eswi, ZERO, fp)/CMPLX(ONE, f*tauw, fp)
  esw = es1 + es2
  esm = ONE + (ess**alpha - ONE)*rhob/rhos + vmc**beta*esw**alpha - vmc

  ! Analytic d(esm)/d(vmc): esm above is the pre-power "inner" value. The water
  ! permittivity esw depends on vmc only through the imaginary part of es1
  ! (= -K/vmc), so d(esw)/d(vmc) = +K/vmc^2 (imaginary). es2 is independent of vmc.
  ! At vmc = 0 the vmc**(beta-1) term is unbounded for beta < 1 (soil types with
  ! high sand fraction), so treat the boundary as a zero-derivative region — the
  ! same convention as the emissivity/input clip handling.
  IF ( PRESENT(desm_dvmc) ) THEN
    IF ( vmc > ZERO ) THEN
      des1_dvmc = CMPLX(ZERO, rhoef*(rhos-rhob)/(TWOPI*f*esof*rhos*vmc*vmc), fp)
      dinner_dvmc = beta*vmc**(beta-ONE)*esw**alpha &
                  + vmc**beta*alpha*esw**(alpha-ONE)*des1_dvmc - ONE
      desm_dvmc   = (ONE/alpha)*esm**(ONE/alpha - ONE)*dinner_dvmc
    ELSE
      desm_dvmc = CMPLX(ZERO, ZERO, fp)
    END IF
  END IF

  esm = esm**(ONE/alpha)

  ! The forward clamps the imaginary part of esm to a constant when it is
  ! non-negative; in that branch the imaginary part of the derivative is zero.
  IF ( PRESENT(desm_dvmc) ) THEN
    IF ( AIMAG(esm) >= ZERO ) desm_dvmc = CMPLX(REAL(desm_dvmc,fp), ZERO, fp)
  END IF

  if(AIMAG(esm) >= ZERO) esm = CMPLX(REAL(esm,fp),-0.0001_fp, fp)

end subroutine Soil_Diel


subroutine Snow_Diel(frequency,ep_real,ep_imag,rad,frac,ep_eff)

!----------------------------------------------------------------------------------
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    Snow_Diel   compute dielectric constant of snow
!
!   prgmmr: Fuzhong Weng and Banghua Yan                 org: nesdis              date: 2000-11-28
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
!----------------------------------------------------------------------------------

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

  if (AIMAG(ep_eff) >= ZERO) ep_eff = CMPLX(REAL(ep_eff), -0.0001_fp, fp)

end subroutine Snow_Diel


subroutine Canopy_Diel(frequency,mg,esv,rhoveg)


  REAL(fp) :: frequency,  mg, en, vf, vb
  REAL(fp) :: rhoveg, vmv
  COMPLEX(fp) :: esv, xx

  vmv = mg*rhoveg/( ONE - mg*(ONE-rhoveg) )
  en = 1.7_fp + (3.2_fp + 6.5_fp*vmv)*vmv

  vf = vmv*(0.82_fp*vmv + 0.166_fp)
  vb = 31.4_fp*vmv*vmv/( ONE + 59.5_fp*vmv*vmv)

  xx = CMPLX(ZERO,ONE,fp)

  esv = en + vf*(4.9_fp + 75.0_fp/(ONE + xx*frequency/18.0_fp)-xx*(18.0_fp/frequency)) + &
        vb*(2.9_fp + 55.0_fp/(ONE + SQRT(xx*frequency/0.18_fp)))

  if (AIMAG(esv) >= ZERO) esv = CMPLX(REAL(esv), -0.0001_fp, fp)

end subroutine Canopy_Diel


subroutine Reflectance(em1, em2, theta_i, theta_t, rv, rh)


  REAL(fp) :: theta_i, theta_t
  REAL(fp) :: rh, rv,cos_i,cos_t
  COMPLEX(fp) :: em1, em2, m1, m2, angle_i, angle_t

  ! compute the refractive index ratio between medium 2 and 1
  ! using dielectric constant (n = SQRT(e))
  cos_i = COS(theta_i)
  cos_t = COS(theta_t)

  angle_i = CMPLX(cos_i, ZERO, fp)
  angle_t = CMPLX(cos_t, ZERO, fp)

  m1 = SQRT(em1)
  m2 = SQRT(em2)

  rv = (ABS((m1*angle_t-m2*angle_i)/(m1*angle_t+m2*angle_i)))**2
  rh = (ABS((m1*angle_i-m2*angle_t)/(m1*angle_i+m2*angle_t)))**2

end subroutine Reflectance


subroutine Transmittance(em1,em2,theta_i,theta_t,tv,th)


  REAL(fp) :: theta_i, theta_t
  REAL(fp) :: th, tv, rr, cos_i,cos_t
  COMPLEX(fp) :: em1, em2, m1, m2, angle_i, angle_t

  ! compute the refractive index ratio between medium 2 and 1
  ! using dielectric constant (n = SQRT(e))
  cos_i = COS(theta_i)
  cos_t = COS(theta_t)

  angle_i = CMPLX(cos_i, ZERO, fp)
  angle_t = CMPLX(cos_t, ZERO, fp)

  m1 = SQRT(em1)
  m2 = SQRT(em2)

  rr = ABS(m2/m1)*cos_t/cos_i
  tv = rr*(ABS(TWO*m1*angle_i/(m1*angle_t + m2*angle_i)))**2
  th = rr*(ABS(TWO*m1*angle_i/(m1*angle_i + m2*angle_t)))**2

end subroutine Transmittance


subroutine Roughness_Reflectance(frequency,sigma,rv,rh)


  REAL(fp) :: frequency
  REAL(fp) :: q, rh, rv, rh_s, rv_s, sigma

  rh_s = 0.3_fp*rh
  rv_s = 0.3_fp*rv
  q = 0.35_fp*(ONE - EXP(-0.60_fp*frequency*sigma**TWO))
  rh = rh_s + q*(rv_s-rh_s)
  rv = rv_s + q*(rh_s-rv_s)

end subroutine Roughness_Reflectance


subroutine Two_Stream_Solution(mu,gv,gh,ssalb_h,ssalb_v,tau_h,tau_v, &
      r21_h,r21_v,r23_h,r23_v,t21_v,t21_h,esv,esh,frequency,t_soil,t_skin, &
      desv_dtauv,desh_dtauh,desv_dr23v,desh_dr23h)


  REAL(fp) :: mu, gv, gh, ssalb_h, ssalb_v, tau_h,tau_v,                 &
              r21_h, r21_v, r23_h, r23_v, t21_v, t21_h, esv, esh
  REAL(fp) :: alfa_v, alfa_h, kk_h, kk_v, gamma_h, gamma_v, beta_v, beta_h
  REAL(fp) :: fact1,fact2
  REAL(fp) :: frequency, t_soil, t_skin
  REAL(fp) :: gsect0, gsect1_h, gsect1_v, gsect2_h, gsect2_v
  ! Optional tangent-linear/adjoint support: d(emissivity)/d(optical depth)
  ! (tau, LAI/vegetation) and d(emissivity)/d(lower-boundary reflectance r23,
  ! soil moisture).
  REAL(fp), OPTIONAL, INTENT(OUT) :: desv_dtauv, desh_dtauh
  REAL(fp), OPTIONAL, INTENT(OUT) :: desv_dr23v, desh_dr23h
  REAL(fp) :: num_h, den_h, num_v, den_v
  REAL(fp) :: dfact1_dtauh, dfact2_dtauv, dgsect2_h_dtauh, dgsect2_v_dtauv
  REAL(fp) :: dnum_h, dden_h, dnum_v, dden_v
  REAL(fp) :: e1h, e2h, e1v, e2v, dgamma_h, dgamma_v, dAcoef_h, dAcoef_v
  REAL(fp) :: dfact1_dr23h, dfact2_dr23v, dgsect1_dr23, dgsect2_h_dr23h, dgsect2_v_dr23v

  alfa_h  = SQRT((ONE - ssalb_h)/(ONE - gh*ssalb_h))
  kk_h    = SQRT((ONE - ssalb_h)*(ONE -  gh*ssalb_h))/mu
  beta_h  = (ONE - alfa_h)/(ONE + alfa_h)
  gamma_h = (beta_h -r23_h)/(ONE-beta_h*r23_h)

  alfa_v  = SQRT((ONE-ssalb_v)/(ONE - gv*ssalb_v))
  kk_v    = SQRT((ONE-ssalb_v)*(ONE - gv*ssalb_v))/mu
  beta_v  = (ONE - alfa_v)/(ONE + alfa_v)
  gamma_v = (beta_v -r23_v)/(ONE-beta_v*r23_v)

  fact1=gamma_h*EXP(-TWO*kk_h*tau_h)
  fact2=gamma_v*EXP(-TWO*kk_v*tau_v)

  gsect0  =(EXP(C_2*frequency/t_skin) -ONE)/(EXP(C_2*frequency/t_soil) -ONE)

  gsect1_h=(ONE-r23_h)*(gsect0-ONE)
  gsect2_h=((ONE-beta_h*beta_h)/(ONE-beta_h*r23_h))*EXP(-kk_h*tau_h)

  gsect1_v=(ONE-r23_v)*(gsect0-ONE)
  gsect2_v=((ONE-beta_v*beta_v)/(ONE-beta_v*r23_v))*EXP(-kk_h*tau_v)

  num_h = (ONE - beta_h)*(ONE + fact1) + gsect1_h*gsect2_h
  den_h = ONE-beta_h*r21_h-(beta_h-r21_h)*fact1
  num_v = (ONE - beta_v)*(ONE + fact2) + gsect1_v*gsect2_v
  den_v = ONE-beta_v*r21_v-(beta_v-r21_v)*fact2

  esh  = t21_h*num_h/den_h
  esv  = t21_v*num_v/den_v

  ! Tangent-linear/adjoint sensitivities of the (pre-clip) emissivities to the
  ! canopy/snow optical depths. Only tau_h and tau_v vary with LAI/vegetation;
  ! beta, kk, gamma, r21, r23, t21, gsect0 and gsect1 are independent of them.
  ! Note gsect2_v carries EXP(-kk_h*tau_v) in the forward, so its tau_v
  ! derivative uses kk_h (preserved here for exact consistency).
  IF ( PRESENT(desh_dtauh) .OR. PRESENT(desv_dtauv) ) THEN
    dfact1_dtauh    = -TWO*kk_h*fact1
    dgsect2_h_dtauh = -kk_h*gsect2_h
    dnum_h          = (ONE - beta_h)*dfact1_dtauh + gsect1_h*dgsect2_h_dtauh
    dden_h          = -(beta_h-r21_h)*dfact1_dtauh
    dfact2_dtauv    = -TWO*kk_v*fact2
    dgsect2_v_dtauv = -kk_h*gsect2_v
    dnum_v          = (ONE - beta_v)*dfact2_dtauv + gsect1_v*dgsect2_v_dtauv
    dden_v          = -(beta_v-r21_v)*dfact2_dtauv
    IF (PRESENT(desh_dtauh)) desh_dtauh = t21_h*(dnum_h*den_h - num_h*dden_h)/(den_h*den_h)
    IF (PRESENT(desv_dtauv)) desv_dtauv = t21_v*(dnum_v*den_v - num_v*dden_v)/(den_v*den_v)
  END IF

  ! Tangent-linear/adjoint sensitivities of the (pre-clip) emissivities to the
  ! lower-boundary reflectances r23 (the soil-moisture path enters here). r23_h
  ! affects gamma_h, gsect1_h and gsect2_h; likewise r23_v for the V component.
  IF ( PRESENT(desh_dr23h) .OR. PRESENT(desv_dr23v) ) THEN
    e1h = EXP(-TWO*kk_h*tau_h); e2h = EXP(-kk_h*tau_h)
    e1v = EXP(-TWO*kk_v*tau_v); e2v = EXP(-kk_h*tau_v)
    dgsect1_dr23 = -(gsect0-ONE)
    ! H component
    dgamma_h        = (beta_h*beta_h - ONE)/(ONE-beta_h*r23_h)**2
    dfact1_dr23h    = e1h*dgamma_h
    dAcoef_h        = (ONE-beta_h*beta_h)*beta_h/(ONE-beta_h*r23_h)**2
    dgsect2_h_dr23h = e2h*dAcoef_h
    dnum_h = (ONE-beta_h)*dfact1_dr23h + dgsect1_dr23*gsect2_h + gsect1_h*dgsect2_h_dr23h
    dden_h = -(beta_h-r21_h)*dfact1_dr23h
    ! V component
    dgamma_v        = (beta_v*beta_v - ONE)/(ONE-beta_v*r23_v)**2
    dfact2_dr23v    = e1v*dgamma_v
    dAcoef_v        = (ONE-beta_v*beta_v)*beta_v/(ONE-beta_v*r23_v)**2
    dgsect2_v_dr23v = e2v*dAcoef_v
    dnum_v = (ONE-beta_v)*dfact2_dr23v + dgsect1_dr23*gsect2_v + gsect1_v*dgsect2_v_dr23v
    dden_v = -(beta_v-r21_v)*dfact2_dr23v
    IF (PRESENT(desh_dr23h)) desh_dr23h = t21_h*(dnum_h*den_h - num_h*dden_h)/(den_h*den_h)
    IF (PRESENT(desv_dr23v)) desv_dr23v = t21_v*(dnum_v*den_v - num_v*dden_v)/(den_v*den_v)
  END IF

  if (esh < EMISSH_DEFAULT) then
    esh = EMISSH_DEFAULT
    if (PRESENT(desh_dtauh)) desh_dtauh = ZERO
    if (PRESENT(desh_dr23h)) desh_dr23h = ZERO
  end if
  if (esv < EMISSV_DEFAULT) then
    esv = EMISSV_DEFAULT
    if (PRESENT(desv_dtauv)) desv_dtauv = ZERO
    if (PRESENT(desv_dr23v)) desv_dr23v = ZERO
  end if

  if (esh > ONE) then
    esh = ONE
    if (PRESENT(desh_dtauh)) desh_dtauh = ZERO
    if (PRESENT(desh_dr23h)) desh_dr23h = ZERO
  end if
  if (esv > ONE) then
    esv = ONE
    if (PRESENT(desv_dtauv)) desv_dtauv = ZERO
    if (PRESENT(desv_dr23v)) desv_dr23v = ZERO
  end if

end subroutine Two_Stream_Solution

END MODULE NESDIS_LandEM_Module
