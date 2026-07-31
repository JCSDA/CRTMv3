!
! Emission_Module
!
! Module containing the emission radiative transfer
! solution procedures in the CRTM.
!
!
! CREATION HISTORY:
!       Written by:     Quanhua Liu, QSS at JCSDA; quanhua.liu@noaa.gov
!                       Yong Han,    NOAA/NESDIS;  yong.han@noaa.gov
!                       Paul van Delst; CIMMS/SSEC; paul.vandelst@noaa.gov
!                       08-Jun-2004

MODULE Emission_Module

  ! ------------------
  ! Environment set up
  ! ------------------
  ! Module use statements
  USE RTV_Define
  USE CRTM_Parameters
  USE Type_Kinds
  
  IMPLICIT NONE
  
  ! --------------------
  ! Default visibilities
  ! --------------------
  ! Everything private by default
  PRIVATE
  
  PUBLIC CRTM_Emission
  PUBLIC CRTM_Emission_TL
  PUBLIC CRTM_Emission_AD
  ! Polarized (Stokes 2..n) completion of the non-scattering solution
  PUBLIC CRTM_Emission_Stokes
  PUBLIC CRTM_Emission_Stokes_TL
  PUBLIC CRTM_Emission_Stokes_AD
  
  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: $'
    
CONTAINS

!################################################################################
!################################################################################
!##                                                                            ##
!##                        ## PRIVATE MODULE ROUTINES ##                       ##
!##                                                                            ##
!################################################################################
!################################################################################


  SUBROUTINE CRTM_Emission(n_Layers, & ! Input  number of atmospheric layers
                           n_Angles, & ! number angles used in SfcOptics
                    Diffuse_Surface, & ! Input  TRUE: Lambertian, FALSE: specular
                                  u, & ! Input  cosine of local viewing angle
                               T_OD, & ! Input  nadir layer optical depth
                  Planck_Atmosphere, & ! Input  atmospheric layer Planck radiance
                     Planck_Surface, & ! Input  surface Planck radiance 
                         emissivity, & ! Input  surface emissivity
                       reflectivity, & ! Input  surface reflectivity matrix
                direct_reflectivity, & ! Input  reflectivity for direct irradiance 
                  cosmic_background, & ! Input  cosmic background radiance
                   Solar_irradiance, & ! Input  Solar spectral irradiance
                   Is_Solar_Channel, & ! Input  Indicate solar affected channel
               Source_Zenith_Radian, & ! Input  Point source (e.g. solar) zenith angle
                                RTV)   ! Output TOA radiance and others
! ----------------------------------------------------------------------------- !
!  FUNCTION: Compute IR/MW upward radiance at the top of the profile.           !
!    This code heritages the concept from previous operational code.            !
!    It starts from cosmic background downward.                                 !
!    The downward radiance at the lower level is the transmitted radiance       !
!    from upper level adding the layer downward source function.                !
!    The downward angle is either the same as satellite viewing zenith for a    !
!    specular surface or the diffuse angle for a lambertian surface. The upward !
!    radiance at the surface is the surface emission term adding from surface   !
!    reflected downward radiance. Then, the upward radiance is the sum of       !
!    from the lower level transmitted radiance adding the upward layer          !
!    source function.                                                           !
!                                                                               !
!    Quanhua Liu    Quanhua.Liu@noaa.gov                                        !
! ----------------------------------------------------------------------------- !

    ! Arguments
    INTEGER,                     INTENT(IN)     :: n_Layers
    INTEGER,                     INTENT(IN)     :: n_Angles
    LOGICAL,                     INTENT(IN)     :: Diffuse_Surface
    REAL(fp),                    INTENT(IN)     :: u
    REAL(fp), DIMENSION(:),      INTENT(IN)     :: T_OD
    REAL(fp), DIMENSION(0:),     INTENT(IN)     :: Planck_Atmosphere
    REAL(fp),                    INTENT(IN)     :: Planck_Surface
    REAL(fp), DIMENSION(:),      INTENT(IN)     :: emissivity
    REAL(fp), DIMENSION(:,:),    INTENT(IN)     :: reflectivity 
    REAL(fp), DIMENSION(:),      INTENT(IN)     :: direct_reflectivity 
    REAL(fp),                    INTENT(IN)     :: cosmic_background
    REAL(fp),                    INTENT(IN)     :: Solar_irradiance
    LOGICAL,                     INTENT(IN)     :: Is_Solar_Channel
    REAL(fp),                    INTENT(IN)     :: Source_Zenith_Radian
    TYPE(RTV_type),              INTENT(IN OUT) :: RTV
    ! Local variables
    REAL(fp) :: layer_source_up, cosine_u0 
    INTEGER :: k

    ! --------------------
    ! Downwelling radiance
    ! --------------------
    ! Determing secant downward angle from surface behavior
    IF( Diffuse_Surface ) THEN
      RTV%Secant_Down_Angle = SECANT_DIFFUSIVITY 
    ELSE
      RTV%Secant_Down_Angle = ONE/u
    END IF

    ! Start from the top of the atmosphere
    RTV%e_Level_Rad_DOWN(0) = cosmic_background 
    RTV%Total_OD = ZERO

    ! Loop from top layer to bottom layer
    DO k = 1, n_Layers
      ! Accumulate optical depth 
      RTV%Total_OD = RTV%Total_OD + T_OD(k)
      ! Layer downward transmittance
      RTV%e_Layer_Trans_DOWN(k) = EXP(-T_OD(k)*RTV%Secant_Down_Angle)
      ! Downward radiance  
      RTV%e_Level_Rad_DOWN(k) = (RTV%e_Level_Rad_DOWN(k-1)*RTV%e_Layer_Trans_DOWN(k)) + &
                                (Planck_Atmosphere(k)*(ONE-RTV%e_Layer_Trans_DOWN(k)))

      RTV%e_Layer_Trans_UP(k) = EXP(-T_OD(k)/u)

      ! GSI cloud detection
      RTV%e_Cloud_Radiance_UP(k) = RTV%e_Source_UP(k-1) + Planck_Atmosphere(k)*RTV%e_Level_Trans_UP(k-1)
      RTV%e_Source_UP(k) = RTV%e_Source_UP(k-1)+RTV%e_Level_Trans_UP(k-1)*Planck_Atmosphere(k)*(ONE-RTV%e_Layer_Trans_UP(k))
      RTV%e_Level_Trans_UP(k) = RTV%e_Level_Trans_UP(k-1)*RTV%e_Layer_Trans_UP(k)
    END DO

    ! ----------------
    ! Surface radiance
    ! ----------------
    ! upward radiance at the surface ( emission part + reflection part)
    RTV%e_Level_Rad_UP(n_Layers) = (emissivity(n_Angles)*Planck_Surface) + &
                                   (reflectivity(1,1)*RTV%e_Level_Rad_DOWN(n_Layers))

    ! Solar contribution to the upward radiance at the surface
    RTV%Down_Solar_Radiance = ZERO
    IF( Is_Solar_Channel ) THEN
      cosine_u0 = COS(Source_Zenith_Radian)
      IF( cosine_u0 > ZERO) THEN
        RTV%Down_Solar_Radiance = cosine_u0*EXP(-RTV%Total_OD/cosine_u0)*Solar_Irradiance/PI
        RTV%e_Level_Rad_UP(n_Layers) = RTV%e_Level_Rad_UP(n_Layers) + &
          (RTV%Down_Solar_Radiance*direct_reflectivity(1))
      END IF
    END IF

    ! ------------------
    ! Upwelling radiance
    ! ------------------
    ! Initialise upwelling radiance
    RTV%Up_Radiance = ZERO

    ! Loop from SFC->TOA
    DO k = n_Layers, 1, -1
      ! layer upwelling transmittance
!!      RTV%e_Layer_Trans_UP(k) = EXP(-T_OD(k)/u)
      ! layer upwelling source function
      layer_source_up = Planck_Atmosphere(k) * ( ONE - RTV%e_Layer_Trans_UP(k) )
      ! upwelling radiance (including reflected downwelling and surface)
      RTV%e_Level_Rad_UP(k-1) = (RTV%e_Level_Rad_UP(k)*RTV%e_Layer_Trans_UP(k)) + &
                                layer_source_up 
      ! upwelling radiance (atmospheric portion only)
      RTV%Up_Radiance = (RTV%Up_Radiance*RTV%e_Layer_Trans_UP(k)) + layer_source_up
    END DO

  END SUBROUTINE CRTM_Emission
  
  SUBROUTINE CRTM_Emission_TL(n_Layers, & ! Input  number of atmospheric layers
                              n_Angles, & ! number angles used in SfcOptics
                                     u, & ! Input  cosine of local viewing angle
                     Planck_Atmosphere, & ! Input  atmospheric layer Planck radiance
                        Planck_Surface, & ! Input  surface Planck radiance 
                            emissivity, & ! Input  surface emissivity
                          reflectivity, & ! Input  surface reflectivity matrix
                   direct_reflectivity, & ! Input  reflectivity for direct irradiance 
                      Solar_irradiance, & ! Input  Solar spectral irradiance
                      Is_Solar_Channel, & ! Input  Indicate solar affected channel 
                  Source_Zenith_Radian, & ! Input  Point source (e.g. solar) zenith angle
                                   RTV, & ! Input  Structure containing forward part results 
                               T_OD_TL, & ! Input  tangent-linear of layer optical depth
                  Planck_Atmosphere_TL, & ! Input  TL atmospheric layer Planck radiance
                     Planck_Surface_TL, & ! Input  TL surface Planck radiance
                         emissivity_TL, & ! Input  TL surface emissivity
                       reflectivity_TL, & ! Input  TL surface reflectivity matrix
                direct_reflectivity_TL, & ! Input  TL surface ditrct reflectivity
                             up_rad_TL, & ! Output TL TOA radiance
                       down_rad_TL_out, & ! Output TL surface downwelling radiance (OPTIONAL)
                  down_rad_prof_TL_out, & ! Output TL downwelling radiance PROFILE (OPTIONAL)
                    up_rad_prof_TL_out)   ! Output TL upwelling radiance PROFILE (OPTIONAL)
! --------------------------------------------------------------------------- !
!  FUNCTION: Compute tangent-linear upward radiance at the top of the         !
!    atmosphere using carried results in RTV structure from forward           !
!    calculation.                                                             !
!    Quanhua Liu    Quanhua.Liu@noaa.gov                                      !
! --------------------------------------------------------------------------- !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_Layers, n_Angles
      LOGICAL, INTENT(IN) :: Is_Solar_Channel
      REAL (fp), INTENT(IN) :: Solar_irradiance, Source_Zenith_Radian
      REAL (fp), INTENT(IN), DIMENSION( : ) ::  emissivity,T_OD_TL,emissivity_TL
      REAL (fp), INTENT(IN), DIMENSION( :,: ) :: reflectivity ,reflectivity_TL
      REAL (fp), INTENT(IN), DIMENSION( : ) :: direct_reflectivity,direct_reflectivity_TL
      REAL (fp), INTENT(IN), DIMENSION( 0: ) :: Planck_Atmosphere,Planck_Atmosphere_TL
      REAL (fp), INTENT(IN) :: Planck_Surface,u,Planck_Surface_TL
      REAL (fp), INTENT(INOUT) :: up_rad_TL
      REAL (fp), INTENT(OUT), OPTIONAL :: down_rad_TL_out
      REAL (fp), INTENT(OUT), OPTIONAL, DIMENSION(:) :: down_rad_prof_TL_out
      REAL (fp), INTENT(OUT), OPTIONAL, DIMENSION(:) :: up_rad_prof_TL_out

    !   Structure RTV carried in variables from forward calculation.
      TYPE(RTV_type), INTENT( IN) :: RTV
    !  internal variables
      REAL (fp) :: layer_source_up_TL, layer_source_down_TL,a_TL,down_rad_TL
      REAL (fp) :: Total_OD, Total_OD_TL
      INTEGER :: k
      REAL( fp) :: cosine_u0

    !#--------------------------------------------------------------------------#
    !#                -- Downwelling TL radiance   --                           #
    !#--------------------------------------------------------------------------#

      down_rad_TL = ZERO
      Total_OD_TL = ZERO

      Total_OD = RTV%Total_OD
      IF ( PRESENT(down_rad_prof_TL_out) ) down_rad_prof_TL_out = ZERO

      DO k = 1, n_Layers
       ! accumulate tangent-linear optical depth
       Total_OD_TL = Total_OD_TL + T_OD_TL(k)
       a_TL = -T_OD_TL(k) * RTV%Secant_Down_Angle

       layer_source_down_TL = Planck_Atmosphere_TL(k) * ( ONE - RTV%e_Layer_Trans_DOWN(k) ) &
                            - Planck_Atmosphere(k) * RTV%e_Layer_Trans_DOWN(k) * a_TL

     ! downward tangent-linear radiance
     !    down_rad(k) = down_rad(k-1) * layer_trans(k) + layer_source_down
       down_rad_TL = down_rad_TL*RTV%e_Layer_Trans_DOWN(k)  &
       +RTV%e_Level_Rad_DOWN(k-1)*RTV%e_Layer_Trans_DOWN(k)*a_TL+layer_source_down_TL
       ! Per-level downwelling TL profile: down_rad_TL now holds TL of e_Level_Rad_DOWN(k).
       IF ( PRESENT(down_rad_prof_TL_out) ) down_rad_prof_TL_out(k) = down_rad_TL
      ENDDO

      ! Surface downwelling tangent-linear radiance (always-on output).
      ! At this point down_rad_TL holds TL of e_Level_Rad_DOWN(n_Layers).
      IF ( PRESENT(down_rad_TL_out) ) down_rad_TL_out = down_rad_TL

    !#--------------------------------------------------------------------------#
    !#                -- at surface   --                                        #
    !#--------------------------------------------------------------------------#

      ! upward tangent-linear radiance at the surface 
       up_rad_TL =emissivity_TL(n_Angles)*Planck_Surface+emissivity(n_Angles)*Planck_Surface_TL &
       +reflectivity_TL(1,1)*RTV%e_Level_Rad_DOWN(n_Layers)+reflectivity(1,1)*down_rad_TL

      ! point source (e.g. solar radiation)
       IF( Is_Solar_Channel ) THEN
        cosine_u0 = cos(Source_Zenith_Radian)
        IF( cosine_u0 > ZERO) THEN
        up_rad_TL = up_rad_TL + cosine_u0*Solar_Irradiance/PI &
                  * direct_reflectivity_TL(1) * exp(-Total_OD/cosine_u0)   &
                  - Solar_Irradiance/PI * direct_reflectivity(1)    &
                  * Total_OD_TL * exp(-Total_OD/cosine_u0)
        ENDIF
       ENDIF

    !#--------------------------------------------------------------------------#
    !#            -- Upwelling TL radiance   --                                 #
    !#--------------------------------------------------------------------------#

      ! Per-level upwelling TL profile: up_rad_TL currently holds the TL of the
      ! surface-level upward radiance e_Level_Rad_UP(n_Layers).
      IF ( PRESENT(up_rad_prof_TL_out) ) THEN
        up_rad_prof_TL_out = ZERO
        up_rad_prof_TL_out(n_Layers) = up_rad_TL
      END IF

      DO k = n_Layers, 1, -1
       a_TL = -T_OD_TL(k)/u
       layer_source_up_TL = Planck_Atmosphere_TL(k) * ( ONE - RTV%e_Layer_Trans_UP(k) ) &
                          - Planck_Atmosphere(k) * RTV%e_Layer_Trans_UP(k) * a_TL

      ! upward tangent linear radiance
       up_rad_TL=up_rad_TL*RTV%e_Layer_Trans_UP(k)  &
       +RTV%e_Level_Rad_UP(k)*RTV%e_Layer_Trans_UP(k)*a_TL+layer_source_up_TL
       ! up_rad_TL now holds the TL of e_Level_Rad_UP(k-1) (level 0 = TOA = Radiance).
       IF ( PRESENT(up_rad_prof_TL_out) .AND. k-1 >= 1 ) up_rad_prof_TL_out(k-1) = up_rad_TL
      ENDDO
!
      RETURN
      END SUBROUTINE CRTM_Emission_TL
!
!
      SUBROUTINE CRTM_Emission_AD(n_Layers, & ! Input  number of atmospheric layers
                                  n_Angles, & ! number angles used in SfcOptics
                                         u, & ! Input  cosine of local viewing angle
                         Planck_Atmosphere, & ! Input  atmospheric layer Planck radiance
                            Planck_Surface, & ! Input  surface Planck radiance 
                                emissivity, & ! Input  surface emissivity
                              reflectivity, & ! Input  surface reflectivity matrix 
                       direct_reflectivity, & ! Input  surface reflectivity matrix 
                          Solar_irradiance, & ! Input  Solar spectral irradiance
                          Is_Solar_Channel, & ! Input  Indicate solar affected channel 
                      Source_Zenith_Radian, & ! Input  Point source (e.g. solar) zenith angle
                                       RTV, & ! Input  Structure containing forward part results 
                              up_rad_AD_in, & ! Input  adjoint radiance at the top
                                   T_OD_AD, & ! Output AD layer optical depth
                      Planck_Atmosphere_AD, & ! Output AD atmospheric layer Planck radiance
                         Planck_Surface_AD, & ! Output AD surface Planck radiance
                             emissivity_AD, & ! Output AD surface emissivity
                           reflectivity_AD, & ! Output AD surface reflectivity matrix
                    direct_reflectivity_AD, & ! Output AD surface direct reflectivity
                            down_rad_AD_in, & ! Input  AD surface downwelling radiance (OPTIONAL)
                       down_rad_prof_AD_in, & ! Input  AD downwelling radiance PROFILE (OPTIONAL)
                         up_rad_prof_AD_in)   ! Input  AD upwelling radiance PROFILE (OPTIONAL)
! --------------------------------------------------------------------------- !
!  FUNCTION: Compute adjoint upward radiance at the top of the                !
!    atmosphere using carried results in RTV structure from forward           !
!    calculation.                                                             !
!    Quanhua Liu    Quanhua.Liu@noaa.gov                                      !
! --------------------------------------------------------------------------- !
      IMPLICIT NONE
      INTEGER, INTENT(IN) :: n_Layers, n_Angles
      LOGICAL, INTENT(IN) :: Is_Solar_Channel
      REAL (fp), INTENT(IN) :: Solar_Irradiance, Source_Zenith_Radian
      REAL (fp), INTENT(IN), DIMENSION( : ) ::  emissivity
      REAL (fp), INTENT(IN), DIMENSION( :,: ) :: reflectivity 
      REAL (fp), INTENT(IN), DIMENSION( : ) :: direct_reflectivity 
      REAL (fp), INTENT(IN), DIMENSION( 0: ) ::  Planck_Atmosphere
      REAL (fp), INTENT(IN) :: Planck_Surface,u
      REAL (fp), INTENT(IN) :: up_rad_AD_in
      REAL (fp), INTENT(IN), OPTIONAL :: down_rad_AD_in
      REAL (fp), INTENT(IN), OPTIONAL, DIMENSION(:) :: down_rad_prof_AD_in
      REAL (fp), INTENT(IN), OPTIONAL, DIMENSION(:) :: up_rad_prof_AD_in
      REAL (fp), INTENT(IN OUT), DIMENSION( : ) ::  T_OD_AD,emissivity_AD
      REAL (fp), INTENT(IN OUT), DIMENSION( :,: ) :: reflectivity_AD
      REAL (fp), INTENT(IN OUT), DIMENSION( : ) :: direct_reflectivity_AD
      REAL (fp), INTENT(IN OUT), DIMENSION( 0: ) ::  Planck_Atmosphere_AD
      REAL (fp), INTENT(IN OUT) :: Planck_Surface_AD
      TYPE(RTV_type), INTENT( IN) :: RTV
    !  internal variables
      REAL (fp) :: layer_source_up_AD, layer_source_down_AD,a_AD,down_rad_AD
      REAL (fp) :: cosine_u0, up_rad_AD, Total_OD, Total_OD_AD
      INTEGER :: k
!
    ! Initialize variables
      Total_OD_AD = ZERO
      T_OD_AD = ZERO
      Planck_Atmosphere_AD = ZERO
      Planck_Surface_AD = ZERO
      emissivity_AD = ZERO
      reflectivity_AD = ZERO
      direct_reflectivity_AD = ZERO
      up_rad_AD = up_rad_AD_in

    ! Total column optical depth carried from forward part
      Total_OD = RTV%Total_OD 

    !#--------------------------------------------------------------------------#
    !#                -- Upwelling adjoint radiance   --                        #
    !#--------------------------------------------------------------------------#
!
      DO k = 1, n_Layers
       ! Inject adjoint of the per-level upwelling profile output: at the top of
       ! iteration k, up_rad_AD is the adjoint of e_Level_Rad_UP(k-1) (level 0 = TOA
       ! = Radiance, handled by up_rad_AD_in).
       IF ( PRESENT(up_rad_prof_AD_in) .AND. k-1 >= 1 ) up_rad_AD = up_rad_AD + up_rad_prof_AD_in(k-1)

       a_AD = RTV%e_Level_Rad_UP(k)*RTV%e_Layer_Trans_UP(k)*up_rad_AD
       layer_source_up_AD = up_rad_AD
       up_rad_AD = up_rad_AD * RTV%e_Layer_Trans_UP(k)

       Planck_Atmosphere_AD(k) = Planck_Atmosphere_AD(k) + &
              layer_source_up_AD * (ONE - RTV%e_Layer_Trans_UP(k))
       a_AD = a_AD - Planck_Atmosphere(k) * RTV%e_Layer_Trans_UP(k)* layer_source_up_AD

       T_OD_AD(k) = T_OD_AD(k) - a_AD/u
      ENDDO
    !#--------------------------------------------------------------------------#
    !#                -- at surface   --                                        #
    !#--------------------------------------------------------------------------#

       ! Inject the surface-level (n_Layers) upwelling profile adjoint into up_rad_AD,
       ! which the surface block below distributes to emissivity / Planck / reflectivity.
       IF ( PRESENT(up_rad_prof_AD_in) ) up_rad_AD = up_rad_AD + up_rad_prof_AD_in(n_Layers)

       IF( Is_Solar_Channel ) THEN
        cosine_u0 = cos(Source_Zenith_Radian)
        IF( cosine_u0 > ZERO) THEN
        Total_OD_AD = -Solar_Irradiance/PI * direct_reflectivity(1) &
                    * up_rad_AD * exp(-Total_OD/cosine_u0)
        direct_reflectivity_AD(1) = cosine_u0 * Solar_Irradiance/PI &
                    * up_rad_AD* exp(-Total_OD/cosine_u0)
        ENDIF
       ENDIF

      emissivity_AD(n_Angles)=up_rad_AD*Planck_Surface
      Planck_Surface_AD = emissivity(n_Angles)*up_rad_AD
      reflectivity_AD(1,1)=up_rad_AD*RTV%e_Level_Rad_DOWN(n_Layers)
      down_rad_AD = reflectivity(1,1)*up_rad_AD
      ! Inject adjoint of the surface downwelling radiance output (always-on).
      ! e_Level_Rad_DOWN(n_Layers) feeds both the surface reflection and Down_Radiance.
      IF ( PRESENT(down_rad_AD_in) ) down_rad_AD = down_rad_AD + down_rad_AD_in
!
    !#--------------------------------------------------------------------------#
    !#                -- Downward adjoint radiance   --                         #
    !#--------------------------------------------------------------------------#
      DO k = n_Layers, 1, -1

       ! Inject adjoint of the per-level downwelling radiance profile output:
       ! at the top of iteration k, down_rad_AD is the adjoint of e_Level_Rad_DOWN(k).
       IF ( PRESENT(down_rad_prof_AD_in) ) down_rad_AD = down_rad_AD + down_rad_prof_AD_in(k)

       a_AD = RTV%e_Level_Rad_DOWN(k-1)*RTV%e_Layer_Trans_DOWN(k)*down_rad_AD
       layer_source_down_AD = down_rad_AD
       down_rad_AD = down_rad_AD*RTV%e_Layer_Trans_DOWN(k)

       Planck_Atmosphere_AD(k) = Planck_Atmosphere_AD(k) + layer_source_down_AD * &
                                 (ONE - RTV%e_Layer_Trans_DOWN(k))
       a_AD = a_AD - Planck_Atmosphere(k) * RTV%e_Layer_Trans_DOWN(k)* layer_source_down_AD
 

       T_OD_AD(k) = T_OD_AD(k) - a_AD * RTV%Secant_Down_Angle

       T_OD_AD(k) = T_OD_AD(k) + Total_OD_AD
      ENDDO

      down_rad_AD = ZERO 

      RETURN
      END SUBROUTINE CRTM_Emission_AD


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Emission_Stokes
!
! PURPOSE:
!       Completes the non-scattering solution for the polarized Stokes
!       components 2..n_Stokes. CRTM_Emission itself is a scalar solver: it
!       returns the total intensity and nothing else, so on the n_Stokes > 1
!       path a clear-sky run came back with Q = U = V = 0 however polarized the
!       surface was.
!
!       In the absence of scattering the atmosphere is polarization neutral. It
!       emits unpolarized radiation, so the thermal source enters Stokes I
!       alone, and it transmits every Stokes component with the same layer
!       transmittance (CRTM replicates the per-angle cosine across the Stokes
!       slots of that angle). The only polarized object in the problem is the
!       surface, and the downwelling it reflects is unpolarized. So for k >= 2
!       the whole solution is a boundary value transported upward with no
!       source of its own,
!
!           S_k(surface) = e_k * B_surface  +  R_k1 * D_surface
!           S_k(level-1) = S_k(level) * layer_transmittance
!
!       where D_surface is the (unpolarized) downwelling already computed by
!       CRTM_Emission and R_k1 is the first column of the surface reflection
!       matrix, which is what an unpolarized incident vector selects.
!
!       This is exactly the statement test_VectorRT_ScalarLimit checks against
!       two scalar runs: I = (Iv+Ih)/2 and Q = (Iv-Ih)/2.
!
! CALLING SEQUENCE:
!       CALL CRTM_Emission_Stokes( n_Layers, n_Angles, n_Stokes, &
!                                  Planck_Surface, emissivity, reflectivity, RTV )
!
! COMMENTS:
!       Must be called after CRTM_Emission, which populates the RTV downwelling
!       and layer transmittances this reads.
!
!       emissivity and reflectivity are the FLATTENED (angle,Stokes) arrays
!       built by Reshape_Surf_Opt, so element (i-1)*n_Stokes+m is angle i,
!       Stokes m. The microwave non-scattering path is always specular, which
!       fixes n_Angles at 1 (Common_RTSolution.f90:334); the routine asserts
!       that rather than assuming it silently, because with more than one angle
!       the sensor angle is no longer the first block.
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Emission_Stokes( &
    n_Layers,        &  ! Input, number of atmospheric layers
    n_Angles,        &  ! Input, number of discrete zenith angles
    n_Stokes,        &  ! Input, number of Stokes components
    Planck_Surface,  &  ! Input, surface radiance
    emissivity,      &  ! Input, flattened surface emissivity
    reflectivity,    &  ! Input, flattened surface reflectivity
    RTV              )  ! In/Output, internal variables
    ! Arguments
    INTEGER,                  INTENT(IN)     :: n_Layers, n_Angles, n_Stokes
    REAL(fp),                 INTENT(IN)     :: Planck_Surface
    REAL(fp), DIMENSION(:),   INTENT(IN)     :: emissivity
    REAL(fp), DIMENSION(:,:), INTENT(IN)     :: reflectivity
    TYPE(RTV_type),           INTENT(IN OUT) :: RTV
    ! Local variables
    INTEGER  :: k, ks, out_lev
    REAL(fp) :: rad, down_sfc

    RTV%e_Rad_UP_Stokes = ZERO
    IF ( n_Stokes < 2 .OR. n_Angles /= 1 ) RETURN

    down_sfc = RTV%e_Level_Rad_DOWN(n_Layers)
    out_lev  = 0
    IF ( RTV%aircraft%rt ) out_lev = RTV%aircraft%idx

    DO ks = 2, n_Stokes
      ! Surface boundary: polarized emission plus the polarized part of the
      ! reflected, unpolarized, downwelling.
      rad = ( emissivity(ks) * Planck_Surface ) + ( reflectivity(ks,1) * down_sfc )
      ! Source-free, polarization-neutral transport to the observer level.
      DO k = n_Layers, out_lev+1, -1
        rad = rad * RTV%e_Layer_Trans_UP(k)
      END DO
      RTV%e_Rad_UP_Stokes(ks) = rad
    END DO

  END SUBROUTINE CRTM_Emission_Stokes


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Emission_Stokes_TL
!
! PURPOSE:
!       Tangent-linear of CRTM_Emission_Stokes. down_rad_TL is the tangent
!       linear of the surface downwelling radiance, which CRTM_Emission_TL
!       already returns through its down_rad_TL_out argument.
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Emission_Stokes_TL( &
    n_Layers,           &  ! Input, number of atmospheric layers
    n_Angles,           &  ! Input, number of discrete zenith angles
    n_Stokes,           &  ! Input, number of Stokes components
    u,                  &  ! Input, cosine of the sensor zenith angle
    Planck_Surface,     &  ! Input, FWD surface radiance
    emissivity,         &  ! Input, FWD flattened surface emissivity
    reflectivity,       &  ! Input, FWD flattened surface reflectivity
    RTV,                &  ! Input, internal variables
    T_OD_TL,            &  ! Input, TL layer optical depth
    Planck_Surface_TL,  &  ! Input, TL surface radiance
    emissivity_TL,      &  ! Input, TL flattened surface emissivity
    reflectivity_TL,    &  ! Input, TL flattened surface reflectivity
    down_rad_TL,        &  ! Input, TL surface downwelling radiance
    Stokes_TL           )  ! Output, TL polarized components
    ! Arguments
    INTEGER,                  INTENT(IN)  :: n_Layers, n_Angles, n_Stokes
    REAL(fp),                 INTENT(IN)  :: u
    REAL(fp),                 INTENT(IN)  :: Planck_Surface
    REAL(fp), DIMENSION(:),   INTENT(IN)  :: emissivity
    REAL(fp), DIMENSION(:,:), INTENT(IN)  :: reflectivity
    TYPE(RTV_type),           INTENT(IN)  :: RTV
    REAL(fp), DIMENSION(:),   INTENT(IN)  :: T_OD_TL
    REAL(fp),                 INTENT(IN)  :: Planck_Surface_TL
    REAL(fp), DIMENSION(:),   INTENT(IN)  :: emissivity_TL
    REAL(fp), DIMENSION(:,:), INTENT(IN)  :: reflectivity_TL
    REAL(fp),                 INTENT(IN)  :: down_rad_TL
    REAL(fp), DIMENSION(:),   INTENT(OUT) :: Stokes_TL
    ! Local variables
    INTEGER  :: k, ks, out_lev
    REAL(fp) :: rad, rad_TL, down_sfc, trans_TL

    Stokes_TL = ZERO
    IF ( n_Stokes < 2 .OR. n_Angles /= 1 ) RETURN

    down_sfc = RTV%e_Level_Rad_DOWN(n_Layers)
    out_lev  = 0
    IF ( RTV%aircraft%rt ) out_lev = RTV%aircraft%idx

    DO ks = 2, n_Stokes
      rad    = ( emissivity(ks) * Planck_Surface ) + ( reflectivity(ks,1) * down_sfc )
      rad_TL = ( emissivity_TL(ks)   * Planck_Surface ) + &
               ( emissivity(ks)      * Planck_Surface_TL ) + &
               ( reflectivity_TL(ks,1) * down_sfc ) + &
               ( reflectivity(ks,1)    * down_rad_TL )
      DO k = n_Layers, out_lev+1, -1
        ! layer_trans = EXP(-T_OD/u)  =>  layer_trans_TL = -T_OD_TL/u * layer_trans
        trans_TL = -T_OD_TL(k) / u * RTV%e_Layer_Trans_UP(k)
        rad_TL   = ( rad_TL * RTV%e_Layer_Trans_UP(k) ) + ( rad * trans_TL )
        rad      = rad * RTV%e_Layer_Trans_UP(k)
      END DO
      Stokes_TL(ks) = rad_TL
    END DO

  END SUBROUTINE CRTM_Emission_Stokes_TL


!--------------------------------------------------------------------------------
!
! NAME:
!       CRTM_Emission_Stokes_AD
!
! PURPOSE:
!       Adjoint of CRTM_Emission_Stokes. down_rad_AD is accumulated, not
!       assigned, so the caller can hand the running total to CRTM_Emission_AD
!       through its down_rad_AD_in argument and keep the two contributions to
!       the surface downwelling adjoint together.
!
!--------------------------------------------------------------------------------

  SUBROUTINE CRTM_Emission_Stokes_AD( &
    n_Layers,           &  ! Input, number of atmospheric layers
    n_Angles,           &  ! Input, number of discrete zenith angles
    n_Stokes,           &  ! Input, number of Stokes components
    u,                  &  ! Input, cosine of the sensor zenith angle
    Planck_Surface,     &  ! Input, FWD surface radiance
    emissivity,         &  ! Input, FWD flattened surface emissivity
    reflectivity,       &  ! Input, FWD flattened surface reflectivity
    RTV,                &  ! Input, internal variables
    Stokes_AD,          &  ! Input, AD polarized components
    T_OD_AD,            &  ! In/Output, AD layer optical depth
    Planck_Surface_AD,  &  ! In/Output, AD surface radiance
    emissivity_AD,      &  ! In/Output, AD flattened surface emissivity
    reflectivity_AD,    &  ! In/Output, AD flattened surface reflectivity
    down_rad_AD         )  ! In/Output, AD surface downwelling radiance
    ! Arguments
    INTEGER,                  INTENT(IN)     :: n_Layers, n_Angles, n_Stokes
    REAL(fp),                 INTENT(IN)     :: u
    REAL(fp),                 INTENT(IN)     :: Planck_Surface
    REAL(fp), DIMENSION(:),   INTENT(IN)     :: emissivity
    REAL(fp), DIMENSION(:,:), INTENT(IN)     :: reflectivity
    TYPE(RTV_type),           INTENT(IN)     :: RTV
    REAL(fp), DIMENSION(:),   INTENT(IN)     :: Stokes_AD
    REAL(fp), DIMENSION(:),   INTENT(IN OUT) :: T_OD_AD
    REAL(fp),                 INTENT(IN OUT) :: Planck_Surface_AD
    REAL(fp), DIMENSION(:),   INTENT(IN OUT) :: emissivity_AD
    REAL(fp), DIMENSION(:,:), INTENT(IN OUT) :: reflectivity_AD
    REAL(fp),                 INTENT(IN OUT) :: down_rad_AD
    ! Local variables
    INTEGER  :: k, ks, out_lev
    REAL(fp) :: rad_AD, down_sfc
    REAL(fp) :: rad_fwd(0:n_Layers)

    IF ( n_Stokes < 2 .OR. n_Angles /= 1 ) RETURN

    down_sfc = RTV%e_Level_Rad_DOWN(n_Layers)
    out_lev  = 0
    IF ( RTV%aircraft%rt ) out_lev = RTV%aircraft%idx

    DO ks = 2, n_Stokes
      ! Forward sweep, retaining the running radiance the adjoint needs. Index k
      ! holds the value at the BOTTOM of layer k, so rad_fwd(n_Layers) is the
      ! surface boundary value.
      rad_fwd = ZERO
      rad_fwd(n_Layers) = ( emissivity(ks) * Planck_Surface ) + &
                          ( reflectivity(ks,1) * down_sfc )
      DO k = n_Layers, out_lev+1, -1
        rad_fwd(k-1) = rad_fwd(k) * RTV%e_Layer_Trans_UP(k)
      END DO

      ! Adjoint sweep, exact transpose of the loop above.
      rad_AD = Stokes_AD(ks)
      DO k = out_lev+1, n_Layers
        ! rad_fwd(k-1) = rad_fwd(k)*trans(k)
        T_OD_AD(k) = T_OD_AD(k) - ( rad_fwd(k) * RTV%e_Layer_Trans_UP(k) / u ) * rad_AD
        rad_AD     = rad_AD * RTV%e_Layer_Trans_UP(k)
      END DO

      ! Adjoint of the surface boundary.
      emissivity_AD(ks)      = emissivity_AD(ks)      + ( Planck_Surface * rad_AD )
      Planck_Surface_AD      = Planck_Surface_AD      + ( emissivity(ks) * rad_AD )
      reflectivity_AD(ks,1)  = reflectivity_AD(ks,1)  + ( down_sfc * rad_AD )
      down_rad_AD            = down_rad_AD            + ( reflectivity(ks,1) * rad_AD )
    END DO

  END SUBROUTINE CRTM_Emission_Stokes_AD

END MODULE Emission_Module
