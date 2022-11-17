!----------------------------------------------------------------------------------
!:sdoc+:
!
! RTTOV_FASTEM_MODULE
!
! Module to provide a general interface to RTTOV FASTEM modules
! 
! This module implementes the FWD, TL and AD functions for variational data 
! assimilation application and retrieval applications. 
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 28-Feb-2018
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!
MODULE RTTOV_FASTEM_MODULE

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE RTTOV_FASTEM6_MODULE,      ONLY: rttov_fastem6 
  USE RTTOV_FASTEM6_TL_MODULE,   ONLY: rttov_fastem6_tl
  USE RTTOV_FASTEM6_AD_MODULE,   ONLY: rttov_fastem6_ad
  
  USE RTTOV_FASTEM5R1_MODULE,    ONLY: rttov_fastem5r1 
  USE RTTOV_FASTEM5R1_TL_MODULE, ONLY: rttov_fastem5r1_tl
  USE RTTOV_FASTEM5R1_AD_MODULE, ONLY: rttov_fastem5r1_ad
 
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  ! Data types
  PUBLIC :: iVar_type
  ! Science routines
  PUBLIC :: Compute_RTTOV_Fastem
  PUBLIC :: Compute_RTTOV_Fastem_TL
  PUBLIC :: Compute_RTTOV_Fastem_AD

  REAL(fp), PARAMETER :: ZERO   = 0.0_fp, ONE = 1.0_fp
  REAL(fp), PARAMETER :: PI     = 3.141592653589793238462643_fp
  
  ! --------------------------------------
  ! Structure definition to hold internal
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE

    ! Forward model input values
    INTEGER  :: fastem_version      = 6
    REAL(fp) :: Frequency           = ZERO
    REAL(fp) :: Zenith_Angle        = ZERO
    REAL(fp) :: Temperature         = ZERO
    REAL(fp) :: Salinity            = ZERO
    REAL(fp) :: Wind_Speed          = ZERO

    REAL(fp) :: Rel_Azimuth         = ZERO
    REAL(fp) :: Transmittance       = ZERO
       
   ! ...Optional
    LOGICAL  :: Supply_Foam_Fraction = .FALSE.
    REAL(fp) :: Foam_Fraction       = ZERO

  END TYPE iVar_type

CONTAINS

  SUBROUTINE Compute_RTTOV_Fastem( &
             fastem_version      , &
             Frequency           , &         ! Input
             Zenith_Angle        , &         ! Input
             Temperature         , &         ! Input
             Salinity            , &         ! Input
             Wind_Speed          , &         ! Input
             Emissivity          , &         ! Output
             Reflectivity        , &         ! Output
             Transmittance       , &         ! Input, may not be used
             Rel_Azimuth         , &         ! Input, may not be used
             Supply_Foam_Fraction, &         ! Optional input
             Foam_Fraction       , &          ! Optional input
             iVar)                            ! Optional input
    ! Arguments
    INTEGER,        INTENT(IN)            :: fastem_version
    REAL(fp),       INTENT(IN)            :: Frequency
    REAL(fp),       INTENT(IN)            :: Zenith_Angle
    REAL(fp),       INTENT(IN)            :: Temperature
    REAL(fp),       INTENT(IN)            :: Salinity
    REAL(fp),       INTENT(IN)            :: Wind_Speed
    REAL(fp),       INTENT(OUT)           :: Emissivity(4), Reflectivity(4)
    REAL(fp),       INTENT(IN)            :: Transmittance
    REAL(fp),       INTENT(IN)            :: Rel_Azimuth
    LOGICAL,        OPTIONAL, INTENT(IN)  :: Supply_Foam_Fraction
    REAL(fp),       OPTIONAL, INTENT(IN)  :: Foam_Fraction
    TYPE(iVar_type),OPTIONAL, INTENT(OUT) :: iVar

    IF(fastem_version == 6) THEN

       CALL rttov_fastem6(fastem_version, &        ! Input
                          Frequency   , &          ! Input
                          Zenith_Angle, &          ! Input
                          Temperature , &          ! Input
                          Salinity    , &          ! Input
                          Wind_Speed  , &          ! Input
                          Emissivity  , &          ! Output
                          Reflectivity, &          ! Output
                          Transmittance,&          ! Input, may not be used
                          Rel_Azimuth  ,&          ! Input, may not be used
                          Supply_Foam_Fraction, &  ! Optional input
                          Foam_Fraction)           ! Optional input
    ELSE
      CALL rttov_fastem5r1(fastem_version, &       ! Input
                           Frequency   , &         ! Input
                           Zenith_Angle, &         ! Input
                           Temperature , &         ! Input
                           Salinity    , &         ! Input
                           Wind_Speed  , &         ! Input
                           Emissivity  , &         ! Output
                           Reflectivity, &         ! Output
                           Transmittance,&         ! Input, may not be used
                           Rel_Azimuth  ,&         ! Input, may not be used
                           Supply_Foam_Fraction, & ! Optional input
                           Foam_Fraction)          ! Optional input

    END IF
    
    ! ...Save forward input variables for TL and AD calculations
    IF(PRESENT(ivar)) THEN
      iVar%fastem_version =  fastem_version
      iVar%Frequency      =  Frequency
      iVar%Zenith_Angle   =  Zenith_Angle
      iVar%Temperature    =  Temperature
      iVar%Salinity       =  Salinity
      iVar%Wind_Speed     =  Wind_Speed
      iVar%Rel_Azimuth    =  Rel_Azimuth
      iVar%Transmittance  =  Transmittance
    
      IF(PRESENT(Supply_Foam_Fraction)) THEN
        iVar%Supply_Foam_Fraction = Supply_Foam_Fraction
        IF(PRESENT(Foam_Fraction))iVar%Foam_Fraction = Foam_Fraction
      ENDIF
    ENDIF 
  END SUBROUTINE Compute_RTTOV_Fastem
  
  SUBROUTINE Compute_RTTOV_Fastem_TL( &
             Temperature_tl,       &  ! Input
             Salinity_tl ,         &  ! Input
             Wind_Speed_tl,        &  ! Input
             Emissivity_tl,        &  ! Output
             Reflectivity_tl,      &  ! Output
             Transmittance_tl,     &  ! Input, may not be used
             Rel_Azimuth_tl,       &  ! Input, may not be used
             Foam_Fraction_tl,     &  ! Input, may not be used
             iVar)                    ! Optional input
    ! Arguments
    REAL(fp), INTENT(IN)  :: Temperature_tl
    REAL(fp), INTENT(IN)  :: Salinity_tl
    REAL(fp), INTENT(IN)  :: Wind_Speed_tl
    REAL(fp), INTENT(OUT) :: Emissivity_tl(4)
    REAL(fp), INTENT(OUT) :: Reflectivity_tl(4)
    REAL(fp), INTENT(IN)  :: Transmittance_tl
    REAL(fp), INTENT(IN)  :: Rel_Azimuth_tl
    REAL(fp), INTENT(IN)  :: Foam_Fraction_tl
    TYPE(iVar_type), INTENT(IN) :: iVar

    REAL(fp) :: Emissivity(4), Reflectivity(4)

    IF(ivar%fastem_version == 6) THEN
                          
       CALL rttov_fastem6_tl(ivar%fastem_version,       &  ! Input
                             ivar%Frequency   ,         &  ! Input
                             ivar%Zenith_Angle,         &  ! Input
                             ivar%Temperature ,         &  ! Input
                             ivar%Salinity    ,         &  ! Input
                             ivar%Wind_Speed  ,         &  ! Input
                             Temperature_tl,            &  ! Input
                             Salinity_tl ,              &  ! Input
                             Wind_Speed_tl,             &  ! Input
                             Emissivity,                &  ! Output
                             Reflectivity,              &  ! Output
                             Emissivity_tl,             &  ! Output
                             Reflectivity_tl,           &  ! Output
                             ivar%Transmittance,        &  ! Input, may not be used
                             ivar%Rel_Azimuth,          &  ! Input, may not be used
                             Transmittance_tl,          &  ! Input, may not be used
                             Rel_Azimuth_tl,            &  ! Input, may not be used
                             ivar%Supply_Foam_Fraction, &  ! Optional input
                             ivar%Foam_Fraction,        &  ! Optional input
                             Foam_Fraction_tl)             ! Optional input
    ELSE

       CALL rttov_fastem5r1_tl(ivar%fastem_version,       &  ! Input
                             ivar%Frequency   ,         &  ! Input
                             ivar%Zenith_Angle,         &  ! Input
                             ivar%Temperature ,         &  ! Input
                             ivar%Salinity    ,         &  ! Input
                             ivar%Wind_Speed  ,         &  ! Input
                             Temperature_tl,            &  ! Input
                             Salinity_tl ,              &  ! Input
                             Wind_Speed_tl,             &  ! Input
                             Emissivity,                &  ! Output
                             Reflectivity,              &  ! Output
                             Emissivity_tl,             &  ! Output
                             Reflectivity_tl,           &  ! Output
                             ivar%Transmittance,        &  ! Input, may not be used
                             ivar%Rel_Azimuth,          &  ! Input, may not be used
                             Transmittance_tl,          &  ! Input, may not be used
                             Rel_Azimuth_tl,            &  ! Input, may not be used
                             ivar%Supply_Foam_Fraction, &  ! Optional input
                             ivar%Foam_Fraction,        &  ! Optional input
                             Foam_Fraction_tl)             ! Optional input

    END IF
    
  END SUBROUTINE Compute_RTTOV_Fastem_TL
  
  
  SUBROUTINE Compute_RTTOV_Fastem_AD( &
             Emissivity_ad,        &  ! Output
             Reflectivity_ad,      &  ! Output
             Temperature_ad,       &  ! Input
             Salinity_ad ,         &  ! Input
             Wind_Speed_ad,        &  ! Input
             Transmittance_ad,     &  ! Input, may not be used
             Rel_Azimuth_ad,       &  ! Input, may not be used
             Foam_Fraction_ad,     &  ! Input, may not be used
             iVar)                    ! Optional input
    ! Arguments
    REAL(fp), INTENT(INOUT)  :: Temperature_ad
    REAL(fp), INTENT(INOUT)  :: Salinity_ad
    REAL(fp), INTENT(INOUT)  :: Wind_Speed_ad
    REAL(fp), INTENT(INOUT)  :: Emissivity_ad(4)
    REAL(fp), INTENT(INOUT)  :: Reflectivity_ad(4)
    REAL(fp), INTENT(INOUT)  :: Transmittance_ad
    REAL(fp), INTENT(INOUT)  :: Rel_Azimuth_ad
    REAL(fp), INTENT(INOUT)  :: Foam_Fraction_ad
    TYPE(iVar_type), INTENT(IN) :: iVar

    REAL(fp) :: Emissivity(4), Reflectivity(4)


    IF(ivar%fastem_version == 6) THEN
                          
       CALL rttov_fastem6_ad(ivar%fastem_version,    &  ! Input
                          ivar%Frequency   ,         &  ! Input
                          ivar%Zenith_Angle,         &  ! Input
                          ivar%Temperature ,         &  ! Input
                          ivar%Salinity    ,         &  ! Input
                          ivar%Wind_Speed  ,         &  ! Input
                          Emissivity_ad,             &  ! Input
                          Reflectivity_ad,           &  ! Input
                          Temperature_ad,            &  ! Output
                          Salinity_ad ,              &  ! Output
                          Wind_Speed_ad,             &  ! Output
                          Emissivity  ,              &  ! Output
                          Reflectivity,              &  ! Output
                          ivar%Transmittance,        &  ! Input, may not be used
                          ivar%Rel_Azimuth  ,        &  ! Input, may not be used
                          Transmittance_ad,          &  ! Output
                          Rel_Azimuth_ad ,           &  ! Output
                          ivar%Supply_Foam_Fraction, &  ! Optional input
                          ivar%Foam_Fraction,        &  ! Optional input
                          Foam_Fraction_ad)             ! Optional output
    ELSE 
       CALL rttov_fastem5r1_ad(ivar%fastem_version,    &  ! Input
                          ivar%Frequency   ,         &  ! Input
                          ivar%Zenith_Angle,         &  ! Input
                          ivar%Temperature ,         &  ! Input
                          ivar%Salinity    ,         &  ! Input
                          ivar%Wind_Speed  ,         &  ! Input
                          Emissivity_ad,             &  ! Input
                          Reflectivity_ad,           &  ! Input
                          Temperature_ad,            &  ! Output
                          Salinity_ad ,              &  ! Output
                          Wind_Speed_ad,             &  ! Output
                          Emissivity  ,              &  ! Output
                          Reflectivity,              &  ! Output
                          ivar%Transmittance,        &  ! Input, may not be used
                          ivar%Rel_Azimuth  ,        &  ! Input, may not be used
                          Transmittance_ad,          &  ! Output
                          Rel_Azimuth_ad ,           &  ! Output
                          ivar%Supply_Foam_Fraction, &  ! Optional input
                          ivar%Foam_Fraction,        &  ! Optional input
                          Foam_Fraction_ad)             ! Optional output
    END IF
  END SUBROUTINE Compute_RTTOV_Fastem_AD
  
END MODULE RTTOV_FASTEM_MODULE
