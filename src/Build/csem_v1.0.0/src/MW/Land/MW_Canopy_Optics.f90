!
! MW_Canopy_Optics
!
! Module to compute the canopy optical properties for LAND surfaces at
! microwave frequencies required for determining the LAND surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to add their canopy
! optical model codes and to simplify integration into
! the CSEM_LandMW_Emiss module.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!



MODULE MW_Canopy_Optics
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use

  USE CSEM_Type_Kinds, ONLY: fp     => CSEM_fp,    &
                             Short  => CSEM_Short, &
                             Double => CSEM_Double
  IMPLICIT NONE
  PRIVATE 
  PUBLIC :: CSEM_CanopyMW_Optics
  PUBLIC :: CSEM_CanopyMW_Optics_TL
  PUBLIC :: CSEM_CanopyMW_Optics_AD
  PUBLIC :: iVar_type
  
  REAL(fp), PARAMETER, PRIVATE :: PI  = 3.141592653589793238462643_fp  
  INTEGER(SHORT),PRIVATE, PARAMETER:: NumOfLeafAngles = 18
  REAL(fp), PARAMETER :: threshold = 0.999_fp

  TYPE OptControl
    LOGICAL geom,lidf
    LOGICAL Coff0,Coff1
  END TYPE OptControl
   
  TYPE CanopyGeometricParams
    REAL(DOUBLE) :: LAI          ! leaf area index 
    REAL(DOUBLE) :: LeafAngle    ! average leaf angle (°)
    REAL(DOUBLE) :: ViewAngle    ! View zenith angle (°)
  END TYPE CanopyGeometricParams
 
  INTERFACE CSEM_CanopyMW_Optics
      MODULE PROCEDURE CRTM_CanopyMW_Optics
  END INTERFACE CSEM_CanopyMW_Optics
  INTERFACE CSEM_CanopyMW_Optics_TL
      MODULE PROCEDURE CRTM_CanopyMW_Optics_TL
  END INTERFACE CSEM_CanopyMW_Optics_TL
  INTERFACE CSEM_CanopyMW_Optics_AD
      MODULE PROCEDURE CRTM_CanopyMW_Optics_AD
  END INTERFACE CSEM_CanopyMW_Optics_AD

  TYPE iVar_type
    PRIVATE
    ! Forward model input values
    REAL(fp) ::  th, tv, rh, rv
    REAL(fp) ::  lai

    INTEGER, PUBLIC  :: Model_Option = 1
   END TYPE iVar_type
 
CONTAINS


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_CanopyMW_Optic
!
! PURPOSE:
!       Subroutine to compute the canopy optical properties at microwave
!       frequencies over a land surface.
!
!       This subroutine is the canopy optical model currently used by NOAA CRTM
!
! CALLING SEQUENCE:
!       CALL CRTM_Canopy_optic(lai, leaf_refl, leaf_trans, g, ssalb, tau, iVar)
!
! INPUTS:
!       leaf_refl :      Leaf  reflectance 
!                        data.
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  2
!                        ATTRIBUTES: INTENT(IN)
!       leaf_trans :     Leaf transmittance
!                        data.
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  2
!                        ATTRIBUTES: INTENT(IN)
!
!       LAI:             Leaf area index
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!
!
! OUTPUTS: 
!       ssalb:         Single scattering albedo
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  2
!                        ATTRIBUTES: INTENT(OUT)
!
!       g:             symetric parameter
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  2
!                        ATTRIBUTES: INTENT(OUT)
!
!       tau :          canopy optical depth
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  2
!                        ATTRIBUTES: INTENT(OUT)
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!
!:sdoc-:
!----------------------------------------------------------------------------------

  SUBROUTINE CRTM_CanopyMW_Optics(lai, leaf_refl, leaf_trans, g, ssalb, tau, iVar)

    REAL(fp) :: lai
    REAL(fp) :: leaf_refl(2), leaf_trans(2)
    REAL(fp) :: ssalb(2), tau(2), g(2)
 
    TYPE(iVar_type) :: iVar
    
    iVar%rh  = leaf_refl(1) 
    iVar%rv  = leaf_refl(2)
    iVar%th  = leaf_trans(1) 
    iVar%tv  = leaf_trans(2)
    iVar%lai = lai
    
    g = 0.5_fp
    tau = 0.5_fp*lai*(2.0_fp-iVar%th-iVar%tv)
    ssalb = MIN((iVar%rv+iVar%rh)/(2.0_fp-iVar%tv-iVar%th),threshold)

  END SUBROUTINE CRTM_CanopyMW_Optics

  SUBROUTINE CRTM_CanopyMW_Optics_TL(LAI_TL,ssalb_TL,tau_TL, iVar)

    REAL(fp) ::  lai_TL
    REAL(fp) ::  ssalb_TL(2), tau_TL(2)
    TYPE(iVar_type) :: iVar

    tau_TL = 0.5_fp*(2.0_fp-iVar%tv-iVar%th)*lai_TL

    ssalb_TL = 0.0_fp

  END SUBROUTINE CRTM_CanopyMW_Optics_TL
  
  SUBROUTINE CRTM_CanopyMW_Optics_AD(LAI_AD,ssalb_AD, tau_AD, iVar)

    REAL(fp) ::  lai_AD
    REAL(fp) ::  ssalb_AD(2),tau_AD(2)

    TYPE(iVar_type) :: iVar

    LAI_AD =  LAI_AD +  0.5_fp*(2.0_fp-iVar%tv-iVar%th)*tau_AD(2)
    LAI_AD =  LAI_AD +  0.5_fp*(2.0_fp-iVar%tv-iVar%th)*tau_AD(1)
    tau_AD = 0.0_fp   
    ssalb_AD = 0.0_fp 

  END SUBROUTINE CRTM_CanopyMW_Optics_AD


END MODULE MW_Canopy_Optics



