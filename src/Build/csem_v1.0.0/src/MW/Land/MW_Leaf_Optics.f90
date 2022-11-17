!
! MW_Leaf_Optics
!
! Module to compute the Leaf optical properties for LAND surfaces at
! microwave frequencies required for determining the LAND surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to add their leaf
! optical model codes and to simplify integration into
! the CSEM_LandMW_Emiss module.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!

MODULE MW_Leaf_Optics
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp

  IMPLICIT NONE
  
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC  :: CSEM_LeafMW_Permittivity
  PUBLIC  :: CSEM_LeafMW_Optics
  PUBLIC  :: iVar_type

  INTEGER, PARAMETER  :: N_Leaf_Models = 3
  CHARACTER(LEN=20), PARAMETER :: Leaf_Models(N_Leaf_Models) = &
         [CHARACTER(LEN=20) :: 'CSEM', 'CRTM2.1.3', 'CRTM2.0.5']
  
  REAL(fp), PARAMETER :: PI     = 3.141592653589793238462643_fp
  
  INTERFACE CSEM_LeafMW_Permittivity
      MODULE PROCEDURE CRTM_LeafMW_Diel
      MODULE PROCEDURE CRTM_LeafMW_Diel_2     
  END INTERFACE CSEM_LeafMW_Permittivity

  TYPE iVar_type
    PRIVATE
    ! Forward model input values
    INTEGER,  PUBLIC  ::  Model_Option = 1
    REAL(fp), PUBLIC  ::  rhoveg = 0.40_fp
    REAL(fp), PUBLIC  ::  lthick = 0.15_fp
  END TYPE iVar_type


CONTAINS
  
  SUBROUTINE CSEM_LeafMW_Optics(frequency, angle, mge, refl, trans, eveg, iVar)
    REAL(fp),   INTENT(IN)  :: frequency,angle,mge
    REAL(fp),   INTENT(OUT) :: refl(2), trans(2)
    COMPLEX(fp),INTENT(OUT) :: eveg
    TYPE(iVar_type)      :: iVar

    SELECT CASE (TRIM(Leaf_Models(iVar%Model_Option)))

      CASE('CRTM2.0.5')
        CALL CRTM_LeafMW_Diel(Frequency, mge, eveg)
        CALL CRTM_LeafMW_Optics(frequency, angle, eveg, iVar%lthick, &
            refl(1),refl(2),trans(1),trans(2))

      CASE('CRTM2.1.3')
        CALL CRTM_LeafMW_Diel_2(Frequency, mge, eveg, iVar%rhoveg)
        CALL CRTM_LeafMW_Optics(frequency, angle, eveg, iVar%lthick, &
             refl(1),refl(2),trans(1),trans(2))

      CASE Default
        CALL CRTM_LeafMW_Diel_2(Frequency, mge, eveg, iVar%rhoveg)
        CALL Mean_LeafMW_Optics(frequency, eveg, iVar%lthick, &
             refl(1),refl(2),trans(1),trans(2))

    END SELECT
    
  
  END SUBROUTINE CSEM_LeafMW_Optics
  
!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!      CRTM_LeafMW_Diel
!
! PURPOSE:
!       Function to compute the bulk leaf dielectric constant at microwave frequency 
!
!       This function is used as an example demonstrating how to add more optional models.
!
!
! CALLING SEQUENCE:
!       Error_Status = LeafMW_Diel_Options(      &
!                         Frequency,             & ! input
!                         mg,                    & ! leaf water content
!                         esv)                     ! leaf bulk dielectric constant
!
!
! INPUTS:
!
!       frequency:       User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!
!       mg:              User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       esv:             leaf dielectric constant
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!
!:sdoc-:
!----------------------------------------------------------------------------------  
  SUBROUTINE CRTM_LeafMW_Diel(frequency, mg, esv)

    REAL(fp),    INTENT(IN)  :: frequency, mg
    COMPLEX(fp), INTENT(OUT) :: esv
    REAL(fp)    :: en, vf, vb
    COMPLEX(fp) :: C1 = CMPLX(0.0,1.0,fp)

    en = 1.7_fp - (0.74_fp - 6.16_fp*mg)*mg

    vf = mg*(0.55_fp*mg - 0.076_fp)
    vb = 4.64_fp*mg*mg/( 1.0_fp + 7.36_fp*mg*mg)

    esv = en + vf*(4.9_fp + 75.0_fp/(1.0_fp + C1 * frequency/18.0_fp) &
        - C1 * (18.0_fp/frequency)) + &
        vb * ( 2.9_fp + 55.0_fp/(1.0_fp + SQRT(C1*frequency/0.18_fp)) )

    if (AIMAG(esv) >= 0.0) esv = CMPLX(REAL(esv), -0.0001_fp, fp)

  END SUBROUTINE CRTM_LeafMW_Diel


!----------------------------------------------------------------------------------
!
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:   canopy_diel compute the bulk leaf dielectric constant 
!
!   prgmmr:  Fuzhong Weng and Banghua Yan                org: nesdis              date: 2000-11-28
!
! abstract: compute the dielectric constant of the vegetation canopy geomatrical optics approximation
!
!           for vegetation canopy work for horizontal leaves
!
! program history log:
!
! input argument list:
!
!      frequency    -  frequency (ghz)
!      mg           -  gravimetric water content
!
! output argument list:
!
!      esv          -  dielectric constant of leaves
!
! remarks:
!
! references:
!
!     ulaby and el-rayer, 1987: microwave dielectric spectrum of vegetation part ii,
!           dual-dispersion model, ieee trans geosci. remote sensing, 25, 550-557
!
! attributes:
!   language: f90
!   machine:  ibm rs/6000 sp
!
!----------------------------------------------------------------------------------

  SUBROUTINE CRTM_LeafMW_Diel_2(frequency,mg,esv,rhoveg)
    REAL(fp),    INTENT(IN)  :: frequency, mg, rhoveg
    COMPLEX(fp), INTENT(OUT) :: esv
    REAL(fp)    :: en, vf, vb, vmv
    COMPLEX(fp) :: C1 = CMPLX(0.0,1.0,fp)

    vmv = mg*rhoveg/( 1.0_fp - mg*(1.0-rhoveg) )
    en = 1.7_fp + (3.2_fp + 6.5_fp*vmv)*vmv

    vf = vmv*(0.82_fp*vmv + 0.166_fp)
    vb = 31.4_fp*vmv*vmv/( 1.0_fp + 59.5_fp*vmv*vmv)

    esv = en + vf*(4.9_fp + 75.0_fp/(1.0_fp + C1*frequency/18.0_fp) &
        - C1*(18.0_fp/frequency) ) + &
        vb*(2.9_fp + 55.0_fp/(1.0_fp + SQRT(C1*frequency/0.18_fp)))

    IF (AIMAG(esv) >= 0.0) esv = CMPLX(REAL(esv), -0.0001_fp, fp)

  END SUBROUTINE CRTM_LeafMW_Diel_2
  


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CRTM_LeafMW_Optic
!
! PURPOSE:
!       Function to calculate  v-pol and h-pol refelectance and trasmittance of one
!       single leaf at microwave frequency 
!
!       This function is the default leaf optic model currently used by NOAA CRTM
!
! CALLING SEQUENCE:
!       Error_Status = CRTM_LeafMW_Optic(        &
!                         Frequency,             & ! input
!                         theta,                 & ! incident angle in degree
!                         d,                     & ! leaf thicknesss  in mm
!                         esv,                   & ! leaf bulk dielectric constant
!                         rh,th,                 & ! h-pol refelectance and trasmittance
!                         rv,tv)                   ! v-pol refelectance and trasmittance

!
!
! INPUTS:
!
!       frequency:       User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       theta:           User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       d:               User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       rh, th:          H-pol refelectance and trasmittance
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!       rv, tv:          V-pol refelectance and trasmittance
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! CREATION HISTORY:
!       Written by:     Fuzhong Weng and Banghua Yan,2000-11-28
!                       fuzhong.weng@noaa.gov
!                       banghua.yan@noaa.gov
!       Recoded by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!
!:sdoc-:
!----------------------------------------------------------------------------------


  SUBROUTINE CRTM_LeafMW_Optics(frequency, theta, esv, d, rh, rv, th, tv)

    REAL(fp),    INTENT(IN)  :: frequency, theta, d 
    COMPLEX(fp), INTENT(IN)  :: esv
    REAL(fp),    INTENT(OUT) :: rh, rv, th, tv
    
    REAL(fp), PARAMETER :: threshold = 0.999_fp
    REAL(fp)    :: mu
    COMPLEX(fp) :: ix,k0,kz0,kz1,rhc,rvc 
    COMPLEX(fp) :: expval1,factt,factrvc,factrhc

    mu = COS(Pi/180.0_fp*theta)
    ix = CMPLX(0.0, 1.0, fp)

    k0  = CMPLX(2.0_fp*Pi*frequency/300.0_fp, 0.0, fp)   ! 1/mm
    kz0 = k0*mu
    kz1 = k0*SQRT(esv - SIN(Pi/180.0_fp*theta)**2)

    rhc = (kz0 - kz1)/(kz0 + kz1)
    rvc = (esv*kz0 - kz1)/(esv*kz0 + kz1)

    expval1 = EXP(-2.0*ix*kz1*d)
    factrvc = 1.0_fp-rvc**2*expval1
    factrhc = 1.0_fp-rhc**2*expval1
    factt   = 4.0_fp*kz0*kz1*EXP(ix*(kz0-kz1)*d)

    rv = ABS(rvc*(1.0_fp - expval1)/factrvc)**2
    rh = ABS(rhc*(1.0_fp - expval1)/factrhc)**2

    th = ABS(factt/((kz1+kz0)**2*factrhc))**2
    tv = ABS(esv*factt/((kz1+esv*kz0)**2*factrvc))**2
   
  END SUBROUTINE CRTM_LeafMW_Optics


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       LeafMW_Mean_Optic
!
! PURPOSE:
!       Function to calculate averaged refelectance and trasmittance of one
!       single leaf at microwave frequency. 
!       Leaves are taken as individual scatters of a canopy. The averaged refelectance 
!       and trasmittance is used by canopy-level scattering model 
!
!
!
! CALLING SEQUENCE:
!       Error_Status = LeafMW_Mean_Optic(        &
!                         Frequency,             & ! input
!                         leaf_thick,            & ! leaf thicknesss  in mm
!                         esv,                   & ! leaf bulk dielectric constant
!                         rh,th,                 & ! h-pol refelectance and trasmittance
!                         rv,tv)                   ! v-pol refelectance and trasmittance

!
!
! INPUTS:
!
!       frequency:       User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!
!       leaf_thick:      User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
! OUTPUTS:
!       rh, th:          H-pol refelectance and trasmittance
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!       rv, tv:          V-pol refelectance and trasmittance
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!:sdoc-:
!----------------------------------------------------------------------------------

  SUBROUTINE Mean_LeafMW_Optics(frequency, eveg, leaf_thick, rh,rv,th,tv)

    REAL(fp), INTENT(IN)    :: frequency
    REAL(fp), INTENT(IN)    :: leaf_thick
    REAL(fp), INTENT(OUT)   :: rh, rv, th, tv
    COMPLEX(fp), INTENT(IN) :: eveg
    
    REAL(fp) :: rh0, rv0, rh1, rv1, rh2, rv2
    REAL(fp) :: th0, tv0, th1, tv1, th2, tv2
    REAL(fp) :: rr

 
    CALL CRTM_LeafMW_Optics(frequency,0.0_fp,eveg,leaf_thick, rh1,rv1,th1,tv1)
    CALL CRTM_LeafMW_Optics(frequency,75.0_fp,eveg,leaf_thick, rh2,rv2,th2,tv2)
    !rh=(rh2+rh1)/2.0 ; th=(th2+th1)/2.0;rv=(rv2+rv1)/2.0 ; tv=(tv2+tv1)/2.0 

    rr=3.0/4.0
    rr=3.0/4.0
  
    rh0=(rh2*(1.0-rr)+rh1*rr) ; th0=(th2*(1.0-rr)+th1*rr)
    rv0=(rv2*(1.0-rr)+rv1*rr) ; tv0=(tv2*(1.0-rr)+tv1*rr)

    rr=0.75
   
    rh=(rh0*rr+rv0*(1.0-rr)); th=(th0*rr+tv0*(1.0-rr))
    rv=(rh0*(1.0-rr)+rv0*rr); tv=(th0*(1.0-rr)+tv0*rr)
 
  END SUBROUTINE Mean_LeafMW_Optics

 
END MODULE MW_Leaf_Optics



