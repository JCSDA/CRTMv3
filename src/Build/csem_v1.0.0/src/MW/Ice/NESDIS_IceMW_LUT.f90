!
! NESDIS_IceEM_Parameters
!
! Module containing the parameters related to microwave Ice emissivity model
!
!
! CREATION HISTORY:
!       Written by:    Ming Chen 
!                      ming.chen@noaa.gov
!                      Fuzhong Weng
!                      fuzhong.weng@noaa.gov
!

MODULE NESDIS_MW_IceEM_LUT


  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Interpolation, ONLY: NPTS        , &
                                CSEM_LPoly_type  , &
                                CSEM_find_index  , &
                                CSEM_interp_1D   , &
                                CSEM_Clear_LPoly , &
                                CSEM_LPoly       
  ! Disable implicit typing
  IMPLICIT NONE
  
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE


  ! -----------------
  ! Module parameters
  ! -----------------
  ! Version Id for the module
  CHARACTER(*), PARAMETER :: MODULE_VERSION_ID = &
  '$Id: NESDIS_IceEM_Parameters.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'

  INTEGER, PUBLIC, PARAMETER :: N_MWICE_TYPES = 13
   
  CHARACTER(LEN=20), PUBLIC, PARAMETER, DIMENSION(N_MWICE_TYPES):: &
    NESDIS_ICE_TYPE_LIST = (/ &
   'RS_ICE_A            ', &         !1
   'RS_ICE_B            ', &         !2
   'MIXED_NEWICE_SNOW   ', &         !3
   'NARE_NEWICE_        ', &         !4
   'BROKEN_ICE_         ', &         !5
   'FIRST_YEAR_ICE_     ', &         !6
   'COMPOSITE_PACK_ICE  ', &         !7
   'RS_ICE_C            ', &         !8
   'FAST_ICE            ', &         !9
   'RS_ICE_D            ', &         !10
   'RS_ICE_E            ', &         !11
   'RS_ICE_F            ', &         !12
   'GREASE_ICE          ' /)         !13

  INTEGER, PUBLIC, PARAMETER :: N_MWICE_FREQUENCY  = 7

  REAL(fp), PARAMETER, DIMENSION(N_MWICE_FREQUENCY) :: FREQUENCY_AMSRE=  &
    (/ 6.925_fp, 10.65_fp, 18.7_fp,23.8_fp, 36.5_fp, 89.0_fp, 157._fp/)
 
 ! Define sixteen MW H-POL emissivity spectra for AMSRE ALGORITHMS
  REAL(fp),PARAMETER, DIMENSION(N_MWICE_FREQUENCY,N_MWICE_TYPES) ::  &
  NESDIS_AMSRE_IceEM_LUT_H = RESHAPE( (/ &
    0.88_fp, 0.92_fp, 0.94_fp, 0.94_fp, 0.95_fp, 0.92_fp, 0.91_fp,& !RS_ICE_A_EMISS
    0.81_fp, 0.82_fp, 0.85_fp, 0.86_fp, 0.87_fp, 0.88_fp, 0.87_fp,& !RS_ICE_B_EMISS
    0.83_fp, 0.84_fp, 0.86_fp, 0.85_fp, 0.84_fp, 0.82_fp, 0.80_fp,& !MIXED_NEWICE_SNOW_EMISS
    0.74_fp, 0.75_fp, 0.76_fp, 0.76_fp, 0.77_fp, 0.73_fp, 0.73_fp,& !NARE_NEWICE_EMISS
    0.71_fp, 0.73_fp, 0.76_fp, 0.77_fp, 0.80_fp, 0.72_fp, 0.69_fp,& !BROKEN_ICE_EMISS
    0.91_fp, 0.90_fp, 0.89_fp, 0.88_fp, 0.86_fp, 0.76_fp, 0.67_fp,& !FIRST_YEAR_ICE_EMISS
    0.85_fp, 0.84_fp, 0.83_fp, 0.82_fp, 0.79_fp, 0.67_fp, 0.57_fp,& !COMPOSITE_PACK_ICE_EMISS
    0.90_fp, 0.87_fp, 0.81_fp, 0.78_fp, 0.69_fp, 0.60_fp, 0.56_fp,& !RS_ICE_C_EMISS
    0.80_fp, 0.80_fp, 0.78_fp, 0.76_fp, 0.72_fp, 0.60_fp, 0.53_fp,& !FAST_ICE_EMISS 
    0.71_fp, 0.71_fp, 0.70_fp, 0.70_fp, 0.70_fp, 0.59_fp, 0.54_fp,& !RS_ICE_D_EMISS
    0.55_fp, 0.59_fp, 0.60_fp, 0.61_fp, 0.62_fp, 0.67_fp, 0.69_fp,& !RS_ICE_E_EMISS
    0.48_fp, 0.51_fp, 0.56_fp, 0.57_fp, 0.60_fp, 0.64_fp, 0.65_fp,& !RS_ICE_F_EMISS
    0.42_fp, 0.42_fp, 0.43_fp, 0.45_fp, 0.49_fp, 0.54_fp, 0.56_fp & !GREASE_ICE_EMISS
  /),(/N_MWICE_FREQUENCY,N_MWICE_TYPES/))


 ! Define sixteen MW V-POL emissivity spectra for AMSRE ALGORITHMS

  REAL(fp),PARAMETER, DIMENSION(N_MWICE_FREQUENCY,N_MWICE_TYPES) ::  &
  NESDIS_AMSRE_IceEM_LUT_V = RESHAPE( (/ &
    0.96_fp, 0.97_fp, 0.99_fp, 0.99_fp, 0.99_fp, 0.98_fp, 0.97_fp,& !RS_ICE_A_EMISS
    0.95_fp, 0.96_fp, 0.99_fp, 0.98_fp, 0.97_fp, 0.94_fp, 0.93_fp,& !RS_ICE_B_EMISS
    0.96_fp, 0.96_fp, 0.95_fp, 0.94_fp, 0.93_fp, 0.88_fp, 0.86_fp,& !MIXED_NEWICE_SNOW_EMISS
    0.88_fp, 0.89_fp, 0.91_fp, 0.91_fp, 0.91_fp, 0.88_fp, 0.88_fp,& !NARE_NEWICE_EMISS
    0.85_fp, 0.87_fp, 0.91_fp, 0.91_fp, 0.91_fp, 0.87_fp, 0.84_fp,& !BROKEN_ICE_EMISS
    0.98_fp, 0.98_fp, 0.98_fp, 0.97_fp, 0.95_fp, 0.84_fp, 0.75_fp,& !FIRST_YEAR_ICE_EMISS
    0.98_fp, 0.97_fp, 0.95_fp, 0.93_fp, 0.89_fp, 0.72_fp, 0.62_fp,& !COMPOSITE_PACK_ICE_EMISS
    0.99_fp, 0.96_fp, 0.90_fp, 0.86_fp, 0.75_fp, 0.66_fp, 0.62_fp,& !RS_ICE_C_EMISS
    0.95_fp, 0.95_fp, 0.94_fp, 0.91_fp, 0.85_fp, 0.69_fp, 0.62_fp,& !FAST_ICE_EMISS 
    0.87_fp, 0.87_fp, 0.88_fp, 0.88_fp, 0.88_fp, 0.77_fp, 0.72_fp,& !RS_ICE_D_EMISS
    0.77_fp, 0.78_fp, 0.81_fp, 0.82_fp, 0.84_fp, 0.86_fp, 0.88_fp,& !RS_ICE_E_EMISS
    0.71_fp, 0.73_fp, 0.77_fp, 0.78_fp, 0.81_fp, 0.86_fp, 0.87_fp,& !RS_ICE_F_EMISS
    0.66_fp, 0.67_fp, 0.70_fp, 0.72_fp, 0.76_fp, 0.82_fp, 0.84_fp & !GREASE_ICE_EMISS
  /),(/N_MWICE_FREQUENCY,N_MWICE_TYPES/))
  
 ! Define sixteen MW weighted emissivity spectra for AMSRE ALGORITHMS

  REAL(fp),PARAMETER, DIMENSION(N_MWICE_FREQUENCY,N_MWICE_TYPES) :: &
  NESDIS_AMSRE_IceEM_LUT_M  =   RESHAPE( (/ &
    0.93_fp, 0.94_fp, 0.96_fp, 0.97_fp, 0.97_fp, 0.94_fp, 0.93_fp,& !RS_ICE_A_EMISS
    0.86_fp, 0.87_fp, 0.90_fp, 0.91_fp, 0.90_fp, 0.90_fp, 0.89_fp,& !RS_ICE_B_EMISS
    0.88_fp, 0.88_fp, 0.89_fp, 0.88_fp, 0.87_fp, 0.84_fp, 0.82_fp,& !MIXED_NEWICE_SNOW_EMISS
    0.80_fp, 0.81_fp, 0.81_fp, 0.81_fp, 0.80_fp, 0.79_fp, 0.79_fp,& !NARE_NEWICE_EMISS
    0.75_fp, 0.78_fp, 0.80_fp, 0.81_fp, 0.80_fp, 0.77_fp, 0.74_fp,& !BROKEN_ICE_EMISS
    0.93_fp, 0.93_fp, 0.92_fp, 0.92_fp, 0.89_fp, 0.78_fp, 0.69_fp,& !FIRST_YEAR_ICE_EMISS
    0.89_fp, 0.88_fp, 0.87_fp, 0.85_fp, 0.82_fp, 0.69_fp, 0.59_fp,& !COMPOSITE_PACK_ICE_EMISS
    0.92_fp, 0.90_fp, 0.83_fp, 0.78_fp, 0.73_fp, 0.62_fp, 0.58_fp,& !RS_ICE_C_EMISS
    0.85_fp, 0.85_fp, 0.84_fp, 0.81_fp, 0.78_fp, 0.63_fp, 0.56_fp,& !FAST_ICE_EMISS 
    0.76_fp, 0.76_fp, 0.76_fp, 0.76_fp, 0.74_fp, 0.65_fp, 0.60_fp,& !RS_ICE_D_EMISS
    0.63_fp, 0.65_fp, 0.67_fp, 0.68_fp, 0.70_fp, 0.74_fp, 0.75_fp,& !RS_ICE_E_EMISS
    0.54_fp, 0.60_fp, 0.64_fp, 0.67_fp, 0.70_fp, 0.71_fp, 0.72_fp,& !RS_ICE_F_EMISS
    0.49_fp, 0.51_fp, 0.53_fp, 0.55_fp, 0.58_fp, 0.65_fp, 0.67_fp & !GREASE_ICE_EMISS
  /),(/N_MWICE_FREQUENCY,N_MWICE_TYPES/))

  ! --------------------------------------
  ! Structure definitions to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! The interpolating polynomials
    TYPE(CSEM_LPoly_type) :: xlp
    ! The LUT interpolation indices
    INTEGER :: i1, i2
    ! The LUT interpolation boundary check
    LOGICAL :: x_outbound
    ! The interpolation input
    REAL(fp) :: x_int
    ! The data to be interpolated
    REAL(fp) :: x(NPTS)
  END TYPE iVar_type
  

  PUBLIC :: Load_MW_IceEM_LUT
  
  INTERFACE Load_MW_IceEM_LUT
    MODULE PROCEDURE Set_AMSRE_IceEM_LUT
  END INTERFACE Load_MW_IceEM_LUT
 
CONTAINS

  
  FUNCTION Set_AMSRE_IceEM_LUT( &
      frequency,           & ! INPUT, wavelength in micrometer
      emissivity_m,        & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_h,        & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_v)        & ! OUTPUT, surface emissivity (0 - 1)
    RESULT (Error_Status)


    REAL(fp), INTENT( IN )  :: frequency(:)
    REAL(fp), INTENT( OUT ) :: emissivity_m(:,:)
    REAL(fp), INTENT( OUT ) :: emissivity_h(:,:)
    REAL(fp), INTENT( OUT ) :: emissivity_v(:,:)
    
    INTEGER :: Ice_type    
    INTEGER :: i,Error_Status 

    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !
    Error_Status = 0
    DO i = 1, size(frequency)
      DO Ice_type = 1, N_MWICE_TYPES
        Error_Status = Interp_AMSRE_IceEM_LUT( &
          frequency(i),                 & ! INPUT, wavelength in micrometer
          emissivity_m(Ice_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          emissivity_h(Ice_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          emissivity_v(Ice_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          Ice_type )       
      ENDDO
    ENDDO
  END FUNCTION  Set_AMSRE_IceEM_LUT


  FUNCTION Interp_AMSRE_IceEM_LUT( &
      frequency,         & ! INPUT, wavelength in micrometer
      emissivity_m,      & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_h,      & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_v,      & ! OUTPUT, surface emissivity (0 - 1)
      Ice_type )         & ! INPUT, surface type (1 - 16)
    RESULT (Error_Status)


    INTEGER , INTENT( IN )  :: Ice_type
    REAL(fp), INTENT( IN )  :: frequency
    REAL(fp), INTENT( OUT ) :: emissivity_m, emissivity_h, emissivity_v
    INTEGER :: Error_Status 
    TYPE(iVar_type) :: iVar
 
    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !
    iVar%x_int = MAX(MIN(FREQUENCY_AMSRE(N_MWICE_FREQUENCY),&
                 frequency), FREQUENCY_AMSRE(1))
    CALL CSEM_find_index(FREQUENCY_AMSRE, &
                 iVar%x_int, iVar%i1, iVar%i2, iVar%x_outbound)
    iVar%x = FREQUENCY_AMSRE(iVar%i1:iVar%i2)

    ! Calculate the interpolating polynomial
    CALL CSEM_LPoly( iVar%x, iVar%x_int, & ! Input
                iVar%xlp            )      ! Output

    ! Perform Interpolation
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_IceEM_LUT_M(iVar%i1:iVar%i2, Ice_Type), &
          iVar%xlp, Emissivity_m)
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_IceEM_LUT_H(iVar%i1:iVar%i2, Ice_Type), &
          iVar%xlp, Emissivity_h)
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_IceEM_LUT_V(iVar%i1:iVar%i2, Ice_Type), &
          iVar%xlp, Emissivity_v)
    Error_Status = 0
  END FUNCTION Interp_AMSRE_IceEM_LUT

END MODULE NESDIS_MW_IceEM_LUT
