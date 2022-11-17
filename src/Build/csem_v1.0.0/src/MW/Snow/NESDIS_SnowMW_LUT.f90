!
! NESDIS_MW_SnowEM_LUT
!
! Module containing the parameters related to microwave snow emissivity model
!
!
! CREATION HISTORY:
!       Written by:     Banghua Yan, 03-Jun-2005 
!                       banghua.yan@noaa.gov
!                       Fuzhong Weng
!                       fuzhong.weng@noaa.gov
!

MODULE NESDIS_MW_SnowEM_LUT


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
  '$Id: NESDIS_SnowEM_Parameters.f90 21141 2012-09-14 17:40:43Z paul.vandelst@noaa.gov $'

  INTEGER, PUBLIC, PARAMETER :: N_MWSNOW_TYPES = 16
 ! Snow types
  INTEGER, PUBLIC, PARAMETER :: INVALID_SNOW_TYPE   = -999
  INTEGER, PUBLIC, PARAMETER :: WET_SNOW            =  1
  INTEGER, PUBLIC, PARAMETER :: GRASS_AFTER_SNOW    =  2
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_A           =  3
  INTEGER, PUBLIC, PARAMETER :: POWDER_SNOW         =  4
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_B           =  5
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_C           =  6
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_D           =  7
  INTEGER, PUBLIC, PARAMETER :: THIN_CRUST_SNOW     =  8
  INTEGER, PUBLIC, PARAMETER :: RS_SNOW_E           =  9
  INTEGER, PUBLIC, PARAMETER :: BOTTOM_CRUST_SNOW_A = 10
  INTEGER, PUBLIC, PARAMETER :: SHALLOW_SNOW        = 11
  INTEGER, PUBLIC, PARAMETER :: DEEP_SNOW           = 12
  INTEGER, PUBLIC, PARAMETER :: CRUST_SNOW          = 13
  INTEGER, PUBLIC, PARAMETER :: MEDIUM_SNOW         = 14
  INTEGER, PUBLIC, PARAMETER :: BOTTOM_CRUST_SNOW_B = 15
  INTEGER, PUBLIC, PARAMETER :: THICK_CRUST_SNOW    = 16
  
  CHARACTER(LEN=20), PUBLIC, PARAMETER, DIMENSION(N_MWSNOW_TYPES):: &
    NESDIS_SNOW_TYPE_LIST = (/ &
   'WET_SNOW            ', &
   'GRASS_AFTER_SNOW    ', &
   'RS_SNOW_A           ', &
   'POWDER_SNOW         ', &
   'RS_SNOW_B           ', &
   'RS_SNOW_C           ', &
   'RS_SNOW_D           ', &
   'THIN_CRUST_SNOW     ', &
   'RS_SNOW_E           ', &
   'BOTTOM_CRUST_SNOW_A ', &
   'SHALLOW_SNOW        ', &
   'DEEP_SNOW           ', &
   'CRUST_SNOW          ', &
   'MEDIUM_SNOW         ', &
   'BOTTOM_CRUST_SNOW_B ', &
   'THICK_CRUST_SNOW    ' /)

  INTEGER, PUBLIC, PARAMETER :: N_MWSNOW_FREQUENCY  = 10
  INTEGER, PUBLIC, PARAMETER :: N_AMSRE_SNOW_FREQ   = 7

  REAL(fp),PARAMETER, DIMENSION(N_AMSRE_SNOW_FREQ):: FREQUENCY_AMSRE = &
    (/ 6.925_fp, 10.65_fp, 18.7_fp,23.8_fp, 36.5_fp, 89.0_fp, 150._fp/)


 ! Define sixteen MW H-POL emissivity spectra for AMSRE ALGORITHMS
  REAL(fp),PARAMETER, DIMENSION(N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES) ::  &
  NESDIS_AMSRE_SnowEM_LUT_H = RESHAPE( (/ &
    0.93_fp, 0.92_fp, 0.93_fp, 0.94_fp, 0.93_fp, 0.93_fp, 0.90_fp, &
    0.91_fp, 0.90_fp, 0.90_fp, 0.90_fp, 0.91_fp, 0.90_fp, 0.85_fp, &
    0.85_fp, 0.85_fp, 0.84_fp, 0.84_fp, 0.82_fp, 0.80_fp, 0.80_fp, &
    0.90_fp, 0.90_fp, 0.92_fp, 0.92_fp, 0.90_fp, 0.80_fp, 0.79_fp, &
    0.82_fp, 0.81_fp, 0.77_fp, 0.76_fp, 0.74_fp, 0.74_fp, 0.74_fp, &
    0.84_fp, 0.83_fp, 0.80_fp, 0.78_fp, 0.77_fp, 0.75_fp, 0.69_fp, &
    0.77_fp, 0.77_fp, 0.76_fp, 0.75_fp, 0.73_fp, 0.71_fp, 0.71_fp, &
    0.95_fp, 0.94_fp, 0.95_fp, 0.94_fp, 0.89_fp, 0.75_fp, 0.65_fp, &
    0.73_fp, 0.73_fp, 0.74_fp, 0.72_fp, 0.71_fp, 0.68_fp, 0.67_fp, &
    0.88_fp, 0.87_fp, 0.86_fp, 0.85_fp, 0.80_fp, 0.68_fp, 0.63_fp, &
    0.86_fp, 0.84_fp, 0.80_fp, 0.78_fp, 0.72_fp, 0.62_fp, 0.57_fp, &
    0.87_fp, 0.85_fp, 0.83_fp, 0.80_fp, 0.77_fp, 0.68_fp, 0.62_fp, &
    0.82_fp, 0.78_fp, 0.74_fp, 0.71_fp, 0.67_fp, 0.64_fp, 0.64_fp, &
    0.90_fp, 0.90_fp, 0.89_fp, 0.88_fp, 0.83_fp, 0.53_fp, 0.48_fp, &
    0.87_fp, 0.85_fp, 0.84_fp, 0.82_fp, 0.74_fp, 0.53_fp, 0.49_fp, &
    0.85_fp, 0.84_fp, 0.83_fp, 0.81_fp, 0.79_fp, 0.51_fp, 0.46_fp  &
  /),(/N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES/))


 ! Define sixteen MW V-POL emissivity spectra for AMSRE ALGORITHMS

  REAL(fp),PARAMETER, DIMENSION(N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES) ::  &
  NESDIS_AMSRE_SnowEM_LUT_V = RESHAPE( (/ &
    0.96_fp, 0.94_fp, 0.96_fp, 0.95_fp, 0.94_fp, 0.94_fp, 0.91_fp, &
    0.96_fp, 0.94_fp, 0.95_fp, 0.96_fp, 0.96_fp, 0.92_fp, 0.87_fp, &
    0.99_fp, 0.97_fp, 0.96_fp, 0.96_fp, 0.93_fp, 0.87_fp, 0.87_fp, &
    0.98_fp, 0.97_fp, 0.99_fp, 0.98_fp, 0.96_fp, 0.84_fp, 0.83_fp, &
    0.97_fp, 0.95_fp, 0.93_fp, 0.92_fp, 0.89_fp, 0.84_fp, 0.84_fp, &
    1.00_fp, 0.97_fp, 0.96_fp, 0.94_fp, 0.91_fp, 0.84_fp, 0.78_fp, &
    0.99_fp, 0.96_fp, 0.93_fp, 0.90_fp, 0.86_fp, 0.80_fp, 0.80_fp, &
    0.98_fp, 0.97_fp, 0.98_fp, 0.97_fp, 0.92_fp, 0.77_fp, 0.67_fp, &
    0.98_fp, 0.95_fp, 0.90_fp, 0.86_fp, 0.82_fp, 0.74_fp, 0.73_fp, &
    0.96_fp, 0.95_fp, 0.95_fp, 0.93_fp, 0.87_fp, 0.71_fp, 0.66_fp, &
    0.97_fp, 0.95_fp, 0.94_fp, 0.90_fp, 0.84_fp, 0.68_fp, 0.63_fp, &
    0.96_fp, 0.94_fp, 0.92_fp, 0.90_fp, 0.85_fp, 0.77_fp, 0.71_fp, &
    0.98_fp, 0.96_fp, 0.93_fp, 0.90_fp, 0.81_fp, 0.71_fp, 0.71_fp, &
    0.99_fp, 0.97_fp, 0.98_fp, 0.96_fp, 0.92_fp, 0.57_fp, 0.52_fp, &
    1.00_fp, 0.97_fp, 0.97_fp, 0.95_fp, 0.86_fp, 0.58_fp, 0.54_fp, &
    0.98_fp, 0.96_fp, 0.96_fp, 0.94_fp, 0.89_fp, 0.56_fp, 0.51_fp  &
  /),(/N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES/))
  
 ! Define sixteen MW weighted emissivity spectra for AMSRE ALGORITHMS

  REAL(fp),PARAMETER, DIMENSION(N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES) :: &
  NESDIS_AMSRE_SnowEM_LUT_M  =   RESHAPE( (/ &
    0.91_fp, 0.93_fp, 0.94_fp, 0.95_fp, 0.95_fp, 0.93_fp, 0.93_fp, &
    0.91_fp, 0.92_fp, 0.91_fp, 0.90_fp, 0.91_fp, 0.91_fp, 0.91_fp, &
    0.90_fp, 0.89_fp, 0.88_fp, 0.87_fp, 0.86_fp, 0.82_fp, 0.82_fp, &
    0.92_fp, 0.93_fp, 0.94_fp, 0.94_fp, 0.92_fp, 0.80_fp, 0.80_fp, &
    0.87_fp, 0.86_fp, 0.83_fp, 0.80_fp, 0.79_fp, 0.77_fp, 0.77_fp, &
    0.89_fp, 0.88_fp, 0.85_fp, 0.84_fp, 0.83_fp, 0.79_fp, 0.79_fp, &
    0.84_fp, 0.83_fp, 0.82_fp, 0.80_fp, 0.78_fp, 0.72_fp, 0.72_fp, &
    0.95_fp, 0.96_fp, 0.96_fp, 0.95_fp, 0.91_fp, 0.75_fp, 0.75_fp, &
    0.80_fp, 0.80_fp, 0.80_fp, 0.79_fp, 0.75_fp, 0.70_fp, 0.70_fp, &
    0.91_fp, 0.90_fp, 0.89_fp, 0.87_fp, 0.82_fp, 0.69_fp, 0.69_fp, &
    0.90_fp, 0.89_fp, 0.85_fp, 0.82_fp, 0.76_fp, 0.65_fp, 0.65_fp, &
    0.89_fp, 0.88_fp, 0.86_fp, 0.83_fp, 0.78_fp, 0.70_fp, 0.70_fp, &
    0.88_fp, 0.86_fp, 0.80_fp, 0.75_fp, 0.69_fp, 0.67_fp, 0.67_fp, &
    0.96_fp, 0.97_fp, 0.92_fp, 0.87_fp, 0.72_fp, 0.50_fp, 0.50_fp, &
    0.93_fp, 0.94_fp, 0.89_fp, 0.85_fp, 0.74_fp, 0.48_fp, 0.48_fp, &
    0.88_fp, 0.88_fp, 0.87_fp, 0.85_fp, 0.77_fp, 0.52_fp, 0.52_fp  &
  /),(/N_AMSRE_SNOW_FREQ,N_MWSNOW_TYPES/))

  
  REAL(fp), PARAMETER, DIMENSION(N_MWSNOW_FREQUENCY) :: NESDIS_MWSNOW_FREQUENCY= (/  &
    4.90_fp,6.93_fp,10.65_fp,18.7_fp,23.8_fp,31.4_fp,50.3_fp,52.5_fp,89._fp,150._fp /)
    
  REAL(fp),PARAMETER :: NESDIS_MWSNOW_Emiss_LUT(N_MWSNOW_FREQUENCY,N_MWSNOW_TYPES) = &
    RESHAPE( (/ &
    0.87_fp,0.89_fp,0.91_fp,0.93_fp,0.94_fp,0.94_fp,0.94_fp,0.93_fp,0.92_fp,0.90_fp, &
    0.91_fp,0.91_fp,0.92_fp,0.91_fp,0.90_fp,0.90_fp,0.91_fp,0.91_fp,0.91_fp,0.86_fp, &
    0.90_fp,0.89_fp,0.88_fp,0.87_fp,0.86_fp,0.86_fp,0.85_fp,0.85_fp,0.82_fp,0.82_fp, &
    0.91_fp,0.91_fp,0.93_fp,0.93_fp,0.93_fp,0.93_fp,0.89_fp,0.88_fp,0.79_fp,0.79_fp, &
    0.90_fp,0.89_fp,0.88_fp,0.85_fp,0.84_fp,0.83_fp,0.83_fp,0.82_fp,0.79_fp,0.73_fp, &
    0.90_fp,0.89_fp,0.86_fp,0.82_fp,0.80_fp,0.79_fp,0.78_fp,0.78_fp,0.77_fp,0.77_fp, &
    0.88_fp,0.86_fp,0.85_fp,0.80_fp,0.78_fp,0.77_fp,0.77_fp,0.76_fp,0.72_fp,0.72_fp, &
    0.93_fp,0.94_fp,0.96_fp,0.96_fp,0.95_fp,0.93_fp,0.87_fp,0.86_fp,0.74_fp,0.65_fp, &
    0.87_fp,0.86_fp,0.84_fp,0.80_fp,0.76_fp,0.76_fp,0.75_fp,0.75_fp,0.70_fp,0.69_fp, &
    0.87_fp,0.86_fp,0.83_fp,0.77_fp,0.73_fp,0.68_fp,0.66_fp,0.66_fp,0.68_fp,0.67_fp, &
    0.89_fp,0.89_fp,0.88_fp,0.87_fp,0.86_fp,0.82_fp,0.77_fp,0.76_fp,0.69_fp,0.64_fp, &
    0.88_fp,0.87_fp,0.86_fp,0.83_fp,0.81_fp,0.77_fp,0.74_fp,0.73_fp,0.69_fp,0.64_fp, &
    0.86_fp,0.86_fp,0.86_fp,0.85_fp,0.82_fp,0.78_fp,0.69_fp,0.68_fp,0.51_fp,0.47_fp, &
    0.89_fp,0.88_fp,0.87_fp,0.83_fp,0.80_fp,0.75_fp,0.70_fp,0.70_fp,0.64_fp,0.60_fp, &
    0.91_fp,0.92_fp,0.93_fp,0.88_fp,0.84_fp,0.76_fp,0.66_fp,0.64_fp,0.48_fp,0.44_fp, &
    0.94_fp,0.95_fp,0.97_fp,0.91_fp,0.86_fp,0.74_fp,0.63_fp,0.63_fp,0.50_fp,0.45_fp  &
  /), (/N_MWSNOW_FREQUENCY, N_MWSNOW_TYPES/))
  
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
  

  PUBLIC :: Load_MW_SnowEM_LUT
  
  INTERFACE Load_MW_SnowEM_LUT
    MODULE PROCEDURE Set_MW_SnowEM_LUT
    MODULE PROCEDURE Set_AMSRE_SnowEM_LUT
  END INTERFACE Load_MW_SnowEM_LUT
 
CONTAINS

  FUNCTION Set_MW_SnowEM_LUT( &
      frequency,         & ! INPUT, wavelength in micrometer
      emissivity)        & ! OUTPUT, surface emissivity (0 - 1)
      RESULT (Error_Status)


    REAL(fp), INTENT( IN )  :: frequency(:)
    REAL(fp), INTENT( OUT ) :: emissivity(:,:)
    
    INTEGER :: snow_type    
    INTEGER :: i,Error_Status 

    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !
    Error_Status = 0
    DO i = 1, size(frequency)
      DO snow_type = 1, N_MWSNOW_TYPES
        Error_Status = Interp_MW_SnowEM_LUT( &
          frequency(i),               & ! INPUT, wavelength in micrometer
          emissivity(snow_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          snow_type )       
      ENDDO
    ENDDO
  END FUNCTION  Set_MW_SnowEM_LUT
  
  FUNCTION Set_AMSRE_SnowEM_LUT( &
      frequency,           & ! INPUT, wavelength in micrometer
      emissivity_m,        & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_h,        & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_v)        & ! OUTPUT, surface emissivity (0 - 1)
    RESULT (Error_Status)


    REAL(fp), INTENT( IN )  :: frequency(:)
    REAL(fp), INTENT( OUT ) :: emissivity_m(:,:)
    REAL(fp), INTENT( OUT ) :: emissivity_h(:,:)
    REAL(fp), INTENT( OUT ) :: emissivity_v(:,:)
    
    INTEGER :: snow_type    
    INTEGER :: i,Error_Status 

    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !
    Error_Status = 0
    
    DO i = 1, size(frequency)
      DO snow_type = 1, N_MWSNOW_TYPES
        Error_Status = Interp_AMSRE_SnowEM_LUT( &
          frequency(i),                 & ! INPUT, wavelength in micrometer
          emissivity_m(snow_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          emissivity_h(snow_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          emissivity_v(snow_type,i),    & ! OUTPUT, surface emissivity (0 - 1)
          snow_type )       
      ENDDO
    ENDDO
  END FUNCTION  Set_AMSRE_SnowEM_LUT

  FUNCTION Interp_MW_SnowEM_LUT( &
      frequency,         & ! INPUT, wavelength in micrometer
      emissivity,        & ! OUTPUT, surface emissivity (0 - 1)
      snow_type )        & ! INPUT, surface type (1 - 16)
    RESULT (Error_Status)


    INTEGER , INTENT( IN )  :: snow_type
    REAL(fp), INTENT( IN )  :: frequency
    REAL(fp), INTENT( OUT ) :: emissivity
    INTEGER :: Error_Status 
    TYPE(iVar_type) :: iVar
 
    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !

    iVar%x_int = MAX(MIN(NESDIS_MWSNOW_FREQUENCY(N_MWSNOW_FREQUENCY),&
                 frequency), NESDIS_MWSNOW_FREQUENCY(1))
    CALL CSEM_find_index(NESDIS_MWSNOW_FREQUENCY, &
                 iVar%x_int, iVar%i1, iVar%i2, iVar%x_outbound)
    iVar%x = NESDIS_MWSNOW_FREQUENCY(iVar%i1:iVar%i2)

    ! Calculate the interpolating polynomial
    CALL CSEM_LPoly( iVar%x, iVar%x_int, & ! Input
                iVar%xlp            )      ! Output

    ! Perform Interpolation
    CALL CSEM_interp_1D( &
          NESDIS_MWSNOW_Emiss_LUT(iVar%i1:iVar%i2, Snow_Type), &
          iVar%xlp, Emissivity)
    Error_Status = 0
  END FUNCTION Interp_MW_SnowEM_LUT


  FUNCTION Interp_AMSRE_SnowEM_LUT( &
      frequency,         & ! INPUT, wavelength in micrometer
      emissivity_m,      & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_h,      & ! OUTPUT, surface emissivity (0 - 1)
      emissivity_v,      & ! OUTPUT, surface emissivity (0 - 1)
      snow_type )        & ! INPUT, surface type (1 - 16)
    RESULT (Error_Status)


    INTEGER , INTENT( IN )  :: snow_type
    REAL(fp), INTENT( IN )  :: frequency
    REAL(fp), INTENT( OUT ) :: emissivity_m, emissivity_h, emissivity_v
    INTEGER :: Error_Status 
    TYPE(iVar_type) :: iVar
 
    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !

    iVar%x_int = MAX(MIN(FREQUENCY_AMSRE(N_AMSRE_SNOW_FREQ),&
                 frequency), FREQUENCY_AMSRE(1))
    CALL CSEM_find_index(FREQUENCY_AMSRE, &
                 iVar%x_int, iVar%i1, iVar%i2, iVar%x_outbound)
    iVar%x = FREQUENCY_AMSRE(iVar%i1:iVar%i2)

    ! Calculate the interpolating polynomial
    CALL CSEM_LPoly( iVar%x, iVar%x_int, & ! Input
                iVar%xlp            )      ! Output

    ! Perform Interpolation
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_SnowEM_LUT_M(iVar%i1:iVar%i2, Snow_Type), &
          iVar%xlp, Emissivity_m)
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_SnowEM_LUT_H(iVar%i1:iVar%i2, Snow_Type), &
          iVar%xlp, Emissivity_h)
    CALL CSEM_interp_1D( &
          NESDIS_AMSRE_SnowEM_LUT_V(iVar%i1:iVar%i2, Snow_Type), &
          iVar%xlp, Emissivity_v)
    Error_Status = 0
  END FUNCTION Interp_AMSRE_SnowEM_LUT

END MODULE NESDIS_MW_SnowEM_LUT
