!
! NESDIS_SnowIR_PhyModel
!
! Module containing the NESDIS Snow emissivity model of infrared bands.
! Currently, this model is not well developed.
!
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, Fuzhong Weng, 02-28-2015
!                       Ming.Chen@noaa.gov
!                       Fuzhong.Weng@noaa.gov



MODULE NESDIS_SnowIR_PhyModel
 
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! IRibilities
  ! ------------
  PRIVATE
  PUBLIC  :: NESDIS_SnowIR_Emiss
  
  ! -----------------
  ! Module parameters
  ! -----------------
  REAL(fp), PARAMETER :: ZERO   = 0.0_fp
  REAL(fp), PARAMETER :: ONE    = 1.0_fp
  REAL(fp), PARAMETER :: TWO    = 2.0_fp
  REAL(fp), PARAMETER :: PI     = 3.141592653589793238462643_fp
  REAL(fp), PARAMETER :: EMISSH_DEFAULT = 0.25_fp
  REAL(fp), PARAMETER :: EMISSV_DEFAULT = 0.30_fp

CONTAINS


!
!--------------------------------------------------------------------------------
!
! NAME:
!       NESDIS_SnowIR_Emiss
!
! PURPOSE:
!       Function to simulate IRible emissivity over  Snow conditions.
!
! REFERENCES:
!       Weng, F., B. Yan, and N. Grody, 2001: "A infrared Snow emissivity model",
!         J. Geophys. Res., 106, 20, 115-20, 123
!
! CALLING SEQUENCE:
!              Error_Status = NESDIS_SnowIR_Emiss(              &
!                                  Frequency,                   &  ! Input
!                                  Angle,                       &  ! Input
!                                  Snow_Temperature,            &  ! Input
!                                  Soil_Temperature,            &  ! Input
!                                  Soil_Moisture_Content,       &  ! Input
!                                  Soil_Type,                   &  ! Input
!                                  Emissivity_H,                &  ! Output
!                                  Emissivity_V)                   ! Output
! 	    
!
! INPUT ARGUMENTS:
!                                  DIMENSION:  Scalar
!
! OUTPUT ARGUMENTS:
!         Emissivity_H:            The surface emissivity at a horizontal
!                                  polarization.
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!
!         Emissivity_V:            The surface emissivity at a vertical polarization.
!                                  UNITS:      N/A
!                                  TYPE:       REAL(fp)
!                                  DIMENSION:  Scalar
!
!
! INTERNAL ARGUMENTS:
!
! CREATION HISTORY:
!       Recoded by:     Ming Chen, Feb 2015
!                       Ming.Chen@noaa.gov
!
!------------------------------------------------------------------------------------------------------------

  FUNCTION NESDIS_SnowIR_Emiss(      &
      & Frequency,                   &  ! Input
      & Angle,                       &  ! Input
      & Snow_Temperature,            &  ! Input
      & Soil_Temperature,            &  ! Input
      & Soil_Moisture_Content,       &  ! Input
      & Soil_Type,                   &  ! Input
      & Emissivity_H,                &  ! Output
      & Emissivity_V)                &  ! Output
    RESULT (Error_Status)
     
    ! Arguments
    REAL(fp), INTENT(IN) :: Frequency
    REAL(fp), INTENT(IN) :: Angle
    REAL(fp), INTENT(IN) :: Soil_Moisture_Content
    REAL(fp), INTENT(IN) :: Soil_Temperature
    REAL(fp), INTENT(IN) :: Snow_Temperature
    INTEGER,  INTENT(IN) :: Soil_Type
    REAL(fp), INTENT(OUT):: Emissivity_V, Emissivity_H
    INTEGER :: Error_Status
    REAL(fp) :: rdummy, idummy
    
    Error_Status = 0
   
    IF(.FALSE.) THEN 
    rdummy =  Frequency
    rdummy =  Angle
    rdummy =  Soil_Moisture_Content
    rdummy =  Soil_Temperature
    rdummy =  Snow_Temperature
    idummy =  Soil_Type
    ENDIF
    
    Emissivity_H = EMISSH_DEFAULT
    Emissivity_V = EMISSV_DEFAULT
  
      
  END FUNCTION NESDIS_SnowIR_Emiss


END MODULE NESDIS_SnowIR_PhyModel
