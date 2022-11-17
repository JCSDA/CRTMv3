!
! NESDIS_SnowVIS_PhyModel
!
! Module containing the NESDIS Snow emissivity model of visible bands.
! Currently, this model is not well developed.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, Fuzhong Weng, 02-28-2015
!                       Ming.Chen@noaa.gov
!                       Fuzhong.Weng@noaa.gov



MODULE NESDIS_SnowVIS_PhyModel
 
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC  :: NESDIS_SnowVIS_Emiss
  
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
!       NESDIS_SnowVIS_Emiss
!
! PURPOSE:
!       Function to simulate visible emissivity over Snow conditions.
!
! REFERENCES:
!
! CALLING SEQUENCE:
!              Error_Status = NESDIS_SnowVIS_Emiss(             &
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

  FUNCTION NESDIS_SnowVIS_Emiss(     &
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
  
      
  END FUNCTION NESDIS_SnowVIS_Emiss


END MODULE NESDIS_SnowVIS_PhyModel
