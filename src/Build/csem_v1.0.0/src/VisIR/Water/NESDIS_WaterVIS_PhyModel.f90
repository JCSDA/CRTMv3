!
! NESDIS_WaterVIS_PhyModel
!
! Module containing the NESDIS physical Water emissivity model of visibal bands.
! Currently, this physical model is not well developed.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, Fuzhong Weng, 02-28-2015
!                       Ming.Chen@noaa.gov
!                       Fuzhong.Weng@noaa.gov



MODULE NESDIS_WaterVIS_PhyModel
 
  ! -----------------
  ! Enviroment set up
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  
  ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! VISibilities
  ! ------------
  PRIVATE
  PUBLIC  :: NESDIS_WaterVIS_Emiss
  
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
!       NESDIS_WaterVIS_Emiss
!
! PURPOSE:
!       Function to simulate VISible emissivity over Water conditions.
!
!!
! CALLING SEQUENCE:
!              Error_Status = NESDIS_WaterVIS_Emiss(              &
!                                    Frequency,                   &  ! Input
!                                    Angle,                       &  ! Input
!                                    Water_Temperature,           &  ! Input
!                                    Salinity,                    &  ! Input
!                                    Wind_Speed,                  &  ! Input
!                                    Wind_Direction,              &  ! Input
!                                    Emissivity_H,                &  ! Output
!                                    Emissivity_V)                &  ! Output
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

  FUNCTION NESDIS_WaterVIS_Emiss(     &
      & Frequency,                   &  ! Input
      & Angle,                       &  ! Input
      & Water_Temperature,           &  ! Input
      & Salinity,                    &  ! Input
      & Wind_Speed,                  &  ! Input
      & Wind_Direction,              &  ! Input
      & Emissivity_H,                &  ! Output
      & Emissivity_V)                &  ! Output
    RESULT (Error_Status)
     
    ! Arguments
    REAL(fp), INTENT(IN) :: Frequency
    REAL(fp), INTENT(IN) :: Angle
    REAL(fp), INTENT(IN) :: Water_Temperature
    REAL(fp), INTENT(IN) :: Wind_Speed
    REAL(fp), INTENT(IN) :: Wind_Direction
    REAL(fp), INTENT(IN) :: Salinity
    REAL(fp), INTENT(OUT):: Emissivity_V, Emissivity_H
    INTEGER :: Error_Status

    REAL(fp) :: rdummy
    
    Error_Status = 0
    
    IF(.FALSE.) THEN 
    rdummy =  Frequency
    rdummy =  Angle
    rdummy =  Water_Temperature
    rdummy =  Wind_Speed
    rdummy =  Wind_Direction
    rdummy =  Salinity
    ENDIF

    Emissivity_H = EMISSH_DEFAULT
    Emissivity_V = EMISSV_DEFAULT
  
      
  END FUNCTION NESDIS_WaterVIS_Emiss


END MODULE NESDIS_WaterVIS_PhyModel
