!
! NESDIS_IceVIS_PhyModel
!
! Module containing the NESDIS physical Ice emissivity model of visible channels.
! Currently, the physical model is not well established. 
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, Fuzhong Weng, 02-28-2015
!                       Ming.Chen@noaa.gov
!                       Fuzhong.Weng@noaa.gov



MODULE NESDIS_IceVIS_PhyModel
 
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
  PUBLIC  :: NESDIS_IceVIS_Emiss
  
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
!       NESDIS_IceVIS_Emiss
!
! PURPOSE:
!       Function to physically simulate visible emissivity over Ice conditions.
!
! REFERENCES:
!
! CALLING SEQUENCE:
!              Error_Status = NESDIS_IceVIS_Emiss(               &
!                                  Frequency,                   &  ! Input
!                                  Angle,                       &  ! Input
!                                  Ice_Temperature,             &  ! Input
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

  FUNCTION NESDIS_IceVIS_Emiss(     &
      & Frequency,                   &  ! Input
      & Angle,                       &  ! Input
      & Ice_Temperature,             &  ! Input
      & Emissivity_H,                &  ! Output
      & Emissivity_V)                &  ! Output
    RESULT (Error_Status)
     
    ! Arguments
    REAL(fp), INTENT(IN) :: Frequency
    REAL(fp), INTENT(IN) :: Angle
    REAL(fp), INTENT(IN) :: Ice_Temperature
    REAL(fp), INTENT(OUT):: Emissivity_V, Emissivity_H
    INTEGER :: Error_Status
    REAL(fp) :: rdummy
    
    Error_Status = 0
    
    IF(.FALSE.) THEN 
    rdummy =  Frequency
    rdummy =  Angle
    rdummy =  Ice_Temperature
    ENDIF
    Emissivity_H = EMISSH_DEFAULT
    Emissivity_V = EMISSV_DEFAULT
  
      
  END FUNCTION NESDIS_IceVIS_Emiss


END MODULE NESDIS_IceVIS_PhyModel
