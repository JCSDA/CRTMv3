!
! CSEM_NPOESS_LUT
!
! This module is provided to allow users to use the LUT of the land IR-VIS surface emissivity/reflectivity 
! spectrum with respect to NPOESS surface types.
! NPOESS LUT includes 24 sets of land surface VIS-IR emissivity/reflectivity spectra, with each spectrum 
! corresponds to one NPOESS surface type. The user provides the surface type or the lat-lon to get the  
! emissivity/reflectivity value at individula wavelength or the spectrum of the multiple wavelengths defined 
! by  the user.
!
! The interfacing follows the general CSEM design where each emissivity model is required 
! to implement two interfaces with one to provide individual emissivity value 
! of a single frequecy and the other to provide the emissivity values of all the channels of
! a specific sensor. 
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!



MODULE NPOESS_LUT_Module

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp, Byte => CSEM_Byte
  USE CSEM_Exception_Handler
  USE NPOESS_LUT_READER
  
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  
  INTERFACE NPOESS_LUT_Emiss
    MODULE PROCEDURE NPOESS_LUT_Emiss_Channel
    MODULE PROCEDURE NPOESS_LUT_Emiss_nChannels
  END INTERFACE NPOESS_LUT_Emiss


  PUBLIC :: NPOESS_LUT_Emiss
 

CONTAINS

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NPOESS_LUT_Emiss_Channel
!
! PURPOSE:
!!      Function to provide the surface  emissivity of a specific 
!       Wavenumber (cm-1) over a NPOESS land surface.
!
!       This function is dedicated to using NPOESS LUT.

! CALLING SEQUENCE:
!       Error_Status = NPOESS_LUT_Emiss(    &
!                         Wavenumber,            & ! input
!                         Emissivity,            & ! output
!                         stype,                 & ! input output
!                         Latitude,              & ! input
!                         Longitude)             & ! input
!
!
! INPUTS:
!       Wavenumber:      User's wavenumber
!                        UNITS:      cm-1
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Latitude:        User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Longitude:       User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       stype  :         User's surface type
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
!
! OUTPUTS:
!       Emissivity  :    Emissivity value
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION NPOESS_LUT_Emiss_Channel(  &
     & Wavenumber,                    & ! input
     & Emissivity,                    & ! output
     & stype,                         & ! input-output
     & Latitude,                      & ! optional input
     & Longitude)                     & ! optional input
    RESULT ( Error_Status )

 
    REAL(fp), INTENT(IN)      :: Wavenumber
    REAL(fp), INTENT(OUT)     :: Emissivity
    INTEGER,  INTENT(INOUT)   :: stype
    REAL(fp), INTENT(IN), OPTIONAL      :: Latitude
    REAL(fp), INTENT(IN), OPTIONAL      :: Longitude
    

    ! Function result
    INTEGER :: Error_Status
    ! local
    ! CHARACTER(LEN=100) :: Message
       
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NPOESS_LUT_Emiss'
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    
    IF(PRESENT(Latitude) .AND. PRESENT(Longitude)) &
       CALL Read_Stype_Map(Latitude,Longitude,stype)
    
    Error_Status =  Read_NPOESS_LUT(Wavenumber,emissivity,stype)
   
            
  END FUNCTION NPOESS_LUT_Emiss_Channel
  

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       NPOESS_LUT_Emiss_nChannels
!
! PURPOSE:
!       Function to provide the surface  emissivity of n_channel specific 
!       Wavenumbers (cm-1) over a NPOESS land surface.
!
!       This function is dedicated to using NPOESS LUT.
!
! CALLING SEQUENCE:
!       Error_Status = NPOESS_LUT_Emiss(    &
!                         Wavenumber,            & ! input
!                         Emissivity,            & ! output
!                         stype,                 & ! input output
!                         Latitude,              & ! input
!                         Longitude)             & ! input
!
!
! INPUTS:
!       Wavenumber:      User's wavenumber
!                        UNITS:      cm-1
!                        TYPE:       float
!                        DIMENSION:  Array
!                        ATTRIBUTES: INTENT(IN)
!
!       Latitude:        User's latitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       Longitude:       User's longitude
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       stype  :         User's surface type
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(INOUT)
!
! OUTPUTS:
!       Emissivity  :    Emissivity value
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Array
!                        ATTRIBUTES: INTENT(OUT)
!
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION NPOESS_LUT_Emiss_nChannels(  &
    & Wavenumber,                    & ! input
    & Emissivity,                    & ! output
    & stype,                         & ! input-output
    & Latitude,                      & ! optional input
    & Longitude )                    & ! optional input
    RESULT ( Error_Status )

 
    REAL(fp), INTENT(IN)      :: Wavenumber(:)
    REAL(fp), INTENT(OUT)     :: Emissivity(:)
    INTEGER,  INTENT(INOUT)   :: stype
    REAL(fp), INTENT(IN), OPTIONAL      :: Latitude
    REAL(fp), INTENT(IN), OPTIONAL      :: Longitude
    INTEGER :: i, n_Channels  

    ! Function result
    INTEGER :: Error_Status
    ! local
    ! CHARACTER(LEN=100) :: Message
       
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'NPOESS_LUT_Emiss'
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    n_Channels = size(Wavenumber)

    IF(PRESENT(Latitude) .AND. PRESENT(Longitude)) &
        CALL Read_Stype_Map(Latitude,Longitude,stype)
    DO i = 1, n_Channels 
       Error_Status =  Read_NPOESS_LUT(Wavenumber(i),emissivity(i),stype)      
    END DO
            
  END FUNCTION NPOESS_LUT_Emiss_nChannels
  
          

END MODULE NPOESS_LUT_MOdule
