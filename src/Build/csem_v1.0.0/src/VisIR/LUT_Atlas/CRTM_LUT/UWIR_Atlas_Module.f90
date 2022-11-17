!
! CSEM_UWIR_Atlas
!
! This module is provided to allow users to use UWIR land surface emissivity data sets by
! CSEM interfaces.
! UWIR includes the monthly land surface emissivity atlas based on the multiple-year
! retrievals from SSMI and some trievals from other sensors. UWIR is a generalized atlas
! which means it may be applicable for different sensors besides SSMI. 
!! 
! The interfacing follows the general CSEM design where each emissivity model is required 
! to implement two interfaces with one to provide the h-pol and v-pol emissivity values 
! of a single frequecy and the other to provide the emissivity values of all the channels of
! a specific sensor. 
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!



MODULE UWIR_Atlas_Module

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler
  USE UWIR_ATLAS_READER
  
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  
  LOGICAL,SAVE :: UWIR_atlas_init
  INTEGER,SAVE :: UWIR_atlas_month
  
  PUBLIC :: UWIR_Atlas_Initialized
  PUBLIC :: UWIR_Atlas_Setup
  PUBLIC :: UWIR_Atlas_Emiss
  PUBLIC :: UWIR_Atlas_Emiss_nChannels
  PUBLIC :: UWIR_Atlas_Close

CONTAINS

  FUNCTION UWIR_Atlas_Setup( &
               &  imonth,       &! in
               &  path,         &! in, optional
               &  mw_atlas_ver) &! in, optional
    RESULT( Error_Status )

    INTEGER, INTENT(IN)           :: imonth       ! month for which atlas data required
    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path         ! path to atlas data (if not the default)
    INTEGER, INTENT(IN), OPTIONAL :: mw_atlas_ver ! version of MW atlas to use
    
    ! Function result
    INTEGER :: Error_Status
     

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_UWIR_Atlas_Setup'
    ! Local variables
    CHARACTER(256) :: Message
    CHARACTER(LEN=300),SAVE :: fpath = './'
  
    ! set path
 
    IF (PRESENT(path)) fpath = TRIM(path)//'/'
    IF (PRESENT(mw_atlas_ver)) print*,mw_atlas_ver
    ! ------  
    ! Set up  
    ! ------  
    Error_Status = SUCCESS
     
     
    IF (UWIR_atlas_init) THEN
        IF(imonth == UWIR_atlas_month) THEN
          Message = 'MW emissivity atlas already initialised.'
          RETURN
        ENDIF
        CALL UWIR_Atlas_Close( )
    ENDIF
    
    Error_Status=crtm_uwiremis_init( &
         &     fpath,             &! in
         &     imonth )        
    

    !PRINT*,'Initialising UWIR Atlas...'
    IF (Error_Status /= SUCCESS) THEN
        Error_Status = FAILURE
       PRINT*, 'Error initialising MW emissivity atlas...'
       RETURN
    END IF
    UWIR_atlas_month = imonth
    UWIR_atlas_init = .TRUE.

  END FUNCTION UWIR_Atlas_Setup

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       UWIR_Atlas_Emiss
!
! PURPOSE:
!       Function to provide the surface emissivity of a specific 
!       IR frequency over a land surface.
!
!       This function is dedicated to using UWIR monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = UWIR_Atlas_Emiss(    &
!                         Frequency,             & ! input
!                         Sensor_Zenith_Angle,   & ! input
!                         Latitude,              & ! input
!                         Longitude,             & ! input
!                         Emissivity_H,          & ! output
!                         Emissivity_V,          & ! output
!                         stype)                   ! out, optional 
!
!
! INPUTS:
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
!
! OUTPUTS:
!       Emissivity_H:    H-pol emissivity value
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!       Emissivity_V:    V-pol emissivity value
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(OUT)
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION UWIR_Atlas_Emiss(  &
     & Wavenumber,              & ! input
     & Latitude,                & ! input
     & Longitude,               & ! input
     & imonth,                  & ! input
     & Emissivity,              & ! output
     & emis_cov,                & ! out, optional
     & stype)                   & ! inout, optional
  RESULT ( Error_Status )

 
    REAL(fp), INTENT(IN)      :: Wavenumber
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
    REAL(fp), INTENT(OUT)     :: Emissivity
    
    REAL(fp), INTENT(OUT), OPTIONAL :: emis_cov
    INTEGER,  INTENT(INOUT), OPTIONAL :: stype

    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER  :: stype1=0
    REAL(fp) :: emis_cov1
    INTEGER  :: emis_flag

    CHARACTER(LEN=100) :: Message
       
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'UWIR_Atlas_Emiss'
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    IF(.NOT. UWIR_Atlas_Initialized(imonth)) THEN
      PRINT*,'UWIR Atlas has not been initilizaed for the month ', imonth
      Error_Status=UWIR_Atlas_Setup(imonth)
    ENDIF
    
    IF (PRESENT(emis_cov)) emis_cov = 0.0_fp
     
    IF (.NOT. UWIR_atlas_init) THEN
        Message='UW atlas not initialised'
        PRINT*,Message
        RETURN
    END IF
    
    IF (PRESENT(stype)) stype1=stype
   
    CALL csem_uwiremis_single( &
        &  wavenumber,           &! in 
        &  latitude,             &! in 
        &  longitude,            &! in
        &  stype1,               &! in
        &  emissivity,           &! out 
        &  emis_cov1,            &! out 
        &  emis_flag)             ! out 
     IF (PRESENT(emis_cov)) emis_cov = emis_cov1

  END FUNCTION UWIR_Atlas_Emiss

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       UWIR_Atlas_Emiss_nChannels
!
! PURPOSE:
!       Function to provide the surface emissivity values of multiple 
!       IR channels of a sensor over a land surface.
!
!       This function is dedicated to using UWIR monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = UWIR_Atlas_Emiss(    &
!                         Frequency,             & ! input
!                         Sensor_Zenith_Angle,   & ! input
!                         Latitude,              & ! input
!                         Longitude,             & ! input
!                         Emissivity,            & ! output
!                         stype)                   ! out, optional 
!
!
! INPUTS:
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
!
! OUTPUTS:
!       Emissivity  :    Emissivity values of sensor channels
!                        UNITS:      N/A
!                        TYPE:       float
!                        DIMENSION:  Array of n_Channels
!                        ATTRIBUTES: INTENT(OUT)
!
!
!:sdoc-:
!----------------------------------------------------------------------------------

  
  FUNCTION UWIR_Atlas_Emiss_nChannels(   &
     & Wavenumber,                         & ! input
     & Latitude,                           & ! input
     & Longitude,                          & ! input
     & imonth,                             & ! input
     & n_Channels,                         & ! input
     & emissivity,                         & ! output
     & emis_cov,                           & ! out, optional
     & stype)                              & ! inout, optional 
  RESULT ( Error_Status )
 
    INTEGER,  INTENT(IN)      :: n_Channels          ! number of channels 
    REAL(fp), INTENT(IN)      :: Wavenumber(n_Channels)
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
    REAL(fp), INTENT(OUT)     :: emissivity(n_Channels)
    REAL(fp), INTENT(OUT),   OPTIONAL :: emis_cov(n_Channels)
    INTEGER, INTENT(INOUT) , OPTIONAL :: stype

    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER  :: emis_flag
    INTEGER  :: stype1=0
    REAL(fp) :: emis_cov1(size(Wavenumber))

    CHARACTER(LEN=100) :: Message
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'UWIR_Atlas_Emiss_nChannels'
      
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    
    IF(.NOT. UWIR_Atlas_Initialized(imonth)) THEN
      PRINT*,'UWIR Atlas has not been initilizaed for the month ', imonth
      Error_Status=UWIR_Atlas_Setup(imonth)
    ENDIF
    
    IF (PRESENT(emis_cov)) emis_cov = 0.0_fp
   
    IF (.NOT. UWIR_atlas_init) THEN
        Message='UWIR atlas not initialised'
        PRINT*,Message
        RETURN
    END IF
    
    IF (PRESENT(stype)) stype1=stype
  
    CALL csem_uwiremis_multi( &
        &  wavenumber,           &! in 
        &  latitude,             &! in 
        &  longitude,            &! in
        &  stype1,               &! in
        &  n_Channels,           &! in 
        &  emissivity,           &! out 
        &  emis_cov1,            &! out 
        &  emis_flag)            ! out 

    IF (PRESENT(emis_cov)) emis_cov = emis_cov1

  END FUNCTION UWIR_Atlas_Emiss_nChannels
            

  SUBROUTINE UWIR_Atlas_Close( )
    IF (uwir_atlas_init) THEN
      !PRINT*,'Clean UWIR Atlas ...'
      CALL crtm_uwiremis_close_atlas()
      uwir_atlas_init = .FALSE. 
    ENDIF
   
  END SUBROUTINE UWIR_Atlas_Close
  
  FUNCTION UWIR_Atlas_Initialized(imonth) RESULT( Atlas_Status )
     INTEGER :: imonth 
     LOGICAL :: Atlas_Status
     Atlas_Status=(uwir_atlas_init .AND. (imonth == uwir_atlas_month)) 
  END FUNCTION UWIR_Atlas_Initialized
END MODULE UWIR_Atlas_Module
