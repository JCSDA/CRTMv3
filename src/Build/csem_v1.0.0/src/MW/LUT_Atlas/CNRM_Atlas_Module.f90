!
! CSEM_CNRM_Atlas
!
! This module is provided to allow users to use CNRM land surface emissivity data sets by
! CSEM interfaces.
! CNRM data includes the monthly land surface emissivity atlas retrieved from AMSU-A,
! AMSU-B, SSMI, SSMIS, TMI and AMSRE (http://www.cnrm.meteo.fr/gmap/mwemis/get_data.html).
! Only the interfaces for the monthly AMSU-A atlas are implemented in this module. Similar 
! interfces may be implemented for the atlas retrieved from other sensors. 
! 
! The interfacing follows the general CSEM design where each emissivity model is required 
! to implement two interfaces with one to provide the h-pol and v-pol emissivity values 
! of a single frequecy and the other to provide the emissivity values of all the channels of
! a specific sensor. 
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!                       ming.chen@noaa.gov
!


    
MODULE CNRM_Atlas_Module

  USE CSEM_Type_Kinds, ONLY: fp   => CSEM_fp, &
                             Byte => CSEM_Byte
  USE CSEM_Exception_Handler,ONLY: &
           SUCCESS, FAILURE, Display_Message
  USE CNRM_AMSUA_Reader
  
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  
  PRIVATE
  
  LOGICAL,SAVE :: cnrm_amsua_init
  INTEGER,SAVE :: cnrm_amsua_month
 
  PUBLIC :: CNRM_Atlas_Initialized
  PUBLIC :: CNRM_Atlas_Setup
  PUBLIC :: CNRM_Atlas_Emiss
  PUBLIC :: CNRM_Atlas_Emiss_nChannels
  PUBLIC :: CNRM_Atlas_Close  

CONTAINS

  FUNCTION CNRM_Atlas_Setup( &
        imonth,       &! in
        path,         &! in, optional
        Atlas_ID,     &
        mw_atlas_ver) &! in, optional
    RESULT( Error_Status )
 
    CHARACTER(LEN=50),INTENT(IN), OPTIONAL :: Atlas_ID
    INTEGER, INTENT(IN)           :: imonth       ! month for which atlas data required
    CHARACTER(LEN=*),INTENT(IN), OPTIONAL :: path ! path to atlas data (if not the default)
    INTEGER, INTENT(IN), OPTIONAL :: mw_atlas_ver ! version of MW atlas to use
    
    ! Function result
    INTEGER :: Error_Status
     
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CNRM_Atlas_Setup'
    ! Local variables
    CHARACTER(256) :: Message
    CHARACTER(LEN=300),SAVE :: fpath
    
    ! set path
    IF (PRESENT(path)) fpath = TRIM(path)//'/'
    
    IF (PRESENT(mw_atlas_ver)) print*,mw_atlas_ver
    IF (PRESENT(Atlas_ID)) print*,Atlas_ID
    ! ------  
    ! Set up  
    ! ------  
    Error_Status = SUCCESS
     
    IF (cnrm_amsua_init) THEN
      IF(imonth == cnrm_amsua_month) THEN
         Message = 'MW emissivity atlas already initialised.'
         RETURN
      ENDIF
    ENDIF
    
    Error_Status=cnrm_amsua_setup(TRIM(fpath),imonth)
    !PRINT*,'Initializing CNRM Atlas ....'
    IF (Error_Status/= SUCCESS) THEN
      Error_Status = FAILURE
      PRINT*,'Error initialising CNRM emissivity atlas.'
      RETURN
    END IF
    cnrm_amsua_month = imonth
    cnrm_amsua_init = .TRUE.
  
 
  END FUNCTION CNRM_Atlas_Setup

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CNRM_Atlas_Emiss
!
! PURPOSE:
!       Function to provide the surface v-pol and h-pol emissivity of a specific 
!       microwave frequency over a land surface.
!
!       This function is dedicated to using CNRM monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = CNRM_Atlas_Emiss(    &
!                         Atlas_ID,              & ! input
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
!       Atlas_ID:        CNRM atlas from a specific sensor retrieval 
!                        UNITS:      N/A
!                        TYPE:       Character
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

  FUNCTION CNRM_Atlas_Emiss(    &
      &  Frequency,             & ! input
      &  Angle,                 & ! input
      &  Latitude,              & ! input
      &  Longitude,             & ! input
      &  imonth,               & ! input
      &  Emissivity_H,          & ! output
      &  Emissivity_V,          & ! output
      &  stype)                 & ! out, optional 
    RESULT ( Error_Status )

    
    REAL(fp), INTENT(IN)      :: Frequency
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
    REAL(fp), INTENT(IN)      :: Angle
    REAL(fp), INTENT(OUT)     :: Emissivity_h, Emissivity_v
    INTEGER,  INTENT(OUT), OPTIONAL :: stype

    ! Function result
    INTEGER :: Error_Status
    ! local
    !CHARACTER(LEN=100) :: Message
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CNRM_Atlas_Emiss'
    REAL(fp) :: pbats_veg   ! vegetation
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    IF(.NOT. CNRM_Atlas_Initialized(imonth)) THEN
      PRINT*,'CNRM Atlas has not been initilizaed for the month ', imonth
      Error_Status=CNRM_Atlas_Setup(imonth)
    ENDIF
    IF (PRESENT(stype))  stype= 0
    
    ! Emissivity errors or covariances requested
    Error_Status = cnrm_amsua_emiss( &
           latitude,        & ! in 
           longitude,       & ! in
           frequency,       & ! in
           Angle,           & ! in
           emissivity_v,    & ! out 
           emissivity_h,    & ! out 
           pbats_veg )      

  END FUNCTION CNRM_Atlas_Emiss



!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CNRM_Atlas_Emiss_nChannels
!
! PURPOSE:
!       Function to provide the surface emissivity values of multiple 
!       microwave channels of a sensor over a land surface.
!
!       This function is dedicated to using CNRM monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = CNRM_Atlas_Emiss(    &
!                         Atlas_ID,              & ! input
!                         Frequency,             & ! input
!                         Sensor_Zenith_Angle,   & ! input
!                         Latitude,              & ! input
!                         Longitude,             & ! input
!                         Emissivity,            & ! output
!                         stype)                   ! out, optional 
!
!
! INPUTS:
!       Atlas_ID:        CNRM atlas from a specific sensor retrieval 
!                        UNITS:      N/A
!                        TYPE:       Character
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


  FUNCTION CNRM_Atlas_Emiss_nChannels(  &
        Frequency,            & ! in
        Angle,                &
        Latitude,             & ! in
        Longitude,            & ! in
        imonth,               & ! input
        n_Channel,            & !
        emissivity,           & ! out
        stype )               & ! out, optional 
    RESULT ( Error_Status )

    INTEGER,  INTENT(IN)      :: n_Channel          ! number of channels 
    REAL(fp), INTENT(IN)      :: Frequency(n_Channel)
    REAL(fp), INTENT(IN)      :: Angle
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    REAL(fp), INTENT(OUT)     :: emissivity(n_Channel)
    INTEGER,  INTENT(IN)      :: imonth
    INTEGER, INTENT(OUT) , OPTIONAL :: stype

    ! Function result
    INTEGER :: Error_Status
    !CHARACTER(LEN=100) :: Message
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_Get_Emis'
    REAL(fp)  ::  pbats_veg 
    
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS

    IF(.NOT. CNRM_Atlas_Initialized(imonth)) THEN
      PRINT*,'CNRM Atlas has not been initilizaed for the month ', imonth
      Error_Status=CNRM_Atlas_Setup(imonth)
    ENDIF
    
    IF (PRESENT(stype))  stype= 0
 
    ! Emissivity errors or covariances requested
    Error_Status = cnrm_amsua_emiss_multi(                   &
        &   latitude,        &! in 
        &   longitude,       &! in
        &   frequency,       &! in
        &   Angle,           &! in
        &   n_Channel,       &! in 
        &   emissivity,      &! out 
        &   pbats_veg)        
  
  END FUNCTION CNRM_Atlas_Emiss_nChannels
  
  FUNCTION CNRM_Atlas_Initialized(imonth) RESULT( Atlas_Status )
     INTEGER :: imonth 
     LOGICAL :: Atlas_Status
     Atlas_Status = cnrm_amsua_init .AND. (imonth == cnrm_amsua_month) 
  END FUNCTION CNRM_Atlas_Initialized
  
  SUBROUTINE CNRM_Atlas_Close( )
    IF (cnrm_amsua_init) THEN
      PRINT*,'Clean CNRM Atlas ...'
      cnrm_amsua_init = .FALSE.
    ENDIF
  END SUBROUTINE CNRM_Atlas_Close
END MODULE CNRM_Atlas_Module
