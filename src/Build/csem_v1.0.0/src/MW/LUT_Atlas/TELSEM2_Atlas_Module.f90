!
! CSEM_TELSEM2_Atlas
!
! This module is provided to allow users to use TELSEM2 land surface emissivity data sets by
! CSEM interfaces.
! TELSEM2 includes the monthly land surface emissivity atlas based on the multiple-year
! retrievals from SSMI and some trievals from other sensors. TELSEM2 is a generalized atlas
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



MODULE TELSEM2_Atlas_Module

  USE CSEM_Type_Kinds, ONLY: fp   => CSEM_fp, &
                             Byte => CSEM_Byte
  USE CSEM_Exception_Handler,ONLY: &
           SUCCESS, FAILURE, Display_Message
  USE TELSEM2_ATLAS_READER
  
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  
  LOGICAL,SAVE :: TELSEM2_atlas_init
  INTEGER,SAVE :: TELSEM2_atlas_month
  
  PUBLIC :: TELSEM2_Atlas_Initialized
  PUBLIC :: TELSEM2_Atlas_Setup
  PUBLIC :: TELSEM2_Atlas_Emiss
  PUBLIC :: TELSEM2_Atlas_Emiss_nChannels
  PUBLIC :: TELSEM2_Atlas_Close

CONTAINS

  FUNCTION TELSEM2_Atlas_Setup( &
               &  imonth,       &! in
               &  path,         &! in, optional
               &  mw_atlas_ver) &! in, optional
    RESULT( Error_Status )

    INTEGER, INTENT(IN)           :: imonth       ! month for which atlas data required
    CHARACTER(LEN=*),   INTENT(IN), OPTIONAL :: path         ! path to atlas data (if not the default)
    INTEGER, INTENT(IN), OPTIONAL :: mw_atlas_ver ! version of MW atlas to use
    
    ! Function result
    INTEGER :: Error_Status
     

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_TELSEM2_Atlas_Setup'
    ! Local variables
    CHARACTER(256) :: Message
    CHARACTER(LEN=300),SAVE :: fpath= './'
    
    ! set path
    IF (PRESENT(path)) fpath = TRIM(path)//'/'
    IF (PRESENT(mw_atlas_ver))print*,mw_atlas_ver
    ! ------  
    ! Set up  
    ! ------  
    Error_Status = SUCCESS
     
     
    IF (TELSEM2_atlas_init) THEN
        IF(imonth == TELSEM2_atlas_month) THEN
          Message = 'MW emissivity atlas already initialised.'
          RETURN
    ENDIF
        CALL TELSEM2_Atlas_Close( )
    ENDIF
    
    Error_Status=load_TELSEM2_atlas(TRIM(fpath),imonth)
    !PRINT*,'Initialising TELSEM2 Atlas of month ', imonth
    IF (Error_Status /= SUCCESS) THEN
        Error_Status = FAILURE
       PRINT*, 'Error initialising MW emissivity atlas TELSEM2 ...'
       RETURN
    END IF
    TELSEM2_atlas_month = imonth
    TELSEM2_atlas_init = .TRUE.

  END FUNCTION TELSEM2_Atlas_Setup

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       TELSEM2_Atlas_Emiss
!
! PURPOSE:
!       Function to provide the surface v-pol and h-pol emissivity of a specific 
!       microwave frequency over a land surface.
!
!       This function is dedicated to using TELSEM2 monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = TELSEM2_Atlas_Emiss(    &
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

  FUNCTION TELSEM2_Atlas_Emiss(  &
     & Frequency,               & ! input
     & Angle,                   & ! input
     & Latitude,                & ! input
     & Longitude,               & ! input
     & imonth,                  & ! input
     & Emissivity_H,            & ! output
     & Emissivity_V,            & ! output
     & resolution,              & ! in,  optional 
     & emis_std_v,              & ! out, optional 
     & emis_std_h,              & ! out, optional 
     & emis_cov,                & ! out, optional
     & stype)                   & ! out, optional 
    RESULT ( Error_Status )

 
    REAL(fp), INTENT(IN)      :: Frequency
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
    REAL(fp), INTENT(IN)      :: Angle
    REAL(fp), INTENT(OUT)     :: Emissivity_H,emissivity_V
    
    REAL(fp), INTENT(IN),  OPTIONAL :: resolution
    REAL(fp), INTENT(OUT), OPTIONAL :: emis_std_h,emis_std_v
    REAL(fp), INTENT(OUT), OPTIONAL :: emis_cov
    INTEGER,  INTENT(OUT), OPTIONAL :: stype

    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER(KIND=1) :: input_type = 0
    CHARACTER(LEN=100) :: Message
       
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'TELSEM2_Atlas_Emiss'
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    IF(.NOT. TELSEM2_Atlas_Initialized(imonth)) THEN
      PRINT*,'TELSEM2 Atlas has not been initilizaed for the month ', imonth
      Error_Status=TELSEM2_Atlas_Setup(imonth)
    ENDIF
    
    IF (PRESENT(emis_std_v)) emis_std_v = 0.0_fp
    IF (PRESENT(emis_std_h)) emis_std_h = 0.0_fp
    IF (PRESENT(emis_cov)) emis_cov = 0.0_fp
    IF (PRESENT(stype))  stype= 0
    
    IF ( PRESENT(resolution) ) input_type=IBSET(input_type, 0)
    IF ( PRESENT(emis_std_h) .AND. PRESENT(emis_std_v) .AND. PRESENT(emis_cov) ) THEN
        input_type=IBSET(input_type, 1)
    ENDIF

      
    IF (.NOT. TELSEM2_atlas_init) THEN
        Message='MW atlas not initialised'
        PRINT*,Message
        RETURN
    END IF
    
    ! Emissivity errors or covariances requested
    SELECT CASE (input_type)
       CASE (0_byte)
           CALL emis_interp_ind_sing(               &
                    latitude,                       &! in 
                    longitude,                      &! in
                    Angle,                          &! in
                    Frequency,                      &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    verb=0)! in
        CASE (1_byte)
            CALL emis_interp_int_sing(               &
                    latitude,                       &! in 
                    longitude,                      &! in
                    resolution,                     &! in 
                    Angle,                          &! in
                    Frequency,                      &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    verb=0)! in
        CASE (2_byte)
           CALL emis_interp_ind_sing(                &
                    latitude,                       &! in 
                    longitude,                      &! in
                    Angle,                          &! in
                    Frequency,                      &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    stdv=emis_std_v,                &! out
                    stdh=emis_std_h,                &! out
                    covvh=emis_cov,verb=0)! in
        CASE (3_byte)
           CALL emis_interp_int_sing(                &
                    latitude,                       &! in 
                    longitude,                      &! in
                    resolution,                     &! in 
                    Angle,                          &! in
                    Frequency,                      &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    stdv=emis_std_v,                &! out
                    stdh=emis_std_h,                &! out
                    covvh=emis_cov,verb=0)! in
  
       CASE DEFAULT
       
           PRINT*,'Wrong input...'
    END SELECT
           
  END FUNCTION TELSEM2_Atlas_Emiss

!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       TELSEM2_Atlas_Emiss_nChannels
!
! PURPOSE:
!       Function to provide the surface emissivity values of multiple 
!       microwave channels of a sensor over a land surface.
!
!       This function is dedicated to using TELSEM2 monthly atals.
!
! CALLING SEQUENCE:
!       Error_Status = TELSEM2_Atlas_Emiss(    &
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


  
  FUNCTION TELSEM2_Atlas_Emiss_nChannels(   &
     & Frequency,                          & ! input
     & Angle,                              & ! input
     & Latitude,                           & ! input
     & Longitude,                          & ! input
     & imonth,                             & ! input
     & crtm_polar_idx,                     & ! input 
     & n_Channels,                         & ! input
     & emissivity,                         & ! output
     & resolution,                         & ! in,  optional 
     & emis_std,                           & ! out, optional 
     & emis_cov,                           & ! out, optional
     & stype)                              & ! out, optional 
    RESULT ( Error_Status )

 
    INTEGER,  INTENT(IN)      :: n_Channels          ! number of channels 
    REAL(fp), INTENT(IN)      :: Frequency(n_Channels)
    REAL(fp), INTENT(IN)      :: Angle
    REAL(fp), INTENT(IN)      :: Latitude, Longitude
    INTEGER,  INTENT(IN)      :: imonth
    REAL(fp), INTENT(OUT)     :: emissivity(n_Channels)
    INTEGER,  INTENT(IN)      :: crtm_polar_idx(n_Channels)
    REAL(fp), INTENT(IN),  OPTIONAL :: resolution
    REAL(fp), INTENT(OUT), OPTIONAL :: emis_std(n_Channels)
    REAL(fp), INTENT(OUT), OPTIONAL :: emis_cov(n_Channels,n_Channels)
    INTEGER, INTENT(OUT) , OPTIONAL :: stype
    
    REAL(fp) :: emis_std2(2*n_Channels,2*n_Channels)
    REAL(fp) :: factor_v(n_Channels),factor_h(n_Channels)
    REAL(fp) :: emissivity_h(n_Channels),emissivity_v(n_Channels)

    ! Function result
    INTEGER :: Error_Status
    ! local
    INTEGER :: i, j, nchan

    !CHARACTER(LEN=100) :: Message
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'TELSEM2_Atlas_Emiss_nChannels'
    INTEGER(KIND=1) :: input_type = 0
     
       
    !-----------------------------
    ! Initialise output arguments
    !-----------------------------
    Error_Status = SUCCESS
    IF(.NOT. TELSEM2_Atlas_Initialized(imonth)) THEN
      PRINT*,'TELSEM2 Atlas has not been initilizaed for the month ', imonth
      Error_Status=TELSEM2_Atlas_Setup(imonth)
    ENDIF
    nchan = n_Channels 
     
    IF (PRESENT(emis_std)) emis_std = 0.0_fp
    IF (PRESENT(emis_cov)) emis_cov = 0.0_fp
    IF (PRESENT(stype))  stype= 0
   
    IF ( PRESENT(resolution) ) input_type=IBSET(input_type, 0)
    IF ( PRESENT(emis_std) .AND. PRESENT(emis_cov) ) THEN
        input_type=IBSET(input_type, 1)
    ENDIF
   
  
    SELECT CASE (input_type)
       CASE (0_byte)
             CALL emis_interp_ind_mult(                &
                    latitude,                       &! in 
                    longitude,                      &! in
                    Angle,                          &! in
                    Frequency,                      &! in
                    n_Channels,                     &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    verb=0)                          ! in
       CASE (1_byte)
           CALL emis_interp_int_mult(               &
                    latitude,                       &! in 
                    longitude,                      &! in
                    resolution,                     &! in 
                    Angle,                          &! in
                    Frequency,                      &! in
                    n_Channels,                     &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    verb=0             )             ! in

       CASE (2_byte)
           CALL emis_interp_ind_mult(                &
                    latitude,                       &! in 
                    longitude,                      &! in
                    Angle,                          &! in
                    Frequency,                      &! in
                    n_Channels,                     &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    std=emis_std2,                  &! out
                    verb=0)                          ! in
        
       CASE (3_byte)
          CALL emis_interp_int_mult(                &
                    latitude,                       &! in 
                    longitude,                      &! in
                    resolution,                     &! in 
                    Angle,                          &! in
                    Frequency,                      &! in
                    n_Channels,                     &! in
                    emissivity_v,                   &! out
                    emissivity_h,                   &! out
                    std=emis_std2,                  &! out
                    verb=0             )             ! in
       CASE DEFAULT
       
           PRINT*,'Wrong input...'
    END SELECT
    
    CALL Polarization_Factor(Angle, crtm_polar_idx, factor_v,factor_h)
    
    emissivity = emissivity_v*factor_v +  emissivity_h*factor_h
   
    IF ( PRESENT(emis_std) .AND. PRESENT(emis_cov) ) THEN
       ! Generate the stdv    
       ! Each dimension of std(:,:) has V-pol values for all channels followed by H-pol values
       ! for all channels.
       ! In this case std(:,:) contains covariances, so diagonal elements are variances (not stddev)
       ! but we want output emis_std(:) to be stddev
              
        DO j = 1, nchan
           emis_std(j) = SQRT(factor_v(j) * factor_v(j) * emis_std2(j,j)  + &
                    factor_v(j) * factor_h(j) * 2 * emis_std2(j,j+nchan)  + &
                    factor_h(j) * factor_h(j) * emis_std2(j+nchan,j+nchan))
        END DO

        ! Generate the covariances
        ! std(:,:) contains covariances. Combine the H- and V-pol values for each pair of
        ! channels to find the covariance of the combined emissivities in the two channels.
              
        DO i = 1, nchan
           DO j = 1, nchan
              emis_cov(i,j) = factor_v(i) * factor_v(j) * emis_std2(i,       j      ) + &
                    factor_v(i) * factor_h(j) * emis_std2(i,       j+nchan) + &
                    factor_h(i) * factor_v(j) * emis_std2(i+nchan, j      ) + &
                    factor_h(i) * factor_h(j) * emis_std2(i+nchan, j+nchan)
            END DO
         END DO
     END IF

  END FUNCTION TELSEM2_Atlas_Emiss_nChannels
            
  SUBROUTINE Polarization_Factor(Sensor_Zenith_Angle, crtm_polar_idx, factor_v,factor_h)
  
     REAL(fp) :: Sensor_Zenith_Angle
     INTEGER  :: crtm_polar_idx(:)
     REAL(fp) :: factor_v(:),factor_h(:)
        
     REAL(fp), PARAMETER :: PI = 3.141592653589793238462643383279_fp
     REAL(fp), PARAMETER :: DEGREES_TO_RADIANS = PI / 180.0_fp
     REAL(fp), PARAMETER :: EARTH_RADIUS     = 6370.0_fp  ! Mean earth radius 
     REAL(fp), PARAMETER :: SATELLITE_HEIGHT = 800.0_fp

    ! pol_v and pol_h give proportion of v and h pol to use in emissivity calculation
    ! pol_s3 adds the 3rd/4th stokes vectors 
     Real(fp), Parameter :: pol_v(3,7) = Reshape( &
       & (/ 0.5_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 1.0_fp, &
       & 0.0_fp, 1.0_fp, 0.0_fp, &
       & 1.0_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp  /), (/3,7/) )
     Real(fp), Parameter :: pol_h(3,7) = Reshape( &
       & (/ 0.5_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 1.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 1.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp, &
       & 1.0_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, 0.0_fp  /), (/3,7/) )
     Real(fp), Parameter :: pol_s3(0:1,7) = Reshape( &
       & (/ 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, &
       & 0.0_fp, 0.0_fp, &
       & 1.0_fp, 0.0_fp, &
       & 0.0_fp, 1.0_fp  /), (/2,7/) )
    ! mapping CRTM polarization flags to RTTOV definiations     
      INTEGER, Dimension(12),Parameter :: pol_map_index = &
    !   & (/ 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12/)  ! CRTM
        & (/ 1, 1, 6, 7, 4, 5, 6, 6, 2,  3,  7,  7/)  ! RTTOV
    
     INTEGER :: nchan,pol_id,j
     REAL(fp) :: sinzen, sinview, sinview_sq, cosview_sq 
   
     nchan=size(crtm_polar_idx)
     sinzen     = SIN(Sensor_Zenith_Angle * DEGREES_TO_RADIANS)
     sinview    = sinzen * EARTH_RADIUS/(EARTH_RADIUS + SATELLITE_HEIGHT) 
     sinview_sq = sinview * sinview
     cosview_sq = 1.0_fp - sinview_sq
            
     DO j = 1, nchan
        pol_id = pol_map_index(crtm_polar_idx(j))
        factor_v(j) = pol_v(1,pol_id) + pol_v(2,pol_id) * sinview_sq + pol_v(3,pol_id) * cosview_sq
        factor_h(j) = pol_h(1,pol_id) + pol_h(2,pol_id) * sinview_sq + pol_h(3,pol_id) * cosview_sq
     END DO

  END SUBROUTINE Polarization_Factor

  SUBROUTINE TELSEM2_Atlas_Close( )
    IF (TELSEM2_atlas_init) THEN
      PRINT*,'Clean TELSEM2 Atlas ...'
      CALL rttov_closemw_atlas
      TELSEM2_atlas_init = .FALSE.
    ENDIF
   
  END SUBROUTINE TELSEM2_Atlas_Close
  
  FUNCTION TELSEM2_Atlas_Initialized(imonth) RESULT( Atlas_Status )
     INTEGER :: imonth 
     LOGICAL :: Atlas_Status
     Atlas_Status=(TELSEM2_atlas_init .AND. (imonth == TELSEM2_atlas_month)) 
  END FUNCTION TELSEM2_Atlas_Initialized
END MODULE TELSEM2_Atlas_Module
