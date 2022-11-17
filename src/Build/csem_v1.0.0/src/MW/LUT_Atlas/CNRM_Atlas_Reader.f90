
MODULE CNRM_AMSUA_Reader

  ! Description:
  !   Data and routines for MW emissivity atlas METEO-FRANCE CNRM.
  ! 
  !   Following hypothesis is done:
  !   Altas is ordered in latitude, longitude:
  !     all longitudes for a given latitude, next latitude etc...
  !
  ! Karbou, F., E. Gérard, and F. Rabier, 2006, 
  ! Microwave land emissivity and skin temperature for AMSU-A and -B assimilation over land,
  ! Q. J. R.Meteorol. Soc., vol. 132, No. 620, Part A, pp.2333-2355(23). 
  !
  ! Copyright:
  !    This software was developed within the context of
  !    the EUMETSAT Satellite Application Facility on
  !    Numerical Weather Prediction (NWP SAF), under the
  !    Cooperation Agreement dated 25 November 1998, between
  !    EUMETSAT and the Met Office, UK, by one or more partners
  !    within the NWP SAF. The partners in the NWP SAF are
  !    the Met Office, ECMWF, KNMI and MeteoFrance.
  !
  !    Copyright 2010, EUMETSAT, All Rights Reserved.
  !
  ! Method:
  !
  ! Current Code Owner: SAF NWP
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  1.0      27/05/2010  Created, based on IFS/Arpege code (J. Hocking P. Brunel)
  !           22/11/2010  Change 30deg zenith angle threshold to 40deg (P. Brunel)
  !
  ! Code Description:
  !   Language:           Fortran 90.
  !   Software Standards: "European Standards for Writing and
  !     Documenting Exchangeable Fortran 90 Code".
  !
  ! Declarations:
  ! Modules used:
  !
 
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler,ONLY: &
           SUCCESS, FAILURE, Display_Message

  ! Disable implicit typing
  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE

  INTEGER :: cnrm_amsua_version=200   ! Version of atlas
  
  PUBLIC :: cnrm_amsua_setup 
  PUBLIC :: cnrm_amsua_emiss
  PUBLIC :: cnrm_amsua_emiss_multi
  PUBLIC :: cnrm_amsua_version 

  ! Atlas constants (extracted from MODULE YOMEMIS of IFS/ARPEGE cy36)
   
  ! Emissivity from an atlas
  ! - low angles (an1) and high angles (an2)
  ! the parameters should be updated in case the emissivity atlases change
 
  INTEGER , PARAMETER :: TOT_EM_AN1 = 352207
  REAL(fp) ,DIMENSION(12,TOT_EM_AN1)   ::  DATA_EMIS_AN1
  INTEGER :: lat_offset_an1(-90:90)
 
  INTEGER , PARAMETER :: TOT_EM_AN2 = 352207
  REAL(fp) ,DIMENSION(12,TOT_EM_AN2)   ::  DATA_EMIS_AN2
  INTEGER  :: lat_offset_an2(-90:90)

  
CONTAINS

!------------------------------------------
! Routines for initialising database
!------------------------------------------
  
  FUNCTION cnrm_amsua_setup( &
       path,     &      ! in
       imonth)   &      ! in
    RESULT( Error_Status )
                    
    IMPLICIT NONE
 
    CHARACTER (LEN=*),  INTENT(IN)  :: path
    INTEGER,            INTENT(IN)  :: imonth
 
     ! Function result
    INTEGER :: Error_Status
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'crtm_cnrmmwemis_init'
 
    INTEGER ::  file_id
    CHARACTER(LEN=300) :: fname
    CHARACTER(LEN=2) :: cmonth

    LOGICAL :: file_exists
    CHARACTER(256) :: Message

    Error_Status = SUCCESS
    WRITE(cmonth,"(i2.2)") imonth
     
    !----------------------------------------------------------------------------

    fname=path//'AMSU_CNRM_'//cmonth//'_atlas_an1.dat'
    Inquire(FILE=fname, EXIST=file_exists)
    IF ( .Not. file_exists )  Error_Status = FAILURE 
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE( Message,'("CNRM Atlas file ",a," not found")') TRIM(fname) 
      PRINT*,Message
      RETURN
    ENDIF
 
    file_id = Get_Lun()
    OPEN(file_id, FILE=TRIM(fname), STATUS='OLD',IOSTAT=Error_Status)
    IF ( Error_Status /= SUCCESS ) THEN
      Error_Status = FAILURE 
      WRITE( Message,'("Error Open CNRM Atlas file ",a)') TRIM(fname) 
      PRINT*,Message
      RETURN
    ENDIF

    READ(file_id,*, iostat=Error_Status) DATA_EMIS_AN1
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'io status while reading DATA_EMIS_AN1' 
      Error_Status = FAILURE 
      PRINT*,Message
      RETURN
    ENDIF
    CLOSE (file_id)
 
    CALL cal_lat_offset( data_emis_an1, lat_offset_an1 ) 

    fname=path//'AMSU_CNRM_'//cmonth//'_atlas_an2.dat'
    Inquire(FILE=fname, EXIST=file_exists)
    IF ( .Not. file_exists )  Error_Status = FAILURE 
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE( Message,'("CNRM Atlas file ",a," not found")') TRIM(fname) 
      PRINT*,Message
      RETURN
    ENDIF
 
    file_id = Get_Lun()
    OPEN(file_id, FILE=TRIM(fname), STATUS='OLD',IOSTAT=Error_Status)
    IF ( Error_Status /= SUCCESS ) THEN
      WRITE( Message,'("Error Open CNRM Atlas file ",a)') TRIM(fname) 
      Error_Status = FAILURE 
      PRINT*,Message
      RETURN
    ENDIF
 
    READ(file_id,*, iostat=Error_Status) DATA_EMIS_AN2
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'io status while reading DATA_EMIS_AN1' 
      Error_Status = FAILURE 
      PRINT*,Message
      RETURN
    ENDIF
    CLOSE (file_id)
 
    CALL cal_lat_offset( data_emis_an2, lat_offset_an2 )

  END FUNCTION cnrm_amsua_setup 
  

  FUNCTION cnrm_amsua_emiss(    &
       & latitude,          &! in  
       & longitude_in,      &! in
       & frequency,         &! in
       & zenangle,          &! in
       & emissivity_v,      &! out 
       & emissivity_h,      &! out 
       & pbats_veg)         & ! out 
    RESULT( Error_Status )
 
    REAL(fp),    INTENT(in)  :: latitude    ! degrees
    REAL(fp),    INTENT(in)  :: longitude_in! degrees
    REAL(fp),    INTENT(in)  :: frequency
    REAL(fp),    INTENT(in)  :: zenangle    ! zenith angle (degrees)
    REAL(fp),    INTENT(out) :: emissivity_v ! emissivity values 
    REAL(fp),    INTENT(out) :: emissivity_h ! emissivity values 
    REAL(fp),    INTENT(out) :: pbats_veg   ! vegetation
    ! Function result
    INTEGER :: Error_Status
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CRTM_cnrmmwemis'

    !CHARACTER(256) :: Message
   
    REAL(fp)    :: longitude ! in [-180,180]
    REAL(fp)    :: ze_ret(4)
    REAL(fp)    :: emissivity
    REAL(fp),  PARAMETER :: freqs(4)=(/23.8,31.4,50.3,89.0/)

    Error_Status = SUCCESS

    IF (longitude_in .GT. 180.0_fp) THEN
      longitude = longitude_in - 360.0_fp
    ELSEIF (longitude_in .LE. -180.0_fp) THEN
      longitude = longitude_in + 360.0_fp
    ELSE
      longitude = longitude_in
    ENDIF

    IF (zenangle > 40) THEN
      CALL LAND_AMSUA_AN2(latitude,longitude,ZE_RET,PBATS_VEG)
    ELSE
      CALL LAND_AMSUA_AN1(latitude,longitude,ZE_RET,PBATS_VEG)
    ENDIF

    IF (frequency <= freqs(1) ) THEN 
      emissivity=ZE_RET(1)
    ELSEIF (frequency > freqs(1) .AND. frequency <= freqs(2) ) THEN
      emissivity=ZE_RET(1) + (ZE_RET(2)-ZE_RET(1))* &
      (frequency-freqs(1))/(freqs(2)-freqs(1))
    ELSEIF (frequency > freqs(2) .AND. frequency <= freqs(3) ) THEN
      emissivity=ZE_RET(2) + (ZE_RET(3)-ZE_RET(2))* &
      (frequency-freqs(2))/(freqs(3)-freqs(2))
    ELSEIF (frequency > freqs(3) .AND. frequency < freqs(4) ) THEN
      emissivity=ZE_RET(3) 
    ELSE 
      emissivity=ZE_RET(4) 
    ENDIF
    emissivity_v = emissivity
    emissivity_h = emissivity

  END FUNCTION cnrm_amsua_emiss

  FUNCTION cnrm_amsua_emiss_multi(    &
      &  latitude,       &! in 
      &  longitude_in,   &! in
      &  frequency,      &
      &  zenangle,       &! in
      &  n_Channel,      &! in 
      &  emissivity,     &! out 
      &  pbats_veg)      & ! out 
    RESULT( Error_Status )
 
    INTEGER, INTENT(in)  :: n_Channel       ! number of channels
    REAL(fp),    INTENT(in)  :: latitude    ! degrees
    REAL(fp),    INTENT(in)  :: longitude_in   ! degrees
    REAL(fp),    INTENT(in)  :: zenangle    ! zenith angle (degrees)
    REAL(fp),    INTENT(in)  :: frequency(n_Channel)   
    REAL(fp),    INTENT(out) :: emissivity(n_Channel) ! emissivity values 
    REAL(fp),    INTENT(out) :: pbats_veg   ! vegetation
    ! Function result
    INTEGER :: Error_Status
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CNRM_AMSUA_Multi'
  
    !CHARACTER(256) :: Message
    REAL(fp)    :: longitude ! in [-180,180]
    REAL(fp)    :: ze_ret(4)
    REAL(fp),  PARAMETER :: freqs(4)=(/23.8,31.4,50.3,89.0/)
    INTEGER :: i


    Error_Status = SUCCESS

    IF (longitude_in .GT. 180.0_fp) THEN
      longitude = longitude_in - 360.0_fp
    ELSEIF (longitude_in .LE. -180.0_fp) THEN
      longitude = longitude_in + 360.0_fp
    ELSE
      longitude = longitude_in
    ENDIF

    IF (zenangle > 40) THEN
      CALL LAND_AMSUA_AN2(latitude,longitude,ZE_RET,PBATS_VEG)
    ELSE
      CALL LAND_AMSUA_AN1(latitude,longitude,ZE_RET,PBATS_VEG)
    ENDIF
    DO I = 1, n_Channel
      IF (frequency(I) <= freqs(1) ) THEN 
        emissivity(I)=ZE_RET(1)
      ELSEIF (frequency(I) > freqs(1) .AND. frequency(I) <= freqs(2) ) THEN
        emissivity(I)=ZE_RET(1) + (ZE_RET(2)-ZE_RET(1))* &
        (frequency(I)-freqs(1))/(freqs(2)-freqs(1))
      ELSEIF (frequency(I) > freqs(2) .AND. frequency(I) <= freqs(3) ) THEN
        emissivity(I)=ZE_RET(2) + (ZE_RET(3)-ZE_RET(2))* &
        (frequency(I)-freqs(2))/(freqs(3)-freqs(2))
      ELSEIF (frequency(I) > freqs(3) .AND. frequency(I) < freqs(4) ) THEN
        emissivity(I)=ZE_RET(3) 
      ELSE 
        emissivity(I)=ZE_RET(4) 
      ENDIF
    ENDDO

 
  END FUNCTION cnrm_amsua_emiss_multi


  SUBROUTINE cal_lat_offset(data_emis, lat_offset )
  ! calculates the offset in the data_emis array
  ! for all interger latitude value
    IMPLICIT NONE
    REAL(fp), INTENT(IN)  :: data_emis(:,:)
    INTEGER,  INTENT(OUT) :: lat_offset(-90:90)
    INTEGER :: last_lat, new_lat, jj

    lat_offset(:) = 1
    last_lat      = -90

    DO jj = 1, SIZE(DATA_EMIS(1,:))
      new_lat = int(DATA_EMIS(2,jj) + 90.0_fp) -90
      IF ( new_lat > last_lat ) THEN
        lat_offset(last_lat+1:new_lat) = jj
        last_lat = new_lat
      ENDIF
    ENDDO

    IF (last_lat < 90 ) THEN
      lat_offset(last_lat+1:) = lat_offset(last_lat) 
    ENDIF

  END SUBROUTINE cal_lat_offset

  !=======================================================================
  !
  ! Following are IFS/ARPEGE routines cleaned for:
  ! Dr_HOOK 
  ! Unused variables and arguments
  ! Commented code
  ! efficient loops with latitude limits (use of lat_offset_xxx)
  !
  !=======================================================================
  !
  SUBROUTINE LAND_AMSUA_AN1 (PLAT,PLON,E_RET,PBATS_VEG)

  ! ROUTINE.
  ! --------
  !    LAND_MF_AN1

  ! PURPOSE.
  ! --------
  !     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
  !     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

  ! INPUT.
  ! ------
  !     
  !     PLAT - observation latitude
  !     PLON - observation longitude


  ! OUTPUT.
  ! -------
  !     E_RET- surface emissivity     
  !     PBATS_VEG- BATS vegetation type

  ! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
  ! should be modified according to the emissivity atlas file format
    
      
    REAL(fp), INTENT(IN)  :: PLAT
    REAL(fp), INTENT(IN)  :: PLON
    REAL(fp), INTENT(OUT) :: E_RET (4)
    REAL(fp), INTENT(OUT) :: PBATS_VEG
 
    INTEGER indice, jj
    INTEGER lo, hi
    REAL(fp) r1,r2,r3,r4,dx,dy
          
  ! find within AMSU mean emissivity map the closet emissivity sets to the observation      
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    indice=0

    lo = 1
    hi = SIZE(DATA_EMIS_AN1(1,:)) 
    IF ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
    IF ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

    DO jj = lo, hi

      dx = abs(DATA_EMIS_AN1(1,jj)-PLON)
      dy = abs(DATA_EMIS_AN1(2,jj)-PLAT)
        
      IF ((dx <= 0.5) .AND. (dy <= 0.5)) THEN
        IF ((DATA_EMIS_AN1(3,jj) > 0) .AND. (DATA_EMIS_AN1(5,jj) > 0).AND. &
        & (DATA_EMIS_AN1(7,jj) > 0) .AND. (DATA_EMIS_AN1(9,jj) > 0)) THEN
          indice = indice+1
          PBATS_VEG=DATA_EMIS_AN1(11,jj)
          IF (PBATS_VEG == 21) PBATS_VEG=14
          r1 = r1 + DATA_EMIS_AN1(3,jj)
          r2 = r2 + DATA_EMIS_AN1(5,jj)
          r3 = r3 + DATA_EMIS_AN1(7,jj)
          r4 = r4 + DATA_EMIS_AN1(9,jj)
        ENDIF
      ENDIF

    ENDDO
    IF (indice > 0) THEN
      E_RET(1)=(r1/indice)
      E_RET(2)=(r2/indice)
      E_RET(3)=(r3/indice)
      E_RET(4)=(r4/indice)   

    ELSE
      PBATS_VEG=-9.09
      IF (abs(PLAT) .GT. 78) THEN
        E_RET(1)=0.75
        E_RET(2)=0.75
        E_RET(3)=0.75
        E_RET(4)=0.75
      ELSE
        E_RET(1)=0.95005
        E_RET(2)=0.95005
        E_RET(3)=0.95005
        E_RET(4)=0.95005
      END IF      
    ENDIF

  END SUBROUTINE LAND_AMSUA_AN1

  SUBROUTINE LAND_AMSUA_AN2 (PLAT,PLON,E_RET,PBATS_VEG)

  ! ROUTINE.
  ! --------
  !    LAND_MF_AN1

  ! PURPOSE.
  ! --------
  !     ROUTINE CALLED BY IFS OBSERVATION SCREENING MODULE TO 
  !     INTERPOLATE SURFACE EMISSIVITY FROM AN ATLAS

  ! INPUT.
  ! ------
  !     
  !     PLAT - observation latitude
  !     PLON - observation longitude


  ! OUTPUT.
  ! -------
  !     E_RET- surface emissivity     
  !     PBATS_VEG- BATS vegetation type

  ! F Karbou, N Bormann   09/05/07 : to caclucate land emissivities at mw frequencies   
  ! should be modified according to the emissivity atlas file format
         
      
    REAL(fp), INTENT(IN)  :: PLAT
    REAL(fp), INTENT(IN)  :: PLON
    REAL(fp), INTENT(OUT) :: PBATS_VEG
    REAL(fp), INTENT(OUT) :: E_RET (4)
     
    INTEGER indice, jj
    INTEGER lo, hi

    REAL(fp) r1,r2,r3,r4,dx,dy


  ! find within AMSU mean emissivity map the closet emissivity sets to the observation 
    r1 = 0.0
    r2 = 0.0
    r3 = 0.0
    r4 = 0.0
    indice=0

    lo = 1
    hi = SIZE(DATA_EMIS_AN2(1,:)) 
    IF ( nint(plat) > -90 ) lo = lat_offset_an1( nint(plat)-1 )
    IF ( nint(plat) <  90 ) hi = lat_offset_an1( nint(plat)+1 )

    DO jj = lo, hi
          
      dx = abs(DATA_EMIS_AN2(1,jj)-PLON)
      dy = abs(DATA_EMIS_AN2(2,jj)-PLAT)

      IF ((dx <= 0.5) .AND. (dy <= 0.5)) THEN
        IF ((DATA_EMIS_AN2(3,jj) > 0) .AND. (DATA_EMIS_AN2(5,jj) > 0).AND. &
        & (DATA_EMIS_AN2(7,jj) > 0) .AND. (DATA_EMIS_AN2(9,jj) > 0)) THEN
          indice = indice+1
          PBATS_VEG=DATA_EMIS_AN2(11,jj)
          IF (PBATS_VEG == 21) PBATS_VEG=14
          r1 = r1 + DATA_EMIS_AN2(3,jj)
          r2 = r2 + DATA_EMIS_AN2(5,jj)
          r3 = r3 + DATA_EMIS_AN2(7,jj)
          r4 = r4 + DATA_EMIS_AN2(9,jj)
        ENDIF
      ENDIF
    ENDDO
    
    IF (indice > 0) THEN
      E_RET(1)=(r1/indice)
      E_RET(2)=(r2/indice)
      E_RET(3)=(r3/indice)
      E_RET(4)=(r4/indice)       
    ELSE
      PBATS_VEG=-9.09
      IF (abs(PLAT) .GT. 78) THEN
        E_RET(1)=0.75
        E_RET(2)=0.75
        E_RET(3)=0.75
        E_RET(4)=0.75

      ELSE
        E_RET(1)=0.95005
        E_RET(2)=0.95005
        E_RET(3)=0.95005
        E_RET(4)=0.95005

      END IF                
    ENDIF         


  END SUBROUTINE LAND_AMSUA_AN2


  FUNCTION Get_Lun() RESULT( Lun )
    INTEGER :: Lun
    LOGICAL :: Is_Open
    LOGICAL :: Existence
    ! Initialise logical unit number
    Lun = 9

    ! Start open loop for Lun Search
    Lun_Search: DO
      Lun = Lun + 1
      INQUIRE( UNIT = Lun, EXIST = Existence )
      IF ( .NOT. Existence ) THEN
        Lun = -1
        EXIT Lun_Search
      END IF
      INQUIRE( UNIT = Lun, OPENED = Is_Open )
      IF ( .NOT. Is_Open  ) EXIT Lun_Search
    END DO Lun_Search

  END FUNCTION Get_Lun

END MODULE CNRM_AMSUA_Reader
