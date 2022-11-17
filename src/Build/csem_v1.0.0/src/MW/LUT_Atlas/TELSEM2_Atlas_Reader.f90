! Description:
!> @file
!!   Subroutines for TELSEM2 MW emissivity atlas and interpolator
!
!> @brief
!!   Subroutines for TELSEM2 MW emissivity atlas and interpolator
!!
!! @details
!!   It is intended that this atlas be used via the RTTOV interface
!!   rather than by calling these subroutines directly.
!!
!!   Surface emissivity at microwaves to millimeter waves over Polar
!!   Regions: parameterization and evaluation with aircraft experiments
!!   D. Wang, C. Prigent, L. Kilic, S. Fox, R. C. Harlow, C. Jimenez,
!!   F. Aires, C. Grassotti, and F. Karbou
!!   Submitted to QJRMS
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
!    Copyright 2016, EUMETSAT, All Rights Reserved.
!
MODULE TELSEM2_ATLAS_READER

  USE CSEM_Type_Kinds, ONLY: jprb => CSEM_fp
  USE CSEM_Exception_Handler,ONLY: &
           SUCCESS, FAILURE, Display_Message


  IMPLICIT NONE


  PRIVATE
  PUBLIC :: telsem2_atlas_data,    &
            rttov_readmw_atlas,    &
            test_inputs,           &
            atlas2,                &
            load_telsem2_atlas,    &
            emis_interp_ind_sing,  &
            emis_interp_ind_mult,  &
            emis_interp_int_sing,  &
            emis_interp_int_mult,  &
            rttov_closemw_atlas

  !> Data type for TELSEM2 atlas data
  TYPE telsem2_atlas_data
    PRIVATE
    ! number of lines in the atlas
    INTEGER :: ndat
    ! number of channels in the atlas
    INTEGER :: nchan!=7
    ! name of the atlas (including version number)
    CHARACTER(len=22) :: name
    ! month of the atlas
    INTEGER :: month
    ! resolution of the atlas (equal-area)
    REAL(jprb) :: dlat!=0.25
    ! number of cells per lat band
    INTEGER, POINTER :: ncells(:) => NULL()
    ! the first cellnumber of lat band
    INTEGER, POINTER :: firstcell(:) => NULL()
    ! limits of the spatial domain (flaged if global)
    REAL(jprb) :: lat1, lat2, lon1, lon2
    ! Emissivities
    REAL(jprb), POINTER :: emis(:,:) => NULL() !emis(ndat,nchan)
    ! Correlations
    REAL(jprb), POINTER :: correl(:,:,:) => NULL() !correl(10,nchan,nchan)
    ! Emissivity uncertainties (std)
    REAL(jprb), POINTER :: emis_err(:,:) => NULL() !emis_err(ndat,nchan)
    ! Surface classe
    INTEGER, POINTER :: class1(:) => NULL()
    INTEGER, POINTER :: class2(:) => NULL()
    ! cellnumber of each of the pixels in the atlas
    INTEGER, POINTER :: cellnum(:) => NULL()
    ! "Correspondance" vector indicating that for the ith element, the j so that EMIS(j,...) is the emissivity of cellnumber i.
    INTEGER, POINTER :: correspondance(:) => NULL() !(660066)
  END TYPE telsem2_atlas_data
  TYPE(telsem2_atlas_data), SAVE :: atlas2


!=======================================================
!ROUTINES ==============================================
!=======================================================
CONTAINS

  !------------------------------------------------------------------
  !> Subroutine to check input variables: not used by RTTOV.
  !! @param[in]       month           month (1-12)
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            channel frequency (GHz)
  SUBROUTINE test_inputs(month,lat,lon,theta,freq)
    !======to test if inputs are correct
    INTEGER, INTENT (IN) :: month
    REAL(jprb), INTENT (IN)    :: lat, lon, theta, freq
    IF (month<1 .OR. month>12) THEN
       PRINT*,'Pb - month=', month
       STOP
    END IF
    IF (lat<-90 .OR. lat>90) THEN
       PRINT*,'Pb - latitude=',lat
       STOP
    END IF
    IF (lon<0 .OR. lon>360) THEN
       PRINT*,'Pb - longitude=',lon
       STOP
    END IF
    IF (theta<0 .OR. theta>90) THEN
       PRINT*,'Pb - angle=', theta
       STOP
    END IF
    IF (freq<10 .OR. freq>700) THEN
       PRINT*,'Pb - frequency=', freq
       STOP
    END IF
  END SUBROUTINE test_inputs


  !------------------------------------------------------------------
  !> Initialise a TELSEM2 atlas data structure. Atlas data may be initialised
  !! for a region of the globe defined by the lower and upper lat/lon limits,
  !! though this feature is not used by RTTOV.
  !! @param[in]       dir             path to atlas data files
  !! @param[in]       month           month of data to read (1-12)
  !! @param[in,out]   atlas           TELSEM2 atlas data structure to initialise
  !! @param[in]       verbose         flag to turn verbose output on/off
  !! @param[in,out]   err             status on exit
  !! @param[in]       lat1            latitude lower bound, optional
  !! @param[in]       lat2            latitude upper bound, optional
  !! @param[in]       lon1            longitude lower bound, optional
  !! @param[in]       lon2            longitude upper bound, optional
  SUBROUTINE rttov_readmw_atlas(dir,month,atlas, verbose,err,lat1,lat2,lon1,lon2)
    !======Read a monthly atlas
    INTEGER, INTENT (IN)                     :: month
    CHARACTER(len=*), INTENT (IN)            :: dir
    REAL(jprb),OPTIONAL, INTENT (IN)         :: lat1, lat2, lon1, lon2
    TYPE(telsem2_atlas_data), INTENT (INOUT) :: atlas
    LOGICAL, INTENT(IN)                :: verbose
    INTEGER, INTENT (INOUT)            :: err
    !----TRANSITORY VARIABLES
    INTEGER    :: ipos
    INTEGER    :: j,i
    INTEGER    :: iiin=21    ! unit for input
    INTEGER    :: ios
    CHARACTER  :: month2*2
    INTEGER    :: cellnum
    REAL(jprb) :: ssmi(2*7)
    REAL(jprb) :: lat=0.0,lon=0.0
    INTEGER    :: cur_class1,cur_class2
    INTEGER    :: take ! flag to take or not the pixel atlas for location constrains
    CHARACTER(LEN=128) :: msg
    INTEGER    :: lerr

    lerr = err ! to avoid "unused" warning
    !initialisation
    !
    IF (ASSOCIATED(atlas%emis)) THEN
       WRITE(*,*)"TELSEM2 atlas data structure already allocated"
    ENDIF
    WRITE(month2,'(I2.2)') month
    atlas%month=month
    atlas%lat1=-777
    atlas%lat2=-777
    atlas%lon1=-777
    atlas%lon2=-777
    IF (PRESENT(lat1)) atlas%lat1=lat1
    IF (PRESENT(lat2)) atlas%lat2=lat2
    IF (PRESENT(lon1)) atlas%lon1=MODULO(lon1,360.0_jprb)
    IF (PRESENT(lon2)) atlas%lon2=MODULO(lon2,360.0_jprb)
    !----ALLOCATION SPECIFIC TO SSMI ATLAS
    atlas%nchan=7
    atlas%name='ssmi_mean_emis_climato'
    atlas%dlat=0.25
    ALLOCATE(atlas%ncells(int(180./atlas%dlat)))
    ALLOCATE(atlas%firstcell(int(180./atlas%dlat)) )
    CALL equare(atlas%dlat,atlas%ncells,atlas%firstcell)
    !----INITIALISATIONS
!     atlas%correspondance(:)=-777
    !----DEFINING NUMBER OF DATA
    OPEN(iiin,file=dir//atlas%name//'_'//month2//'_cov_interpol_M2',status='old',form='formatted',iostat=ios)
    IF( ios /= 0 ) THEN

       WRITE(msg,'(a,i8)') 'Error opening the monthly input file ios= ',ios
    ENDIF
    READ(iiin,*,IOSTAT=ios) j
    atlas%ndat=j
    IF (verbose) THEN
      WRITE(msg,'(a,i10)') 'Nb data=',atlas%ndat
    ENDIF

    !----ALLOCATE VARIABLES
    ALLOCATE(atlas%emis(atlas%ndat,atlas%nchan))
    ALLOCATE(atlas%emis_err(atlas%ndat,atlas%nchan))
    ALLOCATE(atlas%class1(atlas%ndat))
    ALLOCATE(atlas%class2(atlas%ndat))
    ALLOCATE(atlas%cellnum(atlas%ndat))
    ALLOCATE(atlas%correspondance(660066))
    atlas%correspondance(:)=-777

    ipos=0
    DO j=1,atlas%ndat
       READ(iiin,*) cellnum, (ssmi(i),i=1,2*atlas%nchan),cur_class1,cur_class2

       take=1
       IF (PRESENT(lat1)) THEN
          CALL get_coordinates(cellnum,atlas,lat,lon)
          IF ((lat<lat1).OR.(lat>lat2).OR.(lon<lon1).OR.(lon>lon2)) THEN
             take=0
          END IF
       END IF

       IF ((cur_class1 > 0).AND.(cur_class2 > 0).AND.(ipos<atlas%ndat).AND.(take==1)) THEN
          ipos=ipos+1
          DO i=1,atlas%nchan
             atlas%emis(ipos,i)=ssmi(i)
             atlas%emis_err(ipos,i)=sqrt(ssmi(atlas%nchan+i))
          END DO
          atlas%cellnum(ipos)=cellnum
          atlas%class1(ipos)=cur_class1
          atlas%class2(ipos)=cur_class2
          atlas%correspondance(cellnum)=ipos
       END IF
    END DO
    atlas%ndat=ipos;
    CLOSE(iiin)

    !Correlation of uncertainties
    ALLOCATE(atlas%correl(10,atlas%nchan,atlas%nchan))
    OPEN(iiin,file=dir//'correlations', status='old', iostat=ios)
    IF( ios /= 0 ) THEN
       WRITE(msg,'(a,i8)') 'Error opening the correlations input file ios= ',ios
    ENDIF
    ! Corresponds to the class1
    DO i=1,10
       READ(iiin,*)
       DO j=1,7
          READ(iiin,'(7F5.2)') atlas%correl(i,j,:)
       END DO
    END DO
    CLOSE(iiin)

 
  END SUBROUTINE rttov_readmw_atlas


  !------------------------------------------------------------------
  FUNCTION load_telsem2_atlas(dir,month,lat1,lat2,lon1,lon2) RESULT( Error_Status )
    !======Read a monthly atlas
    INTEGER, INTENT (IN)                     :: month
    CHARACTER(len=*), INTENT (IN)            :: dir
    REAL(jprb),OPTIONAL, INTENT (IN)         :: lat1, lat2, lon1, lon2
    INTEGER  :: Error_Status
    !----TRANSITORY VARIABLES
    INTEGER    :: ipos
    INTEGER    :: j,i
    INTEGER    :: iiin=21    ! unit for input
    INTEGER    :: ios
    CHARACTER  :: month2*2
    INTEGER    :: cellnum
    REAL(jprb) :: ssmi(2*7)
    REAL(jprb) :: lat=0.0,lon=0.0
    INTEGER    :: cur_class1,cur_class2
    INTEGER    :: take ! flag to take or not the pixel atlas for location constrains
    CHARACTER(LEN=128) :: msg
    !initialisation
    !
    IF (ASSOCIATED(atlas2%emis)) THEN
       WRITE(*,*)"TELSEM2 atlas data structure already allocated"
    ENDIF
    WRITE(month2,'(I2.2)') month
    atlas2%month=month
    atlas2%lat1=-777
    atlas2%lat2=-777
    atlas2%lon1=-777
    atlas2%lon2=-777
    IF (PRESENT(lat1)) atlas2%lat1=lat1
    IF (PRESENT(lat2)) atlas2%lat2=lat2
    IF (PRESENT(lon1)) atlas2%lon1=MODULO(lon1,360.0_jprb)
    IF (PRESENT(lon2)) atlas2%lon2=MODULO(lon2,360.0_jprb)
    !----ALLOCATION SPECIFIC TO SSMI ATLAS
    atlas2%nchan=7
    atlas2%name='ssmi_mean_emis_climato'
    atlas2%dlat=0.25
    ALLOCATE(atlas2%ncells(int(180./atlas2%dlat)))
    ALLOCATE(atlas2%firstcell(int(180./atlas2%dlat)) )
    CALL equare(atlas2%dlat,atlas2%ncells,atlas2%firstcell)
    !----INITIALISATIONS
!     atlas2%correspondance(:)=-777
    !----DEFINING NUMBER OF DATA
    OPEN(iiin,file=dir//atlas2%name//'_'//month2//'_cov_interpol_M2',status='old',form='formatted',iostat=ios)
    IF( ios /= 0 ) THEN
       WRITE(msg,'(a,i8)') 'Error opening the monthly input file ios= ',ios
    ENDIF
    READ(iiin,*,IOSTAT=ios) j
    atlas2%ndat=j

    !----ALLOCATE VARIABLES
    ALLOCATE(atlas2%emis(atlas2%ndat,atlas2%nchan))
    ALLOCATE(atlas2%emis_err(atlas2%ndat,atlas2%nchan))
    ALLOCATE(atlas2%class1(atlas2%ndat))
    ALLOCATE(atlas2%class2(atlas2%ndat))
    ALLOCATE(atlas2%cellnum(atlas2%ndat))
    ALLOCATE(atlas2%correspondance(660066))
    atlas2%correspondance(:)=-777

    ipos=0
    DO j=1,atlas2%ndat
       READ(iiin,*) cellnum, (ssmi(i),i=1,2*atlas2%nchan),cur_class1,cur_class2

       take=1
       IF (PRESENT(lat1)) THEN
          CALL get_coordinates(cellnum,atlas2,lat,lon)
          IF ((lat<lat1).OR.(lat>lat2).OR.(lon<lon1).OR.(lon>lon2)) THEN
             take=0
          END IF
       END IF

       IF ((cur_class1 > 0).AND.(cur_class2 > 0).AND.(ipos<atlas2%ndat).AND.(take==1)) THEN
          ipos=ipos+1
          DO i=1,atlas2%nchan
             atlas2%emis(ipos,i)=ssmi(i)
             atlas2%emis_err(ipos,i)=sqrt(ssmi(atlas2%nchan+i))
          END DO
          atlas2%cellnum(ipos)=cellnum
          atlas2%class1(ipos)=cur_class1
          atlas2%class2(ipos)=cur_class2
          atlas2%correspondance(cellnum)=ipos
       END IF
    END DO
    atlas2%ndat=ipos;
    CLOSE(iiin)

    !Correlation of uncertainties
    ALLOCATE(atlas2%correl(10,atlas2%nchan,atlas2%nchan))
    OPEN(iiin,file=dir//'correlations', status='old', iostat=ios)
    IF( ios /= 0 ) THEN
       WRITE(msg,'(a,i8)') 'Error opening the correlations input file ios= ',ios
    ENDIF
    ! Corresponds to the class1
    DO i=1,10
       READ(iiin,*)
       DO j=1,7
          READ(iiin,'(7F5.2)') atlas2%correl(i,j,:)
       END DO
    END DO
    CLOSE(iiin)
 
   !Correlation of uncertainties
    ALLOCATE(atlas2%correl(10,atlas2%nchan,atlas2%nchan)) 
    !WRITE(0,*) 'reading classes...'
    OPEN(iiin,file=dir//'correlations', status='old')
    DO i=1,10
       READ(iiin,*) 
       DO j=1,7
          READ(iiin,'(7F5.2)') atlas2%correl(i,j,1),atlas2%correl(i,j,2), &
               & atlas2%correl(i,j,3),atlas2%correl(i,j,4),atlas2%correl(i,j,5),&
               & atlas2%correl(i,j,6),atlas2%correl(i,j,7)
       END DO
    END DO
    CLOSE(iiin)
    Error_Status = 0
  END FUNCTION load_telsem2_atlas

  !------------------------------------------------------------------
  !> Deallocate data in TELSEM2 atlas data structure
  !! @param[in,out]   atlas   TELSEM2 atlas data structure to deallocate
  SUBROUTINE rttov_closemw_atlas()
    !======free a monthly atlas

    IF (ASSOCIATED(atlas2%emis))           DEALLOCATE(atlas2%emis)
    IF (ASSOCIATED(atlas2%emis_err))       DEALLOCATE(atlas2%emis_err)
    IF (ASSOCIATED(atlas2%class1))         DEALLOCATE(atlas2%class1)
    IF (ASSOCIATED(atlas2%class2))         DEALLOCATE(atlas2%class2)
    IF (ASSOCIATED(atlas2%cellnum))        DEALLOCATE(atlas2%cellnum)
    IF (ASSOCIATED(atlas2%ncells))         DEALLOCATE(atlas2%ncells)
    IF (ASSOCIATED(atlas2%firstcell))      DEALLOCATE(atlas2%firstcell)
    IF (ASSOCIATED(atlas2%correl))         DEALLOCATE(atlas2%correl)
    IF (ASSOCIATED(atlas2%correspondance)) DEALLOCATE(atlas2%correspondance)
    NULLIFY(atlas2%emis,      &
            atlas2%emis_err,  &
            atlas2%class1,    &
            atlas2%class2,    &
            atlas2%class1,    &
            atlas2%cellnum,   &
            atlas2%firstcell, &
            atlas2%correl,    &
            atlas2%correspondance)
  END SUBROUTINE rttov_closemw_atlas

  !------------------------------------------------------------------
  !> Compute the number of cells in a latitude band and the first cell number
  !! of a latitude band.
  !! @param[in]       dlat            width of latitude bands
  !! @param[in,out]   ncells          number of cells in each band
  !! @param[in,out]   firstcell       number of first cell in each band
  SUBROUTINE equare(DLAT,NCELLS,FIRSTCELL)
    !======computes the number of cells in a lattitude band
    !======and the first cellnumber of a lattitude band
    IMPLICIT NONE
    ! EQUAL-AREA COMPUTATIONS
    ! I  TOTCEL             TOTAL NUMBER OF E.A. BOXES
    !  NCELLS(720)        NUMBER OF E.A. BOXES PER LAT ZONE
    !  TOCELL(1440,720)   BOX NUMBER OF E.A. BOX (1-659600)
    REAL(jprb), INTENT (IN) :: dlat
    INTEGER,INTENT (INOUT)  :: NCELLS(:)
    INTEGER,INTENT (INOUT)  :: FIRSTCELL(:)
    INTEGER :: maxlat,maxlon
    INTEGER :: TOTCEL
    INTEGER, ALLOCATABLE   :: TOCELL(:,:)
    INTEGER :: maxlt2,lat,icellr,lat1,lat2,numcel,numcls,lon
    REAL(jprb) :: rcells,rearth,pi,aezon,hezon,aecell,xlatb
    REAL(jprb) :: rlatb,rlate,xlate,htb,hte,htzone,rcelat,azone
    !DOUBLE PRECISION    :: dlongr,acell,asq,twopi,halfpi,rcellr
    INTEGER :: I

    maxlat=FLOOR(180./DLAT)
    maxlon=FLOOR(360./DLAT)

    ALLOCATE(TOCELL(maxlon,maxlat))

    !COMMON /EQUCOM/TOCEll,NCELLS
    REARTH = 6371.2d0
    PI     = 2.0d0 * ASIN(1.d0)
    !HALFPI=PI/2.d0
    !TWOPI=2.d0*PI
    RCELAT=(DLAT*PI)/180.d0
    TOTCEL=0d0
    ! CALCULATE HEIGHT AND AREA OF EQUATORIAL ZONE
    HEZON=REARTH*SIN(RCELAT)
    AEZON=2.d0*PI*REARTH*HEZON
    ! CALCULATE AREA OF EQUATORIAL CELL
    AECELL=(AEZON*DLAT)/360.d0
    ! print*,'aecell',aecell
    ! COMPUTE LONGITUDE ZONES FOR EACH LATITUDE ZONE
    MAXLT2=MAXLAT/2
    DO LAT=1,MAXLT2
       XLATB=(LAT-1)*DLAT
       XLATE=XLATB+DLAT
       RLATB=(2.0d0*PI*XLATB)/360.0d0
       RLATE=(2.0d0*PI*XLATE)/360.0d0
       !CALCULATE HEIGHTS OF LATB,LATE,ZONE
       HTB=REARTH*SIN(RLATB)
       HTE=REARTH*SIN(RLATE)
       HTZONE=HTE-HTB
       AZONE=2.d0*PI*REARTH*HTZONE
       !CALCULATE NUMBER OF CELLS
       RCELLS=AZONE/AECELL
       ICELLR=FLOOR(RCELLS+.50d0)
       !AUGMENT TOTAL # GRID CELLS (BOTH HEMISPHERES)
       TOTCEL=TOTCEL+2*ICELLR
       !RCELLR=ICELLR
       !DLONGR=360.0d0/RCELLR
       !ACELL=AZONE/RCELLR
       !ASQ=AZONE/MAXLON
       !CREATE TABLE OF LONGITUDES
       LAT1=LAT+MAXLT2
       LAT2=MAXLT2+1-LAT
       NCELLS(LAT1)=ICELLR
       NCELLS(LAT2)=ICELLR
    END DO
    NUMCEL = 0
    ! THROUGH EACH LAT ZONE
    DO LAT=1,MAXLAT
       NUMCLS = NCELLS(LAT)
       ! if (lat.eq.9) print*,lat,numcls,360./numcls
       ! LOOP THROUGH EACH LON FOR THIS LAT ZONE
       DO LON=1,NUMCLS
          NUMCEL = NUMCEL + 1
          ! FILL TOCELL ARRAY WITH STARTING CELL NUMBER
          TOCELL(LON,LAT) = NUMCEL
       END DO
    END DO
    DEALLOCATE(TOCELL)

    ! search for the first cellnumber in each lat band
    FIRSTCELL(1)=1
    DO I=2,maxlat
       FIRSTCELL(I)=FIRSTCELL(I-1)+NCELLS(I-1)
    END DO
  END SUBROUTINE equare


  !------------------------------------------------------------------
  !> Compute cell number given latitude and longitude
  !! @param[in]       lat         latitude
  !! @param[in]       lon         longitude
  !! @param[in]       atlas       TELSEM2 atlas data
  FUNCTION calc_cellnum(lat,lon,atlas)
    !======computes the cellnumber given the lat and long
    !======using the ncells included in the atlas
    IMPLICIT none
    REAL(jprb), INTENT (IN) :: lat, lon
    type(telsem2_atlas_data), INTENT (IN) :: atlas
    INTEGER :: calc_cellnum
    INTEGER :: ilat,ilon, i
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=min(int((lat+90.)/atlas%dlat)+1,size(atlas%ncells))
    ilon=int(lon/(360./atlas%ncells(ilat)))+1
    do i=1,ilat-1
       calc_cellnum=atlas%ncells(i)+calc_cellnum
    end do
    calc_cellnum=calc_cellnum+ilon
  END FUNCTION calc_cellnum

  !------------------------------------------------------------------
  !> Compute latitude and longitude given cell number
  !! @param[in]        cellnum     cell number
  !! @param[in]        atlas       TELSEM2 atlas data
  !! @param[out]       lat         latitude
  !! @param[out]       lon         longitude
  SUBROUTINE get_coordinates(cellnum,atlas,lat,lon)
    !======computes lat and lon given the cellnumber
    IMPLICIT none
    INTEGER, INTENT (IN) :: cellnum
    type(telsem2_atlas_data), INTENT (IN) :: atlas
    REAL(jprb), INTENT (OUT) :: lat, lon
    !INTEGER :: ilat,ilon
    INTEGER :: i
    REAL :: res_lat ! latitude resolution
    INTEGER :: index_lat_max,index_lat,index_lon

    res_lat = atlas%dlat

    index_lat_max  = int(180/res_lat)

    IF (cellnum >= atlas%firstcell(index_lat_max)) THEN
       index_lat = index_lat_max
       lat =(index_lat - 0.5)*res_lat - 90
       index_lon = cellnum - atlas%firstcell(index_lat_max)+1
       lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
    ELSE
       DO i=1,index_lat_max-1
          IF ( (cellnum>=atlas%firstcell(i)) .AND. (cellnum<atlas%firstcell(i+1)) ) THEN
             index_lat = i
             lat = (index_lat - 0.5)*res_lat- 90
             ! cout << i << "  " << *lat <<endl;
             index_lon = cellnum - atlas%firstcell(i)+1
             lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
          END IF
       END DO
    END IF
  END SUBROUTINE get_coordinates






  !------------------------------------------------------------------
  !> Compute cell numbers given latitude and longitude and resolution (i.e.
  !! returns multiple cells for footprints larger than native atlas
  !! resolution)
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       resol           resolution
  !! @param[in]       atlas           TELSEM2 atlas data
  !! @param[in,out]   cell_num_mult   disassociated pointer, allocated and populated on exit
  !! @param[out]      nb_cell         number of cells in output
  SUBROUTINE calc_cellnum_mult(lat,lon,resol,atlas,cell_num_mult,nb_cell)
    !======computes the cellnumbers given the lat, long and resol
    !======using the ncells and firstcells included in the atlas
    IMPLICIT none
    REAL(jprb), INTENT (IN) :: lat, lon
    ! lat(-90 90)  lon(-180 180) of the pixel center
    ! resol = spatial resolution of the new grid
    type(telsem2_atlas_data), INTENT (IN) :: atlas
    INTEGER, POINTER :: cell_num_mult(:)
    INTEGER, INTENT (OUT)  :: nb_cell
    REAL(jprb), INTENT(IN) :: resol
    INTEGER :: ilat,ilon, i
    INTEGER :: calc_cellnum
    !--------new variables
    INTEGER :: nbcel,cell(200)
    !                 200 is here the maximum number of cells inside the "big" pixel
    INTEGER :: i2lon,i3lon,i2lat,nbreslat,nbreslon,i4lon
    REAL(jprb) :: maxlat
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=MIN(INT((lat+90.)/atlas%dlat)+1,SIZE(atlas%ncells))
    ilon=INT(lon/(360./atlas%ncells(ilat)))+1
    DO i=1,ilat-1
       calc_cellnum=atlas%ncells(i)+calc_cellnum
    END DO
    calc_cellnum=calc_cellnum+ilon

    ! search for the nbcel pixels that are in the box of resolution
    ! 'resol' around lat and lon.
    maxlat=180./atlas%dlat
    !maxlon=360./atlas%dlat

    nbcel=1
    cell(nbcel)=calc_cellnum
    nbreslat=INT(resol/atlas%dlat/2.)
    IF (nbreslat>=1) THEN

       DO i2lat=ilat-nbreslat,ilat+nbreslat
          IF ( (i2lat<1).OR.(i2lat>INT(maxlat)) ) THEN
             !PRINT*,'Congratulation, you reached the pole!'
          ELSE
             IF (ABS((i2lat-.5)*atlas%dlat-90-lat)<=resol/2) THEN
                i2lon=INT(lon/(360./atlas%ncells(i2lat)))+1
                nbreslon=INT(resol/(360./(1.*atlas%ncells(i2lat)))/2.)
                DO i3lon=i2lon-nbreslon,i2lon+nbreslon
                   IF (mod(ABS((i3lon-.5_jprb)*(360._jprb/atlas%ncells(i2lat))-lon),360._jprb)<=resol/2.) THEN
                      nbcel=nbcel+1
                      i4lon=i3lon
                      IF (i3lon<1.) i4lon=atlas%ncells(i2lat)+i3lon
                      IF (i3lon>=atlas%ncells(i2lat)) i4lon=i3lon-atlas%ncells(i2lat)
                      cell(nbcel)=atlas%firstcell(i2lat)+ i4lon-1
                      IF (cell(nbcel)==cell(1)) nbcel=nbcel-1
                   END IF
                END DO
             END IF
          END IF
       END DO

    END IF
    nb_cell=nbcel
    ALLOCATE(cell_num_mult(int(nb_cell)))
    DO i=1,nb_cell
       cell_num_mult(i)=cell(i)
    END DO
  END SUBROUTINE calc_cellnum_mult

  !------------------------------------------------------------------
  !> Return emissivities for a single channel at the native atlas resolution
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            frequency
  !! @param[in]       atlas           TELSEM2 atlas data
  !! @param[out]      ev              V-pol emissivity
  !! @param[out]      eh              H-pol emissivity
  !! @param[out]      stdv            V-pol emissivity standard deviation, optional
  !! @param[out]      stdh            H-pol emissivity standard deviation, optional
  !! @param[out]      covvh           H-/V-pol emissivity covariance, optional
  !! @param[in]       verb            switch for verbose output
  SUBROUTINE emis_interp_ind_sing(lat,lon,theta,freq,ev,eh,stdv,stdh,covvh,verb)
    !======interpolates emissivities
    !          IND: individual cellnumber atlas
    !          SING: singular channel
    Implicit none
    REAL(jprb), INTENT (IN)               :: lat, lon, theta, freq
    REAL(jprb), INTENT (OUT)              :: ev, eh   !resultat de l'interpolation
    REAL(jprb), OPTIONAL, INTENT (OUT)    :: stdv, stdh, covvh
    INTEGER, INTENT (IN)            :: verb
    INTEGER    :: ipos, i, j, i2, j2
    INTEGER    :: cellnum
    REAL(jprb) :: ev_a(3),eh_a(3)      !emissivity in the atlas

    ! Need to calculate the H-/V-pol covariance as well
    REAL(KIND=jprb) :: FIM(2,3*2)           !Frequency Interpolation Matrix
    REAL(KIND=jprb) :: trans_std(3*2,2)     !transpose Frequency Interpolation Matrix
    REAL(KIND=jprb) :: a,b,c                !frequency interpolation coefficients
    REAL(KIND=jprb) :: stdd(3),d            !dummy variables
    REAL(KIND=jprb) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(KIND=jprb) :: std(2,2)             !calculated covariance matrix
    REAL(KIND=jprb) :: new_FIM(3*2,2)       !transpose Frequency Interpolation Matrix

    !Initialisations
    stdd = 0.
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    IF (PRESENT(covvh)) covvh=0

    cellnum=calc_cellnum(lat,MODULO(lon,360.0_jprb),atlas2)
    ipos=atlas2%correspondance(cellnum)
    IF (ipos>0) THEN
       ev_a(1)=atlas2%emis(ipos,1)
       eh_a(1)=atlas2%emis(ipos,2)
       ev_a(2)=atlas2%emis(ipos,4)
       eh_a(2)=atlas2%emis(ipos,5)
       ev_a(3)=atlas2%emis(ipos,6)
       eh_a(3)=atlas2%emis(ipos,7)
       CALL emis_interp(theta,freq,atlas2%class1(ipos),atlas2%class2(ipos),ev_a,eh_a,ev,eh)
       IF (PRESENT(stdv) .OR. PRESENT(stdh) .OR. PRESENT(covvh)) THEN
          !compute covariance matrix on the atlas frequencies
          DO i=1,6
             IF (i>2) THEN
                i2=i+1
             ELSE
                i2=i
             END IF
             DO j=1,6
                IF (j>2) THEN
                   j2=j+1
                ELSE
                   j2=j
                END IF
                cov(i,j)=atlas2%correl(atlas2%class1(ipos),i2,j2)* &
                    & (atlas2%emis_err(ipos,i2)*atlas2%emis_err(ipos,j2))
             END DO
          END DO

          !compute the Frequency Linear Matrix
          CALL interp_freq2(stdd(1),stdd(2),stdd(3),freq,atlas2%class2(ipos),d,a,b,c)
          FIM(1,1)=a
          FIM(1,2)=b
          FIM(1,3)=c

          DO j=1,3
            FIM(1,3+j)=0
            FIM(2,j)=0
            FIM(2,j+3)=FIM(1,j)
          END DO

          new_FIM=TRANSPOSE(FIM)
          trans_std=MATMUL(cov,new_FIM)
          std=MATMUL(FIM,trans_std)
          stdv=SQRT(std(1,1))
          stdh=SQRT(std(2,2))
          covvh=std(1,2)
       END IF
    END IF
    IF (verb==1) THEN
       WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
    END IF
  END SUBROUTINE emis_interp_ind_sing


  !------------------------------------------------------------------
  !> Return emissivities for multiple channels at the native atlas resolution.
  !! Each dimension of the covariance matrix std(:,:) has V-pol values for all
  !! channels followed by H-pol values for all channels.
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            frequencies
  !! @param[in]       n_chan          number of channels
  !! @param[in]       atlas           TELSEM2 atlas data
  !! @param[out]      ev              V-pol emissivities
  !! @param[out]      eh              H-pol emissivities
  !! @param[out]      std             Covariance matrix, optional
  !! @param[in]       verb            switch for verbose output
  SUBROUTINE emis_interp_ind_mult(lat,lon,theta,freq,n_chan,ev,eh,std,verb)
    !======interpolates emissivities
    !          IND: individual cellnumber atlas
    !          MULT: multiple channels
    IMPLICIT none
    INTEGER, INTENT (IN)            :: n_chan
    REAL(jprb), INTENT (IN)               :: lat, lon, theta, freq(:)
    REAL(jprb), INTENT (OUT)              :: ev(:), eh(:)   !resultat de l'interpolation
    REAL(jprb), OPTIONAL, INTENT (OUT)    :: std(:,:)
    INTEGER, INTENT (IN)            :: verb

    INTEGER :: ipos, i, j, i2, j2
    INTEGER :: cellnum
    REAL(jprb) :: ev_a(3),eh_a(3)      !emissivity in the atlas
    REAL(jprb) :: stdv_a(3),stdh_a(3)  !std in the atlas
    REAL(jprb) :: std2
    REAL(jprb) :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL(jprb) :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    REAL(jprb) :: a,b,c                !frequency interpolation coefficients
    REAL(jprb) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas

    REAL(jprb) :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL(jprb) :: correlation(2*n_chan,2*n_chan) !correlation

    !initialisations
    DO i=1,n_chan
       ev(i)=0
       eh(i)=0
    END DO
    IF (PRESENT(std)) THEN
       DO i=1,2*n_chan
          DO j=1,2*n_chan
             std(i,j)=0
          END DO
       END DO
    END IF

    !localization
    cellnum=calc_cellnum(lat,MODULO(lon,360.0_jprb),atlas2)
    ipos=atlas2%correspondance(cellnum)

    !compute emissivities
    IF (ipos>0) THEN
       ev_a(1)=atlas2%emis(ipos,1)
       eh_a(1)=atlas2%emis(ipos,2)
       ev_a(2)=atlas2%emis(ipos,4)
       eh_a(2)=atlas2%emis(ipos,5)
       ev_a(3)=atlas2%emis(ipos,6)
       eh_a(3)=atlas2%emis(ipos,7)
       stdv_a(1)=atlas2%emis_err(ipos,1)  !don't matter to put
       stdh_a(1)=atlas2%emis_err(ipos,2)  !real data here, we
       stdv_a(2)=atlas2%emis_err(ipos,4)  !just need a, b, c
       stdh_a(2)=atlas2%emis_err(ipos,5)  !don't matter to put
       stdv_a(3)=atlas2%emis_err(ipos,6)  !real data here, we
       stdh_a(3)=atlas2%emis_err(ipos,7)  !just need a, b, c

       DO i=1,n_chan
          CALL emis_interp(theta,freq(i),atlas2%class1(ipos),atlas2%class2(ipos),ev_a,eh_a,ev(i),eh(i))
       END DO

       IF (PRESENT(std)) THEN !compute uncertainties on emissivities
          !compute covariance matrix on the atlas frequencies
          DO i=1,6
             IF (i>2) THEN
                i2=i+1
             ELSE
                i2=i
             END IF
             DO j=1,6
                IF (j>2) THEN
                   j2=j+1
                ELSE
                   j2=j
                END IF
                cov(i,j)=atlas2%correl(atlas2%class1(ipos),i2,j2)* &
                    & (atlas2%emis_err(ipos,i2)*atlas2%emis_err(ipos,j2))
             END DO
          END DO

          !compute the Frequency Linear Matrix
          DO i=1,n_chan
             CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),atlas2%class2(ipos),std2,a,b,c)
             FIM(i,1)=a
             FIM(i,2)=b
             FIM(i,3)=c
          END DO
          DO i=1,n_chan
             DO j=1,3
                FIM(i,3+j)=0
                FIM(n_chan+i,j)=0
                FIM(n_chan+i,j+3)=FIM(i,j)
             END DO
          END DO
          new_FIM=TRANSPOSE(FIM)
          trans_std=MATMUL(cov,new_FIM)
          !std=MATMUL(FIM,MATMUL(cov,new_FIM))
          std=MATMUL(FIM,trans_std)
       END IF
       !verbose ?
       IF (verb==1) THEN
          WRITE(0,*) 'Cellnum=',cellnum,' lat=',lat,' lon=',lon,' ipos=',ipos
          WRITE(0,*)    'emis_SSMI_V(1)=',ev_a(1),'/',stdv_a(1)
          WRITE(0,*)    'emis_SSMI_V(2)=',ev_a(2),'/',stdv_a(2)
          WRITE(0,*)    'emis_SSMI_V(3)=',ev_a(3),'/',stdv_a(3)
          WRITE(0,*)    'emis_SSMI_H(1)=',eh_a(1),'/',stdh_a(1)
          WRITE(0,*)    'emis_SSMI_H(2)=',eh_a(2),'/',stdh_a(2)
          WRITE(0,*)    'emis_SSMI_H(3)=',eh_a(3),'/',stdh_a(3)
          DO i=1,n_chan
             WRITE(0,*) 'IND(',freq(i),') EV=', ev(i),' EH=', eh(i)
          END DO
          WRITE(0,*) 'Matrice de covariance 6x6'
          DO i=1,6
             WRITE(0,'(I3,a,6F8.5)') i,': ',(cov(i,j),j=1,6)
          END DO
          WRITE(0,*) 'Matrice de correl'
          DO i=1,6
             WRITE(0,'(I3,a,6F8.5)') i,': ',(std(i,j),j=1,6)
          END DO
          !the Frequency Linear Matrix
          DO i=1,n_chan
             IF (PRESENT(std)) THEN
                WRITE(0,'(3F8.5)') FIM(i,1),FIM(i,2),FIM(i,3)
             END IF
          END DO
          WRITE(0,*) 'FIM'
          DO i=1,2*n_chan
             WRITE(0,'(I3,a,6F8.5)') i,': ',(FIM(i,j),j=1,6)
          END DO
          WRITE(0,*) 'STD'
          DO i=1,2*n_chan
             WRITE(0,'(I3,a,10F8.5)') i,': ',(std(i,j),j=1,2*n_chan)
          END DO
          WRITE(0,*) 'Matrice de correl'
          DO i=1,2*n_chan
             DO j=1,2*n_chan
                correlation(i,j)=std(i,j)/SQRT(std(i,i)*std(j,j))
             END DO
             WRITE(0,'(I3,a,10F8.5)') i,': ',(correlation(i,j),j=1,2*n_chan)
          END DO
       END IF  !verb
    END IF   !ipos
  END SUBROUTINE emis_interp_ind_mult


  !------------------------------------------------------------------
  !> Return emissivities for a single channel at the user-specified resolution
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       resol           resolution
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            frequency
  !! @param[in]       atlas           TELSEM2 atlas data
  !! @param[out]      ev              V-pol emissivity
  !! @param[out]      eh              H-pol emissivity
  !! @param[out]      stdv            V-pol emissivity standard deviation, optional
  !! @param[out]      stdh            H-pol emissivity standard deviation, optional
  !! @param[out]      covvh           H-/V-pol emissivity covariance, optional
  !! @param[in]       verb            switch for verbose output
  SUBROUTINE emis_interp_int_sing(lat,lon,resol,theta,freq,ev,eh,stdv,stdh,covvh,verb)
    !======interpolates emissivities
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          SING: singular channel
    IMPLICIT none
    REAL(jprb), INTENT (IN)               :: lat, lon, theta, freq
    REAL(jprb), INTENT (OUT)              :: ev, eh
    REAL(jprb), OPTIONAL, INTENT (OUT)    :: stdv,stdh,covvh
    INTEGER, INTENT (IN)            :: verb
    REAL(jprb), INTENT (IN)               :: resol

    INTEGER          :: ipos, i, j, i2, j2
    REAL(jprb)       :: ev_a(3),eh_a(3)     !emissivities in the atlas
    INTEGER          :: ii
    INTEGER          :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL(jprb)       :: ev_mean, eh_mean
    INTEGER          :: inumb

    ! Need to calculate the H-/V-pol covariance as well
    REAL(KIND=jprb) :: FIM(2,3*2)           !Frequency Interpolation Matrix
    REAL(KIND=jprb) :: trans_std(3*2,2)     !transpose Frequency Interpolation Matrix
    REAL(KIND=jprb) :: a,b,c                !frequency interpolation coefficients
    REAL(KIND=jprb) :: stdd(3),d            !dummy variables
    REAL(KIND=jprb) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(KIND=jprb) :: std2(2,2)            !calculated covariance matrix
    REAL(KIND=jprb) :: std_mean(2,2)        !calculated average covariance matrix
    REAL(KIND=jprb) :: new_FIM(3*2,2)       !transpose Frequency Interpolation Matrix

    IF(verb /= 0)PRINT*

    !initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    IF (PRESENT(covvh)) covvh=0
    stdd(:) = 0
    !computes the list of cells that need to be integrated in the atlas
    CALL calc_cellnum_mult(lat,MODULO(lon,360.0_jprb),resol,atlas2,cellnum_mult,nb_cell)
    ev_mean=0
    eh_mean=0
    std_mean(:,:)=0
    inumb=0
    DO ii=1,nb_cell
       ipos=atlas2%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          ev_a(1)=atlas2%emis(ipos,1)
          eh_a(1)=atlas2%emis(ipos,2)
          ev_a(2)=atlas2%emis(ipos,4)
          eh_a(2)=atlas2%emis(ipos,5)
          ev_a(3)=atlas2%emis(ipos,6)
          eh_a(3)=atlas2%emis(ipos,7)
          ev=0
          eh=0
          CALL emis_interp(theta,freq,atlas2%class1(ipos),atlas2%class2(ipos),ev_a,eh_a,ev,eh)
          !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas2%correspondance(cellnum_mult(ii)),' ',ev,' ',eh
          ev_mean=ev_mean+ev
          eh_mean=eh_mean+eh
          IF (PRESENT(stdv) .OR. PRESENT(stdh) .OR. PRESENT(covvh)) THEN
             !compute covariance matrix on the atlas frequencies
             DO i=1,6
                IF (i>2) THEN
                   i2=i+1
                ELSE
                   i2=i
                END IF
                DO j=1,6
                   IF (j>2) THEN
                      j2=j+1
                   ELSE
                      j2=j
                   END IF
                   cov(i,j)=atlas2%correl(atlas2%class1(ipos),i2,j2)* &
                       & (atlas2%emis_err(ipos,i2)*atlas2%emis_err(ipos,j2))
                END DO
             END DO

             !compute the Frequency Linear Matrix
             CALL interp_freq2(stdd(1),stdd(2),stdd(3),freq,atlas2%class2(ipos),d,a,b,c)
             FIM(1,1)=a
             FIM(1,2)=b
             FIM(1,3)=c
             DO j=1,3
                FIM(1,3+j)=0
                FIM(2,j)=0
                FIM(2,j+3)=FIM(1,j)
             END DO
             new_FIM=TRANSPOSE(FIM)
             trans_std=MATMUL(cov,new_FIM)
             std2=MATMUL(FIM,trans_std)
             std_mean(:,:)=std_mean(:,:)+std2(:,:)
          END IF
       END IF
    END DO

    IF (inumb>0) THEN
       ev=ev_mean/inumb
       eh=eh_mean/inumb
       IF (PRESENT(stdv)) THEN
          stdv=SQRT(std_mean(1,1)/inumb)
       ENDIF
       IF (PRESENT(stdh)) THEN
          stdh=SQRT(std_mean(2,2)/inumb)
       ENDIF
       IF (PRESENT(covvh)) THEN
          covvh=std_mean(1,2)/inumb
       END IF
    END IF
    DEALLOCATE(cellnum_mult)
  END SUBROUTINE emis_interp_int_sing



  !------------------------------------------------------------------
  !> Return emissivities for multiple channels at the user-specified resolution.
  !! Each dimension of the covariance matrix std(:,:) has V-pol values for all
  !! channels followed by H-pol values for all channels.
  !! @param[in]       lat             latitude
  !! @param[in]       lon             longitude
  !! @param[in]       resol           resolution
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            frequencies
  !! @param[in]       n_chan          number of channels
  !! @param[in]       atlas           TELSEM2 atlas data
  !! @param[out]      ev              V-pol emissivities
  !! @param[out]      eh              H-pol emissivities
  !! @param[out]      std             Covariance matrix, optional
  !! @param[in]       verb            switch for verbose output
  SUBROUTINE emis_interp_int_mult(lat,lon,resol,theta,freq,n_chan,ev,eh,std, verb)
    !======interpolates emissivities
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          MULT: multiple channel
    IMPLICIT none
    INTEGER, INTENT (IN)            :: n_chan
    REAL(jprb), INTENT (IN)               :: lat, lon, theta, freq(:)
    REAL(jprb), INTENT (OUT)              :: ev(:), eh(:)
    REAL(jprb), OPTIONAL, INTENT (OUT)    :: std(:,:)
    INTEGER, INTENT (IN)            :: verb
    REAL(jprb), INTENT(IN)                :: resol

    INTEGER :: ipos
    REAL(jprb) :: ev_a(3),eh_a(3)     !emissivities in the atlas
    REAL(jprb) :: stdv_a(3) !std emissivities in the atlas
    INTEGER :: ii, i2, j2
    INTEGER :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL(jprb) :: ev_mean(n_chan), eh_mean(n_chan)
    REAL(jprb) :: std_mean(2*n_chan,2*n_chan)
    REAL(jprb) :: std2(2*n_chan,2*n_chan)
    REAL(jprb) :: ev2, eh2
    INTEGER :: inumb, i,j

    REAL(jprb) :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL(jprb) :: a,b,c                !frequency interpolation coefficients
    REAL(jprb) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(jprb) :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL(jprb) :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    !REAL    :: correlation(2*n_chan,2*n_chan) !correlation

    IF(verb /= 0)PRINT*
    !initialisations
    ev(:)=0
    eh(:)=0
    ev_mean(:)=0
    eh_mean(:)=0
    IF (PRESENT(std)) THEN
       std_mean(:,:)=0
       std(:,:)=0
    END IF

    !computes the list of cells that need to be integrated in the atlas
    CALL calc_cellnum_mult(lat,MODULO(lon,360.0_jprb),resol,atlas2,cellnum_mult,nb_cell)

    inumb=0
    DO ii=1,nb_cell
       ipos=atlas2%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          ev_a(1)=atlas2%emis(ipos,1)
          eh_a(1)=atlas2%emis(ipos,2)
          ev_a(2)=atlas2%emis(ipos,4)
          eh_a(2)=atlas2%emis(ipos,5)
          ev_a(3)=atlas2%emis(ipos,6)
          eh_a(3)=atlas2%emis(ipos,7)
          stdv_a(1)=atlas2%emis_err(ipos,1)  !don't matter to put
          !stdh_a(1)=atlas2%emis_err(ipos,2)  !real data here, we
          stdv_a(2)=atlas2%emis_err(ipos,4)  !just need a, b, c
          !stdh_a(2)=atlas2%emis_err(ipos,5)  !don't matter to put
          stdv_a(3)=atlas2%emis_err(ipos,6)  !real data here, we
          !stdh_a(3)=atlas2%emis_err(ipos,7)  !just need a, b, c
          DO i=1,n_chan
             !initialisations
             ev2=0
             eh2=0
             CALL emis_interp(theta,freq(i),atlas2%class1(ipos),atlas2%class2(ipos),ev_a,eh_a,ev2,eh2)
             !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas2%correspondance(cellnum_mult(ii)),' ',ev2,' ',eh2
             ev_mean(i)=ev_mean(i)+ev2
             eh_mean(i)=eh_mean(i)+eh2
          END DO

          IF (PRESENT(std)) THEN !compute uncertainties on emissivities
             !compute covariance matrix on the atlas frequencies
             DO i=1,6
                IF (i>2) THEN
                   i2=i+1
                ELSE
                   i2=i
                END IF
                DO j=1,6
                   IF (j>2) THEN
                      j2=j+1
                   ELSE
                      j2=j
                   END IF
                   cov(i,j)=atlas2%correl(atlas2%class1(ipos),i2,j2)* &
                       & (atlas2%emis_err(ipos,i2)*atlas2%emis_err(ipos,j2))
                END DO
             END DO

             !compute the Frequency Linear Matrix
             DO i=1,n_chan
                CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),atlas2%class2(ipos),ev2,a,b,c) !ev2=dummy variable
                FIM(i,1)=a
                FIM(i,2)=b
                FIM(i,3)=c
             END DO
             DO i=1,n_chan
                DO j=1,3
                   FIM(i,3+j)=0
                   FIM(n_chan+i,j)=0
                   FIM(n_chan+i,j+3)=FIM(i,j)
                END DO
             END DO
             new_FIM=TRANSPOSE(FIM)
             trans_std=MATMUL(cov,new_FIM)
             !std2=MATMUL(FIM,MATMUL(cov,new_FIM))
             std2=MATMUL(FIM,trans_std)
             std_mean(:,:)=std_mean(:,:)+std2(:,:)
          END IF !std ou pas
       END IF !ipos
    END DO  !nb_cell

    IF (inumb>0) THEN
       DO i=1,n_chan
          ev(i)=ev_mean(i)/inumb
          eh(i)=eh_mean(i)/inumb
          !WRITE(0,*) 'INT: EV=',ev,' EH=',eh
       END DO
       IF (PRESENT(std)) THEN
          std(:,:)=std_mean(:,:)/inumb
       END IF
    END IF

    DEALLOCATE(cellnum_mult)
  END SUBROUTINE emis_interp_int_mult

  !------------------------------------------------------------------
  !> Linear interpolation of emissivity up to 85GHz. Above 85 GHz, for most
  !! surfaces constant emissivity except when  emiss85>emiss37 for 3 classes
  !! of sea ice (class2 11 12 13) and for surface water (class2 10) (ie when
  !! water like behaviour). For most surface types, very limited frequency
  !! dependence or not enough convincing evidence for a reliable frequency 
  !! dependence.
  !! @param[in]       emiss19         emissivity at 19GHz
  !! @param[in]       emiss37         emissivity at 37GHz
  !! @param[in]       emiss85         emissivity at 85GHz
  !! @param[in]       f               frequency
  !! @param[in]       class2          surface class
  !! @param[out]      emiss           interpolated emissivity
  !! @param[out]      an              interpolation weight at 19GHz, optional
  !! @param[out]      bn              interpolation weight at 37GHz, optional
  !! @param[out]      cn              interpolation weight at 85GHz, optional
  SUBROUTINE interp_freq2(emiss19,emiss37,emiss85,f,class2,emiss,an,bn,cn)
    !======Linear interpolation of emissivity up to 85 GHz.
    !======Above 85 GHz, for most surfaces constant emissivity
    !======except when emiss85>emiss37 for 3 classes of sea ice (class2 11 12 13)
    !======and for surface water (class2 10) (ie when water like behavior)
    !======For most surface types, very limited frequency dependence or not enough
    !======convincing evidence for a reliable frequency dependence.
    IMPLICIT none
    REAL(jprb), INTENT (IN) :: emiss19,emiss37,emiss85,f
    INTEGER, INTENT (IN) :: class2
    REAL(jprb), INTENT (OUT) :: emiss
    REAL(jprb), OPTIONAL, INTENT (OUT) :: an, bn, cn
    REAL(jprb) :: a, b, c,rapport43_32(4),rapport54_43(4)

    DATA rapport43_32 / 0.62_jprb, 0.37_jprb, 0.46_jprb, 0.63_jprb /
    DATA rapport54_43 / 0.30_jprb, 0.60_jprb, 0.47_jprb, 0.35_jprb /

    ! class2 1 to 5 idem as TELSEM
    ! class2 10=water
    ! class2 11 to 16 sea ice
    ! class2 17 to 22 snow and continental ice
    ! Frequency interpolation above 85 GHz only for classes 10 to 13 (water like classes)

    IF (f<=19.35) THEN
       a=1
       b=0
       c=0
       emiss = emiss19
    ELSE IF (f<=37.) THEN
       a=(37.-f  )/(37.-19.35)
       b=(f-19.35)/(37.-19.35)
       c=0
       emiss = a*emiss19+b*emiss37
    ELSE IF (f<85.5) THEN
       a=0
       b=(85.5-f )/(85.5-37)
       c=(f-37   )/(85.5-37)
       emiss = b*emiss37+c*emiss85
    ELSE
        a=0
        b=0
        c=1
        emiss = emiss85
        IF ((class2>9).AND.(class2<14).AND.(emiss85>emiss37)) THEN
           IF (f<=150) THEN
               emiss=emiss85+(f-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
           ELSE IF (f<=190) THEN
               emiss=emiss85+(150.-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
               emiss=emiss+(f-150.)*((emiss-emiss85)/(150.-85.5))*rapport54_43(class2-9)
           ELSE
               emiss=emiss85+(150.-85.5)*((emiss85-emiss37)/(85.5-37.))*rapport43_32(class2-9)
               emiss=emiss+(190.-150.)*((emiss-emiss85)/(150.-85.5))*rapport54_43(class2-9)
           END IF
           IF (emiss>1) emiss=1
        END IF
    END IF

    IF (PRESENT(an)) an=a
    IF (PRESENT(bn)) bn=b
    IF (PRESENT(cn)) cn=c
  END SUBROUTINE interp_freq2

 !------------------------------------------------------------------
  !> Interpolate atlas emissivities in angle and frequency.
  !! @param[in]       theta           zenith angle
  !! @param[in]       freq            frequency
  !! @param[in]       class1          surface class
  !! @param[in]       class2          surface class
  !! @param[in]       ev              V-pol emissivities at 19, 37, 85 GHz from the atlas
  !! @param[in]       eh              H-pol emissivities at 19, 37, 85 GHz from the atlas
  !! @param[out]      emiss_interp_v  interpolated V-pol emissivity
  !! @param[out]      emiss_interp_h  interpolated H-pol emissivity
   SUBROUTINE emis_interp(theta,freq,class1,class2,ev,eh,emiss_interp_v,emiss_interp_h)
    !======interpolation of emissivities for angle, freq
    IMPLICIT none
    INTEGER, INTENT (IN)       :: class1,class2
    REAL(jprb), INTENT (IN)    :: theta, freq, ev(3), eh(3)
    REAL(jprb), INTENT (OUT)   :: emiss_interp_v, emiss_interp_h
    REAL(jprb)    :: e0, theta0, theta53 , emiss_scal_v(3), emiss_scal_h(3)
    REAL(jprb)    :: S1_v, S1_h, S2_v, S2_h, S_v, S_h, a0, a1, a2, a3, b0, b1
    REAL(jprb)    :: b2, b3, em53_v, em53_h, emtheta_v, emtheta_h
    REAL(jprb)    :: a0_k0(3,10),a0_k1(3,10),a0_k2(3,10)
    REAL(jprb)    :: a0_eveh(3,10),a1_eveh(3,10),a2_eveh(3,10),a3_eveh(3,10)
    REAL(jprb)    :: b0_eveh(3,10),b1_eveh(3,10),b2_eveh(3,10),b3_eveh(3,10)
    INTEGER :: j
    !COMMON /EMISSIVITE/emiss_interp_v,emiss_interp_h
    data a0_k0/0.11509_jprb,0.091535_jprb,0.34796_jprb,0.10525_jprb,0.16627_jprb,0.24434_jprb, &
         & 0.29217_jprb,0.23809_jprb,0.28954_jprb,0.17516_jprb,0.19459_jprb,0.28697_jprb, &
         & 0.10521_jprb,0.12126_jprb,0.30278_jprb,0.18212_jprb,0.19625_jprb,0.14551_jprb, &
         & -0.19202_jprb,0.5411_jprb,0.03739_jprb,0.10292_jprb,0.5486_jprb,-0.058937_jprb, &
         & -0.022672_jprb,0.44492_jprb,-0.058448_jprb,-0.33894_jprb,-0.17621_jprb,0.14742_jprb/
    data a0_k1/0.61168_jprb,0.59095_jprb,0.7918_jprb,0.60271_jprb,0.69213_jprb,0.62218_jprb, &
         &  0.32728_jprb,0.34334_jprb,0.37062_jprb,0.51217_jprb,0.4491_jprb,0.50101_jprb, &
         & 0.48913_jprb,0.41932_jprb,0.29734_jprb,0.64474_jprb,0.30637_jprb,0.031107_jprb, &
         & 1.0405_jprb,0.17538_jprb,1.3215_jprb,0.61819_jprb,0.31298_jprb,1.7218_jprb, &
         & 0.87761_jprb,0.47583_jprb,1.2583_jprb,1.0959_jprb,0.92842_jprb,0.51033_jprb/
    data a0_k2/0.26726_jprb,0.32033_jprb,-0.14778_jprb,0.28547_jprb,0.13592_jprb,0.13193_jprb, &
         & 0.37178_jprb,0.41813_jprb,0.33875_jprb,0.30203_jprb,0.35479_jprb,0.20189_jprb, &
         & 0.40663_jprb,0.47493_jprb,0.40668_jprb,0.14811_jprb,0.52382_jprb,0.86634_jprb, &
         & 0.14286_jprb,0.27164_jprb,-0.37947_jprb,0.2737_jprb,0.12001_jprb,-0.67315_jprb, &
         & 0.13492_jprb,0.065463_jprb,-0.19316_jprb,0.24905_jprb,0.25475_jprb,0.34637_jprb/
    data a0_eveh/0.9592599869E+00_jprb,0.9565299749E+00_jprb,0.9511899948E+00_jprb, &
         & 0.9560700059E+00_jprb,0.9541199803E+00_jprb,0.9483199716E+00_jprb, &
         & 0.9461100101E+00_jprb,0.9439799786E+00_jprb,0.9387800097E+00_jprb, &
         & 0.9317600131E+00_jprb,0.9289000034E+00_jprb,0.9236800075E+00_jprb, &
         & 0.9208700061E+00_jprb,0.9190599918E+00_jprb,0.9105200171E+00_jprb, &
         & 0.9162799716E+00_jprb,0.8937299848E+00_jprb,0.8014699817E+00_jprb, &
         & 0.9570500255E+00_jprb,0.9213600159E+00_jprb,0.7893999815E+00_jprb, &
         & 0.9639400244E+00_jprb,0.9530599713E+00_jprb,0.8850200176E+00_jprb, &
         & 0.9685299993E+00_jprb,0.9622600079E+00_jprb,0.9118800163E+00_jprb, &
         & 0.8997200131E+00_jprb,0.9012699723E+00_jprb,0.9107499719E+00_jprb/
    data a1_eveh/0.3627802414E-07_jprb,-0.7778328204E-08_jprb,0.4396108011E-07_jprb, &
         & 0.2503205394E-06_jprb,0.1996262995E-06_jprb,0.2929977541E-06_jprb, &
         & 0.4190530660E-06_jprb,0.3655744649E-06_jprb,0.3519195673E-06_jprb, &
         & 0.5574374313E-06_jprb,0.5273076340E-06_jprb,0.5376484182E-06_jprb, &
         & 0.1026844529E-05_jprb,0.9679998811E-06_jprb,0.8616486866E-06_jprb, &
         & 0.3180800832E-06_jprb,0.2886778532E-06_jprb,0.2310362675E-06_jprb, &
         & -0.1118036366E-06_jprb,-0.1502856577E-06_jprb,0.4842232926E-07_jprb, &
         & -0.8410978580E-08_jprb,-0.3478669441E-07_jprb,0.2209441590E-06_jprb, &
         & 0.2485776633E-06_jprb,0.1800235907E-06_jprb,0.2510202251E-06_jprb, &
         & 0.2687000915E-06_jprb,0.1740325644E-06_jprb,0.3562134339E-06_jprb/
    data a2_eveh/0.3067140824E-05_jprb,0.2520012231E-05_jprb,0.4831396382E-05_jprb, &
         & 0.8213598448E-05_jprb,0.7378375358E-05_jprb,0.1022081960E-04_jprb, &
         & 0.1225889173E-04_jprb,0.1165553113E-04_jprb,0.1188659007E-04_jprb, &
         & 0.1693615741E-04_jprb,0.1648317448E-04_jprb,0.1715818144E-04_jprb, &
         & 0.2744720041E-04_jprb,0.2642072104E-04_jprb,0.2671847506E-04_jprb, &
         & 0.1349592094E-04_jprb,0.1261523357E-04_jprb,0.5447756394E-05_jprb, &
         & 0.2064244654E-05_jprb,0.1919016057E-06_jprb,0.5940860319E-06_jprb, &
         & 0.5334760772E-05_jprb,0.4130339221E-05_jprb,0.4104662821E-05_jprb, &
         & 0.6530796327E-05_jprb,0.5727014013E-05_jprb,0.7451782039E-05_jprb, &
         & 0.1071246970E-04_jprb,0.9539280654E-05_jprb,0.1034286015E-04_jprb/
    data a3_eveh/-0.2004991551E-07_jprb,-0.6895366056E-07_jprb, &
         & -0.2047409282E-06_jprb, &
         & -0.7322448425E-07_jprb,-0.1273002681E-06_jprb,-0.2729916844E-06_jprb, &
         & -0.9421125213E-07_jprb,-0.1683332300E-06_jprb,-0.2726891637E-06_jprb, &
         & -0.1317753799E-06_jprb,-0.2107972250E-06_jprb,-0.3556060904E-06_jprb, &
         & -0.1889465580E-06_jprb,-0.2757958271E-06_jprb,-0.4909850304E-06_jprb, &
         & 0.7339644004E-08_jprb,-0.4058669560E-06_jprb,-0.4146343997E-06_jprb, &
         & 0.6170279931E-07_jprb,-0.1998567996E-06_jprb,-0.4713119139E-07_jprb, &
         & -0.1361754887E-07_jprb,-0.1765622955E-06_jprb,-0.2348146637E-06_jprb, &
         & -0.3901189061E-07_jprb,-0.1305666189E-06_jprb,-0.1533838798E-06_jprb, &
         & -0.2679148992E-07_jprb,-0.4441960044E-07_jprb,-0.1815613899E-06_jprb/
    data b0_eveh/0.9592599869E+00_jprb,0.9565299749E+00_jprb,0.9511899948E+00_jprb, &
         & 0.9560700059E+00_jprb,0.9541199803E+00_jprb,0.9483199716E+00_jprb, &
         & 0.9461100101E+00_jprb,0.9439799786E+00_jprb,0.9387800097E+00_jprb, &
         & 0.9317600131E+00_jprb,0.9289000034E+00_jprb,0.9236800075E+00_jprb, &
         & 0.9208700061E+00_jprb,0.9190599918E+00_jprb,0.9105200171E+00_jprb, &
         & 0.9162799716E+00_jprb,0.8937299848E+00_jprb,0.8014699817E+00_jprb, &
         & 0.9570500255E+00_jprb,0.9213600159E+00_jprb,0.7893999815E+00_jprb, &
         & 0.9639400244E+00_jprb,0.9530599713E+00_jprb,0.8850200176E+00_jprb, &
         & 0.9685299993E+00_jprb,0.9622600079E+00_jprb,0.9118800163E+00_jprb, &
         & 0.8997200131E+00_jprb,0.9012699723E+00_jprb,0.9107499719E+00_jprb/
    data b1_eveh/0.3626608347E-07_jprb,-0.7786279177E-08_jprb,0.4393379172E-07_jprb, &
         & 0.2502746099E-06_jprb,0.1995944388E-06_jprb,0.2929554341E-06_jprb, &
         & 0.4189516289E-06_jprb,0.3655020180E-06_jprb,0.3518483140E-06_jprb, &
         & 0.5572838404E-06_jprb,0.5271903092E-06_jprb,0.5375342766E-06_jprb, &
         & 0.1026605219E-05_jprb,0.9677979733E-06_jprb,0.8614680951E-06_jprb, &
         & 0.3179358714E-06_jprb,0.2884899004E-06_jprb,0.2308632219E-06_jprb, &
         & -0.1118781370E-06_jprb,-0.1503948681E-06_jprb,0.4834672396E-07_jprb, &
         & -0.8455684153E-08_jprb,-0.3485171618E-07_jprb,0.2208606134E-06_jprb, &
         & 0.2485595019E-06_jprb,0.1799959364E-06_jprb,0.2509846695E-06_jprb, &
         & 0.2686167306E-06_jprb,0.1739760478E-06_jprb,0.3561317214E-06_jprb/
    data b2_eveh/0.3065537157E-05_jprb,0.2518960400E-05_jprb,0.4829731552E-05_jprb, &
         & 0.8209894986E-05_jprb,0.7375769655E-05_jprb,0.1021809931E-04_jprb, &
         & 0.1225203869E-04_jprb,0.1165053800E-04_jprb,0.1188218721E-04_jprb, &
         & 0.1692612022E-04_jprb,0.1647546378E-04_jprb,0.1715117833E-04_jprb, &
         & 0.2743142431E-04_jprb,0.2640772436E-04_jprb,0.2670711910E-04_jprb, &
         & 0.1348545720E-04_jprb,0.1260529825E-04_jprb,0.5439695997E-05_jprb, &
         & 0.2058213340E-05_jprb,0.1860650656E-06_jprb,0.5898303925E-06_jprb, &
         & 0.5330772183E-05_jprb,0.4126528893E-05_jprb,0.4100859314E-05_jprb, &
         & 0.6528573977E-05_jprb,0.5725009032E-05_jprb,0.7449450095E-05_jprb, &
         & 0.1070590315E-04_jprb,0.9534271157E-05_jprb,0.1033751869E-04_jprb/
    data b3_eveh/-0.1370247134E-06_jprb,-0.1436897747E-06_jprb, &
         & -0.2954870411E-06_jprb, &
         & -0.3118435643E-06_jprb,-0.2916583242E-06_jprb,-0.4311032171E-06_jprb, &
         & -0.5048401022E-06_jprb,-0.4662823869E-06_jprb,-0.5206445053E-06_jprb, &
         & -0.7210980471E-06_jprb,-0.6662896794E-06_jprb,-0.7548637200E-06_jprb, &
         & -0.1110204039E-05_jprb,-0.1030801400E-05_jprb,-0.1140921199E-05_jprb, &
         & -0.6330818110E-06_jprb,-0.9186441048E-06_jprb,-0.7947813856E-06_jprb, &
         & -0.3242539890E-06_jprb,-0.5027602583E-06_jprb,-0.2777987334E-06_jprb, &
         & -0.2747250676E-06_jprb,-0.3811997260E-06_jprb,-0.4102405455E-06_jprb, &
         & -0.1994112324E-06_jprb,-0.2555484855E-06_jprb,-0.2842682534E-06_jprb, &
         & -0.4413041665E-06_jprb,-0.3717419474E-06_jprb,-0.4975536854E-06_jprb/
    ! Interpolation en angle
    DO j = 1,3
       ! Calcul par regression multilineaire de la valeur e0 en theta=0°
       e0 = a0_k0(j,class1)+a0_k1(j,class1)*ev(j)+a0_k2(j,class1)*eh(j)
       ! Lecture des coefficients des polynomes ev et eh
       a0 = a0_eveh(j,class1)
       a1 = a1_eveh(j,class1)
       a2 = a2_eveh(j,class1)
       a3 = a3_eveh(j,class1)
       b0 = b0_eveh(j,class1)
       b1 = b1_eveh(j,class1)
       b2 = b2_eveh(j,class1)
       b3 = b3_eveh(j,class1)
       theta0 = 0.
       theta53 = 53.
       ! Polarisation verticale
       S1_v = ((theta-theta53)/(theta0-theta53)) * ((e0-a0)/a0)
       em53_v = a3*(theta53**3) + a2*(theta53**2) + a1*theta53 + a0
       S2_v =((theta-theta0)/(theta53-theta0))*((ev(j)-em53_v)/em53_v)
       S_v = 1 + S1_v + S2_v
       emtheta_v = a3*(theta**3) + a2*(theta**2) + a1*theta + a0
       emiss_scal_v(j) = S_v * emtheta_v
       ! Polarisation horizontale
       S1_h = ((theta-theta53)/(theta0-theta53)) * ((e0-b0)/b0)
       em53_h = b3*(theta53**3) + b2*(theta53**2) + b1*theta53 + b0
       S2_h =((theta-theta0)/(theta53-theta0))*((eh(j)-em53_h)/em53_h)
       S_h = 1 + S1_h + S2_h
       emtheta_h = b3*(theta**3) + b2*(theta**2) + b1*theta + b0
       emiss_scal_h(j) = S_h * emtheta_h
    END DO
    ! Interpolation en frequence
    !emiss_interp_v=interp_freq(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq)
    !emiss_interp_h=interp_freq(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq)
    CALL interp_freq2(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq,class2,emiss_interp_v)
    CALL interp_freq2(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq,class2,emiss_interp_h)
    ! Cas ev<eh: on fait la moyenne entre les deux
    IF (emiss_interp_v < emiss_interp_h) THEN
       emiss_interp_v = (emiss_interp_v + emiss_interp_h)/2.
       emiss_interp_h =  emiss_interp_v
    END IF
    emiss_interp_h = MIN(1._jprb, emiss_interp_h)
    emiss_interp_v = MIN(1._jprb, emiss_interp_v)
  END SUBROUTINE emis_interp
END MODULE TELSEM2_ATLAS_READER
