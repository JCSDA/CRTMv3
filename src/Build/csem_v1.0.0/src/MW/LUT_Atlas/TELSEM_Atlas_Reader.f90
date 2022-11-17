MODULE TELSEM_ATLAS_READER
  ! Description:
  !   Data and routines for MW emissivity atlas.
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
  !  1.0      02/06/2010  Created (C. Prigent, E. Joumouille, F. Aires, J. Hocking)
  !  2.0      03/24/2011  Modified to use CRTM modules (Y. Chen)
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

  IMPLICIT none
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  ! The shared data
  PUBLIC :: atlas
  PUBLIC :: load_telsem_atlas
  PUBLIC :: emis_interp_ind_mult 
  PUBLIC :: emis_interp_int_mult
  PUBLIC :: emis_interp_ind_sing 
  PUBLIC :: emis_interp_int_sing
  PUBLIC :: close_telsem_atlas 
  PUBLIC :: telsem_atlas_version
  PUBLIC :: telsem_atlas
  
  TYPE telsem_atlas !  Atlas Data Structure
                                        
    INTEGER :: ndat                    ! number of lines in the atlas
    INTEGER :: nchan=7                 ! number of channels in the atlas
    CHARACTER(len=22) :: name          ! name of the atlas (including version number)
    INTEGER :: month                   ! month of the atlas            
    REAL(fp) :: dlat=0.25              ! resolution of the atlas (equal-area)
    INTEGER, POINTER :: ncells(:)      ! number of cells per lat band
    INTEGER, POINTER :: firstcell(:)   ! the first cellnumber of lat band
    REAL(fp) :: lat1, lat2, lon1, lon2 ! limits of the spatial domain (flaged if global)
    REAL(fp), POINTER :: emis(:,:)     ! Emissivities(ndat,nchan)
    REAL(fp), POINTER :: correl(:,:,:) ! correl(10,nchan,nchan)
    REAL(fp), POINTER :: emis_err(:,:) ! emis_err(ndat,nchan),Emissivity uncertainties (std)
    INTEGER, POINTER :: class(:)       ! Surface classe
    INTEGER, POINTER :: cellnum(:)     ! cellnumber of each of the pixels in the atlas
    INTEGER :: correspondance(660066)  ! "Correspondance" vector indicating that for the ith element, 
                                       ! the j so that EMIS(j,...) is the emissivity of cellnumber i.
  END TYPE telsem_atlas

  TYPE(telsem_atlas), SAVE :: atlas
  
  ! User can specify atlas version at set-up. Version 100 is the default.
  INTEGER :: telsem_atlas_version=100   ! Version of atlas
  
CONTAINS



  !------------------------------------------------------------------
  FUNCTION load_telsem_atlas(dir,month,lat1,lat2,lon1,lon2) RESULT( Error_Status )
    !======Read a monthly atlas
    INTEGER, INTENT (IN)       :: month
    CHARACTER(len=*), INTENT (IN)         :: dir
    REAL(fp),OPTIONAL, INTENT (IN) :: lat1, lat2, lon1, lon2 
    INTEGER    :: error_status 
    !----TRANSITORY VARIABLES
    INTEGER   :: ipos
    INTEGER   :: j,i
    INTEGER   :: iiin=21    ! unit for input
    INTEGER   :: ios
    CHARACTER :: month2*2
    INTEGER   :: cellnum
    REAL(fp) :: ssmi(2*7)
    REAL(fp) :: lat,lon
    INTEGER   :: cur_class
    INTEGER   :: take ! flag to take or not the pixel atlas for location constrains
    !initialisation
    error_status = SUCCESS
    !
    WRITE(month2,'(I2.2)') month
    atlas%month=month
    atlas%lat1=-777
    atlas%lat2=-777
    atlas%lon1=-777
    atlas%lon2=-777
    IF (PRESENT(lat1)) atlas%lat1=lat1
    IF (PRESENT(lat2)) atlas%lat2=lat2
    IF (PRESENT(lon1)) atlas%lon1=modulo(lon1,360.0_fp)
    IF (PRESENT(lon2)) atlas%lon2=modulo(lon2,360.0_fp)
    !----ALLOCATION SPECIFIC TO SSMI ATLAS
    ALLOCATE(atlas%ncells(int(180./atlas%dlat)))
    ALLOCATE(atlas%firstcell(int(180./atlas%dlat)) )
    CALL equare(atlas%dlat,atlas%ncells,atlas%firstcell)
    atlas%nchan=7
    atlas%name='ssmi_mean_emis_climato'
    !----INITIALISATIONS
    DO j=1,660066
       atlas%correspondance(j)=-777
    END DO
    !----DEFINING NUMBER OF DATA
    !WRITE(0,*) 'Reading number of data in atlas...'
    OPEN(iiin,file=dir//atlas%name//'_'//month2//'_cov_interpol',status='old',form='formatted',iostat=ios)
    IF( ios /= 0 ) THEN
       WRITE(0,*) 'Error opening the input file ios= ',ios
       !STOP
       error_status = FAILURE
    ENDIF
    READ(iiin,*,IOSTAT=ios) j
    atlas%ndat=j
    !WRITE(0,*) 'Nb data=',atlas%ndat

    !----ALLOCATE VARIABLES
    ALLOCATE(atlas%emis(atlas%ndat,atlas%nchan)) 
    ALLOCATE(atlas%emis_err(atlas%ndat,atlas%nchan)) 
    ALLOCATE(atlas%class(atlas%ndat)) 
    ALLOCATE(atlas%cellnum(atlas%ndat)) 
    lat = 45.0_fp; lon = 278.0_fp
    ipos=0
    DO j=1,atlas%ndat
       READ(iiin,*) cellnum, (ssmi(i),i=1,2*atlas%nchan),cur_class

       take=1
       IF (PRESENT(lat1)) THEN
          CALL get_coordinates(cellnum,lat,lon)
          IF ((lat<lat1).OR.(lat>lat2).OR.(lon<lon1).OR.(lon>lon2)) THEN
             take=0
          END IF
       END IF

       IF ((cur_class > 0).AND.(ipos<atlas%ndat).AND.(take==1)) THEN
          ipos=ipos+1
          DO i=1,atlas%nchan
             atlas%emis(ipos,i)=ssmi(i)
             atlas%emis_err(ipos,i)=sqrt(ssmi(atlas%nchan+i))
          END DO
          atlas%cellnum(ipos)=cellnum 
          atlas%class(ipos)=cur_class 
          atlas%correspondance(cellnum)=ipos 
       END IF          
    END DO
    atlas%ndat=ipos;
    CLOSE(iiin)

    !Correlation of uncertainties
    ALLOCATE(atlas%correl(10,atlas%nchan,atlas%nchan)) 
    !WRITE(0,*) 'reading classes...'
    OPEN(iiin,file=dir//'correlations', status='old')
    DO i=1,10
       READ(iiin,*) 
       DO j=1,7
          READ(iiin,'(7F5.2)') atlas%correl(i,j,1),atlas%correl(i,j,2), &
               & atlas%correl(i,j,3),atlas%correl(i,j,4),atlas%correl(i,j,5),&
               & atlas%correl(i,j,6),atlas%correl(i,j,7)
       END DO
    END DO
    CLOSE(iiin)

  END FUNCTION load_telsem_atlas

  !------------------------------------------------------------------
  SUBROUTINE close_telsem_atlas 
    !======free a monthly atlas
    IF ( ASSOCIATED(atlas%emis)      ) DEALLOCATE(atlas%emis) 
    IF ( ASSOCIATED(atlas%emis_err)  ) DEALLOCATE(atlas%emis_err) 
    IF ( ASSOCIATED(atlas%class)     ) DEALLOCATE(atlas%class) 
    IF ( ASSOCIATED(atlas%cellnum)   ) DEALLOCATE(atlas%cellnum) 
    IF ( ASSOCIATED(atlas%ncells)    ) DEALLOCATE(atlas%ncells) 
    IF ( ASSOCIATED(atlas%firstcell) ) DEALLOCATE(atlas%firstcell)
    IF ( ASSOCIATED(atlas%correl)    ) DEALLOCATE(atlas%correl) 
  END SUBROUTINE close_telsem_atlas 

  !------------------------------------------------------------------
  SUBROUTINE equare(DLAT,NCELLS,FIRSTCELL)
    !======computes the number of cells in a lattitude band
    !======and the first cellnumber of a lattitude band
    IMPLICIT none
    ! EQUAL-AREA COMPUTATIONS                                           
    ! I  TOTCEL             TOTAL NUMBER OF E.A. BOXES                  
    !  NCELLS(720)        NUMBER OF E.A. BOXES PER LAT ZONE           
    !  TOCELL(1440,720)   BOX NUMBER OF E.A. BOX (1-659600)         
    REAL(fp), INTENT (IN) :: dlat
    INTEGER,INTENT (INOUT) :: NCELLS(:)
    INTEGER,INTENT (INOUT) :: FIRSTCELL(:)
    INTEGER :: maxlat,maxlon
    INTEGER :: TOTCEL
    INTEGER, ALLOCATABLE   :: TOCELL(:,:)
    INTEGER :: maxlt2,lat,icellr,lat1,lat2,numcel,numcls,lon
    DOUBLE PRECISION    :: rcells,rearth,pi,aezon,hezon,aecell,xlatb
    DOUBLE PRECISION    :: rlatb,rlate,xlate,htb,hte,htzone,rcelat,azone
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
    TOTCEL=0                                                         
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
  FUNCTION calc_cellnum(lat,lon)
    !======computes the cellnumber given the lat and long
    !======using the ncells included in the atlas
    IMPLICIT none
    REAL(fp), INTENT (IN) :: lat, lon
    INTEGER :: calc_cellnum
    INTEGER :: ilat,ilon, i
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=int((lat+90.)/atlas%dlat)+1
    ilon=int(lon/(360./atlas%ncells(ilat)))+1
    do i=1,ilat-1
       calc_cellnum=atlas%ncells(i)+calc_cellnum
    end do
    calc_cellnum=calc_cellnum+ilon
  END FUNCTION calc_cellnum

  !------------------------------------------------------------------
  SUBROUTINE get_coordinates(cellnum,lat,lon)
    !======computes lat and lon given the cellnumber
    IMPLICIT none
    INTEGER, INTENT (IN) :: cellnum
    REAL(fp), INTENT (OUT) :: lat, lon
    !INTEGER :: ilat,ilon
    INTEGER :: i
    REAL(fp) :: res_lat ! latitude resolution
    INTEGER :: index_lat_max,index_lat,index_lon

    res_lat = atlas%dlat 

    index_lat_max  = int(180/res_lat)
 
    IF (cellnum >= atlas%firstcell(index_lat_max)) THEN
       index_lat = index_lat_max
       lat =(index_lat - 0.5)*res_lat - 90
       index_lon = cellnum - atlas%firstcell(index_lat_max) 
       lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
    ELSE
       DO i=1,index_lat_max-1 
          IF ( (cellnum>=atlas%firstcell(i)) .AND. (cellnum<atlas%firstcell(i+1)) ) THEN
             index_lat = i
             lat = (index_lat - 0.5)*res_lat- 90
             ! cout << i << "  " << *lat <<endl;
             index_lon = cellnum - atlas%firstcell(i)
             lon = (index_lon - 0.5)*(360.0/atlas%ncells(index_lat))
          END IF
       END DO
    END IF
  END SUBROUTINE get_coordinates

 !------------------------------------------------------------------
  SUBROUTINE calc_cellnum_mult(lat,lon,resol,cell_num_mult,nb_cell)
    !======computes the cellnumbers given the lat, long and resol
    !======using the ncells and firstcells included in the atlas
    IMPLICIT none
    REAL(fp), INTENT (IN) :: lat, lon
    ! lat(-90 90)  lon(-180 180) of the pixel center
    ! resol = spatial resolution of the new grid
    INTEGER, POINTER :: cell_num_mult(:)
    INTEGER, INTENT (OUT) :: nb_cell
    REAL(fp), INTENT(IN) :: resol
    INTEGER :: ilat,ilon, i
    INTEGER :: calc_cellnum
    !--------new variables
    INTEGER :: nbcel,cell(200) 
    !                 200 is here the maximum number of cells inside the "big" pixel
    INTEGER :: i2lon,i3lon,i2lat,nbreslat,nbreslon,i4lon
    REAL(fp) :: maxlat
    ! search for the cellnum in the older file
    calc_cellnum=0
    ilat=INT((lat+90.)/atlas%dlat)+1
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
                   IF (mod(ABS((i3lon-.5_fp)*(360._fp/atlas%ncells(i2lat))-lon),360._fp)<=resol/2.) THEN 
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
  SUBROUTINE emis_interp_ind_sing(lat,lon,theta,freq,ev,eh,stdv,stdh,covvh,verb)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !          SING: singular channel
    Implicit none     
    REAL(fp), INTENT (IN)    :: lat, lon, theta, freq
    REAL(fp), INTENT (OUT)   :: ev, eh   !resultat de l'interpolation
    REAL(fp), OPTIONAL, INTENT (OUT) :: stdv, stdh, covvh
    INTEGER, INTENT (IN) :: verb
    INTEGER         :: ipos, i, j, i2, j2
    INTEGER         :: cellnum
    REAL(fp) :: ev_a(3),eh_a(3)      !emissivity in the atlas

    ! Need to calculate the H-/V-pol covariance as well
    REAL(fp) :: FIM(2,3*2)           !Frequency Interpolation Matrix
    REAL(fp) :: trans_std(3*2,2)     !transpose Frequency Interpolation Matrix
    REAL(fp) :: a,b,c                !frequency interpolation coefficients
    REAL(fp) :: stdd(3),d            !dummy variables
    REAL(fp) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(fp) :: std(2,2)             !calculated covariance matrix
    REAL(fp) :: new_FIM(3*2,2)       !transpose Frequency Interpolation Matrix
    
    !Initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    IF (PRESENT(covvh)) covvh=0

    cellnum=calc_cellnum(lat,modulo(lon,360.0_fp))
    ipos=atlas%correspondance(cellnum)
    IF (ipos>0) THEN
       ev_a(1)=atlas%emis(ipos,1)
       eh_a(1)=atlas%emis(ipos,2)
       ev_a(2)=atlas%emis(ipos,4)
       eh_a(2)=atlas%emis(ipos,5)
       ev_a(3)=atlas%emis(ipos,6)
       eh_a(3)=atlas%emis(ipos,7)
       CALL emis_interp(theta,freq,atlas%class(ipos),ev_a,eh_a,ev,eh)
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
                cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                    & (atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
             END DO
          END DO

          !compute the Frequency Linear Matrix
          CALL interp_freq2(stdd(1),stdd(2),stdd(3),freq,d,a,b,c)
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
  SUBROUTINE emis_interp_ind_mult(lat,lon,theta,freq,n_chan,ev,eh,std,verb,stype)
    !======interpolates emissivities 
    !          IND: individual cellnumber atlas
    !          MULT: multiple channels
    IMPLICIT none     
    INTEGER, INTENT (IN) :: n_chan
    REAL(fp), INTENT (IN)    :: lat, lon, theta, freq(:)
    REAL(fp), INTENT (OUT)   :: ev(:), eh(:)   !resultat de l'interpolation
    REAL(fp), OPTIONAL, INTENT (OUT) :: std(:,:)
    INTEGER, INTENT (IN) :: verb
    INTEGER ,OPTIONAL         :: stype

    INTEGER         :: ipos, i, j, i2, j2
    INTEGER         :: cellnum
    REAL(fp) :: ev_a(3),eh_a(3)      !emissivity in the atlas
    REAL(fp) :: stdv_a(3),stdh_a(3)  !std in the atlas
    REAL(fp) :: std2
    REAL(fp) :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL(fp) :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    REAL(fp) :: a,b,c                !frequency interpolation coefficients
    REAL(fp) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas

    REAL(fp) :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL(fp) :: correlation(2*n_chan,2*n_chan) !correlation

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
    cellnum=calc_cellnum(lat,modulo(lon,360.0_fp))
    ipos=atlas%correspondance(cellnum)
    
    
    IF(PRESENT(stype)) stype=-1

    !compute emissivities
    IF (ipos>0) THEN
       ev_a(1)=atlas%emis(ipos,1)
       eh_a(1)=atlas%emis(ipos,2)
       ev_a(2)=atlas%emis(ipos,4)
       eh_a(2)=atlas%emis(ipos,5)
       ev_a(3)=atlas%emis(ipos,6)
       eh_a(3)=atlas%emis(ipos,7)
       stdv_a(1)=atlas%emis_err(ipos,1)  !don't matter to put
       stdh_a(1)=atlas%emis_err(ipos,2)  !real data here, we 
       stdv_a(2)=atlas%emis_err(ipos,4)  !just need a, b, c  
       stdh_a(2)=atlas%emis_err(ipos,5)  !don't matter to put
       stdv_a(3)=atlas%emis_err(ipos,6)  !real data here, we 
       stdh_a(3)=atlas%emis_err(ipos,7)  !just need a, b, c  
       IF(PRESENT(stype)) stype=atlas%class(ipos)
       DO i=1,n_chan
          !verifier que j'ai ici la bonne interpretation des canaux !!!!
          CALL emis_interp(theta,freq(i),atlas%class(ipos),ev_a,eh_a,ev(i),eh(i))
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
                cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                    & (atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
             END DO
          END DO

          !compute the Frequency Linear Matrix
          DO i=1,n_chan
             CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),std2,a,b,c)
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
  SUBROUTINE emis_interp_int_sing(lat,lon,resol,theta,freq,ev,eh,stdv,stdh,covvh,verb)
    !======interpolates emissivities 
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          SING: singular channel
    IMPLICIT none     
    REAL(fp), INTENT (IN)            :: lat, lon, theta, freq
    REAL(fp), INTENT (OUT)           :: ev, eh
    REAL(fp), OPTIONAL, INTENT (OUT) :: stdv,stdh,covvh
    INTEGER, INTENT (IN)         :: verb
    REAL(fp), INTENT (IN)            :: resol

    INTEGER          :: ipos, i, j, i2, j2
    REAL(fp)  :: ev_a(3),eh_a(3)     !emissivities in the atlas
    INTEGER          :: ii
    INTEGER          :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL(fp)  :: ev_mean, eh_mean
    INTEGER          :: inumb

    ! Need to calculate the H-/V-pol covariance as well
    REAL(fp) :: FIM(2,3*2)           !Frequency Interpolation Matrix
    REAL(fp) :: trans_std(3*2,2)     !transpose Frequency Interpolation Matrix
    REAL(fp) :: a,b,c                !frequency interpolation coefficients
    REAL(fp) :: stdd(3),d            !dummy variables
    REAL(fp) :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(fp) :: std2(2,2)            !calculated covariance matrix
    REAL(fp) :: std_mean(2,2)        !calculated average covariance matrix
    REAL(fp) :: new_FIM(3*2,2)       !transpose Frequency Interpolation Matrix
    INTEGER :: idummy
    idummy = verb
    !initialisations
    ev=0
    eh=0
    IF (PRESENT(stdv)) stdv=0
    IF (PRESENT(stdh)) stdh=0
    IF (PRESENT(covvh)) covvh=0
    stdd(:) = 0
    !computes the list of cells that need to be integrated in the atlas
    CALL calc_cellnum_mult(lat,modulo(lon,360.0_fp),resol,cellnum_mult,nb_cell)
    ev_mean=0
    eh_mean=0
    std_mean(:,:)=0
    inumb=0
    DO ii=1,nb_cell
       ipos=atlas%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          ev_a(1)=atlas%emis(ipos,1)
          eh_a(1)=atlas%emis(ipos,2)
          ev_a(2)=atlas%emis(ipos,4)
          eh_a(2)=atlas%emis(ipos,5)
          ev_a(3)=atlas%emis(ipos,6)
          eh_a(3)=atlas%emis(ipos,7)
          ev=0
          eh=0
          CALL emis_interp(theta,freq,atlas%class(ipos),ev_a,eh_a,ev,eh)
          !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas%correspondance(cellnum_mult(ii)),' ',ev,' ',eh
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
                   cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                       & (atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
                END DO
             END DO

             !compute the Frequency Linear Matrix
             CALL interp_freq2(stdd(1),stdd(2),stdd(3),freq,d,a,b,c)
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
  SUBROUTINE emis_interp_int_mult(lat,lon,resol,theta,freq,n_chan,ev,eh,std,verb,stype)
    !======interpolates emissivities 
    !          INT: integrate atlas, i.e. multiple cellnumber atlas
    !          MULT: multiple channel
    IMPLICIT none     
    INTEGER, INTENT (IN)             :: n_chan
    REAL(fp), INTENT (IN)                :: lat, lon, theta, freq(:)
    REAL(fp), INTENT (OUT)               :: ev(:), eh(:)
    REAL(fp), OPTIONAL, INTENT (OUT)     :: std(:,:)
    INTEGER, INTENT (IN)             :: verb
    REAL(fp), INTENT(IN)                 :: resol
    INTEGER          :: ipos
    REAL(fp)  :: ev_a(3),eh_a(3)     !emissivities in the atlas
    REAL(fp)  :: stdv_a(3) !std emissivities in the atlas
    INTEGER          :: ii, i2, j2
    INTEGER          :: nb_cell
    INTEGER, POINTER :: cellnum_mult(:) => NULL()
    REAL(fp)  :: ev_mean(n_chan), eh_mean(n_chan)
    REAL(fp)  :: std_mean(2*n_chan,2*n_chan)
    REAL(fp)  :: std2(2*n_chan,2*n_chan)
    REAL(fp)  :: ev2, eh2
    INTEGER          :: inumb, i,j
    INTEGER ,OPTIONAL         :: stype
   
    REAL(fp)  :: FIM(2*n_chan,3*2)    !Frequency Interpolation Matrix
    REAL(fp)  :: a,b,c                !frequency interpolation coefficients
    REAL(fp)  :: cov(6,6)             !covariance matrix of emis uncertainties in the atlas
    REAL(fp)  :: new_FIM(3*2,2*n_chan)          !transpose Frequency Interpolation Matrix
    REAL(fp)  :: trans_std(3*2,2*n_chan)        !transpose Frequency Interpolation Matrix
    !REAL(fp)  :: correlation(2*n_chan,2*n_chan) !correlation
    INTEGER :: idummy 
    idummy = verb
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
    CALL calc_cellnum_mult(lat,modulo(lon,360.0_fp),resol,cellnum_mult,nb_cell)
    
    IF(PRESENT(stype)) stype=-1
 
    inumb=0
    DO ii=1,nb_cell
       ipos=atlas%correspondance(cellnum_mult(ii))
       IF (ipos>0) THEN
          inumb=inumb+1
          IF(PRESENT(stype)) stype=atlas%class(ipos)
 
          ev_a(1)=atlas%emis(ipos,1)
          eh_a(1)=atlas%emis(ipos,2)
          ev_a(2)=atlas%emis(ipos,4)
          eh_a(2)=atlas%emis(ipos,5)
          ev_a(3)=atlas%emis(ipos,6)
          eh_a(3)=atlas%emis(ipos,7)
          stdv_a(1)=atlas%emis_err(ipos,1)  !don't matter to put
          !stdh_a(1)=atlas%emis_err(ipos,2)  !real data here, we 
          stdv_a(2)=atlas%emis_err(ipos,4)  !just need a, b, c  
          !stdh_a(2)=atlas%emis_err(ipos,5)  !don't matter to put
          stdv_a(3)=atlas%emis_err(ipos,6)  !real data here, we 
          !stdh_a(3)=atlas%emis_err(ipos,7)  !just need a, b, c  
          DO i=1,n_chan
             !initialisations    
             ev2=0
             eh2=0
             CALL emis_interp(theta,freq(i),atlas%class(ipos),ev_a,eh_a,ev2,eh2)
             !WRITE(0,*) 'Cellnum(',ii,')=',cellnum_mult(ii),' ',atlas%correspondance(cellnum_mult(ii)),' ',ev2,' ',eh2
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
                   cov(i,j)=atlas%correl(atlas%class(ipos),i2,j2)* &
                       & (atlas%emis_err(ipos,i2)*atlas%emis_err(ipos,j2))
                END DO
             END DO

             !compute the Frequency Linear Matrix
             DO i=1,n_chan
                CALL interp_freq2(stdv_a(1),stdv_a(2),stdv_a(3),freq(i),ev2,a,b,c) !ev2=dummy variable
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
  SUBROUTINE interp_freq2(emiss19,emiss37,emiss85,f,emiss,an,bn,cn)  
    !======linear interpolation of emissivity given the freq
    !======and the atlas values on that cellnum
    IMPLICIT none
    REAL(fp), INTENT (IN) :: emiss19,emiss37,emiss85,f
    REAL(fp), INTENT (OUT) :: emiss
    REAL(fp), OPTIONAL, INTENT (OUT) :: an, bn, cn
    REAL(fp) :: a=0, b=0, c=0
    IF (f<=19.35) THEN
       a=1
       b=0
       c=0
       emiss = emiss19
    ELSE IF ((19.35<f).AND.(f<=37.)) THEN
       a=(37.-f  )/(37.-19.35)
       b=(f-19.35)/(37.-19.35)
       c=0
       emiss = a*emiss19+b*emiss37
    ELSE IF ((f>37.).AND.(f<85.5)) THEN
       a=0
       b=(85.5-f )/(85.5-37)
       c=(f-37   )/(85.5-37)
       emiss = b*emiss37+c*emiss85
    ELSE IF(85.5<=f) THEN
       a=0
       b=0
       c=1
       emiss = emiss85
    END IF
    IF (PRESENT(an)) an=a
    IF (PRESENT(bn)) bn=b
    IF (PRESENT(cn)) cn=c
  END SUBROUTINE interp_freq2

 !------------------------------------------------------------------
   SUBROUTINE emis_interp(theta,freq,classe,ev,eh,emiss_interp_v,emiss_interp_h)
    !======interpolation of emissivities for angle, freq
    IMPLICIT none     
    INTEGER, INTENT (IN)          :: classe
    REAL(fp), INTENT (IN)  :: theta, freq, ev(3), eh(3)
    REAL(fp), INTENT (OUT) :: emiss_interp_v, emiss_interp_h
    REAL(fp) :: e0, theta0, theta53 , emiss_scal_v(3), emiss_scal_h(3)
    REAL(fp) :: S1_v, S1_h, S2_v, S2_h, S_v, S_h, a0, a1, a2, a3, b0, b1
    REAL(fp) :: b2, b3, em53_v, em53_h, emtheta_v, emtheta_h
    REAL(fp) :: a0_k0(3,10),a0_k1(3,10),a0_k2(3,10)
    REAL(fp) :: a0_eveh(3,10),a1_eveh(3,10),a2_eveh(3,10),a3_eveh(3,10)
    REAL(fp) :: b0_eveh(3,10),b1_eveh(3,10),b2_eveh(3,10),b3_eveh(3,10)
    INTEGER         :: j
    !COMMON /EMISSIVITE/emiss_interp_v,emiss_interp_h
    data a0_k0/0.11509_fp,0.091535_fp,0.34796_fp,0.10525_fp,0.16627_fp,0.24434_fp, &
         & 0.29217_fp,0.23809_fp,0.28954_fp,0.17516_fp,0.19459_fp,0.28697_fp, &
         & 0.10521_fp,0.12126_fp,0.30278_fp,0.18212_fp,0.19625_fp,0.14551_fp, &
         & -0.19202_fp,0.5411_fp,0.03739_fp,0.10292_fp,0.5486_fp,-0.058937_fp, &
         & -0.022672_fp,0.44492_fp,-0.058448_fp,-0.33894_fp,-0.17621_fp,0.14742_fp/
    data a0_k1/0.61168_fp,0.59095_fp,0.7918_fp,0.60271_fp,0.69213_fp,0.62218_fp, &
         &  0.32728_fp,0.34334_fp,0.37062_fp,0.51217_fp,0.4491_fp,0.50101_fp, &
         & 0.48913_fp,0.41932_fp,0.29734_fp,0.64474_fp,0.30637_fp,0.031107_fp, &
         & 1.0405_fp,0.17538_fp,1.3215_fp,0.61819_fp,0.31298_fp,1.7218_fp, &
         & 0.87761_fp,0.47583_fp,1.2583_fp,1.0959_fp,0.92842_fp,0.51033_fp/
    data a0_k2/0.26726_fp,0.32033_fp,-0.14778_fp,0.28547_fp,0.13592_fp,0.13193_fp, &
         & 0.37178_fp,0.41813_fp,0.33875_fp,0.30203_fp,0.35479_fp,0.20189_fp, &
         & 0.40663_fp,0.47493_fp,0.40668_fp,0.14811_fp,0.52382_fp,0.86634_fp, &
         & 0.14286_fp,0.27164_fp,-0.37947_fp,0.2737_fp,0.12001_fp,-0.67315_fp, &
         & 0.13492_fp,0.065463_fp,-0.19316_fp,0.24905_fp,0.25475_fp,0.34637_fp/
    data a0_eveh/0.9592599869E+00_fp,0.9565299749E+00_fp,0.9511899948E+00_fp, &
         & 0.9560700059E+00_fp,0.9541199803E+00_fp,0.9483199716E+00_fp, &
         & 0.9461100101E+00_fp,0.9439799786E+00_fp,0.9387800097E+00_fp, &
         & 0.9317600131E+00_fp,0.9289000034E+00_fp,0.9236800075E+00_fp, &
         & 0.9208700061E+00_fp,0.9190599918E+00_fp,0.9105200171E+00_fp, &
         & 0.9162799716E+00_fp,0.8937299848E+00_fp,0.8014699817E+00_fp, &
         & 0.9570500255E+00_fp,0.9213600159E+00_fp,0.7893999815E+00_fp, &
         & 0.9639400244E+00_fp,0.9530599713E+00_fp,0.8850200176E+00_fp, &
         & 0.9685299993E+00_fp,0.9622600079E+00_fp,0.9118800163E+00_fp, &
         & 0.8997200131E+00_fp,0.9012699723E+00_fp,0.9107499719E+00_fp/
    data a1_eveh/0.3627802414E-07_fp,-0.7778328204E-08_fp,0.4396108011E-07_fp, &
         & 0.2503205394E-06_fp,0.1996262995E-06_fp,0.2929977541E-06_fp, &
         & 0.4190530660E-06_fp,0.3655744649E-06_fp,0.3519195673E-06_fp, &
         & 0.5574374313E-06_fp,0.5273076340E-06_fp,0.5376484182E-06_fp, &
         & 0.1026844529E-05_fp,0.9679998811E-06_fp,0.8616486866E-06_fp, &
         & 0.3180800832E-06_fp,0.2886778532E-06_fp,0.2310362675E-06_fp, &
         & -0.1118036366E-06_fp,-0.1502856577E-06_fp,0.4842232926E-07_fp, &
         & -0.8410978580E-08_fp,-0.3478669441E-07_fp,0.2209441590E-06_fp, &
         & 0.2485776633E-06_fp,0.1800235907E-06_fp,0.2510202251E-06_fp, &
         & 0.2687000915E-06_fp,0.1740325644E-06_fp,0.3562134339E-06_fp/
    data a2_eveh/0.3067140824E-05_fp,0.2520012231E-05_fp,0.4831396382E-05_fp, &
         & 0.8213598448E-05_fp,0.7378375358E-05_fp,0.1022081960E-04_fp, &
         & 0.1225889173E-04_fp,0.1165553113E-04_fp,0.1188659007E-04_fp, &
         & 0.1693615741E-04_fp,0.1648317448E-04_fp,0.1715818144E-04_fp, &
         & 0.2744720041E-04_fp,0.2642072104E-04_fp,0.2671847506E-04_fp, &
         & 0.1349592094E-04_fp,0.1261523357E-04_fp,0.5447756394E-05_fp, &
         & 0.2064244654E-05_fp,0.1919016057E-06_fp,0.5940860319E-06_fp, &
         & 0.5334760772E-05_fp,0.4130339221E-05_fp,0.4104662821E-05_fp, &
         & 0.6530796327E-05_fp,0.5727014013E-05_fp,0.7451782039E-05_fp, &
         & 0.1071246970E-04_fp,0.9539280654E-05_fp,0.1034286015E-04_fp/
    data a3_eveh/-0.2004991551E-07_fp,-0.6895366056E-07_fp, &
         & -0.2047409282E-06_fp, &
         & -0.7322448425E-07_fp,-0.1273002681E-06_fp,-0.2729916844E-06_fp, &
         & -0.9421125213E-07_fp,-0.1683332300E-06_fp,-0.2726891637E-06_fp, &
         & -0.1317753799E-06_fp,-0.2107972250E-06_fp,-0.3556060904E-06_fp, &
         & -0.1889465580E-06_fp,-0.2757958271E-06_fp,-0.4909850304E-06_fp, &
         & 0.7339644004E-08_fp,-0.4058669560E-06_fp,-0.4146343997E-06_fp, &
         & 0.6170279931E-07_fp,-0.1998567996E-06_fp,-0.4713119139E-07_fp, &
         & -0.1361754887E-07_fp,-0.1765622955E-06_fp,-0.2348146637E-06_fp, &
         & -0.3901189061E-07_fp,-0.1305666189E-06_fp,-0.1533838798E-06_fp, &
         & -0.2679148992E-07_fp,-0.4441960044E-07_fp,-0.1815613899E-06_fp/
    data b0_eveh/0.9592599869E+00_fp,0.9565299749E+00_fp,0.9511899948E+00_fp, &
         & 0.9560700059E+00_fp,0.9541199803E+00_fp,0.9483199716E+00_fp, &
         & 0.9461100101E+00_fp,0.9439799786E+00_fp,0.9387800097E+00_fp, &
         & 0.9317600131E+00_fp,0.9289000034E+00_fp,0.9236800075E+00_fp, &
         & 0.9208700061E+00_fp,0.9190599918E+00_fp,0.9105200171E+00_fp, &
         & 0.9162799716E+00_fp,0.8937299848E+00_fp,0.8014699817E+00_fp, &
         & 0.9570500255E+00_fp,0.9213600159E+00_fp,0.7893999815E+00_fp, &
         & 0.9639400244E+00_fp,0.9530599713E+00_fp,0.8850200176E+00_fp, &
         & 0.9685299993E+00_fp,0.9622600079E+00_fp,0.9118800163E+00_fp, &
         & 0.8997200131E+00_fp,0.9012699723E+00_fp,0.9107499719E+00_fp/
    data b1_eveh/0.3626608347E-07_fp,-0.7786279177E-08_fp,0.4393379172E-07_fp, &
         & 0.2502746099E-06_fp,0.1995944388E-06_fp,0.2929554341E-06_fp, &
         & 0.4189516289E-06_fp,0.3655020180E-06_fp,0.3518483140E-06_fp, &
         & 0.5572838404E-06_fp,0.5271903092E-06_fp,0.5375342766E-06_fp, &
         & 0.1026605219E-05_fp,0.9677979733E-06_fp,0.8614680951E-06_fp, &
         & 0.3179358714E-06_fp,0.2884899004E-06_fp,0.2308632219E-06_fp, &
         & -0.1118781370E-06_fp,-0.1503948681E-06_fp,0.4834672396E-07_fp, &
         & -0.8455684153E-08_fp,-0.3485171618E-07_fp,0.2208606134E-06_fp, &
         & 0.2485595019E-06_fp,0.1799959364E-06_fp,0.2509846695E-06_fp, &
         & 0.2686167306E-06_fp,0.1739760478E-06_fp,0.3561317214E-06_fp/
    data b2_eveh/0.3065537157E-05_fp,0.2518960400E-05_fp,0.4829731552E-05_fp, &
         & 0.8209894986E-05_fp,0.7375769655E-05_fp,0.1021809931E-04_fp, &
         & 0.1225203869E-04_fp,0.1165053800E-04_fp,0.1188218721E-04_fp, &
         & 0.1692612022E-04_fp,0.1647546378E-04_fp,0.1715117833E-04_fp, &
         & 0.2743142431E-04_fp,0.2640772436E-04_fp,0.2670711910E-04_fp, &
         & 0.1348545720E-04_fp,0.1260529825E-04_fp,0.5439695997E-05_fp, &
         & 0.2058213340E-05_fp,0.1860650656E-06_fp,0.5898303925E-06_fp, &
         & 0.5330772183E-05_fp,0.4126528893E-05_fp,0.4100859314E-05_fp, &
         & 0.6528573977E-05_fp,0.5725009032E-05_fp,0.7449450095E-05_fp, &
         & 0.1070590315E-04_fp,0.9534271157E-05_fp,0.1033751869E-04_fp/
    data b3_eveh/-0.1370247134E-06_fp,-0.1436897747E-06_fp, &
         & -0.2954870411E-06_fp, &
         & -0.3118435643E-06_fp,-0.2916583242E-06_fp,-0.4311032171E-06_fp, &
         & -0.5048401022E-06_fp,-0.4662823869E-06_fp,-0.5206445053E-06_fp, &
         & -0.7210980471E-06_fp,-0.6662896794E-06_fp,-0.7548637200E-06_fp, &
         & -0.1110204039E-05_fp,-0.1030801400E-05_fp,-0.1140921199E-05_fp, &
         & -0.6330818110E-06_fp,-0.9186441048E-06_fp,-0.7947813856E-06_fp, &
         & -0.3242539890E-06_fp,-0.5027602583E-06_fp,-0.2777987334E-06_fp, &
         & -0.2747250676E-06_fp,-0.3811997260E-06_fp,-0.4102405455E-06_fp, &
         & -0.1994112324E-06_fp,-0.2555484855E-06_fp,-0.2842682534E-06_fp, &
         & -0.4413041665E-06_fp,-0.3717419474E-06_fp,-0.4975536854E-06_fp/
    ! Interpolation en angle
    DO j = 1,3
       ! Calcul par regression multilineaire de la valeur e0 en theta=0°
       e0 = a0_k0(j,classe)+a0_k1(j,classe)*ev(j)+a0_k2(j,classe)*eh(j)
       ! Lecture des coefficients des polynomes ev et eh
       a0 = a0_eveh(j,classe)
       a1 = a1_eveh(j,classe)
       a2 = a2_eveh(j,classe)
       a3 = a3_eveh(j,classe)
       b0 = b0_eveh(j,classe)
       b1 = b1_eveh(j,classe)
       b2 = b2_eveh(j,classe)
       b3 = b3_eveh(j,classe)
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
    CALL interp_freq2(emiss_scal_v(1),emiss_scal_v(2),emiss_scal_v(3),freq,emiss_interp_v)
    CALL interp_freq2(emiss_scal_h(1),emiss_scal_h(2),emiss_scal_h(3),freq,emiss_interp_h)  
    ! Cas ev<eh: on fait la moyenne entre les deux
    IF (emiss_interp_v < emiss_interp_h) THEN
       emiss_interp_v = (emiss_interp_v + emiss_interp_h)/2.
       emiss_interp_h =  emiss_interp_v
    END IF
  END SUBROUTINE emis_interp


END MODULE TELSEM_ATLAS_READER
