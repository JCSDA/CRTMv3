!
! TELSEM2Atlas_netCDF_IO
!
! Module for reading and writing the CRTM-native netCDF TELSEM2 atlas coefficient
! file. The file holds all twelve monthly atlases stacked along the data
! dimension with per-month offsets. Only the emissivity-bearing data are stored;
! the equal-area grid geometry and reverse lookup are reconstructed on load.
!

MODULE TELSEM2Atlas_netCDF_IO

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds         , ONLY: fp, Double, Long
  USE Message_Handler    , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility       , ONLY: File_Exists
  USE TELSEM2Atlas_Define, ONLY: TELSEM2Atlas_type     , &
                                 TELSEM2Atlas_Associated, &
                                 TELSEM2Atlas_Create    , &
                                 TELSEM2Atlas_Destroy
  USE netcdf
  ! Disable implicit typing
  IMPLICIT NONE


  ! ------------
  ! Visibilities
  ! ------------
  PRIVATE
  PUBLIC :: TELSEM2Atlas_netCDF_InquireFile
  PUBLIC :: TELSEM2Atlas_netCDF_ReadFile
  PUBLIC :: TELSEM2Atlas_netCDF_WriteFile


  ! -----------------
  ! Module parameters
  ! -----------------
  INTEGER, PARAMETER :: ML = 256
  CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'TELSEM2Atlas_netCDF_IO'
  ! Number of angular-correction classes (second dim of the emis_interp
  ! regression tables in TELSEM2_Atlas_Module); the valid Class1 range is 1..N.
  INTEGER, PARAMETER :: N_CLASS1 = 10
  ! Dimension names
  CHARACTER(*), PARAMETER :: CHANNEL_DIMNAME = 'n_channels'
  CHARACTER(*), PARAMETER :: BAND_DIMNAME    = 'n_latitude_bands'
  CHARACTER(*), PARAMETER :: MONTH_DIMNAME   = 'n_months'
  CHARACTER(*), PARAMETER :: DATA_DIMNAME    = 'n_data'
  ! Variable names
  CHARACTER(*), PARAMETER :: MONTHCOUNT_VARNAME = 'month_data_count'
  CHARACTER(*), PARAMETER :: MONTHOFF_VARNAME   = 'month_offset'
  CHARACTER(*), PARAMETER :: CELLNUM_VARNAME    = 'cell_number'
  CHARACTER(*), PARAMETER :: CLASS1_VARNAME     = 'class1'
  CHARACTER(*), PARAMETER :: CLASS2_VARNAME     = 'class2'
  CHARACTER(*), PARAMETER :: EMIS_VARNAME       = 'emissivity'
  ! Global attribute names
  CHARACTER(*), PARAMETER :: RELEASE_GATTNAME    = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_GATTNAME    = 'Version'
  CHARACTER(*), PARAMETER :: RESOLUTION_GATTNAME = 'Resolution'
  CHARACTER(*), PARAMETER :: NCELLS_GATTNAME     = 'n_Cells'


CONTAINS


!--------------------------------------------------------------------------------
! TELSEM2Atlas_netCDF_InquireFile
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_netCDF_InquireFile( &
    Filename        , &
    n_Channels      , &
    n_Latitude_Bands, &
    n_Cells         , &
    n_Months        , &
    n_Data          , &
    Release         , &
    Version         , &
    Resolution      ) RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN)  :: Filename
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: n_Latitude_Bands
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: n_Cells
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: n_Months
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: n_Data
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: Release
    INTEGER(Long), OPTIONAL, INTENT(OUT) :: Version
    REAL(Double),  OPTIONAL, INTENT(OUT) :: Resolution
    ! Function result
    INTEGER :: err_stat
    ! Local variables
    CHARACTER(ML) :: msg
    INTEGER :: ncid
    INTEGER :: n

    err_stat = SUCCESS
    IF ( .NOT. File_Exists(Filename) ) THEN
      err_stat = FAILURE
      CALL Display_Message( ROUTINE_NAME, 'File '//TRIM(Filename)//' not found.', err_stat )
      RETURN
    END IF

    IF ( .NOT. check(NF90_OPEN(Filename, NF90_NOWRITE, ncid), 'open '//TRIM(Filename)) ) THEN
      err_stat = FAILURE; RETURN
    END IF

    IF ( PRESENT(n_Channels) ) THEN
      n = get_dim(ncid, CHANNEL_DIMNAME); IF (n<0) GOTO 900; n_Channels = n
    END IF
    IF ( PRESENT(n_Latitude_Bands) ) THEN
      n = get_dim(ncid, BAND_DIMNAME); IF (n<0) GOTO 900; n_Latitude_Bands = n
    END IF
    IF ( PRESENT(n_Months) ) THEN
      n = get_dim(ncid, MONTH_DIMNAME); IF (n<0) GOTO 900; n_Months = n
    END IF
    IF ( PRESENT(n_Data) ) THEN
      n = get_dim(ncid, DATA_DIMNAME); IF (n<0) GOTO 900; n_Data = n
    END IF
    IF ( PRESENT(Release) ) THEN
      IF (.NOT. check(NF90_GET_ATT(ncid,NF90_GLOBAL,RELEASE_GATTNAME,Release),'get Release')) GOTO 900
    END IF
    IF ( PRESENT(Version) ) THEN
      IF (.NOT. check(NF90_GET_ATT(ncid,NF90_GLOBAL,VERSION_GATTNAME,Version),'get Version')) GOTO 900
    END IF
    IF ( PRESENT(n_Cells) ) THEN
      IF (.NOT. check(NF90_GET_ATT(ncid,NF90_GLOBAL,NCELLS_GATTNAME,n_Cells),'get n_Cells')) GOTO 900
    END IF
    IF ( PRESENT(Resolution) ) THEN
      IF (.NOT. check(NF90_GET_ATT(ncid,NF90_GLOBAL,RESOLUTION_GATTNAME,Resolution),'get Resolution')) GOTO 900
    END IF

    IF ( .NOT. check(NF90_CLOSE(ncid), 'close') ) err_stat = FAILURE
    RETURN
900 CONTINUE
    err_stat = FAILURE
    msg = NF90_STRERROR(NF90_CLOSE(ncid))
  END FUNCTION TELSEM2Atlas_netCDF_InquireFile


!--------------------------------------------------------------------------------
! TELSEM2Atlas_netCDF_WriteFile
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_netCDF_WriteFile( Filename, Atlas ) RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN) :: Filename
    TYPE(TELSEM2Atlas_type), INTENT(IN) :: Atlas
    ! Function result
    INTEGER :: err_stat
    ! Local variables
    INTEGER :: ncid
    INTEGER :: dimid_chan, dimid_band, dimid_month, dimid_data
    INTEGER :: vid_mcount, vid_moff, vid_cell, vid_c1, vid_c2, vid_emis

    err_stat = FAILURE
    IF ( .NOT. TELSEM2Atlas_Associated(Atlas) ) THEN
      CALL Display_Message( ROUTINE_NAME, 'Atlas not allocated', FAILURE ); RETURN
    END IF

    IF (.NOT. check(NF90_CREATE(Filename, IOR(NF90_NETCDF4,NF90_CLOBBER), ncid),'create '//TRIM(Filename))) RETURN

    ! Global attributes
    IF (.NOT. check(NF90_PUT_ATT(ncid,NF90_GLOBAL,RELEASE_GATTNAME,Atlas%Release),'put Release')) GOTO 900
    IF (.NOT. check(NF90_PUT_ATT(ncid,NF90_GLOBAL,VERSION_GATTNAME,Atlas%Version),'put Version')) GOTO 900
    IF (.NOT. check(NF90_PUT_ATT(ncid,NF90_GLOBAL,RESOLUTION_GATTNAME,Atlas%Resolution),'put Resolution')) GOTO 900
    IF (.NOT. check(NF90_PUT_ATT(ncid,NF90_GLOBAL,NCELLS_GATTNAME,Atlas%n_Cells),'put n_Cells')) GOTO 900

    ! Dimensions
    IF (.NOT. check(NF90_DEF_DIM(ncid,CHANNEL_DIMNAME,INT(Atlas%n_Channels),dimid_chan),'def chan')) GOTO 900
    IF (.NOT. check(NF90_DEF_DIM(ncid,BAND_DIMNAME,INT(Atlas%n_Latitude_Bands),dimid_band),'def band')) GOTO 900
    IF (.NOT. check(NF90_DEF_DIM(ncid,MONTH_DIMNAME,INT(Atlas%n_Months),dimid_month),'def month')) GOTO 900
    IF (.NOT. check(NF90_DEF_DIM(ncid,DATA_DIMNAME,INT(Atlas%n_Data),dimid_data),'def data')) GOTO 900

    ! Variable definitions
    IF (.NOT. check(NF90_DEF_VAR(ncid,MONTHCOUNT_VARNAME,NF90_INT,dimid_month,vid_mcount),'def mcount')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,MONTHOFF_VARNAME,NF90_INT,dimid_month,vid_moff),'def moff')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CELLNUM_VARNAME,NF90_INT,dimid_data,vid_cell),'def cell')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CLASS1_VARNAME,NF90_INT,dimid_data,vid_c1),'def c1')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CLASS2_VARNAME,NF90_INT,dimid_data,vid_c2),'def c2')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,EMIS_VARNAME,NF90_DOUBLE,[dimid_data,dimid_chan],vid_emis),'def emis')) GOTO 900

    IF (.NOT. check(NF90_ENDDEF(ncid),'enddef')) GOTO 900

    ! Write data
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_mcount,INT(Atlas%Month_Data_Count)),'put mcount')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_moff,INT(Atlas%Month_Offset)),'put moff')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_cell,INT(Atlas%Cell_Number)),'put cell')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_c1,INT(Atlas%Class1)),'put c1')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_c2,INT(Atlas%Class2)),'put c2')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_emis,Atlas%Emissivity),'put emis')) GOTO 900

    IF (.NOT. check(NF90_CLOSE(ncid),'close')) RETURN
    err_stat = SUCCESS
    RETURN
900 CONTINUE
    err_stat = FAILURE
    IF (check(NF90_CLOSE(ncid),'close')) CONTINUE
  END FUNCTION TELSEM2Atlas_netCDF_WriteFile


!--------------------------------------------------------------------------------
! TELSEM2Atlas_netCDF_ReadFile
!
! Reads the primary atlas data into a created TELSEM2Atlas object. The derived
! grid geometry / reverse lookup are NOT filled here; the caller must call
! TELSEM2_Setup_Grid after a successful read.
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_netCDF_ReadFile( Filename, Atlas ) RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN)  :: Filename
    TYPE(TELSEM2Atlas_type), INTENT(OUT) :: Atlas
    ! Function result
    INTEGER :: err_stat
    ! Local variables
    INTEGER :: ncid, vid
    INTEGER(Long) :: n_Channels, n_Bands, n_Cells, n_Months, n_Data, rel, ver
    REAL(Double)  :: resolution

    err_stat = FAILURE
    IF ( .NOT. File_Exists(Filename) ) THEN
      CALL Display_Message( ROUTINE_NAME, 'File '//TRIM(Filename)//' not found.', FAILURE ); RETURN
    END IF

    ! Pull dimensions and attributes
    err_stat = TELSEM2Atlas_netCDF_InquireFile( Filename, &
                 n_Channels=n_Channels, n_Latitude_Bands=n_Bands, n_Cells=n_Cells, &
                 n_Months=n_Months, n_Data=n_Data, Release=rel, Version=ver, &
                 Resolution=resolution )
    IF ( err_stat /= SUCCESS ) RETURN
    err_stat = FAILURE

    CALL TELSEM2Atlas_Create( Atlas, n_Channels, n_Bands, n_Cells, n_Months, n_Data )
    IF ( .NOT. TELSEM2Atlas_Associated(Atlas) ) THEN
      CALL Display_Message( ROUTINE_NAME, 'Atlas allocation failed', FAILURE ); RETURN
    END IF
    Atlas%Release    = rel
    Atlas%Version    = ver
    Atlas%Resolution = resolution

    IF (.NOT. check(NF90_OPEN(Filename, NF90_NOWRITE, ncid),'open')) THEN
      CALL TELSEM2Atlas_Destroy(Atlas); RETURN
    END IF

    ! Per-month indexing (read into default-int temporaries, then assign)
    IF (.NOT. read_int_vec(ncid,MONTHCOUNT_VARNAME,Atlas%Month_Data_Count)) GOTO 900
    IF (.NOT. read_int_vec(ncid,MONTHOFF_VARNAME,Atlas%Month_Offset)) GOTO 900
    IF (.NOT. read_int_vec(ncid,CELLNUM_VARNAME,Atlas%Cell_Number)) GOTO 900
    IF (.NOT. read_int_vec(ncid,CLASS1_VARNAME,Atlas%Class1)) GOTO 900
    IF (.NOT. read_int_vec(ncid,CLASS2_VARNAME,Atlas%Class2)) GOTO 900
    ! Class1 indexes the (3,N_CLASS1) angular-regression tables in emis_interp;
    ! the stacked arrays hold only populated cells, so every Class1 must be in
    ! range. Reject a corrupt/wrong file here rather than read out of bounds at
    ! interpolation time.
    IF ( ANY( Atlas%Class1 < 1 .OR. Atlas%Class1 > N_CLASS1 ) ) THEN
      CALL Display_Message( ROUTINE_NAME, &
        'class1 values out of range in '//TRIM(Filename), FAILURE )
      GOTO 900
    END IF
    ! Cell_Number indexes Correspondence(n_Cells,:) at grid-setup time -- an
    ! out-of-range value is an out-of-bounds WRITE during CRTM_Init.
    IF ( ANY( Atlas%Cell_Number < 1 .OR. Atlas%Cell_Number > n_Cells ) ) THEN
      CALL Display_Message( ROUTINE_NAME, &
        'cellnum values out of range in '//TRIM(Filename), FAILURE )
      GOTO 900
    END IF
    ! Per-month slices must partition the stacked data arrays.
    IF ( ANY( Atlas%Month_Offset < 0 ) .OR. &
         ANY( Atlas%Month_Offset + Atlas%Month_Data_Count > n_Data ) ) THEN
      CALL Display_Message( ROUTINE_NAME, &
        'month offset/count indexing exceeds n_data in '//TRIM(Filename), FAILURE )
      GOTO 900
    END IF
    ! The band count must agree with the resolution attribute: equare() loops
    ! to FLOOR(180/resolution) writing per-band arrays sized by n_Bands.
    IF ( n_Bands /= FLOOR( 180.0_Double/resolution ) ) THEN
      CALL Display_Message( ROUTINE_NAME, &
        'n_latitude_bands inconsistent with resolution attribute in '//TRIM(Filename), FAILURE )
      GOTO 900
    END IF
    ! Emissivity
    IF (.NOT. check(NF90_INQ_VARID(ncid,EMIS_VARNAME,vid),'inq emis')) GOTO 900
    IF (.NOT. check(NF90_GET_VAR(ncid,vid,Atlas%Emissivity),'get emis')) GOTO 900

    IF (.NOT. check(NF90_CLOSE(ncid),'close')) THEN
      CALL TELSEM2Atlas_Destroy(Atlas); RETURN
    END IF
    err_stat = SUCCESS
    RETURN
900 CONTINUE
    IF (check(NF90_CLOSE(ncid),'close')) CONTINUE
    CALL TELSEM2Atlas_Destroy(Atlas)
    err_stat = FAILURE
  END FUNCTION TELSEM2Atlas_netCDF_ReadFile


!################################################################################
!##                          ## PRIVATE PROCEDURES ##                          ##
!################################################################################

  ! Read an integer vector variable into an INTEGER(Long) array via a temporary.
  FUNCTION read_int_vec( ncid, varname, out ) RESULT( ok )
    INTEGER,       INTENT(IN)  :: ncid
    CHARACTER(*),  INTENT(IN)  :: varname
    INTEGER(Long), INTENT(OUT) :: out(:)
    LOGICAL :: ok
    INTEGER :: vid
    INTEGER, ALLOCATABLE :: tmp(:)
    ok = .FALSE.
    IF (.NOT. check(NF90_INQ_VARID(ncid,varname,vid),'inq '//varname)) RETURN
    ALLOCATE( tmp(SIZE(out)) )
    IF (.NOT. check(NF90_GET_VAR(ncid,vid,tmp),'get '//varname)) THEN
      DEALLOCATE(tmp); RETURN
    END IF
    out = INT(tmp, Long)
    DEALLOCATE(tmp)
    ok = .TRUE.
  END FUNCTION read_int_vec


  ! Return a dimension length, or -1 on error.
  FUNCTION get_dim( ncid, dimname ) RESULT( n )
    INTEGER,      INTENT(IN) :: ncid
    CHARACTER(*), INTENT(IN) :: dimname
    INTEGER :: n
    INTEGER :: dimid
    n = -1
    IF (.NOT. check(NF90_INQ_DIMID(ncid,dimname,dimid),'inq dim '//dimname)) RETURN
    IF (.NOT. check(NF90_INQUIRE_DIMENSION(ncid,dimid,LEN=n),'len dim '//dimname)) n = -1
  END FUNCTION get_dim


  ! netCDF status check: .TRUE. on success, otherwise displays a message.
  FUNCTION check( status, context ) RESULT( ok )
    INTEGER,      INTENT(IN) :: status
    CHARACTER(*), INTENT(IN) :: context
    LOGICAL :: ok
    ok = ( status == NF90_NOERR )
    IF ( .NOT. ok ) &
      CALL Display_Message( ROUTINE_NAME, TRIM(context)//': '//TRIM(NF90_STRERROR(status)), FAILURE )
  END FUNCTION check

END MODULE TELSEM2Atlas_netCDF_IO
