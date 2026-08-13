!
! TELSEM2Atlas_netCDF_IO
!
! Module for reading and writing the CRTM-native netCDF TELSEM2 atlas coefficient
! file. The file holds all twelve monthly atlases stacked along the data
! dimension with per-month offsets. The equal-area grid geometry and reverse
! lookup are reconstructed on load. Release 2 files additionally carry the
! atlas uncertainty content (per-cell emissivity error std and per-class1
! inter-channel correlations); Release 1 files are accepted with the
! uncertainty flagged unavailable.
!

MODULE TELSEM2Atlas_netCDF_IO

  ! -----------------
  ! Environment setup
  ! -----------------
  USE Type_Kinds         , ONLY: fp, Double, Long
  USE Message_Handler    , ONLY: SUCCESS, FAILURE, INFORMATION, Display_Message
  USE File_Utility       , ONLY: File_Exists
  USE TELSEM2Atlas_Define, ONLY: TELSEM2Atlas_type              , &
                                 TELSEM2Atlas_Associated        , &
                                 TELSEM2Atlas_Create            , &
                                 TELSEM2Atlas_Create_Uncertainty, &
                                 TELSEM2Atlas_ValidRelease      , &
                                 TELSEM2Atlas_Destroy           , &
                                 N_CLASS1
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
  ! Dimension names
  CHARACTER(*), PARAMETER :: CHANNEL_DIMNAME = 'n_channels'
  CHARACTER(*), PARAMETER :: BAND_DIMNAME    = 'n_latitude_bands'
  CHARACTER(*), PARAMETER :: MONTH_DIMNAME   = 'n_months'
  CHARACTER(*), PARAMETER :: DATA_DIMNAME    = 'n_data'
  CHARACTER(*), PARAMETER :: CLASS1_DIMNAME  = 'n_class1'          ! Release >= 2
  ! Variable names
  CHARACTER(*), PARAMETER :: MONTHCOUNT_VARNAME = 'month_data_count'
  CHARACTER(*), PARAMETER :: MONTHOFF_VARNAME   = 'month_offset'
  CHARACTER(*), PARAMETER :: CELLNUM_VARNAME    = 'cell_number'
  CHARACTER(*), PARAMETER :: CLASS1_VARNAME     = 'class1'
  CHARACTER(*), PARAMETER :: CLASS2_VARNAME     = 'class2'
  CHARACTER(*), PARAMETER :: EMIS_VARNAME       = 'emissivity'
  CHARACTER(*), PARAMETER :: EMISERR_VARNAME    = 'emissivity_error' ! Release >= 2
  CHARACTER(*), PARAMETER :: CORREL_VARNAME     = 'correlation'      ! Release >= 2
  ! Global attribute names
  CHARACTER(*), PARAMETER :: RELEASE_GATTNAME    = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_GATTNAME    = 'Version'
  CHARACTER(*), PARAMETER :: RESOLUTION_GATTNAME = 'Resolution'
  CHARACTER(*), PARAMETER :: NCELLS_GATTNAME     = 'n_Cells'

  ! Largest number of elements handed to the netCDF Fortran interface in one
  ! call. The atlas is big: n_data is 2,770,889, so cell_number alone is 11 MB
  ! and emissivity is 155 MB. Reading either whole is not safe, because
  ! nf90_get_var takes an assumed-shape dummy and passes it down to an F77
  ! layer that takes an assumed-size one, and the compiler materialises a
  ! contiguous copy-in temporary to bridge the two. That temporary is created
  ! inside the netCDF library's own compiled code, so it obeys the flags netCDF
  ! was built with and not this project's: -heap-arrays here does not reach it.
  ! With Intel it lands on the stack and overflows the 8 MB a stock Linux
  ! gives you, and CRTM_Init dies with a SIGSEGV raised inside libnetcdff.
  ! Reading in bounded slices keeps every temporary small however netCDF was
  ! built. 262144 elements is 1 MB of INTEGER(4) and 2 MB of REAL(8).
  INTEGER, PARAMETER :: MAX_READ_ELEMENTS = 262144


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
    INTEGER :: dimid_chan, dimid_band, dimid_month, dimid_data, dimid_class1
    INTEGER :: vid_mcount, vid_moff, vid_cell, vid_c1, vid_c2, vid_emis
    INTEGER :: vid_err, vid_correl

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
    IF ( Atlas%Has_Uncertainty ) THEN
      IF (.NOT. check(NF90_DEF_DIM(ncid,CLASS1_DIMNAME,INT(N_CLASS1),dimid_class1),'def class1')) GOTO 900
    END IF

    ! Variable definitions
    IF (.NOT. check(NF90_DEF_VAR(ncid,MONTHCOUNT_VARNAME,NF90_INT,dimid_month,vid_mcount),'def mcount')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,MONTHOFF_VARNAME,NF90_INT,dimid_month,vid_moff),'def moff')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CELLNUM_VARNAME,NF90_INT,dimid_data,vid_cell),'def cell')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CLASS1_VARNAME,NF90_INT,dimid_data,vid_c1),'def c1')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,CLASS2_VARNAME,NF90_INT,dimid_data,vid_c2),'def c2')) GOTO 900
    IF (.NOT. check(NF90_DEF_VAR(ncid,EMIS_VARNAME,NF90_DOUBLE,[dimid_data,dimid_chan],vid_emis),'def emis')) GOTO 900
    IF ( Atlas%Has_Uncertainty ) THEN
      IF (.NOT. check(NF90_DEF_VAR(ncid,EMISERR_VARNAME,NF90_DOUBLE,[dimid_data,dimid_chan],vid_err),'def emis_err')) GOTO 900
      IF (.NOT. check(NF90_DEF_VAR(ncid,CORREL_VARNAME,NF90_DOUBLE,[dimid_class1,dimid_chan,dimid_chan],vid_correl), &
                      'def correl')) GOTO 900
    END IF

    IF (.NOT. check(NF90_ENDDEF(ncid),'enddef')) GOTO 900

    ! Write data
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_mcount,INT(Atlas%Month_Data_Count)),'put mcount')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_moff,INT(Atlas%Month_Offset)),'put moff')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_cell,INT(Atlas%Cell_Number)),'put cell')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_c1,INT(Atlas%Class1)),'put c1')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_c2,INT(Atlas%Class2)),'put c2')) GOTO 900
    IF (.NOT. check(NF90_PUT_VAR(ncid,vid_emis,Atlas%Emissivity),'put emis')) GOTO 900
    IF ( Atlas%Has_Uncertainty ) THEN
      IF (.NOT. check(NF90_PUT_VAR(ncid,vid_err,Atlas%Emis_Err),'put emis_err')) GOTO 900
      IF (.NOT. check(NF90_PUT_VAR(ncid,vid_correl,Atlas%Correlation),'put correl')) GOTO 900
    END IF

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
!
! Release handling: Release 2 files must carry the uncertainty variables
! (emissivity_error, correlation), which are read and validated. Release 1
! files are accepted with Has_Uncertainty left .FALSE.. Anything outside the
! supported release range fails.
!--------------------------------------------------------------------------------
  FUNCTION TELSEM2Atlas_netCDF_ReadFile( Filename, Atlas ) RESULT( err_stat )
    ! Arguments
    CHARACTER(*),            INTENT(IN)  :: Filename
    TYPE(TELSEM2Atlas_type), INTENT(OUT) :: Atlas
    ! Function result
    INTEGER :: err_stat
    ! Local variables
    ! Tolerance for correlation-matrix validation: the source data are written
    ! with F5.2 precision, so 1e-3 comfortably covers representation error
    ! while still rejecting a structurally wrong file.
    REAL(Double), PARAMETER :: CORREL_TOL = 1.0E-3_Double
    INTEGER :: ncid, vid
    INTEGER :: j
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

    ! Reject unsupported releases before reading any data
    IF ( .NOT. TELSEM2Atlas_ValidRelease(Atlas) ) THEN
      CALL Display_Message( ROUTINE_NAME, &
        'Unsupported Release in '//TRIM(Filename), FAILURE )
      CALL TELSEM2Atlas_Destroy(Atlas); RETURN
    END IF

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
    IF (.NOT. read_emis(ncid,EMIS_VARNAME,Atlas%Emissivity)) GOTO 900

    ! Uncertainty content (Release >= 2 files must carry it)
    IF ( rel >= 2 ) THEN
      CALL TELSEM2Atlas_Create_Uncertainty( Atlas )
      IF ( .NOT. Atlas%Has_Uncertainty ) THEN
        CALL Display_Message( ROUTINE_NAME, 'Uncertainty allocation failed', FAILURE )
        GOTO 900
      END IF
      IF (.NOT. read_emis(ncid,EMISERR_VARNAME,Atlas%Emis_Err)) GOTO 900
      ! The correlation array is small (N_CLASS1 x n_Channels x n_Channels =
      ! 490 doubles) -- a single read is safe.
      IF (.NOT. check(NF90_INQ_VARID(ncid,CORREL_VARNAME,vid),'inq '//CORREL_VARNAME)) GOTO 900
      IF (.NOT. check(NF90_GET_VAR(ncid,vid,Atlas%Correlation),'get '//CORREL_VARNAME)) GOTO 900
      ! Validate: stds must be non-negative, correlations bounded with a unit
      ! diagonal. Reject a corrupt/wrong file here rather than compute a
      ! non-positive-definite covariance at query time.
      IF ( ANY( Atlas%Emis_Err < 0.0_Double ) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
          'negative emissivity_error values in '//TRIM(Filename), FAILURE )
        GOTO 900
      END IF
      IF ( ANY( ABS(Atlas%Correlation) > 1.0_Double + CORREL_TOL ) ) THEN
        CALL Display_Message( ROUTINE_NAME, &
          'correlation values outside [-1,1] in '//TRIM(Filename), FAILURE )
        GOTO 900
      END IF
      DO j = 1, INT(n_Channels)
        IF ( ANY( ABS(Atlas%Correlation(:,j,j) - 1.0_Double) > CORREL_TOL ) ) THEN
          CALL Display_Message( ROUTINE_NAME, &
            'correlation diagonal is not unity in '//TRIM(Filename), FAILURE )
          GOTO 900
        END IF
      END DO
    END IF

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

  ! Read an integer vector variable into an INTEGER(Long) array via a temporary,
  ! in slices no larger than MAX_READ_ELEMENTS. See that parameter for why the
  ! whole array cannot be read in one call.
  FUNCTION read_int_vec( ncid, varname, out ) RESULT( ok )
    INTEGER,       INTENT(IN)  :: ncid
    CHARACTER(*),  INTENT(IN)  :: varname
    INTEGER(Long), INTENT(OUT) :: out(:)
    LOGICAL :: ok
    INTEGER :: vid, n, i, nslice
    INTEGER, ALLOCATABLE :: tmp(:)
    ok = .FALSE.
    IF (.NOT. check(NF90_INQ_VARID(ncid,varname,vid),'inq '//varname)) RETURN
    n = SIZE(out)
    ALLOCATE( tmp(MIN(n,MAX_READ_ELEMENTS)) )
    i = 1
    DO WHILE ( i <= n )
      nslice = MIN( MAX_READ_ELEMENTS, n-i+1 )
      IF (.NOT. check(NF90_GET_VAR(ncid,vid,tmp(1:nslice), &
                                   start=[i],count=[nslice]),'get '//varname)) THEN
        DEALLOCATE(tmp); RETURN
      END IF
      out(i:i+nslice-1) = INT(tmp(1:nslice), Long)
      i = i + nslice
    END DO
    DEALLOCATE(tmp)
    ok = .TRUE.
  END FUNCTION read_int_vec


  ! Read a (n_Data x n_Channels) double variable in slices, one channel at a
  ! time and each channel in chunks. Same reason as read_int_vec.
  FUNCTION read_emis( ncid, varname, out ) RESULT( ok )
    INTEGER,      INTENT(IN)  :: ncid
    CHARACTER(*), INTENT(IN)  :: varname
    REAL(Double), INTENT(OUT) :: out(:,:)
    LOGICAL :: ok
    INTEGER :: vid, n, nc, i, k, nslice
    ok = .FALSE.
    IF (.NOT. check(NF90_INQ_VARID(ncid,varname,vid),'inq '//varname)) RETURN
    n  = SIZE(out,1)
    nc = SIZE(out,2)
    DO k = 1, nc
      i = 1
      DO WHILE ( i <= n )
        nslice = MIN( MAX_READ_ELEMENTS, n-i+1 )
        IF (.NOT. check(NF90_GET_VAR(ncid,vid,out(i:i+nslice-1,k), &
                                     start=[i,k],count=[nslice,1]),'get '//varname)) RETURN
        i = i + nslice
      END DO
    END DO
    ok = .TRUE.
  END FUNCTION read_emis


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
