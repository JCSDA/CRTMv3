!
! ODSSU_netCDF_IO
!
! Module containing routines to read and write ODSSU TauCoeff data
! files in netCDF format.
!
! ODSSU files are an ODSSU-container holding M = n_TC_CellPressures
! sub-coefficient sets (one per CO2 cell-pressure epoch). For SSU
! this is the ODPS sub-algorithm; the ODAS branch is not exercised
! by the netCDF path.
!
! Schema (flat, single namespace):
!
!   Global attrs : Release, Version, Algorithm (=ODSSU=3), subAlgorithm,
!                  Sensor_Id, WMO_Satellite_Id, WMO_Sensor_Id, title,
!                  history, comment
!   Dimensions   : n_Layers, n_Levels(=n_Layers+1), n_Components,
!                  n_Absorbers, n_Channels, n_Coeffs, n_OPIndex,
!                  n_OCoeffs (only if >0), n_TC_CellPressures,
!                  n_Ref_CellPressures
!   Container    : Sensor_Channel(L), Sensor_Type, Absorber_ID(Jm),
!                  TC_CellPressure(M,L), Ref_Time(N),
!                  Ref_CellPressure(N,L)
!   Per-set ODPS, stacked on M:
!                  Group_Index(M), Component_ID(J,M),
!                  Ref_Level_Pressure(K+1,M), Ref_Pressure(K,M),
!                  Ref_Temperature(K,M), Ref_Absorber(K,Jm,M),
!                  Min_Absorber(K,Jm,M), Max_Absorber(K,Jm,M),
!                  n_Predictors(J,L,M), Pos_Index(J,L,M),
!                  ODPS_Coefficients(Iuse,M)
!   OPTRAN extras (n_OCoeffs > 0 only):
!                  Alpha(M), Alpha_C1(M), Alpha_C2(M),
!                  OComponent_Index(M), OSignificance(L,M), Order(L,M),
!                  OP_Index(OI+1,L,M), OPos_Index(L,M), OC(n_OCoeffs,M)
!
! Asserts that all M sub-ODPS sets within a file share the same ODPS
! dimensions (true for SSU coefficient training, where only cell
! pressure varies between sets).
!

MODULE ODSSU_netCDF_IO

  ! ------------------
  ! Environment set up
  ! ------------------
  USE Type_Kinds      , ONLY: Long, Double, Single, fp
  USE Message_Handler , ONLY: SUCCESS, FAILURE, WARNING, INFORMATION, &
                              Display_Message
  USE File_Utility    , ONLY: File_Exists
  USE ODPS_Define     , ONLY: ODPS_type, &
                              Allocate_ODPS, &
                              Allocate_ODPS_OPTRAN, &
                              Destroy_ODPS, &
                              Associated_ODPS
  USE ODSSU_Define    , ONLY: ODSSU_type, &
                              Allocate_ODSSU, &
                              Destroy_ODSSU, &
                              CheckRelease_ODSSU, &
                              CheckAlgorithm_ODSSU, &
                              Info_ODSSU, &
                              ODAS_ALGORITHM, ODPS_ALGORITHM
  USE CRTM_Parameters , ONLY: ODSSU_ALGORITHM
  USE netcdf
  USE netCDF_Utility  , ONLY: Put_netCDF_Variable, &
                              Get_netCDF_Variable, &
                              Remove_NULL_Characters

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: Inquire_ODSSU_netCDF
  PUBLIC :: Read_ODSSU_netCDF
  PUBLIC :: Write_ODSSU_netCDF

  ! Module parameters
  INTEGER, PARAMETER :: ML = 1024
  INTEGER, PARAMETER :: SET = 1
  REAL(Double), PARAMETER :: ZERO = 0.0_Double

  ! Global attribute names
  CHARACTER(*), PARAMETER :: TITLE_GATTNAME             = 'title'
  CHARACTER(*), PARAMETER :: HISTORY_GATTNAME           = 'history'
  CHARACTER(*), PARAMETER :: COMMENT_GATTNAME           = 'comment'
  CHARACTER(*), PARAMETER :: RELEASE_GATTNAME           = 'Release'
  CHARACTER(*), PARAMETER :: VERSION_GATTNAME           = 'Version'
  CHARACTER(*), PARAMETER :: ALGORITHM_GATTNAME         = 'Algorithm'
  CHARACTER(*), PARAMETER :: SUBALGORITHM_GATTNAME      = 'subAlgorithm'
  CHARACTER(*), PARAMETER :: SENSOR_ID_GATTNAME         = 'Sensor_Id'
  CHARACTER(*), PARAMETER :: WMO_SATELLITE_ID_GATTNAME  = 'WMO_Satellite_Id'
  CHARACTER(*), PARAMETER :: WMO_SENSOR_ID_GATTNAME     = 'WMO_Sensor_Id'
  CHARACTER(*), PARAMETER :: WRITE_MODULE_HISTORY_GATTNAME   = 'write_module_history'
  CHARACTER(*), PARAMETER :: CREATION_DATE_AND_TIME_GATTNAME = 'creation_date_and_time'
  CHARACTER(*), PARAMETER :: MODULE_HISTORY = 'ODSSU_netCDF_IO (CRTM v3 REL-3.2.0)'

  ! Dimension names
  CHARACTER(*), PARAMETER :: LAYER_DIMNAME      = 'n_Layers'
  CHARACTER(*), PARAMETER :: LEVEL_DIMNAME      = 'n_Levels'
  CHARACTER(*), PARAMETER :: COMPONENT_DIMNAME  = 'n_Components'
  CHARACTER(*), PARAMETER :: ABSORBER_DIMNAME   = 'n_Absorbers'
  CHARACTER(*), PARAMETER :: CHANNEL_DIMNAME    = 'n_Channels'
  CHARACTER(*), PARAMETER :: COEFF_DIMNAME      = 'n_Coeffs'
  CHARACTER(*), PARAMETER :: ODASPRED_DIMNAME   = 'n_OPIndex'
  CHARACTER(*), PARAMETER :: ODASCOEFF_DIMNAME  = 'n_OCoeffs'
  CHARACTER(*), PARAMETER :: TC_CP_DIMNAME      = 'n_TC_CellPressures'
  CHARACTER(*), PARAMETER :: REF_CP_DIMNAME     = 'n_Ref_CellPressures'

  ! Container variable names
  CHARACTER(*), PARAMETER :: SENSOR_CHANNEL_VARNAME    = 'Sensor_Channel'
  CHARACTER(*), PARAMETER :: SENSOR_TYPE_VARNAME       = 'Sensor_Type'
  CHARACTER(*), PARAMETER :: ABSORBER_ID_VARNAME       = 'Absorber_ID'
  CHARACTER(*), PARAMETER :: TC_CELLPRESSURE_VARNAME   = 'TC_CellPressure'
  CHARACTER(*), PARAMETER :: REF_TIME_VARNAME          = 'Ref_Time'
  CHARACTER(*), PARAMETER :: REF_CELLPRESSURE_VARNAME  = 'Ref_CellPressure'

  ! Per-set ODPS variable names (stacked on M)
  CHARACTER(*), PARAMETER :: GROUP_INDEX_VARNAME       = 'Group_Index'
  CHARACTER(*), PARAMETER :: COMPONENT_ID_VARNAME      = 'Component_ID'
  CHARACTER(*), PARAMETER :: REF_LEVEL_PRESSURE_VARNAME= 'Ref_Level_Pressure'
  CHARACTER(*), PARAMETER :: REF_PRESSURE_VARNAME      = 'Ref_Pressure'
  CHARACTER(*), PARAMETER :: REF_TEMPERATURE_VARNAME   = 'Ref_Temperature'
  CHARACTER(*), PARAMETER :: REF_ABSORBER_VARNAME      = 'Ref_Absorber'
  CHARACTER(*), PARAMETER :: MIN_ABSORBER_VARNAME      = 'Min_Absorber'
  CHARACTER(*), PARAMETER :: MAX_ABSORBER_VARNAME      = 'Max_Absorber'
  CHARACTER(*), PARAMETER :: N_PREDICTORS_VARNAME      = 'n_Predictors'
  CHARACTER(*), PARAMETER :: POS_INDEX_VARNAME         = 'Pos_Index'
  CHARACTER(*), PARAMETER :: ODPS_COEFFICIENTS_VARNAME = 'ODPS_Coefficients'

  ! OPTRAN per-set names
  CHARACTER(*), PARAMETER :: ALPHA_VARNAME             = 'Alpha'
  CHARACTER(*), PARAMETER :: ALPHA_C1_VARNAME          = 'Alpha_C1'
  CHARACTER(*), PARAMETER :: ALPHA_C2_VARNAME          = 'Alpha_C2'
  CHARACTER(*), PARAMETER :: OCOMPONENT_INDEX_VARNAME  = 'OComponent_Index'
  CHARACTER(*), PARAMETER :: OSIGNIFICANCE_VARNAME     = 'OSignificance'
  CHARACTER(*), PARAMETER :: ORDER_VARNAME             = 'Order'
  CHARACTER(*), PARAMETER :: OP_INDEX_VARNAME          = 'OP_Index'
  CHARACTER(*), PARAMETER :: OPOS_INDEX_VARNAME        = 'OPos_Index'
  CHARACTER(*), PARAMETER :: ODAS_COEFFICIENTS_VARNAME = 'OC'

CONTAINS

!--------------------------------------------------------------------------------
! Inquire_ODSSU_netCDF
!--------------------------------------------------------------------------------
  FUNCTION Inquire_ODSSU_netCDF( NC_Filename        , &
                                 n_Layers           , &
                                 n_Components       , &
                                 n_Absorbers        , &
                                 n_Channels         , &
                                 n_Coeffs           , &
                                 n_OPIndex          , &
                                 n_OCoeffs          , &
                                 n_TC_CellPressures , &
                                 n_Ref_CellPressures, &
                                 Release            , &
                                 Version            , &
                                 subAlgorithm       , &
                                 Sensor_Id          , &
                                 WMO_Satellite_Id   , &
                                 WMO_Sensor_Id      , &
                                 Message_Log        ) RESULT( Error_Status )
    CHARACTER(*),           INTENT(IN)  :: NC_Filename
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Layers
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Components
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Absorbers
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Channels
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Coeffs
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_OPIndex
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_OCoeffs
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_TC_CellPressures
    INTEGER,      OPTIONAL, INTENT(OUT) :: n_Ref_CellPressures
    INTEGER,      OPTIONAL, INTENT(OUT) :: Release
    INTEGER,      OPTIONAL, INTENT(OUT) :: Version
    INTEGER,      OPTIONAL, INTENT(OUT) :: subAlgorithm
    CHARACTER(*), OPTIONAL, INTENT(OUT) :: Sensor_Id
    INTEGER,      OPTIONAL, INTENT(OUT) :: WMO_Satellite_Id
    INTEGER,      OPTIONAL, INTENT(OUT) :: WMO_Sensor_Id
    CHARACTER(*), OPTIONAL, INTENT(IN)  :: Message_Log
    INTEGER :: Error_Status

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Inquire_ODSSU_netCDF'
    CHARACTER(ML) :: Message
    INTEGER :: NC_FileID, NF90_Status, DimID
    INTEGER :: dimlen
    INTEGER :: alg_in
    CHARACTER(5000) :: GAttString

    Error_Status = SUCCESS

    IF ( .NOT. File_Exists( TRIM(NC_Filename) ) ) THEN
      Message = 'File '//TRIM(NC_Filename)//' not found.'
      CALL Inquire_Cleanup( CloseNeeded=.FALSE. ); RETURN
    END IF

    NF90_Status = NF90_OPEN( TRIM(NC_Filename), NF90_NOWRITE, NC_FileID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      Message = 'Error opening '//TRIM(NC_Filename)//' - '//&
                TRIM(NF90_STRERROR(NF90_Status))
      CALL Inquire_Cleanup( CloseNeeded=.FALSE. ); RETURN
    END IF

    ! Mandatory: Algorithm must be ODSSU
    NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, ALGORITHM_GATTNAME, alg_in )
    IF ( NF90_Status /= NF90_NOERR .OR. alg_in /= ODSSU_ALGORITHM ) THEN
      Message = 'Algorithm attribute missing or not ODSSU in '//TRIM(NC_Filename)
      CALL Inquire_Cleanup( CloseNeeded=.TRUE. ); RETURN
    END IF

    ! Dimensions
    IF ( PRESENT(n_Layers) ) THEN
      IF ( .NOT. Get_Dim( LAYER_DIMNAME, n_Layers ) ) RETURN
    END IF
    IF ( PRESENT(n_Components) ) THEN
      IF ( .NOT. Get_Dim( COMPONENT_DIMNAME, n_Components ) ) RETURN
    END IF
    IF ( PRESENT(n_Absorbers) ) THEN
      IF ( .NOT. Get_Dim( ABSORBER_DIMNAME, n_Absorbers ) ) RETURN
    END IF
    IF ( PRESENT(n_Channels) ) THEN
      IF ( .NOT. Get_Dim( CHANNEL_DIMNAME, n_Channels ) ) RETURN
    END IF
    IF ( PRESENT(n_Coeffs) ) THEN
      ! n_Coeffs may be absent when 0
      NF90_Status = NF90_INQ_DIMID( NC_FileID, COEFF_DIMNAME, DimID )
      IF ( NF90_Status == NF90_NOERR ) THEN
        NF90_Status = NF90_INQUIRE_DIMENSION( NC_FileID, DimID, LEN=dimlen )
        n_Coeffs = dimlen
      ELSE
        n_Coeffs = 0
      END IF
    END IF
    IF ( PRESENT(n_OPIndex) ) THEN
      NF90_Status = NF90_INQ_DIMID( NC_FileID, ODASPRED_DIMNAME, DimID )
      IF ( NF90_Status == NF90_NOERR ) THEN
        NF90_Status = NF90_INQUIRE_DIMENSION( NC_FileID, DimID, LEN=dimlen )
        n_OPIndex = dimlen
      ELSE
        n_OPIndex = 0
      END IF
    END IF
    IF ( PRESENT(n_OCoeffs) ) THEN
      NF90_Status = NF90_INQ_DIMID( NC_FileID, ODASCOEFF_DIMNAME, DimID )
      IF ( NF90_Status == NF90_NOERR ) THEN
        NF90_Status = NF90_INQUIRE_DIMENSION( NC_FileID, DimID, LEN=dimlen )
        n_OCoeffs = dimlen
      ELSE
        n_OCoeffs = 0
      END IF
    END IF
    IF ( PRESENT(n_TC_CellPressures) ) THEN
      IF ( .NOT. Get_Dim( TC_CP_DIMNAME, n_TC_CellPressures ) ) RETURN
    END IF
    IF ( PRESENT(n_Ref_CellPressures) ) THEN
      IF ( .NOT. Get_Dim( REF_CP_DIMNAME, n_Ref_CellPressures ) ) RETURN
    END IF

    ! Global attributes (optional)
    IF ( PRESENT(Release) ) THEN
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, RELEASE_GATTNAME, Release )
      IF ( NF90_Status /= NF90_NOERR ) Release = 0
    END IF
    IF ( PRESENT(Version) ) THEN
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, VERSION_GATTNAME, Version )
      IF ( NF90_Status /= NF90_NOERR ) Version = 0
    END IF
    IF ( PRESENT(subAlgorithm) ) THEN
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, SUBALGORITHM_GATTNAME, subAlgorithm )
      IF ( NF90_Status /= NF90_NOERR ) subAlgorithm = 0
    END IF
    IF ( PRESENT(Sensor_Id) ) THEN
      Sensor_Id = ' '
      GAttString = ' '
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, SENSOR_ID_GATTNAME, GAttString )
      IF ( NF90_Status == NF90_NOERR ) THEN
        CALL Remove_NULL_Characters( GAttString )
        Sensor_Id = GAttString(1:MIN( LEN(Sensor_Id), LEN_TRIM(GAttString) ))
      END IF
    END IF
    IF ( PRESENT(WMO_Satellite_Id) ) THEN
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, WMO_SATELLITE_ID_GATTNAME, WMO_Satellite_Id )
      IF ( NF90_Status /= NF90_NOERR ) WMO_Satellite_Id = -1
    END IF
    IF ( PRESENT(WMO_Sensor_Id) ) THEN
      NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, WMO_SENSOR_ID_GATTNAME, WMO_Sensor_Id )
      IF ( NF90_Status /= NF90_NOERR ) WMO_Sensor_Id = -1
    END IF

    NF90_Status = NF90_CLOSE( NC_FileID )

  CONTAINS

    LOGICAL FUNCTION Get_Dim( Dim_Name, Dim_Out )
      CHARACTER(*), INTENT(IN) :: Dim_Name
      INTEGER     , INTENT(OUT):: Dim_Out
      INTEGER :: Stat, ID, L
      Get_Dim = .FALSE.
      Stat = NF90_INQ_DIMID( NC_FileID, Dim_Name, ID )
      IF ( Stat /= NF90_NOERR ) THEN
        Message = 'Dim '//Dim_Name//' missing in '//TRIM(NC_Filename)
        CALL Inquire_Cleanup( CloseNeeded=.TRUE. ); RETURN
      END IF
      Stat = NF90_INQUIRE_DIMENSION( NC_FileID, ID, LEN=L )
      IF ( Stat /= NF90_NOERR ) THEN
        Message = 'Dim length read failed for '//Dim_Name//' in '//TRIM(NC_Filename)
        CALL Inquire_Cleanup( CloseNeeded=.TRUE. ); RETURN
      END IF
      Dim_Out = L
      Get_Dim = .TRUE.
    END FUNCTION Get_Dim

    SUBROUTINE Inquire_Cleanup( CloseNeeded )
      LOGICAL, INTENT(IN) :: CloseNeeded
      INTEGER :: ignore
      IF ( CloseNeeded ) ignore = NF90_CLOSE( NC_FileID )
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            Error_Status, &
                            Message_Log=Message_Log )
    END SUBROUTINE Inquire_Cleanup

  END FUNCTION Inquire_ODSSU_netCDF


!--------------------------------------------------------------------------------
! Write_ODSSU_netCDF
!--------------------------------------------------------------------------------
  FUNCTION Write_ODSSU_netCDF( NC_Filename, &
                               ODSSU      , &
                               Title      , &
                               History    , &
                               Comment    , &
                               Quiet      , &
                               Message_Log) RESULT( Error_Status )
    CHARACTER(*),           INTENT(IN) :: NC_Filename
    TYPE(ODSSU_type),       INTENT(IN) :: ODSSU
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Title
    CHARACTER(*), OPTIONAL, INTENT(IN) :: History
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Comment
    INTEGER     , OPTIONAL, INTENT(IN) :: Quiet
    CHARACTER(*), OPTIONAL, INTENT(IN) :: Message_Log
    INTEGER :: Error_Status

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Write_ODSSU_netCDF'
    CHARACTER(ML) :: Message
    LOGICAL :: Noisy
    INTEGER :: NC_FileID, NF90_Status, ignore
    INTEGER :: i, K, J, Jm, L, Iuse, OI, OC
    INTEGER :: M, N
    INTEGER :: dim_K, dim_Lev, dim_J, dim_Jm, dim_L, dim_Coeff, &
               dim_OI, dim_OC, dim_M, dim_N
    INTEGER :: VarID
    CHARACTER(8)  :: cdate
    CHARACTER(10) :: ctime
    CHARACTER(5)  :: czone

    ! Per-set buffers
    INTEGER(Long), ALLOCATABLE :: Group_Index_Buf(:)
    INTEGER(Long), ALLOCATABLE :: Component_ID_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Level_Pressure_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Pressure_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Temperature_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Absorber_Buf(:,:,:)
    REAL(fp),      ALLOCATABLE :: Min_Absorber_Buf(:,:,:)
    REAL(fp),      ALLOCATABLE :: Max_Absorber_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: n_Predictors_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: Pos_Index_Buf(:,:,:)
    REAL(Single),  ALLOCATABLE :: ODPS_Coefficients_Buf(:,:)
    ! OPTRAN buffers
    REAL(fp),      ALLOCATABLE :: Alpha_Buf(:), Alpha_C1_Buf(:), Alpha_C2_Buf(:)
    INTEGER(Long), ALLOCATABLE :: OComponent_Index_Buf(:)
    INTEGER(Long), ALLOCATABLE :: OSignificance_Buf(:,:)
    INTEGER(Long), ALLOCATABLE :: Order_Buf(:,:)
    INTEGER(Long), ALLOCATABLE :: OP_Index_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: OPos_Index_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: OC_Buf(:,:)

    Error_Status = SUCCESS
    Noisy = .TRUE.
    IF ( PRESENT(Quiet) ) THEN
      IF ( Quiet == SET ) Noisy = .FALSE.
    END IF

    ! Only ODPS sub-algorithm is implemented for the netCDF path
    IF ( ODSSU%subAlgorithm /= ODPS_ALGORITHM ) THEN
      Message = 'Write_ODSSU_netCDF only supports subAlgorithm=ODPS for now.'
      CALL Write_Cleanup( CloseNeeded=.FALSE. ); RETURN
    END IF
    IF ( ODSSU%n_TC_CellPressures < 1 .OR. ODSSU%n_Ref_CellPressures < 1 .OR. &
         ODSSU%n_Channels < 1 .OR. ODSSU%n_Absorbers < 1 ) THEN
      Message = 'One or more ODSSU dimensions are < 1.'
      CALL Write_Cleanup( CloseNeeded=.FALSE. ); RETURN
    END IF

    M  = ODSSU%n_TC_CellPressures
    N  = ODSSU%n_Ref_CellPressures
    L  = ODSSU%n_Channels
    Jm = ODSSU%n_Absorbers

    ! Verify all M sub-sets share dimensions, and grab them from set 1
    K    = ODSSU%ODPS(1)%n_Layers
    J    = ODSSU%ODPS(1)%n_Components
    Iuse = ODSSU%ODPS(1)%n_Coeffs
    OI   = ODSSU%ODPS(1)%n_OPIndex      ! actual array dim is OI+1
    OC   = ODSSU%ODPS(1)%n_OCoeffs
    DO i = 2, M
      IF ( ODSSU%ODPS(i)%n_Layers     /= K    .OR. &
           ODSSU%ODPS(i)%n_Components /= J    .OR. &
           ODSSU%ODPS(i)%n_Absorbers  /= Jm   .OR. &
           ODSSU%ODPS(i)%n_Channels   /= L    .OR. &
           ODSSU%ODPS(i)%n_Coeffs     /= Iuse .OR. &
           ODSSU%ODPS(i)%n_OPIndex    /= OI   .OR. &
           ODSSU%ODPS(i)%n_OCoeffs    /= OC ) THEN
        WRITE(Message,'("ODPS sub-set #",i0," dims differ from sub-set #1 - flat schema requires uniform dims")') i
        CALL Write_Cleanup( CloseNeeded=.FALSE. ); RETURN
      END IF
    END DO

    ! Create file
    NF90_Status = NF90_CREATE( TRIM(NC_Filename), NF90_NETCDF4, NC_FileID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      Message = 'Error creating '//TRIM(NC_Filename)//' - '//TRIM(NF90_STRERROR(NF90_Status))
      CALL Write_Cleanup( CloseNeeded=.FALSE. ); RETURN
    END IF

    ! Define dimensions
    IF ( DefDim( LAYER_DIMNAME    , K   , dim_K    ) /= 0 ) RETURN
    IF ( DefDim( LEVEL_DIMNAME    , K+1 , dim_Lev  ) /= 0 ) RETURN
    IF ( DefDim( COMPONENT_DIMNAME, J   , dim_J    ) /= 0 ) RETURN
    IF ( DefDim( ABSORBER_DIMNAME , Jm  , dim_Jm   ) /= 0 ) RETURN
    IF ( DefDim( CHANNEL_DIMNAME  , L   , dim_L    ) /= 0 ) RETURN
    IF ( DefDim( TC_CP_DIMNAME    , M   , dim_M    ) /= 0 ) RETURN
    IF ( DefDim( REF_CP_DIMNAME   , N   , dim_N    ) /= 0 ) RETURN
    IF ( Iuse > 0 ) THEN
      IF ( DefDim( COEFF_DIMNAME, Iuse, dim_Coeff ) /= 0 ) RETURN
    END IF
    IF ( OC > 0 ) THEN
      IF ( DefDim( ODASPRED_DIMNAME , OI+1, dim_OI ) /= 0 ) RETURN
      IF ( DefDim( ODASCOEFF_DIMNAME, OC  , dim_OC ) /= 0 ) RETURN
    END IF

    ! Global attributes
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, WRITE_MODULE_HISTORY_GATTNAME, MODULE_HISTORY )
    CALL DATE_AND_TIME( cdate, ctime, czone )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, CREATION_DATE_AND_TIME_GATTNAME, &
                                cdate(1:4)//'/'//cdate(5:6)//'/'//cdate(7:8)//', '// &
                                ctime(1:2)//':'//ctime(3:4)//':'//ctime(5:6)//' '// &
                                czone//'UTC' )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, ALGORITHM_GATTNAME   , ODSSU_ALGORITHM )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, RELEASE_GATTNAME     , INT(ODSSU%Release) )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, VERSION_GATTNAME     , INT(ODSSU%Version) )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, SUBALGORITHM_GATTNAME, INT(ODSSU%subAlgorithm) )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, SENSOR_ID_GATTNAME   , TRIM(ODSSU%Sensor_Id) )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, WMO_SATELLITE_ID_GATTNAME, INT(ODSSU%WMO_Satellite_ID) )
    NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, WMO_SENSOR_ID_GATTNAME   , INT(ODSSU%WMO_Sensor_ID) )
    IF ( PRESENT(Title)   ) NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, TITLE_GATTNAME  , Title )
    IF ( PRESENT(History) ) NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, HISTORY_GATTNAME, History )
    IF ( PRESENT(Comment) ) NF90_Status = NF90_PUT_ATT( NC_FileID, NF90_GLOBAL, COMMENT_GATTNAME, Comment )

    ! Define variables: container
    IF ( DefVar( SENSOR_CHANNEL_VARNAME  , NF90_INT   , (/ dim_L /),       VarID ) /= 0 ) RETURN
    IF ( DefVar( SENSOR_TYPE_VARNAME     , NF90_INT   , Empty(),           VarID ) /= 0 ) RETURN
    IF ( DefVar( ABSORBER_ID_VARNAME     , NF90_INT   , (/ dim_Jm /),      VarID ) /= 0 ) RETURN
    IF ( DefVar( TC_CELLPRESSURE_VARNAME , NF90_DOUBLE, (/ dim_M, dim_L /),VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_TIME_VARNAME        , NF90_DOUBLE, (/ dim_N /),       VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_CELLPRESSURE_VARNAME, NF90_DOUBLE, (/ dim_N, dim_L /),VarID ) /= 0 ) RETURN
    ! Per-set ODPS
    IF ( DefVar( GROUP_INDEX_VARNAME       , NF90_INT   , (/ dim_M /),                 VarID ) /= 0 ) RETURN
    IF ( DefVar( COMPONENT_ID_VARNAME      , NF90_INT   , (/ dim_J, dim_M /),          VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_LEVEL_PRESSURE_VARNAME, NF90_DOUBLE, (/ dim_Lev, dim_M /),        VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_PRESSURE_VARNAME      , NF90_DOUBLE, (/ dim_K, dim_M /),          VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_TEMPERATURE_VARNAME   , NF90_DOUBLE, (/ dim_K, dim_M /),          VarID ) /= 0 ) RETURN
    IF ( DefVar( REF_ABSORBER_VARNAME      , NF90_DOUBLE, (/ dim_K, dim_Jm, dim_M /),  VarID ) /= 0 ) RETURN
    IF ( DefVar( MIN_ABSORBER_VARNAME      , NF90_DOUBLE, (/ dim_K, dim_Jm, dim_M /),  VarID ) /= 0 ) RETURN
    IF ( DefVar( MAX_ABSORBER_VARNAME      , NF90_DOUBLE, (/ dim_K, dim_Jm, dim_M /),  VarID ) /= 0 ) RETURN
    IF ( DefVar( N_PREDICTORS_VARNAME      , NF90_INT   , (/ dim_J, dim_L, dim_M /),   VarID ) /= 0 ) RETURN
    IF ( DefVar( POS_INDEX_VARNAME         , NF90_INT   , (/ dim_J, dim_L, dim_M /),   VarID ) /= 0 ) RETURN
    IF ( Iuse > 0 ) THEN
      IF ( DefVar( ODPS_COEFFICIENTS_VARNAME, NF90_FLOAT, (/ dim_Coeff, dim_M /),      VarID ) /= 0 ) RETURN
    END IF
    IF ( OC > 0 ) THEN
      IF ( DefVar( ALPHA_VARNAME             , NF90_DOUBLE, (/ dim_M /),                  VarID ) /= 0 ) RETURN
      IF ( DefVar( ALPHA_C1_VARNAME          , NF90_DOUBLE, (/ dim_M /),                  VarID ) /= 0 ) RETURN
      IF ( DefVar( ALPHA_C2_VARNAME          , NF90_DOUBLE, (/ dim_M /),                  VarID ) /= 0 ) RETURN
      IF ( DefVar( OCOMPONENT_INDEX_VARNAME  , NF90_INT   , (/ dim_M /),                  VarID ) /= 0 ) RETURN
      IF ( DefVar( OSIGNIFICANCE_VARNAME     , NF90_INT   , (/ dim_L, dim_M /),           VarID ) /= 0 ) RETURN
      IF ( DefVar( ORDER_VARNAME             , NF90_INT   , (/ dim_L, dim_M /),           VarID ) /= 0 ) RETURN
      IF ( DefVar( OP_INDEX_VARNAME          , NF90_INT   , (/ dim_OI, dim_L, dim_M /),   VarID ) /= 0 ) RETURN
      IF ( DefVar( OPOS_INDEX_VARNAME        , NF90_INT   , (/ dim_L, dim_M /),           VarID ) /= 0 ) RETURN
      IF ( DefVar( ODAS_COEFFICIENTS_VARNAME , NF90_DOUBLE, (/ dim_OC, dim_M /),          VarID ) /= 0 ) RETURN
    END IF

    NF90_Status = NF90_ENDDEF( NC_FileID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      Message = 'Error ending define mode - '//TRIM(NF90_STRERROR(NF90_Status))
      CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN
    END IF

    ! ---- Write container variables ----
    Error_Status = Put_netCDF_Variable( NC_FileID, SENSOR_CHANNEL_VARNAME, ODSSU%Sensor_Channel )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put Sensor_Channel'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, SENSOR_TYPE_VARNAME, INT(ODSSU%Sensor_Type) )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put Sensor_Type'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, ABSORBER_ID_VARNAME, ODSSU%Absorber_ID )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put Absorber_ID'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, TC_CELLPRESSURE_VARNAME, ODSSU%TC_CellPressure )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put TC_CellPressure'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_TIME_VARNAME, ODSSU%Ref_Time )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put Ref_Time'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_CELLPRESSURE_VARNAME, ODSSU%Ref_CellPressure )
    IF ( Error_Status /= SUCCESS ) THEN; Message = 'put Ref_CellPressure'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF

    ! ---- Pack per-set arrays ----
    ALLOCATE( Group_Index_Buf(M)              , &
              Component_ID_Buf(J,M)           , &
              Ref_Level_Pressure_Buf(K+1, M)  , &
              Ref_Pressure_Buf(K, M)          , &
              Ref_Temperature_Buf(K, M)       , &
              Ref_Absorber_Buf(K, Jm, M)      , &
              Min_Absorber_Buf(K, Jm, M)      , &
              Max_Absorber_Buf(K, Jm, M)      , &
              n_Predictors_Buf(J, L, M)       , &
              Pos_Index_Buf(J, L, M) )
    IF ( Iuse > 0 ) ALLOCATE( ODPS_Coefficients_Buf(Iuse, M) )
    DO i = 1, M
      Group_Index_Buf(i)                = ODSSU%ODPS(i)%Group_Index
      Component_ID_Buf(:, i)            = ODSSU%ODPS(i)%Component_ID
      Ref_Level_Pressure_Buf(:, i)      = ODSSU%ODPS(i)%Ref_Level_Pressure
      Ref_Pressure_Buf(:, i)            = ODSSU%ODPS(i)%Ref_Pressure
      Ref_Temperature_Buf(:, i)         = ODSSU%ODPS(i)%Ref_Temperature
      Ref_Absorber_Buf(:, :, i)         = ODSSU%ODPS(i)%Ref_Absorber
      Min_Absorber_Buf(:, :, i)         = ODSSU%ODPS(i)%Min_Absorber
      Max_Absorber_Buf(:, :, i)         = ODSSU%ODPS(i)%Max_Absorber
      n_Predictors_Buf(:, :, i)         = ODSSU%ODPS(i)%n_Predictors
      Pos_Index_Buf(:, :, i)            = ODSSU%ODPS(i)%Pos_Index
      IF ( Iuse > 0 ) ODPS_Coefficients_Buf(:, i) = ODSSU%ODPS(i)%C
    END DO

    Error_Status = Put_netCDF_Variable( NC_FileID, GROUP_INDEX_VARNAME       , Group_Index_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Group_Index'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, COMPONENT_ID_VARNAME      , Component_ID_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Component_ID'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_LEVEL_PRESSURE_VARNAME, Ref_Level_Pressure_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Ref_Level_Pressure'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_PRESSURE_VARNAME      , Ref_Pressure_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Ref_Pressure'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_TEMPERATURE_VARNAME   , Ref_Temperature_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Ref_Temperature'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, REF_ABSORBER_VARNAME      , Ref_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Ref_Absorber'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, MIN_ABSORBER_VARNAME      , Min_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Min_Absorber'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, MAX_ABSORBER_VARNAME      , Max_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Max_Absorber'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, N_PREDICTORS_VARNAME      , n_Predictors_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put n_Predictors'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    Error_Status = Put_netCDF_Variable( NC_FileID, POS_INDEX_VARNAME         , Pos_Index_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='put Pos_Index'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    IF ( Iuse > 0 ) THEN
      Error_Status = Put_netCDF_Variable( NC_FileID, ODPS_COEFFICIENTS_VARNAME, ODPS_Coefficients_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put ODPS_Coefficients'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    END IF

    IF ( OC > 0 ) THEN
      ALLOCATE( Alpha_Buf(M), Alpha_C1_Buf(M), Alpha_C2_Buf(M), &
                OComponent_Index_Buf(M), &
                OSignificance_Buf(L, M), Order_Buf(L, M), &
                OP_Index_Buf(0:OI, L, M), OPos_Index_Buf(L, M), &
                OC_Buf(OC, M) )
      DO i = 1, M
        Alpha_Buf(i)               = ODSSU%ODPS(i)%Alpha
        Alpha_C1_Buf(i)            = ODSSU%ODPS(i)%Alpha_C1
        Alpha_C2_Buf(i)            = ODSSU%ODPS(i)%Alpha_C2
        OComponent_Index_Buf(i)    = ODSSU%ODPS(i)%OComponent_Index
        OSignificance_Buf(:, i)    = ODSSU%ODPS(i)%OSignificance
        Order_Buf(:, i)            = ODSSU%ODPS(i)%Order
        OP_Index_Buf(:, :, i)      = ODSSU%ODPS(i)%OP_Index
        OPos_Index_Buf(:, i)       = ODSSU%ODPS(i)%OPos_Index
        OC_Buf(:, i)               = ODSSU%ODPS(i)%OC
      END DO
      Error_Status = Put_netCDF_Variable( NC_FileID, ALPHA_VARNAME           , Alpha_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put Alpha'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, ALPHA_C1_VARNAME        , Alpha_C1_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put Alpha_C1'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, ALPHA_C2_VARNAME        , Alpha_C2_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put Alpha_C2'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, OCOMPONENT_INDEX_VARNAME, OComponent_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put OComponent_Index'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, OSIGNIFICANCE_VARNAME   , OSignificance_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put OSignificance'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, ORDER_VARNAME           , Order_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put Order'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, OP_INDEX_VARNAME        , OP_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put OP_Index'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, OPOS_INDEX_VARNAME      , OPos_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put OPos_Index'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
      Error_Status = Put_netCDF_Variable( NC_FileID, ODAS_COEFFICIENTS_VARNAME, OC_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='put OC'; CALL Write_Cleanup( CloseNeeded=.TRUE. ); RETURN; END IF
    END IF

    NF90_Status = NF90_CLOSE( NC_FileID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      CALL Display_Message( ROUTINE_NAME, &
                            'Error closing netCDF ODSSU file '//TRIM(NC_Filename), &
                            WARNING, Message_Log=Message_Log )
    END IF

    IF ( Noisy ) THEN
      CALL Info_ODSSU( ODSSU, Message )
      CALL Display_Message( ROUTINE_NAME, &
                            'FILE: '//TRIM(NC_Filename)//'; '//TRIM(Message), &
                            INFORMATION, Message_Log=Message_Log )
    END IF

  CONTAINS

    PURE FUNCTION Empty() RESULT(arr)
      INTEGER :: arr(0)
    END FUNCTION Empty

    INTEGER FUNCTION DefDim( Name, Sz, ID )
      CHARACTER(*), INTENT(IN)  :: Name
      INTEGER     , INTENT(IN)  :: Sz
      INTEGER     , INTENT(OUT) :: ID
      INTEGER :: stat
      stat = NF90_DEF_DIM( NC_FileID, Name, Sz, ID )
      IF ( stat /= NF90_NOERR ) THEN
        Message = 'DefDim '//Name//' - '//TRIM(NF90_STRERROR(stat))
        CALL Write_Cleanup( CloseNeeded=.TRUE. )
        DefDim = 1
      ELSE
        DefDim = 0
      END IF
    END FUNCTION DefDim

    INTEGER FUNCTION DefVar( Name, NCType, DimIDs, VID )
      CHARACTER(*), INTENT(IN)  :: Name
      INTEGER     , INTENT(IN)  :: NCType
      INTEGER     , INTENT(IN)  :: DimIDs(:)
      INTEGER     , INTENT(OUT) :: VID
      INTEGER :: stat
      IF ( SIZE(DimIDs) == 0 ) THEN
        stat = NF90_DEF_VAR( NC_FileID, Name, NCType, varid=VID )
      ELSE
        stat = NF90_DEF_VAR( NC_FileID, Name, NCType, dimids=DimIDs, varid=VID )
      END IF
      IF ( stat /= NF90_NOERR ) THEN
        Message = 'DefVar '//Name//' - '//TRIM(NF90_STRERROR(stat))
        CALL Write_Cleanup( CloseNeeded=.TRUE. )
        DefVar = 1
      ELSE
        DefVar = 0
      END IF
    END FUNCTION DefVar

    SUBROUTINE Write_Cleanup( CloseNeeded )
      LOGICAL, INTENT(IN) :: CloseNeeded
      INTEGER :: ignore_local
      IF ( CloseNeeded ) ignore_local = NF90_CLOSE( NC_FileID )
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            Error_Status, &
                            Message_Log=Message_Log )
    END SUBROUTINE Write_Cleanup

  END FUNCTION Write_ODSSU_netCDF


!--------------------------------------------------------------------------------
! Read_ODSSU_netCDF
!--------------------------------------------------------------------------------
  FUNCTION Read_ODSSU_netCDF( NC_Filename, &
                              ODSSU      , &
                              Quiet      , &
                              Message_Log) RESULT( Error_Status )
    CHARACTER(*),           INTENT(IN)     :: NC_Filename
    TYPE(ODSSU_type),       INTENT(IN OUT) :: ODSSU
    INTEGER     , OPTIONAL, INTENT(IN)     :: Quiet
    CHARACTER(*), OPTIONAL, INTENT(IN)     :: Message_Log
    INTEGER :: Error_Status

    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'Read_ODSSU_netCDF'
    CHARACTER(ML)   :: Message
    CHARACTER(5000) :: GAttString
    LOGICAL :: Noisy
    INTEGER :: NC_FileID, NF90_Status, Destroy_Status, ignore
    INTEGER :: i, K, J, Jm, L, Iuse, OI_p1, OC
    INTEGER :: M, N
    INTEGER :: Release_in, Version_in, sub_in, Algorithm_in
    INTEGER :: WMO_Sat, WMO_Sens, Sensor_Type
    INTEGER :: tmp_int

    INTEGER(Long), ALLOCATABLE :: Group_Index_Buf(:)
    INTEGER(Long), ALLOCATABLE :: Component_ID_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Level_Pressure_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Pressure_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Temperature_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Ref_Absorber_Buf(:,:,:)
    REAL(fp),      ALLOCATABLE :: Min_Absorber_Buf(:,:,:)
    REAL(fp),      ALLOCATABLE :: Max_Absorber_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: n_Predictors_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: Pos_Index_Buf(:,:,:)
    REAL(Single),  ALLOCATABLE :: ODPS_Coefficients_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: Alpha_Buf(:), Alpha_C1_Buf(:), Alpha_C2_Buf(:)
    INTEGER(Long), ALLOCATABLE :: OComponent_Index_Buf(:)
    INTEGER(Long), ALLOCATABLE :: OSignificance_Buf(:,:)
    INTEGER(Long), ALLOCATABLE :: Order_Buf(:,:)
    INTEGER(Long), ALLOCATABLE :: OP_Index_Buf(:,:,:)
    INTEGER(Long), ALLOCATABLE :: OPos_Index_Buf(:,:)
    REAL(fp),      ALLOCATABLE :: OC_Buf(:,:)

    Error_Status = SUCCESS
    Noisy = .TRUE.
    IF ( PRESENT(Quiet) ) THEN
      IF ( Quiet == SET ) Noisy = .FALSE.
    END IF

    IF ( .NOT. File_Exists( TRIM(NC_Filename) ) ) THEN
      Message = 'File '//TRIM(NC_Filename)//' not found.'; Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
      RETURN
    END IF

    ! Read dimensions, sub-algorithm, and release
    Error_Status = Inquire_ODSSU_netCDF( NC_Filename, &
                                         n_Layers           = K     , &
                                         n_Components       = J     , &
                                         n_Absorbers        = Jm    , &
                                         n_Channels         = L     , &
                                         n_Coeffs           = Iuse  , &
                                         n_OPIndex          = OI_p1 , &
                                         n_OCoeffs          = OC    , &
                                         n_TC_CellPressures = M     , &
                                         n_Ref_CellPressures= N     , &
                                         Release            = Release_in , &
                                         Version            = Version_in , &
                                         subAlgorithm       = sub_in     , &
                                         WMO_Satellite_Id   = WMO_Sat    , &
                                         WMO_Sensor_Id      = WMO_Sens   , &
                                         Message_Log        = Message_Log )
    IF ( Error_Status /= SUCCESS ) RETURN

    IF ( sub_in /= ODPS_ALGORITHM ) THEN
      Message = 'Read_ODSSU_netCDF only supports subAlgorithm=ODPS for now.'
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
      RETURN
    END IF

    ! Set release/sub-alg before allocating (Allocate_ODSSU keys on subAlgorithm)
    ODSSU%Release      = Release_in
    ODSSU%Version      = Version_in
    ODSSU%subAlgorithm = sub_in

    Error_Status = CheckRelease_ODSSU( ODSSU, Message_Log=Message_Log )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'ODSSU Release check failed for '//TRIM(NC_Filename)
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
      RETURN
    END IF

    Error_Status = Allocate_ODSSU( Jm, L, M, N, ODSSU, Message_Log=Message_Log )
    IF ( Error_Status /= SUCCESS ) THEN
      Message = 'Allocate_ODSSU failed for '//TRIM(NC_Filename)
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
      RETURN
    END IF

    ! Allocate each ODPS sub-set
    DO i = 1, M
      Error_Status = Allocate_ODPS( K, J, Jm, L, Iuse, ODSSU%ODPS(i), Message_Log=Message_Log )
      IF ( Error_Status /= SUCCESS ) THEN
        WRITE(Message,'("Allocate_ODPS failed for sub-set ",i0)') i
        CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
        RETURN
      END IF
      IF ( OC > 0 ) THEN
        Error_Status = Allocate_ODPS_OPTRAN( OC, ODSSU%ODPS(i), Message_Log=Message_Log )
        IF ( Error_Status /= SUCCESS ) THEN
          WRITE(Message,'("Allocate_ODPS_OPTRAN failed for sub-set ",i0)') i
          CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
          RETURN
        END IF
      END IF
    END DO

    ! Open file for reading
    NF90_Status = NF90_OPEN( TRIM(NC_Filename), NF90_NOWRITE, NC_FileID )
    IF ( NF90_Status /= NF90_NOERR ) THEN
      Message = 'Error opening '//TRIM(NC_Filename); Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, TRIM(Message), Error_Status, Message_Log=Message_Log )
      RETURN
    END IF

    ! Sensor_Id global attribute (to populate ODSSU and each ODPS)
    GAttString = ' '
    NF90_Status = NF90_GET_ATT( NC_FileID, NF90_GLOBAL, SENSOR_ID_GATTNAME, GAttString )
    IF ( NF90_Status == NF90_NOERR ) THEN
      CALL Remove_NULL_Characters( GAttString )
      ODSSU%Sensor_Id = GAttString(1:MIN( LEN(ODSSU%Sensor_Id), LEN_TRIM(GAttString) ))
    END IF
    ODSSU%WMO_Satellite_ID = WMO_Sat
    ODSSU%WMO_Sensor_ID    = WMO_Sens

    ! Container variables
    Error_Status = Get_netCDF_Variable( NC_FileID, SENSOR_CHANNEL_VARNAME, ODSSU%Sensor_Channel )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Sensor_Channel'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, SENSOR_TYPE_VARNAME, tmp_int )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Sensor_Type'; CALL Read_Cleanup(); RETURN; END IF
    ODSSU%Sensor_Type = tmp_int
    Error_Status = Get_netCDF_Variable( NC_FileID, ABSORBER_ID_VARNAME, ODSSU%Absorber_ID )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Absorber_ID'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, TC_CELLPRESSURE_VARNAME, ODSSU%TC_CellPressure )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get TC_CellPressure'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_TIME_VARNAME, ODSSU%Ref_Time )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_Time'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_CELLPRESSURE_VARNAME, ODSSU%Ref_CellPressure )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_CellPressure'; CALL Read_Cleanup(); RETURN; END IF

    ! Per-set arrays
    ALLOCATE( Group_Index_Buf(M), Component_ID_Buf(J,M), &
              Ref_Level_Pressure_Buf(K+1, M), Ref_Pressure_Buf(K, M), &
              Ref_Temperature_Buf(K, M), Ref_Absorber_Buf(K, Jm, M), &
              Min_Absorber_Buf(K, Jm, M), Max_Absorber_Buf(K, Jm, M), &
              n_Predictors_Buf(J, L, M), Pos_Index_Buf(J, L, M) )
    IF ( Iuse > 0 ) ALLOCATE( ODPS_Coefficients_Buf(Iuse, M) )

    Error_Status = Get_netCDF_Variable( NC_FileID, GROUP_INDEX_VARNAME       , Group_Index_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Group_Index'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, COMPONENT_ID_VARNAME      , Component_ID_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Component_ID'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_LEVEL_PRESSURE_VARNAME, Ref_Level_Pressure_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_Level_Pressure'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_PRESSURE_VARNAME      , Ref_Pressure_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_Pressure'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_TEMPERATURE_VARNAME   , Ref_Temperature_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_Temperature'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, REF_ABSORBER_VARNAME      , Ref_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Ref_Absorber'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, MIN_ABSORBER_VARNAME      , Min_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Min_Absorber'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, MAX_ABSORBER_VARNAME      , Max_Absorber_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Max_Absorber'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, N_PREDICTORS_VARNAME      , n_Predictors_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get n_Predictors'; CALL Read_Cleanup(); RETURN; END IF
    Error_Status = Get_netCDF_Variable( NC_FileID, POS_INDEX_VARNAME         , Pos_Index_Buf )
    IF ( Error_Status /= SUCCESS ) THEN; Message='get Pos_Index'; CALL Read_Cleanup(); RETURN; END IF
    IF ( Iuse > 0 ) THEN
      Error_Status = Get_netCDF_Variable( NC_FileID, ODPS_COEFFICIENTS_VARNAME, ODPS_Coefficients_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get ODPS_Coefficients'; CALL Read_Cleanup(); RETURN; END IF
    END IF

    DO i = 1, M
      ODSSU%ODPS(i)%Group_Index        = Group_Index_Buf(i)
      ODSSU%ODPS(i)%Component_ID       = Component_ID_Buf(:, i)
      ODSSU%ODPS(i)%Ref_Level_Pressure = Ref_Level_Pressure_Buf(:, i)
      ODSSU%ODPS(i)%Ref_Pressure       = Ref_Pressure_Buf(:, i)
      ODSSU%ODPS(i)%Ref_Temperature    = Ref_Temperature_Buf(:, i)
      ODSSU%ODPS(i)%Ref_Absorber       = Ref_Absorber_Buf(:, :, i)
      ODSSU%ODPS(i)%Min_Absorber       = Min_Absorber_Buf(:, :, i)
      ODSSU%ODPS(i)%Max_Absorber       = Max_Absorber_Buf(:, :, i)
      ODSSU%ODPS(i)%n_Predictors       = n_Predictors_Buf(:, :, i)
      ODSSU%ODPS(i)%Pos_Index          = Pos_Index_Buf(:, :, i)
      IF ( Iuse > 0 ) ODSSU%ODPS(i)%C  = ODPS_Coefficients_Buf(:, i)
      ! Stamp scalar/identifying fields onto each sub-set so downstream
      ! code that reads e.g. ODSSU%ODPS(i)%Sensor_Id sees the right value.
      ODSSU%ODPS(i)%Release          = ODSSU%Release
      ODSSU%ODPS(i)%Version          = ODSSU%Version
      ODSSU%ODPS(i)%Sensor_Id        = ODSSU%Sensor_Id
      ODSSU%ODPS(i)%WMO_Satellite_ID = ODSSU%WMO_Satellite_ID
      ODSSU%ODPS(i)%WMO_Sensor_ID    = ODSSU%WMO_Sensor_ID
      ODSSU%ODPS(i)%Sensor_Channel   = ODSSU%Sensor_Channel
      ODSSU%ODPS(i)%Absorber_ID      = ODSSU%Absorber_ID
    END DO

    IF ( OC > 0 ) THEN
      ALLOCATE( Alpha_Buf(M), Alpha_C1_Buf(M), Alpha_C2_Buf(M), &
                OComponent_Index_Buf(M), &
                OSignificance_Buf(L, M), Order_Buf(L, M), &
                OP_Index_Buf(0:OI_p1-1, L, M), OPos_Index_Buf(L, M), &
                OC_Buf(OC, M) )
      Error_Status = Get_netCDF_Variable( NC_FileID, ALPHA_VARNAME           , Alpha_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get Alpha'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, ALPHA_C1_VARNAME        , Alpha_C1_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get Alpha_C1'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, ALPHA_C2_VARNAME        , Alpha_C2_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get Alpha_C2'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, OCOMPONENT_INDEX_VARNAME, OComponent_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get OComponent_Index'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, OSIGNIFICANCE_VARNAME   , OSignificance_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get OSignificance'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, ORDER_VARNAME           , Order_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get Order'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, OP_INDEX_VARNAME        , OP_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get OP_Index'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, OPOS_INDEX_VARNAME      , OPos_Index_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get OPos_Index'; CALL Read_Cleanup(); RETURN; END IF
      Error_Status = Get_netCDF_Variable( NC_FileID, ODAS_COEFFICIENTS_VARNAME, OC_Buf )
      IF ( Error_Status /= SUCCESS ) THEN; Message='get OC'; CALL Read_Cleanup(); RETURN; END IF
      DO i = 1, M
        ODSSU%ODPS(i)%Alpha            = Alpha_Buf(i)
        ODSSU%ODPS(i)%Alpha_C1         = Alpha_C1_Buf(i)
        ODSSU%ODPS(i)%Alpha_C2         = Alpha_C2_Buf(i)
        ODSSU%ODPS(i)%OComponent_Index = OComponent_Index_Buf(i)
        ODSSU%ODPS(i)%OSignificance    = OSignificance_Buf(:, i)
        ODSSU%ODPS(i)%Order            = Order_Buf(:, i)
        ODSSU%ODPS(i)%OP_Index         = OP_Index_Buf(:, :, i)
        ODSSU%ODPS(i)%OPos_Index       = OPos_Index_Buf(:, i)
        ODSSU%ODPS(i)%OC               = OC_Buf(:, i)
      END DO
    END IF

    NF90_Status = NF90_CLOSE( NC_FileID )

    IF ( Noisy ) THEN
      CALL Info_ODSSU( ODSSU, Message )
      CALL Display_Message( ROUTINE_NAME, &
                            'FILE: '//TRIM(NC_Filename)//'; '//TRIM(Message), &
                            INFORMATION, Message_Log=Message_Log )
    END IF

  CONTAINS

    SUBROUTINE Read_Cleanup()
      INTEGER :: stat
      stat = NF90_CLOSE( NC_FileID )
      Destroy_Status = Destroy_ODSSU( ODSSU, Message_Log=Message_Log )
      Error_Status = FAILURE
      CALL Display_Message( ROUTINE_NAME, &
                            TRIM(Message), &
                            Error_Status, &
                            Message_Log=Message_Log )
    END SUBROUTINE Read_Cleanup

  END FUNCTION Read_ODSSU_netCDF

END MODULE ODSSU_netCDF_IO
