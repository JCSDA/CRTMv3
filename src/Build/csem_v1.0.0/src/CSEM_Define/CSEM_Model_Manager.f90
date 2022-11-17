!
! CSEM_Model_Manager
!
! Module containing functions to manage the all model options already implemented in CSEM 
! and registered in the Model_Registor_File.  
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 08-17-2017
!                       ming.chen@noaa.gov
!



MODULE CSEM_Model_Manager

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  
   ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  PUBLIC :: CSEM_Model_ID
  PUBLIC :: Load_Model_Repo
  PUBLIC :: Set_Model_Option
  PUBLIC :: Inq_Model_Option
  PUBLIC :: Get_Data_Path


  TYPE CSEM_Model_ID 
    CHARACTER(LEN=10)  ::  CLASS
    CHARACTER(LEN=30)  ::  NAME
    CHARACTER(LEN=256) ::  DATA_PATH
  END TYPE CSEM_Model_ID 

  INTEGER, PARAMETER   :: N_CSEM_Models = 100
  TYPE(CSEM_Model_ID),SAVE    :: Model_LIST(N_CSEM_Models)
  LOGICAL, SAVE :: IS_INITED = .FALSE.
  ! Default algorithms 
  TYPE(CSEM_Model_ID), PRIVATE :: &
    MW_LAND_Model         =  CSEM_Model_ID("MW_LAND",   "NESDIS_Land_MW",     "./"),      & 
    MW_WATER_Model        =  CSEM_Model_ID("MW_WATER",  "NESDIS_FASTEM_V6",   "./"),      & 
    MW_SNOW_Model         =  CSEM_Model_ID("MW_SNOW",   "NESDIS_Snow_MW",     "./"),      & 
    MW_ICE_Model          =  CSEM_Model_ID("MW_ICE",    "NESDIS_Ice_MW",      "./"),      & 
    IR_LAND_Model         =  CSEM_Model_ID("IR_LAND",   "NPOESS_LUT",         "./"),      &  
    IR_WATER_Model        =  CSEM_Model_ID("IR_WATER",  "NESDIS_IRW_WuSmith", "./"),      & 
    IR_SNOW_Model         =  CSEM_Model_ID("IR_SNOW",   "NPOESS_LUT",         "./"),      & 
    IR_ICE_Model          =  CSEM_Model_ID("IR_ICE",    "NPOESS_LUT",         "./"),      & 
    VIS_LAND_Model        =  CSEM_Model_ID("VIS_LAND",  "NPOESS_LUT",         "./"),      & 
    VIS_WATER_Model       =  CSEM_Model_ID("VIS_WATER", "NPOESS_LUT",         "./"),      & 
    VIS_SNOW_Model        =  CSEM_Model_ID("VIS_SNOW",  "NPOESS_LUT",         "./"),      &      
    VIS_ICE_Model         =  CSEM_Model_ID("VIS_ICE",   "NPOESS_LUT",         "./")
  
CONTAINS

! Load_Model_Repo is to parse the Model_Registor_File
! and to load all the model options that are implemented in CSEM
!  
  FUNCTION Load_Model_Repo( Model_Registor_File )  RESULT ( Error_Status )

    CHARACTER(LEN=*), INTENT(IN) :: Model_Registor_File
    ! Function result
    INTEGER :: Error_Status
    ! local
    CHARACTER(LEN=256) :: LINE 
    CHARACTER(LEN=256) :: Model_CLASS
    CHARACTER(LEN=256) :: Model_REC(3)
    CHARACTER(LEN=1)   :: Model_ON
    TYPE(CSEM_Model_ID)  :: DEFAULT_Model
    INTEGER :: IC, ID, ILS, IRS,  IRE, IDF, IModel=1
    INTEGER :: FUNIT
    
    Error_Status = 0
    IF(IS_INITED)THEN
      RETURN
    ENDIF
    FUNIT=get_lun()
         
    OPEN(FUNIT, FILE=TRIM(Model_Registor_File), STATUS='OLD', ACTION='READ', &
           FORM='FORMATTED', IOSTAT=Error_Status)
    IF (Error_Status /= 0) THEN
       WRITE(*,*) 'Error in opening Model_Registor_File '//TRIM(Model_Registor_File)//' ....'
       RETURN
    END IF

    DO WHILE(.TRUE.)

       READ(FUNIT, '(A)', IOSTAT=Error_Status) LINE
       IF (Error_Status /= 0) THEN
         IF(Error_Status == -1) THEN
           CLOSE(FUNIT)
           IS_INITED=.TRUE.
           Error_Status = 0
        ELSE
         PRINT*,'Error reading Model_Registor_File '//TRIM(Model_Registor_File)
         PRINT*,'FUNIT: ', FUNIT
        ENDIF
         RETURN
       ENDIF
       LINE = StrCompress(LINE)
       IF(LEN(TRIM(ADJUSTL(LINE))) < 1)CYCLE
       IC=SCAN(LINE, "#")
       IF(IC == 1) CYCLE
       IF(IC>1) LINE=LINE(1:IC-1)

       ILS=SCAN(LINE, "[") 
       IRS=SCAN(LINE, "]") 
       IF(ILS >=1 .AND. IRS > ILS+1) THEN
         Model_CLASS = TRIM(LINE(ILS+1:IRS-1))
         !WRITE(*,*)'Algorithm Class: ', TRIM(Model_CLASS)
         CYCLE
       ENDIF 
       ID   = SCAN(LINE, ",") 
       IF(ID <=1) CYCLE

       IRE = 1 ; Model_REC =  ""
       DO WHILE(ID >=1 .AND. IRE <= 3) 
          Model_REC(IRE) = TRIM(ADJUSTL(LINE(1:ID-1)))
          IRE = IRE + 1
          LINE = TRIM(ADJUSTL(LINE(ID+1:)))      
          ID   = SCAN(LINE, ",") 
       ENDDO
       IF(IRE <=3 ) Model_REC(IRE) = TRIM(ADJUSTL(LINE))
    
       Model_LIST(IModel) = CSEM_Model_ID(TRIM(Model_CLASS), TRIM(Model_REC(1)),  TRIM(Model_REC(3)))
       Model_ON = TRIM(ADJUSTL(Model_REC(2))) ;  IDF = 0
       IF(LEN(Model_ON) < 1) Model_ON = '0'
       READ(Model_ON,'(I1)') IDF

       IF(IDF > 0 )THEN
          DEFAULT_Model = Model_LIST(IModel)
          CALL SET_Model_Option( DEFAULT_Model) 
       ENDIF
       IModel = IModel + 1
    END DO
    CLOSE(FUNIT)
    IS_INITED=.TRUE.
  END FUNCTION Load_Model_Repo


  SUBROUTINE SET_Model_Option(Model, verbose)

    TYPE(CSEM_Model_ID), INTENT(INOUT) :: Model
    LOGICAL, OPTIONAL :: verbose

    IF(.NOT. Valid_Check(Model)) THEN
       WRITE(*,*)"Not Valid Model%CLASS: ", TRIM(Model%CLASS), &
          " or Model%NAME: ",  TRIM(Model%NAME)
       RETURN    
    END IF
    
    IF(PRESENT(verbose)) &
        WRITE(*,*)"Setting Model%CLASS : ", TRIM(Model%CLASS), &
          "  Model%NAME: ",  TRIM(Model%NAME)

    SELECT CASE (TRIM(Model%CLASS))
     
        CASE ("MW_LAND")
          MW_LAND_Model   = Model
        CASE ("MW_SNOW")
          MW_SNOW_Model   = Model
        CASE ("MW_WATER")
          MW_WATER_Model  = Model
        CASE ("MW_ICE")
          MW_ICE_Model    = Model
  
        CASE ("IR_LAND")
          IR_LAND_Model   = Model
        CASE ("IR_SNOW")
          IR_SNOW_Model   = Model
        CASE ("IR_WATER")
          IR_WATER_Model  = Model
        CASE ("IR_ICE")
          IR_ICE_Model    = Model

        CASE ("VIS_LAND")
          VIS_LAND_Model  = Model
        CASE ("VIS_SNOW")
          VIS_SNOW_Model  = Model
        CASE ("VIS_WATER")
          VIS_WATER_Model = Model
        CASE ("VIS_ICE")
          VIS_ICE_Model   = Model
  
        CASE DEFAULT
          WRITE(*,*)'Wrong AlgType NAME ',TRIM(Model%CLASS)
    END SELECT 
     
  END SUBROUTINE SET_Model_Option
    
  TYPE(CSEM_Model_ID) FUNCTION Inq_Model_Option(ModelClass) 

    CHARACTER(LEN=*), INTENT(IN)  :: ModelClass
    TYPE(CSEM_Model_ID) :: Model

    SELECT CASE (TRIM(ModelClass))
     
        CASE ("MW_LAND")
          Model = MW_LAND_Model 
        CASE ("MW_SNOW")
          Model = MW_SNOW_Model 
        CASE ("MW_WATER")
          Model = MW_WATER_Model 
        CASE ("MW_ICE")
          Model = MW_ICE_Model 
  
        CASE ("IR_LAND")
          Model = IR_LAND_Model
        CASE ("IR_SNOW")
          Model = IR_SNOW_Model
        CASE ("IR_WATER")
          Model = IR_WATER_Model
        CASE ("IR_ICE")
          Model = IR_ICE_Model

        CASE ("VIS_LAND")
          Model = VIS_LAND_Model 
        CASE ("VIS_SNOW")
          Model = VIS_SNOW_Model  
        CASE ("VIS_WATER")
          Model = VIS_WATER_Model 
        CASE ("VIS_ICE")
          Model = VIS_ICE_Model  
  
        CASE DEFAULT
          WRITE(*,*)'Wrong ModelClass NAME'
          Model = CSEM_Model_ID(ModelClass,"","") 
  
    END SELECT 
    
    Inq_Model_Option = Model

  END FUNCTION Inq_Model_Option

  LOGICAL FUNCTION Valid_Check(Model)

    TYPE(CSEM_Model_ID), INTENT(INOUT) :: Model
    INTEGER :: i
   
    Valid_Check = .FALSE.
    DO i = 1, N_CSEM_Models
      IF(TRIM(ADJUSTL(Model%CLASS)) .EQ. TRIM(ADJUSTL(Model_LIST(i)%CLASS)) .AND. &
         TRIM(ADJUSTL(Model%NAME))  .EQ. TRIM(ADJUSTL(Model_LIST(i)%NAME))) THEN
         Model%DATA_PATH = Model_LIST(i)%DATA_PATH
         Valid_Check = .TRUE.
         RETURN
      ENDIF
    ENDDO
     
  END FUNCTION Valid_Check

  FUNCTION GET_DATA_PATH(ModelClass, ModelName) RESULT(ModelPath)

    CHARACTER(LEN=*), INTENT(IN)  :: ModelClass, ModelName
    CHARACTER(LEN=256) :: ModelPath
    INTEGER :: i
   
    ModelPath='./'
    DO i = 1, N_CSEM_Models
      IF(TRIM(ADJUSTL(ModelClass)) .EQ. TRIM(ADJUSTL(Model_LIST(i)%CLASS)) .AND. &
         TRIM(ADJUSTL(ModelName))  .EQ. TRIM(ADJUSTL(Model_LIST(i)%NAME))) THEN
         ModelPath = TRIM(Model_LIST(i)%DATA_PATH)
         RETURN
      ENDIF
    ENDDO
     
  END FUNCTION GET_DATA_PATH
    
  FUNCTION StrCompress( Input_String, n ) RESULT( Output_String )
    ! Arguments
    CHARACTER(*),      INTENT(IN)  :: Input_String
    INTEGER, OPTIONAL, INTENT(OUT) :: n
    ! Function result
    CHARACTER(LEN(Input_String)) :: Output_String
    ! Local parameters
    INTEGER, PARAMETER :: IACHAR_SPACE = 32
    INTEGER, PARAMETER :: IACHAR_TAB   = 9
    ! Local variables
    INTEGER :: i, j
    INTEGER :: IACHAR_Character

    ! Setup
    ! -----
    ! Initialise output string
    Output_String = ' '
    ! Initialise output string "useful" length counter
    j = 0

    ! Loop over string contents character by character
    ! ------------------------------------------------
    DO i = 1, LEN(Input_String)

      ! Convert the current character to its position
      ! in the ASCII collating sequence
      IACHAR_Character = IACHAR(Input_String(i:i))

      ! If the character is NOT a space ' ' or a tab '-&gt;|'
      ! copy it to the output string.
      IF ( IACHAR_Character /= IACHAR_SPACE .AND. &
           IACHAR_Character /= IACHAR_TAB         ) THEN
        j = j + 1
        Output_String(j:j) = Input_String(i:i)
      END IF

    END DO

    ! Save the non-whitespace count
    ! -----------------------------
    IF ( PRESENT(n) ) n = j

  END FUNCTION StrCompress

  FUNCTION Get_Lun() RESULT( Lun )
    INTEGER :: Lun
    LOGICAL :: Is_Open
    LOGICAL :: Existence

    ! Initialise logical unit number
    Lun = 9

    ! Start open loop for Lun Search
    Lun_Search: DO
      Lun = Lun + 1
      INQUIRE( UNIT = LUN, EXIST = Existence )
      IF ( .NOT. Existence ) THEN
        Lun = -1
        EXIT Lun_Search
      END IF
      INQUIRE( UNIT = Lun, OPENED = Is_Open )
      IF ( .NOT. Is_Open  ) EXIT Lun_Search
    END DO Lun_Search

  END FUNCTION Get_Lun

END MODULE CSEM_Model_Manager
