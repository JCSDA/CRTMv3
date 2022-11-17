!
! CSEM_LifeCycle
!
! This module contains CSEM life cycle functions to initialize and destroy
! the CSEM space. The initalization function is to call the Model manager to
! load all the model options.
! 

!
! CREATION HISTORY:
!       Written by:     Ming Chen, 08-17-2017
!                       ming.chen@noaa.gov
!



MODULE CSEM_LifeCycle

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds,          ONLY:  fp => CSEM_fp

  USE CSEM_Model_Manager,       ONLY:  CSEM_Model_ID, &
                                       Load_Model_Repo, &
                                       Inq_Model_Option
  USE TELSEM_Atlas_Module,      ONLY:  TELSEM_Atlas_Close
  USE CNRM_Atlas_Module,        ONLY:  CNRM_Atlas_Close
  USE CRTM_FASTEM_MODULE,       ONLY:  CRTM_FASTEM_Destroy
  USE UWIR_Atlas_Module,        ONLY:  UWIR_Atlas_Close
  
   ! Disable implicit typing
  IMPLICIT NONE

  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  PUBLIC :: CSEM_INIT
  PUBLIC :: CSEM_Destroy

  CONTAINS

  FUNCTION CSEM_INIT( Model_Registor_File )  RESULT ( Error_Status )

    CHARACTER(LEN=*), INTENT(IN) :: Model_Registor_File
    ! Function result
    INTEGER :: Error_Status
    ! local
    TYPE(CSEM_MODEL_ID) :: MODEL

    Error_Status = Load_Model_Repo(TRIM(Model_Registor_File ))
    WRITE(*,*)"CSEM Water surface models selected: "

    MODEL = Inq_Model_Option("VIS_WATER") 
    WRITE(*,*)"VIS_Water: ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("IR_WATER") 
    WRITE(*,*)"IR_Water : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("MW_WATER") 
    WRITE(*,*)"MW_Water : ", TRIM(MODEL%NAME)


    WRITE(*,*)"CSEM Land surface models selected: "
    MODEL = Inq_Model_Option("VIS_LAND") 
    WRITE(*,*)"VIS_Land : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("IR_LAND") 
    WRITE(*,*)"IR_Land  : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("MW_Land") 
    WRITE(*,*)"MW_Land  : ", TRIM(MODEL%NAME)

    WRITE(*,*)"CSEM Snow surface models selected: "
    MODEL = Inq_Model_Option("VIS_SNOW") 
    WRITE(*,*)"VIS_Snow : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("IR_SNOW") 
    WRITE(*,*)"IR_Snow  : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("MW_SNOW") 
    WRITE(*,*)"MW_Snow  : ", TRIM(MODEL%NAME)

    WRITE(*,*)"CSEM Ice surface models selected: "
    MODEL = Inq_Model_Option("VIS_ICE") 
    WRITE(*,*)"VIS_Ice  : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("IR_ICE") 
    WRITE(*,*)"IR_Ice   : ", TRIM(MODEL%NAME)
    MODEL = Inq_Model_Option("MW_ICE") 
    WRITE(*,*)"MW_Ice   : ", TRIM(MODEL%NAME)
    
  END FUNCTION CSEM_INIT

  SUBROUTINE CSEM_Destroy()
      INTEGER  :: Status
      CALL TELSEM_Atlas_Close
      CALL UWIR_Atlas_Close
      CALL CNRM_Atlas_Close
      Status = CRTM_FASTEM_Destroy()
  END SUBROUTINE CSEM_Destroy

END MODULE CSEM_LifeCycle
