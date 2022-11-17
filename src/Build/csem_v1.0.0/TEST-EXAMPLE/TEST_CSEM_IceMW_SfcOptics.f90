   PROGRAM TEST_CSEM_IceMW_SfcOptics
      USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp, &
                                Short => CSEM_Short
      USE CSEM_IceMW_SfcOptics
      USE CSEM_Define
      USE CSEM_Model_Manager
        
      IMPLICIT NONE
      
      TYPE(CSEM_Ice_Surface)      :: Surface
      TYPE(CSEM_Options_Type)     :: Options
      TYPE(CSEM_SfcOptics_Type)   :: SfcOptics
      CHARACTER(LEN=100 )         :: AlgID
      TYPE(CSEM_Model_ID)         :: CSEM_ALG

      INTEGER(SHORT) :: I
      INTEGER(SHORT) :: Error_Status

      REAL(fp) :: freqs(5)=(/20.0,30.0,50.0,90.0,160.0/)
            
      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("MW_ICE") 
      
      CALL SfcOptics%Init(N_Angles=1)
      SfcOptics%Frequency = 23.8
      SfcOptics%Angle     = 50.0
      Surface%Ice_Temperature  = 285.0
      Surface%Salinity  = 0.005
 
      DO I=1,5
         SfcOptics%Frequency = Freqs(I)
         Error_Status = CSEM_Compute_IceMW_SfcOptics(          &
             Surface                        ,&  ! Input
             SfcOptics                      ,&  ! Output
             Options)                           ! Input
         WRITE(*,FMT='(F5.0,F8.2,F5.2,3F10.3)') &
             SfcOptics%Frequency, &
             SfcOptics%Angle(1),  &
             SfcOptics%Emissivity(1,1)
      ENDDO
       
   
    END PROGRAM TEST_CSEM_IceMW_SfcOptics
