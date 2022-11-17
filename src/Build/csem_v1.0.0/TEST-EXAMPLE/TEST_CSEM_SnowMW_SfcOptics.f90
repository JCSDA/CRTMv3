   PROGRAM TEST_CSEM_SnowMW_SfcOptics

      USE CSEM_Type_Kinds, ONLY: fp    => CSEM_fp, &
                                 Short => CSEM_Short
      USE CSEM_SnowMW_SfcOptics
      USE CSEM_Define
      USE CSEM_Model_Manager
        
      IMPLICIT NONE
      
      TYPE(CSEM_Snow_Surface)    :: Surface
      TYPE(CSEM_SfcOptics_Type)  :: SfcOptics
      TYPE(CSEM_Options_Type)    :: Options

      TYPE(CSEM_Model_ID)        :: CSEM_ALG

      INTEGER(SHORT) :: I, Error_Status

      
      INTEGER, PARAMETER :: N_Channels = 10
      REAL(fp) :: Frequency(N_Channels) = &
         (/5.,10.,20.,30.,50.,60.,90.,110.,160.,200./)
     
      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("MW_SNOW") 
 
      CALL SfcOptics%Init(N_Angles=1)
 
      SfcOptics%Frequency = 23.8
      SfcOptics%Angle(1) = 53.0
      Surface%Top_Soil_Moisture_Content = 0.2
      Surface%Top_Soil_Temperature  = 285.0
      Surface%Snow_Temperature      = 285.0
      Surface%soil_Type             = 3 
      Surface%Snow_depth            = 10.0

     
      DO I=1, N_Channels
         SfcOptics%Frequency = Frequency(I)
         Error_Status = CSEM_Compute_SnowMW_SfcOptics( &
             Surface                        ,&  ! Input
             SfcOptics                      ,&  ! Output
             Options)                           ! Input
        WRITE(*,FMT='(F5.0,F8.2,3F18.10)')SfcOptics%Frequency,SfcOptics%Angle,&
             SfcOptics%Emissivity(1,1),SfcOptics%Emissivity(1,2)
      ENDDO
       
   
    END PROGRAM TEST_CSEM_SnowMW_SfcOptics
