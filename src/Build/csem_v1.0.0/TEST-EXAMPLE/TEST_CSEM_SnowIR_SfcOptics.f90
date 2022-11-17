   PROGRAM TEST_CSEM_SnowIR_SfcOptics

      USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
      USE CSEM_Define
      USE CSEM_Model_Manager
      USE CSEM_SnowIR_SfcOptics

      IMPLICIT NONE
      
      TYPE(CSEM_Snow_Surface)    :: Surface
      TYPE(CSEM_SfcOptics_Type)  :: SfcOptics
      TYPE(CSEM_Options_Type)    :: Options
 
      TYPE(CSEM_Model_ID)        :: CSEM_ALG

      INTEGER :: I, Error_Status
      REAL(fp), PARAMETER:: SPEED_OF_LIGHT = 2.99792458e+08_fp

      !REAL(fp) :: wvs(10)=(/0.35_fp,0.8_fp,1.4_fp,2.6_fp,3.0_fp,4.0_fp,&
      !   5.3_fp,10.3_fp,14.0_fp,18.0_fp/)
      REAL(fp) :: wvs(10)=  (/3.9_fp,6.2_fp,&
                            6.9_fp,7.3_fp,8.6_fp,9.6_fp,10.4_fp,11.2_fp,12.4_fp,13.3_fp/)
      INTEGER :: Snow_Type = 1
      
      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("IR_SNOW") 
      CALL SfcOptics%Init(N_Angles=1)
      Surface%Snow_Type         = Snow_Type
      Options%GeoInfo%Latitude  = 38.79
      Options%GeoInfo%Longitude = -102.73
      Options%GeoInfo%Month     = 5

      DO I=1,size(wvs)
         SfcOptics%Wavenumber = 1.0E+4_fp/wvs(I)
         Error_Status =  CSEM_Compute_SnowIR_SfcOptics( &     
                Surface                  ,&  ! Input
                SfcOptics                ,&  ! Output
                Options                 )    ! Input

         print*,I,wvs(I), SfcOptics%Wavenumber,SfcOptics%Emissivity(1,1),snow_type
     
      END DO
   
    END PROGRAM TEST_CSEM_SnowIR_SfcOptics
