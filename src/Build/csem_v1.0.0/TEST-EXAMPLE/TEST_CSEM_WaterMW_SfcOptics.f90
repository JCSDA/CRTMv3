   PROGRAM TEST_CSEM_WaterMW_SfcOptics
    
 
      USE CSEM_Type_Kinds, ONLY: fp    => CSEM_fp, &
                                 Short => CSEM_Short
      USE CSEM_WaterMW_SfcOptics ,ONLY: yy_type=>ivar_type, &
                              CSEM_Compute_WaterMW_SfcOptics, &
                              CSEM_Compute_WaterMW_SfcOptics_TL,&
                              CSEM_Compute_WaterMW_SfcOptics_AD
      USE CSEM_Define
      USE CSEM_Model_Manager
      USE FASTEM_Coeff_Reader
      USE CSEM_LifeCycle
     
      IMPLICIT NONE
      
      INTEGER, PARAMETER  :: N_Angles = 1 
      INTEGER, PARAMETER  :: N_STOKES = 4 
      INTEGER, PARAMETER  :: SUCCESS  = 0 
      
      TYPE(CSEM_Water_Surface)    :: Surface
      TYPE(CSEM_SfcOptics_Type)   :: SfcOptics
      TYPE(CSEM_Options_Type)     :: Options
      TYPE(yy_type)               :: iVar
      
      TYPE(CSEM_Model_ID)         :: CSEM_ALG
    
      INTEGER(SHORT) :: Error_Status
     
      !Arguments
      REAL(fp) :: Sensor_Azimuth_Angle
      REAL(fp) :: Transmittance
         
      ! REAL(fp) :: Angles(N_Angles) =  (/15.0_fp,30.0_fp,45.0_fp,60.0_fp,75.0_fp/) 
      REAL(fp) :: Angles(N_Angles) =  (/45.0_fp/) 
     
      INTEGER(SHORT) ::j 
      Error_Status =CSEM_INIT("alg.list")
      CSEM_ALG = inq_model_option("MW_WATER") 
      !CSEM_ALG%NAME = "RTTOV_FASTEM_V6"
      !CALL set_model_option(CSEM_ALG) 
      
      CALL SfcOptics%Init(N_Angles=N_Angles)

      SfcOptics%Frequency          =   88.997_fp !23.800904_fp 
      SfcOptics%Angle              =   Angles
      Surface%Wind_Speed           =   11.0_fp
      Surface%Water_Temperature    =   275.0_fp
      Surface%Salinity             =   33.00_fp
      Surface%Wind_Direction       =   50.0_fp
      Sensor_Azimuth_Angle         =   0.0_fp 
      Transmittance                =   0.6_fp
      
      SfcOptics%Sensor_Azimuth_Angle = Sensor_Azimuth_Angle
      Options%Atmos%Transmittance = Transmittance 

      ! forward 
      Error_Status = CSEM_Compute_WaterMW_SfcOptics( &
                 Surface                  ,&  ! Input
                 SfcOptics                ,&  ! Input 
                 Options, ivar)

      WRITE(*,*)TRIM(CSEM_ALG%NAME)
      WRITE(*,*)'Frequency:', SfcOptics%Frequency
      WRITE(*,*)'Geometry:',  SfcOptics%Angle
      WRITE(*,*)'Emissivity(1):', SfcOptics%Emissivity(:,1)
      WRITE(*,*)'Reflectivity(1):', (SfcOptics%Reflectivity(j,1,j,1), j=1,N_Angles)
      CALL CSEM_Destroy

  END PROGRAM TEST_CSEM_WaterMW_SfcOptics
