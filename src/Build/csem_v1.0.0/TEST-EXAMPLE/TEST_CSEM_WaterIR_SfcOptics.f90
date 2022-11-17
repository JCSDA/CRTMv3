   PROGRAM TEST_CSEM_WaterIR_SfcOptics
    
 
      USE CSEM_Type_Kinds, ONLY: fp    => CSEM_fp, &
                                 Short => CSEM_Short
      USE CSEM_WaterIR_SfcOptics ,ONLY: yy_type=>ivar_type, &
                              CSEM_Compute_WaterIR_SfcOptics, &
                              CSEM_Compute_WaterIR_SfcOptics_TL,&
                              CSEM_Compute_WaterIR_SfcOptics_AD
      USE CSEM_Define
      USE CSEM_Model_Manager
      
      IMPLICIT NONE
      
      INTEGER, PARAMETER  :: N_Angles = 6 
      INTEGER, PARAMETER  :: N_STOKES = 4 
      INTEGER, PARAMETER  :: SUCCESS  = 0 
      
      TYPE(CSEM_Water_Surface)    :: Surface
      TYPE(CSEM_Options_Type)     :: Options
      TYPE(CSEM_SfcOptics_Type)   :: SfcOptics
      TYPE(yy_type)               :: iVar
      
      TYPE(CSEM_Model_ID)  :: CSEM_ALG

      INTEGER(SHORT) :: Error_Status
     
      !Arguments
      REAL(fp) :: Sensor_Azimuth_Angle
      REAL(fp) :: Transmittance
          
      REAL(fp) :: Angles(N_Angles) =  (/0.0_fp, 15.0_fp, 30.0_fp, 45.0_fp, 60.0_fp, 75.0_fp/)  
      REAL(fp) :: SST(4) =  (/275.0_fp, 285.0_fp, 295.0_fp, 310.0_fp/)  
      
      INTEGER(SHORT) :: j

      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("IR_WATER") 
      !CSEM_ALG%NAME = "NESDIS_IRW_Nalli" ! "NESDIS_IRW_WuSmith" "NESDIS_IRW_Nalli"
      !CSEM_ALG%NAME = "RTTOV_IRSSEM_V2D1" ! "NESDIS_IRW_WuSmith" "NESDIS_IRW_Nalli"
      CALL set_model_option(CSEM_ALG)

      CALL SfcOptics%Init(N_Angles=N_Angles)
    
      Options%SensorObs%Sensor_Id = 'hirs4_n18'
      SfcOptics%Wavenumber = 1.0E4_fp/3.0_fp      
      SfcOptics%Wavenumber = 668.177369738041_fp ! 668.177369738041  2666.14088520248 (hirs ch1 and 19)

      SfcOptics%Angle            =  Angles
      Surface%Wind_Speed         =  5.0_fp
      Surface%Water_Temperature  =  SST(1)
      Surface%Salinity           =  33.00_fp
      Surface%Wind_Direction     =  150.0_fp
      Sensor_Azimuth_Angle       =  100.0_fp 
      Transmittance              =  0.6_fp
      
      SfcOptics%Is_Solar   = .FALSE.
      SfcOptics%Sensor_Zenith_Angle           =  30.0_fp
      SfcOptics%Sensor_Azimuth_Angle          =  100.0_fp
      SfcOptics%Source_Zenith_Angle           =  50.0_fp
      SfcOptics%Source_Azimuth_Angle          =  80.0_fp

      ! forward 
      Error_Status = CSEM_Compute_WaterIR_SfcOptics( &
                 Surface                  ,&  ! Input
                 SfcOptics                ,&  ! OptionalInput 
                 Options  ,ivar)     
      WRITE(*,*)"IR_Water model selected : ",TRIM(CSEM_ALG%NAME)       
      WRITE(*,'(A, F12.4)')' Frequency:', SfcOptics%Wavenumber
      WRITE(*,'(A, F12.4)')' Wind_Speed:',Surface%Wind_Speed
      WRITE(*,'(A, F12.4)')' Water_Temperature:',Surface%Water_Temperature
      WRITE(*,*)'     Angles    ', ' Emissivity'
      DO j=1, N_Angles
        WRITE(*,'(F12.2, F12.4)')SfcOptics%Angle(j), SfcOptics%Emissivity(j,1)
      ENDDO

      Surface%Water_Temperature  =  SST(2)
      ! forward 
      Error_Status = CSEM_Compute_WaterIR_SfcOptics( &
                 Surface                  ,&  ! Input
                 SfcOptics                ,&  ! OptionalInput 
                 Options  ,ivar)     
      WRITE(*,'(A, F12.4)')' Water_Temperature:',Surface%Water_Temperature
      WRITE(*,*)'     Angles    ', ' Emissivity'
      DO j=1, N_Angles
        WRITE(*,'(F12.2, F12.4)')SfcOptics%Angle(j), SfcOptics%Emissivity(j,1)
      ENDDO

      Surface%Water_Temperature  =  SST(3)
      ! forward 
      Error_Status = CSEM_Compute_WaterIR_SfcOptics( &
                 Surface                  ,&  ! Input
                 SfcOptics                ,&  ! OptionalInput 
                 Options  ,ivar)     
      WRITE(*,'(A, F12.4)')' Water_Temperature:',Surface%Water_Temperature
      WRITE(*,*)'     Angles    ', ' Emissivity'
      DO j=1, N_Angles
        WRITE(*,'(F12.2, F12.4)')SfcOptics%Angle(j), SfcOptics%Emissivity(j,1)
      ENDDO     
      
     
   
  END PROGRAM TEST_CSEM_WaterIR_SfcOptics
