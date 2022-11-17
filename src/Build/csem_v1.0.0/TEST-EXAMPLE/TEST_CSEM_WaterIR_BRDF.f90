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
      
      TYPE(CSEM_Water_Surface)   :: Surface
      TYPE(CSEM_Options_Type)    :: Options
      TYPE(CSEM_SfcOptics_Type)  :: SfcOptics
      TYPE(yy_type)              :: iVar
      
      TYPE(CSEM_Model_ID)        :: CSEM_ALG
     
      !Arguments
      REAL(fp) :: Sensor_Azimuth_Angle
      REAL(fp) :: Transmittance
          
      REAL(fp) :: Angles(N_Angles) =  (/0.0_fp,15.0_fp,30.0_fp,45.0_fp,60.0_fp,75.0_fp/) 
      !REAL(fp) :: Angles(N_Angles) =  (/86.0_fp,71.0_fp,48.0_fp,21.0_fp,30.0_fp/) 
      
      INTEGER(SHORT) :: i, j, Error_Status

      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("IR_WATER") 
      !CSEM_ALG%NAME = "NESDIS_IRW_Nalli" ! "NESDIS_IRW_WuSmith" "NESDIS_IRW_Nalli"
      !CSEM_ALG%NAME = "RTTOV_IRSSEM_V2D1" ! "NESDIS_IRW_WuSmith" "NESDIS_IRW_Nalli"
      CALL set_model_option(CSEM_ALG)


      CALL SfcOptics%Init(N_Angles=N_Angles)
    
      Options%SensorObs%Sensor_Id = 'hirs4_n18'
      SfcOptics%Wavenumber = 1.0E4_fp/2.5_fp      
      !SfcOptics%Wavenumber = 668.177369738041_fp ! 668.177369738041  2666.14088520248 (hirs ch1 and 19)

      SfcOptics%Angle            =  Angles
      Surface%Wind_Speed         =  5.0_fp
      Surface%Water_Temperature  =  275.0_fp
      Surface%Salinity           =  33.00_fp
      Surface%Wind_Direction     =  150.0_fp
      Sensor_Azimuth_Angle       =  100.0_fp 
      Transmittance              =  0.6_fp
      
      SfcOptics%Is_Solar   = .TRUE.
      SfcOptics%Sensor_Zenith_Angle       =  30.0_fp
      SfcOptics%Sensor_Azimuth_Angle      =  0.0_fp
      SfcOptics%Source_Zenith_Angle       =  SfcOptics%Angle(1)
      SfcOptics%Source_Azimuth_Angle      =  180.0_fp

      OPEN (UNIT=15, FILE='IRWATER_BRDF.DAT', STATUS='UNKNOWN',FORM='FORMATTED')
      ! forward 

      DO i=1,N_Angles
        SfcOptics%Sensor_Zenith_Angle          =  SfcOptics%Angle(i)
        SfcOptics%Source_Zenith_Angle          =  SfcOptics%Angle(i)
   
        Error_Status = CSEM_Compute_WaterIR_SfcOptics( &
                 Surface                  ,&  ! Input
                 SfcOptics                ,&  ! OptionalInput 
                 Options  ,ivar)     

        WRITE(*,*)TRIM(CSEM_ALG%NAME)       
        WRITE(*,*)'Frequency:', SfcOptics%Wavenumber
        WRITE(*,*)'Geometry:',  SfcOptics%Angle
        WRITE(*,*)'Emissivity(1):', SfcOptics%Emissivity(:,1)
        DO j=1,N_Angles
          WRITE(15,'(F5.0, 4E15.4)')SfcOptics%Angle(j), SfcOptics%Reflectivity(j,1,j,1), &
           SfcOptics%Direct_Reflectivity(j,1)/3.1415 *cos(3.1415926*SfcOptics%Angle(i)/180.),cos(3.1415926*SfcOptics%Angle(i)/180.0)
        ENDDO
      ENDDO
      CLOSE(15)
  
  END PROGRAM TEST_CSEM_WaterIR_SfcOptics
