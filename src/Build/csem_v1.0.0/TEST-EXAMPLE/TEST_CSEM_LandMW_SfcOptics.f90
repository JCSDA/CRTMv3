PROGRAM TEST_CSEM_LandMW_SfcOptics
  ! Module usage

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Define
  USE CSEM_LandMW_SfcOptics
  USE CSEM_Model_Manager

       
  IMPLICIT NONE
      
  TYPE(CSEM_Land_Surface)     :: Surface
  TYPE(CSEM_SfcOptics_Type)   :: SfcOptics
  TYPE(CSEM_Options_Type)     :: Options
  TYPE(iVar_type)             :: iVar
 
  INTEGER, PARAMETER :: N_Channels = 10
  INTEGER :: Error_Status
 
  REAL(fp) :: Frequency(N_Channels) = (/5.,10.,20.,30.,50.,60.,90.,110., 160.,200./)
  REAL(fp) :: Angles(6)  = (/5.,10.,20.,30.,50.,60./)
 
  INTEGER :: J = 1

  TYPE(CSEM_Model_ID)  :: CSEM_ALG
   
  Error_Status = load_model_repo("alg.list")
  CSEM_ALG = inq_model_option("MW_LAND") 

  CALL SfcOptics%Init(N_Angles=6)
  SfcOptics%Angle = Angles
  
  Options%GeoInfo%Latitude   =   45.25
  Options%GeoInfo%Longitude  =  -89.50
  Options%GeoInfo%Month      =   8

 
  Surface%Top_Soil_Moisture      =   0.2
  Surface%vegetation_Fraction    =   0.5
  Surface%Top_Soil_Temperature   =   285.0
  Surface%Land_Skin_Temperature  =   290.0
  Surface%LAI                    =   3.5
  Surface%soil_Type              =   1
  Surface%Vegetation_Type        =   2
  
  SfcOptics%Frequency = Frequency(J)

  Error_Status = CSEM_Compute_LandMW_SfcOptics( &
      Surface                   ,&  ! Input
      SfcOptics                 ,&  ! Output
      Options                   ,&  ! Input 
      iVar)                         ! Output
  WRITE(*,*)TRIM(CSEM_ALG%NAME)
  WRITE(*,*)SfcOptics%Emissivity(1,:)
  
  CSEM_ALG%NAME = "TELSEM2_ATLAS"
  CALL set_model_option(CSEM_ALG)
  Error_Status = CSEM_Compute_LandMW_SfcOptics( &
      Surface                   ,&  ! Input
      SfcOptics                 ,&  ! Output
      Options                   ,&  ! Input 
      iVar)                         ! Output
  WRITE(*,*)TRIM(CSEM_ALG%NAME)
  WRITE(*,*)SfcOptics%Emissivity(1,:)

END PROGRAM TEST_CSEM_LandMW_SfcOptics
