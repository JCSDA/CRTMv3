   PROGRAM TEST_CSEM_LandVis_SfcOptics

      USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
      USE CSEM_Define
      USE CSEM_Model_Manager
      USE CSEM_LandVis_SfcOptics ,ONLY:  &
                              CSEM_Compute_LandVis_SfcOptics, &
                              CSEM_Compute_LandVis_SfcOptics_TL,&
                              CSEM_Compute_LandVis_SfcOptics_AD
  
      IMPLICIT NONE
      
      TYPE(CSEM_Land_Surface)    :: Surface
      TYPE(CSEM_SfcOptics_Type)  :: SfcOptics
      TYPE(CSEM_Options_Type)    :: Options

      TYPE(CSEM_Model_ID)  :: CSEM_ALG
      INTEGER :: I, Error_Status
      
      REAL(fp) :: wvs(6)=(/0.47_fp,0.51_fp,0.64_fp,0.86_fp,1.6_fp,2.3_fp/)
      REAL(fp) :: angls(6)=(/0.0_fp,15.0_fp,30.0_fp,45.0_fp,60.0_fp,75.0_fp/)
  
      REAL(fp) :: wavelength
      INTEGER(KIND=1) :: surface_type = 6

      CHARACTER(LEN=110) :: IGBP_FILE='/data/jcsda/mchen/CSEM/fix/IGBPa_2006.nc'
      Surface%Land_Cover_Type = surface_type

      Options%GeoInfo%Latitude  = 38.8
      Options%GeoInfo%Longitude = -77.0
      Options%GeoInfo%Month = 3

      IF( Surface%Is_Allocated)THEN
         Surface%Temperature_Profile = (/1.0,3.0,4.0,5.0,9.0/)
         print*,Surface%Temperature_Profile
      ENDIF
      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("VIS_LAND") 

      CALL SfcOptics%Init(N_Angles=1)
      sfcOptics%Source_Zenith_Angle  = 30.0_fp
      sfcOptics%Source_Azimuth_Angle = 30.0_fp
  
      DO I=1,size(wvs)
         wavelength = wvs(I)
         SfcOptics%Wavenumber = 1.0E4_fp/wavelength       
         Error_Status = CSEM_Compute_LandVis_SfcOptics(          &
             Surface                        ,&  ! Input
             SfcOptics                      ,&  ! Output
             Options)                           ! Input
        print*,I,wavelength,SfcOptics%Emissivity(1,1),surface_type
       END DO
 
     
    END PROGRAM TEST_CSEM_LandVis_SfcOptics
