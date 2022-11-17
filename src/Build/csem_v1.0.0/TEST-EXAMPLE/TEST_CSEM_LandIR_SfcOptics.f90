   PROGRAM TEST_CSEM_LandIR_SfcOptics

      USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
      USE CSEM_Define
      USE CSEM_Model_Manager
      USE CSEM_LandIR_SfcOptics ,ONLY:  &
                              CSEM_Compute_LandIR_SfcOptics, &
                              CSEM_Compute_LandIR_SfcOptics_TL,&
                              CSEM_Compute_LandIR_SfcOptics_AD

      IMPLICIT NONE
      
      TYPE(CSEM_Land_Surface)   :: Surface
      TYPE(CSEM_SfcOptics_Type) :: SfcOptics
      TYPE(CSEM_Options_Type)   :: Options

      TYPE(CSEM_Model_ID)  :: CSEM_ALG

      INTEGER :: I, Error_Status
      
      !REAL(fp) :: wvs(10)=(/0.35_fp,0.8_fp,1.4_fp,2.6_fp,3.0_fp,4.0_fp,&
      !                      5.3_fp,10.3_fp,11.2_fp,13.3_fp/)
      REAL(fp) :: wvs(10)=(/3.9_fp,6.2_fp,&
                            6.9_fp,7.3_fp,8.6_fp,9.6_fp,10.4_fp,11.2_fp,12.4_fp,13.3_fp/)
  
      REAL(fp) :: wavelength
      INTEGER  :: surface_type = 5
       REAL(fp) :: NPOS(10), UW(10)
     
      Surface%Land_Cover_Type = surface_type
      Options%GeoInfo%Latitude=18.79
      Options%GeoInfo%Longitude=-102.73
      Options%GeoInfo%Month = 5
      !CALL Surface%Alloc_Profile(5)
      IF( Surface%Is_Allocated)THEN
         Surface%Temperature_Profile = (/1.0,3.0,4.0,5.0,9.0/)
         print*,Surface%Temperature_Profile
      ENDIF
      Error_Status = load_model_repo("alg.list")
      CSEM_ALG = inq_model_option("IR_LAND") 

      CALL SfcOptics%Init(N_Angles=1)
      DO I=1,size(wvs)
         wavelength = wvs(I)
         SfcOptics%Wavenumber = 1.0E4_fp/wavelength       
         Error_Status = CSEM_Compute_LandIR_SfcOptics(          &
             Surface                        ,&  ! Input
             SfcOptics                      ,&  ! Output
             Options)                           ! Input
        print*,I,wavelength,SfcOptics%Emissivity(1,1),surface_type
        NPOS(I)=SfcOptics%Emissivity(1,1)
      END DO
      CSEM_ALG%NAME = "RTTOV_CAMEL_ATLAS"
      CALL set_model_option(CSEM_ALG)

      DO I=1,size(wvs)
         wavelength = wvs(I)
         SfcOptics%Wavenumber = 1.0E4_fp/wavelength      
         Error_Status = CSEM_Compute_LandIR_SfcOptics(          &
             Surface                        ,&  ! Input
             SfcOptics                      ,&  ! Output
             Options)                           ! Input
         UW(I)=SfcOptics%Emissivity(1,1)
        print*,I,wavelength,SfcOptics%Emissivity(1,1),surface_type
       
      END DO
       DO I=1,size(wvs)
       WRITE(*, '(F5.1, 2F10.3)')wvs(I), NPOS(i),UW(I)
       ENDDO
       
    END PROGRAM TEST_CSEM_LandIR_SfcOptics
