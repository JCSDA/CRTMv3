!
! CSEM_MWSensors_IceEmiss_Module
!
! Module to compute the surface optical properties for snow surfaces at
! microwave frequencies required for determining the snow surface
! contribution to the radiative transfer.
!
! This module is provided to allow developers to "wrap" their algorithm (Model)
! as an option into CSEM and to simplify the integration with the upper-level
! module.
!
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2014
!


MODULE NESDIS_Sensors_IceMW_Modules
  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use

  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Exception_Handler
  USE CSEM_Define
  USE NESDIS_AMSU_IceEM_Module,    ONLY: NESDIS_AMSU_ICEEM
  USE NESDIS_AMSRE_IceEM_Module,   ONLY: NESDIS_AMSRE_ICEEM
  USE NESDIS_SSMI_IceEM_Module,    ONLY: NESDIS_SSMI_IceEM
  USE NESDIS_MHS_IceEM_Module,     ONLY: NESDIS_MHS_ICEEM
  USE NESDIS_SSMIS_IceEM_Module,   ONLY: NESDIS_SSMIS_IceEM
  USE NESDIS_ATMS_IceEM_Module,    ONLY: NESDIS_ATMS_SeaICE

  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  REAL, PARAMETER :: ONE = 1.0_fp, ZERO = 0.0_fp
  REAL, PARAMETER :: PI = 3.141592653589793238462643_fp

  PUBLIC :: CRTM_Sensors_IceMW_Emiss
  PUBLIC :: IceMW_SensorName

CONTAINS


!----------------------------------------------------------------------------------
!:sdoc+:
!
! NAME:
!       CSEM_Compute_MW_Ice_SfcOptics
!
! PURPOSE:
!       Function to compute the surface emissivity and reflectivity at microwave
!       frequencies over a Ice surface.
!
!       This function is a wrapper for third party code.
!
! CALLING SEQUENCE:
!       IO_Status = CSEM_Compute_MW_Ice_SfcOptics( &
!                        Surface     , & 
!                        SensorObs,    & 
!                        SfcOptics) 
!
! INPUTS:
!       Surface:         CSEM Ice Surface structure containing the surface state
!                        data.
!                        UNITS:      N/A
!                        TYPE:       CSEM_Ice_Surface
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!       SensorObs:       CSEM input extension structure containing sensor observations,
!                        which may be used by emprirical and semi-empirical models.
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN)
!
!
! OUTPUTS:
!       SfcOptics:       CSEM output structure containing the surface
!                        optical properties required for the radiative
!                        transfer calculation. 
!                        UNITS:      N/A
!                        TYPE:       CSEM_SfcOptics_Type
!                        DIMENSION:  Scalar
!                        ATTRIBUTES: INTENT(IN OUT)
!
! FUNCTION RESULT:
!       IO_Status:    The return value is an integer defining the error status.
!                        The error codes are defined in the CSEM_Exception_Handler module.
!                        If == SUCCESS the computation was sucessful
!                           == FAILURE an unrecoverable error occurred
!                        UNITS:      N/A
!                        TYPE:       INTEGER
!                        DIMENSION:  Scalar
!
! COMMENTS:
!       Note the INTENT on the output SfcOptics argument is IN OUT rather
!       than just OUT as it is assumed to contain some data upon input.
!
!:sdoc-:
!----------------------------------------------------------------------------------

  FUNCTION CRTM_Sensors_IceMW_Emiss( &
      Surface                   ,&  ! Input
      SensorObs                 ,&  ! OptionalInput 
      SfcOptics)                 &  ! Output
    RESULT (IO_Status)
    
    TYPE(CSEM_Ice_Surface),       INTENT(IN)     :: Surface
    TYPE(CSEM_SensorObs_Struct),  INTENT(IN)     :: SensorObs
    TYPE(CSEM_SfcOptics_Type),    INTENT(INOUT)  :: SfcOptics
 
    INTEGER ::  IO_Status 
 
    CHARACTER(*), PARAMETER :: ROUTINE_NAME = 'CSEM_IceMW_Sensors_Emiss'

    !local vars
    REAL(fp),     PARAMETER :: DEFAULT_EMISSIVITY = 0.92_fp
    REAL(fp),     PARAMETER :: NOT_USED(4) = -99.9_fp
    INTEGER,      PARAMETER :: AMSRE_V_INDEX(6) = (/1, 3, 5, 7,  9, 11/) ! AMSRE channels with V pol.
    INTEGER,      PARAMETER :: AMSRE_H_INDEX(6) = (/2, 4, 6, 8, 10, 12/) ! AMSRE channels with H pol.
    INTEGER,      PARAMETER :: AMSUA_INDEX(4)   = (/1, 2, 3, 15/)
    INTEGER,      PARAMETER :: SSMIS_INDEX(8)   = (/13,12,14,16,15,17,18,8/)  ! With swapped polarisations
    INTEGER,      PARAMETER :: ATMS_INDEX(5)    = (/1, 2, 3, 16,17/)          ! With mixed polarisations

    INTEGER :: i
    CHARACTER (LEN=100) :: Sensor_Name = ' '
    
    IO_Status  = SUCCESS

   ! Compute the surface emissivities
    i= INDEX(SensorObs%Sensor_Id,'_')
    IF( i > 1)  Sensor_Name = SensorObs%Sensor_Id(1:I-1)
    
    Sensor_Type: SELECT CASE( TRIM(Sensor_Name) )
      ! ATMSemissivity model
      CASE( 'atms' )    
        DO i = 1, SfcOptics%n_Angles
          CALL NESDIS_ATMS_SeaICE( SfcOptics%Sensor_Zenith_Angle,       &  ! Input, Degree            
                                   SfcOptics%Angle(i),                  &  ! Input, Degree           
                                   SfcOptics%Frequency,                 &  ! Input, GHz                  
                                   Surface%Ice_Temperature,             &  ! Input, K
                                   SensorObs%Tb(ATMS_INDEX),            &  ! Input, ATMS           
                                   SfcOptics%Emissivity(i,2),           &  ! Output, H component      
                                   SfcOptics%Emissivity(i,1)   )               ! Output, V component 
         END DO                                                                                           

      ! AMSU-A emissivity model
      CASE( 'amsua' )
         DO i = 1,SfcOptics% N_Angles  
           CALL NESDIS_AMSU_ICEEM(  SfcOptics%Sensor_Zenith_Angle,      &  ! Input, Degree 
                                   SfcOptics%Angle(i),                  &  ! Input, Degree
                                   SfcOptics%Frequency,                 &  ! Input, GHz
                                   Surface%Ice_Temperature,             &  ! Input, K
                                   SensorObs%Tb(AMSUA_INDEX),           &  ! Input, AMSUA
                                   NOT_USED(1:2),                       &  ! Input, AMSUB  *** NO AMSU-B DATA ***
                                   SfcOptics%Emissivity(i,2),           &  ! Output, H component
                                   SfcOptics%Emissivity(i,1)            )  ! Output, V component
         ENDDO
 
      ! AMSU-B emissivity model
      CASE( 'amsub')
         DO i = 1,SfcOptics%N_Angles  
            CALL NESDIS_AMSU_ICEEM( SfcOptics%Sensor_Zenith_Angle,      &  ! Input, Degree 
                                   SfcOptics%Angle(i),                  &  ! Input, Degree
                                   SfcOptics%Frequency,                 &  ! Input, GHz
                                   Surface%Ice_Temperature,             &  ! Input, K
                                   NOT_USED,                            &  ! Input  AMSUA  *** NO AMSU-A DATA ***
                                   SensorObs%Tb(1:2),                   &  ! Input,, AMSUB
                                   SfcOptics%Emissivity(i,2),           &  ! Output, H component
                                   SfcOptics%Emissivity(i,1)            )  ! Output, V component
          ENDDO
 
      ! MHS emissivity model
      CASE ('mhs')
         DO i = 1,SfcOptics% N_Angles  
           CALL NESDIS_MHS_ICEEM( SfcOptics%Sensor_Zenith_Angle,        &  ! Input, Degree 
                                  SfcOptics%Angle(i),                   &  ! Input, Degree
                                  SfcOptics%Frequency,                  &  ! Input, GHz
                                  Surface%Ice_Temperature,              &  ! Input, K
                                  SensorObs%Tb(1:2),                    &  ! Input,, AMSUB
                                  SfcOptics%Emissivity(i,2),            &  ! Output, H component
                                  SfcOptics%Emissivity(i,1)             )  ! Output, V component
          ENDDO
 
      ! AMSR-E emissivity model
      CASE( 'amsre' )
         DO i = 1,SfcOptics%N_Angles  
           CALL NESDIS_AMSRE_ICEEM(SfcOptics%Frequency,                 &  ! Input, GHz
                                  SfcOptics%Angle(i),                   &  ! Input, Degree
                                  SensorObs%Tb(AMSRE_V_INDEX),          &  ! Input,, AMSUB
                                  SensorObs%Tb(AMSRE_H_INDEX),          &  ! Input,, AMSUB
                                  Surface%Ice_Temperature,              &  ! Input, Ts, K
                                  Surface%Ice_Temperature,              &  ! Input, Tsnow, K
                                  SfcOptics%Emissivity(i,2),            &  ! Output, H component
                                  SfcOptics%Emissivity(i,1)             )  ! Output, V component
         ENDDO

      ! SSM/I emissivity model
      CASE( 'ssmi' )
         DO i = 1,SfcOptics%N_Angles  
           CALL NESDIS_SSMI_IceEM(SfcOptics%Frequency,                  &  ! Input, GHz
                                  SfcOptics%Angle(i),                   &  ! Input, Degree
                                  Surface%Ice_Temperature,              &  ! Input, K
                                  SensorObs%Tb,                         &  ! Input, K
                                  Surface%Ice_Thickness,                &  ! Input, mm
                                  SfcOptics%Emissivity(i,2),            &  ! Output, H component
                                  SfcOptics%Emissivity(i,1)             )  ! Output, V component
         ENDDO

      ! SSMIS emissivity model
      CASE( 'ssmis' )
         DO i = 1,SfcOptics%N_Angles  
           CALL NESDIS_SSMIS_IceEM(SfcOptics%Frequency,                 &  ! Input, GHz
                                   SfcOptics%Angle(i),                  &  ! Input, Degree
                                   Surface%Ice_Temperature,             &  ! Input, K
                                   SensorObs%Tb(SSMIS_INDEX),           &  ! Input, K
                                   Surface%Ice_Thickness,               &  ! Input, mm
                                   SfcOptics%Emissivity(i,2),           &  ! Output, H component
                                   SfcOptics%Emissivity(i,1)            )  ! Output, V component

         ENDDO

      ! Default physical model
      CASE DEFAULT
           
          IO_Status  = FAILURE
          SfcOptics%Emissivity(:,1) = DEFAULT_EMISSIVITY
          SfcOptics%Emissivity(:,2) = DEFAULT_EMISSIVITY
       
 
    END SELECT Sensor_Type

    ! Compute the surface reflectivities,
    ! assuming a specular surface
    DO i = 1,SfcOptics% N_Angles  
      SfcOptics%Reflectivity(i,1,i,1) = ONE-SfcOptics%Emissivity(i,1)
      SfcOptics%Reflectivity(i,2,i,2) = ONE-SfcOptics%Emissivity(i,2)
    ENDDO
     
  END FUNCTION CRTM_Sensors_IceMW_Emiss

  FUNCTION IceMW_SensorName ( sensor_id ) RESULT ( sensor_name )
    USE CSEM_String_Utility
    ! -- Argument and result
    CHARACTER( * ), INTENT( IN )     :: sensor_id
    CHARACTER( LEN( sensor_id ) )    :: sensor_name
    
    IF(INDEX(Str2LowCase(TRIM(sensor_id)),  'amsua') > 0) THEN
      sensor_name = 'amsua'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'amsub') > 0) THEN
      sensor_name = 'amsub'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'mhs') > 0) THEN
      sensor_name = 'mhs'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'amsre') > 0) THEN
      sensor_name = 'amsre'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'atms') > 0) THEN
      sensor_name = 'atms'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'ssmis') > 0) THEN
      sensor_name = 'ssmis'
    ELSE IF(INDEX(Str2LowCase(TRIM(sensor_id)),'ssmi') > 0) THEN
      sensor_name = 'ssmi'
    ELSE
      sensor_name = 'unknown'
    ENDIF   
  END FUNCTION IceMW_SensorName
  
END MODULE NESDIS_Sensors_IceMW_Modules
