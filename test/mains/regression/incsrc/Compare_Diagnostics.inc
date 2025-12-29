  !
  ! Include file containing internal subprograms to report comparison diagnostics
  !
  SUBROUTINE Report_RTSolution_DiffStats_R1( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_RTSolution_Define, ONLY: CRTM_RTSolution_type, CRTM_RTSolution_Associated
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: x(:), y(:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i

    IF ( SIZE(x) /= SIZE(y) ) THEN
      Message = TRIM(Label)//' diff stats: size mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO i = 1, SIZE(x)
      IF ( .NOT. CRTM_RTSolution_Associated(x(i)) .OR. &
           .NOT. CRTM_RTSolution_Associated(y(i)) ) CYCLE

      CALL Update_Diff_Scalar( x(i)%SOD, y(i)%SOD, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Surface_Emissivity, y(i)%Surface_Emissivity, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Surface_Reflectivity, y(i)%Surface_Reflectivity, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Up_Radiance, y(i)%Up_Radiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Down_Radiance, y(i)%Down_Radiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Down_Solar_Radiance, y(i)%Down_Solar_Radiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Surface_Planck_Radiance, y(i)%Surface_Planck_Radiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Total_Cloud_Cover, y(i)%Total_Cloud_Cover, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%R_clear, y(i)%R_clear, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Tb_clear, y(i)%Tb_clear, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Reflectance_clear, y(i)%Reflectance_clear, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Radiance, y(i)%Radiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Brightness_Temperature, y(i)%Brightness_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Solar_Irradiance, y(i)%Solar_Irradiance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Reflectance, y(i)%Reflectance, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_1D( x(i)%Stokes, y(i)%Stokes, &
                           min_diff, max_diff, max_rel, has_value )

      IF ( ALLOCATED(x(i)%Upwelling_Overcast_Radiance) .AND. ALLOCATED(y(i)%Upwelling_Overcast_Radiance) ) THEN
        CALL Update_Diff_1D( x(i)%Upwelling_Overcast_Radiance, y(i)%Upwelling_Overcast_Radiance, &
                             min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Upwelling_Radiance) .AND. ALLOCATED(y(i)%Upwelling_Radiance) ) THEN
        CALL Update_Diff_1D( x(i)%Upwelling_Radiance, y(i)%Upwelling_Radiance, &
                             min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Layer_Optical_Depth) .AND. ALLOCATED(y(i)%Layer_Optical_Depth) ) THEN
        CALL Update_Diff_1D( x(i)%Layer_Optical_Depth, y(i)%Layer_Optical_Depth, &
                             min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Reflectivity) .AND. ALLOCATED(y(i)%Reflectivity) ) THEN
        CALL Update_Diff_1D( x(i)%Reflectivity, y(i)%Reflectivity, &
                             min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Reflectivity_Attenuated) .AND. ALLOCATED(y(i)%Reflectivity_Attenuated) ) THEN
        CALL Update_Diff_1D( x(i)%Reflectivity_Attenuated, y(i)%Reflectivity_Attenuated, &
                             min_diff, max_diff, max_rel, has_value )
      END IF
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_RTSolution_DiffStats_R1


  SUBROUTINE Report_RTSolution_DiffStats( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_RTSolution_Define, ONLY: CRTM_RTSolution_type, CRTM_RTSolution_Associated
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_RTSolution_type), INTENT(IN) :: x(:,:), y(:,:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i, j

    IF ( ANY(SHAPE(x) /= SHAPE(y)) ) THEN
      Message = TRIM(Label)//' diff stats: shape mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO j = 1, SIZE(x, DIM=2)
      DO i = 1, SIZE(x, DIM=1)
        IF ( .NOT. CRTM_RTSolution_Associated(x(i,j)) .OR. &
             .NOT. CRTM_RTSolution_Associated(y(i,j)) ) CYCLE

        CALL Update_Diff_Scalar( x(i,j)%SOD, y(i,j)%SOD, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Surface_Emissivity, y(i,j)%Surface_Emissivity, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Surface_Reflectivity, y(i,j)%Surface_Reflectivity, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Up_Radiance, y(i,j)%Up_Radiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Down_Radiance, y(i,j)%Down_Radiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Down_Solar_Radiance, y(i,j)%Down_Solar_Radiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Surface_Planck_Radiance, y(i,j)%Surface_Planck_Radiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Total_Cloud_Cover, y(i,j)%Total_Cloud_Cover, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%R_clear, y(i,j)%R_clear, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Tb_clear, y(i,j)%Tb_clear, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Reflectance_clear, y(i,j)%Reflectance_clear, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Radiance, y(i,j)%Radiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Brightness_Temperature, y(i,j)%Brightness_Temperature, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Solar_Irradiance, y(i,j)%Solar_Irradiance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Reflectance, y(i,j)%Reflectance, &
                                 min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_1D( x(i,j)%Stokes, y(i,j)%Stokes, &
                             min_diff, max_diff, max_rel, has_value )

        IF ( ALLOCATED(x(i,j)%Upwelling_Overcast_Radiance) .AND. ALLOCATED(y(i,j)%Upwelling_Overcast_Radiance) ) THEN
          CALL Update_Diff_1D( x(i,j)%Upwelling_Overcast_Radiance, y(i,j)%Upwelling_Overcast_Radiance, &
                               min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Upwelling_Radiance) .AND. ALLOCATED(y(i,j)%Upwelling_Radiance) ) THEN
          CALL Update_Diff_1D( x(i,j)%Upwelling_Radiance, y(i,j)%Upwelling_Radiance, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Layer_Optical_Depth) .AND. ALLOCATED(y(i,j)%Layer_Optical_Depth) ) THEN
          CALL Update_Diff_1D( x(i,j)%Layer_Optical_Depth, y(i,j)%Layer_Optical_Depth, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Reflectivity) .AND. ALLOCATED(y(i,j)%Reflectivity) ) THEN
          CALL Update_Diff_1D( x(i,j)%Reflectivity, y(i,j)%Reflectivity, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Reflectivity_Attenuated) .AND. ALLOCATED(y(i,j)%Reflectivity_Attenuated) ) THEN
          CALL Update_Diff_1D( x(i,j)%Reflectivity_Attenuated, y(i,j)%Reflectivity_Attenuated, &
                               min_diff, max_diff, max_rel, has_value )
        END IF
      END DO
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_RTSolution_DiffStats


  SUBROUTINE Report_Atmosphere_DiffStats_R1( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, CRTM_Atmosphere_Associated
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: x(:), y(:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i

    IF ( SIZE(x) /= SIZE(y) ) THEN
      Message = TRIM(Label)//' diff stats: size mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO i = 1, SIZE(x)
      IF ( .NOT. CRTM_Atmosphere_Associated(x(i)) .OR. &
           .NOT. CRTM_Atmosphere_Associated(y(i)) ) CYCLE

      IF ( ALLOCATED(x(i)%Level_Pressure) .AND. ALLOCATED(y(i)%Level_Pressure) ) THEN
        CALL Update_Diff_1D( x(i)%Level_Pressure, y(i)%Level_Pressure, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Height) .AND. ALLOCATED(y(i)%Height) ) THEN
        CALL Update_Diff_1D( x(i)%Height, y(i)%Height, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Pressure) .AND. ALLOCATED(y(i)%Pressure) ) THEN
        CALL Update_Diff_1D( x(i)%Pressure, y(i)%Pressure, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Temperature) .AND. ALLOCATED(y(i)%Temperature) ) THEN
        CALL Update_Diff_1D( x(i)%Temperature, y(i)%Temperature, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Relative_Humidity) .AND. ALLOCATED(y(i)%Relative_Humidity) ) THEN
        CALL Update_Diff_1D( x(i)%Relative_Humidity, y(i)%Relative_Humidity, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Absorber) .AND. ALLOCATED(y(i)%Absorber) ) THEN
        CALL Update_Diff_2D( x(i)%Absorber, y(i)%Absorber, min_diff, max_diff, max_rel, has_value )
      END IF
      IF ( ALLOCATED(x(i)%Cloud_Fraction) .AND. ALLOCATED(y(i)%Cloud_Fraction) ) THEN
        CALL Update_Diff_1D( x(i)%Cloud_Fraction, y(i)%Cloud_Fraction, min_diff, max_diff, max_rel, has_value )
      END IF
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_Atmosphere_DiffStats_R1


  SUBROUTINE Report_Atmosphere_DiffStats( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_Atmosphere_Define, ONLY: CRTM_Atmosphere_type, CRTM_Atmosphere_Associated
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_Atmosphere_type), INTENT(IN) :: x(:,:), y(:,:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i, j

    IF ( ANY(SHAPE(x) /= SHAPE(y)) ) THEN
      Message = TRIM(Label)//' diff stats: shape mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO j = 1, SIZE(x, DIM=2)
      DO i = 1, SIZE(x, DIM=1)
        IF ( .NOT. CRTM_Atmosphere_Associated(x(i,j)) .OR. &
             .NOT. CRTM_Atmosphere_Associated(y(i,j)) ) CYCLE

        IF ( ALLOCATED(x(i,j)%Level_Pressure) .AND. ALLOCATED(y(i,j)%Level_Pressure) ) THEN
          CALL Update_Diff_1D( x(i,j)%Level_Pressure, y(i,j)%Level_Pressure, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Height) .AND. ALLOCATED(y(i,j)%Height) ) THEN
          CALL Update_Diff_1D( x(i,j)%Height, y(i,j)%Height, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Pressure) .AND. ALLOCATED(y(i,j)%Pressure) ) THEN
          CALL Update_Diff_1D( x(i,j)%Pressure, y(i,j)%Pressure, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Temperature) .AND. ALLOCATED(y(i,j)%Temperature) ) THEN
          CALL Update_Diff_1D( x(i,j)%Temperature, y(i,j)%Temperature, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Relative_Humidity) .AND. ALLOCATED(y(i,j)%Relative_Humidity) ) THEN
          CALL Update_Diff_1D( x(i,j)%Relative_Humidity, y(i,j)%Relative_Humidity, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Absorber) .AND. ALLOCATED(y(i,j)%Absorber) ) THEN
          CALL Update_Diff_2D( x(i,j)%Absorber, y(i,j)%Absorber, min_diff, max_diff, max_rel, has_value )
        END IF
        IF ( ALLOCATED(x(i,j)%Cloud_Fraction) .AND. ALLOCATED(y(i,j)%Cloud_Fraction) ) THEN
          CALL Update_Diff_1D( x(i,j)%Cloud_Fraction, y(i,j)%Cloud_Fraction, min_diff, max_diff, max_rel, has_value )
        END IF
      END DO
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_Atmosphere_DiffStats


  SUBROUTINE Report_Surface_DiffStats_R1( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_Surface_Define, ONLY: CRTM_Surface_type
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_Surface_type), INTENT(IN) :: x(:), y(:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i

    IF ( SIZE(x) /= SIZE(y) ) THEN
      Message = TRIM(Label)//' diff stats: size mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO i = 1, SIZE(x)
      CALL Update_Diff_Scalar( x(i)%Land_Coverage, y(i)%Land_Coverage, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Water_Coverage, y(i)%Water_Coverage, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Snow_Coverage, y(i)%Snow_Coverage, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Ice_Coverage, y(i)%Ice_Coverage, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Land_Temperature, y(i)%Land_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Soil_Moisture_Content, y(i)%Soil_Moisture_Content, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Canopy_Water_Content, y(i)%Canopy_Water_Content, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Vegetation_Fraction, y(i)%Vegetation_Fraction, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Soil_Temperature, y(i)%Soil_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%LAI, y(i)%LAI, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Water_Temperature, y(i)%Water_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Wind_Speed, y(i)%Wind_Speed, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Wind_Direction, y(i)%Wind_Direction, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Salinity, y(i)%Salinity, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Snow_Temperature, y(i)%Snow_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Snow_Depth, y(i)%Snow_Depth, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Snow_Density, y(i)%Snow_Density, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Snow_Grain_Size, y(i)%Snow_Grain_Size, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Ice_Temperature, y(i)%Ice_Temperature, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Ice_Thickness, y(i)%Ice_Thickness, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Ice_Density, y(i)%Ice_Density, &
                               min_diff, max_diff, max_rel, has_value )
      CALL Update_Diff_Scalar( x(i)%Ice_Roughness, y(i)%Ice_Roughness, &
                               min_diff, max_diff, max_rel, has_value )
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_Surface_DiffStats_R1


  SUBROUTINE Report_Surface_DiffStats( Program_Name, Label, x, y )
    USE Type_Kinds, ONLY: fp
    USE Message_Handler, ONLY: INFORMATION, Display_Message
    USE CRTM_Surface_Define, ONLY: CRTM_Surface_type
    CHARACTER(*), INTENT(IN) :: Program_Name
    CHARACTER(*), INTENT(IN) :: Label
    TYPE(CRTM_Surface_type), INTENT(IN) :: x(:,:), y(:,:)
    REAL(fp) :: min_diff, max_diff, max_rel
    LOGICAL :: has_value
    CHARACTER(256) :: Message
    INTEGER :: i, j

    IF ( ANY(SHAPE(x) /= SHAPE(y)) ) THEN
      Message = TRIM(Label)//' diff stats: shape mismatch'
      CALL Display_Message( Program_Name, Message, INFORMATION )
      RETURN
    END IF

    has_value = .FALSE.

    DO j = 1, SIZE(x, DIM=2)
      DO i = 1, SIZE(x, DIM=1)
        CALL Update_Diff_Scalar( x(i,j)%Land_Coverage,        y(i,j)%Land_Coverage,        min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Water_Coverage,       y(i,j)%Water_Coverage,       min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Snow_Coverage,        y(i,j)%Snow_Coverage,        min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Ice_Coverage,         y(i,j)%Ice_Coverage,         min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Land_Temperature,     y(i,j)%Land_Temperature,     min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Soil_Moisture_Content,y(i,j)%Soil_Moisture_Content,min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Canopy_Water_Content, y(i,j)%Canopy_Water_Content, min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Vegetation_Fraction,  y(i,j)%Vegetation_Fraction,  min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Soil_Temperature,     y(i,j)%Soil_Temperature,     min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%LAI,                  y(i,j)%LAI,                  min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Water_Temperature,    y(i,j)%Water_Temperature,    min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Wind_Speed,           y(i,j)%Wind_Speed,           min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Wind_Direction,       y(i,j)%Wind_Direction,       min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Salinity,             y(i,j)%Salinity,             min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Snow_Temperature,     y(i,j)%Snow_Temperature,     min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Snow_Depth,           y(i,j)%Snow_Depth,           min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Snow_Density,         y(i,j)%Snow_Density,         min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Snow_Grain_Size,      y(i,j)%Snow_Grain_Size,      min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Ice_Temperature,      y(i,j)%Ice_Temperature,      min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Ice_Thickness,        y(i,j)%Ice_Thickness,        min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Ice_Density,          y(i,j)%Ice_Density,          min_diff, max_diff, max_rel, has_value )
        CALL Update_Diff_Scalar( x(i,j)%Ice_Roughness,        y(i,j)%Ice_Roughness,        min_diff, max_diff, max_rel, has_value )
      END DO
    END DO

    IF ( .NOT. has_value ) THEN
      Message = TRIM(Label)//' diff stats: no floating-point components compared'
    ELSE
      WRITE( Message, '(A,": min abs=",es13.6," max abs=",es13.6," max rel=",es13.6)' ) &
        TRIM(Label)//' diff stats', min_diff, max_diff, max_rel
    END IF
    CALL Display_Message( Program_Name, Message, INFORMATION )
  END SUBROUTINE Report_Surface_DiffStats


  SUBROUTINE Update_Diff_Scalar( x, y, min_diff, max_diff, max_rel, has_value )
    USE Type_Kinds, ONLY: fp
    REAL(fp), INTENT(IN) :: x
    REAL(fp), INTENT(IN) :: y
    REAL(fp), INTENT(IN OUT) :: min_diff
    REAL(fp), INTENT(IN OUT) :: max_diff
    REAL(fp), INTENT(IN OUT) :: max_rel
    LOGICAL, INTENT(IN OUT) :: has_value
    REAL(fp) :: diff
    REAL(fp) :: rel
    REAL(fp) :: denom

    diff = ABS(x - y)
    denom = MAX(ABS(y), TINY(1.0_fp))
    rel = diff / denom
    IF ( .NOT. has_value ) THEN
      min_diff = diff
      max_diff = diff
      max_rel = rel
      has_value = .TRUE.
    ELSE
      min_diff = MIN(min_diff, diff)
      max_diff = MAX(max_diff, diff)
      max_rel = MAX(max_rel, rel)
    END IF
  END SUBROUTINE Update_Diff_Scalar


  SUBROUTINE Update_Diff_1D( x, y, min_diff, max_diff, max_rel, has_value )
    USE Type_Kinds, ONLY: fp
    REAL(fp), INTENT(IN) :: x(:)
    REAL(fp), INTENT(IN) :: y(:)
    REAL(fp), INTENT(IN OUT) :: min_diff
    REAL(fp), INTENT(IN OUT) :: max_diff
    REAL(fp), INTENT(IN OUT) :: max_rel
    LOGICAL, INTENT(IN OUT) :: has_value
    REAL(fp) :: min_local
    REAL(fp) :: max_local
    REAL(fp) :: max_rel_local
    REAL(fp) :: denom_floor

    IF ( SIZE(x) == 0 ) RETURN
    min_local = MINVAL( ABS(x - y) )
    max_local = MAXVAL( ABS(x - y) )
    denom_floor = TINY(1.0_fp)
    max_rel_local = MAXVAL( ABS(x - y) / MAX(ABS(y), denom_floor) )

    IF ( .NOT. has_value ) THEN
      min_diff = min_local
      max_diff = max_local
      max_rel = max_rel_local
      has_value = .TRUE.
    ELSE
      min_diff = MIN(min_diff, min_local)
      max_diff = MAX(max_diff, max_local)
      max_rel = MAX(max_rel, max_rel_local)
    END IF
  END SUBROUTINE Update_Diff_1D


  SUBROUTINE Update_Diff_2D( x, y, min_diff, max_diff, max_rel, has_value )
    USE Type_Kinds, ONLY: fp
    REAL(fp), INTENT(IN) :: x(:,:)
    REAL(fp), INTENT(IN) :: y(:,:)
    REAL(fp), INTENT(IN OUT) :: min_diff
    REAL(fp), INTENT(IN OUT) :: max_diff
    REAL(fp), INTENT(IN OUT) :: max_rel
    LOGICAL, INTENT(IN OUT) :: has_value
    REAL(fp) :: min_local
    REAL(fp) :: max_local
    REAL(fp) :: max_rel_local
    REAL(fp) :: denom_floor

    IF ( SIZE(x) == 0 ) RETURN
    min_local = MINVAL( ABS(x - y) )
    max_local = MAXVAL( ABS(x - y) )
    denom_floor = TINY(1.0_fp)
    max_rel_local = MAXVAL( ABS(x - y) / MAX(ABS(y), denom_floor) )

    IF ( .NOT. has_value ) THEN
      min_diff = min_local
      max_diff = max_local
      max_rel = max_rel_local
      has_value = .TRUE.
    ELSE
      min_diff = MIN(min_diff, min_local)
      max_diff = MAX(max_diff, max_local)
      max_rel = MAX(max_rel, max_rel_local)
    END IF
  END SUBROUTINE Update_Diff_2D
