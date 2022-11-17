!----------------------------------------------------------------------------------
!:sdoc+:
!
!
! NPOESS_LUT_READER
!
! Module containing the load/destruction routines to handel
! the shared NPOESS LUT. NPOESS LUT is a type-based reflectance
! spectral LUT. The LUT is composed of the reflectivity values at 74 discrete &
! wavenumber points for each of the 24 NPOESS(UMD) surface types. 
! The discrete spectrum ranges from 0.2um to 15.0um. Interpolated value if assigned
! if the wavenumber is not of the 74 points.
!
!
!
! CREATION HISTORY:
!       Written by:     Ming  Chen, 12-Nov-2015
!                       ming.chen@noaa.gov
!:sdoc-:
!-----------------------------------------------------------------------------------------------
!
MODULE NPOESS_LUT_READER
 
  !
  ! History:
  ! Version   Date     Comment
  ! -------   ----     -------
  !  
  ! Code Description:
  !   Language:           Fortran 90.
  !
  ! Declarations:
  ! Modules used:
  !
  
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  USE CSEM_Interpolation, ONLY: NPTS        , &
                                CSEM_LPoly_type  , &
                                CSEM_find_index  , &
                                CSEM_interp_1D   , &
                                CSEM_Clear_LPoly , &
                                CSEM_LPoly       

  IMPLICIT NONE
  ! ------------
  ! Visibilities
  ! ------------
  ! Everything private by default
  PRIVATE
  INTEGER, PARAMETER :: N_NPOESS_TYPES = 24
  INTEGER, PARAMETER :: N_NPOESS_WAVELENGTHS = 74

  INTEGER, PARAMETER :: NPOESS_COMPACTED_SOIL           =  1
  INTEGER, PARAMETER :: NPOESS_TILLED_SOIL              =  2
  INTEGER, PARAMETER :: NPOESS_SAND                     =  3
  INTEGER, PARAMETER :: NPOESS_ROCK                     =  4
  INTEGER, PARAMETER :: NPOESS_IRRIGATED_LOW_VEGETATION =  5
  INTEGER, PARAMETER :: NPOESS_MEADOW_GRASS             =  6
  INTEGER, PARAMETER :: NPOESS_SCRUB                    =  7
  INTEGER, PARAMETER :: NPOESS_BROADLEAF_FOREST         =  8
  INTEGER, PARAMETER :: NPOESS_PINE_FOREST              =  9
  INTEGER, PARAMETER :: NPOESS_TUNDRA                   = 10
  INTEGER, PARAMETER :: NPOESS_GRASS_SOIL               = 11
  INTEGER, PARAMETER :: NPOESS_BROADLEAF_PINE_FOREST    = 12
  INTEGER, PARAMETER :: NPOESS_GRASS_SCRUB              = 13
  INTEGER, PARAMETER :: NPOESS_SOIL_GRASS_SCRUB         = 14
  INTEGER, PARAMETER :: NPOESS_URBAN_CONCRETE           = 15
  INTEGER, PARAMETER :: NPOESS_PINE_BRUSH               = 16
  INTEGER, PARAMETER :: NPOESS_BROADLEAF_BRUSH          = 17
  INTEGER, PARAMETER :: NPOESS_WET_SOIL                 = 18
  INTEGER, PARAMETER :: NPOESS_SCRUB_SOIL               = 19
  INTEGER, PARAMETER :: NPOESS_BROADLEAF70_PINE30       = 20
  INTEGER, PARAMETER :: NPOESS_WATER                    = 21
  INTEGER, PARAMETER :: NPOESS_OLD_SNOW                 = 22
  INTEGER, PARAMETER :: NPOESS_FRESH_SNOW               = 23
  INTEGER, PARAMETER :: NPOESS_NEW_ICE                  = 24

  REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS), SAVE :: npoess_lut_wavenumbers
  REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS,N_NPOESS_TYPES), SAVE :: npoess_lut_reflectivity
  LOGICAL, SAVE :: NPOESS_LUT_INIT = .TRUE.
  
  ! --------------------------------------
  ! Structure definitions to hold forward
  ! variables across FWD, TL, and AD calls
  ! --------------------------------------
  TYPE :: iVar_type
    PRIVATE
    ! The interpolating polynomials
    TYPE(CSEM_LPoly_type) :: xlp
    ! The LUT interpolation indices
    INTEGER :: i1, i2
    ! The LUT interpolation boundary check
    LOGICAL :: x_outbound
    ! The interpolation input
    REAL(fp) :: x_int
    ! The data to be interpolated
    REAL(fp) :: x(NPTS)
  END TYPE iVar_type
  
  ! this is for the precision testing of using wavelength as input
  REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS,N_NPOESS_TYPES), SAVE :: npoess_lut_reflectivity_0
  PUBLIC :: Read_NPOESS_LUT_0

  PUBLIC :: Read_NPOESS_LUT
  PUBLIC :: Load_NPOESS_LUT
  PUBLIC :: Read_Stype_Map
  
CONTAINS

  FUNCTION Read_NPOESS_LUT( &
      wavenumber,        & ! INPUT, wavenumer in cm-1
      emissivity,        & ! OUTPUT, surface emissivity (0 - 1)
      surface_type)      & ! OPTIONAL INPUT, surface type (1 - 24)
    RESULT (Error_Status)


    INTEGER , INTENT( IN )  :: surface_type
    REAL(fp), INTENT( IN )  :: wavenumber
    REAL(fp), INTENT( OUT ) :: emissivity
    INTEGER :: Error_Status 
    TYPE(iVar_type) :: iVar
 
    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !

    REAL(fp) :: reflectivity
    
    IF (NPOESS_LUT_INIT) THEN
        CALL Load_NPOESS_LUT()
        NPOESS_LUT_INIT = .FALSE.
    ENDIF
    
    iVar%x_int = MAX(MIN(npoess_lut_wavenumbers(N_NPOESS_WAVELENGTHS),&
                 wavenumber), npoess_lut_wavenumbers(1))
    CALL CSEM_find_index(npoess_lut_wavenumbers, &
                 iVar%x_int, iVar%i1, iVar%i2, iVar%x_outbound)
    iVar%x = npoess_lut_wavenumbers(iVar%i1:iVar%i2)

    ! Calculate the interpolating polynomial
    CALL CSEM_LPoly( iVar%x, iVar%x_int, & ! Input
                iVar%xlp            )      ! Output


    ! Perform Interpolation
    CALL CSEM_interp_1D( &
          npoess_lut_reflectivity(iVar%i1:iVar%i2, Surface_Type), &
          iVar%xlp, reflectivity )
    Emissivity = 1.0_fp - reflectivity
    Error_Status = 0
  END FUNCTION Read_NPOESS_LUT


  FUNCTION Read_NPOESS_LUT_0( &
      wavelength,        & ! INPUT, wavelength in micrometer
      emissivity,        & ! OUTPUT, surface emissivity (0 - 1)
      surface_type)      & ! OPTIONAL INPUT, surface type (1 - 24)
    RESULT (Error_Status)


    INTEGER , INTENT( IN )  :: surface_type
    REAL(fp), INTENT( IN )  :: wavelength
    REAL(fp), INTENT( OUT ) :: emissivity
    INTEGER :: Error_Status 
    TYPE(iVar_type) :: iVar
 
    ! --------------------------------------------------  !
    !       internal variables                            !
    ! --------------------------------------------------  !

    REAL(fp) :: reflectivity
    REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS):: wavelengthS = (/ &
    0.200_fp, 0.225_fp, 0.250_fp, 0.275_fp, 0.300_fp, 0.325_fp, 0.350_fp, 0.375_fp, &
    0.400_fp, 0.425_fp, 0.450_fp, 0.475_fp, 0.500_fp, 0.525_fp, 0.550_fp, 0.575_fp, &
    0.600_fp, 0.625_fp, 0.650_fp, 0.675_fp, 0.700_fp, 0.725_fp, 0.750_fp, 0.775_fp, &
    0.800_fp, 0.825_fp, 0.850_fp, 0.875_fp, 0.900_fp, 0.925_fp, 0.950_fp, 0.975_fp, &
    1.000_fp, 1.050_fp, 1.100_fp, 1.150_fp, 1.200_fp, 1.250_fp, 1.300_fp, 1.350_fp, &
    1.400_fp, 1.450_fp, 1.500_fp, 1.550_fp, 1.600_fp, 1.650_fp, 1.700_fp, 1.750_fp, &
    1.800_fp, 1.850_fp, 1.900_fp, 1.950_fp, 2.000_fp, 2.500_fp, 3.000_fp, 3.500_fp, &
    4.000_fp, 4.500_fp, 5.000_fp, 5.500_fp, 6.000_fp, 6.500_fp, 7.000_fp, 7.500_fp, &
    8.000_fp, 8.500_fp, 9.000_fp, 9.500_fp,10.000_fp,11.000_fp,12.000_fp,13.000_fp, &
    14.000_fp,15.000_fp /)
    
    IF (NPOESS_LUT_INIT) THEN
        CALL Load_NPOESS_LUT()
        NPOESS_LUT_INIT = .FALSE.
    ENDIF
    
    iVar%x_int = MAX(MIN(WAVELENGTHS(N_NPOESS_WAVELENGTHS),&
                 wavelength), wavelengths(1))
    CALL CSEM_find_index(WAVELENGTHS, &
                 iVar%x_int, iVar%i1, iVar%i2, iVar%x_outbound)
    iVar%x = WAVELENGTHS(iVar%i1:iVar%i2)

    ! Calculate the interpolating polynomial
    CALL CSEM_LPoly( iVar%x, iVar%x_int, & ! Input
                iVar%xlp            )      ! Output


    ! Perform Interpolation
    CALL CSEM_interp_1D( &
          npoess_lut_reflectivity_0(iVar%i1:iVar%i2, Surface_Type), &
          iVar%xlp, reflectivity )
    Emissivity = 1.0_fp - reflectivity
    Error_Status = 0
  END FUNCTION Read_NPOESS_LUT_0


  SUBROUTINE Read_Stype_Map(alat,alon,stype)


    REAL(fp), INTENT( IN ) :: alat, alon
    REAL(fp) :: alon1
    INTEGER, PARAMETER :: Nlat = 1080, NLon = 2160
    INTEGER i,j,k,stype,init_topo,index_lat,index_lon,FileID
    INTEGER(KIND=1), DIMENSION( NLon, Nlat) :: Surface_Type
    INTEGER(KIND=1) :: I3      ! One byte
    INTEGER(KIND=2) :: I2         ! Two bytes
    DATA init_topo/0/
    SAVE init_topo, Surface_Type

    IF (init_topo .EQ. 0) THEN 
    ! data contains surface topography parameters
    !  FileID = Get_Lun()
      FileID = 16
      PRINT *,' READ TOPOGRAPHY DATA '
      OPEN(FileID,file='topography.bin.Big_Endian',form='unformatted', &
      access='direct',RECL=4)
      DO i = 1, Nlat
        DO j = 1, Nlon
          k = (i-1) * Nlon + j
          READ(FileID,rec=k) Surface_Type(j,i), I3, I2
        ENDDO
      ! alt=I2*1.0
      ! scover = I3 * 0.01
      ENDDO
      CLOSE(FileID)
      init_topo=1
    ENDIF 

    IF (alon .GT. 180.0_fp) THEN
      alon1 = alon - 360.0_fp
    ELSE
      alon1=alon
    ENDIF
    index_lat = INT((90.0_fp - alat))*6+1
    index_lon = INT((180_fp + alon1))*6+1
    IF (index_lat .GT. Nlat) index_lat=Nlat
    IF (index_lon .GT. Nlon) index_lon=Nlon
    stype = Surface_Type(index_lon, index_lat)
 
    RETURN
  END SUBROUTINE Read_Stype_Map

   
  SUBROUTINE Load_NPOESS_LUT()

    REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS):: wavelength
    REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS):: wavenumber
    REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS,N_NPOESS_TYPES):: reflectivity
    REAL(fp), DIMENSION(N_NPOESS_WAVELENGTHS,N_NPOESS_TYPES):: reflectance
    INTEGER :: i
    
    ! ----------------------------------------------------------------------------------------
    ! The 24 surface types are: original order
    !    1.     water             2.  old snow                  3.   fresh snow
    !    4.  compacted soil       5.  tilled soil               6.   sand
    !    7.     rock              8.  irrigated low vegetation  9.   meadow grass
    !   10.    scrub             11.  broadleaf forest         12.   pine forest
    !   13.    tundra            14.  grass soil               15.   broadleaf pine forest
    !   16.    grass scrub       17.  oil grass                18.   urban concrete  
    !   19.    pine brush        20.  broadleaf brush          21.   wet soil
    !   22.    scrub soil        23.  broadleaf 70-pine 30     24.   new ice
    ! ----------------------------------------------------------------------------------------
    !
    
    ! ----------------------------------------------------------------------------------------
    ! The 24 surface types are: 1-20(Land) 21 (Water) 22-23(Snow) 24 Ice
    !    1.  compacted soil       2.  tilled soil               3.   sand
    !    4.    rock               5.  irrigated low vegetation  6.   meadow grass
    !    7.    scrub              8.  broadleaf forest          9.   pine forest
    !   10.    tundra            11.  grass soil               12.   broadleaf pine forest
    !   13.    grass scrub       14.  oil grass                15.   urban concrete  
    !   16.    pine brush        17.  broadleaf brush          18.   wet soil
    !   19.    scrub soil        20.  broadleaf 70-pine 30     21.   Water
    !   22.    old snow          23.  fresh snow               24.    new ice 
    ! ----------------------------------------------------------------------------------------
    !

 
   ! wavelength um

    DATA (WAVELENGTH(i), i=1,N_NPOESS_WAVELENGTHS)/ &
    0.200_fp, 0.225_fp, 0.250_fp, 0.275_fp, 0.300_fp, 0.325_fp, 0.350_fp, 0.375_fp, &
    0.400_fp, 0.425_fp, 0.450_fp, 0.475_fp, 0.500_fp, 0.525_fp, 0.550_fp, 0.575_fp, &
    0.600_fp, 0.625_fp, 0.650_fp, 0.675_fp, 0.700_fp, 0.725_fp, 0.750_fp, 0.775_fp, &
    0.800_fp, 0.825_fp, 0.850_fp, 0.875_fp, 0.900_fp, 0.925_fp, 0.950_fp, 0.975_fp, &
    1.000_fp, 1.050_fp, 1.100_fp, 1.150_fp, 1.200_fp, 1.250_fp, 1.300_fp, 1.350_fp, &
    1.400_fp, 1.450_fp, 1.500_fp, 1.550_fp, 1.600_fp, 1.650_fp, 1.700_fp, 1.750_fp, &
    1.800_fp, 1.850_fp, 1.900_fp, 1.950_fp, 2.000_fp, 2.500_fp, 3.000_fp, 3.500_fp, &
    4.000_fp, 4.500_fp, 5.000_fp, 5.500_fp, 6.000_fp, 6.500_fp, 7.000_fp, 7.500_fp, &
    8.000_fp, 8.500_fp, 9.000_fp, 9.500_fp,10.000_fp,11.000_fp,12.000_fp,13.000_fp, &
    14.000_fp,15.000_fp /

    ! ...The reflectivity data
    DATA (reflectance(i,NPOESS_WATER),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.033_fp, 0.033_fp, 0.033_fp, 0.033_fp, 0.027_fp, 0.027_fp, 0.027_fp, 0.027_fp, &
    0.046_fp, 0.050_fp, 0.052_fp, 0.045_fp, 0.025_fp, 0.024_fp, 0.024_fp, 0.024_fp, &
    0.024_fp, 0.005_fp, 0.005_fp, 0.005_fp, 0.005_fp, 0.006_fp, 0.007_fp, 0.007_fp, &
    0.007_fp, 0.007_fp, 0.007_fp, 0.007_fp, 0.007_fp, 0.024_fp, 0.024_fp, 0.024_fp, &
    0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, &
    0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, &
    0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.022_fp, 0.001_fp, 0.001_fp, 0.020_fp, &
    0.020_fp, 0.020_fp, 0.020_fp, 0.010_fp, 0.010_fp, 0.020_fp, 0.020_fp, 0.020_fp, &
    0.020_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.020_fp, &
    0.030_fp, 0.030_fp/
    ! ...old snow
    DATA (reflectance(i,NPOESS_OLD_SNOW),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.882_fp, 0.892_fp, 0.903_fp, 0.913_fp, 0.924_fp, 0.934_fp, 0.946_fp, 0.955_fp,  &
    0.958_fp, 0.959_fp, 0.954_fp, 0.951_fp, 0.944_fp, 0.940_fp, 0.930_fp, 0.917_fp,  &
    0.908_fp, 0.900_fp, 0.889_fp, 0.876_fp, 0.856_fp, 0.833_fp, 0.807_fp, 0.779_fp,  &
    0.754_fp, 0.727_fp, 0.677_fp, 0.642_fp, 0.583_fp, 0.537_fp, 0.473_fp, 0.419_fp,  &
    0.351_fp, 0.367_fp, 0.360_fp, 0.287_fp, 0.164_fp, 0.124_fp, 0.132_fp, 0.125_fp,  &
    0.094_fp, 0.021_fp, 0.002_fp, 0.000_fp, 0.000_fp, 0.000_fp, 0.008_fp, 0.009_fp,  &
    0.015_fp, 0.037_fp, 0.025_fp, 0.025_fp, 0.000_fp, 0.000_fp, 0.031_fp, 0.013_fp,  &
    0.015_fp, 0.010_fp, 0.014_fp, 0.010_fp, 0.008_fp, 0.009_fp, 0.009_fp, 0.009_fp,  &
    0.009_fp, 0.009_fp, 0.008_fp, 0.008_fp, 0.006_fp, 0.010_fp, 0.018_fp, 0.020_fp,  &
    0.025_fp, 0.020_fp/
    ! ...fresh snow
    DATA (reflectance(i,NPOESS_FRESH_SNOW),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.961_fp, 0.964_fp, 0.966_fp, 0.969_fp, 0.972_fp, 0.974_fp, 0.977_fp, 0.979_fp,  &
    0.979_fp, 0.979_fp, 0.979_fp, 0.979_fp, 0.979_fp, 0.977_fp, 0.977_fp, 0.974_fp,  &
    0.972_fp, 0.970_fp, 0.969_fp, 0.969_fp, 0.964_fp, 0.960_fp, 0.953_fp, 0.946_fp,  &
    0.937_fp, 0.923_fp, 0.912_fp, 0.893_fp, 0.874_fp, 0.858_fp, 0.840_fp, 0.813_fp,  &
    0.790_fp, 0.779_fp, 0.786_fp, 0.758_fp, 0.654_fp, 0.611_fp, 0.612_fp, 0.607_fp,  &
    0.516_fp, 0.291_fp, 0.145_fp, 0.115_fp, 0.158_fp, 0.183_fp, 0.251_fp, 0.291_fp,  &
    0.306_fp, 0.425_fp, 0.150_fp, 0.070_fp, 0.044_fp, 0.092_fp, 0.031_fp, 0.013_fp,  &
    0.023_fp, 0.009_fp, 0.015_fp, 0.008_fp, 0.007_fp, 0.008_fp, 0.008_fp, 0.009_fp,  &
    0.008_fp, 0.009_fp, 0.006_fp, 0.006_fp, 0.005_fp, 0.008_fp, 0.018_fp, 0.020_fp,  &
    0.025_fp, 0.020_fp/
    ! ...new ice
    DATA (reflectance(i,NPOESS_NEW_ICE),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.032_fp, 0.030_fp, 0.029_fp, 0.027_fp, 0.025_fp, 0.025_fp, 0.024_fp, 0.023_fp,  &
    0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp, 0.023_fp,  &
    0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp,  &
    0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp, 0.022_fp,  &
    0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp,  &
    0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.020_fp, 0.020_fp, 0.020_fp, 0.020_fp,  &
    0.019_fp, 0.019_fp, 0.018_fp, 0.018_fp, 0.018_fp, 0.010_fp, 0.050_fp, 0.030_fp,  &
    0.020_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp,  &
    0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.020_fp, 0.040_fp, 0.050_fp,  &
    0.050_fp, 0.050_fp/
 
    ! ...compacted_soil
    DATA (reflectance(i,NPOESS_COMPACTED_SOIL),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, &
    0.024_fp, 0.025_fp, 0.030_fp, 0.037_fp, 0.050_fp, 0.089_fp, 0.102_fp, 0.122_fp, &
    0.141_fp, 0.158_fp, 0.174_fp, 0.198_fp, 0.206_fp, 0.200_fp, 0.219_fp, 0.237_fp, &
    0.248_fp, 0.256_fp, 0.263_fp, 0.264_fp, 0.273_fp, 0.190_fp, 0.200_fp, 0.205_fp, &
    0.205_fp, 0.210_fp, 0.225_fp, 0.225_fp, 0.232_fp, 0.237_fp, 0.237_fp, 0.225_fp, &
    0.187_fp, 0.149_fp, 0.162_fp, 0.190_fp, 0.212_fp, 0.225_fp, 0.225_fp, 0.212_fp, &
    0.205_fp, 0.187_fp, 0.095_fp, 0.062_fp, 0.075_fp, 0.050_fp, 0.040_fp, 0.070_fp, &
    0.140_fp, 0.125_fp, 0.118_fp, 0.100_fp, 0.060_fp, 0.050_fp, 0.040_fp, 0.030_fp, &
    0.030_fp, 0.040_fp, 0.050_fp, 0.050_fp, 0.045_fp, 0.040_fp, 0.025_fp, 0.030_fp, &
    0.020_fp, 0.020_fp/
    ! ...tilled_soil
    DATA (reflectance(i,NPOESS_TILLED_SOIL),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp,  &
    0.011_fp, 0.012_fp, 0.013_fp, 0.014_fp, 0.016_fp, 0.023_fp, 0.025_fp, 0.025_fp,  &
    0.028_fp, 0.031_fp, 0.034_fp, 0.035_fp, 0.035_fp, 0.035_fp, 0.034_fp, 0.034_fp,  &
    0.038_fp, 0.042_fp, 0.045_fp, 0.049_fp, 0.051_fp, 0.040_fp, 0.043_fp, 0.046_fp,  &
    0.049_fp, 0.051_fp, 0.054_fp, 0.057_fp, 0.065_fp, 0.070_fp, 0.076_fp, 0.076_fp,  &
    0.074_fp, 0.074_fp, 0.075_fp, 0.090_fp, 0.098_fp, 0.100_fp, 0.101_fp, 0.101_fp,  &
    0.100_fp, 0.100_fp, 0.050_fp, 0.037_fp, 0.040_fp, 0.035_fp, 0.028_fp, 0.049_fp,  &
    0.098_fp, 0.088_fp, 0.083_fp, 0.070_fp, 0.042_fp, 0.035_fp, 0.028_fp, 0.021_fp,  &
    0.021_fp, 0.028_fp, 0.035_fp, 0.035_fp, 0.031_fp, 0.028_fp, 0.018_fp, 0.030_fp,  &
    0.020_fp, 0.020_fp/
    ! ...sand
    DATA (reflectance(i,NPOESS_SAND),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.200_fp, 0.210_fp, 0.220_fp, 0.230_fp, 0.250_fp, 0.260_fp, 0.270_fp, 0.287_fp,  &
    0.300_fp, 0.312_fp, 0.325_fp, 0.350_fp, 0.375_fp, 0.387_fp, 0.400_fp, 0.423_fp,  &
    0.450_fp, 0.455_fp, 0.460_fp, 0.470_fp, 0.475_fp, 0.483_fp, 0.487_fp, 0.495_fp,  &
    0.500_fp, 0.501_fp, 0.505_fp, 0.508_fp, 0.510_fp, 0.515_fp, 0.517_fp, 0.520_fp,  &
    0.525_fp, 0.535_fp, 0.545_fp, 0.555_fp, 0.565_fp, 0.575_fp, 0.585_fp, 0.595_fp,  &
    0.605_fp, 0.615_fp, 0.625_fp, 0.630_fp, 0.636_fp, 0.641_fp, 0.647_fp, 0.652_fp,  &
    0.658_fp, 0.663_fp, 0.669_fp, 0.674_fp, 0.680_fp, 0.325_fp, 0.315_fp, 0.500_fp,  &
    0.400_fp, 0.150_fp, 0.050_fp, 0.050_fp, 0.150_fp, 0.100_fp, 0.100_fp, 0.100_fp,  &
    0.100_fp, 0.100_fp, 0.100_fp, 0.090_fp, 0.080_fp, 0.020_fp, 0.020_fp, 0.020_fp,  &
    0.020_fp, 0.020_fp/
    ! ...rock
    DATA (reflectance(i,NPOESS_ROCK),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.025_fp, 0.028_fp, 0.037_fp, 0.049_fp, 0.072_fp, 0.078_fp, 0.088_fp, 0.100_fp,  &
    0.118_fp, 0.132_fp, 0.151_fp, 0.151_fp, 0.150_fp, 0.150_fp, 0.152_fp, 0.170_fp,  &
    0.176_fp, 0.175_fp, 0.175_fp, 0.175_fp, 0.177_fp, 0.199_fp, 0.218_fp, 0.215_fp,  &
    0.208_fp, 0.202_fp, 0.210_fp, 0.212_fp, 0.213_fp, 0.217_fp, 0.225_fp, 0.234_fp,  &
    0.240_fp, 0.265_fp, 0.285_fp, 0.279_fp, 0.270_fp, 0.265_fp, 0.252_fp, 0.247_fp,  &
    0.228_fp, 0.227_fp, 0.227_fp, 0.227_fp, 0.226_fp, 0.225_fp, 0.225_fp, 0.225_fp,  &
    0.225_fp, 0.220_fp, 0.212_fp, 0.205_fp, 0.200_fp, 0.050_fp, 0.100_fp, 0.200_fp,  &
    0.120_fp, 0.180_fp, 0.070_fp, 0.050_fp, 0.070_fp, 0.080_fp, 0.050_fp, 0.040_fp,  &
    0.100_fp, 0.110_fp, 0.130_fp, 0.140_fp, 0.120_fp, 0.050_fp, 0.030_fp, 0.020_fp,  &
    0.020_fp, 0.020_fp/
    ! ...irrigated_low_vegetation
    DATA (reflectance(i,NPOESS_IRRIGATED_LOW_VEGETATION),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.025_fp, 0.029_fp, 0.034_fp,  &
    0.038_fp, 0.040_fp, 0.046_fp, 0.050_fp, 0.048_fp, 0.045_fp, 0.046_fp, 0.048_fp,  &
    0.045_fp, 0.034_fp, 0.040_fp, 0.061_fp, 0.080_fp, 0.268_fp, 0.338_fp, 0.429_fp,  &
    0.478_fp, 0.510_fp, 0.511_fp, 0.577_fp, 0.598_fp, 0.654_fp, 0.659_fp, 0.662_fp,  &
    0.663_fp, 0.665_fp, 0.677_fp, 0.593_fp, 0.570_fp, 0.546_fp, 0.456_fp, 0.366_fp,  &
    0.275_fp, 0.185_fp, 0.100_fp, 0.138_fp, 0.180_fp, 0.223_fp, 0.266_fp, 0.308_fp,  &
    0.252_fp, 0.196_fp, 0.141_fp, 0.085_fp, 0.029_fp, 0.018_fp, 0.034_fp, 0.038_fp,  &
    0.043_fp, 0.039_fp, 0.034_fp, 0.032_fp, 0.027_fp, 0.034_fp, 0.036_fp, 0.036_fp,  &
    0.036_fp, 0.036_fp, 0.036_fp, 0.036_fp, 0.036_fp, 0.036_fp, 0.045_fp, 0.036_fp,  &
    0.018_fp, 0.018_fp/
    ! ...meadow_grass
    DATA (reflectance(i,NPOESS_MEADOW_GRASS),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.025_fp, 0.025_fp, 0.025_fp, 0.025_fp, 0.025_fp, 0.025_fp, 0.025_fp, 0.030_fp,  &
    0.035_fp, 0.040_fp, 0.065_fp, 0.065_fp, 0.065_fp, 0.085_fp, 0.100_fp, 0.085_fp,  &
    0.075_fp, 0.070_fp, 0.065_fp, 0.140_fp, 0.225_fp, 0.325_fp, 0.500_fp, 0.600_fp,  &
    0.635_fp, 0.640_fp, 0.641_fp, 0.645_fp, 0.645_fp, 0.648_fp, 0.650_fp, 0.660_fp,  &
    0.675_fp, 0.670_fp, 0.665_fp, 0.660_fp, 0.655_fp, 0.650_fp, 0.546_fp, 0.442_fp,  &
    0.338_fp, 0.234_fp, 0.130_fp, 0.169_fp, 0.208_fp, 0.247_fp, 0.286_fp, 0.325_fp,  &
    0.272_fp, 0.219_fp, 0.166_fp, 0.113_fp, 0.060_fp, 0.150_fp, 0.050_fp, 0.125_fp,  &
    0.200_fp, 0.260_fp, 0.285_fp, 0.295_fp, 0.060_fp, 0.105_fp, 0.060_fp, 0.045_fp,  &
    0.050_fp, 0.060_fp, 0.068_fp, 0.075_fp, 0.100_fp, 0.165_fp, 0.150_fp, 0.125_fp,  &
    0.105_fp, 0.090_fp/
    ! ...scrub
    DATA (reflectance(i,NPOESS_SCRUB),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.012_fp, 0.012_fp, 0.012_fp, 0.012_fp, 0.012_fp, 0.012_fp, 0.012_fp, 0.012_fp,  &
    0.015_fp, 0.025_fp, 0.045_fp, 0.050_fp, 0.067_fp, 0.040_fp, 0.054_fp, 0.052_fp,  &
    0.043_fp, 0.030_fp, 0.031_fp, 0.033_fp, 0.035_fp, 0.058_fp, 0.073_fp, 0.114_fp,  &
    0.166_fp, 0.180_fp, 0.190_fp, 0.199_fp, 0.210_fp, 0.355_fp, 0.350_fp, 0.348_fp,  &
    0.343_fp, 0.334_fp, 0.326_fp, 0.317_fp, 0.308_fp, 0.300_fp, 0.269_fp, 0.238_fp,  &
    0.207_fp, 0.176_fp, 0.145_fp, 0.163_fp, 0.181_fp, 0.200_fp, 0.180_fp, 0.160_fp,  &
    0.147_fp, 0.135_fp, 0.123_fp, 0.110_fp, 0.098_fp, 0.042_fp, 0.050_fp, 0.065_fp,  &
    0.100_fp, 0.150_fp, 0.130_fp, 0.120_fp, 0.030_fp, 0.040_fp, 0.060_fp, 0.055_fp,  &
    0.050_fp, 0.040_fp, 0.030_fp, 0.035_fp, 0.040_fp, 0.050_fp, 0.060_fp, 0.050_fp,  &
    0.040_fp, 0.070_fp/
    ! ...broadleaf_forest
    DATA (reflectance(i,NPOESS_BROADLEAF_FOREST),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.030_fp, 0.030_fp, 0.032_fp, 0.030_fp, 0.030_fp, 0.029_fp, 0.027_fp, 0.034_fp,  &
    0.035_fp, 0.034_fp, 0.036_fp, 0.038_fp, 0.038_fp, 0.023_fp, 0.056_fp, 0.042_fp,  &
    0.035_fp, 0.012_fp, 0.014_fp, 0.013_fp, 0.013_fp, 0.074_fp, 0.255_fp, 0.337_fp,  &
    0.318_fp, 0.314_fp, 0.316_fp, 0.315_fp, 0.314_fp, 0.465_fp, 0.456_fp, 0.453_fp,  &
    0.451_fp, 0.448_fp, 0.446_fp, 0.445_fp, 0.444_fp, 0.441_fp, 0.391_fp, 0.341_fp,  &
    0.290_fp, 0.240_fp, 0.190_fp, 0.222_fp, 0.255_fp, 0.287_fp, 0.319_fp, 0.351_fp,  &
    0.288_fp, 0.226_fp, 0.163_fp, 0.101_fp, 0.038_fp, 0.020_fp, 0.038_fp, 0.048_fp,  &
    0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp,  &
    0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.038_fp, 0.050_fp, 0.038_fp,  &
    0.019_fp, 0.019_fp/
    ! ...pine_forest
    DATA (reflectance(i,NPOESS_PINE_FOREST),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.012_fp, 0.013_fp, 0.014_fp, 0.015_fp, 0.016_fp, 0.017_fp, 0.018_fp, 0.019_fp,  &
    0.020_fp, 0.023_fp, 0.025_fp, 0.027_fp, 0.030_fp, 0.035_fp, 0.040_fp, 0.037_fp,  &
    0.035_fp, 0.033_fp, 0.031_fp, 0.030_fp, 0.075_fp, 0.150_fp, 0.250_fp, 0.260_fp,  &
    0.270_fp, 0.280_fp, 0.290_fp, 0.295_fp, 0.300_fp, 0.310_fp, 0.305_fp, 0.300_fp,  &
    0.280_fp, 0.295_fp, 0.300_fp, 0.270_fp, 0.230_fp, 0.250_fp, 0.210_fp, 0.160_fp,  &
    0.110_fp, 0.075_fp, 0.095_fp, 0.115_fp, 0.125_fp, 0.130_fp, 0.140_fp, 0.150_fp,  &
    0.120_fp, 0.090_fp, 0.070_fp, 0.090_fp, 0.105_fp, 0.120_fp, 0.100_fp, 0.095_fp,  &
    0.090_fp, 0.070_fp, 0.050_fp, 0.030_fp, 0.025_fp, 0.022_fp, 0.015_fp, 0.012_fp,  &
    0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp, 0.010_fp,  &
    0.010_fp, 0.010_fp/
    ! ...tundra
    DATA (reflectance(i,NPOESS_TUNDRA),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.032_fp, 0.032_fp, 0.034_fp, 0.037_fp, 0.039_fp, 0.041_fp, 0.042_fp, 0.046_fp,  &
    0.058_fp, 0.063_fp, 0.075_fp, 0.073_fp, 0.070_fp, 0.075_fp, 0.081_fp, 0.083_fp,  &
    0.084_fp, 0.084_fp, 0.086_fp, 0.103_fp, 0.123_fp, 0.154_fp, 0.201_fp, 0.228_fp,  &
    0.242_fp, 0.245_fp, 0.250_fp, 0.254_fp, 0.257_fp, 0.256_fp, 0.257_fp, 0.258_fp,  &
    0.270_fp, 0.274_fp, 0.279_fp, 0.276_fp, 0.274_fp, 0.273_fp, 0.244_fp, 0.216_fp,  &
    0.180_fp, 0.150_fp, 0.128_fp, 0.143_fp, 0.156_fp, 0.167_fp, 0.174_fp, 0.179_fp,  &
    0.164_fp, 0.147_fp, 0.118_fp, 0.098_fp, 0.087_fp, 0.065_fp, 0.050_fp, 0.087_fp,  &
    0.099_fp, 0.124_fp, 0.104_fp, 0.095_fp, 0.041_fp, 0.055_fp, 0.039_fp, 0.032_fp,  &
    0.045_fp, 0.046_fp, 0.053_fp, 0.056_fp, 0.057_fp, 0.056_fp, 0.047_fp, 0.043_fp,  &
    0.041_fp, 0.040_fp/
    ! ...grass_soil
    DATA (reflectance(i,NPOESS_GRASS_SOIL),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.024_fp, 0.026_fp,  &
    0.028_fp, 0.031_fp, 0.044_fp, 0.048_fp, 0.056_fp, 0.087_fp, 0.101_fp, 0.107_fp,  &
    0.115_fp, 0.123_fp, 0.130_fp, 0.175_fp, 0.214_fp, 0.250_fp, 0.331_fp, 0.382_fp,  &
    0.403_fp, 0.410_fp, 0.414_fp, 0.416_fp, 0.422_fp, 0.373_fp, 0.380_fp, 0.387_fp,  &
    0.393_fp, 0.394_fp, 0.401_fp, 0.399_fp, 0.401_fp, 0.402_fp, 0.361_fp, 0.312_fp,  &
    0.247_fp, 0.183_fp, 0.149_fp, 0.182_fp, 0.210_fp, 0.234_fp, 0.249_fp, 0.257_fp,  &
    0.232_fp, 0.200_fp, 0.123_fp, 0.082_fp, 0.069_fp, 0.090_fp, 0.044_fp, 0.092_fp,  &
    0.164_fp, 0.179_fp, 0.185_fp, 0.178_fp, 0.060_fp, 0.072_fp, 0.048_fp, 0.036_fp,  &
    0.038_fp, 0.048_fp, 0.057_fp, 0.060_fp, 0.067_fp, 0.090_fp, 0.075_fp, 0.068_fp,  &
    0.054_fp, 0.048_fp/
    ! ...broadleaf_pine_forest
    DATA (reflectance(i,NPOESS_BROADLEAF_PINE_FOREST),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.017_fp, 0.018_fp, 0.019_fp, 0.020_fp, 0.020_fp, 0.021_fp, 0.021_fp, 0.024_fp,  &
    0.025_fp, 0.026_fp, 0.028_fp, 0.030_fp, 0.032_fp, 0.031_fp, 0.045_fp, 0.039_fp,  &
    0.035_fp, 0.027_fp, 0.026_fp, 0.025_fp, 0.056_fp, 0.127_fp, 0.252_fp, 0.283_fp,  &
    0.284_fp, 0.290_fp, 0.298_fp, 0.301_fp, 0.304_fp, 0.357_fp, 0.350_fp, 0.346_fp,  &
    0.331_fp, 0.341_fp, 0.344_fp, 0.323_fp, 0.294_fp, 0.307_fp, 0.264_fp, 0.214_fp,  &
    0.164_fp, 0.124_fp, 0.124_fp, 0.147_fp, 0.164_fp, 0.177_fp, 0.194_fp, 0.280_fp,  &
    0.170_fp, 0.131_fp, 0.098_fp, 0.093_fp, 0.085_fp, 0.090_fp, 0.081_fp, 0.081_fp,  &
    0.074_fp, 0.060_fp, 0.046_fp, 0.032_fp, 0.029_fp, 0.027_fp, 0.022_fp, 0.020_fp,  &
    0.018_fp, 0.018_fp, 0.018_fp, 0.018_fp, 0.018_fp, 0.018_fp, 0.022_fp, 0.018_fp,  &
    0.013_fp, 0.013_fp/
    ! ...grass_scrub
    DATA (reflectance(i,NPOESS_GRASS_SCRUB),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.021_fp,  &
    0.025_fp, 0.033_fp, 0.055_fp, 0.058_fp, 0.066_fp, 0.063_fp, 0.077_fp, 0.069_fp,  &
    0.059_fp, 0.050_fp, 0.048_fp, 0.087_fp, 0.130_fp, 0.192_fp, 0.287_fp, 0.357_fp,  &
    0.401_fp, 0.410_fp, 0.416_fp, 0.422_fp, 0.428_fp, 0.502_fp, 0.500_fp, 0.504_fp,  &
    0.509_fp, 0.502_fp, 0.496_fp, 0.489_fp, 0.482_fp, 0.475_fp, 0.408_fp, 0.340_fp,  &
    0.273_fp, 0.205_fp, 0.138_fp, 0.166_fp, 0.195_fp, 0.224_fp, 0.233_fp, 0.243_fp,  &
    0.210_fp, 0.177_fp, 0.145_fp, 0.112_fp, 0.079_fp, 0.096_fp, 0.050_fp, 0.095_fp,  &
    0.150_fp, 0.205_fp, 0.208_fp, 0.208_fp, 0.045_fp, 0.073_fp, 0.060_fp, 0.050_fp,  &
    0.050_fp, 0.050_fp, 0.049_fp, 0.055_fp, 0.070_fp, 0.108_fp, 0.105_fp, 0.088_fp,  &
    0.073_fp, 0.080_fp/
    ! ...soil_grass_scrub
    DATA (reflectance(i,NPOESS_SOIL_GRASS_SCRUB),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.021_fp, 0.022_fp,  &
    0.025_fp, 0.030_fp, 0.045_fp, 0.049_fp, 0.060_fp, 0.073_fp, 0.087_fp, 0.090_fp,  &
    0.092_fp, 0.093_fp, 0.098_fp, 0.131_fp, 0.160_fp, 0.195_fp, 0.260_fp, 0.309_fp,  &
    0.340_fp, 0.348_fp, 0.355_fp, 0.359_fp, 0.366_fp, 0.377_fp, 0.380_fp, 0.384_fp,  &
    0.387_fp, 0.385_fp, 0.387_fp, 0.383_fp, 0.382_fp, 0.380_fp, 0.339_fp, 0.294_fp,  &
    0.238_fp, 0.183_fp, 0.147_fp, 0.176_fp, 0.201_fp, 0.224_fp, 0.230_fp, 0.230_fp,  &
    0.208_fp, 0.181_fp, 0.125_fp, 0.092_fp, 0.077_fp, 0.078_fp, 0.046_fp, 0.085_fp,  &
    0.146_fp, 0.173_fp, 0.172_fp, 0.165_fp, 0.051_fp, 0.064_fp, 0.052_fp, 0.042_fp,  &
    0.042_fp, 0.046_fp, 0.049_fp, 0.053_fp, 0.060_fp, 0.081_fp, 0.073_fp, 0.065_fp,  &
    0.052_fp, 0.056_fp/
    ! ...urban_concrete
    DATA (reflectance(i,NPOESS_URBAN_CONCRETE),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.026_fp, 0.026_fp, 0.026_fp, 0.026_fp, 0.026_fp, 0.026_fp, 0.027_fp, 0.028_fp,  &
    0.028_fp, 0.029_fp, 0.031_fp, 0.032_fp, 0.035_fp, 0.037_fp, 0.043_fp, 0.047_fp,  &
    0.051_fp, 0.052_fp, 0.055_fp, 0.060_fp, 0.063_fp, 0.084_fp, 0.091_fp, 0.099_fp,  &
    0.100_fp, 0.104_fp, 0.109_fp, 0.116_fp, 0.119_fp, 0.126_fp, 0.129_fp, 0.130_fp,  &
    0.133_fp, 0.134_fp, 0.137_fp, 0.129_fp, 0.127_fp, 0.125_fp, 0.119_fp, 0.109_fp,  &
    0.094_fp, 0.083_fp, 0.078_fp, 0.087_fp, 0.094_fp, 0.099_fp, 0.102_fp, 0.107_fp,  &
    0.102_fp, 0.096_fp, 0.079_fp, 0.069_fp, 0.065_fp, 0.032_fp, 0.029_fp, 0.039_fp,  &
    0.033_fp, 0.074_fp, 0.072_fp, 0.058_fp, 0.026_fp, 0.023_fp, 0.023_fp, 0.021_fp,  &
    0.027_fp, 0.029_fp, 0.032_fp, 0.033_fp, 0.030_fp, 0.023_fp, 0.020_fp, 0.019_fp,  &
    0.016_fp, 0.016_fp/
    ! ...pine_brush
    DATA (reflectance(i,NPOESS_PINE_BRUSH),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.012_fp, 0.013_fp, 0.013_fp, 0.014_fp, 0.014_fp, 0.015_fp, 0.015_fp, 0.016_fp,  &
    0.018_fp, 0.024_fp, 0.035_fp, 0.039_fp, 0.049_fp, 0.038_fp, 0.047_fp, 0.045_fp,  &
    0.039_fp, 0.032_fp, 0.031_fp, 0.032_fp, 0.055_fp, 0.104_fp, 0.162_fp, 0.187_fp,  &
    0.218_fp, 0.230_fp, 0.240_fp, 0.247_fp, 0.255_fp, 0.333_fp, 0.328_fp, 0.324_fp,  &
    0.312_fp, 0.315_fp, 0.313_fp, 0.294_fp, 0.269_fp, 0.275_fp, 0.240_fp, 0.199_fp,  &
    0.159_fp, 0.126_fp, 0.120_fp, 0.139_fp, 0.153_fp, 0.165_fp, 0.160_fp, 0.205_fp,  &
    0.134_fp, 0.113_fp, 0.097_fp, 0.100_fp, 0.102_fp, 0.081_fp, 0.075_fp, 0.080_fp,  &
    0.095_fp, 0.110_fp, 0.090_fp, 0.075_fp, 0.028_fp, 0.031_fp, 0.038_fp, 0.034_fp,  &
    0.030_fp, 0.025_fp, 0.020_fp, 0.023_fp, 0.025_fp, 0.030_fp, 0.035_fp, 0.030_fp,  &
    0.025_fp, 0.040_fp/
    ! ...broadleaf_brush
    DATA (reflectance(i,NPOESS_BROADLEAF_BRUSH),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.019_fp, 0.019_fp, 0.020_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.018_fp, 0.021_fp,  &
    0.023_fp, 0.029_fp, 0.041_fp, 0.045_fp, 0.055_fp, 0.033_fp, 0.055_fp, 0.048_fp,  &
    0.040_fp, 0.023_fp, 0.024_fp, 0.025_fp, 0.026_fp, 0.064_fp, 0.146_fp, 0.203_fp,  &
    0.227_fp, 0.234_fp, 0.240_fp, 0.245_fp, 0.252_fp, 0.399_fp, 0.392_fp, 0.390_fp,  &
    0.386_fp, 0.380_fp, 0.374_fp, 0.368_fp, 0.362_fp, 0.356_fp, 0.318_fp, 0.279_fp,  &
    0.240_fp, 0.202_fp, 0.163_fp, 0.187_fp, 0.211_fp, 0.235_fp, 0.236_fp, 0.236_fp,  &
    0.203_fp, 0.171_fp, 0.139_fp, 0.106_fp, 0.074_fp, 0.033_fp, 0.045_fp, 0.058_fp,  &
    0.075_fp, 0.105_fp, 0.093_fp, 0.087_fp, 0.033_fp, 0.039_fp, 0.051_fp, 0.048_fp,  &
    0.045_fp, 0.039_fp, 0.033_fp, 0.036_fp, 0.039_fp, 0.045_fp, 0.056_fp, 0.045_fp,  &
    0.032_fp, 0.050_fp/
    ! ...wet_soil
    DATA (reflectance(i,NPOESS_WET_SOIL),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.029_fp, 0.029_fp, 0.029_fp, 0.029_fp, 0.026_fp, 0.026_fp, 0.026_fp, 0.026_fp,  &
    0.035_fp, 0.038_fp, 0.041_fp, 0.041_fp, 0.038_fp, 0.057_fp, 0.063_fp, 0.073_fp,  &
    0.082_fp, 0.082_fp, 0.090_fp, 0.102_fp, 0.106_fp, 0.103_fp, 0.113_fp, 0.122_fp,  &
    0.128_fp, 0.132_fp, 0.135_fp, 0.136_fp, 0.140_fp, 0.107_fp, 0.112_fp, 0.114_fp,  &
    0.114_fp, 0.117_fp, 0.125_fp, 0.125_fp, 0.128_fp, 0.131_fp, 0.131_fp, 0.125_fp,  &
    0.106_fp, 0.087_fp, 0.093_fp, 0.107_fp, 0.118_fp, 0.124_fp, 0.124_fp, 0.118_fp,  &
    0.114_fp, 0.105_fp, 0.059_fp, 0.043_fp, 0.049_fp, 0.026_fp, 0.021_fp, 0.045_fp,  &
    0.080_fp, 0.073_fp, 0.069_fp, 0.055_fp, 0.035_fp, 0.035_fp, 0.030_fp, 0.025_fp,  &
    0.025_fp, 0.025_fp, 0.030_fp, 0.030_fp, 0.028_fp, 0.025_fp, 0.018_fp, 0.025_fp,  &
    0.025_fp, 0.025_fp/
    ! ...scrub_soil
    DATA (reflectance(i,NPOESS_SCRUB_SOIL),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp, 0.019_fp,  &
    0.020_fp, 0.025_fp, 0.036_fp, 0.042_fp, 0.057_fp, 0.069_fp, 0.083_fp, 0.094_fp,  &
    0.102_fp, 0.107_fp, 0.117_fp, 0.132_fp, 0.138_fp, 0.143_fp, 0.161_fp, 0.188_fp,  &
    0.215_fp, 0.226_fp, 0.234_fp, 0.238_fp, 0.248_fp, 0.256_fp, 0.260_fp, 0.262_fp,  &
    0.260_fp, 0.260_fp, 0.265_fp, 0.262_fp, 0.262_fp, 0.262_fp, 0.250_fp, 0.230_fp,  &
    0.195_fp, 0.160_fp, 0.155_fp, 0.179_fp, 0.200_fp, 0.215_fp, 0.207_fp, 0.191_fp,  &
    0.182_fp, 0.166_fp, 0.106_fp, 0.081_fp, 0.084_fp, 0.047_fp, 0.044_fp, 0.068_fp,  &
    0.124_fp, 0.135_fp, 0.123_fp, 0.108_fp, 0.048_fp, 0.046_fp, 0.048_fp, 0.040_fp,  &
    0.038_fp, 0.040_fp, 0.042_fp, 0.044_fp, 0.043_fp, 0.044_fp, 0.039_fp, 0.038_fp,  &
    0.028_fp, 0.040_fp/
    ! ...broadleaf70_pine30
    DATA (reflectance(i,NPOESS_BROADLEAF70_PINE30),i=1,N_NPOESS_WAVELENGTHS)/ &
    0.025_fp, 0.025_fp, 0.027_fp, 0.026_fp, 0.026_fp, 0.025_fp, 0.024_fp, 0.030_fp,  &
    0.031_fp, 0.031_fp, 0.033_fp, 0.035_fp, 0.036_fp, 0.027_fp, 0.051_fp, 0.041_fp,  &
    0.035_fp, 0.018_fp, 0.019_fp, 0.018_fp, 0.032_fp, 0.097_fp, 0.254_fp, 0.314_fp,  &
    0.304_fp, 0.304_fp, 0.308_fp, 0.309_fp, 0.310_fp, 0.419_fp, 0.411_fp, 0.407_fp,  &
    0.400_fp, 0.402_fp, 0.402_fp, 0.393_fp, 0.380_fp, 0.384_fp, 0.337_fp, 0.287_fp,  &
    0.236_fp, 0.191_fp, 0.162_fp, 0.190_fp, 0.216_fp, 0.240_fp, 0.265_fp, 0.321_fp,  &
    0.238_fp, 0.185_fp, 0.135_fp, 0.098_fp, 0.058_fp, 0.050_fp, 0.057_fp, 0.062_fp,  &
    0.054_fp, 0.048_fp, 0.042_fp, 0.036_fp, 0.034_fp, 0.033_fp, 0.031_fp, 0.030_fp,  &
    0.030_fp, 0.030_fp, 0.030_fp, 0.030_fp, 0.030_fp, 0.030_fp, 0.038_fp, 0.030_fp,  &
    0.016_fp, 0.016_fp/
    
    DO i=1,N_NPOESS_WAVELENGTHS
        wavenumber(i)     =  1.0E4_fp/wavelength(N_NPOESS_WAVELENGTHS-i+1)
        reflectivity(i,:) =  reflectance(N_NPOESS_WAVELENGTHS-i+1,:)
    ENDDO
    npoess_lut_wavenumbers  = wavenumber
    npoess_lut_reflectivity = reflectivity
    npoess_lut_reflectivity_0 = reflectance
   END SUBROUTINE Load_NPOESS_LUT


END MODULE NPOESS_LUT_READER
