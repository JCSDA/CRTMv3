
!
! CSEM_Define
!
! Module to define the general CSEM  data structures. These data structures 
! are used to develop generic interfaces handling the top-down and bottomo-up 
! data streams. 
!
!  
! Some FORTRAN 2003 features are used in the type class constructor for dynamic 
! array allocation and in the destructor for array deallocation
! 
!
! CREATION HISTORY:
!       Written by:     Ming Chen, 09-Jul-2016
!                       ming.chen@noaa.gov
!




MODULE CSEM_Define

  ! -----------------
  ! Environment setup
  ! -----------------
  ! Module use
  USE CSEM_Type_Kinds, ONLY: fp => CSEM_fp
  
  ! ---------------------------
  ! Disable all implicit typing
  ! ---------------------------
  IMPLICIT NONE
  

  ! ----------------------------
  ! Land Surface structure definition
  ! ----------------------------
  TYPE :: CSEM_Land_Surface
    INTEGER  :: Land_Cover_Type        =  1         
    INTEGER  :: Vegetation_Type        =  1                   
    INTEGER  :: Soil_Type              =  1          
    REAL(fp) :: Vegetation_Fraction    =  0.3_fp    ! 30%
    REAL(fp) :: Land_Skin_Temperature  =  283.0_fp  ! K
    REAL(fp) :: Top_Soil_Temperature   =  283.0_fp  ! K  
    REAL(fp) :: Top_Soil_Moisture      =  0.05_fp   ! g/cm^3  
    REAL(fp) :: LAI                    =  3.5 
    REAL(fp) :: Canopy_Water_Content   =  0.05_fp   ! g/cm^3
    INTEGER  :: n_Soil_Layers          =  0          
    LOGICAL  :: Is_Allocated           = .FALSE.         

    REAL(fp), ALLOCATABLE :: Temperature_Profile(:)    
    REAL(fp), ALLOCATABLE :: Moisture_Profile(:)     
    REAL(fp), ALLOCATABLE :: Soil_Depth(:)  
    
  CONTAINS
        PROCEDURE :: init=>Alloc_Soil_Profile
        FINAL     :: Clean_Land
             
  END TYPE CSEM_Land_Surface
  
  ! ----------------------------
  ! Water Surface structure definition
  ! ----------------------------
  TYPE :: CSEM_Water_Surface
    ! Water type data
    INTEGER  :: Water_Type             =  1              
    REAL(fp) :: Water_Temperature      =  283.0_fp   ! K
    REAL(fp) :: Wind_Speed             =  5.0_fp     ! m/s
    REAL(fp) :: Wind_Direction         =  0.0_fp     ! North
    REAL(fp) :: Salinity               =  33.0_fp    ! ppmv
    REAL(fp) :: Foam_Fraction          =  0.0_fp     ! N/A
  END TYPE CSEM_Water_Surface

  ! ----------------------------
  ! Snow Surface structure definition
  ! ----------------------------  
  TYPE :: CSEM_Snow_Surface
    ! Snow surface type data
    INTEGER  :: Snow_Type              =  1          ! First item in list  
    REAL(fp) :: Snow_Temperature       =  263.0_fp   ! K
    REAL(fp) :: Snow_Depth             =  50.0_fp    ! mm
    REAL(fp) :: Snow_Density           =  0.2_fp     ! g/cm^3
    REAL(fp) :: Snow_Grain_Size        =  2.0_fp     ! mm
    INTEGER  :: Soil_Type              =  1          ! First item in list 
    REAL(fp) :: Top_Soil_Temperature   =  283.0_fp   ! K  
    REAL(fp) :: Top_Soil_Moisture_Content = 0.05_fp  ! g/cm^3 
    INTEGER  :: Vegetation_Type        =  1         ! First item in list          
    REAL(fp) :: LAI                    =  3.5 
  END TYPE CSEM_Snow_Surface
  
  ! ----------------------------
  ! Sea ice Surface structure definition
  ! ----------------------------  
  TYPE :: CSEM_Ice_Surface
    ! Ice surface type data
    INTEGER  :: Ice_Type               =  1          ! First item in list  
    REAL(fp) :: Ice_Temperature        =  263.0_fp  ! K
    REAL(fp) :: Ice_Thickness          =  10.0_fp   ! mm
    REAL(fp) :: Ice_Density            =  0.9_fp    ! g/cm^3
    REAL(fp) :: Ice_Roughness          =  0.0_fp
    REAL(fp) :: Salinity               =  33.0_fp    ! ppmv
  END TYPE CSEM_Ice_Surface
  
 
  ! ----------------------------
  ! Surface Optic output structure 
  ! ----------------------------
  TYPE :: CSEM_SfcOptics_Type
    LOGICAL  :: Is_Allocated  =  .FALSE.
    LOGICAL  :: Is_Solar      =  .FALSE.
    LOGICAL  :: Is_Spectral   =  .FALSE.
    REAL(fp) :: Frequency
    REAL(fp) :: Wavenumber
    REAL(fp) :: Source_Zenith_Angle    = 0.0_fp
    REAL(fp) :: Source_Azimuth_Angle   = 0.0_fp  
    REAL(fp) :: Sensor_Zenith_Angle    = 0.0_fp
    REAL(fp) :: Sensor_Scan_Angle      = 0.0_fp
    REAL(fp) :: Sensor_Azimuth_Angle   = 0.0_fp  
    REAL(fp) :: Relative_Azimuth_Angle = 0.0_fp  
    
    INTEGER  :: n_Angles      = 1 
    INTEGER  :: n_Stokes      = 4 

    ! The counter for the m'th component of the Fourier exapnsion of
    ! the radiance for azimuth angle
    INTEGER  :: mth_Azi = 0

    REAL(fp), ALLOCATABLE :: Angle(:)
    REAL(fp), ALLOCATABLE :: Weight(:)                ! I
    REAL(fp), ALLOCATABLE :: Emissivity(:,:) 
    REAL(fp), ALLOCATABLE :: Direct_Reflectivity(:,:) 
    REAL(fp), ALLOCATABLE :: Reflectivity(:,:,:,:)    ! I x Ls x I x Ls

  CONTAINS
        PROCEDURE, PASS(self) :: init  => Init_SfcOptics
        FINAL     :: Clean_SfcOptics
  END TYPE CSEM_SfcOptics_Type
  
  ! ---------------------------------------------------
  ! Extended Data structure definitions
  ! These data structures are used to provide
  ! additional inputs for empirical and semi-empirical
  ! model development  
  ! --------------------------------------------------
  
  ! -----------------------------------------  
  ! Sensor Observation-based structure definition
  !--------------------------------------------  
  TYPE :: CSEM_SensorObs_Struct
    ! The data sensor IDs
    CHARACTER(LEN=100) :: Sensor_Id        = ' '
    LOGICAL :: Is_Allocated = .FALSE.
    ! Dimension values
    INTEGER :: n_Channels       =   0 
    
    ! Geometric
    !REAL(fp) :: Scan_Angle      =   0.0_fp
    !REAL(fp) :: Zenith_Angle    =   0.0_fp
    !REAL(fp) :: Azimuth_Angle   =   0.0_fp

    ! The sensor channels and brightness temperatures
    REAL(fp), ALLOCATABLE :: Channel_Frequency(:) 
    INTEGER,  ALLOCATABLE :: Channel_Polarization(:)            
    REAL(fp), ALLOCATABLE :: Tb(:)     
  CONTAINS
        PROCEDURE, PASS(self) :: init  => Alloc_SensorObs
        FINAL     :: Clean_SensorObs
  END TYPE CSEM_SensorObs_Struct


  ! -----------------------------------------  
  ! Spacetime metadata structure definition
  !--------------------------------------------   
  TYPE :: CSEM_GeoInfo_Struct
    ! Geolocation and time
    REAL(fp) :: Latitude        =   0.0_fp
    REAL(fp) :: Longitude       =   0.0_fp
    INTEGER  :: Year            =   2001
    INTEGER  :: Month           =   1
    INTEGER  :: Day             =   1
    INTEGER  :: Hour            =   1

  END TYPE CSEM_GeoInfo_Struct
  
  
  ! -------------------------------------------  
  ! AtmosphericInfo structure definition
  !--------------------------------------------
  TYPE :: CSEM_Atmosphere_Parameters
    REAL(fp) :: Downward_Atm_Radiance = 0.0_fp
    REAL(fp) :: Transmittance         = 0.0_fp
    REAL(fp) :: Downward_Solar_Irradiance  = 0.0_fp
  END TYPE CSEM_Atmosphere_Parameters
  
  ! ---------------------------------------------- 
  ! Composite Optional data structure definition
  !-----------------------------------------------
  TYPE :: CSEM_Options_Type
   TYPE(CSEM_SensorObs_Struct)          :: SensorObs
   TYPE(CSEM_Atmosphere_Parameters)     :: Atmos 
   TYPE(CSEM_GeoInfo_Struct)            :: GeoInfo
  END TYPE CSEM_Options_Type
  
  
CONTAINS

 
  ! -------------------------------------------
  ! Type-bined procedure for the allocation
  ! of soile profiles in the land data structure
  ! --------------------------------------------   
    SUBROUTINE Alloc_Soil_Profile(land, n_Layers)
      CLASS(CSEM_Land_Surface) :: land  
      INTEGER :: n_Layers
      INTEGER :: Alloc_Status
      ALLOCATE(land%Temperature_Profile(n_Layers), &
               land%Moisture_Profile(n_Layers),    &
               land%Soil_Depth(n_Layers), STAT = Alloc_Status) 
      IF ( Alloc_Status /= 0 ) RETURN
      land%n_Soil_Layers = n_Layers
      land%Is_Allocated = .TRUE.
    END SUBROUTINE Alloc_Soil_Profile
    
  ! -------------------------------------------
  ! Type-bined procedure for the deallocation
  ! of soile profiles in the land data structure
  ! --------------------------------------------  
    ELEMENTAL SUBROUTINE Clean_Land(land)
      TYPE(CSEM_Land_Surface), INTENT(INOUT) :: land
      land%n_Soil_Layers = 0
      land%Is_Allocated = .FALSE.
    END SUBROUTINE Clean_Land
    

  ! --------------------------------------------
  ! Type-bined procedure for the allocation of the
  ! dynamic arrays in the surface optics structure
  ! --------------------------------------------  
    SUBROUTINE Init_SfcOptics(self, n_Angles)
      CLASS(CSEM_SfcOptics_Type) :: self
      INTEGER, OPTIONAL :: n_Angles
      INTEGER :: Alloc_Status
      self%Is_Allocated = .FALSE.
      IF(PRESENT(n_Angles))self%n_Angles   = n_Angles
      ALLOCATE(self%Emissivity(self%n_Angles,self%n_Stokes), &
               self%Direct_Reflectivity(self%n_Angles,self%n_Stokes), &
               self%Reflectivity(self%n_Angles, self%n_Stokes, &
               self%n_Angles, self%n_Stokes), &
               STAT = Alloc_Status)
      IF ( Alloc_Status /= 0 )  RETURN
      ALLOCATE(self%Angle(self%n_Angles),  &
               self%Weight(self%n_Angles), &
               STAT = Alloc_Status)
      IF ( Alloc_Status /= 0 )  RETURN
      self%Is_Allocated = .TRUE.
    END SUBROUTINE Init_SfcOptics

  ! --------------------------------------------
  ! Type-bined procedure for the allocation of the
  ! dynamic arrays in the sensor obs data structure
  ! --------------------------------------------  
    SUBROUTINE Alloc_SensorObs(self, n_Channels)
      CLASS(CSEM_SensorObs_Struct) :: self
      INTEGER :: n_Channels
      INTEGER :: Alloc_Status
      self%Is_Allocated = .FALSE.
      IF(n_Channels < 1) RETURN
      ALLOCATE(self%Channel_Frequency(n_Channels),    &
               self%Channel_Polarization(n_Channels), &
               self%TB(n_Channels), STAT = Alloc_Status)
      IF ( Alloc_Status /= 0 )  RETURN
      self%n_Channels = n_Channels
      self%Is_Allocated = .TRUE.
    END SUBROUTINE Alloc_SensorObs
    

  ! --------------------------------------------
  ! Type-bined procedure for the deallocation of the
  ! dynamic arrays in the sensor obs data structure
  ! --------------------------------------------  
    ELEMENTAL SUBROUTINE Clean_SensorObs(sensor)
      TYPE(CSEM_SensorObs_Struct), INTENT(INOUT) :: sensor
      sensor%n_Channels = 0
      sensor%Is_Allocated = .FALSE.
      
    END SUBROUTINE Clean_SensorObs
  ! --------------------------------------------
  ! Type-bined procedure for the deallocation of the
  ! dynamic arrays in the sensor obs data structure
  ! --------------------------------------------  
    ELEMENTAL SUBROUTINE Clean_SfcOptics(sfcOptics)
      TYPE(CSEM_SfcOptics_Type), INTENT(INOUT) :: sfcOptics
      sfcOptics%Is_Allocated = .FALSE.
       
    END SUBROUTINE Clean_SfcOptics
    
END MODULE CSEM_Define
