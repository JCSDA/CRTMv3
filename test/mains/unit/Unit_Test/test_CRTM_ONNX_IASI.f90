PROGRAM test_crtm_onnx_iasi
  USE CRTM_Module
  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: SENSOR_ID = 'iasi_metop-a'
  INTEGER,      PARAMETER :: N_PROFILES = 1
  INTEGER,      PARAMETER :: N_LAYERS   = 92
  INTEGER,      PARAMETER :: N_ABSORBERS = 2
  
  TYPE(CRTM_ChannelInfo_type) :: chinfo(1)
  TYPE(CRTM_Atmosphere_type)  :: atm(1)
  TYPE(CRTM_Surface_type)     :: sfc(1)
  TYPE(CRTM_Geometry_type)    :: geo(1)
  TYPE(CRTM_RTSolution_type), ALLOCATABLE :: rts(:,:)
  
  INTEGER :: status, n_channels, l
  
  PRINT *, "Initializing CRTM with ONNX support for ", SENSOR_ID
  
  ! 1. Initialize
  status = CRTM_Init([CHARACTER(len=20) :: SENSOR_ID], chinfo, Use_ONNX=.TRUE., File_Path="testinput/")
  IF (status /= SUCCESS) THEN
    PRINT *, "CRTM_Init failed"
    STOP
  END IF
  
  n_channels = chinfo(1)%n_Channels
  ALLOCATE(rts(n_channels, N_PROFILES))
  
  ! 2. Create and Load Data
  CALL CRTM_Atmosphere_Create(atm, N_LAYERS, N_ABSORBERS, 0, 0)
  CALL CRTM_Surface_Create(sfc, n_channels)
  CALL CRTM_Geometry_Create(geo)
  
  CALL Load_Atm_Data()
  CALL Load_Sfc_Data()
  
  CALL CRTM_Geometry_SetValue(geo, Sensor_Zenith_Angle = 30.0_fp)
  
  DO l = 1, n_channels
    rts(l, 1)%Sensor_Id = chinfo(1)%Sensor_Id
    rts(l, 1)%Sensor_Channel = chinfo(1)%Sensor_Channel(l)
       rts(l, 1)%WMO_Satellite_Id = chinfo(1)%WMO_Satellite_Id
       rts(l, 1)%WMO_Sensor_Id = chinfo(1)%WMO_Sensor_Id
  END DO

  ! 3. Call Forward
  PRINT *, "Calling CRTM_Forward with ONNX..."
  status = CRTM_Forward(atm, sfc, geo, chinfo, rts)
  
  IF (status == SUCCESS) THEN
    PRINT *, "CRTM_Forward successful"
    PRINT *, "First 5 radiances: ", rts(1:5, 1)%Radiance
  ELSE
    PRINT *, "CRTM_Forward failed"
  END IF

  ! 4. Cleanup
  status = CRTM_Destroy(chinfo)
  CALL CRTM_Atmosphere_Destroy(atm)
  CALL CRTM_Surface_Destroy(sfc)
  CALL CRTM_Geometry_Destroy(geo)
  DEALLOCATE(rts)

CONTAINS

  SUBROUTINE Load_Atm_Data()
    atm(1)%Climatology         = US_STANDARD_ATMOSPHERE
    atm(1)%Absorber_Id(1:2)    = (/ H2O_ID                 , O3_ID /)
    atm(1)%Absorber_Units(1:2) = (/ MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS /)
    atm(1)%Level_Pressure = &
    (/0.714_fp,   0.975_fp,   1.297_fp,   1.687_fp,   2.153_fp,   2.701_fp,   3.340_fp,   4.077_fp, &
      4.920_fp,   5.878_fp,   6.957_fp,   8.165_fp,   9.512_fp,  11.004_fp,  12.649_fp,  14.456_fp, &
     16.432_fp,  18.585_fp,  20.922_fp,  23.453_fp,  26.183_fp,  29.121_fp,  32.274_fp,  35.650_fp, &
     39.257_fp,  43.100_fp,  47.188_fp,  51.528_fp,  56.126_fp,  60.990_fp,  66.125_fp,  71.540_fp, &
     77.240_fp,  83.231_fp,  89.520_fp,  96.114_fp, 103.017_fp, 110.237_fp, 117.777_fp, 125.646_fp, &
    133.846_fp, 142.385_fp, 151.266_fp, 160.496_fp, 170.078_fp, 180.018_fp, 190.320_fp, 200.989_fp, &
    212.028_fp, 223.441_fp, 235.234_fp, 247.409_fp, 259.969_fp, 272.919_fp, 286.262_fp, 300.000_fp, &
    314.137_fp, 328.675_fp, 343.618_fp, 358.967_fp, 374.724_fp, 390.893_fp, 407.474_fp, 424.470_fp, &
    441.882_fp, 459.712_fp, 477.961_fp, 496.630_fp, 515.720_fp, 535.232_fp, 555.167_fp, 575.525_fp, &
    596.306_fp, 617.511_fp, 639.140_fp, 661.192_fp, 683.667_fp, 706.565_fp, 729.886_fp, 753.627_fp, &
    777.790_fp, 802.371_fp, 827.371_fp, 852.788_fp, 878.620_fp, 904.866_fp, 931.524_fp, 958.591_fp, &
    986.067_fp,1013.948_fp,1042.232_fp,1070.917_fp,1100.000_fp/)
    atm(1)%Pressure = &
    (/0.838_fp,   1.129_fp,   1.484_fp,   1.910_fp,   2.416_fp,   3.009_fp,   3.696_fp,   4.485_fp, &
      5.385_fp,   6.402_fp,   7.545_fp,   8.822_fp,  10.240_fp,  11.807_fp,  13.532_fp,  15.423_fp, &
     17.486_fp,  19.730_fp,  22.163_fp,  24.793_fp,  27.626_fp,  30.671_fp,  33.934_fp,  37.425_fp, &
     41.148_fp,  45.113_fp,  49.326_fp,  53.794_fp,  58.524_fp,  63.523_fp,  68.797_fp,  74.353_fp, &
     80.198_fp,  86.338_fp,  92.778_fp,  99.526_fp, 106.586_fp, 113.965_fp, 121.669_fp, 129.703_fp, &
    138.072_fp, 146.781_fp, 155.836_fp, 165.241_fp, 175.001_fp, 185.121_fp, 195.606_fp, 206.459_fp, &
    217.685_fp, 229.287_fp, 241.270_fp, 253.637_fp, 266.392_fp, 279.537_fp, 293.077_fp, 307.014_fp, &
    321.351_fp, 336.091_fp, 351.236_fp, 366.789_fp, 382.751_fp, 399.126_fp, 415.914_fp, 433.118_fp, &
    450.738_fp, 468.777_fp, 487.236_fp, 506.115_fp, 525.416_fp, 545.139_fp, 565.285_fp, 585.854_fp, &
    606.847_fp, 628.263_fp, 650.104_fp, 672.367_fp, 695.054_fp, 718.163_fp, 741.693_fp, 765.645_fp, &
    790.017_fp, 814.807_fp, 840.016_fp, 865.640_fp, 891.679_fp, 918.130_fp, 944.993_fp, 972.264_fp, &
    999.942_fp,1028.025_fp,1056.510_fp,1085.394_fp/)
    atm(1)%Temperature = &
    (/256.186_fp, 252.608_fp, 247.762_fp, 243.314_fp, 239.018_fp, 235.282_fp, 233.777_fp, 234.909_fp, &
      237.889_fp, 241.238_fp, 243.194_fp, 243.304_fp, 242.977_fp, 243.133_fp, 242.920_fp, 242.026_fp, &
      240.695_fp, 239.379_fp, 238.252_fp, 236.928_fp, 235.452_fp, 234.561_fp, 234.192_fp, 233.774_fp, &
      233.305_fp, 233.053_fp, 233.103_fp, 233.307_fp, 233.702_fp, 234.219_fp, 234.959_fp, 235.940_fp, &
      236.744_fp, 237.155_fp, 237.374_fp, 238.244_fp, 239.736_fp, 240.672_fp, 240.688_fp, 240.318_fp, &
      239.888_fp, 239.411_fp, 238.512_fp, 237.048_fp, 235.388_fp, 233.551_fp, 231.620_fp, 230.418_fp, &
      229.927_fp, 229.511_fp, 229.197_fp, 228.947_fp, 228.772_fp, 228.649_fp, 228.567_fp, 228.517_fp, &
      228.614_fp, 228.861_fp, 229.376_fp, 230.223_fp, 231.291_fp, 232.591_fp, 234.013_fp, 235.508_fp, &
      237.041_fp, 238.589_fp, 240.165_fp, 241.781_fp, 243.399_fp, 244.985_fp, 246.495_fp, 247.918_fp, &
      249.228_fp, 250.418_fp, 251.488_fp, 252.438_fp, 253.268_fp, 253.978_fp, 254.568_fp, 255.038_fp, &
      255.388_fp, 255.618_fp, 255.728_fp, 255.718_fp, 255.588_fp, 255.338_fp, 254.968_fp, 254.478_fp, &
      253.868_fp, 253.138_fp, 252.288_fp, 251.318_fp/)
    atm(1)%Absorber(1:92, 1) = 0.01_fp ! dummy H2O
    atm(1)%Absorber(1:92, 2) = 0.00001_fp ! dummy O3
  END SUBROUTINE Load_Atm_Data

  SUBROUTINE Load_Sfc_Data()
    sfc(1)%Water_Coverage    = 1.0_fp
    sfc(1)%Water_Type        = 1
    sfc(1)%Water_Temperature = 275.0_fp
  END SUBROUTINE Load_Sfc_Data

END PROGRAM test_crtm_onnx_iasi
