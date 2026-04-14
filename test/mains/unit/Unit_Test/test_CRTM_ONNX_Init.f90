program test_crtm_onnx_init
    use CRTM_Module
    implicit none

    integer :: status
    type(CRTM_ChannelInfo_type) :: chinfo(1)
    type(CRTM_Atmosphere_type) :: atm(1)
    type(CRTM_Surface_type) :: sfc(1)
    type(CRTM_Geometry_type) :: geo(1)
    type(CRTM_RTSolution_type), allocatable :: rts(:,:)
    integer :: n_layers = 92
    character(len=255) :: sensor_id = "iasi_metop-a"
    integer :: l

    print *, "Initializing CRTM with ONNX support for ", trim(sensor_id)
    
    ! 1. Initialize with Use_ONNX=.TRUE.
    status = CRTM_Init([sensor_id], chinfo, Use_ONNX=.TRUE., File_Path="testinput/")
    if (status /= SUCCESS) then
        print *, "CRTM_Init failed"
        stop
    end if

    print *, "CRTM_Init successful. Use_ONNX = ", chinfo(1)%Use_ONNX

    ! 2. Setup atmosphere (using 92 layers to match typical IASI configs)
    call CRTM_Atmosphere_Create(atm(1), n_layers, 2, 0, 0)
    atm(1)%Absorber_ID(1:2) = [H2O_ID, O3_ID]
    atm(1)%Absorber_Units(1:2) = [MASS_MIXING_RATIO_UNITS, VOLUME_MIXING_RATIO_UNITS]
    atm(1)%Level_Pressure = (/(l*10.0_fp, l=0,n_layers)/)
    atm(1)%Pressure = (/(l*10.0_fp + 5.0_fp, l=0,n_layers-1)/)
    atm(1)%Temperature = 250.0_fp
    atm(1)%Absorber = 0.001_fp

    ! 3. Setup surface and geometry
    call CRTM_Surface_Create(sfc(1), chinfo(1)%n_Channels)
    sfc(1)%Water_Coverage = 1.0_fp
    sfc(1)%Water_Type = 1
    
    ! Fill SensorData to pass validation
    sfc(1)%SensorData%Sensor_Id = chinfo(1)%Sensor_Id
    sfc(1)%SensorData%Sensor_Channel = chinfo(1)%Sensor_Channel
    sfc(1)%SensorData%Tb = 280.0_fp
    
    call CRTM_Geometry_Create(geo(1))
    geo(1)%Sensor_Zenith_Angle = 0.0

    ! 4. Allocate RTSolution
    allocate(rts(chinfo(1)%n_Channels, 1))
    do l = 1, chinfo(1)%n_Channels
       rts(l, 1)%Sensor_Id = chinfo(1)%Sensor_Id
       rts(l, 1)%Sensor_Channel = chinfo(1)%Sensor_Channel(l)
       rts(l, 1)%WMO_Satellite_Id = chinfo(1)%WMO_Satellite_Id
       rts(l, 1)%WMO_Sensor_Id = chinfo(1)%WMO_Sensor_Id
    end do

    ! 5. Call CRTM_Forward (this will trigger the ONNX path internally)
    print *, "Calling CRTM_Forward..."
    status = CRTM_Forward(atm, sfc, geo, chinfo, rts)
    
    if (status == SUCCESS) then
        print *, "CRTM_Forward successful with ONNX"
        print *, "First 5 radiances: ", rts(1:5, 1)%Radiance
    else
        print *, "CRTM_Forward failed with ONNX"
    end if

    ! 6. Cleanup
    status = CRTM_Destroy(chinfo)
    call CRTM_Atmosphere_Destroy(atm(1))
    call CRTM_Surface_Destroy(sfc(1))
    call CRTM_Geometry_Destroy(geo(1))
    deallocate(rts)

end program test_crtm_onnx_init
