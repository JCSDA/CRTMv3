program test_crtm_onnx_init
    use CRTM_Module
    implicit none

    integer :: status
    type(CRTM_ChannelInfo_type) :: chinfo(1)
    type(CRTM_Atmosphere_type) :: atm(1)
    type(CRTM_Surface_type) :: sfc(1)
    type(CRTM_Geometry_type) :: geo(1)
    type(CRTM_RTSolution_type), allocatable :: rts(:,:)
    integer :: n_layers = 100
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

    ! 2. Setup dummy atmosphere
    call CRTM_Atmosphere_Create(atm(1), n_layers, 3, 0, 0)
    atm(1)%Absorber_ID(1:3) = [H2O_ID, CO2_ID, O3_ID]
    atm(1)%Absorber_Units(1:3) = MASS_MIXING_RATIO_UNITS
    atm(1)%Pressure(1:n_layers) = 1000.0 ! dummy
    atm(1)%Temperature(1:n_layers) = 300.0 ! dummy
    atm(1)%Absorber(1:n_layers, 1) = 0.01 ! dummy H2O
    atm(1)%Absorber(1:n_layers, 2) = 0.0003 ! dummy CO2
    atm(1)%Absorber(1:n_layers, 3) = 0.00001 ! dummy O3

    ! 3. Setup dummy surface and geometry
    call CRTM_Surface_Create(sfc(1), chinfo(1)%n_Channels)
    sfc(1)%Land_Coverage = 1.0
    sfc(1)%Land_Type = 1
    
    call CRTM_Geometry_Create(geo(1))
    geo(1)%Sensor_Zenith_Angle = 0.0

    ! 4. Allocate RTSolution
    allocate(rts(chinfo(1)%n_Channels, 1))
    do l = 1, chinfo(1)%n_Channels
       rts(l, 1)%Sensor_Id = sensor_id
       rts(l, 1)%Sensor_Channel = chinfo(1)%Sensor_Channel(l)
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
