program test_onnx_bridge
    use crtm_onnx_interface
    use iso_c_binding
    implicit none

    integer(c_int) :: status
    real(c_float), dimension(6) :: features  ! P, T, H2O, CO2, O3, Angle
    real(c_float), dimension(8461) :: transmittances
    character(kind=c_char), dimension(:), allocatable :: model_path_c
    character(len=255) :: model_path_fortran
    character(len=255) :: sensor_id
    integer :: i

    ! Get Sensor_ID from command line
    if (command_argument_count() < 1) then
        print *, "Usage: test_ONNX_Interface <sensor_id>"
        print *, "Example: test_ONNX_Interface iasi_metop-a"
        stop
    end if
    call get_command_argument(1, sensor_id)
    sensor_id = trim(adjustl(sensor_id))

    ! Path: ./testinput/ONNX/<sensor_id>/model.onnx
    model_path_fortran = "./testinput/ONNX/" // trim(sensor_id) // "/model.onnx"
    print *, "Loading model: ", trim(model_path_fortran)

    allocate(model_path_c(len_trim(model_path_fortran) + 1))
    
    do i = 1, len_trim(model_path_fortran)
        model_path_c(i) = model_path_fortran(i:i)
    end do
    model_path_c(len_trim(model_path_fortran) + 1) = c_null_char

    status = crtm_onnx_init(model_path_c)
    if (status /= 0) then
        print *, "Failed to initialize ONNX bridge"
        stop
    end if

    ! Dummy features (normalized)
    features = [0.0, 0.0, 0.0, 0.0, 0.0, 1.0]

    status = crtm_onnx_predict(features, 1_c_size_t, &
                               size(features, kind=c_size_t), &
                               transmittances, &
                               size(transmittances, kind=c_size_t))
    
    if (status == 0) then
        print *, "Successfully predicted transmittances"
        print *, "First 5 channels: ", transmittances(1:5)
    else
        print *, "Prediction failed"
    end if

    call crtm_onnx_cleanup()

end program test_onnx_bridge
