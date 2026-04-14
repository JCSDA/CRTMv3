module crtm_onnx_interface
    use iso_c_binding
    implicit none

    interface
        function crtm_onnx_init(model_path) bind(c, name="crtm_onnx_init")
            import :: c_int, c_char
            character(kind=c_char), intent(in) :: model_path(*)
            integer(c_int) :: crtm_onnx_init
        end function crtm_onnx_init

        function crtm_onnx_predict(input_data, input_dim, output_data, output_dim) &
                bind(c, name="crtm_onnx_predict")
            import :: c_int, c_float, c_size_t
            real(c_float), intent(in) :: input_data(*)
            integer(c_size_t), value :: input_dim
            real(c_float), intent(out) :: output_data(*)
            integer(c_size_t), value :: output_dim
            integer(c_int) :: crtm_onnx_predict
        end function crtm_onnx_predict

        subroutine crtm_onnx_cleanup() bind(c, name="crtm_onnx_cleanup")
        end subroutine crtm_onnx_cleanup
    end interface

end module crtm_onnx_interface

program test_onnx_bridge
    use crtm_onnx_interface
    use iso_c_binding
    implicit none

    integer(c_int) :: status
    real(c_float), dimension(6) :: features  ! P, T, H2O, CO2, O3, Angle
    real(c_float), dimension(8461) :: transmittances
    character(kind=c_char), dimension(:), allocatable :: model_path_c
    character(len=255) :: model_path_fortran
    integer :: i

    model_path_fortran = "/home/ben/CRTM/LBL/ml_emulator/output_iasi_resnet_dynamic_v2/model.onnx"
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

    status = crtm_onnx_predict(features, size(features, kind=c_size_t), &
                               transmittances, size(transmittances, kind=c_size_t))
    
    if (status == 0) then
        print *, "Successfully predicted transmittances"
        print *, "First 5 channels: ", transmittances(1:5)
    else
        print *, "Prediction failed"
    end if

    call crtm_onnx_cleanup()

end program test_onnx_bridge
