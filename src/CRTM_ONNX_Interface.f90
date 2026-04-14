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
