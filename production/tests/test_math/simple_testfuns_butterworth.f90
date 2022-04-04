
module simple_testfuns_butterworth
    use simple_defs
    implicit none
    
    !>  \brief  defines the test function interface
    abstract interface
        function testfun_butterworth( fun_self, vec, D ) result( cost )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: D
            real,     intent(in)    :: vec(D)
            real                    :: cost
        end function
    end interface
    
    contains
        function butterworth_cost( fun_self, x, d ) result( r )
            use simple_testfuns_constants, only: target_img, obj_img, ker_img
            use simple_math,               only: butterworth, butterworth_kernel
            use simple_image,              only: image
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(in)    :: x(d)
            real                    :: r
            real,                          pointer :: rmat_target(:,:,:), rmat_ker(:,:,:)
            complex(kind=c_float_complex), pointer :: cmat_conv(:,:,:)  , cmat_obj(:,:,:)

            call ker_img%get_rmat_ptr(rmat_ker)
            call butterworth_kernel(rmat_ker, 202, 8, x(1))  ! WARNING: fix the constants here
            call ker_img%set_rmat(rmat_ker, .false.)
            call target_img%get_rmat_ptr(rmat_target)
            call obj_img%get_cmat_ptr(cmat_obj)
            call ker_img%fft()
            call ker_img%get_cmat_ptr(cmat_conv)
            cmat_conv = cmat_conv*cmat_obj
            call ker_img%ifft()
            call ker_img%get_rmat_ptr(rmat_ker)
            r = sum(abs(rmat_ker - rmat_target))
        end function
end module
    