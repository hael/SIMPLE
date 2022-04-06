
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
            use simple_testfuns_constants, only: target_img, obj_img, ker_img, ker_der_img
            use simple_math,               only: butterworth_kernel
            use simple_image,              only: image
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(in)    :: x(d)
            real                    :: r
            real                          :: rmat_target(202,202,1), rmat_ker(202,202,1), rmat_ker_der(202,202,1)
            complex(kind=c_float_complex) :: cmat_conv(202,202,1)  , cmat_obj(202,202,1)

            call butterworth_kernel(rmat_ker, rmat_ker_der, 202, 8, x(1))  ! WARNING: fix the constants here
            call ker_img%set_rmat(rmat_ker, .false.)
            call ker_der_img%set_rmat(rmat_ker_der, .false.)
            
            call target_img%get_rmat_sub(rmat_target)
            call obj_img%fft()
            call obj_img%get_cmat_sub(cmat_obj)
            call ker_img%fft()
            call ker_img%get_cmat_sub(cmat_conv)
            cmat_conv = cmat_conv*cmat_obj
            call ker_img%set_cmat(cmat_conv)
            call ker_img%ifft()
            call ker_img%get_rmat_sub(rmat_ker)
            call obj_img%ifft()
            rmat_ker = rmat_ker*202*202;

            r = sum(abs(rmat_ker - rmat_target)**2)
            !write(*, *) 'cost: ', x(1), r
        end function

        subroutine butterworth_gcost( fun_self, x, grad, d )
            use simple_testfuns_constants, only: target_img, obj_img, ker_img, ker_der_img
            use simple_math,               only: butterworth_kernel
            use simple_image,              only: image
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(inout) :: x(d)
            real,     intent(out)   :: grad(d)
            real,                          pointer :: rmat_target(:,:,:), rmat_ker(:,:,:), rmat_ker_der(:,:,:)
            complex(kind=c_float_complex), pointer :: cmat_conv(:,:,:)  , cmat_obj(:,:,:), cmat_der_conv(:,:,:)

            allocate(rmat_target(  202, 202, 1))
            allocate(rmat_ker(     202, 202, 1))
            allocate(rmat_ker_der( 202, 202, 1))
            allocate(cmat_conv(    202, 202, 1))
            allocate(cmat_obj(     202, 202, 1))
            allocate(cmat_der_conv(202, 202, 1))
            
            call butterworth_kernel(rmat_ker, rmat_ker_der, 202, 8, x(1))  ! WARNING: fix the constants here
            call ker_img%set_rmat(rmat_ker, .false.)
            call ker_der_img%set_rmat(rmat_ker_der, .false.)
            
            call target_img%get_rmat_sub(rmat_target)
            call obj_img%fft()
            call obj_img%get_cmat_sub(cmat_obj)
            call ker_img%fft()
            call ker_img%get_cmat_sub(cmat_conv)
            

            call ker_der_img%fft()
            call ker_der_img%get_cmat_sub(cmat_der_conv)
            cmat_conv     = cmat_conv    *cmat_obj
            cmat_der_conv = cmat_der_conv*cmat_obj
            call ker_img%set_cmat(cmat_conv)
            call ker_der_img%set_cmat(cmat_der_conv)
            call ker_der_img%ifft()
            call ker_img%ifft()       
            call obj_img%ifft()
            call ker_img%get_rmat_sub(rmat_ker)
            call ker_der_img%get_rmat_sub(rmat_ker_der)
            rmat_ker     = rmat_ker*202*202;
            rmat_ker_der = rmat_ker_der*202*202;

            grad(1) = 2*sum((rmat_ker - rmat_target)*rmat_ker_der)
            !write(*, *) 'grad: ', x(1), grad(1)
        end subroutine
end module
    