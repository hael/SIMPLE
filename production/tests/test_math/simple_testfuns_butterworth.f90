
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
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(in)    :: x(d)
            real                    :: r
            real                         , pointer :: rmat_target(:,:,:), rmat_ker(:,:,:), rmat_ker_der(:,:,:)
            complex(kind=c_float_complex), pointer :: cmat_conv(:,:,:)  , cmat_obj(:,:,:), cmat_der_conv(:,:,:)
            integer :: ldim(3)

            ldim = ker_img%get_ldim()

            allocate(rmat_target(  ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der( ldim(1), ldim(2), ldim(3)))
            allocate(cmat_conv(    ldim(1), ldim(2), ldim(3)))
            allocate(cmat_obj(     ldim(1), ldim(2), ldim(3)))

            call ker_img%get_rmat_ptr(rmat_ker)
            call ker_der_img%get_rmat_ptr(rmat_ker_der)
            call butterworth_kernel(rmat_ker, rmat_ker_der, ldim(1), 8, x(1))  ! WARNING: fix the constants here
            
            call target_img%get_rmat_ptr(rmat_target)
            call obj_img%fft()
            call obj_img%get_cmat_ptr(cmat_obj)
            call ker_img%fft()
            call ker_img%get_cmat_ptr(cmat_conv)
            call ker_der_img%fft()
            call ker_der_img%get_cmat_ptr(cmat_der_conv)
            cmat_conv     = cmat_conv    *cmat_obj
            cmat_der_conv = cmat_der_conv*cmat_obj
            call ker_img%ifft()
            call ker_der_img%ifft()
            rmat_ker     = rmat_ker    *ldim(1)*ldim(2)   ! TODO: check why scaling here
            rmat_ker_der = rmat_ker_der*ldim(1)*ldim(2)   ! TODO: check why scaling here

            r = sum(abs(rmat_ker - rmat_target)**2)
            write(*, *) 'x = ', x(1), 'cost = ', r
        end function

        subroutine butterworth_gcost( fun_self, x, grad, d )
            use simple_testfuns_constants, only: target_img, obj_img, ker_img, ker_der_img
            use simple_math,               only: butterworth_kernel
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(inout) :: x(d)
            real,     intent(out)   :: grad(d)
            real,     pointer       :: rmat_target(:,:,:), rmat_ker(:,:,:), rmat_ker_der(:,:,:)
            integer                 :: ldim(3)

            ldim = ker_img%get_ldim()

            allocate(rmat_target(  ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der( ldim(1), ldim(2), ldim(3)))
            
            call ker_img    %get_rmat_ptr(rmat_ker)
            call ker_der_img%get_rmat_ptr(rmat_ker_der)
            call target_img %get_rmat_ptr(rmat_target)
        
            grad(1) = 2*sum((rmat_ker - rmat_target)*rmat_ker_der)
        end subroutine
end module
    