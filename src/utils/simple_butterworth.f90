
module simple_butterworth
    use simple_defs
    use simple_image, only: image
    use simple_math,  only: hyp
    implicit none
    type(image) :: target_img, obj_img, ker_img, ker_der_img
    
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
        ! Compute the value of the Butterworth transfer function of order n(th)
        ! at a given frequency s, with the cut-off frequency fc
        ! SOURCE :
        ! https://en.wikipedia.org/wiki/Butterworth_filter
        function butterworth(s, n, fc) result(val)
            real   , intent(in)  :: s
            integer, intent(in)  :: n
            real   , intent(in)  :: fc
            real                 :: val(2)

            real, parameter :: an(9) = (/ 1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371, 5.1258, 1./)
            complex :: Bn, dBn, Kn, dKn  ! Normalized Butterworth polynomial, its derivative and its reciprocal
            complex :: j = (0, 1)        ! Complex identity: j = sqrt(-1)
            complex :: js                ! frequency is multiplied by the complex identity j
            integer :: k
            val = [0., 0.]
            Bn  = (0., 0.)
            dBn = (0., 0.)
            js  = j*s/fc
            do k = 0, n
                Bn  = Bn  +   an(k+1)*js**k
                dBn = dBn + k*an(k+1)*js**k
            end do
            dBn = -dBn/fc
            Kn  = 1/Bn
            dKn = -dBn/Bn**2
            val(1) = sqrt(real(Kn)**2 + aimag(Kn)**2)
            val(2) = real( Kn*conjg(dKn) )/val(1)
        end function butterworth

        ! Compute the Butterworth kernel of the order n-th of width w
        ! with the cut-off frequency fc
        ! https://en.wikipedia.org/wiki/Butterworth_filter
        subroutine butterworth_kernel(ker, ker_der, w, n, fc)
            real,    intent(inout) :: ker    (:, :, :)    ! assuming 2D kernel for now!!!
            real,    intent(inout) :: ker_der(:, :, :)    ! assuming 2D kernel for now!!!
            integer, intent(in)    :: w
            integer, intent(in)    :: n
            real   , intent(in)    :: fc
            integer :: k, l, half_w
            real    :: freq_val, val(2)    ! current frequency value

            freq_val = 0
            half_w   = int(w/2)
            do k = 1, w
                do l = 1, w
                    freq_val = hyp(real(k-half_w), real(l-half_w))

                    ! compute the value of Butterworth transfer function at current frequency value
                    val = butterworth(freq_val, n, fc)
                    ker(k,l,1)     = val(1)
                    ker_der(k,l,1) = val(2)
                end do
            end do
        end subroutine butterworth_kernel

        function butterworth_cost( fun_self, x, d ) result( r )
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
    