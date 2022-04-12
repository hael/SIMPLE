
module simple_butterworth
    use simple_defs
    use simple_image, only: image
    use simple_math,  only: hyp
    implicit none
    type(image) :: odd_img, even_img, ker_odd_img, ker_even_img, ker_der_odd_img, ker_der_even_img
    
    !>  \brief  defines the function interface
    abstract interface
        function fun_butterworth( fun_self, vec, D ) result( cost )
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
            val(2) = real( Kn*conjg(dKn) )/Kn
        end function butterworth

        ! Compute the Butterworth kernel of the order n-th of width w
        ! with the cut-off frequency fc
        ! https://en.wikipedia.org/wiki/Butterworth_filter
        subroutine butterworth_kernel(ker, ker_der, w, n, fc)
            real,    intent(inout) :: ker    (:, :, :)
            real,    intent(inout) :: ker_der(:, :, :)
            integer, intent(in)    :: w
            integer, intent(in)    :: n
            real   , intent(in)    :: fc
            integer :: k, l, j, half_w
            real    :: freq_val, val(2)    ! current frequency value

            freq_val = 0
            half_w   = int(w/2)
            do k = 1, w
                do l = 1, w
                    do j = 1, w
                        freq_val = hyp(real(k-half_w), real(l-half_w), real(j-half_w))

                        ! compute the value of Butterworth transfer function at current frequency value
                        val = butterworth(freq_val, n, fc)
                        ker(k,l,j)     = val(1)
                        ker_der(k,l,j) = val(2)
                    end do
                end do
            end do
        end subroutine butterworth_kernel

        ! use even as the target and odd as the object to be convolved with Butterworth Kernel
        function butterworth_cost( fun_self, x, d ) result( r )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(in)    :: x(d)
            real                    :: r
            real                         , pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker(:,:,:),  rmat_ker_der(:,:,:), orig_ker(:,:,:), orig_ker_der(:,:,:)
            complex(kind=c_float_complex), pointer :: cmat_odd(:,:,:), cmat_even(:,:,:), cmat_conv(:,:,:), cmat_der_conv(:,:,:)
            integer :: ldim(3)

            ldim = ker_odd_img%get_ldim()

            allocate(rmat_odd(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_even(    ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der( ldim(1), ldim(2), ldim(3)))
            allocate(orig_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(orig_ker_der( ldim(1), ldim(2), ldim(3)))
            allocate(cmat_odd(     ldim(1), ldim(2), ldim(3)))
            allocate(cmat_even(    ldim(1), ldim(2), ldim(3)))

            allocate(cmat_conv(    int(ldim(1)/2)+1, ldim(2), ldim(3))) ! FT change the first dimenstion into n/2+1
            allocate(cmat_der_conv(int(ldim(1)/2)+1, ldim(2), ldim(3)))

            call butterworth_kernel(orig_ker, orig_ker_der, ldim(1), 8, x(1))  ! WARNING: fix the constants here
            
            call odd_img%get_rmat_ptr(rmat_odd)
            call odd_img%fft()
            call odd_img%get_cmat_ptr(cmat_odd)
            call even_img%get_rmat_ptr(rmat_even)
            call even_img%fft()
            call even_img%get_cmat_ptr(cmat_even)

            ! compute cost of ||ker*odd - even||
            call ker_odd_img%get_rmat_ptr(rmat_ker)
            call ker_odd_img%set_rmat(orig_ker, .false.)
            call ker_odd_img%fft()
            call ker_odd_img%get_cmat_sub(cmat_conv)        ! no use of pointer since cmat_conv is used again below

            call ker_der_odd_img%get_rmat_ptr(rmat_ker_der)
            call ker_der_odd_img%set_rmat(orig_ker_der, .false.)
            call ker_der_odd_img%fft()
            call ker_der_odd_img%get_cmat_sub(cmat_der_conv)
            cmat_conv     = cmat_conv    *cmat_odd
            cmat_der_conv = cmat_der_conv*cmat_odd
            call ker_odd_img%set_cmat(cmat_conv)
            call ker_der_odd_img%set_cmat(cmat_der_conv)
            call ker_odd_img%ifft()
            call ker_der_odd_img%ifft()
            call odd_img%ifft()
            call even_img%ifft()
            rmat_ker     = rmat_ker    *ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here
            rmat_ker_der = rmat_ker_der*ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here

            r = sum(abs(rmat_ker - rmat_even)**2)
            write(*, *) 'x = ', x(1), 'cost = ', r
        end function

        subroutine butterworth_gcost( fun_self, x, grad, d )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(inout) :: x(d)
            real,     intent(out)   :: grad(d)
            real,     pointer       :: rmat_even(:,:,:), rmat_ker_odd(:,:,:), rmat_ker_der_odd(:,:,:)
            integer                 :: ldim(3)

            ldim = ker_odd_img%get_ldim()
            allocate(rmat_even(ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_odd( ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der_odd( ldim(1), ldim(2), ldim(3)))
            
            call even_img%get_rmat_ptr(rmat_even)
            call ker_odd_img%get_rmat_ptr(rmat_ker_odd)
            call ker_der_odd_img%get_rmat_ptr(rmat_ker_der_odd)
            
            grad(1) = 2*sum((rmat_ker_odd - rmat_even)*rmat_ker_der_odd)
            write(*, *) 'x = ', x(1), 'dcost = ', grad(1)
        end subroutine

        function butterworth_evenodd_cost( fun_self, x, d ) result( r )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(in)    :: x(d)
            real                    :: r
            real                         , pointer :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker(:,:,:),  rmat_ker_der(:,:,:), orig_ker(:,:,:), orig_ker_der(:,:,:)
            complex(kind=c_float_complex), pointer :: cmat_odd(:,:,:), cmat_even(:,:,:), cmat_conv(:,:,:), cmat_der_conv(:,:,:)
            integer :: ldim(3)

            ldim = ker_odd_img%get_ldim()

            allocate(rmat_odd(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_even(    ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der( ldim(1), ldim(2), ldim(3)))
            allocate(orig_ker(     ldim(1), ldim(2), ldim(3)))
            allocate(orig_ker_der( ldim(1), ldim(2), ldim(3)))
            allocate(cmat_odd(     ldim(1), ldim(2), ldim(3)))
            allocate(cmat_even(    ldim(1), ldim(2), ldim(3)))

            allocate(cmat_conv(    int(ldim(1)/2)+1, ldim(2), ldim(3))) ! FT change the first dimension into n/2+1
            allocate(cmat_der_conv(int(ldim(1)/2)+1, ldim(2), ldim(3)))

            call butterworth_kernel(orig_ker, orig_ker_der, ldim(1), 8, x(1))  ! WARNING: fix the constants here
            
            call odd_img%get_rmat_ptr(rmat_odd)
            call odd_img%fft()
            call odd_img%get_cmat_ptr(cmat_odd)
            call even_img%get_rmat_ptr(rmat_even)
            call even_img%fft()
            call even_img%get_cmat_ptr(cmat_even)

            ! compute cost of ||ker*odd - even||
            call ker_odd_img%get_rmat_ptr(rmat_ker)
            call ker_odd_img%set_rmat(orig_ker, .false.)
            call ker_odd_img%fft()
            call ker_odd_img%get_cmat_sub(cmat_conv)        ! no use of pointer since cmat_conv is used again below

            call ker_der_odd_img%get_rmat_ptr(rmat_ker_der)
            call ker_der_odd_img%set_rmat(orig_ker_der, .false.)
            call ker_der_odd_img%fft()
            call ker_der_odd_img%get_cmat_sub(cmat_der_conv)
            cmat_conv     = cmat_conv    *cmat_odd
            cmat_der_conv = cmat_der_conv*cmat_odd
            call ker_odd_img%set_cmat(cmat_conv)
            call ker_der_odd_img%set_cmat(cmat_der_conv)
            call ker_odd_img%ifft()
            call ker_der_odd_img%ifft()
            call even_img%ifft()
            rmat_ker     = rmat_ker    *ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here
            rmat_ker_der = rmat_ker_der*ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here

            r = sum(abs(rmat_ker - rmat_even)**2)

            ! compute cost of ||ker*even - odd||
            call ker_even_img%get_rmat_ptr(rmat_ker)
            call ker_even_img%set_rmat(orig_ker, .false.)
            call ker_even_img%fft()
            call ker_even_img%get_cmat_ptr(cmat_conv)        ! it's safe to get the pointer here since cmat_conv is not used after

            call ker_der_even_img%get_rmat_ptr(rmat_ker_der)
            call ker_der_even_img%set_rmat(orig_ker_der, .false.)
            call ker_der_even_img%fft()
            call ker_der_even_img%get_cmat_ptr(cmat_der_conv)
            cmat_conv     = cmat_conv    *cmat_even
            cmat_der_conv = cmat_der_conv*cmat_even
            call ker_even_img%ifft()
            call ker_der_even_img%ifft()
            call odd_img%ifft()
            rmat_ker     = rmat_ker    *ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here
            rmat_ker_der = rmat_ker_der*ldim(1)*ldim(2)*ldim(3)   ! TODO: check why scaling here

            r = r + sum(abs(rmat_ker - rmat_odd)**2)
            write(*, *) 'x = ', x(1), 'cost = ', r
        end function

        subroutine butterworth_evenodd_gcost( fun_self, x, grad, d )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: d
            real,     intent(inout) :: x(d)
            real,     intent(out)   :: grad(d)
            real,     pointer       :: rmat_odd(:,:,:), rmat_even(:,:,:), rmat_ker_odd(:,:,:), rmat_ker_even(:,:,:), rmat_ker_der_odd(:,:,:), rmat_ker_der_even(:,:,:)
            integer                 :: ldim(3)

            ldim = ker_odd_img%get_ldim()

            allocate(rmat_odd( ldim(1), ldim(2), ldim(3)))
            allocate(rmat_even(ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_odd( ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_even(ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der_odd( ldim(1), ldim(2), ldim(3)))
            allocate(rmat_ker_der_even(ldim(1), ldim(2), ldim(3)))
            
            call odd_img %get_rmat_ptr(rmat_odd)
            call even_img%get_rmat_ptr(rmat_even)

            call ker_odd_img %get_rmat_ptr(rmat_ker_odd)
            call ker_even_img%get_rmat_ptr(rmat_ker_even)

            call ker_der_odd_img %get_rmat_ptr(rmat_ker_der_odd)
            call ker_der_even_img%get_rmat_ptr(rmat_ker_der_even)
            
        
            grad(1) = 2*sum((rmat_ker_odd - rmat_even)*rmat_ker_der_odd) + 2*sum((rmat_ker_even - rmat_odd)*rmat_ker_der_even)
            write(*, *) 'x = ', x(1), 'dcost = ', grad(1)
        end subroutine
end module
    