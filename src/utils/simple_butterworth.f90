
module simple_butterworth
!$ use omp_lib
!$ use omp_lib_kinds
use simple_defs
use simple_image, only: image
use simple_math,  only: hyp
implicit none

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
        real,    parameter :: AN(9) = (/ 1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371, 5.1258, 1./)
        complex, parameter :: J = (0, 1) ! Complex identity: j = sqrt(-1)
        complex :: Bn, dBn, Kn, dKn      ! Normalized Butterworth polynomial, its derivative and its reciprocal
        complex :: js                    ! frequency is multiplied by the complex identity j
        integer :: k
        val = [0., 0.]
        Bn  = (0., 0.)
        dBn = (0., 0.)
        if (s/fc < 100) then
            js  = j*s/fc
            do k = 0, n
                Bn  = Bn  +   AN(k+1)*js**k
                dBn = dBn + k*AN(k+1)*js**k
            end do
            dBn = -dBn/fc
            Kn  = 1/Bn
            dKn = -dBn/Bn/Bn
            val(1) = sqrt(real(Kn)**2 + aimag(Kn)**2)
            val(2) = real( Kn*conjg(dKn) )/val(1)
        else
            val(1) = epsilon(val(1))
            val(2) = epsilon(val(2))
        endif
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
        integer :: k, l, j, half_w, ldim3
        real    :: freq_val, val(2) ! current frequency value
        freq_val = 0
        half_w   = int(w/2)
        ldim3    = size(ker, 3)
        ! loop over pixels
        if( ldim3 == 1 )then
            !$omp parallel do collapse(2) default(shared) private(l,k,freq_val,val) schedule(static) proc_bind(close)
            do k = 1, w
                do l = 1, w
                    freq_val = hyp(real(k-half_w), real(l-half_w))
                    ! compute the value of Butterworth transfer function at current frequency value
                    val            = butterworth(freq_val, n, fc)
                    ker(k,l,1)     = val(1)
                    ker_der(k,l,1) = val(2)
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) default(shared) private(l,j,k,freq_val,val) schedule(static) proc_bind(close)
            do k = 1, w
                do l = 1, w
                    do j = 1, ldim3
                        freq_val = hyp(real(k-half_w), real(l-half_w), real(j-half_w))
                        ! compute the value of Butterworth transfer function at current frequency value
                        val            = butterworth(freq_val, n, fc)
                        ker(k,l,j)     = val(1)
                        ker_der(k,l,j) = val(2)
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    end subroutine butterworth_kernel

end module simple_butterworth
    