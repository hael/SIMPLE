! Butterworth kernel
module simple_butterworth
include 'simple_lib.f08'
implicit none
#include "simple_local_flags.inc"

interface butterworth_filter
    module procedure butterworth_filter_1, butterworth_filter_2, butterworth_filter_3, butterworth_filter_4
end interface butterworth_filter

integer, parameter :: BW_ORDER = 8

contains

    ! Compute the value of the Butterworth transfer function of order n(th)
    ! at a given frequency s, with the cut-off frequency fc
    ! SOURCE :
    ! https://en.wikipedia.org/wiki/Butterworth_filter
    pure function butterworth(s, n, fc) result(val)
        real   , intent(in)  :: s
        integer, intent(in)  :: n
        real   , intent(in)  :: fc
        real                 :: val
        real,    parameter :: AN(11,10) = reshape((/ 1., 1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 1.4142,  1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 2.    ,  2.    ,  1.    ,  0.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 2.6131,  3.4142,  2.6131,  1.    ,  0.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 3.2361,  5.2361,  5.2361,  3.2361,  1.    ,  0.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 3.8637,  7.4641,  9.1416,  7.4641,  3.8637,  1.    ,  0.    ,  0.    , 0.    , 0.,&
                                                    &1., 4.4940, 10.0978, 14.5918, 14.5918, 10.0978,  4.4940,  1.    ,  0.    , 0.    , 0.,&
                                                    &1., 5.1258, 13.1371, 21.8462, 25.6884, 21.8462, 13.1371,  5.1258,  1.    , 0.    , 0.,&
                                                    &1., 5.7588, 16.5817, 31.1634, 41.9864, 41.9864, 31.1634, 16.5817,  5.7588, 1.    , 0.,&
                                                    &1., 6.3925, 20.4317, 42.8021, 64.8824, 74.2334, 64.8824, 42.8021, 20.4317, 6.3925, 1. /),&
                                                    &(/11,10/))
        complex, parameter :: J = (0, 1) ! Complex identity: j = sqrt(-1)
        complex :: Bn, Kn                ! Normalized Butterworth polynomial, its derivative and its reciprocal
        complex :: js                    ! frequency is multiplied by the complex identity j
        integer :: k
        Bn  = (0., 0.)
        if (s/fc < 100) then
            js  = J*s/fc
            do k = 0, n
                Bn  = Bn + AN(k+1,n)*js**k
            end do
            Kn  = 1/Bn
            val = cabs(Kn)
        else
            val = epsilon(val)
        endif
        val = min(1.,max(0.,val))
    end function butterworth

    subroutine butterworth_filter_1(img, cutoff_find, cur_fil)
        use simple_image, only: image
        class(image), intent(inout) :: img
        integer,      intent(in)    :: cutoff_find
        real,         intent(inout) :: cur_fil(:)
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
        call img%apply_filter(cur_fil)
    end subroutine butterworth_filter_1

    subroutine butterworth_filter_2(cutoff_find, cur_fil)
        integer,      intent(in)    :: cutoff_find
        real,         intent(inout) :: cur_fil(:)
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
    end subroutine butterworth_filter_2

    !< notch filter squeezed between two cut_off frequencies
    subroutine butterworth_filter_3(c1, c2, cur_fil)
        integer,      intent(in)    :: c1, c2           ! two cut-off frequencies
        real,         intent(inout) :: cur_fil(:)
        integer :: freq_val
        do freq_val = 1, size(cur_fil)
            cur_fil(freq_val) =        butterworth(real(freq_val-1), BW_ORDER, real(c2)) * &
                                &( 1 - butterworth(real(freq_val-1), BW_ORDER, real(c1)) )
        enddo
    end subroutine butterworth_filter_3

    subroutine butterworth_filter_4(cutoff_find, kfromto, cur_fil)
        integer,      intent(in)    :: cutoff_find
        integer,      intent(in)    :: kfromto(2)
        real,         intent(inout) :: cur_fil(kfromto(1):kfromto(2))
        integer :: freq_val
        do freq_val = kfromto(1), kfromto(2)
            cur_fil(freq_val) = butterworth(real(freq_val-1), BW_ORDER, real(cutoff_find))
        enddo
    end subroutine butterworth_filter_4

end module simple_butterworth
