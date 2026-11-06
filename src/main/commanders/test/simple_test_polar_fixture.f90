!@descr: polar-ring fixtures shared by the shift-map testers: the calculator's shift arguments, a half ring rotated on the grid, a ladder of shift magnitudes
module simple_test_polar_fixture
use simple_defs, only: sp, dp, DTWOPI
implicit none
private
public :: shift_args, rotate_half_ring, magnitudes

contains

    !> the shift arguments of the polar shift matrix, as the calculator forms them (argtransf):
    !! ax = k sin(theta_p) 2 pi / box, ay = -k cos(theta_p) 2 pi / box, theta_p = (p-1) 2 pi / nrots
    pure subroutine shift_args( pftsz, kfromto, box, ax, ay )
        integer,  intent(in)  :: pftsz, kfromto(2), box
        real(dp), intent(out) :: ax(pftsz,kfromto(1):kfromto(2)), ay(pftsz,kfromto(1):kfromto(2))
        real(dp) :: theta
        integer  :: p, k
        do p = 1, pftsz
            theta = real(p-1,dp) * DTWOPI / real(2*pftsz,dp)
            do k = kfromto(1), kfromto(2)
                ax(p,k) =  real(k,dp) * sin(theta) * DTWOPI / real(box,dp)
                ay(p,k) = -real(k,dp) * cos(theta) * DTWOPI / real(box,dp)
            enddo
        enddo
    end subroutine shift_args

    !> a half ring rotated by grid rotation irot, as rotate_ref_8 does it: entry p of the result is entry
    !! p - irot + 1 (modulo the turn) of the Friedel-extended ring, whose second half is the conjugate of the first
    pure subroutine rotate_half_ring( pftsz, nk, pft, irot, pft_rot )
        integer,     intent(in)  :: pftsz, nk, irot
        complex(sp), intent(in)  :: pft(pftsz,nk)
        complex(sp), intent(out) :: pft_rot(pftsz,nk)
        integer :: p, q
        do p = 1, pftsz
            q = modulo(p - irot, 2*pftsz) + 1
            if( q <= pftsz )then
                pft_rot(p,:) = pft(q,:)
            else
                pft_rot(p,:) = conjg(pft(q-pftsz,:))
            endif
        enddo
    end subroutine rotate_half_ring

    !> the magnitudes 0, step, 2 step, ... up to radius
    subroutine magnitudes( radius, step, rhos )
        real(dp),              intent(in)  :: radius, step
        real(dp), allocatable, intent(out) :: rhos(:)
        integer :: n, i
        n = floor(radius / step) + 1
        allocate(rhos(n))
        do i = 1, n
            rhos(i) = real(i-1,dp) * step
        enddo
    end subroutine magnitudes

end module simple_test_polar_fixture
