!@descr: unit tests of wave_field: transforms, round trips, transposes, embedding and crop
! The forward transform equals a double-precision direct DFT (12-point grid) to 1e-5; bwd(fwd(x)) = x to 1e-6;
! the transposes of F, embed and crop hold to 1e-5 relative; crop(embed(S)) = S, and bwd(embed(S)) on the even
! fine-grid points is the native image of S to 1e-5; on the native grid (Nyquist pairs sharing a bin) embed and
! crop stay transposes; the real-field transforms equal their complex counterparts to 1e-5 on both grids.
module simple_wave_field_tester
use simple_defs,       only: dp, PI
use simple_wave_field, only: wave_field, wave_freq
use simple_test_utils, only: assert_true, assert_int, fix_test_seed
implicit none
private
public :: run_all_wave_field_tests

integer, parameter :: NFINE = 24, RNAT = 6     ! fine grid 2*nat with nat = 2*RNAT = 12

contains

    subroutine run_all_wave_field_tests()
        write(*,'(A)') '**** running all wave field tests ****'
        call fix_test_seed()
        call test_freq()
        call test_dft_and_round_trip()
        call test_transform_adjoints()
        call test_disk_transposes()
        call test_embedded_amplitude()
        call test_native_grid()
        call test_real_transforms(NFINE)
        call test_real_transforms(2*RNAT)
    end subroutine run_all_wave_field_tests

    subroutine test_freq()
        call assert_int(0,  wave_freq(1, 8),  'index 1 is the zero frequency')
        call assert_int(4,  wave_freq(5, 8),  'index n/2+1 is the Nyquist frequency, positive')
        call assert_int(-3, wave_freq(6, 8),  'indices past Nyquist are negative frequencies')
        call assert_int(-1, wave_freq(8, 8),  'the last index is frequency -1')
        call assert_int(-2, wave_freq(4, 5),  'odd grids have no Nyquist bin')
    end subroutine test_freq

    subroutine test_dft_and_round_trip()
        integer, parameter :: N = 12
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        complex          :: x(N,N)
        complex(dp)      :: ref
        real(dp)         :: err, scale
        integer          :: i, j, u, v
        call f%new(N)
        x = random_field(N)
        p => f%ptr()
        p = x
        call f%fwd
        err = 0.0_dp; scale = 0.0_dp
        do v = 1, N
            do u = 1, N
                ref = (0.0_dp, 0.0_dp)
                do j = 1, N
                    do i = 1, N
                        ref = ref + cmplx(x(i,j), kind=dp) * exp(cmplx(0.0_dp, &
                            &-2.0_dp*PI*real((i-1)*(u-1) + (j-1)*(v-1),dp)/real(N,dp), kind=dp))
                    end do
                end do
                err   = max(err, abs(cmplx(p(u,v), kind=dp) - ref))
                scale = max(scale, abs(ref))
            end do
        end do
        call assert_true(err <= 1.0e-5_dp * scale, 'the forward transform is the unnormalised DFT of sign -1')
        call f%bwd
        call assert_true(maxval(abs(p - x)) <= 1.0e-6 * maxval(abs(x)), 'bwd(fwd(x)) recovers x')
        call f%kill
        call f%kill
        call assert_int(0, f%get_n(), 'kill is idempotent and empties the field')
    end subroutine test_dft_and_round_trip

    !> Re<F x, y> = Re<x, n^2 F^-1 y>
    subroutine test_transform_adjoints()
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        complex          :: x(NFINE,NFINE), y(NFINE,NFINE), fx(NFINE,NFINE), fhy(NFINE,NFINE)
        real(dp)         :: lhs, rhs
        call f%new(NFINE)
        p => f%ptr()
        x = random_field(NFINE)
        y = random_field(NFINE)
        p = x
        call f%fwd
        fx = p
        p = y
        call f%bwd
        fhy = p * real(NFINE)**2
        lhs = redot(fx, y)
        rhs = redot(x, fhy)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), &
            &'the transpose of the forward transform is n^2 times the inverse')
        call f%kill
    end subroutine test_transform_adjoints

    !> Re<embed S, u> = Re<S, embed_t u>, Re<crop u, S> = Re<u, crop_t S>, crop(embed(S)) = S
    subroutine test_disk_transposes()
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        complex          :: s(-RNAT:RNAT,-RNAT:RNAT), t(-RNAT:RNAT,-RNAT:RNAT)
        complex          :: u(NFINE,NFINE), es(NFINE,NFINE)
        real(dp)         :: lhs, rhs
        call f%new(NFINE)
        p => f%ptr()
        s = random_disk()
        u = random_field(NFINE)
        call f%embed(s, RNAT)
        es = p
        p  = u
        call f%embed_t(t, RNAT)
        lhs = redot(es, u)
        rhs = redot_disk(s, t)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), 'embed_t is the transpose of embed')
        p = u
        call f%crop(t, RNAT)
        lhs = redot_disk(t, s)
        call f%crop_t(s, RNAT)
        rhs = redot(u, p)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), 'crop_t is the transpose of crop')
        call f%embed(s, RNAT)
        call f%crop(t, RNAT)
        call assert_true(maxval(abs(t - s)) <= 1.0e-6 * maxval(abs(s)), 'crop undoes embed on the disk')
        call assert_true(all(p(RNAT+2:NFINE-RNAT,:) == cmplx(0.,0.)), 'embed leaves the field beyond the disk zero')
        call f%kill
    end subroutine test_disk_transposes

    !> bwd(embed(S)) at the fine-grid points 2i-1 is the native image of S, the both-sign Nyquist
    !! pair of the disk folding onto one native bin
    subroutine test_embedded_amplitude()
        integer, parameter :: NAT = 2*RNAT
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        complex          :: s(-RNAT:RNAT,-RNAT:RNAT)
        complex(dp)      :: ref
        real(dp)         :: err, scale
        integer          :: i, j, h, k
        call f%new(NFINE)
        p => f%ptr()
        s = random_disk()
        call f%embed(s, RNAT)
        call f%bwd
        err = 0.0_dp; scale = 0.0_dp
        do j = 1, NAT
            do i = 1, NAT
                ref = (0.0_dp, 0.0_dp)
                do k = -RNAT, RNAT
                    do h = -RNAT, RNAT
                        if( h*h + k*k > RNAT*RNAT ) cycle
                        ref = ref + cmplx(s(h,k), kind=dp) * exp(cmplx(0.0_dp, &
                            &2.0_dp*PI*real(h*(i-1) + k*(j-1),dp)/real(NAT,dp), kind=dp))
                    end do
                end do
                ref   = ref / real(NAT*NAT,dp)
                err   = max(err, abs(cmplx(p(2*i-1,2*j-1), kind=dp) - ref))
                scale = max(scale, abs(ref))
            end do
        end do
        call assert_true(err <= 1.0e-5_dp * scale, 'the embedded spectrum is the native image at its own amplitude')
        call f%kill
    end subroutine test_embedded_amplitude

    !> on a grid of the native box itself (n = nat) the disk's Nyquist pairs share a bin: embed adds them,
    !! embed_t gives the bin to both, and the pairs stay transposes
    subroutine test_native_grid()
        integer, parameter :: N = 2*RNAT
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        complex  :: s(-RNAT:RNAT,-RNAT:RNAT), t(-RNAT:RNAT,-RNAT:RNAT), u(N,N), es(N,N)
        real(dp) :: lhs, rhs
        call f%new(N)
        p => f%ptr()
        s = random_disk()
        u = random_field(N)
        call f%embed(s, RNAT)
        es = p
        p  = u
        call f%embed_t(t, RNAT)
        lhs = redot(es, u)
        rhs = redot_disk(s, t)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), 'on the native grid embed_t is the transpose of embed')
        call f%crop(t, RNAT)
        lhs = redot_disk(t, s)
        call f%crop_t(s, RNAT)
        rhs = redot(u, p)
        call assert_true(abs(lhs - rhs) <= 1.0e-5_dp * max(abs(lhs), abs(rhs)), 'on the native grid crop_t is the transpose of crop')
        call f%kill
    end subroutine test_native_grid

    !> fwd_real_crop is crop(fwd) of the real field and crop_t_real is n^2 Re(bwd(crop_t)), on the grid of n
    !! points (n = 2 RNAT: the native grid, the disk's Nyquist pairs sharing a bin)
    subroutine test_real_transforms( n )
        integer, intent(in) :: n
        type(wave_field) :: f
        complex, pointer :: p(:,:)
        real,    pointer :: r(:,:)
        complex :: s(-RNAT:RNAT,-RNAT:RNAT), t1(-RNAT:RNAT,-RNAT:RNAT), t2(-RNAT:RNAT,-RNAT:RNAT), u(n,n)
        real    :: rref(n,n)
        call f%new(n, with_real=.true.)
        p => f%ptr()
        r => f%rptr()
        u = random_field(n)
        r = real(u)
        call f%fwd_real_crop(t1, RNAT)
        p = cmplx(real(u), 0.)
        call f%fwd
        call f%crop(t2, RNAT)
        call assert_true(maxval(abs(t1 - t2)) <= 1.0e-5 * maxval(abs(t2)), 'fwd_real_crop is crop(fwd) of the real field')
        s = random_disk()
        call f%crop_t_real(s, RNAT)
        call f%crop_t(s, RNAT)
        call f%bwd
        rref = real(n)**2 * real(p)
        call assert_true(maxval(abs(r - rref)) <= 1.0e-5 * maxval(abs(rref)), 'crop_t_real is n^2 Re(bwd(crop_t))')
        call f%kill
    end subroutine test_real_transforms

    function random_field( n ) result( x )
        integer, intent(in) :: n
        complex :: x(n,n)
        real    :: re(n,n), im(n,n)
        call random_number(re)
        call random_number(im)
        x = cmplx(re - 0.5, im - 0.5)
    end function random_field

    function random_disk() result( s )
        complex :: s(-RNAT:RNAT,-RNAT:RNAT)
        real    :: re(-RNAT:RNAT,-RNAT:RNAT), im(-RNAT:RNAT,-RNAT:RNAT)
        integer :: h, k
        call random_number(re)
        call random_number(im)
        s = cmplx(re - 0.5, im - 0.5)
        do k = -RNAT, RNAT
            do h = -RNAT, RNAT
                if( h*h + k*k > RNAT*RNAT ) s(h,k) = cmplx(0.,0.)
            end do
        end do
    end function random_disk

    real(dp) function redot( a, b )
        complex, intent(in) :: a(:,:), b(:,:)
        redot = sum(real(a,dp)*real(b,dp) + real(aimag(a),dp)*real(aimag(b),dp))
    end function redot

    real(dp) function redot_disk( a, b )
        complex, intent(in) :: a(-RNAT:RNAT,-RNAT:RNAT), b(-RNAT:RNAT,-RNAT:RNAT)
        redot_disk = sum(real(a,dp)*real(b,dp) + real(aimag(a),dp)*real(aimag(b),dp))
    end function redot_disk

end module simple_wave_field_tester
