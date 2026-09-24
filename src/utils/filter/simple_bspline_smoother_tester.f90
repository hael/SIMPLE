!@descr: unit tests for the quadratic B-spline smoother (simple_bspline_smoother): closed-form transfer function in 2D and 3D
! Replaces the in-module test_bspline_smoother and test_bspline_smoother_3d (linearity, reproducibility
! and monotone smoothing only; two MRC files left behind). The smoother is a Tikhonov solve whose
! transfer function is H = B**2 / (B**2 + lambda*R), with B the DFT of the sampled quadratic B-spline
! (3/4 + cos(w)/4 per dimension) and R that of its gradient Gram kernel, a0 (11/20, 13/60, 1/120) in
! one dimension times a2 (1, -1/3, -1/6) in another, summed over the dimensions. A cosine comes out
! scaled by H at its frequency; a numpy emulation of the code agrees with the closed form to 2e-7,
! for square and non-square (micrograph) boxes and volumes.
module simple_bspline_smoother_tester
use simple_test_utils
use simple_defs
use simple_string_utils,     only: int2str, real2str
use simple_image,            only: image
use simple_bspline_smoother, only: bspline_smoother
implicit none
private
public :: run_all_bspline_smoother_tests

real,    parameter :: LAMBDAS(2) = [0.2, 5.]
real(dp), parameter :: TWOPI_DP = 2.d0 * acos(-1.d0)

contains

    subroutine run_all_bspline_smoother_tests()
        write(*,'(A)') '**** running all bspline_smoother tests ****'
        call test_transfer_2d()
        call test_transfer_3d()
        call test_ft_state()
    end subroutine run_all_bspline_smoother_tests

    real(dp) function beta( w )
        real(dp), intent(in) :: w
        beta = 0.75d0 + 0.25d0 * cos(w)
    end function beta

    real(dp) function a0( w )
        real(dp), intent(in) :: w
        a0 = 11.d0/20.d0 + 13.d0/30.d0 * cos(w) + 1.d0/60.d0 * cos(2.d0*w)
    end function a0

    real(dp) function a2( w )
        real(dp), intent(in) :: w
        a2 = 1.d0 - 2.d0/3.d0 * cos(w) - 1.d0/3.d0 * cos(2.d0*w)
    end function a2

    ! transfer function at frequency k of an n box (k(3) = 0 and n(3) = 1 in 2D)
    real(dp) function transfer( k, n, lambda )
        integer, intent(in) :: k(3), n(3)
        real,    intent(in) :: lambda
        real(dp) :: w(3), b, r
        w = TWOPI_DP * real(k, dp) / real(n, dp)
        if( n(3) == 1 )then
            b = beta(w(1)) * beta(w(2))
            r = a0(w(1)) * a2(w(2)) + a2(w(1)) * a0(w(2))
        else
            b = beta(w(1)) * beta(w(2)) * beta(w(3))
            r = a2(w(1)) * a0(w(2)) * a0(w(3)) + a0(w(1)) * a2(w(2)) * a0(w(3)) + a0(w(1)) * a0(w(2)) * a2(w(3))
        endif
        transfer = b * b / (b * b + real(lambda, dp) * r)
    end function transfer

    ! cos(2 pi (k . (x-1)) / n) on the pixel grid
    subroutine make_cosine( img, n, k )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: n(3), k(3)
        real, allocatable :: rmat(:,:,:)
        integer :: i, j, l
        allocate(rmat(n(1),n(2),n(3)))
        do l = 1,n(3)
            do j = 1,n(2)
                do i = 1,n(1)
                    rmat(i,j,l) = real(cos(TWOPI_DP * (real(k(1)*(i-1),dp)/n(1) + real(k(2)*(j-1),dp)/n(2) + real(k(3)*(l-1),dp)/n(3))))
                end do
            end do
        end do
        call img%new(n, 1.)
        call img%set_rmat(rmat, .false.)
    end subroutine make_cosine

    ! one smoother for both boxes: the kernels are remade when the box changes
    subroutine test_transfer_2d()
        integer, parameter :: NBOX = 2, NFREQ = 5
        integer, parameter :: BOXES(3,NBOX)  = reshape([64,64,1, 64,48,1], [3,NBOX])
        integer, parameter :: FREQS(3,NFREQ) = reshape([0,0,0, 3,0,0, 5,7,0, 12,9,0, 0,11,0], [3,NFREQ])
        type(bspline_smoother) :: bs
        type(image)            :: img
        real, allocatable      :: before(:,:,:), after(:,:,:)
        real    :: h
        integer :: ib, il, ik
        write(*,'(A)') 'test_transfer_2d'
        call bs%new
        do ib = 1,NBOX
            do il = 1,size(LAMBDAS)
                do ik = 1,NFREQ
                    call make_cosine(img, BOXES(:,ib), FREQS(:,ik))
                    before = img%get_rmat()
                    call bs%smooth(img, LAMBDAS(il))
                    after  = img%get_rmat()
                    h      = real(transfer(FREQS(:,ik), BOXES(:,ib), LAMBDAS(il)))
                    call assert_real(0., maxval(abs(after - h * before)), 1.e-4, &
                        &int2str(BOXES(1,ib))//'x'//int2str(BOXES(2,ib))//', lambda '//trim(real2str(LAMBDAS(il)))//&
                        &', frequency ('//int2str(FREQS(1,ik))//','//int2str(FREQS(2,ik))//'): scaled by H = '//trim(real2str(h)))
                end do
            end do
        end do
        call bs%kill
        call img%kill
    end subroutine test_transfer_2d

    subroutine test_transfer_3d()
        integer, parameter :: NFREQ = 3, N(3) = [32,32,32]
        integer, parameter :: FREQS(3,NFREQ) = reshape([0,0,0, 3,2,1, 5,0,4], [3,NFREQ])
        type(bspline_smoother) :: bs
        type(image)            :: img
        real, allocatable      :: before(:,:,:), after(:,:,:)
        real    :: h
        integer :: il, ik
        write(*,'(A)') 'test_transfer_3d'
        call bs%new
        do il = 1,size(LAMBDAS)
            do ik = 1,NFREQ
                call make_cosine(img, N, FREQS(:,ik))
                before = img%get_rmat()
                call bs%smooth_3d(img, LAMBDAS(il))
                after  = img%get_rmat()
                h      = real(transfer(FREQS(:,ik), N, LAMBDAS(il)))
                call assert_real(0., maxval(abs(after - h * before)), 1.e-4, &
                    &'32^3, lambda '//trim(real2str(LAMBDAS(il)))//', frequency ('//int2str(FREQS(1,ik))//','//&
                    &int2str(FREQS(2,ik))//','//int2str(FREQS(3,ik))//'): scaled by H = '//trim(real2str(h)))
            end do
        end do
        call bs%kill
        call img%kill
    end subroutine test_transfer_3d

    ! Fourier input comes back in Fourier space, filtered the same as real-space input
    subroutine test_ft_state()
        type(bspline_smoother) :: bs
        type(image)            :: img_r, img_f
        write(*,'(A)') 'test_ft_state'
        call img_r%new([64,64,1], 1.)
        call img_r%gauran(0., 1.)
        call img_f%copy(img_r)
        call bs%new
        call bs%smooth(img_r, 0.2)
        call img_f%fft()
        call bs%smooth(img_f, 0.2)
        call assert_true(img_f%is_ft(), 'Fourier input stays in Fourier space')
        call img_f%ifft()
        call assert_real(0., maxval(abs(img_f%get_rmat() - img_r%get_rmat())), 1.e-5, 'real and Fourier input are filtered alike')
        call bs%kill
        call img_r%kill
        call img_f%kill
    end subroutine test_ft_state

end module simple_bspline_smoother_tester
