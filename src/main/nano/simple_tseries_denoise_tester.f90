!@descr: rectangular Fourier invariants, noise reduction and windowed thread identity of tracking denoising
module simple_tseries_denoise_tester
use simple_defs,           only: PI, logfhandle
use simple_image,          only: image
use simple_tseries_denoise, only: flex_denoise_imgarr
use simple_imgarr_utils,   only: dealloc_imgarr
use simple_rnd,            only: gasdev
use simple_test_utils,     only: assert_true, assert_false, assert_int, assert_real, set_fixed_seed
!$ use omp_lib, only: omp_get_max_threads, omp_set_num_threads, omp_get_dynamic, omp_set_dynamic
implicit none
private
public :: run_all_tseries_denoise_tests

integer, parameter :: NX = 24, NY = 16, NFRAMES = 24
! Two single-precision transforms and the kernel reductions, with amplitude bounded by five.
real, parameter :: FFT_TOL = 64.*epsilon(1.)*5., SMPD = 1.
real, parameter :: NOISE_SD = 0.2, ERROR_RATIO_MAX = 0.95

contains

    subroutine run_all_tseries_denoise_tests()
        integer :: saved_log, scratch
        write(*,'(A)') '**** running all tseries denoise tests ****'
        saved_log = logfhandle
        open(newunit=scratch, status='scratch', action='readwrite')
        logfhandle = scratch
        call test_noiseless_rectangle()
        call test_static_noise_and_windows()
        logfhandle = saved_log
        close(scratch)
    end subroutine run_all_tseries_denoise_tests

    !> Preserve DC, low spatial frequencies and the Nyquist checkerboard on a rectangular grid.
    subroutine test_noiseless_rectangle()
        type(image), allocatable :: imgs(:), den(:)
        real, allocatable :: actual(:,:,:)
        real :: truth(NX,NY,1)
        integer :: n, i, x, y
        write(*,'(A)') 'test_noiseless_rectangle'
        do y = 1, NY
            do x = 1, NX
                truth(x,y,1) = 3. + cos(2.*PI*real(x-1)/real(NX)) + &
                    &0.3*sin(2.*PI*real(y-1)/real(NY)) + 0.2*real((-1)**(x+y))
            enddo
        enddo
        do n = 2, 12, 10
            allocate(imgs(n))
            do i = 1, n
                call imgs(i)%new([NX,NY,1], SMPD, wthreads=.false.)
                call imgs(i)%set_rmat(truth, .false.)
            enddo
            call flex_denoise_imgarr(imgs, den, 4, 0, .false.)
            call assert_int(n, size(den), 'tracking denoise preserves frame count')
            do i = 1, n
                actual = den(i)%get_rmat()
                call assert_true(all(den(i)%get_ldim() == [NX,NY,1]), 'tracking denoise preserves rectangular geometry')
                call assert_real(SMPD, den(i)%get_smpd(), 0., 'tracking denoise preserves sampling')
                call assert_true(maxval(abs(actual-truth)) < FFT_TOL, 'noiseless DC and spatial frequencies are preserved')
                actual = imgs(i)%get_rmat()
                call assert_true(all(actual == truth), 'tracking denoise leaves source pixels unchanged')
                call assert_false(imgs(i)%is_ft(), 'tracking denoise leaves the source in real space')
            enddo
            if( n == 2 )then
                call den(1)%zero
                actual = imgs(1)%get_rmat()
                call assert_true(all(actual == truth), 'short-batch copied frames own independent buffers')
            endif
            call dealloc_imgarr(den)
            call dealloc_imgarr(imgs)
        enddo
    end subroutine test_noiseless_rectangle

    !> Static signal plus independent white noise has a known clean target; narrow windows have too few frames to fit.
    subroutine test_static_noise_and_windows()
        type(image), allocatable :: imgs(:), den(:), serial(:)
        real, allocatable :: actual(:,:,:), expected(:,:,:)
        real :: truth(NX,NY,1), noisy(NX,NY,1), raw_error, den_error
        integer :: t, x, y, mode, halfwin
        !$ integer :: threads_saved
        !$ logical :: dynamic_saved
        write(*,'(A)') 'test_static_noise_and_windows'
        do y = 1, NY
            do x = 1, NX
                truth(x,y,1) = 2. + exp(-real((x-12)**2+(y-8)**2)/8.)
            enddo
        enddo
        allocate(imgs(NFRAMES))
        call set_fixed_seed(20261106)
        raw_error = 0.
        do t = 1, NFRAMES
            do y = 1, NY
                do x = 1, NX
                    noisy(x,y,1) = truth(x,y,1) + NOISE_SD*gasdev()
                enddo
            enddo
            raw_error = raw_error + sum((noisy-truth)**2)
            call imgs(t)%new([NX,NY,1], SMPD, wthreads=.false.)
            call imgs(t)%set_rmat(noisy, .false.)
        enddo
        call assert_true(raw_error > 0., 'independent noise supplies a nonzero raw-error control')
        do mode = 1, 2
            call flex_denoise_imgarr(imgs, den, 4, merge(0,8,mode == 1), mode == 2)
            den_error = 0.
            do t = 1, NFRAMES
                actual = den(t)%get_rmat()
                den_error = den_error + sum((actual-truth)**2)
            enddo
            call assert_true(den_error < ERROR_RATIO_MAX*raw_error, 'kernel pooling reduces error against the clean planted image')
            call dealloc_imgarr(den)
        enddo
        !$ threads_saved = omp_get_max_threads()
        !$ dynamic_saved = omp_get_dynamic()
        !$ call omp_set_dynamic(.false.)
        do halfwin = 1, 4, 3
            !$ call omp_set_num_threads(1)
            call flex_denoise_imgarr(imgs, serial, 4, halfwin, .false.)
            !$ call omp_set_num_threads(3)
            call flex_denoise_imgarr(imgs, den, 4, halfwin, .false.)
            do t = 1, NFRAMES
                expected = imgs(t)%get_rmat()
                actual = serial(t)%get_rmat()
                call assert_true(maxval(abs(actual-expected)) < FFT_TOL, 'undersized windows retain every target frame')
                expected = den(t)%get_rmat()
                call assert_true(all(actual == expected), 'narrow-window outputs are identical on one and three threads')
            enddo
            call dealloc_imgarr(serial)
            call dealloc_imgarr(den)
        enddo
        !$ call omp_set_num_threads(threads_saved)
        !$ call omp_set_dynamic(dynamic_saved)
        call dealloc_imgarr(imgs)
    end subroutine test_static_noise_and_windows

end module simple_tseries_denoise_tester
