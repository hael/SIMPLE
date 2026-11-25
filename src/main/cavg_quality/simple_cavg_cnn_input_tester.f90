!@descr: unit tests of simple_cavg_cnn_input: the resampling against band-limited images sampled exactly at the new pixel size (odd and even sizes, down and up), the centre pixel, the scaling by the average's statistics, empty classes and preparation in a team
submodule (simple_cavg_cnn_input) simple_cavg_cnn_input_tester
use simple_test_utils
use simple_string_utils, only: int2str
!$ use omp_lib,          only: omp_set_num_threads, omp_get_max_threads, omp_get_num_threads, omp_get_thread_limit
implicit none
#include "simple_local_flags.inc"

! the operator in double precision on Gaussians whose tails at the box edge and spectra at Nyquist are below 1e-10
real(dp), parameter :: TOL_RESAMPLE = 1.d-9
integer,  parameter :: TBOX         = 48       !< the output box: the 35-37 px images are padded, the 49 px cropped

contains

    module subroutine run_all_cavg_cnn_input_tests()
        write(*,'(A)') '**** running all cavg_cnn_input tests ****'
        call test_resampling_exact(56, 1.3, 5.d0)   ! even to even, padded
        call test_resampling_exact(54, 1.3, 5.d0)   ! even to odd
        call test_resampling_exact(57, 1.3, 5.d0)   ! odd to odd
        call test_resampling_exact(38, 2.6, 7.d0)   ! up to an odd size, cropped
        call test_centre_pixel()
        call test_scaling()
        call test_empty()
        call test_team()
    end subroutine run_all_cavg_cnn_input_tests

    !> the band-limited test function at u, v (A): a Gaussian of sigma (A) times a slope along u, not symmetric
    real(dp) function test_fun( u, v, sigma )
        real(dp), intent(in) :: u, v, sigma
        test_fun = exp(-(u**2 + v**2) / (2.d0 * sigma**2)) * (1.d0 + 0.5d0 * sin(u / 7.d0))
    end function test_fun

    !> the resampled image is the function sampled at the new pixel size n smpd / m about the centre pixel, times
    !! (n / m)^2 (numpy's convention: fft2 unnormalised, ifft2 over m^2)
    subroutine test_resampling_exact( n, smpd, sigma )
        integer,  intent(in) :: n
        real,     intent(in) :: smpd
        real(dp), intent(in) :: sigma
        integer,  parameter   :: BOX = TBOX
        type(cavg_cnn_input)  :: inp
        real(dp), allocatable :: img(:,:), out(:,:)
        real(dp)              :: eff, err, expected
        integer               :: i, j, m
        write(*,'(A,I0,A,F4.2,A)') 'test_resampling_exact, ', n, ' px at ', smpd, ' A'
        m   = nint(real(n, dp) * real(smpd, dp) / 2.d0)
        eff = real(n, dp) * real(smpd, dp) / real(m, dp)
        allocate(img(n,n), out(BOX,BOX))
        do j = 1,n
            do i = 1,n
                img(i,j) = test_fun(real(i - 1 - n/2, dp) * smpd, real(j - 1 - n/2, dp) * smpd, sigma)
            enddo
        enddo
        call inp%new(n, smpd, BOX, 2.0)
        call inp%resample(img, out)
        err = 0.d0
        do j = 1,BOX
            do i = 1,BOX
                if( abs(i - 1 - BOX/2) > m/2 - 1 .or. abs(j - 1 - BOX/2) > m/2 - 1 ) cycle   ! the pad's zeros
                expected = test_fun(real(i - 1 - BOX/2, dp) * eff, real(j - 1 - BOX/2, dp) * eff, sigma) * (real(n, dp) / real(m, dp))**2
                err = max(err, abs(out(i,j) - expected))
            enddo
        enddo
        call assert_true(err <= TOL_RESAMPLE, 'resampling '//int2str(n)//' to '//int2str(m)//&
            &' px samples the band-limited image exactly at the new pixel size about the centre pixel')
        call inp%kill
    end subroutine test_resampling_exact

    !> a centred symmetric image keeps its maximum on the box's centre pixel (index box/2 from 0) at an odd size
    subroutine test_centre_pixel()
        integer, parameter   :: N = 54, BOX = TBOX
        type(cavg_cnn_input) :: inp
        real(dp)             :: img(N,N), out(BOX,BOX)
        integer              :: i, j, loc(2)
        write(*,'(A)') 'test_centre_pixel'
        do j = 1,N
            do i = 1,N
                img(i,j) = exp(-(real(i - 1 - N/2, dp)**2 + real(j - 1 - N/2, dp)**2) * 1.3d0**2 / (2.d0 * 5.d0**2))
            enddo
        enddo
        call inp%new(N, 1.3, BOX, 2.0)
        call inp%resample(img, out)
        loc = maxloc(out)
        call assert_true(all(loc == BOX/2 + 1), 'the centre of an odd-sized image stays on the centre pixel')
        call inp%kill
    end subroutine test_centre_pixel

    !> without resampling the average leaves at mean 0 and standard deviation 1 inside the inscribed circle, and the
    !! halves go through the same affine map
    subroutine test_scaling()
        integer, parameter   :: N = 24
        type(cavg_cnn_input) :: inp
        real(sp)             :: raw(N,N,3), x(N,N,3)
        real(dp)             :: mean, sdev, s1, s2, sd_ref
        integer              :: i, j, nin
        logical              :: empty
        write(*,'(A)') 'test_scaling'
        do j = 1,N
            do i = 1,N
                raw(i,j,1) = 5. + sin(0.3 * real(i)) + cos(0.7 * real(j))
                raw(i,j,2) = raw(i,j,1) + 0.1 * sin(real(i * j))
                raw(i,j,3) = raw(i,j,1) - 0.1 * sin(real(i * j))
            enddo
        enddo
        call inp%new(N, 2.0, N, 2.0)
        call inp%prepare(raw, x, empty)
        call assert_false(empty, 'a varying average is not empty')
        s1  = 0.d0
        s2  = 0.d0
        nin = 0
        do j = 1,N
            do i = 1,N
                if( real((i - 1 - N/2)**2 + (j - 1 - N/2)**2, dp) < (real(N, dp) / 2.d0)**2 )then
                    s1  = s1 + real(x(i,j,1), dp)
                    s2  = s2 + real(x(i,j,1), dp)**2
                    nin = nin + 1
                endif
            enddo
        enddo
        mean = s1 / real(nin, dp)
        sdev = sqrt(s2 / real(nin, dp) - mean**2)
        call assert_true(abs(mean) <= 1.d-6 .and. abs(sdev - 1.d0) <= 1.d-6, &
            &'the average has mean 0 and standard deviation 1 inside the inscribed circle')
        ! the halves minus the average scale by the average's standard deviation alone
        sd_ref = (real(raw(3,5,2), dp) - real(raw(3,5,1), dp)) / (real(x(3,5,2), dp) - real(x(3,5,1), dp))
        call assert_true(maxval(abs((real(raw(:,:,2), dp) - real(raw(:,:,1), dp)) - sd_ref * (real(x(:,:,2), dp) - real(x(:,:,1), dp)))) <= 1.d-5, &
            &'the halves are scaled by the average''s statistics, not their own')
        call inp%kill
    end subroutine test_scaling

    !> a constant average is empty and its input zero
    subroutine test_empty()
        integer, parameter   :: N = 16
        type(cavg_cnn_input) :: inp
        real(sp)             :: raw(N,N,3), x(N,N,3)
        logical              :: empty
        write(*,'(A)') 'test_empty'
        raw = 0.
        raw(:,:,2) = 1.
        x = 7.
        call inp%new(N, 2.0, N, 2.0)
        call inp%prepare(raw, x, empty)
        call assert_true(empty, 'a constant average is empty')
        call assert_true(all(x == 0.), 'an empty class''s input is zero')
        call inp%kill
    end subroutine test_empty

    !> preparing images concurrently from one object gives the serial result bit for bit
    subroutine test_team()
        integer, parameter   :: N = 54, BOX = TBOX, NIMG = 6, NTEAM = 3
        type(cavg_cnn_input) :: inp
        real(sp)             :: raw(N,N,1,NIMG), x1(BOX,BOX,1,NIMG), xn(BOX,BOX,1,NIMG)
        logical              :: empty(NIMG)
        integer              :: i, j, k, nthr_prev, limit, nteam_run
        write(*,'(A)') 'test_team'
        do k = 1,NIMG
            do j = 1,N
                do i = 1,N
                    raw(i,j,1,k) = real(test_fun(real(i - N/2 - k, dp), real(j - N/2 + k, dp), 5.d0))
                enddo
            enddo
        enddo
        call inp%new(N, 1.3, BOX, 2.0)
        do k = 1,NIMG
            call inp%prepare(raw(:,:,:,k), x1(:,:,:,k), empty(k))
        enddo
        nthr_prev = 1
        !$ nthr_prev = omp_get_max_threads()
        limit = NTEAM
        !$ limit = omp_get_thread_limit()
        if( limit >= NTEAM )then
            nteam_run = 1
            !$omp parallel do num_threads(NTEAM) default(shared) private(k) schedule(static,1)
            do k = 1,NIMG
                call inp%prepare(raw(:,:,:,k), xn(:,:,:,k), empty(k))
                !$ if( k == 1 ) nteam_run = omp_get_num_threads()
            enddo
            !$omp end parallel do
            call assert_int(NTEAM, nteam_run, 'the team of three formed')
            call assert_true(all(transfer(x1, 0, size(x1)) == transfer(xn, 0, size(xn))), 'images prepared in a team are the serial ones bit for bit')
        else
            write(*,'(A)') 'test_team: the team of three skipped under the thread limit'
        endif
        !$ call omp_set_num_threads(nthr_prev)
        call inp%kill
    end subroutine test_team

end submodule simple_cavg_cnn_input_tester
