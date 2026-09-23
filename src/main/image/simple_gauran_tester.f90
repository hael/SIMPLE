!@descr: unit test routines for the Gaussian noise generators of simple_image (gauran, add_gauran)
! The observation-noise contracts reconstruction relies on: gauran draws N(mean, sdev^2)
! per voxel; add_gauran adds white noise with variance var(signal)/snr so the realised
! SNR is the requested one; draws are reproducible after reseeding and independent
! between calls. Tolerances are statistical (24^3 = 13824 samples), never exact pins,
! because the stream behind random_number is compiler-specific.
module simple_gauran_tester
use simple_defs,       only: dp
use simple_image,      only: image
use simple_test_utils
use ieee_arithmetic,   only: ieee_is_finite
implicit none
private
public :: run_all_gauran_tests

integer,  parameter :: BOX  = 24
real,     parameter :: SMPD = 1.5
integer,  parameter :: SEED = 20260812
! standard error of the mean of N unit-variance samples is 1/sqrt(N) = 0.0085, of the
! variance sqrt(2/N) = 0.012: the tolerances below are 5-6 sigma
real(dp), parameter :: MEAN_TOL = 0.05_dp
real(dp), parameter :: VAR_TOL  = 0.06_dp
real(dp), parameter :: CORR_TOL = 0.05_dp

contains

    subroutine run_all_gauran_tests()
        write(*,'(A)') '**** running all gauran tests ****'
        call test_gauran_moments()
        call test_add_gauran_snr()
        call test_replay_and_independence()
    end subroutine run_all_gauran_tests

    !> a smooth signal with non-trivial variance: one off-centre Gaussian blob
    subroutine make_signal( img )
        type(image), intent(inout) :: img
        real, allocatable :: rmat(:,:,:)
        real    :: ctr, dx, dy, dz
        integer :: i, j, k
        allocate(rmat(BOX,BOX,BOX))
        ctr = real(BOX)/2. + 0.5
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    dx = real(i) - ctr - 3.; dy = real(j) - ctr + 2.; dz = real(k) - ctr - 1.
                    rmat(i,j,k) = exp(-(dx*dx + dy*dy + dz*dz) / (2. * 2.5**2))
                enddo
            enddo
        enddo
        call img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call img%set_rmat(rmat, .false.)
    end subroutine make_signal

    pure real(dp) function pmean( x )
        real, intent(in) :: x(:,:,:)
        pmean = sum(real(x,dp)) / real(size(x),dp)
    end function pmean

    pure real(dp) function pvar( x )
        real, intent(in) :: x(:,:,:)
        real(dp) :: m
        m = pmean(x)
        pvar = sum((real(x,dp) - m)**2) / real(size(x),dp)
    end function pvar

    pure real(dp) function pcorr( a, b )
        real, intent(in) :: a(:,:,:), b(:,:,:)
        real(dp) :: am, bm, den
        am  = pmean(a)
        bm  = pmean(b)
        den = sqrt(sum((real(a,dp)-am)**2) * sum((real(b,dp)-bm)**2))
        pcorr = sum((real(a,dp)-am) * (real(b,dp)-bm)) / max(den, tiny(den))
    end function pcorr

    !> gauran fills the image with N(mean, sdev^2) samples
    subroutine test_gauran_moments()
        type(image) :: img
        real, allocatable :: x(:,:,:)
        write(*,'(A)') 'test_gauran_moments'
        call img%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call set_fixed_seed(SEED)
        call img%gauran(0., 1.)
        x = img%get_rmat()
        call assert_true(all(ieee_is_finite(x)), 'unit Gaussian samples are finite')
        call assert_true(abs(pmean(x)) < MEAN_TOL, 'unit Gaussian mean is 0 (|mean| < 0.05)')
        call assert_true(abs(pvar(x) - 1._dp) < VAR_TOL, 'unit Gaussian variance is 1 (within 0.06)')
        call assert_false(img%is_ft(), 'gauran leaves a real-space image')
        call img%gauran(3., 0.5)
        x = img%get_rmat()
        call assert_true(abs(pmean(x) - 3._dp) < 0.5_dp*MEAN_TOL, 'gauran(3, 0.5) mean is 3')
        call assert_true(abs(pvar(x) - 0.25_dp) < 0.25_dp*VAR_TOL, 'gauran(3, 0.5) variance is 0.25')
        call img%kill()
    end subroutine test_gauran_moments

    !> add_gauran adds zero-mean white noise of variance var(signal)/snr
    subroutine test_add_gauran_snr()
        type(image) :: clean_img, noisy_img
        real, allocatable :: clean(:,:,:), noisy(:,:,:), noise(:,:,:)
        real(dp) :: snr_realised
        real     :: snr
        integer  :: i
        write(*,'(A)') 'test_add_gauran_snr'
        call make_signal(clean_img)
        clean = clean_img%get_rmat()
        call set_fixed_seed(SEED + 1)
        do i = 1, 3
            snr = snr_of(i)
            call noisy_img%copy(clean_img)
            call noisy_img%add_gauran(snr)
            noisy = noisy_img%get_rmat()
            noise = noisy - clean
            snr_realised = pvar(clean) / pvar(noise)
            call assert_true(all(ieee_is_finite(noisy)), 'noisy volume is finite')
            call assert_true(any(noise /= 0.), 'add_gauran changed the volume')
            call assert_true(abs(snr_realised - real(snr,dp)) / real(snr,dp) < VAR_TOL, &
                &'realised SNR var(signal)/var(noise) matches the request within 6%')
            call assert_true(abs(pmean(noise)) / sqrt(pvar(noise)) < MEAN_TOL, 'added noise has zero mean')
            call assert_true(abs(pcorr(noise, clean)) < CORR_TOL, 'added noise is uncorrelated with the signal')
        enddo
        call assert_true(all(clean_img%get_rmat() == clean), 'the source image is untouched')
        call clean_img%kill()
        call noisy_img%kill()

    contains

        real function snr_of( idx )
            integer, intent(in) :: idx
            real, parameter :: SNRS(3) = [0.1, 0.5, 2.0]
            snr_of = SNRS(idx)
        end function snr_of

    end subroutine test_add_gauran_snr

    !> the same seed reproduces a draw; consecutive draws are independent
    subroutine test_replay_and_independence()
        type(image) :: clean_img, a_img, b_img, replay_img
        real, allocatable :: clean(:,:,:), a(:,:,:), b(:,:,:), replay(:,:,:)
        write(*,'(A)') 'test_replay_and_independence'
        call make_signal(clean_img)
        clean = clean_img%get_rmat()
        call a_img%copy(clean_img)
        call b_img%copy(clean_img)
        call replay_img%copy(clean_img)
        call set_fixed_seed(SEED + 2)
        call a_img%add_gauran(0.5)
        call b_img%add_gauran(0.5)
        call set_fixed_seed(SEED + 2)
        call replay_img%add_gauran(0.5)
        a      = a_img%get_rmat() - clean
        b      = b_img%get_rmat() - clean
        replay = replay_img%get_rmat() - clean
        call assert_true(all(a == replay), 'reseeding reproduces the draw exactly')
        call assert_true(any(a /= b), 'consecutive draws differ')
        call assert_true(abs(pcorr(a, b)) < CORR_TOL, 'consecutive draws are uncorrelated (|r| < 0.05)')
        call clean_img%kill()
        call a_img%kill()
        call b_img%kill()
        call replay_img%kill()
    end subroutine test_replay_and_independence

end module simple_gauran_tester
