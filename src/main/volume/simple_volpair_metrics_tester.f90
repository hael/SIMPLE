!@descr: unit tests for comparing two maps inside a soft mask (simple_volpair_metrics): FSC and band correlation
! Synthetic maps at box 32: a map of three Gaussian blobs over weak white noise
! against itself (FSC and correlation of one), against independent white noise
! (neither correlates), against itself under white noise of a comparable power
! (the band up to 10 A correlates better than the full band, and the FSC falls
! with frequency),
! and against a map of another box (refused).
module simple_volpair_metrics_tester
use simple_core_module_api
use simple_image,           only: image
use simple_volpair_metrics, only: compare_volpair
use simple_test_utils
implicit none

public :: run_all_volpair_metrics_tests
private
#include "simple_local_flags.inc"

integer,          parameter :: BOX     = 32
real,             parameter :: SMPD    = 2.
real,             parameter :: MSKDIAM = 50.
integer,          parameter :: SEED    = 20260927
character(len=*), parameter :: F_MAP   = 'tmp_volpair_metrics_map.mrc'
character(len=*), parameter :: F_NOISE = 'tmp_volpair_metrics_noise.mrc'
character(len=*), parameter :: F_NOISY = 'tmp_volpair_metrics_noisy.mrc'
character(len=*), parameter :: F_SMALL = 'tmp_volpair_metrics_small.mrc'

contains

    subroutine run_all_volpair_metrics_tests()
        write(logfhandle,'(A)') '**** running all volume pair metrics tests ****'
        call write_fixtures()
        call test_identical_maps()
        call test_independent_maps()
        call test_band_correlation()
        call test_mismatched_maps()
        call del_file(F_MAP)
        call del_file(F_NOISE)
        call del_file(F_NOISY)
        call del_file(F_SMALL)
    end subroutine run_all_volpair_metrics_tests

    ! ---- fixtures -------------------------------------------------------------

    !> three blobs over weak noise; independent noise; the blobs under strong
    !! noise; and a smaller map
    subroutine write_fixtures()
        real, allocatable :: blobs(:,:,:), noise(:,:,:)
        call set_fixed_seed(SEED)
        blobs = blob_map()
        noise = white_noise(BOX)
        call write_map(F_MAP,   blobs + 0.02 * noise, BOX)
        noise = white_noise(BOX)
        call write_map(F_NOISE, noise, BOX)
        noise = white_noise(BOX)
        call write_map(F_NOISY, blobs + noise, BOX)
        call write_map(F_SMALL, white_noise(BOX - 8), BOX - 8)
    end subroutine write_fixtures

    function blob_map() result( rmat )
        real, allocatable :: rmat(:,:,:)
        real, parameter :: CENTRES(3,3) = reshape([-4.,2.,1., 5.,-3.,0., 1.,5.,-4.], [3,3])
        real    :: d2
        integer :: i, j, k, b
        allocate(rmat(BOX,BOX,BOX), source=0.)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do b = 1, 3
                        d2 = (real(i - BOX/2 - 1) - CENTRES(1,b))**2 + (real(j - BOX/2 - 1) - CENTRES(2,b))**2 + &
                            &(real(k - BOX/2 - 1) - CENTRES(3,b))**2
                        rmat(i,j,k) = rmat(i,j,k) + exp(-0.5 * d2 / 4.)
                    enddo
                enddo
            enddo
        enddo
    end function blob_map

    function white_noise( n ) result( rmat )
        integer, intent(in) :: n
        real, allocatable :: rmat(:,:,:)
        allocate(rmat(n,n,n))
        call random_number(rmat)
        rmat = rmat - 0.5
    end function white_noise

    subroutine write_map( fname, rmat, n )
        character(len=*), intent(in) :: fname
        real,             intent(in) :: rmat(:,:,:)
        integer,          intent(in) :: n
        type(image) :: vol
        call vol%new([n,n,n], SMPD, wthreads=.false.)
        call vol%set_rmat(rmat, .false.)
        call vol%write(string(fname), del_if_exists=.true.)
        call vol%kill
    end subroutine write_map

    ! ---- tests ----------------------------------------------------------------

    subroutine test_identical_maps()
        real, allocatable :: fsc(:), res(:)
        real    :: corr
        logical :: ok
        write(logfhandle,'(A)') '-- test_identical_maps'
        call compare_volpair(string(F_MAP), string(F_MAP), MSKDIAM, 0., corr, fsc, res, ok)
        call assert_true(ok, 'a map against itself is compared')
        call assert_real(1., corr, 1.e-4, 'a map correlates fully with itself')
        call assert_true(minval(fsc(2:)) > 0.999, 'the FSC of a map with itself is one in every shell')
        call assert_int(size(fsc), size(res), 'one resolution per FSC shell')
        call assert_real(2. * SMPD, res(size(res)), 0.2, 'the last shell lies at Nyquist')
    end subroutine test_identical_maps

    subroutine test_independent_maps()
        real, allocatable :: fsc(:), res(:)
        real    :: corr
        logical :: ok
        write(logfhandle,'(A)') '-- test_independent_maps'
        call compare_volpair(string(F_MAP), string(F_NOISE), MSKDIAM, 0., corr, fsc, res, ok)
        call assert_true(ok, 'a map against independent noise is compared')
        call assert_true(abs(corr) < 0.2, 'a map does not correlate with independent noise')
        call assert_true(sum(abs(fsc(3:))) / real(size(fsc) - 2) < 0.2, 'nor does its FSC')
    end subroutine test_independent_maps

    subroutine test_band_correlation()
        real, allocatable :: fsc(:), res(:)
        real    :: corr_full, corr_band
        logical :: ok
        integer :: n
        write(logfhandle,'(A)') '-- test_band_correlation'
        call compare_volpair(string(F_MAP), string(F_NOISY), MSKDIAM, 0.,  corr_full, fsc, res, ok)
        call assert_true(ok, 'the full band is compared')
        call compare_volpair(string(F_MAP), string(F_NOISY), MSKDIAM, 10., corr_band, fsc, res, ok)
        call assert_true(ok, 'the band to 10 A is compared')
        write(logfhandle,'(A,F7.4,A,F7.4)') '   correlation to 10 A ', corr_band, ', full band ', corr_full
        call assert_true(corr_band > corr_full + 0.1, 'white noise hurts the full band more than the band to 10 A')
        n = size(fsc)
        call assert_true(sum(fsc(2:4)) > sum(fsc(n-2:n)), 'the FSC falls where the noise dominates')
    end subroutine test_band_correlation

    subroutine test_mismatched_maps()
        real, allocatable :: fsc(:), res(:)
        real    :: corr
        logical :: ok
        write(logfhandle,'(A)') '-- test_mismatched_maps'
        call compare_volpair(string(F_MAP), string(F_SMALL), MSKDIAM, 0., corr, fsc, res, ok)
        call assert_false(ok, 'maps of different boxes are refused')
        call assert_real(0., corr, 0., 'a refused pair has no correlation')
        call compare_volpair(string(F_MAP), string('tmp_volpair_metrics_absent.mrc'), MSKDIAM, 0., corr, fsc, res, ok)
        call assert_false(ok, 'a missing map is refused')
    end subroutine test_mismatched_maps

end module simple_volpair_metrics_tester
