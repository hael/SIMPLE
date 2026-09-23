!@descr: unit test routines for segmentation: thresholding and peak detection
!! The connected-component contract of image_bin lives in simple_image_msk_tester
!! (sub-suite 'binary image'); this module is for the routines of simple_segmentation
!! and the thresholding utilities of simple_math that feed them.
module simple_segmentation_tester
use simple_test_utils        ! assertions etc.
use simple_defs
use simple_math,             only: otsu
use simple_rnd,              only: gasdev
use simple_segmentation,     only: detect_peak_thres_fdr
implicit none
private
public :: run_all_segmentation_tests

contains

    subroutine run_all_segmentation_tests()
        write(*,'(A)') '**** running all segmentation tests ****'
        call test_otsu()
        call test_peak_thres_fdr()
    end subroutine run_all_segmentation_tests

    !---------------- Otsu thresholding ----------------

    ! two well-separated Gaussian classes: the threshold falls between them, the
    ! classes come back the right size, and the three overloads agree
    subroutine test_otsu()
        integer, parameter :: NPER = 500, N = 2*NPER
        real,    parameter :: MEAN_LO = 0.0, MEAN_HI = 10.0, SDEV = 1.0
        real    :: x(N), x_copy(N), x_out(N), thresh, thresh2
        logical :: mask(N)
        integer :: i, nhi
        write(*,'(A)') 'test_otsu'
        do i = 1,NPER
            x(i)      = gasdev(MEAN_LO, SDEV)
            x(NPER+i) = gasdev(MEAN_HI, SDEV)
        end do
        call otsu(N, x, thresh)
        call assert_true(thresh > MEAN_LO + 2.0*SDEV .and. thresh < MEAN_HI - 2.0*SDEV,&
            &'otsu: threshold lies between the two classes')
        x_copy = x
        call otsu(N, x_copy, mask)
        nhi = count(mask)
        call assert_true(abs(nhi - NPER) <= 5,               'otsu: the upper class has the expected size (within 1%)')
        call assert_true(all(mask(NPER+1:) .or. x(NPER+1:) <= thresh),&
            &'otsu (mask): upper-class samples above the threshold are selected')
        call assert_true(.not. any(mask(1:NPER) .and. x(1:NPER) <= thresh),&
            &'otsu (mask): lower-class samples below the threshold are not')
        x_copy = x
        call otsu(N, x_copy, x_out, thresh2)
        call assert_real(thresh, thresh2, 1.0e-6,                    'otsu: the binarising overload reports the same threshold')
        call assert_true(all((x_out > 0.5) .eqv. mask),               'otsu: the binarising overload agrees with the mask overload')
    end subroutine test_otsu

    !---------------- FDR peak threshold ----------------

    ! synthetic scores with clear outliers: the upper and lower tails are found, the
    ! threshold and the count agree, and the min/max peak bounds are honoured
    subroutine test_peak_thres_fdr()
        integer, parameter :: N = 20
        real    :: x_up(N), x_low(N), t
        integer :: npeaks
        write(*,'(A)') 'test_peak_thres_fdr'
        ! upper tail: two clear high outliers
        x_up = [ -0.35, -0.25, -0.20, -0.12, -0.08, -0.03, 0.00, 0.01, 0.02, 0.03, &
                  0.04,  0.05,  0.06,  0.07,  0.09,  0.11, 0.14, 0.18, 5.00, 6.00 ]
        call detect_peak_thres_fdr(N, x_up, 0.25, 0, N, t, npeaks)
        call assert_true(npeaks >= 1,                'FDR upper tail detects at least one peak')
        call assert_int(count(x_up >= t), npeaks,    'FDR upper tail: threshold and count agree')
        call detect_peak_thres_fdr(N, x_up, 0.25, 0, 1, t, npeaks)
        call assert_int(1, npeaks,                   'FDR upper tail: max_peaks cap is respected')
        ! lower tail: three clear low outliers
        x_low = [ -6.00, -5.00, -4.00, -0.20, -0.14, -0.09, -0.05, -0.02, 0.00, 0.01, &
                   0.02,  0.03,  0.04,  0.05,  0.07,  0.10,  0.13,  0.16, 0.22, 0.30 ]
        call detect_peak_thres_fdr(N, x_low, 0.25, 0, N, t, npeaks, lower_tail=.true.)
        call assert_true(npeaks >= 1,                'FDR lower tail detects at least one peak')
        call assert_int(count(x_low <= t), npeaks,   'FDR lower tail: threshold and count agree')
        call detect_peak_thres_fdr(N, x_low, 0.25, 3, 3, t, npeaks, lower_tail=.true.)
        call assert_int(3, npeaks,                   'FDR lower tail: min/max peak bounds are respected')
    end subroutine test_peak_thres_fdr

end module simple_segmentation_tester
