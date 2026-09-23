!@descr: unit test routines for segmentation: thresholding and peak detection
!! The connected-component contract of image_bin lives in simple_image_msk_tester
!! (sub-suite 'binary image'); this module is for the routines of simple_segmentation
!! and the thresholding utilities of simple_math that feed them.
module simple_segmentation_tester
use simple_test_utils        ! assertions etc.
use simple_defs
use simple_math,             only: otsu
use simple_rnd,              only: gasdev
use simple_image,            only: image
use simple_segmentation,     only: detect_peak_thres_fdr, detect_peak_thres, detect_peak_thres_for_npeaks,&
                                  &refine_peak_thres_sortmeans, otsu_img, otsu_robust_fast, sauvola, sobel, canny
implicit none
private
public :: run_all_segmentation_tests

contains

    subroutine run_all_segmentation_tests()
        write(*,'(A)') '**** running all segmentation tests ****'
        call test_otsu()
        call test_peak_thres_fdr()
        call test_detect_peak_thres()
        call test_binarize()
        call test_otsu_img()
        call test_otsu_robust_fast()
        call test_sauvola()
        call test_sobel()
        call test_canny()
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

    !==================================================================
    ! thresholding of images
    !==================================================================

    !---------------- binarize ----------------

    ! a ramp image has one pixel per value, so every count is known exactly
    subroutine test_binarize()
        integer, parameter :: BOX = 32
        type(image) :: img, out
        logical     :: mask(BOX,BOX,1)
        integer     :: i, j
        write(*,'(A)') 'test_binarize'
        call img%new([BOX,BOX,1], 1.0)
        do j = 1,BOX
            do i = 1,BOX
                call img%set([i,j,1], real((j-1)*BOX + i))   ! 1..BOX*BOX, distinct
            end do
        end do
        ! threshold form, in place: >= thres is foreground
        call out%copy(img)
        call out%binarize(512.5)
        call assert_int(512, out%nforeground(), 'binarize(thres): pixels at or above the threshold are foreground')
        call assert_int(512, out%nbackground(), 'binarize(thres): the rest are background')
        call assert_real(1.0, out%get([BOX,BOX,1]), 1.0e-6, 'binarize(thres): largest value is 1')
        call assert_real(0.0, out%get([1,1,1]),     1.0e-6, 'binarize(thres): smallest value is 0')
        ! threshold form with a separate output leaves the input untouched
        call out%new([BOX,BOX,1], 1.0)
        call img%binarize(1000.0, out)
        call assert_int(BOX*BOX - 999, out%nforeground(),        'binarize(thres, out): foreground count')
        call assert_real(1.0, img%get([1,1,1]), 1.0e-6,           'binarize(thres, out): input is untouched')
        ! mask form
        call img%binarize(1000.0, mask)
        call assert_int(BOX*BOX - 999, count(mask),               'binarize(thres, mask): mask count')
        call assert_true(mask(BOX,BOX,1) .and. .not. mask(1,1,1), 'binarize(thres, mask): mask orientation')
        ! count form: exactly npix foreground pixels, the largest ones
        call out%copy(img)
        call out%binarize(100)
        call assert_int(100, out%nforeground(),                   'binarize(npix): exactly npix foreground pixels')
        call assert_real(1.0, out%get([BOX,BOX,1]), 1.0e-6,       'binarize(npix): the largest pixel is foreground')
        call img%kill
        call out%kill
    end subroutine test_binarize

    !---------------- otsu_img ----------------

    ! two-level images with a little noise: the threshold falls between the levels
    ! and the binarised image is exactly the object
    subroutine test_otsu_img()
        integer, parameter :: BOX = 64
        real,    parameter :: RAD = 12.0, RAD_CORE = 6.0, NOISE = 0.01
        type(image) :: img
        real    :: thresh
        integer :: ndisc, ncore
        write(*,'(A)') 'test_otsu_img'
        ! background 0, disc 1
        call two_level(img, BOX, RAD, 0.0, 1.0, NOISE, ndisc)
        call otsu_img(img, thresh)
        call assert_true(thresh > 0.0 .and. thresh < 1.0,         'otsu_img: threshold lies between the two levels')
        call assert_int(ndisc, img%nforeground(),                 'otsu_img: the binarised image is the disc')
        call assert_true(all_binary(img),                         'otsu_img: output is 0/1')
        ! zero-mean background: positive=.true. builds the histogram from the positive pixels only
        call two_level(img, BOX, RAD, 0.0, 1.0, 2.0*NOISE, ndisc)
        call otsu_img(img, thresh, positive=.true.)
        call assert_true(thresh > 0.0 .and. thresh < 1.0,         'otsu_img(positive): threshold between 0 and the object level')
        call assert_int(ndisc, img%nforeground(),                 'otsu_img(positive): the binarised image is the disc')
        ! tight: background 0, halo 0.5, core 1: a second Otsu on the foreground keeps the core
        call three_level(img, BOX, RAD, RAD_CORE, [0.0, 0.5, 1.0], NOISE, ncore)
        call otsu_img(img, thresh, tight=.true.)
        call assert_true(thresh > 0.5 .and. thresh < 1.0,         'otsu_img(tight): threshold between halo and core')
        call assert_int(ncore, img%nforeground(),                 'otsu_img(tight): the binarised image is the core')
        ! tighter: four levels (0, 0.5, 0.8, 1 at radii 26, 18, 12), a third Otsu keeps the core
        call four_level(img, BOX, 26.0, 18.0, 12.0, [0.0, 0.5, 0.8, 1.0], NOISE, ncore)
        call otsu_img(img, thresh, tighter=.true.)
        call assert_true(thresh > 0.8 .and. thresh < 1.0,&
            &'otsu_img(tighter): threshold between the third level and the core')
        call assert_int(ncore, img%nforeground(),                 'otsu_img(tighter): the binarised image is the core')
        call img%kill
    end subroutine test_otsu_img

    !---------------- otsu_robust_fast ----------------

    ! the majority vote of the original, the local mean and the local median removes
    ! salt-and-pepper noise that a single Otsu would keep
    subroutine test_otsu_robust_fast()
        integer, parameter :: BOX = 64
        real,    parameter :: RAD = 14.0
        type(image) :: img, clean
        real    :: thresh(3), d
        integer :: ndisc, i, j, c, nflip, nwrong_far, nwrong_in
        logical :: is_in
        write(*,'(A)') 'test_otsu_robust_fast'
        call two_level(clean, BOX, RAD, 0.0, 1.0, 0.01, ndisc)
        call img%copy(clean)
        ! salt and pepper: every 13th pixel flipped, well inside or well outside the disc
        c = BOX/2 + 1
        nflip = 0
        do j = 1,BOX
            do i = 1,BOX
                if( mod((j-1)*BOX + i, 13) /= 0 ) cycle
                d = sqrt(real((i-c)**2 + (j-c)**2))
                if( d < RAD - 3.0 )then
                    call img%set([i,j,1], 0.0); nflip = nflip + 1
                else if( d > RAD + 3.0 )then
                    call img%set([i,j,1], 1.0); nflip = nflip + 1
                endif
            end do
        end do
        call assert_true(nflip > 100, 'fixture: enough salt-and-pepper pixels')
        call otsu_robust_fast(img, is2D=.true., noneg=.false., thresh=thresh)
        call assert_true(all_binary(img), 'otsu_robust_fast: output is 0/1')
        ! every flipped pixel away from the edge is repaired; the edge may move by a pixel
        nwrong_far = 0
        nwrong_in  = 0
        do j = 1,BOX
            do i = 1,BOX
                d = sqrt(real((i-c)**2 + (j-c)**2))
                is_in = img%get([i,j,1]) > 0.5
                if( d > RAD + 2.0 .and. is_in )       nwrong_far = nwrong_far + 1
                if( d < RAD - 2.0 .and. .not. is_in ) nwrong_in  = nwrong_in  + 1
            end do
        end do
        call assert_int(0, nwrong_far, 'otsu_robust_fast: no foreground left in the background (salt removed)')
        call assert_int(0, nwrong_in,  'otsu_robust_fast: no holes left in the object (pepper removed)')
        call assert_true(abs(img%nforeground() - ndisc) < 0.05 * ndisc, 'otsu_robust_fast: object area within 5% of the disc')
        call img%kill
        call clean%kill
    end subroutine test_otsu_robust_fast

    !---------------- sauvola ----------------

    ! the local standard deviations are those of the clipped (2*winsz+1)^2 window, and the
    ! binarisation follows t = m * (1 + bias * (s/s_max - 1)) pixel by pixel
    subroutine test_sauvola()
        integer, parameter :: BOX = 32, WINSZ = 2
        real,    parameter :: BIAS = 0.34
        type(image) :: img, sdevs
        real    :: ref(BOX,BOX), avgs(BOX,BOX), sd(BOX,BOX), t, smax, maxdiff
        integer :: i, j, i1, i2, j1, j2, np, nmismatch
        write(*,'(A)') 'test_sauvola'
        ! a step in x plus a slow ramp in y, so windows see both flat and edge regions
        do j = 1,BOX
            do i = 1,BOX
                ref(i,j) = merge(1.0, 0.0, i > BOX/2) + 0.01*real(j)
            end do
        end do
        call img%new([BOX,BOX,1], 1.0)
        do j = 1,BOX
            do i = 1,BOX
                call img%set([i,j,1], ref(i,j))
            end do
        end do
        ! reference: brute-force window statistics
        do j = 1,BOX
            do i = 1,BOX
                i1 = max(1,i-WINSZ); i2 = min(BOX,i+WINSZ)
                j1 = max(1,j-WINSZ); j2 = min(BOX,j+WINSZ)
                np = (i2-i1+1)*(j2-j1+1)
                avgs(i,j) = sum(ref(i1:i2,j1:j2)) / real(np)
                sd(i,j)   = sqrt(sum((ref(i1:i2,j1:j2) - avgs(i,j))**2) / real(np-1))
            end do
        end do
        smax = maxval(sd)
        call sauvola(img, WINSZ, sdevs, BIAS)
        maxdiff = 0.0
        do j = 1,BOX
            do i = 1,BOX
                maxdiff = max(maxdiff, abs(sdevs%get([i,j,1]) - sd(i,j)))
            end do
        end do
        call assert_real(0.0, maxdiff, 1.0e-5,                'sauvola: local standard deviations match the window statistics')
        call assert_true(sdevs%get([BOX/2,BOX/2,1]) > 0.4,    'sauvola: the standard deviation peaks at the step')
        call assert_true(sdevs%get([4,BOX/2,1]) < 0.02,       'sauvola: the standard deviation is small in a flat region')
        nmismatch = 0
        do j = 1,BOX
            do i = 1,BOX
                t = avgs(i,j) * (1.0 + BIAS * (sd(i,j)/smax - 1.0))
                if( abs(ref(i,j) - t) < 1.0e-4 ) cycle   ! on the threshold: rounding decides
                if( (ref(i,j) >= t) .neqv. (img%get([i,j,1]) > 0.5) ) nmismatch = nmismatch + 1
            end do
        end do
        call assert_int(0, nmismatch, 'sauvola: binarisation follows the Sauvola threshold pixel by pixel')
        call assert_true(all_binary(img), 'sauvola: output is 0/1')
        call img%kill
        call sdevs%kill
    end subroutine test_sauvola

    !==================================================================
    ! edge detection
    !==================================================================

    !---------------- calc_gradient and sobel ----------------

    ! the Sobel gradient of a unit ramp is 1 away from the borders; on a square the
    ! thresholded gradient is a ring along the edge
    subroutine test_sobel()
        integer, parameter :: BOX = 64, HALF = 12
        type(image) :: img
        real    :: grad(BOX,BOX,1), v
        integer :: i, j, c, nedge, nfar, d
        write(*,'(A)') 'test_sobel'
        call img%new([BOX,BOX,1], 1.0)
        do j = 1,BOX
            do i = 1,BOX
                call img%set([i,j,1], real(i))
            end do
        end do
        call img%calc_gradient(grad)
        call assert_real(1.0, grad(BOX/2,BOX/2,1), 1.0e-5,          'calc_gradient: unit ramp has gradient 1 in the interior')
        call assert_real(1.0, maxval(grad(3:BOX-2,3:BOX-2,1)), 1.0e-5, 'calc_gradient: unit ramp gradient is 1 everywhere inside')
        call assert_real(1.0, minval(grad(3:BOX-2,3:BOX-2,1)), 1.0e-5, 'calc_gradient: unit ramp gradient has no dips inside')
        ! square: sobel with a mid threshold marks the edge only
        c = BOX/2 + 1
        img = 0.0
        do j = c-HALF, c+HALF-1
            do i = c-HALF, c+HALF-1
                call img%set([i,j,1], 1.0)
            end do
        end do
        call sobel(img, [0.25])
        call assert_true(all_binary(img), 'sobel: output is 0/1')
        nedge = 0
        nfar  = 0
        do j = 1,BOX
            do i = 1,BOX
                v = img%get([i,j,1])
                if( v < 0.5 ) cycle
                ! distance to the nearest square edge line
                d = min(abs(i-(c-HALF)), abs(i-(c+HALF-1)), abs(j-(c-HALF)), abs(j-(c+HALF-1)))
                if( d <= 2 )then
                    nedge = nedge + 1
                else
                    nfar = nfar + 1
                endif
            end do
        end do
        call assert_true(nedge >= 4*(2*HALF) - 8, 'sobel: the square edge is detected all the way round')
        call assert_int(0, nfar,                 'sobel: nothing is detected away from the edge')
        call assert_real(0.0, img%get([c,c,1]),   1.0e-6, 'sobel: the interior of the square is 0')
        call assert_real(0.0, img%get([2,2,1]),   1.0e-6, 'sobel: the exterior is 0')
        call img%kill
    end subroutine test_sobel

    !---------------- canny ----------------

    ! with explicit thresholds canny returns a thin closed edge along the square
    subroutine test_canny()
        integer, parameter :: BOX = 64, HALF = 14
        type(image) :: img, edges
        real    :: v
        integer :: i, j, c, nedge, nfar, d, nfg
        write(*,'(A)') 'test_canny'
        c = BOX/2 + 1
        call img%new([BOX,BOX,1], 1.0)
        img = 0.0
        do j = c-HALF, c+HALF-1
            do i = c-HALF, c+HALF-1
                call img%set([i,j,1], 1.0)
            end do
        end do
        call canny(img, edges, thresh=[0.02, 0.06])
        call assert_true(all_binary(edges),          'canny: output is 0/1')
        call assert_real(1.0, img%get([c,c,1]), 1.0e-6, 'canny: the input is untouched when an output image is given')
        nfg   = edges%nforeground()
        nedge = 0
        nfar  = 0
        do j = 1,BOX
            do i = 1,BOX
                v = edges%get([i,j,1])
                if( v < 0.5 ) cycle
                d = min(abs(i-(c-HALF)), abs(i-(c+HALF-1)), abs(j-(c-HALF)), abs(j-(c+HALF-1)))
                if( d <= 3 )then
                    nedge = nedge + 1
                else
                    nfar = nfar + 1
                endif
            end do
        end do
        call assert_true(nfg > 0,                    'canny: an edge is found')
        call assert_int(0, nfar,                     'canny: nothing is detected away from the edge')
        call assert_true(nfg <= 3*4*(2*HALF),        'canny: the edge is thin (at most three pixels wide)')
        call assert_true(nedge >= 2*(2*HALF),        'canny: the edge runs around at least half the square')
        call assert_real(0.0, edges%get([c,c,1]), 1.0e-6, 'canny: the interior is 0')
        call img%kill
        call edges%kill
    end subroutine test_canny

    !==================================================================
    ! peak thresholds
    !==================================================================

    ! 200 background scores in [0,1) and 20 peaks in [5,6): every variant must
    ! separate the peaks, and the count-based variants return exactly their count
    subroutine test_detect_peak_thres()
        integer, parameter :: NBG = 200, NPK = 20, N = NBG + NPK
        real    :: x(N), t, t1, t2, t3
        integer :: i
        write(*,'(A)') 'test_detect_peak_thres'
        do i = 1,NBG
            x(i) = real(i-1) / real(NBG)              ! 0 .. 0.995, distinct
        end do
        do i = 1,NPK
            x(NBG+i) = 5.0 + real(i-1) / real(NPK)    ! 5 .. 5.95, distinct
        end do
        ! upper bound only (level 0): the n_ub largest
        call detect_peak_thres(N, 50, 0, x, t)
        call assert_int(50, count(x >= t),      'detect_peak_thres(level 0): threshold keeps n_ub scores')
        ! upper bound then Otsu (level 1): the peaks
        call detect_peak_thres(N, 50, 1, x, t)
        call assert_int(NPK, count(x >= t),     'detect_peak_thres(level 1): Otsu on the upper bound keeps the peaks')
        ! no upper bound, Otsu (level 1) and twice Otsu (level 2)
        call detect_peak_thres(N, 1, x, t)
        call assert_int(NPK, count(x >= t),     'detect_peak_thres(Otsu): keeps the peaks')
        call detect_peak_thres(N, 2, x, t)
        call assert_true(count(x >= t) <= NPK .and. count(x >= t) >= 1,&
            &'detect_peak_thres(twice Otsu): splits the peaks, keeps some')
        ! a fixed number of peaks
        call detect_peak_thres_for_npeaks(N, 7, x, t)
        call assert_int(7, count(x >= t),       'detect_peak_thres_for_npeaks: exactly npeaks scores at or above')
        ! sortmeans refinement: levels 1..3 admit non-decreasing numbers of peaks, thresholds are data values
        t1 = t; call refine_peak_thres_sortmeans(N, 1, x, t1)
        t2 = t; call refine_peak_thres_sortmeans(N, 2, x, t2)
        t3 = t; call refine_peak_thres_sortmeans(N, 3, x, t3)
        call assert_true(t1 >= t2 .and. t2 >= t3,  'refine_peak_thres_sortmeans: level 1 to 3 admits more peaks')
        call assert_true(any(abs(x - t1) < 1.0e-6) .and. any(abs(x - t3) < 1.0e-6),&
            &'refine_peak_thres_sortmeans: thresholds are data values')
    end subroutine test_detect_peak_thres

    !==================================================================
    ! fixtures and helpers
    !==================================================================

    ! background level plus a disc at the object level, with Gaussian noise
    subroutine two_level( img, box, rad, bg, fg, noise, ndisc )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: box
        real,        intent(in)    :: rad, bg, fg, noise
        integer,     intent(out)   :: ndisc
        integer :: i, j, c
        real    :: d, v
        c = box/2 + 1
        call img%new([box,box,1], 1.0)
        ndisc = 0
        do j = 1,box
            do i = 1,box
                d = sqrt(real((i-c)**2 + (j-c)**2))
                if( d <= rad )then
                    v = fg; ndisc = ndisc + 1
                else
                    v = bg
                endif
                if( noise > 0.0 ) v = v + gasdev(0.0, noise)
                call img%set([i,j,1], v)
            end do
        end do
    end subroutine two_level

    subroutine three_level( img, box, rad, rad_core, levels, noise, ncore )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: box
        real,        intent(in)    :: rad, rad_core, levels(3), noise
        integer,     intent(out)   :: ncore
        integer :: i, j, c
        real    :: d, v
        c = box/2 + 1
        call img%new([box,box,1], 1.0)
        ncore = 0
        do j = 1,box
            do i = 1,box
                d = sqrt(real((i-c)**2 + (j-c)**2))
                if( d <= rad_core )then
                    v = levels(3); ncore = ncore + 1
                else if( d <= rad )then
                    v = levels(2)
                else
                    v = levels(1)
                endif
                if( noise > 0.0 ) v = v + gasdev(0.0, noise)
                call img%set([i,j,1], v)
            end do
        end do
    end subroutine three_level

    subroutine four_level( img, box, rad, rad_mid, rad_core, levels, noise, ncore )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: box
        real,        intent(in)    :: rad, rad_mid, rad_core, levels(4), noise
        integer,     intent(out)   :: ncore
        integer :: i, j, c
        real    :: d, v
        c = box/2 + 1
        call img%new([box,box,1], 1.0)
        ncore = 0
        do j = 1,box
            do i = 1,box
                d = sqrt(real((i-c)**2 + (j-c)**2))
                if( d <= rad_core )then
                    v = levels(4); ncore = ncore + 1
                else if( d <= rad_mid )then
                    v = levels(3)
                else if( d <= rad )then
                    v = levels(2)
                else
                    v = levels(1)
                endif
                if( noise > 0.0 ) v = v + gasdev(0.0, noise)
                call img%set([i,j,1], v)
            end do
        end do
    end subroutine four_level

    logical function all_binary( img )
        type(image), intent(in) :: img
        integer :: i, j, ldim(3)
        real    :: v
        ldim = img%get_ldim()
        all_binary = .true.
        do j = 1,ldim(2)
            do i = 1,ldim(1)
                v = img%get([i,j,1])
                if( abs(v) > 1.0e-6 .and. abs(v - 1.0) > 1.0e-6 ) all_binary = .false.
            end do
        end do
    end function all_binary

end module simple_segmentation_tester
