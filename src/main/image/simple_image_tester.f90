!@descr: unit tests for the image class and the Fourier projector (simple_image, simple_projector)
! Replaces the self-test test_image that lived in simple_image (THROW_HARD checks, print-only
! filter, mask, rotation and binarisation parts, six image files left behind) and takes over the
! image basics that only the deleted ptcl_center test touched: get_nyq, masscen, shift2Dserial,
! roavg, power_spectrum and fproject (plan, section 9.7, the open items, 2026-09-25). Expected values
! are closed forms or exact array operations written here; the image files are removed at the end.
module simple_image_tester
use, intrinsic :: ieee_exceptions, only: ieee_get_flag, ieee_set_flag, ieee_divide_by_zero
use simple_test_utils
use simple_defs
use simple_linalg,       only: hyp
use simple_string,       only: string
use simple_string_utils, only: int2str
use simple_syslib,       only: del_file
use simple_image,        only: image
use simple_projector,    only: projector
use simple_ori,          only: ori
implicit none
private
public :: run_all_image_tests

integer, parameter :: BOX  = 64
real,    parameter :: SMPD = 1.0
character(len=*), parameter :: STK_SPI  = 'tmp_image_tester_squares.spi'
character(len=*), parameter :: STK_MRC  = 'tmp_image_tester_squares.mrc'
character(len=*), parameter :: STK_SPI2 = 'tmp_image_tester_squares_converted.spi'
character(len=*), parameter :: STK_MRC2 = 'tmp_image_tester_squares_converted.mrc'
character(len=*), parameter :: VOL_SPI  = 'tmp_image_tester_cube.spi'
character(len=*), parameter :: VOL_MRC  = 'tmp_image_tester_cube.mrc'

contains

    subroutine run_all_image_tests()
        write(*,'(A)') '**** running all image tests ****'
        call test_construct_and_access()
        call test_checkups_and_stats()
        call test_fft_roundtrip_and_nyquist()
        call test_bandpass()
        call test_apply_filter()
        call test_power_spectrum()
        call test_shifts()
        call test_acf()
        call test_rotations_and_roavg()
        call test_masscen()
        call test_corr()
        call test_file_roundtrip()
        call test_file_roundtrip_sizes()
        call test_fproject()
    end subroutine run_all_image_tests

    !---------------- construction and access ----------------

    subroutine test_construct_and_access()
        integer, parameter :: B3 = 24
        type(image) :: img, vol
        integer :: i, j, k, cnt
        logical :: ok2, ok3
        write(*,'(A)') 'test_construct_and_access'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call assert_true(img%exists(), 'a new 2D image exists')
        call assert_true(vol%exists(), 'a new 3D image exists')
        ok2 = .true.
        ok3 = .true.
        cnt = 0
        do i = 1,BOX
            do j = 1,BOX
                cnt = cnt + 1
                call img%set([i,j,1], real(cnt))
                if( img%get([i,j,1]) /= real(cnt) ) ok2 = .false.
            end do
        end do
        cnt = 0
        do i = 1,B3
            do j = 1,B3
                do k = 1,B3
                    cnt = cnt + 1
                    call vol%set([i,j,k], real(cnt))
                    if( vol%get([i,j,k]) /= real(cnt) ) ok3 = .false.
                end do
            end do
        end do
        call assert_true(ok2, 'set/get round trip of every pixel of a 2D image')
        call assert_true(ok3, 'set/get round trip of every voxel of a 3D image')
        call assert_real(real(BOX*BOX), maxval(img%get_rmat()), 0., 'get_rmat returns what set wrote')
        call img%kill
        call vol%kill
        call assert_false(img%exists(), 'a killed 2D image does not exist')
        call assert_false(vol%exists(), 'a killed 3D image does not exist')
    end subroutine test_construct_and_access

    !---------------- dimension checks and statistics ----------------

    ! foreground statistics of N(5,15) noise inside a disc of radius 40: about 5000 pixels, so the
    ! mean and median are within 1 of 5 and the standard deviation within 1 of 15 (5 standard errors)
    subroutine test_checkups_and_stats()
        integer, parameter :: BOXS = 100
        type(image) :: img, img2, rect
        real :: ave, sdev, maxv, minv, med
        write(*,'(A)') 'test_checkups_and_stats'
        call img%new([BOXS,BOXS,1], SMPD, wthreads=.false.)
        call img2%copy(img)
        call rect%new([BOXS,BOXS-20,1], SMPD, wthreads=.false.)
        call assert_true(img%even_dims(),          'a 100x100 image has even dimensions')
        call assert_true(img%square_dims(),        'a 100x100 image is square')
        call assert_false(rect%square_dims(),      'a 100x80 image is not square')
        call assert_true(img .eqdims. img2,        'a copy has the same dimensions')
        call assert_false(img .eqdims. rect,       'a 100x80 image differs in dimensions')
        call assert_true(img%is_2d(),              'a single-slice image is 2D')
        call assert_false(img%is_3d(),             'a single-slice image is not 3D')
        call set_fixed_seed(20260925)
        call img%gauran(5., 15.)
        call img%stats('foreground', ave, sdev, maxv, minv, 40., med)
        call assert_real(5.,  ave,  1., 'foreground mean of N(5,15) noise')
        call assert_real(15., sdev, 1., 'foreground standard deviation of N(5,15) noise')
        call assert_real(5.,  med,  1., 'foreground median of N(5,15) noise')
        call assert_true(maxv >= ave .and. minv <= ave, 'foreground extremes bracket the mean')
        call img%kill
        call img2%kill
        call rect%kill
    end subroutine test_checkups_and_stats

    !---------------- FFT round trip, Nyquist ----------------

    subroutine test_fft_roundtrip_and_nyquist()
        type(image) :: img, odd
        real, allocatable :: before(:,:,:), after(:,:,:)
        write(*,'(A)') 'test_fft_roundtrip_and_nyquist'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call odd%new([BOX-1,BOX-1,1], SMPD, wthreads=.false.)
        call assert_int(BOX/2,     img%get_nyq(), 'get_nyq of an even box is box/2')
        call assert_int(31,        odd%get_nyq(), 'get_nyq of an odd box (63) is 31')
        call set_fixed_seed(20260925)
        call img%gauran(0., 1.)
        before = img%get_rmat()
        call img%fft
        call img%ifft
        after = img%get_rmat()
        call assert_true(maxval(abs(after - before)) < 1.e-5 * maxval(abs(before)), 'fft followed by ifft returns the image')
        call img%kill
        call odd%kill
    end subroutine test_fft_roundtrip_and_nyquist

    !---------------- band-pass ----------------

    ! a flat spectrum (every component 1) through bp: the low-pass zeroes every component beyond
    ! its shell and keeps every one below the cosine edge, the high-pass the reverse; a limit of 0
    ! switches its side off without the division by zero that get_find(0.) was (fixed 2026-09-25)
    subroutine test_bandpass()
        real,    parameter :: LP = 4., HP = 8., WIDTH = 10.
        type(image) :: img
        logical     :: flag
        integer     :: lpf, hpf
        write(*,'(A)') 'test_bandpass'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        lpf = img%get_find(LP)
        hpf = img%get_find(HP)
        call flat_spectrum(img)
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        call img%bp(0., LP)
        call ieee_get_flag(ieee_divide_by_zero, flag)
        call assert_false(flag, 'bp(0., lp) raises no IEEE division by zero')
        call assert_true(band_equals(img, real(lpf) + 0.5, 1.e6, 0.), 'low-pass: every component beyond the limit is zero')
        call assert_true(band_equals(img, 0., real(lpf) - WIDTH - 0.5, 1.), 'low-pass: every component below the edge is kept')
        call flat_spectrum(img)
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        call img%bp(HP, 0.)
        call ieee_get_flag(ieee_divide_by_zero, flag)
        call assert_false(flag, 'bp(hp, 0.) raises no IEEE division by zero')
        call assert_true(band_equals(img, 0., real(hpf) - 0.5, 0.), 'high-pass: every component below the limit is zero')
        call assert_true(band_equals(img, real(hpf) + WIDTH + 0.5, 1.e6, 1.), 'high-pass: every component beyond the edge is kept')
        call ieee_set_flag(ieee_divide_by_zero, .false.)
        call img%kill
    end subroutine test_bandpass

    !---------------- apply_filter ----------------

    ! every component 1, filtered by 1/shell: DC stays 1, shell s becomes 1/s, beyond the filter 0
    subroutine test_apply_filter()
        integer, parameter :: B3 = 32
        type(image)       :: vol
        real, allocatable :: filter(:)
        integer :: h, k, l, sh, sz, lims(3,2), phys(3)
        complex :: expected
        logical :: ok
        write(*,'(A)') 'test_apply_filter'
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call vol%set_ft(.true.)
        call vol%set_cmat(cmplx(1.,0.))
        sz = vol%get_filtsz()
        allocate(filter(sz))
        do sh = 1,sz
            filter(sh) = 1. / real(sh)
        end do
        call vol%apply_filter(filter)
        lims = vol%loop_lims(2)
        ok = .true.
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                do l = lims(3,1),lims(3,2)
                    sh   = nint(hyp(h,k,l))
                    phys = vol%comp_addr_phys(h,k,l)
                    if( sh == 0 )then
                        expected = cmplx(1.,0.)
                    else if( sh > sz )then
                        expected = cmplx(0.,0.)
                    else
                        expected = cmplx(filter(sh),0.)
                    endif
                    if( abs(vol%get_cmat_at(phys(1),phys(2),phys(3)) - expected) > 1.e-6 ) ok = .false.
                end do
            end do
        end do
        call assert_true(ok, 'apply_filter scales every component by the filter of its shell')
        call vol%kill
    end subroutine test_apply_filter

    !---------------- power spectrum ----------------

    ! a plane wave cos(2 pi K0 x / box) has its power in shell K0 only
    subroutine test_power_spectrum()
        integer, parameter :: K0 = 5
        type(image)       :: img
        real, allocatable :: rmat(:,:,:)
        real    :: spec(BOX/2)
        integer :: i
        write(*,'(A)') 'test_power_spectrum'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        allocate(rmat(BOX,BOX,1))
        do i = 1,BOX
            rmat(i,:,1) = cos(TWOPI * real(K0 * (i - 1)) / real(BOX))
        end do
        call img%set_rmat(rmat, .false.)
        call img%fft
        call img%power_spectrum(spec)
        call assert_true(spec(K0) > 0.,                             'power spectrum: the plane-wave shell carries power')
        call assert_true(maxval(spec, mask=[(i /= K0, i=1,BOX/2)]) < 1.e-6 * spec(K0), &
            &'power spectrum: every other shell is empty')
        call img%kill
    end subroutine test_power_spectrum

    !---------------- shifts ----------------

    ! the in-place and out-of-place serial Fourier shifts agree with the general image%shift; an
    ! integer shift is an exact circular shift of the pixels
    subroutine test_shifts()
        real, parameter :: SH(2) = [2.5, -1.25]
        type(image)       :: orig, a, b, c
        real, allocatable :: r0(:,:,:), ra(:,:,:), rb(:,:,:), rc(:,:,:)
        real    :: err_plus, err_minus, amax
        write(*,'(A)') 'test_shifts'
        call smooth_noise(orig, 20260925)
        r0    = orig%get_rmat()
        amax = maxval(abs(r0))
        call a%copy(orig)
        call a%shift([SH(1), SH(2), 0.])
        ra = a%get_rmat()
        call b%copy(orig)
        call b%fft
        call b%shift2Dserial(SH)
        call b%ifft
        rb = b%get_rmat()
        call assert_true(maxval(abs(rb - ra)) < 1.e-4 * amax, 'shift2Dserial (in place) agrees with shift')
        call b%copy(orig)
        call b%fft
        call c%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call c%set_ft(.true.)
        call b%shift2Dserial(SH, c)
        call c%ifft
        rc = c%get_rmat()
        call assert_true(maxval(abs(rc - ra)) < 1.e-4 * amax, 'shift2Dserial (out of place) agrees with shift')
        ! integer shift: shift(s) gives out(x) = in(x + s), circularly, i.e. the content moves by -s
        ! (the opposite sign misses by 0.5 of the maximum; first run, 2026-09-25)
        call a%copy(orig)
        call a%shift([3., -5., 0.])
        ra = a%get_rmat()
        err_plus  = maxval(abs(ra(:,:,1) - cshift(cshift(r0(:,:,1), -3, dim=1),  5, dim=2)))
        err_minus = maxval(abs(ra(:,:,1) - cshift(cshift(r0(:,:,1),  3, dim=1), -5, dim=2)))
        write(logfhandle,'(A,2ES10.2)') 'integer shift, error against the content moved by +s and by -s: ', err_plus, err_minus
        call assert_true(err_minus < 1.e-4 * amax, 'shift by an integer s is the circular shift out(x) = in(x + s)')
        call orig%kill
        call a%kill
        call b%kill
        call c%kill
    end subroutine test_shifts

    !---------------- autocorrelation ----------------

    ! the autocorrelation peaks at the origin, the centre pixel, and does not see a circular shift
    subroutine test_acf()
        type(image)       :: a, b
        real, allocatable :: ra(:,:,:), rb(:,:,:)
        integer :: loc(3)
        write(*,'(A)') 'test_acf'
        call a%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call a%square(10)
        call b%copy(a)
        call b%shift([5., -5., 0.])
        call a%acf
        call b%acf
        ra  = a%get_rmat()
        rb  = b%get_rmat()
        loc = maxloc(ra)
        call assert_int(BOX/2+1, loc(1), 'acf peaks at the centre (x)')
        call assert_int(BOX/2+1, loc(2), 'acf peaks at the centre (y)')
        call assert_true(maxval(abs(ra - rb)) < 1.e-4 * maxval(abs(ra)), 'acf is invariant to a circular shift')
        call a%kill
        call b%kill
    end subroutine test_acf

    !---------------- rotations, rotational average ----------------

    subroutine test_rotations_and_roavg()
        integer, parameter :: C = BOX/2 + 1, RAD = 24
        type(image)       :: orig, rot, back, gau, avg, sq
        real, allocatable :: r0(:,:,:), r1(:,:,:), r2(:,:,:)
        real    :: err_ccw, err_cw, amax
        integer :: i, j, r
        logical :: sym_ok
        write(*,'(A)') 'test_rotations_and_roavg'
        ! a smooth image confined to a disc, so that a rotation keeps it inside the box
        call smooth_noise(orig, 20260926)
        r0 = orig%get_rmat()
        do i = 1,BOX
            do j = 1,BOX
                if( (i-C)**2 + (j-C)**2 > RAD**2 ) r0(i,j,1) = 0.
            end do
        end do
        call orig%set_rmat(r0, .false.)
        amax = maxval(abs(r0))
        call rot%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call orig%rtsq(0., 0., 0., rot)
        r1 = rot%get_rmat()
        call assert_true(maxval(abs(r1 - r0)) < 1.e-4 * amax, 'rtsq by 0 degrees is the identity')
        ! a quarter turn about the centre pixel lands on the grid: exact up to the interpolation
        call orig%rtsq(90., 0., 0., rot)
        r1 = rot%get_rmat()
        err_ccw = 0.
        err_cw  = 0.
        do i = C-RAD,C+RAD
            do j = C-RAD,C+RAD
                err_ccw = max(err_ccw, abs(r1(i,j,1) - r0(j, 2*C-i, 1)))
                err_cw  = max(err_cw,  abs(r1(i,j,1) - r0(2*C-j, i, 1)))
            end do
        end do
        write(logfhandle,'(A,2ES10.2)') 'rtsq by 90 degrees, error against the two quarter turns: ', err_ccw, err_cw
        call assert_true(err_cw < 1.e-3 * amax, 'rtsq by 90 degrees is the quarter turn out(i,j) = in(2c-j,i) about the centre c')
        ! a rotation and its inverse recover the image away from the sharp disc edge, which the
        ! quadratic interpolation smears (a numpy emulation of rtsq gives 0.9999 inside radius 20,
        ! 0.978 over the whole box)
        call back%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call orig%rtsq(37., 0., 0., rot)
        call rot%rtsq(-37., 0., 0., back)
        r2 = back%get_rmat()
        call assert_true(pearson(pack(r0(:,:,1), in_disc(20)), pack(r2(:,:,1), in_disc(20))) > 0.999, &
            &'rtsq by 37 and -37 degrees recovers the image inside radius 20')
        ! rotational average: an isotropic Gaussian is its own average inside the inscribed circle
        ! (the corners take in the circular wrap of the rotations: 0.022 at the corner pixel,
        ! 3e-5 inside radius 31 in the emulation); a square becomes isotropic
        call gau%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call gau%gauimg(10)
        r0 = gau%get_rmat()
        call gau%roavg(10, avg)
        r1 = avg%get_rmat()
        call assert_true(maxval(abs(r1(:,:,1) - r0(:,:,1)), mask=in_disc(30)) < 1.e-3 * maxval(r0), &
            &'roavg of an isotropic Gaussian is the Gaussian inside radius 30')
        call sq%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call sq%square(10)
        r0 = sq%get_rmat()
        call sq%roavg(5, avg)
        r1 = avg%get_rmat()
        sym_ok = .true.
        do r = 3,12,9 ! inside the square's inscribed circle about the centre, and where the circle cuts it
            if( maxval([r1(C+r,C,1), r1(C,C+r,1), r1(C-r,C,1), r1(C,C-r,1)]) - &
               &minval([r1(C+r,C,1), r1(C,C+r,1), r1(C-r,C,1), r1(C,C-r,1)]) > 2.e-2 * maxval(r0) ) sym_ok = .false.
        end do
        call assert_true(sym_ok, 'roavg of a square is the same at radius r along both axes, both signs')
        call assert_real(sum(r0), sum(r1), 0.02 * sum(r0), 'roavg keeps the integral of a square inside the box')
        call orig%kill
        call rot%kill
        call back%kill
        call gau%kill
        call avg%kill
        call sq%kill

      contains

        ! pixels within radius rad of the centre pixel
        function in_disc( rad ) result( msk )
            integer, intent(in) :: rad
            logical :: msk(BOX,BOX)
            integer :: ii, jj
            do jj = 1,BOX
                do ii = 1,BOX
                    msk(ii,jj) = (ii-C)**2 + (jj-C)**2 <= rad**2
                end do
            end do
        end function in_disc

    end subroutine test_rotations_and_roavg

    !---------------- centre of mass ----------------

    ! masscen reports positions relative to the centre pixel box/2+1: pixel i sits at i-1-box/2
    subroutine test_masscen()
        integer, parameter :: B3 = 24
        type(image)          :: img, vol
        real,    allocatable :: rmat(:,:,:), vmat(:,:,:)
        logical, allocatable :: msk(:,:,:)
        real :: xyz(3)
        write(*,'(A)') 'test_masscen'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        allocate(rmat(BOX,BOX,1), source=0.)
        call img%set_rmat(rmat, .false.)
        call img%masscen(xyz)
        call assert_true(all(xyz == 0.), 'masscen of an empty image is the origin')
        rmat(40,20,1) = 1.
        call img%set_rmat(rmat, .false.)
        call img%masscen(xyz)
        call assert_true(all(abs(xyz - [7., -13., 0.]) < 1.e-5), 'masscen of one pixel is its offset from the centre')
        rmat(24,44,1) = 3.
        call img%set_rmat(rmat, .false.)
        call img%masscen(xyz)
        call assert_true(all(abs(xyz - [-5., 5., 0.]) < 1.e-5), 'masscen of two pixels is their weighted mean offset')
        allocate(msk(BOX,BOX,1), source=.true.)
        msk(24,44,1) = .false.
        call img%masscen(xyz, msk)
        call assert_true(all(abs(xyz - [7., -13., 0.]) < 1.e-5), 'masscen ignores pixels outside mask_in')
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        allocate(vmat(B3,B3,B3), source=0.)
        vmat(10,20,5) = 2.
        call vol%set_rmat(vmat, .false.)
        call vol%masscen(xyz)
        call assert_true(all(abs(xyz - [-3., 7., -8.]) < 1.e-5), 'masscen of one voxel is its offset from the centre')
        call img%kill
        call vol%kill
    end subroutine test_masscen

    !---------------- correlation ----------------

    ! two centred Gaussians of widths in the ratio 10:13 correlate at 2*10*13/(10**2+13**2) = 0.9665
    ! (continuous closed form). corr leaves out the lowest Fourier indices (|h|**2 < 2), which weigh
    ! more in a small box: a numpy emulation of it gives 0.9672 in the box of 100 of the original
    ! test_image and 0.9592 in a box of 64 (first run, 2026-09-25); low-passed at 20 A the value
    ! depends on where loop_lims starts h (0.9636-0.9672), so that one keeps the old range
    subroutine test_corr()
        integer, parameter :: BOXC = 100
        type(image) :: g1, g2, g3
        real :: cc, cc_lp, cc_sym
        write(*,'(A)') 'test_corr'
        call g1%new([BOXC,BOXC,1], 2., wthreads=.false.)
        call g2%new([BOXC,BOXC,1], 2., wthreads=.false.)
        call g1%gauimg(10)
        call g2%gauimg(13)
        call g1%fft
        call g2%fft
        cc     = g1%corr(g2)
        cc_lp  = g1%corr(g2, 20.)
        cc_sym = g2%corr(g1)
        call assert_real(0.9672, cc, 0.002,               'corr of Gaussians of widths 10 and 13')
        call assert_true(cc_lp > 0.96 .and. cc_lp < 0.98, 'corr of Gaussians of widths 10 and 13, low-passed at 20 A')
        call assert_real(cc, cc_sym, 1.e-6,               'corr is symmetric')
        call g3%copy(g1)
        call assert_real(1., g1%corr(g3), 1.e-5,          'corr of an image with its copy is 1')
        call g1%kill
        call g2%kill
        call g3%kill
    end subroutine test_corr

    !---------------- file I/O ----------------

    ! five images through SPIDER and MRC stacks and their conversions, and a volume through both
    ! formats, come back bit for bit
    subroutine test_file_roundtrip()
        integer, parameter :: NIMGS = 5, B3 = 32
        type(image)       :: img, back, vol, vback
        real, allocatable :: ref(:,:,:), got(:,:,:)
        integer :: i
        logical :: ok
        write(*,'(A)') 'test_file_roundtrip'
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call back%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call del_all
        do i = 1,NIMGS
            call img%square(4 + 2*i)
            call img%write(string(STK_SPI), i)
            call img%write(string(STK_MRC), i)
        end do
        do i = 1,NIMGS
            call back%read(string(STK_SPI), i)
            call back%write(string(STK_MRC2), i)
            call back%read(string(STK_MRC), i)
            call back%write(string(STK_SPI2), i)
        end do
        ok = .true.
        do i = 1,NIMGS
            call img%square(4 + 2*i)
            ref = img%get_rmat()
            call back%read(string(STK_SPI), i)
            got = back%get_rmat()
            if( any(got /= ref) ) ok = .false.
            call back%read(string(STK_MRC), i)
            got = back%get_rmat()
            if( any(got /= ref) ) ok = .false.
            call back%read(string(STK_SPI2), i)
            got = back%get_rmat()
            if( any(got /= ref) ) ok = .false.
            call back%read(string(STK_MRC2), i)
            got = back%get_rmat()
            if( any(got /= ref) ) ok = .false.
        end do
        call assert_true(ok, 'SPIDER and MRC stacks and their conversions read back every image exactly')
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call vback%new([B3,B3,B3], SMPD, wthreads=.false.)
        call vol%square(6)
        ref = vol%get_rmat()
        call vol%write(string(VOL_SPI))
        call vol%write(string(VOL_MRC))
        call vback%read(string(VOL_SPI))
        got = vback%get_rmat()
        ok  = all(got == ref)
        call vback%read(string(VOL_MRC))
        got = vback%get_rmat()
        ok  = ok .and. all(got == ref)
        call assert_true(ok, 'a volume reads back exactly from SPIDER and MRC')
        call del_all
        call img%kill
        call back%kill
        call vol%kill
        call vback%kill

      contains

        subroutine del_all
            call del_file(STK_SPI)
            call del_file(STK_MRC)
            call del_file(STK_SPI2)
            call del_file(STK_MRC2)
            call del_file(VOL_SPI)
            call del_file(VOL_MRC)
        end subroutine del_all

    end subroutine test_file_roundtrip

    ! round trips at the sizes where a file layout changes (policy, section 4.5): a SPIDER header is
    ! labrec records of 4*nx bytes and at least 1024 bytes, so boxes below 13 and below 43 and boxes
    ! above 256 lay it out differently (it was wrong for every box below 43 until 2026-09-25, while
    ! the tests used 64); odd and non-square images and volumes, stacks of three, in both formats.
    ! The pixel values are integers below 2**24, exact in single precision and distinct per pixel
    ! and per image, so a transposition or a swapped image shows
    subroutine test_file_roundtrip_sizes()
        integer, parameter :: N2 = 6, N3 = 3, NSTK = 3
        integer, parameter :: DIMS2(2,N2) = reshape([5,5, 13,13, 42,42, 43,43, 43,17, 257,9], [2,N2])
        integer, parameter :: DIMS3(3,N3) = reshape([5,5,5, 13,17,5, 43,43,43], [3,N3])
        character(len=4), parameter :: EXTS(2) = ['.mrc', '.spi']
        type(image) :: img, back
        character(len=:), allocatable :: fname, tag
        integer :: id, ie, i, ldim(3)
        logical :: ok
        write(*,'(A)') 'test_file_roundtrip_sizes'
        do ie = 1,size(EXTS)
            fname = 'tmp_image_tester_sizes'//EXTS(ie)
            do id = 1,N2
                ldim = [DIMS2(1,id), DIMS2(2,id), 1]
                tag  = EXTS(ie)//' '//int2str(ldim(1))//'x'//int2str(ldim(2))
                call del_file(fname)
                call img%new(ldim, SMPD, wthreads=.false.)
                call back%new(ldim, SMPD, wthreads=.false.)
                do i = 1,NSTK
                    call img%set_rmat(pattern(ldim, i), .false.)
                    call img%write(string(fname), i)
                end do
                ok = .true.
                do i = 1,NSTK
                    call back%read(string(fname), i)
                    if( any(back%get_rmat() /= pattern(ldim, i)) ) ok = .false.
                end do
                call assert_true(ok, tag//': a stack of three images reads back exactly')
            end do
            do id = 1,N3
                ldim = DIMS3(:,id)
                tag  = EXTS(ie)//' '//int2str(ldim(1))//'x'//int2str(ldim(2))//'x'//int2str(ldim(3))
                call del_file(fname)
                call img%new(ldim, SMPD, wthreads=.false.)
                call back%new(ldim, SMPD, wthreads=.false.)
                call img%set_rmat(pattern(ldim, 1), .false.)
                call img%write(string(fname))
                call back%read(string(fname))
                call assert_true(all(back%get_rmat() == pattern(ldim, 1)), tag//': a volume reads back exactly')
            end do
            call del_file(fname)
        end do
        call img%kill
        call back%kill

      contains

        function pattern( ldim, iimg ) result( p )
            integer, intent(in) :: ldim(3), iimg
            real    :: p(ldim(1),ldim(2),ldim(3))
            integer :: i, j, k
            do k = 1,ldim(3)
                do j = 1,ldim(2)
                    do i = 1,ldim(1)
                        p(i,j,k) = real(i + 1000*j + 100000*k + 1000000*iimg)
                    end do
                end do
            end do
        end function pattern

    end subroutine test_file_roundtrip_sizes

    !---------------- Fourier projector ----------------

    ! projections of a centred isotropic Gaussian: the threaded fproject and fproject_serial give the
    ! same plane, every orientation gives the same projection, and the projection is centred
    subroutine test_fproject()
        integer, parameter :: B3 = 32, NORIS = 4
        real,    parameter :: SIGMA = 3., EULS(3,NORIS) = reshape([0.,0.,0., 30.,40.,50., 123.,77.,200., 250.,160.,15.], [3,NORIS])
        type(image)       :: vol, plane_ser, plane_par, proj2d
        type(projector)   :: vproj
        type(ori)         :: e
        real, allocatable :: vmat(:,:,:), p0(:,:,:), p(:,:,:)
        real    :: xyz(3), dsq
        integer :: i, j, k, iori, c
        logical :: same_plane
        write(*,'(A)') 'test_fproject'
        c = B3/2 + 1
        allocate(vmat(B3,B3,B3))
        do i = 1,B3
            do j = 1,B3
                do k = 1,B3
                    dsq = real((i-c)**2 + (j-c)**2 + (k-c)**2)
                    vmat(i,j,k) = exp(-dsq / (2. * SIGMA**2))
                end do
            end do
        end do
        call vol%new([B3,B3,B3], SMPD, wthreads=.false.)
        call vol%set_rmat(vmat, .false.)
        call vproj%new(OSMPL_PAD_FAC*[B3,B3,B3], SMPD, wthreads=.false.)
        call vol%pad(vproj)
        call vproj%fft()
        call vproj%expand_cmat()
        call plane_ser%new([OSMPL_PAD_FAC*B3,OSMPL_PAD_FAC*B3,1], SMPD, wthreads=.false.)
        call plane_par%new([OSMPL_PAD_FAC*B3,OSMPL_PAD_FAC*B3,1], SMPD, wthreads=.false.)
        call proj2d%new([B3,B3,1], SMPD, wthreads=.false.)
        call e%new(.false.)
        same_plane = .true.
        do iori = 1,NORIS
            call e%set_euler(EULS(:,iori))
            call vproj%fproject_serial(e, plane_ser)
            call vproj%fproject(e, plane_par)
            if( maxval(abs(plane_ser%get_cmat() - plane_par%get_cmat())) > 1.e-6 ) same_plane = .false.
            call plane_ser%ifft()
            call plane_ser%clip(proj2d)
            p = proj2d%get_rmat()
            if( iori == 1 )then
                p0 = p
                call proj2d%masscen(xyz)
                call assert_true(all(abs(xyz) < 0.05), 'the projection of a centred Gaussian is centred')
            else
                call assert_true(pearson(pack(p0, .true.), pack(p, .true.)) > 0.999, &
                    &'an isotropic volume projects the same at every orientation (orientation '//int2str(iori)//')')
            endif
        end do
        call assert_true(same_plane, 'fproject (threaded) and fproject_serial give the same plane')
        call e%kill
        call vproj%kill_expanded()
        call vproj%kill()
        call vol%kill
        call plane_ser%kill
        call plane_par%kill
        call proj2d%kill
    end subroutine test_fproject

    !---------------- helpers ----------------

    ! zero-mean Gaussian noise, Gaussian low-passed: smooth enough for interpolation and shifts
    subroutine smooth_noise( img, seed )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: seed
        call img%new([BOX,BOX,1], SMPD, wthreads=.false.)
        call set_fixed_seed(seed)
        call img%gauran(0., 1.)
        call img%fft
        call img%bpgau2D(0., 6.)
        call img%ifft
    end subroutine smooth_noise

    ! every Fourier component of img set to 1
    subroutine flat_spectrum( img )
        type(image), intent(inout) :: img
        call img%set_ft(.true.)
        call img%set_cmat(cmplx(1.,0.))
    end subroutine flat_spectrum

    ! .true. when every component with lo <= |(h,k)| <= hi has modulus val
    logical function band_equals( img, lo, hi, val )
        type(image), intent(inout) :: img
        real,        intent(in)    :: lo, hi, val
        integer :: h, k, lims(3,2), phys(3)
        real    :: r
        band_equals = .true.
        lims = img%loop_lims(2)
        do h = lims(1,1),lims(1,2)
            do k = lims(2,1),lims(2,2)
                r = hyp(real(h), real(k))
                if( r < lo .or. r > hi ) cycle
                phys = img%comp_addr_phys(h,k,0)
                if( abs(abs(img%get_cmat_at(phys(1),phys(2),phys(3))) - val) > 1.e-6 ) band_equals = .false.
            end do
        end do
    end function band_equals

    real function pearson( a, b )
        real, intent(in) :: a(:), b(:)
        real(dp) :: ma, mb
        ma = sum(real(a,dp)) / real(size(a),dp)
        mb = sum(real(b,dp)) / real(size(b),dp)
        pearson = real(sum((a-ma)*(b-mb)) / sqrt(sum((a-ma)**2) * sum((b-mb)**2)))
    end function pearson

end module simple_image_tester
