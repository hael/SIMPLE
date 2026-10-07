!@descr: unit tests for the nonuniform filter's cross-half competition: the squared-error unary and its noise level
! (image::nu_objective, nu_objective_noise_scale), the ranking's independence of that level, the label field of
! setup_nu_dmats/optimize_nu_cutoff_finds on a two-resolution phantom (a white Gaussian field smoothed finely on
! the left and coarsely on the right, unit sdev, plus independent per-half noise, more of it on the right), and
! the evidence state and envelope on the neutral fixture (a sphere of band-limited common signal in a noise
! support). Expected values: closed forms, the fixtures' construction, floors from a numpy emulation.
module simple_nu_filter_tester
use simple_core_module_api
use simple_image,     only: image
use simple_nu_filter, only: setup_nu_dmats, optimize_nu_cutoff_finds, pack_filtmap_lowpass_limits, &
                           &cleanup_nu_filter, set_nu_filter_report, build_nu_evidence_state, nu_evidence_state, &
                           &nu_evidence_summary, get_nu_evidence_summary, unpack_nu_evidence_state, &
                           &nu_evidence_state_is_valid, calc_nu_evidence_margin, nu_evidence_envelope, &
                           &nu_envmask_params, nu_envmask_stats, NU_EVIDENCE_NBANDS, NU_EVIDENCE_SOURCE_BASE, &
                           &NU_EVIDENCE_MIN_NULL_FRAC, NU_EVIDENCE_MAX_NULL_FRAC
use simple_test_utils
implicit none

public :: run_all_nu_filter_tests
private
#include "simple_local_flags.inc"

integer, parameter :: BOX        = 64
real,    parameter :: SMPD       = 2.
real,    parameter :: MSKDIAM    = 100.  !< A: a support sphere of radius 25 px
integer, parameter :: SEED       = 20261006
real,    parameter :: SIG_FINE   = 0.7   !< px: Gaussian smoothing of the white field, fine half
real,    parameter :: SIG_COARSE = 2.5   !< px: coarse half
real,    parameter :: NOISE      = 0.5   !< per-half noise sdev, unit signal sdev
real,    parameter :: EXTRA      = 0.7   !< additional per-half disagreement sdev, coarse half only
integer, parameter :: MARGIN     = 8     !< px either side of the split left out of the region statistics
! the neutral fixture: a molecule of band-limited common signal (period about five pixels, so an
! intermediate rung wins inside it) in a generous noise support; a hard radius emulates a PCG solve support
integer, parameter :: EBOX        = 48
real,    parameter :: MOL_RAD_PX  = 15.0
real,    parameter :: SUPP_RAD_PX = 22.0
real,    parameter :: HARD_RAD_PX = 18.0
real,    parameter :: ENV_NOISE   = 1.0  !< amplitude of the uniform [-0.5,0.5] per-half noise
real,    parameter :: ENV_LP      = 4.0  !< envelope scale in A, small so the floors test the statistic

contains

    subroutine run_all_nu_filter_tests()
        write(logfhandle,'(A)') '**** running all nonuniform filter tests ****'
        call test_objective_closed_form()
        call test_noise_scale()
        call test_ranking_scale_free()
        call test_coarser_region()
        call test_evidence_envelope()
    end subroutine run_all_nu_filter_tests

    ! ---- fixtures -------------------------------------------------------------

    subroutine fill_gaussian( rmat )
        real, intent(inout) :: rmat(:,:,:)
        integer :: i, j, k
        do k = 1, size(rmat,3)
            do j = 1, size(rmat,2)
                do i = 1, size(rmat,1)
                    rmat(i,j,k) = gasdev()
                end do
            end do
        end do
    end subroutine fill_gaussian

    !> white field smoothed by a separable Gaussian of sigma_px (clamp-to-edge), zero mean and unit sdev
    function smoothed_field( white, sigma_px ) result( rmat )
        real, intent(in) :: white(BOX,BOX,BOX)
        real, intent(in) :: sigma_px
        real, allocatable :: rmat(:,:,:), tmp(:,:,:), kernel(:)
        integer :: r, m
        real    :: avg, sdev
        r = ceiling(3.5 * sigma_px)
        allocate(kernel(-r:r))
        do m = -r, r
            kernel(m) = exp(-0.5 * (real(m) / sigma_px)**2)
        end do
        kernel = kernel / sum(kernel)
        allocate(rmat(BOX,BOX,BOX), tmp(BOX,BOX,BOX))
        call convolve_axis(white, tmp,  kernel, r, 1)
        call convolve_axis(tmp,   rmat, kernel, r, 2)
        call convolve_axis(rmat,  tmp,  kernel, r, 3)
        avg  = sum(tmp) / real(BOX**3)
        sdev = sqrt(sum((tmp - avg)**2) / real(BOX**3))
        rmat = (tmp - avg) / sdev
    end function smoothed_field

    subroutine convolve_axis( src, dst, kernel, r, axis )
        integer, intent(in)  :: r, axis
        real,    intent(in)  :: src(BOX,BOX,BOX), kernel(-r:r)
        real,    intent(out) :: dst(BOX,BOX,BOX)
        integer :: i, j, k, m, n
        dst = 0.
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    do m = -r, r
                        select case( axis )
                            case(1)
                                n = min(BOX, max(1, i + m))
                                dst(i,j,k) = dst(i,j,k) + kernel(m) * src(n,j,k)
                            case(2)
                                n = min(BOX, max(1, j + m))
                                dst(i,j,k) = dst(i,j,k) + kernel(m) * src(i,n,k)
                            case default
                                n = min(BOX, max(1, k + m))
                                dst(i,j,k) = dst(i,j,k) + kernel(m) * src(i,j,n)
                        end select
                    end do
                end do
            end do
        end do
    end subroutine convolve_axis

    !> even/odd = signal + independent noise. l_split: the signal is the fine field for i <= BOX/2 and the
    !! coarse field beyond, where each half also carries EXTRA independent disagreement; else the fine
    !! field everywhere under stationary noise.
    subroutine make_pair( l_split, even, odd )
        logical,     intent(in)    :: l_split
        type(image), intent(inout) :: even, odd
        real, allocatable :: white(:,:,:), s_fine(:,:,:), s_coarse(:,:,:), e(:,:,:), o(:,:,:)
        integer :: i, j, k
        call set_fixed_seed(SEED)
        allocate(white(BOX,BOX,BOX), e(BOX,BOX,BOX), o(BOX,BOX,BOX))
        call fill_gaussian(white)
        s_fine = smoothed_field(white, SIG_FINE)
        if( l_split ) s_coarse = smoothed_field(white, SIG_COARSE)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    if( l_split .and. i > BOX/2 )then
                        e(i,j,k) = s_coarse(i,j,k) + NOISE * gasdev() + EXTRA * gasdev()
                        o(i,j,k) = s_coarse(i,j,k) + NOISE * gasdev() + EXTRA * gasdev()
                    else
                        e(i,j,k) = s_fine(i,j,k) + NOISE * gasdev()
                        o(i,j,k) = s_fine(i,j,k) + NOISE * gasdev()
                    endif
                end do
            end do
        end do
        call even%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call odd%new( [BOX,BOX,BOX], SMPD, wthreads=.false.)
        call even%set_rmat(e, .false.)
        call odd%set_rmat(o,  .false.)
    end subroutine make_pair

    subroutine support_sphere( lmask )
        logical, allocatable, intent(out) :: lmask(:,:,:)
        type(image) :: vol
        call vol%disc([BOX,BOX,BOX], SMPD, 0.5 * MSKDIAM / SMPD, lmask)
        call vol%kill
    end subroutine support_sphere

    ! ---- tests ----------------------------------------------------------------

    !> the unary is ((E - O_c)^2 + (E_c - O)^2) / (2 sigma^2) on the support and zero outside it
    subroutine test_objective_closed_form()
        integer, parameter :: N     = 16
        real,    parameter :: SIGMA = 0.8
        type(image) :: e_raw, e_filt, o_raw, o_filt
        real,    allocatable :: er(:,:,:), ef(:,:,:), orw(:,:,:), of(:,:,:), diff(:,:,:)
        logical, allocatable :: lmask(:,:,:)
        real(dp) :: expected, maxerr, maxexp
        integer  :: i, j, k, n_outside
        write(logfhandle,'(A)') '-- test_objective_closed_form'
        call set_fixed_seed(SEED + 1)
        allocate(er(N,N,N), ef(N,N,N), orw(N,N,N), of(N,N,N), diff(N,N,N), lmask(N,N,N))
        call random_number(er)
        call random_number(ef)
        call random_number(orw)
        call random_number(of)
        er  = er  - 0.5
        ef  = ef  - 0.5
        orw = orw - 0.5
        of  = of  - 0.5
        do k = 1, N
            do j = 1, N
                do i = 1, N
                    lmask(i,j,k) = (i - N/2 - 1)**2 + (j - N/2 - 1)**2 + (k - N/2 - 1)**2 <= 36
                end do
            end do
        end do
        call new_from(e_raw,  er)
        call new_from(e_filt, ef)
        call new_from(o_raw,  orw)
        call new_from(o_filt, of)
        call e_raw%nu_objective(e_filt, o_raw, o_filt, diff, lmask, SIGMA)
        maxerr    = 0.d0
        maxexp    = 0.d0
        n_outside = 0
        do k = 1, N
            do j = 1, N
                do i = 1, N
                    if( lmask(i,j,k) )then
                        expected = (real(er(i,j,k) - of(i,j,k),dp)**2 + real(ef(i,j,k) - orw(i,j,k),dp)**2) / &
                            &(2.d0 * real(SIGMA,dp)**2)
                        maxerr = max(maxerr, abs(real(diff(i,j,k),dp) - expected))
                        maxexp = max(maxexp, expected)
                    else if( diff(i,j,k) /= 0. )then
                        n_outside = n_outside + 1
                    endif
                end do
            end do
        end do
        call assert_true(maxerr <= 1.d-5 * maxexp, 'the unary is the half sum of the squared cross-half residuals at the noise level')
        call assert_int(0, n_outside, 'outside the support the unary is the neutral zero')
        ! one residual of exactly sigma in each term costs one
        er  = 0.;  of  = 0.;  ef  = 0.;  orw = 0.
        er(N/2,N/2,N/2)  = SIGMA
        orw(N/2,N/2,N/2) = -SIGMA
        call e_raw%set_rmat(er, .false.)
        call e_filt%set_rmat(ef, .false.)
        call o_raw%set_rmat(orw, .false.)
        call o_filt%set_rmat(of, .false.)
        call e_raw%nu_objective(e_filt, o_raw, o_filt, diff, lmask, SIGMA)
        call assert_real(1., diff(N/2,N/2,N/2), 1.e-6, 'residuals of one noise level in both terms cost one')
        call e_raw%kill
        call e_filt%kill
        call o_raw%kill
        call o_filt%kill

    contains

        subroutine new_from( img, rmat )
            type(image), intent(inout) :: img
            real,        intent(in)    :: rmat(:,:,:)
            call img%new([N,N,N], SMPD, wthreads=.false.)
            call img%set_rmat(rmat, .false.)
        end subroutine new_from

    end subroutine test_objective_closed_form

    !> sigma_0 recovers the sdev of E - O; exact zero/zero voxels are left out; an all-zero pair gives one
    subroutine test_noise_scale()
        real, parameter :: SDEV = 0.3
        type(image) :: even, odd
        real,    allocatable :: white(:,:,:), s(:,:,:), e(:,:,:), o(:,:,:)
        logical, allocatable :: lmask(:,:,:)
        real    :: truth, est
        integer :: i, j, k
        write(logfhandle,'(A)') '-- test_noise_scale'
        call set_fixed_seed(SEED + 2)
        allocate(white(BOX,BOX,BOX), e(BOX,BOX,BOX), o(BOX,BOX,BOX))
        call fill_gaussian(white)
        s = smoothed_field(white, SIG_FINE)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    e(i,j,k) = s(i,j,k) + SDEV * gasdev()
                    o(i,j,k) = s(i,j,k) + SDEV * gasdev()
                end do
            end do
        end do
        call support_sphere(lmask)
        call even%new([BOX,BOX,BOX], SMPD, wthreads=.false.)
        call odd%new( [BOX,BOX,BOX], SMPD, wthreads=.false.)
        call even%set_rmat(e, .false.)
        call odd%set_rmat(o,  .false.)
        truth = sqrt(2.) * SDEV
        est   = even%nu_objective_noise_scale(odd, lmask)
        write(logfhandle,'(A,F8.5,A,F8.5)') '   sigma_0 ', est, ', sdev of E - O ', truth
        call assert_real(truth, est, 0.03 * truth, 'sigma_0 is the sdev of the half-map difference')
        ! half of the support exactly zero in both halves, as a density-constrained solve leaves it
        e(BOX/2+1:,:,:) = 0.
        o(BOX/2+1:,:,:) = 0.
        call even%set_rmat(e, .false.)
        call odd%set_rmat(o,  .false.)
        est = even%nu_objective_noise_scale(odd, lmask)
        call assert_real(truth, est, 0.03 * truth, 'exact zero/zero voxels do not enter sigma_0')
        e = 0.
        o = 0.
        call even%set_rmat(e, .false.)
        call odd%set_rmat(o,  .false.)
        est = even%nu_objective_noise_scale(odd, lmask)
        call assert_real(1., est, 0., 'an all-zero pair has the unit level')
        call even%kill
        call odd%kill
    end subroutine test_noise_scale

    !> under stationary noise, scaling sigma_0 by a constant scales every unary by its inverse square
    !! and leaves the per-voxel winner among the candidates unchanged
    subroutine test_ranking_scale_free()
        real, parameter :: LPS(4)    = [12., 8., 6., 4.5]
        real, parameter :: SCALES(2) = [3.7, 0.2]
        type(image) :: even, odd
        type(image) :: even_bank(size(LPS)), odd_bank(size(LPS))
        real,    allocatable :: cost(:,:,:,:), cost_s(:,:,:,:)
        logical, allocatable :: lmask(:,:,:)
        real    :: sigma0, best, second, ref
        integer :: i, j, k, c, is, n_flips, n_off, n_resolved, nmask
        write(logfhandle,'(A)') '-- test_ranking_scale_free'
        call make_pair(.false., even, odd)
        call support_sphere(lmask)
        nmask  = count(lmask)
        sigma0 = even%nu_objective_noise_scale(odd, lmask)
        allocate(cost(BOX,BOX,BOX,size(LPS)), cost_s(BOX,BOX,BOX,size(LPS)))
        do c = 1, size(LPS)
            call even_bank(c)%copy(even)
            call odd_bank(c)%copy(odd)
            call even_bank(c)%bp(0., LPS(c))
            call odd_bank(c)%bp(0., LPS(c))
            call even%nu_objective(even_bank(c), odd, odd_bank(c), cost(:,:,:,c), lmask, sigma0)
        end do
        do is = 1, size(SCALES)
            do c = 1, size(LPS)
                call even%nu_objective(even_bank(c), odd, odd_bank(c), cost_s(:,:,:,c), lmask, SCALES(is) * sigma0)
            end do
            n_flips    = 0
            n_off      = 0
            n_resolved = 0
            do k = 1, BOX
                do j = 1, BOX
                    do i = 1, BOX
                        if( .not.lmask(i,j,k) ) cycle
                        do c = 1, size(LPS)
                            ref = cost(i,j,k,c) / SCALES(is)**2
                            if( abs(cost_s(i,j,k,c) - ref) > 1.e-5 * max(ref, TINY) ) n_off = n_off + 1
                        end do
                        best   = minval(cost(i,j,k,:))
                        second = minval(cost(i,j,k,:), mask=cost(i,j,k,:) > best)
                        ! a tie within single-precision rounding has no winner to preserve
                        if( second - best <= 1.e-4 * max(best, TINY) ) cycle
                        n_resolved = n_resolved + 1
                        if( minloc(cost_s(i,j,k,:), dim=1) /= minloc(cost(i,j,k,:), dim=1) ) n_flips = n_flips + 1
                    end do
                end do
            end do
            write(logfhandle,'(A,F5.2,A,I0,A,I0)') '   scale ', SCALES(is), ': resolved voxels ', n_resolved, &
                &' of ', nmask
            call assert_int(0, n_off, 'every unary scales by the inverse square of the noise level')
            call assert_true(n_resolved > nint(0.99 * real(nmask)), 'nearly every support voxel has a resolved winner')
            call assert_int(0, n_flips, 'the per-voxel winner does not depend on the noise level')
        end do
        do c = 1, size(LPS)
            call even_bank(c)%kill
            call odd_bank(c)%kill
        end do
        call even%kill
        call odd%kill
    end subroutine test_ranking_scale_free

    !> the full competition (bank, smoothing, like-for-like selection, Potts prior) on the split phantom:
    !! the fine half takes labels at or finer than 6 A, the coarse disagreeing half at or coarser than 10 A
    subroutine test_coarser_region()
        type(image) :: even, odd
        real,    allocatable :: lp_fine(:), lp_coarse(:)
        logical, allocatable :: region(:,:,:)
        real :: med_fine, med_coarse, frac_fine, frac_coarse
        write(logfhandle,'(A)') '-- test_coarser_region'
        call make_pair(.true., even, odd)
        call set_nu_filter_report(.false.)
        call setup_nu_dmats(even, odd, MSKDIAM, [real ::])
        call optimize_nu_cutoff_finds()
        allocate(region(BOX,BOX,BOX), source=.false.)
        region(1:BOX/2-MARGIN,:,:) = .true.
        call pack_filtmap_lowpass_limits(lp_fine, region)
        region = .false.
        region(BOX/2+1+MARGIN:BOX,:,:) = .true.
        call pack_filtmap_lowpass_limits(lp_coarse, region)
        call cleanup_nu_filter()
        call set_nu_filter_report(.true.)
        call even%kill
        call odd%kill
        call assert_true(size(lp_fine) > 1000 .and. size(lp_coarse) > 1000, 'both halves hold support voxels')
        frac_fine   = real(count(lp_fine   <= 6.5)) / real(max(1, size(lp_fine)))
        frac_coarse = real(count(lp_coarse >= 9.5)) / real(max(1, size(lp_coarse)))
        med_fine    = median_nocopy(lp_fine)
        med_coarse  = median_nocopy(lp_coarse)
        write(logfhandle,'(A,F6.2,A,F6.2,A)') '   median selected low-pass: fine half ', med_fine, &
            &' A, coarse half ', med_coarse, ' A'
        call assert_true(med_fine   <= 6.5, 'the fine half selects labels at or finer than 6 A')
        call assert_true(med_coarse >= 9.5, 'the coarse, disagreeing half selects labels at or coarser than 10 A')
        call assert_true(frac_fine   >= 0.9, 'at least 90% of the fine half is at or finer than 6 A')
        call assert_true(frac_coarse >= 0.9, 'at least 90% of the coarse half is at or coarser than 10 A')
    end subroutine test_coarser_region

    !> the neutral fixture through the evidence path: the compact state favours signal inside the molecule
    !! and the null in the solvent, the margin separates them, the envelope recovers the molecule; a pair
    !! zeroed outside a hard support keeps its noise level and freezes the unobserved voxels at the null
    subroutine test_evidence_envelope()
        type(image) :: even, odd, vol_supp
        type(nu_evidence_state)   :: evstate
        type(nu_evidence_summary) :: evsummary
        type(nu_envmask_params)   :: envp
        type(nu_envmask_stats)    :: envstats
        real,    allocatable :: rmat_even(:,:,:), rmat_odd(:,:,:), evmargin(:), ev_cutoff(:), ev_uncertainty(:)
        real,    allocatable :: ev_band_support(:,:)
        integer, allocatable :: ev_label(:)
        logical, allocatable :: l_true(:,:,:), l_supp(:,:,:), l_hard(:,:,:), l_env(:,:,:)
        integer(kind=8) :: s1, s2
        real    :: cen, dsq, sig, sigma_full, sigma_hard, support_in, support_out, null_in, null_out
        real    :: mean_in, mean_out, recall, fpr
        integer :: i, j, k, imask, n_in, n_out, n_unobserved, n_bad, iband
        write(logfhandle,'(A)') '-- test_evidence_envelope'
        cen = real(EBOX) / 2. + 1.
        s1  = 1234567_8
        s2  = 7654321_8
        allocate(rmat_even(EBOX,EBOX,EBOX), rmat_odd(EBOX,EBOX,EBOX), l_true(EBOX,EBOX,EBOX), &
            &l_hard(EBOX,EBOX,EBOX))
        do k = 1, EBOX
            do j = 1, EBOX
                do i = 1, EBOX
                    dsq = (real(i) - cen)**2 + (real(j) - cen)**2 + (real(k) - cen)**2
                    l_true(i,j,k) = dsq <= MOL_RAD_PX**2
                    l_hard(i,j,k) = dsq <= HARD_RAD_PX**2
                    sig = 0.
                    if( l_true(i,j,k) ) sig = 3.0 * sin(1.26 * real(i) + 0.63 * real(j)) + &
                        &2.0 * cos(0.94 * real(j) - 1.10 * real(k)) + 2.0 * sin(1.05 * real(k) + 0.42 * real(i))
                    rmat_even(i,j,k) = sig + ENV_NOISE * lcg_uniform(s1)
                    rmat_odd(i,j,k)  = sig + ENV_NOISE * lcg_uniform(s2)
                end do
            end do
        end do
        call even%new([EBOX,EBOX,EBOX], SMPD, wthreads=.false.)
        call odd%new( [EBOX,EBOX,EBOX], SMPD, wthreads=.false.)
        call even%set_rmat(rmat_even, .false.)
        call odd%set_rmat(rmat_odd,   .false.)
        call vol_supp%disc([EBOX,EBOX,EBOX], SMPD, SUPP_RAD_PX, l_supp)
        call vol_supp%kill
        call assert_true(all(l_supp .or. .not.l_true), 'the molecule lies inside the support')
        n_in  = count(l_supp .and. l_true)
        n_out = count(l_supp .and. .not.l_true)
        sigma_full = even%nu_objective_noise_scale(odd, l_supp)
        ! ---- the compact evidence state of the full pair ----
        call set_nu_filter_report(.false.)
        call setup_nu_dmats(even, odd, 2. * SUPP_RAD_PX * SMPD, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE)
        call optimize_nu_cutoff_finds()
        call build_nu_evidence_state(even, odd, evstate)
        call assert_true(nu_evidence_state_is_valid(evstate), 'the compact evidence state is valid')
        call get_nu_evidence_summary(evstate, evsummary)
        call assert_true(evsummary%null_fraction >= NU_EVIDENCE_MIN_NULL_FRAC .and. &
            &evsummary%null_fraction <= NU_EVIDENCE_MAX_NULL_FRAC, 'the explicit-null fraction is within its readiness bounds')
        call assert_int(9, evsummary%n_candidates, 'the null plus the eight-rung static bank')
        call assert_int(NU_EVIDENCE_NBANDS, evsummary%n_bands, 'the static evidence bands')
        call assert_int(count(l_supp), evsummary%n_support, 'the state covers the spherical support')
        call assert_true(index(evsummary%provenance, 'algorithm=nu_evidence_v2') > 0, 'the provenance names the squared-error algorithm')
        call assert_true(index(evsummary%provenance, ';noise_scale=') > 0, 'the provenance carries the noise level')
        call assert_real(1., evsummary%observed_fraction, 0., 'a pair without exact zeros is fully observed')
        call unpack_nu_evidence_state(evstate, ev_label, ev_cutoff, ev_uncertainty, ev_band_support)
        support_in  = 0.
        support_out = 0.
        null_in     = 0.
        null_out    = 0.
        imask = 0
        do k = 1, EBOX
            do j = 1, EBOX
                do i = 1, EBOX
                    if( .not.l_supp(i,j,k) ) cycle
                    imask = imask + 1
                    if( l_true(i,j,k) )then
                        support_in = support_in + ev_band_support(imask,1)
                        if( ev_label(imask) == 0 ) null_in = null_in + 1.
                    else
                        support_out = support_out + ev_band_support(imask,1)
                        if( ev_label(imask) == 0 ) null_out = null_out + 1.
                    endif
                end do
            end do
        end do
        support_in  = support_in  / real(n_in)
        support_out = support_out / real(n_out)
        null_in     = null_in     / real(n_in)
        null_out    = null_out    / real(n_out)
        write(logfhandle,'(A,F6.3,A,F6.3,A,F6.3,A,F6.3)') '   coarse-band support in/out ', support_in, '/', &
            &support_out, ', null selections in/out ', null_in, '/', null_out
        call assert_true(support_in > support_out, 'coarse-band support is higher inside the molecule')
        call assert_true(null_out > null_in, 'the explicit null is selected preferentially in the solvent')
        ! ---- the margin and the envelope ----
        call calc_nu_evidence_margin(evmargin, ENV_LP, .false.)
        call assert_int(count(l_supp), size(evmargin), 'the margin is packed over the support')
        call assert_true(all(evmargin >= 0.), 'the margin is non-negative by construction')
        mean_in  = 0.
        mean_out = 0.
        imask    = 0
        do k = 1, EBOX
            do j = 1, EBOX
                do i = 1, EBOX
                    if( .not.l_supp(i,j,k) ) cycle
                    imask = imask + 1
                    if( l_true(i,j,k) )then
                        mean_in = mean_in + evmargin(imask)
                    else
                        mean_out = mean_out + evmargin(imask)
                    endif
                end do
            end do
        end do
        mean_in  = mean_in  / real(n_in)
        mean_out = mean_out / real(n_out)
        write(logfhandle,'(A,ES10.3,A,ES10.3)') '   mean evidence margin in/out ', mean_in, '/', mean_out
        call assert_true(mean_in > 2. * mean_out, 'the margin separates the molecule from the solvent')
        envp%nsigma      = 3.0
        envp%beta        = 1.0
        envp%dens_weight = 0.0
        envp%lp_smooth   = ENV_LP
        envp%l_relative  = .false.
        envp%maxits      = 6
        call nu_evidence_envelope(envp, l_env, envstats)
        call assert_true(envstats%n_signal > 0, 'the envelope is not empty')
        call assert_true(envstats%l_null_valid, 'the solvent holds the majority of the support, so the null is valid')
        call assert_false(any(l_env .and. .not.l_supp), 'the envelope stays inside the support')
        recall = real(count(l_env .and. l_true)) / real(count(l_true))
        fpr    = real(count(l_env .and. l_supp .and. .not.l_true)) / real(n_out)
        write(logfhandle,'(A,F6.3,A,F6.3)') '   envelope recall of true density ', recall, &
            &', solvent false-positive rate ', fpr
        call assert_true(recall >= 0.90, 'the envelope recovers at least 90% of the molecule')
        call assert_true(fpr <= 0.15,    'the envelope admits at most 15% of the solvent')
        call cleanup_nu_filter()
        deallocate(ev_label, ev_cutoff, ev_uncertainty, ev_band_support, evmargin, l_env)
        ! ---- the pair zeroed outside a hard support ----
        where( .not.l_hard )
            rmat_even = 0.
            rmat_odd  = 0.
        end where
        call even%set_rmat(rmat_even, .false.)
        call odd%set_rmat(rmat_odd,   .false.)
        sigma_hard = even%nu_objective_noise_scale(odd, l_supp)
        write(logfhandle,'(A,F8.5,A,F8.5)') '   sigma_0 of the full pair ', sigma_full, ', of the hard-supported pair ', sigma_hard
        call assert_real(sigma_full, sigma_hard, 0.05 * sigma_full, 'exact zero/zero voxels leave the noise level unchanged')
        call setup_nu_dmats(even, odd, 2. * SUPP_RAD_PX * SMPD, [real ::], evidence_source=NU_EVIDENCE_SOURCE_BASE)
        call optimize_nu_cutoff_finds()
        call build_nu_evidence_state(even, odd, evstate)
        call assert_true(nu_evidence_state_is_valid(evstate), 'the hard-supported pair gives a valid evidence state')
        call get_nu_evidence_summary(evstate, evsummary)
        call assert_true(evsummary%null_fraction >= NU_EVIDENCE_MIN_NULL_FRAC .and. &
            &evsummary%null_fraction <= NU_EVIDENCE_MAX_NULL_FRAC, 'its explicit-null fraction is within the readiness bounds')
        call assert_real(real(count(l_supp .and. l_hard)) / real(count(l_supp)), evsummary%observed_fraction, 1.e-6, &
            &'the observed fraction is the hard support inside the sphere')
        call unpack_nu_evidence_state(evstate, ev_label, ev_cutoff, ev_uncertainty, ev_band_support)
        n_unobserved = 0
        n_bad        = 0
        imask        = 0
        do k = 1, EBOX
            do j = 1, EBOX
                do i = 1, EBOX
                    if( .not.l_supp(i,j,k) ) cycle
                    imask = imask + 1
                    if( l_hard(i,j,k) ) cycle
                    n_unobserved = n_unobserved + 1
                    if( ev_label(imask) /= 0 ) n_bad = n_bad + 1
                    if( abs(ev_uncertainty(imask) - 1.) > 1.e-6 ) n_bad = n_bad + 1
                    do iband = 1, size(ev_band_support,2)
                        if( abs(ev_band_support(imask,iband)) > 1.e-6 ) n_bad = n_bad + 1
                    end do
                end do
            end do
        end do
        call assert_int(count(l_supp .and. .not.l_hard), n_unobserved, 'every voxel outside the hard support is unobserved')
        call assert_int(0, n_bad, 'unobserved voxels are frozen at the null, maximally uncertain, with no band support')
        call cleanup_nu_filter()
        call set_nu_filter_report(.true.)
        call even%kill
        call odd%kill
        deallocate(rmat_even, rmat_odd, l_true, l_supp, l_hard, ev_label, ev_cutoff, ev_uncertainty, ev_band_support)

    contains

        !> deterministic uniform deviate on [-0.5,0.5]; two streams give the halves independent noise
        real function lcg_uniform( lcg_state )
            integer(kind=8), intent(inout) :: lcg_state
            lcg_state = mod(1103515245_8 * lcg_state + 12345_8, 2147483648_8)
            lcg_uniform = real(lcg_state) / 2147483648.0 - 0.5
        end function lcg_uniform

    end subroutine test_evidence_envelope

end module simple_nu_filter_tester
