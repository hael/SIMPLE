!@descr: order prior of the PCG base solve (pcg_order=yes): a ridge on density without NU evidence of order
! Penalty lambda h (1 - w_ord) x^2, w_ord a logistic on the NU-evidence margin at ordered fraction f, h the data
! precision per voxel. f and lambda in closed form (octave bands, each shrunk per voxel by D_b/(D_b + lambda h
! (1 - w_ord))), keeping 95 % of the protein core's cross-half signal: the most removed from the dense density around
! it, or with none around it the strongest lambda. Design: doc/implementation_notes/planned/pcg_order_prior.md.
module simple_pcg_order_prior
use, intrinsic :: iso_c_binding, only: c_float
use simple_defs,              only: dp, TINY, logfhandle
use simple_defs_fname,        only: MRC_EXT, NUEVIDENCE_SUFFIX, NULOCRES_SUFFIX, PCG_ORDER_WEIGHT_FBODY, &
    &PCG_ORDER_REGIONS_FBODY
use simple_error,             only: simple_exception
use simple_math,              only: otsu
use simple_math_ft,           only: fdim, get_resarr
use simple_estimate_ssnr,     only: get_resolution
use simple_srch_sort_loc,     only: hpsort
use simple_string,            only: string
use simple_string_utils,      only: int2str_pad
use simple_fileio,            only: add2fbody
use simple_syslib,            only: file_exists
use simple_imghead,           only: find_ldim_nptcls
use simple_image,             only: image
use simple_parameters,        only: parameters
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_refine3D_fnames,   only: refine3D_state_halfvol_fname, refine3D_state_vol_fname
implicit none

public :: build_density_weight, order_regions, build_order_weights, select_order_prior, order_prior_stats
public :: choose_order_setting, ORDER_FRAC_GRID, ORDER_LAMBDA_GRID, CORE_SIGNAL_MIN, CORE_MIN_VOX
public :: prepare_order_prior_on_pair, report_order_prior, write_pcg_order_pair, order_prior_provenance
public :: run_all_pcg_order_prior_tests
private
#include "simple_local_flags.inc"

!> ordered fractions of the support scanned by the evidence weight
real,    parameter :: ORDER_FRAC_GRID(8)   = [0.01, 0.02, 0.04, 0.08, 0.16, 0.32, 0.48, 0.64]
!> ridge strengths relative to the data precision per voxel, ascending: a fully penalized voxel shrinks by 1/(1+lambda)
real,    parameter :: ORDER_LAMBDA_GRID(8) = [0.1, 0.2, 0.5, 1., 2., 5., 10., 20.]
real,    parameter :: CORE_RES_FAC    = 1.5   !< protein core: dense, NU label at or finer than this multiple of FSC=0.143
real,    parameter :: CORE_SIGNAL_MIN = 0.95  !< the core keeps at least this fraction of its cross-half signal
integer, parameter :: CORE_MIN_VOX    = 1000  !< fewer core (or surrounding) voxels than this: nothing to check
!> smoothing scale of the density statistic, a multiple of the pair's FSC=0.143 (Wang: about twice), floored in A
real,    parameter :: DENSITY_LP_FACTOR = 2.0
real,    parameter :: DENSITY_LP_MIN_A  = 8.0
real,    parameter :: SUPPORT_MIN       = 0.5  !< production-support value taken as the evaluation domain

type :: order_prior_stats
    real    :: frac        = 0.    !< chosen ordered fraction
    real    :: lambda_rel  = 0.
    real    :: removed     = 0.    !< fraction of the dense non-core cross-half signal removed
    real    :: core_kept   = 1.    !< fraction of the core's cross-half signal kept
    real    :: pct_core    = 0.    !< core, % of the evaluation domain
    real    :: pct_peri    = 0.    !< dense outside the core, % of the evaluation domain
    logical :: l_off       = .true. !< no setting keeps the core: no prior
    ! diagnostics (report_order_prior): what a guard on every band to FSC=0.143 would choose
    real    :: guard_frac = 0., guard_lambda = 0., guard_removed = 0., guard_core_kept = 1.
    logical :: l_guard_off = .true.
    real    :: thresh = 0., width = 0.    !< evidence threshold and logistic width at the chosen fraction
    real    :: pen(4,3) = 0.              !< 1 - w_ord: mean, > 0.05, > 0.2, > 0.5 (rows) in core, peri, solvent
    real    :: smpd = 0., res0143 = 0.
    integer :: shell0143 = 0
    integer, allocatable :: band_lo(:), band_hi(:), idx(:,:)  !< octave bands; the dense domain voxels
    logical, allocatable :: pcore(:)                          !< dense voxel in the core
    real,    allocatable :: dband_rel(:)  !< band data diagonal over the per-voxel precision h, mean of the halves
    real,    allocatable :: share(:)      !< band's share of the core's base cross-half signal
    real,    allocatable :: kept_b(:), rem_b(:), guard_kept_b(:) !< predicted per band: chosen; band guard's choice
end type order_prior_stats

interface
    !> unit tests (submodule simple_pcg_order_prior_tester, built only with BUILD_TESTS=ON)
    module subroutine run_all_pcg_order_prior_tests()
    end subroutine run_all_pcg_order_prior_tests
end interface

contains

    !> Density weight w(r) in [0,1] of ONE prior-free half (gold standard kept): |x| smoothed at
    !! max(8 A, 2 x res0143), Otsu threshold over support >= 1/2, logistic of the solvent class's spread
    subroutine build_density_weight( base_half, res0143, support, weight, solvent_frac )
        class(image),         intent(in)    :: base_half
        real,                 intent(in)    :: res0143
        class(image), target, intent(in)    :: support
        type(image),          intent(inout) :: weight
        real,                 intent(out)   :: solvent_frac
        real(c_float), pointer :: rmat_s(:,:,:), rmat_sup(:,:,:)
        real,    allocatable :: vals(:)
        integer :: ldim(3), i, j, k, n, n_solv, n_below
        real    :: smpd, lp, thresh, width, arg
        real(dp) :: s_solv, ss_solv
        ldim = base_half%get_ldim()
        smpd = base_half%get_smpd()
        if( any(support%get_ldim() /= ldim) ) THROW_HARD('support and base half on different grids; build_density_weight')
        if( base_half%is_ft() .or. support%is_ft() ) THROW_HARD('inputs must be in real space; build_density_weight')
        if( res0143 <= 0. ) THROW_HARD('requires a positive base-pair resolution; build_density_weight')
        lp = max(DENSITY_LP_MIN_A, DENSITY_LP_FACTOR * res0143)
        lp = min(lp, real(ldim(1)) * smpd / 4.)
        call weight%copy(base_half)
        call weight%get_rmat_ptr(rmat_s)
        rmat_s(:ldim(1),:ldim(2),:ldim(3)) = abs(rmat_s(:ldim(1),:ldim(2),:ldim(3)))
        nullify(rmat_s)
        call weight%bp(0., lp)
        call weight%get_rmat_ptr(rmat_s)
        call support%get_rmat_ptr(rmat_sup)
        n = count(rmat_sup(:ldim(1),:ldim(2),:ldim(3)) >= SUPPORT_MIN)
        if( n < 1 ) THROW_HARD('the production support is empty; build_density_weight')
        vals = pack(rmat_s(:ldim(1),:ldim(2),:ldim(3)), rmat_sup(:ldim(1),:ldim(2),:ldim(3)) >= SUPPORT_MIN)
        thresh = 0.
        call otsu(n, vals, thresh)
        n_solv = count(vals < thresh)
        if( n_solv < 2 .or. n_solv >= n ) THROW_HARD('Otsu partition degenerate inside the support; build_density_weight')
        s_solv  = sum(real(vals,dp), mask=vals < thresh) / real(n_solv,dp)
        ss_solv = sum(real(vals,dp)**2, mask=vals < thresh) / real(n_solv,dp)
        width   = real(sqrt(max(ss_solv - s_solv**2, 0.0_dp)))
        width   = max(width, 1.e-3 * max(abs(thresh), TINY))
        deallocate(vals)
        n_below = 0
        !$omp parallel do collapse(3) default(shared) private(i,j,k,arg) reduction(+:n_below) proc_bind(close) schedule(static)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    arg = min(40., max(-40., (rmat_s(i,j,k) - thresh) / width))
                    rmat_s(i,j,k) = 1. / (1. + exp(-arg))
                    if( rmat_sup(i,j,k) >= SUPPORT_MIN .and. rmat_s(i,j,k) < 0.5 ) n_below = n_below + 1
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rmat_s, rmat_sup)
        solvent_frac = real(n_below) / real(n)
    end subroutine build_density_weight

    !> Ordered fraction and strength in closed form, and the order weight w_ord they give (both halves).
    !! base, w_dens: prior-free halves and density weights (regions only); evmap, locres: NU evidence margin and local
    !! resolution (A); lsup: evaluation domain; dshell(0:,2), hdiag(2): per half, the operator's diagonal per shell and
    !! its mean (the data precision per voxel)
    subroutine select_order_prior( base, w_dens, evmap, locres, res0143, lsup, dshell, hdiag, lambda_fixed, &
            &w_order, st, l_tables )
        class(image),            intent(in)    :: base(2), w_dens(2), evmap, locres
        real,                    intent(in)    :: res0143
        logical,                 intent(in)    :: lsup(:,:,:)
        real,                    intent(in)    :: dshell(0:,:), hdiag(2), lambda_fixed
        type(image),             intent(inout) :: w_order(2)
        type(order_prior_stats), intent(out)   :: st
        logical,                 intent(in)    :: l_tables
        integer, parameter :: NF = size(ORDER_FRAC_GRID), NLMAX = size(ORDER_LAMBDA_GRID), NBMAX = 16
        real(c_float), pointer :: rw1(:,:,:), rw2(:,:,:), rev(:,:,:), rlr(:,:,:)
        real,     allocatable  :: vals(:), lambdas(:), xb(:,:,:), dband(:,:), pev(:), kept(:,:), rem(:,:), kept_g(:,:)
        integer,  allocatable  :: band_lo(:), band_hi(:), idx(:,:)
        logical,  allocatable  :: pcore(:)
        real(dp) :: core_kept(NF,NLMAX), peri_kept(NF,NLMAX), core_base, peri_base, p
        real(dp) :: core_kept_b(NBMAX,NF,NLMAX), peri_kept_b(NBMAX,NF,NLMAX), core_base_b(NBMAX), peri_base_b(NBMAX)
        real     :: thresh(NF), width(NF), lim, wo, lam1, lam2, ye, yo, se, so, pen
        logical  :: l_frac(NF), l_valid, l_guard(NBMAX)
        integer  :: ldim(3), n, nl, nb, nd, i, j, k, iv, ib, ifr, il, ncore, nperi, ndom, fsel, lsel, gsel, lgsel
        ldim = base(1)%get_ldim()
        n    = ldim(1)
        if( lambda_fixed > 0. )then
            lambdas = [lambda_fixed]
        else
            lambdas = ORDER_LAMBDA_GRID
        endif
        nl = size(lambdas)
        call w_dens(1)%get_rmat_ptr(rw1)
        call w_dens(2)%get_rmat_ptr(rw2)
        call evmap%get_rmat_ptr(rev)
        call locres%get_rmat_ptr(rlr)
        ! evidence thresholds at the (1-f) quantiles over the domain
        vals = pack(rev(:n,:n,:n), lsup)
        call hpsort(vals)
        call evidence_thresholds(vals, ORDER_FRAC_GRID, thresh, width, l_frac)
        ndom = size(vals)
        deallocate(vals)
        ! the dense domain voxels (core and peri), packed with their weights and evidence
        lim = CORE_RES_FAC * res0143
        nd  = count(lsup .and. 0.5 * (rw1(:n,:n,:n) + rw2(:n,:n,:n)) >= 0.5)
        allocate(idx(3,nd), pev(nd), pcore(nd))
        iv = 0
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if( .not. lsup(i,j,k) ) cycle
                    if( 0.5 * (rw1(i,j,k) + rw2(i,j,k)) < 0.5 ) cycle
                    iv = iv + 1
                    idx(:,iv) = [i,j,k]
                    pev(iv)   = rev(i,j,k)
                    pcore(iv) = rlr(i,j,k) > 0. .and. rlr(i,j,k) <= lim
                end do
            end do
        end do
        ncore = count(pcore)
        nperi = nd - ncore
        ! octave bands of Fourier shells, each shrunk by its own share of the operator's diagonal
        call octave_bands(n, band_lo, band_hi)
        nb = size(band_lo)
        if( nb > NBMAX ) THROW_HARD('more octave bands than NBMAX; select_order_prior')
        allocate(dband(nb,2), xb(nb,nd,2))
        do ib = 1, nb
            dband(ib,1) = band_diagonal(dshell(:,1), band_lo(ib), band_hi(ib))
            dband(ib,2) = band_diagonal(dshell(:,2), band_lo(ib), band_hi(ib))
        end do
        call split_bands(base(1), band_lo, band_hi, idx, xb(:,:,1))
        call split_bands(base(2), band_lo, band_hi, idx, xb(:,:,2))
        ! cross-half signal x_e x_o (noise cancels) of the band-shrunk halves, in the core and outside it
        core_kept   = 0._dp
        peri_kept   = 0._dp
        core_base   = 0._dp
        peri_base   = 0._dp
        core_kept_b = 0._dp
        peri_kept_b = 0._dp
        core_base_b = 0._dp
        peri_base_b = 0._dp
        ! the per-band terms (same band in both halves) are diagnostics: the rule reads the totals
        !$omp parallel do default(shared) private(iv,ifr,il,ib,p,wo,pen,lam1,lam2,ye,yo,se,so) &
        !$omp reduction(+:core_kept,peri_kept,core_base,peri_base,core_kept_b,peri_kept_b,core_base_b,peri_base_b) &
        !$omp schedule(static) proc_bind(close)
        do iv = 1, nd
            p = real(sum(xb(:,iv,1)),dp) * real(sum(xb(:,iv,2)),dp)
            if( pcore(iv) )then
                core_base = core_base + p
                core_base_b(:nb) = core_base_b(:nb) + real(xb(:,iv,1),dp) * real(xb(:,iv,2),dp)
            else
                peri_base = peri_base + p
                peri_base_b(:nb) = peri_base_b(:nb) + real(xb(:,iv,1),dp) * real(xb(:,iv,2),dp)
            endif
            do ifr = 1, NF
                if( .not. l_frac(ifr) ) cycle
                wo  = 1. / (1. + exp(-min(40., max(-40., (pev(iv) - thresh(ifr)) / width(ifr)))))
                pen = 1. - wo
                do il = 1, nl
                    lam1 = lambdas(il) * hdiag(1) * pen
                    lam2 = lambdas(il) * hdiag(2) * pen
                    ye = 0.
                    yo = 0.
                    do ib = 1, nb
                        se = xb(ib,iv,1) * shrink(dband(ib,1), lam1)
                        so = xb(ib,iv,2) * shrink(dband(ib,2), lam2)
                        ye = ye + se
                        yo = yo + so
                        if( pcore(iv) )then
                            core_kept_b(ib,ifr,il) = core_kept_b(ib,ifr,il) + real(se,dp) * real(so,dp)
                        else
                            peri_kept_b(ib,ifr,il) = peri_kept_b(ib,ifr,il) + real(se,dp) * real(so,dp)
                        endif
                    end do
                    p = real(ye,dp) * real(yo,dp)
                    if( pcore(iv) )then
                        core_kept(ifr,il) = core_kept(ifr,il) + p
                    else
                        peri_kept(ifr,il) = peri_kept(ifr,il) + p
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(xb)
        st%pct_core = 100. * real(ncore) / real(max(1,ndom))
        st%pct_peri = 100. * real(nperi) / real(max(1,ndom))
        ! fractions of the core's signal kept and of the surrounding signal removed, and the rule on them
        allocate(kept(NF,nl), rem(NF,nl), source=0.)
        if( core_base > 0._dp ) kept = real(core_kept(:,:nl) / core_base)
        if( peri_base > 0._dp ) rem  = real(1._dp - peri_kept(:,:nl) / peri_base)
        call choose_order_setting(kept, rem, l_frac, ncore >= CORE_MIN_VOX .and. core_base > 0._dp, &
            &nperi >= CORE_MIN_VOX .and. peri_base > 0._dp, lambda_fixed > 0., fsel, lsel)
        st%l_off      = fsel == 0
        st%frac       = 0.
        st%lambda_rel = 0.
        st%core_kept  = 1.
        st%removed    = 0.
        if( .not. st%l_off )then
            st%frac       = ORDER_FRAC_GRID(fsel)
            st%lambda_rel = lambdas(lsel)
            st%core_kept  = kept(fsel,lsel)
            st%removed    = rem(fsel,lsel)
        endif
        call band_diagnostics
        if( l_tables ) call print_tables
        ! the chosen weight, the same for both halves (none: no prior)
        if( .not. st%l_off ) call build_order_weights(evmap, lsup, st%frac, w_order, l_valid)
        nullify(rw1, rw2, rev, rlr)
        call move_alloc(idx, st%idx)
        call move_alloc(pcore, st%pcore)
        call move_alloc(band_lo, st%band_lo)
        call move_alloc(band_hi, st%band_hi)
        deallocate(lambdas, dband, kept, rem, kept_g, pev)

    contains

        !> per-band kept and removed at the chosen setting, the choice of a guard on every band that starts below
        !! the FSC=0.143 shell, and the spread of the penalty over the regions
        subroutine band_diagnostics
            real, parameter :: PEN_LEVELS(3) = [0.05, 0.2, 0.5]
            real(dp) :: psum(3)
            real     :: pen_v
            integer  :: ireg, cnt(3), nabv(3,3)
            st%smpd      = base(1)%get_smpd()
            st%res0143   = res0143
            st%shell0143 = nint(real(n) * st%smpd / res0143)
            allocate(st%dband_rel(nb), st%share(nb), st%kept_b(nb), st%rem_b(nb), st%guard_kept_b(nb), source=0.)
            do ib = 1, nb
                st%dband_rel(ib) = 0.5 * (dband(ib,1) / max(hdiag(1),TINY) + dband(ib,2) / max(hdiag(2),TINY))
            end do
            if( core_base > 0._dp ) st%share = real(core_base_b(:nb) / core_base)
            l_guard = .false.
            l_guard(:nb) = band_lo <= st%shell0143 .and. core_base_b(:nb) > 0._dp
            ! the band guard: the worst guarded band in place of the total, the same rule
            allocate(kept_g(NF,nl))
            kept_g = kept
            do ib = 1, nb
                if( l_guard(ib) ) kept_g = min(kept_g, real(core_kept_b(ib,:,:nl) / core_base_b(ib)))
            end do
            call choose_order_setting(kept_g, rem, l_frac, ncore >= CORE_MIN_VOX .and. core_base > 0._dp, &
                &nperi >= CORE_MIN_VOX .and. peri_base > 0._dp, lambda_fixed > 0., gsel, lgsel)
            st%l_guard_off = gsel == 0
            if( .not. st%l_guard_off )then
                st%guard_frac      = ORDER_FRAC_GRID(gsel)
                st%guard_lambda    = lambdas(lgsel)
                st%guard_removed   = rem(gsel,lgsel)
                st%guard_core_kept = kept(gsel,lgsel)
                call per_band(gsel, lgsel, st%guard_kept_b)
            endif
            if( st%l_off ) return
            call per_band(fsel, lsel, st%kept_b)
            do ib = 1, nb
                if( peri_base_b(ib) > 0._dp ) st%rem_b(ib) = real(1._dp - peri_kept_b(ib,fsel,lsel) / peri_base_b(ib))
            end do
            ! 1 - w_ord at the chosen fraction over core, dense around it, and the rest of the domain
            st%thresh = thresh(fsel)
            st%width  = width(fsel)
            cnt  = 0
            psum = 0._dp
            nabv = 0
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        if( .not. lsup(i,j,k) ) cycle
                        pen_v = 1. - 1. / (1. + exp(-min(40., max(-40., (rev(i,j,k) - st%thresh) / st%width))))
                        ireg = 3
                        if( 0.5 * (rw1(i,j,k) + rw2(i,j,k)) >= 0.5 )then
                            ireg = 2
                            if( rlr(i,j,k) > 0. .and. rlr(i,j,k) <= lim ) ireg = 1
                        endif
                        cnt(ireg)  = cnt(ireg) + 1
                        psum(ireg) = psum(ireg) + real(pen_v,dp)
                        where( pen_v > PEN_LEVELS ) nabv(:,ireg) = nabv(:,ireg) + 1
                    end do
                end do
            end do
            do ireg = 1, 3
                st%pen(1,ireg)   = real(psum(ireg) / real(max(1,cnt(ireg)),dp))
                st%pen(2:4,ireg) = 100. * real(nabv(:,ireg)) / real(max(1,cnt(ireg)))
            end do
        end subroutine band_diagnostics

        subroutine per_band( jf, jl, kb )
            integer, intent(in)  :: jf, jl
            real,    intent(out) :: kb(:)
            kb = 0.
            do ib = 1, nb
                if( core_base_b(ib) > 0._dp ) kb(ib) = real(core_kept_b(ib,jf,jl) / core_base_b(ib))
            end do
        end subroutine per_band

        subroutine print_tables
            integer :: jf, jl
            write(logfhandle,'(A)') '>>> ORDER PRIOR: closed form, % of the dense cross-half signal around the core removed,'//&
                &' rows ordered fraction, columns lambda_rel; * keeps 95 % of the core, < chosen'
            write(logfhandle,'(A10,*(F9.2))') 'frac', lambdas
            do jf = 1, NF
                if( .not. l_frac(jf) ) cycle
                write(logfhandle,'(F10.3,*(F8.1,A1))') ORDER_FRAC_GRID(jf), (100. * rem(jf,jl), mark(jf,jl), jl=1,nl)
            end do
            write(logfhandle,'(A)') '>>> ORDER PRIOR: closed form, % of the core cross-half signal kept'
            write(logfhandle,'(A10,*(F9.2))') 'frac', lambdas
            do jf = 1, NF
                if( .not. l_frac(jf) ) cycle
                write(logfhandle,'(F10.3,*(F8.1,A1))') ORDER_FRAC_GRID(jf), (100. * kept(jf,jl), mark(jf,jl), jl=1,nl)
            end do
        end subroutine print_tables

        character(len=1) function mark( jf, jl )
            integer, intent(in) :: jf, jl
            mark = ' '
            if( kept(jf,jl) >= CORE_SIGNAL_MIN ) mark = '*'
            if( jf == fsel .and. jl == lsel ) mark = '<'
        end function mark

    end subroutine select_order_prior

    !> The rule on the fractions of the core's signal kept, kept(f,lambda), and of the dense signal around the core
    !! removed, rem(f,lambda): among the settings that keep CORE_SIGNAL_MIN of the core, the most removal (l_peri), or
    !! with nothing around the core the strongest lambda at its smallest fraction; none (fsel = 0): no prior. A set
    !! strength (l_fixed, one column) is applied at the smallest fraction that keeps the core, else the largest
    subroutine choose_order_setting( kept, rem, l_frac, l_core, l_peri, l_fixed, fsel, lsel )
        real,    intent(in)  :: kept(:,:), rem(:,:)
        logical, intent(in)  :: l_frac(:), l_core, l_peri, l_fixed
        integer, intent(out) :: fsel, lsel
        real    :: best
        integer :: ifr, il
        fsel = 0
        lsel = 0
        if( l_fixed )then
            lsel = 1
            do ifr = 1, size(kept,1)
                if( .not. l_frac(ifr) ) cycle
                fsel = ifr
                if( l_core .and. kept(ifr,1) >= CORE_SIGNAL_MIN ) exit
            end do
            if( fsel == 0 ) lsel = 0
            return
        endif
        if( .not. l_core ) return
        best = -huge(best)
        do il = 1, size(kept,2)
            do ifr = 1, size(kept,1)
                if( .not. l_frac(ifr) ) cycle
                if( kept(ifr,il) < CORE_SIGNAL_MIN ) cycle
                if( l_peri )then
                    if( rem(ifr,il) <= best ) cycle
                    best = rem(ifr,il)
                else if( il == lsel )then
                    cycle
                endif
                fsel = ifr
                lsel = il
            end do
        end do
    end subroutine choose_order_setting

    !> a band's data diagonal: the mean of the shell diagonal over its shells, weighted by shell size; past the
    !! last shell given, the last
    real function band_diagonal( dsh, lo, hi ) result( d )
        real,    intent(in) :: dsh(0:)
        integer, intent(in) :: lo, hi
        real(dp) :: s, wsum, w
        integer  :: sh, nmax
        nmax = ubound(dsh,1)
        s    = 0._dp
        wsum = 0._dp
        do sh = lo, hi
            w    = max(1._dp, real(sh,dp)**2)
            s    = s + w * real(dsh(min(sh,nmax)),dp)
            wsum = wsum + w
        end do
        d = real(s / wsum)
    end function band_diagonal

    !> per-voxel shrink of one band under a ridge: d/(d + lam); a band without data stays as it is (it is zero)
    elemental real function shrink( d, lam )
        real, intent(in) :: d, lam
        shrink = 1.
        if( d + lam > 0. ) shrink = d / (d + lam)
    end function shrink

    !> Fourier shells 0..corner in octaves: [0,1], [2,3], [4,7], ...; the last reaches the box corner
    subroutine octave_bands( n, lo, hi )
        integer,              intent(in)    :: n
        integer, allocatable, intent(inout) :: lo(:), hi(:)
        integer :: nsh, nb, ib
        nsh = ceiling(sqrt(3.) * real(n/2)) + 1
        nb  = 1
        do while( 2**nb <= nsh )
            nb = nb + 1
        end do
        if( allocated(lo) ) deallocate(lo)
        if( allocated(hi) ) deallocate(hi)
        allocate(lo(nb), hi(nb))
        lo(1) = 0
        hi(1) = 1
        do ib = 2, nb
            lo(ib) = 2**(ib-1)
            hi(ib) = min(2**ib - 1, nsh)
        end do
        hi(nb) = nsh
    end subroutine octave_bands

    !> one half split into the octave bands at the packed voxels: xb(band, voxel); the bands sum to the half
    subroutine split_bands( half, lo, hi, idx, xb )
        class(image), intent(in)    :: half
        integer,      intent(in)    :: lo(:), hi(:), idx(:,:)
        real,         intent(inout) :: xb(:,:)
        type(image)            :: xft, work
        real(c_float), pointer :: rp(:,:,:)
        real,    allocatable   :: filt(:)
        real    :: xmean
        integer :: ldim(3), ib, iv, nsh
        ldim = half%get_ldim()
        nsh  = hi(size(hi))
        ! the mean is taken out first: apply_filter keeps the DC term at the filter's maximum, so every band would carry it
        call xft%copy(half)
        call xft%get_rmat_ptr(rp)
        xmean = real(sum(real(rp(:ldim(1),:ldim(2),:ldim(3)),dp)) / real(product(ldim),dp))
        rp(:ldim(1),:ldim(2),:ldim(3)) = rp(:ldim(1),:ldim(2),:ldim(3)) - xmean
        nullify(rp)
        call xft%fft
        allocate(filt(nsh))
        do ib = 1, size(lo)
            filt = 0.
            filt(max(1,lo(ib)):hi(ib)) = 1.
            call work%copy(xft)
            call work%apply_filter(filt)
            call work%ifft
            call work%get_rmat_ptr(rp)
            do iv = 1, size(idx,2)
                xb(ib,iv) = rp(idx(1,iv),idx(2,iv),idx(3,iv))
            end do
            nullify(rp)
            if( lo(ib) == 0 ) xb(ib,:) = xb(ib,:) + xmean
        end do
        call xft%kill
        call work%kill
        deallocate(filt)
    end subroutine split_bands

    !> the logistic order weight's threshold at the (1-f) quantile of the sorted evidence values and its width,
    !! the spread below it; a fraction of 1 orders everything; invalid where the threshold margin is zero
    subroutine evidence_thresholds( vals, fracs, thresh, width, l_valid )
        real,    intent(in)  :: vals(:), fracs(:)
        real,    intent(out) :: thresh(:), width(:)
        logical, intent(out) :: l_valid(:)
        real(dp) :: s, ss
        integer  :: ifr, it
        thresh  = 0.
        width   = 1.
        l_valid = .false.
        do ifr = 1, size(fracs)
            if( fracs(ifr) >= 1. )then
                l_valid(ifr) = .true.
                cycle
            endif
            it = min(size(vals), max(2, nint((1. - fracs(ifr)) * real(size(vals)))))
            thresh(ifr) = vals(it)
            if( thresh(ifr) <= TINY ) cycle
            l_valid(ifr) = .true.
            s  = sum(real(vals(1:it-1),dp)) / real(it-1,dp)
            ss = sum(real(vals(1:it-1),dp)**2) / real(it-1,dp)
            width(ifr) = max(real(sqrt(max(ss - s**2, 0._dp))), 1.e-3 * thresh(ifr))
        end do
    end subroutine evidence_thresholds

    !> protein core: domain voxels dense in both halves on average (mean density weight >= 1/2) with an NU label
    !! at or finer than CORE_RES_FAC x res0143; lperi: the dense domain voxels outside it
    subroutine order_regions( w_dens, locres, res0143, lsup, lcore, lperi )
        class(image),         intent(in)    :: w_dens(2), locres
        real,                 intent(in)    :: res0143
        logical,              intent(in)    :: lsup(:,:,:)
        logical, allocatable, intent(inout) :: lcore(:,:,:), lperi(:,:,:)
        real(c_float), pointer :: rw1(:,:,:), rw2(:,:,:), rlr(:,:,:)
        integer :: ldim(3), n, i, j, k
        real    :: lim
        ldim = locres%get_ldim()
        n    = ldim(1)
        lim  = CORE_RES_FAC * res0143
        call w_dens(1)%get_rmat_ptr(rw1)
        call w_dens(2)%get_rmat_ptr(rw2)
        call locres%get_rmat_ptr(rlr)
        if( allocated(lcore) ) deallocate(lcore)
        if( allocated(lperi) ) deallocate(lperi)
        allocate(lcore(n,n,n), lperi(n,n,n), source=.false.)
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    if( .not. lsup(i,j,k) ) cycle
                    if( 0.5 * (rw1(i,j,k) + rw2(i,j,k)) < 0.5 ) cycle
                    if( rlr(i,j,k) > 0. .and. rlr(i,j,k) <= lim )then
                        lcore(i,j,k) = .true.
                    else
                        lperi(i,j,k) = .true.
                    endif
                end do
            end do
        end do
        nullify(rw1, rw2, rlr)
    end subroutine order_regions

    !> the order weight w_ord at ordered fraction frac, a logistic on the evidence margin, for both halves;
    !! l_valid false (and w_ord = 1, no penalty) when the evidence is zero at that fraction
    subroutine build_order_weights( evmap, lsup, frac, w_order, l_valid )
        class(image), intent(in)    :: evmap
        logical,      intent(in)    :: lsup(:,:,:)
        real,         intent(in)    :: frac
        type(image),  intent(inout) :: w_order(2)
        logical,      intent(out)   :: l_valid
        real(c_float), pointer :: rev(:,:,:), rwo(:,:,:)
        real, allocatable :: vals(:)
        real     :: thresh(1), width(1)
        logical  :: lv(1)
        integer  :: ldim(3), n, i, j, k
        ldim = evmap%get_ldim()
        n    = ldim(1)
        call evmap%get_rmat_ptr(rev)
        vals = pack(rev(:n,:n,:n), lsup)
        call hpsort(vals)
        call evidence_thresholds(vals, [frac], thresh, width, lv)
        deallocate(vals)
        l_valid = lv(1)
        call w_order(1)%new(ldim, evmap%get_smpd())
        call w_order(1)%get_rmat_ptr(rwo)
        if( l_valid )then
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, n
                do j = 1, n
                    do i = 1, n
                        rwo(i,j,k) = 1. / (1. + exp(-min(40., max(-40., (rev(i,j,k) - thresh(1)) / width(1)))))
                    end do
                end do
            end do
            !$omp end parallel do
        else
            rwo(:n,:n,:n) = 1.
        endif
        nullify(rwo, rev)
        call w_order(2)%copy(w_order(1))
    end subroutine build_order_weights

    !> pcg_order=yes, between the prior-free base solve and its re-solve: the pair's FSC, the density weights (the
    !! regions), the previous iteration's NU evidence and local resolution of this state (lag one: the first
    !! iteration at a box has none and runs without the prior), the closed-form choice, the ridge on both operators;
    !! st keeps the choice's diagnostics for report_order_prior after the re-solve
    subroutine prepare_order_prior_on_pair( params, state_here, pcgop_even, pcgop_odd, x_even, x_odd, &
            &state_support, l_state_support, weight, l_weight, lambda_rel_out, st )
        class(parameters),        intent(in)    :: params
        integer,                  intent(in)    :: state_here
        class(reconstructor_pcg), intent(inout) :: pcgop_even, pcgop_odd
        real,                     intent(in)    :: x_even(:,:,:), x_odd(:,:,:)
        class(image),             intent(in)    :: state_support
        logical,                  intent(in)    :: l_state_support
        type(image),              intent(inout) :: weight(2)
        logical,                  intent(out)   :: l_weight
        real,                     intent(out)   :: lambda_rel_out
        type(order_prior_stats),  intent(out)   :: st
        type(image)             :: base(2), w_dens(2), support, evmap, locres
        type(string)            :: fname_vol, fname_ev, fname_lr, fbody
        real(c_float), pointer  :: rsup(:,:,:)
        logical, allocatable    :: lsup(:,:,:)
        real,    allocatable    :: corrs(:), res(:), ones(:,:,:), dsh(:,:), dsh_e(:), dsh_o(:)
        real    :: fsc05, res0143, solv(2), lambda_fixed
        integer :: n, ldim_ev(3), ldim_lr(3), nptcls
        l_weight       = .false.
        lambda_rel_out = 0.
        if( .not. params%l_pcg_order ) return
        n = params%box_crop
        fname_vol = refine3D_state_vol_fname(state_here)
        fname_ev  = add2fbody(fname_vol, MRC_EXT, NUEVIDENCE_SUFFIX)
        fname_lr  = add2fbody(fname_vol, MRC_EXT, NULOCRES_SUFFIX)
        if( .not. (file_exists(fname_ev) .and. file_exists(fname_lr)) )then
            write(logfhandle,'(A,I0,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
                &', no NU evidence of a previous iteration; prior off this iteration'
            call kill_locals
            return
        endif
        call find_ldim_nptcls(fname_ev, ldim_ev, nptcls)
        call find_ldim_nptcls(fname_lr, ldim_lr, nptcls)
        if( any(ldim_ev /= n) .or. any(ldim_lr /= n) )then
            write(logfhandle,'(A,I0,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
                &', the previous NU evidence is at another box; prior off this iteration'
            call kill_locals
            return
        endif
        call base(1)%new([n,n,n], params%smpd_crop)
        call base(2)%new([n,n,n], params%smpd_crop)
        call base(1)%set_rmat(x_even, .false.)
        call base(2)%set_rmat(x_odd,  .false.)
        ! the prior-free pair's FSC sets the density smoothing scale and the core's resolution limit
        allocate(corrs(fdim(n)-1), source=0.)
        call base(1)%fft
        call base(2)%fft
        call base(1)%fsc(base(2), corrs)
        call base(1)%ifft
        call base(2)%ifft
        res = get_resarr(n, params%smpd_crop)
        call get_resolution(corrs, res, fsc05, res0143)
        if( res0143 <= 0. )then
            write(logfhandle,'(A,I0,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, ', no base-pair resolution; prior off'
            call kill_locals
            return
        endif
        ! the production support: explicit pcg_mskfile, the state support, or the soft sphere
        if( params%pcg_mskfile%is_allocated() .and. len_trim(params%pcg_mskfile%to_char()) > 0 )then
            call support%read_and_crop(params%pcg_mskfile, params%smpd, params%box_crop, params%smpd_crop)
        else if( l_state_support )then
            call support%copy(state_support)
        else
            allocate(ones(n,n,n), source=1.0)
            call support%new([n,n,n], params%smpd_crop)
            call support%set_rmat(ones, .false.)
            call support%mask3D_soft(params%msk_crop, backgr=0.)
            deallocate(ones)
        endif
        call support%get_rmat_ptr(rsup)
        allocate(lsup(n,n,n))
        lsup = rsup(:n,:n,:n) >= SUPPORT_MIN
        nullify(rsup)
        call build_density_weight(base(1), res0143, support, w_dens(1), solv(1))
        call build_density_weight(base(2), res0143, support, w_dens(2), solv(2))
        call evmap%new([n,n,n], params%smpd_crop)
        call locres%new([n,n,n], params%smpd_crop)
        call evmap%read(fname_ev)
        call locres%read(fname_lr)
        lambda_fixed = -1.
        if( .not. params%l_pcg_order_lambda_auto ) lambda_fixed = params%pcg_order_lambda
        dsh_e = pcgop_even%get_shell_diagonal()
        dsh_o = pcgop_odd%get_shell_diagonal()
        allocate(dsh(0:size(dsh_e)-1,2))
        dsh(:,1) = dsh_e
        dsh(:,2) = dsh_o
        call select_order_prior(base, w_dens, evmap, locres, res0143, lsup, dsh, &
            &[pcgop_even%get_realspace_diagonal(), pcgop_odd%get_realspace_diagonal()], lambda_fixed, weight, st, .true.)
        deallocate(dsh, dsh_e, dsh_o)
        call write_order_regions(st, n, params%smpd_crop, string(PCG_ORDER_REGIONS_FBODY)//int2str_pad(state_here,2)//MRC_EXT)
        if( st%l_guard_off )then
            write(logfhandle,'(A,I0,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
                &', a guard on every band to FSC=0.143 would turn the prior off'
        else
            write(logfhandle,'(A,I0,A,F6.2,A,F7.2,A,F6.1,A,F6.1,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
                &', a guard on every band to FSC=0.143 would choose fraction ', 100.*st%guard_frac, ' %, lambda_rel ', &
                &st%guard_lambda, ': removes ', 100.*st%guard_removed, ' %, keeps ', 100.*st%guard_core_kept, ' % of the core'
        endif
        if( st%l_off )then
            write(logfhandle,'(A,I0,A,F6.1,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
                &', no setting keeps 95 % of the core (', st%pct_core, ' % of the support); prior off this iteration'
            call kill_locals
            return
        endif
        write(logfhandle,'(A,I0,A,F6.2,A,F7.2,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
            &', ordered fraction ', 100.*st%frac, ' %, lambda_rel ', st%lambda_rel, ': removes ', 100.*st%removed, &
            &' % of the dense signal around the core (', st%pct_peri, ' % of the support), keeps ', 100.*st%core_kept, &
            &' % of the core (', st%pct_core, ' %)'
        write(logfhandle,'(A,I0,A,ES10.3,A,ES10.3,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
            &', evidence threshold ', st%thresh, ', logistic width ', st%width, &
            &'; 1 - w_ord: mean, % of voxels > 0.05, > 0.2, > 0.5'
        write(logfhandle,'(A,F8.3,3F8.1)') '>>>   core            ', st%pen(1,1), st%pen(2:4,1)
        write(logfhandle,'(A,F8.3,3F8.1)') '>>>   dense around it ', st%pen(1,2), st%pen(2:4,2)
        write(logfhandle,'(A,F8.3,3F8.1)') '>>>   rest of domain  ', st%pen(1,3), st%pen(2:4,3)
        call pcgop_even%set_order_prior(weight(1), st%lambda_rel)
        call pcgop_odd%set_order_prior( weight(2), st%lambda_rel)
        lambda_rel_out = st%lambda_rel
        l_weight       = .true.
        ! pcg_order_weight_stateNN.mrc (one weight for both halves), overwritten every iteration
        fbody = string(PCG_ORDER_WEIGHT_FBODY)//int2str_pad(state_here,2)
        call weight(1)%write(fbody//MRC_EXT, del_if_exists=.true.)
        call fbody%kill
        call kill_locals

    contains

        subroutine kill_locals
            call base(1)%kill
            call base(2)%kill
            call w_dens(1)%kill
            call w_dens(2)%kill
            call support%kill
            call evmap%kill
            call locres%kill
            call fname_vol%kill
            call fname_ev%kill
            call fname_lr%kill
            if( allocated(lsup)  ) deallocate(lsup)
            if( allocated(corrs) ) deallocate(corrs)
            if( allocated(res)   ) deallocate(res)
        end subroutine kill_locals

    end subroutine prepare_order_prior_on_pair

    !> the regions of the choice as a map: 2 core, 1 dense around it, 0 elsewhere (inspection only)
    subroutine write_order_regions( st, n, smpd, fname )
        type(order_prior_stats), intent(in) :: st
        integer,                 intent(in) :: n
        real,                    intent(in) :: smpd
        type(string),            intent(in) :: fname
        type(image)            :: reg
        real(c_float), pointer :: rp(:,:,:)
        integer :: iv
        if( .not. allocated(st%idx) ) return
        call reg%new([n,n,n], smpd)
        call reg%get_rmat_ptr(rp)
        rp = 0.
        do iv = 1, size(st%idx,2)
            rp(st%idx(1,iv),st%idx(2,iv),st%idx(3,iv)) = merge(2., 1., st%pcore(iv))
        end do
        nullify(rp)
        call reg%write(fname, del_if_exists=.true.)
        call reg%kill
    end subroutine write_order_regions

    !> After the re-solve: per octave band, the core's cross-half signal kept and the dense signal around it removed,
    !! predicted by the closed form (st) and measured on the order pair against the base pair
    subroutine report_order_prior( st, state_here, base_even, base_odd, order_even, order_odd )
        type(order_prior_stats), intent(in) :: st
        integer,                 intent(in) :: state_here
        class(image),            intent(in) :: base_even, base_odd, order_even, order_odd
        real,     allocatable :: xb(:,:,:)
        real(dp), allocatable :: cb(:,:), pb(:,:)
        real(dp) :: ct(2), pt(2), pv
        real, parameter :: SHARE_MIN = 1.e-3
        real     :: rlo, rhi, kept_act, rem_act, sh_c, sh_p
        integer  :: nb, nd, ib, iv, ip, n, ldim(3)
        logical  :: l_c, l_p
        character(len=1) :: g
        if( st%l_off .or. .not. allocated(st%idx) ) return
        nb = size(st%band_lo)
        nd = size(st%idx,2)
        ldim = base_even%get_ldim()
        n    = ldim(1)
        allocate(xb(nb,nd,4), cb(nb,2), pb(nb,2))
        call split_bands(base_even,  st%band_lo, st%band_hi, st%idx, xb(:,:,1))
        call split_bands(base_odd,   st%band_lo, st%band_hi, st%idx, xb(:,:,2))
        call split_bands(order_even, st%band_lo, st%band_hi, st%idx, xb(:,:,3))
        call split_bands(order_odd,  st%band_lo, st%band_hi, st%idx, xb(:,:,4))
        ! (band, pair): cross-half signal of the base (1) and order (2) pair in the core (cb) and around it (pb)
        cb = 0._dp
        pb = 0._dp
        ct = 0._dp
        pt = 0._dp
        do iv = 1, nd
            do ip = 1, 2
                pv = real(sum(xb(:,iv,2*ip-1)),dp) * real(sum(xb(:,iv,2*ip)),dp)
                if( st%pcore(iv) )then
                    cb(:,ip) = cb(:,ip) + real(xb(:,iv,2*ip-1),dp) * real(xb(:,iv,2*ip),dp)
                    ct(ip)   = ct(ip) + pv
                else
                    pb(:,ip) = pb(:,ip) + real(xb(:,iv,2*ip-1),dp) * real(xb(:,iv,2*ip),dp)
                    pt(ip)   = pt(ip) + pv
                endif
            end do
        end do
        write(logfhandle,'(A,I0,A,F5.1,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
            &', per octave band (* starts below FSC=0.143, guarded by the band guard; core: NU label <= ', &
            &CORE_RES_FAC * st%res0143, ' A)'
        write(logfhandle,'(A)') '>>>   band  res (A)        share  around   D_b/h  ridge/D_b  core kept: pred  guard actual'//&
            &'  around removed: pred actual'
        do ib = 1, nb
            rhi = huge(rhi)
            if( st%band_lo(ib) > 0 ) rhi = real(n) * st%smpd / real(st%band_lo(ib))
            rlo = real(n) * st%smpd / real(max(1,st%band_hi(ib)))
            g = ' '
            if( st%band_lo(ib) <= st%shell0143 ) g = '*'
            ! a band with (almost) no base signal in a region has no meaningful ratio there
            sh_c = 0.
            sh_p = 0.
            if( ct(1) > 0._dp ) sh_c = real(cb(ib,1) / ct(1))
            if( pt(1) > 0._dp ) sh_p = real(pb(ib,1) / pt(1))
            l_c = sh_c > SHARE_MIN
            l_p = sh_p > SHARE_MIN
            kept_act = 0.
            rem_act  = 0.
            if( l_c ) kept_act = real(cb(ib,2) / cb(ib,1))
            if( l_p ) rem_act  = real(1._dp - pb(ib,2) / pb(ib,1))
            write(logfhandle,'(A,I4,A1,F7.1,A,F6.1,2F8.3,F8.3,ES11.2,3A7,4X,2A7)') '>>> ', ib, g, &
                &min(rhi, 9999.), '-', rlo, sh_c, sh_p, st%dband_rel(ib), &
                &st%lambda_rel / max(st%dband_rel(ib), TINY), pct(st%kept_b(ib), l_c), pct(st%guard_kept_b(ib), l_c), &
                &pct(kept_act, l_c), pct(st%rem_b(ib), l_p), pct(rem_act, l_p)
        end do
        kept_act = 0.
        rem_act  = 0.
        if( ct(1) > 0._dp ) kept_act = real(ct(2) / ct(1))
        if( pt(1) > 0._dp ) rem_act  = real(1._dp - pt(2) / pt(1))
        write(logfhandle,'(A,I0,A,F6.1,A,F6.1,A,F6.1,A,F6.1,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
            &', total core kept ', 100. * st%core_kept, ' % predicted, ', 100. * kept_act, ' % actual; around removed ', &
            &100. * st%removed, ' % predicted, ', 100. * rem_act, ' % actual'
        if( ct(1) > 0._dp .and. pt(1) > 0._dp ) write(logfhandle,'(A,F6.1,A,F6.1,A)') &
            &'>>>   same-band terms carry ', 100. * sum(cb(:,1)) / ct(1), ' % of the core signal and ', &
            &100. * sum(pb(:,1)) / pt(1), ' % of the signal around it (shares: of the totals; - : no signal)'
        deallocate(xb, cb, pb)

    contains

        character(len=7) function pct( x, l_ok )
            real,    intent(in) :: x
            logical, intent(in) :: l_ok
            pct = '      -'
            if( l_ok ) write(pct,'(F7.1)') 100. * x
        end function pct

    end subroutine report_order_prior

    !> pcg_order=yes: the ridge re-solve written beside the base pair as <half>_order, a diagnostic
    subroutine write_pcg_order_pair( params, state_here, order_even, order_odd )
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: state_here
        class(image),      intent(inout) :: order_even, order_odd
        type(string) :: fname
        fname = add2fbody(refine3D_state_halfvol_fname(state_here, 'even'), MRC_EXT, '_order')
        call order_even%write(fname, del_if_exists=.true.)
        fname = add2fbody(refine3D_state_halfvol_fname(state_here, 'odd'), MRC_EXT, '_order')
        call order_odd%write(fname, del_if_exists=.true.)
        if( params%part == 1 ) write(logfhandle,'(A,I0,A)') '>>> PCG ORDER PRIOR: STATE ', state_here, &
            &', prior-free pair written as _unfil (FSC, NU competition, evidence), prior''d pair as _order (references)'
        call fname%kill
    end subroutine write_pcg_order_pair

    !> the order_prior= line of the support provenance sidecar; lambda_eff is the ridge in force
    function order_prior_provenance( params, lambda_eff ) result( str )
        class(parameters), intent(in) :: params
        real,              intent(in) :: lambda_eff
        character(len=64) :: str
        if( params%l_pcg_order_lambda_auto )then
            write(str,'(A,F0.3,A)') 'shared_weight base_pair lambda_rel=', lambda_eff, ' auto'
        else
            write(str,'(A,F0.3,A)') 'shared_weight base_pair lambda_rel=', params%pcg_order_lambda, ' set'
        endif
    end function order_prior_provenance

end module simple_pcg_order_prior
