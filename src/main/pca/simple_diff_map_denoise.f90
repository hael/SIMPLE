!@descr: CSR graph coefficient-projection denoising helpers
module simple_diff_map_denoise
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, graph_matvec
use simple_image,           only: image
use simple_imgarr_utils,    only: copy_imgarr, dealloc_imgarr
use simple_linalg,          only: sparse_eigh
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc
implicit none
#include "simple_local_flags.inc"

private
public :: graph_coeffproj_denoise
public :: estimate_diffmap_denoise_rank
public :: graph_nystrom_residual_preimage
public :: calc_diffmap_reconstruction_error
public :: calc_diffmap_residual_energy_ratio
public :: select_spectral_rank_icm

integer, parameter :: DIFFMAP_DENOISE_ICM_RANK_MAXITS          = 16
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_BETA_FRAC       = 0.35
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_COMPLEXITY_FRAC = 0.10
real,    parameter :: DIFFMAP_DENOISE_ICM_RANK_LOWER_SEED_FRAC = 0.50

contains

    subroutine graph_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, rank_keep_override, verbose)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:)
        type(diffmap_graph), target, intent(in) :: graph
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        integer, optional,        intent(in)  :: rank_keep_override
        logical, optional,        intent(in)  :: verbose
        call so2_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, force_zero_theta=.true., &
                                   rank_keep_override=rank_keep_override, verbose=verbose)
    end subroutine graph_coeffproj_denoise

    subroutine so2_coeffproj_denoise(params, imgs, avg, graph, coeff_ptcls, force_zero_theta, rank_keep_override, verbose)
        type(parameters),         intent(in)  :: params
        class(image),             intent(in)  :: imgs(:)
        real,                     intent(in)  :: avg(:)
        type(diffmap_graph), target, intent(in) :: graph
        type(image), allocatable, intent(out) :: coeff_ptcls(:)
        logical, optional,        intent(in)  :: force_zero_theta
        integer, optional,        intent(in)  :: rank_keep_override
        logical, optional,        intent(in)  :: verbose
        type(parameters) :: params_pft
        type(polarft_calc) :: pftc
        type(image), allocatable :: imgs_work(:)
        type(image) :: avg_img, synth_img
        complex(sp), allocatable :: pfts(:,:,:), pfts_hat(:,:,:), pft_tmp(:,:)
        complex(sp), allocatable :: mode_coeff(:,:), mode_hat(:,:)
        integer :: nptcls, ldim(3), kfromto(2), pdim_srch(3), pftsz, nrots, nk
        integer :: i, k, mode, mode_max, rank_keep
        real :: smpd
        logical :: zero_theta, l_verbose
        nptcls = size(imgs)
        if( nptcls < 2 ) return
        if( graph%n /= nptcls ) return
        zero_theta = .false.
        if( present(force_zero_theta) ) zero_theta = force_zero_theta
        l_verbose = .true.
        if( present(verbose) ) l_verbose = verbose
        if( .not. zero_theta .and. .not. graph%n > 0 ) return
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        if( size(avg) /= product(ldim) ) return
        call prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        kfromto(1) = max(2, calc_fourier_index(params_pft%hp, params_pft%box, params_pft%smpd))
        kfromto(2) =        calc_fourier_index(params_pft%lp, params_pft%box, params_pft%smpd)
        if( kfromto(2) < kfromto(1) ) return
        call avg_img%new(ldim, smpd)
        call avg_img%unserialize(avg)
        call pftc%new(params_pft, 1, [1,nptcls], kfromto)
        pdim_srch = pftc%get_pdim_srch()
        pftsz     = pftc%get_pftsz()
        nrots     = pftc%get_nrots()
        nk        = kfromto(2) - kfromto(1) + 1
        mode_max  = 0
        if( present(rank_keep_override) )then
            rank_keep = min(max(1, rank_keep_override), nptcls)
        else
            rank_keep = min(max(2, params%k_nn), nptcls)
        endif
        allocate(pfts(pftsz,nk,nptcls), pfts_hat(pftsz,nk,nptcls))
        pfts     = cmplx(0.,0.,kind=sp)
        pfts_hat = cmplx(0.,0.,kind=sp)
        allocate(pft_tmp(pftsz,kfromto(1):kfromto(2)))
        imgs_work = copy_imgarr(imgs)
        call imgs_work(1)%memoize4polarize(pdim_srch)
        do i = 1,nptcls
            call imgs_work(i)%subtr(avg_img)
            call imgs_work(i)%fft()
            call pftc%polarize_ptcl_pft(imgs_work(i), i, pdim_srch, oversamp=.false.)
            call pftc%get_ptcl_pft(i, pft_tmp)
            do k = 1,nk
                pfts(:,k,i) = pft_tmp(:,kfromto(1) + k - 1)
            end do
        end do
        call pftc%kill
        call dealloc_imgarr(imgs_work)
        allocate(mode_coeff(nptcls,nk), mode_hat(nptcls,nk))
        do mode = 0,mode_max
            call extract_pft_angular_mode(pfts, kfromto, nrots, mode, mode_coeff)
            call project_pft_angular_mode(graph, mode, mode_coeff, mode_hat, rank_keep, zero_theta)
            call accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, mode, mode_hat)
        end do
        coeff_ptcls = copy_imgarr(imgs)
        call synth_img%new(ldim, smpd, wthreads=.false.)
        do i = 1,nptcls
            call synthesize_pft_residual(pfts_hat(:,:,i), kfromto, nrots, ldim, smpd, synth_img)
            call coeff_ptcls(i)%copy_fast(avg_img)
            call coeff_ptcls(i)%add(synth_img)
        end do
        call synth_img%kill
        if( l_verbose )then
            write(logfhandle,'(A,I8,A,I8,A,I8,A,I8,A,I8,A,I8)') 'Graph coefficient denoising: n=', nptcls, &
                ' modes=', mode_max, ' rank=', rank_keep, ' pftsz=', pftsz, ' kmin=', kfromto(1), ' kmax=', kfromto(2)
            call flush(logfhandle)
        endif
        call avg_img%kill
        deallocate(pfts, pfts_hat, pft_tmp, mode_coeff, mode_hat)
    end subroutine so2_coeffproj_denoise

    subroutine prepare_steerable_pft_params(params, ldim, smpd, params_pft)
        type(parameters), intent(in)  :: params
        integer,          intent(in)  :: ldim(3)
        real,             intent(in)  :: smpd
        type(parameters), intent(out) :: params_pft
        params_pft             = params
        params_pft%objfun      = 'cc'
        params_pft%cc_objfun   = OBJFUN_CC
        params_pft%ctf         = 'no'
        params_pft%ldim        = ldim
        params_pft%box         = ldim(1)
        params_pft%box_crop    = ldim(1)
        params_pft%smpd        = smpd
        params_pft%smpd_crop   = smpd
        if( params_pft%pftsz < 1 ) params_pft%pftsz = magic_pftsz((real(ldim(1)) - COSMSKHALFWIDTH) / 2.)
    end subroutine prepare_steerable_pft_params

    subroutine project_pft_angular_mode(graph, mode, mode_coeff, mode_hat, rank_keep, zero_theta)
        type(diffmap_graph), target, intent(in) :: graph
        integer,     intent(in)  :: mode, rank_keep
        complex(sp), intent(in)  :: mode_coeff(:,:)
        complex(sp), intent(out) :: mode_hat(:,:)
        logical,     intent(in)  :: zero_theta
        real, allocatable :: evals(:), evecs(:,:)
        complex(sp) :: coeff_ck
        integer :: n, nk, nkeep, a, i, k, info, max_basis
        n  = size(mode_coeff,1)
        nk = size(mode_coeff,2)
        mode_hat = cmplx(0.,0.,kind=sp)
        if( n < 1 .or. nk < 1 ) return
        nkeep = min(max(1, rank_keep), n)
        allocate(evals(nkeep), evecs(n,nkeep))
        max_basis = min(n, max(160, 8 * nkeep + 80))
        call sparse_eigh(graph_matvec, graph, n, nkeep, evals, evecs, tol=1.e-5, max_basis=max_basis, info=info)
        !$omp parallel do default(shared) private(k,a,i,coeff_ck) schedule(static) proc_bind(close)
        do k = 1,nk
            do a = 1,nkeep
                coeff_ck = cmplx(0.,0.,kind=sp)
                do i = 1,n
                    coeff_ck = coeff_ck + real(evecs(i,a), kind=sp) * mode_coeff(i,k)
                end do
                do i = 1,n
                    mode_hat(i,k) = mode_hat(i,k) + real(evecs(i,a), kind=sp) * coeff_ck
                end do
            end do
        end do
        !$omp end parallel do
        deallocate(evals, evecs)
    end subroutine project_pft_angular_mode

    subroutine extract_pft_angular_mode(pfts, kfromto, nrots, mode, mode_coeff)
        complex(sp), intent(in)  :: pfts(:,:,:)
        integer,     intent(in)  :: kfromto(2), nrots, mode
        complex(sp), intent(out) :: mode_coeff(:,:)
        complex(sp) :: acc, phase, val
        real :: theta, phase_arg, parity
        integer :: nptcls, pftsz, nk, i, irot, k
        nptcls = size(pfts,3)
        pftsz  = size(pfts,1)
        nk     = kfromto(2) - kfromto(1) + 1
        mode_coeff = cmplx(0.,0.,kind=sp)
        parity = merge(1., -1., mod(abs(mode), 2) == 0)
        !$omp parallel do collapse(2) default(shared) private(i,k,irot,acc,theta,phase_arg,phase,val) schedule(static) proc_bind(close)
        do i = 1,nptcls
            do k = 1,nk
                acc = cmplx(0.,0.,kind=sp)
                do irot = 1,pftsz
                    theta     = twopi * real(irot - 1) / real(nrots)
                    phase_arg = -real(mode) * theta
                    phase     = cmplx(cos(phase_arg), sin(phase_arg), kind=sp)
                    val       = pfts(irot,k,i)
                    acc       = acc + val * phase + conjg(val) * phase * real(parity, kind=sp)
                end do
                mode_coeff(i,k) = acc / real(nrots, kind=sp)
            end do
        end do
        !$omp end parallel do
    end subroutine extract_pft_angular_mode

    subroutine accumulate_pft_angular_mode(pfts_hat, kfromto, nrots, mode, mode_hat)
        complex(sp), intent(inout) :: pfts_hat(:,:,:)
        integer,     intent(in)    :: kfromto(2), nrots, mode
        complex(sp), intent(in)    :: mode_hat(:,:)
        complex(sp) :: phase
        real :: theta, phase_arg
        integer :: nptcls, pftsz, nk, i, irot, k
        nptcls = size(pfts_hat,3)
        pftsz  = size(pfts_hat,1)
        nk     = kfromto(2) - kfromto(1) + 1
        !$omp parallel do collapse(2) default(shared) private(i,k,irot,theta,phase_arg,phase) schedule(static) proc_bind(close)
        do i = 1,nptcls
            do k = 1,nk
                do irot = 1,pftsz
                    theta     = twopi * real(irot - 1) / real(nrots)
                    phase_arg = real(mode) * theta
                    phase     = cmplx(cos(phase_arg), sin(phase_arg), kind=sp)
                    pfts_hat(irot,k,i) = pfts_hat(irot,k,i) + mode_hat(i,k) * phase
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine accumulate_pft_angular_mode

    subroutine synthesize_pft_residual(pft_half, kfromto, nrots, ldim, smpd, img_out)
        complex(sp), intent(in)    :: pft_half(:,:)
        integer,     intent(in)    :: kfromto(2), nrots, ldim(3)
        real,        intent(in)    :: smpd
        type(image), intent(inout) :: img_out
        complex, allocatable :: cmat(:,:,:)
        real,    allocatable :: rho(:,:,:)
        real :: theta, x, y
        integer :: ashape(3), lims(3,2), pftsz, nk, irot, k, krad, h, l
        if( .not. img_out%exists() ) THROW_HARD('synthesize_pft_residual requires preconstructed output image')
        if( any(img_out%get_ldim() /= ldim) ) THROW_HARD('synthesize_pft_residual output image dimension mismatch')
        call img_out%zero_and_flag_ft
        ashape = img_out%get_array_shape()
        allocate(rho(ashape(1),ashape(2),ashape(3)), source=0.)
        lims  = img_out%loop_lims(3)
        pftsz = size(pft_half,1)
        nk    = kfromto(2) - kfromto(1) + 1
        do irot = 1,pftsz
            theta = twopi * real(irot - 1) / real(nrots)
            do k = 1,nk
                krad = kfromto(1) + k - 1
                x =  sin(theta) * real(krad)
                y = -cos(theta) * real(krad)
                call scatter_fourier_sample(img_out, rho, lims, x, y, pft_half(irot,k))
            end do
        end do
        cmat = img_out%get_cmat()
        do l = 1,ashape(3)
            do k = 1,ashape(2)
                do h = 1,ashape(1)
                    if( rho(h,k,l) > 1.e-6 )then
                        cmat(h,k,l) = cmat(h,k,l) / rho(h,k,l)
                    else
                        cmat(h,k,l) = cmplx(0.,0.)
                    endif
                end do
            end do
        end do
        call img_out%set_cmat(cmat)
        call img_out%ifft()
        deallocate(cmat, rho)
    end subroutine synthesize_pft_residual

    subroutine scatter_fourier_sample(img, rho, lims, x, y, comp)
        type(image), intent(inout) :: img
        real,        intent(inout) :: rho(:,:,:)
        integer,     intent(in)    :: lims(3,2)
        real,        intent(in)    :: x, y
        complex(sp), intent(in)    :: comp
        complex :: comp_here
        integer :: h0, h1, k0, k1, h, k, ih, ik, phys(3)
        real :: dh, dk, wh(2), wk(2), w
        if( .not. ieee_is_finite(real(comp)) .or. .not. ieee_is_finite(aimag(comp)) ) return
        h0 = floor(x)
        k0 = floor(y)
        h1 = h0 + 1
        k1 = k0 + 1
        dh = x - real(h0)
        dk = y - real(k0)
        wh = [1. - dh, dh]
        wk = [1. - dk, dk]
        do ih = 1,2
            h = merge(h0, h1, ih == 1)
            if( h < lims(1,1) .or. h > lims(1,2) ) cycle
            do ik = 1,2
                k = merge(k0, k1, ik == 1)
                if( k < lims(2,1) .or. k > lims(2,2) ) cycle
                w = wh(ih) * wk(ik)
                if( w <= 1.e-6 ) cycle
                comp_here = cmplx(real(comp) * w, aimag(comp) * w)
                call img%add([h,k,0], comp_here, phys_out=phys)
                rho(phys(1),phys(2),phys(3)) = rho(phys(1),phys(2),phys(3)) + w
            end do
        end do
    end subroutine scatter_fourier_sample

    subroutine estimate_diffmap_denoise_rank(params, graph, set_id, nptcls, den_rank, icm_converged, icm_iters, icm_score)
        use simple_diffusion_maps,  only: embed_graph
        type(parameters),    intent(in)  :: params
        type(diffmap_graph), intent(in)  :: graph
        integer,             intent(in)  :: set_id, nptcls
        integer,             intent(out) :: den_rank, icm_iters
        logical,             intent(out) :: icm_converged
        real,                intent(out) :: icm_score
        real, allocatable :: coords(:,:), eigvals(:)
        integer :: rank_scan, max_rank
        max_rank = max(1, nptcls - 2)
        if( params%neigs > 0 )then
            rank_scan = min(max(1, params%neigs), max_rank)
        else
            rank_scan = diffmap_denoise_auto_neigs_scan(nptcls)
        endif
        rank_scan = min(max(1, rank_scan), max_rank)
        call embed_graph(graph, rank_scan, coords, eigvals)
        if( allocated(eigvals) .and. size(eigvals) > 0 )then
            call select_spectral_rank_icm(eigvals, size(eigvals), den_rank, &
                                                 icm_converged, icm_iters, icm_score)
        else
            den_rank      = rank_scan
            icm_converged = .false.
            icm_iters     = 0
            icm_score     = huge(icm_score)
            write(logfhandle,'(A,I8,A,I8)') 'Diffmap denoise rank warning: no eigenspectrum; set=', set_id, &
                ' fallback_features=', den_rank
            call flush(logfhandle)
        endif
        den_rank = min(max(1, den_rank), nptcls)
        if( allocated(coords) ) deallocate(coords)
        if( allocated(eigvals) ) deallocate(eigvals)
    end subroutine estimate_diffmap_denoise_rank

    integer function diffmap_denoise_auto_neigs_scan(nptcls) result(neigs_scan)
        integer, intent(in) :: nptcls
        neigs_scan = min(DIFFMAP_NEIGS_AUTO_SCAN_MAX, max(1, nptcls - 2))
        if( nptcls > 3 ) neigs_scan = max(2, neigs_scan)
    end function diffmap_denoise_auto_neigs_scan

    subroutine select_spectral_rank_icm(eigvals, max_neigs, nkeep, converged, niter, score, min_rank)
        real,    intent(in)  :: eigvals(:)
        integer, intent(in)  :: max_neigs
        integer, intent(out) :: nkeep, niter
        logical, intent(out) :: converged
        real,    intent(out) :: score
        integer, optional, intent(in) :: min_rank
        real, allocatable :: spec(:)
        real :: smin, smax, delta, beta, complexity, trial_score, best_score
        integer :: i, n, nmin_rank, upper_rank, lower_rank, k_trial, trial_iter
        integer :: seeds(2)
        logical :: trial_converged, best_converged
        n = min(size(eigvals), max(1, max_neigs))
        nmin_rank = diffmap_denoise_min_neigs(n)
        if( present(min_rank) ) nmin_rank = max(1, min(n, min_rank))
        if( n <= 0 )then
            nkeep     = 1
            niter     = 0
            score     = 0.
            converged = .false.
            return
        endif
        allocate(spec(n), source=0.)
        do i = 1, n
            if( ieee_is_finite(eigvals(i)) .and. eigvals(i) > real(DTINY) )then
                spec(i) = log(eigvals(i))
            else
                spec(i) = log(real(DTINY))
            endif
        end do
        smin  = minval(spec)
        smax  = maxval(spec)
        delta = smax - smin
        if( .not. ieee_is_finite(delta) .or. delta <= 1.e-6 )then
            nkeep     = nmin_rank
            niter     = 0
            score     = 0.
            converged = .true.
            deallocate(spec)
            return
        endif
        spec = (spec - smin) / delta
        upper_rank = diffmap_denoise_initial_neigs_from_gap(spec, nmin_rank)
        upper_rank = max(nmin_rank, min(n, upper_rank))
        lower_rank = nint(DIFFMAP_DENOISE_ICM_RANK_LOWER_SEED_FRAC * real(upper_rank))
        lower_rank = max(nmin_rank, min(upper_rank, lower_rank))
        seeds(1) = upper_rank
        seeds(2) = lower_rank
        beta       = DIFFMAP_DENOISE_ICM_RANK_BETA_FRAC
        complexity = DIFFMAP_DENOISE_ICM_RANK_COMPLEXITY_FRAC
        best_score     = huge(best_score)
        nkeep          = upper_rank
        niter          = 0
        best_converged = .false.
        do i = 1, size(seeds)
            call run_diffmap_denoise_icm_rank_seed(spec, seeds(i), nmin_rank, upper_rank, &
                beta, complexity, k_trial, trial_score, trial_iter, trial_converged)
            if( trial_score < best_score )then
                best_score     = trial_score
                nkeep          = k_trial
                niter          = trial_iter
                best_converged = trial_converged
            endif
        end do
        score     = best_score
        converged = best_converged
        deallocate(spec)
    end subroutine select_spectral_rank_icm

    subroutine run_diffmap_denoise_icm_rank_seed(spec, seed_rank, nmin_rank, nmax_rank, beta, alpha, &
                                                 nkeep, score, niter, converged)
        real,             intent(in)  :: spec(:), beta, alpha
        integer,          intent(in)  :: seed_rank, nmin_rank, nmax_rank
        integer,          intent(out) :: nkeep, niter
        real,             intent(out) :: score
        logical,          intent(out) :: converged
        integer, allocatable :: labels(:), prev_labels(:)
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: iter, i, n, nchanged, maxits
        n = size(spec)
        allocate(labels(n), prev_labels(n), source=0)
        labels = 0
        labels(1:max(nmin_rank, min(nmax_rank, seed_rank))) = 1
        if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
        converged = .false.
        niter     = 0
        maxits    = max(DIFFMAP_DENOISE_ICM_RANK_MAXITS, n)
        do iter = 1, maxits
            prev_labels = labels
            call estimate_diffmap_denoise_icm_rank_stats(spec, prev_labels, mu_drop, mu_keep, var_drop, var_keep)
            do i = 1, nmax_rank
                call update_diffmap_denoise_icm_rank_site(spec, prev_labels, i, beta, alpha, nmin_rank, nmax_rank, &
                    mu_drop, mu_keep, var_drop, var_keep, labels(i))
            end do
            labels(1:nmin_rank) = 1
            if( nmax_rank < n ) labels(nmax_rank+1:n) = 0
            nchanged = count(labels /= prev_labels)
            score = score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
            niter = iter
            if( nchanged == 0 )then
                converged = .true.
                exit
            endif
        end do
        nkeep = diffmap_denoise_rank_prefix(labels, nmin_rank, nmax_rank)
        score = score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank)
        deallocate(labels, prev_labels)
    end subroutine run_diffmap_denoise_icm_rank_seed

    integer function diffmap_denoise_rank_prefix(labels, nmin_rank, nmax_rank) result(nkeep)
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        integer :: i
        nkeep = 0
        do i = 1, min(size(labels), nmax_rank)
            if( labels(i) == 1 ) nkeep = i
        end do
        nkeep = max(nmin_rank, min(nmax_rank, nkeep))
    end function diffmap_denoise_rank_prefix

    integer function diffmap_denoise_min_neigs(max_neigs) result(nmin_rank)
        integer, intent(in) :: max_neigs
        nmin_rank = 1
        if( max_neigs >= 2 ) nmin_rank = 2
    end function diffmap_denoise_min_neigs

    integer function diffmap_denoise_initial_neigs_from_gap(spec, nmin_rank) result(nkeep)
        real,    intent(in) :: spec(:)
        integer, intent(in) :: nmin_rank
        real :: best_gap, gap
        integer :: i
        if( size(spec) <= 1 )then
            nkeep = max(1, nmin_rank)
            return
        endif
        best_gap = -huge(best_gap)
        nkeep    = nmin_rank
        do i = nmin_rank, size(spec) - 1
            gap = spec(i) - spec(i+1)
            if( gap > best_gap )then
                best_gap = gap
                nkeep = i
            endif
        end do
    end function diffmap_denoise_initial_neigs_from_gap

    subroutine estimate_diffmap_denoise_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
        real,    intent(in)  :: spec(:)
        integer, intent(in)  :: labels(:)
        real,    intent(out) :: mu_drop, mu_keep, var_drop, var_keep
        real :: sum_drop, sum_keep, ssq_drop, ssq_keep
        integer :: i, n_drop, n_keep
        sum_drop = 0.; sum_keep = 0.; ssq_drop = 0.; ssq_keep = 0.
        n_drop = 0; n_keep = 0
        do i = 1, size(spec)
            if( labels(i) == 1 )then
                sum_keep = sum_keep + spec(i)
                ssq_keep = ssq_keep + spec(i) * spec(i)
                n_keep = n_keep + 1
            else
                sum_drop = sum_drop + spec(i)
                ssq_drop = ssq_drop + spec(i) * spec(i)
                n_drop = n_drop + 1
            endif
        end do
        if( n_keep > 0 )then
            mu_keep = sum_keep / real(n_keep)
            var_keep = max(1.e-4, ssq_keep / real(n_keep) - mu_keep * mu_keep)
        else
            mu_keep = 0.; var_keep = 1.
        endif
        if( n_drop > 0 )then
            mu_drop = sum_drop / real(n_drop)
            var_drop = max(1.e-4, ssq_drop / real(n_drop) - mu_drop * mu_drop)
        else
            mu_drop = 0.; var_drop = 1.
        endif
    end subroutine estimate_diffmap_denoise_icm_rank_stats

    subroutine update_diffmap_denoise_icm_rank_site(spec, labels, ind, beta, alpha, nmin_rank, nmax_rank, &
                                                    mu_drop, mu_keep, var_drop, var_keep, label_new)
        real,    intent(in)  :: spec(:), beta, alpha, mu_drop, mu_keep, var_drop, var_keep
        integer, intent(in)  :: labels(:), ind, nmin_rank, nmax_rank
        integer, intent(out) :: label_new
        real :: cost_drop, cost_keep
        integer :: kfree
        if( ind <= nmin_rank )then
            label_new = 1
            return
        endif
        if( ind > nmax_rank )then
            label_new = 0
            return
        endif
        kfree = min(max(nmin_rank, 4), nmax_rank)
        cost_keep = (spec(ind) - mu_keep)**2 / var_keep + alpha * real(max(0, ind - kfree)) / real(max(1, nmax_rank - kfree))
        cost_drop = (spec(ind) - mu_drop)**2 / var_drop
        if( ind > 1 )then
            if( labels(ind-1) /= 0 ) cost_drop = cost_drop + beta
            if( labels(ind-1) /= 1 ) cost_keep = cost_keep + beta
        endif
        if( ind < size(labels) )then
            if( labels(ind+1) /= 0 ) cost_drop = cost_drop + beta
            if( labels(ind+1) /= 1 ) cost_keep = cost_keep + beta
        endif
        if( ind > 1 )then
            if( labels(ind-1) == 0 ) cost_keep = cost_keep + 2. * beta
        endif
        if( ind < size(labels) )then
            if( labels(ind+1) == 1 ) cost_drop = cost_drop + 2. * beta
        endif
        if( cost_keep < cost_drop )then
            label_new = 1
        else
            label_new = 0
        endif
    end subroutine update_diffmap_denoise_icm_rank_site

    real function score_diffmap_denoise_icm_rank_solution(spec, labels, beta, alpha, nmin_rank, nmax_rank) result(score)
        real,    intent(in) :: spec(:), beta, alpha
        integer, intent(in) :: labels(:), nmin_rank, nmax_rank
        real :: mu_drop, mu_keep, var_drop, var_keep
        integer :: i, kfree
        call estimate_diffmap_denoise_icm_rank_stats(spec, labels, mu_drop, mu_keep, var_drop, var_keep)
        kfree = min(max(nmin_rank, 4), nmax_rank)
        score = 0.
        do i = 1, size(spec)
            if( labels(i) == 1 )then
                score = score + (spec(i) - mu_keep)**2 / var_keep
                score = score + alpha * real(max(0, i - kfree)) / real(max(1, nmax_rank - kfree))
            else
                score = score + (spec(i) - mu_drop)**2 / var_drop
            endif
        end do
        do i = 1, size(spec) - 1
            if( labels(i) /= labels(i+1) ) score = score + beta
            if( labels(i) == 0 .and. labels(i+1) == 1 ) score = score + 2. * beta
        end do
    end function score_diffmap_denoise_icm_rank_solution

    subroutine graph_nystrom_residual_preimage(raw_imgs, avg_img, graph, den_imgs, rank_keep)
        type(image),         intent(inout) :: raw_imgs(:)
        type(image),         intent(inout) :: avg_img
        type(diffmap_graph), intent(in)    :: graph
        type(image), allocatable, intent(out) :: den_imgs(:)
        integer,             intent(in)    :: rank_keep
        type(image), allocatable :: mode_imgs(:)
        type(image) :: resid_img
        real, allocatable :: evals(:), evecs(:,:), phi_ext(:,:)
        real :: coeff
        real :: smpd
        integer :: i, j, p, k, n, ldim(3), rank_used, nev, eig_idx, eig_info, max_basis
        n = size(raw_imgs)
        if( graph%n /= n ) THROW_HARD('nystrom preimage graph/image count mismatch; diffusion-map denoise')
        if( n < 3 ) THROW_HARD('nystrom preimage requires at least three graph nodes; diffusion-map denoise')
        ldim = avg_img%get_ldim()
        smpd = avg_img%get_smpd()
        rank_used = min(max(1, rank_keep), max(1, n - 2))
        nev = rank_used + 1
        allocate(evals(nev), evecs(n,nev), phi_ext(n,rank_used), source=0.)
        max_basis = min(n, max(160, 8 * nev + 80))
        call sparse_eigh(graph_matvec, graph, n, nev, evals, evecs, tol=1.e-5, max_basis=max_basis, info=eig_info)
        if( eig_info /= 0 ) THROW_HARD('sparse eigensolve failed in nystrom preimage; diffusion-map denoise')
        if( allocated(graph%wnorm) )then
            do k = 1, rank_used
                eig_idx = nev - k
                if( abs(evals(eig_idx)) <= real(DTINY) ) cycle
                !$omp parallel do default(shared) private(i,p,j,coeff) schedule(static) proc_bind(close)
                do i = 1, n
                    coeff = 0.
                    do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                        j = graph%colind(p)
                        if( j < 1 .or. j > n ) cycle
                        coeff = coeff + graph%wnorm(p) * evecs(j,eig_idx)
                    end do
                    phi_ext(i,k) = coeff / evals(eig_idx)
                end do
                !$omp end parallel do
            end do
        else
            do k = 1, rank_used
                eig_idx = nev - k
                if( abs(evals(eig_idx)) <= real(DTINY) ) cycle
                !$omp parallel do default(shared) private(i,p,j,coeff) schedule(static) proc_bind(close)
                do i = 1, n
                    coeff = 0.
                    do p = graph%rowptr(i), graph%rowptr(i+1) - 1
                        j = graph%colind(p)
                        if( j < 1 .or. j > n ) cycle
                        coeff = coeff + graph%w(p) * evecs(j,eig_idx)
                    end do
                    phi_ext(i,k) = coeff / evals(eig_idx)
                end do
                !$omp end parallel do
            end do
        endif
        allocate(den_imgs(n))
        allocate(mode_imgs(rank_used))
        do i = 1, n
            call den_imgs(i)%new(ldim, smpd, wthreads=.false.)
        end do
        do k = 1, rank_used
            call mode_imgs(k)%new(ldim, smpd, wthreads=.false.)
            call mode_imgs(k)%zero()
        end do
        call resid_img%new(ldim, smpd, wthreads=.false.)
        do i = 1, n
            call resid_img%copy_fast(raw_imgs(i))
            call resid_img%subtr(avg_img)
            do k = 1, rank_used
                eig_idx = nev - k
                coeff = evecs(i,eig_idx)
                if( abs(coeff) <= real(DTINY) ) cycle
                call mode_imgs(k)%add(resid_img, coeff)
            end do
        end do
        !$omp parallel do default(shared) private(i,k,coeff) schedule(dynamic) proc_bind(close)
        do i = 1, n
            call den_imgs(i)%copy_fast(avg_img)
            do k = 1, rank_used
                coeff = phi_ext(i,k)
                if( abs(coeff) <= real(DTINY) ) cycle
                call den_imgs(i)%add(mode_imgs(k), coeff)
            end do
        end do
        !$omp end parallel do
        write(logfhandle,'(A,I8,A,I8,A,I8,A,F8.4)') 'Diffmap Nystrom preimage: n=', n, ' rank=', rank_used, &
            ' nnz=', graph%nnz, ' lambda_min=', minval(evals(max(1,nev-rank_used):nev-1))
        call flush(logfhandle)
        do k = 1, rank_used
            if( mode_imgs(k)%exists() ) call mode_imgs(k)%kill
        end do
        if( resid_img%exists() ) call resid_img%kill
        deallocate(mode_imgs, evals, evecs, phi_ext)
    end subroutine graph_nystrom_residual_preimage

    subroutine calc_diffmap_residual_energy_ratio(raw_imgs, den_imgs, avg_img, ratio)
        type(image), intent(inout) :: raw_imgs(:), den_imgs(:), avg_img
        real,        intent(out)   :: ratio
        real, allocatable :: raw_vec(:), den_vec(:), avg_vec(:)
        real(kind=8) :: raw_ssq, den_ssq
        integer :: i
        if( size(raw_imgs) /= size(den_imgs) ) THROW_HARD('residual energy image count mismatch; diffusion-map denoise')
        avg_vec = avg_img%serialize()
        raw_ssq = 0.d0
        den_ssq = 0.d0
        !$omp parallel do default(shared) private(i,raw_vec,den_vec) schedule(static) proc_bind(close) &
        !$omp& reduction(+:raw_ssq,den_ssq)
        do i = 1, size(raw_imgs)
            raw_vec = raw_imgs(i)%serialize()
            den_vec = den_imgs(i)%serialize()
            if( size(raw_vec) /= size(avg_vec) .or. size(den_vec) /= size(avg_vec) )then
                THROW_HARD('residual energy image size mismatch; diffusion-map denoise')
            endif
            raw_ssq = raw_ssq + sum(real(raw_vec - avg_vec, kind=8) * real(raw_vec - avg_vec, kind=8))
            den_ssq = den_ssq + sum(real(den_vec - avg_vec, kind=8) * real(den_vec - avg_vec, kind=8))
            deallocate(raw_vec, den_vec)
        end do
        !$omp end parallel do
        ratio = real(sqrt(den_ssq / max(raw_ssq, real(DTINY, kind=8))))
        deallocate(avg_vec)
    end subroutine calc_diffmap_residual_energy_ratio

    subroutine calc_diffmap_reconstruction_error(raw_imgs, den_imgs, rmse, rel_rmse)
        type(image), intent(inout) :: raw_imgs(:), den_imgs(:)
        real,        intent(out)   :: rmse, rel_rmse
        real, allocatable :: raw_vec(:), den_vec(:)
        real(kind=8) :: ssq, raw_ssq
        integer(kind=8) :: npix_total
        integer :: i
        if( size(raw_imgs) /= size(den_imgs) ) THROW_HARD('reconstruction error image count mismatch; diffusion-map denoise')
        ssq        = 0.d0
        raw_ssq    = 0.d0
        npix_total = 0_8
        !$omp parallel do default(shared) private(i,raw_vec,den_vec) schedule(static) proc_bind(close) &
        !$omp& reduction(+:ssq,raw_ssq,npix_total)
        do i = 1, size(raw_imgs)
            raw_vec = raw_imgs(i)%serialize()
            den_vec = den_imgs(i)%serialize()
            if( size(raw_vec) /= size(den_vec) ) THROW_HARD('reconstruction error image size mismatch; diffusion-map denoise')
            ssq        = ssq + sum(real(raw_vec - den_vec, kind=8) * real(raw_vec - den_vec, kind=8))
            raw_ssq    = raw_ssq + sum(real(raw_vec, kind=8) * real(raw_vec, kind=8))
            npix_total = npix_total + int(size(raw_vec), kind=8)
            deallocate(raw_vec, den_vec)
        end do
        !$omp end parallel do
        if( npix_total > 0_8 )then
            rmse = real(sqrt(ssq / real(npix_total, kind=8)))
        else
            rmse = 0.
        endif
        rel_rmse = real(sqrt(ssq / max(raw_ssq, real(DTINY, kind=8))))
    end subroutine calc_diffmap_reconstruction_error

end module simple_diff_map_denoise
