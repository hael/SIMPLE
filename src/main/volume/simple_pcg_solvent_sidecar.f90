!@descr: opt-in soft solvent prior of the PCG base solve (pcg_solvent=yes)
!
!  What it does. Neither the spherical support nor the density envelope (drawn
!  at envmsklp, 20 A) sees solvent finer than the scale it was drawn at: the
!  dilation ring and skirt, and every cavity, channel and gap below ~20 A.
!  With pcg_solvent=yes the base system gains a real-space, position-dependent
!  ridge, (H + lambda I + Lambda_s) x = b with Lambda_s = lambda_s (1 - w(r)):
!  a Gaussian prior with position-dependent variance, the real-space twin of
!  the replay's P_tau. Each half is first solved prior-free; that pair's
!  FSC=0.143 sets the smoothing scale and each half's own prior-free map
!  yields its protein weight w(r) in [0,1] (this module): a Wang-type solvent
!  statistic (smoothed absolute density), an Otsu threshold inside the
!  production support, and a logistic of the statistic around that threshold
!  whose width is the spread of the statistic in the solvent class. Both
!  halves are then solved again, cold, with the same budget, ridge installed.
!  The prior is half-independent, so the pair stays gold standard; the FSC is
!  solvent-flattened and reported as such. Nothing is zeroed and nothing is
!  masked: where the data term is strong the prior is irrelevant, where it is
!  weak solvent is pulled toward zero; a misassigned voxel is over-regularized,
!  not deleted, and the partition is redrawn from prior-free maps every
!  iteration. The prior-free pair stays the base pair (FSC, NU competition
!  and calibration, evidence, _unfil); the prior'd pair is the base of the
!  closed-form replay and the pair the NU label field is applied to
!  (2026-09-21). The solve support is untouched.
!
!  Reporting. One PCG SOLVENT PRIOR line per half per reconstruction (the
!  smoothing scale, the Otsu threshold and the logistic width, the fraction of
!  the production support with w < 1/2, the mean weight, the relative ridge
!  coefficient), the prior-free pair's FSC, the even/odd weight agreement, and
!  the weight volumes pcg_solvent_weight_stateNN_even|odd.mrc (overwritten
!  each iteration). Provenance: solvent_prior=soft per_half base_pair lambda_rel=<x> auto|set.
!
!  Strength (2026-09-22). pcg_solvent_lambda not given = estimated per state
!  and iteration by cross-validation with the NU objective over the production
!  support, in closed form on the prior-free pair (estimate_solvent_prior_lambda):
!  x(lambda) ~ h/(h + lambda data_scale (1-w)) x_pre with h the real-space
!  diagonal of the data operator; J(lambda) = whitened Huber cross-half
!  prediction error of the shrunk halves against the prior-free other halves;
!  grid, argmin, one parabolic step in log lambda. Table and verdict in the
!  log (PCG SOLVENT PRIOR LAMBDA), edge and flat curves flagged. The real
!  re-solve then runs once at the chosen strength. pcg_solvent_check=yes
!  (shared-memory and distributed paths) prints the same grid by real
!  re-solves, with their residuals, beside it, to validate the closed form
!  on a data set.
module simple_pcg_solvent_sidecar
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image
implicit none

public :: build_solvent_prior_weight, pcg_solvent_stats, estimate_solvent_prior_lambda, pcg_solvent_lambda_stats, &
    &solvent_prior_cross_half_objective, PCG_SOLVENT_LAMBDA_GRID
private
#include "simple_local_flags.inc"

!> smoothing scale of the solvent statistic, as a multiple of the base pair's
!! FSC=0.143 resolution (Wang: about twice the resolution), floored in A
real,    parameter :: PCG_SOLVENT_LP_FACTOR = 2.0
real,    parameter :: PCG_SOLVENT_LP_MIN_A  = 8.0
!> voxels of the production support (window >= this) that take part in the
!! Otsu partition; the weight is still evaluated everywhere
real,    parameter :: PCG_SOLVENT_SUPPORT_MIN = 0.5
!> strength estimate (2026-09-22): the ridge coefficient is chosen by
!! cross-validation with the NU objective over the production support, in
!! closed form on the prior-free pair: x(lambda) ~ s(r) x_pre with
!! s = h / (h + lambda_rel data_scale (1-w(r))), h the real-space diagonal of
!! the data operator; J(lambda) = the whitened Huber cross-half prediction
!! error of s x_even against y_odd and s x_odd against y_even (image%nu_objective,
!! profile from the prior-free pair). Grid, argmin, one parabolic step in
!! log lambda between the neighbours. An edge argmin and a flat curve are
!! flagged: flat means the weight map, not the strength, is the limit.
real,    parameter :: PCG_SOLVENT_LAMBDA_GRID(8) = [0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
real,    parameter :: PCG_SOLVENT_LAMBDA_FLAT_TOL = 1.e-3 !< relative J range below which the curve is flat
type :: pcg_solvent_lambda_stats
    real    :: lambda_opt = 0.   !< chosen relative coefficient
    real    :: j_ref      = 0.   !< objective without the prior (lambda=0)
    real    :: j_min      = 0.   !< objective at the grid minimum
    real    :: j_gain     = 0.   !< 1 - j_min/j_ref
    logical :: l_edge     = .false. !< argmin on the grid edge
    logical :: l_flat     = .false. !< curve flat within PCG_SOLVENT_LAMBDA_FLAT_TOL
end type pcg_solvent_lambda_stats
type :: pcg_solvent_stats
    real    :: lp_smooth    = 0. !< smoothing scale used (A)
    real    :: thresh       = 0. !< Otsu threshold on the smoothed absolute density
    real    :: width        = 0. !< logistic width (spread of the statistic in the solvent class)
    real    :: solvent_frac = 0. !< fraction of the production support with w < 1/2
    real    :: weight_mean  = 0. !< mean protein weight over the production support
    integer :: n_support    = 0
end type pcg_solvent_stats

contains

    !> Protein weight w(r) in [0,1] from ONE prior-free base half, so that the
    !! prior of each half depends on its own data only (gold standard kept).
    !! base_half: the prior-free base half map at the crop box, real space.
    !! res0143: the prior-free base pair's FSC=0.143 resolution (A).
    !! base_support: the production support of the solve (density envelope,
    !! explicit pcg_mskfile, or the soft sphere the caller builds); it only
    !! selects the voxels the threshold is estimated on.
    !! weight: w(r), real space, same grid.
    subroutine build_solvent_prior_weight( state, half, base_half, res0143, base_support, lambda_rel, weight, stats )
        integer,                  intent(in)    :: state
        character(len=*),         intent(in)    :: half      !< 'even' | 'odd' (reporting only)
        class(image),             intent(in)    :: base_half
        real,                     intent(in)    :: res0143
        class(image), target,     intent(in)    :: base_support
        real,                     intent(in)    :: lambda_rel
        type(image),              intent(inout) :: weight
        type(pcg_solvent_stats),  intent(out)   :: stats
        real(kind=c_float), pointer :: rmat_s(:,:,:), rmat_sup(:,:,:)
        real,    allocatable :: vals(:)
        integer :: ldim(3), i, j, k, n, n_sup, n_solv, n_below
        real    :: smpd, lp, thresh, width, arg
        real(dp) :: s_solv, ss_solv, s_w
        ldim = base_half%get_ldim()
        smpd = base_half%get_smpd()
        if( any(base_support%get_ldim() /= ldim) ) THROW_HARD('solvent prior: base support and base pair are on different grids')
        if( base_half%is_ft() .or. base_support%is_ft() ) THROW_HARD('solvent prior: inputs must be in real space')
        if( res0143 <= 0. ) THROW_HARD('solvent prior: requires a positive base-pair resolution')
        ! smoothed absolute density at a multiple of the working resolution
        lp = max(PCG_SOLVENT_LP_MIN_A, PCG_SOLVENT_LP_FACTOR * res0143)
        lp = min(lp, real(ldim(1)) * smpd / 4.)   ! never coarser than a quarter of the box
        call weight%copy(base_half)
        call weight%get_rmat_ptr(rmat_s)
        rmat_s(:ldim(1),:ldim(2),:ldim(3)) = abs(rmat_s(:ldim(1),:ldim(2),:ldim(3)))
        nullify(rmat_s)
        call weight%bp(0., lp)   ! transforms in, filters, and transforms back to real space
        call weight%get_rmat_ptr(rmat_s)
        ! Otsu threshold of the statistic inside the production support
        call base_support%get_rmat_ptr(rmat_sup)
        n = count(rmat_sup(:ldim(1),:ldim(2),:ldim(3)) >= PCG_SOLVENT_SUPPORT_MIN)
        if( n < 1 ) THROW_HARD('solvent prior: the production support is empty')
        allocate(vals(n))
        n = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( rmat_sup(i,j,k) >= PCG_SOLVENT_SUPPORT_MIN )then
                        n = n + 1
                        vals(n) = rmat_s(i,j,k)
                    endif
                end do
            end do
        end do
        thresh = 0.
        call otsu(n, vals, thresh)
        ! logistic width: the spread of the statistic in the solvent class
        n_sup   = n
        n_solv  = count(vals < thresh)
        if( n_solv < 2 .or. n_solv >= n_sup ) THROW_HARD('solvent prior: Otsu partition degenerate inside the support')
        s_solv  = sum(real(vals,dp), mask=vals < thresh)
        ss_solv = sum(real(vals,dp)**2, mask=vals < thresh)
        s_solv  = s_solv / real(n_solv,dp)
        width   = real(sqrt(max(ss_solv / real(n_solv,dp) - s_solv**2, 0.0_dp)))
        width   = max(width, 1.e-3 * max(abs(thresh), TINY))
        deallocate(vals)
        ! w(r) = logistic((s - t) / width), everywhere on the grid
        n_below = 0
        s_w     = 0.0_dp
        !$omp parallel do collapse(3) default(shared) private(i,j,k,arg) reduction(+:n_below,s_w) proc_bind(close) schedule(static)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    arg = (rmat_s(i,j,k) - thresh) / width
                    arg = min(40., max(-40., arg))
                    rmat_s(i,j,k) = 1. / (1. + exp(-arg))
                    if( rmat_sup(i,j,k) >= PCG_SOLVENT_SUPPORT_MIN )then
                        s_w = s_w + real(rmat_s(i,j,k),dp)
                        if( rmat_s(i,j,k) < 0.5 ) n_below = n_below + 1
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        nullify(rmat_s, rmat_sup)
        stats%lp_smooth    = lp
        stats%thresh       = thresh
        stats%width        = width
        stats%n_support    = n_sup
        stats%solvent_frac = real(n_below) / real(max(1,n_sup))
        stats%weight_mean  = real(s_w / real(max(1,n_sup),dp))
        if( lambda_rel > 0. )then
            write(logfhandle,'(A,I0,A,F6.1,A,ES10.3,A,ES10.3,A,F6.2,A,F5.2,A,F6.2,A)') &
                &'>>> PCG SOLVENT PRIOR: STATE ', state, ' '//trim(half)//', smoothing ', lp, ' A, threshold ', thresh, &
                &', width ', width, ', solvent ', 100.*stats%solvent_frac, ' % of the production support, mean weight ', &
                &stats%weight_mean, ', ridge lambda_rel ', lambda_rel, ' (soft; per-half base re-solve)'
        else
            write(logfhandle,'(A,I0,A,F6.1,A,ES10.3,A,ES10.3,A,F6.2,A,F5.2,A)') &
                &'>>> PCG SOLVENT PRIOR: STATE ', state, ' '//trim(half)//', smoothing ', lp, ' A, threshold ', thresh, &
                &', width ', width, ', solvent ', 100.*stats%solvent_frac, ' % of the production support, mean weight ', &
                &stats%weight_mean, ', ridge lambda_rel auto (soft; per-half base re-solve)'
        endif
    end subroutine build_solvent_prior_weight

    !> Relative ridge coefficient by cross-validation, closed form on the
    !! prior-free pair (see the type comment above). x_even/odd: prior-free
    !! base halves (real space); weight_even/odd: their protein weights;
    !! support: the production support of the solve (window >= 1/2 takes part);
    !! h: real-space diagonal of the data operator; data_scale: the reference
    !! the relative coefficient multiplies. No operator, no solve.
    subroutine estimate_solvent_prior_lambda( state, x_even, x_odd, weight_even, weight_odd, support, &
            &h_even, h_odd, data_scale_even, data_scale_odd, lambda_opt, stats )
        integer,                       intent(in)    :: state
        class(image),                  intent(in)    :: x_even, x_odd, weight_even, weight_odd
        class(image), target,          intent(in)    :: support
        real,                          intent(in)    :: h_even, h_odd, data_scale_even, data_scale_odd !< per operator
        real,                          intent(out)   :: lambda_opt
        type(pcg_solvent_lambda_stats),intent(out)   :: stats
        type(image) :: s_even, s_odd
        real(kind=c_float), pointer :: rmat_sup(:,:,:)
        real,    allocatable :: sigma_r(:), diff(:,:,:), jgrid(:)
        logical, allocatable :: lmask(:,:,:)
        real    :: rmax, jref, y0, y1, y2, xa, xb, xc, denom, xopt
        integer :: ldim(3), ig, ng, imin, nmask
        ldim = x_even%get_ldim()
        if( any(x_odd%get_ldim() /= ldim) .or. any(weight_even%get_ldim() /= ldim) .or. &
            &any(weight_odd%get_ldim() /= ldim) .or. any(support%get_ldim() /= ldim) ) &
            &THROW_HARD('dimension mismatch; estimate_solvent_prior_lambda')
        if( min(h_even, h_odd, data_scale_even, data_scale_odd) <= 0. ) &
            &THROW_HARD('the data operator scales must be positive; estimate_solvent_prior_lambda')
        call support%get_rmat_ptr(rmat_sup)
        allocate(lmask(ldim(1),ldim(2),ldim(3)))
        lmask = rmat_sup(:ldim(1),:ldim(2),:ldim(3)) >= PCG_SOLVENT_SUPPORT_MIN
        nullify(rmat_sup)
        nmask = count(lmask)
        if( nmask < 1 ) THROW_HARD('empty production support; estimate_solvent_prior_lambda')
        ! whitening of the prior-free pair, the same profile the NU unary uses
        call x_even%nu_objective_noise_profile(x_odd, lmask, sigma_r, rmax)
        allocate(diff(ldim(1),ldim(2),ldim(3)), source=0.)
        ! reference: no prior (the shrunk halves are the prior-free halves)
        call x_even%nu_objective(x_even, x_odd, x_odd, diff, lmask, sigma_r, rmax)
        jref = sum(diff, mask=lmask) / real(nmask)
        ng = size(PCG_SOLVENT_LAMBDA_GRID)
        allocate(jgrid(ng), source=0.)
        call s_even%new(ldim, x_even%get_smpd())
        call s_odd%new( ldim, x_even%get_smpd())
        write(logfhandle,'(A,I0,A)') '>>> PCG SOLVENT PRIOR LAMBDA: STATE ', state, &
            &', cross-half objective over the production support, closed-form shrink of the prior-free pair'
        write(logfhandle,'(A,ES10.3,A,ES10.3,A,ES10.3,A,ES10.3)') '    h (real-space data diagonal) even/odd ', h_even, '/', h_odd, &
            &', data_scale even/odd ', data_scale_even, '/', data_scale_odd
        write(logfhandle,'(A)') '    lambda_rel     J/J(0)'
        write(logfhandle,'(A,F10.4)') '        0.000  ', 1.0
        do ig = 1, ng
            call shrink_half(x_even, weight_even, h_even, PCG_SOLVENT_LAMBDA_GRID(ig) * data_scale_even, s_even)
            call shrink_half(x_odd,  weight_odd,  h_odd,  PCG_SOLVENT_LAMBDA_GRID(ig) * data_scale_odd,  s_odd)
            call x_even%nu_objective(s_even, x_odd, s_odd, diff, lmask, sigma_r, rmax)
            jgrid(ig) = sum(diff, mask=lmask) / real(nmask)
            write(logfhandle,'(A,F8.3,A,F10.4)') '     ', PCG_SOLVENT_LAMBDA_GRID(ig), '  ', jgrid(ig) / max(TINY, jref)
        enddo
        imin = minloc(jgrid, dim=1)
        stats%j_ref  = jref
        stats%j_min  = jgrid(imin)
        stats%j_gain = 1. - jgrid(imin) / max(TINY, jref)
        stats%l_flat = (maxval(jgrid) - minval(jgrid)) / max(TINY, jref) < PCG_SOLVENT_LAMBDA_FLAT_TOL
        stats%l_edge = imin == 1 .or. imin == ng
        lambda_opt = PCG_SOLVENT_LAMBDA_GRID(imin)
        if( jref <= jgrid(imin) )then
            ! no prior beats every strength on the grid: the weakest grid
            ! value, flagged; the prior is kept so that the run stays
            ! comparable with an explicit setting
            lambda_opt   = PCG_SOLVENT_LAMBDA_GRID(1)
            stats%l_edge = .true.
        else if( .not. stats%l_edge )then
            ! one parabolic step in log lambda through the three points
            xa = log(PCG_SOLVENT_LAMBDA_GRID(imin-1)); y0 = jgrid(imin-1)
            xb = log(PCG_SOLVENT_LAMBDA_GRID(imin));   y1 = jgrid(imin)
            xc = log(PCG_SOLVENT_LAMBDA_GRID(imin+1)); y2 = jgrid(imin+1)
            denom = (xb-xa)*(y1-y2) - (xb-xc)*(y1-y0)
            if( abs(denom) > TINY )then
                xopt = xb - 0.5 * ((xb-xa)**2*(y1-y2) - (xb-xc)**2*(y1-y0)) / denom
                xopt = min(max(xopt, xa), xc)
                lambda_opt = exp(xopt)
            endif
        endif
        stats%lambda_opt = lambda_opt
        write(logfhandle,'(A,F8.3,A,F6.1,A)') '>>> PCG SOLVENT PRIOR LAMBDA: chosen lambda_rel ', lambda_opt, &
            &' (objective ', 100.*stats%j_gain, ' % below no prior)'
        if( stats%l_edge ) write(logfhandle,'(A)') '>>> PCG SOLVENT PRIOR LAMBDA: minimum on the grid edge'
        if( stats%l_flat ) write(logfhandle,'(A)') &
            &'>>> PCG SOLVENT PRIOR LAMBDA: flat curve; the weight map, not the strength, is the limit'
        call s_even%kill
        call s_odd%kill
        deallocate(sigma_r, diff, jgrid, lmask)

    contains

        !> closed-form action of the ridge on one half: s = h / (h + lam_s (1-w))
        subroutine shrink_half( x, w, h, lam_s, sx )
            class(image), intent(in)    :: x, w
            real,         intent(in)    :: h, lam_s
            type(image),  intent(inout) :: sx
            real(kind=c_float), pointer :: rx(:,:,:), rw(:,:,:), rs(:,:,:)
            integer :: i, j, k
            call x%get_rmat_ptr(rx)
            call w%get_rmat_ptr(rw)
            call sx%get_rmat_ptr(rs)
            !$omp parallel do collapse(3) default(shared) private(i,j,k) schedule(static) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        rs(i,j,k) = rx(i,j,k) * h / (h + lam_s * (1. - min(1., max(0., rw(i,j,k)))))
                    end do
                end do
            end do
            !$omp end parallel do
            nullify(rx, rw, rs)
        end subroutine shrink_half

    end subroutine estimate_solvent_prior_lambda

    !> The same objective on any candidate pair (the dev check compares the
    !! closed form with real re-solves): whitened Huber cross-half prediction
    !! error of the candidate halves against the prior-free other halves, mean
    !! over the production support (window >= 1/2).
    real function solvent_prior_cross_half_objective( x_even, x_odd, cand_even, cand_odd, support ) result( j )
        class(image),         intent(in) :: x_even, x_odd, cand_even, cand_odd
        class(image), target, intent(in) :: support
        real(kind=c_float), pointer :: rmat_sup(:,:,:)
        real,    allocatable :: sigma_r(:), diff(:,:,:)
        logical, allocatable :: lmask(:,:,:)
        real    :: rmax
        integer :: ldim(3), nmask
        ldim = x_even%get_ldim()
        call support%get_rmat_ptr(rmat_sup)
        allocate(lmask(ldim(1),ldim(2),ldim(3)))
        lmask = rmat_sup(:ldim(1),:ldim(2),:ldim(3)) >= PCG_SOLVENT_SUPPORT_MIN
        nullify(rmat_sup)
        nmask = count(lmask)
        if( nmask < 1 ) THROW_HARD('empty production support; solvent_prior_cross_half_objective')
        call x_even%nu_objective_noise_profile(x_odd, lmask, sigma_r, rmax)
        allocate(diff(ldim(1),ldim(2),ldim(3)), source=0.)
        call x_even%nu_objective(cand_even, x_odd, cand_odd, diff, lmask, sigma_r, rmax)
        j = sum(diff, mask=lmask) / real(nmask)
        deallocate(sigma_r, diff, lmask)
    end function solvent_prior_cross_half_objective

end module simple_pcg_solvent_sidecar
