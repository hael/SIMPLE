!@descr: opt-in soft solvent prior of the PCG regularized solve (pcg_solvent=yes)
!
!  What it does. The regularized (shipped) PCG solve is (H + P_tau + lambda I) x = b
!  on the production support. Neither the spherical support nor the density
!  envelope (drawn at envmsklp, 20 A) sees solvent finer than the scale it was
!  drawn at: the dilation ring and skirt, and every cavity, channel and gap
!  below ~20 A. With pcg_solvent=yes the regularized system gains a real-space,
!  position-dependent ridge, (H + P_tau + lambda I + Lambda_s) x = b with
!  Lambda_s = lambda_s (1 - w(r)): a Gaussian prior with position-dependent
!  variance, the real-space twin of P_tau. The protein weight w(r) in [0,1] is
!  built here per half from that half's own base map at the working
!  resolution, so the prior is half-independent (gold standard kept): a
!  Wang-type solvent statistic (smoothed absolute density, smoothing scale a
!  multiple of the base pair's FSC=0.143 resolution), an Otsu threshold inside
!  the production support, and a logistic of the statistic around that
!  threshold whose width is the spread of the statistic in the solvent class.
!  Nothing is zeroed and nothing is masked: where the data term is strong the
!  prior is irrelevant, where it is weak solvent is pulled toward zero. A
!  misassigned voxel is over-regularized, not deleted, and the partition is
!  redrawn from the base pair every iteration.
!
!  What it does not touch. The solve support, the base solve, the base pair,
!  the FSC/cFAR and the NU candidate bank are unchanged: the FSC oracle carries
!  no extra mask. With pcg_solvent=no (the default) no code in this module runs.
!
!  Reporting. One PCG SOLVENT PRIOR line per state per reconstruction: the
!  smoothing scale, the Otsu threshold and the logistic width, the fraction of
!  the production support with w < 1/2, the mean weight, and the relative
!  ridge coefficient. The shipped pair's provenance sidecar records
!  solvent_prior=soft lambda_rel=<x>.
module simple_pcg_solvent_sidecar
use simple_core_module_api
use simple_parameters, only: parameters
use simple_image,      only: image
implicit none

public :: build_solvent_prior_weight, pcg_solvent_stats
private
#include "simple_local_flags.inc"

!> smoothing scale of the solvent statistic, as a multiple of the base pair's
!! FSC=0.143 resolution (Wang: about twice the resolution), floored in A
real,    parameter :: PCG_SOLVENT_LP_FACTOR = 2.0
real,    parameter :: PCG_SOLVENT_LP_MIN_A  = 8.0
!> voxels of the production support (window >= this) that take part in the
!! Otsu partition; the weight is still evaluated everywhere
real,    parameter :: PCG_SOLVENT_SUPPORT_MIN = 0.5
type :: pcg_solvent_stats
    real    :: lp_smooth    = 0. !< smoothing scale used (A)
    real    :: thresh       = 0. !< Otsu threshold on the smoothed absolute density
    real    :: width        = 0. !< logistic width (spread of the statistic in the solvent class)
    real    :: solvent_frac = 0. !< fraction of the production support with w < 1/2
    real    :: weight_mean  = 0. !< mean protein weight over the production support
    integer :: n_support    = 0
end type pcg_solvent_stats

contains

    !> Protein weight w(r) in [0,1] from ONE base half, so that the prior of
    !! each regularized half depends on its own data only (gold standard kept).
    !! base_half: the base (unfil) half map at the crop box, real space.
    !! res0143: the base pair's FSC=0.143 resolution (A).
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
        write(logfhandle,'(A,I0,A,F6.1,A,ES10.3,A,ES10.3,A,F6.2,A,F5.2,A,F6.2,A)') &
            &'>>> PCG SOLVENT PRIOR: STATE ', state, ' '//trim(half)//', smoothing ', lp, ' A, threshold ', thresh, &
            &', width ', width, ', solvent ', 100.*stats%solvent_frac, ' % of the production support, mean weight ', &
            &stats%weight_mean, ', ridge lambda_rel ', lambda_rel, ' (soft; support, base pair and FSC untouched)'
    end subroutine build_solvent_prior_weight

end module simple_pcg_solvent_sidecar
