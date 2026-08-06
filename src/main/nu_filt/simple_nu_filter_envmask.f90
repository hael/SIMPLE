!@descr: NU-evidence-driven envelope masking for volume-domain nonuniform filtering
!
! The NU unary is a cross-half prediction error. In solvent the two half maps are
! uncorrelated at every bandwidth, so no candidate can beat the coarsest one and
! the per-voxel improvement over that baseline concentrates near zero. Inside real
! density the objective has a genuine minimum at the local SNR crossover, so the
! improvement is positive and large. That improvement, not the selected label, is
! the statistic this submodule segments on:
!
!     margin(v) = dmats_mask(v, coarsest) - min over candidates of dmats_mask(v, c)
!
! The optional scale-free form is baseline / best - 1, with a robust cost
! floor. Unlike a fractional reduction, this ratio remains compatible with an
! unbounded median-plus-MAD decision threshold.
!
! The selected label is a poor substitute: the coarsest bank member is a saturating
! bin that absorbs both solvent and genuinely coarse density, and taking the argmin
! discards the confidence.
!
! Segmentation is a binary MRF solved with the same parallel 8-color ICM schedule
! the ordered-label Potts prior uses. beta regularizes boundary area only; it does
! not enforce a connected result. Topology is the caller's job, via the connected
! component and morphology tail in simple_image_msk.
!
submodule (simple_nu_filter) simple_nu_filter_envmask
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine calc_nu_evidence_margin( margin, lp_smooth, l_relative )
        real, allocatable, intent(inout) :: margin(:)
        real,    optional, intent(in)    :: lp_smooth
        logical, optional, intent(in)    :: l_relative
        real, allocatable :: full(:,:,:), base_full(:,:,:), tmp(:,:,:)
        real    :: lp, floor_val, base_val, best_val
        logical :: l_rel
        integer :: imask, i, j, k
        if( .not.allocated(nu_ev_base) .or. .not.allocated(nu_ev_best) ) &
            &THROW_HARD('raw evidence not allocated; run setup_nu_dmats before calc_nu_evidence_margin')
        if( n_nu_mask < 1 ) THROW_HARD('empty NU support mask; calc_nu_evidence_margin')
        lp    = 8.0
        l_rel = .false.
        if( present(lp_smooth) )then
            if( lp_smooth > TINY ) lp = lp_smooth
        endif
        if( present(l_relative) ) l_rel = l_relative
        if( allocated(margin) ) deallocate(margin)
        allocate(margin(n_nu_mask), source=0.)
        allocate(full(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(tmp(ldim(1),ldim(2),ldim(3)),  source=0.)
        ! The baseline grid is only needed by the scale-free form. Allocating it
        ! unconditionally costs a whole box^3 real per state for nothing.
        if( l_rel ) allocate(base_full(ldim(1),ldim(2),ldim(3)), source=0.)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            full(i,j,k) = max(0., nu_ev_base(imask) - nu_ev_best(imask))
        end do
        !$omp end parallel do
        if( l_rel )then
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                base_full(i,j,k) = nu_ev_base(imask)
            end do
            !$omp end parallel do
        endif
        ! Smoothing the difference is equivalent to smoothing both terms with the
        ! same kernel, which is exactly what the candidate-scale path fails to do.
        call smooth_nu_objective(full, tmp, lp)
        if( l_rel ) call smooth_nu_objective(base_full, tmp, lp)
        floor_val = 0.
        if( l_rel ) floor_val = nu_evidence_baseline_floor(base_full)
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(imask,i,j,k,base_val,best_val) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( full(i,j,k) >= NU_EVIDENCE_INVALID )then
                margin(imask) = 0.
                cycle
            endif
            margin(imask) = max(0., full(i,j,k))
            if( l_rel )then
                ! The cost-improvement ratio is scale-free without the hard
                ! upper bound of a fractional reduction. That bound is
                ! incompatible with a median-plus-MAD threshold, which can
                ! legitimately exceed one.
                base_val = base_full(i,j,k)
                if( base_val >= NU_EVIDENCE_INVALID ) base_val = 0.
                base_val = max(base_val, floor_val)
                best_val = max(base_val - margin(imask), floor_val)
                margin(imask) = base_val / best_val - 1.
            endif
        end do
        !$omp end parallel do
        deallocate(full, tmp)
        if( allocated(base_full) ) deallocate(base_full)
    end subroutine calc_nu_evidence_margin

    !>  Robust floor for the relative cost ratio, so that near-zero baseline or
    !!  best-candidate costs cannot manufacture arbitrarily large evidence.
    module real function nu_evidence_baseline_floor( base_full ) result( floor_val )
        real, intent(in) :: base_full(:,:,:)
        real, allocatable :: work(:)
        integer :: imask, i, j, k, n
        real    :: val
        allocate(work(n_nu_mask), source=0.)
        n = 0
        do imask = 1, n_nu_mask
            i   = nu_mask_vox(1,imask)
            j   = nu_mask_vox(2,imask)
            k   = nu_mask_vox(3,imask)
            val = base_full(i,j,k)
            if( val >= NU_EVIDENCE_INVALID ) cycle
            n       = n + 1
            work(n) = val
        end do
        floor_val = TINY
        if( n > 0 ) floor_val = max(TINY, 0.1 * median_nocopy(work(:n)))
        deallocate(work)
    end function nu_evidence_baseline_floor

    module subroutine calc_nu_evidence_score( margin, nsigma, score, stats )
        real,                    intent(in)    :: margin(:)
        real,                    intent(in)    :: nsigma
        real, allocatable,       intent(inout) :: score(:)
        type(nu_envmask_stats),  intent(inout) :: stats
        real, allocatable :: work(:)
        real    :: med, mad_val, denom
        integer :: n
        n = size(margin)
        if( n < 1 ) THROW_HARD('empty margin vector; calc_nu_evidence_score')
        ! The null is estimated from the support itself rather than from geometry.
        ! With a generous support the solvent is the majority population, so the
        ! median and MAD of the margin describe the no-evidence distribution. The
        ! reported signal fraction is what tells the caller whether that held.
        allocate(work(n), source=margin)
        med     = median_nocopy(work)
        mad_val = mad_gau(margin, med)
        deallocate(work)
        denom = max(mad_val, TINY)
        stats%null_med = med
        stats%null_mad = mad_val
        stats%thres    = med + nsigma * denom
        if( allocated(score) ) deallocate(score)
        allocate(score(n), source=0.)
        score = (margin - stats%thres) / denom
    end subroutine calc_nu_evidence_score

    module subroutine add_nu_evidence_density( vol_dens, weight, score, stats )
        class(image), target,    intent(in)    :: vol_dens
        real,                    intent(in)    :: weight
        real, allocatable,       intent(inout) :: score(:)
        type(nu_envmask_stats),  intent(inout) :: stats
        real(kind=c_float), pointer :: rmat(:,:,:) => null()
        real, allocatable :: dens(:), work(:)
        real    :: med, mad_val, denom
        integer :: imask, i, j, k
        if( abs(weight) <= TINY ) return
        if( .not.allocated(score) ) THROW_HARD('score not allocated; add_nu_evidence_density')
        if( any(vol_dens%get_ldim() /= ldim) ) &
            &THROW_HARD('density volume dimension mismatch; add_nu_evidence_density')
        call vol_dens%get_rmat_ptr(rmat)
        allocate(dens(n_nu_mask), source=0.)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            dens(imask) = rmat(i,j,k)
        end do
        !$omp end parallel do
        allocate(work(n_nu_mask), source=dens)
        med     = median_nocopy(work)
        mad_val = mad_gau(dens, med)
        deallocate(work)
        denom = max(mad_val, TINY)
        stats%dens_med    = med
        stats%dens_mad    = mad_val
        stats%dens_weight = weight
        ! Additive so that strong density can hold in a poorly ordered region that
        ! the resolution evidence alone would carve out. This is the term that
        ! protects flexible periphery, which is otherwise indistinguishable from
        ! solvent by cross-half consistency.
        score = score + weight * (dens - med) / denom
        deallocate(dens)
    end subroutine add_nu_evidence_density

    module subroutine segment_nu_evidence( score, p, lmask, stats )
        real,                    intent(in)    :: score(:)
        type(nu_envmask_params), intent(in)    :: p
        logical, allocatable,    intent(inout) :: lmask(:,:,:)
        type(nu_envmask_stats),  intent(inout) :: stats
        logical, allocatable :: lab(:,:,:)
        integer :: iter, color, imask, i, j, k, ineigh, ni, nj, nk
        integer :: n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, deg, nsig, nchanged
        real    :: e_sig, e_sol, beta
        logical :: newlab
        if( size(score) /= n_nu_mask ) THROW_HARD('score size mismatch; segment_nu_evidence')
        beta = max(0., p%beta)
        if( allocated(lmask) ) deallocate(lmask)
        allocate(lab(ldim(1),ldim(2),ldim(3)), source=.false.)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            lab(i,j,k) = score(imask) > 0.
        end do
        !$omp end parallel do
        stats%n_seed = count(lab)
        if( beta > TINY )then
            ! Same 8-color schedule as the ordered-label prior: two voxels sharing
            ! a color differ by at least two in some coordinate, so they are never
            ! 26-neighbors and can be updated concurrently.
            do iter = 1, max(1,p%maxits)
                nchanged = 0
                do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                    !$omp parallel do schedule(static) default(shared) &
                    !$omp private(imask,i,j,k,ineigh,ni,nj,nk,n_full,nsz,deg,nsig,e_sig,e_sol,newlab) &
                    !$omp reduction(+:nchanged) proc_bind(close)
                    do imask = 1, n_nu_mask
                        i = nu_mask_vox(1,imask)
                        j = nu_mask_vox(2,imask)
                        k = nu_mask_vox(3,imask)
                        if( nu_label_smooth_color(i,j,k) /= color ) cycle
                        call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                        deg  = 0
                        nsig = 0
                        do ineigh = 1, nsz
                            ni = n_full(1,ineigh)
                            nj = n_full(2,ineigh)
                            nk = n_full(3,ineigh)
                            if( .not.nu_lmask(ni,nj,nk) ) cycle
                            deg = deg + 1
                            if( lab(ni,nj,nk) ) nsig = nsig + 1
                        end do
                        e_sig = -score(imask)
                        e_sol =  score(imask)
                        if( deg > 0 )then
                            ! Degree-normalized so support-boundary voxels are not
                            ! penalized for their missing neighbors.
                            e_sig = e_sig + beta * real(deg - nsig) / real(deg)
                            e_sol = e_sol + beta * real(nsig)       / real(deg)
                        endif
                        newlab = e_sig < e_sol
                        if( newlab .neqv. lab(i,j,k) )then
                            lab(i,j,k) = newlab
                            nchanged   = nchanged + 1
                        endif
                    end do
                    !$omp end parallel do
                end do
                stats%nits = iter
                if( nu_l_report ) write(logfhandle,'(A,I2,A,I10)') &
                    &'>>> NU envelope ICM iteration ', iter, ' changed voxels: ', nchanged
                if( nchanged == 0 ) exit
            end do
        endif
        stats%n_signal = count(lab)
        call move_alloc(lab, lmask)
    end subroutine segment_nu_evidence

    module subroutine nu_evidence_envelope( p, lmask, stats, vol_dens )
        type(nu_envmask_params), intent(in)    :: p
        logical, allocatable,    intent(inout) :: lmask(:,:,:)
        type(nu_envmask_stats),  intent(inout) :: stats
        class(image), optional, target, intent(in) :: vol_dens
        type(nu_envmask_stats) :: fresh
        real, allocatable :: margin(:), score(:)
        stats           = fresh
        stats%n_support = n_nu_mask
        stats%beta_used = max(0., p%beta)
        stats%nsigma    = p%nsigma
        stats%lp_smooth = p%lp_smooth
        stats%l_relative = p%l_relative
        call calc_nu_evidence_margin(margin, p%lp_smooth, p%l_relative)
        call calc_nu_evidence_score(margin, p%nsigma, score, stats)
        if( present(vol_dens) ) call add_nu_evidence_density(vol_dens, p%dens_weight, score, stats)
        call segment_nu_evidence(score, p, lmask, stats)
        if( stats%n_support > 0 )then
            stats%pct_seed   = 100. * real(stats%n_seed)   / real(stats%n_support)
            stats%pct_signal = 100. * real(stats%n_signal) / real(stats%n_support)
        endif
        if( allocated(margin) ) deallocate(margin)
        if( allocated(score)  ) deallocate(score)
    end subroutine nu_evidence_envelope

    module subroutine write_nu_evidence_map( fname, lp_smooth, l_relative )
        class(string),     intent(in) :: fname
        real,    optional, intent(in) :: lp_smooth
        logical, optional, intent(in) :: l_relative
        type(image) :: evmap
        real(kind=c_float), pointer :: rmat(:,:,:) => null()
        real, allocatable :: margin(:)
        logical :: l_rel
        integer :: imask, i, j, k
        l_rel = .false.
        if( present(l_relative) ) l_rel = l_relative
        call calc_nu_evidence_margin(margin, lp_smooth, l_relative)
        call evmap%new(ldim, smpd, wthreads=.false.)
        call evmap%get_rmat_ptr(rmat)
        rmat(:ldim(1),:ldim(2),:ldim(3)) = 0.
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            rmat(i,j,k) = margin(imask)
        end do
        !$omp end parallel do
        call evmap%write(fname, del_if_exists=.true.)
        write(logfhandle,'(A,A)') '>>> WROTE NU EVIDENCE MARGIN MAP: ', fname%to_char()
        if( l_rel )then
            write(logfhandle,'(A)') '    Values are dimensionless baseline-to-best Huber cost-improvement ratios.'
        else
            write(logfhandle,'(A)') '    Units are normalized Huber objective improvement over the coarsest candidate.'
        endif
        call evmap%kill
        deallocate(margin)
    end subroutine write_nu_evidence_map

    module subroutine print_nu_envmask_stats( stats )
        type(nu_envmask_stats), intent(in) :: stats
        write(logfhandle,'(A)')            '>>> NU EVIDENCE ENVELOPE MASK'
        write(logfhandle,'(A,I12)')        '    Support voxels             : ', stats%n_support
        write(logfhandle,'(A,ES12.4)')     '    Null median margin         : ', stats%null_med
        write(logfhandle,'(A,ES12.4)')     '    Null MAD (Gaussian-scaled) : ', stats%null_mad
        write(logfhandle,'(A,F12.3)')      '    Envelope scale (A)         : ', stats%lp_smooth
        write(logfhandle,'(A,L12)')        '    Scale-free margin          : ', stats%l_relative
        write(logfhandle,'(A,F12.3)')      '    Threshold (n MADs)         : ', stats%nsigma
        write(logfhandle,'(A,ES12.4)')     '    Threshold margin           : ', stats%thres
        if( abs(stats%dens_weight) > TINY )then
            write(logfhandle,'(A,F12.3)')  '    Density term weight        : ', stats%dens_weight
            write(logfhandle,'(A,ES12.4)') '    Density median             : ', stats%dens_med
            write(logfhandle,'(A,ES12.4)') '    Density MAD                : ', stats%dens_mad
        endif
        write(logfhandle,'(A,F12.3)')      '    ICM beta                   : ', stats%beta_used
        write(logfhandle,'(A,I12)')        '    ICM iterations             : ', stats%nits
        write(logfhandle,'(A,I12,F9.2,A)') '    Seed voxels (raw threshold): ', stats%n_seed,   stats%pct_seed,   '%'
        write(logfhandle,'(A,I12,F9.2,A)') '    Signal voxels (after ICM)  : ', stats%n_signal, stats%pct_signal, '%'
        if( stats%pct_signal > 50. )then
            write(logfhandle,'(A)') '    WARNING: signal occupies more than half the support, so the median/MAD null'
            write(logfhandle,'(A)') '             estimate is not trustworthy. Widen mskdiam or raise nu_msk_sig.'
        endif
    end subroutine print_nu_envmask_stats

end submodule simple_nu_filter_envmask
