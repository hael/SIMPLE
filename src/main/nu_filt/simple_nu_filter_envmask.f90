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
! not enforce a connected result. Topology is applied by the connected component
! and morphology tail in simple_image_msk: piecewise callers invoke it themselves
! after nu_evidence_envelope, while write_nu_evidence_envmask runs the complete
! evidence -> segmentation -> topology -> artifact chain as the one shared
! producer for every workflow that regenerates the envelope from live evidence.
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

    !> The packed voxel set on which evidence labels are free: the observed
    !! density envelope when a Euclidean null shell is set (voxels outside it
    !! are fixed solvent), else the observed support, else all of it. Never
    !! the raw support once a pair carries exact zero/zero voxels: those are
    !! a degenerate spike at margin 0 that pins the median and collapses the
    !! MAD.
    subroutine nu_evidence_calibration_mask( calib )
        logical, allocatable, intent(inout) :: calib(:)
        if( allocated(calib) ) deallocate(calib)
        if( allocated(nu_calib_lmask) )then
            if( size(nu_calib_lmask) /= n_nu_mask ) THROW_HARD('NU evidence label domain size mismatch')
            allocate(calib(n_nu_mask), source=nu_calib_lmask)
        else if( allocated(nu_observed_mask) )then
            if( size(nu_observed_mask) /= n_nu_mask ) THROW_HARD('NU observation mask size mismatch')
            allocate(calib(n_nu_mask), source=nu_observed_mask)
        else
            allocate(calib(n_nu_mask), source=.true.)
        endif
        if( count(calib) < 1 ) THROW_HARD('empty NU evidence label domain')
    end subroutine nu_evidence_calibration_mask

    !> The packed voxel set the null statistics (margin median/MAD, density
    !! median/MAD) are estimated on: the Euclidean null shell when one is set
    !! (envelope-constrained base pair), else the label domain itself, where
    !! the robust statistics resolve the solvent-majority mixture (spherical
    !! base pair). l_shell reports which regime applies.
    subroutine nu_evidence_null_mask( nullm, l_shell )
        logical, allocatable, intent(inout) :: nullm(:)
        logical,              intent(out)   :: l_shell
        if( allocated(nullm) ) deallocate(nullm)
        l_shell = allocated(nu_null_lmask)
        if( l_shell )then
            if( size(nu_null_lmask) /= n_nu_mask ) THROW_HARD('NU evidence null shell size mismatch')
            allocate(nullm(n_nu_mask), source=nu_null_lmask)
        else
            call nu_evidence_calibration_mask(nullm)
        endif
    end subroutine nu_evidence_null_mask

    !>  Robust floor for the relative cost ratio, so that near-zero baseline or
    !!  best-candidate costs cannot manufacture arbitrarily large evidence.
    !!  Estimated on the calibration domain.
    module real function nu_evidence_baseline_floor( base_full ) result( floor_val )
        real, intent(in) :: base_full(:,:,:)
        real,    allocatable :: work(:)
        logical, allocatable :: calib(:)
        integer :: imask, i, j, k, n
        real    :: val
        call nu_evidence_calibration_mask(calib)
        allocate(work(n_nu_mask), source=0.)
        n = 0
        do imask = 1, n_nu_mask
            if( .not.calib(imask) ) cycle
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
        deallocate(work, calib)
    end function nu_evidence_baseline_floor

    module subroutine calc_nu_evidence_score( margin, nsigma, score, stats )
        real,                    intent(in)    :: margin(:)
        real,                    intent(in)    :: nsigma
        real, allocatable,       intent(inout) :: score(:)
        type(nu_envmask_stats),  intent(inout) :: stats
        real,    allocatable :: work(:), margin_null(:)
        logical, allocatable :: calib(:), nullm(:)
        real    :: med, mad_val, denom
        integer :: n
        logical :: l_shell
        n = size(margin)
        if( n < 1 ) THROW_HARD('empty margin vector; calc_nu_evidence_score')
        if( n /= n_nu_mask ) THROW_HARD('margin size mismatch; calc_nu_evidence_score')
        ! Two regimes. Spherical base pair: the null is the robust median/MAD
        ! of the margin over the observed support, a solvent-majority mixture
        ! by construction of the generous sphere; l_null_majority reports
        ! whether that held. Envelope-constrained base pair: the null is
        ! designated by Euclidean geometry (the dilation ring of the density
        ! envelope) and the median/MAD are taken there, where residual weak
        ! density is the only contamination and the robust pair absorbs it.
        call nu_evidence_calibration_mask(calib)
        call nu_evidence_null_mask(nullm, l_shell)
        stats%n_calib      = count(calib)
        stats%n_null       = count(nullm)
        stats%l_null_shell = l_shell
        if( stats%n_null < 1 )then
            ! an empty null set (a shell the base support carries nowhere at
            ! full weight) cannot be calibrated: an unattainable threshold
            ! yields an empty envelope, and the caller's fallback applies
            stats%null_med = 0.
            stats%null_mad = 0.
            stats%thres    = huge(1.)
            if( allocated(score) ) deallocate(score)
            allocate(score(n), source=NU_ENVMASK_EXCLUDED_SCORE)
            deallocate(calib, nullm)
            return
        endif
        margin_null = pack(margin, nullm)
        allocate(work(size(margin_null)), source=margin_null)
        med     = median_nocopy(work)
        mad_val = mad_gau(margin_null, med)
        deallocate(work, margin_null)
        denom = max(mad_val, TINY)
        stats%null_med = med
        stats%null_mad = mad_val
        stats%thres    = med + nsigma * denom
        if( l_shell )then
            ! separation diagnostic: the median margin inside the core
            ! (domain minus shell) against the shell's null median
            work = pack(margin, calib .and. .not.nullm)
            if( size(work) > 0 ) stats%core_med = median_nocopy(work)
            deallocate(work)
        endif
        if( allocated(score) ) deallocate(score)
        allocate(score(n), source=0.)
        ! outside the domain the label is solvent by construction, so the
        ! evidence envelope is nested inside the density envelope
        where( calib )
            score = (margin - stats%thres) / denom
        elsewhere
            score = NU_ENVMASK_EXCLUDED_SCORE
        end where
        deallocate(calib, nullm)
    end subroutine calc_nu_evidence_score

    module subroutine add_nu_evidence_density( vol_dens, weight, score, stats )
        class(image), target,    intent(in)    :: vol_dens
        real,                    intent(in)    :: weight
        real, allocatable,       intent(inout) :: score(:)
        type(nu_envmask_stats),  intent(inout) :: stats
        real(kind=c_float), pointer :: rmat(:,:,:) => null()
        real,    allocatable :: dens(:), work(:), dens_null(:)
        logical, allocatable :: calib(:), nullm(:)
        real    :: med, mad_val, denom
        integer :: imask, i, j, k
        logical :: l_shell
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
        call nu_evidence_calibration_mask(calib)
        call nu_evidence_null_mask(nullm, l_shell)
        if( count(nullm) < 1 )then
            deallocate(dens, calib, nullm)
            return
        endif
        dens_null = pack(dens, nullm)
        allocate(work(size(dens_null)), source=dens_null)
        med     = median_nocopy(work)
        mad_val = mad_gau(dens_null, med)
        deallocate(work, dens_null, nullm)
        denom = max(mad_val, TINY)
        stats%dens_med    = med
        stats%dens_mad    = mad_val
        stats%dens_weight = weight
        ! Additive so that strong density can hold in a poorly ordered region that
        ! the resolution evidence alone would carve out. This is the term that
        ! protects flexible periphery, which is otherwise indistinguishable from
        ! solvent by cross-half consistency. Voxels outside the calibration
        ! domain keep their fixed solvent score.
        where( calib ) score = score + weight * (dens - med) / denom
        deallocate(dens, calib)
    end subroutine add_nu_evidence_density

    module subroutine segment_nu_evidence( score, p, lmask, stats )
        real,                    intent(in)    :: score(:)
        type(nu_envmask_params), intent(in)    :: p
        logical, allocatable,    intent(inout) :: lmask(:,:,:)
        type(nu_envmask_stats),  intent(inout) :: stats
        logical, allocatable :: lab(:,:,:)
        integer :: iter, color, imask, i, j, k, ineigh, ni, nj, nk, nsig
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
                if( NU_DEV_OUTPUT .and. nu_l_report ) write(logfhandle,'(A,I2,A,I10)') &
                    &'>>> NU envelope ICM iteration ', iter, ' changed voxels: ', nchanged
                if( nchanged == 0 ) exit
            end do
        endif
        stats%n_signal = count(lab)
        ! how much of the Euclidean null shell (the density envelope's dilation
        ! ring) the evidence labels signal: the empirical answer to whether the
        ! dilation is capturing density or is pure margin
        stats%pct_signal_null = 0.
        if( allocated(nu_null_lmask) )then
            if( size(nu_null_lmask) == n_nu_mask .and. stats%n_null > 0 )then
                nsig = 0
                !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) reduction(+:nsig) proc_bind(close)
                do imask = 1, n_nu_mask
                    if( .not.nu_null_lmask(imask) ) cycle
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( lab(i,j,k) ) nsig = nsig + 1
                end do
                !$omp end parallel do
                stats%pct_signal_null = 100. * real(nsig) / real(stats%n_null)
            endif
        endif
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
            stats%pct_calib  = 100. * real(stats%n_calib)  / real(stats%n_support)
        endif
        ! validity of the null, per regime: a mixture null (spherical base
        ! pair) is only meaningful if solvent held the majority of the domain;
        ! a Euclidean shell (envelope-constrained base pair) needs enough
        ! voxels for its median/MAD to be stable. Signal cannot exceed the
        ! domain since everything outside it is fixed solvent.
        if( stats%n_calib > 0 )then
            stats%pct_signal_calib = 100. * real(stats%n_signal) / real(stats%n_calib)
            stats%pct_null         = 100. * real(stats%n_null)   / real(stats%n_calib)
        endif
        stats%l_null_majority = stats%pct_signal_calib <= 50.
        if( stats%l_null_shell )then
            stats%l_null_valid = stats%n_null >= NU_ENVMASK_MIN_NULL_VOX .and. &
                &real(stats%n_null) >= NU_ENVMASK_MIN_NULL_FRAC * real(max(1,stats%n_calib))
        else
            stats%l_null_valid = stats%l_null_majority
        endif
        if( allocated(margin) ) deallocate(margin)
        if( allocated(score)  ) deallocate(score)
    end subroutine nu_evidence_envelope

    !> The one shared NU-evidence envelope producer: build the envelope from
    !! the LIVE raw evidence (setup_nu_dmats plus a completed candidate
    !! evaluation must have run), apply the connected-component/morphology
    !! topology tail, and write the artifact to the explicit filename. The
    !! filename is passed in by the caller, which regenerates the envelope
    !! every cycle it runs the competition under automsk=yes.
    !! An empty evidence field warns and writes nothing; every envelope
    !! consumer handles absence.
    module subroutine write_nu_evidence_envmask( nsigma, lp_smooth, smpd, state, fname, l_arm_background, l_armed )
        use simple_image_msk, only: image_msk
        real,              intent(in)  :: nsigma, lp_smooth, smpd
        integer,           intent(in)  :: state
        class(string),     intent(in)  :: fname
        logical, optional, intent(in)  :: l_arm_background
        !! l_armed: armed from the evidence envelope (non-empty, valid null)
        logical, optional, intent(out) :: l_armed
        type(nu_envmask_params) :: envp
        type(nu_envmask_stats)  :: envstats
        type(image_msk)         :: envmask
        logical, allocatable    :: l_env(:,:,:)
        integer :: grow_px, edge_px, n_ccs, n_ccs_kept
        logical :: l_arm
        l_arm = .false.
        if( present(l_arm_background) ) l_arm = l_arm_background
        if( present(l_armed) ) l_armed = .false.
        envp%nsigma      = nsigma
        envp%beta        = NU_ENVMASK_BETA
        envp%dens_weight = NU_ENVMASK_DENS_WEIGHT
        envp%lp_smooth   = lp_smooth
        envp%l_relative  = NU_ENVMASK_RELATIVE
        call nu_evidence_envelope(envp, l_env, envstats)
        call print_nu_envmask_stats(envstats)
        grow_px = max(1, nint(NU_ENVMASK_GROW_A / smpd))
        edge_px = max(1, nint(NU_ENVMASK_EDGE_A / smpd))
        call envmask%envmask3D_from_lmask(l_env, smpd, grow_px, edge_px, &
            &NU_ENVMASK_MINVOL_FRAC, .true., n_ccs, n_ccs_kept)
        ! an empty evidence field returns without constructing the mask image;
        ! writing it would abort on invalid MRC dimensions. Skip the write --
        ! every envelope consumer handles absence.
        if( n_ccs_kept < 1 )then
            THROW_WARN('NU evidence envelope is empty; no envelope mask written this iteration')
        else
            call envmask%write(fname, del_if_exists=.true.)
            call wait_for_closure(fname)
            write(logfhandle,'(A,I0,A,1X,A)') &
                &'>>> NU EVIDENCE ENVELOPE: STATE ', state, ', MASK', fname%to_char()
            if( .not. envstats%l_null_valid )then
                if( envstats%l_null_shell )then
                    write(logfhandle,'(A,I0,A,I0,A,F5.1,A)') '>>> NU EVIDENCE ENVELOPE: STATE ', state, &
                        &', Euclidean null shell too thin (', envstats%n_null, ' voxels, ', envstats%pct_null, &
                        &' % of the domain); widen the density envelope dilation (binwidth); the envelope is not armed'
                else
                    write(logfhandle,'(A,I0,A,F6.1,A)') '>>> NU EVIDENCE ENVELOPE: STATE ', state, &
                        &', signal occupies ', envstats%pct_signal_calib, &
                        &' % of the support; the median/MAD null is not trustworthy and the envelope is not armed'
                endif
            endif
            ! automsk=yes background policy: the filter-field background is the
            ! complement of this envelope, derived from the SAME evidence pass
            ! (no second compute). Voxels outside it take the coarsest bank
            ! candidate -- a heavy background low-pass (cisTEM-style) that
            ! down-weights the excluded density's contribution to alignment
            ! without removing it from the reference. The PCG SOLVE support
            ! stays on the conservative density envelope (automsk=yes only),
            ! never on this evidence mask. Armed only on a valid null;
            ! the caller owns the fallback to the density envelope.
            if( l_arm .and. envstats%l_null_valid )then
                call set_nu_solvent_envelope(envmask, source='nu_evidence_envelope')
                if( present(l_armed) ) l_armed = .true.
                write(logfhandle,'(A,I0)') &
                    &'>>> NU BACKGROUND: FILTER-FIELD BACKGROUND ARMED FROM THE EVIDENCE ENVELOPE, STATE ', state
            endif
        endif
        ! One greppable line per state per cycle. The envelope is allowed to
        ! shrink as resolution improves, but reference masking suppresses
        ! whatever it excludes, so a monotonically falling occupancy is the
        ! signal that it is clipping rather than tightening.
        write(logfhandle,'(A,I0,A,F8.3,A,I0,A,I0)') &
            &'>>> NU ENVELOPE OCCUPANCY: STATE ', state, ', SUPPORT FRACTION ', &
            &envstats%pct_signal, ' %, COMPONENTS KEPT ', n_ccs_kept, ' OF ', n_ccs
        ! the dilation-ring occupancy says whether the density envelope's
        ! dilation is capturing density (signal in the ring) or is pure
        ! margin (ring null throughout): the number to consult before
        ! tightening binwidth / ENVMSKWIDTH_A_MIN
        if( envstats%l_null_shell ) write(logfhandle,'(A,I0,A,F8.3,A)') &
            &'>>> NU DILATION RING OCCUPANCY: STATE ', state, ', SIGNAL FRACTION OF THE RING ', &
            &envstats%pct_signal_null, ' %'
        call envmask%kill_bimg
        if( allocated(l_env) ) deallocate(l_env)
    end subroutine write_nu_evidence_envmask

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
        if( .not. (NU_DEV_OUTPUT .and. nu_l_report) )then
            ! the validity warning must never be silenced
            call warn_if_invalid
            return
        endif
        write(logfhandle,'(A)')            '>>> NU EVIDENCE ENVELOPE MASK'
        write(logfhandle,'(A,I12)')        '    Support voxels             : ', stats%n_support
        write(logfhandle,'(A,I12,F9.2,A)') '    Label domain voxels        : ', stats%n_calib, stats%pct_calib, '%'
        if( stats%l_null_shell )then
            write(logfhandle,'(A)')        '    Null model                 :   Euclidean shell (dilation ring)'
            write(logfhandle,'(A,I12,F9.2,A)') '    Null shell voxels          : ', stats%n_null, stats%pct_null, '%'
            write(logfhandle,'(A,ES12.4)') '    Core median margin         : ', stats%core_med
            write(logfhandle,'(A,F12.2,A)') '    Signal within null shell   : ', stats%pct_signal_null, '%'
        else
            write(logfhandle,'(A)')        '    Null model                 :   robust mixture over the support'
        endif
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
        write(logfhandle,'(A,F12.2,A)')    '    Signal within domain       : ', stats%pct_signal_calib, '%'
        call warn_if_invalid

    contains

        subroutine warn_if_invalid
            if( stats%l_null_valid ) return
            if( stats%l_null_shell )then
                write(logfhandle,'(A)') '    WARNING: the Euclidean null shell is too thin for a stable median/MAD;'
                write(logfhandle,'(A)') '             widen the density envelope dilation (binwidth).'
            else
                write(logfhandle,'(A)') '    WARNING: signal occupies more than half the support, so the median/MAD null'
                write(logfhandle,'(A)') '             estimate is not trustworthy. Widen mskdiam or raise nu_msk_sig.'
            endif
        end subroutine warn_if_invalid

    end subroutine print_nu_envmask_stats

end submodule simple_nu_filter_envmask
