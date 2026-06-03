!@descr: simple nu filter potts implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_potts
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine refine_nu_candidate_map_ordered_labels( candmap, n_candidates, histogram_potts )
        integer(kind=NU_LABEL_KIND), intent(inout) :: candmap(:,:,:)
        integer, intent(in)    :: n_candidates
        logical, optional, intent(in) :: histogram_potts
        integer :: iter, color, i, j, k, imask, icand, cur_icand, best_icand, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: nchanged
        integer :: n_base
        real    :: beta, e, best_e, site_energy
        logical :: l_histogram_potts
        if( n_candidates < 2 ) return
        if( .not. allocated(candidate_coords) ) THROW_HARD('candidate_coords not allocated; refine_nu_candidate_map_ordered_labels')
        if( size(candidate_coords) /= n_candidates ) &
            &THROW_HARD('candidate_coords size mismatch; refine_nu_candidate_map_ordered_labels')
        l_histogram_potts = .false.
        if( present(histogram_potts) ) l_histogram_potts = histogram_potts
        n_base = size(cutoff_finds)
        beta = estimate_nu_label_smooth_beta(n_candidates)
        write(logfhandle,'(A,ES12.4,A,I0,A,I0,A,I0,A,I0)') '>>> NU ordered-label smoothing: beta=', beta, &
            &', max iterations=', NU_LABEL_SMOOTH_MAXITS, ', candidates=', n_candidates, &
            &', base labels=', n_base, ', step tolerance=', NU_LABEL_SMOOTH_STEP_TOL
        write(logfhandle,'(A,I0,A,I0)') '>>> NU ordered-label smoothing neighborhood: ', &
            &NU_LABEL_SMOOTH_NNEIGH, '-connected, color passes=', NU_LABEL_SMOOTH_NCOLORS
        write(logfhandle,'(A,F6.3)') '>>> NU ordered-label smoothing quadratic jump fraction: ', &
            &NU_LABEL_SMOOTH_QUAD_FRAC
        call log_nu_candidate_coords
        if( beta <= TINY )then
            write(logfhandle,'(A)') '>>> NU ordered-label smoothing skipped: beta <= TINY'
            return
        endif
        site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
        write(logfhandle,'(A,F12.5)') '>>> NU ordered-label smoothing initial mean site energy: ', site_energy
        if( l_histogram_potts )then
            call refine_nu_candidate_map_histogram_ordered_labels(candmap, n_candidates, beta)
            return
        endif
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do schedule(static) default(shared) &
                !$omp private(i,j,k,imask,icand,cur_icand,best_icand,n_full,nsz,e,best_e) &
                !$omp reduction(+:nchanged) proc_bind(close)
                do imask = 1, n_nu_mask
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( nu_label_smooth_color(i,j,k) /= color ) cycle
                    call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                    cur_icand  = int(candmap(i,j,k))
                    best_icand = cur_icand
                    best_e     = dmats_mask(imask,cur_icand) + beta * &
                        &nu_label_smooth_neighborhood_cost(cur_icand, candmap, n_full, nsz)
                    do icand = 1, n_candidates
                        if( icand == cur_icand ) cycle
                        e = dmats_mask(imask,icand) + beta * &
                            &nu_label_smooth_neighborhood_cost(icand, candmap, n_full, nsz)
                        if( nu_label_smooth_is_better(e, best_e) )then
                            best_e     = e
                            best_icand = icand
                        endif
                    end do
                    if( best_icand /= cur_icand )then
                        nchanged = nchanged + 1
                        candmap(i,j,k) = int(best_icand, kind=NU_LABEL_KIND)
                    endif
                end do
                !$omp end parallel do
            end do
            site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
            write(logfhandle,'(A,I0,A,I0,A,F12.5)') '>>> NU ordered-label smoothing iteration ', iter, &
                &' changed voxels: ', nchanged, ', mean site energy: ', site_energy
            if( nchanged == 0 ) exit
        end do
    end subroutine refine_nu_candidate_map_ordered_labels

    subroutine refine_nu_candidate_map_histogram_ordered_labels( candmap, n_candidates, beta )
        integer(kind=NU_LABEL_KIND), intent(inout) :: candmap(:,:,:)
        integer, intent(in) :: n_candidates
        real,    intent(in) :: beta
        integer, allocatable :: counts(:), target_counts(:), slack_counts(:)
        integer :: iter, color, i, j, k, imask, icand, cur_icand, best_icand
        integer :: n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nchanged, label_drift
        real    :: e, best_e, site_energy, hist_scale, hist_energy
        allocate(counts(n_candidates),        source=0)
        allocate(target_counts(n_candidates), source=0)
        allocate(slack_counts(n_candidates),  source=0)
        call count_nu_candidate_labels(candmap, n_candidates, target_counts)
        counts = target_counts
        slack_counts = max(1, nint(NU_LABEL_HIST_SLACK_FRAC * real(max(1, n_nu_mask))))
        hist_scale = NU_LABEL_HIST_BETA_FRAC * beta / real(max(1, n_nu_mask))
        hist_energy = calc_nu_label_histogram_energy(counts, target_counts, slack_counts, hist_scale)
        write(logfhandle,'(A,F6.3,A,ES12.4,A,F6.3,A,I0)') &
            &'>>> NU slack-histogram ordered-label smoothing: beta fraction=', &
            &NU_LABEL_HIST_BETA_FRAC, ', scale=', hist_scale, ', slack fraction=', &
            &NU_LABEL_HIST_SLACK_FRAC, ', slack voxels/bin=', slack_counts(1)
        write(logfhandle,'(A,ES12.4)') &
            &'>>> NU slack-histogram ordered-label smoothing initial mean histogram energy: ', &
            &hist_energy / real(max(1, n_nu_mask))
        call log_nu_label_histogram_counts('target', target_counts, n_candidates)
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do schedule(static) default(shared) &
                !$omp private(i,j,k,imask,icand,cur_icand,best_icand,n_full,nsz,e,best_e) &
                !$omp reduction(+:nchanged) proc_bind(close)
                do imask = 1, n_nu_mask
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( nu_label_smooth_color(i,j,k) /= color ) cycle
                    call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                    cur_icand  = int(candmap(i,j,k))
                    best_icand = cur_icand
                    best_e     = dmats_mask(imask,cur_icand) + beta * &
                        &nu_label_smooth_neighborhood_cost(cur_icand, candmap, n_full, nsz)
                    do icand = 1, n_candidates
                        if( icand == cur_icand ) cycle
                        e = dmats_mask(imask,icand) + beta * &
                            &nu_label_smooth_neighborhood_cost(icand, candmap, n_full, nsz) + &
                            &nu_label_histogram_delta(icand, cur_icand, counts, target_counts, &
                            &slack_counts, hist_scale)
                        if( nu_label_smooth_is_better(e, best_e) )then
                            best_e     = e
                            best_icand = icand
                        endif
                    end do
                    if( best_icand /= cur_icand )then
                        nchanged = nchanged + 1
                        candmap(i,j,k) = int(best_icand, kind=NU_LABEL_KIND)
                    endif
                end do
                !$omp end parallel do
                call count_nu_candidate_labels(candmap, n_candidates, counts)
            end do
            site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
            hist_energy = calc_nu_label_histogram_energy(counts, target_counts, slack_counts, hist_scale)
            write(logfhandle,'(A,I0,A,I0,A,F12.5,A,ES12.4)') &
                &'>>> NU histogram-constrained ordered-label smoothing iteration ', iter, &
                &' changed voxels: ', nchanged, ', mean site energy: ', site_energy, &
                &', mean histogram energy: ', hist_energy / real(max(1, n_nu_mask))
            if( nchanged == 0 ) exit
        end do
        label_drift = sum(abs(counts - target_counts))
        hist_energy = calc_nu_label_histogram_energy(counts, target_counts, slack_counts, hist_scale)
        write(logfhandle,'(A,I0,A,ES12.4)') &
            &'>>> NU histogram-constrained ordered-label smoothing final label-count drift: ', &
            &label_drift, ', mean histogram energy: ', hist_energy / real(max(1, n_nu_mask))
        call log_nu_label_histogram_counts('final', counts, n_candidates)
        deallocate(counts, target_counts, slack_counts)
    end subroutine refine_nu_candidate_map_histogram_ordered_labels

    module real function estimate_nu_label_smooth_beta( n_candidates )
        integer, intent(in) :: n_candidates
        integer :: imask, icand, nvox
        real    :: best_e, second_e, cur_e, beta_sum
        estimate_nu_label_smooth_beta = 0.
        beta_sum = 0.
        nvox = 0
        if( n_candidates < 2 ) return
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(imask,icand,best_e,second_e,cur_e) reduction(+:beta_sum,nvox) proc_bind(close)
        do imask = 1, n_nu_mask
            best_e   = huge(best_e)
            second_e = huge(second_e)
            do icand = 1, n_candidates
                cur_e = dmats_mask(imask,icand)
                if( cur_e < best_e )then
                    second_e = best_e
                    best_e   = cur_e
                else if( cur_e < second_e )then
                    second_e = cur_e
                endif
            end do
            if( second_e < huge(second_e) )then
                beta_sum = beta_sum + max(0., second_e - best_e)
                nvox = nvox + 1
            endif
        end do
        !$omp end parallel do
        if( nvox > 0 ) estimate_nu_label_smooth_beta = &
            &NU_LABEL_SMOOTH_BETA_FRAC * beta_sum / real(nvox)
    end function estimate_nu_label_smooth_beta

    module real function nu_label_smooth_neighborhood_cost( icand, candmap, neigh, nsz )
        integer, intent(in) :: icand, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer :: ineigh, ni, nj, nk, degree
        nu_label_smooth_neighborhood_cost = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            nu_label_smooth_neighborhood_cost = nu_label_smooth_neighborhood_cost + &
                &nu_label_smooth_pair_cost(icand, int(candmap(ni,nj,nk)))
        end do
        if( degree > 0 ) nu_label_smooth_neighborhood_cost = nu_label_smooth_neighborhood_cost / real(degree)
    end function nu_label_smooth_neighborhood_cost

    subroutine count_nu_candidate_labels( candmap, n_candidates, counts )
        integer(kind=NU_LABEL_KIND), intent(in)  :: candmap(:,:,:)
        integer,                     intent(in)  :: n_candidates
        integer,                     intent(out) :: counts(:)
        integer :: i, j, k, imask, ilabel, nbad
        if( size(counts) /= n_candidates ) THROW_HARD('counts size mismatch; count_nu_candidate_labels')
        counts = 0
        nbad = 0
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,ilabel) &
        !$omp reduction(+:counts,nbad) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            ilabel = int(candmap(i,j,k))
            if( ilabel < 1 .or. ilabel > n_candidates )then
                nbad = nbad + 1
            else
                counts(ilabel) = counts(ilabel) + 1
            endif
        end do
        !$omp end parallel do
        if( nbad > 0 ) THROW_HARD('candidate label out of range; count_nu_candidate_labels')
    end subroutine count_nu_candidate_labels

    real function nu_label_histogram_delta( icand, cur_icand, counts, target_counts, slack_counts, hist_scale )
        integer, intent(in) :: icand, cur_icand, counts(:), target_counts(:), slack_counts(:)
        real,    intent(in) :: hist_scale
        real :: cur_before, cand_before
        if( icand == cur_icand )then
            nu_label_histogram_delta = 0.
            return
        endif
        cur_before  = real(counts(cur_icand) - target_counts(cur_icand))
        cand_before = real(counts(icand)     - target_counts(icand))
        nu_label_histogram_delta = hist_scale * (nu_label_histogram_bin_penalty(cur_before - 1., slack_counts(cur_icand)) - &
            &nu_label_histogram_bin_penalty(cur_before, slack_counts(cur_icand)) + &
            &nu_label_histogram_bin_penalty(cand_before + 1., slack_counts(icand)) - &
            &nu_label_histogram_bin_penalty(cand_before, slack_counts(icand)))
    end function nu_label_histogram_delta

    real function calc_nu_label_histogram_energy( counts, target_counts, slack_counts, hist_scale )
        integer, intent(in) :: counts(:), target_counts(:), slack_counts(:)
        real,    intent(in) :: hist_scale
        integer :: ilabel
        if( size(counts) /= size(target_counts) ) THROW_HARD('histogram size mismatch; calc_nu_label_histogram_energy')
        if( size(counts) /= size(slack_counts)  ) THROW_HARD('slack size mismatch; calc_nu_label_histogram_energy')
        calc_nu_label_histogram_energy = 0.
        do ilabel = 1, size(counts)
            calc_nu_label_histogram_energy = calc_nu_label_histogram_energy + &
                &nu_label_histogram_bin_penalty(real(counts(ilabel) - target_counts(ilabel)), slack_counts(ilabel))
        end do
        calc_nu_label_histogram_energy = hist_scale * calc_nu_label_histogram_energy
    end function calc_nu_label_histogram_energy

    real function nu_label_histogram_bin_penalty( drift, slack_count )
        real,    intent(in) :: drift
        integer, intent(in) :: slack_count
        real :: excess
        excess = max(0., abs(drift) - real(slack_count))
        nu_label_histogram_bin_penalty = excess * excess
    end function nu_label_histogram_bin_penalty

    subroutine log_nu_label_histogram_counts( stage, counts, n_candidates )
        character(len=*), intent(in) :: stage
        integer,          intent(in) :: counts(:), n_candidates
        integer :: ilabel
        if( size(counts) /= n_candidates ) THROW_HARD('counts size mismatch; log_nu_label_histogram_counts')
        write(logfhandle,'(A,1X,A,A)') '>>> NU histogram-constrained ordered-label smoothing ', trim(stage), ' counts'
        do ilabel = 1, n_candidates
            write(logfhandle,'(A,I0,A,F8.3,A,I0)') '    label ', ilabel, ' (', &
                &nu_label_lowpass_limit(ilabel), ' A): ', counts(ilabel)
        end do
    end subroutine log_nu_label_histogram_counts

    module real function nu_label_smooth_pair_cost( icand, jcand )
        integer, intent(in) :: icand, jcand
        nu_label_smooth_pair_cost = nu_label_smooth_coord_pair_cost(candidate_coords(icand), &
            &candidate_coords(jcand))
    end function nu_label_smooth_pair_cost

    module real function nu_label_smooth_coord_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
        real :: step_jump
        if( abs(icoord - jcoord) <= TINY )then
            nu_label_smooth_coord_pair_cost = 0.
        else
            step_jump = max(0., abs(icoord - jcoord) - real(NU_LABEL_SMOOTH_STEP_TOL))
            nu_label_smooth_coord_pair_cost = step_jump + NU_LABEL_SMOOTH_QUAD_FRAC * step_jump * step_jump
        endif
    end function nu_label_smooth_coord_pair_cost

    module integer function nu_label_smooth_color( i, j, k )
        integer, intent(in) :: i, j, k
        nu_label_smooth_color = mod(i,2) + 2 * mod(j,2) + 4 * mod(k,2)
    end function nu_label_smooth_color

    module logical function nu_label_smooth_is_better( e, best_e )
        real, intent(in) :: e, best_e
        real :: tol
        tol = NU_LABEL_SMOOTH_TIE_EPS * max(1., abs(best_e))
        nu_label_smooth_is_better = e < best_e - tol
    end function nu_label_smooth_is_better

    module real function calc_nu_label_smooth_site_energy( candmap, beta )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        real,    intent(in) :: beta
        integer :: i, j, k, imask, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nvox
        real :: energy_sum
        calc_nu_label_smooth_site_energy = 0.
        energy_sum = 0.
        nvox = 0
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(i,j,k,imask,n_full,nsz) reduction(+:energy_sum,nvox) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
            energy_sum = energy_sum + dmats_mask(imask,int(candmap(i,j,k))) + beta * &
                &nu_label_smooth_neighborhood_cost(int(candmap(i,j,k)), candmap, n_full, nsz)
            nvox = nvox + 1
        end do
        !$omp end parallel do
        if( nvox > 0 ) calc_nu_label_smooth_site_energy = energy_sum / real(nvox)
    end function calc_nu_label_smooth_site_energy

end submodule simple_nu_filter_potts
