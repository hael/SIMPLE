!@descr: simple nu filter extend implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_extend
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, new_limit, stats, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: new_limit        ! Angstrom limit for the proposed shell
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
        type(image)       :: vol_even_filt_new, vol_odd_filt_new
        type(string)      :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_new(:,:,:), dmat_tmp(:,:,:), dmat_finest_mask(:)
        integer, allocatable :: frontier_vox(:)
        integer           :: new_find, n_finest, n_total, n_extended, sz_old, old_label
        integer           :: old_radius_px, new_radius_px, n_seed_min
        real              :: accept_pct_eff, pct_finest, x, noise_sigma
        real              :: old_radius_angstrom, new_radius_angstrom
        integer(kind=NU_LABEL_KIND), allocatable :: extend_choice(:)
        type(nu_highres_extension_stats) :: local_stats
        integer           :: i
        logical           :: l_permissive_accept
        local_stats%new_limit = new_limit
        accept_pct_eff = NU_HIGHRES_EXTENSION_ACCEPT_PCT
        if( present(accept_pct) ) accept_pct_eff = max(0., accept_pct)
        l_permissive_accept = present(accept_pct) .and. accept_pct_eff <= TINY
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats first')
        if( nu_aux_replacement_label > 0 ) &
            &THROW_HARD('NU high-resolution extension cannot run with an auxiliary replacement candidate')
        sz_old    = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        n_total   = n_nu_mask
        local_stats%n_mask = n_total
        if( n_total == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        old_label = 0
        n_finest  = 0
        do i = sz_old, 1, -1
            n_finest = count_nu_base_label_voxels(i)
            if( n_finest > 0 )then
                old_label = i
                exit
            endif
        end do
        if( old_label == 0 )then
            local_stats%n_tested = 0
            if( present(stats) ) stats = local_stats
            return
        endif
        if( old_label < sz_old )then
            write(logfhandle,'(A,I0,A,I0,A,F8.3,A)') &
                &'>>> NU high-resolution extension frontier reset from empty bank label ', &
                &sz_old, ' to populated label ', old_label, ' (', cutoff_find_to_lowpass_limit(old_label), ' A)'
            call prune_nu_highres_extension_bank(old_label)
            sz_old = old_label
        endif
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        local_stats%n_tested = n_finest
        pct_finest = 100. * real(n_finest) / real(n_total)
        local_stats%pct_tested_mask = pct_finest
        if( n_finest == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        if( pct_finest < threshold_pct )then
            if( present(stats) ) stats = local_stats
            return   ! trigger not met, nothing to do
        endif
        call collect_nu_base_label_voxels(sz_old, frontier_vox)
        if( size(frontier_vox) /= n_finest )then
            n_finest = size(frontier_vox)
            local_stats%n_tested = n_finest
            if( n_finest == 0 )then
                if( present(stats) ) stats = local_stats
                return
            endif
        endif
        new_find = calc_fourier_index(new_limit, box, smpd)
        if( new_find <= cutoff_finds(sz_old) )then
            if( present(stats) ) stats = local_stats
            return
        endif
        if( any(cutoff_finds == new_find) )then
            if( present(stats) ) stats = local_stats
            return
        endif
        local_stats%new_find = new_find
        local_stats%attempted = .true.
        ! --- build the new Butterworth filter and cache filtered vols ---
        even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), new_find)
        odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  new_find)
        call delete_cached_filtered_pair(new_find)
        ! Generate just this one new filtered pair after deleting any stale
        ! scratch products left by a prior interrupted run.
        call generate_single_filtered_pair(vol_even, vol_odd, new_find, even_cache_fname, odd_cache_fname)
        ! --- evaluate the new objective only within the mask ---
        allocate(dmat_new(ldim(1),ldim(2),ldim(3)), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)), source=0.)
        call vol_even_filt_new%new(ldim, smpd)
        call vol_odd_filt_new%new(ldim, smpd)
        call vol_even_filt_new%read(even_cache_fname)
        call vol_odd_filt_new%read(odd_cache_fname)
        if( nu_noise_sigma_cached > TINY )then
            noise_sigma = nu_noise_sigma_cached
        else
            noise_sigma = vol_even%nu_objective_noise_scale(vol_odd, nu_lmask)
            nu_noise_sigma_cached = noise_sigma
        endif
        call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_new, &
            &nu_lmask, noise_sigma)
        call vol_even_filt_new%kill
        call vol_odd_filt_new%kill
        call smooth_nu_objective(dmat_new, dmat_tmp, new_limit)
        allocate(dmat_finest_mask(n_nu_mask), source=huge(x))
        if( allocated(dmat_finest_cached) .and. size(dmat_finest_cached) == n_nu_mask ) then
            dmat_finest_mask = dmat_finest_cached
        else
            call fill_nu_frontier_dmat_from_bank(frontier_vox, sz_old, dmat_finest_mask)
        end if
        deallocate(dmat_tmp)
        ! Keep the smoothed mask-normalization volume across adjacent shell
        ! challenges; nu_filter_vols/cleanup releases it before output synthesis.
        ! --- update filtmap in place for the masked voxels ---
        n_extended = 0
        allocate(extend_choice(n_finest), source=0_NU_LABEL_KIND)
        call init_nu_highres_extension_selection(frontier_vox, dmat_finest_mask, dmat_new, &
            &extend_choice, n_extended)
        ! The extension experiment is already constrained to the current
        ! finest frontier and one Fourier-shell step, so we deliberately do not
        ! apply the full-bank Potts prior here. This lets any local unary
        ! evidence seed the next refinement shell.
        local_stats%n_unary_wins = n_extended
        if( n_finest > 0 )then
            local_stats%pct_unary_wins_tested = 100. * real(local_stats%n_unary_wins) / real(n_finest)
        endif
        if( n_total > 0 )then
            local_stats%pct_unary_wins_mask = 100. * real(local_stats%n_unary_wins) / real(n_total)
        endif
        n_seed_min = min(n_total, NU_HIGHRES_EXTENSION_MIN_SEED_VOXELS)
        local_stats%n_seed_min = n_seed_min
        write(logfhandle,'(A,F8.3,A,I0,A,F8.3,A,I0,A,I12,A,I12,A,F8.3,A,I12,A,F8.3,A)') &
            &'>>> NU high-resolution extension challenge ', local_stats%old_limit, ' A (k=', &
            &local_stats%old_find, ') -> ', local_stats%new_limit, ' A (k=', local_stats%new_find, &
            &'); frontier ', n_finest, '/', n_total, ' (', local_stats%pct_tested_mask, &
            &'% mask), unary wins ', local_stats%n_unary_wins, ' (', local_stats%pct_unary_wins_tested, '%)'
        old_radius_angstrom = nu_objective_smooth_radius_angstrom(local_stats%old_limit)
        new_radius_angstrom = nu_objective_smooth_radius_angstrom(local_stats%new_limit)
        old_radius_px = nu_objective_smooth_radius_pixels(local_stats%old_limit)
        new_radius_px = nu_objective_smooth_radius_pixels(local_stats%new_limit)
        write(logfhandle,'(A,F8.3,A,I0,A,F8.3,A,I0,A)') &
            &'>>> NU high-resolution extension AWF radii old/new: ', &
            &old_radius_angstrom, ' A (px=', old_radius_px, ') / ', &
            &new_radius_angstrom, ' A (px=', new_radius_px, ')'
        if( l_permissive_accept )then
            local_stats%applied = n_extended > 0
        else
            local_stats%accepted_by_frontier = n_extended >= n_seed_min .and. &
                &local_stats%pct_unary_wins_tested >= accept_pct_eff
            local_stats%applied = n_extended > 0 .and. local_stats%accepted_by_frontier
        endif
        local_stats%promote_next = local_stats%applied
        if( .not. local_stats%applied ) then
            if( n_extended > 0 )then
                write(logfhandle,'(A,F8.3,A,F8.3,A,I0,A,I0,A)') &
                    &'>>> NU high-resolution extension rejected: challenger wins ', &
                    &local_stats%pct_unary_wins_tested, '% below ', accept_pct_eff, &
                    &'% of tested frontier or lacks absolute support ', n_extended, '/', n_seed_min, ' voxels'
            else
                write(logfhandle,'(A,F8.3,A)') &
                    &'>>> NU high-resolution extension stopped: no unary wins for challenger ', new_limit, ' A'
            endif
            call delete_cached_filtered_pair(new_find)
            if( allocated(extend_choice) ) deallocate(extend_choice)
            deallocate(dmat_new, dmat_finest_mask)
            if( present(stats) ) stats = local_stats
            return
        end if
        call compact_nu_highres_dmat_bank_for_capacity()
        if( allocated(dmats_mask) )then
            if( size(dmats_mask,2) >= NU_DMAT_CANDIDATE_CAP )then
                local_stats%applied = .false.
                local_stats%promote_next = .false.
                local_stats%memory_limited = .true.
                write(logfhandle,'(A,I0,A)') &
                    &'>>> NU high-resolution extension stopped: distance-matrix cap full at ', &
                    &NU_DMAT_CANDIDATE_CAP, ' candidates'
                call delete_cached_filtered_pair(new_find)
                if( allocated(extend_choice) ) deallocate(extend_choice)
                deallocate(dmat_new, dmat_finest_mask)
                if( present(stats) ) stats = local_stats
                return
            endif
        endif
        sz_old = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        local_stats%n_extended = n_extended
        if( n_finest > 0 ) local_stats%pct_extended_tested = 100. * real(n_extended) / real(n_finest)
        call apply_nu_highres_extension_selection(frontier_vox, extend_choice, sz_old, sz_old + 1)
        call append_and_thin_nu_highres_candidate(dmat_new, sz_old, new_find)
        call cache_nu_highres_extension_frontier_after_selection(dmat_new, frontier_vox, extend_choice)
        write(logfhandle,'(A,I12,A,F8.2,A)') '>>> Extended ', n_extended, ' voxels to ', new_limit, ' A'
        if( allocated(extend_choice) ) deallocate(extend_choice)
        deallocate(dmat_new, dmat_finest_mask)
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres

    module subroutine extend_nu_filter_highres_shell_next( vol_even, vol_odd, stats, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        real, optional, intent(in) :: accept_pct
        type(nu_highres_extension_stats) :: local_stats
        integer :: new_find, sz_old, frontier_label, n_frontier, ilabel
        real    :: new_limit
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; extend_nu_filter_highres_shell_next')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; extend_nu_filter_highres_shell_next')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; extend_nu_filter_highres_shell_next')
        sz_old   = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        local_stats%n_mask = n_nu_mask
        frontier_label = 0
        n_frontier = 0
        do ilabel = sz_old, 1, -1
            n_frontier = count_nu_base_label_voxels(ilabel)
            if( n_frontier > 0 )then
                frontier_label = ilabel
                exit
            endif
        end do
        if( frontier_label == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        local_stats%old_find  = cutoff_finds(frontier_label)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(frontier_label)
        local_stats%n_tested  = n_frontier
        new_find = cutoff_finds(frontier_label) + 1
        do while( new_find <= box/2 )
            if( .not.any(cutoff_finds(:frontier_label) == new_find) ) exit
            new_find = new_find + 1
        end do
        if( new_find > box/2 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        new_limit = calc_lowpass_lim(new_find, box, smpd)
        if( present(accept_pct) )then
            call extend_nu_filter_highres(vol_even, vol_odd, NU_HIGHRES_EXTENSION_THRESHOLD_PCT, &
                &new_limit, stats=local_stats, accept_pct=accept_pct)
        else
            call extend_nu_filter_highres(vol_even, vol_odd, NU_HIGHRES_EXTENSION_THRESHOLD_PCT, &
                &new_limit, stats=local_stats)
        endif
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres_shell_next

    module subroutine extend_nu_filter_highres_shells( vol_even, vol_odd, nsteps, accept_pct )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(out) :: nsteps
        real, optional, intent(in) :: accept_pct
        type(nu_highres_extension_stats) :: step_stats
        integer :: nsteps_local
        nsteps_local = 0
        do
            if( present(accept_pct) )then
                call extend_nu_filter_highres_shell_next(vol_even, vol_odd, stats=step_stats, accept_pct=accept_pct)
            else
                call extend_nu_filter_highres_shell_next(vol_even, vol_odd, stats=step_stats)
            endif
            if( .not. step_stats%attempted ) exit
            if( .not. step_stats%applied   ) exit
            nsteps_local = nsteps_local + 1
            if( .not. step_stats%promote_next ) exit
        end do
        if( nsteps_local > 0 ) call refine_nu_extension_filtmap_ordered_labels
        if( present(nsteps) ) nsteps = nsteps_local
    end subroutine extend_nu_filter_highres_shells

    module subroutine refine_nu_extension_filtmap_ordered_labels
        integer :: n_base, n_candidates
        if( .not.allocated(filtmap)          ) THROW_HARD('filtmap not allocated; refine_nu_extension_filtmap_ordered_labels')
        if( .not.allocated(dmats_mask)       ) THROW_HARD('dmats_mask not allocated; refine_nu_extension_filtmap_ordered_labels')
        if( .not.allocated(candidate_coords) ) &
            &THROW_HARD('candidate_coords not allocated; refine_nu_extension_filtmap_ordered_labels')
        n_base = size(cutoff_finds)
        n_candidates = size(dmats_mask, 2)
        if( n_candidates /= size(candidate_coords) ) &
            &THROW_HARD('candidate/unary size mismatch; refine_nu_extension_filtmap_ordered_labels')
        if( n_candidates < 2 ) return
        write(logfhandle,'(A)') '>>> NU post-extension ordered-label cleanup'
        call clamp_nu_filtmap_labels(n_base)
        call log_nu_candidate_selection_counts(filtmap, n_base, 'before post-extension ordered-label cleanup')
        call refine_nu_candidate_map_ordered_labels(filtmap, n_candidates)
        call clamp_nu_filtmap_labels(n_base)
        call log_nu_candidate_selection_counts(filtmap, n_base, 'after post-extension ordered-label cleanup')
        call compact_nu_highres_dmat_bank_for_capacity()
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
    end subroutine refine_nu_extension_filtmap_ordered_labels

    subroutine append_and_thin_nu_highres_candidate( dmat_new, old_n_base, new_find )
        real,    intent(in)  :: dmat_new(:,:,:)
        integer, intent(in)  :: old_n_base, new_find
        logical, allocatable :: keep(:)
        integer, allocatable :: old_to_new(:), drop_to_new(:), new_to_old(:), new_cutoff_finds(:)
        real,    allocatable :: new_coords(:), new_bwfilters(:,:), new_dmats(:,:)
        integer :: new_n_base, n_aux, base_keep_n, n_keep, i, j, k, ikeep, old_label, new_label
        integer :: base_find, finest_step, istep, jlabel, src_find, src_label, imask, ibin, iaux
        logical :: l_thinned
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; append_and_thin_nu_highres_candidate')
        if( .not.allocated(dmats_mask)   ) THROW_HARD('dmats_mask not allocated; append_and_thin_nu_highres_candidate')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; append_and_thin_nu_highres_candidate')
        if( any(shape(dmat_new) /= ldim) ) THROW_HARD('dmat_new shape mismatch; append_and_thin_nu_highres_candidate')
        if( old_n_base /= size(cutoff_finds) ) THROW_HARD('base-bank size mismatch; append_and_thin_nu_highres_candidate')
        if( old_n_base < 1 .or. old_n_base > size(dmats_mask,2) ) &
            &THROW_HARD('invalid base-bank size; append_and_thin_nu_highres_candidate')
        new_n_base = old_n_base + 1
        n_aux = max(0, size(dmats_mask,2) - old_n_base)
        allocate(keep(new_n_base), source=.true.)
        base_keep_n = min(size(lowpass_limits), new_n_base)
        if( NU_HIGHRES_EXTENSION_RETAIN_STRIDE > 1 .and. new_n_base > base_keep_n )then
            keep = .false.
            keep(:base_keep_n) = .true.
            base_find = cutoff_finds(min(base_keep_n, old_n_base))
            finest_step = new_find - base_find
            do i = base_keep_n + 1, new_n_base
                if( i <= old_n_base )then
                    src_find = cutoff_finds(i)
                else
                    src_find = new_find
                endif
                istep = src_find - base_find
                keep(i) = keep_nu_highres_extension_step(istep, finest_step)
            end do
        endif
        n_keep = count(keep)
        if( n_keep < 1 ) THROW_HARD('empty retained bank; append_and_thin_nu_highres_candidate')
        allocate(old_to_new(new_n_base), source=0)
        allocate(drop_to_new(new_n_base), source=1)
        allocate(new_to_old(n_keep), source=0)
        allocate(new_cutoff_finds(n_keep))
        ikeep = 0
        do i = 1, new_n_base
            if( keep(i) )then
                ikeep = ikeep + 1
                old_to_new(i) = ikeep
                new_to_old(ikeep) = i
                if( i <= old_n_base )then
                    new_cutoff_finds(ikeep) = cutoff_finds(i)
                else
                    new_cutoff_finds(ikeep) = new_find
                endif
            else if( i <= old_n_base )then
                call delete_cached_filtered_pair(cutoff_finds(i))
            else
                call delete_cached_filtered_pair(new_find)
            endif
        end do
        new_label = old_to_new(new_n_base)
        if( new_label == 0 ) THROW_HARD('accepted high-resolution shell missing after bank thinning')
        do i = 1, new_n_base
            if( old_to_new(i) > 0 )then
                drop_to_new(i) = old_to_new(i)
            else
                do jlabel = i - 1, 1, -1
                    if( old_to_new(jlabel) > 0 )then
                        drop_to_new(i) = old_to_new(jlabel)
                        exit
                    endif
                end do
            endif
        end do
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,old_label) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            old_label = int(filtmap(i,j,k))
            if( old_label >= 1 .and. old_label <= new_n_base )then
                filtmap(i,j,k) = int(drop_to_new(old_label), kind=NU_LABEL_KIND)
            else
                filtmap(i,j,k) = 1_NU_LABEL_KIND
            endif
        end do
        !$omp end parallel do
        call move_alloc(new_cutoff_finds, cutoff_finds)
        if( allocated(bwfilters) )then
            allocate(new_bwfilters(size(bwfilters,1), n_keep), source=0.)
            !$omp parallel do collapse(2) schedule(static) default(shared) private(ikeep,ibin,src_label) proc_bind(close)
            do ikeep = 1, n_keep
                do ibin = 1, size(bwfilters,1)
                    src_label = new_to_old(ikeep)
                    if( src_label <= old_n_base .and. src_label <= size(bwfilters,2) )then
                        new_bwfilters(ibin,ikeep) = bwfilters(ibin,src_label)
                    endif
                end do
            end do
            !$omp end parallel do
            do ikeep = 1, n_keep
                if( new_to_old(ikeep) > old_n_base .or. new_to_old(ikeep) > size(bwfilters,2) )then
                    call butterworth_filter(cutoff_finds(ikeep), new_bwfilters(:,ikeep))
                endif
            end do
            call move_alloc(new_bwfilters, bwfilters)
        endif
        if( allocated(candidate_coords) )then
            allocate(new_coords(n_keep + n_aux), source=0.)
            do ikeep = 1, n_keep
                new_coords(ikeep) = real(ikeep)
            end do
            call move_alloc(new_coords, candidate_coords)
        endif
        allocate(new_dmats(n_nu_mask, n_keep + n_aux))
        !$omp parallel do collapse(2) schedule(static) default(shared) &
        !$omp private(ikeep,imask,src_label,i,j,k) proc_bind(close)
        do ikeep = 1, n_keep
            do imask = 1, n_nu_mask
                src_label = new_to_old(ikeep)
                if( src_label <= old_n_base )then
                    new_dmats(imask,ikeep) = dmats_mask(imask,src_label)
                else
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    new_dmats(imask,ikeep) = dmat_new(i,j,k)
                endif
            end do
        end do
        !$omp end parallel do
        if( n_aux > 0 )then
            !$omp parallel do collapse(2) schedule(static) default(shared) private(iaux,imask) proc_bind(close)
            do iaux = 1, n_aux
                do imask = 1, n_nu_mask
                    new_dmats(imask,n_keep + iaux) = dmats_mask(imask,old_n_base + iaux)
                end do
            end do
            !$omp end parallel do
        endif
        call move_alloc(new_dmats, dmats_mask)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        l_thinned = n_keep < new_n_base
        if( l_thinned )then
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') &
                &'>>> NU high-resolution extension bank thinned: ', new_n_base, &
                &' -> ', n_keep, ' base labels; retaining every ', &
                &NU_HIGHRES_EXTENSION_RETAIN_STRIDE, ' shell step(s)'
        endif
        deallocate(keep, old_to_new, drop_to_new, new_to_old)
    end subroutine append_and_thin_nu_highres_candidate

    subroutine compact_nu_highres_dmat_bank_for_capacity()
        logical, allocatable :: keep(:)
        integer, allocatable :: old_to_new(:), new_cutoff_finds(:)
        real,    allocatable :: new_coords(:), new_bwfilters(:,:), new_dmats(:,:)
        integer :: old_n_base, n_aux, base_keep_n, n_keep, i, j, k, ikeep, old_label, imask, iaux
        if( .not.allocated(dmats_mask)   ) return
        if( .not.allocated(cutoff_finds) ) return
        if( .not.allocated(filtmap)      ) return
        if( .not.allocated(nu_lmask)     ) return
        if( size(dmats_mask,2) < NU_DMAT_CANDIDATE_CAP ) return
        old_n_base = size(cutoff_finds)
        if( old_n_base <= size(lowpass_limits) ) return
        allocate(keep(old_n_base), source=.false.)
        base_keep_n = min(size(lowpass_limits), old_n_base)
        keep(:base_keep_n) = .true.
        do i = base_keep_n + 1, old_n_base
            keep(i) = count_nu_base_label_voxels(i) > 0
        end do
        n_keep = count(keep)
        if( n_keep >= old_n_base )then
            deallocate(keep)
            return
        endif
        allocate(old_to_new(old_n_base), source=0)
        allocate(new_cutoff_finds(n_keep))
        ikeep = 0
        do i = 1, old_n_base
            if( keep(i) )then
                ikeep = ikeep + 1
                old_to_new(i) = ikeep
                new_cutoff_finds(ikeep) = cutoff_finds(i)
            else
                call delete_cached_filtered_pair(cutoff_finds(i))
            endif
        end do
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,old_label)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            old_label = filtmap(i,j,k)
            if( old_label >= 1 .and. old_label <= old_n_base )then
                if( old_to_new(old_label) > 0 )then
                    filtmap(i,j,k) = int(old_to_new(old_label), kind=NU_LABEL_KIND)
                else
                    filtmap(i,j,k) = 1
                endif
            else
                filtmap(i,j,k) = 1
            endif
        end do
        !$omp end parallel do
        call move_alloc(new_cutoff_finds, cutoff_finds)
        if( allocated(bwfilters) )then
            if( size(bwfilters,2) >= old_n_base )then
                allocate(new_bwfilters(size(bwfilters,1), n_keep))
                ikeep = 0
                do i = 1, old_n_base
                    if( .not.keep(i) ) cycle
                    ikeep = ikeep + 1
                    new_bwfilters(:,ikeep) = bwfilters(:,i)
                end do
                call move_alloc(new_bwfilters, bwfilters)
            endif
        endif
        if( allocated(candidate_coords) )then
            n_aux = max(0, size(candidate_coords) - old_n_base)
            allocate(new_coords(n_keep + n_aux), source=0.)
            ikeep = 0
            do i = 1, old_n_base
                if( .not.keep(i) ) cycle
                ikeep = ikeep + 1
                new_coords(ikeep) = real(ikeep)
            end do
            call move_alloc(new_coords, candidate_coords)
        endif
        n_aux = max(0, size(dmats_mask,2) - old_n_base)
        allocate(new_dmats(size(dmats_mask,1), n_keep + n_aux))
        !$omp parallel do collapse(2) schedule(static) default(shared) private(i,imask,ikeep) proc_bind(close)
        do i = 1, old_n_base
            do imask = 1, n_nu_mask
                ikeep = old_to_new(i)
                if( ikeep > 0 ) new_dmats(imask,ikeep) = dmats_mask(imask,i)
            end do
        end do
        !$omp end parallel do
        if( n_aux > 0 )then
            !$omp parallel do collapse(2) schedule(static) default(shared) private(iaux,imask) proc_bind(close)
            do iaux = 1, n_aux
                do imask = 1, n_nu_mask
                    new_dmats(imask,n_keep + iaux) = dmats_mask(imask,old_n_base + iaux)
                end do
            end do
            !$omp end parallel do
        endif
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') &
            &'>>> NU distance-matrix bank compacted for memory: ', size(dmats_mask,2), &
            &' -> ', n_keep + n_aux, ' candidates (cap ', NU_DMAT_CANDIDATE_CAP, ')'
        call move_alloc(new_dmats, dmats_mask)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        deallocate(keep, old_to_new)
    end subroutine compact_nu_highres_dmat_bank_for_capacity

    subroutine pack_nu_dmat_to_mask_vector( dmat_full, dmat_mask )
        real, intent(in)    :: dmat_full(:,:,:)
        real, intent(inout) :: dmat_mask(:)
        integer :: imask, i, j, k
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; pack_nu_dmat_to_mask_vector')
        if( size(dmat_mask) /= n_nu_mask ) THROW_HARD('mask-vector size mismatch; pack_nu_dmat_to_mask_vector')
        if( any(shape(dmat_full) /= ldim) ) THROW_HARD('dmat shape mismatch; pack_nu_dmat_to_mask_vector')
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            dmat_mask(imask) = dmat_full(i,j,k)
        end do
        !$omp end parallel do
    end subroutine pack_nu_dmat_to_mask_vector

    integer function count_nu_base_label_voxels( label ) result(nvox)
        integer, intent(in) :: label
        integer :: imask, i, j, k
        nvox = 0
        if( label < 1 ) return
        if( .not.allocated(nu_mask_vox) ) return
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) reduction(+:nvox) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( int(filtmap(i,j,k)) == label ) nvox = nvox + 1
        end do
        !$omp end parallel do
    end function count_nu_base_label_voxels

    subroutine collect_nu_base_label_voxels( label, frontier_vox )
        integer, intent(in) :: label
        integer, allocatable, intent(inout) :: frontier_vox(:)
        integer :: imask, i, j, k, ifront, nfront, idx
        if( allocated(frontier_vox) ) deallocate(frontier_vox)
        nfront = count_nu_base_label_voxels(label)
        allocate(frontier_vox(nfront))
        if( nfront == 0 ) return
        ifront = 0
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,idx) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( int(filtmap(i,j,k)) /= label ) cycle
            !$omp atomic capture
            ifront = ifront + 1
            idx = ifront
            !$omp end atomic
            frontier_vox(idx) = imask
        end do
        !$omp end parallel do
        if( ifront /= nfront ) THROW_HARD('frontier voxel count changed during collection')
    end subroutine collect_nu_base_label_voxels

    subroutine fill_nu_frontier_dmat_from_bank( frontier_vox, old_label, dmat_finest_mask )
        integer, intent(in) :: frontier_vox(:), old_label
        real,    intent(inout) :: dmat_finest_mask(:)
        integer :: ifront, imask
        real :: x
        if( .not.allocated(dmats_mask) ) THROW_HARD('dmats_mask not allocated; fill_nu_frontier_dmat_from_bank')
        if( old_label < 1 .or. old_label > size(dmats_mask,2) ) &
            &THROW_HARD('old label out of range; fill_nu_frontier_dmat_from_bank')
        if( size(dmat_finest_mask) /= n_nu_mask ) THROW_HARD('mask-vector size mismatch; fill_nu_frontier_dmat_from_bank')
        dmat_finest_mask = huge(x)
        !$omp parallel do schedule(static) default(shared) private(ifront,imask) proc_bind(close)
        do ifront = 1, size(frontier_vox)
            imask = frontier_vox(ifront)
            dmat_finest_mask(imask) = dmats_mask(imask,old_label)
        end do
        !$omp end parallel do
    end subroutine fill_nu_frontier_dmat_from_bank

    module subroutine init_nu_highres_extension_selection( frontier_vox, dmat_old, dmat_new, &
            &extend_choice, n_extended )
        integer, intent(in)    :: frontier_vox(:)
        real,    intent(in)    :: dmat_old(:), dmat_new(:,:,:)
        integer(kind=NU_LABEL_KIND), intent(inout) :: extend_choice(:)
        integer, intent(out)   :: n_extended
        integer :: i, j, k, imask, ifront
        if( size(extend_choice) /= size(frontier_vox) ) &
            &THROW_HARD('extension choice/frontier size mismatch; init_nu_highres_extension_selection')
        extend_choice = 0_NU_LABEL_KIND
        n_extended = 0
        !$omp parallel do schedule(static) default(shared) private(ifront,i,j,k,imask) reduction(+:n_extended) proc_bind(close)
        do ifront = 1, size(frontier_vox)
            imask = frontier_vox(ifront)
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( dmat_new(i,j,k) < dmat_old(imask) )then
                extend_choice(ifront) = 1_NU_LABEL_KIND
                n_extended = n_extended + 1
            endif
        end do
        !$omp end parallel do
    end subroutine init_nu_highres_extension_selection

    module subroutine apply_nu_highres_extension_selection( frontier_vox, extend_choice, old_label, new_label )
        integer, intent(in) :: frontier_vox(:)
        integer(kind=NU_LABEL_KIND), intent(in) :: extend_choice(:)
        integer, intent(in) :: old_label, new_label
        integer :: i, j, k, imask, ifront, choice
        if( size(extend_choice) /= size(frontier_vox) ) &
            &THROW_HARD('extension choice/frontier size mismatch; apply_nu_highres_extension_selection')
        if( old_label < 1 ) THROW_HARD('invalid old label; apply_nu_highres_extension_selection')
        !$omp parallel do schedule(static) default(shared) private(ifront,imask,i,j,k,choice) proc_bind(close)
        do ifront = 1, size(frontier_vox)
            choice = int(extend_choice(ifront))
            if( choice /= 1 ) cycle
            imask = frontier_vox(ifront)
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            filtmap(i,j,k) = int(new_label, kind=NU_LABEL_KIND)
        end do
        !$omp end parallel do
    end subroutine apply_nu_highres_extension_selection

    subroutine cache_nu_highres_extension_frontier_after_selection( dmat_new, frontier_vox, extend_choice )
        real,    intent(in) :: dmat_new(:,:,:)
        integer, intent(in) :: frontier_vox(:)
        integer(kind=NU_LABEL_KIND), intent(in) :: extend_choice(:)
        integer :: i, j, k, imask, ifront
        if( size(extend_choice) /= size(frontier_vox) ) &
            &THROW_HARD('extension choice/frontier size mismatch; cache_nu_highres_extension_frontier_after_selection')
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(n_nu_mask), source=huge(0.))
        !$omp parallel do schedule(static) default(shared) private(ifront,imask,i,j,k) proc_bind(close)
        do ifront = 1, size(frontier_vox)
            if( int(extend_choice(ifront)) /= 1 ) cycle
            imask = frontier_vox(ifront)
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            dmat_finest_cached(imask) = dmat_new(i,j,k)
        end do
        !$omp end parallel do
    end subroutine cache_nu_highres_extension_frontier_after_selection

    subroutine prune_nu_highres_extension_bank( active_label )
        integer, intent(in) :: active_label
        integer, allocatable :: new_cutoff_finds(:)
        real,    allocatable :: new_coords(:), new_bwfilters(:,:), new_dmats(:,:)
        integer :: old_n_base, n_aux, i
        if( .not. allocated(cutoff_finds) ) return
        old_n_base = size(cutoff_finds)
        if( active_label < 1 .or. active_label >= old_n_base ) return
        do i = active_label + 1, old_n_base
            call delete_cached_filtered_pair(cutoff_finds(i))
        end do
        allocate(new_cutoff_finds(active_label))
        new_cutoff_finds = cutoff_finds(:active_label)
        call move_alloc(new_cutoff_finds, cutoff_finds)
        if( allocated(bwfilters) )then
            if( size(bwfilters,2) >= old_n_base )then
                allocate(new_bwfilters(size(bwfilters,1), active_label))
                new_bwfilters = bwfilters(:,1:active_label)
                call move_alloc(new_bwfilters, bwfilters)
            endif
        endif
        if( allocated(candidate_coords) )then
            n_aux = max(0, size(candidate_coords) - old_n_base)
            allocate(new_coords(active_label + n_aux), source=0.)
            do i = 1, active_label
                new_coords(i) = real(i)
            end do
            call move_alloc(new_coords, candidate_coords)
        endif
        if( allocated(dmats_mask) )then
            n_aux = max(0, size(dmats_mask, 2) - old_n_base)
            allocate(new_dmats(size(dmats_mask,1), active_label + n_aux))
            new_dmats(:,:active_label) = dmats_mask(:,:active_label)
            if( n_aux > 0 ) new_dmats(:,active_label + 1:) = dmats_mask(:,old_n_base + 1:)
            call move_alloc(new_dmats, dmats_mask)
        endif
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
    end subroutine prune_nu_highres_extension_bank

end submodule simple_nu_filter_extend
