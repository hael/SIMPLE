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
        real, allocatable :: dmat_new(:,:,:), dmat_tmp(:,:,:), dmat_finest(:,:,:)
        integer, allocatable :: cutoff_finds_new(:)
        integer           :: new_find, n_finest, n_total, n_extended, sz_old, old_label
        integer           :: old_radius_px, new_radius_px, n_seed_min
        real              :: accept_pct_eff, pct_finest, x, noise_sigma
        real              :: old_radius_angstrom, new_radius_angstrom
        logical, allocatable :: extend_mask(:,:,:), extend_to_new(:,:,:)
        integer, allocatable :: extend_choice(:,:,:)
        type(nu_highres_extension_stats) :: local_stats
        integer           :: i, j, k
        logical           :: l_use_aux_extension, l_permissive_accept
        local_stats%new_limit = new_limit
        accept_pct_eff = NU_HIGHRES_EXTENSION_ACCEPT_PCT
        if( present(accept_pct) ) accept_pct_eff = max(0., accept_pct)
        l_permissive_accept = present(accept_pct) .and. accept_pct_eff <= TINY
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats first')
        sz_old    = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        n_total   = count(nu_lmask)
        local_stats%n_mask = n_total
        if( n_total == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        old_label = 0
        n_finest  = 0
        do i = sz_old, 1, -1
            n_finest = count(nu_lmask .and. srcmap == 1 .and. filtmap == i)
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
        if( .not.file_exists(even_cache_fname) .or. .not.file_exists(odd_cache_fname) ) then
            ! generate just this one new filtered pair — cheap, one FFT+filter+IFFT pass
            call generate_single_filtered_pair(vol_even, vol_odd, new_find, even_cache_fname, odd_cache_fname)
        end if
        ! --- build the extend mask: voxels currently at the finest limit ---
        allocate(extend_mask(ldim(1),ldim(2),ldim(3)), source=.false.)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k)
        do k = 1, ldim(3)
        do j = 1, ldim(2)
            do i = 1, ldim(1)
                extend_mask(i,j,k) = (nu_lmask(i,j,k) .and. srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == sz_old)
            end do
        end do
        end do
        !$omp end parallel do
        ! --- evaluate the new objective only within the mask ---
        allocate(dmat_new(ldim(1),ldim(2),ldim(3)), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)), source=0.)
        call vol_even_filt_new%new(ldim, smpd)
        call vol_odd_filt_new%new(ldim, smpd)
        call vol_even_filt_new%read(even_cache_fname)
        call vol_odd_filt_new%read(odd_cache_fname)
        noise_sigma = vol_even%nu_objective_noise_scale(vol_odd, nu_lmask)
        call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_new, &
            &nu_lmask, noise_sigma)
        call smooth_nu_objective(dmat_new, dmat_tmp, new_limit)
        ! dmat_tmp is only a work buffer; smooth_nu_objective updates dmat_new in place.
        allocate(dmat_finest(ldim(1),ldim(2),ldim(3)), source=huge(x))
        if( allocated(dmat_finest_cached) ) then
            if( all(shape(dmat_finest_cached) == ldim) ) then
                dmat_finest = dmat_finest_cached
            else
                call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
                call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
                call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, &
                    &dmat_finest, nu_lmask, noise_sigma)
                call smooth_nu_objective(dmat_finest, dmat_tmp, local_stats%old_limit)
            end if
        else
            call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
            call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
            call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, &
                &dmat_finest, nu_lmask, noise_sigma)
            call smooth_nu_objective(dmat_finest, dmat_tmp, local_stats%old_limit)
        end if
        ! --- update filtmap in place for the masked voxels ---
        l_use_aux_extension = l_aux_source_unordered_potts .and. allocated(dmats_aux_mask)
        n_extended = 0
        if( l_use_aux_extension )then
            allocate(extend_choice(ldim(1),ldim(2),ldim(3)), source=0)
            call init_nu_highres_extension_selection_aux(extend_mask, dmat_finest, dmat_new, &
                &extend_choice, n_extended)
        else
            allocate(extend_to_new(ldim(1),ldim(2),ldim(3)), source=.false.)
            call init_nu_highres_extension_selection(extend_mask, dmat_finest, dmat_new, extend_to_new, n_extended)
        endif
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
            call vol_even_filt_new%kill
            call vol_odd_filt_new%kill
            if( allocated(extend_to_new) ) deallocate(extend_to_new)
            if( allocated(extend_choice) ) deallocate(extend_choice)
            deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
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
                call vol_even_filt_new%kill
                call vol_odd_filt_new%kill
                if( allocated(extend_to_new) ) deallocate(extend_to_new)
                if( allocated(extend_choice) ) deallocate(extend_choice)
                deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
                if( present(stats) ) stats = local_stats
                return
            endif
        endif
        sz_old = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        local_stats%n_extended = n_extended
        if( n_finest > 0 ) local_stats%pct_extended_tested = 100. * real(n_extended) / real(n_finest)
        if( l_use_aux_extension )then
            call apply_nu_highres_extension_selection_aux(extend_choice, sz_old, sz_old + 1)
        else
            call apply_nu_highres_extension_selection(extend_to_new, sz_old + 1)
        endif
        ! --- grow cutoff_finds to include the new level ---
        allocate(cutoff_finds_new(sz_old + 1))
        cutoff_finds_new(:sz_old)  = cutoff_finds
        cutoff_finds_new(sz_old+1) = new_find
        call move_alloc(cutoff_finds_new, cutoff_finds)
        call append_nu_highres_dmat_candidate(dmat_new, sz_old)
        if( allocated(candidate_coords) .and. size(candidate_coords) >= sz_old )then
            call append_nu_highres_candidate_coord(sz_old, candidate_coords(sz_old) + &
                &real(new_find - cutoff_finds(sz_old)))
        else
            call append_nu_highres_candidate_coord(sz_old, real(sz_old + 1))
        endif
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( l_use_aux_extension )then
            call cache_nu_highres_extension_frontier_after_aux_selection(dmat_new, sz_old + 1)
        else
            allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=dmat_new)
        endif
        write(logfhandle,'(A,I12,A,F8.2,A)') '>>> Extended ', n_extended, ' voxels to ', new_limit, ' A'
        call vol_even_filt_new%kill
        call vol_odd_filt_new%kill
        if( allocated(extend_to_new) ) deallocate(extend_to_new)
        if( allocated(extend_choice) ) deallocate(extend_choice)
        deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
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
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; extend_nu_filter_highres_shell_next')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; extend_nu_filter_highres_shell_next')
        sz_old   = size(cutoff_finds)
        local_stats%old_find  = cutoff_finds(sz_old)
        local_stats%old_limit = cutoff_find_to_lowpass_limit(sz_old)
        local_stats%n_mask = count(nu_lmask)
        frontier_label = 0
        n_frontier = 0
        do ilabel = sz_old, 1, -1
            n_frontier = count(nu_lmask .and. srcmap == 1 .and. filtmap == ilabel)
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
        if( nsteps_local > 0 ) call refine_nu_extension_filtmap_ordered_labels()
        if( present(nsteps) ) nsteps = nsteps_local
    end subroutine extend_nu_filter_highres_shells

    module subroutine refine_nu_extension_filtmap_ordered_labels
        integer, allocatable :: candmap(:,:,:)
        integer :: i, j, k, n_base, n_candidates, aux_icand
        if( .not.allocated(filtmap)          ) THROW_HARD('filtmap not allocated; refine_nu_extension_filtmap_ordered_labels')
        if( .not.allocated(srcmap)           ) THROW_HARD('srcmap not allocated; refine_nu_extension_filtmap_ordered_labels')
        if( .not.allocated(dmats_mask)       ) THROW_HARD('dmats_mask not allocated; refine_nu_extension_filtmap_ordered_labels')
        if( .not.allocated(candidate_coords) ) &
            &THROW_HARD('candidate_coords not allocated; refine_nu_extension_filtmap_ordered_labels')
        n_base = size(cutoff_finds)
        n_candidates = size(dmats_mask, 2)
        if( n_candidates /= size(candidate_coords) ) &
            &THROW_HARD('candidate/unary size mismatch; refine_nu_extension_filtmap_ordered_labels')
        if( n_candidates < 2 ) return
        write(logfhandle,'(A)') '>>> NU post-extension ordered-label cleanup'
        allocate(candmap(ldim(1),ldim(2),ldim(3)), source=1)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,aux_icand) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    if( srcmap(i,j,k) == 1 )then
                        candmap(i,j,k) = max(1, min(n_base, filtmap(i,j,k)))
                    else
                        aux_icand = n_base + srcmap(i,j,k) - 1
                        candmap(i,j,k) = max(1, min(n_candidates, aux_icand))
                    endif
                end do
            end do
        end do
        !$omp end parallel do
        call log_nu_candidate_selection_counts(candmap, n_base, 'before post-extension ordered-label cleanup')
        call refine_nu_candidate_map_ordered_labels(candmap, n_candidates)
        call log_nu_candidate_selection_counts(candmap, n_base, 'after post-extension ordered-label cleanup')
        call candidate_map_to_filt_and_src(candmap, n_base)
        call compact_nu_highres_dmat_bank_for_capacity()
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        deallocate(candmap)
    end subroutine refine_nu_extension_filtmap_ordered_labels

    subroutine compact_nu_highres_dmat_bank_for_capacity()
        logical, allocatable :: keep(:)
        integer, allocatable :: old_to_new(:), new_cutoff_finds(:)
        real,    allocatable :: new_coords(:), new_bwfilters(:,:), new_dmats(:,:)
        integer :: old_n_base, n_aux, base_keep_n, n_keep, i, j, k, ikeep, old_label
        if( .not.allocated(dmats_mask)   ) return
        if( .not.allocated(cutoff_finds) ) return
        if( .not.allocated(filtmap)      ) return
        if( .not.allocated(srcmap)       ) return
        if( .not.allocated(nu_lmask)     ) return
        if( size(dmats_mask,2) < NU_DMAT_CANDIDATE_CAP ) return
        old_n_base = size(cutoff_finds)
        if( old_n_base <= size(lowpass_limits) ) return
        allocate(keep(old_n_base), source=.false.)
        base_keep_n = min(size(lowpass_limits), old_n_base)
        keep(:base_keep_n) = .true.
        do i = base_keep_n + 1, old_n_base
            keep(i) = any(nu_lmask .and. srcmap == 1 .and. filtmap == i)
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
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,old_label)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( srcmap(i,j,k) /= 1 ) cycle
                    old_label = filtmap(i,j,k)
                    if( old_label >= 1 .and. old_label <= old_n_base )then
                        if( old_to_new(old_label) > 0 )then
                            filtmap(i,j,k) = old_to_new(old_label)
                        else
                            filtmap(i,j,k) = 1
                        endif
                    else
                        filtmap(i,j,k) = 1
                    endif
                end do
            end do
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
                new_coords(ikeep) = candidate_coords(i)
            end do
            if( n_aux > 0 ) new_coords(n_keep + 1:) = candidate_coords(old_n_base + 1:)
            call move_alloc(new_coords, candidate_coords)
        endif
        n_aux = max(0, size(dmats_mask,2) - old_n_base)
        allocate(new_dmats(size(dmats_mask,1), n_keep + n_aux), source=huge(0.))
        ikeep = 0
        do i = 1, old_n_base
            if( .not.keep(i) ) cycle
            ikeep = ikeep + 1
            new_dmats(:,ikeep) = dmats_mask(:,i)
        end do
        if( n_aux > 0 ) new_dmats(:,n_keep + 1:) = dmats_mask(:,old_n_base + 1:)
        write(logfhandle,'(A,I0,A,I0,A,I0,A)') &
            &'>>> NU distance-matrix bank compacted for memory: ', size(dmats_mask,2), &
            &' -> ', n_keep + n_aux, ' candidates (cap ', NU_DMAT_CANDIDATE_CAP, ')'
        call move_alloc(new_dmats, dmats_mask)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        deallocate(keep, old_to_new)
    end subroutine compact_nu_highres_dmat_bank_for_capacity

    module subroutine init_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(out)   :: n_extended
        integer :: i, j, k
        extend_to_new = .false.
        n_extended = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:n_extended)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    if( dmat_new(i,j,k) < dmat_old(i,j,k) )then
                        extend_to_new(i,j,k) = .true.
                        n_extended = n_extended + 1
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine init_nu_highres_extension_selection

    module subroutine init_nu_highres_extension_selection_aux( extend_mask, dmat_old, dmat_new, &
            &extend_choice, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer, intent(inout) :: extend_choice(:,:,:)
        integer, intent(out)   :: n_extended
        integer :: i, j, k, iaux, imask, best_choice
        real    :: best_dmat, cur_dmat
        extend_choice = 0
        n_extended = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,iaux,imask,best_choice,best_dmat,cur_dmat) reduction(+:n_extended)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    imask = nu_mask_index(i,j,k)
                    best_choice = 0
                    best_dmat   = dmat_old(i,j,k)
                    if( dmat_new(i,j,k) < best_dmat )then
                        best_choice = 1
                        best_dmat   = dmat_new(i,j,k)
                    endif
                    do iaux = 1, size(dmats_aux_mask,2)
                        cur_dmat = dmats_aux_mask(imask,iaux)
                        if( cur_dmat < best_dmat )then
                            best_choice = iaux + 1
                            best_dmat   = cur_dmat
                        endif
                    end do
                    extend_choice(i,j,k) = best_choice
                    if( best_choice == 1 ) n_extended = n_extended + 1
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine init_nu_highres_extension_selection_aux

    module subroutine apply_nu_highres_extension_selection( extend_to_new, new_label )
        logical, intent(in) :: extend_to_new(:,:,:)
        integer, intent(in) :: new_label
        integer :: i, j, k
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_to_new(i,j,k) ) cycle
                    srcmap(i,j,k)  = 1
                    filtmap(i,j,k) = new_label
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine apply_nu_highres_extension_selection

    module subroutine apply_nu_highres_extension_selection_aux( extend_choice, old_label, new_label )
        integer, intent(in) :: extend_choice(:,:,:)
        integer, intent(in) :: old_label, new_label
        integer :: i, j, k, choice, aux_icand
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,choice,aux_icand)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    choice = extend_choice(i,j,k)
                    if( choice == 0 ) cycle
                    if( choice == 1 )then
                        srcmap(i,j,k)  = 1
                        filtmap(i,j,k) = new_label
                    else
                        srcmap(i,j,k) = choice
                        aux_icand = old_label + choice - 1
                        filtmap(i,j,k) = nu_effective_base_label_for_candidate(aux_icand, old_label)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine apply_nu_highres_extension_selection_aux

    subroutine cache_nu_highres_extension_frontier_after_aux_selection( dmat_new, new_label )
        real,    intent(in) :: dmat_new(:,:,:)
        integer, intent(in) :: new_label
        integer :: i, j, k, imask, iaux
        allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=huge(0.))
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,imask,iaux)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    if( filtmap(i,j,k) /= new_label ) cycle
                    if( srcmap(i,j,k) == 1 )then
                        dmat_finest_cached(i,j,k) = dmat_new(i,j,k)
                    else
                        iaux = srcmap(i,j,k) - 1
                        if( iaux < 1 .or. iaux > size(dmats_aux_mask,2) ) cycle
                        imask = nu_mask_index(i,j,k)
                        dmat_finest_cached(i,j,k) = dmats_aux_mask(imask,iaux)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cache_nu_highres_extension_frontier_after_aux_selection

    module subroutine append_nu_highres_candidate_coord( old_n_base, new_coord )
        integer, intent(in) :: old_n_base
        real,    intent(in) :: new_coord
        real, allocatable :: new_coords(:)
        integer :: n_aux
        if( .not.allocated(candidate_coords) ) return
        if( size(candidate_coords) < old_n_base ) &
            &THROW_HARD('candidate_coords shorter than base bank in append_nu_highres_candidate_coord')
        n_aux = size(candidate_coords) - old_n_base
        allocate(new_coords(old_n_base + 1 + n_aux), source=0.)
        new_coords(:old_n_base) = candidate_coords(:old_n_base)
        new_coords(old_n_base + 1) = new_coord
        if( n_aux > 0 ) new_coords(old_n_base + 2:) = candidate_coords(old_n_base + 1:)
        call move_alloc(new_coords, candidate_coords)
    end subroutine append_nu_highres_candidate_coord

    subroutine append_nu_highres_dmat_candidate( dmat_new, old_n_base )
        real,    intent(in) :: dmat_new(:,:,:)
        integer, intent(in) :: old_n_base
        real, allocatable :: new_dmats(:,:)
        integer :: i, j, k, imask, n_aux, n_candidates_old
        if( .not.allocated(dmats_mask)  ) THROW_HARD('dmats_mask not allocated; append_nu_highres_dmat_candidate')
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; append_nu_highres_dmat_candidate')
        if( any(shape(dmat_new) /= ldim) ) THROW_HARD('dmat_new shape mismatch; append_nu_highres_dmat_candidate')
        n_candidates_old = size(dmats_mask, 2)
        if( old_n_base < 1 .or. old_n_base > n_candidates_old ) &
            &THROW_HARD('base-bank size mismatch; append_nu_highres_dmat_candidate')
        n_aux = n_candidates_old - old_n_base
        allocate(new_dmats(n_nu_mask, old_n_base + 1 + n_aux), source=huge(0.))
        new_dmats(:,:old_n_base) = dmats_mask(:,:old_n_base)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            new_dmats(imask,old_n_base + 1) = dmat_new(i,j,k)
        end do
        !$omp end parallel do
        if( n_aux > 0 ) new_dmats(:,old_n_base + 2:) = dmats_mask(:,old_n_base + 1:)
        call move_alloc(new_dmats, dmats_mask)
    end subroutine append_nu_highres_dmat_candidate

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
            allocate(new_coords(active_label + n_aux))
            new_coords(:active_label) = candidate_coords(:active_label)
            if( n_aux > 0 ) new_coords(active_label + 1:) = candidate_coords(old_n_base + 1:)
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
