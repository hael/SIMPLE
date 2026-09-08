!@descr: simple nu filter bank implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_bank
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine setup_nu_dmats( vol_even, vol_odd, mskdiam, aux_resolutions, aux_even, aux_odd, &
            &n_highres_steps, evidence_source, fsc_res )
        class(image),          intent(in) :: vol_even, vol_odd
        real,                  intent(in) :: mskdiam
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
        character(len=*), optional, intent(in) :: evidence_source
        real,        optional, intent(in) :: fsc_res
        type(image) :: vol_even_filt, vol_odd_filt, vol_support
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:), dmat_cand(:,:,:)
        real, allocatable :: noise_profile(:)
        real :: noise_rmax, finest_lp
        integer :: i, n_candidates, aux_replacement_idx
        real    :: x
        call init_nu_filter(vol_even, vol_odd, n_highres_steps, fsc_res)
        if( nu_l_report .and. nu_bank_cap_find > 0 )then
            write(logfhandle,'(A,F8.3,A,F8.3,A,I0,A,I0,A)') '>>> NU BANK CAP: FSC=0.143 ', fsc_res, &
                &' A; candidates finer than ', calc_lowpass_lim(nu_bank_cap_find, box, smpd), &
                &' A dropped; ', size(cutoff_finds), ' candidates retained of ', size(lowpass_limits), ' static labels'
        endif
        if( mskdiam <= TINY ) THROW_HARD('mskdiam must be positive in setup_nu_dmats')
        nu_support_mskdiam = mskdiam
        if( allocated(nu_lmask) ) deallocate(nu_lmask)
        call vol_support%disc(ldim, smpd, 0.5 * mskdiam / smpd, nu_lmask)
        call vol_support%kill
        if( .not. any(nu_lmask) ) THROW_HARD('spherical support has no true voxels in setup_nu_dmats')
        call setup_nu_mask_voxels
        nu_evidence_requested = present(evidence_source)
        nu_evidence_source = ''
        nu_evidence_source_fingerprint = 0.d0
        if( nu_evidence_requested )then
            if( trim(evidence_source) /= NU_EVIDENCE_SOURCE_BASE .and. &
                &trim(evidence_source) /= NU_EVIDENCE_SOURCE_PREV )then
                THROW_HARD('NU replay evidence source must be base_unfil or previous_shipped')
            endif
            nu_evidence_source = trim(evidence_source)
        endif
        aux_replacement_idx = 0
        if( present(aux_even) ) then
            if( .not. present(aux_odd) ) THROW_HARD('Auxiliary odd bank missing; setup_nu_dmats')
            if( size(aux_resolutions) /= size(aux_even) ) THROW_HARD('Auxiliary resolutions size mismatch; setup_nu_dmats')
            call validate_aux_volumes(aux_even, aux_odd)
            finest_lp = cutoff_find_to_lowpass_limit(size(cutoff_finds))
            do i = 1, size(aux_resolutions)
                if( aux_resolutions(i) <= TINY ) THROW_HARD('Auxiliary resolution must be positive; setup_nu_dmats')
                if( aux_resolutions(i) < finest_lp - TINY )then
                    if( aux_replacement_idx == 0 )then
                        aux_replacement_idx = i
                    else if( aux_resolutions(i) < aux_resolutions(aux_replacement_idx) )then
                        aux_replacement_idx = i
                    endif
                endif
            end do
            if( aux_replacement_idx > 0 )then
                call stash_aux_volumes(aux_even(aux_replacement_idx:aux_replacement_idx), &
                    &aux_odd(aux_replacement_idx:aux_replacement_idx))
                nu_aux_replacement_label = size(cutoff_finds)
                nu_aux_replacement_resolution = aux_resolutions(aux_replacement_idx)
                if( NU_DEV_OUTPUT .and. nu_l_report )then
                    write(logfhandle,'(A,I0,A,F8.3,A,F8.3,A)') &
                        &'>>> NU auxiliary replacement: auxiliary pair ', aux_replacement_idx, &
                        &' replaces finest discrete label at ', finest_lp, ' A with effective ', &
                        &nu_aux_replacement_resolution, ' A'
                endif
            else
                call cleanup_aux_bank
                if( NU_DEV_OUTPUT .and. nu_l_report .and. size(aux_resolutions) > 0 )then
                    write(logfhandle,'(A,F8.3,A,F8.3,A)') &
                        &'>>> NU auxiliary ignored: finest supplied effective resolution ', minval(aux_resolutions), &
                        &' A does not extend beyond finest discrete label ', finest_lp, ' A'
                endif
            endif
        else
            if( present(aux_odd) ) THROW_HARD('Auxiliary odd bank supplied without even bank; setup_nu_dmats')
            if( size(aux_resolutions) /= 0 ) THROW_HARD('Auxiliary resolutions supplied without auxiliary volumes; setup_nu_dmats')
            call cleanup_aux_bank
        end if
        if( nu_evidence_requested .and. nu_aux_replacement_label > 0 ) &
            &THROW_HARD('NU replay evidence cannot include an auxiliary replacement pair')
        if( nu_evidence_requested ) &
            &call calculate_nu_source_fingerprint(vol_even, vol_odd, nu_evidence_source_fingerprint)
        ! Filter caches are local scratch products. Rebuild the current bank
        ! after setup so a prior interrupted run cannot satisfy existence-only
        ! cache checks with stale volumes.
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_EVEN))
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_ODD))
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call cache_filtered_vols(vol_even, vol_odd)
        call vol_even%nu_objective_noise_profile(vol_odd, nu_lmask, noise_profile, noise_rmax)
        ! Cache for reuse during high-resolution shell extension; the raw
        ! E/O noise profile is candidate-independent so it does not need to be
        ! recomputed per shell challenge.
        if( allocated(nu_noise_profile_cached) ) deallocate(nu_noise_profile_cached)
        nu_noise_profile_cached = noise_profile
        nu_noise_rmax_cached    = noise_rmax
        call setup_nu_observed_mask(vol_even, vol_odd)
        if( NU_DEV_OUTPUT .and. nu_l_report ) &
            &write(logfhandle,'(A,I0,A,ES11.4,A,ES11.4,A,F6.3)') '>>> NU WHITENING PROFILE: ', &
            &size(noise_profile), ' shells, sigma(r) min ', minval(noise_profile), ' max ', &
            &maxval(noise_profile), ' edge/centre ', noise_profile(size(noise_profile))/max(noise_profile(1),TINY)
        if( allocated(dmats_mask) ) deallocate(dmats_mask)
        n_candidates = size(cutoff_finds)
        if( n_candidates > NU_DMAT_CANDIDATE_CAP )then
            THROW_HARD('NU distance-matrix candidate cap exceeded in setup_nu_dmats')
        endif
        call setup_nu_candidate_coords(n_candidates)
        if( NU_DEV_OUTPUT .and. nu_l_report ) call log_nu_objective_smoothing_bank()
        allocate(dmats_mask(n_nu_mask,n_candidates), source=huge(x))
        if( allocated(raw_dmats_mask) ) deallocate(raw_dmats_mask)
        allocate(raw_dmats_mask(n_nu_mask,n_candidates), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),  source=0.)
        allocate(dmat_cand(ldim(1),ldim(2),ldim(3)), source=huge(x))
        if( allocated(nu_ev_base) ) deallocate(nu_ev_base)
        if( allocated(nu_ev_best) ) deallocate(nu_ev_best)
        allocate(nu_ev_base(n_nu_mask), source=0.)
        allocate(nu_ev_best(n_nu_mask), source=huge(x))
        do i = 1, size(cutoff_finds)
            dmat_cand = huge(x)
            if( nu_label_is_aux_replacement(i) )then
                call vol_even%nu_objective(aux_even_bank(1), vol_odd, aux_odd_bank(1), dmat_cand, &
                    &nu_lmask, noise_profile, noise_rmax)
            else
                even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
                odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
                if( .not.file_exists(even_cache_fname) ) THROW_HARD('Missing filtered volume cache: '//even_cache_fname%to_char())
                if( .not.file_exists(odd_cache_fname)  ) THROW_HARD('Missing filtered volume cache: '//odd_cache_fname%to_char())
                call vol_even_filt%read(even_cache_fname)
                call vol_odd_filt%read(odd_cache_fname)
                call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmat_cand, &
                    &nu_lmask, noise_profile, noise_rmax)
            endif
            ! Snapshot the raw cost before candidate-scale smoothing; the envelope
            ! needs terms that were blurred identically, not per-candidate.
            call accumulate_nu_evidence_raw(dmat_cand, i)
            call pack_nu_raw_candidate(dmat_cand, i)
            call smooth_nu_objective(dmat_cand, dmat_tmp, nu_label_lowpass_limit(i))
            call pack_nu_dmat_candidate(dmat_cand, i)
        end do
        call vol_even_filt%kill
        call vol_odd_filt%kill
        deallocate(dmat_tmp, dmat_cand)
        call release_nu_smooth_norm()
    end subroutine setup_nu_dmats


    module subroutine setup_nu_candidate_coords( n_candidates )
        integer, intent(in) :: n_candidates
        integer :: i, n_base
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; setup_nu_candidate_coords')
        n_base = size(cutoff_finds)
        if( n_candidates /= n_base ) THROW_HARD('candidate count must match base bank in setup_nu_candidate_coords')
        if( allocated(candidate_coords) ) deallocate(candidate_coords)
        allocate(candidate_coords(n_candidates), source=0.)
        do i = 1, n_base
            candidate_coords(i) = real(i)
        end do
    end subroutine setup_nu_candidate_coords

    module real function get_nu_filter_bank_finest_lp()
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; get_nu_filter_bank_finest_lp')
        if( size(cutoff_finds) < 1 ) THROW_HARD('empty filter bank; get_nu_filter_bank_finest_lp')
        get_nu_filter_bank_finest_lp = nu_label_lowpass_limit(size(cutoff_finds))
    end function get_nu_filter_bank_finest_lp

    module integer function get_nu_filtmap_highres_shell_depth()
        integer :: base_n, finest_label, i, j, k, imask
        get_nu_filtmap_highres_shell_depth = 0
        if( .not.allocated(cutoff_finds) ) return
        if( .not.allocated(filtmap)      ) return
        if( .not.allocated(nu_lmask)     ) return
        base_n = min(size(lowpass_limits), size(cutoff_finds))
        if( base_n < 1 ) return
        finest_label = 0
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            finest_label = max(finest_label, int(filtmap(i,j,k)))
        end do
        if( finest_label == 0 ) return
        if( finest_label <= base_n ) return
        finest_label = min(finest_label, size(cutoff_finds))
        get_nu_filtmap_highres_shell_depth = max(0, cutoff_finds(finest_label) - cutoff_finds(base_n))
    end function get_nu_filtmap_highres_shell_depth

    module subroutine optimize_nu_cutoff_finds()
        integer :: nx, ny, nz, i, j, k, icand, n_base, n_candidates, imask, n_clamped, ilevel
        integer, allocatable :: sel(:)
        real,    allocatable :: lvl(:,:), full(:,:,:), tmp(:,:,:)
        real    :: lp_level
        if( .not.allocated(dmats_mask) ) THROW_HARD('dmats_mask not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        n_base       = size(cutoff_finds)
        ! dmats_mask has one column per retained label. If an auxiliary pair is
        ! eligible, it backs the finest label rather than appending a new one.
        n_candidates = size(dmats_mask, 2)
        if( .not.allocated(raw_dmats_mask) ) THROW_HARD('raw_dmats_mask not allocated; run setup_nu_dmats before optimize_nu_cutoff_finds')
        if( allocated(filtmap) ) deallocate(filtmap)
        allocate(filtmap(nx,ny,nz), source=1_NU_LABEL_KIND)
        ! Coarse-to-fine selection with like-for-like smoothing (2026-09-08).
        ! dmats_mask holds each candidate smoothed at its own radius (1.5 x LP),
        ! so an argmin over it compares differently smoothed fields: two
        ! candidates with near-identical raw unaries do not tie, the smaller
        ! radius wins at local minima of the unary field, the larger at maxima,
        ! an intermediate one almost never. An honest gridding pair never
        ! exposes this (adjacent fine candidates differ by the admitted noise
        ! band); a regularized pair does, and the populated fine label then
        ! follows the radius table (PfCRT 2026-09-07: box 140, radii 4/3/3 px
        ! for 5.97/5.0/4.44 A, 7.2% at 5.0 and 0.08% at 4.44; box 150, 4/4/3,
        ! 0.1% at 5.0 and 3% at 4.14). Here each finer candidate replaces the
        ! incumbent only if it wins at ITS scale with both smoothed alike, so
        ! identical unaries tie exactly and the coarser label keeps (strict <).
        ! The cost is the smoothing passes: n(n+1)/2 - 1 instead of n.
        allocate(sel(n_nu_mask), source=1)
        allocate(lvl(n_nu_mask,n_candidates), source=0.)
        allocate(full(nx,ny,nz), source=0.)
        allocate(tmp(nx,ny,nz),  source=0.)
        do ilevel = 2, n_candidates
            lp_level = nu_label_lowpass_limit(ilevel)
            do icand = 1, ilevel
                call unpack_nu_raw_candidate(icand, full)
                call smooth_nu_objective(full, tmp, lp_level)
                call pack_nu_full_to_mask(full, lvl(:,icand))
            end do
            !$omp parallel do schedule(static) default(shared) private(imask) proc_bind(close)
            do imask = 1, n_nu_mask
                if( lvl(imask,ilevel) < lvl(imask,sel(imask)) ) sel(imask) = ilevel
            end do
            !$omp end parallel do
        end do
        deallocate(lvl, full, tmp)
        call release_nu_smooth_norm()
        !$omp parallel do schedule(static) default(shared) private(i,j,k,imask) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            filtmap(i,j,k) = int(sel(imask), kind=NU_LABEL_KIND)
        end do
        !$omp end parallel do
        deallocate(sel)
        ! solvent-constraint clamp: outside the envelope the label is the
        ! coarsest candidate. Applied before Potts smoothing as the intended
        ! initialization AND re-applied after: the smoothing re-optimizes on
        ! the unaries at every sweep, so a pre-smoothing override alone is
        ! eroded wherever solvent unaries prefer fine labels
        if( nu_l_solvent_clamp )then
            call apply_nu_solvent_clamp(n_clamped)
            if( nu_l_report ) write(logfhandle,'(A,I0,A)') &
                &'>>> NU SOLVENT CLAMP: ', n_clamped, ' support voxels outside the envelope set to the coarsest candidate'
        endif
        if( NU_DEV_OUTPUT .and. nu_l_report ) call log_nu_aux_replacement_margin_stats()
        if( NU_DEV_OUTPUT .and. nu_l_report ) &
            &call log_nu_candidate_selection_counts(filtmap, n_base, 'before ordered-label smoothing')
        call refine_nu_candidate_map_ordered_labels(filtmap, n_candidates)
        if( nu_l_solvent_clamp ) call apply_nu_solvent_clamp()
        if( NU_DEV_OUTPUT .and. nu_l_report ) &
            &call log_nu_candidate_selection_counts(filtmap, n_base, 'after ordered-label smoothing')
        call clamp_nu_filtmap_labels(n_base)
        call cache_nu_extension_frontier_dmats(filtmap, n_base)
        ! Keep the mask-packed unary bank. High-resolution extension appends
        ! accepted challenger unaries and can then run a final ordered-label
        ! cleanup over the expanded label field.
    end subroutine optimize_nu_cutoff_finds

    !> Enforce the solvent-constraint clamp on the module label field:
    !! support voxels outside the envelope take the coarsest candidate.
    !! Called before AND after every ordered-label smoothing pass, because
    !! the smoothing re-optimizes labels from the unaries.
    module subroutine apply_nu_solvent_clamp( n_clamped )
        integer, optional, intent(out) :: n_clamped
        integer :: i, j, k, imask, n_here
        if( .not. nu_l_solvent_clamp ) then
            if( present(n_clamped) ) n_clamped = 0
            return
        endif
        if( .not. allocated(filtmap) ) THROW_HARD('filtmap not allocated; apply_nu_solvent_clamp')
        n_here = 0
        !$omp parallel do schedule(static) default(shared) private(i,j,k,imask) &
        !$omp reduction(+:n_here) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( nu_solvent_lmask(i,j,k) )then
                if( filtmap(i,j,k) /= 1_NU_LABEL_KIND )then
                    filtmap(i,j,k) = 1_NU_LABEL_KIND
                    n_here = n_here + 1
                endif
            endif
        end do
        !$omp end parallel do
        if( present(n_clamped) ) n_clamped = n_here
    end subroutine apply_nu_solvent_clamp

    subroutine cache_nu_extension_frontier_dmats( candmap, n_base )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer, intent(in) :: n_base
        integer :: i, j, k, imask, icand
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(n_nu_mask), source=huge(0.))
        !$omp parallel do schedule(static) default(shared) private(i,j,k,imask,icand) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            icand = int(candmap(i,j,k))
            if( nu_effective_base_label_for_candidate(icand, n_base) /= n_base ) cycle
            dmat_finest_cached(imask) = dmats_mask(imask,icand)
        end do
        !$omp end parallel do
    end subroutine cache_nu_extension_frontier_dmats

    module subroutine clamp_nu_filtmap_labels( n_base )
        integer, intent(in) :: n_base
        integer :: i, j, k, icand, imask
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; clamp_nu_filtmap_labels')
        if( n_base < 1 ) THROW_HARD('empty base bank; clamp_nu_filtmap_labels')
        !$omp parallel do schedule(static) default(shared) private(i,j,k,icand,imask) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            icand = int(filtmap(i,j,k))
            filtmap(i,j,k) = int(max(1, min(n_base, icand)), kind=NU_LABEL_KIND)
        end do
        !$omp end parallel do
    end subroutine clamp_nu_filtmap_labels

    module subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
        integer(kind=NU_LABEL_KIND), intent(in) :: candmap(:,:,:)
        integer,          intent(in) :: n_base
        character(len=*), intent(in) :: stage
        integer, allocatable :: cand_counts(:)
        integer :: icand, n_candidates, nmask, nvox
        integer :: i, j, k, imask
        real    :: pct
        character(len=16) :: source_tag
        if( .not.allocated(nu_lmask) ) return
        if( .not.allocated(candidate_coords) ) return
        n_candidates = size(candidate_coords)
        if( n_candidates /= n_base ) return
        nmask = n_nu_mask
        if( nmask == 0 ) return
        allocate(cand_counts(n_candidates), source=0)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,icand) reduction(+:cand_counts) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            icand = int(candmap(i,j,k))
            if( icand >= 1 .and. icand <= n_candidates ) cand_counts(icand) = cand_counts(icand) + 1
        end do
        !$omp end parallel do
        write(logfhandle,'(A,A)') '>>> NU candidate label assignments ', trim(stage)
        write(logfhandle,'(A,I12)') '    Mask voxels: ', nmask
        write(logfhandle,'(A)') '    Source      Bank  Coord      LP(A)        Voxels    Pct mask'
        do icand = 1, n_candidates
            nvox = cand_counts(icand)
            pct = 100. * real(nvox) / real(nmask)
            if( nu_label_is_aux_replacement(icand) )then
                source_tag = 'AuxReplace'
            else
                source_tag = 'Base'
            endif
            write(logfhandle,'(4X,A10,2X,I4,2X,F7.2,2X,F9.3,2X,I12,2X,F8.2,A)') &
                &source_tag, icand, candidate_coords(icand), nu_label_lowpass_limit(icand), nvox, pct, '%'
        end do
        deallocate(cand_counts)
    end subroutine log_nu_candidate_selection_counts

    subroutine log_nu_aux_replacement_margin_stats()
        integer :: ibase, n_base, nmask, nwins
        integer :: imask
        real    :: best_base, margin, margin_sum, win_margin_sum, avg_margin, avg_win_margin, pct
        if( .not.allocated(dmats_mask) ) return
        if( .not.allocated(nu_lmask) ) return
        if( .not.allocated(candidate_coords) ) return
        if( .not.allocated(nu_mask_vox) ) return
        n_base = size(cutoff_finds)
        if( n_base < 2 ) return
        if( .not.nu_label_is_aux_replacement(n_base) ) return
        nmask = n_nu_mask
        if( nmask == 0 ) return
        nwins = 0
        margin_sum = 0.
        win_margin_sum = 0.
        !$omp parallel do schedule(static) default(shared) &
        !$omp private(imask,ibase,best_base,margin) reduction(+:nwins,margin_sum,win_margin_sum) proc_bind(close)
        do imask = 1, n_nu_mask
            best_base = dmats_mask(imask,1)
            do ibase = 2, n_base - 1
                best_base = min(best_base, dmats_mask(imask,ibase))
            end do
            margin = best_base - dmats_mask(imask,n_base)
            margin_sum = margin_sum + margin
            if( margin > 0. )then
                nwins = nwins + 1
                win_margin_sum = win_margin_sum + margin
            endif
        end do
        !$omp end parallel do
        avg_margin = margin_sum / real(nmask)
        avg_win_margin = 0.
        if( nwins > 0 ) avg_win_margin = win_margin_sum / real(nwins)
        pct = 100. * real(nwins) / real(nmask)
        write(logfhandle,'(A)') '>>> NU auxiliary replacement unary margins versus best coarser retained label'
        write(logfhandle,'(A)') '    Positive margin means the replacement label has the lower unary objective.'
        write(logfhandle,'(A,F8.3,A,I8,2X,F8.2,A,2X,F12.5,2X,F12.5)') &
            &'    Replacement LP(A): ', nu_aux_replacement_resolution, '; wins: ', nwins, pct, '%', &
            &avg_margin, avg_win_margin
    end subroutine log_nu_aux_replacement_margin_stats

    module subroutine log_nu_candidate_coords
        integer :: icand
        if( .not.allocated(candidate_coords) ) return
        write(logfhandle,'(A)', advance='no') '>>> NU ordered-label smoothing candidate coordinates:'
        do icand = 1, size(candidate_coords)
            write(logfhandle,'(1X,F6.2)', advance='no') candidate_coords(icand)
        end do
        write(logfhandle,*)
    end subroutine log_nu_candidate_coords

    subroutine log_nu_objective_smoothing_bank()
        integer :: i
        real    :: lp_angstrom, radius_angstrom
        character(len=10) :: source_tag
        write(logfhandle,'(A,F6.2,A,F6.2,A,F7.2,A)') &
            &'>>> NU objective AWF smoothing: radius_A=', NU_OBJECTIVE_SMOOTH_RADIUS_FRAC, &
            &' * AWF * LP(A), AWF=', NU_OBJECTIVE_SMOOTH_AWF, &
            &', cap=', NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A, ' A'
        write(logfhandle,'(A)') '    Source  Bank  Fourier k    LP(A)  Radius(A)  Radius(px)'
        do i = 1, size(cutoff_finds)
            lp_angstrom = nu_label_lowpass_limit(i)
            radius_angstrom = nu_objective_smooth_radius_angstrom(lp_angstrom)
            if( nu_label_is_aux_replacement(i) )then
                source_tag = 'AuxReplace'
            else
                source_tag = 'Base'
            endif
            write(logfhandle,'(4X,A10,2X,I4,2X,I9,2X,F7.3,2X,F9.3,2X,I10)') &
                &source_tag, i, cutoff_finds(i), lp_angstrom, radius_angstrom, &
                &nu_objective_smooth_radius_pixels(lp_angstrom)
        end do
    end subroutine log_nu_objective_smoothing_bank

    module real function nu_candidate_coord_for_label( ilabel )
        integer, intent(in) :: ilabel
        if( allocated(candidate_coords) )then
            if( ilabel >= 1 .and. ilabel <= size(candidate_coords) )then
                nu_candidate_coord_for_label = candidate_coords(ilabel)
                return
            endif
        endif
        nu_candidate_coord_for_label = real(ilabel)
    end function nu_candidate_coord_for_label

    module integer function nu_effective_base_label_for_candidate( icand, n_base )
        integer, intent(in) :: icand, n_base
        if( n_base < 1 ) THROW_HARD('empty base bank; nu_effective_base_label_for_candidate')
        nu_effective_base_label_for_candidate = nint(nu_candidate_coord_for_label(icand))
        nu_effective_base_label_for_candidate = max(1, min(n_base, nu_effective_base_label_for_candidate))
    end function nu_effective_base_label_for_candidate

end submodule simple_nu_filter_bank
