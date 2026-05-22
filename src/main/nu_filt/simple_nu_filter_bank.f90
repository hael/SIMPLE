!@descr: simple nu filter bank implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_bank
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine setup_nu_dmats( vol_even, vol_odd, l_mask, aux_resolutions, aux_even, aux_odd, &
            &n_highres_steps, l_aux_source_unordered )
        class(image),          intent(in) :: vol_even, vol_odd
        logical,               intent(in) :: l_mask(:,:,:)
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
        logical,     optional, intent(in) :: l_aux_source_unordered
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:), dmat_cand(:,:,:)
        real :: noise_sigma
        integer :: i, n_candidates
        real    :: x
        call init_nu_filter(vol_even, vol_odd, n_highres_steps)
        l_aux_source_unordered_potts = .false.
        if( present(l_aux_source_unordered) ) l_aux_source_unordered_potts = l_aux_source_unordered
        if( any(shape(l_mask) /= ldim) ) THROW_HARD('l_mask shape mismatch in setup_nu_dmats')
        if( allocated(nu_lmask) ) deallocate(nu_lmask)
        allocate(nu_lmask(ldim(1),ldim(2),ldim(3)), source=l_mask)
        if( .not. any(nu_lmask) ) THROW_HARD('l_mask has no true voxels in setup_nu_dmats')
        call setup_nu_mask_index
        if( present(aux_even) ) then
            if( .not. present(aux_odd) ) THROW_HARD('Auxiliary odd bank missing; setup_nu_dmats')
            if( size(aux_resolutions) /= size(aux_even) ) THROW_HARD('Auxiliary resolutions size mismatch; setup_nu_dmats')
            call stash_aux_volumes(aux_even, aux_odd)
        else
            if( size(aux_resolutions) /= 0 ) THROW_HARD('Auxiliary resolutions supplied without auxiliary volumes; setup_nu_dmats')
            call cleanup_aux_bank
        end if
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call cache_filtered_vols(vol_even, vol_odd)
        noise_sigma = vol_even%nu_objective_noise_scale(vol_odd, nu_lmask)
        write(logfhandle,'(A,ES12.4)') '>>> NU normalized Huber objective raw E/O scale: ', noise_sigma
        if( allocated(dmats_mask) ) deallocate(dmats_mask)
        n_candidates = size(cutoff_finds)
        if( allocated(aux_even_bank) ) n_candidates = n_candidates + size(aux_even_bank)
        call setup_nu_candidate_coords(n_candidates, aux_resolutions)
        call log_nu_objective_smoothing_bank(aux_resolutions)
        allocate(dmats_mask(n_nu_mask,n_candidates), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),  source=0.)
        allocate(dmat_cand(ldim(1),ldim(2),ldim(3)), source=huge(x))
        do i = 1, size(cutoff_finds)
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
            if( .not.file_exists(even_cache_fname) ) THROW_HARD('Missing filtered volume cache: '//even_cache_fname%to_char())
            if( .not.file_exists(odd_cache_fname)  ) THROW_HARD('Missing filtered volume cache: '//odd_cache_fname%to_char())
            call vol_even_filt%read(even_cache_fname)
            call vol_odd_filt%read(odd_cache_fname)
            dmat_cand = huge(x)
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmat_cand, &
                &nu_lmask, noise_sigma)
            call smooth_nu_objective(dmat_cand, dmat_tmp, cutoff_find_to_lowpass_limit(i))
            call pack_nu_dmat_candidate(dmat_cand, i)
        end do
        if( allocated(aux_even_bank) ) then
            do i = 1, size(aux_even_bank)
                dmat_cand = huge(x)
                call vol_even%nu_objective(aux_even_bank(i), vol_odd, aux_odd_bank(i), &
                    &dmat_cand, nu_lmask, noise_sigma)
                call smooth_nu_objective(dmat_cand, dmat_tmp, aux_resolutions(i))
                call pack_nu_dmat_candidate(dmat_cand, size(cutoff_finds)+i)
            end do
        end if
        call vol_even_filt%kill
        call vol_odd_filt%kill
        deallocate(dmat_tmp, dmat_cand)
    end subroutine setup_nu_dmats

    module subroutine setup_nu_candidate_coords( n_candidates, aux_resolutions )
        integer, intent(in) :: n_candidates
        real,    intent(in) :: aux_resolutions(:)
        integer :: i, n_base, iaux
        n_base = size(cutoff_finds)
        if( allocated(candidate_coords) ) deallocate(candidate_coords)
        allocate(candidate_coords(n_candidates), source=0.)
        do i = 1, n_base
            candidate_coords(i) = real(i)
        end do
        if( n_candidates == n_base ) return
        if( size(aux_resolutions) /= n_candidates - n_base )then
            THROW_HARD('aux_resolutions size mismatch in setup_nu_candidate_coords')
        endif
        do iaux = 1, size(aux_resolutions)
            candidate_coords(n_base + iaux) = lowpass_limit_to_candidate_coord(aux_resolutions(iaux))
        end do
    end subroutine setup_nu_candidate_coords

    module real function lowpass_limit_to_candidate_coord( lp_angstrom )
        real, intent(in) :: lp_angstrom
        integer :: i, n_base
        real :: denom, lp_hi, lp_lo
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; lowpass_limit_to_candidate_coord')
        n_base = size(cutoff_finds)
        lp_hi = cutoff_find_to_lowpass_limit(1)
        if( lp_angstrom >= lp_hi )then
            lowpass_limit_to_candidate_coord = 1.
            return
        endif
        lp_lo = cutoff_find_to_lowpass_limit(n_base)
        if( lp_angstrom <= lp_lo )then
            lowpass_limit_to_candidate_coord = real(n_base)
            return
        endif
        do i = 1, n_base - 1
            lp_hi = cutoff_find_to_lowpass_limit(i)
            lp_lo = cutoff_find_to_lowpass_limit(i+1)
            if( lp_angstrom <= lp_hi .and. lp_angstrom >= lp_lo )then
                denom = lp_hi - lp_lo
                if( denom <= TINY ) cycle
                lowpass_limit_to_candidate_coord = real(i) + (lp_hi - lp_angstrom) / denom
                return
            endif
        end do
        lowpass_limit_to_candidate_coord = real(n_base)
    end function lowpass_limit_to_candidate_coord

    module real function get_nu_filter_bank_finest_lp()
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; get_nu_filter_bank_finest_lp')
        if( size(cutoff_finds) < 1 ) THROW_HARD('empty filter bank; get_nu_filter_bank_finest_lp')
        get_nu_filter_bank_finest_lp = cutoff_find_to_lowpass_limit(size(cutoff_finds))
    end function get_nu_filter_bank_finest_lp

    module subroutine optimize_nu_cutoff_finds()
        integer, allocatable :: candmap(:,:,:)
        integer :: nx, ny, nz, i, j, k, icand, best_icand, n_base, n_candidates, imask
        real    :: best_dmat
        if( .not.allocated(dmats_mask) ) THROW_HARD('dmats_mask not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        if( .not.allocated(nu_mask_index) ) THROW_HARD('nu_mask_index not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        n_base       = size(cutoff_finds)
        ! dmats_mask is laid out as [base low-pass bank | auxiliary pre-filtered pairs].
        ! The first n_base entries correspond to cutoff_finds(:); any remaining
        ! entries are caller-supplied auxiliary sources.
        n_candidates = size(dmats_mask, 2)
        if( allocated(dmats_aux_mask) ) deallocate(dmats_aux_mask)
        if( l_aux_source_unordered_potts .and. n_candidates > n_base )then
            allocate(dmats_aux_mask(n_nu_mask,n_candidates-n_base))
            dmats_aux_mask = dmats_mask(:,n_base+1:n_candidates)
        endif
        if( allocated(filtmap) ) deallocate(filtmap)
        if( allocated(srcmap)  ) deallocate(srcmap)
        allocate(candmap(nx,ny,nz), source=1)
        allocate(filtmap(nx,ny,nz), source=1)
        allocate(srcmap(nx,ny,nz),  source=1)
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,icand,best_icand,best_dmat,imask) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not. nu_lmask(i,j,k) )then
                        filtmap(i,j,k) = 1
                        srcmap(i,j,k)  = 1
                        cycle
                    endif
                    imask = nu_mask_index(i,j,k)
                    best_icand = 1
                    best_dmat = dmats_mask(imask,1)
                    do icand = 2, n_candidates
                        if( dmats_mask(imask,icand) < best_dmat ) then
                            best_dmat = dmats_mask(imask,icand)
                            best_icand = icand
                        end if
                    end do
                    candmap(i,j,k) = best_icand
                end do
            end do
        end do
        !$omp end parallel do
        call log_nu_aux_unary_margin_stats(n_base)
        call log_nu_candidate_selection_counts(candmap, n_base, 'before ordered-label smoothing')
        call refine_nu_candidate_map_ordered_labels(candmap, n_candidates)
        call log_nu_candidate_selection_counts(candmap, n_base, 'after ordered-label smoothing')
        call candidate_map_to_filt_and_src(candmap, n_base)
        call cache_nu_extension_frontier_dmats(candmap, n_base)
        ! Keep the mask-packed unary bank. High-resolution extension appends
        ! accepted challenger unaries and can then run a final ordered-label
        ! cleanup over the expanded label field.
        deallocate(candmap)
    end subroutine optimize_nu_cutoff_finds

    subroutine cache_nu_extension_frontier_dmats( candmap, n_base )
        integer, intent(in) :: candmap(:,:,:), n_base
        integer :: i, j, k, imask, icand, nx, ny, nz
        nx = size(candmap, 1)
        ny = size(candmap, 2)
        nz = size(candmap, 3)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(nx,ny,nz), source=huge(0.))
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,imask,icand) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not.nu_lmask(i,j,k) ) cycle
                    icand = candmap(i,j,k)
                    if( nu_effective_base_label_for_candidate(icand, n_base) /= n_base ) cycle
                    imask = nu_mask_index(i,j,k)
                    dmat_finest_cached(i,j,k) = dmats_mask(imask,icand)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine cache_nu_extension_frontier_dmats

    module subroutine candidate_map_to_filt_and_src( candmap, n_base )
        integer, intent(in) :: candmap(:,:,:), n_base
        integer :: i, j, k, icand, nx, ny, nz
        nx = size(candmap, 1)
        ny = size(candmap, 2)
        nz = size(candmap, 3)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,icand) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    icand = candmap(i,j,k)
                    ! Base-bank winners preserve their low-pass index in filtmap.
                    ! Auxiliary winners preserve their provenance in srcmap while
                    ! filtmap stores the nearest effective base-bank label for
                    ! diagnostics and map visualization.
                    if( icand <= n_base ) then
                        srcmap(i,j,k)  = 1
                        filtmap(i,j,k) = icand
                    else
                        ! srcmap numbering:
                        !   1   -> base low-pass bank
                        !   2+  -> auxiliary pair 1, 2, ...
                        srcmap(i,j,k)  = icand - n_base + 1
                        filtmap(i,j,k) = nu_effective_base_label_for_candidate(icand, n_base)
                    end if
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine candidate_map_to_filt_and_src

    module subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
        integer,          intent(in) :: candmap(:,:,:), n_base
        character(len=*), intent(in) :: stage
        integer :: icand, iaux, n_candidates, nmask, nvox, nbasevox, nauxvox, eff_label
        real    :: pct
        character(len=16) :: auxtag
        if( .not.allocated(nu_lmask) ) return
        if( .not.allocated(candidate_coords) ) return
        n_candidates = size(candidate_coords)
        if( n_candidates <= n_base ) return
        nmask = count(nu_lmask)
        if( nmask == 0 ) return
        nbasevox = 0
        do icand = 1, n_base
            nbasevox = nbasevox + count(nu_lmask .and. candmap == icand)
        end do
        nauxvox = 0
        do icand = n_base + 1, n_candidates
            nauxvox = nauxvox + count(nu_lmask .and. candmap == icand)
        end do
        write(logfhandle,'(A,A)') '>>> NU candidate source assignments ', trim(stage)
        write(logfhandle,'(A,I12)') '    Mask voxels:      ', nmask
        pct = 100. * real(nbasevox) / real(nmask)
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Base-bank voxels: ', nbasevox, ' (', pct, '%)'
        pct = 100. * real(nauxvox) / real(nmask)
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Auxiliary voxels: ', nauxvox, ' (', pct, '%)'
        write(logfhandle,'(A)') '    Source      Coord  Nearest LP(A)        Voxels    Pct mask'
        do iaux = 1, n_candidates - n_base
            icand = n_base + iaux
            eff_label = nu_effective_base_label_for_candidate(icand, n_base)
            nvox = count(nu_lmask .and. candmap == icand)
            pct = 100. * real(nvox) / real(nmask)
            write(auxtag,'(A,I0)') 'Aux', iaux
            write(logfhandle,'(4X,A8,2X,F7.2,2X,F13.3,2X,I12,2X,F8.2,A)') &
                &auxtag, candidate_coords(icand), cutoff_find_to_lowpass_limit(eff_label), nvox, pct, '%'
        end do
    end subroutine log_nu_candidate_selection_counts

    module subroutine log_nu_aux_unary_margin_stats( n_base )
        integer, intent(in) :: n_base
        integer :: iaux, ibase, icand, n_candidates, nmask, nwins, eff_label
        integer :: i, j, k, imask
        real    :: best_base, margin, margin_sum, win_margin_sum, avg_margin, avg_win_margin, pct
        character(len=16) :: auxtag
        if( .not.allocated(dmats_mask) ) return
        if( .not.allocated(nu_lmask) ) return
        if( .not.allocated(candidate_coords) ) return
        if( .not.allocated(nu_mask_index) ) return
        n_candidates = size(dmats_mask, 2)
        if( n_candidates <= n_base ) return
        nmask = count(nu_lmask)
        if( nmask == 0 ) return
        write(logfhandle,'(A)') '>>> NU auxiliary unary margins versus best base-bank candidate'
        write(logfhandle,'(A)') '    Positive margin means the auxiliary candidate has the lower unary objective.'
        write(logfhandle,'(A)') '    Source      Coord  Nearest LP(A)      Wins    Pct mask   Mean margin    Win margin'
        do iaux = 1, n_candidates - n_base
            icand = n_base + iaux
            eff_label = nu_effective_base_label_for_candidate(icand, n_base)
            nwins = 0
            margin_sum = 0.
            win_margin_sum = 0.
            !$omp parallel do collapse(3) schedule(static) default(shared) &
            !$omp private(i,j,k,imask,ibase,best_base,margin) reduction(+:nwins,margin_sum,win_margin_sum) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.nu_lmask(i,j,k) ) cycle
                        imask = nu_mask_index(i,j,k)
                        best_base = dmats_mask(imask,1)
                        do ibase = 2, n_base
                            best_base = min(best_base, dmats_mask(imask,ibase))
                        end do
                        margin = best_base - dmats_mask(imask,icand)
                        margin_sum = margin_sum + margin
                        if( margin > 0. )then
                            nwins = nwins + 1
                            win_margin_sum = win_margin_sum + margin
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
            avg_margin = margin_sum / real(nmask)
            avg_win_margin = 0.
            if( nwins > 0 ) avg_win_margin = win_margin_sum / real(nwins)
            pct = 100. * real(nwins) / real(nmask)
            write(auxtag,'(A,I0)') 'Aux', iaux
            write(logfhandle,'(4X,A8,2X,F7.2,2X,F13.3,2X,I8,2X,F8.2,A,2X,F12.5,2X,F12.5)') &
                &auxtag, candidate_coords(icand), cutoff_find_to_lowpass_limit(eff_label), &
                &nwins, pct, '%', avg_margin, avg_win_margin
        end do
    end subroutine log_nu_aux_unary_margin_stats

    module subroutine log_nu_candidate_coords
        integer :: icand
        if( .not.allocated(candidate_coords) ) return
        write(logfhandle,'(A)', advance='no') '>>> NU ordered-label smoothing candidate coordinates:'
        do icand = 1, size(candidate_coords)
            write(logfhandle,'(1X,F6.2)', advance='no') candidate_coords(icand)
        end do
        write(logfhandle,*)
    end subroutine log_nu_candidate_coords

    subroutine log_nu_objective_smoothing_bank( aux_resolutions )
        real, intent(in) :: aux_resolutions(:)
        integer :: i
        real    :: lp_angstrom, radius_angstrom
        write(logfhandle,'(A,F6.2,A,F6.2,A,F7.2,A)') &
            &'>>> NU objective AWF smoothing: radius_A=', NU_OBJECTIVE_SMOOTH_RADIUS_FRAC, &
            &' * AWF * LP(A), AWF=', NU_OBJECTIVE_SMOOTH_AWF, &
            &', cap=', NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A, ' A'
        write(logfhandle,'(A)') '    Source  Bank  Fourier k    LP(A)  Radius(A)  Radius(px)'
        do i = 1, size(cutoff_finds)
            lp_angstrom = cutoff_find_to_lowpass_limit(i)
            radius_angstrom = nu_objective_smooth_radius_angstrom(lp_angstrom)
            write(logfhandle,'(4X,A6,2X,I4,2X,I9,2X,F7.3,2X,F9.3,2X,I10)') &
                &'Base', i, cutoff_finds(i), lp_angstrom, radius_angstrom, &
                &nu_objective_smooth_radius_pixels(lp_angstrom)
        end do
        do i = 1, size(aux_resolutions)
            lp_angstrom = aux_resolutions(i)
            radius_angstrom = nu_objective_smooth_radius_angstrom(lp_angstrom)
            write(logfhandle,'(4X,A6,2X,I4,2X,I9,2X,F7.3,2X,F9.3,2X,I10)') &
                &'Aux', size(cutoff_finds) + i, 0, lp_angstrom, radius_angstrom, &
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
