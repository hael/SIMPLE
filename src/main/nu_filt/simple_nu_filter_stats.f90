!@descr: simple nu filter stats implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_stats
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
        real, allocatable, intent(inout) :: lowpass_vals(:)
        logical, optional, intent(in)    :: mask(:,:,:)
        integer :: i, j, k, imask, nvals, ival
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before pack_filtmap_lowpass_limits')
        call require_valid_stats_mask(mask, 'pack_filtmap_lowpass_limits')
        nvals = count_active_nu_mask(mask, 1)
        if( allocated(lowpass_vals) ) deallocate(lowpass_vals)
        allocate(lowpass_vals(nvals), source=0.)
        ival = 0
        if( present(mask) )then
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                        if( allocated(srcmap) ) then
                            if( srcmap(i,j,k) /= 1 ) cycle
                        end if
                        ival = ival + 1
                        lowpass_vals(ival) = cutoff_find_to_lowpass_limit(int(filtmap(i,j,k)))
                    end do
                end do
            end do
        else
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( allocated(srcmap) ) then
                    if( srcmap(i,j,k) /= 1 ) cycle
                end if
                ival = ival + 1
                lowpass_vals(ival) = cutoff_find_to_lowpass_limit(int(filtmap(i,j,k)))
            end do
        endif
    end subroutine pack_filtmap_lowpass_limits

    module subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, optional, intent(in)   :: mask(:,:,:)
        integer, allocatable :: counts(:)
        real,    allocatable :: percentages(:)
        integer :: icut, nselected, pos1, pos2
        real(dp) :: lp_dp, nr_dp, sum_dp, sumsq_dp, var_dp
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_stats')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before calc_filtmap_lowpass_stats')
        allocate(counts(size(cutoff_finds)), percentages(size(cutoff_finds)))
        call calc_filtmap_lowpass_histogram(counts, percentages, mask)
        nselected = sum(counts)
        if( nselected == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_stats')
        sum_dp = 0._dp
        sumsq_dp = 0._dp
        do icut = 1, size(cutoff_finds)
            if( counts(icut) == 0 ) cycle
            lp_dp = real(cutoff_find_to_lowpass_limit(icut), dp)
            sum_dp   = sum_dp   + real(counts(icut), dp) * lp_dp
            sumsq_dp = sumsq_dp + real(counts(icut), dp) * lp_dp * lp_dp
        end do
        nr_dp = real(nselected, dp)
        statvars%avg = real(sum_dp / nr_dp)
        var_dp = 0._dp
        if( nselected > 1 ) var_dp = max(0._dp, (sumsq_dp - sum_dp * sum_dp / nr_dp) / (nr_dp - 1._dp))
        statvars%sdev = real(sqrt(var_dp))
        pos1 = (nselected + 1) / 2
        pos2 = (nselected + 2) / 2
        statvars%med = 0.5 * (lowpass_histogram_value_at(counts, pos1) + &
            &lowpass_histogram_value_at(counts, pos2))
        statvars%minv = lowpass_histogram_value_at(counts, 1)
        statvars%maxv = lowpass_histogram_value_at(counts, nselected)
        deallocate(counts, percentages)
    end subroutine calc_filtmap_lowpass_stats

    module subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
        integer, intent(out) :: counts(:)
        real,    intent(out) :: percentages(:)
        logical, optional, intent(in) :: mask(:,:,:)
        integer :: i, j, k, imask, icut, nselected
        if( .not.allocated(filtmap)                 ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_histogram')
        if( .not.allocated(cutoff_finds)            ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before calc_filtmap_lowpass_histogram')
        if( size(counts) /= size(cutoff_finds)      ) THROW_HARD('counts size mismatch in calc_filtmap_lowpass_histogram')
        if( size(percentages) /= size(cutoff_finds) ) THROW_HARD('percentages size mismatch in calc_filtmap_lowpass_histogram')
        call require_valid_stats_mask(mask, 'calc_filtmap_lowpass_histogram')
        nselected = count_active_nu_mask(mask, 1)
        counts       = 0
        percentages  = 0.
        if( nselected == 0 ) return
        if( present(mask) )then
            !$omp parallel do collapse(3) schedule(static) default(shared) &
            !$omp private(i,j,k,icut) reduction(+:counts) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                        if( allocated(srcmap) )then
                            if( srcmap(i,j,k) /= 1 ) cycle
                        endif
                        icut = int(filtmap(i,j,k))
                        if( icut < 1 .or. icut > size(cutoff_finds) ) cycle
                        counts(icut) = counts(icut) + 1
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) &
            !$omp private(imask,i,j,k,icut) reduction(+:counts) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( allocated(srcmap) )then
                    if( srcmap(i,j,k) /= 1 ) cycle
                endif
                icut = int(filtmap(i,j,k))
                if( icut < 1 .or. icut > size(cutoff_finds) ) cycle
                counts(icut) = counts(icut) + 1
            end do
            !$omp end parallel do
        endif
        do icut = 1, size(cutoff_finds)
            percentages(icut) = 100. * real(counts(icut)) / real(nselected)
        end do
    end subroutine calc_filtmap_lowpass_histogram

    module real function get_nu_filtmap_finest_selected_lp( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
        integer :: icut, iaux, ncur
        real    :: selected_lp
        if( .not.allocated(filtmap) )then
            THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before get_nu_filtmap_finest_selected_lp')
        endif
        if( .not.allocated(cutoff_finds) )then
            THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before get_nu_filtmap_finest_selected_lp')
        endif
        call require_valid_stats_mask(mask, 'get_nu_filtmap_finest_selected_lp')
        get_nu_filtmap_finest_selected_lp = 0.
        do icut = size(cutoff_finds), 1, -1
            ncur = count_active_nu_label(mask, icut, 1)
            if( ncur == 0 ) cycle
            get_nu_filtmap_finest_selected_lp = cutoff_find_to_lowpass_limit(icut)
            exit
        end do
        if( present(aux_resolutions) )then
            if( .not. allocated(srcmap) ) return
            if( allocated(aux_even_bank) )then
                if( size(aux_resolutions) /= size(aux_even_bank) ) &
                    &THROW_HARD('aux_resolutions size mismatch in get_nu_filtmap_finest_selected_lp')
            endif
            do iaux = 1, size(aux_resolutions)
                ncur = count_active_nu_mask(mask, iaux + 1)
                if( ncur == 0 ) cycle
                selected_lp = aux_resolutions(iaux)
                if( selected_lp <= TINY ) cycle
                if( get_nu_filtmap_finest_selected_lp <= TINY .or. &
                    &selected_lp < get_nu_filtmap_finest_selected_lp )then
                    get_nu_filtmap_finest_selected_lp = selected_lp
                endif
            enddo
        endif
    end function get_nu_filtmap_finest_selected_lp

    module subroutine print_filtmap_lowpass_histogram( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
        integer, allocatable :: counts(:)
        real,    allocatable :: percentages(:)
        integer :: icut, iaux, nselected, nvox
        real    :: pct
        character(len=8) :: auxtag
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before print_filtmap_lowpass_histogram')
        allocate(counts(size(cutoff_finds)), percentages(size(cutoff_finds)))
        call calc_filtmap_lowpass_histogram(counts, percentages, mask)
        nselected = count_active_nu_mask(mask, 1)
        write(logfhandle,'(A)') '>>> NU LOW-PASS ASSIGNMENTS (base filter bank)'
        write(logfhandle,'(A,I12)') '    Base-bank voxels: ', nselected
        write(logfhandle,'(A)')     '    Bank  Fourier k  LP limit (A)        Voxels    Pct base'
        do icut = 1, size(cutoff_finds)
            write(logfhandle,'(4X,I4,2X,I9,2X,F12.3,2X,I12,2X,F8.2,A)') &
                &icut, cutoff_finds(icut), cutoff_find_to_lowpass_limit(icut), counts(icut), percentages(icut), '%'
        end do
        ! Print auxiliary pair assignments if present
        if( allocated(aux_even_bank) .and. present(aux_resolutions) ) then
            if( size(aux_resolutions) /= size(aux_even_bank) ) THROW_HARD('aux_resolutions size mismatch in print_filtmap_lowpass_histogram')
            nselected = count_active_nu_mask(mask, 0)
            write(logfhandle,'(A)') ''
            write(logfhandle,'(A)')     '>>> NU AUXILIARY SOURCE ASSIGNMENTS'
            write(logfhandle,'(A,I12)') '    Mask voxels:      ', nselected
            write(logfhandle,'(A)')     '    Source    Resolution (A)        Voxels    Pct mask'
            do iaux = 1, size(aux_even_bank)
                nvox = count_active_nu_mask(mask, iaux + 1)
                pct = 0.
                if( nselected > 0 ) pct = 100. * real(nvox) / real(nselected)
                write(auxtag,'(A,I0,A)') 'Aux', iaux, '@'
                write(logfhandle,'(4X,A8,2X,F14.3,2X,I12,2X,F8.2,A)') auxtag, aux_resolutions(iaux), nvox, pct, '%'
            end do
        end if
        deallocate(counts, percentages)
    end subroutine print_filtmap_lowpass_histogram

    module subroutine print_nu_filtmap_lowpass_stats( mask, aux_resolutions )
        logical, optional, intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
        type(stats_struct) :: statvars
        integer :: nbase
        if( allocated(srcmap) ) then
            nbase = count_active_nu_mask(mask, 1)
            if( nbase == 0 ) then
                write(logfhandle,'(A)') '>>> No base low-pass selections remain after auxiliary-source optimization'
                call print_filtmap_lowpass_histogram(mask, aux_resolutions)
                return
            end if
        end if
        call calc_filtmap_lowpass_stats(statvars, mask)
        nbase = count_active_nu_mask(mask, 1)
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NU FILTER LOCAL RESOLUTION SUMMARY'
        write(logfhandle,'(A,I12)') '    Voxels analyzed: ', nbase
        write(logfhandle,'(A)')     '              Mean    Median     Sigma       Min       Max'
        write(logfhandle,'(A,5F10.3)') '    Angstrom ', statvars%avg, statvars%med, statvars%sdev, statvars%minv, statvars%maxv
        write(logfhandle,'(A)') ''
        call print_filtmap_lowpass_histogram(mask, aux_resolutions)
    end subroutine print_nu_filtmap_lowpass_stats

    module subroutine analyze_filtmap_neighbor_continuity( mask )
        logical, optional, intent(in) :: mask(:,:,:)
        integer :: i, j, k, di, dj, dk, ni, nj, nk
        integer :: lp_i, lp_j, step_diff, n_discontinuous_neighbors
        integer :: n_total_neighbor_pairs, n_discontinuous_pairs, n_voxels_with_discontinuity, nx, ny, nz, ii
        integer :: thresh, large_thresh, n_analyzed, n_identical_pairs, n_one_step_pairs, n_large_jump_pairs
        integer :: max_step_diff
        integer, allocatable :: stepdiff_counts(:)
        real :: coord_i, coord_j, pct, pct_vox, pct_pairs, pct_large
        logical :: l_count_pair
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before analyze_filtmap_neighbor_continuity')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before analyze_filtmap_neighbor_continuity')
        call require_valid_stats_mask(mask, 'analyze_filtmap_neighbor_continuity')
        thresh = DISCONT_STEP_THRESH
        large_thresh = NU_CONTINUITY_LARGE_STEP_THRESH
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        max_step_diff = max(1, size(cutoff_finds) - 1)
        allocate(stepdiff_counts(max_step_diff), source=0)
        n_voxels_with_discontinuity = 0
        n_total_neighbor_pairs      = 0
        n_discontinuous_pairs       = 0
        n_identical_pairs           = 0
        n_one_step_pairs            = 0
        n_large_jump_pairs          = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,di,dj,dk,ni,nj,nk,lp_i,lp_j,step_diff,n_discontinuous_neighbors, &
        !$omp coord_i,coord_j,l_count_pair) &
        !$omp reduction(+:n_total_neighbor_pairs,n_discontinuous_pairs,n_voxels_with_discontinuity, &
        !$omp n_identical_pairs,n_one_step_pairs,n_large_jump_pairs,stepdiff_counts)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                    if( allocated(srcmap) ) then
                        if( srcmap(i,j,k) /= 1 ) cycle
                    end if
                    lp_i = int(filtmap(i,j,k))
                    coord_i = nu_candidate_coord_for_label(lp_i)
                    n_discontinuous_neighbors = 0
                    do dk = -1, 1
                        do dj = -1, 1
                            do di = -1, 1
                                if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
                                ni = i + di
                                nj = j + dj
                                nk = k + dk
                                if( ni < 1 .or. ni > nx ) cycle
                                if( nj < 1 .or. nj > ny ) cycle
                                if( nk < 1 .or. nk > nz ) cycle
                                if( .not.active_nu_mask_at(mask, ni, nj, nk) ) cycle
                                if( allocated(srcmap) ) then
                                    if( srcmap(ni,nj,nk) /= 1 ) cycle
                                end if
                                lp_j = int(filtmap(ni,nj,nk))
                                coord_j = nu_candidate_coord_for_label(lp_j)
                                step_diff = nint(abs(coord_i - coord_j))
                                if( step_diff > thresh ) n_discontinuous_neighbors = n_discontinuous_neighbors + 1
                                l_count_pair = dk > 0 .or. (dk == 0 .and. dj > 0) .or. &
                                    &(dk == 0 .and. dj == 0 .and. di > 0)
                                if( .not.l_count_pair ) cycle
                                n_total_neighbor_pairs = n_total_neighbor_pairs + 1
                                if( step_diff == 0 )then
                                    n_identical_pairs = n_identical_pairs + 1
                                else
                                    n_discontinuous_pairs = n_discontinuous_pairs + 1
                                    if( step_diff <= large_thresh )then
                                        n_one_step_pairs = n_one_step_pairs + 1
                                    else
                                        n_large_jump_pairs = n_large_jump_pairs + 1
                                    endif
                                    if( step_diff >= 1 .and. step_diff <= max_step_diff ) &
                                        &stepdiff_counts(step_diff) = stepdiff_counts(step_diff) + 1
                                endif
                            end do
                        end do
                    end do
                    if( n_discontinuous_neighbors > 0 ) n_voxels_with_discontinuity = &
                        &n_voxels_with_discontinuity + 1
                end do
            end do
        end do
        !$omp end parallel do
        n_analyzed = count_active_nu_mask(mask, 1)
        pct_vox = 0.
        if( n_analyzed > 0 ) pct_vox = 100. * real(n_voxels_with_discontinuity) / real(n_analyzed)
        pct_pairs = 0.
        pct_large = 0.
        if( n_total_neighbor_pairs > 0 )then
            pct_pairs = 100. * real(n_discontinuous_pairs) / real(n_total_neighbor_pairs)
            pct_large = 100. * real(n_large_jump_pairs) / real(n_total_neighbor_pairs)
        endif
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NU NEIGHBOR CONTINUITY'
        write(logfhandle,'(A,I0,A,I0,A)') '    Retained-bank boundary: > ', thresh, &
            &'; large jump: > ', large_thresh, ' coordinate step(s)'
        write(logfhandle,'(A,I12)') '    Voxels analyzed: ', n_analyzed
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Voxels touching a retained-bank boundary: ', &
            &n_voxels_with_discontinuity, ' (', pct_vox, '%)'
        if( n_total_neighbor_pairs > 0 ) then
            write(logfhandle,'(A,I14)') '    Neighbor pairs examined (unique links): ', n_total_neighbor_pairs
            pct = 100. * real(n_identical_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(A,I14,A,F8.2,A)') '      identical:          ', n_identical_pairs, ' (', pct, '%)'
            pct = 100. * real(n_one_step_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(A,I14,A,F8.2,A)') '      one-step boundary:  ', n_one_step_pairs, ' (', pct, '%)'
            write(logfhandle,'(A,I14,A,F8.2,A)') '      larger jump:        ', n_large_jump_pairs, ' (', pct_large, '%)'
            write(logfhandle,'(A,I14,A,F8.2,A)') '      boundary total:     ', n_discontinuous_pairs, ' (', pct_pairs, '%)'
        endif
        write(logfhandle,'(A)') '    Retained-bank coordinate-step difference distribution:'
        write(logfhandle,'(A)') '      Step           Pairs       Pct    Class'
        if( n_total_neighbor_pairs > 0 )then
            pct = 100. * real(n_identical_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') 0, n_identical_pairs, pct, 'identical'
        endif
        if( n_total_neighbor_pairs > 0 )then
            do ii = 1, max_step_diff
                if( stepdiff_counts(ii) > 0 ) then
                    pct = 100. * real(stepdiff_counts(ii)) / real(n_total_neighbor_pairs)
                    if( ii <= large_thresh )then
                        write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') &
                            &ii, stepdiff_counts(ii), pct, 'one-step boundary'
                    else
                        write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') &
                            &ii, stepdiff_counts(ii), pct, 'larger jump'
                    endif
                end if
            end do
        endif
        if( n_total_neighbor_pairs > 0 ) then
            if( pct_pairs <= 5. .and. pct_large <= 1. ) then
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: low boundary rate; local resolution map is spatially smooth'
            else if( pct_pairs <= 10. .and. pct_large <= 5. ) then
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: moderate boundary rate; inspect one-step and larger-jump structure'
            else
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: high boundary rate; inspect mask support, objective maps, and label-prior convergence'
            end if
        else
            write(logfhandle,'(A)') '    Continuity assessment: no neighbor pairs found in mask; analysis is inconclusive'
        end if
        call log_nu_source_neighbor_continuity(mask)
        write(logfhandle,'(A)') ''
        deallocate(stepdiff_counts)
    end subroutine analyze_filtmap_neighbor_continuity

    subroutine log_nu_source_neighbor_continuity( mask )
        logical, optional, intent(in) :: mask(:,:,:)
        integer :: i, j, k, di, dj, dk, ni, nj, nk, nx, ny, nz
        integer :: n_analyzed, n_base_selected, n_aux_selected
        integer :: n_source_pairs, n_source_boundary_pairs, n_source_identical_pairs
        integer :: n_voxels_with_source_boundary, n_boundary_neighbors
        real    :: pct_aux, pct_pairs, pct_vox
        logical :: l_count_pair
        if( .not.allocated(srcmap) ) return
        call require_valid_stats_mask(mask, 'log_nu_source_neighbor_continuity')
        n_analyzed      = count_active_nu_mask(mask, 0)
        n_base_selected = count_active_nu_mask(mask, 1)
        n_aux_selected  = n_analyzed - n_base_selected
        if( n_aux_selected == 0 ) return
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        n_source_pairs                   = 0
        n_source_boundary_pairs          = 0
        n_source_identical_pairs         = 0
        n_voxels_with_source_boundary    = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,di,dj,dk,ni,nj,nk,n_boundary_neighbors,l_count_pair) &
        !$omp reduction(+:n_source_pairs,n_source_boundary_pairs,n_source_identical_pairs, &
        !$omp n_voxels_with_source_boundary)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                    n_boundary_neighbors = 0
                    do dk = -1, 1
                        do dj = -1, 1
                            do di = -1, 1
                                if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle
                                ni = i + di
                                nj = j + dj
                                nk = k + dk
                                if( ni < 1 .or. ni > nx ) cycle
                                if( nj < 1 .or. nj > ny ) cycle
                                if( nk < 1 .or. nk > nz ) cycle
                                if( .not.active_nu_mask_at(mask, ni, nj, nk) ) cycle
                                if( srcmap(i,j,k) == srcmap(ni,nj,nk) )then
                                    l_count_pair = dk > 0 .or. (dk == 0 .and. dj > 0) .or. &
                                        &(dk == 0 .and. dj == 0 .and. di > 0)
                                    if( l_count_pair )then
                                        n_source_pairs = n_source_pairs + 1
                                        n_source_identical_pairs = n_source_identical_pairs + 1
                                    endif
                                else
                                    n_boundary_neighbors = n_boundary_neighbors + 1
                                    l_count_pair = dk > 0 .or. (dk == 0 .and. dj > 0) .or. &
                                        &(dk == 0 .and. dj == 0 .and. di > 0)
                                    if( l_count_pair )then
                                        n_source_pairs = n_source_pairs + 1
                                        n_source_boundary_pairs = n_source_boundary_pairs + 1
                                    endif
                                endif
                            end do
                        end do
                    end do
                    if( n_boundary_neighbors > 0 ) n_voxels_with_source_boundary = &
                        &n_voxels_with_source_boundary + 1
                end do
            end do
        end do
        !$omp end parallel do
        pct_aux = 0.
        if( n_analyzed > 0 ) pct_aux = 100. * real(n_aux_selected) / real(n_analyzed)
        pct_pairs = 0.
        if( n_source_pairs > 0 ) pct_pairs = 100. * real(n_source_boundary_pairs) / real(n_source_pairs)
        pct_vox = 0.
        if( n_analyzed > 0 ) pct_vox = 100. * real(n_voxels_with_source_boundary) / real(n_analyzed)
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NU SOURCE CONTINUITY'
        write(logfhandle,'(A,I12)') '    Mask voxels:      ', n_analyzed
        write(logfhandle,'(A,I12)') '    Base-bank voxels: ', n_base_selected
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Auxiliary voxels: ', n_aux_selected, ' (', pct_aux, '%)'
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Voxels touching a source boundary: ', &
            &n_voxels_with_source_boundary, ' (', pct_vox, '%)'
        if( n_source_pairs > 0 )then
            write(logfhandle,'(A,I14)') '    Neighbor pairs examined (unique links): ', n_source_pairs
            write(logfhandle,'(A,I14)') '      same source:      ', n_source_identical_pairs
            write(logfhandle,'(A,I14,A,F8.2,A)') '      source boundary:  ', &
                &n_source_boundary_pairs, ' (', pct_pairs, '%)'
        endif
        if( n_source_pairs > 0 )then
            if( pct_pairs <= 5. )then
                write(logfhandle,'(A)') &
                    &'    Source assessment: low aux/base boundary rate'
            else if( pct_pairs <= 10. )then
                write(logfhandle,'(A)') &
                    &'    Source assessment: moderate aux/base boundary rate; inspect auxiliary assignments'
            else
                write(logfhandle,'(A)') &
                    &'    Source assessment: high aux/base boundary rate; auxiliary source may be patchy'
            endif
        endif
    end subroutine log_nu_source_neighbor_continuity

    subroutine require_valid_stats_mask( mask, caller )
        logical, optional, intent(in) :: mask(:,:,:)
        character(len=*), intent(in) :: caller
        if( present(mask) )then
            if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; '//caller)
            if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in '//caller)
        else
            if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; '//caller)
            if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; '//caller)
            if( any(shape(nu_lmask) /= ldim) ) THROW_HARD('internal mask shape mismatch in '//caller)
        endif
    end subroutine require_valid_stats_mask

    logical function active_nu_mask_at( mask, i, j, k )
        logical, optional, intent(in) :: mask(:,:,:)
        integer, intent(in) :: i, j, k
        if( present(mask) )then
            active_nu_mask_at = mask(i,j,k)
            if( allocated(nu_lmask) ) active_nu_mask_at = active_nu_mask_at .and. nu_lmask(i,j,k)
        else
            active_nu_mask_at = nu_lmask(i,j,k)
        endif
    end function active_nu_mask_at

    integer function count_active_nu_mask( mask, source_id ) result(nactive)
        logical, optional, intent(in) :: mask(:,:,:)
        integer, intent(in) :: source_id
        integer :: i, j, k, imask
        nactive = 0
        if( present(mask) )then
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:nactive) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                        if( source_id > 0 )then
                            if( allocated(srcmap) )then
                                if( int(srcmap(i,j,k)) /= source_id ) cycle
                            else if( source_id /= 1 )then
                                cycle
                            endif
                        endif
                        nactive = nactive + 1
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) reduction(+:nactive) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( source_id > 0 )then
                    if( allocated(srcmap) )then
                        if( int(srcmap(i,j,k)) /= source_id ) cycle
                    else if( source_id /= 1 )then
                        cycle
                    endif
                endif
                nactive = nactive + 1
            end do
            !$omp end parallel do
        endif
    end function count_active_nu_mask

    integer function count_active_nu_label( mask, label_id, source_id ) result(nactive)
        logical, optional, intent(in) :: mask(:,:,:)
        integer, intent(in) :: label_id, source_id
        integer :: i, j, k, imask
        nactive = 0
        if( present(mask) )then
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:nactive) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.active_nu_mask_at(mask, i, j, k) ) cycle
                        if( int(filtmap(i,j,k)) /= label_id ) cycle
                        if( source_id > 0 )then
                            if( allocated(srcmap) )then
                                if( int(srcmap(i,j,k)) /= source_id ) cycle
                            else if( source_id /= 1 )then
                                cycle
                            endif
                        endif
                        nactive = nactive + 1
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) reduction(+:nactive) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( int(filtmap(i,j,k)) /= label_id ) cycle
                if( source_id > 0 )then
                    if( allocated(srcmap) )then
                        if( int(srcmap(i,j,k)) /= source_id ) cycle
                    else if( source_id /= 1 )then
                        cycle
                    endif
                endif
                nactive = nactive + 1
            end do
            !$omp end parallel do
        endif
    end function count_active_nu_label

    real function lowpass_histogram_value_at( counts, pos )
        integer, intent(in) :: counts(:), pos
        integer :: icut, nseen
        if( pos < 1 ) THROW_HARD('invalid order statistic position; lowpass_histogram_value_at')
        nseen = 0
        do icut = size(counts), 1, -1
            nseen = nseen + counts(icut)
            if( nseen >= pos )then
                lowpass_histogram_value_at = cutoff_find_to_lowpass_limit(icut)
                return
            endif
        end do
        THROW_HARD('order statistic exceeds histogram count; lowpass_histogram_value_at')
    end function lowpass_histogram_value_at

end submodule simple_nu_filter_stats
