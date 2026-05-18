!@descr: simple nu filter stats implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_stats
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
        real, allocatable, intent(inout) :: lowpass_vals(:)
        logical, intent(in)              :: mask(:,:,:)
        integer :: i, j, k, nvals, ival
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before pack_filtmap_lowpass_limits')
        if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in pack_filtmap_lowpass_limits')
        if( allocated(srcmap) ) then
            nvals = count(mask .and. srcmap == 1)
        else
            nvals = count(mask)
        end if
        if( allocated(lowpass_vals) ) deallocate(lowpass_vals)
        allocate(lowpass_vals(nvals), source=0.)
        ival = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not. mask(i,j,k) ) cycle
                    if( allocated(srcmap) ) then
                        if( srcmap(i,j,k) /= 1 ) cycle
                    end if
                    ival = ival + 1
                    lowpass_vals(ival) = cutoff_find_to_lowpass_limit(filtmap(i,j,k))
                end do
            end do
        end do
    end subroutine pack_filtmap_lowpass_limits

    module subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, intent(in)             :: mask(:,:,:)
        real, allocatable :: lowpass_vals(:)
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_stats')
        call pack_filtmap_lowpass_limits(lowpass_vals, mask)
        if( size(lowpass_vals) == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_stats')
        call calc_stats(lowpass_vals, statvars)
        deallocate(lowpass_vals)
    end subroutine calc_filtmap_lowpass_stats

    module subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
        integer, intent(out) :: counts(:)
        real,    intent(out) :: percentages(:)
        logical, intent(in)  :: mask(:,:,:)
        integer :: icut, nselected
        if( .not.allocated(filtmap)                 ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_histogram')
        if( .not.allocated(cutoff_finds)            ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before calc_filtmap_lowpass_histogram')
        if( size(counts) /= size(cutoff_finds)      ) THROW_HARD('counts size mismatch in calc_filtmap_lowpass_histogram')
        if( size(percentages) /= size(cutoff_finds) ) THROW_HARD('percentages size mismatch in calc_filtmap_lowpass_histogram')
        if( any(shape(mask) /= shape(filtmap))      ) THROW_HARD('mask shape mismatch in calc_filtmap_lowpass_histogram')
        if( allocated(srcmap) ) then
            nselected = count(mask .and. srcmap == 1)
        else
            nselected = count(mask)
        end if
        counts       = 0
        percentages  = 0.
        if( nselected == 0 ) return
        do icut = 1, size(cutoff_finds)
            if( allocated(srcmap) ) then
                counts(icut) = count(filtmap == icut .and. srcmap == 1 .and. mask)
            else
                counts(icut) = count(filtmap == icut .and. mask)
            end if
            percentages(icut) = 100. * real(counts(icut)) / real(nselected)
        end do
    end subroutine calc_filtmap_lowpass_histogram

    module subroutine print_filtmap_lowpass_histogram( mask, aux_resolutions )
        logical,        intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
        integer, allocatable :: counts(:)
        real,    allocatable :: percentages(:)
        integer :: icut, iaux, nselected, nvox
        real    :: pct
        character(len=8) :: auxtag
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before print_filtmap_lowpass_histogram')
        allocate(counts(size(cutoff_finds)), percentages(size(cutoff_finds)))
        call calc_filtmap_lowpass_histogram(counts, percentages, mask)
        if( allocated(srcmap) )then
            nselected = count(mask .and. srcmap == 1)
        else
            nselected = count(mask)
        endif
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
            nselected = count(mask)
            write(logfhandle,'(A)') ''
            write(logfhandle,'(A)')     '>>> NU AUXILIARY SOURCE ASSIGNMENTS'
            write(logfhandle,'(A,I12)') '    Mask voxels:      ', nselected
            write(logfhandle,'(A)')     '    Source    Resolution (A)        Voxels    Pct mask'
            do iaux = 1, size(aux_even_bank)
                nvox = count(srcmap == iaux + 1 .and. mask)
                pct = 0.
                if( nselected > 0 ) pct = 100. * real(nvox) / real(nselected)
                write(auxtag,'(A,I0,A)') 'Aux', iaux, '@'
                write(logfhandle,'(4X,A8,2X,F14.3,2X,I12,2X,F8.2,A)') auxtag, aux_resolutions(iaux), nvox, pct, '%'
            end do
        end if
        deallocate(counts, percentages)
    end subroutine print_filtmap_lowpass_histogram

    module subroutine print_nu_filtmap_lowpass_stats( mask, aux_resolutions )
        logical,        intent(in) :: mask(:,:,:)
        real, optional, intent(in) :: aux_resolutions(:)
        type(stats_struct) :: statvars
        integer :: nbase
        if( allocated(srcmap) ) then
            nbase = count(srcmap == 1 .and. mask)
            if( nbase == 0 ) then
                write(logfhandle,'(A)') '>>> No base low-pass selections remain after auxiliary-source optimization'
                call print_filtmap_lowpass_histogram(mask, aux_resolutions)
                return
            end if
        end if
        call calc_filtmap_lowpass_stats(statvars, mask)
        if( allocated(srcmap) )then
            nbase = count(srcmap == 1 .and. mask)
        else
            nbase = count(mask)
        endif
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NU FILTER LOCAL RESOLUTION SUMMARY'
        write(logfhandle,'(A,I12)') '    Voxels analyzed: ', nbase
        write(logfhandle,'(A)')     '              Mean    Median     Sigma       Min       Max'
        write(logfhandle,'(A,5F10.3)') '    Angstrom ', statvars%avg, statvars%med, statvars%sdev, statvars%minv, statvars%maxv
        write(logfhandle,'(A)') ''
        call print_filtmap_lowpass_histogram(mask, aux_resolutions)
    end subroutine print_nu_filtmap_lowpass_stats

    module subroutine analyze_filtmap_neighbor_continuity( mask )
        logical, intent(in) :: mask(:,:,:)
        integer :: i, j, k, di, dj, dk, ni, nj, nk
        integer :: lp_i, lp_j, lp_diff, max_diff, n_neighbors, n_discontinuous_neighbors
        integer :: n_total_neighbor_pairs, n_discontinuous_pairs, n_voxels_with_discontinuity, nx, ny, nz, ii, thresh, n_analyzed
        integer :: n_identical_pairs, n_tolerated_pairs, max_step_diff
        integer, allocatable :: stepdiff_counts(:)
        real :: pct, pct_vox, pct_pairs
        if( .not.allocated(filtmap)            ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before analyze_filtmap_neighbor_continuity')
        if( .not.allocated(cutoff_finds)       ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before analyze_filtmap_neighbor_continuity')
        if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in analyze_filtmap_neighbor_continuity')
        ! One low-pass step is tolerated by the ordered-label prior, so only
        ! neighbor pairs beyond this threshold are counted as discontinuities.
        thresh = DISCONT_STEP_THRESH
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        max_step_diff = max(1, size(cutoff_finds) - 1)
        allocate(stepdiff_counts(max_step_diff), source=0)
        n_voxels_with_discontinuity = 0
        n_total_neighbor_pairs      = 0
        n_discontinuous_pairs       = 0
        n_identical_pairs           = 0
        n_tolerated_pairs           = 0
        ! Iterate through all voxels
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,di,dj,dk,ni,nj,nk,lp_i,lp_j,lp_diff,n_neighbors,n_discontinuous_neighbors,max_diff) &
        !$omp reduction(+:n_total_neighbor_pairs, n_discontinuous_pairs, n_voxels_with_discontinuity, &
        !$omp n_identical_pairs, n_tolerated_pairs, stepdiff_counts)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not.mask(i,j,k) ) cycle
                    if( allocated(srcmap) ) then
                        if( srcmap(i,j,k) /= 1 ) cycle  ! only check base bank voxels
                    end if
                    lp_i = filtmap(i,j,k)
                    n_neighbors = 0
                    n_discontinuous_neighbors = 0
                    max_diff = 0
                    ! Check 26-connected neighborhood (3x3x3 cube centered at (i,j,k))
                    do dk = -1, 1
                        do dj = -1, 1
                            do di = -1, 1
                                if( di == 0 .and. dj == 0 .and. dk == 0 ) cycle  ! skip self
                                ni = i + di
                                nj = j + dj
                                nk = k + dk
                                ! Check bounds
                                if( ni < 1 .or. ni > nx ) cycle
                                if( nj < 1 .or. nj > ny ) cycle
                                if( nk < 1 .or. nk > nz ) cycle
                                ! Check mask
                                if( .not.mask(ni,nj,nk) ) cycle
                                if( allocated(srcmap) ) then
                                    if( srcmap(ni,nj,nk) /= 1 ) cycle
                                end if
                                lp_j = filtmap(ni,nj,nk)
                                lp_diff = abs(lp_i - lp_j)
                                n_neighbors = n_neighbors + 1
                                n_total_neighbor_pairs = n_total_neighbor_pairs + 1
                                if( lp_diff == 0 )then
                                    n_identical_pairs = n_identical_pairs + 1
                                else if( lp_diff <= thresh )then
                                    n_tolerated_pairs = n_tolerated_pairs + 1
                                endif
                                if( lp_diff > thresh ) then
                                    n_discontinuous_neighbors = n_discontinuous_neighbors + 1
                                    n_discontinuous_pairs     = n_discontinuous_pairs + 1
                                end if
                                max_diff = max(max_diff, lp_diff)
                                if( lp_diff >= 1 .and. lp_diff <= max_step_diff ) then
                                    stepdiff_counts(lp_diff) = stepdiff_counts(lp_diff) + 1
                                end if
                            end do
                        end do
                    end do
                    ! Check if this voxel has any discontinuous neighbors
                    if( n_discontinuous_neighbors > 0 ) then
                        n_voxels_with_discontinuity = n_voxels_with_discontinuity + 1
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        ! Count total masked voxels
        if( allocated(srcmap) ) then
            n_analyzed = count(mask .and. srcmap == 1)
        else
            n_analyzed = count(mask)
        end if
        pct_vox = 0.
        if( n_analyzed > 0 ) pct_vox = 100. * real(n_voxels_with_discontinuity) / real(n_analyzed)
        pct_pairs = 0.
        if( n_total_neighbor_pairs > 0 ) pct_pairs = 100. * real(n_discontinuous_pairs) / real(n_total_neighbor_pairs)
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NU NEIGHBOR CONTINUITY'
        write(logfhandle,'(A,I0,A,I0,A)') '    LP-step tolerance: <= ', thresh, '; discontinuity: > ', thresh, ' step(s)'
        write(logfhandle,'(A,I12)') '    Voxels analyzed: ', n_analyzed
        write(logfhandle,'(A,I12,A,F8.2,A)') '    Voxels with discontinuous neighbors: ', &
            &n_voxels_with_discontinuity, ' (', pct_vox, '%)'
        if( n_total_neighbor_pairs > 0 ) then
            write(logfhandle,'(A,I14)') '    Neighbor pairs examined: ', n_total_neighbor_pairs
            pct = 100. * real(n_identical_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(A,I14,A,F8.2,A)') '      identical:      ', n_identical_pairs, ' (', pct, '%)'
            pct = 100. * real(n_tolerated_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(A,I14,A,F8.2,A)') '      tolerated:      ', n_tolerated_pairs, ' (', pct, '%)'
            write(logfhandle,'(A,I14,A,F8.2,A)') '      discontinuous:  ', n_discontinuous_pairs, ' (', pct_pairs, '%)'
        endif
        ! Print distribution of all pair step-differences without implying that
        ! tolerated one-step differences are discontinuities.
        write(logfhandle,'(A)') '    LP-step difference distribution:'
        write(logfhandle,'(A)') '      Step           Pairs       Pct    Class'
        if( n_total_neighbor_pairs > 0 )then
            pct = 100. * real(n_identical_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') 0, n_identical_pairs, pct, 'identical'
        endif
        do ii = 1, max_step_diff
            if( stepdiff_counts(ii) > 0 ) then
                pct = 100. * real(stepdiff_counts(ii)) / real(n_total_neighbor_pairs)
                if( ii <= thresh )then
                    write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') ii, stepdiff_counts(ii), pct, 'tolerated'
                else
                    write(logfhandle,'(6X,I4,2X,I14,2X,F8.3,4X,A)') ii, stepdiff_counts(ii), pct, 'discontinuous'
                endif
            end if
        end do
        if( n_total_neighbor_pairs > 0 ) then
            if( pct_pairs <= 5. ) then
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: low discontinuity rate; local resolution map is spatially smooth'
            else if( pct_pairs <= 10. ) then
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: moderate discontinuity rate; label smoothing may benefit from additional convergence'
            else
                write(logfhandle,'(A)') &
                    &'    Continuity assessment: high discontinuity rate; inspect mask support, objective maps, and label-prior convergence'
            end if
        else
            write(logfhandle,'(A)') '    Continuity assessment: no neighbor pairs found in mask; analysis is inconclusive'
        end if
        write(logfhandle,'(A)') ''
        deallocate(stepdiff_counts)
    end subroutine analyze_filtmap_neighbor_continuity

end submodule simple_nu_filter_stats
