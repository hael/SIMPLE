!@descr: nonuniform filtering of even/odd volumes
!
! A typical call sequence would be:
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::])
!    call optimize_nu_cutoff_finds()
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! Updated call sequence with iterative high-resolution refinement
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::])
!    call optimize_nu_cutoff_finds()
!    call extend_nu_filter_highres_shell_next(vol_even, vol_odd) ! optional shell step
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! Postprocess workflows can instead build the full available Fourier-shell bank
! up front and then apply the resulting local filter map to a merged/sharpened map:
!    call setup_nu_dmats(vol_even, vol_odd, l_mask, [real ::], l_full_shell_bank=.true.)
!    call optimize_nu_cutoff_finds()
!    call nu_filter_vol(vol, vol_filt)
!
! supports auxiliary candidate pairs that can compete with the
! base low-pass bank during voxelwise optimization
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
use simple_neighs,      only: neigh_8_3D
implicit none
#include "simple_local_flags.inc"

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, nu_filter_vol, cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, extend_nu_filter_highres_shell_next, extend_nu_filter_highres_shells,&
          analyze_filtmap_neighbor_continuity, nu_highres_extension_stats
private

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
! Minimum finest-frontier fraction of the NU mask required before testing a
! finer shell. Zero means challenge whenever at least one frontier voxel exists.
real,             parameter   :: NU_HIGHRES_EXTENSION_THRESHOLD_PCT = 0.
real,             parameter   :: NU_HIGHRES_PROMOTION_TESTED_PCT    = 5.  ! promotion threshold among tested frontier voxels
! Physical half-width of the tent regularization kernel. The smoother consumes
! this as an integer pixel radius, so the full tent base spans 2*radius + 1
! voxels along each axis; 8 A at 1 A/px gives radius=8 and a 17-voxel base.
real,             parameter   :: WINSZ_TENT_ANGSTROM       = 8.
integer,          parameter   :: DISCONT_STEP_THRESH       = 1
integer,          parameter   :: NU_LABEL_SMOOTH_MAXITS    = 3
integer,          parameter   :: NU_LABEL_SMOOTH_STEP_TOL  = 1
integer,          parameter   :: NU_LABEL_SMOOTH_NNEIGH    = 26
integer,          parameter   :: NU_LABEL_SMOOTH_NCOLORS   = 8
integer,          parameter   :: NU_FULL_SHELL_ZERO_GUARD   = 3
real,             parameter   :: NU_LABEL_SMOOTH_BETA_FRAC = 2.0
real,             parameter   :: NU_LABEL_SMOOTH_QUAD_FRAC = 1.0
real,             parameter   :: NU_LABEL_SMOOTH_TIE_EPS   = 1.e-6
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN      = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD       = 'nu_filter_cache_odd'
real,             allocatable :: dmats(:,:,:,:)
real,             allocatable :: bwfilters(:,:)
real,             allocatable :: candidate_coords(:)
integer,          allocatable :: filtmap(:,:,:)
integer,          allocatable :: srcmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
real,             allocatable :: dmat_finest_cached(:,:,:)
logical,          allocatable :: nu_lmask(:,:,:)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
integer :: winsz_tent
real    :: smpd
logical :: l_full_shell_bank_active = .false.

type :: nu_highres_extension_stats
    logical :: attempted    = .false.
    logical :: applied      = .false.
    logical :: promote_next = .false.
    real    :: new_limit = 0.
    integer :: n_mask     = 0
    integer :: n_tested   = 0
    integer :: n_extended = 0
    real    :: pct_tested_mask     = 0.
    real    :: pct_extended_tested = 0.
end type nu_highres_extension_stats

contains

    real function cutoff_find_to_lowpass_limit( icut )
        integer, intent(in) :: icut
        if( .not.allocated(cutoff_finds)            ) THROW_HARD('cutoff_finds not allocated; cutoff_find_to_lowpass_limit')
        if( icut < 1 .or. icut > size(cutoff_finds) ) THROW_HARD('cutoff index out of range; cutoff_find_to_lowpass_limit')
        cutoff_find_to_lowpass_limit = calc_lowpass_lim(cutoff_finds(icut), box, smpd)
    end function cutoff_find_to_lowpass_limit

    subroutine init_nu_filter( vol_even, vol_odd, n_highres_steps, l_full_shell_bank )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(in) :: n_highres_steps
        logical, optional, intent(in) :: l_full_shell_bank
        integer, allocatable :: cutoff_finds_tmp(:)
        integer :: i, n_extra, n_valid, new_find, n_added
        logical :: l_full_bank
        ldim = vol_even%get_ldim()
        smpd = vol_even%get_smpd()
        box  = ldim(1)
        if( any(vol_odd%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; init_nu_filter')
        if( abs(vol_odd%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; init_nu_filter')
        if( smpd <= TINY ) THROW_HARD('Input volume smpd must be positive; init_nu_filter')
        ! Convert the physical half-width to pixel radius for tent_smooth_3d.
        ! The resulting support is 2*winsz_tent + 1 voxels.
        winsz_tent = max(1, nint(WINSZ_TENT_ANGSTROM / smpd))
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        l_full_bank = .false.
        if( present(l_full_shell_bank) ) l_full_bank = l_full_shell_bank
        l_full_shell_bank_active = l_full_bank
        n_extra = 0
        if( l_full_bank )then
            n_extra = max(0, box / 2)
        else if( present(n_highres_steps) )then
            n_extra = min(max(0, n_highres_steps), max(0, box / 2))
        endif
        allocate(cutoff_finds_tmp(size(lowpass_limits) + n_extra))
        do i = 1, size(lowpass_limits)
            cutoff_finds_tmp(i) = calc_fourier_index(lowpass_limits(i), box, smpd)
        end do
        n_valid = size(lowpass_limits)
        if( l_full_bank )then
            do new_find = cutoff_finds_tmp(n_valid) + 1, box / 2
                if( any(cutoff_finds_tmp(:n_valid) == new_find) ) cycle
                n_valid = n_valid + 1
                cutoff_finds_tmp(n_valid) = new_find
            end do
        else
            n_added = 0
            new_find = cutoff_finds_tmp(n_valid) + 1
            do while( n_added < n_extra .and. new_find <= box / 2 )
                if( .not.any(cutoff_finds_tmp(:n_valid) == new_find) )then
                    n_valid = n_valid + 1
                    cutoff_finds_tmp(n_valid) = new_find
                    n_added = n_added + 1
                endif
                new_find = new_find + 1
            end do
        endif
        allocate(cutoff_finds(n_valid))
        cutoff_finds = cutoff_finds_tmp(:n_valid)
        deallocate(cutoff_finds_tmp)
        if( allocated(bwfilters) ) deallocate(bwfilters)
        allocate(bwfilters(box,size(cutoff_finds)), source=0.)
        do i = 1, size(cutoff_finds)
            call butterworth_filter(cutoff_finds(i), bwfilters(:,i))
        end do
    end subroutine init_nu_filter

    function filtered_vol_fname( cache_prefix, cutoff_find ) result( fname )
        class(string), intent(in) :: cache_prefix
        integer,       intent(in) :: cutoff_find
        type(string) :: fname
        fname = cache_prefix%to_char()//'_k_'//int2str(cutoff_find)//'.mrc'
    end function filtered_vol_fname

    logical function filtered_vols_cached( cache_prefix )
        class(string), intent(in) :: cache_prefix
        integer :: i
        filtered_vols_cached = .false.
        if( .not.allocated(cutoff_finds) ) return
        do i = 1, size(cutoff_finds)
            if( .not.file_exists(filtered_vol_fname(cache_prefix, cutoff_finds(i))) ) return
        end do
        filtered_vols_cached = .true.
    end function filtered_vols_cached

    subroutine delete_cached_filtered_vols( cache_prefix )
        class(string), intent(in) :: cache_prefix
        integer :: i
        if( allocated(cutoff_finds) ) then
            do i = 1, size(cutoff_finds)
                if( file_exists(filtered_vol_fname(cache_prefix, cutoff_finds(i))) ) then
                    call del_file(filtered_vol_fname(cache_prefix, cutoff_finds(i)))
                end if
            end do
        end if
    end subroutine delete_cached_filtered_vols

    subroutine cleanup_nu_filter
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_EVEN))
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_ODD))
        if( allocated(dmats)              ) deallocate(dmats)
        if( allocated(bwfilters)          ) deallocate(bwfilters)
        if( allocated(candidate_coords)   ) deallocate(candidate_coords)
        if( allocated(filtmap)            ) deallocate(filtmap)
        if( allocated(srcmap)             ) deallocate(srcmap)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( allocated(nu_lmask)           ) deallocate(nu_lmask)
        call cleanup_aux_bank
        ldim = 0
        box  = 0
        winsz_tent = 0
        smpd = 0.
        l_full_shell_bank_active = .false.
    end subroutine cleanup_nu_filter

    subroutine cleanup_aux_bank
        integer :: i
        if( allocated(aux_even_bank) ) then
            do i = 1, size(aux_even_bank)
                call aux_even_bank(i)%kill
            end do
            deallocate(aux_even_bank)
        end if
        if( allocated(aux_odd_bank) ) then
            do i = 1, size(aux_odd_bank)
                call aux_odd_bank(i)%kill
            end do
            deallocate(aux_odd_bank)
        end if
    end subroutine cleanup_aux_bank

    subroutine validate_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
        integer :: i
        if( size(aux_even) /= size(aux_odd) ) THROW_HARD('Auxiliary even/odd banks differ in size; validate_aux_volumes')
        do i = 1, size(aux_even)
            if( any(aux_even(i)%get_ldim() /= ldim)       ) THROW_HARD('Auxiliary even volume dimensions differ; validate_aux_volumes')
            if( any(aux_odd(i)%get_ldim()  /= ldim)       ) THROW_HARD('Auxiliary odd volume dimensions differ; validate_aux_volumes')
            if( abs(aux_even(i)%get_smpd() - smpd) > TINY ) THROW_HARD('Auxiliary even volume smpd differs; validate_aux_volumes')
            if( abs(aux_odd(i)%get_smpd()  - smpd) > TINY ) THROW_HARD('Auxiliary odd volume smpd differs; validate_aux_volumes')
        end do
    end subroutine validate_aux_volumes

    subroutine stash_aux_volumes( aux_even, aux_odd )
        type(image), intent(in) :: aux_even(:), aux_odd(:)
        integer :: i
        call cleanup_aux_bank
        call validate_aux_volumes(aux_even, aux_odd)
        allocate(aux_even_bank(size(aux_even)))
        allocate(aux_odd_bank(size(aux_odd)))
        do i = 1, size(aux_even)
            call aux_even_bank(i)%copy(aux_even(i))
            call aux_odd_bank(i)%copy(aux_odd(i))
        end do
    end subroutine stash_aux_volumes

    subroutine cache_filtered_vols( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
        type(string) :: even_cache_fname, odd_cache_fname
        integer :: i, winsz
        real    :: edge_mean
        logical :: even_cached, odd_cached
        winsz = nint(COSMSKHALFWIDTH)
        even_cached = filtered_vols_cached(string(NU_FILTER_CACHE_EVEN))
        odd_cached  = filtered_vols_cached(string(NU_FILTER_CACHE_ODD))
        if( even_cached .and. odd_cached ) return
        call vol_even_copy_cmat%copy(vol_even)
        call vol_odd_copy_cmat%copy(vol_odd)
        call vol_even_copy_cmat%set_wthreads(.true.)
        call vol_odd_copy_cmat%set_wthreads(.true.)
        call vol_even_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_odd_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_even_copy_cmat%fft
        call vol_odd_copy_cmat%fft
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_ft(.true.)
        call vol_odd_filt%set_ft(.true.)
        call vol_even_filt%set_wthreads(.true.)
        call vol_odd_filt%set_wthreads(.true.)
        do i = 1, size(cutoff_finds)
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
            if( file_exists(even_cache_fname) .and. file_exists(odd_cache_fname) ) cycle
            call vol_even_filt%copy_fast(vol_even_copy_cmat)
            call vol_odd_filt%copy_fast(vol_odd_copy_cmat)
            call vol_even_filt%apply_filter(vol_odd_filt, bwfilters(:,i)) ! double filter application
            call vol_even_filt%ifft
            call vol_odd_filt%ifft
            call vol_even_filt%write(even_cache_fname, del_if_exists=.true.)
            call vol_odd_filt%write(odd_cache_fname,  del_if_exists=.true.)
        end do
        call vol_even_copy_cmat%kill
        call vol_odd_copy_cmat%kill
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine cache_filtered_vols

    subroutine generate_single_filtered_pair( vol_even, vol_odd, cutoff_find, even_cache_fname, odd_cache_fname )
        class(image),  intent(in) :: vol_even, vol_odd
        integer,       intent(in) :: cutoff_find
        class(string), intent(in) :: even_cache_fname, odd_cache_fname
        type(image) :: vol_even_filt, vol_odd_filt
        type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
        integer :: winsz
        real    :: edge_mean
        real, allocatable :: bwfilter(:)
        if( .not.allocated(cutoff_finds) ) call init_nu_filter(vol_even, vol_odd)
        winsz = nint(COSMSKHALFWIDTH)
        call vol_even_copy_cmat%copy(vol_even)
        call vol_odd_copy_cmat%copy(vol_odd)
        call vol_even_copy_cmat%set_wthreads(.true.)
        call vol_odd_copy_cmat%set_wthreads(.true.)
        call vol_even_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_odd_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_even_copy_cmat%fft
        call vol_odd_copy_cmat%fft
        allocate(bwfilter(box), source=0.)
        call butterworth_filter(cutoff_find, bwfilter)
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_ft(.true.)
        call vol_odd_filt%set_ft(.true.)
        call vol_even_filt%set_wthreads(.true.)
        call vol_odd_filt%set_wthreads(.true.)
        call vol_even_filt%copy_fast(vol_even_copy_cmat)
        call vol_odd_filt%copy_fast(vol_odd_copy_cmat)
        call vol_even_filt%apply_filter(vol_odd_filt, bwfilter)
        call vol_even_filt%ifft
        call vol_odd_filt%ifft
        call vol_even_filt%write(even_cache_fname, del_if_exists=.true.)
        call vol_odd_filt%write(odd_cache_fname,  del_if_exists=.true.)
        call vol_even_copy_cmat%kill
        call vol_odd_copy_cmat%kill
        call vol_even_filt%kill
        call vol_odd_filt%kill
        deallocate(bwfilter)
    end subroutine generate_single_filtered_pair

    subroutine delete_cached_filtered_pair( cutoff_find )
        integer, intent(in) :: cutoff_find
        type(string) :: cache_fname
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_find)
        if( file_exists(cache_fname) ) call del_file(cache_fname)
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_find)
        if( file_exists(cache_fname) ) call del_file(cache_fname)
    end subroutine delete_cached_filtered_pair

    subroutine trim_nu_base_bank( n_keep )
        integer, intent(in) :: n_keep
        integer, allocatable :: cutoff_finds_new(:)
        real,    allocatable :: bwfilters_new(:,:), dmats_new(:,:,:,:)
        integer :: i, n_old
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; trim_nu_base_bank')
        n_old = size(cutoff_finds)
        if( n_keep < 1 .or. n_keep > n_old ) THROW_HARD('invalid n_keep; trim_nu_base_bank')
        if( n_keep == n_old ) return
        do i = n_keep + 1, n_old
            call delete_cached_filtered_pair(cutoff_finds(i))
        end do
        allocate(cutoff_finds_new(n_keep))
        cutoff_finds_new = cutoff_finds(:n_keep)
        call move_alloc(cutoff_finds_new, cutoff_finds)
        if( allocated(bwfilters) )then
            allocate(bwfilters_new(box,n_keep), source=0.)
            bwfilters_new = bwfilters(:,:n_keep)
            call move_alloc(bwfilters_new, bwfilters)
        endif
        if( allocated(dmats) )then
            if( size(dmats,4) < n_keep ) THROW_HARD('dmats too short; trim_nu_base_bank')
            allocate(dmats_new(ldim(1),ldim(2),ldim(3),n_keep), source=0.)
            dmats_new = dmats(:,:,:,:n_keep)
            call move_alloc(dmats_new, dmats)
        endif
    end subroutine trim_nu_base_bank

    integer function update_full_shell_best_counts( dmat_cand, best_dmat, l_first ) result( nwins )
        real,    intent(in)    :: dmat_cand(:,:,:)
        real,    intent(inout) :: best_dmat(:,:,:)
        logical, intent(in)    :: l_first
        integer :: i, j, k
        nwins = 0
        if( l_first )then
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:nwins) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.nu_lmask(i,j,k) ) cycle
                        best_dmat(i,j,k) = dmat_cand(i,j,k)
                        nwins = nwins + 1
                    end do
                end do
            end do
            !$omp end parallel do
        else
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:nwins) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.nu_lmask(i,j,k) ) cycle
                        if( dmat_cand(i,j,k) < best_dmat(i,j,k) )then
                            best_dmat(i,j,k) = dmat_cand(i,j,k)
                            nwins = nwins + 1
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        endif
    end function update_full_shell_best_counts

    subroutine setup_nu_dmats_full_shell_guarded( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:), best_dmat(:,:,:), no_aux_resolutions(:)
        integer :: i, n_base_max, n_eval, zero_run, nwins, winsz
        real    :: edge_mean, x
        logical :: stopped
        n_base_max = size(cutoff_finds)
        n_eval     = n_base_max
        zero_run   = 0
        stopped    = .false.
        if( n_base_max < 1 ) THROW_HARD('empty cutoff bank; setup_nu_dmats_full_shell_guarded')
        call vol_even_copy_cmat%copy(vol_even)
        call vol_odd_copy_cmat%copy(vol_odd)
        call vol_even_copy_cmat%set_wthreads(.true.)
        call vol_odd_copy_cmat%set_wthreads(.true.)
        winsz = nint(COSMSKHALFWIDTH)
        call vol_even_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_odd_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_even_copy_cmat%fft
        call vol_odd_copy_cmat%fft
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_ft(.true.)
        call vol_odd_filt%set_ft(.true.)
        call vol_even_filt%set_wthreads(.true.)
        call vol_odd_filt%set_wthreads(.true.)
        if( allocated(dmats) ) deallocate(dmats)
        allocate(dmats(ldim(1),ldim(2),ldim(3),n_base_max), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),  source=0.)
        allocate(best_dmat(ldim(1),ldim(2),ldim(3)), source=huge(x))
        do i = 1, n_base_max
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
            if( file_exists(even_cache_fname) .and. file_exists(odd_cache_fname) )then
                call vol_even_filt%read(even_cache_fname)
                call vol_odd_filt%read(odd_cache_fname)
            else
                call vol_even_filt%copy_fast(vol_even_copy_cmat)
                call vol_odd_filt%copy_fast(vol_odd_copy_cmat)
                call vol_even_filt%apply_filter(vol_odd_filt, bwfilters(:,i))
                call vol_even_filt%ifft
                call vol_odd_filt%ifft
                call vol_even_filt%write(even_cache_fname, del_if_exists=.true.)
                call vol_odd_filt%write(odd_cache_fname,  del_if_exists=.true.)
            endif
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmats(:,:,:,i), nu_lmask)
            call tent_smooth_3d(dmats(:,:,:,i), dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
            nwins = update_full_shell_best_counts(dmats(:,:,:,i), best_dmat, i == 1)
            if( i > size(lowpass_limits) )then
                if( nwins == 0 )then
                    zero_run = zero_run + 1
                else
                    zero_run = 0
                endif
                if( zero_run >= NU_FULL_SHELL_ZERO_GUARD )then
                    n_eval  = i
                    stopped = .true.
                    exit
                endif
            endif
        end do
        if( stopped )then
            write(logfhandle,'(A,I0,A,I0,A,I0,A,F6.2,A)') '>>> NU full-shell candidate bank stopped after ', &
                &n_eval, '/', n_base_max, ' candidates (zero guard=', NU_FULL_SHELL_ZERO_GUARD, &
                &', finest LP=', cutoff_find_to_lowpass_limit(n_eval), ' A)'
        endif
        if( n_eval < n_base_max ) call trim_nu_base_bank(n_eval)
        allocate(no_aux_resolutions(0))
        call setup_nu_candidate_coords(size(cutoff_finds), no_aux_resolutions)
        deallocate(no_aux_resolutions)
        call vol_even_copy_cmat%kill
        call vol_odd_copy_cmat%kill
        call vol_even_filt%kill
        call vol_odd_filt%kill
        deallocate(dmat_tmp, best_dmat)
    end subroutine setup_nu_dmats_full_shell_guarded

    subroutine setup_nu_dmats( vol_even, vol_odd, l_mask, aux_resolutions, aux_even, aux_odd, n_highres_steps, &
        &l_full_shell_bank )
        class(image),          intent(in) :: vol_even, vol_odd
        logical,               intent(in) :: l_mask(:,:,:)
        real,                  intent(in) :: aux_resolutions(:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        integer,     optional, intent(in) :: n_highres_steps
        logical,     optional, intent(in) :: l_full_shell_bank
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:)
        integer :: i, n_candidates
        real    :: x
        call init_nu_filter(vol_even, vol_odd, n_highres_steps, l_full_shell_bank)
        if( any(shape(l_mask) /= ldim) ) THROW_HARD('l_mask shape mismatch in setup_nu_dmats')
        if( allocated(nu_lmask) ) deallocate(nu_lmask)
        allocate(nu_lmask(ldim(1),ldim(2),ldim(3)), source=l_mask)
        if( .not. any(nu_lmask) ) THROW_HARD('l_mask has no true voxels in setup_nu_dmats')
        if( present(aux_even) ) then
            if( .not. present(aux_odd) ) THROW_HARD('Auxiliary odd bank missing; setup_nu_dmats')
            if( size(aux_resolutions) /= size(aux_even) ) THROW_HARD('Auxiliary resolutions size mismatch; setup_nu_dmats')
            call stash_aux_volumes(aux_even, aux_odd)
        else
            if( size(aux_resolutions) /= 0 ) THROW_HARD('Auxiliary resolutions supplied without auxiliary volumes; setup_nu_dmats')
            call cleanup_aux_bank
        end if
        if( l_full_shell_bank_active .and. .not.allocated(aux_even_bank) )then
            call setup_nu_dmats_full_shell_guarded(vol_even, vol_odd)
            return
        endif
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call cache_filtered_vols(vol_even, vol_odd)
        if( allocated(dmats) ) deallocate(dmats)
        n_candidates = size(cutoff_finds)
        if( allocated(aux_even_bank) ) n_candidates = n_candidates + size(aux_even_bank)
        call setup_nu_candidate_coords(n_candidates, aux_resolutions)
        allocate(dmats(ldim(1),ldim(2),ldim(3),n_candidates), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),           source=0.)
        do i = 1, size(cutoff_finds)
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
            if( .not.file_exists(even_cache_fname) ) THROW_HARD('Missing filtered volume cache: '//even_cache_fname%to_char())
            if( .not.file_exists(odd_cache_fname)  ) THROW_HARD('Missing filtered volume cache: '//odd_cache_fname%to_char())
            call vol_even_filt%read(even_cache_fname)
            call vol_odd_filt%read(odd_cache_fname)
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmats(:,:,:,i), nu_lmask)
            call tent_smooth_3d(dmats(:,:,:,i), dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
            ! dmat_tmp is never just a temporary buffer, the result is in dmats(:,:,:,i)
        end do
        if( allocated(aux_even_bank) ) then
            do i = 1, size(aux_even_bank)
                call vol_even%nu_objective(aux_even_bank(i), vol_odd, aux_odd_bank(i), &
                    &dmats(:,:,:,size(cutoff_finds)+i), nu_lmask)
                call tent_smooth_3d(dmats(:,:,:,size(cutoff_finds)+i), dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
            end do
        end if
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine setup_nu_dmats

    subroutine setup_nu_candidate_coords( n_candidates, aux_resolutions )
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

    real function lowpass_limit_to_candidate_coord( lp_angstrom )
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

    subroutine optimize_nu_cutoff_finds()
        integer, allocatable :: candmap(:,:,:)
        integer :: nx, ny, nz, i, j, k, icand, best_icand, n_base, n_candidates
        real    :: best_dmat
        if( .not.allocated(dmats) ) THROW_HARD('dmats not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        n_base       = size(cutoff_finds)
        ! dmats is laid out as [base low-pass bank | auxiliary pre-filtered pairs].
        ! The first n_base entries correspond to cutoff_finds(:); any remaining
        ! entries are caller-supplied auxiliary sources.
        n_candidates = size(dmats, 4)
        if( allocated(filtmap) ) deallocate(filtmap)
        if( allocated(srcmap)  ) deallocate(srcmap)
        allocate(candmap(nx,ny,nz), source=1)
        allocate(filtmap(nx,ny,nz), source=1)
        allocate(srcmap(nx,ny,nz),  source=1)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,icand,best_icand,best_dmat) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    if( .not. nu_lmask(i,j,k) )then
                        filtmap(i,j,k) = 1
                        srcmap(i,j,k)  = 1
                        cycle
                    endif
                    best_icand = 1
                    best_dmat = dmats(i,j,k,1)
                    do icand = 2, n_candidates
                        if( dmats(i,j,k,icand) < best_dmat ) then
                            best_dmat = dmats(i,j,k,icand)
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
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(nx,ny,nz), source=dmats(:,:,:,n_base))
        ! this is the big memory consumer, so deallocate it here
        if( allocated(dmats) ) deallocate(dmats)
        deallocate(candmap)
    end subroutine optimize_nu_cutoff_finds

    subroutine candidate_map_to_filt_and_src( candmap, n_base )
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

    subroutine log_nu_candidate_selection_counts( candmap, n_base, stage )
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
            write(logfhandle,'(4X,A8,2X,F7.2,2X,F13.1,2X,I12,2X,F8.2,A)') &
                &auxtag, candidate_coords(icand), cutoff_find_to_lowpass_limit(eff_label), nvox, pct, '%'
        end do
    end subroutine log_nu_candidate_selection_counts

    subroutine log_nu_aux_unary_margin_stats( n_base )
        integer, intent(in) :: n_base
        integer :: iaux, ibase, icand, n_candidates, nmask, nwins, eff_label
        integer :: i, j, k
        real    :: best_base, margin, margin_sum, win_margin_sum, avg_margin, avg_win_margin, pct
        character(len=16) :: auxtag
        if( .not.allocated(dmats) ) return
        if( .not.allocated(nu_lmask) ) return
        if( .not.allocated(candidate_coords) ) return
        n_candidates = size(dmats, 4)
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
            !$omp private(i,j,k,ibase,best_base,margin) reduction(+:nwins,margin_sum,win_margin_sum) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( .not.nu_lmask(i,j,k) ) cycle
                        best_base = dmats(i,j,k,1)
                        do ibase = 2, n_base
                            best_base = min(best_base, dmats(i,j,k,ibase))
                        end do
                        margin = best_base - dmats(i,j,k,icand)
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
            write(logfhandle,'(4X,A8,2X,F7.2,2X,F13.1,2X,I8,2X,F8.2,A,2X,F12.5,2X,F12.5)') &
                &auxtag, candidate_coords(icand), cutoff_find_to_lowpass_limit(eff_label), &
                &nwins, pct, '%', avg_margin, avg_win_margin
        end do
    end subroutine log_nu_aux_unary_margin_stats

    subroutine refine_nu_candidate_map_ordered_labels( candmap, n_candidates )
        integer, intent(inout) :: candmap(:,:,:)
        integer, intent(in)    :: n_candidates
        integer :: iter, color, i, j, k, icand, cur_icand, best_icand, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: nchanged
        integer :: n_base, n_aux
        real    :: beta, e, best_e, site_energy
        if( n_candidates < 2 ) return
        if( .not. allocated(candidate_coords) ) THROW_HARD('candidate_coords not allocated; refine_nu_candidate_map_ordered_labels')
        if( size(candidate_coords) /= n_candidates ) &
            &THROW_HARD('candidate_coords size mismatch; refine_nu_candidate_map_ordered_labels')
        n_base = size(cutoff_finds)
        n_aux  = max(0, n_candidates - n_base)
        beta = estimate_nu_label_smooth_beta(n_candidates)
        write(logfhandle,'(A,ES12.4,A,I0,A,I0,A,I0,A,I0)') '>>> NU ordered-label smoothing: beta=', beta, &
            &', max iterations=', NU_LABEL_SMOOTH_MAXITS, ', candidates=', n_candidates, &
            &', auxiliary=', n_aux, ', step tolerance=', NU_LABEL_SMOOTH_STEP_TOL
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
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do collapse(3) schedule(static) default(shared) &
                !$omp private(i,j,k,icand,cur_icand,best_icand,n_full,nsz,e,best_e) &
                !$omp reduction(+:nchanged) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( .not.nu_lmask(i,j,k) ) cycle
                            if( nu_label_smooth_color(i,j,k) /= color ) cycle
                            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                            cur_icand  = candmap(i,j,k)
                            best_icand = cur_icand
                            best_e     = dmats(i,j,k,cur_icand) + beta * &
                                &nu_label_smooth_neighborhood_cost(cur_icand, candmap, n_full, nsz)
                            do icand = 1, n_candidates
                                if( icand == cur_icand ) cycle
                                e = dmats(i,j,k,icand) + beta * &
                                    &nu_label_smooth_neighborhood_cost(icand, candmap, n_full, nsz)
                                if( nu_label_smooth_is_better(e, best_e) )then
                                    best_e     = e
                                    best_icand = icand
                                endif
                            end do
                            if( best_icand /= cur_icand )then
                                nchanged = nchanged + 1
                                candmap(i,j,k) = best_icand
                            endif
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            site_energy = calc_nu_label_smooth_site_energy(candmap, beta)
            write(logfhandle,'(A,I0,A,I0,A,F12.5)') '>>> NU ordered-label smoothing iteration ', iter, &
                &' changed voxels: ', nchanged, ', mean site energy: ', site_energy
            if( nchanged == 0 ) exit
        end do
    end subroutine refine_nu_candidate_map_ordered_labels

    real function estimate_nu_label_smooth_beta( n_candidates )
        integer, intent(in) :: n_candidates
        integer :: i, j, k, icand, nvox
        real    :: best_e, second_e, cur_e
        estimate_nu_label_smooth_beta = 0.
        nvox = 0
        if( n_candidates < 2 ) return
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    best_e   = huge(best_e)
                    second_e = huge(second_e)
                    do icand = 1, n_candidates
                        cur_e = dmats(i,j,k,icand)
                        if( cur_e < best_e )then
                            second_e = best_e
                            best_e   = cur_e
                        else if( cur_e < second_e )then
                            second_e = cur_e
                        endif
                    end do
                    if( second_e < huge(second_e) )then
                        estimate_nu_label_smooth_beta = estimate_nu_label_smooth_beta + max(0., second_e - best_e)
                        nvox = nvox + 1
                    endif
                end do
            end do
        end do
        if( nvox > 0 ) estimate_nu_label_smooth_beta = &
            &NU_LABEL_SMOOTH_BETA_FRAC * estimate_nu_label_smooth_beta / real(nvox)
    end function estimate_nu_label_smooth_beta

    real function nu_label_smooth_neighborhood_cost( icand, candmap, neigh, nsz )
        integer, intent(in) :: icand, candmap(:,:,:), neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
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
                &nu_label_smooth_pair_cost(icand, candmap(ni,nj,nk))
        end do
        if( degree > 0 ) nu_label_smooth_neighborhood_cost = nu_label_smooth_neighborhood_cost / real(degree)
    end function nu_label_smooth_neighborhood_cost

    real function nu_label_smooth_pair_cost( icand, jcand )
        integer, intent(in) :: icand, jcand
        nu_label_smooth_pair_cost = nu_label_smooth_coord_pair_cost(candidate_coords(icand), candidate_coords(jcand))
    end function nu_label_smooth_pair_cost

    real function nu_label_smooth_coord_pair_cost( icoord, jcoord )
        real, intent(in) :: icoord, jcoord
        real :: step_jump
        if( abs(icoord - jcoord) <= TINY )then
            nu_label_smooth_coord_pair_cost = 0.
        else
            step_jump = max(0., abs(icoord - jcoord) - real(NU_LABEL_SMOOTH_STEP_TOL))
            nu_label_smooth_coord_pair_cost = step_jump + NU_LABEL_SMOOTH_QUAD_FRAC * step_jump * step_jump
        endif
    end function nu_label_smooth_coord_pair_cost

    integer function nu_label_smooth_color( i, j, k )
        integer, intent(in) :: i, j, k
        nu_label_smooth_color = mod(i,2) + 2 * mod(j,2) + 4 * mod(k,2)
    end function nu_label_smooth_color

    logical function nu_label_smooth_is_better( e, best_e )
        real, intent(in) :: e, best_e
        real :: tol
        tol = NU_LABEL_SMOOTH_TIE_EPS * max(1., abs(best_e))
        nu_label_smooth_is_better = e < best_e - tol
    end function nu_label_smooth_is_better

    real function calc_nu_label_smooth_site_energy( candmap, beta )
        integer, intent(in) :: candmap(:,:,:)
        real,    intent(in) :: beta
        integer :: i, j, k, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nvox
        real :: energy_sum
        calc_nu_label_smooth_site_energy = 0.
        energy_sum = 0.
        nvox = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,n_full,nsz) reduction(+:energy_sum,nvox) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                    energy_sum = energy_sum + dmats(i,j,k,candmap(i,j,k)) + beta * &
                        &nu_label_smooth_neighborhood_cost(candmap(i,j,k), candmap, n_full, nsz)
                    nvox = nvox + 1
                end do
            end do
        end do
        !$omp end parallel do
        if( nvox > 0 ) calc_nu_label_smooth_site_energy = energy_sum / real(nvox)
    end function calc_nu_label_smooth_site_energy

    subroutine log_nu_candidate_coords
        integer :: icand
        if( .not.allocated(candidate_coords) ) return
        write(logfhandle,'(A)', advance='no') '>>> NU ordered-label smoothing candidate coordinates:'
        do icand = 1, size(candidate_coords)
            write(logfhandle,'(1X,F6.2)', advance='no') candidate_coords(icand)
        end do
        write(logfhandle,*)
    end subroutine log_nu_candidate_coords

    subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, new_limit, stats )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: new_limit        ! Angstrom limit for the proposed shell
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        type(image)       :: vol_even_filt_new, vol_odd_filt_new
        type(string)      :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_new(:,:,:), dmat_tmp(:,:,:), dmat_finest(:,:,:)
        integer, allocatable :: cutoff_finds_new(:)
        integer           :: new_find, n_finest, n_total, n_extended, sz_old
        real              :: pct_finest, x
        logical, allocatable :: extend_mask(:,:,:), extend_to_new(:,:,:)
        type(nu_highres_extension_stats) :: local_stats
        integer           :: i, j, k
        local_stats%new_limit = new_limit
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats first')
        sz_old    = size(cutoff_finds)
        n_total   = count(nu_lmask)
        local_stats%n_mask = n_total
        if( n_total == 0 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        n_finest  = count(nu_lmask .and. srcmap == 1 .and. filtmap == sz_old)
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
        local_stats%attempted = .true.
        write(logfhandle,'(A,F6.2,A)') '>>> Extending NU filter to ', new_limit, ' A'
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
        call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_new, nu_lmask)
        call tent_smooth_3d(dmat_new, dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
        ! dmat_tmp is never just a temporary buffer, the result is in dmats(:,:,:,i)
        allocate(dmat_finest(ldim(1),ldim(2),ldim(3)), source=huge(x))
        if( allocated(dmat_finest_cached) ) then
            if( all(shape(dmat_finest_cached) == ldim) ) then
                dmat_finest = dmat_finest_cached
            else
                call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
                call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
                call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
                call tent_smooth_3d(dmat_finest, dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
            end if
        else
            call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
            call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
            call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
            call tent_smooth_3d(dmat_finest, dmat_tmp, ldim(1), ldim(2), ldim(3), winsz_tent)
        end if
        ! --- update filtmap in place for the masked voxels ---
        allocate(extend_to_new(ldim(1),ldim(2),ldim(3)), source=.false.)
        n_extended = 0
        call init_nu_highres_extension_selection(extend_mask, dmat_finest, dmat_new, extend_to_new, n_extended)
        call refine_nu_highres_extension_selection(extend_mask, dmat_finest, dmat_new, &
            &extend_to_new, sz_old, n_extended)
        local_stats%n_extended = n_extended
        if( n_finest > 0 ) local_stats%pct_extended_tested = 100. * real(n_extended) / real(n_finest)
        local_stats%applied      = n_extended > 0
        local_stats%promote_next = local_stats%applied .and. &
            &local_stats%pct_extended_tested > NU_HIGHRES_PROMOTION_TESTED_PCT
        write(logfhandle,'(A,I0,A,I0,A,F6.2,A)') '>>> NU high-resolution extension changed ', &
            &n_extended, '/', n_finest, ' tested voxels (', local_stats%pct_extended_tested, '%)'
        if( n_extended == 0 ) then
            call vol_even_filt_new%kill
            call vol_odd_filt_new%kill
            deallocate(extend_mask, extend_to_new, dmat_new, dmat_finest, dmat_tmp)
            if( present(stats) ) stats = local_stats
            return
        end if
        call apply_nu_highres_extension_selection(extend_to_new, sz_old + 1)
        ! --- grow cutoff_finds to include the new level ---
        allocate(cutoff_finds_new(sz_old + 1))
        cutoff_finds_new(:sz_old)  = cutoff_finds
        cutoff_finds_new(sz_old+1) = new_find
        call move_alloc(cutoff_finds_new, cutoff_finds)
        call append_nu_highres_candidate_coord(sz_old, real(sz_old + 1))
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=dmat_new)
        write(logfhandle,'(A,I0,A,F6.2,A)') '>>> Extended ', n_extended, ' voxels to ', new_limit, ' A'
        call vol_even_filt_new%kill
        call vol_odd_filt_new%kill
        deallocate(extend_mask, extend_to_new, dmat_new, dmat_finest, dmat_tmp)
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres

    subroutine extend_nu_filter_highres_shell_next( vol_even, vol_odd, stats )
        class(image), intent(in) :: vol_even, vol_odd
        type(nu_highres_extension_stats), optional, intent(out) :: stats
        type(nu_highres_extension_stats) :: local_stats
        integer :: new_find, sz_old
        real    :: new_limit
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; extend_nu_filter_highres_shell_next')
        sz_old   = size(cutoff_finds)
        new_find = cutoff_finds(sz_old) + 1
        do while( new_find <= box/2 )
            if( .not.any(cutoff_finds == new_find) ) exit
            new_find = new_find + 1
        end do
        if( new_find > box/2 )then
            if( present(stats) ) stats = local_stats
            return
        endif
        new_limit = calc_lowpass_lim(new_find, box, smpd)
        call extend_nu_filter_highres(vol_even, vol_odd, NU_HIGHRES_EXTENSION_THRESHOLD_PCT, &
            &new_limit, stats=local_stats)
        if( present(stats) ) stats = local_stats
    end subroutine extend_nu_filter_highres_shell_next

    subroutine extend_nu_filter_highres_shells( vol_even, vol_odd, nsteps )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(out) :: nsteps
        type(nu_highres_extension_stats) :: step_stats
        integer :: nsteps_local
        nsteps_local = 0
        do
            call extend_nu_filter_highres_shell_next(vol_even, vol_odd, stats=step_stats)
            if( .not. step_stats%attempted ) exit
            if( .not. step_stats%applied   ) exit
            nsteps_local = nsteps_local + 1
            if( .not. step_stats%promote_next ) exit
        end do
        if( present(nsteps) ) nsteps = nsteps_local
    end subroutine extend_nu_filter_highres_shells

    subroutine init_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, n_extended )
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

    subroutine refine_nu_highres_extension_selection( extend_mask, dmat_old, dmat_new, extend_to_new, old_label, n_extended )
        logical, intent(in)    :: extend_mask(:,:,:)
        real,    intent(in)    :: dmat_old(:,:,:), dmat_new(:,:,:)
        logical, intent(inout) :: extend_to_new(:,:,:)
        integer, intent(in)    :: old_label
        integer, intent(out)   :: n_extended
        integer :: iter, color, i, j, k, n_full(3,NU_LABEL_SMOOTH_NNEIGH), nsz, nchanged
        real    :: beta, old_coord, new_coord, e_old, e_new
        beta = estimate_nu_highres_extension_beta(extend_mask, dmat_old, dmat_new)
        old_coord = nu_candidate_coord_for_label(old_label)
        new_coord = real(old_label + 1)
        write(logfhandle,'(A,F10.4,A,I0,A,I0,A,F6.2,A,F6.2)') &
            &'>>> NU high-resolution extension smoothing: beta=', beta, &
            &', max iterations=', NU_LABEL_SMOOTH_MAXITS, ', step tolerance=', NU_LABEL_SMOOTH_STEP_TOL, &
            &', old/new coords=', old_coord, '/', new_coord
        if( beta <= TINY )then
            write(logfhandle,'(A)') '>>> NU high-resolution extension smoothing skipped: beta <= TINY'
            n_extended = count(extend_to_new)
            return
        endif
        do iter = 1, NU_LABEL_SMOOTH_MAXITS
            nchanged = 0
            do color = 0, NU_LABEL_SMOOTH_NCOLORS - 1
                !$omp parallel do collapse(3) schedule(static) default(shared) &
                !$omp private(i,j,k,n_full,nsz,e_old,e_new) reduction(+:nchanged) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( .not.extend_mask(i,j,k) ) cycle
                            if( nu_label_smooth_color(i,j,k) /= color ) cycle
                            call neigh_8_3D(ldim, [i,j,k], n_full, nsz)
                            e_old = dmat_old(i,j,k) + beta * &
                                &nu_highres_extension_neighborhood_cost(old_coord, extend_mask, &
                                &extend_to_new, old_label, new_coord, n_full, nsz)
                            e_new = dmat_new(i,j,k) + beta * &
                                &nu_highres_extension_neighborhood_cost(new_coord, extend_mask, &
                                &extend_to_new, old_label, new_coord, n_full, nsz)
                            if( extend_to_new(i,j,k) )then
                                if( nu_label_smooth_is_better(e_old, e_new) )then
                                    extend_to_new(i,j,k) = .false.
                                    nchanged = nchanged + 1
                                endif
                            else
                                if( nu_label_smooth_is_better(e_new, e_old) )then
                                    extend_to_new(i,j,k) = .true.
                                    nchanged = nchanged + 1
                                endif
                            endif
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
            write(logfhandle,'(A,I0,A,I0)') '>>> NU high-resolution extension smoothing iteration ', iter, &
                &' changed voxels: ', nchanged
            if( nchanged == 0 ) exit
        end do
        n_extended = count(extend_to_new)
    end subroutine refine_nu_highres_extension_selection

    real function estimate_nu_highres_extension_beta( extend_mask, dmat_old, dmat_new )
        logical, intent(in) :: extend_mask(:,:,:)
        real,    intent(in) :: dmat_old(:,:,:), dmat_new(:,:,:)
        integer :: i, j, k, nvox
        estimate_nu_highres_extension_beta = 0.
        nvox = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    estimate_nu_highres_extension_beta = estimate_nu_highres_extension_beta + &
                        &abs(dmat_new(i,j,k) - dmat_old(i,j,k))
                    nvox = nvox + 1
                end do
            end do
        end do
        if( nvox > 0 ) estimate_nu_highres_extension_beta = &
            &NU_LABEL_SMOOTH_BETA_FRAC * estimate_nu_highres_extension_beta / real(nvox)
    end function estimate_nu_highres_extension_beta

    real function nu_highres_extension_neighborhood_cost( icoord, extend_mask, extend_to_new, old_label, &
        &new_coord, neigh, nsz )
        real,    intent(in) :: icoord, new_coord
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        integer, intent(in) :: old_label, neigh(3,NU_LABEL_SMOOTH_NNEIGH), nsz
        integer :: ineigh, ni, nj, nk, degree
        real    :: jcoord
        nu_highres_extension_neighborhood_cost = 0.
        degree = 0
        do ineigh = 1, nsz
            ni = neigh(1,ineigh)
            nj = neigh(2,ineigh)
            nk = neigh(3,ineigh)
            if( .not.nu_lmask(ni,nj,nk) ) cycle
            degree = degree + 1
            jcoord = nu_highres_extension_current_coord(ni, nj, nk, extend_mask, extend_to_new, old_label, new_coord)
            nu_highres_extension_neighborhood_cost = nu_highres_extension_neighborhood_cost + &
                &nu_label_smooth_coord_pair_cost(icoord, jcoord)
        end do
        if( degree > 0 ) nu_highres_extension_neighborhood_cost = &
            &nu_highres_extension_neighborhood_cost / real(degree)
    end function nu_highres_extension_neighborhood_cost

    real function nu_highres_extension_current_coord( i, j, k, extend_mask, extend_to_new, old_label, new_coord )
        integer, intent(in) :: i, j, k, old_label
        logical, intent(in) :: extend_mask(:,:,:), extend_to_new(:,:,:)
        real,    intent(in) :: new_coord
        integer :: aux_label
        if( extend_mask(i,j,k) )then
            if( extend_to_new(i,j,k) )then
                nu_highres_extension_current_coord = new_coord
            else
                nu_highres_extension_current_coord = nu_candidate_coord_for_label(old_label)
            endif
        else
            if( allocated(srcmap) )then
                if( srcmap(i,j,k) > 1 )then
                    aux_label = old_label + srcmap(i,j,k) - 1
                    nu_highres_extension_current_coord = nu_candidate_coord_for_label(aux_label)
                    return
                endif
            endif
            nu_highres_extension_current_coord = nu_candidate_coord_for_label(filtmap(i,j,k))
        endif
    end function nu_highres_extension_current_coord

    real function nu_candidate_coord_for_label( ilabel )
        integer, intent(in) :: ilabel
        if( allocated(candidate_coords) )then
            if( ilabel >= 1 .and. ilabel <= size(candidate_coords) )then
                nu_candidate_coord_for_label = candidate_coords(ilabel)
                return
            endif
        endif
        nu_candidate_coord_for_label = real(ilabel)
    end function nu_candidate_coord_for_label

    integer function nu_effective_base_label_for_candidate( icand, n_base )
        integer, intent(in) :: icand, n_base
        if( n_base < 1 ) THROW_HARD('empty base bank; nu_effective_base_label_for_candidate')
        nu_effective_base_label_for_candidate = nint(nu_candidate_coord_for_label(icand))
        nu_effective_base_label_for_candidate = max(1, min(n_base, nu_effective_base_label_for_candidate))
    end function nu_effective_base_label_for_candidate

    subroutine apply_nu_highres_extension_selection( extend_to_new, new_label )
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

    subroutine append_nu_highres_candidate_coord( old_n_base, new_coord )
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

    subroutine nu_filter_vols( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real(kind=c_float), pointer :: rmat_even_filt(:,:,:), rmat_odd_filt(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_even(:,:,:), rmat_aux_odd(:,:,:)
        integer :: i, j, k, icut, iaux
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_wthreads(.false.)
        call vol_odd_filt%set_wthreads(.false.)
        call vol_even%new(ldim, smpd, wthreads=.false.)
        call vol_odd%new(ldim, smpd, wthreads=.false.)
        call vol_even%get_rmat_ptr(rmat_even_out)
        call vol_odd%get_rmat_ptr(rmat_odd_out)
        rmat_even_out(:ldim(1),:ldim(2),:ldim(3)) = 0.
        rmat_odd_out(:ldim(1),:ldim(2),:ldim(3))  = 0.
        do icut = 1, size(cutoff_finds)
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(icut))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(icut))
            if( .not.file_exists(even_cache_fname) ) THROW_HARD('Missing filtered volume cache: '//even_cache_fname%to_char()//' ; run setup_nu_dmats first')
            if( .not.file_exists(odd_cache_fname)  ) THROW_HARD('Missing filtered volume cache: '//odd_cache_fname%to_char()//' ; run setup_nu_dmats first')
            call vol_even_filt%read(even_cache_fname)
            call vol_odd_filt%read(odd_cache_fname)
            call vol_even_filt%get_rmat_ptr(rmat_even_filt)
            call vol_odd_filt%get_rmat_ptr(rmat_odd_filt)
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == icut ) then
                            rmat_even_out(i,j,k) = rmat_even_filt(i,j,k)
                            rmat_odd_out(i,j,k)  = rmat_odd_filt(i,j,k)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        if( allocated(aux_even_bank) ) then
            do iaux = 1, size(aux_even_bank)
                call aux_even_bank(iaux)%get_rmat_ptr(rmat_aux_even)
                call aux_odd_bank(iaux)%get_rmat_ptr(rmat_aux_odd)
                !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
                do k = 1, ldim(3)
                    do j = 1, ldim(2)
                        do i = 1, ldim(1)
                            if( srcmap(i,j,k) == iaux + 1 ) then
                                rmat_even_out(i,j,k) = rmat_aux_even(i,j,k)
                                rmat_odd_out(i,j,k)  = rmat_aux_odd(i,j,k)
                            end if
                        end do
                    end do
                end do
                !$omp end parallel do
            end do
        end if
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine nu_filter_vols

    subroutine nu_filter_vol( vol_in, vol_out )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
        type(image) :: vol_in_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_out(:,:,:)
        real, allocatable :: bwfilter(:)
        integer :: i, j, k, icut, winsz
        real    :: edge_mean
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vol')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vol')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vol')
        if( any(vol_in%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; nu_filter_vol')
        if( abs(vol_in%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; nu_filter_vol')
        if( any(nu_lmask .and. srcmap /= 1) )then
            THROW_HARD('single-map NU filtering requires a base-bank-only filter map; nu_filter_vol')
        endif
        call vol_in_ft%copy(vol_in)
        call vol_in_ft%set_wthreads(.true.)
        if( .not. vol_in_ft%is_ft() )then
            winsz = nint(COSMSKHALFWIDTH)
            call vol_in_ft%taper_edges_vol(winsz, edge_mean)
            call vol_in_ft%fft
        endif
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_ft(.true.)
        call vol_filt%set_wthreads(.true.)
        call vol_out%new(ldim, smpd, wthreads=.false.)
        call vol_out%get_rmat_ptr(rmat_out)
        rmat_out(:ldim(1),:ldim(2),:ldim(3)) = 0.
        allocate(bwfilter(box), source=0.)
        do icut = 1, size(cutoff_finds)
            call butterworth_filter(cutoff_finds(icut), bwfilter)
            call vol_filt%copy_fast(vol_in_ft)
            call vol_filt%apply_filter(bwfilter)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == icut )then
                            rmat_out(i,j,k) = rmat_filt(i,j,k)
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        deallocate(bwfilter)
        call vol_in_ft%kill
        call vol_filt%kill
    end subroutine nu_filter_vol

    ! statistics about local resolution

     subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
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

    subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, intent(in)             :: mask(:,:,:)
        real, allocatable :: lowpass_vals(:)
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_stats')
        call pack_filtmap_lowpass_limits(lowpass_vals, mask)
        if( size(lowpass_vals) == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_stats')
        call calc_stats(lowpass_vals, statvars)
        deallocate(lowpass_vals)
    end subroutine calc_filtmap_lowpass_stats

    subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
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

    subroutine print_filtmap_lowpass_histogram( mask, aux_resolutions )
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
        write(logfhandle,'(A)')     '    LP limit (A)        Voxels    Pct base'
        do icut = 1, size(cutoff_finds)
            write(logfhandle,'(4X,F10.1,2X,I12,2X,F8.2,A)') cutoff_find_to_lowpass_limit(icut), counts(icut), percentages(icut), '%'
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
                write(logfhandle,'(4X,A8,2X,F14.1,2X,I12,2X,F8.2,A)') auxtag, aux_resolutions(iaux), nvox, pct, '%'
            end do
        end if
        deallocate(counts, percentages)
    end subroutine print_filtmap_lowpass_histogram

    subroutine print_nu_filtmap_lowpass_stats( mask, aux_resolutions )
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

    subroutine analyze_filtmap_neighbor_continuity( mask )
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

end module simple_nu_filter
