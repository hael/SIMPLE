!@descr: nonuniform filtering of even/odd volumes
!
! A typical call sequence would be:
!    call setup_nu_dmats(vol_even, vol_odd, l_mask)
!    call optimize_nu_cutoff_finds()
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! Updated call sequence with iterative high-resolution refinement
!    call setup_nu_dmats(vol_even, vol_odd, l_mask)
!    call optimize_nu_cutoff_finds()
!    call extend_nu_filter_highres_iterative(vol_even, vol_odd)  ! optional
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
! supports auxiliary candidate pairs that can compete with the
! base low-pass bank during voxelwise optimization
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
implicit none
#include "simple_local_flags.inc"

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_nu_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram, extend_nu_filter_highres, extend_nu_filter_highres_iterative,&
          analyze_filtmap_neighbor_continuity
private

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
real,             parameter   :: EXTRA_LIMITS(3)   = [3.5, 3.0, 2.5]
integer,          parameter   :: WINSZ_TENT = 3
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD  = 'nu_filter_cache_odd'
real,             allocatable :: dmats(:,:,:,:)
real,             allocatable :: bwfilters(:,:)
integer,          allocatable :: filtmap(:,:,:)
integer,          allocatable :: srcmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
real,             allocatable :: dmat_finest_cached(:,:,:)
logical,          allocatable :: nu_lmask(:,:,:)
type(image),      allocatable :: aux_even_bank(:), aux_odd_bank(:)
integer :: ldim(3), box
real    :: smpd

contains

    real function cutoff_find_to_lowpass_limit( icut )
        integer, intent(in) :: icut
        if( .not.allocated(cutoff_finds)            ) THROW_HARD('cutoff_finds not allocated; cutoff_find_to_lowpass_limit')
        if( icut < 1 .or. icut > size(cutoff_finds) ) THROW_HARD('cutoff index out of range; cutoff_find_to_lowpass_limit')
        cutoff_find_to_lowpass_limit = calc_lowpass_lim(cutoff_finds(icut), box, smpd)
    end function cutoff_find_to_lowpass_limit

    subroutine init_nu_filter( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        integer :: i
        ldim = vol_even%get_ldim()
        smpd = vol_even%get_smpd()
        box  = ldim(1)
        if( any(vol_odd%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; init_nu_filter')
        if( abs(vol_odd%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; init_nu_filter')
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        allocate(cutoff_finds(size(lowpass_limits)))
        do i = 1, size(lowpass_limits)
            cutoff_finds(i) = calc_fourier_index(lowpass_limits(i), box, smpd)
        end do
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
        if( allocated(filtmap)            ) deallocate(filtmap)
        if( allocated(srcmap)             ) deallocate(srcmap)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( allocated(nu_lmask)           ) deallocate(nu_lmask)
        call cleanup_aux_bank
        ldim = 0
        box  = 0
        smpd = 0.
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

    subroutine setup_nu_dmats( vol_even, vol_odd, l_mask, aux_even, aux_odd )
        class(image),          intent(in) :: vol_even, vol_odd
        logical,               intent(in) :: l_mask(:,:,:)
        type(image), optional, intent(in) :: aux_even(:), aux_odd(:)
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:)
        integer :: i, n_candidates
        real    :: x
        call init_nu_filter(vol_even, vol_odd)
        if( any(shape(l_mask) /= ldim) ) THROW_HARD('l_mask shape mismatch in setup_nu_dmats')
        if( allocated(nu_lmask) ) deallocate(nu_lmask)
        allocate(nu_lmask(ldim(1),ldim(2),ldim(3)), source=l_mask)
        if( .not. any(nu_lmask) ) THROW_HARD('l_mask has no true voxels in setup_nu_dmats')
        if( present(aux_even) ) then
            if( .not. present(aux_odd) ) THROW_HARD('Auxiliary odd bank missing; setup_nu_dmats')
            call stash_aux_volumes(aux_even, aux_odd)
        else
            call cleanup_aux_bank
        end if
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call cache_filtered_vols(vol_even, vol_odd)
        if( allocated(dmats) ) deallocate(dmats)
        n_candidates = size(cutoff_finds)
        if( allocated(aux_even_bank) ) n_candidates = n_candidates + size(aux_even_bank)
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
            call tent_smooth_3d(dmats(:,:,:,i), dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
            ! dmat_tmp is never just a temporary buffer, the result is in dmats(:,:,:,i)
        end do
        if( allocated(aux_even_bank) ) then
            do i = 1, size(aux_even_bank)
                call vol_even%nu_objective(aux_even_bank(i), vol_odd, aux_odd_bank(i), &
                    &dmats(:,:,:,size(cutoff_finds)+i), nu_lmask)
                call tent_smooth_3d(dmats(:,:,:,size(cutoff_finds)+i), dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
            end do
        end if
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine setup_nu_dmats

    subroutine optimize_nu_cutoff_finds
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
                    ! Base-bank winners preserve their low-pass index in filtmap.
                    ! Auxiliary winners are tracked through srcmap only, because
                    ! they do not correspond to one of cutoff_finds(:).
                    if( best_icand <= n_base ) then
                        srcmap(i,j,k)  = 1
                        filtmap(i,j,k) = best_icand
                    else
                        ! srcmap numbering:
                        !   1   -> base low-pass bank
                        !   2+  -> auxiliary pair 1, 2, ...
                        srcmap(i,j,k)  = best_icand - n_base + 1
                        filtmap(i,j,k) = 1
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(nx,ny,nz), source=dmats(:,:,:,n_base))
        ! this is the big memory consumer, so deallocate it here
        if( allocated(dmats) ) deallocate(dmats)
    end subroutine optimize_nu_cutoff_finds

    subroutine extend_nu_filter_highres( vol_even, vol_odd, threshold_pct, new_limit )
        class(image), intent(in) :: vol_even, vol_odd
        real,         intent(in) :: threshold_pct   ! e.g. 10.0
        real,         intent(in) :: new_limit        ! e.g. 3.5 Angstroms
        type(image)       :: vol_even_filt_new, vol_odd_filt_new
        type(string)      :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_new(:,:,:), dmat_tmp(:,:,:), dmat_finest(:,:,:)
        integer, allocatable :: cutoff_finds_new(:)
        integer           :: new_find, n_finest, n_total, n_extended, sz_old
        real              :: pct_finest, x
        logical, allocatable :: extend_mask(:,:,:)
        integer           :: i, j, k
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds first')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated')
        sz_old    = size(cutoff_finds)
        n_total   = size(filtmap)
        n_finest  = count(srcmap == 1 .and. filtmap == sz_old)
        pct_finest = 100. * real(n_finest) / real(n_total)
        if( pct_finest < threshold_pct ) return   ! trigger not met, nothing to do
        new_find = calc_fourier_index(new_limit, box, smpd)
        if( new_find <= cutoff_finds(sz_old) ) return
        if( any(cutoff_finds == new_find) ) return
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
                extend_mask(i,j,k) = (srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == sz_old)
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
        call tent_smooth_3d(dmat_new, dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
        ! dmat_tmp is never just a temporary buffer, the result is in dmats(:,:,:,i)
        allocate(dmat_finest(ldim(1),ldim(2),ldim(3)), source=huge(x))
        if( allocated(dmat_finest_cached) ) then
            if( all(shape(dmat_finest_cached) == ldim) ) then
                dmat_finest = dmat_finest_cached
            else
                call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
                call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
                call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
                call tent_smooth_3d(dmat_finest, dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
            end if
        else
            call vol_even_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(sz_old)))
            call vol_odd_filt_new%read(filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(sz_old)))
            call vol_even%nu_objective(vol_even_filt_new, vol_odd, vol_odd_filt_new, dmat_finest, nu_lmask)
            call tent_smooth_3d(dmat_finest, dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
        end if
        ! --- update filtmap in place for the masked voxels ---
        n_extended = 0
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) reduction(+:n_extended)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.extend_mask(i,j,k) ) cycle
                    if( dmat_new(i,j,k) < dmat_finest(i,j,k) ) then
                        srcmap(i,j,k)  = 1
                        filtmap(i,j,k) = sz_old + 1
                        n_extended = n_extended + 1
                    end if
                end do
            end do
        end do
        !$omp end parallel do
        if( n_extended == 0 ) then
            call vol_even_filt_new%kill
            call vol_odd_filt_new%kill
            deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
            return
        end if
        ! --- grow cutoff_finds to include the new level ---
        allocate(cutoff_finds_new(sz_old + 1))
        cutoff_finds_new(:sz_old)  = cutoff_finds
        cutoff_finds_new(sz_old+1) = new_find
        call move_alloc(cutoff_finds_new, cutoff_finds)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        allocate(dmat_finest_cached(ldim(1),ldim(2),ldim(3)), source=dmat_new)
        write(logfhandle,'(A,I0,A,F6.2,A)') '>>> Extended ', n_extended, ' voxels to ', new_limit, ' A'
        call vol_even_filt_new%kill
        call vol_odd_filt_new%kill
        deallocate(extend_mask, dmat_new, dmat_finest, dmat_tmp)
    end subroutine extend_nu_filter_highres

    subroutine extend_nu_filter_highres_iterative( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        integer :: i
        do i = 1, size(EXTRA_LIMITS)
            call extend_nu_filter_highres(vol_even, vol_odd, 10.0, EXTRA_LIMITS(i))
        end do
    end subroutine extend_nu_filter_highres_iterative

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
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_histogram')
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before calc_filtmap_lowpass_histogram')
        if( size(counts) /= size(cutoff_finds) ) THROW_HARD('counts size mismatch in calc_filtmap_lowpass_histogram')
        if( size(percentages) /= size(cutoff_finds) ) THROW_HARD('percentages size mismatch in calc_filtmap_lowpass_histogram')
        if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in calc_filtmap_lowpass_histogram')
        if( allocated(srcmap) ) then
            nselected = count(mask .and. srcmap == 1)
        else
            nselected = count(mask)
        end if
        if( nselected == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_histogram')
        counts       = 0
        percentages  = 0.
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
        write(logfhandle,'(A)') '>>> LOCAL RESOLUTION HISTOGRAM'
        do icut = 1, size(cutoff_finds)
            write(logfhandle,'(A8,1X,F8.1,A,I12,A,F8.1,A)') '', cutoff_find_to_lowpass_limit(icut), ' A : ', counts(icut), ' voxels, ', percentages(icut), '%'
        end do
        ! Print auxiliary pair assignments if present
        if( allocated(aux_even_bank) .and. present(aux_resolutions) ) then
            if( size(aux_resolutions) /= size(aux_even_bank) ) THROW_HARD('aux_resolutions size mismatch in print_filtmap_lowpass_histogram')
            ! Determine nselected for percentage calculation
            nselected = count(mask)
            do iaux = 1, size(aux_even_bank)
                nvox = count(srcmap == iaux + 1 .and. mask)
                pct = 100. * real(nvox) / real(nselected)
                write(auxtag,'(A,I0,A)') 'Aux', iaux, '@'
                write(logfhandle,'(A8,1X,F8.1,A,I12,A,F8.1,A)') auxtag, aux_resolutions(iaux), ' A : ', nvox, ' voxels, ', pct, '%'
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
                return
            end if
        end if
        call calc_filtmap_lowpass_stats(statvars, mask)
        write(logfhandle,'(A)') '>>> LOCAL RESOLUTION STATS'
        write(logfhandle,'(A,F8.4)') 'Average: ', statvars%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', statvars%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', statvars%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', statvars%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', statvars%minv
        call print_filtmap_lowpass_histogram(mask, aux_resolutions)
    end subroutine print_nu_filtmap_lowpass_stats

    subroutine analyze_filtmap_neighbor_continuity( mask, discontinuity_threshold )
        logical, intent(in) :: mask(:,:,:)
        integer, optional, intent(in) :: discontinuity_threshold
        integer :: i, j, k, di, dj, dk, ni, nj, nk
        integer :: lp_i, lp_j, lp_diff, max_diff, n_neighbors, n_discontinuous_neighbors
        integer :: n_total_neighbor_pairs, n_discontinuous_pairs, n_voxels_with_discontinuity, nx, ny, nz, ii, thresh
        integer, allocatable :: discontinuity_counts(:)
        real :: pct
        if( .not.allocated(filtmap)            ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before analyze_filtmap_neighbor_continuity')
        if( .not.allocated(cutoff_finds)       ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before analyze_filtmap_neighbor_continuity')
        if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in analyze_filtmap_neighbor_continuity')
        ! Use threshold of 1 by default (neighbors differ by more than 1 step)
        thresh = 1
        if( present(discontinuity_threshold) ) thresh = discontinuity_threshold
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        allocate(discontinuity_counts(8), source=0)
        n_voxels_with_discontinuity = 0
        n_total_neighbor_pairs      = 0
        n_discontinuous_pairs       = 0
        ! Iterate through all voxels
        !$omp parallel do collapse(3) schedule(static) default(shared) &
        !$omp private(i,j,k,di,dj,dk,ni,nj,nk,lp_i,lp_j,lp_diff,n_neighbors,n_discontinuous_neighbors,max_diff) &
        !$omp reduction(+:n_total_neighbor_pairs, n_discontinuous_pairs, n_voxels_with_discontinuity, discontinuity_counts)
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
                                if( lp_diff > thresh ) then
                                    n_discontinuous_neighbors = n_discontinuous_neighbors + 1
                                    n_discontinuous_pairs     = n_discontinuous_pairs + 1
                                end if
                                max_diff = max(max_diff, lp_diff)
                                if( lp_diff >= 1 .and. lp_diff <= 8 ) then
                                    discontinuity_counts(lp_diff) = discontinuity_counts(lp_diff) + 1
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
        ! Print diagnostics
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> NONUNIFORM FILTER NEIGHBOR CONTINUITY ANALYSIS'
        write(logfhandle,'(A,I0,A)') '>>> Threshold: neighbors differing by > ', thresh, ' step(s) in LP ordering'
        write(logfhandle,'(A)') ''
        ! Count total masked voxels
        if( allocated(srcmap) ) then
            ii = count(mask .and. srcmap == 1)
        else
            ii = count(mask)
        end if
        write(logfhandle,'(A,I12)') 'Total voxels analyzed:                   ', ii
        write(logfhandle,'(A,I12)') 'Voxels with discontinuous neighbors:     ', n_voxels_with_discontinuity
        if( ii > 0 ) then
            pct = 100. * real(n_voxels_with_discontinuity) / real(ii)
            write(logfhandle,'(A,F8.2,A)') 'Percentage of voxels with discontinuity: ', pct, '%'
        end if
        write(logfhandle,'(A)') ''
        if( n_total_neighbor_pairs > 0 ) then
            pct = 100. * real(n_discontinuous_pairs) / real(n_total_neighbor_pairs)
            write(logfhandle,'(A,I14)') 'Total neighbor pairs examined:            ', n_total_neighbor_pairs
            write(logfhandle,'(A,I0,A,I14)') 'Neighbor pairs with discontinuity (>', thresh, ' step): ', n_discontinuous_pairs
            write(logfhandle,'(A,F8.2,A)') 'Percentage of discontinuous pairs:        ', pct, '%'
        end if
        write(logfhandle,'(A)') ''
        ! Print distribution of all pair step-differences
        write(logfhandle,'(A)') '>>> DISCONTINUITY MAGNITUDE DISTRIBUTION (all LP step differences)'
        do ii = 1, 8
            if( discontinuity_counts(ii) > 0 ) then
                pct = 100. * real(discontinuity_counts(ii)) / real(n_total_neighbor_pairs)
                write(logfhandle,'(A,I1,A,I14,A,F7.3,A)') &
                    'LP step diff = ', ii, ': ', discontinuity_counts(ii), ' pairs (', pct, '%)'
            end if
        end do
        write(logfhandle,'(A)') ''
        write(logfhandle,'(A)') '>>> INTERPRETATION'
        if( n_voxels_with_discontinuity > 0 ) then
            pct = 100. * real(n_voxels_with_discontinuity) / real(max(ii,1))
            if( pct > 50. ) then
                write(logfhandle,'(A)') '>>> Very high discontinuity rate — tent regularization kernel likely too narrow'
                write(logfhandle,'(A,I0,A)') '    Current WINSZ_TENT = ', WINSZ_TENT, '; consider increasing to reduce spatial noise in LP map'
            else if( pct > 20. ) then
                write(logfhandle,'(A)') '>>> Moderate discontinuity rate — tent regularization may benefit from widening'
            else
                write(logfhandle,'(A)') '>>> Low discontinuity rate — local resolution map is spatially smooth'
            end if
        else
            write(logfhandle,'(A)') '>>> No discontinuities found — local resolution map is perfectly smooth'
        end if
        write(logfhandle,'(A)') ''
        deallocate(discontinuity_counts)
    end subroutine analyze_filtmap_neighbor_continuity

end module simple_nu_filter
