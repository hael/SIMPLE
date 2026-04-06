!@descr: nonuniform filtering of even/odd volumes
!
! A typical call sequence would be:
!    call setup_nu_dmats(vol_even, vol_odd)
!    call optimize_nu_cutoff_finds()
!    call nu_filter_vols(vol_even_filt, vol_odd_filt)
!    call cleanup_nu_filter()
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
implicit none
#include "simple_local_flags.inc"

public :: setup_nu_dmats, optimize_nu_cutoff_finds, nu_filter_vols, cleanup_nu_filter, pack_filtmap_lowpass_limits,&
          calc_filtmap_lowpass_stats, print_filtmap_lowpass_stats, calc_filtmap_lowpass_histogram,&
          print_filtmap_lowpass_histogram
private

real,             parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
integer,          parameter   :: WINSZ_TENT = 1
character(len=*), parameter   :: NU_FILTER_CACHE_EVEN = 'nu_filter_cache_even'
character(len=*), parameter   :: NU_FILTER_CACHE_ODD  = 'nu_filter_cache_odd'
real,             allocatable :: dmats(:,:,:,:)
real,             allocatable :: bwfilters(:,:)
integer,          allocatable :: filtmap(:,:,:)
integer,          allocatable :: cutoff_finds(:)
integer :: ldim(3), box
real    :: smpd

contains

    subroutine init_nu_filter( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        integer :: i
        ldim = vol_even%get_ldim()
        smpd = vol_even%get_smpd()
        box  = ldim(1)
        if( any(vol_odd%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; init_nu_filter')
        if( abs(vol_odd%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; init_nu_filter')
        if( allocated(cutoff_finds) ) deallocate(cutoff_finds)
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
        if( allocated(dmats) )        deallocate(dmats)
        if( allocated(bwfilters) )    deallocate(bwfilters)
        if( allocated(filtmap) )      deallocate(filtmap)
        if( allocated(cutoff_finds) ) deallocate(cutoff_finds)
        ldim = 0
        box  = 0
        smpd = 0.
    end subroutine cleanup_nu_filter

    subroutine cache_filtered_vols( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
        type(string) :: even_cache_fname, odd_cache_fname
        integer :: i, winsz
        real    :: edge_mean
        logical :: even_cached, odd_cached
        call init_nu_filter(vol_even, vol_odd)
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

    subroutine setup_nu_dmats( vol_even, vol_odd )
        class(image),  intent(in) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real, allocatable :: dmat_tmp(:,:,:)
        integer :: i
        real    :: x
        call init_nu_filter(vol_even, vol_odd)
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call cache_filtered_vols(vol_even, vol_odd)
        if( allocated(dmats) ) deallocate(dmats)
        allocate(dmats(ldim(1),ldim(2),ldim(3),size(cutoff_finds)), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),                 source=0.)
        do i = 1, size(cutoff_finds)
            even_cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(i))
            odd_cache_fname  = filtered_vol_fname(string(NU_FILTER_CACHE_ODD),  cutoff_finds(i))
            if( .not.file_exists(even_cache_fname) ) THROW_HARD('Missing filtered volume cache: '//even_cache_fname%to_char())
            if( .not.file_exists(odd_cache_fname)  ) THROW_HARD('Missing filtered volume cache: '//odd_cache_fname%to_char())
            call vol_even_filt%read(even_cache_fname)
            call vol_odd_filt%read(odd_cache_fname)
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmats(:,:,:,i))
            call tent_smooth_3d(dmats(:,:,:,i), dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
        end do
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine setup_nu_dmats

    ! this is where the mask goes in
    subroutine optimize_nu_cutoff_finds
        integer :: nx, ny, nz, i, j, k, icut, best_icut, sz
        real    :: best_dmat
        if( .not.allocated(dmats) ) THROW_HARD('dmats not allocated; run setup_nu_dmats before nonuniform_filter_vol')
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        sz  = size(cutoff_finds)
        if( allocated(filtmap) ) deallocate(filtmap)
        allocate(filtmap(nx,ny,nz), source=1)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k,icut,best_icut,best_dmat) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    best_icut = 1
                    best_dmat = dmats(i,j,k,1)
                    do icut = 2, sz
                        if( dmats(i,j,k,icut) < best_dmat ) then
                            best_dmat = dmats(i,j,k,icut)
                            best_icut = icut
                        end if
                    end do
                    filtmap(i,j,k) = best_icut
                end do
            end do
        end do
        !$omp end parallel do
        ! this is the big memory consumer, so deallocate it here
        if( allocated(dmats) ) deallocate(dmats)
    end subroutine optimize_nu_cutoff_finds

    subroutine nu_filter_vols( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: even_cache_fname, odd_cache_fname
        real(kind=c_float), pointer :: rmat_even_filt(:,:,:), rmat_odd_filt(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        integer :: i, j, k, icut
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
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
                        if( filtmap(i,j,k) == icut ) then
                            rmat_even_out(i,j,k) = rmat_even_filt(i,j,k)
                            rmat_odd_out(i,j,k)  = rmat_odd_filt(i,j,k)
                        end if
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine nu_filter_vols

    ! statistics about local resolution

     subroutine pack_filtmap_lowpass_limits( lowpass_vals, mask )
        real, allocatable, intent(inout) :: lowpass_vals(:)
        logical, optional, intent(in)    :: mask(:,:,:)
        integer :: i, j, k, nvals, ival
        logical :: l_mask_present
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before pack_filtmap_lowpass_limits')
        if( present(mask) ) then
            if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in pack_filtmap_lowpass_limits')
            nvals = count(mask)
            l_mask_present = .true.
        else
            nvals = size(filtmap)
            l_mask_present = .false.
        end if
        if( allocated(lowpass_vals) ) deallocate(lowpass_vals)
        allocate(lowpass_vals(nvals), source=0.)
        ival = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( l_mask_present ) then
                        if( .not.mask(i,j,k) ) cycle
                    end if
                    ival = ival + 1
                    lowpass_vals(ival) = lowpass_limits(filtmap(i,j,k))
                end do
            end do
        end do
    end subroutine pack_filtmap_lowpass_limits

    subroutine calc_filtmap_lowpass_stats( statvars, mask )
        type(stats_struct), intent(out) :: statvars
        logical, optional,  intent(in)  :: mask(:,:,:)
        real, allocatable :: lowpass_vals(:)
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_stats')
        call pack_filtmap_lowpass_limits(lowpass_vals, mask)
        if( size(lowpass_vals) == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_stats')
        call calc_stats(lowpass_vals, statvars)
        deallocate(lowpass_vals)
    end subroutine calc_filtmap_lowpass_stats

    subroutine calc_filtmap_lowpass_histogram( counts, percentages, mask )
        integer,          intent(out) :: counts(:)
        real,             intent(out) :: percentages(:)
        logical, optional, intent(in) :: mask(:,:,:)
        integer :: icut, nselected
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before calc_filtmap_lowpass_histogram')
        if( size(counts) /= size(lowpass_limits) ) THROW_HARD('counts size mismatch in calc_filtmap_lowpass_histogram')
        if( size(percentages) /= size(lowpass_limits) ) THROW_HARD('percentages size mismatch in calc_filtmap_lowpass_histogram')
        if( present(mask) ) then
            if( any(shape(mask) /= shape(filtmap)) ) THROW_HARD('mask shape mismatch in calc_filtmap_lowpass_histogram')
            nselected = count(mask)
        else
            nselected = size(filtmap)
        end if
        if( nselected == 0 ) THROW_HARD('No local resolution values selected in calc_filtmap_lowpass_histogram')
        counts       = 0
        percentages  = 0.
        do icut = 1, size(lowpass_limits)
            if( present(mask) ) then
                counts(icut) = count(filtmap == icut .and. mask)
            else
                counts(icut) = count(filtmap == icut)
            end if
            percentages(icut) = 100. * real(counts(icut)) / real(nselected)
        end do
    end subroutine calc_filtmap_lowpass_histogram

    subroutine print_filtmap_lowpass_histogram( mask, title )
        logical,          optional, intent(in) :: mask(:,:,:)
        character(len=*), optional, intent(in) :: title
        integer :: counts(size(lowpass_limits)), icut
        real    :: percentages(size(lowpass_limits))
        call calc_filtmap_lowpass_histogram(counts, percentages, mask)
        if( present(title) ) then
            write(logfhandle,'(A)') trim(title)
        else
            write(logfhandle,'(A)') '>>> LOCAL RESOLUTION HISTOGRAM'
        end if
        do icut = 1, size(lowpass_limits)
            write(logfhandle,'(F8.3,A,I12,A,F8.3,A)') lowpass_limits(icut), ' A : ', counts(icut), ' voxels, ', percentages(icut), '%'
        end do
    end subroutine print_filtmap_lowpass_histogram

    subroutine print_filtmap_lowpass_stats( mask, title )
        logical, optional, intent(in) :: mask(:,:,:)
        character(len=*), optional, intent(in) :: title
        type(stats_struct) :: statvars
        call calc_filtmap_lowpass_stats(statvars, mask)
        if( present(title) ) then
            write(logfhandle,'(A)') trim(title)
        else
            write(logfhandle,'(A)') '>>> LOCAL RESOLUTION STATS'
        end if
        write(logfhandle,'(A,F8.4)') 'Average: ', statvars%avg
        write(logfhandle,'(A,F8.4)') 'Median : ', statvars%med
        write(logfhandle,'(A,F8.4)') 'Sigma  : ', statvars%sdev
        write(logfhandle,'(A,F8.4)') 'Max    : ', statvars%maxv
        write(logfhandle,'(A,F8.4)') 'Min    : ', statvars%minv
        call print_filtmap_lowpass_histogram(mask, '>>> LOCAL RESOLUTION HISTOGRAM')
    end subroutine print_filtmap_lowpass_stats

end module simple_nu_filter
