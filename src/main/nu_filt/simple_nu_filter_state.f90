!@descr: simple nu filter state implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_state
implicit none
#include "simple_local_flags.inc"

contains

    module real function cutoff_find_to_lowpass_limit( icut )
        integer, intent(in) :: icut
        if( .not.allocated(cutoff_finds)            ) THROW_HARD('cutoff_finds not allocated; cutoff_find_to_lowpass_limit')
        if( icut < 1 .or. icut > size(cutoff_finds) ) THROW_HARD('cutoff index out of range; cutoff_find_to_lowpass_limit')
        cutoff_find_to_lowpass_limit = calc_lowpass_lim(cutoff_finds(icut), box, smpd)
    end function cutoff_find_to_lowpass_limit

    module real function nu_label_lowpass_limit( ilabel )
        integer, intent(in) :: ilabel
        if( nu_label_is_aux_replacement(ilabel) )then
            nu_label_lowpass_limit = nu_aux_replacement_resolution
        else
            nu_label_lowpass_limit = cutoff_find_to_lowpass_limit(ilabel)
        endif
    end function nu_label_lowpass_limit

    module logical function nu_label_is_aux_replacement( ilabel )
        integer, intent(in) :: ilabel
        nu_label_is_aux_replacement = nu_aux_replacement_label > 0 .and. &
            &ilabel == nu_aux_replacement_label .and. nu_aux_replacement_resolution > TINY .and. &
            &allocated(aux_even_bank) .and. allocated(aux_odd_bank)
    end function nu_label_is_aux_replacement

    module subroutine init_nu_filter( vol_even, vol_odd, n_highres_steps )
        class(image), intent(in) :: vol_even, vol_odd
        integer, optional, intent(in) :: n_highres_steps
        integer, allocatable :: cutoff_finds_tmp(:)
        integer :: i, n_extra, n_extra_requested, n_valid, max_extra, base_find
        integer :: istep, n_extra_retained_requested, n_extra_skip, n_kept_seen
        ldim = vol_even%get_ldim()
        smpd = vol_even%get_smpd()
        box  = ldim(1)
        if( any(vol_odd%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; init_nu_filter')
        if( abs(vol_odd%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; init_nu_filter')
        if( smpd <= TINY ) THROW_HARD('Input volume smpd must be positive; init_nu_filter')
        nu_smooth_norm_radius = -1
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        base_find = calc_fourier_index(lowpass_limits(size(lowpass_limits)), box, smpd)
        n_extra_requested = 0
        if( present(n_highres_steps) )then
            n_extra_requested = min(max(0, n_highres_steps), max(0, box / 2 - base_find))
        endif
        max_extra = max(0, NU_DMAT_CANDIDATE_CAP - size(lowpass_limits) - NU_DMAT_CANDIDATE_HEADROOM)
        n_extra_retained_requested = count_nu_highres_extension_retained_steps(n_extra_requested)
        n_extra = min(n_extra_retained_requested, max_extra)
        n_extra_skip = max(0, n_extra_retained_requested - n_extra)
        allocate(cutoff_finds_tmp(size(lowpass_limits) + n_extra))
        do i = 1, size(lowpass_limits)
            cutoff_finds_tmp(i) = calc_fourier_index(lowpass_limits(i), box, smpd)
        end do
        n_valid = size(lowpass_limits)
        if( NU_DEV_OUTPUT .and. nu_l_report .and. n_extra_retained_requested > n_extra )then
            write(logfhandle,'(A,I0,A,I0,A,I0,A)') &
                &'>>> NU high-resolution depth ', n_extra_requested, &
                &' exceeds distance-matrix memory window; using finest ', n_extra, &
                &' retained shell step(s) within cap ', NU_DMAT_CANDIDATE_CAP, ' candidates'
        endif
        if( NU_DEV_OUTPUT .and. nu_l_report .and. n_extra_requested > 0 .and. &
            &NU_HIGHRES_EXTENSION_RETAIN_STRIDE > 1 )then
            write(logfhandle,'(A,I0,A,I0,A)') &
                &'>>> NU high-resolution extension bank retention: every ', &
                &NU_HIGHRES_EXTENSION_RETAIN_STRIDE, &
                &' shell step(s), plus current terminal shell; requested depth ', &
                &n_extra_requested
        endif
        n_kept_seen = 0
        do istep = 1, n_extra_requested
            if( .not.keep_nu_highres_extension_step(istep, n_extra_requested) ) cycle
            n_kept_seen = n_kept_seen + 1
            if( n_kept_seen <= n_extra_skip ) cycle
            if( .not.any(cutoff_finds_tmp(:n_valid) == base_find + istep) )then
                n_valid = n_valid + 1
                cutoff_finds_tmp(n_valid) = base_find + istep
            endif
        end do
        allocate(cutoff_finds(n_valid))
        cutoff_finds = cutoff_finds_tmp(:n_valid)
        deallocate(cutoff_finds_tmp)
        if( allocated(bwfilters) ) deallocate(bwfilters)
        allocate(bwfilters(box,size(cutoff_finds)), source=0.)
        do i = 1, size(cutoff_finds)
            call butterworth_filter(cutoff_finds(i), bwfilters(:,i))
        end do
    end subroutine init_nu_filter

    module subroutine set_nu_filter_report( l_report )
        logical, intent(in) :: l_report
        nu_l_report = l_report
    end subroutine set_nu_filter_report

    module logical function keep_nu_highres_extension_step( istep, finest_step )
        integer, intent(in) :: istep, finest_step
        keep_nu_highres_extension_step = .false.
        if( istep <= 0 .or. finest_step <= 0 ) return
        if( NU_HIGHRES_EXTENSION_RETAIN_STRIDE <= 1 )then
            keep_nu_highres_extension_step = .true.
        else
            keep_nu_highres_extension_step = &
                &mod(istep, NU_HIGHRES_EXTENSION_RETAIN_STRIDE) == 0 .or. istep == finest_step
        endif
    end function keep_nu_highres_extension_step

    module integer function count_nu_highres_extension_retained_steps( nsteps )
        integer, intent(in) :: nsteps
        integer :: istep
        count_nu_highres_extension_retained_steps = 0
        do istep = 1, max(0, nsteps)
            if( keep_nu_highres_extension_step(istep, nsteps) ) &
                &count_nu_highres_extension_retained_steps = count_nu_highres_extension_retained_steps + 1
        end do
    end function count_nu_highres_extension_retained_steps

    module function filtered_vol_fname( cache_prefix, cutoff_find ) result( fname )
        class(string), intent(in) :: cache_prefix
        integer,       intent(in) :: cutoff_find
        type(string) :: fname
        fname = cache_prefix%to_char()//'_k_'//int2str(cutoff_find)//'.mrc'
    end function filtered_vol_fname

    module logical function filtered_vols_cached( cache_prefix )
        class(string), intent(in) :: cache_prefix
        integer :: i
        filtered_vols_cached = .false.
        if( .not.allocated(cutoff_finds) ) return
        do i = 1, size(cutoff_finds)
            if( nu_label_is_aux_replacement(i) ) cycle
            if( .not.file_exists(filtered_vol_fname(cache_prefix, cutoff_finds(i))) ) return
        end do
        filtered_vols_cached = .true.
    end function filtered_vols_cached

    module subroutine delete_cached_filtered_vols( cache_prefix )
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

    module subroutine release_nu_filter_unary_storage
        if( allocated(dmats_mask)         ) deallocate(dmats_mask)
        if( allocated(dmat_finest_cached) ) deallocate(dmat_finest_cached)
        call release_nu_smooth_norm
    end subroutine release_nu_filter_unary_storage

    module subroutine release_nu_smooth_norm
        if( allocated(nu_smooth_norm)     ) deallocate(nu_smooth_norm)
        nu_smooth_norm_radius = -1
    end subroutine release_nu_smooth_norm

    module subroutine cleanup_nu_filter
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_EVEN))
        call delete_cached_filtered_vols(string(NU_FILTER_CACHE_ODD))
        call release_nu_filter_unary_storage
        if( allocated(bwfilters)          ) deallocate(bwfilters)
        if( allocated(candidate_coords)   ) deallocate(candidate_coords)
        if( allocated(filtmap)            ) deallocate(filtmap)
        if( allocated(cutoff_finds)       ) deallocate(cutoff_finds)
        if( allocated(nu_lmask)           ) deallocate(nu_lmask)
        if( allocated(nu_mask_vox)        ) deallocate(nu_mask_vox)
        call cleanup_aux_bank
        ldim = 0
        box  = 0
        n_nu_mask = 0
        smpd = 0.
        nu_noise_sigma_cached = 0.
    end subroutine cleanup_nu_filter

    module subroutine cleanup_aux_bank
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
        nu_aux_replacement_label = 0
        nu_aux_replacement_resolution = 0.
    end subroutine cleanup_aux_bank

    module subroutine validate_aux_volumes( aux_even, aux_odd )
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

    module subroutine stash_aux_volumes( aux_even, aux_odd )
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

    module subroutine setup_nu_mask_voxels
        integer :: i, j, k, imask
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; setup_nu_mask_voxels')
        if( allocated(nu_mask_vox)    ) deallocate(nu_mask_vox)
        call release_nu_smooth_norm
        n_nu_mask = count(nu_lmask)
        if( n_nu_mask < 1 ) THROW_HARD('l_mask has no true voxels; setup_nu_mask_voxels')
        allocate(nu_mask_vox(3,n_nu_mask), source=0)
        imask = 0
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) cycle
                    imask = imask + 1
                    nu_mask_vox(1,imask) = i
                    nu_mask_vox(2,imask) = j
                    nu_mask_vox(3,imask) = k
                end do
            end do
        end do
        if( imask /= n_nu_mask ) THROW_HARD('mask voxel count mismatch; setup_nu_mask_voxels')
    end subroutine setup_nu_mask_voxels

    module real function nu_objective_smooth_radius_angstrom( lp_angstrom )
        real, intent(in) :: lp_angstrom
        real :: radius_angstrom
        if( smpd <= TINY ) THROW_HARD('invalid smpd; nu_objective_smooth_radius_angstrom')
        radius_angstrom = NU_OBJECTIVE_SMOOTH_RADIUS_FRAC * NU_OBJECTIVE_SMOOTH_AWF * max(lp_angstrom, smpd)
        if( NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A > TINY )then
            radius_angstrom = min(radius_angstrom, NU_OBJECTIVE_SMOOTH_MAX_RADIUS_A)
        endif
        nu_objective_smooth_radius_angstrom = max(smpd, radius_angstrom)
    end function nu_objective_smooth_radius_angstrom

    module integer function nu_objective_smooth_radius_pixels( lp_angstrom )
        real, intent(in) :: lp_angstrom
        nu_objective_smooth_radius_pixels = max(1, nint(nu_objective_smooth_radius_angstrom(lp_angstrom) / smpd))
    end function nu_objective_smooth_radius_pixels

    module subroutine smooth_nu_objective( dmat, tmp, lp_angstrom )
        real, intent(inout) :: dmat(:,:,:)
        real, intent(inout) :: tmp(:,:,:)
        real, intent(in)    :: lp_angstrom
        integer :: i, j, k, radius_px
        real :: x
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; smooth_nu_objective')
        if( any(shape(dmat) /= ldim) ) THROW_HARD('dmat shape mismatch; smooth_nu_objective')
        if( any(shape(tmp)  /= ldim) ) THROW_HARD('tmp shape mismatch; smooth_nu_objective')
        radius_px = nu_objective_smooth_radius_pixels(lp_angstrom)
        call prepare_nu_smooth_norm(radius_px, tmp)
        ! Outside support is a neutral fill for normalized convolution, not a
        ! candidate preference. optimize_nu_cutoff_finds assigns those voxels
        ! to the coarsest filter-bank label after the in-mask competition.
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( .not.nu_lmask(i,j,k) ) dmat(i,j,k) = 0.
                end do
            end do
        end do
        !$omp end parallel do
        call tent_smooth_3d(dmat, tmp, ldim(1), ldim(2), ldim(3), radius_px)
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( nu_lmask(i,j,k) .and. nu_smooth_norm(i,j,k) > TINY )then
                        dmat(i,j,k) = dmat(i,j,k) / nu_smooth_norm(i,j,k)
                    else
                        dmat(i,j,k) = huge(x)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine smooth_nu_objective

    subroutine prepare_nu_smooth_norm( radius_px, tmp )
        integer, intent(in)    :: radius_px
        real,    intent(inout) :: tmp(:,:,:)
        integer :: i, j, k
        if( radius_px < 1 ) THROW_HARD('invalid radius; prepare_nu_smooth_norm')
        if( nu_smooth_norm_radius == radius_px ) return
        if( .not.allocated(nu_smooth_norm) )then
            allocate(nu_smooth_norm(ldim(1),ldim(2),ldim(3)), source=0.)
        endif
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    nu_smooth_norm(i,j,k) = merge(1., 0., nu_lmask(i,j,k))
                end do
            end do
        end do
        !$omp end parallel do
        call tent_smooth_3d(nu_smooth_norm, tmp, ldim(1), ldim(2), ldim(3), radius_px)
        nu_smooth_norm_radius = radius_px
    end subroutine prepare_nu_smooth_norm

    module subroutine pack_nu_dmat_candidate( dmat_full, icand )
        real,    intent(in) :: dmat_full(:,:,:)
        integer, intent(in) :: icand
        integer :: i, j, k, imask
        if( .not.allocated(dmats_mask) ) THROW_HARD('dmats_mask not allocated; pack_nu_dmat_candidate')
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; pack_nu_dmat_candidate')
        if( size(dmats_mask,1) /= n_nu_mask ) THROW_HARD('dmats_mask mask dimension mismatch; pack_nu_dmat_candidate')
        if( icand < 1 .or. icand > size(dmats_mask,2) ) THROW_HARD('candidate index out of range; pack_nu_dmat_candidate')
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            dmats_mask(imask,icand) = dmat_full(i,j,k)
        end do
        !$omp end parallel do
    end subroutine pack_nu_dmat_candidate

    module subroutine unpack_nu_dmat_candidate( icand, dmat_full )
        integer, intent(in)  :: icand
        real,    intent(out) :: dmat_full(:,:,:)
        integer :: i, j, k, imask
        real :: x
        if( .not.allocated(dmats_mask) ) THROW_HARD('dmats_mask not allocated; unpack_nu_dmat_candidate')
        if( .not.allocated(nu_mask_vox) ) THROW_HARD('nu_mask_vox not allocated; unpack_nu_dmat_candidate')
        if( icand < 1 .or. icand > size(dmats_mask,2) ) THROW_HARD('candidate index out of range; unpack_nu_dmat_candidate')
        dmat_full = huge(x)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            dmat_full(i,j,k) = dmats_mask(imask,icand)
        end do
        !$omp end parallel do
    end subroutine unpack_nu_dmat_candidate

    module subroutine cache_filtered_vols( vol_even, vol_odd )
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
            if( nu_label_is_aux_replacement(i) ) cycle
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

    module subroutine generate_single_filtered_pair( vol_even, vol_odd, cutoff_find, even_cache_fname, odd_cache_fname )
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

    module subroutine delete_cached_filtered_pair( cutoff_find )
        integer, intent(in) :: cutoff_find
        type(string) :: cache_fname
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_find)
        if( file_exists(cache_fname) ) call del_file(cache_fname)
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_find)
        if( file_exists(cache_fname) ) call del_file(cache_fname)
        call cache_fname%kill
    end subroutine delete_cached_filtered_pair

end submodule simple_nu_filter_state
