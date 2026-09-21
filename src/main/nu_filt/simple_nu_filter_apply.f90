!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine nu_filter_vols( vol_even, vol_odd, vol_apply_even, vol_apply_odd )
        class(image),           intent(out) :: vol_even, vol_odd
        class(image), optional, intent(in)  :: vol_apply_even, vol_apply_odd
        type(image) :: vol_filt
        type(string) :: cache_fname
        real(kind=c_float), pointer :: rmat_filt(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_even(:,:,:), rmat_aux_odd(:,:,:)
        integer :: i, j, k, icut, imask
        logical :: l_apply
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nu_filter_vols')
        l_apply = present(vol_apply_even) .or. present(vol_apply_odd)
        if( l_apply )then
            if( .not.(present(vol_apply_even) .and. present(vol_apply_odd)) ) &
                &THROW_HARD('an apply pair needs both halves; nu_filter_vols')
            ! the label field of the competition pair applied to another pair
            ! (the solvent-prior'd base pair): per-label Butterworth of the
            ! apply halves, scattered by the label field
            call release_nu_filter_unary_storage
            call compose_nu_labels_from_volume(vol_apply_even, vol_even)
            call compose_nu_labels_from_volume(vol_apply_odd,  vol_odd)
            call vol_even%get_rmat_ptr(rmat_even_out)
            call vol_odd%get_rmat_ptr(rmat_odd_out)
            call overlay_aux_label
            return
        endif
        call release_nu_filter_unary_storage
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_wthreads(.false.)
        call vol_even%new(ldim, smpd, wthreads=.false.)
        call vol_odd%new(ldim, smpd, wthreads=.false.)
        call vol_even%get_rmat_ptr(rmat_even_out)
        call vol_odd%get_rmat_ptr(rmat_odd_out)
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(1))
        if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
        call vol_filt%read(cache_fname)
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_even_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        do icut = 2, size(cutoff_finds)
            if( nu_label_is_aux_replacement(icut) ) cycle
            cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(icut))
            if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
            call vol_filt%read(cache_fname)
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( filtmap(i,j,k) == icut ) rmat_even_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(1))
        if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
        call vol_filt%read(cache_fname)
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_odd_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        do icut = 2, size(cutoff_finds)
            if( nu_label_is_aux_replacement(icut) ) cycle
            cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(icut))
            if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
            call vol_filt%read(cache_fname)
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( filtmap(i,j,k) == icut ) rmat_odd_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        call overlay_aux_label
        call vol_filt%kill

    contains

        subroutine overlay_aux_label
            if( nu_aux_replacement_label > 0 ) then
                if( .not.allocated(aux_even_bank) .or. .not.allocated(aux_odd_bank) ) &
                    &THROW_HARD('missing NU auxiliary replacement volumes; nu_filter_vols')
                call aux_even_bank(1)%get_rmat_ptr(rmat_aux_even)
                call aux_odd_bank(1)%get_rmat_ptr(rmat_aux_odd)
                !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
                do imask = 1, n_nu_mask
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( int(filtmap(i,j,k)) == nu_aux_replacement_label ) then
                        rmat_even_out(i,j,k) = rmat_aux_even(i,j,k)
                        rmat_odd_out(i,j,k)  = rmat_aux_odd(i,j,k)
                    end if
                end do
                !$omp end parallel do
            end if
        end subroutine overlay_aux_label

    end subroutine nu_filter_vols

    !> One volume filtered per label of the current field (no auxiliary
    !! fill): the coarsest filter seeds the whole grid, the mask-packed voxels
    !! then take their own label's Butterworth. The same preparation as
    !! nu_filter_vol (edge taper, FFT).
    subroutine compose_nu_labels_from_volume( vol_in, vol_out )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
        type(image) :: vol_in_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_out(:,:,:)
        real, allocatable :: bwfilter(:)
        integer :: i, j, k, icut, imask, winsz
        real    :: edge_mean
        if( any(vol_in%get_ldim() /= ldim)       ) THROW_HARD('apply volume dimensions differ; compose_nu_labels_from_volume')
        if( abs(vol_in%get_smpd() - smpd) > TINY ) THROW_HARD('apply volume smpd differs; compose_nu_labels_from_volume')
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
        allocate(bwfilter(box), source=0.)
        call butterworth_filter(cutoff_finds(1), bwfilter)
        call vol_filt%copy_fast(vol_in_ft)
        call vol_filt%apply_filter(bwfilter)
        call vol_filt%ifft
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        do icut = 2, size(cutoff_finds)
            if( nu_label_is_aux_replacement(icut) ) cycle
            call butterworth_filter(cutoff_finds(icut), bwfilter)
            call vol_filt%copy_fast(vol_in_ft)
            call vol_filt%apply_filter(bwfilter)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( filtmap(i,j,k) == icut ) rmat_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        deallocate(bwfilter)
        call vol_in_ft%kill
        call vol_filt%kill
    end subroutine compose_nu_labels_from_volume

    module subroutine nu_filter_vol( vol_in, vol_out )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
        type(image) :: vol_in_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_out(:,:,:)
        real, allocatable :: bwfilter(:)
        integer :: i, j, k, icut, imask, winsz
        real    :: edge_mean
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vol')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vol')
        if( .not.allocated(nu_lmask)     ) THROW_HARD('nu_lmask not allocated; run setup_nu_dmats before nu_filter_vol')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nu_filter_vol')
        if( any(vol_in%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; nu_filter_vol')
        if( abs(vol_in%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; nu_filter_vol')
        if( nu_aux_replacement_label > 0 )then
            if( any(nu_lmask .and. filtmap == int(nu_aux_replacement_label, kind=NU_LABEL_KIND)) )then
                THROW_HARD('single-map NU filtering cannot synthesize an auxiliary replacement label; nu_filter_vol')
            endif
        endif
        call release_nu_filter_unary_storage
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
        allocate(bwfilter(box), source=0.)
        ! Seed the output (including outside-mask voxels) with the coarsest
        ! filter, matching nu_filter_vols semantics, then scatter only the
        ! mask-packed voxels for the remaining bank entries.
        call butterworth_filter(cutoff_finds(1), bwfilter)
        call vol_filt%copy_fast(vol_in_ft)
        call vol_filt%apply_filter(bwfilter)
        call vol_filt%ifft
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        do icut = 2, size(cutoff_finds)
            if( nu_label_is_aux_replacement(icut) ) cycle
            call butterworth_filter(cutoff_finds(icut), bwfilter)
            call vol_filt%copy_fast(vol_in_ft)
            call vol_filt%apply_filter(bwfilter)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( filtmap(i,j,k) == icut ) rmat_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        deallocate(bwfilter)
        call vol_in_ft%kill
        call vol_filt%kill
    end subroutine nu_filter_vol

end submodule simple_nu_filter_apply
