!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

contains

    !> the filtered pair of the label field: each half of the source pair (the competition pair itself, or another
    !! pair the field is applied to) Butterworth-filtered per label and scattered by the field, then the auxiliary label
    module subroutine nu_filter_vols( vol_even, vol_odd, vol_src_even, vol_src_odd )
        class(image), intent(out) :: vol_even, vol_odd
        class(image), intent(in)  :: vol_src_even, vol_src_odd
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_even(:,:,:), rmat_aux_odd(:,:,:)
        integer(timer_int_kind) :: t0
        integer :: i, j, k, imask
        if( L_BENCH_GLOB ) t0 = tic()
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nu_filter_vols')
        call release_nu_filter_unary_storage
        call compose_nu_labels_from_volume(vol_src_even, vol_even)
        call compose_nu_labels_from_volume(vol_src_odd,  vol_odd)
        call vol_even%get_rmat_ptr(rmat_even_out)
        call vol_odd%get_rmat_ptr(rmat_odd_out)
        call overlay_aux_label
        if( L_BENCH_GLOB ) nu_rt_compose = nu_rt_compose + toc(t0)

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
    !! then take their own label's Butterworth. The same preparation and filter as
    !! the candidates of setup_nu_dmats (edge taper, FFT), so the voxels are theirs.
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

end submodule simple_nu_filter_apply
