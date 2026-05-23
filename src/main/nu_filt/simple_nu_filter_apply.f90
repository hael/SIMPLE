!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine nu_filter_vols( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
        type(image) :: vol_filt
        type(string) :: cache_fname
        real(kind=c_float), pointer :: rmat_filt(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_even(:,:,:), rmat_aux_odd(:,:,:)
        integer :: i, j, k, icut, iaux, imask
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(srcmap)       ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nu_filter_vols')
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
            cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(icut))
            if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
            call vol_filt%read(cache_fname)
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == icut ) rmat_even_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(1))
        if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
        call vol_filt%read(cache_fname)
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_odd_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        do icut = 2, size(cutoff_finds)
            cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(icut))
            if( .not.file_exists(cache_fname) ) THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
            call vol_filt%read(cache_fname)
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                if( srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == icut ) rmat_odd_out(i,j,k) = rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        if( allocated(aux_even_bank) ) then
            do iaux = 1, size(aux_even_bank)
                call aux_even_bank(iaux)%get_rmat_ptr(rmat_aux_even)
                call aux_odd_bank(iaux)%get_rmat_ptr(rmat_aux_odd)
                !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
                do imask = 1, n_nu_mask
                    i = nu_mask_vox(1,imask)
                    j = nu_mask_vox(2,imask)
                    k = nu_mask_vox(3,imask)
                    if( srcmap(i,j,k) == iaux + 1 ) then
                        rmat_even_out(i,j,k) = rmat_aux_even(i,j,k)
                        rmat_odd_out(i,j,k)  = rmat_aux_odd(i,j,k)
                    end if
                end do
                !$omp end parallel do
            end do
        end if
        call vol_filt%kill
    end subroutine nu_filter_vols

    module subroutine nu_filter_vol( vol_in, vol_out )
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

end submodule simple_nu_filter_apply
