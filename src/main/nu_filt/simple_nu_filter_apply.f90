!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine nu_filter_vols( vol_even, vol_odd, soft_synthesis )
        class(image), intent(out) :: vol_even, vol_odd
        logical, optional, intent(in) :: soft_synthesis
        type(image) :: vol_filt
        type(string) :: cache_fname
        real(kind=c_float), pointer :: rmat_filt(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:),  rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_even(:,:,:), rmat_aux_odd(:,:,:)
        integer :: i, j, k, icut, imask
        logical :: l_soft_synthesis
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_filter_vols')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_filter_vols')
        if( .not.allocated(nu_mask_vox)  ) THROW_HARD('nu_mask_vox not allocated; run setup_nu_dmats before nu_filter_vols')
        call release_nu_filter_unary_storage
        l_soft_synthesis = .false.
        if( present(soft_synthesis) ) l_soft_synthesis = soft_synthesis
        if( l_soft_synthesis )then
            call synthesize_nu_filter_vols_soft(vol_even, vol_odd)
            return
        endif
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
        call vol_filt%kill
    end subroutine nu_filter_vols

    module subroutine nu_filter_vol( vol_in, vol_out, soft_synthesis )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_out
        logical, optional, intent(in) :: soft_synthesis
        type(image) :: vol_in_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_out(:,:,:)
        real, allocatable :: bwfilter(:)
        integer :: i, j, k, icut, imask, winsz
        real    :: edge_mean
        logical :: l_soft_synthesis
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
        l_soft_synthesis = .false.
        if( present(soft_synthesis) ) l_soft_synthesis = soft_synthesis
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
        if( l_soft_synthesis )then
            call synthesize_nu_filter_vol_soft(vol_in_ft, vol_filt, bwfilter, vol_out)
            deallocate(bwfilter)
            call vol_in_ft%kill
            call vol_filt%kill
            return
        endif
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

    subroutine synthesize_nu_filter_vols_soft( vol_even, vol_odd )
        class(image), intent(out) :: vol_even, vol_odd
        type(image) :: vol_even_filt, vol_odd_filt
        type(string) :: cache_fname
        real, allocatable :: synth_coord(:,:,:), synth_tmp(:,:,:), synth_norm(:,:,:)
        real(kind=c_float), pointer :: rmat_even_out(:,:,:), rmat_odd_out(:,:,:)
        real(kind=c_float), pointer :: rmat_even_src(:,:,:), rmat_odd_src(:,:,:)
        integer :: i, j, k, icut, imask, nlabels
        real    :: weight
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; synthesize_nu_filter_vols_soft')
        nlabels = size(cutoff_finds)
        if( nlabels < 1 ) THROW_HARD('empty NU filter bank; synthesize_nu_filter_vols_soft')
        call vol_even%new(ldim, smpd, wthreads=.false.)
        call vol_odd%new(ldim, smpd, wthreads=.false.)
        call vol_even%get_rmat_ptr(rmat_even_out)
        call vol_odd%get_rmat_ptr(rmat_odd_out)
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_wthreads(.false.)
        call vol_odd_filt%set_wthreads(.false.)
        allocate(synth_coord(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(synth_tmp(ldim(1),ldim(2),ldim(3)),   source=0.)
        allocate(synth_norm(ldim(1),ldim(2),ldim(3)),  source=0.)
        call build_nu_synthesis_coord(synth_coord, synth_tmp, synth_norm)
        call log_nu_soft_synthesis()
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(1))
        if( .not.file_exists(cache_fname) ) &
            &THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
        call vol_even_filt%read(cache_fname)
        call vol_even_filt%get_rmat_ptr(rmat_even_src)
        rmat_even_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_even_src(:ldim(1),:ldim(2),:ldim(3))
        cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(1))
        if( .not.file_exists(cache_fname) ) &
            &THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
        call vol_odd_filt%read(cache_fname)
        call vol_odd_filt%get_rmat_ptr(rmat_odd_src)
        rmat_odd_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_odd_src(:ldim(1),:ldim(2),:ldim(3))
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            rmat_even_out(i,j,k) = 0.
            rmat_odd_out(i,j,k)  = 0.
        end do
        !$omp end parallel do
        do icut = 1, nlabels
            if( nu_label_is_aux_replacement(icut) )then
                call aux_even_bank(1)%get_rmat_ptr(rmat_even_src)
                call aux_odd_bank(1)%get_rmat_ptr(rmat_odd_src)
            else
                cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_EVEN), cutoff_finds(icut))
                if( .not.file_exists(cache_fname) ) &
                    &THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
                call vol_even_filt%read(cache_fname)
                call vol_even_filt%get_rmat_ptr(rmat_even_src)
                cache_fname = filtered_vol_fname(string(NU_FILTER_CACHE_ODD), cutoff_finds(icut))
                if( .not.file_exists(cache_fname) ) &
                    &THROW_HARD('Missing filtered volume cache: '//cache_fname%to_char()//' ; run setup_nu_dmats first')
                call vol_odd_filt%read(cache_fname)
                call vol_odd_filt%get_rmat_ptr(rmat_odd_src)
            endif
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,weight) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                weight = nu_synthesis_weight(synth_coord(i,j,k), icut)
                if( weight > TINY )then
                    rmat_even_out(i,j,k) = rmat_even_out(i,j,k) + weight * rmat_even_src(i,j,k)
                    rmat_odd_out(i,j,k)  = rmat_odd_out(i,j,k)  + weight * rmat_odd_src(i,j,k)
                endif
            end do
            !$omp end parallel do
        end do
        deallocate(synth_coord, synth_tmp, synth_norm)
        call vol_even_filt%kill
        call vol_odd_filt%kill
    end subroutine synthesize_nu_filter_vols_soft

    subroutine synthesize_nu_filter_vol_soft( vol_in_ft, vol_filt, bwfilter, vol_out )
        class(image), intent(in)    :: vol_in_ft
        class(image), intent(inout) :: vol_filt, vol_out
        real,         intent(inout) :: bwfilter(:)
        real, allocatable :: synth_coord(:,:,:), synth_tmp(:,:,:), synth_norm(:,:,:)
        real(kind=c_float), pointer :: rmat_out(:,:,:), rmat_filt(:,:,:)
        integer :: i, j, k, icut, imask, nlabels
        real    :: weight
        if( .not.allocated(nu_lmask) ) THROW_HARD('nu_lmask not allocated; synthesize_nu_filter_vol_soft')
        nlabels = size(cutoff_finds)
        if( nlabels < 1 ) THROW_HARD('empty NU filter bank; synthesize_nu_filter_vol_soft')
        if( size(bwfilter) /= box ) THROW_HARD('filter size mismatch; synthesize_nu_filter_vol_soft')
        call vol_out%get_rmat_ptr(rmat_out)
        allocate(synth_coord(ldim(1),ldim(2),ldim(3)), source=0.)
        allocate(synth_tmp(ldim(1),ldim(2),ldim(3)),   source=0.)
        allocate(synth_norm(ldim(1),ldim(2),ldim(3)),  source=0.)
        call build_nu_synthesis_coord(synth_coord, synth_tmp, synth_norm)
        call log_nu_soft_synthesis()
        call butterworth_filter(cutoff_finds(1), bwfilter)
        call vol_filt%copy_fast(vol_in_ft)
        call vol_filt%apply_filter(bwfilter)
        call vol_filt%ifft
        call vol_filt%get_rmat_ptr(rmat_filt)
        rmat_out(:ldim(1),:ldim(2),:ldim(3)) = rmat_filt(:ldim(1),:ldim(2),:ldim(3))
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            rmat_out(i,j,k) = 0.
        end do
        !$omp end parallel do
        do icut = 1, nlabels
            if( nu_label_is_aux_replacement(icut) ) cycle
            call butterworth_filter(cutoff_finds(icut), bwfilter)
            call vol_filt%copy_fast(vol_in_ft)
            call vol_filt%apply_filter(bwfilter)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,weight) proc_bind(close)
            do imask = 1, n_nu_mask
                i = nu_mask_vox(1,imask)
                j = nu_mask_vox(2,imask)
                k = nu_mask_vox(3,imask)
                weight = nu_synthesis_weight(synth_coord(i,j,k), icut)
                if( weight > TINY ) rmat_out(i,j,k) = rmat_out(i,j,k) + weight * rmat_filt(i,j,k)
            end do
            !$omp end parallel do
        end do
        deallocate(synth_coord, synth_tmp, synth_norm)
    end subroutine synthesize_nu_filter_vol_soft

    subroutine build_nu_synthesis_coord( synth_coord, work, norm )
        real, intent(inout) :: synth_coord(:,:,:), work(:,:,:), norm(:,:,:)
        integer :: i, j, k, ilabel, imask, nlabels, radius_px
        real    :: coord_min, coord_max
        if( .not.allocated(cutoff_finds)      ) THROW_HARD('cutoff_finds not allocated; build_nu_synthesis_coord')
        if( .not.allocated(candidate_coords)  ) THROW_HARD('candidate_coords not allocated; build_nu_synthesis_coord')
        if( .not.allocated(filtmap)           ) THROW_HARD('filtmap not allocated; build_nu_synthesis_coord')
        if( .not.allocated(nu_lmask)          ) THROW_HARD('nu_lmask not allocated; build_nu_synthesis_coord')
        if( .not.allocated(nu_mask_vox)       ) THROW_HARD('nu_mask_vox not allocated; build_nu_synthesis_coord')
        if( any(shape(synth_coord) /= ldim)   ) THROW_HARD('synth_coord shape mismatch; build_nu_synthesis_coord')
        if( any(shape(work)        /= ldim)   ) THROW_HARD('work shape mismatch; build_nu_synthesis_coord')
        if( any(shape(norm)        /= ldim)   ) THROW_HARD('norm shape mismatch; build_nu_synthesis_coord')
        nlabels = size(cutoff_finds)
        if( size(candidate_coords) < nlabels ) THROW_HARD('candidate_coords size mismatch; build_nu_synthesis_coord')
        radius_px = nu_synthesis_radius_pixels()
        coord_min = minval(candidate_coords(:nlabels))
        coord_max = maxval(candidate_coords(:nlabels))
        synth_coord = 0.
        norm        = 0.
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,ilabel) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            ilabel = max(1, min(nlabels, int(filtmap(i,j,k))))
            synth_coord(i,j,k) = nu_candidate_coord_for_label(ilabel)
            norm(i,j,k) = 1.
        end do
        !$omp end parallel do
        call tent_smooth_3d(synth_coord, work, ldim(1), ldim(2), ldim(3), radius_px)
        call tent_smooth_3d(norm,        work, ldim(1), ldim(2), ldim(3), radius_px)
        !$omp parallel do schedule(static) default(shared) private(imask,i,j,k,ilabel) proc_bind(close)
        do imask = 1, n_nu_mask
            i = nu_mask_vox(1,imask)
            j = nu_mask_vox(2,imask)
            k = nu_mask_vox(3,imask)
            if( norm(i,j,k) > TINY )then
                synth_coord(i,j,k) = synth_coord(i,j,k) / norm(i,j,k)
            else
                ilabel = max(1, min(nlabels, int(filtmap(i,j,k))))
                synth_coord(i,j,k) = nu_candidate_coord_for_label(ilabel)
            endif
            synth_coord(i,j,k) = max(coord_min, min(coord_max, synth_coord(i,j,k)))
        end do
        !$omp end parallel do
    end subroutine build_nu_synthesis_coord

    integer function nu_synthesis_radius_pixels()
        if( smpd <= TINY ) THROW_HARD('invalid smpd; nu_synthesis_radius_pixels')
        nu_synthesis_radius_pixels = max(1, nint(NU_SYNTH_LABEL_SMOOTH_RADIUS_A / smpd))
    end function nu_synthesis_radius_pixels

    real function nu_synthesis_weight( coord, ilabel )
        real,    intent(in) :: coord
        integer, intent(in) :: ilabel
        if( .not.allocated(candidate_coords) ) THROW_HARD('candidate_coords not allocated; nu_synthesis_weight')
        if( ilabel < 1 .or. ilabel > size(candidate_coords) ) &
            &THROW_HARD('label out of range; nu_synthesis_weight')
        nu_synthesis_weight = max(0., 1. - abs(coord - candidate_coords(ilabel)))
    end function nu_synthesis_weight

    subroutine log_nu_soft_synthesis()
        write(logfhandle,'(A,F5.1,A,I0,A)') &
            &'>>> NU filter synthesis: blended bank members over a ', &
            &NU_SYNTH_LABEL_SMOOTH_RADIUS_A, &
            &' A masked tent-smoothed label field (radius ', &
            &nu_synthesis_radius_pixels(), ' px)'
    end subroutine log_nu_soft_synthesis

end submodule simple_nu_filter_apply
