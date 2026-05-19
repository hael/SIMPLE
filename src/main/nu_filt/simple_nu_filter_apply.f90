!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

contains

    module subroutine nu_filter_vols( vol_even, vol_odd )
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

    module subroutine nu_postprocess_vol( vol_in, vol_lp, vol_pproc, global_lp, global_bfac )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_lp, vol_pproc
        real,         intent(in)  :: global_lp, global_bfac
        type(image) :: vol_in_ft, vol_sharp_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_lp(:,:,:), rmat_pproc(:,:,:)
        integer :: icut, winsz
        real    :: edge_mean, local_lp
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_postprocess_vol')
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_postprocess_vol')
        if( .not.allocated(srcmap)  ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds before nu_postprocess_vol')
        if( any(vol_in%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; nu_postprocess_vol')
        if( abs(vol_in%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; nu_postprocess_vol')
        if( global_lp <= TINY ) THROW_HARD('Global low-pass limit must be positive; nu_postprocess_vol')
        if( any(nu_lmask .and. srcmap /= 1) )then
            THROW_HARD('NU postprocess transfer functions require a base-bank-only filter map; nu_postprocess_vol')
        endif
        call vol_in_ft%copy(vol_in)
        call vol_in_ft%set_wthreads(.true.)
        if( .not. vol_in_ft%is_ft() )then
            winsz = nint(COSMSKHALFWIDTH)
            call vol_in_ft%taper_edges_vol(winsz, edge_mean)
            call vol_in_ft%fft
        endif
        call vol_sharp_ft%copy(vol_in_ft)
        call vol_sharp_ft%set_wthreads(.true.)
        call vol_sharp_ft%apply_bfac(global_bfac)
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_ft(.true.)
        call vol_filt%set_wthreads(.true.)
        call vol_lp%new(ldim, smpd, wthreads=.false.)
        call vol_pproc%new(ldim, smpd, wthreads=.false.)
        call vol_lp%get_rmat_ptr(rmat_lp)
        call vol_pproc%get_rmat_ptr(rmat_pproc)
        rmat_lp(:ldim(1),:ldim(2),:ldim(3))    = 0.
        rmat_pproc(:ldim(1),:ldim(2),:ldim(3)) = 0.
        call log_nu_postprocess_transfer_bank(global_lp, global_bfac)
        do icut = 1, size(cutoff_finds)
            local_lp = cutoff_find_to_lowpass_limit(icut)
            call vol_filt%copy_fast(vol_in_ft)
            call vol_filt%bp(0., local_lp, NU_POSTPROCESS_HANN_WIDTH)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            call copy_nu_postprocess_voxels(icut, rmat_filt, rmat_lp)
            call vol_filt%copy_fast(vol_sharp_ft)
            call vol_filt%bp(0., local_lp, NU_POSTPROCESS_HANN_WIDTH)
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            call copy_nu_postprocess_voxels(icut, rmat_filt, rmat_pproc)
        end do
        call vol_in_ft%kill
        call vol_sharp_ft%kill
        call vol_filt%kill
    end subroutine nu_postprocess_vol

    subroutine log_nu_postprocess_transfer_bank( global_lp, global_bfac )
        real, intent(in) :: global_lp, global_bfac
        integer :: icut
        real :: local_lp
        write(logfhandle,'(A)') '>>> NU postprocess transfer-function bank'
        write(logfhandle,'(4X,A,F8.3,A,F9.2,A,F6.1)') &
            &'Global FSC LP(A): ', global_lp, '  Global B: ', global_bfac, &
            &'  Hann width(k): ', NU_POSTPROCESS_HANN_WIDTH
        write(logfhandle,'(4X,A)') 'LP limit (A)    Fourier k'
        do icut = 1, size(cutoff_finds)
            local_lp   = cutoff_find_to_lowpass_limit(icut)
            write(logfhandle,'(4X,F10.3,4X,I9)') local_lp, cutoff_finds(icut)
        end do
    end subroutine log_nu_postprocess_transfer_bank

    subroutine copy_nu_postprocess_voxels( icut, rmat_src, rmat_dst )
        integer,            intent(in)    :: icut
        real(kind=c_float), intent(in)    :: rmat_src(:,:,:)
        real(kind=c_float), intent(inout) :: rmat_dst(:,:,:)
        integer :: i, j, k
        !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
        do k = 1, ldim(3)
            do j = 1, ldim(2)
                do i = 1, ldim(1)
                    if( srcmap(i,j,k) == 1 .and. filtmap(i,j,k) == icut )then
                        rmat_dst(i,j,k) = rmat_src(i,j,k)
                    endif
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine copy_nu_postprocess_voxels

end submodule simple_nu_filter_apply
