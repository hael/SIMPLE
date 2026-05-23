!@descr: simple nu filter apply implementation for volume-domain nonuniform filtering
submodule (simple_nu_filter) simple_nu_filter_apply
implicit none
#include "simple_local_flags.inc"

! NU postprocess B-factor model tuning. These constants live in the applying
! submodule so local tuning edits force recompilation of the code that computes
! the transfer-function bank.
! - Bins at or better than the global FSC resolution use the global B factor.
!   This preserves classical sharpening wherever the global FSC supports it.
! - DAMPING_BFAC_REF fixes the positive-B endpoint of the low-resolution branch
!   independently of the fitted global B factor. Increase the magnitude if
!   low-resolution regions remain too sharp/noisy; decrease it if those regions
!   become too blurred.
! - ALPHA scales that positive-B endpoint. With the defaults, maximally damped
!   bins approach the same effective B factor regardless of the fitted global
!   sharpening B.
! - SIGMOID_MID is the resolution where damping turns on most rapidly.
!   Increase it if mid-resolution density is damped too early; decrease it if
!   damping should begin closer to the global FSC limit.
! - SIGMOID_WIDTH controls transition sharpness. Increase it for a gentler
!   change across local-resolution bins; decrease it for a sharper transition.
! - CLASSICAL_LP_WINDOW_A protects the immediate FSC-resolution neighborhood:
!   bins within this distance of the global FSC low-pass reuse the exact
!   classical postprocess transfer map supplied by the caller.
real, parameter :: NU_POSTPROCESS_BFAC_ALPHA           = 0.75
real, parameter :: NU_POSTPROCESS_DAMPING_BFAC_REF     = -133.333
real, parameter :: NU_POSTPROCESS_BFAC_SIGMOID_MID     = 6.
real, parameter :: NU_POSTPROCESS_BFAC_SIGMOID_WIDTH   = 1.5
real, parameter :: NU_POSTPROCESS_ANTIALIAS_HANN_WIDTH = 4.
real, parameter :: NU_POSTPROCESS_CLASSICAL_LP_WINDOW_A = 0.5

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
        call release_nu_filter_unary_storage
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

    module subroutine nu_postprocess_vol( vol_in, vol_lp, vol_pproc, global_lp, global_bfac, aux_vols, &
            &global_filter )
        class(image), intent(in)  :: vol_in
        class(image), intent(out) :: vol_lp, vol_pproc
        real,         intent(in)  :: global_lp, global_bfac
        type(image), optional, intent(in) :: aux_vols(:)
        real, optional, intent(in) :: global_filter(:)
        type(image) :: vol_in_ft, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_pproc(:,:,:)
        real(kind=c_float), pointer :: rmat_aux_classical(:,:,:) => null()
        integer :: icut, winsz
        real    :: edge_mean, local_lp, local_bfac, antialias_lp
        logical :: l_have_classical_aux, l_classical_window, l_have_global_filter
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_nu_dmats before nu_postprocess_vol')
        if( .not.allocated(filtmap) ) THROW_HARD('filtmap not allocated; run optimize_nu_cutoff_finds before nu_postprocess_vol')
        if( .not.allocated(srcmap)  ) THROW_HARD('srcmap not allocated; run optimize_nu_cutoff_finds before nu_postprocess_vol')
        if( any(vol_in%get_ldim() /= ldim)       ) THROW_HARD('Input volume dimensions differ; nu_postprocess_vol')
        if( abs(vol_in%get_smpd() - smpd) > TINY ) THROW_HARD('Input volume smpd differs; nu_postprocess_vol')
        if( global_lp <= TINY ) THROW_HARD('Global low-pass limit must be positive; nu_postprocess_vol')
        l_have_global_filter = present(global_filter)
        call validate_nu_postprocess_aux_vols(aux_vols)
        l_have_classical_aux = present(aux_vols)
        if( l_have_classical_aux ) l_have_classical_aux = size(aux_vols) >= 1
        if( l_have_classical_aux ) call aux_vols(1)%get_rmat_ptr(rmat_aux_classical)
        call release_nu_filter_unary_storage
        call vol_in_ft%copy(vol_in)
        call vol_in_ft%set_wthreads(.true.)
        if( .not. vol_in_ft%is_ft() )then
            winsz = nint(COSMSKHALFWIDTH)
            call vol_in_ft%taper_edges_vol(winsz, edge_mean)
            call vol_in_ft%fft
        endif
        antialias_lp = global_lp
        call vol_lp%copy(vol_in_ft)
        call vol_lp%set_wthreads(.true.)
        call apply_nu_postprocess_global_filter(vol_lp)
        call vol_lp%ifft
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_ft(.true.)
        call vol_filt%set_wthreads(.true.)
        call vol_pproc%new(ldim, smpd, wthreads=.false.)
        call vol_pproc%get_rmat_ptr(rmat_pproc)
        rmat_pproc(:ldim(1),:ldim(2),:ldim(3)) = 0.
        call log_nu_postprocess_transfer_bank(global_lp, global_bfac, antialias_lp, &
            &l_have_classical_aux, l_have_global_filter)
        do icut = 1, size(cutoff_finds)
            local_lp = cutoff_find_to_lowpass_limit(icut)
            l_classical_window = l_have_classical_aux .and. &
                &abs(local_lp - global_lp) <= NU_POSTPROCESS_CLASSICAL_LP_WINDOW_A
            if( l_classical_window )then
                call copy_nu_postprocess_voxels(icut, rmat_aux_classical, rmat_pproc)
            else
                local_bfac = nu_postprocess_resolution_bfac(local_lp, global_lp, global_bfac)
                call vol_filt%copy_fast(vol_in_ft)
                call vol_filt%apply_bfac(local_bfac)
                call apply_nu_postprocess_global_filter(vol_filt)
                call vol_filt%ifft
                call vol_filt%get_rmat_ptr(rmat_filt)
                call copy_nu_postprocess_voxels(icut, rmat_filt, rmat_pproc)
            endif
        end do
        if( present(aux_vols) ) call copy_nu_postprocess_aux_voxels(aux_vols, rmat_pproc)
        call vol_in_ft%kill
        call vol_filt%kill
    contains
        subroutine apply_nu_postprocess_global_filter( vol )
            type(image), intent(inout) :: vol
            if( l_have_global_filter )then
                call vol%apply_filter(global_filter)
            else
                call vol%bp(0., antialias_lp, NU_POSTPROCESS_ANTIALIAS_HANN_WIDTH)
            endif
        end subroutine apply_nu_postprocess_global_filter
    end subroutine nu_postprocess_vol

    subroutine validate_nu_postprocess_aux_vols( aux_vols )
        type(image), optional, intent(in) :: aux_vols(:)
        integer :: iaux, n_aux_needed, n_aux_to_check
        n_aux_needed = max(0, maxval(srcmap) - 1)
        if( n_aux_needed == 0 .and. .not. present(aux_vols) ) return
        if( .not. present(aux_vols) )then
            THROW_HARD('NU postprocess selected aux source but no aux output volume; nu_postprocess_vol')
        endif
        if( size(aux_vols) < n_aux_needed )then
            THROW_HARD('NU postprocess aux output volume count too small; nu_postprocess_vol')
        endif
        n_aux_to_check = size(aux_vols)
        do iaux = 1, n_aux_to_check
            if( any(aux_vols(iaux)%get_ldim() /= ldim) )then
                THROW_HARD('NU postprocess auxiliary output volume dimensions differ; nu_postprocess_vol')
            endif
            if( abs(aux_vols(iaux)%get_smpd() - smpd) > TINY )then
                THROW_HARD('NU postprocess auxiliary output volume smpd differs; nu_postprocess_vol')
            endif
            if( aux_vols(iaux)%is_ft() )then
                THROW_HARD('NU postprocess auxiliary output volumes must be real-space maps; nu_postprocess_vol')
            endif
        end do
    end subroutine validate_nu_postprocess_aux_vols

    subroutine log_nu_postprocess_transfer_bank( global_lp, global_bfac, antialias_lp, &
            &l_have_classical_aux, l_have_global_filter )
        real, intent(in) :: global_lp, global_bfac, antialias_lp
        logical, intent(in) :: l_have_classical_aux, l_have_global_filter
        integer :: icut
        real :: local_lp, damp_frac, local_bfac, damping_plateau
        character(len=9) :: mode
        damping_plateau = nu_postprocess_bfac_damping_plateau()
        write(logfhandle,'(A)') '>>> NU postprocess transfer-function bank'
        write(logfhandle,'(4X,A,F8.3,A,F9.2)') &
            &'Damping reference LP(A): ', global_lp, '  Global B: ', global_bfac
        write(logfhandle,'(4X,A,F9.2,A,F9.2)') &
            &'Fixed damping B ref: ', NU_POSTPROCESS_DAMPING_BFAC_REF, &
            &'  B floor: ', damping_plateau
        write(logfhandle,'(4X,A,F5.2,A,F6.2,A,F5.2)') &
            &'B-factor model alpha/mid/width: ', NU_POSTPROCESS_BFAC_ALPHA, '/', &
            &NU_POSTPROCESS_BFAC_SIGMOID_MID, '/', NU_POSTPROCESS_BFAC_SIGMOID_WIDTH
        if( l_have_global_filter )then
            write(logfhandle,'(4X,A)') 'Antialias filter: global FSC-derived filter'
        else
            write(logfhandle,'(4X,A,F8.3,A,F5.1,A)') &
                &'Antialias Hann LP(A): ', antialias_lp, '  Width: ', &
                &NU_POSTPROCESS_ANTIALIAS_HANN_WIDTH, ' Fourier pixels'
        endif
        if( l_have_classical_aux )then
            write(logfhandle,'(4X,A,F5.2,A)') &
                &'Classical transfer window: +/- ', NU_POSTPROCESS_CLASSICAL_LP_WINDOW_A, &
                &' A around global FSC LP'
        endif
        write(logfhandle,'(4X,A)') 'LP limit (A)    Fourier k    Damp frac      B_eff       Mode'
        do icut = 1, size(cutoff_finds)
            local_lp   = cutoff_find_to_lowpass_limit(icut)
            damp_frac  = nu_postprocess_damping_frac(local_lp, global_lp)
            local_bfac = nu_postprocess_resolution_bfac(local_lp, global_lp, global_bfac)
            if( l_have_classical_aux .and. &
                &abs(local_lp - global_lp) <= NU_POSTPROCESS_CLASSICAL_LP_WINDOW_A )then
                mode = 'classical'
            else
                mode = 'NU'
            endif
            write(logfhandle,'(4X,F10.3,4X,I9,4X,F9.3,4X,F9.2,4X,A)') &
                &local_lp, cutoff_finds(icut), damp_frac, local_bfac, trim(mode)
        end do
    end subroutine log_nu_postprocess_transfer_bank

    real function nu_postprocess_resolution_bfac( local_lp, global_lp, global_bfac ) result( local_bfac )
        real, intent(in) :: local_lp, global_lp, global_bfac
        real :: damp_frac
        ! The global B factor still defines classical sharpening at and better
        ! than the global FSC resolution. Worse local-resolution bins
        ! interpolate toward a fixed positive-B endpoint derived from the
        ! reference damping B factor, so weakly sharpened maps do not also get a
        ! weak low-resolution penalty.
        damp_frac  = nu_postprocess_damping_frac(local_lp, global_lp)
        local_bfac = global_bfac + damp_frac * (nu_postprocess_bfac_damping_plateau() - global_bfac)
    end function nu_postprocess_resolution_bfac

    real function nu_postprocess_bfac_damping_plateau() result( damping_plateau )
        damping_plateau = NU_POSTPROCESS_BFAC_ALPHA * abs(NU_POSTPROCESS_DAMPING_BFAC_REF)
    end function nu_postprocess_bfac_damping_plateau

    real function nu_postprocess_damping_frac( local_lp, global_lp ) result( frac )
        real, intent(in) :: local_lp, global_lp
        real :: sig, sig_global, width
        if( local_lp <= global_lp )then
            frac = 0.
            return
        endif
        width      = max(TINY, NU_POSTPROCESS_BFAC_SIGMOID_WIDTH)
        sig        = nu_postprocess_sigmoid(local_lp,  width)
        sig_global = nu_postprocess_sigmoid(global_lp, width)
        if( sig_global >= 1. - TINY )then
            frac = 0.
        else
            frac = (sig - sig_global) / (1. - sig_global)
        endif
        frac = min(1., max(0., frac))
    end function nu_postprocess_damping_frac

    real function nu_postprocess_sigmoid( lp_angstrom, width ) result( sig )
        real, intent(in) :: lp_angstrom, width
        sig = 1. / (1. + exp(-(lp_angstrom - NU_POSTPROCESS_BFAC_SIGMOID_MID) / width))
    end function nu_postprocess_sigmoid

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

    subroutine copy_nu_postprocess_aux_voxels( aux_vols, rmat_dst )
        type(image),        intent(in)    :: aux_vols(:)
        real(kind=c_float), intent(inout) :: rmat_dst(:,:,:)
        real(kind=c_float), pointer :: rmat_aux(:,:,:)
        integer :: iaux, i, j, k
        do iaux = 1, size(aux_vols)
            call aux_vols(iaux)%get_rmat_ptr(rmat_aux)
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( srcmap(i,j,k) == iaux + 1 )then
                            rmat_dst(i,j,k) = rmat_aux(i,j,k)
                        endif
                    end do
                end do
            end do
            !$omp end parallel do
        end do
    end subroutine copy_nu_postprocess_aux_voxels

end submodule simple_nu_filter_apply
