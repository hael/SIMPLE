!@descr: nonuniform filtering of even/odd volumes
!
! A typical call sequence would be:
!    call setup_dmats(vol_even, vol_odd)
!    call optimize_cutoff_finds()
!    call nonuniform_filter_vol(vol_in, vol_out)
!
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
use simple_tent_smooth, only: tent_smooth_3d
implicit none
#include "simple_local_flags.inc"

public :: write_filtered_vols, setup_dmats, optimize_cutoff_finds, nonuniform_filter_vol
private

real,    parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
integer, parameter   :: WINSZ_TENT = 1

real,    allocatable :: dmats(:,:,:,:)
real,    allocatable :: bwfilters(:,:)
integer, allocatable :: filtmap(:,:,:)
integer, allocatable :: cutoff_finds(:)
integer :: ldim(3), box
real    :: smpd

contains

    subroutine write_filtered_vols( vol_in )
        class(image), intent(in) :: vol_in
        real,    allocatable :: cur_fil(:)
        type(image) :: vol_filt, vol_copy_cmat
        integer :: i, ldim(3), box
        real    :: smpd
        ldim = vol_in%get_ldim()
        smpd = vol_in%get_smpd()
        box  = ldim(1)
        allocate(cur_fil(box), source=0.)
        call vol_copy_cmat%copy(vol_in)
        call vol_copy_cmat%fft
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_ft(.true.)
        call vol_filt%set_wthreads(.true.)
        call vol_copy_cmat%set_wthreads(.true.)
        cutoff_finds = gen_lp_cutoff_finds( box, smpd )
        do i = 1, size(cutoff_finds)
            call vol_filt%copy_fast(vol_copy_cmat)
            call butterworth_filter(vol_filt, cutoff_finds(i), cur_fil)
            call vol_filt%ifft
            call vol_filt%write( string('vol_filt_k_'//int2str(cutoff_finds(i))//'.mrc'))
        end do
        call vol_copy_cmat%kill
        call vol_filt%kill
    end subroutine write_filtered_vols

    function gen_lp_cutoff_finds( box, smpd ) result( cutoff_finds )
        integer, intent(in) :: box
        real,    intent(in) :: smpd
        integer, allocatable :: cutoff_finds(:)
        integer :: i
        if( allocated(cutoff_finds) ) deallocate(cutoff_finds)
        allocate(cutoff_finds(size(lowpass_limits)))
        do i = 1, size(lowpass_limits)
            cutoff_finds(i) = calc_fourier_index(lowpass_limits(i), box, smpd )
        end do
    end function gen_lp_cutoff_finds

    ! subroutine setup_dmats( vol_even, vol_odd, vol_msk )
    !     class(image), intent(in) :: vol_even, vol_odd, vol_msk
    !     logical, allocatable :: lmsk(:,:,:)
    !     type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
    !     type(image) :: vol_even_filt, vol_odd_filt
    !     integer :: ldim(3), box, winsz
    !     real    :: smpd, edge_mean

    !     lmsk = vol_msk%bin2logical()
    !     ldim  = vol_even%get_ldim()
    !     smpd  = vol_even%get_smpd()
    !     box   = ldim(1)
    !     winsz = nint(COSMSKHALFWIDTH)

    !     ! prep cmat copies
    !     call vol_even_copy_cmat%copy(vol_even)
    !     call vol_odd_copy_cmat%copy(vol_odd)
    !     call vol_even_copy_cmat%fft
    !     call vol_odd_copy_cmat%fft
    !     call vol_even_copy_cmat%set_wthreads(.true.)
    !     call vol_odd_copy_cmat%set_wthreads(.true.)
    !     call vol_even_copy_cmat%taper_edges_vol(winsz, edge_mean)
    !     call vol_odd_copy_cmat%taper_edges_vol(winsz, edge_mean)

    !     ! prep filt vols
    !     call vol_even_filt%new(ldim, smpd)
    !     call vol_odd_filt%new(ldim, smpd)
    !     call vol_even_filt%set_ft(.true.)
    !     call vol_odd_filt%set_ft(.true.)
    !     call vol_even_filt%set_wthreads(.true.)
    !     call vol_odd_filt%set_wthreads(.true.)

    ! end subroutine setup_dmats

    subroutine setup_dmats( vol_even, vol_odd )
        class(image), intent(in) :: vol_even, vol_odd
        type(image) :: vol_even_copy_cmat, vol_odd_copy_cmat
        type(image) :: vol_even_filt, vol_odd_filt
        real, allocatable :: dmat_tmp(:,:,:)
        integer :: winsz, i
        real    :: edge_mean, x
        ldim  = vol_even%get_ldim()
        smpd  = vol_even%get_smpd()
        box   = ldim(1)
        winsz = nint(COSMSKHALFWIDTH)

        ! prep cmat copies
        call vol_even_copy_cmat%copy(vol_even)
        call vol_odd_copy_cmat%copy(vol_odd)
        call vol_even_copy_cmat%set_wthreads(.true.)
        call vol_odd_copy_cmat%set_wthreads(.true.)
        call vol_even_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_odd_copy_cmat%taper_edges_vol(winsz, edge_mean)
        call vol_even_copy_cmat%fft
        call vol_odd_copy_cmat%fft

        ! prep filt vols
        call vol_even_filt%new(ldim, smpd)
        call vol_odd_filt%new(ldim, smpd)
        call vol_even_filt%set_ft(.true.)
        call vol_odd_filt%set_ft(.true.)
        call vol_even_filt%set_wthreads(.false.) ! faster with threads off on mac
        call vol_odd_filt%set_wthreads(.false.)  ! faster with threads off on mac

        cutoff_finds = gen_lp_cutoff_finds( box, smpd )
        allocate(bwfilters(box,size(cutoff_finds)),                 source=0.)
        allocate(dmats(size(cutoff_finds),ldim(1),ldim(2),ldim(3)), source=huge(x))
        allocate(dmat_tmp(ldim(1),ldim(2),ldim(3)),                 source=0.)

        do i = 1, size(cutoff_finds)
            call vol_even_filt%copy_fast(vol_even_copy_cmat)
            call vol_odd_filt%copy_fast(vol_odd_copy_cmat)
            call butterworth_filter(vol_even_filt, cutoff_finds(i), bwfilters(:,i))
            call vol_odd_filt%apply_filter(bwfilters(:,i))
            call vol_even_filt%ifft
            call vol_odd_filt%ifft
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmats(i,:,:,:))
            call tent_smooth_3d(dmats(i,:,:,:), dmat_tmp, ldim(1), ldim(2), ldim(3), WINSZ_TENT)
        end do
    end subroutine setup_dmats

    ! this is where the mask goes in
    subroutine optimize_cutoff_finds
        integer :: nx, ny, nz, loc(1), i, j, k
        if( .not.allocated(dmats) ) THROW_HARD('dmats not allocated; run setup_dmats before nonuniform_filter_vol')
        nx = ldim(1)
        ny = ldim(2)
        nz = ldim(3)
        if( allocated(filtmap) ) deallocate(filtmap)
        allocate(filtmap(nx,ny,nz), source=1)
        !$omp parallel do private(loc) collapse(3) schedule(static) default(shared) proc_bind(close)
        do k = 1, nz
            do j = 1, ny
                do i = 1, nx
                    loc = minloc(dmats(:,i,j,k))
                    filtmap(i,j,k) = loc(1)
                end do
            end do
        end do
        !$omp end parallel do
    end subroutine optimize_cutoff_finds

    subroutine nonuniform_filter_vol( vol_in, vol_out )
        class(image), intent(in)    :: vol_in
        class(image), intent(inout) :: vol_out
        type(image) :: vol_copy_cmat, vol_filt
        real(kind=c_float), pointer :: rmat_filt(:,:,:), rmat_out(:,:,:)
        integer :: i, j, k, icut, ldim_in(3)
        real    :: smpd_in
        if( .not.allocated(cutoff_finds) ) THROW_HARD('cutoff_finds not allocated; run setup_dmats before nonuniform_filter_vol')
        if( .not.allocated(filtmap)      ) THROW_HARD('filtmap not allocated; run optimize_cutoff_finds before nonuniform_filter_vol')
        if( .not.allocated(bwfilters)    ) THROW_HARD('bwfilters not allocated; run setup_dmats before nonuniform_filter_vol')
        ldim_in = vol_in%get_ldim()
        smpd_in = vol_in%get_smpd()
        if( any(ldim_in /= ldim)       ) THROW_HARD('Input volume dimensions differ from setup_dmats; nonuniform_filter_vol')
        if( abs(smpd_in - smpd) > TINY ) THROW_HARD('Input volume smpd differs from setup_dmats; nonuniform_filter_vol')
        call vol_copy_cmat%copy(vol_in)
        call vol_copy_cmat%set_wthreads(.false.)
        call vol_copy_cmat%fft
        call vol_filt%new(ldim, smpd)
        call vol_filt%set_ft(.true.)
        call vol_filt%set_wthreads(.false.)
        call vol_out%new(ldim, smpd, wthreads=.false.)
        call vol_out%get_rmat_ptr(rmat_out)
        rmat_out(:ldim(1),:ldim(2),:ldim(3)) = 0.
        do icut = 1, size(cutoff_finds)
            call vol_filt%copy_fast(vol_copy_cmat)
            call vol_filt%apply_filter_serial(bwfilters(:,icut))
            call vol_filt%ifft
            call vol_filt%get_rmat_ptr(rmat_filt)
            !$omp parallel do collapse(3) schedule(static) default(shared) private(i,j,k) proc_bind(close)
            do k = 1, ldim(3)
                do j = 1, ldim(2)
                    do i = 1, ldim(1)
                        if( filtmap(i,j,k) == icut ) rmat_out(i,j,k) = rmat_filt(i,j,k)
                    end do
                end do
            end do
            !$omp end parallel do
        end do
        call vol_copy_cmat%kill
        call vol_filt%kill
    end subroutine nonuniform_filter_vol

end module simple_nu_filter
