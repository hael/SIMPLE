!@descr: nonuniform filtering of even/odd volumes
module simple_nu_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
implicit none
#include "simple_local_flags.inc"

public :: write_filtered_vols, setup_dmats
private

real, parameter   :: lowpass_limits(8) = [20.,15.,12.,10.,8.,6.,5.,4.]
real, allocatable :: dmats(:,:,:,:)

contains

    subroutine write_filtered_vols( vol_in )
        class(image), intent(in) :: vol_in
        integer, allocatable :: cutoff_finds(:)
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
        integer, allocatable :: cutoff_finds(:)
        real,    allocatable :: cur_fil(:)
        integer :: ldim(3), box, winsz, i
        real    :: smpd, edge_mean, x

        ldim  = vol_even%get_ldim()
        smpd  = vol_even%get_smpd()
        box   = ldim(1)
        winsz = nint(COSMSKHALFWIDTH)
        allocate(cur_fil(box), source=0.)

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
        call vol_even_filt%set_wthreads(.true.)
        call vol_odd_filt%set_wthreads(.true.)

        cutoff_finds = gen_lp_cutoff_finds( box, smpd )

        allocate(dmats(ldim(1),ldim(2),ldim(3),size(cutoff_finds)), source=huge(x))

        do i = 1, size(cutoff_finds)
            call vol_even_filt%copy_fast(vol_even_copy_cmat)
            call vol_odd_filt%copy_fast(vol_odd_copy_cmat)
            call butterworth_filter(vol_even_filt, cutoff_finds(i), cur_fil)
            call butterworth_filter(vol_odd_filt,  cutoff_finds(i), cur_fil)
            call vol_even_filt%ifft
            call vol_odd_filt%ifft
            call vol_even%nu_objective(vol_even_filt, vol_odd, vol_odd_filt, dmats(:,:,:,i))
        end do
    end subroutine setup_dmats


end module simple_nu_filter
