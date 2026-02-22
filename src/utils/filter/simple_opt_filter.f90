!@descr: optimization(search)-based filtering
module simple_opt_filter
use simple_core_module_api
use simple_image, only: image
use simple_butterworth
implicit none
#include "simple_local_flags.inc"

public :: estimate_lplim, estimate_lplims2D
private

interface estimate_lplim
    module procedure estimate_lplim_1, estimate_lplim_2
end interface estimate_lplim


contains

    subroutine estimate_lplim_1( odd, even, mskimg, kfromto, best_ind, odd_filt_out )
        class(image),           intent(inout) :: odd, even, mskimg
        integer,                intent(in)    :: kfromto(2)
        integer,                intent(out)   :: best_ind
        class(image), optional, intent(out)   :: odd_filt_out
        real,    allocatable :: cur_fil(:)
        logical, allocatable :: lmask(:,:,:)
        type(image) :: odd_copy_cmat, odd_filt
        integer     :: box, ldim(3), cutoff_find
        real        :: smpd, cur_cost, best_cost
        if( kfromto(1) == kfromto(2) )then  
            best_ind = kfromto(1)
            return
        endif
        if( odd%is_ft() .or. even%is_ft() ) THROW_HARD('Input even & odd has to be in real-space representation')
        if( mskimg%is_ft() ) THROW_HARD('Input mskimg has to be in real-space representation')
        ldim  = odd%get_ldim()
        box   = ldim(1)
        smpd  = odd%get_smpd()
        lmask = mskimg%bin2logical()
        call odd_copy_cmat%copy(odd)
        call odd_copy_cmat%fft
        call odd_filt%new(ldim, smpd)
        call odd_filt%set_ft(.true.)
        allocate(cur_fil(box), source=0.)
        best_ind  = kfromto(1)
        best_cost = huge(best_cost)
        if( ldim(3) > 1 )then ! 3D, so parallelize
            call odd_filt%set_wthreads(.true.)
            call even%set_wthreads(.true.)
        endif
        do cutoff_find = kfromto(1), kfromto(2)
            ! apply BW kernel to odd
            call odd_filt%copy_fast(odd_copy_cmat)
            call butterworth_filter(odd_filt, cutoff_find, cur_fil)
            call odd_filt%ifft
            ! calculate squared Euclidean distance
            cur_cost = odd_filt%sqeuclid(even, lmask)
            if( cur_cost <= best_cost )then
                best_cost = cur_cost
                best_ind  = cutoff_find
            endif
        enddo
        if( present(odd_filt_out) )then
            call odd_filt_out%new(ldim, smpd)
            call odd_filt_out%set_ft(.true.)
            call odd_filt_out%copy_fast(odd_copy_cmat)
            call butterworth_filter(odd_filt_out, best_ind, cur_fil)
            call odd_filt_out%ifft
        endif
        call odd_copy_cmat%kill
        call odd_filt%kill
    end subroutine estimate_lplim_1

    subroutine estimate_lplim_2( odd, even, mskimg, lprange, lpopt, odd_filt_out )
        class(image),           intent(inout) :: odd, even, mskimg
        real,                   intent(in)    :: lprange(2)
        real,                   intent(out)   :: lpopt
        class(image), optional, intent(out)   :: odd_filt_out
        integer :: ldim(3), box, kfromto(2), best_ind
        real    :: smpd
        ldim       = odd%get_ldim()
        box        = ldim(1)
        smpd       = odd%get_smpd()
        kfromto(1) = calc_fourier_index(lprange(1), box, smpd)
        kfromto(2) = calc_fourier_index(lprange(2), box, smpd)
        call estimate_lplim_1( odd, even, mskimg, kfromto, best_ind, odd_filt_out )
        lpopt = calc_lowpass_lim(best_ind, box, smpd)
    end subroutine estimate_lplim_2

    subroutine estimate_lplims2D( odd, even, mskrad_px, lprange, lpsopt, odd_filt_out )
        class(image),                       intent(inout) :: odd(:), even(:)
        real,                               intent(in)    :: mskrad_px, lprange(2)
        real,                  allocatable, intent(inout) :: lpsopt(:)
        type(image), optional, allocatable, intent(inout) :: odd_filt_out(:)
        type(image), allocatable :: masks(:)
        real    :: smpd
        integer :: iptcl, box, ldim(3), nptcls, npix
        write(logfhandle,'(A)') '>>> 2D UNIFORM FILTERING FOR LP ESTIMATION'
        ! init
        ldim    = odd(1)%get_ldim()
        ldim(3) = 1 ! because we operate on stacks
        box     = ldim(1)
        smpd    = odd(1)%get_smpd()
        nptcls  = size(odd)
        if( allocated(lpsopt) ) deallocate(lpsopt)
        allocate(masks(nptcls), lpsopt(nptcls))
        lpsopt = 0.
        call masks(1)%disc(ldim, smpd, mskrad_px, npix)
        do iptcl = 2, nptcls
            call masks(iptcl)%copy(masks(1))
        end do
        if( present(odd_filt_out) )then
            do iptcl = 1, size(odd_filt_out)
                call odd_filt_out(iptcl)%kill
            end do
            deallocate(odd_filt_out)
            allocate(odd_filt_out(nptcls))
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                call estimate_lplim(odd(iptcl), even(iptcl), masks(iptcl), lprange, lpsopt(iptcl), odd_filt_out(iptcl))
            enddo
            !$omp end parallel do
        else
            !$omp parallel do default(shared) private(iptcl) schedule(static) proc_bind(close)
            do iptcl = 1, nptcls
                call estimate_lplim(odd(iptcl), even(iptcl), masks(iptcl), lprange, lpsopt(iptcl))
            enddo
            !$omp end parallel do
        endif
        do iptcl = 1, nptcls
            call masks(iptcl)%kill
        end do
        deallocate(masks)
    end subroutine estimate_lplims2D

end module simple_opt_filter
