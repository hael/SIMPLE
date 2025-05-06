module simple_strategy2D_utils
include 'simple_lib.f08'
use simple_image,  only: image
use simple_masker, only: density_outside_mask
implicit none

public :: flag_non_junk_cavgs
private
#include "simple_local_flags.inc"

contains

    subroutine flag_non_junk_cavgs( cavgs, os_cls2D, lp_bin, msk, l_non_junk )
        class(image),         intent(inout) :: cavgs(:)
        class(oris),          intent(in)    :: os_cls2D
        real,                 intent(in)    :: lp_bin, msk
        logical, allocatable, intent(inout) :: l_non_junk(:)
        real,        parameter   :: DYNRANGE_THRES = 1e-6
        real,        parameter   :: HP_SPEC        = 20.
        real,        parameter   :: LP_SPEC        = 6.
        integer,     parameter   :: MINPOP         = 20
        type(image), allocatable :: cavg_threads(:)
        real,        allocatable :: pspec(:)
        integer :: ncls, icls, ldim(3), kfromto(2), ithr
        real    :: dynrange, smpd
        logical :: l_dens_outside
        ncls = size(cavgs)
        if( os_cls2D%get_noris() /= ncls ) THROW_HARD('# cavgs /= # entries in os_cls2D')
        ldim = cavgs(1)%get_ldim()
        smpd = cavgs(1)%get_smpd()
        if( allocated(l_non_junk) ) deallocate(l_non_junk)
        allocate(l_non_junk(ncls), source=.false.)
        kfromto(1) = calc_fourier_index(HP_SPEC, ldim(1), smpd)
        kfromto(2) = calc_fourier_index(LP_SPEC, ldim(1), smpd)
        allocate(cavg_threads(nthr_glob))
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%new(ldim, smpd)
        end do
        !$omp parallel do default(shared) private(icls,ithr,pspec,dynrange,l_dens_outside) proc_bind(close) schedule(static)
        do icls = 1, ncls
            ithr = omp_get_thread_num() + 1
            call cavg_threads(ithr)%copy(cavgs(icls))
            call cavg_threads(ithr)%norm
            l_dens_outside = density_outside_mask(cavg_threads(ithr), lp_bin, msk)
            call cavg_threads(ithr)%mask(msk, 'soft')
            call cavg_threads(ithr)%spectrum('sqrt', pspec)
            dynrange = pspec(kfromto(1)) - pspec(kfromto(2))
            if( dynrange > DYNRANGE_THRES .and. os_cls2D%get_int(icls, 'pop') >= MINPOP )then
                if( .not. l_dens_outside ) l_non_junk(icls) = .true.
            endif
        enddo
        !$omp end parallel do
        do ithr = 1, nthr_glob
            call cavg_threads(ithr)%kill
        end do
        deallocate(cavg_threads)
    end subroutine flag_non_junk_cavgs

end module simple_strategy2D_utils
