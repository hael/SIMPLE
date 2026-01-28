!@descr: Light-weight module to memoize logical to physical address/spatial frequency mapping and avoid re-computing them repeatedly
module simple_memoize_ft_maps
use simple_core_module_api
use simple_ftiter, only: ftiter
implicit none

public :: memoize_ft_maps, forget_ft_maps
public :: ft_map_lims, ft_map_lims_nr
public :: ft_map_phys_addrh, ft_map_phys_addrk, ft_map_spafreqsq, ft_map_astigang
public :: ft_map_get_lfny, ft_map_get_farray_shape
private
#include "simple_local_flags.inc"

type(ftiter), private :: ftiterator
integer               :: ft_map_lims(3,2)    = 0     ! Redundant Fourier limits
integer               :: ft_map_lims_nr(3,2) = 0     ! Non-redudant Fourier limits
integer,  allocatable :: ft_map_phys_addrh(:,:)      ! Physical address along h-axis
integer,  allocatable :: ft_map_phys_addrk(:,:)      ! Physical address along k-axis
real,     allocatable :: ft_map_spafreqsq(:,:)       ! Spatial frequency squared
real,     allocatable :: ft_map_astigang(:,:)        ! Angle of astigmatism

contains

    subroutine memoize_ft_maps( ldim, smpd )
        integer, intent(in) :: ldim(2)
        real,    intent(in) :: smpd
        real         :: spafreqh, spafreqk, spafreqsqk
        integer      :: phys(2), h,k
        if( OMP_IN_PARALLEL() ) THROW_HARD('No memoization inside OpenMP regions')
        ! Friedel symmetry
        ! init & alloc
        call forget_ft_maps
        ftiterator     = ftiter([ldim(1), ldim(2), 1], smpd)
        ft_map_lims_nr = ftiterator%loop_lims(3)
        ft_map_lims    = ftiterator%loop_lims(2)
        allocate(ft_map_phys_addrh(ft_map_lims_nr(1,1):ft_map_lims_nr(1,2), ft_map_lims_nr(2,1):ft_map_lims_nr(2,2)),&
                &ft_map_phys_addrk(ft_map_lims_nr(1,1):ft_map_lims_nr(1,2), ft_map_lims_nr(2,1):ft_map_lims_nr(2,2)),&
                &ft_map_spafreqsq( ft_map_lims_nr(1,1):ft_map_lims_nr(1,2), ft_map_lims_nr(2,1):ft_map_lims_nr(2,2)),&
                &ft_map_astigang(  ft_map_lims_nr(1,1):ft_map_lims_nr(1,2), ft_map_lims_nr(2,1):ft_map_lims_nr(2,2)))
        ! fill matrix
        do k = ft_map_lims_nr(2,1),ft_map_lims_nr(2,2)
            spafreqk   = real(k) / real(ldim(2))
            spafreqsqk = spafreqk * spafreqk
            do h = ft_map_lims_nr(1,1),ft_map_lims_nr(1,2)
                phys     = ftiterator%comp_addr_phys(h,k)
                spafreqh = real(h) / real(ldim(1))
                ft_map_phys_addrh(h,k) = phys(1)
                ft_map_phys_addrk(h,k) = phys(2)
                ft_map_spafreqsq(h,k)  = spafreqh*spafreqh + spafreqsqk
                ft_map_astigang(h,k)   = atan2(real(k), real(h))
            end do
        end do
    end subroutine memoize_ft_maps

    subroutine forget_ft_maps
        if( allocated(ft_map_phys_addrh) )then
            ft_map_lims    = 0
            ft_map_lims_nr = 0
            deallocate(ft_map_phys_addrh,ft_map_phys_addrk,ft_map_spafreqsq,ft_map_astigang)
        endif
    end subroutine forget_ft_maps

    pure integer function ft_map_get_lfny()
        ft_map_get_lfny = ftiterator%get_lfny(1)
    end function ft_map_get_lfny

    function ft_map_get_farray_shape( )result( cshape )
        integer :: cshape(3)
        if( allocated(ft_map_phys_addrh) )then
            cshape = [fdim(size(ft_map_phys_addrh,1)), size(ft_map_phys_addrh,2), 1]
        else
            cshape = -1
        endif
    end function ft_map_get_farray_shape

end module simple_memoize_ft_maps
