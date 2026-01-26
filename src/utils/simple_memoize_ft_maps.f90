!@descr: Light-weight module to memoize logical to physical address/spatial frequency mapping and avoid re-computing them repeatedly
module simple_memoize_ft_maps
use simple_core_module_api
implicit none

public :: memoize_ft_maps, forget_ft_maps
public :: ft_map_lims, ft_map_lims_nr
public :: ft_map_phys_addrh, ft_map_phys_addrk, ft_map_spafreqsq, ft_map_astigang
private
#include "simple_local_flags.inc"

integer              :: ft_map_lims(3,2)    = 0     ! Redundant Fourier limits
integer              :: ft_map_lims_nr(3,2) = 0     ! Non-redudant Fourier limits
integer, allocatable :: ft_map_phys_addrh(:,:)      ! Physical address along h-axis
integer, allocatable :: ft_map_phys_addrk(:,:)      ! Physical address along k-axis
real,    allocatable :: ft_map_spafreqsq(:,:)       ! Spatial frequency squared
real,    allocatable :: ft_map_astigang(:,:)        ! Angle of astigmatism

contains

    subroutine memoize_ft_maps( ldim, mode )
        use simple_ftiter, only: ftiter
        integer,           intent(in) :: ldim(2)
        integer, optional, intent(in) :: mode   ! 2:non-redundant | 3: redudant adresses
        type(ftiter) :: ftiterator
        real         :: spafreqh, spafreqk, spafreqsqk
        integer      :: lims(3,2), phys(2), h,k, imode
        if( OMP_IN_PARALLEL() )     THROW_HARD('No memoization inside OpenMP regions')
        ! Friedel symmetry
        if( present(mode) )then
            if( (mode<2) .or. (mode>3) )THROW_HARD('Incorrect mode')
            imode = mode
            ! when mode=3 the Fourier component corresponding
            ! to the address must be complex conjugated when h<0
        else
            imode = 2
        endif
        ! init & alloc
        call forget_ft_maps
        ftiterator     = ftiter([ldim(1), ldim(2), 1], 1.0)
        ft_map_lims_nr = ftiterator%loop_lims(3)
        ft_map_lims    = ftiterator%loop_lims(2)
        if( imode==2 )then
            lims = ft_map_lims
        else
            lims = ft_map_lims_nr
        endif
        allocate(ft_map_phys_addrh(lims(1,1):lims(1,2), lims(2,1):lims(2,2)),&
                &ft_map_phys_addrk(lims(1,1):lims(1,2), lims(2,1):lims(2,2)),&
                &ft_map_spafreqsq(lims(1,1):lims(1,2), lims(2,1):lims(2,2)),&
                &ft_map_astigang(lims(1,1):lims(1,2), lims(2,1):lims(2,2)))
        ! fill matrix
        do k = lims(2,1),lims(2,2)
            spafreqk   = real(k) / real(ldim(2))
            spafreqsqk = spafreqk * spafreqk
            do h = lims(1,1),lims(1,2)
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

end module simple_memoize_ft_maps
