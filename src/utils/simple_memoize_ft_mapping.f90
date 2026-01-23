! Light-weight module to memoize logical to physical address/spatial
! frequency mapping and avoid re-computing them repeatedly
module simple_memoize_ft_mapping
use simple_core_module_api
implicit none

public :: memoize_ft_map, forget_ft_map
public :: get_phys_addr, get_spaFreqSq, get_astigang
private
#include "simple_local_flags.inc"

type ftmapping
    integer :: phys(2)      ! Physical h/k addresses
    real    :: spaFreqSq    ! Spatial frequency
    real    :: ang          ! angle on unit circle (astigmatism, radians)
end type ftmapping

type(ftmapping), allocatable :: ft_map(:,:)

contains

    subroutine memoize_ft_map( ldim, mode )
        use simple_ftiter, only: ftiter
        integer,           intent(in) :: ldim(2)
        integer, optional, intent(in) :: mode   ! 2:non-redundant | 3: redudant adresses
        type(ftiter) :: ftiterator
        real         :: spafreqh, spafreqk, spafreqsqk
        integer      :: lims(3,2), h,k, imode
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
        ftiterator = ftiter([ldim(1), ldim(2), 1], 1.0)
        lims       = ftiterator%loop_lims(imode)
        if( allocated(ft_map) ) deallocate(ft_map)
        allocate(ft_map(lims(1,1):lims(1,2), lims(2,1):lims(2,2)))
        ! fill matrix
        do k = lims(2,1),lims(2,2)
            spafreqk   = real(k) / real(ldim(2))
            spafreqsqk = spafreqk * spafreqk
            do h = lims(1,1),lims(1,2)
                spafreqh = real(h) / real(ldim(1))
                ft_map(h,k)%phys      = ftiterator%comp_addr_phys(h,k)
                ft_map(h,k)%spafreqsq = spafreqh*spafreqh + spafreqsqk
                ft_map(h,k)%ang       = atan2(spafreqk,spafreqh)
            end do
        end do
    end subroutine memoize_ft_map

    subroutine forget_ft_map
        if( allocated(ft_map) ) deallocate(ft_map)
    end subroutine forget_ft_map

    pure function get_phys_addr(h, k) result (phys)
        integer, intent(in) :: h,k
        integer :: phys(2)
        phys = ft_map(h,k)%phys
    end function get_phys_addr

    pure real function get_spafreqsq(h, k)
        integer, intent(in) :: h,k
        get_spafreqsq = ft_map(h,k)%spaFreqSq
    end function get_spafreqsq

    pure real function get_astigang(h, k)
        integer, intent(in) :: h,k
        get_astigang = ft_map(h,k)%ang
    end function get_astigang

end module simple_memoize_ft_mapping