!@descr: projector-to-polarft adapters kept outside simple_projector to avoid module cycles
module simple_projector_pft
use simple_core_module_api
use simple_projector,    only: projector
use simple_polarft_calc, only: polarft_calc
implicit none

public :: fproject_polar, fproject_polar_oversamp
private
#include "simple_local_flags.inc"

contains

    subroutine fproject_polar( self, iref, e, pftc, iseven )
        class(projector),    intent(inout) :: self
        integer,             intent(in)    :: iref
        class(ori),          intent(in)    :: e
        class(polarft_calc), intent(inout) :: pftc
        logical,             intent(in)    :: iseven
        integer :: pdim(3), irot, k
        real    :: loc(3), e_rotmat(3,3), hk(2)
        pdim     = pftc%get_pdim_srch()
        e_rotmat = e%get_mat()
        do irot = 1, pdim(1)
            do k = pdim(2), pdim(3)
                hk  = pftc%get_coord(irot,k)
                loc = matmul([hk(1), hk(2), 0.0], e_rotmat)
                call pftc%set_ref_fcomp(iref, irot, k, self%interp_fcomp(loc), iseven)
            enddo
        enddo
    end subroutine fproject_polar

    subroutine fproject_polar_oversamp( self, iref, e, pftc, iseven )
        class(projector),    intent(inout) :: self
        integer,             intent(in)    :: iref
        class(ori),          intent(in)    :: e
        class(polarft_calc), intent(inout) :: pftc
        logical,             intent(in)    :: iseven
        integer :: pdim(3), irot, k
        real    :: loc(3), e_rotmat(3,3), hk(2)
        pdim     = pftc%get_pdim_srch()
        e_rotmat = e%get_mat()
        do irot = 1, pdim(1)
            do k = pdim(2), pdim(3)
                hk  = pftc%get_coord(irot,k)
                loc = matmul([hk(1), hk(2), 0.0], e_rotmat)
                call pftc%set_ref_fcomp(iref, irot, k, self%interp_fcomp_oversamp(loc), iseven)
            enddo
        enddo
    end subroutine fproject_polar_oversamp

end module simple_projector_pft
