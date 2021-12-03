module simple_tvlam_opt
!$ use omp_lib
!$ use omp_lib_kinds
use simple_opt_spec,    only: opt_spec
use simple_opt_simplex, only: opt_simplex
use simple_image,       only: image
implicit none

public :: tvlam_opt
private
#include "simple_local_flags.inc"

type tvlam_opt
    private
    type(opt_spec)            :: ospec          !< optimizer specification object
    type(opt_simplex)         :: nlopt          !< optimizer object
    type(image),      pointer :: img1 => null() !< img1 to compare with img2
    type(image),      pointer :: img2 => null() !< img2 to compare with img1
    logical, allocatable      :: lmsk(:,:,:)    !< logical mask for distance calc
contains
    procedure :: new
end type tvlam_opt

contains

    subroutine new( self, ldim, smpd, msk )
        class(tvlam_opt), intent(inout) :: self
        integer, intent(in) :: ldim(3)
        real,    intent(in) :: smpd, msk
        type(image) :: mskimg
        real        :: lims(1,2)
        ! logical mask for distance calculation
        call mskimg%disc(ldim, smpd, msk, self%lmsk)
        call mskimg%kill
        ! optimizer specification
        lims(1,1) = 0.5
        lims(1,2) = 5.0
        call self%ospec%specify('simplex', 1, ftol=1e-4, gtol=1e-4, limits=lims, nrestarts=3, maxits=100)
        call self%ospec%set_costfun(tvlam_costfun)
        call self%nlopt%new(self%ospec)
    end subroutine new

    function tvlam_costfun( self, vec, D ) result( cost )
        class(*), intent(inout) :: self
        integer,  intent(in)    :: D
        real,     intent(in)    :: vec(D)
        real                    :: cost
        select type( self )
            class is (tvlam_opt)
                cost = self%img1%sqeuclid(self%img2, mask=self%lmsk)
            class default
                THROW_HARD('error in tvlam_costfun: unknown type')
        end select
    end function tvlam_costfun

end module simple_tvlam_opt
