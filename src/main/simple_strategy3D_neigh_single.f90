! concrete strategy3D: neigh single-state refinement
module simple_strategy3D_neigh_single
include 'simple_lib.f08'
use simple_strategy3D_alloc        ! use all in there
use simple_strategy3D_utils        ! use all in there
use simple_strategy3D_neigh_multi, only: strategy3D_neigh_multi
use simple_parameters,             only: params_glob
use simple_builder,                only: build_glob
implicit none

public :: strategy3D_neigh_single
private
#include "simple_local_flags.inc"

type, extends(strategy3D_neigh_multi) :: strategy3D_neigh_single

contains
    procedure :: oris_assign => oris_assign_neigh_single
end type strategy3D_neigh_single

contains

    subroutine oris_assign_neigh_single( self )
        class(strategy3D_neigh_single), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_neigh_single

end module simple_strategy3D_neigh_single
