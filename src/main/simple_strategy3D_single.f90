! concrete strategy3D: probabilistic single-state refinement
module simple_strategy3D_single
include 'simple_lib.f08'
use simple_strategy3D_utils
use simple_strategy3D_multi, only: strategy3D_multi
implicit none

public :: strategy3D_single
private

#include "simple_local_flags.inc"

type, extends(strategy3D_multi) :: strategy3D_single
contains
    procedure :: oris_assign => oris_assign_single
end type strategy3D_single

contains

    subroutine oris_assign_single( self )
        class(strategy3D_single), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_single

end module simple_strategy3D_single
