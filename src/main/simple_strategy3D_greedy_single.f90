! concrete strategy3D: greedy single-state refinement
module simple_strategy3D_greedy_single
include 'simple_lib.f08'
use simple_strategy3D_utils
use simple_strategy3D_greedy_multi, only: strategy3D_greedy_multi
implicit none

public :: strategy3D_greedy_single
private
#include "simple_local_flags.inc"

type, extends(strategy3D_greedy_multi) :: strategy3D_greedy_single

contains
    procedure :: oris_assign => oris_assign_greedy_single
end type strategy3D_greedy_single

contains

    subroutine oris_assign_greedy_single( self )
        class(strategy3D_greedy_single), intent(inout) :: self
        call extract_peak_ori(self%s)
    end subroutine oris_assign_greedy_single

end module simple_strategy3D_greedy_single
