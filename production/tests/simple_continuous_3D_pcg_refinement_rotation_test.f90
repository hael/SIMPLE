module continuous_3D_pcg_refinement_rotation_test
use continuous_3D_pcg_refinement_test_helpers, only: skip_unimplemented_case
implicit none
private
public :: run_rotation_gradient

contains

subroutine run_rotation_gradient()
    call skip_unimplemented_case('rotation_gradient', 'Stage 3 -- tangent-space rotation derivatives')
end subroutine run_rotation_gradient

end module continuous_3D_pcg_refinement_rotation_test
