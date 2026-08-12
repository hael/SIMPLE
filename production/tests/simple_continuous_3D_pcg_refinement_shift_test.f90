module continuous_3D_pcg_refinement_shift_test
use continuous_3D_pcg_refinement_test_helpers, only: skip_unimplemented_case
implicit none
private
public :: run_shift_gradient

contains

subroutine run_shift_gradient()
    call skip_unimplemented_case('shift_gradient', 'Stage 1 -- analytic shift derivatives')
end subroutine run_shift_gradient

end module continuous_3D_pcg_refinement_shift_test
