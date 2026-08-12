module continuous_3D_pcg_refinement_recovery_test
use continuous_3D_pcg_refinement_test_helpers, only: skip_unimplemented_case
implicit none
private
public :: run_pose_recovery

contains

subroutine run_pose_recovery()
    call skip_unimplemented_case('pose_recovery', 'Stage 3/4 -- local pose and alternating-loop recovery')
end subroutine run_pose_recovery

end module continuous_3D_pcg_refinement_recovery_test
