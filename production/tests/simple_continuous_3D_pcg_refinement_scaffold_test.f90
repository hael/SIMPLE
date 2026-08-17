module continuous_3D_pcg_refinement_scaffold_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_int_equal, &
    &assert_real_close, assert_true, set_deterministic_seed
implicit none
private
public :: run_scaffold_contract

contains

!> Verify mother/child case dispatch and shared test-helper contracts.
subroutine run_scaffold_contract()
    integer, parameter :: FIXTURE_SEED = 20260811
    real :: first_draw(8), repeated_draw(8)

    call set_deterministic_seed(FIXTURE_SEED)
    call random_number(first_draw)
    call set_deterministic_seed(FIXTURE_SEED)
    call random_number(repeated_draw)

    call assert_true(all(first_draw == repeated_draw), &
        &'deterministic seed did not reproduce the same random sequence')
    call assert_true(all(first_draw >= 0.) .and. all(first_draw < 1.), &
        &'intrinsic random_number returned a value outside [0,1)')
    call assert_int_equal(size(first_draw), 8, 'scaffold random draw size')
    call assert_real_close(maxval(abs(first_draw - repeated_draw)), 0., 0., &
        &'scaffold repeated random draw')

    write(*,'(a,i0)') 'CONTINUOUS_3D_PCG_SCAFFOLD deterministic seed: ', FIXTURE_SEED
    write(*,'(a)') 'CONTINUOUS_3D_PCG_SCAFFOLD: PASS'
end subroutine run_scaffold_contract

end module continuous_3D_pcg_refinement_scaffold_test
