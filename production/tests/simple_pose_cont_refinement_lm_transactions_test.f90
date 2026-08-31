! Acceptance test for the scaled pose LM proposal and complete-pose transaction.
module pose_cont_refinement_lm_transactions_test
use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use pose_cont_refinement_calibration_helpers, only: calibration_fixture, FROZEN_ABSOLUTE_TOLERANCES, &
    &FROZEN_RELATIVE_TOLERANCES, build_acceptance_fixture, combined_real_passes
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp, sp, PI, OSMPL_PAD_FAC
use simple_linalg, only: hermitian_solve
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, right_increment_rotation, &
    &pose_lm_system_test, pose_lm_transaction_test, POSE_LM_ACCEPTED_IMPROVEMENT, &
    &POSE_LM_FINITE_NO_IMPROVEMENT, POSE_LM_NO_RELIABLE_UPDATE, &
    &POSE_LM_STEP_BOUND_REJECTED, POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
implicit none
private
public :: run_lm_transactions

integer, parameter :: LM_TOLERANCE = 2
real(dp), parameter :: ROTATION_SCALE = 0.10_dp
real(dp), parameter :: MU = 1.e-3_dp

contains

!> Compare the shared LM mechanism with pivoted dense and LAPACK solves.
subroutine run_lm_transactions()
    real(dp) :: hessian(5,5), gradient(5), singular_hessian(5,5), invalid_gradient(5)
    real(dp) :: input_rotmat(3,3), trial_rotmat(3,3), input_shift(2), trial_shift(2)
    character(len=:), allocatable :: evidence_directory
    integer :: system_unit, transaction_unit

    evidence_directory = required_argument('evidence_dir')
    open(newunit=system_unit,file=evidence_directory//'/lm_systems.tsv', &
        &status='replace',action='write')
    open(newunit=transaction_unit,file=evidence_directory//'/pose_transactions.tsv', &
        &status='replace',action='write')
    write(system_unit,'(a)') 'case'//achar(9)//'target_ratio'//achar(9)//'actual_ratio'// &
        &achar(9)//'pivoted_backward_error'//achar(9)//'lapack_backward_error'// &
        &achar(9)//'lapack_info'//achar(9)//'bounded'//achar(9)//'accepted'
    write(transaction_unit,'(a)') 'status'//achar(9)//'rotation_exact'//achar(9)//'shift_exact'

    call build_reliable_system(hessian,gradient)
    call test_reliable_system('accepted',hessian,gradient,0.80_dp,.false.,.true.,system_unit)
    call test_reliable_system('rejected',hessian,gradient,0.10_dp,.false.,.false.,system_unit)
    call test_reliable_system('step_bound',hessian,100._dp*gradient,0.10_dp,.true.,.false.,system_unit)

    singular_hessian = 0._dp
    call test_unreliable_system('singular',singular_hessian,gradient,system_unit)
    invalid_gradient = gradient
    invalid_gradient(3) = ieee_value(0._dp,ieee_quiet_nan)
    call test_unreliable_system('invalid',hessian,invalid_gradient,system_unit)
    call test_executed_transactions()

    input_rotmat = reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    trial_rotmat = reshape([0._dp,1._dp,0._dp,-1._dp,0._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    input_shift = [0.25_dp,-0.50_dp]
    trial_shift = [-0.75_dp,0.80_dp]
    call test_pose_transaction('rejected',0,input_rotmat,input_shift,trial_rotmat,trial_shift, &
        &.false.,transaction_unit)
    call test_pose_transaction('finite_no_improvement',POSE_LM_FINITE_NO_IMPROVEMENT, &
        &input_rotmat,input_shift,trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('singular_unreliable',POSE_LM_NO_RELIABLE_UPDATE, &
        &input_rotmat,input_shift,trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('invalid',POSE_LM_INVALID_NUMERICS,input_rotmat,input_shift, &
        &trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('cumulative_bound',POSE_LM_STEP_BOUND_REJECTED,input_rotmat, &
        &input_shift,trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('step_bound',POSE_LM_STEP_BOUND_REJECTED,input_rotmat,input_shift, &
        &trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('iteration_limit',POSE_LM_ITERATION_LIMIT,input_rotmat,input_shift, &
        &trial_rotmat,trial_shift,.false.,transaction_unit)
    call test_pose_transaction('accepted',POSE_LM_ACCEPTED_IMPROVEMENT,input_rotmat,input_shift, &
        &trial_rotmat,trial_shift, &
        &.true.,transaction_unit)
    close(system_unit)
    close(transaction_unit)

    write(*,'(a)') 'POSE_CONT_REFINEMENT_LM_TRANSACTIONS ratios: accepted=0.80 rejected=0.10'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_LM_TRANSACTIONS solves: cholesky pivoted_gaussian lapack_dposv'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_LM_TRANSACTIONS outcomes: finite singular invalid bounds limit'
    write(*,'(a)') 'POSE_CONT_REFINEMENT_LM_TRANSACTIONS: PASS'
end subroutine run_lm_transactions

!> Exercise accepted traces and production terminal rollback on a physical fixture.
subroutine test_executed_transactions()
    integer, parameter :: MAX_ITERATIONS = 20
    type(calibration_fixture) :: fixture
    type(cartesian_pose_refiner) :: workspace
    complex, allocatable :: observed(:,:), invalid_observed(:,:), zero_transfer(:,:)
    real(dp) :: identity(3,3), truth_rotmat(3,3), seed_rotmat(3,3), rotmat(3,3)
    real(dp) :: truth_shift(2), seed_shift(2), shift(2), input_rotmat(3,3), input_shift(2)
    real(dp) :: objectives(0:MAX_ITERATIONS), max_rotation_step, max_shift_step, nan_value
    integer :: naccepted, nattempted, nstencil_switches, status

    call build_acceptance_fixture(10,fixture)
    call workspace%new(real(fixture%volume,sp))
    call workspace%set_shell_range([0,fixture%box/2])
    identity = reshape([1._dp,0._dp,0._dp,0._dp,1._dp,0._dp,0._dp,0._dp,1._dp],[3,3])
    truth_rotmat = right_increment_rotation(identity,fixture%exact_pose(1:3))
    truth_shift = fixture%exact_pose(4:5)
    seed_rotmat = right_increment_rotation(identity,fixture%nonstationary_pose(1:3))
    seed_shift = fixture%nonstationary_pose(4:5)
    call build_observation(workspace,fixture%box,truth_rotmat,truth_shift,observed)

    rotmat = seed_rotmat
    shift = seed_shift
    call workspace%refine_pose_lm(rotmat,observed,shift,ROTATION_SCALE,MAX_ITERATIONS, &
        &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)
    call assert_true(status == POSE_LM_ACCEPTED_IMPROVEMENT .and. naccepted > 0, &
        &'executed LM fixture produced no accepted trace')
    call assert_true(all(objectives(1:naccepted) < objectives(0:naccepted-1)), &
        &'executed accepted-objective trace is not strictly decreasing')

    rotmat = truth_rotmat
    shift = truth_shift
    input_rotmat = rotmat
    input_shift = shift
    call workspace%refine_pose_lm(rotmat,observed,shift,ROTATION_SCALE,MAX_ITERATIONS, &
        &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)
    call assert_true(status == POSE_LM_FINITE_NO_IMPROVEMENT .and. naccepted == 0, &
        &'exact pose did not return finite no improvement')
    call assert_true(all(rotmat == input_rotmat) .and. all(shift == input_shift), &
        &'finite no-improvement changed the complete input pose')

    allocate(zero_transfer,mold=observed)
    zero_transfer = cmplx(0.,0.)
    rotmat = seed_rotmat
    shift = seed_shift
    input_rotmat = rotmat
    input_shift = shift
    call workspace%refine_pose_lm(rotmat,observed,shift,ROTATION_SCALE,MAX_ITERATIONS, &
        &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nstencil_switches, &
        &transfer=zero_transfer)
    call assert_true(status == POSE_LM_NO_RELIABLE_UPDATE .and. naccepted == 0, &
        &'singular pose did not return no reliable update')
    call assert_true(all(rotmat == input_rotmat) .and. all(shift == input_shift), &
        &'singular pose transaction changed the complete input pose')

    allocate(invalid_observed,source=observed)
    nan_value = ieee_value(0._dp,ieee_quiet_nan)
    invalid_observed(0,0) = cmplx(real(nan_value,sp),0._sp,kind=sp)
    rotmat = seed_rotmat
    shift = seed_shift
    input_rotmat = rotmat
    input_shift = shift
    call workspace%refine_pose_lm(rotmat,invalid_observed,shift,ROTATION_SCALE,MAX_ITERATIONS, &
        &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)
    call assert_true(status == POSE_LM_INVALID_NUMERICS .and. naccepted == 0, &
        &'non-finite pose did not return invalid numerics')
    call assert_true(all(rotmat == input_rotmat) .and. all(shift == input_shift), &
        &'invalid pose transaction changed the complete input pose')

    rotmat = seed_rotmat
    shift = seed_shift
    input_rotmat = rotmat
    input_shift = shift
    call workspace%refine_pose_lm(rotmat,observed,shift,ROTATION_SCALE,1,objectives, &
        &naccepted,status,nattempted,max_rotation_step,max_shift_step,nstencil_switches, &
        &anchor_rotmat=input_rotmat,anchor_shift=input_shift,max_total_rotation=1.e-12_dp, &
        &max_total_shift=1.e-12_dp)
    call assert_true(status == POSE_LM_STEP_BOUND_REJECTED .and. naccepted == 0, &
        &'cumulative-bound pose did not return step-bound rejected')
    call assert_true(all(rotmat == input_rotmat) .and. all(shift == input_shift), &
        &'cumulative-bound transaction changed the complete input pose')
    call workspace%kill
end subroutine test_executed_transactions

!> Build the full redundant-disk observation at one known pose.
subroutine build_observation(workspace,box,rotmat,shift,observed)
    type(cartesian_pose_refiner), intent(in) :: workspace
    integer, intent(in) :: box
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    complex, allocatable, intent(out) :: observed(:,:)
    complex :: value, derivative(3), phase
    real :: margin(3)
    real(sp) :: location(3)
    real(dp) :: argument
    integer :: h, k, lims2(2,2)

    lims2 = workspace%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)),source=cmplx(0.,0.))
    do k = lims2(2,1), lims2(2,2)
        do h = lims2(1,1), lims2(1,2)
            if( h*h+k*k > (box/2)**2 ) cycle
            location = real(OSMPL_PAD_FAC,sp)*real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(location,value,derivative,margin)
            argument = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(box,dp)
            phase = cmplx(cos(argument),sin(argument),kind=sp)
            observed(h,k) = phase*value
        enddo
    enddo
end subroutine build_observation

!> Define a deterministic symmetric positive-definite five-parameter block.
subroutine build_reliable_system(hessian,gradient)
    real(dp), intent(out) :: hessian(5,5), gradient(5)

    hessian = 0._dp
    hessian(1,1) = 12._dp
    hessian(2,2) = 10._dp
    hessian(3,3) = 9._dp
    hessian(4,4) = 8._dp
    hessian(5,5) = 7._dp
    hessian(1,2) = 0.7_dp
    hessian(2,1) = hessian(1,2)
    hessian(2,3) = -0.4_dp
    hessian(3,2) = hessian(2,3)
    hessian(3,4) = 0.5_dp
    hessian(4,3) = hessian(3,4)
    hessian(4,5) = -0.3_dp
    hessian(5,4) = hessian(4,5)
    gradient = [0.20_dp,-0.15_dp,0.10_dp,-0.18_dp,0.12_dp]
end subroutine build_reliable_system

!> Compare scaling, damping, dense solution, step, reductions, ratio, and decision.
subroutine test_reliable_system(name,hessian,gradient,target_ratio,expect_bounded, &
    &expect_accepted,unit)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: hessian(5,5), gradient(5), target_ratio
    logical, intent(in) :: expect_bounded, expect_accepted
    integer, intent(in) :: unit
    real(dp) :: actual_gradient(5), actual_hessian(5,5), actual_diagonal(5)
    real(dp) :: actual_matrix(5,5), actual_scaled_step(5), actual_step(5)
    real(dp) :: expected_gradient(5), expected_hessian(5,5), expected_diagonal(5)
    real(dp) :: expected_matrix(5,5), expected_scaled_step(5), expected_step(5)
    real(dp) :: lapack_scaled_step(5), lapack_step(5), lapack_residual(5)
    real(dp) :: objective, trial_objective, expected_predicted, actual_predicted
    real(dp) :: expected_actual, actual_reduction, actual_ratio, backward_error
    real(dp) :: lapack_backward_error, lapack_denominator
    integer :: lapack_info
    logical :: identifiable, stationary, reliable, bounded, accepted
    call independent_lm_system(hessian,gradient,expected_gradient,expected_hessian, &
        &expected_diagonal,expected_matrix,expected_scaled_step,expected_step,backward_error)
    call hermitian_solve(expected_matrix,-expected_gradient,lapack_scaled_step,lapack_info)
    call assert_true(lapack_info == 0,'LAPACK DPOSV rejected a reliable LM system')
    call bound_physical_step(lapack_scaled_step,lapack_step)
    lapack_residual = matmul(expected_matrix,lapack_scaled_step)+expected_gradient
    lapack_denominator = maxval(sum(abs(expected_matrix),dim=2))* &
        &maxval(abs(lapack_scaled_step))+maxval(abs(expected_gradient))
    lapack_backward_error = maxval(abs(lapack_residual))/lapack_denominator
    expected_predicted = -dot_product(gradient,expected_step)-0.5_dp* &
        &dot_product(expected_step,matmul(hessian,expected_step))
    objective = 3._dp
    trial_objective = objective-target_ratio*expected_predicted
    expected_actual = objective-trial_objective
    call pose_lm_system_test(gradient,hessian,ROTATION_SCALE,MU,objective,trial_objective, &
        &[.true.,.true.,.true.,.true.,.true.],actual_gradient,actual_hessian,actual_diagonal, &
        &actual_matrix,actual_scaled_step,actual_step,actual_predicted,actual_reduction, &
        &actual_ratio,identifiable,stationary,reliable,bounded,accepted)

    call assert_true(identifiable .and. .not. stationary .and. reliable, &
        &'reliable LM fixture was not solvable')
    call assert_true((bounded .eqv. expect_bounded),'LM step-bound decision differs')
    call assert_true((accepted .eqv. expect_accepted),'LM gain-ratio decision differs')
    call assert_combined(actual_gradient,expected_gradient,'scaled gradient')
    call assert_combined(reshape(actual_hessian,[25]),reshape(expected_hessian,[25]), &
        &'scaled Hessian')
    call assert_combined(actual_diagonal,expected_diagonal,'damping diagonal')
    call assert_combined(reshape(actual_matrix,[25]),reshape(expected_matrix,[25]), &
        &'damped matrix')
    call assert_combined(actual_scaled_step,expected_scaled_step,'scaled dense step')
    call assert_combined(actual_step,expected_step,'physical pose step')
    call assert_combined(lapack_scaled_step,expected_scaled_step,'LAPACK scaled step')
    call assert_combined(lapack_step,expected_step,'LAPACK physical pose step')
    call assert_combined([actual_predicted,actual_reduction,actual_ratio], &
        &[expected_predicted,expected_actual,target_ratio],'LM reductions and ratio')
    call assert_true(backward_error <= FROZEN_RELATIVE_TOLERANCES(LM_TOLERANCE), &
        &'independent dense solve has excessive backward error')
    call assert_true(lapack_backward_error <= FROZEN_RELATIVE_TOLERANCES(LM_TOLERANCE), &
        &'LAPACK DPOSV solve has excessive backward error')
    write(unit,'(a,4(a,es24.16),a,i0,2(a,l1))') trim(name),achar(9),target_ratio,achar(9), &
        &actual_ratio,achar(9),backward_error,achar(9),lapack_backward_error,achar(9), &
        &lapack_info,achar(9),bounded,achar(9),accepted
end subroutine test_reliable_system

!> Verify that singular and non-finite systems do not produce a reliable step.
subroutine test_unreliable_system(name,hessian,gradient,unit)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: hessian(5,5), gradient(5)
    integer, intent(in) :: unit
    real(dp) :: scaled_gradient(5), scaled_hessian(5,5), diagonal(5), matrix(5,5)
    real(dp) :: scaled_step(5), step(5), predicted, actual, ratio
    logical :: identifiable, stationary, reliable, bounded, accepted

    call pose_lm_system_test(gradient,hessian,ROTATION_SCALE,MU,3._dp,2._dp, &
        &[.true.,.true.,.true.,.true.,.true.],scaled_gradient,scaled_hessian,diagonal, &
        &matrix,scaled_step,step,predicted,actual,ratio,identifiable,stationary,reliable, &
        &bounded,accepted)
    call assert_true(.not. reliable .and. .not. accepted, &
        &'singular or invalid LM fixture produced an accepted step')
    write(unit,'(a,4(a,es24.16),a,i0,2(a,l1))') trim(name),achar(9),0._dp,achar(9),ratio, &
        &achar(9),0._dp,achar(9),0._dp,achar(9),-1,achar(9),bounded,achar(9),accepted
end subroutine test_unreliable_system

!> Solve the scaled damped system with pivoted Gaussian elimination.
subroutine independent_lm_system(hessian,gradient,scaled_gradient,scaled_hessian,diagonal, &
    &matrix,scaled_step,physical_step,backward_error)
    real(dp), intent(in) :: hessian(5,5), gradient(5)
    real(dp), intent(out) :: scaled_gradient(5), scaled_hessian(5,5), diagonal(5)
    real(dp), intent(out) :: matrix(5,5), scaled_step(5), physical_step(5), backward_error
    real(dp) :: coordinate_scale(5), hessian_scale, residual(5), denominator
    integer :: axis, jaxis

    coordinate_scale = [ROTATION_SCALE,ROTATION_SCALE,ROTATION_SCALE,1._dp,1._dp]
    do axis = 1, 5
        scaled_gradient(axis) = coordinate_scale(axis)*gradient(axis)
        do jaxis = 1, 5
            scaled_hessian(axis,jaxis) = &
                &coordinate_scale(axis)*hessian(axis,jaxis)*coordinate_scale(jaxis)
        enddo
    enddo
    hessian_scale = maxval(abs(scaled_hessian))
    do axis = 1, 5
        diagonal(axis) = max(scaled_hessian(axis,axis), &
            &sqrt(epsilon(1._dp))*hessian_scale,epsilon(1._dp)**2)
    enddo
    matrix = scaled_hessian
    do axis = 1, 5
        matrix(axis,axis) = matrix(axis,axis)+MU*diagonal(axis)
    enddo
    call gaussian_solve(matrix,-scaled_gradient,scaled_step)
    residual = matmul(matrix,scaled_step)+scaled_gradient
    denominator = maxval(sum(abs(matrix),dim=2))*maxval(abs(scaled_step))+ &
        &maxval(abs(scaled_gradient))
    backward_error = maxval(abs(residual))/denominator
    call bound_physical_step(scaled_step,physical_step)
end subroutine independent_lm_system

!> Convert a scaled solution to the physical bounded pose step.
subroutine bound_physical_step(scaled_step,physical_step)
    real(dp), intent(in) :: scaled_step(5)
    real(dp), intent(out) :: physical_step(5)
    real(dp) :: coordinate_scale(5), rotation_norm, shift_norm

    coordinate_scale = [ROTATION_SCALE,ROTATION_SCALE,ROTATION_SCALE,1._dp,1._dp]
    physical_step = coordinate_scale*scaled_step
    rotation_norm = sqrt(sum(physical_step(1:3)**2))
    if( rotation_norm > ROTATION_SCALE ) &
        &physical_step(1:3) = physical_step(1:3)*ROTATION_SCALE/rotation_norm
    shift_norm = sqrt(sum(physical_step(4:5)**2))
    if( shift_norm > 1._dp ) physical_step(4:5) = physical_step(4:5)/shift_norm
end subroutine bound_physical_step

!> Dense partial-pivot elimination independent of the production Cholesky path.
subroutine gaussian_solve(matrix,rhs,solution)
    real(dp), intent(in) :: matrix(5,5), rhs(5)
    real(dp), intent(out) :: solution(5)
    real(dp) :: work(5,5), work_rhs(5), factor, row(5), rhs_value
    integer :: column, pivot, row_index

    work = matrix
    work_rhs = rhs
    do column = 1, 4
        pivot = column
        do row_index = column+1, 5
            if( abs(work(row_index,column)) > abs(work(pivot,column)) ) pivot = row_index
        enddo
        call assert_true(abs(work(pivot,column)) > epsilon(1._dp), &
            &'independent dense solve found a zero pivot')
        if( pivot /= column )then
            row = work(column,:)
            work(column,:) = work(pivot,:)
            work(pivot,:) = row
            rhs_value = work_rhs(column)
            work_rhs(column) = work_rhs(pivot)
            work_rhs(pivot) = rhs_value
        endif
        do row_index = column+1, 5
            factor = work(row_index,column)/work(column,column)
            work(row_index,column:5) = work(row_index,column:5)- &
                &factor*work(column,column:5)
            work_rhs(row_index) = work_rhs(row_index)-factor*work_rhs(column)
        enddo
    enddo
    do row_index = 5, 1, -1
        solution(row_index) = (work_rhs(row_index)- &
            &dot_product(work(row_index,row_index+1:5),solution(row_index+1:5)))/ &
            &work(row_index,row_index)
    enddo
end subroutine gaussian_solve

!> Apply the frozen LM componentwise combined comparison rule.
subroutine assert_combined(actual,expected,label)
    real(dp), intent(in) :: actual(:), expected(:)
    character(len=*), intent(in) :: label

    call assert_true(combined_real_passes(actual,expected, &
        &FROZEN_ABSOLUTE_TOLERANCES(LM_TOLERANCE), &
        &FROZEN_RELATIVE_TOLERANCES(LM_TOLERANCE)), &
        &trim(label)//' differs from the independent LM oracle')
end subroutine assert_combined

!> Verify exact complete-pose selection for one terminal transaction class.
subroutine test_pose_transaction(name,status,input_rotmat,input_shift,trial_rotmat,trial_shift, &
    &accepted,unit)
    character(len=*), intent(in) :: name
    integer, intent(in) :: status, unit
    real(dp), intent(in) :: input_rotmat(3,3), input_shift(2), trial_rotmat(3,3), trial_shift(2)
    logical, intent(in) :: accepted
    real(dp) :: output_rotmat(3,3), output_shift(2), expected_rotmat(3,3), expected_shift(2)

    call pose_lm_transaction_test(input_rotmat,input_shift,trial_rotmat,trial_shift,accepted, &
        &output_rotmat,output_shift)
    expected_rotmat = merge(trial_rotmat,input_rotmat,accepted)
    expected_shift = merge(trial_shift,input_shift,accepted)
    call assert_true(all(output_rotmat == expected_rotmat) .and. all(output_shift == expected_shift), &
        &'pose transaction did not retain or commit both pose fields exactly')
    write(unit,'(a,a,i0,2(a,l1))') trim(name),achar(9),status,achar(9), &
        &all(output_rotmat == expected_rotmat),achar(9),all(output_shift == expected_shift)
end subroutine test_pose_transaction

!> Read one required key=value command argument.
function required_argument(key) result(value)
    character(len=*), intent(in) :: key
    character(len=:), allocatable :: value
    character(len=4096) :: argument
    integer :: argument_index, argument_status, separator

    value = ''
    do argument_index = 1, command_argument_count()
        call get_command_argument(argument_index,argument,status=argument_status)
        if( argument_status /= 0 ) error stop 'could not read LM transaction argument'
        separator = index(argument,'=')
        if( separator <= 0 ) cycle
        if( trim(argument(:separator-1)) /= key ) cycle
        value = trim(argument(separator+1:))
    enddo
    if( len(value) == 0 ) &
        &error stop 'LM transaction test requires evidence_dir=<existing-directory>'
end function required_argument

end module pose_cont_refinement_lm_transactions_test
