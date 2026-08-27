module pose_cont_refinement_shift_polish_test
use pose_cont_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use simple_defs, only: dp, sp, DPI
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, &
    &SHIFT_LM_ACCEPTED_IMPROVEMENT, SHIFT_LM_FINITE_NO_IMPROVEMENT, SHIFT_LM_NO_RELIABLE_UPDATE
use pose_cont_refinement_batch_helpers, only: pose_refinement_summary, refine_fixed_volume_shifts
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_shift_polish

integer, parameter :: NPARTICLES = 8, NHALF = NPARTICLES/2
integer, parameter :: MAX_POLISH_ITERATIONS = 20
real(dp), parameter :: RECOVERY_TOL = 2.e-3_dp

contains

!> Verify bounded shift-only LM, half isolation, rollback, and terminal accounting.
subroutine run_shift_polish()
    type(reconstructor_pcg) :: even_operator, odd_operator, weak_operator
    type(cartesian_pose_refiner) :: even_workspace, odd_workspace, weak_workspace
    type(pose_refinement_summary) :: even_summary, odd_summary, weak_summary
    real, allocatable :: even_volume(:,:,:), odd_volume(:,:,:), weak_volume(:,:,:)
    complex, allocatable :: observed(:,:,:), observed_before(:,:,:), zero_plane(:,:)
    real(dp) :: rotmats(3,3,NPARTICLES), rotations_before(3,3,NPARTICLES)
    real(dp) :: true_shifts(2,NPARTICLES), shifts(2,NPARTICLES), shifts_before(2,NPARTICLES)
    real(dp) :: even_shifts(2,NHALF), odd_shifts(2,NHALF), ignored_objective
    real(dp) :: weak_rotation(3,3,1), weak_shift(2,1), weak_shift_before(2,1)
    integer :: half_ids(NPARTICLES), half_ids_before(NPARTICLES)
    integer :: even_ids(NHALF), odd_ids(NHALF), even_statuses(NHALF), odd_statuses(NHALF)
    integer :: weak_status(1)
    integer :: iparticle, lims2(2,2)

    call build_truth_volume(even_volume)
    allocate(odd_volume, source=even_volume(TRUTH_VOLUME_BOX:1:-1,:,:))
    call assert_true(any(even_volume /= odd_volume), 'shift-polish half volumes are identical')
    call even_operator%new(TRUTH_VOLUME_BOX,1._sp)
    call odd_operator%new(TRUTH_VOLUME_BOX,1._sp)
    call even_operator%set_volume(even_volume)
    call odd_operator%set_volume(odd_volume)
    call even_workspace%new_prepared_test(even_volume)
    call odd_workspace%new_prepared_test(odd_volume)
    lims2 = even_operator%get_lims2()
    call assert_true(all(lims2 == odd_operator%get_lims2()), &
        &'even and odd shift-polish workspaces have different Fourier bounds')
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),NPARTICLES))
    allocate(zero_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    zero_plane = cmplx(0.,0.)

    half_ids = [(mod(iparticle-1,2),iparticle=1,NPARTICLES)]
    even_ids = pack([(iparticle,iparticle=1,NPARTICLES)],half_ids==0)
    odd_ids = pack([(iparticle,iparticle=1,NPARTICLES)],half_ids==1)
    call set_fixture_poses(rotmats,true_shifts,shifts)
    do iparticle = 1, NPARTICLES
        if( half_ids(iparticle) == 0 )then
            call even_workspace%shift_residual(rotmats(:,:,iparticle),true_shifts(:,iparticle),&
                &zero_plane,observed(:,:,iparticle),ignored_objective)
        else
            call odd_workspace%shift_residual(rotmats(:,:,iparticle),true_shifts(:,iparticle),&
                &zero_plane,observed(:,:,iparticle),ignored_objective)
        endif
    enddo
    allocate(observed_before, source=observed)
    rotations_before = rotmats
    half_ids_before = half_ids
    shifts_before = shifts

    ! Each half uses its own fixed-volume Fourier workspace.
    even_shifts = shifts(:,even_ids)
    call refine_fixed_volume_shifts(even_workspace,rotmats(:,:,even_ids),observed(:,:,even_ids),&
        &even_shifts,MAX_POLISH_ITERATIONS,even_statuses,even_summary)
    shifts(:,even_ids) = even_shifts
    call assert_true(all(shifts(:,odd_ids) == shifts_before(:,odd_ids)), &
        &'even shift polishing changed odd particle shifts')

    odd_shifts = shifts(:,odd_ids)
    call refine_fixed_volume_shifts(odd_workspace,rotmats(:,:,odd_ids),observed(:,:,odd_ids),&
        &odd_shifts,MAX_POLISH_ITERATIONS,odd_statuses,odd_summary)
    shifts(:,odd_ids) = odd_shifts
    call assert_true(all(shifts(:,even_ids) == even_shifts), &
        &'odd shift polishing changed even particle shifts')

    call assert_half_summary(even_summary,'even')
    call assert_half_summary(odd_summary,'odd')
    call assert_half_statuses(even_statuses,'even')
    call assert_half_statuses(odd_statuses,'odd')
    call assert_true(maxval(sqrt(sum(&
        &(shifts(:,1:NPARTICLES-2)-true_shifts(:,1:NPARTICLES-2))**2,dim=1))) < RECOVERY_TOL, &
        &'shift-polish batch did not recover every injected shift')
    call assert_true(all(shifts(:,NPARTICLES-1:NPARTICLES) == shifts_before(:,NPARTICLES-1:NPARTICLES)), &
        &'exact-solution particles did not preserve their input shifts')
    call assert_true(all(rotmats == rotations_before), 'shift polishing changed particle rotations')
    call assert_true(all(half_ids == half_ids_before), 'shift polishing changed particle half ownership')
    call assert_true(all(observed == observed_before), 'shift polishing changed observed Fourier planes')

    ! A particle with no shift information must retain its incoming metadata.
    allocate(weak_volume, mold=even_volume)
    weak_volume = 0.
    call weak_operator%new(TRUTH_VOLUME_BOX,1._sp)
    call weak_operator%set_volume(weak_volume)
    call weak_workspace%new_prepared_test(weak_volume)
    weak_rotation = 0._dp
    weak_rotation(1,1,1) = 1._dp
    weak_rotation(2,2,1) = 1._dp
    weak_rotation(3,3,1) = 1._dp
    weak_shift(:,1) = [0.23_dp,-0.17_dp]
    weak_shift_before = weak_shift
    call refine_fixed_volume_shifts(weak_workspace,weak_rotation,observed(:,:,1:1)*0.,weak_shift,&
        &MAX_POLISH_ITERATIONS,weak_status,weak_summary)
    call assert_true(weak_status(1) == SHIFT_LM_NO_RELIABLE_UPDATE, &
        &'zero-information shift polish returned the wrong status')
    call assert_true(all(weak_shift == weak_shift_before), &
        &'zero-information shift polish changed the incoming shift')
    call assert_true(weak_summary%nunreliable == 1 .and. weak_summary%nparticles == 1, &
        &'zero-information shift polish was not accounted as unreliable')

    write(*,'(a,4(i0,1x))') 'POSE_CONT_REFINEMENT_SHIFT_POLISH particles/improved/unchanged/steps: ', &
        &even_summary%nparticles+odd_summary%nparticles, even_summary%nimproved+odd_summary%nimproved, &
        &even_summary%nunchanged+odd_summary%nunchanged, &
        &even_summary%naccepted_steps+odd_summary%naccepted_steps
    write(*,'(a,2(es14.6,1x))') 'POSE_CONT_REFINEMENT_SHIFT_POLISH objective before/after: ', &
        &even_summary%objective_before+odd_summary%objective_before, &
        &even_summary%objective_after+odd_summary%objective_after
    write(*,'(a,es14.6)') 'POSE_CONT_REFINEMENT_SHIFT_POLISH maximum trial step: ', &
        &max(even_summary%max_trial_step,odd_summary%max_trial_step)
    write(*,'(a)') 'POSE_CONT_REFINEMENT_SHIFT_POLISH: PASS'

    deallocate(observed,observed_before,zero_plane,even_volume,odd_volume,weak_volume)
    call even_workspace%kill
    call odd_workspace%kill
    call weak_workspace%kill
    call even_operator%kill
    call odd_operator%kill
    call weak_operator%kill
end subroutine run_shift_polish

!> Build fixed rotations and known initial shift errors for both halves.
subroutine set_fixture_poses(rotmats,true_shifts,shifts)
    real(dp), intent(out) :: rotmats(3,3,NPARTICLES)
    real(dp), intent(out) :: true_shifts(2,NPARTICLES), shifts(2,NPARTICLES)
    real(dp), parameter :: angles(NPARTICLES) = &
        &[0._dp,0.13_dp,-0.21_dp,0.31_dp,-0.38_dp,0.47_dp,-0.16_dp,0.26_dp]
    integer :: iparticle

    true_shifts = reshape([ &
        &0.35_dp,-0.27_dp, -0.42_dp,0.18_dp, 0.24_dp,0.39_dp, -0.31_dp,-0.22_dp, &
        &0.46_dp,0.11_dp, -0.19_dp,0.44_dp, 0.28_dp,-0.33_dp, -0.37_dp,0.29_dp], &
        &shape(true_shifts))
    shifts = true_shifts
    shifts(:,1:NPARTICLES-2) = shifts(:,1:NPARTICLES-2) + &
        &reshape([0.22_dp,-0.17_dp, -0.18_dp,0.21_dp, 0.16_dp,0.19_dp, &
        &-0.20_dp,-0.14_dp, 0.17_dp,-0.22_dp, -0.15_dp,0.18_dp], [2,NPARTICLES-2])
    do iparticle = 1, NPARTICLES
        rotmats(:,:,iparticle) = 0._dp
        rotmats(1,1,iparticle) = cos(angles(iparticle)*DPI)
        rotmats(1,2,iparticle) = -sin(angles(iparticle)*DPI)
        rotmats(2,1,iparticle) = sin(angles(iparticle)*DPI)
        rotmats(2,2,iparticle) = cos(angles(iparticle)*DPI)
        rotmats(3,3,iparticle) = 1._dp
    enddo
end subroutine set_fixture_poses

!> Require balanced shift-polish outcomes and an aggregate objective reduction.
subroutine assert_half_summary(summary,label)
    type(pose_refinement_summary), intent(in) :: summary
    character(len=*), intent(in) :: label
    integer :: accounted

    accounted = summary%nimproved + summary%nunchanged + summary%nunreliable + &
        &summary%nstep_bound + summary%ninvalid + summary%niteration_limit
    call assert_true(accounted == summary%nparticles, trim(label)//' shift-polish summary lost a particle')
    call assert_true(summary%nimproved == NHALF-1, trim(label)//' shift-polish improvement count is wrong')
    call assert_true(summary%nunchanged == 1, trim(label)//' exact-solution count is wrong')
    call assert_true(summary%objective_after < summary%objective_before, &
        &trim(label)//' shift-polish batch did not reduce its aggregate objective')
    call assert_true(summary%max_trial_step <= 1._dp+epsilon(1._dp), &
        &trim(label)//' shift-polish batch exceeded the one-pixel trial bound')
end subroutine assert_half_summary

!> Verify the expected accepted and finite-no-improvement terminal statuses.
subroutine assert_half_statuses(statuses,label)
    integer, intent(in) :: statuses(NHALF)
    character(len=*), intent(in) :: label

    call assert_true(count(statuses == SHIFT_LM_ACCEPTED_IMPROVEMENT) == NHALF-1, &
        &trim(label)//' shift-polish accepted-status count is wrong')
    call assert_true(count(statuses == SHIFT_LM_FINITE_NO_IMPROVEMENT) == 1, &
        &trim(label)//' shift-polish no-improvement status count is wrong')
end subroutine assert_half_statuses

end module pose_cont_refinement_shift_polish_test
