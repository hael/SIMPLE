module continuous_3D_pcg_refinement_pose_contract_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp
use simple_core_module_api, only: euler2m
use simple_pcg_pose_polisher, only: pcg_pose_polish_summary, polish_fixed_volume_poses
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, &
    &POSE_LM_ACCEPTED_IMPROVEMENT
implicit none
private
public :: run_pose_contract

integer, parameter :: NPARTICLES = 6, NHALF = NPARTICLES/2
integer, parameter :: MAX_ITERATIONS = 20
real(dp), parameter :: ROTATION_SCALE = 0.10_dp
real(dp), parameter :: ROTATION_TOL = 2.e-3_dp
real(dp), parameter :: SHIFT_TOL = 2.e-3_dp
real(dp), parameter :: FRAME_TOL = 1.e-3_dp

contains

subroutine run_pose_contract()
    type(reconstructor_pcg) :: even_operator, odd_operator
    type(pcg_fourier_workspace) :: even_workspace, odd_workspace
    type(pcg_pose_polish_summary) :: even_summary, odd_summary
    real, allocatable :: even_volume(:,:,:), odd_volume(:,:,:)
    complex, allocatable :: observed(:,:,:), zero_plane(:,:)
    real(dp) :: truth_rotmats(3,3,NPARTICLES), estimate_rotmats(3,3,NPARTICLES)
    real(dp) :: truth_shifts(2,NPARTICLES), estimate_shifts(2,NPARTICLES)
    real(dp) :: even_rotmats(3,3,NHALF), odd_rotmats(3,3,NHALF)
    real(dp) :: even_shifts(2,NHALF), odd_shifts(2,NHALF)
    real(dp) :: odd_rotmats_before(3,3,NHALF), odd_shifts_before(2,NHALF)
    real(dp) :: max_rotation_error, max_shift_error, ignored_objective
    integer :: even_statuses(NHALF), odd_statuses(NHALF), lims2(2,2), iparticle

    call build_truth_volume(even_volume)
    allocate(odd_volume,source=even_volume(TRUTH_VOLUME_BOX:1:-1,:,:))
    call assert_true(any(even_volume /= odd_volume), &
        &'joint pose contract requires distinct even and odd half volumes')
    call even_operator%new(TRUTH_VOLUME_BOX,1._sp)
    call odd_operator%new(TRUTH_VOLUME_BOX,1._sp)
    call even_operator%set_volume(even_volume)
    call odd_operator%set_volume(odd_volume)
    call even_operator%begin_fourier_workspace(even_workspace)
    call odd_operator%begin_fourier_workspace(odd_workspace)
    call even_workspace%set_shell_range([2,TRUTH_VOLUME_BOX/2])
    call odd_workspace%set_shell_range([2,TRUTH_VOLUME_BOX/2])
    lims2 = even_operator%get_lims2()
    call assert_true(all(lims2 == odd_operator%get_lims2()), &
        &'joint pose half workspaces have different Fourier bounds')
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),NPARTICLES))
    allocate(zero_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    zero_plane = cmplx(0.,0.)

    call build_pose_fixture(truth_rotmats,estimate_rotmats,truth_shifts,estimate_shifts)
    do iparticle = 1, NPARTICLES
        if( iparticle <= NHALF )then
            call even_workspace%shift_residual(truth_rotmats(:,:,iparticle), &
                &truth_shifts(:,iparticle),zero_plane,observed(:,:,iparticle),ignored_objective)
        else
            call odd_workspace%shift_residual(truth_rotmats(:,:,iparticle), &
                &truth_shifts(:,iparticle),zero_plane,observed(:,:,iparticle),ignored_objective)
        endif
    enddo

    even_rotmats = estimate_rotmats(:,:,1:NHALF)
    even_shifts = estimate_shifts(:,1:NHALF)
    odd_rotmats = estimate_rotmats(:,:,NHALF+1:NPARTICLES)
    odd_shifts = estimate_shifts(:,NHALF+1:NPARTICLES)
    odd_rotmats_before = odd_rotmats
    odd_shifts_before = odd_shifts

    ! Each call owns one fixed half-map and only its matching particles.
    call polish_fixed_volume_poses(even_workspace,even_rotmats,observed(:,:,1:NHALF), &
        &even_shifts,ROTATION_SCALE,MAX_ITERATIONS,even_statuses,even_summary)
    call assert_true(all(odd_rotmats == odd_rotmats_before) .and. all(odd_shifts == odd_shifts_before), &
        &'even pose polishing changed odd particle metadata')
    call polish_fixed_volume_poses(odd_workspace,odd_rotmats,observed(:,:,NHALF+1:NPARTICLES), &
        &odd_shifts,ROTATION_SCALE,MAX_ITERATIONS,odd_statuses,odd_summary)
    estimate_rotmats(:,:,1:NHALF) = even_rotmats
    estimate_rotmats(:,:,NHALF+1:NPARTICLES) = odd_rotmats
    estimate_shifts(:,1:NHALF) = even_shifts
    estimate_shifts(:,NHALF+1:NPARTICLES) = odd_shifts

    call assert_summary(even_summary,even_statuses,'even')
    call assert_summary(odd_summary,odd_statuses,'odd')
    max_rotation_error = 0._dp
    max_shift_error = 0._dp
    do iparticle = 1, NPARTICLES
        max_rotation_error = max(max_rotation_error, &
            &rotation_distance(estimate_rotmats(:,:,iparticle),truth_rotmats(:,:,iparticle)))
        max_shift_error = max(max_shift_error,sqrt(sum( &
            &(estimate_shifts(:,iparticle)-truth_shifts(:,iparticle))**2)))
    enddo
    call assert_true(max_rotation_error < ROTATION_TOL, &
        &'multiple-orientation joint polishing did not recover every rotation')
    call assert_true(max_shift_error < SHIFT_TOL, &
        &'multiple-orientation joint polishing did not recover every shift')
    call test_symmetry_and_global_gauge(truth_rotmats)
    call assert_true(all(ieee_is_finite([max_rotation_error,max_shift_error])), &
        &'joint pose contract produced non-finite recovery errors')

    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_POSE_CONTRACT max rotation/shift error: ', &
        &max_rotation_error,max_shift_error
    write(*,'(a,4(i0,1x))') 'CONTINUOUS_3D_PCG_POSE_CONTRACT particles/improved/accepted/switches: ', &
        &even_summary%nparticles+odd_summary%nparticles,even_summary%nimproved+odd_summary%nimproved, &
        &even_summary%naccepted_steps+odd_summary%naccepted_steps, &
        &even_summary%nstencil_switches+odd_summary%nstencil_switches
    write(*,'(a)') 'CONTINUOUS_3D_PCG_POSE_CONTRACT: PASS'

    deallocate(even_volume,odd_volume,observed,zero_plane)
    call even_workspace%kill
    call odd_workspace%kill
    call even_operator%kill
    call odd_operator%kill
end subroutine run_pose_contract

subroutine build_pose_fixture(truth_rotmats,estimate_rotmats,truth_shifts,estimate_shifts)
    real(dp), intent(out) :: truth_rotmats(3,3,NPARTICLES), estimate_rotmats(3,3,NPARTICLES)
    real(dp), intent(out) :: truth_shifts(2,NPARTICLES), estimate_shifts(2,NPARTICLES)
    real, parameter :: truth_eulers(3,NPARTICLES) = reshape([ &
        &11.,28.,17., 47.,39.,63., 83.,54.,22., &
        &129.,31.,74., 168.,67.,141., 233.,46.,197.],[3,NPARTICLES])
    real, parameter :: estimate_eulers(3,NPARTICLES) = reshape([ &
        &11.4,27.7,17.3, 46.6,39.3,62.7, 83.3,53.6,22.4, &
        &128.7,31.4,73.6, 168.4,66.7,141.3, 232.6,46.3,196.7],[3,NPARTICLES])
    integer :: iparticle

    truth_shifts = reshape([0.31_dp,-0.24_dp, -0.27_dp,0.18_dp, 0.22_dp,0.35_dp, &
        &-0.34_dp,-0.16_dp, 0.29_dp,0.13_dp, -0.19_dp,0.32_dp],[2,NPARTICLES])
    estimate_shifts = truth_shifts+reshape([0.16_dp,-0.12_dp, -0.14_dp,0.15_dp, &
        &0.13_dp,0.17_dp, -0.15_dp,-0.11_dp, 0.12_dp,-0.16_dp, -0.13_dp,0.14_dp], &
        &[2,NPARTICLES])
    do iparticle = 1, NPARTICLES
        truth_rotmats(:,:,iparticle) = real(euler2m(truth_eulers(:,iparticle)),dp)
        estimate_rotmats(:,:,iparticle) = real(euler2m(estimate_eulers(:,iparticle)),dp)
    enddo
end subroutine build_pose_fixture

subroutine assert_summary(summary,statuses,label)
    type(pcg_pose_polish_summary), intent(in) :: summary
    integer, intent(in) :: statuses(NHALF)
    character(len=*), intent(in) :: label
    integer :: accounted

    accounted = summary%nimproved+summary%nunchanged+summary%nunreliable+ &
        &summary%nstep_bound+summary%ninvalid+summary%niteration_limit
    call assert_true(accounted == summary%nparticles,trim(label)//' pose summary lost a particle')
    call assert_true(all(statuses == POSE_LM_ACCEPTED_IMPROVEMENT), &
        &trim(label)//' pose batch did not improve every perturbed particle')
    call assert_true(summary%objective_after < summary%objective_before, &
        &trim(label)//' pose batch did not reduce its aggregate objective')
    call assert_true(summary%max_rotation_step <= ROTATION_SCALE+epsilon(1._dp), &
        &trim(label)//' pose batch exceeded its rotation-step bound')
    call assert_true(summary%max_shift_step <= 1._dp+epsilon(1._dp), &
        &trim(label)//' pose batch exceeded its shift-step bound')
end subroutine assert_summary

subroutine test_symmetry_and_global_gauge(truth_rotmats)
    real(dp), intent(in) :: truth_rotmats(3,3,NPARTICLES)
    real(dp) :: symmetries(3,3,4), equivalent(3,3), gauge(3,3), gauge_estimate(3,3)
    real(dp) :: gauged(3,3,NPARTICLES), aligned(3,3), symmetry_error, gauge_error
    integer :: iparticle

    symmetries = d2_symmetries()
    equivalent = matmul(truth_rotmats(:,:,1),symmetries(:,:,2))
    symmetry_error = point_group_distance(truth_rotmats(:,:,1),equivalent,symmetries)
    call assert_true(symmetry_error < FRAME_TOL, &
        &'D2-equivalent rotations were reported as different poses')

    ! Fix one common left gauge before comparing the remaining particle poses.
    gauge = real(euler2m([4.,3.,2.]),dp)
    do iparticle = 1, NPARTICLES
        gauged(:,:,iparticle) = matmul(gauge,truth_rotmats(:,:,iparticle))
    enddo
    gauge_estimate = matmul(gauged(:,:,1),transpose(truth_rotmats(:,:,1)))
    gauge_error = 0._dp
    do iparticle = 1, NPARTICLES
        aligned = matmul(transpose(gauge_estimate),gauged(:,:,iparticle))
        gauge_error = max(gauge_error, &
            &point_group_distance(truth_rotmats(:,:,iparticle),aligned,symmetries))
    enddo
    call assert_true(gauge_error < FRAME_TOL, &
        &'global-gauge removal did not restore symmetry-aware pose agreement')
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_PCG_POSE_CONTRACT symmetry/gauge errors: ', &
        &symmetry_error,gauge_error
end subroutine test_symmetry_and_global_gauge

pure function d2_symmetries() result(symmetries)
    real(dp) :: symmetries(3,3,4)
    symmetries = 0._dp
    symmetries(1,1,:) = [1._dp,1._dp,-1._dp,-1._dp]
    symmetries(2,2,:) = [1._dp,-1._dp,1._dp,-1._dp]
    symmetries(3,3,:) = [1._dp,-1._dp,-1._dp,1._dp]
end function d2_symmetries

pure function point_group_distance(left,right,symmetries) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3), symmetries(:,:,:)
    real(dp) :: distance
    integer :: isym
    distance = huge(0._dp)
    do isym = 1, size(symmetries,3)
        distance = min(distance,rotation_distance(left,matmul(right,symmetries(:,:,isym))))
    enddo
end function point_group_distance

pure function rotation_distance(left,right) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3)
    real(dp) :: distance, cosine
    cosine = 0.5_dp*(trace3(matmul(transpose(left),right))-1._dp)
    distance = acos(max(-1._dp,min(1._dp,cosine)))
end function rotation_distance

pure function trace3(matrix) result(trace)
    real(dp), intent(in) :: matrix(3,3)
    real(dp) :: trace
    trace = matrix(1,1)+matrix(2,2)+matrix(3,3)
end function trace3

end module continuous_3D_pcg_refinement_pose_contract_test
