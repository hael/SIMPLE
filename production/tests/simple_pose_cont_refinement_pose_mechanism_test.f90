module pose_cont_refinement_pose_mechanism_test
use pose_cont_refinement_test_helpers, only: build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use simple_core_module_api, only: euler2m
use simple_defs, only: dp, PI
use simple_image, only: image
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, right_increment_rotation, &
    &POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT, &
    &POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED, &
    &POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_string, only: string
implicit none
private
public :: run_pose_capture_mechanism

integer, parameter :: MAX_ITERATIONS = 40
integer, parameter :: NVOLUMES = 3
integer, parameter :: NCASES = 5
integer, parameter :: NPATH_POINTS = 41
real, parameter :: CAPTURE_SMPD = 1.
real(dp), parameter :: STEP_ROTATION_BOUND = 0.10_dp
real(dp), parameter :: TOTAL_ROTATION_BOUND = 15._dp*real(PI,dp)/180._dp
real(dp), parameter :: TOTAL_SHIFT_BOUND = 5._dp
real(dp), parameter :: CASE_ROTATION_DEG(NCASES) = [10._dp,12.5_dp,15._dp,15._dp,-15._dp]
real(dp), parameter :: CASE_SHIFT_PX(NCASES) = [5._dp,5._dp,5._dp,4._dp,5._dp]
character(len=24), parameter :: VOLUME_NAMES(NVOLUMES) = [character(len=24) :: &
    &'gaussian_blobs', 'shifted_mixture', 'permuted_texture']
character(len=24), parameter :: CASE_NAMES(NCASES) = [character(len=24) :: &
    &'joint_p10_p5', 'joint_p12p5_p5', 'joint_p15_p5', 'joint_p15_p4', 'joint_m15_p5']

contains

!> Diagnose the local basin with multiple volumes, trajectories and staged solves.
subroutine run_pose_capture_mechanism()
    type(cartesian_pose_refiner) :: workspace
    type(reconstructor_pcg) :: pcgop
    real, allocatable :: volume(:,:,:)
    complex, allocatable :: observed(:,:), zero_plane(:,:), residual(:,:)
    real(dp) :: truth_rotmat(3,3), truth_shift(2), seed_rotmat(3,3), seed_shift(2)
    real(dp) :: direct_rotmat(3,3), direct_shift(2), ignored_objective
    real(dp) :: rotation_vector(3)
    character(len=4096) :: output_dir
    integer :: case_index, env_length, env_status, lims2(2,2)
    integer :: path_unit, summary_unit, trajectory_unit, volume_index

    call get_environment_variable('CONTINUOUS_3D_POSE_CAPTURE_DIR',output_dir, &
        &length=env_length,status=env_status)
    if( env_status /= 0 .or. env_length < 1 ) &
        &error stop 'CONTINUOUS_3D_POSE_CAPTURE_DIR is required'
    output_dir = output_dir(:env_length)
    call open_evidence_files(trim(output_dir),trajectory_unit,summary_unit,path_unit)

    truth_rotmat = real(euler2m([19.,37.,28.]),dp)
    truth_shift = [0.31_dp,-0.24_dp]
    do volume_index = 1, NVOLUMES
        call build_mechanism_volume(volume_index,volume)
        call pcgop%new(BOX,CAPTURE_SMPD)
        call pcgop%set_volume(volume)
        call workspace%new_prepared_test(volume)
        call workspace%set_shell_range([2,BOX/2])
        lims2 = workspace%get_lims2()
        allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
        allocate(zero_plane,mold=observed)
        allocate(residual,mold=observed)
        zero_plane = cmplx(0.,0.)
        call workspace%shift_residual(truth_rotmat,truth_shift,zero_plane,observed,ignored_objective)
        call write_real_volume(volume,trim(output_dir)//'/mechanism_'// &
            &trim(VOLUME_NAMES(volume_index))//'_volume.mrc')
        call write_fourier_image(observed,trim(output_dir)//'/mechanism_'// &
            &trim(VOLUME_NAMES(volume_index))//'_observation.mrc')

        do case_index = 1, NCASES
            rotation_vector = [CASE_ROTATION_DEG(case_index),0._dp,0._dp]*real(PI,dp)/180._dp
            seed_rotmat = right_increment_rotation(truth_rotmat,rotation_vector)
            seed_shift = truth_shift+[CASE_SHIFT_PX(case_index),0._dp]

            call run_route(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
                &trim(VOLUME_NAMES(volume_index)),trim(CASE_NAMES(case_index)),'joint', &
                &trajectory_unit,summary_unit,direct_rotmat,direct_shift)
            call run_route(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
                &trim(VOLUME_NAMES(volume_index)),trim(CASE_NAMES(case_index)),'guarded_joint', &
                &trajectory_unit,summary_unit)
            call run_route(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
                &trim(VOLUME_NAMES(volume_index)),trim(CASE_NAMES(case_index)),'shift_then_joint', &
                &trajectory_unit,summary_unit)
            call run_route(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
                &trim(VOLUME_NAMES(volume_index)),trim(CASE_NAMES(case_index)),'rotation_then_joint', &
                &trajectory_unit,summary_unit)
            call write_objective_path(workspace,observed,residual,truth_rotmat,truth_shift, &
                &seed_rotmat,seed_shift,truth_rotmat,truth_shift,trim(VOLUME_NAMES(volume_index)), &
                &trim(CASE_NAMES(case_index)),'truth',path_unit)
            call write_objective_path(workspace,observed,residual,truth_rotmat,truth_shift, &
                &seed_rotmat,seed_shift,direct_rotmat,direct_shift,trim(VOLUME_NAMES(volume_index)), &
                &trim(CASE_NAMES(case_index)),'joint_endpoint',path_unit)
        enddo

        deallocate(volume,observed,zero_plane,residual)
        call workspace%kill()
        call pcgop%kill()
    enddo
    close(trajectory_unit)
    close(summary_unit)
    close(path_unit)
    write(*,'(a,i0)') 'CONTINUOUS_3D_POSE_MECHANISM volumes: ',NVOLUMES
    write(*,'(a,i0)') 'CONTINUOUS_3D_POSE_MECHANISM seed cases per volume: ',NCASES
    write(*,'(a)') 'CONTINUOUS_3D_POSE_MECHANISM: EVIDENCE COMPLETE'
end subroutine run_pose_capture_mechanism

!> Run one direct or staged route and preserve every accepted pose.
subroutine run_route(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
    &volume_name,case_name,route_name,trajectory_unit,summary_unit,final_rotmat_out,final_shift_out)
    type(cartesian_pose_refiner), intent(in) :: workspace
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    real(dp), intent(in) :: truth_rotmat(3,3), truth_shift(2), seed_rotmat(3,3), seed_shift(2)
    character(len=*), intent(in) :: volume_name, case_name, route_name
    integer, intent(in) :: trajectory_unit, summary_unit
    real(dp), optional, intent(out) :: final_rotmat_out(3,3), final_shift_out(2)
    real(dp) :: current_rotmat(3,3), current_shift(2), objective_before, objective_after
    real(dp) :: stage_before, stage_after, rotation_vector(3), shift_vector(2)
    integer :: accepted_total, attempted_total, stage_accepted, stage_attempted, status
    logical :: rotation_mask(5), shift_mask(5), joint_mask(5)

    rotation_mask = [.true.,.true.,.true.,.false.,.false.]
    shift_mask = [.false.,.false.,.false.,.true.,.true.]
    joint_mask = .true.
    current_rotmat = seed_rotmat
    current_shift = seed_shift
    accepted_total = 0
    attempted_total = 0
    objective_before = -1._dp

    select case(route_name)
    case('joint')
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'joint',joint_mask,.false.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_before = stage_before
        objective_after = stage_after
        accepted_total = stage_accepted
        attempted_total = stage_attempted
    case('guarded_joint')
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'joint',joint_mask,.true.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_before = stage_before
        objective_after = stage_after
        accepted_total = stage_accepted
        attempted_total = stage_attempted
    case('shift_then_joint')
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'shift_only',shift_mask,.false.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_before = stage_before
        accepted_total = stage_accepted
        attempted_total = stage_attempted
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'joint',joint_mask,.false.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_after = stage_after
        accepted_total = accepted_total+stage_accepted
        attempted_total = attempted_total+stage_attempted
    case('rotation_then_joint')
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'rotation_only',rotation_mask,.false.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_before = stage_before
        accepted_total = stage_accepted
        attempted_total = stage_attempted
        call run_stage(workspace,observed,truth_rotmat,truth_shift,seed_rotmat,seed_shift, &
            &volume_name,case_name,route_name,'joint',joint_mask,.false.,trajectory_unit, &
            &current_rotmat,current_shift,stage_before,stage_after,stage_accepted,stage_attempted,status)
        objective_after = stage_after
        accepted_total = accepted_total+stage_accepted
        attempted_total = attempted_total+stage_attempted
    case default
        error stop 'unknown pose-mechanism route'
    end select

    rotation_vector = rotation_log_vector(truth_rotmat,current_rotmat)*180._dp/real(PI,dp)
    shift_vector = current_shift-truth_shift
    write(summary_unit,'(*(g0,:,","))') trim(volume_name),trim(case_name),trim(route_name), &
        &rotation_vector,shift_vector,sqrt(sum(rotation_vector**2)),sqrt(sum(shift_vector**2)), &
        &objective_before,objective_after,accepted_total,attempted_total,trim(pose_status_name(status)), &
        &merge(1,0,rotation_distance(seed_rotmat,current_rotmat) <= TOTAL_ROTATION_BOUND+1.e-10_dp .and. &
        &sqrt(sum((current_shift-seed_shift)**2)) <= TOTAL_SHIFT_BOUND+1.e-10_dp)
    if( present(final_rotmat_out) ) final_rotmat_out = current_rotmat
    if( present(final_shift_out) ) final_shift_out = current_shift
end subroutine run_route

!> Run one masked LM stage and write its accepted five-parameter trajectory.
subroutine run_stage(workspace,observed,truth_rotmat,truth_shift,anchor_rotmat,anchor_shift, &
    &volume_name,case_name,route_name,stage_name,active,guarded,trajectory_unit, &
    &rotmat,shift,objective_before,objective_after,naccepted,nattempted,status)
    type(cartesian_pose_refiner), intent(in) :: workspace
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    real(dp), intent(in) :: truth_rotmat(3,3), truth_shift(2), anchor_rotmat(3,3), anchor_shift(2)
    character(len=*), intent(in) :: volume_name, case_name, route_name, stage_name
    logical, intent(in) :: active(5), guarded
    integer, intent(in) :: trajectory_unit
    real(dp), intent(inout) :: rotmat(3,3), shift(2)
    real(dp), intent(out) :: objective_before, objective_after
    integer, intent(out) :: naccepted, nattempted, status
    real(dp) :: objectives(0:MAX_ITERATIONS), rotmats(3,3,0:MAX_ITERATIONS)
    real(dp) :: shifts(2,0:MAX_ITERATIONS), rotation_vector(3), shift_vector(2)
    real(dp) :: max_rotation_step, max_shift_step
    integer :: accepted_index, nswitches

    if( guarded )then
        call workspace%refine_pose_lm(rotmat,observed,shift,STEP_ROTATION_BOUND,MAX_ITERATIONS, &
            &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nswitches, &
            &accepted_rotmats=rotmats,accepted_shifts=shifts,active_parameters=active, &
            &anchor_rotmat=anchor_rotmat,anchor_shift=anchor_shift, &
            &max_total_rotation=TOTAL_ROTATION_BOUND,max_total_shift=TOTAL_SHIFT_BOUND)
    else
        call workspace%refine_pose_lm(rotmat,observed,shift,STEP_ROTATION_BOUND,MAX_ITERATIONS, &
            &objectives,naccepted,status,nattempted,max_rotation_step,max_shift_step,nswitches, &
            &accepted_rotmats=rotmats,accepted_shifts=shifts,active_parameters=active)
    endif
    objective_before = objectives(0)
    objective_after = objectives(naccepted)
    do accepted_index = 0, naccepted
        rotation_vector = rotation_log_vector(truth_rotmat,rotmats(:,:,accepted_index))* &
            &180._dp/real(PI,dp)
        shift_vector = shifts(:,accepted_index)-truth_shift
        write(trajectory_unit,'(*(g0,:,","))') trim(volume_name),trim(case_name), &
            &trim(route_name),trim(stage_name),accepted_index,rotation_vector,shift_vector, &
            &sqrt(sum(rotation_vector**2)),sqrt(sum(shift_vector**2)),objectives(accepted_index)
    enddo
end subroutine run_stage

!> Sample a straight tangent/shift path from the seed to one target pose.
subroutine write_objective_path(workspace,observed,residual,truth_rotmat,truth_shift, &
    &seed_rotmat,seed_shift,target_rotmat,target_shift,volume_name,case_name,target_name,path_unit)
    type(cartesian_pose_refiner), intent(in) :: workspace
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    complex, intent(inout) :: residual(-BOX/2:,-BOX/2:)
    real(dp), intent(in) :: truth_rotmat(3,3), truth_shift(2), seed_rotmat(3,3), seed_shift(2)
    real(dp), intent(in) :: target_rotmat(3,3), target_shift(2)
    character(len=*), intent(in) :: volume_name, case_name, target_name
    integer, intent(in) :: path_unit
    real(dp) :: alpha, objective, path_omega(3), path_rotmat(3,3), path_shift(2)
    real(dp) :: rotation_vector(3), shift_vector(2)
    integer :: point

    path_omega = rotation_log_vector(seed_rotmat,target_rotmat)
    do point = 0, NPATH_POINTS-1
        alpha = real(point,dp)/real(NPATH_POINTS-1,dp)
        path_rotmat = right_increment_rotation(seed_rotmat,alpha*path_omega)
        path_shift = seed_shift+alpha*(target_shift-seed_shift)
        call workspace%shift_residual(path_rotmat,path_shift,observed,residual,objective)
        rotation_vector = rotation_log_vector(truth_rotmat,path_rotmat)*180._dp/real(PI,dp)
        shift_vector = path_shift-truth_shift
        write(path_unit,'(*(g0,:,","))') trim(volume_name),trim(case_name),trim(target_name), &
            &alpha,rotation_vector,shift_vector,sqrt(sum(rotation_vector**2)), &
            &sqrt(sum(shift_vector**2)),objective
    enddo
end subroutine write_objective_path

!> Construct three deterministic asymmetric morphologies at one fixed box size.
subroutine build_mechanism_volume(volume_index,volume)
    integer, intent(in) :: volume_index
    real, allocatable, intent(out) :: volume(:,:,:)
    real, allocatable :: base(:,:,:), shifted(:,:,:)
    integer :: i, j, k

    call build_truth_volume(base)
    allocate(volume(BOX,BOX,BOX))
    select case(volume_index)
    case(1)
        volume = base
    case(2)
        shifted = cshift(cshift(cshift(base,3,dim=1),-2,dim=2),1,dim=3)
        volume = base+0.45*shifted
        shifted = cshift(cshift(cshift(base,-4,dim=1),3,dim=2),-2,dim=3)
        volume = volume+0.25*shifted
        deallocate(shifted)
    case(3)
        do k = 1, BOX
            do j = 1, BOX
                do i = 1, BOX
                    volume(i,j,k) = base(i,j,k)**1.5+0.35*base(k,i,j)+ &
                        &0.20*base(BOX-j+1,k,i)
                enddo
            enddo
        enddo
    case default
        error stop 'unknown pose-mechanism volume'
    end select
    volume = volume/maxval(volume)
    deallocate(base)
end subroutine build_mechanism_volume

!> Return the right tangent vector that maps the first rotation to the second.
pure function rotation_log_vector(from_rotmat,to_rotmat) result(omega)
    real(dp), intent(in) :: from_rotmat(3,3), to_rotmat(3,3)
    real(dp) :: omega(3), relative(3,3), cosine, theta, scale

    relative = matmul(transpose(from_rotmat),to_rotmat)
    cosine = max(-1._dp,min(1._dp,(relative(1,1)+relative(2,2)+relative(3,3)-1._dp)/2._dp))
    theta = acos(cosine)
    if( theta < 1.e-8_dp )then
        scale = 0.5_dp
    else
        scale = theta/(2._dp*sin(theta))
    endif
    omega = scale*[relative(3,2)-relative(2,3),relative(1,3)-relative(3,1), &
        &relative(2,1)-relative(1,2)]
end function rotation_log_vector

!> Return the geodesic rotation distance in radians.
pure function rotation_distance(left,right) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3)
    real(dp) :: distance, sine_half

    sine_half = sqrt(sum((left-right)**2))/(2._dp*sqrt(2._dp))
    distance = 2._dp*asin(max(0._dp,min(1._dp,sine_half)))
end function rotation_distance

!> Open stable CSV schemas for trajectory, endpoint and objective-path evidence.
subroutine open_evidence_files(output_dir,trajectory_unit,summary_unit,path_unit)
    character(len=*), intent(in) :: output_dir
    integer, intent(out) :: trajectory_unit, summary_unit, path_unit

    open(newunit=trajectory_unit,file=trim(output_dir)//'/pose_mechanism_trajectory.csv', &
        &status='replace',action='write',form='formatted')
    write(trajectory_unit,'(a)') 'volume,case,route,stage,accepted_step,rotation_x_deg,'// &
        &'rotation_y_deg,rotation_z_deg,shift_x_px,shift_y_px,rotation_norm_deg,'// &
        &'shift_norm_px,objective'
    open(newunit=summary_unit,file=trim(output_dir)//'/pose_mechanism_summary.csv', &
        &status='replace',action='write',form='formatted')
    write(summary_unit,'(a)') 'volume,case,route,final_rotation_x_deg,final_rotation_y_deg,'// &
        &'final_rotation_z_deg,final_shift_x_px,final_shift_y_px,final_rotation_norm_deg,'// &
        &'final_shift_norm_px,objective_before,objective_after,accepted_steps,attempted_steps,'// &
        &'status_name,within_15deg_5px_from_seed'
    open(newunit=path_unit,file=trim(output_dir)//'/pose_objective_paths.csv', &
        &status='replace',action='write',form='formatted')
    write(path_unit,'(a)') 'volume,case,target,alpha,rotation_x_deg,rotation_y_deg,'// &
        &'rotation_z_deg,shift_x_px,shift_y_px,rotation_norm_deg,shift_norm_px,objective'
end subroutine open_evidence_files

!> Write one Hermitian Fourier plane as a real-space MRC image.
subroutine write_fourier_image(plane,filename)
    complex, intent(in) :: plane(-BOX/2:,-BOX/2:)
    character(len=*), intent(in) :: filename
    type(image) :: plane_image
    integer :: h, k, phys(3)

    call plane_image%new([BOX,BOX,1],CAPTURE_SMPD,wthreads=.false.)
    call plane_image%zero_and_flag_ft()
    do k = -BOX/2, BOX/2-1
        do h = 0, BOX/2
            if( h*h+k*k > (BOX/2)**2 ) cycle
            phys = [h+1,k+1+merge(BOX,0,k < 0),1]
            call plane_image%set_cmat_at(phys(1),phys(2),phys(3),plane(h,k))
        enddo
    enddo
    call plane_image%ifft()
    call plane_image%write(string(filename),del_if_exists=.true.)
    call plane_image%kill()
end subroutine write_fourier_image

!> Write one deterministic mechanism-test volume.
subroutine write_real_volume(volume,filename)
    real, intent(in) :: volume(BOX,BOX,BOX)
    character(len=*), intent(in) :: filename
    type(image) :: volume_image

    call volume_image%new([BOX,BOX,BOX],CAPTURE_SMPD,wthreads=.false.)
    call volume_image%set_rmat(volume,.false.)
    call volume_image%write(string(filename),del_if_exists=.true.)
    call volume_image%kill()
end subroutine write_real_volume

!> Return the LM terminal status as a stable report label.
pure function pose_status_name(status) result(name)
    integer, intent(in) :: status
    character(len=24) :: name

    select case(status)
    case(POSE_LM_ACCEPTED_IMPROVEMENT)
        name = 'accepted_improvement'
    case(POSE_LM_FINITE_NO_IMPROVEMENT)
        name = 'finite_no_improvement'
    case(POSE_LM_NO_RELIABLE_UPDATE)
        name = 'no_reliable_update'
    case(POSE_LM_STEP_BOUND_REJECTED)
        name = 'step_bound_rejected'
    case(POSE_LM_INVALID_NUMERICS)
        name = 'invalid_numerics'
    case(POSE_LM_ITERATION_LIMIT)
        name = 'iteration_limit'
    case default
        name = 'unknown'
    end select
end function pose_status_name

end module pose_cont_refinement_pose_mechanism_test
