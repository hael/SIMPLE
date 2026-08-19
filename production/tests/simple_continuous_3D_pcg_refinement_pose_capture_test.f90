module continuous_3D_pcg_refinement_pose_capture_test
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api, only: euler2m
use simple_defs, only: dp, PI
use simple_image, only: image
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, &
    &POSE_LM_ACCEPTED_IMPROVEMENT, POSE_LM_FINITE_NO_IMPROVEMENT, &
    &POSE_LM_NO_RELIABLE_UPDATE, POSE_LM_STEP_BOUND_REJECTED, &
    &POSE_LM_INVALID_NUMERICS, POSE_LM_ITERATION_LIMIT, right_increment_rotation
use simple_string, only: string
implicit none
private
public :: run_pose_capture_range

integer, parameter :: MAX_CAPTURE_ITERATIONS = 40
integer, parameter :: NSHIFT_MAGNITUDES = 15
integer, parameter :: NROTATION_MAGNITUDES = 10
integer, parameter :: NSMALL_JOINT_CASES = 2
integer, parameter :: NJOINT_ROTATIONS = 3
integer, parameter :: NJOINT_SHIFTS = 3
integer, parameter :: NMIXED_SIGN_CASES = 3
integer, parameter :: NMULTI_AXIS_CASES = 3
integer, parameter :: NJOINT_CASES = NSMALL_JOINT_CASES+NJOINT_ROTATIONS*NJOINT_SHIFTS+ &
    &NMIXED_SIGN_CASES+NMULTI_AXIS_CASES
real, parameter :: CAPTURE_SMPD = 1.
real(dp), parameter :: ROTATION_STEP_BOUND = 0.10_dp
real(dp), parameter :: SHIFT_MAGNITUDES(NSHIFT_MAGNITUDES) = &
    &[0.25_dp,0.50_dp,1._dp,2._dp,3._dp,4._dp,5._dp,5.5_dp,6._dp,6.5_dp,7._dp, &
    &7.5_dp,10._dp,12.5_dp,15._dp]
real(dp), parameter :: ROTATION_MAGNITUDES_DEG(NROTATION_MAGNITUDES) = &
    &[0.25_dp,0.50_dp,1._dp,2._dp,3._dp,5._dp,7.5_dp,10._dp,12.5_dp,15._dp]
real(dp), parameter :: SMALL_JOINT_ROTATION_DEG(NSMALL_JOINT_CASES) = &
    &[2._dp,5._dp]
real(dp), parameter :: SMALL_JOINT_SHIFT_PX(NSMALL_JOINT_CASES) = &
    &[0.5_dp,1._dp]
real(dp), parameter :: JOINT_ROTATION_DEG(NJOINT_ROTATIONS) = &
    &[10._dp,12.5_dp,15._dp]
real(dp), parameter :: JOINT_SHIFT_PX(NJOINT_SHIFTS) = &
    &[3._dp,4._dp,5._dp]
real(dp), parameter :: MIXED_ROTATION_DEG(3,NMIXED_SIGN_CASES) = reshape( &
    &[5._dp,0._dp,0._dp, -10._dp,0._dp,0._dp, -15._dp,0._dp,0._dp], &
    &[3,NMIXED_SIGN_CASES])
real(dp), parameter :: MIXED_SHIFT_PX(2,NMIXED_SIGN_CASES) = reshape( &
    &[-1._dp,0._dp, 3._dp,0._dp, 5._dp,0._dp], [2,NMIXED_SIGN_CASES])
real(dp), parameter :: MULTI_ROTATION_DEG(3,NMULTI_AXIS_CASES) = reshape( &
    &[2._dp,-3._dp,4._dp, 5._dp,-7.5_dp,10._dp, -5._dp,7.5_dp,-10._dp], &
    &[3,NMULTI_AXIS_CASES])
real(dp), parameter :: MULTI_SHIFT_PX(2,NMULTI_AXIS_CASES) = reshape( &
    &[0.5_dp,-1._dp, 2._dp,-3._dp, -2._dp,3._dp], [2,NMULTI_AXIS_CASES])

contains

!> Characterize the matched PCG Cartesian five-parameter LM capture range.
!! Separate sweeps bracket each coordinate. Joint grids, mixed signs and
!! simultaneous multi-axis rotations then expose cross-family coupling.
subroutine run_pose_capture_range()
    type(pcg_fourier_workspace) :: workspace
    type(reconstructor_pcg) :: pcgop
    real, allocatable :: phantom(:,:,:)
    complex, allocatable :: observed(:,:), zero_plane(:,:)
    real(dp) :: truth_rotmat(3,3), truth_shift(2), ignored_objective
    character(len=4096) :: output_dir
    integer :: axis, csv_unit, env_length, env_status, integrity_failures
    integer :: joint_index, magnitude_index, rotation_index, shift_index
    integer :: sign_value, lims2(2,2)
    real(dp) :: rotation_vector_deg(3), shift_vector_px(2)

    call get_environment_variable('CONTINUOUS_3D_POSE_CAPTURE_DIR',output_dir, &
        &length=env_length,status=env_status)
    if( env_status /= 0 .or. env_length < 1 ) &
        &error stop 'CONTINUOUS_3D_POSE_CAPTURE_DIR is required'
    output_dir = output_dir(:env_length)

    call build_truth_volume(phantom)
    call pcgop%new(BOX,CAPTURE_SMPD)
    call pcgop%set_volume(phantom)
    call pcgop%begin_fourier_workspace(workspace)
    call workspace%set_shell_range([2,BOX/2])
    lims2 = workspace%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(zero_plane,mold=observed)
    zero_plane = cmplx(0.,0.)

    truth_rotmat = real(euler2m([19.,37.,28.]),dp)
    truth_shift = [0.31_dp,-0.24_dp]
    ! With y=0, the residual routine returns the executed shifted prediction.
    call workspace%shift_residual(truth_rotmat,truth_shift,zero_plane,observed,ignored_objective)
    call write_real_volume(phantom,trim(output_dir)//'/truth_volume.mrc')
    call write_fourier_image(observed,trim(output_dir)//'/truth_observation.mrc')
    call write_configuration(trim(output_dir))

    open(newunit=csv_unit,file=trim(output_dir)//'/pose_capture.csv',status='replace', &
        &action='write',form='formatted')
    write(csv_unit,'(a)') 'trial_id,family,axis,sign,injected_magnitude,'// &
        &'injected_rotation_deg,injected_shift_px,injected_rotation_x_deg,'// &
        &'injected_rotation_y_deg,injected_rotation_z_deg,injected_shift_x_px,'// &
        &'injected_shift_y_px,initial_rotation_deg,'// &
        &'final_rotation_deg,initial_shift_px,final_shift_px,objective_before,'// &
        &'objective_after,accepted_steps,attempted_steps,status,status_name,'// &
        &'max_rotation_step_deg,max_shift_step_px,stencil_switches,'// &
        &'monotone,rotation_improved,shift_improved'

    integrity_failures = 0
    call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
        &trim(output_dir),csv_unit,'exact',0,0,0,[0._dp,0._dp,0._dp], &
        &[0._dp,0._dp],integrity_failures)

    do axis = 1, 2
        do magnitude_index = 1, NSHIFT_MAGNITUDES
            do sign_value = -1, 1, 2
                rotation_vector_deg = 0._dp
                shift_vector_px = 0._dp
                shift_vector_px(axis) = real(sign_value,dp)*SHIFT_MAGNITUDES(magnitude_index)
                call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
                    &trim(output_dir),csv_unit,'shift',axis,sign_value,0,rotation_vector_deg, &
                    &shift_vector_px,integrity_failures)
            enddo
        enddo
    enddo
    do axis = 1, 3
        do magnitude_index = 1, NROTATION_MAGNITUDES
            do sign_value = -1, 1, 2
                rotation_vector_deg = 0._dp
                shift_vector_px = 0._dp
                rotation_vector_deg(axis) = real(sign_value,dp)* &
                    &ROTATION_MAGNITUDES_DEG(magnitude_index)
                call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
                    &trim(output_dir),csv_unit,'rotation',axis,sign_value,0,rotation_vector_deg, &
                    &shift_vector_px,integrity_failures)
            enddo
        enddo
    enddo
    joint_index = 0
    do magnitude_index = 1, NSMALL_JOINT_CASES
        joint_index = joint_index+1
        rotation_vector_deg = [SMALL_JOINT_ROTATION_DEG(magnitude_index),0._dp,0._dp]
        shift_vector_px = [SMALL_JOINT_SHIFT_PX(magnitude_index),0._dp]
        call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
            &trim(output_dir),csv_unit,'joint',1,1,joint_index,rotation_vector_deg, &
            &shift_vector_px,integrity_failures)
    enddo
    ! Bracket the coupled boundary on a complete positive x/x grid.
    do rotation_index = 1, NJOINT_ROTATIONS
        do shift_index = 1, NJOINT_SHIFTS
            joint_index = joint_index+1
            rotation_vector_deg = [JOINT_ROTATION_DEG(rotation_index),0._dp,0._dp]
            shift_vector_px = [JOINT_SHIFT_PX(shift_index),0._dp]
            call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
                &trim(output_dir),csv_unit,'joint',1,1,joint_index,rotation_vector_deg, &
                &shift_vector_px,integrity_failures)
        enddo
    enddo
    do magnitude_index = 1, NMIXED_SIGN_CASES
        joint_index = joint_index+1
        call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
            &trim(output_dir),csv_unit,'joint',1,0,joint_index, &
            &MIXED_ROTATION_DEG(:,magnitude_index),MIXED_SHIFT_PX(:,magnitude_index), &
            &integrity_failures)
    enddo
    do magnitude_index = 1, NMULTI_AXIS_CASES
        joint_index = joint_index+1
        call run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
            &trim(output_dir),csv_unit,'joint',0,0,joint_index, &
            &MULTI_ROTATION_DEG(:,magnitude_index),MULTI_SHIFT_PX(:,magnitude_index), &
            &integrity_failures)
    enddo
    call assert_true(joint_index == NJOINT_CASES, 'joint pose-capture matrix is incomplete')
    close(csv_unit)

    write(*,'(a,i0)') 'CONTINUOUS_3D_POSE_CAPTURE trials: ', &
        &1+2*2*NSHIFT_MAGNITUDES+2*3*NROTATION_MAGNITUDES+NJOINT_CASES
    write(*,'(a,i0)') 'CONTINUOUS_3D_POSE_CAPTURE integrity failures: ',integrity_failures
    call assert_true(integrity_failures == 0, &
        &'pose-capture diagnostic violated a numerical integrity gate')
    write(*,'(a)') 'CONTINUOUS_3D_POSE_CAPTURE: EVIDENCE COMPLETE'

    deallocate(phantom,observed,zero_plane)
    call workspace%kill()
    call pcgop%kill()
end subroutine run_pose_capture_range

!> Run one isolated injected error and retain the complete five-parameter solve.
subroutine run_capture_trial(workspace,truth_rotmat,truth_shift,observed,zero_plane, &
    &output_dir,csv_unit,family,axis,sign_value,scenario_index,rotation_vector_deg, &
    &shift_vector_px,integrity_failures)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: truth_rotmat(3,3), truth_shift(2)
    real(dp), intent(in) :: rotation_vector_deg(3), shift_vector_px(2)
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:), zero_plane(-BOX/2:,-BOX/2:)
    character(len=*), intent(in) :: output_dir, family
    integer, intent(in) :: csv_unit, axis, sign_value, scenario_index
    integer, intent(inout) :: integrity_failures
    complex, allocatable :: initial_prediction(:,:), final_prediction(:,:)
    complex, allocatable :: initial_residual(:,:), final_residual(:,:)
    real(dp) :: estimate_rotmat(3,3), estimate_shift(2), omega(3)
    real(dp) :: accepted_objectives(0:MAX_CAPTURE_ITERATIONS)
    real(dp) :: initial_rotation_error, final_rotation_error
    real(dp) :: initial_shift_error, final_shift_error
    real(dp) :: objective_before, objective_after, ignored_objective
    real(dp) :: max_rotation_step, max_shift_step
    real(dp) :: injected_magnitude, injected_rotation_deg, injected_shift_px
    character(len=64) :: trial_id
    logical :: finite_metrics, monotone, rotation_improved, shift_improved
    integer :: naccepted, nattempted, nstencil_switches, status

    allocate(initial_prediction,mold=observed)
    allocate(final_prediction,mold=observed)
    allocate(initial_residual,mold=observed)
    allocate(final_residual,mold=observed)
    omega = rotation_vector_deg*real(PI,dp)/180._dp
    estimate_rotmat = right_increment_rotation(truth_rotmat,omega)
    estimate_shift = truth_shift+shift_vector_px
    injected_magnitude = 0._dp
    injected_rotation_deg = sqrt(sum(rotation_vector_deg**2))
    injected_shift_px = sqrt(sum(shift_vector_px**2))

    select case(family)
    case('exact')
        trial_id = 'exact_pose'
    case('shift')
        injected_magnitude = injected_shift_px
        trial_id = capture_trial_id(family,axis,sign_value,scenario_index, &
            &rotation_vector_deg,shift_vector_px)
    case('rotation')
        injected_magnitude = injected_rotation_deg
        trial_id = capture_trial_id(family,axis,sign_value,scenario_index, &
            &rotation_vector_deg,shift_vector_px)
    case('joint')
        trial_id = capture_trial_id(family,axis,sign_value,scenario_index, &
            &rotation_vector_deg,shift_vector_px)
    case default
        error stop 'unknown pose-capture trial family'
    end select

    initial_rotation_error = rotation_distance(estimate_rotmat,truth_rotmat)
    initial_shift_error = sqrt(sum((estimate_shift-truth_shift)**2))
    call workspace%shift_residual(estimate_rotmat,estimate_shift,zero_plane, &
        &initial_prediction,ignored_objective)
    initial_residual = initial_prediction-observed

    call workspace%refine_pose_lm(estimate_rotmat,observed,estimate_shift, &
        &ROTATION_STEP_BOUND,MAX_CAPTURE_ITERATIONS,accepted_objectives,naccepted, &
        &status,nattempted,max_rotation_step,max_shift_step,nstencil_switches)

    final_rotation_error = rotation_distance(estimate_rotmat,truth_rotmat)
    final_shift_error = sqrt(sum((estimate_shift-truth_shift)**2))
    objective_before = accepted_objectives(0)
    call workspace%shift_residual(estimate_rotmat,estimate_shift,observed, &
        &final_residual,objective_after)
    call workspace%shift_residual(estimate_rotmat,estimate_shift,zero_plane, &
        &final_prediction,ignored_objective)

    monotone = .true.
    if( naccepted > 0 ) &
        &monotone = all(accepted_objectives(1:naccepted) < accepted_objectives(0:naccepted-1))
    rotation_improved = final_rotation_error < initial_rotation_error
    shift_improved = final_shift_error < initial_shift_error
    finite_metrics = all(ieee_is_finite([initial_rotation_error,final_rotation_error, &
        &initial_shift_error,final_shift_error,objective_before,objective_after, &
        &max_rotation_step,max_shift_step]))

    call record_integrity(family,truth_rotmat,truth_shift,estimate_rotmat,estimate_shift, &
        &finite_metrics,monotone,objective_before,objective_after,naccepted,nattempted, &
        &status,max_rotation_step,max_shift_step,integrity_failures)
    call write_capture_csv(csv_unit,trim(trial_id),family,axis,sign_value,injected_magnitude, &
        &injected_rotation_deg,injected_shift_px,rotation_vector_deg,shift_vector_px, &
        &initial_rotation_error, &
        &final_rotation_error,initial_shift_error,final_shift_error,objective_before, &
        &objective_after,naccepted,nattempted,status,max_rotation_step,max_shift_step, &
        &nstencil_switches,monotone,rotation_improved,shift_improved)

    call write_fourier_image(initial_prediction, &
        &trim(output_dir)//'/'//trim(trial_id)//'_initial_prediction.mrc')
    call write_fourier_image(initial_residual, &
        &trim(output_dir)//'/'//trim(trial_id)//'_initial_residual.mrc')
    call write_fourier_image(final_prediction, &
        &trim(output_dir)//'/'//trim(trial_id)//'_final_prediction.mrc')
    call write_fourier_image(final_residual, &
        &trim(output_dir)//'/'//trim(trial_id)//'_final_residual.mrc')

    write(*,'(a,1x,a,4(1x,es12.4),2(1x,i0),1x,a)') &
        &'CONTINUOUS_3D_POSE_CAPTURE',trim(trial_id),initial_rotation_error*180._dp/real(PI,dp), &
        &final_rotation_error*180._dp/real(PI,dp),initial_shift_error,final_shift_error, &
        &naccepted,nattempted,trim(pose_status_name(status))
    deallocate(initial_prediction,final_prediction,initial_residual,final_residual)
end subroutine run_capture_trial

!> Count only harness corruption; capture-range failures remain reportable evidence.
subroutine record_integrity(family,truth_rotmat,truth_shift,estimate_rotmat,estimate_shift, &
    &finite_metrics,monotone,objective_before,objective_after,naccepted,nattempted,status, &
    &max_rotation_step,max_shift_step,integrity_failures)
    character(len=*), intent(in) :: family
    real(dp), intent(in) :: truth_rotmat(3,3), truth_shift(2)
    real(dp), intent(in) :: estimate_rotmat(3,3), estimate_shift(2)
    real(dp), intent(in) :: objective_before, objective_after
    real(dp), intent(in) :: max_rotation_step, max_shift_step
    logical, intent(in) :: finite_metrics, monotone
    integer, intent(in) :: naccepted, nattempted, status
    integer, intent(inout) :: integrity_failures
    logical :: valid_status

    valid_status = status >= POSE_LM_ACCEPTED_IMPROVEMENT .and. &
        &status <= POSE_LM_ITERATION_LIMIT
    if( .not. finite_metrics .or. .not. monotone .or. nattempted < naccepted .or. &
        &.not. valid_status .or. status == POSE_LM_INVALID_NUMERICS .or. &
        &max_rotation_step > ROTATION_STEP_BOUND+epsilon(1._dp) .or. &
        &max_shift_step > 1._dp+epsilon(1._dp) .or. &
        &objective_after > objective_before+epsilon(1._dp)*max(1._dp,abs(objective_before)) )then
        integrity_failures = integrity_failures+1
    endif
    if( family == 'exact' )then
        if( status /= POSE_LM_FINITE_NO_IMPROVEMENT .or. naccepted /= 0 .or. &
            &maxval(abs(estimate_rotmat-truth_rotmat)) > 10._dp*epsilon(1._dp) .or. &
            &maxval(abs(estimate_shift-truth_shift)) > 10._dp*epsilon(1._dp) ) &
            &integrity_failures = integrity_failures+1
    endif
end subroutine record_integrity

!> Write one machine-readable result row without assigning a capture threshold.
subroutine write_capture_csv(csv_unit,trial_id,family,axis,sign_value,magnitude, &
    &injected_rotation_deg,injected_shift_px,rotation_vector_deg,shift_vector_px, &
    &initial_rotation_error,final_rotation_error,initial_shift_error,final_shift_error, &
    &objective_before,objective_after,naccepted,nattempted,status,max_rotation_step, &
    &max_shift_step,nstencil_switches,monotone,rotation_improved,shift_improved)
    integer, intent(in) :: csv_unit, axis, sign_value, naccepted, nattempted, status
    integer, intent(in) :: nstencil_switches
    character(len=*), intent(in) :: trial_id, family
    real(dp), intent(in) :: magnitude, injected_rotation_deg, injected_shift_px
    real(dp), intent(in) :: rotation_vector_deg(3), shift_vector_px(2)
    real(dp), intent(in) :: initial_rotation_error, final_rotation_error
    real(dp), intent(in) :: initial_shift_error, final_shift_error
    real(dp), intent(in) :: objective_before, objective_after
    real(dp), intent(in) :: max_rotation_step, max_shift_step
    logical, intent(in) :: monotone, rotation_improved, shift_improved

    write(csv_unit,'(*(g0,:,","))') trim(trial_id),trim(family),axis,sign_value,magnitude, &
        &injected_rotation_deg,injected_shift_px,rotation_vector_deg,shift_vector_px, &
        &initial_rotation_error*180._dp/real(PI,dp), &
        &final_rotation_error*180._dp/real(PI,dp),initial_shift_error,final_shift_error, &
        &objective_before,objective_after,naccepted,nattempted,status, &
        &trim(pose_status_name(status)),max_rotation_step*180._dp/real(PI,dp), &
        &max_shift_step,nstencil_switches,merge(1,0,monotone), &
        &merge(1,0,rotation_improved),merge(1,0,shift_improved)
end subroutine write_capture_csv

!> Create a stable artifact name from one signed perturbation magnitude.
function capture_trial_id(family,axis,sign_value,scenario_index,rotation_vector_deg, &
    &shift_vector_px) result(trial_id)
    character(len=*), intent(in) :: family
    integer, intent(in) :: axis, sign_value, scenario_index
    real(dp), intent(in) :: rotation_vector_deg(3), shift_vector_px(2)
    character(len=64) :: trial_id
    character(len=6) :: rotation_digits, shift_digits
    character(len=3) :: index_digits
    character(len=1) :: axis_name, sign_name
    real(dp) :: rotation_magnitude_deg, shift_magnitude_px

    axis_name = 'x'
    if( axis >= 1 .and. axis <= 3 ) axis_name = achar(iachar('x')+axis-1)
    sign_name = merge('p','m',sign_value > 0)
    rotation_magnitude_deg = sqrt(sum(rotation_vector_deg**2))
    shift_magnitude_px = sqrt(sum(shift_vector_px**2))
    write(rotation_digits,'(i6.6)') nint(1000._dp*rotation_magnitude_deg)
    write(shift_digits,'(i6.6)') nint(1000._dp*shift_magnitude_px)
    write(index_digits,'(i3.3)') scenario_index
    if( family == 'rotation' )then
        trial_id = 'rotation_'//axis_name//'_'//sign_name//rotation_digits//'mdeg'
    else if( family == 'joint' )then
        if( axis == 1 .and. sign_value == 1 )then
            trial_id = 'joint_x_p'//rotation_digits//'mdeg_p'//shift_digits//'mpix'
        else if( axis == 1 )then
            trial_id = 'joint_mixed_'//index_digits
        else
            trial_id = 'joint_multi_'//index_digits
        endif
    else
        trial_id = 'shift_'//axis_name//'_'//sign_name//shift_digits//'mpix'
    endif
end function capture_trial_id

!> Return the optimizer terminal status as a stable report label.
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

!> Write the fixed input volume used by every capture trial.
subroutine write_real_volume(volume,filename)
    real, intent(in) :: volume(BOX,BOX,BOX)
    character(len=*), intent(in) :: filename
    type(image) :: volume_image

    call volume_image%new([BOX,BOX,BOX],CAPTURE_SMPD,wthreads=.false.)
    call volume_image%set_rmat(volume,.false.)
    call volume_image%write(string(filename),del_if_exists=.true.)
    call volume_image%kill()
end subroutine write_real_volume

!> Record the fixed numerical conditions needed to interpret the sweep.
subroutine write_configuration(output_dir)
    character(len=*), intent(in) :: output_dir
    integer :: unit

    open(newunit=unit,file=trim(output_dir)//'/configuration.txt',status='replace', &
        &action='write',form='formatted')
    write(unit,'(a,i0)') 'box=',BOX
    write(unit,'(a,f0.6)') 'sampling_distance_A=',CAPTURE_SMPD
    write(unit,'(a,i0)') 'shell_min=',2
    write(unit,'(a,i0)') 'shell_max=',BOX/2
    write(unit,'(a,f0.6)') 'nominal_lowpass_A=',real(BOX,dp)*real(CAPTURE_SMPD,dp)/real(BOX/2,dp)
    write(unit,'(a,f0.6)') 'rotation_step_bound_deg=', &
        &ROTATION_STEP_BOUND*180._dp/real(PI,dp)
    write(unit,'(a,f0.6)') 'shift_step_bound_px=',1._dp
    write(unit,'(a,i0)') 'max_lm_iterations=',MAX_CAPTURE_ITERATIONS
    write(unit,'(a)') 'truth_euler_deg=19,37,28'
    write(unit,'(a)') 'truth_shift_px=0.31,-0.24'
    write(unit,'(a)') 'observation_model=matched_clean_pcg_cartesian'
    write(unit,'(a)') 'scenario_version=4'
    write(unit,'(a)') 'shift_magnitudes_px=0.25,0.5,1,2,3,4,5,5.5,6,6.5,7,7.5,10,12.5,15'
    write(unit,'(a)') 'joint_boundary_rotation_deg=10,12.5,15'
    write(unit,'(a)') 'joint_boundary_shift_px=3,4,5'
    write(unit,'(a)') 'perturbation_policy=separate_boundary_joint_grid_mixed_sign_multi_axis'
    write(unit,'(a)') 'mechanism_volumes=gaussian_blobs,shifted_mixture,permuted_texture'
    write(unit,'(a)') 'mechanism_routes=joint,guarded_joint,shift_then_joint,rotation_then_joint'
    write(unit,'(a)') 'cumulative_guard_from_seed=15deg,5px'
    close(unit)
end subroutine write_configuration

!> Return the chordal form of the rotation geodesic, stable at equal matrices.
pure function rotation_distance(left,right) result(distance)
    real(dp), intent(in) :: left(3,3), right(3,3)
    real(dp) :: distance, sine_half

    sine_half = sqrt(sum((left-right)**2))/(2._dp*sqrt(2._dp))
    distance = 2._dp*asin(max(0._dp,min(1._dp,sine_half)))
end function rotation_distance

end module continuous_3D_pcg_refinement_pose_capture_test
