module pose_cont_refinement_matched_window_test
use pose_cont_refinement_forward_path_test, only: build_forward_path_observations
use pose_cont_refinement_reference_support, only: build_disjoint_half_orientations, &
    &HALFSET_LAMBDA, HALFSET_SMPD
use pose_cont_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, sp, PI, OSMPL_PAD_FAC
use simple_image, only: image
use simple_ori, only: ori
use simple_oris, only: oris
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, right_increment_rotation
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_matched_window_diagnostic
public :: matched_reference_metrics, evaluate_matched_reference_batch
public :: gauge_corrected_errors

integer, parameter :: NPARTICLES = 48
integer, parameter :: NCOORDS = 5
integer, parameter :: MAX_DIAGNOSTIC_ITERATIONS = 8
integer, parameter :: NFD = 3
integer, parameter :: PADDED_BOX = OSMPL_PAD_FAC*BOX
integer, parameter :: PADDED_RADIUS = PADDED_BOX/2
real(dp), parameter :: FD_STEPS(NFD) = [2.e-4_dp,1.e-4_dp,5.e-5_dp]
real(dp), parameter :: POSE_DIRECTION(NCOORDS) = [0.30_dp,-0.40_dp,0.50_dp,0.20_dp,-0.10_dp]
real(dp), parameter :: GRADIENT_FD_TOLERANCE = 5.e-3_dp
real(dp), parameter :: EXACT_ROTATION_TOLERANCE = 0.001745_dp
real(dp), parameter :: EXACT_SHIFT_TOLERANCE = 0.05_dp
real(dp), parameter :: MAX_ROTATION_STEP = 0.05_dp
real(dp), parameter :: MAX_SHIFT_STEP = 0.25_dp
real(dp), parameter :: NUMERIC_FLOOR = epsilon(1._dp)**2
real(dp), parameter :: DIAGNOSTIC_ROTATION_OFFSET(3) = [0.010_dp,-0.008_dp,0.006_dp]
real(dp), parameter :: DIAGNOSTIC_SHIFT_OFFSET(2) = [0.20_dp,-0.15_dp]

type :: matched_reference_metrics
    real(dp) :: reference_truth_correlation = 0._dp
    real(dp) :: fitted_residual = 0._dp
    real(dp) :: amplitude_scale = 0._dp
    real(dp) :: exact_objective = 0._dp
    real(dp) :: rotation_gradient_rms = 0._dp
    real(dp) :: shift_gradient_rms = 0._dp
    real(dp) :: exact_rotation_rms = 0._dp
    real(dp) :: exact_shift_rms = 0._dp
    real(dp) :: gauge_rotation_rms = 0._dp
    real(dp) :: gauge_shift_rms = 0._dp
    real(dp) :: recovery_rotation_rms = 0._dp
    real(dp) :: recovery_shift_rms = 0._dp
    real(dp) :: exact_objective_before = 0._dp
    real(dp) :: exact_objective_after = 0._dp
    real(dp) :: recovery_objective_before = 0._dp
    real(dp) :: recovery_objective_after = 0._dp
    integer :: exact_accepted = 0
    integer :: recovery_accepted = 0
end type matched_reference_metrics

contains

!> Test an explicit simulator clip transform on PCG predictions and pose derivatives.
subroutine run_matched_window_diagnostic()
    type(image) :: padded_model, native_model
    type(image) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    type(ori) :: orientation
    type(oris) :: even_oris, odd_oris
    type(reconstructor_pcg) :: model_operator
    type(cartesian_pose_refiner) :: workspace
    complex, allocatable :: simulator_stages(:,:,:,:), observed(:,:,:)
    complex, allocatable :: model_planes(:,:,:), jacobian_planes(:,:,:,:)
    real, allocatable :: truth_volume(:,:,:)
    real(dp), allocatable :: exact_rotmats(:,:,:)
    integer, allocatable :: even_ids(:), odd_ids(:)
    real(dp) :: amplitude_scale, fitted_residual, exact_objective
    real(dp) :: exact_rotation_gradient, exact_shift_gradient
    real(dp) :: exact_final_rotmat(3,3), exact_final_shift(2)
    real(dp) :: perturbed_rotmat(3,3), perturbed_shift(2)
    real(dp) :: recovered_rotmat(3,3), recovered_shift(2)
    real(dp) :: exact_before, exact_after, recovery_before, recovery_after
    real(dp) :: exact_rotation_error, exact_shift_error
    real(dp) :: initial_rotation_error, initial_shift_error
    real(dp) :: recovered_rotation_error, recovered_shift_error
    real(dp) :: fd_errors(NFD)
    integer :: axis, exact_accepted, i, recovery_accepted

    call build_disjoint_half_orientations(NPARTICLES,even_oris,odd_oris,even_ids,odd_ids)
    call build_truth_volume(truth_volume)
    call model_operator%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call model_operator%set_volume(truth_volume)
    call workspace%new_prepared_test(truth_volume)
    call workspace%set_shell_range([2,BOX/2])
    call initialize_window_images(padded_model,native_model,padded_jacobian,native_jacobian)

    allocate(exact_rotmats(3,3,NPARTICLES))
    call orientation%new(.false.)
    do i = 1, NPARTICLES
        call even_oris%get_ori(i,orientation)
        exact_rotmats(:,:,i) = real(orientation%get_mat(),dp)
    enddo
    call build_forward_path_observations(model_operator,even_oris,truth_volume,simulator_stages)
    allocate(observed,mold=simulator_stages(:,:,:,3))
    observed = simulator_stages(:,:,:,3)
    allocate(model_planes(-BOX/2:BOX/2,-BOX/2:BOX/2,NPARTICLES),source=cmplx(0.,0.))
    allocate(jacobian_planes(-BOX/2:BOX/2,-BOX/2:BOX/2,NCOORDS,NPARTICLES), &
        &source=cmplx(0.,0.))

    ! Apply the identical linear clip transform to each PCG value and Jacobian plane.
    do i = 1, NPARTICLES
        call evaluate_windowed_prediction(workspace,model_operator,exact_rotmats(:,:,i), &
            &[0._dp,0._dp],padded_model,native_model,padded_jacobian,native_jacobian, &
            &model_planes(:,:,i),jacobian_planes(:,:,:,i))
    enddo
    call fit_global_amplitude(model_planes,observed,amplitude_scale,fitted_residual)
    call measure_exact_batch(model_planes,jacobian_planes,observed,amplitude_scale, &
        &exact_objective,exact_rotation_gradient,exact_shift_gradient)

    perturbed_rotmat = right_increment_rotation(exact_rotmats(:,:,1),DIAGNOSTIC_ROTATION_OFFSET)
    perturbed_shift = DIAGNOSTIC_SHIFT_OFFSET
    call measure_windowed_gradient_fd(workspace,model_operator,perturbed_rotmat,perturbed_shift, &
        &observed(:,:,1),amplitude_scale,padded_model,native_model,padded_jacobian, &
        &native_jacobian,fd_errors)

    exact_final_rotmat = exact_rotmats(:,:,1)
    exact_final_shift = 0._dp
    call refine_windowed_pose(workspace,model_operator,observed(:,:,1),amplitude_scale, &
        &padded_model,native_model,padded_jacobian,native_jacobian,exact_final_rotmat, &
        &exact_final_shift,exact_before,exact_after,exact_accepted)
    exact_rotation_error = rotation_error(exact_rotmats(:,:,1),exact_final_rotmat)
    exact_shift_error = sqrt(sum(exact_final_shift**2))

    recovered_rotmat = perturbed_rotmat
    recovered_shift = perturbed_shift
    initial_rotation_error = rotation_error(exact_rotmats(:,:,1),recovered_rotmat)
    initial_shift_error = sqrt(sum(recovered_shift**2))
    call refine_windowed_pose(workspace,model_operator,observed(:,:,1),amplitude_scale, &
        &padded_model,native_model,padded_jacobian,native_jacobian,recovered_rotmat, &
        &recovered_shift,recovery_before,recovery_after,recovery_accepted)
    recovered_rotation_error = rotation_error(exact_rotmats(:,:,1),recovered_rotmat)
    recovered_shift_error = sqrt(sum(recovered_shift**2))

    write(*,'(a,2(1x,es14.6))') &
        &'CONTINUOUS_3D_MATCHED_WINDOW fitted residual/scale',fitted_residual,amplitude_scale
    write(*,'(a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_MATCHED_WINDOW exact objective/rotation-gradient/shift-gradient', &
        &exact_objective,exact_rotation_gradient,exact_shift_gradient
    write(*,'(a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_MATCHED_WINDOW gradient finite-difference errors',fd_errors
    write(*,'(a,4(1x,es14.6),1x,i0)') &
        &'CONTINUOUS_3D_MATCHED_WINDOW exact rot/shift/objective before/after/accepted', &
        &exact_rotation_error,exact_shift_error,exact_before,exact_after,exact_accepted
    write(*,'(a,6(1x,es14.6),1x,i0)') &
        &'CONTINUOUS_3D_MATCHED_WINDOW recovery initial/final rot/shift objective before/after/accepted', &
        &initial_rotation_error,initial_shift_error,recovered_rotation_error, &
        &recovered_shift_error,recovery_before,recovery_after,recovery_accepted
    call print_matched_window_decision(fitted_residual,exact_rotation_error,exact_shift_error, &
        &initial_rotation_error,initial_shift_error,recovered_rotation_error,recovered_shift_error)

    call assert_true(all(ieee_is_finite([amplitude_scale,fitted_residual,exact_objective, &
        &exact_rotation_gradient,exact_shift_gradient,fd_errors,exact_rotation_error, &
        &exact_shift_error,exact_before,exact_after,initial_rotation_error,initial_shift_error, &
        &recovered_rotation_error,recovered_shift_error,recovery_before,recovery_after])), &
        &'matched-window diagnostic produced non-finite evidence')
    call assert_true(minval(fd_errors) < GRADIENT_FD_TOLERANCE, &
        &'matched-window Jacobian disagrees with centered finite differences')
    write(*,'(a)') 'CONTINUOUS_3D_MATCHED_WINDOW: EVIDENCE COMPLETE'

    do axis = 1, NCOORDS
        call native_jacobian(axis)%kill()
        call padded_jacobian(axis)%kill()
    enddo
    call native_model%kill()
    call padded_model%kill()
    call orientation%kill()
    call even_oris%kill()
    call odd_oris%kill()
    call workspace%kill()
    call model_operator%kill()
end subroutine run_matched_window_diagnostic

!> Measure exact stationarity and local recovery against one window-matched reference.
subroutine evaluate_matched_reference_batch(reference_volume,exact_rotmats,observed,shell_range,metrics)
    real, intent(in) :: reference_volume(BOX,BOX,BOX)
    real(dp), intent(in) :: exact_rotmats(:,:,:)
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:,:)
    integer, intent(in) :: shell_range(2)
    type(matched_reference_metrics), intent(out) :: metrics
    type(image) :: padded_model, native_model
    type(image) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    type(reconstructor_pcg) :: model_operator
    type(cartesian_pose_refiner) :: workspace
    complex, allocatable :: model_planes(:,:,:), jacobian_planes(:,:,:,:)
    real(dp), allocatable :: exact_estimates(:,:,:), recovery_estimates(:,:,:)
    real(dp), allocatable :: exact_shifts(:,:), recovery_shifts(:,:)
    real(dp) :: objective_before, objective_after
    integer :: accepted, axis, i, nparticles

    metrics = matched_reference_metrics()
    nparticles = size(exact_rotmats,3)
    if( size(observed,3) /= nparticles ) error stop 'matched reference batch size mismatch'
    call model_operator%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call model_operator%set_volume(reference_volume)
    call workspace%new_prepared_test(reference_volume)
    call workspace%set_shell_range(shell_range)
    call initialize_window_images(padded_model,native_model,padded_jacobian,native_jacobian)
    allocate(model_planes(-BOX/2:BOX/2,-BOX/2:BOX/2,nparticles),source=cmplx(0.,0.))
    allocate(jacobian_planes(-BOX/2:BOX/2,-BOX/2:BOX/2,NCOORDS,nparticles), &
        &source=cmplx(0.,0.))
    allocate(exact_estimates,source=exact_rotmats)
    allocate(recovery_estimates,source=exact_rotmats)
    allocate(exact_shifts(2,nparticles),source=0._dp)
    allocate(recovery_shifts(2,nparticles),source=0._dp)

    do i = 1, nparticles
        call evaluate_windowed_prediction(workspace,model_operator,exact_rotmats(:,:,i), &
            &[0._dp,0._dp],padded_model,native_model,padded_jacobian,native_jacobian, &
            &model_planes(:,:,i),jacobian_planes(:,:,:,i))
    enddo
    call fit_global_amplitude(model_planes,observed,metrics%amplitude_scale, &
        &metrics%fitted_residual,shell_range)
    call measure_exact_batch(model_planes,jacobian_planes,observed,metrics%amplitude_scale, &
        &metrics%exact_objective,metrics%rotation_gradient_rms,metrics%shift_gradient_rms, &
        &shell_range)

    ! Run identical exact and perturbed batches so stationarity and attraction are separate.
    do i = 1, nparticles
        call refine_windowed_pose(workspace,model_operator,observed(:,:,i),metrics%amplitude_scale, &
            &padded_model,native_model,padded_jacobian,native_jacobian,exact_estimates(:,:,i), &
            &exact_shifts(:,i),objective_before,objective_after,accepted,shell_range)
        metrics%exact_objective_before = metrics%exact_objective_before+objective_before
        metrics%exact_objective_after = metrics%exact_objective_after+objective_after
        metrics%exact_accepted = metrics%exact_accepted+accepted
        recovery_estimates(:,:,i) = right_increment_rotation( &
            &exact_rotmats(:,:,i),DIAGNOSTIC_ROTATION_OFFSET)
        recovery_shifts(:,i) = DIAGNOSTIC_SHIFT_OFFSET
        call refine_windowed_pose(workspace,model_operator,observed(:,:,i),metrics%amplitude_scale, &
            &padded_model,native_model,padded_jacobian,native_jacobian,recovery_estimates(:,:,i), &
            &recovery_shifts(:,i),objective_before,objective_after,accepted,shell_range)
        metrics%recovery_objective_before = metrics%recovery_objective_before+objective_before
        metrics%recovery_objective_after = metrics%recovery_objective_after+objective_after
        metrics%recovery_accepted = metrics%recovery_accepted+accepted
    enddo
    metrics%exact_rotation_rms = batch_rotation_rms(exact_rotmats,exact_estimates)
    metrics%exact_shift_rms = sqrt(sum(exact_shifts**2)/real(nparticles,dp))
    metrics%recovery_rotation_rms = batch_rotation_rms(exact_rotmats,recovery_estimates)
    metrics%recovery_shift_rms = sqrt(sum(recovery_shifts**2)/real(nparticles,dp))
    call gauge_corrected_errors(exact_rotmats,exact_estimates,exact_shifts, &
        &metrics%gauge_rotation_rms,metrics%gauge_shift_rms)

    do axis = 1, NCOORDS
        call native_jacobian(axis)%kill()
        call padded_jacobian(axis)%kill()
    enddo
    call native_model%kill()
    call padded_model%kill()
    call workspace%kill()
    call model_operator%kill()
end subroutine evaluate_matched_reference_batch

!> Allocate reusable padded and native images for one model and five derivatives.
subroutine initialize_window_images(padded_model,native_model,padded_jacobian,native_jacobian)
    type(image), intent(inout) :: padded_model, native_model
    type(image), intent(inout) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    integer :: axis

    call padded_model%new([PADDED_BOX,PADDED_BOX,1],HALFSET_SMPD,wthreads=.false.)
    call native_model%new([BOX,BOX,1],HALFSET_SMPD,wthreads=.false.)
    do axis = 1, NCOORDS
        call padded_jacobian(axis)%new([PADDED_BOX,PADDED_BOX,1], &
            &HALFSET_SMPD,wthreads=.false.)
        call native_jacobian(axis)%new([BOX,BOX,1],HALFSET_SMPD,wthreads=.false.)
    enddo
end subroutine initialize_window_images

!> Evaluate PCG values and derivatives on the padded grid, then clip them identically.
subroutine evaluate_windowed_prediction(workspace,sampler,rotmat,shift,padded_model,native_model, &
    &padded_jacobian,native_jacobian,model_plane,jacobian_plane)
    type(cartesian_pose_refiner), intent(in) :: workspace
    type(reconstructor_pcg), intent(in) :: sampler
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    type(image), intent(inout) :: padded_model, native_model
    type(image), intent(inout) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    complex, intent(out) :: model_plane(-BOX/2:,-BOX/2:)
    complex, intent(out) :: jacobian_plane(-BOX/2:,-BOX/2:,:)
    complex :: dvalue_dloc(3), phase, value
    complex(dp) :: jacobian(NCOORDS), model
    real(sp) :: loc(3), switch_margin(3)
    real(dp) :: arg, dloc(3,NCOORDS), frequency(2)
    integer :: axis, h, k, phys(3)

    call padded_model%zero_and_flag_ft()
    do axis = 1, NCOORDS
        call padded_jacobian(axis)%zero_and_flag_ft()
    enddo
    do k = -PADDED_RADIUS, PADDED_RADIUS-1
        do h = 0, PADDED_RADIUS
            if( h*h+k*k > PADDED_RADIUS**2 ) cycle
            loc = real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
            arg = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/ &
                &real(PADDED_BOX,dp)
            phase = cmplx(cos(arg),sin(arg),kind=sp)
            model = cmplx(phase*value,kind=dp)
            dloc(:,1) = [0._dp,real(loc(3),dp),-real(loc(2),dp)]
            dloc(:,2) = [-real(loc(3),dp),0._dp,real(loc(1),dp)]
            dloc(:,3) = [real(loc(2),dp),-real(loc(1),dp),0._dp]
            do axis = 1, 3
                jacobian(axis) = cmplx(phase,kind=dp)* &
                    &sum(cmplx(dvalue_dloc,kind=dp)*dloc(:,axis))
            enddo
            frequency = 2._dp*real(PI,dp)*real([h,k],dp)/real(PADDED_BOX,dp)
            jacobian(4:5) = cmplx(0._dp,frequency,kind=dp)*model
            phys = [h+1,k+1+merge(PADDED_BOX,0,k<0),1]
            call padded_model%set_cmat_at(phys(1),phys(2),phys(3),cmplx(model,kind=sp))
            do axis = 1, NCOORDS
                call padded_jacobian(axis)%set_cmat_at(phys(1),phys(2),phys(3), &
                    &cmplx(jacobian(axis),kind=sp))
            enddo
        enddo
    enddo

    ! The linear simulator window must act on the value and every Jacobian column.
    call padded_model%ifft()
    call padded_model%fft()
    call padded_model%ifft()
    call padded_model%clip(native_model)
    call native_model%fft()
    model_plane = sampler%extract_native_plane(native_model)
    do axis = 1, NCOORDS
        call padded_jacobian(axis)%ifft()
        call padded_jacobian(axis)%fft()
        call padded_jacobian(axis)%ifft()
        call padded_jacobian(axis)%clip(native_jacobian(axis))
        call native_jacobian(axis)%fft()
        jacobian_plane(:,:,axis) = sampler%extract_native_plane(native_jacobian(axis))
    enddo
end subroutine evaluate_windowed_prediction

!> Fit one fixed global real amplitude between the matched predictions and observations.
subroutine fit_global_amplitude(model,observed,scale,residual,shell_range)
    complex, intent(in) :: model(-BOX/2:,-BOX/2:,:), observed(-BOX/2:,-BOX/2:,:)
    real(dp), intent(out) :: scale, residual
    integer, optional, intent(in) :: shell_range(2)
    complex(dp) :: model_value, observed_value
    real(dp) :: cross_term, data_norm, fitted_squared, model_norm
    integer :: active_shells(2), h, i, k

    active_shells = [2,BOX/2]
    if( present(shell_range) ) active_shells = shell_range
    cross_term = 0._dp
    data_norm = 0._dp
    model_norm = 0._dp
    do i = 1, size(model,3)
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                if( h*h+k*k < active_shells(1)**2 .or. &
                    &h*h+k*k > active_shells(2)**2 ) cycle
                model_value = cmplx(model(h,k,i),kind=dp)
                observed_value = cmplx(observed(h,k,i),kind=dp)
                model_norm = model_norm+abs(model_value)**2
                data_norm = data_norm+abs(observed_value)**2
                cross_term = cross_term+real(conjg(model_value)*observed_value,dp)
            enddo
        enddo
    enddo
    scale = cross_term/max(model_norm,tiny(model_norm))
    fitted_squared = max(0._dp,data_norm+scale**2*model_norm-2._dp*scale*cross_term)
    residual = sqrt(fitted_squared/max(data_norm,tiny(data_norm)))
end subroutine fit_global_amplitude

!> Measure the exact-pose objective and particle-wise gradient RMS after window matching.
subroutine measure_exact_batch(model,jacobian,observed,scale,objective,rotation_gradient, &
    &shift_gradient,shell_range)
    complex, intent(in) :: model(-BOX/2:,-BOX/2:,:)
    complex, intent(in) :: jacobian(-BOX/2:,-BOX/2:,:,:)
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:,:)
    real(dp), intent(in) :: scale
    real(dp), intent(out) :: objective, rotation_gradient, shift_gradient
    integer, optional, intent(in) :: shell_range(2)
    real(dp) :: gradient(NCOORDS), hessian(NCOORDS,NCOORDS), particle_objective
    integer :: i

    objective = 0._dp
    rotation_gradient = 0._dp
    shift_gradient = 0._dp
    do i = 1, size(model,3)
        call plane_normal_terms(model(:,:,i),jacobian(:,:,:,i),observed(:,:,i),scale, &
            &particle_objective,gradient,hessian,shell_range)
        objective = objective+particle_objective
        rotation_gradient = rotation_gradient+sum(gradient(1:3)**2)
        shift_gradient = shift_gradient+sum(gradient(4:5)**2)
    enddo
    rotation_gradient = sqrt(rotation_gradient/real(size(model,3),dp))
    shift_gradient = sqrt(shift_gradient/real(size(model,3),dp))
end subroutine measure_exact_batch

!> Form one five-coordinate Gauss-Newton block from already windowed planes.
subroutine plane_normal_terms(model,jacobian,observed,scale,objective,gradient,hessian,shell_range)
    complex, intent(in) :: model(-BOX/2:,-BOX/2:)
    complex, intent(in) :: jacobian(-BOX/2:,-BOX/2:,:)
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    real(dp), intent(in) :: scale
    real(dp), intent(out) :: objective, gradient(NCOORDS), hessian(NCOORDS,NCOORDS)
    integer, optional, intent(in) :: shell_range(2)
    complex(dp) :: residual, scaled_jacobian(NCOORDS)
    integer :: active_shells(2), axis, h, jaxis, k

    active_shells = [2,BOX/2]
    if( present(shell_range) ) active_shells = shell_range
    objective = 0._dp
    gradient = 0._dp
    hessian = 0._dp
    do k = -BOX/2, BOX/2
        do h = -BOX/2, BOX/2
            if( h*h+k*k < active_shells(1)**2 .or. &
                &h*h+k*k > active_shells(2)**2 ) cycle
            residual = scale*cmplx(model(h,k),kind=dp)-cmplx(observed(h,k),kind=dp)
            scaled_jacobian = scale*cmplx(jacobian(h,k,:),kind=dp)
            objective = objective+0.5_dp*real(conjg(residual)*residual,dp)
            do axis = 1, NCOORDS
                gradient(axis) = gradient(axis)+real(conjg(scaled_jacobian(axis))*residual,dp)
                do jaxis = 1, NCOORDS
                    hessian(axis,jaxis) = hessian(axis,jaxis)+ &
                        &real(conjg(scaled_jacobian(axis))*scaled_jacobian(jaxis),dp)
                enddo
            enddo
        enddo
    enddo
end subroutine plane_normal_terms

!> Compare the transformed five-vector gradient with centered objective differences.
subroutine measure_windowed_gradient_fd(workspace,sampler,rotmat,shift,observed,scale, &
    &padded_model,native_model,padded_jacobian,native_jacobian,errors)
    type(cartesian_pose_refiner), intent(in) :: workspace
    type(reconstructor_pcg), intent(in) :: sampler
    real(dp), intent(in) :: rotmat(3,3), shift(2), scale
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    type(image), intent(inout) :: padded_model, native_model
    type(image), intent(inout) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    real(dp), intent(out) :: errors(NFD)
    real(dp) :: gradient(NCOORDS), hessian(NCOORDS,NCOORDS)
    real(dp) :: minus_objective, plus_objective, center_objective
    real(dp) :: minus_rotmat(3,3), plus_rotmat(3,3), finite_difference, directional_derivative
    integer :: istep

    call evaluate_pose_terms(workspace,sampler,rotmat,shift,observed,scale,padded_model, &
        &native_model,padded_jacobian,native_jacobian,center_objective,gradient,hessian)
    directional_derivative = dot_product(gradient,POSE_DIRECTION)
    if( .not. ieee_is_finite(center_objective) )then
        errors = huge(0._dp)
        return
    endif
    do istep = 1, NFD
        minus_rotmat = right_increment_rotation(rotmat,-FD_STEPS(istep)*POSE_DIRECTION(1:3))
        plus_rotmat = right_increment_rotation(rotmat,FD_STEPS(istep)*POSE_DIRECTION(1:3))
        call evaluate_pose_terms(workspace,sampler,minus_rotmat, &
            &shift-FD_STEPS(istep)*POSE_DIRECTION(4:5),observed,scale,padded_model, &
            &native_model,padded_jacobian,native_jacobian,minus_objective,gradient,hessian)
        call evaluate_pose_terms(workspace,sampler,plus_rotmat, &
            &shift+FD_STEPS(istep)*POSE_DIRECTION(4:5),observed,scale,padded_model, &
            &native_model,padded_jacobian,native_jacobian,plus_objective,gradient,hessian)
        finite_difference = (plus_objective-minus_objective)/(2._dp*FD_STEPS(istep))
        errors(istep) = abs(finite_difference-directional_derivative)/ &
            &max(1.e-8_dp,abs(finite_difference),abs(directional_derivative))
    enddo
end subroutine measure_windowed_gradient_fd

!> Evaluate one window-matched objective, gradient, and Gauss-Newton block.
subroutine evaluate_pose_terms(workspace,sampler,rotmat,shift,observed,scale,padded_model, &
    &native_model,padded_jacobian,native_jacobian,objective,gradient,hessian,shell_range)
    type(cartesian_pose_refiner), intent(in) :: workspace
    type(reconstructor_pcg), intent(in) :: sampler
    real(dp), intent(in) :: rotmat(3,3), shift(2), scale
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    type(image), intent(inout) :: padded_model, native_model
    type(image), intent(inout) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    real(dp), intent(out) :: objective, gradient(NCOORDS), hessian(NCOORDS,NCOORDS)
    integer, optional, intent(in) :: shell_range(2)
    complex :: model(-BOX/2:BOX/2,-BOX/2:BOX/2)
    complex :: jacobian(-BOX/2:BOX/2,-BOX/2:BOX/2,NCOORDS)

    call evaluate_windowed_prediction(workspace,sampler,rotmat,shift,padded_model,native_model, &
        &padded_jacobian,native_jacobian,model,jacobian)
    call plane_normal_terms(model,jacobian,observed,scale,objective,gradient,hessian,shell_range)
end subroutine evaluate_pose_terms

!> Run a small diagnostic LM solve with fixed amplitude and bounded pose steps.
subroutine refine_windowed_pose(workspace,sampler,observed,scale,padded_model,native_model, &
    &padded_jacobian,native_jacobian,rotmat,shift,objective_before,objective_after,naccepted, &
    &shell_range)
    type(cartesian_pose_refiner), intent(in) :: workspace
    type(reconstructor_pcg), intent(in) :: sampler
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:)
    real(dp), intent(in) :: scale
    type(image), intent(inout) :: padded_model, native_model
    type(image), intent(inout) :: padded_jacobian(NCOORDS), native_jacobian(NCOORDS)
    real(dp), intent(inout) :: rotmat(3,3), shift(2)
    real(dp), intent(out) :: objective_before, objective_after
    integer, intent(out) :: naccepted
    integer, optional, intent(in) :: shell_range(2)
    real(dp) :: gradient(NCOORDS), hessian(NCOORDS,NCOORDS), solve_matrix(NCOORDS,NCOORDS)
    real(dp) :: direction(NCOORDS), diagonal(NCOORDS), trial_gradient(NCOORDS)
    real(dp) :: trial_hessian(NCOORDS,NCOORDS), trial_rotmat(3,3), trial_shift(2)
    real(dp) :: objective, trial_objective, mu, maxdiag, rotation_norm, shift_norm, step_scale
    integer :: axis, iteration
    logical :: reliable

    call evaluate_pose_terms(workspace,sampler,rotmat,shift,observed,scale,padded_model, &
        &native_model,padded_jacobian,native_jacobian,objective,gradient,hessian,shell_range)
    objective_before = objective
    mu = 1.e-3_dp
    naccepted = 0
    do iteration = 1, MAX_DIAGNOSTIC_ITERATIONS
        maxdiag = max(maxval([(hessian(axis,axis),axis=1,NCOORDS)]),NUMERIC_FLOOR)
        diagonal = [(max(hessian(axis,axis),1.e-6_dp*maxdiag),axis=1,NCOORDS)]
        solve_matrix = hessian
        do axis = 1, NCOORDS
            solve_matrix(axis,axis) = solve_matrix(axis,axis)+mu*diagonal(axis)
        enddo
        call solve_cholesky(solve_matrix,-gradient,direction,reliable)
        if( .not. reliable )then
            mu = 10._dp*mu
            cycle
        endif
        rotation_norm = sqrt(sum(direction(1:3)**2))
        shift_norm = sqrt(sum(direction(4:5)**2))
        step_scale = 1._dp
        if( rotation_norm > MAX_ROTATION_STEP ) step_scale = min(step_scale,MAX_ROTATION_STEP/rotation_norm)
        if( shift_norm > MAX_SHIFT_STEP ) step_scale = min(step_scale,MAX_SHIFT_STEP/shift_norm)
        direction = step_scale*direction
        trial_rotmat = right_increment_rotation(rotmat,direction(1:3))
        trial_shift = shift+direction(4:5)
        call evaluate_pose_terms(workspace,sampler,trial_rotmat,trial_shift,observed,scale, &
            &padded_model,native_model,padded_jacobian,native_jacobian,trial_objective, &
            &trial_gradient,trial_hessian,shell_range)
        if( trial_objective < objective )then
            rotmat = trial_rotmat
            shift = trial_shift
            objective = trial_objective
            gradient = trial_gradient
            hessian = trial_hessian
            naccepted = naccepted+1
            mu = max(1.e-8_dp,0.25_dp*mu)
            if( max(rotation_norm,shift_norm)*step_scale < 1.e-8_dp ) exit
        else
            mu = 4._dp*mu
        endif
    enddo
    objective_after = objective
end subroutine refine_windowed_pose

!> Solve one symmetric positive-definite five-coordinate diagnostic system.
pure subroutine solve_cholesky(matrix,rhs,solution,reliable)
    real(dp), intent(in) :: matrix(NCOORDS,NCOORDS), rhs(NCOORDS)
    real(dp), intent(out) :: solution(NCOORDS)
    logical, intent(out) :: reliable
    real(dp) :: lower(NCOORDS,NCOORDS), intermediate(NCOORDS)
    real(dp) :: matrix_scale, pivot, pivot_floor
    integer :: i, j

    solution = 0._dp
    intermediate = 0._dp
    lower = 0._dp
    reliable = .false.
    if( any(.not. ieee_is_finite(matrix)) .or. any(.not. ieee_is_finite(rhs)) ) return
    matrix_scale = maxval(abs(matrix))
    if( matrix_scale <= NUMERIC_FLOOR ) return
    pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
    do i = 1, NCOORDS
        do j = 1, i-1
            lower(i,j) = (matrix(i,j)-dot_product(lower(i,1:j-1),lower(j,1:j-1)))/lower(j,j)
        enddo
        pivot = matrix(i,i)-dot_product(lower(i,1:i-1),lower(i,1:i-1))
        if( .not. ieee_is_finite(pivot) .or. pivot <= pivot_floor ) return
        lower(i,i) = sqrt(pivot)
    enddo
    do i = 1, NCOORDS
        intermediate(i) = (rhs(i)-dot_product(lower(i,1:i-1),intermediate(1:i-1)))/lower(i,i)
    enddo
    do i = NCOORDS, 1, -1
        solution(i) = (intermediate(i)-dot_product(lower(i+1:NCOORDS,i), &
            &solution(i+1:NCOORDS)))/lower(i,i)
    enddo
    reliable = all(ieee_is_finite(solution))
end subroutine solve_cholesky

!> Remove the best common left rotation and projected three-dimensional translation.
subroutine gauge_corrected_errors(reference,estimate,shifts,rotation_rms_value,shift_rms_value)
    real(dp), intent(in) :: reference(:,:,:), estimate(:,:,:), shifts(:,:)
    real(dp), intent(out) :: rotation_rms_value, shift_rms_value
    real(dp) :: gauge(3,3), identity(3,3), mean_rotation(3), relative(3,3)
    real(dp) :: normal(3,3), rhs(3), translation(3), row(3), predicted_shift(2)
    real(dp) :: aligned(3,3), squared_rotation, squared_shift
    integer :: axis, i
    logical :: reliable

    mean_rotation = 0._dp
    do i = 1, size(reference,3)
        relative = matmul(estimate(:,:,i),transpose(reference(:,:,i)))
        mean_rotation = mean_rotation+rotation_vector(relative)
    enddo
    mean_rotation = mean_rotation/real(size(reference,3),dp)
    identity = 0._dp
    identity(1,1) = 1._dp
    identity(2,2) = 1._dp
    identity(3,3) = 1._dp
    gauge = right_increment_rotation(identity,mean_rotation)
    squared_rotation = 0._dp
    do i = 1, size(reference,3)
        aligned = matmul(transpose(gauge),estimate(:,:,i))
        squared_rotation = squared_rotation+rotation_error(reference(:,:,i),aligned)**2
    enddo
    rotation_rms_value = sqrt(squared_rotation/real(size(reference,3),dp))

    ! A global volume translation t induces image shifts -[R t]_1:2.
    normal = 0._dp
    rhs = 0._dp
    do i = 1, size(reference,3)
        do axis = 1, 2
            row = reference(axis,:,i)
            normal = normal+spread(row,2,3)*spread(row,1,3)
            rhs = rhs-row*shifts(axis,i)
        enddo
    enddo
    call solve_cholesky3(normal,rhs,translation,reliable)
    if( .not. reliable )then
        shift_rms_value = huge(0._dp)
        return
    endif
    squared_shift = 0._dp
    do i = 1, size(reference,3)
        predicted_shift = -matmul(reference(1:2,:,i),translation)
        squared_shift = squared_shift+sum((shifts(:,i)-predicted_shift)**2)
    enddo
    shift_rms_value = sqrt(squared_shift/real(size(reference,3),dp))
end subroutine gauge_corrected_errors

!> Return the axis-angle logarithm of one proper rotation matrix.
pure function rotation_vector(matrix) result(vector)
    real(dp), intent(in) :: matrix(3,3)
    real(dp) :: vector(3), cosine, scale, sine, theta

    cosine = 0.5_dp*(matrix(1,1)+matrix(2,2)+matrix(3,3)-1._dp)
    theta = acos(max(-1._dp,min(1._dp,cosine)))
    vector = [matrix(3,2)-matrix(2,3),matrix(1,3)-matrix(3,1), &
        &matrix(2,1)-matrix(1,2)]
    if( theta < 1.e-8_dp )then
        vector = 0.5_dp*vector
    else
        sine = sin(theta)
        scale = theta/(2._dp*sign(max(abs(sine),NUMERIC_FLOOR),sine))
        vector = scale*vector
    endif
end function rotation_vector

!> Solve the three-coordinate translation-gauge normal system.
pure subroutine solve_cholesky3(matrix,rhs,solution,reliable)
    real(dp), intent(in) :: matrix(3,3), rhs(3)
    real(dp), intent(out) :: solution(3)
    logical, intent(out) :: reliable
    real(dp) :: lower(3,3), intermediate(3), matrix_scale, pivot, pivot_floor
    integer :: i, j

    lower = 0._dp
    intermediate = 0._dp
    solution = 0._dp
    reliable = .false.
    matrix_scale = maxval(abs(matrix))
    if( matrix_scale <= NUMERIC_FLOOR ) return
    pivot_floor = sqrt(epsilon(1._dp))*matrix_scale
    do i = 1, 3
        do j = 1, i-1
            lower(i,j) = (matrix(i,j)-dot_product(lower(i,1:j-1),lower(j,1:j-1)))/lower(j,j)
        enddo
        pivot = matrix(i,i)-dot_product(lower(i,1:i-1),lower(i,1:i-1))
        if( pivot <= pivot_floor ) return
        lower(i,i) = sqrt(pivot)
    enddo
    do i = 1, 3
        intermediate(i) = (rhs(i)-dot_product(lower(i,1:i-1),intermediate(1:i-1)))/lower(i,i)
    enddo
    do i = 3, 1, -1
        solution(i) = (intermediate(i)-dot_product(lower(i+1:3,i),solution(i+1:3)))/lower(i,i)
    enddo
    reliable = all(ieee_is_finite(solution))
end subroutine solve_cholesky3

!> Return rotation RMS for two pose batches in one fixed gauge.
pure real(dp) function batch_rotation_rms(reference,estimate) result(error)
    real(dp), intent(in) :: reference(:,:,:), estimate(:,:,:)
    real(dp) :: squared_error
    integer :: i

    squared_error = 0._dp
    do i = 1, size(reference,3)
        squared_error = squared_error+rotation_error(reference(:,:,i),estimate(:,:,i))**2
    enddo
    error = sqrt(squared_error/real(size(reference,3),dp))
end function batch_rotation_rms

!> Report whether matched clipping restores exact stationarity and local recovery.
subroutine print_matched_window_decision(residual,exact_rotation,exact_shift,initial_rotation, &
    &initial_shift,recovered_rotation,recovered_shift)
    real(dp), intent(in) :: residual, exact_rotation, exact_shift
    real(dp), intent(in) :: initial_rotation, initial_shift, recovered_rotation, recovered_shift

    if( residual > 1.e-3_dp )then
        write(*,'(a)') 'CONTINUOUS_3D_MATCHED_WINDOW DIAGNOSIS: WINDOWED_FORWARD_MISMATCH'
    else if( exact_rotation > EXACT_ROTATION_TOLERANCE .or. exact_shift > EXACT_SHIFT_TOLERANCE )then
        write(*,'(a)') 'CONTINUOUS_3D_MATCHED_WINDOW DIAGNOSIS: EXACT_POSE_NOT_STATIONARY'
    else if( recovered_rotation >= initial_rotation .or. recovered_shift >= initial_shift )then
        write(*,'(a)') 'CONTINUOUS_3D_MATCHED_WINDOW DIAGNOSIS: LOCAL_RECOVERY_FAILURE'
    else
        write(*,'(a)') 'CONTINUOUS_3D_MATCHED_WINDOW DIAGNOSIS: MATCHED_WINDOW_SUPPORTED'
    endif
end subroutine print_matched_window_decision

!> Return the geodesic rotation error for one matrix pair.
pure real(dp) function rotation_error(reference,estimate) result(error)
    real(dp), intent(in) :: reference(3,3), estimate(3,3)
    real(dp) :: cosine

    cosine = 0.5_dp*(sum(reference*estimate)-1._dp)
    error = acos(max(-1._dp,min(1._dp,cosine)))
end function rotation_error

end module pose_cont_refinement_matched_window_test
