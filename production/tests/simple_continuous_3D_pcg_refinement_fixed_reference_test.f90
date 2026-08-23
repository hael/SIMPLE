module continuous_3D_pcg_refinement_fixed_reference_test
use continuous_3D_pcg_refinement_halfset_support, only: build_disjoint_half_orientations, &
    &build_independent_observations, HALFSET_LAMBDA, HALFSET_MASK_RADIUS, HALFSET_SMPD, &
    &reconstruct_half_pair_fixed
use continuous_3D_pcg_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, OSMPL_PAD_FAC
use simple_image, only: image
use simple_ori, only: ori
use simple_oris, only: oris
use simple_pcg_pose_polisher, only: pcg_pose_polish_summary, polish_fixed_volume_poses
use simple_projector, only: projector
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg
implicit none
private
public :: run_fixed_reference_diagnostic, measure_exact_pose_fit

integer, parameter :: NPARTICLES = 48
integer, parameter :: MAX_POLISH_ITERATIONS = 8
integer, parameter :: NOISE_SEED = 20260818
integer, parameter :: NBANDS = 3
integer, parameter :: NARMS = 11
integer, parameter :: ARM_NATIVE_RAW = 1
integer, parameter :: ARM_NATIVE_CORRECTED = 2
integer, parameter :: ARM_SIMULATOR_RAW = 3
integer, parameter :: ARM_SIMULATOR_SCALED_RAW = 4
integer, parameter :: ARM_SIMULATOR_CORRECTED = 5
integer, parameter :: ARM_SIMULATOR_SCALED_CORRECTED = 6
integer, parameter :: ARM_EXACT_KB_CORRECTED = 7
integer, parameter :: ARM_EXACT_KB_SCALED = 8
integer, parameter :: ARM_MATCHED_TRUTH = 9
integer, parameter :: ARM_SIMULATOR_RAW_RECONSTRUCTION = 10
integer, parameter :: ARM_SIMULATOR_CORRECTED_RECONSTRUCTION = 11
real(dp), parameter :: EXACT_ROTATION_TOLERANCE = 0.001745_dp
real(dp), parameter :: EXACT_SHIFT_TOLERANCE = 0.05_dp
character(len=36), parameter :: ARM_NAMES(NARMS) = [character(len=36) :: &
    &'native_raw_truth', 'native_corrected_truth', 'simulator_raw_truth', &
    &'simulator_scaled_raw_truth', 'simulator_corrected_truth', &
    &'simulator_scaled_corrected_truth', &
    &'exact_kb_corrected_truth', 'exact_kb_scaled_corrected_truth', &
    &'pcg_matched_corrected_truth', 'simulator_raw_recon', 'simulator_corrected_recon']

contains

!> Collect one matrix that isolates padding, scale, KB, LM, and map-reference drift.
subroutine run_fixed_reference_diagnostic()
    type(ori) :: orientation
    type(oris) :: even_oris, odd_oris
    type(reconstructor_pcg) :: sampler, model_operator
    type(pcg_fourier_workspace) :: raw_truth_workspace, corrected_truth_workspace
    type(pcg_fourier_workspace) :: simulator_scaled_raw_workspace
    type(pcg_fourier_workspace) :: simulator_scaled_corrected_workspace, exact_kb_scaled_workspace
    type(pcg_fourier_workspace) :: reconstructed_workspace, corrected_reconstructed_workspace
    type(pcg_pose_polish_summary) :: summaries(NARMS)
    complex, allocatable :: noisy_planes(:,:,:), native_planes(:,:,:), matched_planes(:,:,:)
    complex, allocatable :: simulator_planes(:,:,:), exact_kb_planes(:,:,:)
    real, allocatable :: clean_images(:,:,:), noisy_images(:,:,:), noise(:,:,:)
    real, allocatable :: truth_volume(:,:,:), reconstructed_volume(:,:,:), duplicate_volume(:,:,:)
    real, allocatable :: inverse_envelope(:,:,:), corrected_truth_volume(:,:,:)
    real, allocatable :: simulator_scaled_raw_volume(:,:,:)
    real, allocatable :: simulator_scaled_corrected_volume(:,:,:), exact_kb_scaled_volume(:,:,:)
    real, allocatable :: corrected_reconstructed_volume(:,:,:)
    real(dp), allocatable :: exact_rotmats(:,:,:), arm_rotmats(:,:,:,:), arm_shifts(:,:,:)
    integer, allocatable :: even_ids(:), odd_ids(:), statuses(:,:)
    integer :: arm, i, niters(2)
    real(dp) :: ignored_snr, residuals(2), rotation_scale, ignored_objective
    real(dp) :: ignored_residual, ignored_fitted_residual, ignored_rotation_gradient
    real(dp) :: ignored_shift_gradient, ignored_band_residuals(NBANDS), ignored_band_correlations(NBANDS)
    real(dp) :: simulator_raw_scale, simulator_corrected_scale, exact_kb_scale
    real(dp) :: rotation_rms_values(NARMS), shift_rms_values(NARMS)
    real(dp) :: exact_objectives(NARMS), forward_residuals(NARMS), fitted_residuals(NARMS)
    real(dp) :: amplitude_scales(NARMS), band_residuals(NBANDS,NARMS)
    real(dp) :: band_correlations(NBANDS,NARMS)
    real(dp) :: rotation_gradient_rms(NARMS), shift_gradient_rms(NARMS)

    ! Keep the legacy native projector as one explicit numerical convention.
    call build_disjoint_half_orientations(NPARTICLES,even_oris,odd_oris,even_ids,odd_ids)
    call sampler%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call build_independent_observations(sampler,even_oris,NOISE_SEED,1.0, &
        &noisy_planes,simulator_planes,clean_images,noisy_images,noise,ignored_snr)
    call build_truth_volume(truth_volume)
    ! Retain the earlier native-box generator as a labeled comparison arm.
    call build_native_observations(sampler,even_oris,truth_volume,native_planes)

    ! Reconstruct the same simulator-path observations used by the production-like arms.
    call reconstruct_half_pair_fixed(even_oris,simulator_planes,simulator_planes,HALFSET_LAMBDA,2, &
        &reconstructed_volume,duplicate_volume,niters,residuals)
    call model_operator%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call model_operator%set_volume(truth_volume)
    call model_operator%begin_fourier_workspace(raw_truth_workspace)
    call raw_truth_workspace%set_shell_range([2,BOX/2])

    ! Real-data modeling requires A x = A_tilde E^-1 x.
    inverse_envelope = model_operator%get_invenv()
    corrected_truth_volume = truth_volume*inverse_envelope
    call model_operator%set_volume(corrected_truth_volume)
    call model_operator%begin_fourier_workspace(corrected_truth_workspace)
    call corrected_truth_workspace%set_shell_range([2,BOX/2])
    call model_operator%set_volume(reconstructed_volume)
    call model_operator%begin_fourier_workspace(reconstructed_workspace)
    call reconstructed_workspace%set_shell_range([2,BOX/2])
    corrected_reconstructed_volume = reconstructed_volume*inverse_envelope
    call model_operator%set_volume(corrected_reconstructed_volume)
    call model_operator%begin_fourier_workspace(corrected_reconstructed_workspace)
    call corrected_reconstructed_workspace%set_shell_range([2,BOX/2])

    allocate(exact_rotmats(3,3,NPARTICLES))
    allocate(arm_rotmats(3,3,NPARTICLES,NARMS),arm_shifts(2,NPARTICLES,NARMS),source=0._dp)
    allocate(statuses(NPARTICLES,NARMS))
    call orientation%new(.false.)
    do i = 1, NPARTICLES
        call even_oris%get_ori(i,orientation)
        exact_rotmats(:,:,i) = real(orientation%get_mat(),dp)
    enddo
    do arm = 1, NARMS
        arm_rotmats(:,:,:,arm) = exact_rotmats
    enddo
    rotation_scale = 1._dp/real(HALFSET_MASK_RADIUS,dp)

    allocate(matched_planes,mold=native_planes)
    call build_matched_observations(corrected_truth_workspace,exact_rotmats,matched_planes)

    ! Compare the exact KB padded slice with the fast KB PCG gather at one volume.
    call build_exact_kb_observations(corrected_truth_volume,even_oris,exact_kb_planes)

    ! Fit only a global volume amplitude; pose, shells, and interpolation stay fixed.
    call measure_exact_pose_fit(raw_truth_workspace,exact_rotmats,simulator_planes, &
        &ignored_objective,ignored_residual,ignored_fitted_residual,simulator_raw_scale, &
        &ignored_rotation_gradient,ignored_shift_gradient,ignored_band_residuals, &
        &ignored_band_correlations)
    simulator_scaled_raw_volume = simulator_raw_scale*truth_volume
    call model_operator%set_volume(simulator_scaled_raw_volume)
    call model_operator%begin_fourier_workspace(simulator_scaled_raw_workspace)
    call simulator_scaled_raw_workspace%set_shell_range([2,BOX/2])
    call measure_exact_pose_fit(corrected_truth_workspace,exact_rotmats,simulator_planes, &
        &ignored_objective,ignored_residual,ignored_fitted_residual,simulator_corrected_scale, &
        &ignored_rotation_gradient,ignored_shift_gradient,ignored_band_residuals, &
        &ignored_band_correlations)
    simulator_scaled_corrected_volume = simulator_corrected_scale*corrected_truth_volume
    call model_operator%set_volume(simulator_scaled_corrected_volume)
    call model_operator%begin_fourier_workspace(simulator_scaled_corrected_workspace)
    call simulator_scaled_corrected_workspace%set_shell_range([2,BOX/2])
    call measure_exact_pose_fit(corrected_truth_workspace,exact_rotmats,exact_kb_planes, &
        &ignored_objective,ignored_residual,ignored_fitted_residual,exact_kb_scale, &
        &ignored_rotation_gradient,ignored_shift_gradient,ignored_band_residuals, &
        &ignored_band_correlations)
    exact_kb_scaled_volume = exact_kb_scale*corrected_truth_volume
    call model_operator%set_volume(exact_kb_scaled_volume)
    call model_operator%begin_fourier_workspace(exact_kb_scaled_workspace)
    call exact_kb_scaled_workspace%set_shell_range([2,BOX/2])

    call run_diagnostic_arm(raw_truth_workspace,native_planes,ARM_NATIVE_RAW)
    call run_diagnostic_arm(corrected_truth_workspace,native_planes,ARM_NATIVE_CORRECTED)
    call run_diagnostic_arm(raw_truth_workspace,simulator_planes,ARM_SIMULATOR_RAW)
    call run_diagnostic_arm(simulator_scaled_raw_workspace,simulator_planes,ARM_SIMULATOR_SCALED_RAW)
    call run_diagnostic_arm(corrected_truth_workspace,simulator_planes,ARM_SIMULATOR_CORRECTED)
    call run_diagnostic_arm(simulator_scaled_corrected_workspace,simulator_planes, &
        &ARM_SIMULATOR_SCALED_CORRECTED)
    call run_diagnostic_arm(corrected_truth_workspace,exact_kb_planes,ARM_EXACT_KB_CORRECTED)
    call run_diagnostic_arm(exact_kb_scaled_workspace,exact_kb_planes,ARM_EXACT_KB_SCALED)
    call run_diagnostic_arm(corrected_truth_workspace,matched_planes,ARM_MATCHED_TRUTH)
    call run_diagnostic_arm(reconstructed_workspace,simulator_planes, &
        &ARM_SIMULATOR_RAW_RECONSTRUCTION)
    call run_diagnostic_arm(corrected_reconstructed_workspace,simulator_planes, &
        &ARM_SIMULATOR_CORRECTED_RECONSTRUCTION)

    do arm = 1, NARMS
        call print_diagnostic_arm(ARM_NAMES(arm),exact_objectives(arm),forward_residuals(arm), &
            &fitted_residuals(arm),amplitude_scales(arm),rotation_gradient_rms(arm), &
            &shift_gradient_rms(arm),band_residuals(:,arm),band_correlations(:,arm), &
            &rotation_rms_values(arm),shift_rms_values(arm),summaries(arm))
    enddo
    write(*,'(a,2(es14.6,1x))') 'CONTINUOUS_3D_FIXED_REFERENCE PCG residuals: ', residuals
    call print_diagnostic_decision(rotation_rms_values,shift_rms_values,fitted_residuals)

    call assert_true(niters(1) == 2 .and. niters(2) == 2, &
        &'fixed-reference diagnostic did not run two PCG iterations')
    call assert_true(all(summaries%ninvalid == 0), &
        &'fixed-reference diagnostic produced invalid numerics')
    do arm = 1, NARMS
        call assert_true(terminal_count(summaries(arm)) == NPARTICLES, &
            &'fixed-reference diagnostic did not account for every particle')
    enddo
    call assert_true(all(ieee_is_finite([rotation_rms_values,shift_rms_values,exact_objectives, &
        &forward_residuals,fitted_residuals,amplitude_scales,rotation_gradient_rms, &
        &shift_gradient_rms,reshape(band_residuals,[NBANDS*NARMS]), &
        &reshape(band_correlations,[NBANDS*NARMS])])), &
        &'fixed-reference diagnostic produced a non-finite pose error')
    call assert_true(rotation_rms_values(ARM_MATCHED_TRUTH) <= EXACT_ROTATION_TOLERANCE .and. &
        &shift_rms_values(ARM_MATCHED_TRUTH) <= EXACT_SHIFT_TOLERANCE, &
        &'operator-matched exact poses drift in the LM or Jacobian path')
    write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE: EVIDENCE COMPLETE'

    call orientation%kill()
    call even_oris%kill()
    call odd_oris%kill()
    call raw_truth_workspace%kill()
    call corrected_truth_workspace%kill()
    call simulator_scaled_raw_workspace%kill()
    call simulator_scaled_corrected_workspace%kill()
    call exact_kb_scaled_workspace%kill()
    call reconstructed_workspace%kill()
    call corrected_reconstructed_workspace%kill()
    call model_operator%kill()
    call sampler%kill()

contains

    !> Run one collect-all arm without applying unresolved scientific gates.
    subroutine run_diagnostic_arm(workspace,observed,arm_index)
        type(pcg_fourier_workspace), intent(in) :: workspace
        complex, intent(in) :: observed(-BOX/2:,-BOX/2:,:)
        integer, intent(in) :: arm_index

        call measure_and_polish(workspace,observed,exact_rotmats,rotation_scale, &
            &arm_rotmats(:,:,:,arm_index),arm_shifts(:,:,arm_index),statuses(:,arm_index), &
            &summaries(arm_index),exact_objectives(arm_index),forward_residuals(arm_index), &
            &fitted_residuals(arm_index),amplitude_scales(arm_index), &
            &rotation_gradient_rms(arm_index),shift_gradient_rms(arm_index), &
            &band_residuals(:,arm_index),band_correlations(:,arm_index), &
            &rotation_rms_values(arm_index),shift_rms_values(arm_index))
    end subroutine run_diagnostic_arm
end subroutine run_fixed_reference_diagnostic

!> Retain the earlier native-box projector convention as a diagnostic control.
subroutine build_native_observations(sampler,orientations,volume,observed)
    type(reconstructor_pcg), intent(in) :: sampler
    type(oris), intent(inout) :: orientations
    real, intent(in) :: volume(BOX,BOX,BOX)
    complex, allocatable, intent(out) :: observed(:,:,:)
    type(image) :: projection
    type(projector) :: native_projector
    type(ori) :: orientation
    integer :: i, lims2(2,2), nprojs

    nprojs = orientations%get_noris()
    lims2 = sampler%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojs))
    call native_projector%new([BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call native_projector%set_rmat(volume,.false.)
    call native_projector%fft()
    call native_projector%expand_cmat()
    call projection%new([BOX,BOX,1],HALFSET_SMPD,wthreads=.false.)
    call orientation%new(.false.)
    do i = 1, nprojs
        call orientations%get_ori(i,orientation)
        call native_projector%fproject_serial(orientation,projection)
        observed(:,:,i) = sampler%extract_native_plane(projection)
    enddo
    call orientation%kill()
    call projection%kill()
    call native_projector%kill_expanded()
    call native_projector%kill()
end subroutine build_native_observations

!> Generate padded exact-KB Fourier slices without the simulator's real-space clip.
subroutine build_exact_kb_observations(volume,orientations,observed)
    real, intent(in) :: volume(BOX,BOX,BOX)
    type(oris), intent(inout) :: orientations
    complex, allocatable, intent(out) :: observed(:,:,:)
    type(image) :: source_volume
    type(projector) :: padded_projector
    type(ori) :: orientation
    real :: loc(3), rotmat(3,3)
    integer :: h, i, k, nprojs

    nprojs = orientations%get_noris()
    allocate(observed(-BOX/2:BOX/2,-BOX/2:BOX/2,nprojs),source=cmplx(0.,0.))
    call source_volume%new([BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%set_rmat(volume,.false.)
    call padded_projector%new(OSMPL_PAD_FAC*[BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%pad(padded_projector)
    call padded_projector%fft()
    call padded_projector%expand_cmat()
    call orientation%new(.false.)
    do i = 1, nprojs
        call orientations%get_ori(i,orientation)
        rotmat = orientation%get_mat()
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                if( h*h+k*k > (BOX/2)**2 ) cycle
                loc = matmul(real([h,k,0]),rotmat)
                observed(h,k,i) = padded_projector%interp_fcomp_oversamp(loc)
            enddo
        enddo
    enddo
    call orientation%kill()
    call padded_projector%kill_expanded()
    call padded_projector%kill()
    call source_volume%kill()
end subroutine build_exact_kb_observations

!> Generate an inverse-crime control with the corrected PCG model itself.
subroutine build_matched_observations(workspace,rotmats,observed)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:)
    complex, intent(out) :: observed(-BOX/2:,-BOX/2:,:)
    complex, allocatable :: zero_plane(:,:)
    real(dp) :: ignored_objective
    integer :: i

    allocate(zero_plane(lbound(observed,1):ubound(observed,1), &
        &lbound(observed,2):ubound(observed,2)),source=cmplx(0.,0.))
    do i = 1, size(rotmats,3)
        call workspace%shift_residual(rotmats(:,:,i),[0._dp,0._dp],zero_plane, &
            &observed(:,:,i),ignored_objective)
    enddo
end subroutine build_matched_observations

!> Measure exact-pose fit, then run the unchanged five-coordinate LM policy.
subroutine measure_and_polish(workspace,observed,exact_rotmats,rotation_scale,rotmats,shifts, &
    &statuses,summary,exact_objective,forward_residual,fitted_residual,amplitude_scale, &
    &rotation_gradient,shift_gradient,band_residuals,band_correlations,final_rotation_rms, &
    &final_shift_rms)
    type(pcg_fourier_workspace), intent(in) :: workspace
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:,:)
    real(dp), intent(in) :: exact_rotmats(:,:,:), rotation_scale
    real(dp), intent(inout) :: rotmats(:,:,:), shifts(:,:)
    integer, intent(out) :: statuses(:)
    type(pcg_pose_polish_summary), intent(out) :: summary
    real(dp), intent(out) :: exact_objective, forward_residual, fitted_residual, amplitude_scale
    real(dp), intent(out) :: rotation_gradient, shift_gradient
    real(dp), intent(out) :: band_residuals(NBANDS), band_correlations(NBANDS)
    real(dp), intent(out) :: final_rotation_rms, final_shift_rms

    call measure_exact_pose_fit(workspace,exact_rotmats,observed,exact_objective, &
        &forward_residual,fitted_residual,amplitude_scale,rotation_gradient,shift_gradient, &
        &band_residuals,band_correlations)
    call polish_fixed_volume_poses(workspace,rotmats,observed,shifts,rotation_scale, &
        &MAX_POLISH_ITERATIONS,statuses,summary)
    final_rotation_rms = rotation_rms(exact_rotmats,rotmats)
    final_shift_rms = sqrt(sum(shifts**2)/real(size(shifts,2),dp))
end subroutine measure_and_polish

!> Report raw and amplitude-fitted residuals, shell fits, and pose gradients at truth.
subroutine measure_exact_pose_fit(workspace,rotmats,observed,objective,relative_residual, &
    &fitted_residual,amplitude_scale,rotation_gradient_rms,shift_gradient_rms, &
    &band_residuals,band_correlations)
    type(pcg_fourier_workspace), intent(in) :: workspace
    real(dp), intent(in) :: rotmats(:,:,:)
    complex, intent(in) :: observed(-BOX/2:,-BOX/2:,:)
    real(dp), intent(out) :: objective, relative_residual, fitted_residual, amplitude_scale
    real(dp), intent(out) :: rotation_gradient_rms, shift_gradient_rms
    real(dp), intent(out) :: band_residuals(NBANDS), band_correlations(NBANDS)
    complex, allocatable :: model(:,:), residual(:,:), zero_plane(:,:)
    complex(dp) :: model_value, observed_value
    real(dp) :: data_norm, model_norm, cross_term, residual_norm
    real(dp) :: band_data(NBANDS), band_model(NBANDS), band_cross(NBANDS)
    real(dp) :: fitted_squared
    real(dp) :: particle_objective, ignored_objective, gradient(5)
    integer :: band, h, i, k, radius

    allocate(model(lbound(observed,1):ubound(observed,1),lbound(observed,2):ubound(observed,2)))
    allocate(residual,mold=model)
    allocate(zero_plane,mold=model)
    zero_plane = cmplx(0.,0.)
    objective = 0._dp
    data_norm = 0._dp
    model_norm = 0._dp
    cross_term = 0._dp
    residual_norm = 0._dp
    band_data = 0._dp
    band_model = 0._dp
    band_cross = 0._dp
    rotation_gradient_rms = 0._dp
    shift_gradient_rms = 0._dp
    do i = 1, size(rotmats,3)
        call workspace%shift_residual(rotmats(:,:,i),[0._dp,0._dp],zero_plane,model,ignored_objective)
        call workspace%shift_residual(rotmats(:,:,i),[0._dp,0._dp],observed(:,:,i), &
            &residual,ignored_objective)
        call workspace%pose_objective_gradient(rotmats(:,:,i),[0._dp,0._dp],observed(:,:,i), &
            &particle_objective,gradient)
        objective = objective+particle_objective
        rotation_gradient_rms = rotation_gradient_rms+sum(gradient(1:3)**2)
        shift_gradient_rms = shift_gradient_rms+sum(gradient(4:5)**2)
        do k = lbound(observed,2), ubound(observed,2)
            do h = lbound(observed,1), ubound(observed,1)
                if( h*h+k*k < 4 .or. h*h+k*k > (BOX/2)**2 ) cycle
                model_value = cmplx(model(h,k),kind=dp)
                observed_value = cmplx(observed(h,k,i),kind=dp)
                model_norm = model_norm+abs(model_value)**2
                data_norm = data_norm+abs(observed_value)**2
                cross_term = cross_term+real(conjg(model_value)*observed_value,dp)
                residual_norm = residual_norm+abs(cmplx(residual(h,k),kind=dp))**2
                radius = nint(sqrt(real(h*h+k*k,dp)))
                band = min(NBANDS,max(1,(NBANDS*radius-1)/(BOX/2)+1))
                band_model(band) = band_model(band)+abs(model_value)**2
                band_data(band) = band_data(band)+abs(observed_value)**2
                band_cross(band) = band_cross(band)+real(conjg(model_value)*observed_value,dp)
            enddo
        enddo
    enddo
    relative_residual = sqrt(residual_norm/max(data_norm,tiny(data_norm)))
    amplitude_scale = cross_term/max(model_norm,tiny(model_norm))
    fitted_squared = max(0._dp,data_norm+amplitude_scale**2*model_norm- &
        &2._dp*amplitude_scale*cross_term)
    fitted_residual = sqrt(fitted_squared/max(data_norm,tiny(data_norm)))
    do band = 1, NBANDS
        fitted_squared = max(0._dp,band_data(band)+amplitude_scale**2*band_model(band)- &
            &2._dp*amplitude_scale*band_cross(band))
        band_residuals(band) = sqrt(fitted_squared/max(band_data(band),tiny(band_data(band))))
        band_correlations(band) = band_cross(band)/sqrt(max(band_model(band)*band_data(band), &
            &tiny(band_model(band))))
    enddo
    rotation_gradient_rms = sqrt(rotation_gradient_rms/real(size(rotmats,3),dp))
    shift_gradient_rms = sqrt(shift_gradient_rms/real(size(rotmats,3),dp))
end subroutine measure_exact_pose_fit

!> Print all evidence for one fixed-reference arm before any scientific gate.
subroutine print_diagnostic_arm(name,exact_objective,forward_residual,fitted_residual, &
    &amplitude_scale,rotation_gradient,shift_gradient,band_residuals,band_correlations, &
    &final_rotation_rms,final_shift_rms,summary)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: exact_objective, forward_residual, fitted_residual, amplitude_scale
    real(dp), intent(in) :: rotation_gradient, shift_gradient
    real(dp), intent(in) :: band_residuals(NBANDS), band_correlations(NBANDS)
    real(dp), intent(in) :: final_rotation_rms, final_shift_rms
    type(pcg_pose_polish_summary), intent(in) :: summary

    write(*,'(a,1x,a,4(1x,es14.6))') &
        &'CONTINUOUS_3D_FIXED_REFERENCE exact objective/raw-residual/fitted-residual/scale', &
        &trim(name),exact_objective,forward_residual,fitted_residual,amplitude_scale
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FIXED_REFERENCE fitted low/mid/high residual',trim(name),band_residuals
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FIXED_REFERENCE low/mid/high correlation',trim(name),band_correlations
    write(*,'(a,1x,a,2(1x,es14.6))') 'CONTINUOUS_3D_FIXED_REFERENCE exact rotation/shift gradient RMS', &
        &trim(name),rotation_gradient,shift_gradient
    write(*,'(a,1x,a,4(1x,es14.6))') 'CONTINUOUS_3D_FIXED_REFERENCE final rot/shift/objective before/after', &
        &trim(name),final_rotation_rms,final_shift_rms,summary%objective_before,summary%objective_after
    write(*,'(a,1x,a,6(1x,i0))') 'CONTINUOUS_3D_FIXED_REFERENCE terminal improved/unchanged/unreliable/bound/invalid/limit', &
        &trim(name),summary%nimproved,summary%nunchanged,summary%nunreliable,summary%nstep_bound, &
        &summary%ninvalid,summary%niteration_limit
end subroutine print_diagnostic_arm

!> Print a deterministic diagnosis ladder without turning open hypotheses into pass gates.
subroutine print_diagnostic_decision(rotation_errors,shift_errors,fitted_residuals)
    real(dp), intent(in) :: rotation_errors(NARMS), shift_errors(NARMS), fitted_residuals(NARMS)
    logical :: exact_kb_stationary, matched_stationary
    logical :: simulator_scaled_raw_stationary, simulator_scaled_corrected_stationary

    matched_stationary = rotation_errors(ARM_MATCHED_TRUTH) <= EXACT_ROTATION_TOLERANCE .and. &
        &shift_errors(ARM_MATCHED_TRUTH) <= EXACT_SHIFT_TOLERANCE
    exact_kb_stationary = rotation_errors(ARM_EXACT_KB_SCALED) <= EXACT_ROTATION_TOLERANCE .and. &
        &shift_errors(ARM_EXACT_KB_SCALED) <= EXACT_SHIFT_TOLERANCE
    simulator_scaled_raw_stationary = &
        &rotation_errors(ARM_SIMULATOR_SCALED_RAW) <= EXACT_ROTATION_TOLERANCE .and. &
        &shift_errors(ARM_SIMULATOR_SCALED_RAW) <= EXACT_SHIFT_TOLERANCE
    simulator_scaled_corrected_stationary = &
        &rotation_errors(ARM_SIMULATOR_SCALED_CORRECTED) <= EXACT_ROTATION_TOLERANCE .and. &
        &shift_errors(ARM_SIMULATOR_SCALED_CORRECTED) <= EXACT_SHIFT_TOLERANCE
    write(*,'(a,4(1x,l1))') &
        &'CONTINUOUS_3D_FIXED_REFERENCE stationary matched/exact-kb/simulator-raw/simulator-corrected', &
        &matched_stationary,exact_kb_stationary,simulator_scaled_raw_stationary, &
        &simulator_scaled_corrected_stationary
    write(*,'(a,5(1x,es14.6))') &
        &'CONTINUOUS_3D_FIXED_REFERENCE fitted residual exact-kb/sim-raw/sim-corrected/raw-recon/corrected-recon', &
        &fitted_residuals(ARM_EXACT_KB_SCALED),fitted_residuals(ARM_SIMULATOR_SCALED_RAW), &
        &fitted_residuals(ARM_SIMULATOR_SCALED_CORRECTED), &
        &fitted_residuals(ARM_SIMULATOR_RAW_RECONSTRUCTION), &
        &fitted_residuals(ARM_SIMULATOR_CORRECTED_RECONSTRUCTION)
    if( .not. matched_stationary )then
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: LM_OR_JACOBIAN_MISMATCH'
    else if( .not. exact_kb_stationary )then
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: EXACT_VERSUS_FAST_PADDED_GATHER_MISMATCH'
    else if( simulator_scaled_raw_stationary .and. .not. simulator_scaled_corrected_stationary )then
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: SIMULATOR_MATCHES_RAW_MODEL_NOT_INVERSE_ENVELOPE'
    else if( .not. simulator_scaled_raw_stationary )then
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: SIMULATOR_PADDING_OR_CLIP_MISMATCH'
    else if( rotation_errors(ARM_SIMULATOR_RAW_RECONSTRUCTION) > EXACT_ROTATION_TOLERANCE .or. &
        &shift_errors(ARM_SIMULATOR_RAW_RECONSTRUCTION) > EXACT_SHIFT_TOLERANCE )then
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: RECONSTRUCTION_REFERENCE_BIAS'
    else
        write(*,'(a)') 'CONTINUOUS_3D_FIXED_REFERENCE DIAGNOSIS: NO_FIXED_REFERENCE_DRIFT'
    endif
end subroutine print_diagnostic_decision

!> Count all mutually exclusive pose-polishing terminal outcomes.
pure integer function terminal_count(summary) result(total)
    type(pcg_pose_polish_summary), intent(in) :: summary

    total = summary%nimproved+summary%nunchanged+summary%nunreliable+summary%nstep_bound+ &
        &summary%ninvalid+summary%niteration_limit
end function terminal_count

!> Return the rotation-error RMS over one matched batch.
pure real(dp) function rotation_rms(reference,estimate) result(error)
    real(dp), intent(in) :: reference(:,:,:), estimate(:,:,:)
    real(dp) :: cosine, squared_error
    integer :: i

    squared_error = 0._dp
    do i = 1, size(reference,3)
        cosine = 0.5_dp*(trace3(matmul(transpose(reference(:,:,i)),estimate(:,:,i)))-1._dp)
        squared_error = squared_error+acos(max(-1._dp,min(1._dp,cosine)))**2
    enddo
    error = sqrt(squared_error/real(size(reference,3),dp))
end function rotation_rms

!> Return the trace of one three-by-three matrix.
pure real(dp) function trace3(matrix) result(trace)
    real(dp), intent(in) :: matrix(3,3)

    trace = matrix(1,1)+matrix(2,2)+matrix(3,3)
end function trace3

end module continuous_3D_pcg_refinement_fixed_reference_test
