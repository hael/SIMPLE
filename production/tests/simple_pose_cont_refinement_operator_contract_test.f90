module pose_cont_refinement_operator_contract_test
use pose_cont_refinement_operator_contract_support, only: &
    &CONTRACT_MODEL_BARE, CONTRACT_MODEL_DEAPOD, CONTRACT_MODEL_NAMES, CONTRACT_MODEL_WINDOW, &
    &CONTRACT_MODEL_WINDOW_DEAPOD, &
    &CONTRACT_NMODELS, CONTRACT_SUPPORT_RADIUS, build_contract_observations, build_contract_orientations, &
    &build_contract_truth, build_envelope_free_observations, evaluate_contract_model, metrics_are_finite, &
    &operator_contract_metrics, reconstruct_contract_reference
use pose_cont_refinement_test_helpers, only: assert_true
use simple_defs, only: dp
use simple_oris, only: oris
implicit none
private

integer, parameter :: NBOXES = 3
integer, parameter :: NITERATIONS = 3
integer, parameter :: NPROJECTIONS = 48
integer, parameter :: NPROBES = 8
integer, parameter :: BOXES(NBOXES) = [24,32,48]
integer, parameter :: ITERATION_COUNTS(NITERATIONS) = [2,4,8]
real(dp), parameter :: MATCHED_TRUTH_TOLERANCE = 1.e-3_dp
real(dp), parameter :: ENVELOPE_FREE_TOLERANCE = 1.e-5_dp
real(dp), parameter :: ENVELOPE_EFFECT_MINIMUM = 1.e-3_dp

public :: run_operator_contract_diagnostic

contains

!> Separate inverse-envelope, finite-box, taper, iteration, and support-margin effects.
subroutine run_operator_contract_diagnostic()
    type(operator_contract_metrics) :: clip_truth(NBOXES,CONTRACT_NMODELS)
    type(operator_contract_metrics) :: envelope_truth(NBOXES,2)
    type(operator_contract_metrics) :: tapered_truth(NBOXES,CONTRACT_NMODELS)
    type(operator_contract_metrics) :: references(NBOXES,NITERATIONS,CONTRACT_NMODELS)
    type(oris) :: orientations
    complex, allocatable :: clip_observed(:,:,:), envelope_observed(:,:,:), tapered_observed(:,:,:)
    real, allocatable :: truth(:,:,:), reference(:,:,:)
    integer :: ibox, iteration, model, niters
    logical :: deapod_contract_passes, deapod_improves_reference, margin_reduces_bare
    logical :: window_matches_clip, window_matches_taper

    call build_contract_orientations(NPROJECTIONS,orientations)
    do ibox = 1, NBOXES
        call build_contract_truth(BOXES(ibox),truth)
        call build_contract_observations(truth,orientations,.false.,clip_observed)
        call build_contract_observations(truth,orientations,.true.,tapered_observed)
        call build_envelope_free_observations(truth,orientations,envelope_observed)

        call evaluate_contract_model(truth,orientations,envelope_observed,NPROBES, &
            &CONTRACT_MODEL_BARE,.false.,envelope_truth(ibox,1))
        call print_arm(BOXES(ibox),'truth_envelope_free',0,CONTRACT_MODEL_NAMES(CONTRACT_MODEL_BARE), &
            &envelope_truth(ibox,1))
        call evaluate_contract_model(truth,orientations,envelope_observed,NPROBES, &
            &CONTRACT_MODEL_DEAPOD,.false.,envelope_truth(ibox,2))
        call print_arm(BOXES(ibox),'truth_envelope_free',0,CONTRACT_MODEL_NAMES(CONTRACT_MODEL_DEAPOD), &
            &envelope_truth(ibox,2))

        do model = 1, CONTRACT_NMODELS
            call evaluate_contract_model(truth,orientations,clip_observed,NPROBES,model,.false., &
                &clip_truth(ibox,model))
            call print_arm(BOXES(ibox),'truth_clip',0,CONTRACT_MODEL_NAMES(model), &
                &clip_truth(ibox,model))
            call evaluate_contract_model(truth,orientations,tapered_observed,NPROBES,model,.true., &
                &tapered_truth(ibox,model))
            call print_arm(BOXES(ibox),'truth_taper',0,CONTRACT_MODEL_NAMES(model), &
                &tapered_truth(ibox,model))
        enddo

        do iteration = 1, NITERATIONS
            call reconstruct_contract_reference(orientations,tapered_observed,CONTRACT_SUPPORT_RADIUS, &
                &ITERATION_COUNTS(iteration),reference,niters)
            write(*,'(a,3(1x,i0),1x,f6.2)') &
                &'CONTINUOUS_3D_OPERATOR_CONTRACT reconstruction box/requested/completed/support', &
                &BOXES(ibox),ITERATION_COUNTS(iteration),niters,CONTRACT_SUPPORT_RADIUS
            do model = 1, CONTRACT_NMODELS
                call evaluate_contract_model(reference,orientations,tapered_observed,NPROBES, &
                    &model,.true.,references(ibox,iteration,model))
                call print_arm(BOXES(ibox),'reconstructed',ITERATION_COUNTS(iteration), &
                    &CONTRACT_MODEL_NAMES(model),references(ibox,iteration,model))
            enddo
            deallocate(reference)
        enddo

        deallocate(clip_observed,envelope_observed,tapered_observed,truth)
    enddo
    call orientations%kill()

    do ibox = 1, NBOXES
        do model = 1, CONTRACT_NMODELS
            call assert_true(metrics_are_finite(clip_truth(ibox,model)), &
                &'operator contract clip-truth arm is not finite')
            call assert_true(metrics_are_finite(tapered_truth(ibox,model)), &
                &'operator contract tapered-truth arm is not finite')
            do iteration = 1, NITERATIONS
                call assert_true(metrics_are_finite(references(ibox,iteration,model)), &
                    &'operator contract reconstructed-reference arm is not finite')
            enddo
        enddo
        call assert_true(clip_truth(ibox,CONTRACT_MODEL_WINDOW)%fitted_residual < &
            &MATCHED_TRUTH_TOLERANCE, &
            &'finite-box model does not match the independent clipped truth')
        call assert_true(tapered_truth(ibox,CONTRACT_MODEL_WINDOW)%fitted_residual < &
            &MATCHED_TRUTH_TOLERANCE, &
            &'finite-box model does not match the tapered independent truth')
        call assert_true(metrics_are_finite(envelope_truth(ibox,1)) .and. &
            &metrics_are_finite(envelope_truth(ibox,2)), &
            &'operator contract envelope-free arm is not finite')
        call assert_true(envelope_truth(ibox,2)%fitted_residual < ENVELOPE_FREE_TOLERANCE, &
            &'inverse-envelope model does not match envelope-free observations')
        call assert_true(envelope_truth(ibox,1)%fitted_residual > ENVELOPE_EFFECT_MINIMUM, &
            &'envelope-free control does not distinguish the bare gather')
    enddo

    window_matches_clip = all(clip_truth(:,CONTRACT_MODEL_WINDOW)%fitted_residual < &
        &MATCHED_TRUTH_TOLERANCE)
    window_matches_taper = all(tapered_truth(:,CONTRACT_MODEL_WINDOW)%fitted_residual < &
        &MATCHED_TRUTH_TOLERANCE)
    deapod_contract_passes = all(envelope_truth(:,2)%fitted_residual < &
        &ENVELOPE_FREE_TOLERANCE) .and. all(envelope_truth(:,1)%fitted_residual > &
        &ENVELOPE_EFFECT_MINIMUM)
    deapod_improves_reference = minval(references(:,:,CONTRACT_MODEL_WINDOW_DEAPOD)%fitted_residual) < &
        &0.9_dp*minval(references(:,:,CONTRACT_MODEL_WINDOW)%fitted_residual)
    margin_reduces_bare = clip_truth(NBOXES,CONTRACT_MODEL_BARE)%fitted_residual < &
        &0.5_dp*clip_truth(1,CONTRACT_MODEL_BARE)%fitted_residual
    write(*,'(a,4(1x,l1))') &
        &'CONTINUOUS_3D_OPERATOR_CONTRACT components clip/taper/envelope-control/margin', &
        &window_matches_clip,window_matches_taper,deapod_contract_passes,margin_reduces_bare
    write(*,'(a,1x,l1)') &
        &'CONTINUOUS_3D_OPERATOR_CONTRACT reconstructed-reference deapod-improves-10pct', &
        &deapod_improves_reference
    write(*,'(a)') 'CONTINUOUS_3D_OPERATOR_CONTRACT DIAGNOSIS: FOUR_MODEL_MATRIX_COMPLETE'
    write(*,'(a)') 'CONTINUOUS_3D_OPERATOR_CONTRACT: EVIDENCE COMPLETE'
end subroutine run_operator_contract_diagnostic

!> Print one complete model/reference arm in a machine-searchable record.
subroutine print_arm(box,reference_name,maxits,model_name,metrics)
    integer, intent(in) :: box, maxits
    character(len=*), intent(in) :: reference_name, model_name
    type(operator_contract_metrics), intent(in) :: metrics

    write(*,'(a,1x,i0,1x,a,1x,i0,1x,a,5(1x,es14.6))') &
        &'CONTINUOUS_3D_OPERATOR_CONTRACT arm',box,trim(reference_name),maxits,trim(model_name), &
        &metrics%fitted_residual,metrics%amplitude_scale,metrics%objective, &
        &metrics%rotation_gradient_rms,metrics%shift_gradient_rms
end subroutine print_arm

end module pose_cont_refinement_operator_contract_test
