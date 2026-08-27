module pose_cont_refinement_forward_path_test
use pose_cont_refinement_fixed_reference_test, only: measure_exact_pose_fit
use pose_cont_refinement_reference_support, only: build_disjoint_half_orientations, &
    &HALFSET_LAMBDA, HALFSET_SMPD
use pose_cont_refinement_test_helpers, only: assert_true, build_truth_volume, &
    &BOX => TRUTH_VOLUME_BOX
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: dp, OSMPL_PAD_FAC
use simple_image, only: image
use simple_ori, only: ori
use simple_oris, only: oris
use simple_projector, only: projector
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner
use simple_reconstructor_pcg, only: reconstructor_pcg
implicit none
private
public :: run_forward_path_diagnostic, build_forward_path_observations

integer, parameter :: NPARTICLES = 48
integer, parameter :: NBANDS = 3
integer, parameter :: NSTAGES = 3
integer, parameter :: NTRANSITIONS = NSTAGES-1
real(dp), parameter :: BOUNDARY_RESIDUAL_TOLERANCE = 1.e-3_dp
character(len=24), parameter :: STAGE_NAMES(NSTAGES) = [character(len=24) :: &
    &'padded_projected', 'padded_roundtrip', 'clipped_native']
character(len=32), parameter :: TRANSITION_NAMES(NTRANSITIONS) = [character(len=32) :: &
    &'projected_to_roundtrip', 'roundtrip_to_clipped']

contains

!> Locate the first simulator operation that disagrees with the raw-truth PCG gather.
subroutine run_forward_path_diagnostic()
    type(ori) :: orientation
    type(oris) :: even_oris, odd_oris
    type(reconstructor_pcg) :: model_operator
    type(cartesian_pose_refiner) :: raw_truth_workspace
    complex, allocatable :: stage_planes(:,:,:,:)
    real, allocatable :: truth_volume(:,:,:)
    real(dp), allocatable :: exact_rotmats(:,:,:)
    integer, allocatable :: even_ids(:), odd_ids(:)
    real(dp) :: objectives(NSTAGES), raw_residuals(NSTAGES), fitted_residuals(NSTAGES)
    real(dp) :: amplitude_scales(NSTAGES), rotation_gradients(NSTAGES), shift_gradients(NSTAGES)
    real(dp) :: band_residuals(NBANDS,NSTAGES), band_correlations(NBANDS,NSTAGES)
    real(dp) :: transition_residuals(NTRANSITIONS), transition_scales(NTRANSITIONS)
    real(dp) :: transition_band_residuals(NBANDS,NTRANSITIONS)
    real(dp) :: transition_band_correlations(NBANDS,NTRANSITIONS)
    integer :: i, stage, transition

    call build_disjoint_half_orientations(NPARTICLES,even_oris,odd_oris,even_ids,odd_ids)
    call build_truth_volume(truth_volume)
    call model_operator%new(BOX,HALFSET_SMPD,HALFSET_LAMBDA)
    call model_operator%set_volume(truth_volume)
    call raw_truth_workspace%new_prepared_test(truth_volume)
    call raw_truth_workspace%set_shell_range([2,BOX/2])

    allocate(exact_rotmats(3,3,NPARTICLES))
    call orientation%new(.false.)
    do i = 1, NPARTICLES
        call even_oris%get_ori(i,orientation)
        exact_rotmats(:,:,i) = real(orientation%get_mat(),dp)
    enddo
    call build_forward_path_observations(model_operator,even_oris,truth_volume,stage_planes)

    do stage = 1, NSTAGES
        call measure_exact_pose_fit(raw_truth_workspace,exact_rotmats,stage_planes(:,:,:,stage), &
            &objectives(stage),raw_residuals(stage),fitted_residuals(stage), &
            &amplitude_scales(stage),rotation_gradients(stage),shift_gradients(stage), &
            &band_residuals(:,stage),band_correlations(:,stage))
        call print_stage_metrics(STAGE_NAMES(stage),objectives(stage),raw_residuals(stage), &
            &fitted_residuals(stage),amplitude_scales(stage),rotation_gradients(stage), &
            &shift_gradients(stage),band_residuals(:,stage),band_correlations(:,stage))
    enddo

    do transition = 1, NTRANSITIONS
        call measure_plane_transition(stage_planes(:,:,:,transition), &
            &stage_planes(:,:,:,transition+1),transition_residuals(transition), &
            &transition_scales(transition),transition_band_residuals(:,transition), &
            &transition_band_correlations(:,transition))
        call print_transition_metrics(TRANSITION_NAMES(transition), &
            &transition_residuals(transition),transition_scales(transition), &
            &transition_band_residuals(:,transition), &
            &transition_band_correlations(:,transition))
    enddo
    call print_forward_path_decision(fitted_residuals,transition_residuals)

    call assert_true(all(ieee_is_finite([objectives,raw_residuals,fitted_residuals, &
        &amplitude_scales,rotation_gradients,shift_gradients, &
        &reshape(band_residuals,[NBANDS*NSTAGES]), &
        &reshape(band_correlations,[NBANDS*NSTAGES]),transition_residuals, &
        &transition_scales,reshape(transition_band_residuals,[NBANDS*NTRANSITIONS]), &
        &reshape(transition_band_correlations,[NBANDS*NTRANSITIONS])])), &
        &'forward-path diagnostic produced non-finite evidence')
    write(*,'(a)') 'CONTINUOUS_3D_FORWARD_PATH: EVIDENCE COMPLETE'

    call orientation%kill()
    call even_oris%kill()
    call odd_oris%kill()
    call raw_truth_workspace%kill()
    call model_operator%kill()
end subroutine run_forward_path_diagnostic

!> Capture simulator observations before round trip, after round trip, and after clipping.
subroutine build_forward_path_observations(sampler,orientations,volume,stage_planes)
    type(reconstructor_pcg), intent(in) :: sampler
    type(oris), intent(inout) :: orientations
    real, intent(in) :: volume(BOX,BOX,BOX)
    complex, allocatable, intent(out) :: stage_planes(:,:,:,:)
    type(image) :: source_volume, padded_projection, native_projection
    type(projector) :: padded_projector
    type(ori) :: orientation
    integer :: i, lims2(2,2), nprojs

    nprojs = orientations%get_noris()
    lims2 = sampler%get_lims2()
    allocate(stage_planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2), &
        &nprojs,NSTAGES),source=cmplx(0.,0.))
    call source_volume%new([BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%set_rmat(volume,.false.)
    call padded_projector%new(OSMPL_PAD_FAC*[BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%pad(padded_projector)
    call padded_projector%fft()
    call padded_projector%expand_cmat()
    call padded_projection%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], &
        &HALFSET_SMPD,wthreads=.false.)
    call native_projection%new([BOX,BOX,1],HALFSET_SMPD,wthreads=.false.)
    call orientation%new(.false.)
    do i = 1, nprojs
        call orientations%get_ori(i,orientation)
        call padded_projector%fproject_serial(orientation,padded_projection)
        ! Native frequencies occupy every pad-factor-th sample on the padded plane.
        call extract_decimated_plane(padded_projection,stage_planes(:,:,i,1))
        ! simimg enforces a real-image Fourier round trip before detector operations.
        call padded_projection%ifft()
        call padded_projection%fft()
        call extract_decimated_plane(padded_projection,stage_planes(:,:,i,2))
        ! Production then returns to real space, clips, and later FFTs the native image.
        call padded_projection%ifft()
        call padded_projection%clip(native_projection)
        call native_projection%fft()
        stage_planes(:,:,i,3) = sampler%extract_native_plane(native_projection)
    enddo
    call orientation%kill()
    call native_projection%kill()
    call padded_projection%kill()
    call padded_projector%kill_expanded()
    call padded_projector%kill()
    call source_volume%kill()
end subroutine build_forward_path_observations

!> Read native Fourier frequencies from their coincident padded-plane locations.
subroutine extract_decimated_plane(padded_projection,plane)
    type(image), intent(in) :: padded_projection
    complex, intent(out) :: plane(-BOX/2:,-BOX/2:)
    integer :: h, k

    plane = cmplx(0.,0.)
    do k = -BOX/2, BOX/2
        do h = -BOX/2, BOX/2
            if( h*h+k*k > (BOX/2)**2 ) cycle
            plane(h,k) = padded_projection%get_fcomp2D(OSMPL_PAD_FAC*h,OSMPL_PAD_FAC*k)
        enddo
    enddo
end subroutine extract_decimated_plane

!> Fit one global amplitude between adjacent observation stages and report shell evidence.
subroutine measure_plane_transition(reference,candidate,residual,scale,band_residuals,band_correlations)
    complex, intent(in) :: reference(-BOX/2:,-BOX/2:,:), candidate(-BOX/2:,-BOX/2:,:)
    real(dp), intent(out) :: residual, scale
    real(dp), intent(out) :: band_residuals(NBANDS), band_correlations(NBANDS)
    complex(dp) :: candidate_value, reference_value
    real(dp) :: candidate_norm, cross_term, fitted_squared, reference_norm
    real(dp) :: band_candidate(NBANDS), band_cross(NBANDS), band_reference(NBANDS)
    integer :: band, h, i, k, radius

    candidate_norm = 0._dp
    cross_term = 0._dp
    reference_norm = 0._dp
    band_candidate = 0._dp
    band_cross = 0._dp
    band_reference = 0._dp
    do i = 1, size(reference,3)
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                if( h*h+k*k < 4 .or. h*h+k*k > (BOX/2)**2 ) cycle
                reference_value = cmplx(reference(h,k,i),kind=dp)
                candidate_value = cmplx(candidate(h,k,i),kind=dp)
                reference_norm = reference_norm+abs(reference_value)**2
                candidate_norm = candidate_norm+abs(candidate_value)**2
                cross_term = cross_term+real(conjg(candidate_value)*reference_value,dp)
                radius = nint(sqrt(real(h*h+k*k,dp)))
                band = min(NBANDS,max(1,(NBANDS*radius-1)/(BOX/2)+1))
                band_reference(band) = band_reference(band)+abs(reference_value)**2
                band_candidate(band) = band_candidate(band)+abs(candidate_value)**2
                band_cross(band) = band_cross(band)+ &
                    &real(conjg(candidate_value)*reference_value,dp)
            enddo
        enddo
    enddo
    scale = cross_term/max(candidate_norm,tiny(candidate_norm))
    fitted_squared = max(0._dp,reference_norm+scale**2*candidate_norm-2._dp*scale*cross_term)
    residual = sqrt(fitted_squared/max(reference_norm,tiny(reference_norm)))
    do band = 1, NBANDS
        fitted_squared = max(0._dp,band_reference(band)+scale**2*band_candidate(band)- &
            &2._dp*scale*band_cross(band))
        band_residuals(band) = sqrt(fitted_squared/ &
            &max(band_reference(band),tiny(band_reference(band))))
        band_correlations(band) = band_cross(band)/sqrt(max( &
            &band_reference(band)*band_candidate(band),tiny(band_candidate(band))))
    enddo
end subroutine measure_plane_transition

!> Print PCG-fit evidence for one forward-path observation stage.
subroutine print_stage_metrics(name,objective,raw_residual,fitted_residual,scale, &
    &rotation_gradient,shift_gradient,band_residuals,band_correlations)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: objective, raw_residual, fitted_residual, scale
    real(dp), intent(in) :: rotation_gradient, shift_gradient
    real(dp), intent(in) :: band_residuals(NBANDS), band_correlations(NBANDS)

    write(*,'(a,1x,a,4(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH objective/raw-residual/fitted-residual/scale', &
        &trim(name),objective,raw_residual,fitted_residual,scale
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH fitted low/mid/high residual',trim(name),band_residuals
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH low/mid/high correlation',trim(name),band_correlations
    write(*,'(a,1x,a,2(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH rotation/shift gradient RMS', &
        &trim(name),rotation_gradient,shift_gradient
end subroutine print_stage_metrics

!> Print the fitted change introduced by one adjacent simulator operation.
subroutine print_transition_metrics(name,residual,scale,band_residuals,band_correlations)
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: residual, scale
    real(dp), intent(in) :: band_residuals(NBANDS), band_correlations(NBANDS)

    write(*,'(a,1x,a,2(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH transition fitted-residual/scale',trim(name),residual,scale
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH transition low/mid/high residual',trim(name),band_residuals
    write(*,'(a,1x,a,3(1x,es14.6))') &
        &'CONTINUOUS_3D_FORWARD_PATH transition low/mid/high correlation', &
        &trim(name),band_correlations
end subroutine print_transition_metrics

!> Label the first boundary above the predeclared fitted-residual tolerance.
subroutine print_forward_path_decision(stage_residuals,transition_residuals)
    real(dp), intent(in) :: stage_residuals(NSTAGES), transition_residuals(NTRANSITIONS)

    write(*,'(a,5(1x,es14.6))') 'CONTINUOUS_3D_FORWARD_PATH residuals direct/roundtrip/clipped/transitions', &
        &stage_residuals,transition_residuals
    if( stage_residuals(1) > BOUNDARY_RESIDUAL_TOLERANCE )then
        write(*,'(a)') 'CONTINUOUS_3D_FORWARD_PATH DIAGNOSIS: PADDED_PROJECTOR_VERSUS_PCG_GATHER'
    else if( transition_residuals(1) > BOUNDARY_RESIDUAL_TOLERANCE )then
        write(*,'(a)') 'CONTINUOUS_3D_FORWARD_PATH DIAGNOSIS: FOURIER_ROUNDTRIP'
    else if( transition_residuals(2) > BOUNDARY_RESIDUAL_TOLERANCE )then
        write(*,'(a)') 'CONTINUOUS_3D_FORWARD_PATH DIAGNOSIS: REALSPACE_CLIP_NATIVE_FFT'
    else
        write(*,'(a)') 'CONTINUOUS_3D_FORWARD_PATH DIAGNOSIS: NO_MEASURABLE_BOUNDARY_MISMATCH'
    endif
end subroutine print_forward_path_decision

end module pose_cont_refinement_forward_path_test
