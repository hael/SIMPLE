module continuous_3D_pcg_refinement_operator_contract_support
use ieee_arithmetic, only: ieee_is_finite
use simple_defs, only: COSMSKHALFWIDTH, OSMPL_PAD_FAC, PI, dp, sp
use simple_image, only: image
use simple_ori, only: ori
use simple_oris, only: oris
use simple_projector, only: projector
use simple_reconstructor_pcg, only: pcg_fourier_workspace, reconstructor_pcg, &
    &right_increment_rotation
implicit none
private

integer, parameter, public :: CONTRACT_MODEL_BARE = 1
integer, parameter, public :: CONTRACT_MODEL_DEAPOD = 2
integer, parameter, public :: CONTRACT_MODEL_WINDOW = 3
integer, parameter, public :: CONTRACT_MODEL_WINDOW_DEAPOD = 4
integer, parameter, public :: CONTRACT_NMODELS = 4
character(len=24), parameter, public :: CONTRACT_MODEL_NAMES(CONTRACT_NMODELS) = &
    &[character(len=24) :: 'bare', 'inverse_envelope', 'finite_box', &
      &'finite_box_envelope']
real, parameter :: CONTRACT_SMPD = 1.5
real, parameter :: CONTRACT_LAMBDA = 1.e-3
real, parameter, public :: CONTRACT_SUPPORT_RADIUS = 10.
real(dp), parameter :: ROTATION_FD_STEP = 5.e-4_dp
real(dp), parameter :: SHIFT_FD_STEP = 5.e-3_dp

type, public :: operator_contract_metrics
    real(dp) :: fitted_residual = 0._dp
    real(dp) :: amplitude_scale = 0._dp
    real(dp) :: objective = 0._dp
    real(dp) :: rotation_gradient_rms = 0._dp
    real(dp) :: shift_gradient_rms = 0._dp
end type operator_contract_metrics

type :: prediction_buffers
    type(image) :: padded
    type(image) :: native
    logical :: exists = .false.
end type prediction_buffers

public :: build_contract_observations
public :: build_contract_orientations
public :: build_contract_truth
public :: build_envelope_free_observations
public :: evaluate_contract_model
public :: metrics_are_finite
public :: reconstruct_contract_reference

contains

!> Build one deterministic orientation set shared by every operator arm.
subroutine build_contract_orientations(nprojections,orientations)
    integer, intent(in) :: nprojections
    type(oris), intent(inout) :: orientations

    call orientations%new(nprojections,.false.)
    call orientations%spiral()
end subroutine build_contract_orientations

!> Embed the established asymmetric Gaussian truth in a requested even box.
subroutine build_contract_truth(box,volume)
    integer, intent(in) :: box
    real, allocatable, intent(out) :: volume(:,:,:)
    real, parameter :: centres(3,4) = reshape([&
        &-5.0,-3.0, 2.0, 4.0, 5.0,-3.0, 0.0,-6.0,-5.0, 3.0,-2.0, 6.0],[3,4])
    real, parameter :: sigmas(4) = [2.0,2.5,1.8,2.2]
    real, parameter :: amplitudes(4) = [1.0,0.8,0.6,0.5]
    type(image) :: support_image
    real :: centre, dx, dy, dz
    integer :: blob, i, j, k

    if( box < 24 .or. mod(box,2) /= 0 ) error stop 'operator contract requires an even box of at least 24'
    allocate(volume(box,box,box),source=0.)
    centre = real(box)/2.+0.5
    do k = 1, box
        do j = 1, box
            do i = 1, box
                do blob = 1, size(amplitudes)
                    dx = real(i)-centre-centres(1,blob)
                    dy = real(j)-centre-centres(2,blob)
                    dz = real(k)-centre-centres(3,blob)
                    volume(i,j,k) = volume(i,j,k)+amplitudes(blob)* &
                        &exp(-(dx*dx+dy*dy+dz*dz)/(2.*sigmas(blob)**2))
                enddo
            enddo
        enddo
    enddo
    ! Keep object support fixed while the distance to the particle-box edge changes.
    call support_image%new([box,box,box],CONTRACT_SMPD,wthreads=.false.)
    call support_image%set_rmat(volume,.false.)
    call support_image%mask3D_soft(CONTRACT_SUPPORT_RADIUS,backgr=0.)
    volume = support_image%get_rmat()
    call support_image%kill()
end subroutine build_contract_truth

!> Generate independent padded-projector observations, with optional production edge taper.
subroutine build_contract_observations(volume,orientations,apply_taper,observed)
    real, intent(in) :: volume(:,:,:)
    type(oris), intent(inout) :: orientations
    logical, intent(in) :: apply_taper
    complex, allocatable, intent(out) :: observed(:,:,:)
    type(image) :: source_volume, padded_projection, native_projection
    type(projector) :: padded_projector
    type(reconstructor_pcg) :: sampler
    type(ori) :: orientation
    real :: edge_mean
    integer :: box, i, lims2(2,2), nprojections, padded_box

    box = size(volume,1)
    if( any(shape(volume) /= [box,box,box]) ) error stop 'operator contract truth must be cubic'
    padded_box = OSMPL_PAD_FAC*box
    nprojections = orientations%get_noris()
    call sampler%new(box,CONTRACT_SMPD,CONTRACT_LAMBDA)
    lims2 = sampler%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojections), &
        &source=cmplx(0.,0.))
    call source_volume%new([box,box,box],CONTRACT_SMPD,wthreads=.false.)
    call source_volume%set_rmat(volume,.false.)
    call padded_projector%new([padded_box,padded_box,padded_box],CONTRACT_SMPD,wthreads=.false.)
    call source_volume%pad(padded_projector)
    call padded_projector%fft()
    call padded_projector%expand_cmat()
    call padded_projection%new([padded_box,padded_box,1],CONTRACT_SMPD,wthreads=.false.)
    call native_projection%new([box,box,1],CONTRACT_SMPD,wthreads=.false.)
    call orientation%new(.false.)
    do i = 1, nprojections
        call orientations%get_ori(i,orientation)
        call padded_projector%fproject_serial(orientation,padded_projection)
        ! Match simulate_particles through the padded round trip and native clip.
        call padded_projection%ifft()
        call padded_projection%fft()
        call padded_projection%ifft()
        call padded_projection%clip(native_projection)
        if( apply_taper )then
            edge_mean = 0.
            call native_projection%taper_edges_particle(nint(COSMSKHALFWIDTH),edge_mean)
        endif
        call native_projection%fft()
        observed(:,:,i) = sampler%extract_native_plane(native_projection)
    enddo
    call orientation%kill()
    call native_projection%kill()
    call padded_projection%kill()
    call padded_projector%kill_expanded()
    call padded_projector%kill()
    call source_volume%kill()
    call sampler%kill()
end subroutine build_contract_observations

!> Build the established PCG Stage 8 control by cancelling the gather envelope.
subroutine build_envelope_free_observations(volume,orientations,observed)
    real, intent(in) :: volume(:,:,:)
    type(oris), intent(inout) :: orientations
    complex, allocatable, intent(out) :: observed(:,:,:)
    type(reconstructor_pcg) :: sampler
    type(ori) :: orientation
    real, allocatable :: deapodized_volume(:,:,:)
    integer :: box, i, lims2(2,2), nprojections

    box = size(volume,1)
    if( any(shape(volume) /= [box,box,box]) ) error stop 'envelope-free truth must be cubic'
    nprojections = orientations%get_noris()
    call sampler%new(box,CONTRACT_SMPD,CONTRACT_LAMBDA)
    allocate(deapodized_volume,source=volume*sampler%get_invenv())
    call sampler%set_volume(deapodized_volume)
    lims2 = sampler%get_lims2()
    allocate(observed(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojections), &
        &source=cmplx(0.,0.))
    call orientation%new(.false.)
    do i = 1, nprojections
        call orientations%get_ori(i,orientation)
        ! G(E^-1 V) is the policy's envelope-free target for this isolation arm.
        call sampler%forward_plane(orientation,observed(:,:,i))
    enddo
    call orientation%kill()
    call sampler%kill()
end subroutine build_envelope_free_observations

!> Reconstruct one fixed reference with the current matrix-free PCG contract.
subroutine reconstruct_contract_reference(orientations,observed,mask_radius,maxits,reference,niters)
    type(oris), intent(inout) :: orientations
    complex, intent(in) :: observed(:,:,:)
    real, intent(in) :: mask_radius
    integer, intent(in) :: maxits
    real, allocatable, intent(out) :: reference(:,:,:)
    integer, intent(out) :: niters
    type(reconstructor_pcg) :: pcgop
    integer :: box

    box = size(observed,1)-1
    call pcgop%new(box,CONTRACT_SMPD,CONTRACT_LAMBDA)
    if( mask_radius > 0. ) call pcgop%set_mask(mask_radius)
    call pcgop%prep_particles(orientations,use_ctf=.false.)
    call pcgop%build_operators(.false.)
    allocate(reference(box,box,box),source=0.)
    call pcgop%solve(observed,reference,maxits=maxits,rtol=0.,niters=niters)
    call pcgop%kill()
end subroutine reconstruct_contract_reference

!> Measure fit and exact-pose stationarity for one declared forward model.
subroutine evaluate_contract_model(reference,orientations,observed,nprobes,model_kind,apply_taper,metrics)
    real, intent(in) :: reference(:,:,:)
    type(oris), intent(inout) :: orientations
    complex, intent(in) :: observed(:,:,:)
    integer, intent(in) :: nprobes, model_kind
    logical, intent(in) :: apply_taper
    type(operator_contract_metrics), intent(out) :: metrics
    type(reconstructor_pcg) :: model_operator
    type(pcg_fourier_workspace) :: workspace
    type(prediction_buffers) :: buffers
    type(ori) :: orientation
    complex, allocatable :: models(:,:,:), minus_plane(:,:), plus_plane(:,:)
    real, allocatable :: model_volume(:,:,:)
    real(dp), allocatable :: rotmats(:,:,:)
    real(dp) :: data_norm, gradient, minus_objective, plus_objective
    real(dp) :: minus_rotmat(3,3), plus_rotmat(3,3), shift(2)
    integer :: axis, box, i, lims2(2,2)

    box = size(reference,1)
    if( any(shape(reference) /= [box,box,box]) ) error stop 'operator contract reference must be cubic'
    if( nprobes < 1 .or. nprobes > size(observed,3) ) error stop 'invalid operator contract probe count'
    if( model_kind < 1 .or. model_kind > CONTRACT_NMODELS ) error stop 'invalid operator contract model'
    call model_operator%new(box,CONTRACT_SMPD,CONTRACT_LAMBDA)
    allocate(model_volume,source=reference)
    if( model_uses_deapod(model_kind) ) model_volume = model_volume*model_operator%get_invenv()
    call model_operator%set_volume(model_volume)
    call model_operator%begin_fourier_workspace(workspace)
    call workspace%set_shell_range([2,box/2])
    lims2 = model_operator%get_lims2()
    allocate(models(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprobes),source=cmplx(0.,0.))
    allocate(minus_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(plus_plane(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2)))
    allocate(rotmats(3,3,nprobes))
    if( model_uses_window(model_kind) ) call initialize_buffers(box,buffers)
    call orientation%new(.false.)
    do i = 1, nprobes
        call orientations%get_ori(i,orientation)
        rotmats(:,:,i) = real(orientation%get_mat(),dp)
        call predict_plane(box,workspace,model_operator,rotmats(:,:,i),[0._dp,0._dp], &
            &model_uses_window(model_kind),apply_taper,buffers,models(:,:,i))
    enddo
    call fit_amplitude(box,models,observed(:,:,1:nprobes),metrics%amplitude_scale, &
        &metrics%fitted_residual,metrics%objective,data_norm)

    metrics%rotation_gradient_rms = 0._dp
    metrics%shift_gradient_rms = 0._dp
    do i = 1, nprobes
        do axis = 1, 5
            shift = 0._dp
            minus_rotmat = rotmats(:,:,i)
            plus_rotmat = rotmats(:,:,i)
            if( axis <= 3 )then
                minus_rotmat = right_increment_rotation(rotmats(:,:,i), &
                    &axis_vector(axis,-ROTATION_FD_STEP))
                plus_rotmat = right_increment_rotation(rotmats(:,:,i), &
                    &axis_vector(axis, ROTATION_FD_STEP))
            else
                shift(axis-3) = -SHIFT_FD_STEP
            endif
            call predict_plane(box,workspace,model_operator,minus_rotmat,shift, &
                &model_uses_window(model_kind),apply_taper,buffers,minus_plane)
            if( axis > 3 ) shift(axis-3) = SHIFT_FD_STEP
            call predict_plane(box,workspace,model_operator,plus_rotmat,shift, &
                &model_uses_window(model_kind),apply_taper,buffers,plus_plane)
            minus_objective = plane_objective(box,minus_plane,observed(:,:,i), &
                &metrics%amplitude_scale,data_norm)
            plus_objective = plane_objective(box,plus_plane,observed(:,:,i), &
                &metrics%amplitude_scale,data_norm)
            if( axis <= 3 )then
                gradient = (plus_objective-minus_objective)/(2._dp*ROTATION_FD_STEP)
                metrics%rotation_gradient_rms = metrics%rotation_gradient_rms+gradient*gradient
            else
                gradient = (plus_objective-minus_objective)/(2._dp*SHIFT_FD_STEP)
                metrics%shift_gradient_rms = metrics%shift_gradient_rms+gradient*gradient
            endif
        enddo
    enddo
    metrics%rotation_gradient_rms = sqrt(metrics%rotation_gradient_rms/real(3*nprobes,dp))
    metrics%shift_gradient_rms = sqrt(metrics%shift_gradient_rms/real(2*nprobes,dp))

    call orientation%kill()
    call kill_buffers(buffers)
    call workspace%kill()
    call model_operator%kill()
end subroutine evaluate_contract_model

!> Return true only when every reported diagnostic scalar is finite.
pure logical function metrics_are_finite(metrics) result(is_finite)
    type(operator_contract_metrics), intent(in) :: metrics

    is_finite = all(ieee_is_finite([metrics%fitted_residual,metrics%amplitude_scale, &
        &metrics%objective,metrics%rotation_gradient_rms,metrics%shift_gradient_rms]))
end function metrics_are_finite

!> Generate one raw or finite-box prediction at a requested local pose.
subroutine predict_plane(box,workspace,sampler,rotmat,shift,use_window,apply_taper,buffers,plane)
    integer, intent(in) :: box
    type(pcg_fourier_workspace), intent(in) :: workspace
    type(reconstructor_pcg), intent(in) :: sampler
    real(dp), intent(in) :: rotmat(3,3), shift(2)
    logical, intent(in) :: use_window, apply_taper
    type(prediction_buffers), intent(inout) :: buffers
    complex, intent(out) :: plane(-box/2:,-box/2:)
    complex, allocatable :: zero_plane(:,:)
    complex :: dvalue_dloc(3), phase, value
    real(sp) :: loc(3), switch_margin(3)
    real(dp) :: arg, ignored_objective
    real :: edge_mean
    integer :: h, k, padded_box, padded_radius, phys(3)

    if( .not. use_window )then
        allocate(zero_plane(lbound(plane,1):ubound(plane,1),lbound(plane,2):ubound(plane,2)), &
            &source=cmplx(0.,0.))
        call workspace%shift_residual(rotmat,shift,zero_plane,plane,ignored_objective)
        return
    endif
    if( .not. buffers%exists ) error stop 'finite-box prediction buffers are not initialized'
    padded_box = OSMPL_PAD_FAC*box
    padded_radius = padded_box/2
    call buffers%padded%zero_and_flag_ft()
    do k = -padded_radius, padded_radius-1
        do h = 0, padded_radius
            if( h*h+k*k > padded_radius**2 ) cycle
            loc = real(matmul(real([h,k,0],dp),rotmat),sp)
            call workspace%sample_with_grad(loc,value,dvalue_dloc,switch_margin)
            arg = 2._dp*real(PI,dp)*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(padded_box,dp)
            phase = cmplx(cos(arg),sin(arg),kind=sp)
            phys = [h+1,k+1+merge(padded_box,0,k<0),1]
            call buffers%padded%set_cmat_at(phys(1),phys(2),phys(3),phase*value)
        enddo
    enddo
    ! Apply the simulator round trip, finite native box, and selected edge treatment.
    call buffers%padded%ifft()
    call buffers%padded%fft()
    call buffers%padded%ifft()
    call buffers%padded%clip(buffers%native)
    if( apply_taper )then
        edge_mean = 0.
        call buffers%native%taper_edges_particle(nint(COSMSKHALFWIDTH),edge_mean)
    endif
    call buffers%native%fft()
    plane = sampler%extract_native_plane(buffers%native)
end subroutine predict_plane

!> Fit one real nuisance amplitude and return normalized exact-pose evidence.
subroutine fit_amplitude(box,model,observed,scale,residual,objective,data_norm)
    integer, intent(in) :: box
    complex, intent(in) :: model(-box/2:,-box/2:,:), observed(-box/2:,-box/2:,:)
    real(dp), intent(out) :: scale, residual, objective, data_norm
    complex(dp) :: model_value, observed_value
    real(dp) :: cross_term, fitted_squared, model_norm
    integer :: h, i, k

    cross_term = 0._dp
    data_norm = 0._dp
    model_norm = 0._dp
    do i = 1, size(model,3)
        do k = lbound(model,2), ubound(model,2)
            do h = lbound(model,1), ubound(model,1)
                if( h*h+k*k < 4 .or. h*h+k*k > (box/2)**2 ) cycle
                model_value = cmplx(model(h,k,i),kind=dp)
                observed_value = cmplx(observed(h,k,i),kind=dp)
                model_norm = model_norm+abs(model_value)**2
                data_norm = data_norm+abs(observed_value)**2
                cross_term = cross_term+real(conjg(model_value)*observed_value,dp)
            enddo
        enddo
    enddo
    scale = cross_term/max(model_norm,tiny(model_norm))
    fitted_squared = max(0._dp,data_norm+scale*scale*model_norm-2._dp*scale*cross_term)
    residual = sqrt(fitted_squared/max(data_norm,tiny(data_norm)))
    objective = 0.5_dp*fitted_squared/max(data_norm,tiny(data_norm))
end subroutine fit_amplitude

!> Evaluate one particle objective with the arm's amplitude held fixed.
pure real(dp) function plane_objective(box,model,observed,scale,total_data_norm) result(objective)
    integer, intent(in) :: box
    complex, intent(in) :: model(-box/2:,-box/2:), observed(-box/2:,-box/2:)
    real(dp), intent(in) :: scale, total_data_norm
    complex(dp) :: residual
    integer :: h, k

    objective = 0._dp
    do k = lbound(model,2), ubound(model,2)
        do h = lbound(model,1), ubound(model,1)
            if( h*h+k*k < 4 .or. h*h+k*k > (box/2)**2 ) cycle
            residual = scale*cmplx(model(h,k),kind=dp)-cmplx(observed(h,k),kind=dp)
            objective = objective+0.5_dp*abs(residual)**2/max(total_data_norm,tiny(total_data_norm))
        enddo
    enddo
end function plane_objective

!> Construct a Cartesian tangent perturbation along one rotation axis.
pure function axis_vector(axis,magnitude) result(omega)
    integer, intent(in) :: axis
    real(dp), intent(in) :: magnitude
    real(dp) :: omega(3)

    omega = 0._dp
    omega(axis) = magnitude
end function axis_vector

!> Allocate reusable FFT plans for the finite-box prediction branch.
subroutine initialize_buffers(box,buffers)
    integer, intent(in) :: box
    type(prediction_buffers), intent(inout) :: buffers

    call buffers%padded%new([OSMPL_PAD_FAC*box,OSMPL_PAD_FAC*box,1], &
        &CONTRACT_SMPD,wthreads=.false.)
    call buffers%native%new([box,box,1],CONTRACT_SMPD,wthreads=.false.)
    buffers%exists = .true.
end subroutine initialize_buffers

!> Release optional finite-box prediction buffers.
subroutine kill_buffers(buffers)
    type(prediction_buffers), intent(inout) :: buffers

    if( .not. buffers%exists ) return
    call buffers%native%kill()
    call buffers%padded%kill()
    buffers%exists = .false.
end subroutine kill_buffers

!> Report whether one model includes the PCG inverse gather envelope.
pure logical function model_uses_deapod(model_kind) result(uses_deapod)
    integer, intent(in) :: model_kind

    uses_deapod = model_kind == CONTRACT_MODEL_DEAPOD .or. &
        &model_kind == CONTRACT_MODEL_WINDOW_DEAPOD
end function model_uses_deapod

!> Report whether one model includes finite-box and edge preprocessing.
pure logical function model_uses_window(model_kind) result(uses_window)
    integer, intent(in) :: model_kind

    uses_window = model_kind == CONTRACT_MODEL_WINDOW .or. &
        &model_kind == CONTRACT_MODEL_WINDOW_DEAPOD
end function model_uses_window

end module continuous_3D_pcg_refinement_operator_contract_support
