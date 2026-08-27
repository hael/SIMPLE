! Pose-diagnostic reference construction, including explicit PCG-reference controls.
! The production pose owner does not depend on this test-only module.
module pose_cont_refinement_reference_support
use pose_cont_refinement_test_helpers, only: build_truth_volume, set_deterministic_seed, &
    &BOX => TRUTH_VOLUME_BOX
use simple_defs, only: dp, OSMPL_PAD_FAC
use simple_image, only: image
use simple_ori, only: ori
use simple_oris, only: oris
use simple_projector, only: projector
use simple_reconstructor_pcg, only: reconstructor_pcg
use simple_string, only: string
implicit none
private

real, parameter :: HALFSET_SMPD = 1.5
real, parameter :: HALFSET_LAMBDA = 1.e-3
real, parameter :: HALFSET_MASK_RADIUS = 10.

public :: build_disjoint_half_orientations
public :: build_independent_observations
public :: apply_fixed_support
public :: array_l2_norm
public :: calculate_fsc
public :: centered_array_correlation
public :: HALFSET_LAMBDA
public :: HALFSET_MASK_RADIUS
public :: HALFSET_SMPD
public :: reconstruct_half
public :: reconstruct_half_fixed
public :: reconstruct_half_pair_fixed
public :: reconstruct_half_trajectory
public :: relative_forward_residual
public :: write_mrc_volume

contains

!> Build deterministic, nonoverlapping even and odd orientation sets.
subroutine build_disjoint_half_orientations(nhalf, even_oris, odd_oris, even_ids, odd_ids)
    integer,          intent(in)    :: nhalf
    type(oris),       intent(inout) :: even_oris, odd_oris
    integer, allocatable, intent(out) :: even_ids(:), odd_ids(:)
    type(ori) :: orientation
    type(oris) :: all_oris
    integer :: i

    call all_oris%new(2*nhalf, .false.)
    call all_oris%spiral()
    call even_oris%new(nhalf, .false.)
    call odd_oris%new(nhalf, .false.)
    call orientation%new(.false.)
    allocate(even_ids(nhalf), odd_ids(nhalf))
    do i = 1, nhalf
        ! Fixed parity gives disjoint ownership of the common spiral grid.
        odd_ids(i) = 2*i - 1
        even_ids(i) = 2*i
        call all_oris%get_ori(odd_ids(i), orientation)
        call odd_oris%set_ori(i, orientation)
        call all_oris%get_ori(even_ids(i), orientation)
        call even_oris%set_ori(i, orientation)
    enddo

    call orientation%kill()
    call all_oris%kill()
end subroutine build_disjoint_half_orientations

!> Generate simulator-path clean projections and independent seeded noise.
subroutine build_independent_observations(sampler, orientations, noise_seed, requested_snr, &
                                          &planes, clean_planes, clean_images, noisy_images, noise, realized_snr)
    class(reconstructor_pcg), intent(in)    :: sampler
    type(oris),                intent(inout) :: orientations
    integer,                   intent(in)    :: noise_seed
    real,                      intent(in)    :: requested_snr
    complex, allocatable,      intent(out)   :: planes(:,:,:)
    complex, allocatable,      intent(out)   :: clean_planes(:,:,:)
    real, allocatable,         intent(out)   :: clean_images(:,:,:)
    real, allocatable,         intent(out)   :: noisy_images(:,:,:)
    real, allocatable,         intent(out)   :: noise(:,:,:)
    real(dp),                  intent(out)   :: realized_snr
    type(image) :: source_volume, padded_projection, projection
    type(ori) :: orientation
    type(projector) :: truth_projector
    real, allocatable :: clean(:,:,:), noisy(:,:,:), phantom(:,:,:)
    integer :: i, lims2(2,2), nprojs
    real(dp) :: clean_mean, noise_mean, noise_power, signal_power

    nprojs = orientations%get_noris()
    lims2 = sampler%get_lims2()
    allocate(planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojs))
    allocate(clean_planes(lims2(1,1):lims2(1,2),lims2(2,1):lims2(2,2),nprojs))
    allocate(clean_images(BOX,BOX,nprojs), noisy_images(BOX,BOX,nprojs))
    allocate(noise(BOX,BOX,nprojs), source=0.)
    allocate(clean(BOX,BOX,1), noisy(BOX,BOX,1))
    call build_truth_volume(phantom)
    call source_volume%new([BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%set_rmat(phantom,.false.)
    call truth_projector%new(OSMPL_PAD_FAC*[BOX,BOX,BOX],HALFSET_SMPD,wthreads=.false.)
    call source_volume%pad(truth_projector)
    call truth_projector%fft()
    call truth_projector%expand_cmat()
    call padded_projection%new([OSMPL_PAD_FAC*BOX,OSMPL_PAD_FAC*BOX,1], &
        &HALFSET_SMPD,wthreads=.false.)
    call projection%new([BOX,BOX,1], HALFSET_SMPD, wthreads=.false.)
    call orientation%new(.false.)
    call set_deterministic_seed(noise_seed)
    signal_power = 0._dp
    noise_power = 0._dp
    do i = 1, nprojs
        call orientations%get_ori(i, orientation)
        ! Match simulate_particles: padded projection followed by a real-space clip.
        call truth_projector%fproject_serial(orientation,padded_projection)
        ! simimg performs this no-CTF Fourier round trip before detector noise.
        call padded_projection%ifft()
        call padded_projection%fft()
        call padded_projection%ifft()
        call padded_projection%clip(projection)
        call projection%get_rmat_sub(clean)
        clean_images(:,:,i) = clean(:,:,1)
        call projection%fft()
        clean_planes(:,:,i) = sampler%extract_native_plane(projection)
        ! Calibrate white noise on the tested native observations so their SNR is explicit.
        call projection%ifft()
        call projection%add_gauran(requested_snr)
        call projection%get_rmat_sub(noisy)
        noisy_images(:,:,i) = noisy(:,:,1)
        noise(:,:,i) = noisy(:,:,1) - clean(:,:,1)
        clean_mean = sum(real(clean(:,:,1),dp)) / real(BOX*BOX,dp)
        noise_mean = sum(real(noise(:,:,i),dp)) / real(BOX*BOX,dp)
        signal_power = signal_power + sum((real(clean(:,:,1),dp)-clean_mean)**2)
        noise_power = noise_power + sum((real(noise(:,:,i),dp)-noise_mean)**2)
        call projection%fft()
        planes(:,:,i) = sampler%extract_native_plane(projection)
    enddo
    ! SNR = total centred signal power / total centred noise power.
    realized_snr = signal_power / noise_power

    call orientation%kill()
    call projection%kill()
    call padded_projection%kill()
    call truth_projector%kill_expanded()
    call truth_projector%kill()
    call source_volume%kill()
end subroutine build_independent_observations

!> Reconstruct one fixed half-set with the matrix-free PCG operator.
subroutine reconstruct_half(orientations, planes, reconstruction, niters)
    type(oris),           intent(inout) :: orientations
    complex,              intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
    real, allocatable,    intent(out)   :: reconstruction(:,:,:)
    integer,              intent(out)   :: niters
    type(reconstructor_pcg) :: pcgop

    call pcgop%new(BOX, HALFSET_SMPD, HALFSET_LAMBDA)
    ! x = P u: both halves use the same predetermined soft spherical support.
    call pcgop%set_mask(HALFSET_MASK_RADIUS)
    call pcgop%prep_particles(orientations, use_ctf=.false.)
    call pcgop%build_operators(.false.)
    allocate(reconstruction(BOX,BOX,BOX), source=0.)
    call pcgop%solve(planes, reconstruction, maxits=40, rtol=1.e-3, niters=niters)
    call pcgop%kill()
end subroutine reconstruct_half

!> Reconstruct one half with explicit regularization, iteration, and support controls.
subroutine reconstruct_half_fixed(orientations,planes,lambda,maxits,mask_radius,reconstruction, &
    &niters,data_residual)
    type(oris), intent(inout) :: orientations
    complex, intent(in) :: planes(-BOX/2:,-BOX/2:,:)
    real, intent(in) :: lambda, mask_radius
    integer, intent(in) :: maxits
    real, allocatable, intent(out) :: reconstruction(:,:,:)
    integer, intent(out) :: niters
    real(dp), intent(out) :: data_residual
    type(reconstructor_pcg) :: pcgop

    call pcgop%new(BOX,HALFSET_SMPD,lambda)
    if( mask_radius > 0. ) call pcgop%set_mask(mask_radius)
    call pcgop%prep_particles(orientations,use_ctf=.false.)
    call pcgop%build_operators(.false.)
    allocate(reconstruction(BOX,BOX,BOX),source=0.)
    call pcgop%solve(planes,reconstruction,maxits=maxits,rtol=0.,niters=niters)
    call relative_forward_residual(pcgop,orientations,planes,reconstruction,data_residual)
    call pcgop%kill()
end subroutine reconstruct_half_fixed

!> Reconstruct matched clean and noisy half data with one fixed PCG policy.
subroutine reconstruct_half_pair_fixed(orientations, clean_planes, noisy_planes, lambda, maxits, &
                                       &clean_reconstruction, noisy_reconstruction, niters, data_residuals)
    type(oris),        intent(inout) :: orientations
    complex,           intent(in)    :: clean_planes(-BOX/2:,-BOX/2:,:), noisy_planes(-BOX/2:,-BOX/2:,:)
    real,              intent(in)    :: lambda
    integer,           intent(in)    :: maxits
    real, allocatable, intent(out)   :: clean_reconstruction(:,:,:), noisy_reconstruction(:,:,:)
    integer,           intent(out)   :: niters(2)
    real(dp),          intent(out)   :: data_residuals(2)
    type(reconstructor_pcg) :: pcgop

    call pcgop%new(BOX, HALFSET_SMPD, lambda)
    call pcgop%set_mask(HALFSET_MASK_RADIUS)
    call pcgop%prep_particles(orientations, use_ctf=.false.)
    call pcgop%build_operators(.false.)
    allocate(clean_reconstruction(BOX,BOX,BOX), source=0.)
    allocate(noisy_reconstruction(BOX,BOX,BOX), source=0.)
    ! Fixed-count paired solves isolate lambda under an identical iteration budget.
    call pcgop%solve(clean_planes, clean_reconstruction, maxits=maxits, rtol=0., niters=niters(1))
    call pcgop%solve(noisy_planes, noisy_reconstruction, maxits=maxits, rtol=0., niters=niters(2))
    call relative_forward_residual(pcgop, orientations, clean_planes, clean_reconstruction, data_residuals(1))
    call relative_forward_residual(pcgop, orientations, noisy_planes, noisy_reconstruction, data_residuals(2))
    call pcgop%kill()
end subroutine reconstruct_half_pair_fixed

!> Measure ||A x-y||/||y|| with the same matrix-free forward operator.
subroutine relative_forward_residual(pcgop, orientations, planes, volume, relative_residual)
    class(reconstructor_pcg), intent(inout) :: pcgop
    type(oris),               intent(inout) :: orientations
    complex,                  intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
    real,                     intent(in)    :: volume(BOX,BOX,BOX)
    real(dp),                 intent(out)   :: relative_residual
    type(ori) :: orientation
    complex, allocatable :: predicted(:,:)
    real(dp) :: denominator, numerator
    integer :: h, i, k

    allocate(predicted(-BOX/2:BOX/2,-BOX/2:BOX/2))
    call orientation%new(.false.)
    call pcgop%set_volume(volume)
    numerator = 0._dp
    denominator = 0._dp
    do i = 1, orientations%get_noris()
        call orientations%get_ori(i, orientation)
        call pcgop%forward_plane(orientation, predicted)
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                if( h*h+k*k > (BOX/2)**2 ) cycle
                numerator = numerator + abs(cmplx(predicted(h,k)-planes(h,k,i),kind=dp))**2
                denominator = denominator + abs(cmplx(planes(h,k,i),kind=dp))**2
            enddo
        enddo
    enddo
    relative_residual = sqrt(numerator/max(denominator,tiny(denominator)))
    call orientation%kill()
end subroutine relative_forward_residual

!> Save PCG solutions at predefined iteration counts for trajectory analysis.
subroutine reconstruct_half_trajectory(orientations, planes, iteration_counts, mask_radius, reconstructions, niters)
    type(oris),        intent(inout) :: orientations
    complex,           intent(in)    :: planes(-BOX/2:,-BOX/2:,:)
    integer,           intent(in)    :: iteration_counts(:)
    real,              intent(in)    :: mask_radius
    real, allocatable, intent(out)   :: reconstructions(:,:,:,:)
    integer, allocatable, intent(out) :: niters(:)
    type(reconstructor_pcg) :: pcgop
    integer :: i

    call pcgop%new(BOX, HALFSET_SMPD, HALFSET_LAMBDA)
    if( mask_radius > 0. ) call pcgop%set_mask(mask_radius)
    call pcgop%prep_particles(orientations, use_ctf=.false.)
    call pcgop%build_operators(.false.)
    allocate(reconstructions(BOX,BOX,BOX,size(iteration_counts)), source=0.)
    allocate(niters(size(iteration_counts)), source=0)
    do i = 1, size(iteration_counts)
        ! rtol=0 forces each solve to restart at zero and run exactly the requested count.
        call pcgop%solve(planes, reconstructions(:,:,:,i), maxits=iteration_counts(i), &
            &rtol=0., niters=niters(i))
    enddo
    call pcgop%kill()
end subroutine reconstruct_half_trajectory

!> Apply the common spherical support used by all reconstruction comparisons.
subroutine apply_fixed_support(volume)
    real, intent(inout) :: volume(:,:,:)
    type(image) :: volume_image

    call volume_image%new([BOX,BOX,BOX], HALFSET_SMPD, wthreads=.false.)
    call volume_image%set_rmat(volume, .false.)
    ! P x_truth uses the exact soft support applied inside the PCG operator.
    call volume_image%mask3D_soft(HALFSET_MASK_RADIUS, backgr=0.)
    volume = volume_image%get_rmat()
    call volume_image%kill()
end subroutine apply_fixed_support

!> Calculate the shell FSC between two independently reconstructed volumes.
subroutine calculate_fsc(even_volume, odd_volume, fsc)
    real,              intent(in)  :: even_volume(:,:,:), odd_volume(:,:,:)
    real, allocatable, intent(out) :: fsc(:)
    type(image) :: even_image, odd_image

    call even_image%new([BOX,BOX,BOX], HALFSET_SMPD, wthreads=.false.)
    call odd_image%new([BOX,BOX,BOX], HALFSET_SMPD, wthreads=.false.)
    call even_image%set_rmat(even_volume, .false.)
    call odd_image%set_rmat(odd_volume, .false.)
    call even_image%fft()
    call odd_image%fft()
    allocate(fsc(even_image%get_filtsz()))
    ! FSC(s) = Re sum_s F_even conjg(F_odd) / sqrt(sum_s|F_even|^2 sum_s|F_odd|^2).
    call even_image%fsc(odd_image, fsc)
    call even_image%kill()
    call odd_image%kill()
end subroutine calculate_fsc

!> Return the centered voxel correlation between two real volumes.
pure real(dp) function centered_array_correlation(a, b) result(correlation)
    real, intent(in) :: a(:,:,:), b(:,:,:)
    real(dp) :: amean, bmean, denominator

    amean = sum(real(a,dp)) / real(size(a),dp)
    bmean = sum(real(b,dp)) / real(size(b),dp)
    denominator = sqrt(sum((real(a,dp)-amean)**2) * sum((real(b,dp)-bmean)**2))
    correlation = sum((real(a,dp)-amean) * (real(b,dp)-bmean)) / &
        &max(denominator, tiny(denominator))
end function centered_array_correlation

!> Return the Euclidean norm of a real volume array.
pure real(dp) function array_l2_norm(volume) result(norm)
    real, intent(in) :: volume(:,:,:)

    ! ||x||_2 = sqrt(sum_j x_j^2).
    norm = sqrt(sum(real(volume,dp)**2))
end function array_l2_norm

!> Write a real test array as a reviewable MRC volume.
subroutine write_mrc_volume(volume, filename)
    real, intent(in) :: volume(:,:,:)
    character(len=*), intent(in) :: filename
    type(image) :: volume_image

    call volume_image%new([BOX,BOX,BOX], HALFSET_SMPD, wthreads=.false.)
    call volume_image%set_rmat(volume, .false.)
    call volume_image%write(string(filename), del_if_exists=.true.)
    call volume_image%kill()
    write(*,'(a,a)') 'CONTINUOUS_3D_MATRIX wrote volume: ', trim(filename)
end subroutine write_mrc_volume

end module pose_cont_refinement_reference_support
