!@descr: library test of Cartesian pose refinement on simulated 1JYX particles (simple_pose_cont_refine3D_adapter)
! A long-running quality gate for the five-parameter LM: 1JYX from the embedded
! coordinates at box 144, 5 000 reproducible projections with varying CTF and finite
! noise, every starting pose perturbed by exactly 15 degrees and two pixels, no PFTC
! search. Refines every particle through the adapter, reconstructs the truth, perturbed
! and refined pose sets and scores them by FSC and truth-map correlation. Pinned: the
! aggregate objective, rotation error and shift error fall, and the refined
! reconstruction correlates better with the truth than the perturbed one. The run
! directory keeps 1JYX.mrc, the three reconstructions, pose_metrics.tsv and
! reconstruction_fsc.tsv as the reviewable record. Nightly (lib_cart_align3D).
module simple_pose_cont_1jyx_tester
use ieee_arithmetic, only: ieee_is_finite
use iso_fortran_env, only: int64
!$ use omp_lib, only: omp_get_max_threads, omp_get_thread_num
use simple_atoms, only: atoms
use simple_cartesian_pose_refiner, only: right_increment_rotation
use simple_cmdline, only: cmdline
use simple_commanders_sim, only: commander_simulate_particles
use simple_core_module_api
use simple_image, only: image
use simple_memoize_ft_maps, only: forget_ft_maps, memoize_ft_maps
use simple_molecule_data, only: betagal_1jyx, molecule_data
use simple_ori, only: ori
use simple_oris, only: oris
use simple_parameters, only: parameters
use simple_pose_cont_refine3D_adapter, only: cartesian_pose_data, &
    &LM_ACCEPTED_IMPROVEMENT, pose_cont_config, pose_cont_limits, pose_cont_pose, &
    &pose_cont_reference_workspace, pose_cont_transaction_result, &
    &POSE_CONT_ROUTE_SHIFT_THEN_JOINT, prepare_pose_cont_observation, &
    &remove_pose_cont_reference_artifacts, write_pose_cont_reference_artifact
use simple_reconstructor, only: reconstructor
use simple_sp_project, only: sp_project
use simple_sym, only: sym
use simple_ui, only: make_ui
use simple_test_utils
implicit none
private

#include "simple_local_flags.inc"

public :: run_all_pose_cont_1jyx_tests

integer, parameter :: TEST_BOX = 144
integer, parameter :: TEST_PARTICLES = 5000
integer, parameter :: SIMULATION_SEED = 20260914
integer, parameter :: PERTURBATION_SEED = 20260915
integer, parameter :: NOISE_SEED = 20260916
real, parameter :: TEST_SMPD = 1.3
real, parameter :: TEST_MSKDIAM = 120.
real, parameter :: TEST_MSKRAD = 46.
real, parameter :: TEST_SNR = 1.
real, parameter :: TEST_LP = 8.
real, parameter :: TEST_KV = 300.
real, parameter :: TEST_CS = 2.7
real, parameter :: TEST_FRACA = 0.1
real, parameter :: TEST_DEFOCUS = 1.5
real, parameter :: TEST_DFERR = 0.5
real, parameter :: TEST_ASTIGERR = 0.1
real(dp), parameter :: ROTATION_ERROR = 15._dp*real(PI,dp)/180._dp
real(dp), parameter :: SHIFT_ERROR = 2._dp
character(len=*), parameter :: ATOM_VOLUME_FILE = '1JYX_atoms.mrc'
character(len=*), parameter :: TRUTH_VOLUME_FILE = '1JYX.mrc'
character(len=*), parameter :: PARTICLE_FILE = '1JYX_particles.mrcs'
character(len=*), parameter :: SIMULATION_ORIENTATION_FILE = '1JYX_simulation_orientations.txt'
character(len=*), parameter :: TRUTH_ORIENTATION_FILE = '1JYX_truth_orientations.txt'

contains

subroutine run_all_pose_cont_1jyx_tests()
    write(*,'(A)') '**** running all pose_cont 1JYX tests ****'
    write(*,'(A)') 'test_pose_cont_1jyx_reconstruction'
    call run_pose_cont_1jyx_reconstruction()
end subroutine run_all_pose_cont_1jyx_tests

subroutine run_pose_cont_1jyx_reconstruction()
    type(atoms) :: molecule
    type(commander_simulate_particles) :: simulator
    type(image) :: mask_image, particle_reader, truth_image
    type(image), allocatable :: raw_work(:), fourier_work(:)
    type(molecule_data) :: molecule_record
    type(ori) :: truth_orientation
    type(oris) :: truth_orientations
    type(pose_cont_config) :: config
    type(pose_cont_limits) :: limits
    type(pose_cont_reference_workspace) :: workspace
    type(cmdline) :: simulation_command
    type(ctfparams), allocatable :: ctf_parameters(:)
    type(string) :: result_directory
    real, allocatable :: particles(:,:,:), sigma2(:)
    real(dp), allocatable :: truth_rotations(:,:,:), initial_rotations(:,:,:), terminal_rotations(:,:,:)
    real(dp), allocatable :: truth_shifts(:,:), initial_shifts(:,:), terminal_shifts(:,:)
    real(dp), allocatable :: initial_rotation_errors(:), terminal_rotation_errors(:)
    real(dp), allocatable :: initial_shift_errors(:), terminal_shift_errors(:)
    real(dp), allocatable :: objectives_before(:), objectives_after(:)
    integer, allocatable :: statuses(:)
    logical, allocatable :: noise_mask(:,:,:)
    logical :: reconstruction_control_valid, reconstruction_improved
    integer :: i, ldim(3), nsections, nworkers

    call make_ui
    call simple_getcwd(result_directory)
    write(logfhandle,'(a)') 'POSE_CONT_1JYX_ROOT: '//result_directory%to_char()
    nworkers = 1
!$  nworkers = omp_get_max_threads()

    ! Stage 1: build the reference workspace from a centered atomic truth.
    ! The simulator masks ATOM_VOLUME_FILE once;
    ! the matching masked reference is persisted as the user-facing 1JYX.mrc.
    molecule_record = betagal_1jyx()
    call molecule%pdb2mrc(volfile=string(ATOM_VOLUME_FILE),smpd=TEST_SMPD, &
        &mol=molecule_record,center_pdb=.true.,vol_dim=[TEST_BOX,TEST_BOX,TEST_BOX])
    call molecule%kill
    call truth_image%new([TEST_BOX,TEST_BOX,TEST_BOX],TEST_SMPD,wthreads=.false.)
    call truth_image%read(string(ATOM_VOLUME_FILE))
    call truth_image%mask3D_soft(TEST_MSKRAD,backgr=0.)
    call truth_image%write(string(TRUTH_VOLUME_FILE),del_if_exists=.true.)

    ! The standalone adapter needs no discrete orientation bank. Both half slots
    ! intentionally contain the same known truth reference for this C1 fixture.
    ! Remove artifacts left by an interrupted rerun before materializing this truth.
    call remove_pose_cont_reference_artifacts(1)
    call write_pose_cont_reference_artifact(truth_image,1,'even')
    call write_pose_cont_reference_artifact(truth_image,1,'odd')
    call workspace%new_from_artifacts(1,TEST_BOX,TEST_SMPD)

    ! Stage 2: freeze orientations and CTFs, then ask SIMPLE's simulator to
    ! generate the corresponding clean observations. Noise is added below from
    ! a separate deterministic stream because parameters%new reseeds SIMPLE.
    call truth_orientations%new(TEST_PARTICLES,is_ptcl=.true.)
    call set_deterministic_seed(SIMULATION_SEED)
    call truth_orientations%rnd_oris(0.)
    call truth_orientations%rnd_ctf(TEST_KV,TEST_CS,TEST_FRACA,TEST_DEFOCUS, &
        &TEST_DFERR,TEST_ASTIGERR)
    call truth_orientations%set_all2single('state',1.)
    call truth_orientations%write(string(SIMULATION_ORIENTATION_FILE),[1,TEST_PARTICLES])

    call simulation_command%set('prg','simulate_particles')
    call simulation_command%set('mkdir','no')
    call simulation_command%set('vol1',ATOM_VOLUME_FILE)
    call simulation_command%set('outstk',PARTICLE_FILE)
    call simulation_command%set('oritab',SIMULATION_ORIENTATION_FILE)
    call simulation_command%set('outfile',TRUTH_ORIENTATION_FILE)
    call simulation_command%set('nptcls',TEST_PARTICLES)
    call simulation_command%set('nthr',nworkers)
    call simulation_command%set('smpd',TEST_SMPD)
    call simulation_command%set('mskdiam',TEST_MSKDIAM)
    call simulation_command%set('pgrp','c1')
    call simulation_command%set('ctf','yes')
    ! snr>=5 suppresses simulator noise; deterministic noise is added on read.
    call simulation_command%set('snr',5.)
    call simulation_command%set('kv',TEST_KV)
    call simulation_command%set('cs',TEST_CS)
    call simulation_command%set('fraca',TEST_FRACA)
    call simulation_command%set('defocus',TEST_DEFOCUS)
    call simulation_command%set('dferr',TEST_DFERR)
    call simulation_command%set('astigerr',TEST_ASTIGERR)
    call simulation_command%set('bfac',0.)
    call simulation_command%set('bfacerr',0.)
    call simulation_command%set('sherr',0.)
    call simulator%execute(simulation_command)
    call simulation_command%kill

    call find_ldim_nptcls(string(PARTICLE_FILE),ldim,nsections)
    if( any(ldim(1:2) /= [TEST_BOX,TEST_BOX]) .or. nsections /= TEST_PARTICLES ) &
        &THROW_HARD('1JYX simulator returned an incompatible particle stack')

    ! oris%read fills existing records; it does not allocate the destination.
    call truth_orientations%kill
    call truth_orientations%new(TEST_PARTICLES,is_ptcl=.true.)
    call truth_orientations%read(string(TRUTH_ORIENTATION_FILE),[1,TEST_PARTICLES])
    if( truth_orientations%get_noris() /= TEST_PARTICLES ) &
        &THROW_HARD('1JYX simulator did not return 5,000 truth orientations')

    ! Stage 3: retain the finite particle stack and construct controlled seeds.
    allocate(particles(TEST_BOX,TEST_BOX,TEST_PARTICLES))
    call particle_reader%new([TEST_BOX,TEST_BOX,1],TEST_SMPD,wthreads=.false.)
    call set_deterministic_seed(NOISE_SEED)
    do i = 1, TEST_PARTICLES
        call particle_reader%read(string(PARTICLE_FILE),i)
        call particle_reader%add_gauran(TEST_SNR)
        ! Copy only the logical image region; image%rmat may include FFT padding.
        call particle_reader%get_rmat_sub(particles(:,:,i:i))
    enddo
    call particle_reader%kill

    allocate(ctf_parameters(TEST_PARTICLES))
    allocate(truth_rotations(3,3,TEST_PARTICLES),truth_shifts(2,TEST_PARTICLES))
    allocate(initial_rotations(3,3,TEST_PARTICLES),terminal_rotations(3,3,TEST_PARTICLES))
    allocate(initial_shifts(2,TEST_PARTICLES),terminal_shifts(2,TEST_PARTICLES))
    call set_deterministic_seed(PERTURBATION_SEED)
    do i = 1, TEST_PARTICLES
        call truth_orientations%get_ori(i,truth_orientation)
        ctf_parameters(i) = truth_orientation%get_ctfvars()
        ! A standalone orientation table has no stack-level sampling metadata.
        ! Production obtains it from the project; this fixture owns it directly.
        ctf_parameters(i)%smpd = TEST_SMPD
        if( ctf_parameters(i)%ctfflag /= CTFFLAG_YES ) &
            &THROW_HARD('1JYX simulator did not preserve enabled CTF metadata')
        truth_rotations(:,:,i) = real(truth_orientation%get_mat(),dp)
        truth_shifts(:,i) = real(truth_orientation%get_2Dshift(),dp)
        call make_perturbed_pose(truth_rotations(:,:,i),truth_shifts(:,i), &
            &initial_rotations(:,:,i),initial_shifts(:,i))
    enddo
    call truth_orientation%kill

    call mask_image%disc([TEST_BOX,TEST_BOX,1],TEST_SMPD,TEST_MSKRAD,noise_mask)
    call mask_image%kill
    allocate(sigma2(0:TEST_BOX/2),source=1.)
    allocate(initial_rotation_errors(TEST_PARTICLES),terminal_rotation_errors(TEST_PARTICLES))
    allocate(initial_shift_errors(TEST_PARTICLES),terminal_shift_errors(TEST_PARTICLES))
    allocate(objectives_before(TEST_PARTICLES),objectives_after(TEST_PARTICLES))
    allocate(statuses(TEST_PARTICLES))
    config = pose_cont_config(route=POSE_CONT_ROUTE_SHIFT_THEN_JOINT)
    limits = pose_cont_limits(shift_step_bound=1._dp,max_total_shift=5._dp)

    allocate(raw_work(nworkers),fourier_work(nworkers))
    do i = 1, nworkers
        call raw_work(i)%new([TEST_BOX,TEST_BOX,1],TEST_SMPD,wthreads=.false.)
        call fourier_work(i)%new([TEST_BOX,TEST_BOX,1],TEST_SMPD,wthreads=.false.)
        call fourier_work(i)%memoize_mask_coords()
    enddo

    ! Stage 4: refine all particles. Each OpenMP worker shares the immutable
    ! Cartesian reference and owns its
    ! image/LM temporaries. There is no PFTC search or stored-pose handoff here.
!$omp parallel default(shared)
    block
        type(cartesian_pose_data) :: data
        type(pose_cont_pose) :: seed
        type(pose_cont_transaction_result) :: result
        type(ctfparams) :: cropped_ctf
        type(ori) :: orientation_truth
        complex, allocatable :: observed(:,:)
        real(dp) :: truth_rotation(3,3), truth_shift(2)
        integer :: iparticle, thread_index

        thread_index = 1
!$      thread_index = omp_get_thread_num()+1
!$omp do schedule(dynamic,1)
        do iparticle = 1, TEST_PARTICLES
            call truth_orientations%get_ori(iparticle,orientation_truth)
            truth_rotation = real(orientation_truth%get_mat(),dp)
            truth_shift = real(orientation_truth%get_2Dshift(),dp)
            call raw_work(thread_index)%set_rmat(particles(:,:,iparticle:iparticle),.false.)
            call prepare_pose_cont_observation(raw_work(thread_index),noise_mask, &
                &fourier_work(thread_index),TEST_MSKRAD,TEST_SMPD,ctf_parameters(iparticle), &
                &observed,cropped_ctf)
            call workspace%prepare_particle(1,.true.,observed,cropped_ctf,sigma2, &
                &[2,active_upper_shell()],data)
            seed = pose_cont_pose(rotmat=initial_rotations(:,:,iparticle), &
                &shift=initial_shifts(:,iparticle))
            call workspace%refine_particle(1,.true.,seed,data,config,limits,result)
            terminal_rotations(:,:,iparticle) = result%pose%rotmat
            terminal_shifts(:,iparticle) = result%pose%shift
            objectives_before(iparticle) = result%objective_before
            objectives_after(iparticle) = result%objective_after
            statuses(iparticle) = result%status
            initial_rotation_errors(iparticle) = rotation_distance(seed%rotmat,truth_rotation)
            terminal_rotation_errors(iparticle) = rotation_distance(result%pose%rotmat,truth_rotation)
            initial_shift_errors(iparticle) = norm2(seed%shift-truth_shift)
            terminal_shift_errors(iparticle) = norm2(result%pose%shift-truth_shift)
        enddo
!$omp end do
        call orientation_truth%kill
    end block
!$omp end parallel

    do i = 1, nworkers
        call fourier_work(i)%kill
        call raw_work(i)%kill
    enddo
    deallocate(fourier_work,raw_work)

    ! Stage 5: persist particle metrics, reconstruct both pose sets, and apply
    ! aggregate quantitative acceptance checks.
    call write_pose_metrics(statuses,objectives_before,objectives_after, &
        &initial_rotation_errors,terminal_rotation_errors,initial_shift_errors,terminal_shift_errors)
    call reconstruct_and_score(particles,ctf_parameters,truth_image,noise_mask, &
        &truth_rotations,truth_shifts,initial_rotations,initial_shifts, &
        &terminal_rotations,terminal_shifts,reconstruction_control_valid,reconstruction_improved)
    call assert_pose_improvement(statuses,objectives_before,objectives_after, &
        &initial_rotation_errors,terminal_rotation_errors,initial_shift_errors,terminal_shift_errors)
    call assert_true(reconstruction_control_valid, 'the exact-pose reconstruction control is valid')
    call assert_true(reconstruction_improved, 'the refined reconstruction correlates better with the truth than the perturbed one')

    call workspace%kill
    call remove_pose_cont_reference_artifacts(1)
    call truth_orientations%kill
    call truth_image%kill
    write(logfhandle,'(a)') 'POSE_CONT_1JYX_RESULTS: '//result_directory%to_char()
end subroutine run_pose_cont_1jyx_reconstruction

subroutine make_perturbed_pose(truth_rotation,truth_shift,rotation,shift)
    real(dp), intent(in) :: truth_rotation(3,3), truth_shift(2)
    real(dp), intent(out) :: rotation(3,3), shift(2)
    real(dp) :: axis(3), azimuth, radial, uniform(3), z

    call random_number(uniform)
    z = 2._dp*uniform(1)-1._dp
    azimuth = 2._dp*real(PI,dp)*uniform(2)
    radial = sqrt(max(0._dp,1._dp-z*z))
    axis = [radial*cos(azimuth),radial*sin(azimuth),z]
    rotation = right_increment_rotation(truth_rotation,ROTATION_ERROR*axis)
    azimuth = 2._dp*real(PI,dp)*uniform(3)
    shift = truth_shift+SHIFT_ERROR*[cos(azimuth),sin(azimuth)]
end subroutine make_perturbed_pose

integer pure function active_upper_shell() result(shell)
    shell = min(TEST_BOX/2-1,int(real(TEST_BOX)*TEST_SMPD/TEST_LP))
end function active_upper_shell

pure real(dp) function rotation_distance(rotation_a,rotation_b) result(distance)
    real(dp), intent(in) :: rotation_a(3,3), rotation_b(3,3)
    real(dp) :: relative_rotation(3,3), cosine

    relative_rotation = matmul(transpose(rotation_a),rotation_b)
    cosine = 0.5_dp*(relative_rotation(1,1)+relative_rotation(2,2)+ &
        &relative_rotation(3,3)-1._dp)
    distance = acos(max(-1._dp,min(1._dp,cosine)))
end function rotation_distance

subroutine reconstruct_and_score(particles,ctf_parameters,truth_image,noise_mask, &
    &truth_rotations,truth_shifts,initial_rotations,initial_shifts, &
    &terminal_rotations,terminal_shifts,control_valid,refinement_improved)
    real, intent(in) :: particles(:,:,:)
    type(ctfparams), intent(in) :: ctf_parameters(:)
    class(image), intent(in) :: truth_image
    logical, intent(in) :: noise_mask(:,:,:)
    real(dp), intent(in) :: truth_rotations(:,:,:), truth_shifts(:,:)
    real(dp), intent(in) :: initial_rotations(:,:,:), initial_shifts(:,:)
    real(dp), intent(in) :: terminal_rotations(:,:,:), terminal_shifts(:,:)
    logical, intent(out) :: control_valid, refinement_improved
    type(fplane_type) :: plane
    type(image) :: particle, padded_particle, control_map, perturbed_map, refined_map
    type(image) :: truth_metric, control_metric, perturbed_metric, refined_metric
    type(ori) :: orientation
    type(parameters), target :: reconstruction_parameters
    type(reconstructor) :: control_reconstructor, perturbed_reconstructor, refined_reconstructor
    type(sp_project) :: empty_project
    type(sym) :: c1
    real, allocatable :: control_array(:,:,:), perturbed_array(:,:,:), refined_array(:,:,:), truth_array(:,:,:)
    real, allocatable :: control_fsc(:), perturbed_fsc(:), refined_fsc(:), resolutions(:)
    real(dp) :: control_correlation, perturbed_correlation, refined_correlation
    real :: control_fsc05, control_fsc0143, perturbed_fsc05, perturbed_fsc0143
    real :: refined_fsc05, refined_fsc0143
    integer :: i

    control_valid = .false.
    refinement_improved = .false.

    ! Step 1: create matched gridding accumulators for the exact, perturbed,
    ! and LM poses. The exact-pose arm validates this reconstruction harness.
    reconstruction_parameters%box = TEST_BOX
    reconstruction_parameters%box_crop = TEST_BOX
    reconstruction_parameters%box_croppd = OSMPL_PAD_FAC*TEST_BOX
    reconstruction_parameters%smpd_crop = TEST_SMPD
    reconstruction_parameters%nstates = 1
    reconstruction_parameters%numlen = 1
    reconstruction_parameters%oritype = 'cls3D'
    call control_reconstructor%new_accumulator(reconstruction_parameters,empty_project, &
        &expand=.true.,wthreads=.false.)
    call perturbed_reconstructor%new_accumulator(reconstruction_parameters,empty_project, &
        &expand=.true.,wthreads=.false.)
    call refined_reconstructor%new_accumulator(reconstruction_parameters,empty_project, &
        &expand=.true.,wthreads=.false.)
    call control_reconstructor%set_ft(.true.)
    call perturbed_reconstructor%set_ft(.true.)
    call refined_reconstructor%set_ft(.true.)
    call particle%new([TEST_BOX,TEST_BOX,1],TEST_SMPD,wthreads=.false.)
    call padded_particle%new([OSMPL_PAD_FAC*TEST_BOX,OSMPL_PAD_FAC*TEST_BOX,1], &
        &TEST_SMPD,wthreads=.false.)
    call control_map%new([TEST_BOX,TEST_BOX,TEST_BOX],TEST_SMPD,wthreads=.false.)
    call perturbed_map%new([TEST_BOX,TEST_BOX,TEST_BOX],TEST_SMPD,wthreads=.false.)
    call refined_map%new([TEST_BOX,TEST_BOX,TEST_BOX],TEST_SMPD,wthreads=.false.)
    call orientation%new(.false.)
    call c1%new('c1')
    call memoize_ft_maps([OSMPL_PAD_FAC*TEST_BOX,OSMPL_PAD_FAC*TEST_BOX,1],TEST_SMPD)

    ! Step 2: insert each noisy observation into all three accumulators; only
    ! the supplied pose differs.
    do i = 1, size(particles,3)
        call particle%set_rmat(particles(:,:,i:i),.false.)
        call particle%norm_noise_taper_edge_pad_fft(noise_mask,padded_particle)

        call orientation%set_euler(real(dm2euler(truth_rotations(:,:,i))))
        call padded_particle%gen_fplane4rec([2,active_upper_shell()],TEST_SMPD, &
            &ctf_parameters(i),real(truth_shifts(:,i)),plane)
        call control_reconstructor%insert_plane_oversamp(c1,orientation,plane)

        call orientation%set_euler(real(dm2euler(initial_rotations(:,:,i))))
        call padded_particle%gen_fplane4rec([2,active_upper_shell()],TEST_SMPD, &
            &ctf_parameters(i),real(initial_shifts(:,i)),plane)
        call perturbed_reconstructor%insert_plane_oversamp(c1,orientation,plane)

        call orientation%set_euler(real(dm2euler(terminal_rotations(:,:,i))))
        call padded_particle%gen_fplane4rec([2,active_upper_shell()],TEST_SMPD, &
            &ctf_parameters(i),real(terminal_shifts(:,i)),plane)
        call refined_reconstructor%insert_plane_oversamp(c1,orientation,plane)
    enddo

    ! Step 3: finalize and retain all maps for qualitative inspection.
    call control_reconstructor%compress_exp()
    call perturbed_reconstructor%compress_exp()
    call refined_reconstructor%compress_exp()
    call control_reconstructor%restore_final(control_map)
    call perturbed_reconstructor%restore_final(perturbed_map)
    call refined_reconstructor%restore_final(refined_map)
    call control_map%write(string('reconstruction_truth_pose.mrc'),del_if_exists=.true.)
    call perturbed_map%write(string('reconstruction_perturbed.mrc'),del_if_exists=.true.)
    call refined_map%write(string('reconstruction_pose_cont.mrc'),del_if_exists=.true.)

    ! Step 4: compare every map with the same masked truth by correlation and FSC.
    call truth_metric%copy(truth_image)
    call control_metric%copy(control_map)
    call perturbed_metric%copy(perturbed_map)
    call refined_metric%copy(refined_map)
    call control_metric%mask3D_soft(TEST_MSKRAD,backgr=0.)
    call perturbed_metric%mask3D_soft(TEST_MSKRAD,backgr=0.)
    call refined_metric%mask3D_soft(TEST_MSKRAD,backgr=0.)
    truth_array = truth_metric%get_rmat()
    control_array = control_metric%get_rmat()
    perturbed_array = perturbed_metric%get_rmat()
    refined_array = refined_metric%get_rmat()
    control_correlation = centered_array_correlation(control_array,truth_array)
    perturbed_correlation = centered_array_correlation(perturbed_array,truth_array)
    refined_correlation = centered_array_correlation(refined_array,truth_array)

    call truth_metric%fft()
    call control_metric%fft()
    call perturbed_metric%fft()
    call refined_metric%fft()
    allocate(control_fsc(truth_metric%get_filtsz()),source=0.)
    allocate(perturbed_fsc(truth_metric%get_filtsz()),source=0.)
    allocate(refined_fsc(truth_metric%get_filtsz()),source=0.)
    call truth_metric%fsc(control_metric,control_fsc)
    call truth_metric%fsc(perturbed_metric,perturbed_fsc)
    call truth_metric%fsc(refined_metric,refined_fsc)
    allocate(resolutions(size(perturbed_fsc)))
    do i = 1, size(resolutions)
        resolutions(i) = real(TEST_BOX)*TEST_SMPD/real(i)
    enddo
    call get_resolution(control_fsc,resolutions,control_fsc05,control_fsc0143)
    call get_resolution(perturbed_fsc,resolutions,perturbed_fsc05,perturbed_fsc0143)
    call get_resolution(refined_fsc,resolutions,refined_fsc05,refined_fsc0143)
    call write_reconstruction_metrics(resolutions,control_fsc,perturbed_fsc,refined_fsc, &
        &control_correlation,perturbed_correlation,refined_correlation, &
        &control_fsc05,control_fsc0143,perturbed_fsc05,perturbed_fsc0143, &
        &refined_fsc05,refined_fsc0143)

    ! Step 5: return scientific gates only after every diagnostic is written.
    control_valid = ieee_is_finite(control_correlation) .and. control_correlation > 0._dp .and. &
        &all(ieee_is_finite(control_fsc)) .and. control_fsc0143 > 0.
    refinement_improved = ieee_is_finite(perturbed_correlation) .and. &
        &ieee_is_finite(refined_correlation) .and. all(ieee_is_finite(perturbed_fsc)) .and. &
        &all(ieee_is_finite(refined_fsc)) .and. refined_correlation > perturbed_correlation

    if( allocated(plane%cmplx_plane) ) deallocate(plane%cmplx_plane)
    if( allocated(plane%ctfsq_plane) ) deallocate(plane%ctfsq_plane)
    if( allocated(plane%transfer_plane) ) deallocate(plane%transfer_plane)
    call forget_ft_maps()
    call c1%kill
    call orientation%kill
    call control_metric%kill
    call refined_metric%kill
    call perturbed_metric%kill
    call truth_metric%kill
    call control_map%kill
    call refined_map%kill
    call perturbed_map%kill
    call padded_particle%kill
    call particle%kill
    call control_reconstructor%kill
    call refined_reconstructor%kill
    call perturbed_reconstructor%kill
    call empty_project%kill
end subroutine reconstruct_and_score

subroutine write_pose_metrics(statuses,objectives_before,objectives_after, &
    &rotation_before,rotation_after,shift_before,shift_after)
    integer, intent(in) :: statuses(:)
    real(dp), intent(in) :: objectives_before(:), objectives_after(:)
    real(dp), intent(in) :: rotation_before(:), rotation_after(:), shift_before(:), shift_after(:)
    integer :: file_unit, i

    open(newunit=file_unit,file='pose_metrics.tsv',status='replace',action='write')
    write(file_unit,'(a)') 'particle'//char(9)//'status'//char(9)// &
        &'objective_before'//char(9)//'objective_after'//char(9)// &
        &'rotation_before_deg'//char(9)//'rotation_after_deg'//char(9)// &
        &'shift_before_px'//char(9)//'shift_after_px'
    do i = 1, size(statuses)
        write(file_unit,'(i0,a,i0,6(a,es16.8))') i,char(9),statuses(i), &
            &char(9),objectives_before(i),char(9),objectives_after(i), &
            &char(9),rotation_before(i)*180._dp/real(PI,dp), &
            &char(9),rotation_after(i)*180._dp/real(PI,dp), &
            &char(9),shift_before(i),char(9),shift_after(i)
    enddo
    close(file_unit)
end subroutine write_pose_metrics

subroutine write_reconstruction_metrics(resolutions,control_fsc,perturbed_fsc,refined_fsc, &
    &control_correlation,perturbed_correlation,refined_correlation,control_fsc05, &
    &control_fsc0143,perturbed_fsc05,perturbed_fsc0143,refined_fsc05,refined_fsc0143)
    real, intent(in) :: resolutions(:), control_fsc(:), perturbed_fsc(:), refined_fsc(:)
    real(dp), intent(in) :: control_correlation, perturbed_correlation, refined_correlation
    real, intent(in) :: control_fsc05, control_fsc0143, perturbed_fsc05, perturbed_fsc0143
    real, intent(in) :: refined_fsc05, refined_fsc0143
    integer :: file_unit, i

    open(newunit=file_unit,file='reconstruction_fsc.tsv',status='replace',action='write')
    write(file_unit,'(a)') 'shell'//char(9)//'resolution_A'//char(9)// &
        &'truth_vs_exact_pose'//char(9)//'truth_vs_perturbed'//char(9)//'truth_vs_pose_cont'
    do i = 1, size(resolutions)
        write(file_unit,'(i0,4(a,es16.8))') i,char(9),resolutions(i), &
            &char(9),control_fsc(i),char(9),perturbed_fsc(i),char(9),refined_fsc(i)
    enddo
    close(file_unit)
    write(logfhandle,'(a,3(1x,es12.4))') 'POSE_CONT_1JYX map correlation exact/perturbed/refined:', &
        &control_correlation,perturbed_correlation,refined_correlation
    write(logfhandle,'(a,3(1x,f8.3))') 'POSE_CONT_1JYX FSC=0.5 resolution exact/perturbed/refined (A):', &
        &control_fsc05,perturbed_fsc05,refined_fsc05
    write(logfhandle,'(a,3(1x,f8.3))') 'POSE_CONT_1JYX FSC=0.143 resolution exact/perturbed/refined (A):', &
        &control_fsc0143,perturbed_fsc0143,refined_fsc0143
end subroutine write_reconstruction_metrics

subroutine assert_pose_improvement(statuses,objectives_before,objectives_after, &
    &rotation_before,rotation_after,shift_before,shift_after)
    integer, intent(in) :: statuses(:)
    real(dp), intent(in) :: objectives_before(:), objectives_after(:)
    real(dp), intent(in) :: rotation_before(:), rotation_after(:), shift_before(:), shift_after(:)
    logical, allocatable :: finite_objective(:)
    real(dp) :: objective_before_mean, objective_after_mean
    real(dp) :: rotation_before_rms, rotation_after_rms, shift_before_rms, shift_after_rms
    integer :: accepted

    finite_objective = ieee_is_finite(objectives_before) .and. ieee_is_finite(objectives_after) .and. &
        &objectives_before >= 0._dp .and. objectives_after >= 0._dp
    call assert_int(size(statuses), count(finite_objective), 'every particle has finite, non-negative objectives before and after')
    accepted = count(statuses == LM_ACCEPTED_IMPROVEMENT)
    objective_before_mean = sum(objectives_before)/real(size(statuses),dp)
    objective_after_mean = sum(objectives_after)/real(size(statuses),dp)
    rotation_before_rms = sqrt(sum(rotation_before**2)/real(size(statuses),dp))
    rotation_after_rms = sqrt(sum(rotation_after**2)/real(size(statuses),dp))
    shift_before_rms = sqrt(sum(shift_before**2)/real(size(statuses),dp))
    shift_after_rms = sqrt(sum(shift_after**2)/real(size(statuses),dp))

    write(logfhandle,'(a,i0,a,i0)') 'POSE_CONT_1JYX accepted particles: ',accepted,' / ',size(statuses)
    write(logfhandle,'(a,2(1x,es12.4))') 'POSE_CONT_1JYX mean objective before/after:', &
        &objective_before_mean,objective_after_mean
    write(logfhandle,'(a,2(1x,f8.3))') 'POSE_CONT_1JYX rotation RMS before/after (deg):', &
        &rotation_before_rms*180._dp/real(PI,dp),rotation_after_rms*180._dp/real(PI,dp)
    write(logfhandle,'(a,2(1x,f8.3))') 'POSE_CONT_1JYX shift RMS before/after (pixels):', &
        &shift_before_rms,shift_after_rms

    call assert_true(accepted >= 1, 'at least one particle was accepted')
    call assert_true(objective_after_mean < objective_before_mean, 'the aggregate Cartesian objective falls')
    call assert_true(rotation_after_rms < rotation_before_rms, 'the aggregate rotation error falls')
    call assert_true(shift_after_rms < shift_before_rms, 'the aggregate shift error falls')
end subroutine assert_pose_improvement

pure real(dp) function centered_array_correlation(array_a,array_b) result(correlation)
    real, intent(in) :: array_a(:,:,:), array_b(:,:,:)
    real(dp) :: mean_a, mean_b, denominator

    mean_a = sum(real(array_a,dp))/real(size(array_a),dp)
    mean_b = sum(real(array_b,dp))/real(size(array_b),dp)
    denominator = sqrt(sum((real(array_a,dp)-mean_a)**2)* &
        &sum((real(array_b,dp)-mean_b)**2))
    correlation = sum((real(array_a,dp)-mean_a)*(real(array_b,dp)-mean_b))/ &
        &max(denominator,epsilon(denominator))
end function centered_array_correlation

subroutine set_deterministic_seed(base_seed)
    integer, intent(in) :: base_seed
    integer, allocatable :: seed(:)
    integer(int64) :: candidate, modulus
    integer :: i, seed_size

    call random_seed(size=seed_size)
    allocate(seed(seed_size))
    modulus = int(huge(0),int64)-1_int64
    do i = 1, seed_size
        candidate = int(base_seed,int64)+104729_int64*int(i-1,int64)
        seed(i) = int(modulo(candidate,modulus))+1
    enddo
    call random_seed(put=seed)
    deallocate(seed)
end subroutine set_deterministic_seed

end module simple_pose_cont_1jyx_tester
