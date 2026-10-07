!@descr: library test of Cartesian pose refinement on simulated 1JYX particles (simple_strategy3D_cont)
! Nightly gate (lib_cart_align3D): 1 000 simulated 1JYX particles at box 144 with varying CTF and
! noise, every pose perturbed by 15 deg and 2 px, refined through strategy3D_cont as the matcher
! runs it (shift_then_joint, trs 5 px, athres_cont 15 deg) under objfun=cc and euclid at lp 8 A,
! then at 4 A from the 8 A poses. Pinned: objective, rotation and shift errors fall, the refined
! map beats the perturbed one, and the floors declared below. Fixture files are removed.
module simple_strategy3D_cont_1jyx_tester
use ieee_arithmetic, only: ieee_is_finite
!$  use omp_lib, only: omp_get_max_threads, omp_get_thread_num
use simple_atoms, only: atoms
use simple_cartft_calc, only: cartft_calc
use simple_builder, only: builder
use simple_cartft_pose_opt, only: right_increment_rotation, CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT
use simple_cmdline, only: cmdline
use simple_commanders_sim, only: commander_simulate_particles
use simple_core_module_api
use simple_image, only: image
use simple_memoize_ft_maps, only: forget_ft_maps, memoize_ft_maps
use simple_molecule_data, only: betagal_1jyx, molecule_data
use simple_ori, only: ori
use simple_oris, only: oris
use simple_parameters, only: parameters
use simple_matcher_2Dprep, only: prepimg4align_cart
use simple_strategy3D, only: strategy3D
use simple_strategy3D_cont, only: strategy3D_cont
use simple_strategy3D_srch, only: strategy3D_spec
use simple_reconstructor, only: reconstructor
use simple_sp_project, only: sp_project
use simple_sym, only: sym
use simple_syslib, only: del_file
use simple_ui, only: make_ui
use simple_test_utils
implicit none
private

#include "simple_local_flags.inc"

public :: run_all_strategy3D_cont_1jyx_tests

integer, parameter :: TEST_BOX = 144
integer, parameter :: TEST_PARTICLES = 1000
integer, parameter :: SIMULATION_SEED = 20260914
integer, parameter :: PERTURBATION_SEED = 20260915
integer, parameter :: NOISE_SEED = 20260916
real, parameter :: TEST_SMPD = 1.3
real, parameter :: TEST_MSKDIAM = 120.
real, parameter :: TEST_MSKRAD = 46.
real, parameter :: TEST_SNR = 1.
real, parameter :: TEST_LPS(2) = [8., 4.]             !< the E25 band and its high-resolution case
real, parameter :: TEST_KV = 300.
real, parameter :: TEST_CS = 2.7
real, parameter :: TEST_FRACA = 0.1
real, parameter :: TEST_DEFOCUS = 1.5
real, parameter :: TEST_DFERR = 0.5
real, parameter :: TEST_ASTIGERR = 0.1
real(dp), parameter :: ROTATION_ERROR = 15._dp*real(PI, dp)/180._dp
real(dp), parameter :: SHIFT_ERROR = 2._dp
character(len=*), parameter :: ATOM_VOLUME_FILE = '1JYX_atoms.mrc'
character(len=*), parameter :: TRUTH_VOLUME_FILE = '1JYX.mrc'
character(len=*), parameter :: PARTICLE_FILE = '1JYX_particles.mrcs'
character(len=*), parameter :: SIMULATION_ORIENTATION_FILE = '1JYX_simulation_orientations.txt'
character(len=*), parameter :: TRUTH_ORIENTATION_FILE = '1JYX_truth_orientations.txt'
! bounds of the Phase 0 baseline: trs (pixels, box_crop == box) and athres_cont (degrees)
real, parameter :: TEST_TRS = 5., TEST_ATHRES_CONT = 15.
integer, parameter :: NOBJFUNS = 2
integer, parameter :: TEST_OBJFUNS(NOBJFUNS) = [OBJFUN_CC, OBJFUN_EUCLID]
character(len=6), parameter :: OBJFUN_NAMES(NOBJFUNS) = ['cc    ', 'euclid']
! floors (section 10.1), declared before the first run of the Phase 3 code
real(dp), parameter :: TRUTH_FACTOR = 0.1_dp                       !< median error / injected error
real(dp), parameter :: BASIN_FRACTION_FLOOR(2) = [0.95_dp, 0.90_dp] !< per band
real(dp), parameter :: MAP_CORR_MARGIN = 0.01_dp                    !< below the exact-pose map
real(dp), parameter :: PHASE0_MEDIAN_ROTATION_DEG = 0.194_dp        !< cc at 8 A, Phase 0 at 1 000 particles (R2)
real(dp), parameter :: PHASE0_BASIN_FRACTION = 0.9870_dp
real(dp), parameter :: PHASE0_MAP_CORR = 0.96888_dp

contains

    subroutine run_all_strategy3D_cont_1jyx_tests()
        write (*, '(A)') '**** running all pose_cont 1JYX tests ****'
        write (*, '(A)') 'test_pose_cont_1jyx_reconstruction'
        call run_pose_cont_1jyx_reconstruction()
    end subroutine run_all_strategy3D_cont_1jyx_tests

    subroutine run_pose_cont_1jyx_reconstruction()
        type(atoms) :: molecule
        type(commander_simulate_particles) :: simulator
        type(image) :: mask_image, particle_reader, truth_image
        type(molecule_data) :: molecule_record
        type(ori) :: truth_orientation
        type(oris) :: truth_orientations
        type(builder), target :: b
        type(parameters), target :: p
        type(cartft_calc), pointer :: calc
        type(cmdline) :: simulation_command
        type(ctfparams), allocatable :: ctf_parameters(:)
        real, allocatable :: particles(:, :, :), sigma2(:), truth_volume(:, :, :), taper(:)
        real(dp), allocatable :: truth_rotations(:, :, :), initial_rotations(:, :, :), terminal_rotations(:, :, :, :)
        real(dp), allocatable :: truth_shifts(:, :), initial_shifts(:, :), terminal_shifts(:, :, :)
        real(dp), allocatable :: start_rotations(:, :, :), start_shifts(:, :)
        real(dp), allocatable :: initial_rotation_errors(:), terminal_rotation_errors(:)
        real(dp), allocatable :: initial_shift_errors(:), terminal_shift_errors(:)
        real(dp), allocatable :: objectives_before(:), objectives_after(:), objectives_truth(:)
        real(dp) :: median_rotation(NOBJFUNS, size(TEST_LPS)), gain, basin_width
        integer, allocatable :: statuses(:)
        logical, allocatable :: noise_mask(:, :, :)
        integer :: i, iband, iobj, ldim(3), nsections, nworkers, upper_shell

        call make_ui
        call cleanup_1jyx_fixture()
        nworkers = 1
!$      nworkers = omp_get_max_threads()

        ! Stage 1: build the Cartesian references from a centered atomic truth.
        ! The simulator masks ATOM_VOLUME_FILE once;
        ! the matching masked reference is persisted as the user-facing 1JYX.mrc.
        molecule_record = betagal_1jyx()
        call molecule%pdb2mrc(volfile=string(ATOM_VOLUME_FILE), smpd=TEST_SMPD, &
            &mol=molecule_record, center_pdb=.true., vol_dim=[TEST_BOX, TEST_BOX, TEST_BOX])
        call molecule%kill
        call truth_image%new([TEST_BOX, TEST_BOX, TEST_BOX], TEST_SMPD, wthreads=.false.)
        call truth_image%read(string(ATOM_VOLUME_FILE))
        call truth_image%mask3D_soft(TEST_MSKRAD, backgr=0.)
        call truth_image%write(string(TRUTH_VOLUME_FILE), del_if_exists=.true.)

        ! The builder's calculator needs no discrete orientation bank. Both halves intentionally
        ! hold the same known truth reference for this C1 fixture (set per case below);
        ! one particle slot per particle, particle i in slot i (the batch is the whole set).
        truth_volume = truth_image%get_rmat()
        calc => b%cftc
        call calc%new(1, TEST_BOX, TEST_PARTICLES)
        call calc%set_ptcl_inds([(i, i=1, TEST_PARTICLES)])
        taper = calc%get_ptcl_taper()

        ! Stage 2: freeze orientations and CTFs, then ask SIMPLE's simulator to
        ! generate the corresponding clean observations. Noise is added below from
        ! a separate deterministic stream because parameters%new reseeds SIMPLE.
        call truth_orientations%new(TEST_PARTICLES, is_ptcl=.true.)
        call set_fixed_seed(SIMULATION_SEED)
        call truth_orientations%rnd_oris(0.)
        call truth_orientations%rnd_ctf(TEST_KV, TEST_CS, TEST_FRACA, TEST_DEFOCUS, &
            &TEST_DFERR, TEST_ASTIGERR)
        call truth_orientations%set_all2single('state', 1.)
        call truth_orientations%write(string(SIMULATION_ORIENTATION_FILE), [1, TEST_PARTICLES])

        call simulation_command%set('prg', 'simulate_particles')
        call simulation_command%set('mkdir', 'no')
        call simulation_command%set('vol1', ATOM_VOLUME_FILE)
        call simulation_command%set('outstk', PARTICLE_FILE)
        call simulation_command%set('oritab', SIMULATION_ORIENTATION_FILE)
        call simulation_command%set('outfile', TRUTH_ORIENTATION_FILE)
        call simulation_command%set('nptcls', TEST_PARTICLES)
        call simulation_command%set('nthr', nworkers)
        call simulation_command%set('smpd', TEST_SMPD)
        call simulation_command%set('mskdiam', TEST_MSKDIAM)
        call simulation_command%set('pgrp', 'c1')
        call simulation_command%set('ctf', 'yes')
        ! snr>=5 suppresses simulator noise; deterministic noise is added on read.
        call simulation_command%set('snr', 5.)
        call simulation_command%set('kv', TEST_KV)
        call simulation_command%set('cs', TEST_CS)
        call simulation_command%set('fraca', TEST_FRACA)
        call simulation_command%set('defocus', TEST_DEFOCUS)
        call simulation_command%set('dferr', TEST_DFERR)
        call simulation_command%set('astigerr', TEST_ASTIGERR)
        call simulation_command%set('bfac', 0.)
        call simulation_command%set('bfacerr', 0.)
        call simulation_command%set('sherr', 0.)
        call simulator%execute(simulation_command)
        call simulation_command%kill

        call find_ldim_nptcls(string(PARTICLE_FILE), ldim, nsections)
        if (any(ldim(1:2) /= [TEST_BOX, TEST_BOX]) .or. nsections /= TEST_PARTICLES) &
            &THROW_HARD('1JYX simulator returned an incompatible particle stack')

        ! oris%read fills existing records; it does not allocate the destination.
        call truth_orientations%kill
        call truth_orientations%new(TEST_PARTICLES, is_ptcl=.true.)
        call truth_orientations%read(string(TRUTH_ORIENTATION_FILE), [1, TEST_PARTICLES])
        if (truth_orientations%get_noris() /= TEST_PARTICLES) &
            &THROW_HARD('1JYX simulator did not return the truth orientations')

        ! Stage 3: retain the finite particle stack and construct controlled seeds.
        allocate (particles(TEST_BOX, TEST_BOX, TEST_PARTICLES))
        call particle_reader%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD, wthreads=.false.)
        call set_fixed_seed(NOISE_SEED)
        do i = 1, TEST_PARTICLES
            call particle_reader%read(string(PARTICLE_FILE), i)
            call particle_reader%add_gauran(TEST_SNR)
            ! Copy only the logical image region; image%rmat may include FFT padding.
            call particle_reader%get_rmat_sub(particles(:, :, i:i))
        end do
        call particle_reader%kill

        allocate (ctf_parameters(TEST_PARTICLES))
        allocate (truth_rotations(3, 3, TEST_PARTICLES), truth_shifts(2, TEST_PARTICLES))
        allocate (initial_rotations(3, 3, TEST_PARTICLES), initial_shifts(2, TEST_PARTICLES))
        call set_fixed_seed(PERTURBATION_SEED)
        do i = 1, TEST_PARTICLES
            call truth_orientations%get_ori(i, truth_orientation)
            ctf_parameters(i) = truth_orientation%get_ctfvars()
            ! A standalone orientation table has no stack-level sampling metadata.
            ! Production obtains it from the project; this fixture owns it directly.
            ctf_parameters(i)%smpd = TEST_SMPD
            if (ctf_parameters(i)%ctfflag /= CTFFLAG_YES) &
                &THROW_HARD('1JYX simulator did not preserve enabled CTF metadata')
            truth_rotations(:, :, i) = real(truth_orientation%get_mat(), dp)
            truth_shifts(:, i) = real(truth_orientation%get_2Dshift(), dp)
            call make_perturbed_pose(truth_rotations(:, :, i), truth_shifts(:, i), &
                &initial_rotations(:, :, i), initial_shifts(:, i))
        end do
        call truth_orientation%kill

        call mask_image%disc([TEST_BOX, TEST_BOX, 1], TEST_SMPD, TEST_MSKRAD, noise_mask)
        call mask_image%kill
        allocate (initial_rotation_errors(TEST_PARTICLES), terminal_rotation_errors(TEST_PARTICLES))
        allocate (initial_shift_errors(TEST_PARTICLES), terminal_shift_errors(TEST_PARTICLES))
        allocate (objectives_before(TEST_PARTICLES), objectives_after(TEST_PARTICLES), objectives_truth(TEST_PARTICLES))
        allocate (statuses(TEST_PARTICLES))
        allocate (terminal_rotations(3, 3, TEST_PARTICLES, NOBJFUNS), terminal_shifts(2, TEST_PARTICLES, NOBJFUNS))

        ! Stage 4: the Euclidean noise model and reference amplitude, as production would have
        ! them at convergence: one least-squares gain of the truth reference to the data and
        ! the per-shell sigma2 (mean squared residual over two per sample) at the truth poses,
        ! over the widest band of the test.
        call calc%set_ref(1, .true., truth_volume)
        call calc%set_ref(1, .false., truth_volume)
        call estimate_noise_model(calc, particles, ctf_parameters, noise_mask, taper, truth_rotations, &
            &truth_shifts, initial_shifts, active_upper_shell(minval(TEST_LPS)), nworkers, gain, sigma2)
        write (logfhandle, '(a,1x,es12.4)') 'POSE_CONT_1JYX euclid reference gain:', gain

        ! Stage 5: refine all particles under each objective at each band, then reconstruct the
        ! exact, perturbed and refined pose sets of the band and score them.
        median_rotation = 0._dp
        do iband = 1, size(TEST_LPS)
            upper_shell = active_upper_shell(TEST_LPS(iband))
            do iobj = 1, NOBJFUNS
                ! coarse to fine (R1): the 8 A case starts from the injected perturbation, the 4 A
                ! case from the poses the same objective refined at 8 A
                if (iband == 1) then
                    start_rotations = initial_rotations
                    start_shifts = initial_shifts
                else
                    start_rotations = terminal_rotations(:, :, :, iobj)
                    start_shifts = terminal_shifts(:, :, iobj)
                end if
                if (TEST_OBJFUNS(iobj) == OBJFUN_EUCLID) then
                    call calc%set_ref(1, .true., real(gain)*truth_volume)
                    call calc%set_ref(1, .false., real(gain)*truth_volume)
                else
                    call calc%set_ref(1, .true., truth_volume)
                    call calc%set_ref(1, .false., truth_volume)
                end if
                call refine_all_particles(b, p, particles, ctf_parameters, noise_mask, taper, TEST_OBJFUNS(iobj), &
                    &sigma2, upper_shell, nworkers, truth_rotations, truth_shifts, start_rotations, start_shifts, &
                    &terminal_rotations(:, :, :, iobj), terminal_shifts(:, :, iobj), objectives_before, &
                    &objectives_after, objectives_truth, statuses)
                do i = 1, TEST_PARTICLES
                    initial_rotation_errors(i) = rotation_distance(start_rotations(:, :, i), truth_rotations(:, :, i))
                    terminal_rotation_errors(i) = rotation_distance(terminal_rotations(:, :, i, iobj), truth_rotations(:, :, i))
                    initial_shift_errors(i) = norm2(start_shifts(:, i) - truth_shifts(:, i))
                    terminal_shift_errors(i) = norm2(terminal_shifts(:, i, iobj) - truth_shifts(:, i))
                end do
                ! diagnostic only: are the particles that end outside the basin width in a local
                ! minimum (objective above its truth value) or below the truth objective?
                basin_width = real(TEST_LPS(iband), dp)/(0.5_dp*real(TEST_MSKDIAM, dp))
                write (logfhandle, '(a,1x,a,a,i0,a,i0,a,f8.4)') 'POSE_CONT_1JYX', '['//trim(case_label(iobj, iband))//']', &
                    &' outside the basin: ', count(terminal_rotation_errors > basin_width), &
                    &'; of these with objective below the truth objective: ', &
                    &count(terminal_rotation_errors > basin_width .and. objectives_after < objectives_truth), &
                    &'; mean truth objective: ', sum(objectives_truth)/real(TEST_PARTICLES, dp)
                call assert_pose_improvement(case_label(iobj, iband), iband, TEST_OBJFUNS(iobj), statuses, &
                    &objectives_before, objectives_after, initial_rotation_errors, terminal_rotation_errors, &
                    &initial_shift_errors, terminal_shift_errors, median_rotation(iobj, iband))
            end do
            call reconstruct_and_score(iband, particles, ctf_parameters, truth_image, noise_mask, &
                &truth_rotations, truth_shifts, initial_rotations, initial_shifts, terminal_rotations, terminal_shifts)
        end do
        ! F5: a finer band refines the rotation more precisely, which needs the high-resolution
        ! signal that the centring and phase flip before the mask keep (O4)
        do iobj = 1, NOBJFUNS
            write (logfhandle, '(a,1x,a,2(1x,f8.3))') 'POSE_CONT_1JYX median rotation error lp 8 / 4 A (deg):', &
                &trim(OBJFUN_NAMES(iobj)), median_rotation(iobj, :)*180._dp/real(PI, dp)
            call assert_true(median_rotation(iobj, 2) < median_rotation(iobj, 1), &
                &trim(OBJFUN_NAMES(iobj))//': the 4 A median rotation error is below the 8 A one')
        end do

        deallocate (truth_volume)
        call calc%kill
        call b%spproj%os_ptcl3D%kill
        nullify (b%spproj_field, calc)
        call truth_orientations%kill
        call truth_image%kill
        call cleanup_1jyx_fixture()
    end subroutine run_pose_cont_1jyx_reconstruction

    subroutine make_perturbed_pose(truth_rotation, truth_shift, rotation, shift)
        real(dp), intent(in) :: truth_rotation(3, 3), truth_shift(2)
        real(dp), intent(out) :: rotation(3, 3), shift(2)
        real(dp) :: axis(3), azimuth, radial, uniform(3), z

        call random_number(uniform)
        z = 2._dp*uniform(1) - 1._dp
        azimuth = 2._dp*real(PI, dp)*uniform(2)
        radial = sqrt(max(0._dp, 1._dp - z*z))
        axis = [radial*cos(azimuth), radial*sin(azimuth), z]
        rotation = right_increment_rotation(truth_rotation, ROTATION_ERROR*axis)
        azimuth = 2._dp*real(PI, dp)*uniform(3)
        shift = truth_shift + SHIFT_ERROR*[cos(azimuth), sin(azimuth)]
    end subroutine make_perturbed_pose

    integer pure function active_upper_shell(lp) result(shell)
        real, intent(in) :: lp
        shell = min(TEST_BOX/2 - 1, int(real(TEST_BOX)*TEST_SMPD/lp))
    end function active_upper_shell

    pure function case_label(iobj, iband) result(label)
        integer, intent(in) :: iobj, iband
        character(len=24) :: label
        write (label, '(a,a,i0,a)') trim(OBJFUN_NAMES(iobj)), ' lp ', nint(TEST_LPS(iband)), ' A'
    end function case_label

    !> Prepare particle iparticle into the calculator slot islot as production prepares it: the
    !! observation centred on the stored (perturbed) shift, phase-flipped, masked and tapered
    !! (prepimg4align_cart, O4 and O5 (a)), for objfun over [2, upper_shell].
    subroutine prepare_particle(calc, particle, ctf_parameters, noise_mask, taper, stored_shift, objfun, sigma2, &
        &upper_shell, raw_work, fourier_work, islot)
        class(cartft_calc), intent(inout) :: calc
        real, intent(in) :: particle(:, :, :), taper(:), sigma2(0:)
        type(ctfparams), intent(in) :: ctf_parameters
        logical, intent(in) :: noise_mask(:, :, :)
        real(dp), intent(in) :: stored_shift(2)
        integer, intent(in) :: objfun, upper_shell, islot
        type(image), intent(inout) :: raw_work, fourier_work
        type(ctfparams) :: cropped_ctf
        complex, allocatable :: observed(:, :)

        call raw_work%set_rmat(particle, .false.)
        call prepimg4align_cart(raw_work, noise_mask, fourier_work, TEST_MSKRAD, TEST_SMPD, &
            &real(stored_shift), taper, ctf_parameters, observed, cropped_ctf)
        if (objfun == OBJFUN_EUCLID) then
            call calc%set_ptcl(islot, observed, cropped_ctf, sigma2, [2, upper_shell])
        else
            call calc%set_ptcl(islot, observed, cropped_ctf, [2, upper_shell])
        end if
    end subroutine prepare_particle

    !> The soft-mask coordinates of a TEST_BOX image (process-global, memoized serially).
    subroutine memoize_mask_coords_of_box()
        type(image) :: mask_coords_image
        call mask_coords_image%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD, wthreads=.false.)
        call mask_coords_image%memoize_mask_coords()
        call mask_coords_image%kill
    end subroutine memoize_mask_coords_of_box

    !> The Euclidean fixture inputs at the truth poses (the state production reaches at
    !! convergence): the least-squares gain g of the truth reference to the observations,
    !! g = sum Re(X* C M)/sum |C M|^2 over the band, and the per-shell sigma2 of the residual
    !! X - g C M, mean squared residual over two per sample. Uses the calculator's unwhitened
    !! per-shell accounting (sigma_contribution) of a unit-sigma2 Euclidean slot.
    subroutine estimate_noise_model(calc, particles, ctf_parameters, noise_mask, taper, truth_rotations, &
        &truth_shifts, initial_shifts, upper_shell, nworkers, gain, sigma2)
        class(cartft_calc), intent(inout) :: calc
        real, intent(in) :: particles(:, :, :), taper(:)
        type(ctfparams), intent(in) :: ctf_parameters(:)
        logical, intent(in) :: noise_mask(:, :, :)
        real(dp), intent(in) :: truth_rotations(:, :, :), truth_shifts(:, :), initial_shifts(:, :)
        integer, intent(in) :: upper_shell, nworkers
        real(dp), intent(out) :: gain
        real, allocatable, intent(out) :: sigma2(:)
        type(image), allocatable :: raw_work(:), fourier_work(:)
        real, allocatable :: unit_sigma2(:)
        real(dp), allocatable :: residual_sum(:, :), ref_sum(:, :), ptcl_sum(:, :)
        real(dp) :: cross(2:upper_shell), ref_mean(2:upper_shell), ptcl_mean(2:upper_shell), npix(2:upper_shell)
        integer :: h, k, shell

        allocate (unit_sigma2(0:TEST_BOX/2), source=1.)
        allocate (residual_sum(2:upper_shell, nworkers), ref_sum(2:upper_shell, nworkers), &
            &ptcl_sum(2:upper_shell, nworkers), source=0._dp)
        ! the fused phase flip and the mask of the observation read the memoized Fourier maps and
        ! mask coordinates of the box, memoized outside the parallel region
        call memoize_ft_maps([TEST_BOX, TEST_BOX, 1], TEST_SMPD)
        call memoize_mask_coords_of_box()
        call new_work_images(nworkers, raw_work, fourier_work)
!$omp parallel default(shared)
        block
            real, allocatable :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
            real :: v
            integer :: iparticle, thread_index

            thread_index = 1
!$          thread_index = omp_get_thread_num() + 1
!$omp do schedule(dynamic,1)
            do iparticle = 1, size(particles, 3)
                call prepare_particle(calc, particles(:, :, iparticle:iparticle), ctf_parameters(iparticle), &
                    &noise_mask, taper, initial_shifts(:, iparticle), OBJFUN_EUCLID, unit_sigma2, upper_shell, &
                    &raw_work(thread_index), fourier_work(thread_index), thread_index)
                ! the truth pose relative to the observation centred on the stored shift
                call calc%sigma_contribution(1, .true., thread_index, truth_rotations(:, :, iparticle), &
                    &truth_shifts(:, iparticle) - initial_shifts(:, iparticle), sigma_contrib, ref_pow, ptcl_pow, v)
                residual_sum(:, thread_index) = residual_sum(:, thread_index) + real(sigma_contrib, dp)
                ref_sum(:, thread_index) = ref_sum(:, thread_index) + real(ref_pow, dp)
                ptcl_sum(:, thread_index) = ptcl_sum(:, thread_index) + real(ptcl_pow, dp)
            end do
!$omp end do
        end block
!$omp end parallel
        call kill_work_images(raw_work, fourier_work)
        call forget_ft_maps()
        ! per-shell means per sample over all particles: |X - CM|^2 = 2 residual, so
        ! Re(X* CM) = (|X|^2 + |CM|^2 - 2 residual)/2
        ref_mean = sum(ref_sum, dim=2)/real(size(particles, 3), dp)
        ptcl_mean = sum(ptcl_sum, dim=2)/real(size(particles, 3), dp)
        cross = (ptcl_mean + ref_mean - 2._dp*sum(residual_sum, dim=2)/real(size(particles, 3), dp))/2._dp
        npix = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if (shell >= 2 .and. shell <= upper_shell) npix(shell) = npix(shell) + 1._dp
            end do
        end do
        gain = sum(npix*cross)/sum(npix*ref_mean)
        allocate (sigma2(0:TEST_BOX/2), source=1.)
        sigma2(2:upper_shell) = real((ptcl_mean - 2._dp*gain*cross + gain*gain*ref_mean)/2._dp)
        if (.not. (ieee_is_finite(gain) .and. gain > 0._dp) .or. any(sigma2(2:upper_shell) <= 0.)) &
            &THROW_HARD('1JYX Euclidean noise model is not positive')
    end subroutine estimate_noise_model

    !> Refine every particle from its seed through strategy3D_cont, as the matcher runs a
    !! Cartesian pass: the seeds in the ptcl3D field, every observation centred on its stored
    !! shift in its slot, one strategy per particle committing the pose when its solve is
    !! accepted. The objectives before and after are re-evaluated on the slot at the seed and
    !! at the committed pose (relative to the stored shift); the status is CARTFT_ACCEPTED for
    !! an improved particle and CARTFT_NO_IMPROVEMENT for a kept seed (the improved flag).
    subroutine refine_all_particles(b, p, particles, ctf_parameters, noise_mask, taper, objfun, sigma2, &
        &upper_shell, nworkers, truth_rotations, truth_shifts, initial_rotations, initial_shifts, terminal_rotations, &
        &terminal_shifts, objectives_before, objectives_after, objectives_truth, statuses)
        type(builder), target, intent(inout) :: b
        type(parameters), target, intent(inout) :: p
        real, intent(in) :: particles(:, :, :), taper(:), sigma2(0:)
        type(ctfparams), intent(in) :: ctf_parameters(:)
        logical, intent(in) :: noise_mask(:, :, :)
        integer, intent(in) :: objfun, upper_shell, nworkers
        real(dp), intent(in) :: truth_rotations(:, :, :), truth_shifts(:, :)
        real(dp), intent(in) :: initial_rotations(:, :, :), initial_shifts(:, :)
        real(dp), intent(out) :: terminal_rotations(:, :, :), terminal_shifts(:, :)
        real(dp), intent(out) :: objectives_before(:), objectives_after(:), objectives_truth(:)
        integer, intent(out) :: statuses(:)
        type(image), allocatable :: raw_work(:), fourier_work(:)
        integer :: i

        ! the pass's policy: the Phase 0 bounds (O6), shift before joint
        p%oritype = 'ptcl3D'
        p%cc_objfun = objfun
        p%inpl_cont = 'no'
        p%box = TEST_BOX
        p%box_crop = TEST_BOX
        p%trs = TEST_TRS
        p%athres_cont = TEST_ATHRES_CONT
        ! the route of the Phase 0 floors; production runs the joint stage alone
        p%cont_route = 'shift_then_joint'
        p%l_cont_shift_first = .true.
        p%kfromto = [2, upper_shell]
        p%fromp = 1
        p%top = size(particles, 3)
        ! the seeds: the C1 fixture refines every particle against the even reference
        call b%spproj%os_ptcl3D%kill
        call b%spproj%os_ptcl3D%new(size(particles, 3), is_ptcl=.true.)
        b%spproj_field => b%spproj%os_ptcl3D
        call b%pgrpsyms%new('c1')
        do i = 1, size(particles, 3)
            call b%spproj_field%set_euler(i, real(dm2euler(initial_rotations(:, :, i))))
            call b%spproj_field%set_shift(i, real(initial_shifts(:, i)))
            call b%spproj_field%set_state(i, 1)
            call b%spproj_field%set(i, 'eo', 0.)
            call b%spproj_field%set(i, 'proj', 1.)
        end do
        ! sigma2 belongs to the Euclidean objective (C5): the strategy records the residual
        if (objfun == OBJFUN_EUCLID) then
            call b%esig%new(p, string('tmp_strategy3D_cont_1jyx_sigma2.bin'), TEST_BOX)
            call b%esig%allocate_ptcls
        end if
        call memoize_ft_maps([TEST_BOX, TEST_BOX, 1], TEST_SMPD)
        call memoize_mask_coords_of_box()
        call new_work_images(nworkers, raw_work, fourier_work)
!$omp parallel default(shared)
        block
            class(strategy3D), pointer :: strat
            type(strategy3D_spec) :: spec
            type(ori) :: seed
            real(dp) :: gradient(5)
            integer :: iparticle, thread_index

            thread_index = 1
!$          thread_index = omp_get_thread_num() + 1
!$omp do schedule(dynamic,1)
            do iparticle = 1, size(particles, 3)
                call prepare_particle(b%cftc, particles(:, :, iparticle:iparticle), ctf_parameters(iparticle), &
                    &noise_mask, taper, initial_shifts(:, iparticle), objfun, sigma2, upper_shell, &
                    &raw_work(thread_index), fourier_work(thread_index), iparticle)
                call b%spproj_field%get_ori(iparticle, seed)
                call b%cftc%objective_gradient(1, .true., iparticle, real(seed%get_mat(), dp), [0._dp, 0._dp], &
                    &objectives_before(iparticle), gradient)
                spec%iptcl = iparticle
                spec%iptcl_map = iparticle
                allocate (strategy3D_cont :: strat)
                call strat%new(p, spec, b)
                call strat%srch(b%spproj_field, thread_index)
                call strat%kill
                deallocate (strat)
                terminal_rotations(:, :, iparticle) = real(b%spproj_field%get_mat(iparticle), dp)
                terminal_shifts(:, iparticle) = real(b%spproj_field%get_2Dshift(iparticle), dp)
                call b%cftc%objective_gradient(1, .true., iparticle, terminal_rotations(:, :, iparticle), &
                    &terminal_shifts(:, iparticle) - real(seed%get_2Dshift(), dp), objectives_after(iparticle), gradient)
                statuses(iparticle) = merge(CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT, &
                    &b%spproj_field%get(iparticle, 'pose_cont_improved') > 0.5)
                ! diagnostic: the objective at the truth pose, relative to the centred observation
                call b%cftc%objective_gradient(1, .true., iparticle, truth_rotations(:, :, iparticle), &
                    &truth_shifts(:, iparticle) - real(seed%get_2Dshift(), dp), objectives_truth(iparticle), gradient)
            end do
!$omp end do
            call seed%kill
        end block
!$omp end parallel
        call kill_work_images(raw_work, fourier_work)
        call forget_ft_maps()
        if (objfun == OBJFUN_EUCLID) call b%esig%kill
        call b%pgrpsyms%kill
    end subroutine refine_all_particles

    !> One raw and one Fourier work image per OpenMP worker, created serially (image%new
    !! makes FFTW plans, which is not thread-safe).
    subroutine new_work_images(nworkers, raw_work, fourier_work)
        integer, intent(in) :: nworkers
        type(image), allocatable, intent(out) :: raw_work(:), fourier_work(:)
        integer :: i
        allocate (raw_work(nworkers), fourier_work(nworkers))
        do i = 1, nworkers
            call raw_work(i)%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD, wthreads=.false.)
            call fourier_work(i)%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD, wthreads=.false.)
        end do
    end subroutine new_work_images

    subroutine kill_work_images(raw_work, fourier_work)
        type(image), allocatable, intent(inout) :: raw_work(:), fourier_work(:)
        integer :: i
        do i = 1, size(raw_work)
            call raw_work(i)%kill
            call fourier_work(i)%kill
        end do
        deallocate (raw_work, fourier_work)
    end subroutine kill_work_images

    pure real(dp) function rotation_distance(rotation_a, rotation_b) result(distance)
        real(dp), intent(in) :: rotation_a(3, 3), rotation_b(3, 3)
        real(dp) :: relative_rotation(3, 3), cosine

        relative_rotation = matmul(transpose(rotation_a), rotation_b)
        cosine = 0.5_dp*(relative_rotation(1, 1) + relative_rotation(2, 2) + &
            &relative_rotation(3, 3) - 1._dp)
        distance = acos(max(-1._dp, min(1._dp, cosine)))
    end function rotation_distance

    !> Reconstruct the exact, perturbed and (per objective) refined pose sets of band iband
    !! with matched gridding accumulators, compare every map with the same masked truth by
    !! correlation and FSC, and apply the map gates: the exact-pose control is valid, each
    !! refined map beats the perturbed one and lies within MAP_CORR_MARGIN of the exact-pose map
    !! (F4), and for cc at 8 A it meets the Phase 0 correlation with its margin (F3).
    subroutine reconstruct_and_score(iband, particles, ctf_parameters, truth_image, noise_mask, &
        &truth_rotations, truth_shifts, initial_rotations, initial_shifts, terminal_rotations, terminal_shifts)
        integer, intent(in) :: iband
        real, intent(in) :: particles(:, :, :)
        type(ctfparams), intent(in) :: ctf_parameters(:)
        class(image), intent(in) :: truth_image
        logical, intent(in) :: noise_mask(:, :, :)
        real(dp), intent(in) :: truth_rotations(:, :, :), truth_shifts(:, :)
        real(dp), intent(in) :: initial_rotations(:, :, :), initial_shifts(:, :)
        real(dp), intent(in) :: terminal_rotations(:, :, :, :), terminal_shifts(:, :, :)
        integer, parameter :: NMAPS = 2 + NOBJFUNS !< exact, perturbed, refined per objective
        type(fplane_type) :: plane
        type(image) :: particle, padded_particle, maps(NMAPS), truth_metric, metric
        type(ori) :: orientation
        type(parameters), target :: reconstruction_parameters
        type(reconstructor) :: reconstructors(NMAPS)
        type(sp_project) :: empty_project
        type(sym) :: c1
        real, allocatable :: map_array(:, :, :), truth_array(:, :, :), fsc(:), resolutions(:)
        real(dp) :: correlations(NMAPS)
        real :: fsc05(NMAPS), fsc0143(NMAPS)
        integer :: i, imap, upper_shell
        logical :: finite_fsc(NMAPS)
        character(len=16) :: map_names(NMAPS)

        upper_shell = active_upper_shell(TEST_LPS(iband))
        map_names(1) = 'exact'
        map_names(2) = 'perturbed'
        do imap = 1, NOBJFUNS
            map_names(2 + imap) = 'refined '//trim(OBJFUN_NAMES(imap))
        end do
        ! Step 1: matched gridding accumulators; the exact-pose arm validates the harness.
        reconstruction_parameters%box = TEST_BOX
        reconstruction_parameters%box_crop = TEST_BOX
        reconstruction_parameters%box_croppd = OSMPL_PAD_FAC*TEST_BOX
        reconstruction_parameters%smpd_crop = TEST_SMPD
        reconstruction_parameters%nstates = 1
        reconstruction_parameters%numlen = 1
        reconstruction_parameters%oritype = 'cls3D'
        do imap = 1, NMAPS
            call reconstructors(imap)%new_accumulator(reconstruction_parameters, empty_project, &
                &expand=.true., wthreads=.false.)
            call reconstructors(imap)%set_ft(.true.)
            call maps(imap)%new([TEST_BOX, TEST_BOX, TEST_BOX], TEST_SMPD, wthreads=.false.)
        end do
        call particle%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD, wthreads=.false.)
        call padded_particle%new([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], &
            &TEST_SMPD, wthreads=.false.)
        call orientation%new(.false.)
        call c1%new('c1')
        call memoize_ft_maps([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, 1], TEST_SMPD)

        ! Step 2: insert each noisy observation into every accumulator; only the pose differs.
        do i = 1, size(particles, 3)
            call particle%set_rmat(particles(:, :, i:i), .false.)
            call particle%norm_noise_taper_edge_pad_fft(noise_mask, padded_particle)
            call insert(1, truth_rotations(:, :, i), truth_shifts(:, i))
            call insert(2, initial_rotations(:, :, i), initial_shifts(:, i))
            do imap = 1, NOBJFUNS
                call insert(2 + imap, terminal_rotations(:, :, i, imap), terminal_shifts(:, i, imap))
            end do
        end do

        ! Step 3: finalize the maps and compare each with the same masked truth.
        call truth_metric%copy(truth_image)
        call truth_metric%mask3D_soft(TEST_MSKRAD, backgr=0.)
        truth_array = truth_metric%get_rmat()
        call truth_metric%fft()
        allocate (fsc(truth_metric%get_filtsz()), resolutions(truth_metric%get_filtsz()))
        do i = 1, size(resolutions)
            resolutions(i) = real(TEST_BOX)*TEST_SMPD/real(i)
        end do
        do imap = 1, NMAPS
            call reconstructors(imap)%compress_exp()
            call reconstructors(imap)%restore_final(maps(imap))
            call metric%copy(maps(imap))
            call metric%mask3D_soft(TEST_MSKRAD, backgr=0.)
            map_array = metric%get_rmat()
            correlations(imap) = centered_array_correlation(map_array, truth_array)
            call metric%fft()
            fsc = 0.
            call truth_metric%fsc(metric, fsc)
            finite_fsc(imap) = all(ieee_is_finite(fsc))
            call get_resolution(fsc, resolutions, fsc05(imap), fsc0143(imap))
            call metric%kill
            write (logfhandle, '(a,i0,a,a16,a,es12.4,a,2(1x,f8.3))') 'POSE_CONT_1JYX [lp ', nint(TEST_LPS(iband)), &
                &' A] map ', map_names(imap), ' correlation', correlations(imap), '  FSC=0.5/0.143 (A)', &
                &fsc05(imap), fsc0143(imap)
        end do

        ! Step 4: the gates, after every diagnostic is written.
        call assert_true(ieee_is_finite(correlations(1)) .and. correlations(1) > 0._dp .and. finite_fsc(1) .and. &
            &fsc0143(1) > 0., 'the exact-pose reconstruction control is valid')
        do imap = 1, NOBJFUNS
            call assert_true(ieee_is_finite(correlations(2 + imap)) .and. finite_fsc(2) .and. finite_fsc(2 + imap) .and. &
                &correlations(2 + imap) > correlations(2), &
                &trim(case_label(imap, iband))//': the refined map correlates better with the truth than the perturbed one')
            call assert_true(correlations(2 + imap) >= correlations(1) - MAP_CORR_MARGIN, &
                &trim(case_label(imap, iband))//': the refined map is within the margin of the exact-pose map (F4)')
            if (TEST_OBJFUNS(imap) == OBJFUN_CC .and. iband == 1) &
                &call assert_true(correlations(2 + imap) >= PHASE0_MAP_CORR - 0.005_dp, &
                &trim(case_label(imap, iband))//': the refined map meets the Phase 0 correlation with its margin (F3)')
        end do

        if (allocated(plane%cmplx_plane)) deallocate (plane%cmplx_plane)
        if (allocated(plane%ctfsq_plane)) deallocate (plane%ctfsq_plane)
        if (allocated(plane%transfer_plane)) deallocate (plane%transfer_plane)
        call forget_ft_maps()
        call c1%kill
        call orientation%kill
        call truth_metric%kill
        do imap = 1, NMAPS
            call maps(imap)%kill
            call reconstructors(imap)%kill
        end do
        call padded_particle%kill
        call particle%kill
        call empty_project%kill

    contains

        subroutine insert(imap_insert, rotation, shift)
            integer, intent(in) :: imap_insert
            real(dp), intent(in) :: rotation(3, 3), shift(2)
            call orientation%set_euler(real(dm2euler(rotation)))
            call padded_particle%gen_fplane4rec([2, upper_shell], TEST_SMPD, ctf_parameters(i), real(shift), plane)
            call reconstructors(imap_insert)%insert_plane_oversamp(c1, orientation, plane)
        end subroutine insert

    end subroutine reconstruct_and_score

    subroutine cleanup_1jyx_fixture()
        call del_file(ATOM_VOLUME_FILE)
        call del_file(TRUTH_VOLUME_FILE)
        call del_file(PARTICLE_FILE)
        call del_file(SIMULATION_ORIENTATION_FILE)
        call del_file(TRUTH_ORIENTATION_FILE)
    end subroutine cleanup_1jyx_fixture

    !> Print the metrics of one case and apply its pose gates: the existing ones (an accepted
    !! particle; the mean objective, the RMS rotation and the RMS shift error fall) and the
    !! floors of section 10.1 (F1 truth, F2 analytic, F3 Phase 0 for cc at 8 A). Returns the
    !! median terminal rotation error for F5.
    subroutine assert_pose_improvement(label, iband, objfun, statuses, objectives_before, objectives_after, &
        &rotation_before, rotation_after, shift_before, shift_after, rotation_median)
        character(len=*), intent(in) :: label
        integer, intent(in) :: iband, objfun
        integer, intent(in) :: statuses(:)
        real(dp), intent(in) :: objectives_before(:), objectives_after(:)
        real(dp), intent(in) :: rotation_before(:), rotation_after(:), shift_before(:), shift_after(:)
        real(dp), intent(out) :: rotation_median
        logical, allocatable :: finite_objective(:)
        real(dp) :: objective_before_mean, objective_after_mean
        real(dp) :: rotation_before_rms, rotation_after_rms, shift_before_rms, shift_after_rms
        real(dp) :: basin_width, shift_median, basin_fraction
        integer :: accepted, status_code
        character(len=40) :: tag

        tag = 'POSE_CONT_1JYX ['//trim(label)//']'
        finite_objective = ieee_is_finite(objectives_before) .and. ieee_is_finite(objectives_after) .and. &
            &objectives_before >= 0._dp .and. objectives_after >= 0._dp
        call assert_int(size(statuses), count(finite_objective), &
            &trim(label)//': every particle has finite, non-negative objectives before and after')
        accepted = count(statuses == CARTFT_ACCEPTED)
        objective_before_mean = sum(objectives_before)/real(size(statuses), dp)
        objective_after_mean = sum(objectives_after)/real(size(statuses), dp)
        rotation_before_rms = sqrt(sum(rotation_before**2)/real(size(statuses), dp))
        rotation_after_rms = sqrt(sum(rotation_after**2)/real(size(statuses), dp))
        shift_before_rms = sqrt(sum(shift_before**2)/real(size(statuses), dp))
        shift_after_rms = sqrt(sum(shift_after**2)/real(size(statuses), dp))
        rotation_median = real(median(real(rotation_after)), dp)
        shift_median = real(median(real(shift_after)), dp)
        basin_width = real(TEST_LPS(iband), dp)/(0.5_dp*real(TEST_MSKDIAM, dp))
        basin_fraction = real(count(rotation_after <= basin_width), dp)/real(size(statuses), dp)

        write (logfhandle, '(a,1x,a,i0,a,i0)') trim(tag), 'accepted particles: ', accepted, ' / ', size(statuses)
        write (logfhandle, '(a,1x,a,2(1x,es12.4))') trim(tag), 'mean objective before/after:', &
            &objective_before_mean, objective_after_mean
        write (logfhandle, '(a,1x,a,2(1x,f8.3))') trim(tag), 'rotation RMS before/after (deg):', &
            &rotation_before_rms*180._dp/real(PI, dp), rotation_after_rms*180._dp/real(PI, dp)
        write (logfhandle, '(a,1x,a,2(1x,f8.3))') trim(tag), 'shift RMS before/after (pixels):', &
            &shift_before_rms, shift_after_rms
        write (logfhandle, '(a,1x,a,2(1x,f8.3))') trim(tag), 'rotation median before/after (deg):', &
            &median(real(rotation_before))*180./PI, rotation_median*180._dp/real(PI, dp)
        write (logfhandle, '(a,1x,a,2(1x,f8.3))') trim(tag), 'shift median before/after (pixels):', &
            &median(real(shift_before)), shift_median
        write (logfhandle, '(a,1x,a,2(1x,f8.4))') trim(tag), 'fraction with lower rotation/shift error:', &
            &real(count(rotation_after < rotation_before))/real(size(statuses)), &
            &real(count(shift_after < shift_before))/real(size(statuses))
        write (logfhandle, '(a,1x,a,f8.3,a,2(1x,f8.4))') trim(tag), 'fraction within basin width ', &
            &basin_width*180._dp/real(PI, dp), ' deg before/after:', &
            &real(count(rotation_before <= basin_width))/real(size(statuses)), basin_fraction
        do status_code = minval(statuses), maxval(statuses)
            if (count(statuses == status_code) == 0) cycle
            write (logfhandle, '(a,1x,a,i0,a,i0)') trim(tag), 'LM status ', status_code, ' count: ', &
                &count(statuses == status_code)
        end do

        call assert_true(accepted >= 1, trim(label)//': at least one particle was accepted')
        call assert_true(objective_after_mean < objective_before_mean, trim(label)//': the aggregate objective falls')
        call assert_true(rotation_after_rms < rotation_before_rms, trim(label)//': the aggregate rotation error falls')
        call assert_true(shift_after_rms < shift_before_rms, trim(label)//': the aggregate shift error falls')
        ! F1: the median errors fall to a tenth of the injected 15 degrees and 2 pixels
        call assert_true(rotation_median <= TRUTH_FACTOR*ROTATION_ERROR, &
            &trim(label)//': median rotation error at most a tenth of the injected error (F1)')
        call assert_true(shift_median <= TRUTH_FACTOR*SHIFT_ERROR, &
            &trim(label)//': median shift error at most a tenth of the injected error (F1)')
        ! F2: the particles end inside the basin width lp/(mskdiam/2) at the band limit
        call assert_true(basin_fraction >= BASIN_FRACTION_FLOOR(iband), &
            &trim(label)//': fraction inside the basin width meets its floor (F2)')
        ! F3: the pre-refactor implementation on the same fixture, with margins
        if (objfun == OBJFUN_CC .and. iband == 1) then
            call assert_true(rotation_median <= 2._dp*PHASE0_MEDIAN_ROTATION_DEG*real(PI, dp)/180._dp, &
                &trim(label)//': median rotation error within twice the Phase 0 value (F3)')
            call assert_true(basin_fraction >= PHASE0_BASIN_FRACTION - 0.01_dp, &
                &trim(label)//': basin fraction within 0.01 of the Phase 0 value (F3)')
        end if
    end subroutine assert_pose_improvement

    pure real(dp) function centered_array_correlation(array_a, array_b) result(correlation)
        real, intent(in) :: array_a(:, :, :), array_b(:, :, :)
        real(dp) :: mean_a, mean_b, denominator

        mean_a = sum(real(array_a, dp))/real(size(array_a), dp)
        mean_b = sum(real(array_b, dp))/real(size(array_b), dp)
        denominator = sqrt(sum((real(array_a, dp) - mean_a)**2)* &
            &sum((real(array_b, dp) - mean_b)**2))
        correlation = sum((real(array_a, dp) - mean_a)*(real(array_b, dp) - mean_b))/ &
            &max(denominator, epsilon(denominator))
    end function centered_array_correlation

end module simple_strategy3D_cont_1jyx_tester
