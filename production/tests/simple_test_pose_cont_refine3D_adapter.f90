program simple_test_pose_cont_refine3D_adapter
use simple_core_module_api, only: CTFFLAG_NO, ctfparams, dp, euler2m
use simple_image, only: image
use simple_ori, only: ori
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner
use simple_strategy3D_pose_cont, only: pose_cont_seed_is_valid
use simple_pose_cont_refine3D_adapter, only: pose_cont_reference_workspace, &
    &pose_cont_particle_workspace, &
    &pose_cont_pose, pose_cont_limits, pose_cont_config, pose_cont_particle_spec, &
    &pose_cont_observation_spec, pose_cont_observation, &
    &pose_cont_transaction_result, &
    &cartesian_pose_data, &
    &write_pose_cont_reference_artifact, &
    &remove_pose_cont_reference_artifacts, prepare_pose_cont_observation, &
    &pose_cont_seed_from_orientation, pose_cont_pose_to_orientation, &
    &shift_native_to_crop, shift_crop_to_native, pose_cont_config_from_route, &
    &pose_cont_limits_from_boxes, &
    &POSE_CONT_INVALID_PREPARATION, LM_ACCEPTED_IMPROVEMENT, &
    &LM_FINITE_NO_IMPROVEMENT, LM_STEP_BOUND_REJECTED, POSE_CONT_NOT_ATTEMPTED, &
    &POSE_CONT_ROUTE_SHIFT_THEN_JOINT, POSE_CONT_ROUTE_JOINT
use pose_cont_refine3D_adapter_1jyx_test, only: run_pose_cont_1jyx_reconstruction
implicit none

integer, parameter :: TEST_BOX = 16
real, parameter :: TEST_SMPD = 1.5
character(len=32) :: selected_case
integer :: occurrences

call find_selected_case(selected_case, occurrences)
if (occurrences > 1) error stop 'pose-cont adapter suite accepts only one case= argument'
if (occurrences == 0) then
    call run_adapter_contracts()
else
    select case (trim(selected_case))
    case ('adapter')
        call run_adapter_contracts()
    case ('1jyx_reconstruction')
        call run_pose_cont_1jyx_reconstruction()
    case default
        error stop 'pose-cont adapter suite requires case=adapter or case=1jyx_reconstruction'
    end select
end if

contains

    subroutine run_adapter_contracts()
        call test_reference_workspace_lifecycle()
        call test_observation_and_coordinate_adapters()
        call test_inpl_pose_cont_handoff()
        call test_transaction_contracts()
        call test_strategy_seed_contract()
        write (*, '(a)') 'POSE_CONT_REFINE3D_ADAPTER: PASS'
    end subroutine run_adapter_contracts

    subroutine find_selected_case(case_name, count)
        character(len=*), intent(out) :: case_name
        integer, intent(out) :: count
        character(len=256) :: argument
        integer :: iarg, separator, status

        case_name = ''
        count = 0
        do iarg = 1, command_argument_count()
            call get_command_argument(iarg, argument, status=status)
            if (status /= 0) error stop 'could not read pose-cont adapter test argument'
            separator = index(argument, '=')
            if (separator <= 1) cycle
            if (trim(argument(:separator - 1)) /= 'case') cycle
            count = count + 1
            case_name = trim(argument(separator + 1:))
        end do
    end subroutine find_selected_case

    subroutine test_reference_workspace_lifecycle()
        type(pose_cont_reference_workspace) :: workspace
        real, allocatable :: even_volume(:, :, :), odd_volume(:, :, :)

        call build_test_volume(even_volume)
        odd_volume = -0.5*even_volume
        call write_reference(even_volume, .true.)
        call write_reference(odd_volume, .false.)

        call workspace%new_from_artifacts(1, TEST_BOX, TEST_SMPD)
        call assert_true(workspace%is_ready(1, .true.) .and. workspace%is_ready(1, .false.), &
            &'reference workspace did not load both half-set slots')
        call assert_true(.not. workspace%is_ready(0, .true.) .and. &
            &.not. workspace%is_ready(2, .false.), 'reference workspace accepted an invalid state')

        ! Rebuilding must safely replace the existing slot arrays and refiners.
        call workspace%new_from_artifacts(1, TEST_BOX, TEST_SMPD)
        call assert_true(workspace%is_ready(1, .true.) .and. workspace%is_ready(1, .false.), &
            &'reference workspace could not be rebuilt')
        call workspace%kill()
        call assert_true(.not. workspace%is_ready(1, .true.) .and. &
            &.not. workspace%is_ready(1, .false.), 'killed reference workspace remained ready')
        call workspace%kill()
        call remove_pose_cont_reference_artifacts(1)
    end subroutine test_reference_workspace_lifecycle

    subroutine test_observation_and_coordinate_adapters()
        integer, parameter :: NATIVE_BOX = 32
        type(image) :: raw, oracle, work, oracle_work, preserved_work
        type(pose_cont_particle_workspace) :: particles
        type(ctfparams) :: input_ctf, output_ctf
        type(pose_cont_observation_spec) :: observation_spec
        type(pose_cont_observation) :: preserved_observation
        logical, allocatable :: noise_mask(:, :, :)
        complex, allocatable :: observed(:, :), expected(:, :)
        real :: native_shift(2), crop_shift(2)
        integer :: i, j

        call raw%new([NATIVE_BOX, NATIVE_BOX, 1], TEST_SMPD*real(TEST_BOX)/real(NATIVE_BOX))
        do j = 1, NATIVE_BOX
            do i = 1, NATIVE_BOX
                call raw%set_rmat_at(i, j, 1, sin(0.21*real(i)) + cos(0.16*real(j)) + 0.01*real(i*j))
            end do
        end do
        call oracle%copy(raw)
        call work%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD)
        call oracle_work%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD)
        call preserved_work%new([TEST_BOX, TEST_BOX, 1], TEST_SMPD)
        call work%memoize_mask_coords()
        call preserved_work%memoize_mask_coords()
        allocate (noise_mask(NATIVE_BOX, NATIVE_BOX, 1), source=.false.)
        input_ctf%smpd = raw%get_smpd()
        input_ctf%ctfflag = CTFFLAG_NO

        ! Preserve the real-space particle before the established preparation
        ! mutates the caller-owned image for PFTC matching.
        call particles%new(1, NATIVE_BOX, raw%get_smpd())
        call particles%capture(1, raw)
        call prepare_pose_cont_observation(raw, noise_mask, work, 6., TEST_SMPD, &
            &input_ctf, observed, output_ctf)
        observation_spec = pose_cont_observation_spec(batch_index=1, mask_radius=6., &
            &smpd_crop=TEST_SMPD, ctfparms=input_ctf)
        call particles%prepare_observation(observation_spec, noise_mask, preserved_work, &
            &preserved_observation)
        call oracle%norm_noise_fft_clip_shift(noise_mask, oracle_work, [0., 0.])
        call oracle_work%ifft_mask_fft(6.)
        expected = oracle_work%expand_ft()
        call assert_true(maxval(abs(observed - expected)) <= 2.e-5, &
            &'observation adapter disagrees with the established particle path')
        call assert_true(all(lbound(observed) == [-TEST_BOX/2, -TEST_BOX/2]) .and. &
            &all(ubound(observed) == [TEST_BOX/2, TEST_BOX/2]), &
            &'observation adapter did not return the full redundant disk')
        call assert_true(abs(output_ctf%smpd - TEST_SMPD) <= epsilon(TEST_SMPD), &
            &'observation adapter did not update CTF sampling for the cropped box')
        call assert_true(maxval(abs(preserved_observation%samples - expected)) <= 2.e-5, &
            &'preserved particle workspace changed the Cartesian observation')
        call assert_true(abs(preserved_observation%ctfparms%smpd - TEST_SMPD) <= &
            &epsilon(TEST_SMPD), &
            &'preserved particle workspace did not update cropped CTF sampling')

        native_shift = [2.25, -1.75]
        crop_shift = shift_native_to_crop(native_shift, NATIVE_BOX, TEST_BOX)
        call assert_true(maxval(abs(crop_shift - [1.125, -0.875])) <= epsilon(1.), &
            &'native-to-cropped shift conversion used the wrong scale')
        call assert_true(maxval(abs(shift_crop_to_native(crop_shift, NATIVE_BOX, TEST_BOX) - native_shift)) &
            &<= epsilon(1.), 'cropped-to-native shift conversion is not reversible')
        call raw%kill()
        call oracle%kill()
        call work%kill()
        call oracle_work%kill()
        call preserved_work%kill()
        call particles%kill()
        call particles%kill()
    end subroutine test_observation_and_coordinate_adapters

    ! The inpl_cont winner is stored as a fractional SIMPLE Euler pose and a
    ! native-pixel shift. Verify the exact conversion into Cartesian LM
    ! coordinates and the accepted-pose round trip back to project units.
    subroutine test_inpl_pose_cont_handoff()
        integer, parameter :: NATIVE_BOX = 256, CROP_BOX = 144
        type(ori) :: inpl_winner, recovered
        type(pose_cont_pose) :: seed
        real(dp) :: expected_rotation(3, 3)
        real :: eulers(3, 2), native_shifts(2, 2), expected_crop_shift(2)
        integer :: icase

        eulers(:, 1) = [23.5, 67.25, 14.375]
        eulers(:, 2) = [201.125, 88.75, 359.625]
        native_shifts(:, 1) = [2.25, -1.75]
        native_shifts(:, 2) = [-3.5, 0.625]

        do icase = 1, 2
            call inpl_winner%set_euler(eulers(:, icase))
            call inpl_winner%set_shift(native_shifts(:, icase))
            call inpl_winner%set('state', 2.)
            call inpl_winner%set('eo', 1.)
            call inpl_winner%set('corr', 0.42)
            call inpl_winner%set('proj', 7.)
            call inpl_winner%set('inpl', 13.)
            expected_rotation = real(inpl_winner%get_mat(), dp)
            expected_crop_shift = shift_native_to_crop(native_shifts(:, icase), &
                &NATIVE_BOX, CROP_BOX)

            call pose_cont_seed_from_orientation(inpl_winner, NATIVE_BOX, CROP_BOX, seed)
            call assert_true(maxval(abs(seed%rotmat - expected_rotation)) <= 1.e-12_dp, &
                &'inpl_cont Euler pose changed while constructing the pose_cont seed')
            call assert_true(maxval(abs(seed%shift - real(expected_crop_shift, dp))) <= 1.e-12_dp, &
                &'inpl_cont native shift was not converted to cropped pixels')

            call recovered%set('state', 2.)
            call recovered%set('eo', 1.)
            call recovered%set('corr', 0.42)
            call recovered%set('proj', 7.)
            call recovered%set('inpl', 13.)
            call pose_cont_pose_to_orientation(seed, NATIVE_BOX, CROP_BOX, recovered)
            call assert_true(maxval(abs(real(recovered%get_mat(), dp) - expected_rotation)) <= 2.e-6_dp, &
                &'pose_cont rotation did not round-trip to the inpl_cont winner')
            call assert_true(maxval(abs(recovered%get_2Dshift() - native_shifts(:, icase))) <= 2.e-6, &
                &'pose_cont shift did not round-trip to native pixels')
            call assert_true(recovered%get_state() == 2 .and. recovered%get_eo() == 1 .and. &
                &abs(recovered%get('corr') - 0.42) <= epsilon(1.) .and. &
                &nint(recovered%get('proj')) == 7 .and. nint(recovered%get('inpl')) == 13, &
                &'pose handoff changed authoritative non-pose metadata')
        end do
        call inpl_winner%kill()
        call recovered%kill()
    end subroutine test_inpl_pose_cont_handoff

    subroutine test_transaction_contracts()
        type(cartesian_pose_refiner) :: generator
        type(pose_cont_reference_workspace) :: workspace
        type(cartesian_pose_data) :: data, invalid_data
        type(pose_cont_transaction_result) :: result
        type(pose_cont_pose) :: seed
        type(pose_cont_config) :: config
        type(pose_cont_limits) :: limits
        type(pose_cont_particle_spec) :: particle_spec
        type(ctfparams) :: no_ctf
        real, allocatable :: volume(:, :, :), sigma2_noise(:, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: truth_rotation(3, 3), initial_rotation(3, 3), truth_shift(2), initial_shift(2)

        call build_test_volume(volume)
        call write_reference(volume, .true.)
        call write_reference(volume, .false.)
        call workspace%new_from_artifacts(1, TEST_BOX, TEST_SMPD)
        call generator%new_physical_reference(volume)
        truth_rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift = [0.31_dp, -0.24_dp]
        call generator%predict_unweighted(truth_rotation, truth_shift, observed)
        allocate (sigma2_noise(2:TEST_BOX/2 - 1, 5:5), source=1.)
        no_ctf%smpd = TEST_SMPD
        no_ctf%ctfflag = CTFFLAG_NO
        particle_spec = pose_cont_particle_spec(state=1, particle=5, &
            &shell_range=[2, TEST_BOX/2 - 1], even=.true., ctfparms=no_ctf)
        call workspace%prepare_particle_from_sigma_noise(observed, sigma2_noise, &
            &particle_spec, data)
        call assert_true(data%is_valid(), &
            &'adapter did not prepare a valid particle')
        call assert_true(all(data%get_shell_range() == [2, TEST_BOX/2 - 1]), &
            &'adapter changed the requested active shell range')

        initial_rotation = real(euler2m([20., 36.2, 28.7]), dp)
        initial_shift = [-0.08_dp, 0.06_dp]
        seed = pose_cont_pose(rotmat=initial_rotation, shift=initial_shift)
        config = pose_cont_config_from_route('shift_then_joint')
        call assert_true(config%route == POSE_CONT_ROUTE_SHIFT_THEN_JOINT, &
            &'adapter route constructor did not select shift-then-joint')
        limits = pose_cont_limits_from_boxes(2*TEST_BOX, TEST_BOX)
        call assert_true(abs(limits%shift_step_bound - 0.5_dp) < epsilon(1._dp) .and. &
            &abs(limits%max_total_shift - 2.5_dp) < epsilon(1._dp), &
            &'adapter limits constructor changed cropped-grid bounds')
        limits = pose_cont_limits_from_boxes(TEST_BOX, TEST_BOX)
        call assert_true(abs(limits%shift_step_bound - 1._dp) < epsilon(1._dp) .and. &
            &abs(limits%max_total_shift - 5._dp) < epsilon(1._dp), &
            &'adapter limits constructor changed native-grid bounds')
        call workspace%refine_particle(1, .true., seed, data, config, limits, result)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &result%objective_after < result%objective_before, &
            &'adapter transaction did not commit an improving pose')
        call assert_true(result%shift_stage%status >= LM_ACCEPTED_IMPROVEMENT .and. &
            &result%shift_stage%status <= LM_FINITE_NO_IMPROVEMENT .and. &
            &result%joint_stage%status >= LM_ACCEPTED_IMPROVEMENT .and. &
            &result%joint_stage%status <= LM_FINITE_NO_IMPROVEMENT, &
            &'adapter transaction did not complete both LM stages')
        call assert_true(result%attempts == result%shift_stage%attempts + result%joint_stage%attempts .and. &
            &result%accepts == result%shift_stage%accepts + result%joint_stage%accepts .and. &
            &result%bound_hits == result%shift_stage%bound_hits + result%joint_stage%bound_hits, &
            &'adapter stage accounting does not balance')

        config = pose_cont_config_from_route('joint')
        call workspace%refine_particle(1, .true., seed, data, config, limits, result)
        call assert_true(result%status == LM_ACCEPTED_IMPROVEMENT .and. &
            &result%objective_after < result%objective_before, &
            &'direct-joint adapter route did not commit an improving pose')
        call assert_true(result%shift_stage%status == POSE_CONT_NOT_ATTEMPTED .and. &
            &result%shift_stage%attempts == 0 .and. result%joint_stage%attempts > 0, &
            &'direct-joint adapter route executed or accounted for a shift-only stage')

        seed = pose_cont_pose(rotmat=truth_rotation, shift=initial_shift)
        config = pose_cont_config_from_route('shift_then_joint')
        limits%max_total_shift = 1.e-30_dp
        call workspace%refine_particle(1, .true., seed, data, config, limits, result)
        call assert_true(result%status == LM_STEP_BOUND_REJECTED .and. result%bound_hits > 0, &
            &'adapter did not report a cumulative shift-bound rejection')
        call assert_true(all(result%pose%rotmat == truth_rotation) .and. &
            &all(result%pose%shift == initial_shift), &
            &'bound-rejected adapter transaction did not preserve the input pose')

        seed = pose_cont_pose(rotmat=truth_rotation, shift=truth_shift)
        limits%max_total_shift = 5._dp
        call workspace%refine_particle(1, .true., seed, data, config, limits, result)
        write (*, '(a,3(1x,i0),4(1x,es12.4))') 'POSE_CONT_ADAPTER_EXACT', result%status, &
            &result%shift_stage%status, result%joint_stage%status, result%objective_before, &
            &result%objective_after, maxval(abs(result%pose%rotmat - truth_rotation)), &
            &maxval(abs(result%pose%shift - truth_shift))
        call assert_true(result%status == LM_FINITE_NO_IMPROVEMENT, &
            &'exact-pose adapter transaction did not report finite no-improvement')
        call assert_true(all(result%pose%rotmat == truth_rotation) .and. &
            &all(result%pose%shift == truth_shift), &
            &'non-improving adapter transaction changed the input pose')
        config = pose_cont_config_from_route('joint')
        call workspace%refine_particle(1, .true., seed, data, config, limits, result)
        call assert_true(result%status == LM_FINITE_NO_IMPROVEMENT .and. &
            &result%shift_stage%status == POSE_CONT_NOT_ATTEMPTED, &
            &'direct-joint exact-pose transaction did not report finite no-improvement')
        call assert_true(all(result%pose%rotmat == truth_rotation) .and. &
            &all(result%pose%shift == truth_shift), &
            &'direct-joint non-improving transaction changed the input pose')
        call workspace%refine_particle(1, .true., seed, invalid_data, config, limits, result)
        call assert_true(result%status == POSE_CONT_INVALID_PREPARATION .and. &
            &all(result%pose%rotmat == truth_rotation) .and. &
            &all(result%pose%shift == truth_shift), &
            &'invalid particle preparation changed the input pose')

        call generator%kill()
        call workspace%kill()
        call remove_pose_cont_reference_artifacts(1)
    end subroutine test_transaction_contracts

    ! Identity is a valid initialized pose; explicit state/half metadata, not
    ! nonzero Euler coordinates, defines readiness for the standalone class.
    subroutine test_strategy_seed_contract()
        type(ori) :: seed

        call seed%set_euler([0., 0., 0.])
        call seed%set_shift([1.25, -0.75])
        call seed%set('state', 1.)
        call seed%set('eo', 0.)
        call seed%set('proj', 1.)
        call assert_true(pose_cont_seed_is_valid(seed), &
            &'standalone pose strategy rejected a valid identity seed')
        call seed%set('proj', 0.)
        call assert_true(.not. pose_cont_seed_is_valid(seed), &
            &'standalone pose strategy accepted a missing projection seed')
        call seed%kill()
    end subroutine test_strategy_seed_contract

    subroutine write_reference(volume, even)
        real, intent(in) :: volume(:, :, :)
        logical, intent(in) :: even
        type(image) :: reference

        call reference%new([TEST_BOX, TEST_BOX, TEST_BOX], TEST_SMPD)
        call reference%set_rmat(volume, .false.)
        call write_pose_cont_reference_artifact(reference, 1, merge('even', 'odd ', even))
        call reference%kill()
    end subroutine write_reference

    subroutine build_test_volume(volume)
        real, allocatable, intent(out) :: volume(:, :, :)
        integer :: i, j, k

        allocate (volume(TEST_BOX, TEST_BOX, TEST_BOX))
        do k = 1, TEST_BOX
            do j = 1, TEST_BOX
                do i = 1, TEST_BOX
                    volume(i, j, k) = sin(0.11*real(i) + 0.17*real(j) - 0.07*real(k)) + &
                        &0.02*real(i*j - k) + 0.003*real(i*k)
                end do
            end do
        end do
    end subroutine build_test_volume

    subroutine assert_true(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) error stop trim(message)
    end subroutine assert_true

end program simple_test_pose_cont_refine3D_adapter
