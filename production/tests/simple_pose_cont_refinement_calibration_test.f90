module pose_cont_refinement_calibration_test
    use pose_cont_refinement_calibration_helpers
    use pose_cont_refinement_test_helpers, only: assert_true
    use simple_defs, only: dp, sp, OSMPL_PAD_FAC
    use simple_image, only: image
    use simple_cartesian_fourier, only: extract_native_fourier_plane
    use simple_cartesian_pose_refiner, only: cartesian_pose_refiner
    use simple_ori_utils, only: euler2m
    implicit none
    private
    public :: run_tolerance_calibration

    integer, parameter :: EXECUTED_DFT_FAMILY = 5
    integer, parameter :: SLOW_GATHER_FAMILY = 6
    integer, parameter :: FINITE_PROJECTION_FAMILY = 7
    integer, parameter :: N_FORWARD_ROTATIONS = 4
    integer, parameter :: N_CALIBRATION_FIXTURES = size(CALIBRATION_BOXES)*N_CALIBRATION_VARIANTS
    integer, parameter :: N_FORWARD_METRICS = 3*size(CALIBRATION_BOXES)* &
        &N_CALIBRATION_VARIANTS*N_FORWARD_ROTATIONS

    type :: forward_calibration_metric
        character(len=24) :: family = ''
        integer :: box = 0
        integer :: fixture_id = 0
        integer :: rotation = 0
        integer :: samples = 0
        real(dp) :: maximum_absolute_error = 0._dp
        real(dp) :: maximum_scaled_relative_error = 0._dp
        real(dp) :: clipping_fraction = 0._dp
        real(dp) :: original_volume_error = 0._dp
        real(dp) :: minimum_switch_margin = 0._dp
        integer :: stencil_switches = 0
    end type forward_calibration_metric

contains

!> Calibrate every tolerance family without constructing an acceptance fixture.
!! Calibration observes deterministic implementation-oracle disagreement twice,
!! applies the predeclared safety rule, and freezes only the approved proposals.
!! Holdout boxes and volumes are named here but are never sampled by this test.
    subroutine run_tolerance_calibration()
        type(calibration_fixture) :: first(N_CALIBRATION_FIXTURES)
        type(calibration_fixture) :: repeated(N_CALIBRATION_FIXTURES)
        type(tolerance_record) :: records(N_TOLERANCE_FAMILIES), repeated_records(N_TOLERANCE_FAMILIES)
        type(forward_calibration_metric) :: metrics(N_FORWARD_METRICS), repeated_metrics(N_FORWARD_METRICS)
        character(len=:), allocatable :: evidence_directory
        integer :: box_index, fixture_index, metric_count, repeated_metric_count, variant
        logical :: scientifically_loose, repeated_scientifically_loose

        evidence_directory = required_argument('evidence_dir')
        call initialize_tolerances(records)
        call initialize_tolerances(repeated_records)
        fixture_index = 0
        metric_count = 0
        repeated_metric_count = 0
        do box_index = 1, size(CALIBRATION_BOXES)
            do variant = 1, N_CALIBRATION_VARIANTS
                fixture_index = fixture_index + 1
                call build_calibration_fixture(CALIBRATION_BOXES(box_index), variant, first(fixture_index))
                call build_calibration_fixture(CALIBRATION_BOXES(box_index), variant, repeated(fixture_index))
                call assert_true(fixtures_are_identical(first(fixture_index), repeated(fixture_index)), &
                    &'calibration fixture is not deterministic')
                if (variant == 1) then
                    call collect_retained_observations(first(fixture_index), records)
                    call collect_retained_observations(repeated(fixture_index), repeated_records)
                end if
                call collect_forward_observations(first(fixture_index), records, metrics, metric_count)
                call collect_forward_observations(repeated(fixture_index), repeated_records, &
                    &repeated_metrics, repeated_metric_count)
            end do
        end do
        ! Derive proposals only after all calibration observations are collected.
        call finalize_tolerances(records, scientifically_loose)
        call finalize_tolerances(repeated_records, repeated_scientifically_loose)
        call require_repeatable_records(records, repeated_records)
        call require_calibration_coverage(records)
        call require_frozen_forward_tolerances(records)
        call require_repeatable_forward_metrics(metrics, metric_count, repeated_metrics, repeated_metric_count)
        call check_combined_metric()
        call assert_true(.not. any(CALIBRATION_BOXES(1) == ACCEPTANCE_BOXES) .and. &
            &.not. any(CALIBRATION_BOXES(2) == ACCEPTANCE_BOXES), &
            &'calibration and acceptance box sets overlap')
        call write_calibration_artifacts(evidence_directory, records, first)
        call write_forward_calibration_artifacts(evidence_directory, metrics(:metric_count))
        call write_holdout_static_review(evidence_directory)
        call assert_true(scientifically_loose .eqv. repeated_scientifically_loose, &
            &'repeated calibration changed the scientific-review outcome')
        call assert_true(.not. scientifically_loose, &
            &'forward calibration produced a scientifically loose proposal')

        write (*, '(a)') 'POSE_CONT_REFINEMENT_CALIBRATION boxes/variants: 8 12 / 2'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CALIBRATION full-disk rotations: 4'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CALIBRATION observed diagnostics not sampled: 10 16'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CALIBRATION frozen holdouts not sampled: 14 18'
        write (*, '(a,es14.6)') 'POSE_CONT_REFINEMENT_CALIBRATION safety factor: ', CALIBRATION_SAFETY_FACTOR
        write (*, '(a)') 'POSE_CONT_REFINEMENT_FORWARD_AMENDMENT: FROZEN'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_TOLERANCE_CALIBRATION: PASS'
    end subroutine run_tolerance_calibration

!> Verify the componentwise absolute-plus-relative comparison boundary.
    subroutine check_combined_metric()
        real(dp), parameter :: absolute_tolerance = 1.e-4_dp, relative_tolerance = 2.e-3_dp

        call assert_true(combined_real_passes([1._dp, 0._dp], [1.001_dp, 5.e-5_dp], &
            &absolute_tolerance, relative_tolerance), &
            'combined real metric rejected an in-bound vector')
        call assert_true(.not. combined_real_passes([1._dp], [1.01_dp], absolute_tolerance, relative_tolerance), &
            &'combined real metric accepted an out-of-bound component')
        call assert_true(combined_complex_passes([cmplx(1._dp, 1._dp, kind=dp)], &
            &[cmplx(1.001_dp, 1._dp, kind=dp)], absolute_tolerance, relative_tolerance), &
            &'combined complex metric rejected an in-bound value')
    end subroutine check_combined_metric

!> Reproduce the four tolerance families that remain frozen.
    subroutine collect_retained_observations(fixture, records)
        type(calibration_fixture), intent(in) :: fixture
        type(tolerance_record), intent(inout) :: records(N_TOLERANCE_FAMILIES)
        complex(dp) :: analytic_values(4), brute_values(4)
        real(dp) :: exact_values(5), rounded_real(5), derivative_exact(5), derivative_fd(5)
        real(dp) :: lm_matrix(5, 5), lm_rhs(5), lm_solution(5), lm_residual(5)
        real(dp) :: volume_sum_dp, volume_sum_sp
        integer :: axis, mode

        volume_sum_dp = sum(fixture%volume)
        volume_sum_sp = real(sum(real(fixture%volume, sp)), dp)
        call observe_real_comparison(records(1), [volume_sum_sp], [volume_sum_dp])

        do axis = 1, 5
            exact_values(axis) = fixture%nonstationary_pose(axis) + 0.03_dp*real(axis, dp)
            rounded_real(axis) = real(real(exact_values(axis), sp), dp)
        end do
        lm_matrix = 0._dp
        do axis = 1, 5
            lm_matrix(axis, axis) = 1._dp + 0.2_dp*real(axis, dp)
        end do
        lm_rhs = matmul(lm_matrix, exact_values)
        lm_solution = rounded_real
        lm_residual = matmul(lm_matrix, lm_solution) - lm_rhs
        call observe_real_comparison(records(2), lm_solution, exact_values)
        call observe_real_comparison(records(2), lm_residual, [0._dp, 0._dp, 0._dp, 0._dp, 0._dp])

        do axis = 1, 5
            derivative_exact(axis) = cos(fixture%nonstationary_pose(axis)) + &
                &0.2_dp*fixture%nonstationary_pose(axis)
            derivative_fd(axis) = centred_probe_derivative(fixture%nonstationary_pose(axis))
        end do
        call observe_real_comparison(records(3), derivative_fd, derivative_exact)

        do mode = 1, 4
            analytic_values(mode) = analytic_point_dft(fixture%box, mode)
            brute_values(mode) = brute_point_dft(fixture%box, mode)
        end do
        call observe_complex_comparison(records(4), brute_values, analytic_values)
    end subroutine collect_retained_observations

!> Calibrate the three reopened forward families over the full redundant disk.
    subroutine collect_forward_observations(fixture, records, metrics, metric_count)
        type(calibration_fixture), intent(in) :: fixture
        type(tolerance_record), intent(inout) :: records(N_TOLERANCE_FAMILIES)
        type(forward_calibration_metric), intent(inout) :: metrics(:)
        integer, intent(inout) :: metric_count
        type(cartesian_pose_refiner) :: workspace
        complex(dp), allocatable :: direct(:), executed(:), fast(:), slow(:)
        complex(dp), allocatable :: projected(:), matched_direct(:), original_direct(:)
        complex, allocatable :: zero_plane(:, :), executed_plane(:, :), projection_plane(:, :)
        real, allocatable :: physical_volume(:, :, :)
        real(dp), allocatable :: rotated_volume(:, :, :)
        real(dp) :: rotations(3, 3, N_FORWARD_ROTATIONS), ignored_objective
        real(dp) :: clipping_fraction, minimum_switch_margin
        integer :: active_samples, lims2(2, 2), rotation_index, stencil_switches

        call build_forward_rotations(rotations)
        allocate (physical_volume, source=real(fixture%volume, sp))
        call workspace%new(physical_volume)
        call workspace%set_shell_range([0, fixture%box/2])
        lims2 = workspace%get_lims2()
        allocate (zero_plane(lims2(1, 1):lims2(1, 2), lims2(2, 1):lims2(2, 2)), source=cmplx(0., 0.))
        allocate (executed_plane, mold=zero_plane)
        active_samples = full_disk_sample_count(fixture%box)
        allocate (direct(active_samples), executed(active_samples), fast(active_samples), slow(active_samples))
        allocate (projected(active_samples), matched_direct(active_samples), original_direct(active_samples))
        allocate (rotated_volume(fixture%box, fixture%box, fixture%box))

        do rotation_index = 1, N_FORWARD_ROTATIONS
            call workspace%shift_residual(rotations(:, :, rotation_index), [0._dp, 0._dp], &
                &zero_plane, executed_plane, ignored_objective)
            call collect_full_disk_planes(fixture%volume, workspace, rotations(:, :, rotation_index), &
                &executed_plane, direct, executed, fast, slow, minimum_switch_margin)
            ! Contract 1: executed gather G(E^-1 X) versus the direct DFT of X.
            call observe_complex_comparison(records(EXECUTED_DFT_FAMILY), executed, direct)
            call append_forward_metric(records(EXECUTED_DFT_FAMILY), fixture, rotation_index, &
                &'executed_dft', executed, direct, minimum_switch_margin, metrics, metric_count)

            ! Contract 2: fast packed gather versus a structurally slower traversal.
            call observe_complex_comparison(records(SLOW_GATHER_FAMILY), fast, slow)
            stencil_switches = workspace%count_stencil_switches(rotations(:, :, 1), &
                &rotations(:, :, rotation_index), [0, fixture%box/2])
            call append_forward_metric(records(SLOW_GATHER_FAMILY), fixture, rotation_index, &
                &'slow_gather', fast, slow, minimum_switch_margin, metrics, metric_count, &
                &stencil_switches=stencil_switches)

            call rotate_finite_volume(fixture%volume, rotations(:, :, rotation_index), &
                &rotated_volume, clipping_fraction)
            call finite_projection_fft(rotated_volume, lims2, projection_plane)
            call flatten_full_disk(projection_plane, fixture%box, projected)
            call collect_direct_planes(fixture%volume, rotated_volume, rotations(:, :, rotation_index), &
                &original_direct, matched_direct)
            ! Contract 3: projection-FFT versus the direct DFT of the same finite X_R.
            call observe_complex_comparison(records(FINITE_PROJECTION_FAMILY), projected, matched_direct)
            call append_forward_metric(records(FINITE_PROJECTION_FAMILY), fixture, rotation_index, &
                &'finite_projection', projected, matched_direct, minimum_switch_margin, metrics, &
                &metric_count, clipping_fraction=clipping_fraction, &
                &original_volume_error=maxval(abs(matched_direct - original_direct)))
            deallocate (projection_plane)
        end do

        deallocate (physical_volume, zero_plane, executed_plane, direct, executed, fast, slow)
        deallocate (projected, matched_direct, original_direct, rotated_volume)
        call workspace%kill
    end subroutine collect_forward_observations

!> Collect direct, executed, fast, and slow samples in one deterministic order.
    subroutine collect_full_disk_planes(volume, workspace, rotation, executed_plane, direct, &
        &executed, fast, slow, minimum_switch_margin)
        real(dp), intent(in) :: volume(:, :, :), rotation(3, 3)
        type(cartesian_pose_refiner), intent(in) :: workspace
        complex, intent(in) :: executed_plane(-size(volume, 1)/2:, -size(volume, 1)/2:)
        complex(dp), intent(out) :: direct(:), executed(:), fast(:), slow(:)
        real(dp), intent(out) :: minimum_switch_margin
        complex :: fast_value, ignored_gradient(3)
        real :: margin(3)
        real(sp) :: location(3)
        integer :: box, h, k, sample

        box = size(volume, 1)
        sample = 0
        minimum_switch_margin = huge(0._dp)
        do k = -box/2, box/2
            do h = -box/2, box/2
                if (h*h + k*k > (box/2)**2) cycle
                sample = sample + 1
                location = real(OSMPL_PAD_FAC, sp)* &
                    &real(matmul(real([h, k, 0], dp), rotation), sp)
                call workspace%sample_with_grad(location, fast_value, ignored_gradient, margin)
                call workspace%sample_slow_test(location, slow(sample))
                direct(sample) = direct_volume_dft(volume, rotation, h, k)
                executed(sample) = cmplx(executed_plane(h, k), kind=dp)
                fast(sample) = cmplx(fast_value, kind=dp)
                minimum_switch_margin = min(minimum_switch_margin, real(minval(margin), dp))
            end do
        end do
        call assert_true(sample == size(direct), 'forward calibration omitted full-disk samples')
    end subroutine collect_full_disk_planes

!> Evaluate original-frequency and matched finite-object DFT planes.
    subroutine collect_direct_planes(original, rotated, rotation, original_direct, matched_direct)
        real(dp), intent(in) :: original(:, :, :), rotated(:, :, :), rotation(3, 3)
        complex(dp), intent(out) :: original_direct(:), matched_direct(:)
        real(dp) :: identity(3, 3)
        integer :: box, h, k, sample

        identity = 0._dp
        identity(1, 1) = 1._dp
        identity(2, 2) = 1._dp
        identity(3, 3) = 1._dp
        box = size(original, 1)
        sample = 0
        do k = -box/2, box/2
            do h = -box/2, box/2
                if (h*h + k*k > (box/2)**2) cycle
                sample = sample + 1
                original_direct(sample) = direct_volume_dft(original, rotation, h, k)
                matched_direct(sample) = direct_volume_dft(rotated, identity, h, k)
            end do
        end do
        call assert_true(sample == size(original_direct), 'direct calibration plane is incomplete')
    end subroutine collect_direct_planes

!> Evaluate the normalized direct three-dimensional DFT at one rotated frequency.
    function direct_volume_dft(volume, rotation, h, k) result(value)
        real(dp), intent(in) :: volume(:, :, :), rotation(3, 3)
        integer, intent(in) :: h, k
        complex(dp) :: value
        real(dp) :: angle, location(3), twopi
        integer :: box, i, j, l, x, y, z

        box = size(volume, 1)
        twopi = 2._dp*acos(-1._dp)
        location = matmul(real([h, k, 0], dp), rotation)
        value = cmplx(0._dp, 0._dp, kind=dp)
        do l = 1, box
            z = l - (box/2 + 1)
            do j = 1, box
                y = j - (box/2 + 1)
                do i = 1, box
                    x = i - (box/2 + 1)
                    angle = -twopi*dot_product(location, real([x, y, z], dp))/real(box, dp)
                    value = value + volume(i, j, l)*cmplx(cos(angle), sin(angle), kind=dp)
                end do
            end do
        end do
        value = value/real(box**3, dp)
    end function direct_volume_dft

!> Resample the physical volume once so both finite-projection paths use X_R.
    subroutine rotate_finite_volume(volume, rotation, rotated, clipping_fraction)
        real(dp), intent(in) :: volume(:, :, :), rotation(3, 3)
        real(dp), intent(out) :: rotated(:, :, :), clipping_fraction
        real(dp) :: coordinate(3)
        integer :: box, clipped, i, j, k, u, v, w

        box = size(volume, 1)
        clipped = 0
        do k = 1, box
            w = k - (box/2 + 1)
            do j = 1, box
                v = j - (box/2 + 1)
                do i = 1, box
                    u = i - (box/2 + 1)
                    coordinate = matmul(real([u, v, w], dp), rotation)
                    call trilinear_finite_sample(volume, coordinate, rotated(i, j, k), clipped)
                end do
            end do
        end do
        clipping_fraction = real(clipped, dp)/real(box**3, dp)
    end subroutine rotate_finite_volume

!> Sample a centered finite box with trilinear interpolation and zero outside.
    subroutine trilinear_finite_sample(volume, coordinate, value, clipped)
        real(dp), intent(in) :: volume(:, :, :), coordinate(3)
        real(dp), intent(out) :: value
        integer, intent(inout) :: clipped
        real(dp) :: fraction(3), weight
        integer :: box, lower(3), upper(3), ix, iy, iz, xindex, yindex, zindex

        box = size(volume, 1)
        lower = floor(coordinate)
        fraction = coordinate - real(lower, dp)
        upper = lower + 1
        where (abs(fraction) <= 8._dp*epsilon(1._dp)) upper = lower
        if (any(lower < -box/2) .or. any(upper > box/2 - 1)) then
            value = 0._dp
            clipped = clipped + 1
            return
        end if
        value = 0._dp
        do iz = 0, 1
            zindex = merge(lower(3), upper(3), iz == 0) + box/2 + 1
            do iy = 0, 1
                yindex = merge(lower(2), upper(2), iy == 0) + box/2 + 1
                do ix = 0, 1
                    xindex = merge(lower(1), upper(1), ix == 0) + box/2 + 1
                    weight = merge(1._dp - fraction(1), fraction(1), ix == 0)* &
                        &merge(1._dp - fraction(2), fraction(2), iy == 0)* &
                        &merge(1._dp - fraction(3), fraction(3), iz == 0)
                    value = value + weight*volume(xindex, yindex, zindex)
                end do
            end do
        end do
    end subroutine trilinear_finite_sample

!> Form a line-sum projection and the executed normalized two-dimensional FFT.
    subroutine finite_projection_fft(volume, lims2, plane)
        real(dp), intent(in) :: volume(:, :, :)
        integer, intent(in) :: lims2(2, 2)
        complex, allocatable, intent(out) :: plane(:, :)
        type(image) :: projection_image
        real :: projection(size(volume, 1), size(volume, 1), 1)
        integer :: box, i, j

        box = size(volume, 1)
        do j = 1, box
            do i = 1, box
                projection(i, j, 1) = real(sum(volume(i, j, :)), sp)
            end do
        end do
        call projection_image%new([box, box, 1], 1.0)
        call projection_image%set_rmat(projection, .false.)
        call projection_image%fft()
        plane = extract_native_fourier_plane(projection_image, lims2, (box/2)**2)/real(box, sp)
        call projection_image%kill
    end subroutine finite_projection_fft

!> Flatten the full redundant disk without discarding Nyquist mates.
    subroutine flatten_full_disk(plane, box, values)
        integer, intent(in) :: box
        complex, intent(in) :: plane(-box/2:, -box/2:)
        complex(dp), intent(out) :: values(:)
        integer :: h, k, sample

        sample = 0
        do k = -box/2, box/2
            do h = -box/2, box/2
                if (h*h + k*k > (box/2)**2) cycle
                sample = sample + 1
                values(sample) = cmplx(plane(h, k), kind=dp)
            end do
        end do
        call assert_true(sample == size(values), 'projection calibration omitted full-disk samples')
    end subroutine flatten_full_disk

!> Append one aggregate row while raw component errors remain in the family record.
    subroutine append_forward_metric(record, fixture, rotation_index, family, actual, expected, &
        &minimum_switch_margin, metrics, metric_count, clipping_fraction, original_volume_error, &
        &stencil_switches)
        type(tolerance_record), intent(in) :: record
        type(calibration_fixture), intent(in) :: fixture
        integer, intent(in) :: rotation_index
        character(len=*), intent(in) :: family
        complex(dp), intent(in) :: actual(:), expected(:)
        real(dp), intent(in) :: minimum_switch_margin
        type(forward_calibration_metric), intent(inout) :: metrics(:)
        integer, intent(inout) :: metric_count
        real(dp), optional, intent(in) :: clipping_fraction, original_volume_error
        integer, optional, intent(in) :: stencil_switches
        real(dp) :: absolute_error, scale
        integer :: sample

        metric_count = metric_count + 1
        if (metric_count > size(metrics)) error stop 'forward calibration metric capacity exceeded'
        metrics(metric_count)%family = family
        metrics(metric_count)%box = fixture%box
        metrics(metric_count)%fixture_id = fixture%fixture_id
        metrics(metric_count)%rotation = rotation_index
        metrics(metric_count)%samples = size(actual)
        metrics(metric_count)%minimum_switch_margin = minimum_switch_margin
        if (present(clipping_fraction)) metrics(metric_count)%clipping_fraction = clipping_fraction
        if (present(original_volume_error)) &
            &metrics(metric_count)%original_volume_error = original_volume_error
        if (present(stencil_switches)) metrics(metric_count)%stencil_switches = stencil_switches
        do sample = 1, size(actual)
            absolute_error = abs(actual(sample) - expected(sample))
            scale = max(record%relative_scale_floor, abs(actual(sample)), abs(expected(sample)))
            metrics(metric_count)%maximum_absolute_error = max( &
                &metrics(metric_count)%maximum_absolute_error, absolute_error)
            metrics(metric_count)%maximum_scaled_relative_error = max( &
                &metrics(metric_count)%maximum_scaled_relative_error, absolute_error/scale)
        end do
    end subroutine append_forward_metric

!> Predeclare identity and three arbitrary rotations for amended calibration.
    subroutine build_forward_rotations(rotations)
        real(dp), intent(out) :: rotations(3, 3, N_FORWARD_ROTATIONS)

        rotations = 0._dp
        rotations(1, 1, 1) = 1._dp
        rotations(2, 2, 1) = 1._dp
        rotations(3, 3, 1) = 1._dp
        rotations(:, :, 2) = real(euler2m([17._sp, 31._sp, -23._sp]), dp)
        rotations(:, :, 3) = real(euler2m([-29._sp, 13._sp, 41._sp]), dp)
        rotations(:, :, 4) = real(euler2m([38._sp, -27._sp, 9._sp]), dp)
    end subroutine build_forward_rotations

    integer function full_disk_sample_count(box) result(count)
        integer, intent(in) :: box
        integer :: h, k

        count = 0
        do k = -box/2, box/2
            do h = -box/2, box/2
                if (h*h + k*k <= (box/2)**2) count = count + 1
            end do
        end do
    end function full_disk_sample_count

!> Return the finite-difference derivative of the independent scalar probe.
    pure real(dp) function centred_probe_derivative(value) result(derivative)
        real(dp), intent(in) :: value
        real(dp), parameter :: step = 2.e-5_dp
        derivative = (probe_value(value + step) - probe_value(value - step))/(2._dp*step)
    end function centred_probe_derivative

!> Return a smooth independent scalar probe for derivative calibration.
    pure real(dp) function probe_value(value) result(probe)
        real(dp), intent(in) :: value
        probe = sin(value) + 0.1_dp*value*value
    end function probe_value

!> Evaluate the closed-form DFT of a non-collinear unequal-amplitude point set.
    function analytic_point_dft(box, mode) result(value)
        integer, intent(in) :: box, mode
        complex(dp) :: value
        integer, parameter :: points(3, 4) = reshape([0, 0, 0, 1, -2, 1, -2, 1, -1, 2, 2, -2], [3, 4])
        integer, parameter :: h_modes(4) = [0, 1, -2, 1]
        integer, parameter :: k_modes(4) = [0, -1, 1, 2]
        integer, parameter :: l_modes(4) = [0, 1, 1, -1]
        real(dp), parameter :: amplitudes(4) = [1._dp, 0.61_dp, -0.37_dp, 0.23_dp]
        real(dp) :: angle, twopi
        integer :: point

        twopi = 2._dp*acos(-1._dp)
        value = cmplx(0._dp, 0._dp, kind=dp)
        do point = 1, 4
            angle = -twopi*real(h_modes(mode)*points(1, point) + k_modes(mode)*points(2, point) + &
                &l_modes(mode)*points(3, point), dp)/real(box, dp)
            value = value + amplitudes(point)*cmplx(cos(angle), sin(angle), kind=dp)
        end do
        value = value/real(box**3, dp)
    end function analytic_point_dft

!> Evaluate the same point set through a full voxel-loop brute-force DFT.
    function brute_point_dft(box, mode) result(value)
        integer, intent(in) :: box, mode
        complex(dp) :: value
        integer, parameter :: points(3, 4) = reshape([0, 0, 0, 1, -2, 1, -2, 1, -1, 2, 2, -2], [3, 4])
        integer, parameter :: h_modes(4) = [0, 1, -2, 1]
        integer, parameter :: k_modes(4) = [0, -1, 1, 2]
        integer, parameter :: l_modes(4) = [0, 1, 1, -1]
        real(dp), parameter :: amplitudes(4) = [1._dp, 0.61_dp, -0.37_dp, 0.23_dp]
        real(dp) :: angle, sample, twopi
        integer :: h, k, l, point, x, y, z

        h = h_modes(mode)
        k = k_modes(mode)
        l = l_modes(mode)
        twopi = 2._dp*acos(-1._dp)
        value = cmplx(0._dp, 0._dp, kind=dp)
        do z = -box/2, box/2 - 1
            do y = -box/2, box/2 - 1
                do x = -box/2, box/2 - 1
                    sample = 0._dp
                    do point = 1, 4
                        if (all([x, y, z] == points(:, point))) sample = sample + amplitudes(point)
                    end do
                    angle = -twopi*real(h*x + k*y + l*z, dp)/real(box, dp)
                    value = value + sample*cmplx(cos(angle), sin(angle), kind=dp)
                end do
            end do
        end do
        value = value/real(box**3, dp)
    end function brute_point_dft

!> Write aggregate forward errors and diagnostics without freezing proposals.
    subroutine write_forward_calibration_artifacts(output_directory, metrics)
        character(len=*), intent(in) :: output_directory
        type(forward_calibration_metric), intent(in) :: metrics(:)
        integer :: item, unit

        open (newunit=unit, file=trim(output_directory)//'/forward_calibration.tsv', &
            &status='replace', action='write')
        write (unit, '(a)') 'family'//achar(9)//'box'//achar(9)//'fixture_id'//achar(9)// &
            &'rotation'//achar(9)//'samples'//achar(9)//'max_absolute'//achar(9)// &
            &'max_scaled_relative'//achar(9)//'clipping_fraction'//achar(9)// &
            &'original_volume_error'//achar(9)//'minimum_switch_margin'//achar(9)// &
            &'stencil_switches'
        do item = 1, size(metrics)
            write (unit, '(a,4(a,i0),5(a,es24.16),a,i0)') trim(metrics(item)%family), achar(9), &
                &metrics(item)%box, achar(9), metrics(item)%fixture_id, achar(9), &
                &metrics(item)%rotation, achar(9), metrics(item)%samples, achar(9), &
                &metrics(item)%maximum_absolute_error, achar(9), &
                &metrics(item)%maximum_scaled_relative_error, achar(9), &
                &metrics(item)%clipping_fraction, achar(9), metrics(item)%original_volume_error, &
                &achar(9), metrics(item)%minimum_switch_margin, achar(9), metrics(item)%stencil_switches
        end do
        close (unit)
    end subroutine write_forward_calibration_artifacts

!> Record frozen holdout cost without constructing or sampling a holdout fixture.
    subroutine write_holdout_static_review(output_directory)
        character(len=*), intent(in) :: output_directory
        integer :: box, box_index, disk_samples, unit

        open (newunit=unit, file=trim(output_directory)//'/forward_holdout_static_review.tsv', &
            &status='replace', action='write')
        write (unit, '(a)') 'box'//achar(9)//'full_disk_samples'//achar(9)//'voxels'//achar(9)// &
            &'direct_dft_terms_per_plane'//achar(9)//'status'
        do box_index = 1, size(FORWARD_HOLDOUT_BOXES)
            box = FORWARD_HOLDOUT_BOXES(box_index)
            disk_samples = full_disk_sample_count(box)
            write (unit, '(i0,a,i0,a,i0,a,i0,a,a)') box, achar(9), disk_samples, achar(9), &
                &box**3, achar(9), disk_samples*box**3, achar(9), 'frozen_holdout_not_sampled'
        end do
        close (unit)
    end subroutine write_holdout_static_review

!> Require deterministic aggregate forward metrics from the repeated pass.
    subroutine require_repeatable_forward_metrics(first, nfirst, repeated, nrepeated)
        type(forward_calibration_metric), intent(in) :: first(:), repeated(:)
        integer, intent(in) :: nfirst, nrepeated
        integer :: item

        call assert_true(nfirst == nrepeated .and. nfirst == N_FORWARD_METRICS, &
            &'repeated forward calibration changed its metric count')
        do item = 1, nfirst
            call assert_true(first(item)%family == repeated(item)%family .and. &
                &first(item)%box == repeated(item)%box .and. &
                &first(item)%fixture_id == repeated(item)%fixture_id .and. &
                &first(item)%rotation == repeated(item)%rotation .and. &
                &first(item)%samples == repeated(item)%samples, &
                &'repeated forward calibration changed metric identity')
            call assert_true(first(item)%maximum_absolute_error == &
                &repeated(item)%maximum_absolute_error .and. &
                &first(item)%maximum_scaled_relative_error == &
                &repeated(item)%maximum_scaled_relative_error .and. &
                &first(item)%clipping_fraction == repeated(item)%clipping_fraction .and. &
                &first(item)%original_volume_error == repeated(item)%original_volume_error .and. &
                &first(item)%minimum_switch_margin == repeated(item)%minimum_switch_margin .and. &
                &first(item)%stencil_switches == repeated(item)%stencil_switches, &
                &'repeated forward calibration changed metric values')
        end do
    end subroutine require_repeatable_forward_metrics

!> Verify that each proposal covers its observations through the combined rule.
    subroutine require_calibration_coverage(records)
        type(tolerance_record), intent(in) :: records(N_TOLERANCE_FAMILIES)
        real(dp) :: combined_limit
        integer :: family, observation

        do family = 1, N_TOLERANCE_FAMILIES
            if (family >= EXECUTED_DFT_FAMILY) then
                call assert_true(records(family)%relative_tolerance == records(family)%relative_floor, &
                    &'forward proposal changed its predeclared relative floor')
            end if
            do observation = 1, records(family)%observations
                combined_limit = records(family)%absolute_tolerance + records(family)%relative_tolerance* &
                    &records(family)%comparison_scales(observation)
                call assert_true(records(family)%absolute_errors(observation) <= combined_limit, &
                    &'combined tolerance proposal does not cover a calibration observation')
            end do
        end do
    end subroutine require_calibration_coverage

!> Require the reviewed forward derivation to reproduce the frozen constants exactly.
    subroutine require_frozen_forward_tolerances(records)
        type(tolerance_record), intent(in) :: records(N_TOLERANCE_FAMILIES)
        integer :: family

        do family = EXECUTED_DFT_FAMILY, N_TOLERANCE_FAMILIES
            call assert_true(records(family)%absolute_tolerance == &
                &FROZEN_ABSOLUTE_TOLERANCES(family), &
                &'derived forward absolute tolerance differs from its frozen value')
            call assert_true(records(family)%relative_tolerance == &
                &FROZEN_RELATIVE_TOLERANCES(family), &
                &'derived forward relative tolerance differs from its frozen value')
        end do
    end subroutine require_frozen_forward_tolerances

!> Require bitwise-identical calibration records from the repeated run.
    subroutine require_repeatable_records(first, repeated)
        type(tolerance_record), intent(in) :: first(N_TOLERANCE_FAMILIES)
        type(tolerance_record), intent(in) :: repeated(N_TOLERANCE_FAMILIES)
        integer :: family

        do family = 1, N_TOLERANCE_FAMILIES
            call assert_true(first(family)%name == repeated(family)%name, &
                &'repeated calibration changed a family name')
            call assert_true(first(family)%observations == repeated(family)%observations, &
                &'repeated calibration changed an observation count')
            call assert_true(first(family)%maximum_absolute_error == &
                &repeated(family)%maximum_absolute_error, &
                'repeated calibration changed an absolute maximum')
            call assert_true(first(family)%maximum_scaled_relative_error == &
                &repeated(family)%maximum_scaled_relative_error, &
                &'repeated calibration changed a relative maximum')
            call assert_true(all(first(family)%absolute_errors == repeated(family)%absolute_errors) .and. &
                &all(first(family)%scaled_relative_errors == repeated(family)%scaled_relative_errors), &
                &'repeated calibration changed a raw error table')
            call assert_true(all(first(family)%comparison_scales == repeated(family)%comparison_scales) .and. &
                &first(family)%maximum_combined_excess == repeated(family)%maximum_combined_excess, &
                &'repeated calibration changed a combined comparison scale')
            call assert_true(first(family)%absolute_tolerance == repeated(family)%absolute_tolerance .and. &
                &first(family)%relative_tolerance == repeated(family)%relative_tolerance, &
                &'repeated calibration changed a derived tolerance')
        end do
    end subroutine require_repeatable_records

!> Read one required key=value command argument.
    function required_argument(key) result(value)
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: value
        character(len=4096) :: argument
        integer :: argument_index, argument_status, separator

        value = ''
        do argument_index = 1, command_argument_count()
            call get_command_argument(argument_index, argument, status=argument_status)
            if (argument_status /= 0) error stop 'could not read calibration argument'
            separator = index(argument, '=')
            if (separator <= 0) cycle
            if (trim(argument(:separator - 1)) /= key) cycle
            value = trim(argument(separator + 1:))
        end do
        if (len(value) == 0) error stop 'tolerance calibration requires evidence_dir=<existing-directory>'
    end function required_argument

end module pose_cont_refinement_calibration_test
