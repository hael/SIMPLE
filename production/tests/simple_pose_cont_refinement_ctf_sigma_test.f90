! Independent CTF, whitening, and variance-range evidence for Cartesian pose data.
! Raw CTF-distorted observations are formed before whitening so the oracle does
! not reproduce the pose owner's operation ordering by construction.
module pose_cont_refinement_ctf_sigma_test
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use pose_cont_refinement_test_helpers, only: assert_true
    use simple_defs, only: dp, sp
    use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
        &POSE_DATA_VALID, POSE_DATA_INVALID_NOISE_RANGE
    use simple_reconstructor_pcg, only: reconstructor_pcg
    use simple_ctf, only: ctf
    use simple_type_defs, only: ctfparams, CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP
    implicit none
    private
    public :: run_ctf_sigma

    integer, parameter :: BOX = 10
    real(dp), parameter :: COMPONENT_TOLERANCE = 4.e-6_dp

contains

!> Verify CTF mode, whitening order, variance limits, and distinct PCG fallback.
!! The pose transfer is T_pose=C/sqrt(sigma2) and contains no shift phase.
!! Short variance caps the pose range, while the existing PCG fallback is retained.
    subroutine run_ctf_sigma()
        type(cartesian_pose_refiner) :: workspace
        type(reconstructor_pcg) :: pcgop
        type(ctfparams) :: ordinary_ctf, zero_ctf, no_ctf
        real, allocatable :: volume(:, :, :), constant_sigma(:), varying_sigma(:), short_sigma(:)
        real, allocatable :: invalid_sigma(:), unused_invalid_sigma(:)
        character(len=:), allocatable :: evidence_directory
        integer :: case_unit, fallback_unit, outcome_unit, shell

        evidence_directory = required_argument('evidence_dir')
        open (newunit=case_unit, file=evidence_directory//'/ctf_sigma_cases.tsv', &
            &status='replace', action='write')
        open (newunit=outcome_unit, file=evidence_directory//'/variance_outcomes.tsv', &
            &status='replace', action='write')
        open (newunit=fallback_unit, file=evidence_directory//'/pcg_fallback.tsv', &
            &status='replace', action='write')
        write (case_unit, '(a)') 'case'//achar(9)//'requested_low'//achar(9)//'requested_high'// &
            &achar(9)//'effective_low'//achar(9)//'effective_high'//achar(9)//'active_samples'// &
            &achar(9)//'max_transfer_error'//achar(9)//'max_observed_error'
        write (outcome_unit, '(a)') 'case'//achar(9)//'status'//achar(9)//'valid'//achar(9)// &
            &'requested_low'//achar(9)//'requested_high'//achar(9)//'effective_low'//achar(9)// &
            &'effective_high'//achar(9)//'active_samples'
        write (fallback_unit, '(a)') 'pose_effective_high'//achar(9)//'pose_tail_abs'//achar(9)// &
            &'pcg_tail_abs'//achar(9)//'pcg_extended_max_error'

        allocate (volume(BOX, BOX, BOX), source=0.0)
        call workspace%new(volume)
        call pcgop%new(BOX, 1._sp)
        call make_ctf_parameters(ordinary_ctf, zero_ctf, no_ctf)
        allocate (constant_sigma(0:BOX/2), source=2.25)
        allocate (varying_sigma(0:BOX/2))
        varying_sigma = [(1.0 + 0.17*real(shell), shell=0, BOX/2)]
        allocate (short_sigma(0:BOX/2 - 2))
        short_sigma = [(1.2 + 0.11*real(shell), shell=0, BOX/2 - 2)]

        call test_valid_case(workspace, 'ctf_no_constant', no_ctf, constant_sigma, [0, BOX/2], &
            &case_unit, .false.)
        call test_valid_case(workspace, 'ctf_signed_constant', ordinary_ctf, constant_sigma, &
            &[0, BOX/2], case_unit, .false.)
        call test_valid_case(workspace, 'ctf_signed_varying', ordinary_ctf, varying_sigma, &
            &[1, BOX/2], case_unit, .false.)
        ordinary_ctf%ctfflag = CTFFLAG_FLIP
        call test_valid_case(workspace, 'ctf_flip_varying', ordinary_ctf, varying_sigma, &
            &[1, BOX/2], case_unit, .false.)
        call test_valid_case(workspace, 'physical_ctf_zero', zero_ctf, constant_sigma, &
            &[0, BOX/2], case_unit, .true.)
        ordinary_ctf%ctfflag = CTFFLAG_YES
        call test_valid_case(workspace, 'short_variance', ordinary_ctf, short_sigma, &
            &[1, BOX/2], case_unit, .false.)

        allocate (invalid_sigma, source=varying_sigma)
        invalid_sigma(2) = -1.0
        call test_invalid_case(workspace, 'invalid_used_negative', ordinary_ctf, invalid_sigma, &
            &[1, BOX/2], [1, BOX/2], outcome_unit)
        invalid_sigma = varying_sigma
        invalid_sigma(3) = ieee_value(0.0, ieee_quiet_nan)
        call test_invalid_case(workspace, 'invalid_used_nan', ordinary_ctf, invalid_sigma, &
            &[1, BOX/2], [1, BOX/2], outcome_unit)
        call test_invalid_case(workspace, 'no_usable_overlap', ordinary_ctf, varying_sigma, &
            &[BOX/2 + 1, BOX/2 + 2], [BOX/2 + 1, BOX/2], outcome_unit)

        allocate (unused_invalid_sigma, source=varying_sigma)
        unused_invalid_sigma(BOX/2) = -1.0
        call test_valid_case(workspace, 'invalid_unused_shell', ordinary_ctf, unused_invalid_sigma, &
            &[1, 2], case_unit, .false.)
        call test_pcg_last_shell_fallback(workspace, pcgop, no_ctf, short_sigma, fallback_unit)

        close (case_unit)
        close (outcome_unit)
        close (fallback_unit)
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CTF_SIGMA modes: no signed flip physical_zero'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CTF_SIGMA variance: constant varying short invalid no_overlap'
        write (*, '(a)') 'POSE_CONT_REFINEMENT_CTF_SIGMA: PASS'
        call workspace%kill
        call pcgop%kill
    end subroutine run_ctf_sigma

!> Define ordinary, zero-valued, and disabled CTF states.
    subroutine make_ctf_parameters(ordinary_ctf, zero_ctf, no_ctf)
        type(ctfparams), intent(out) :: ordinary_ctf, zero_ctf, no_ctf

        ordinary_ctf%smpd = 1.0
        ordinary_ctf%kv = 300.0
        ordinary_ctf%cs = 2.7
        ordinary_ctf%fraca = 0.1
        ordinary_ctf%dfx = 1.4
        ordinary_ctf%dfy = 1.65
        ordinary_ctf%angast = 23.0
        ordinary_ctf%phshift = 0.37
        ordinary_ctf%ctfflag = CTFFLAG_YES
        zero_ctf = ordinary_ctf
        zero_ctf%cs = 0.0
        zero_ctf%fraca = 0.0
        zero_ctf%dfx = 0.0
        zero_ctf%dfy = 0.0
        zero_ctf%phshift = 0.0
        no_ctf = ordinary_ctf
        no_ctf%ctfflag = CTFFLAG_NO
    end subroutine make_ctf_parameters

!> Compare prepared components with a scalar CTF and separate whitening oracle.
    subroutine test_valid_case(workspace, name, ctfparms, sigma2, requested_range, unit, expect_zero_ctf)
        type(cartesian_pose_refiner), intent(in) :: workspace
        character(len=*), intent(in) :: name
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: requested_range(2), unit
        logical, intent(in) :: expect_zero_ctf
        type(cartesian_pose_data) :: data
        complex, allocatable :: raw_observed(:, :), actual_observed(:, :), actual_transfer(:, :)
        complex, allocatable :: expected_observed(:, :), expected_transfer(:, :)
        real(dp) :: max_observed_error, max_transfer_error, max_active_transfer
        integer :: effective_range(2), expected_count, status

        ! Pose data may use only shells represented by both the request and sigma2.
        effective_range = [max(0, requested_range(1)), &
            &min(requested_range(2), ubound(sigma2, 1), BOX/2)]
        call build_independent_components(ctfparms, sigma2, effective_range, raw_observed, &
            &expected_observed, expected_transfer)
        call workspace%prepare_particle(raw_observed, ctfparms, sigma2, requested_range, data, status)
        call assert_true(status == POSE_DATA_VALID .and. data%is_valid(), &
            &trim(name)//' did not produce valid prepared pose data')
        call assert_true(all(data%get_requested_shell_range() == requested_range), &
            &trim(name)//' did not retain the requested shell range')
        call assert_true(all(data%get_shell_range() == effective_range), &
            &trim(name)//' reported the wrong effective shell range')
        expected_count = count_active_samples(effective_range)
        call assert_true(data%get_active_sample_count() == expected_count, &
            &trim(name)//' reported the wrong full-disk sample count')
        call data%copy_components_test(actual_observed, actual_transfer)
        max_observed_error = maxval(abs(cmplx(actual_observed, kind=dp) - &
            &cmplx(expected_observed, kind=dp)))
        max_transfer_error = maxval(abs(cmplx(actual_transfer, kind=dp) - &
            &cmplx(expected_transfer, kind=dp)))
        call assert_true(max_observed_error <= COMPONENT_TOLERANCE, &
            &trim(name)//' observation whitening disagrees with the independent oracle')
        call assert_true(max_transfer_error <= COMPONENT_TOLERANCE, &
            &trim(name)//' transfer disagrees with the independent CTF/noise oracle')
        if (expect_zero_ctf) then
            max_active_transfer = max_active_abs(actual_transfer, effective_range)
            call assert_true(max_active_transfer <= COMPONENT_TOLERANCE, &
                &'a physical CTF zero did not remain a valid zero transfer')
        end if
        write (unit, '(a,5(a,i0),2(a,es24.16))') trim(name), achar(9), requested_range(1), &
            &achar(9), requested_range(2), achar(9), effective_range(1), achar(9), effective_range(2), &
            &achar(9), expected_count, achar(9), max_transfer_error, achar(9), max_observed_error
    end subroutine test_valid_case

!> Verify a nonfatal invalid-noise outcome and its retained range diagnostics.
    subroutine test_invalid_case(workspace, name, ctfparms, sigma2, requested_range, effective_range, unit)
        type(cartesian_pose_refiner), intent(in) :: workspace
        character(len=*), intent(in) :: name
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: requested_range(2), effective_range(2), unit
        type(cartesian_pose_data) :: data
        complex, allocatable :: raw_observed(:, :)
        integer :: lims2(2, 2), status

        lims2 = workspace%get_lims2()
        allocate (raw_observed(lims2(1, 1):lims2(1, 2), lims2(2, 1):lims2(2, 2)), &
            &source=cmplx(0., 0.))
        call workspace%prepare_particle(raw_observed, ctfparms, sigma2, requested_range, data, status)
        call assert_true(status == POSE_DATA_INVALID_NOISE_RANGE .and. .not. data%is_valid(), &
            &trim(name)//' did not return the nonfatal invalid-noise outcome')
        call assert_true(all(data%get_requested_shell_range() == requested_range), &
            &trim(name)//' did not retain its requested range')
        call assert_true(all(data%get_shell_range() == effective_range), &
            &trim(name)//' reported the wrong attempted effective range')
        call assert_true(data%get_active_sample_count() == 0, &
            &trim(name)//' retained active samples after invalid preparation')
        write (unit, '(a,a,i0,a,l1,5(a,i0))') trim(name), achar(9), status, achar(9), data%is_valid(), &
            &achar(9), requested_range(1), achar(9), requested_range(2), achar(9), effective_range(1), &
            &achar(9), effective_range(2), achar(9), data%get_active_sample_count()
    end subroutine test_invalid_case

!> Build raw CTF-distorted data first, then whiten it independently.
    subroutine build_independent_components(ctfparms, sigma2, effective_range, raw_observed, &
        &expected_observed, expected_transfer)
        type(ctfparams), intent(in) :: ctfparms
        real, intent(in) :: sigma2(0:)
        integer, intent(in) :: effective_range(2)
        complex, allocatable, intent(out) :: raw_observed(:, :), expected_observed(:, :)
        complex, allocatable, intent(out) :: expected_transfer(:, :)
        type(ctf) :: tfun
        complex :: signal
        real :: cval, inverse_sigma
        integer :: h, k, radius_squared, shell
        logical :: use_ctf

        allocate (raw_observed(-BOX/2:BOX/2, -BOX/2:BOX/2), source=cmplx(0., 0.))
        allocate (expected_observed(-BOX/2:BOX/2, -BOX/2:BOX/2), source=cmplx(0., 0.))
        allocate (expected_transfer(-BOX/2:BOX/2, -BOX/2:BOX/2), source=cmplx(0., 0.))
        use_ctf = ctfparms%ctfflag /= CTFFLAG_NO
        if (use_ctf) then
            tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
            call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        end if
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                radius_squared = h*h + k*k
                if (radius_squared > (BOX/2)**2) cycle
                cval = 1.0
                if (use_ctf) cval = scalar_ctf_value(tfun, ctfparms, h, k)
                signal = cmplx(0.7 + 0.03*real(h) - 0.02*real(k), &
                    &-0.4 + 0.01*real(h*k), kind=sp)
                raw_observed(h, k) = cval*signal
                if (radius_squared < effective_range(1)**2 .or. &
                    &radius_squared > effective_range(2)**2) cycle
                shell = nint(sqrt(real(radius_squared)))
                inverse_sigma = 1.0/sqrt(sigma2(shell))
                expected_transfer(h, k) = cval*inverse_sigma
                expected_observed(h, k) = raw_observed(h, k)*inverse_sigma
            end do
        end do
    end subroutine build_independent_components

!> Evaluate the scalar CTF API without calling the pose transfer routine.
    real function scalar_ctf_value(tfun, ctfparms, h, k) result(cval)
        type(ctf), intent(in) :: tfun
        type(ctfparams), intent(in) :: ctfparms
        integer, intent(in) :: h, k
        real :: angle, spatial_frequency_squared

        spatial_frequency_squared = real(h*h + k*k)/(real(BOX)*ctfparms%smpd)**2
        angle = atan2(real(k), real(h))
        cval = tfun%eval(spatial_frequency_squared, angle, ctfparms%phshift)
        if (ctfparms%ctfflag == CTFFLAG_FLIP) cval = abs(cval)
    end function scalar_ctf_value

!> Demonstrate the pose cap and unchanged PCG last-shell extension separately.
    subroutine test_pcg_last_shell_fallback(workspace, pcgop, no_ctf, short_sigma, unit)
        type(cartesian_pose_refiner), intent(in) :: workspace
        type(reconstructor_pcg), intent(in) :: pcgop
        type(ctfparams), intent(in) :: no_ctf
        real, intent(in) :: short_sigma(0:)
        integer, intent(in) :: unit
        type(cartesian_pose_data) :: data
        complex, allocatable :: raw_observed(:, :), expected_observed(:, :), expected_transfer(:, :)
        complex, allocatable :: pose_observed(:, :), pose_transfer(:, :)
        complex, allocatable :: pcg_short(:, :), pcg_extended(:, :)
        real, allocatable :: extended_sigma(:)
        real(dp) :: pcg_error, pcg_tail, pose_tail
        integer :: effective_range(2), lims2(2, 2), status

        call build_independent_components(no_ctf, short_sigma, [0, ubound(short_sigma, 1)], &
            &raw_observed, expected_observed, expected_transfer)
        call workspace%prepare_particle(raw_observed, no_ctf, short_sigma, [0, BOX/2], data, status)
        call assert_true(status == POSE_DATA_VALID, 'short pose variance did not remain valid')
        effective_range = data%get_shell_range()
        call data%copy_components_test(pose_observed, pose_transfer)
        allocate (extended_sigma(0:BOX/2), source=short_sigma(ubound(short_sigma, 1)))
        extended_sigma(0:ubound(short_sigma, 1)) = short_sigma
        lims2 = pcgop%get_lims2()
        allocate (pcg_short(lims2(1, 1):lims2(1, 2), lims2(2, 1):lims2(2, 2)))
        allocate (pcg_extended(lims2(1, 1):lims2(1, 2), lims2(2, 1):lims2(2, 2)))
        pcg_short = pcgop%build_transfer(no_ctf, [0., 0.], short_sigma)
        pcg_extended = pcgop%build_transfer(no_ctf, [0., 0.], extended_sigma)
        pcg_error = maxval(abs(cmplx(pcg_short, kind=dp) - cmplx(pcg_extended, kind=dp)))
        pose_tail = abs(cmplx(pose_transfer(BOX/2, 0), kind=dp))
        pcg_tail = abs(cmplx(pcg_short(BOX/2, 0), kind=dp))
        call assert_true(effective_range(2) == ubound(short_sigma, 1), &
            &'pose variance did not cap the effective upper shell')
        call assert_true(pose_tail == 0._dp, 'pose transfer extended its final known variance shell')
        call assert_true(pcg_tail > 0._dp .and. pcg_error <= COMPONENT_TOLERANCE, &
            &'PCG last-shell fallback no longer matches explicit tail extension')
        write (unit, '(i0,3(a,es24.16))') effective_range(2), achar(9), pose_tail, &
            &achar(9), pcg_tail, achar(9), pcg_error
    end subroutine test_pcg_last_shell_fallback

!> Count the exact full redundant disk for one effective radial range.
    pure integer function count_active_samples(shell_range) result(count)
        integer, intent(in) :: shell_range(2)
        integer :: h, k, radius_squared

        count = 0
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                radius_squared = h*h + k*k
                if (radius_squared < shell_range(1)**2 .or. &
                    &radius_squared > shell_range(2)**2) cycle
                count = count + 1
            end do
        end do
    end function count_active_samples

!> Return the largest active transfer magnitude in one radial range.
    pure real(dp) function max_active_abs(values, shell_range) result(maximum)
        complex, intent(in) :: values(-BOX/2:, -BOX/2:)
        integer, intent(in) :: shell_range(2)
        integer :: h, k, radius_squared

        maximum = 0._dp
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                radius_squared = h*h + k*k
                if (radius_squared < shell_range(1)**2 .or. &
                    &radius_squared > shell_range(2)**2) cycle
                maximum = max(maximum, abs(cmplx(values(h, k), kind=dp)))
            end do
        end do
    end function max_active_abs

!> Read one required old-style key=value argument.
    function required_argument(key) result(value)
        character(len=*), intent(in) :: key
        character(len=:), allocatable :: value
        character(len=4096) :: argument
        integer :: argument_status, iarg, separator

        value = ''
        do iarg = 1, command_argument_count()
            call get_command_argument(iarg, argument, status=argument_status)
            if (argument_status /= 0) error stop 'could not read CTF/sigma test argument'
            separator = index(argument, '=')
            if (separator <= 0) cycle
            if (trim(argument(:separator - 1)) /= key) cycle
            value = trim(argument(separator + 1:))
        end do
        if (len(value) == 0) error stop 'CTF/sigma test requires evidence_dir=<existing-directory>'
    end function required_argument

end module pose_cont_refinement_ctf_sigma_test
