module pose_cont_refinement_numerics_test
use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use pose_cont_refinement_test_helpers, only: assert_true, build_test_volume, &
    &identity_rotation, prepare_unweighted_particle, TEST_BOX
use simple_defs, only: dp, sp, DPI, KBALPHA, KBWINSZ, OSMPL_PAD_FAC
use simple_core_module_api, only: euler2m
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &right_increment_rotation, POSE_CONT_OBJECTIVE_CART_NCC, &
    &POSE_CONT_OBJECTIVE_CART_EUCLID
use simple_ctf, only: ctf
use simple_gridding, only: kb_stencil_centered_crop_inv_envelope_1d
use simple_image, only: image
use simple_kbinterpol, only: kbinterpol
use simple_projector, only: projector
use simple_type_defs, only: ctfparams, ctfvars, CTFFLAG_FLIP, CTFFLAG_NO, CTFFLAG_YES
implicit none
private
public :: run_pose_cont_numerics

real(dp), parameter :: GRADIENT_TOL = 3.e-2_dp
real(dp), parameter :: ORTHOGONAL_TOL = 2.e-12_dp
real(dp), parameter :: NCC_FORMULA_TOL = 20._dp*real(epsilon(1.), dp)

contains

    subroutine run_pose_cont_numerics()
        call test_prepared_particle_contract()
        call test_shift_phase_sign()
        call test_ncc_objective_formula()
        call test_five_parameter_gradient()
        call test_matched_projector_boundary()
        call test_rotation_increment()
        write (*, '(a)') 'POSE_CONT_REFINEMENT_NUMERICS: PASS'
    end subroutine run_pose_cont_numerics

    subroutine test_prepared_particle_contract()
        type(cartesian_pose_refiner) :: workspace, corrected_workspace
        type(cartesian_pose_data) :: data
        type(ctfparams) :: ctfparms
        type(ctf) :: tfun
        type(ctfvars) :: ctfvals
        type(kbinterpol) :: kbwin
        real, allocatable :: volume(:, :, :), corrected_volume(:, :, :), inv1d(:)
        real, allocatable :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real :: sigma2(0:TEST_BOX/2), short_sigma2(0:TEST_BOX/2 - 2), angle, cval, v
        complex :: prediction(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex :: corrected_prediction(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex :: raw_observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3, 3), shift(2), objective, gradient(5)
        integer :: effective_range(2), h, i, j, k, shell

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = real(euler2m([17., 31., 23.]), dp)
        shift = [0.23_dp, -0.17_dp]
        call workspace%predict_unweighted(rotation, shift, prediction)

        call prepare_unweighted_particle(workspace, prediction, data, [2, TEST_BOX/2 - 1])
        call assert_true(data%is_valid(), 'valid unweighted particle was rejected')
        call assert_true(all(data%get_shell_range() == [2, TEST_BOX/2 - 1]), &
            &'prepared particle retained the wrong shell range')
        call workspace%prepared_objective_gradient(rotation, shift, data, objective, gradient)
        call assert_true(objective < 1.e-10_dp .and. maxval(abs(gradient)) < 1.e-8_dp, &
            &'exact unweighted particle has nonzero objective or gradient')
        call workspace%prepared_sigma_contribution(rotation, shift, data, sigma_contrib, ref_pow, ptcl_pow, v)
        call assert_true(maxval(abs(sigma_contrib)) < 1.e-8 .and. abs(v) < 1.e-8, &
            &'exact particle has nonzero unwhitened residual accounting')

        sigma2 = [(1.+0.1*real(shell), shell=0, TEST_BOX/2)]
        ctfparms%smpd = 1.
        ctfparms%kv = 300.
        ctfparms%cs = 2.7
        ctfparms%fraca = 0.1
        ctfparms%dfx = 1.4
        ctfparms%dfy = 1.65
        ctfparms%angast = 23.
        ctfparms%phshift = 0.37
        ctfparms%ctfflag = CTFFLAG_YES
        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        ctfvals = tfun%get_ctfvars(ctfparms%phshift)
        raw_observed = cmplx(0., 0.)
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if (h*h + k*k > (TEST_BOX/2)**2) cycle
                angle = 0.
                if (h /= 0 .or. k /= 0) angle = atan2(real(k), real(h))
                cval = tfun%eval_canonical(real(h*h + k*k)/real(TEST_BOX*TEST_BOX), &
                    &angle, ctfvals%phshift)
                raw_observed(h, k) = cval*prediction(h, k)
            end do
        end do
        call workspace%prepare_particle(raw_observed, ctfparms, sigma2, [2, TEST_BOX/2 - 1], data)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective, gradient)
        call assert_true(objective < 1.e-9_dp .and. maxval(abs(gradient)) < 1.e-7_dp, &
            &'CTF and shell whitening changed an exact particle match')

        ctfparms%ctfflag = CTFFLAG_FLIP
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if (h*h + k*k > (TEST_BOX/2)**2) cycle
                angle = 0.
                if (h /= 0 .or. k /= 0) angle = atan2(real(k), real(h))
                cval = abs(tfun%eval_canonical(real(h*h + k*k)/real(TEST_BOX*TEST_BOX), &
                    &angle, ctfvals%phshift))
                raw_observed(h, k) = cval*prediction(h, k)
            end do
        end do
        call workspace%prepare_particle(raw_observed, ctfparms, sigma2, [2, TEST_BOX/2 - 1], data)
        call workspace%prepared_objective_gradient(rotation, shift, data, objective, gradient)
        call assert_true(objective < 1.e-9_dp .and. maxval(abs(gradient)) < 1.e-7_dp, &
            &'phase-flipped CTF changed an exact particle match')

        short_sigma2 = 1.
        call workspace%prepare_particle(prediction, ctfparms, short_sigma2, [2, TEST_BOX/2], data)
        effective_range = data%get_shell_range()
        call assert_true(data%is_valid() .and. effective_range(2) == TEST_BOX/2 - 2, &
            &'short noise spectrum did not cap the active shell range')

        sigma2(TEST_BOX/2) = -1.
        call workspace%prepare_particle(prediction, ctfparms, sigma2, [2, TEST_BOX/2], data)
        call assert_true(.not. data%is_valid(), 'invalid active noise variance was accepted')
        sigma2 = 1.
        sigma2(1) = ieee_value(0., ieee_quiet_nan)
        ctfparms%ctfflag = CTFFLAG_NO
        call workspace%prepare_particle(prediction, ctfparms, sigma2, [2, TEST_BOX/2], data)
        call assert_true(data%is_valid(), 'invalid variance outside the active range rejected the particle')

        kbwin = kbinterpol(KBWINSZ, KBALPHA)
        call kb_stencil_centered_crop_inv_envelope_1d(kbwin, OSMPL_PAD_FAC*TEST_BOX, TEST_BOX, inv1d)
        allocate (corrected_volume, source=volume)
        do k = 1, TEST_BOX
            do j = 1, TEST_BOX
                do i = 1, TEST_BOX
                    corrected_volume(i, j, k) = corrected_volume(i, j, k)*inv1d(i)*inv1d(j)*inv1d(k)
                end do
            end do
        end do
        call workspace%new_inverse_envelope_reference(volume)
        call corrected_workspace%new_physical_reference(corrected_volume)
        call workspace%predict_unweighted(identity_rotation(), [0._dp, 0._dp], prediction)
        call corrected_workspace%predict_unweighted(identity_rotation(), [0._dp, 0._dp], corrected_prediction)
        call assert_true(maxval(abs(prediction - corrected_prediction)) < 1.e-6, &
            &'inverse-envelope constructor did not apply the correction exactly once')
        call workspace%kill
        call corrected_workspace%kill
    end subroutine test_prepared_particle_contract

    subroutine test_shift_phase_sign()
        type(cartesian_pose_refiner) :: workspace
        real, allocatable :: volume(:, :, :)
        complex :: base(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex :: shifted(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: expected, ratio
        real(dp) :: rotation(3, 3), shift(2), argument
        integer, parameter :: modes(2, 3) = reshape([1, 0, 0, 1, 2, -1], [2, 3])
        integer :: h, imode, k

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = identity_rotation()
        shift = [0.35_dp, -0.27_dp]
        call workspace%predict_unweighted(rotation, [0._dp, 0._dp], base)
        call workspace%predict_unweighted(rotation, shift, shifted)
        do imode = 1, size(modes, 2)
            h = modes(1, imode)
            k = modes(2, imode)
            call assert_true(abs(base(h, k)) > 1.e-6, 'shift-sign test selected an absent mode')
            argument = 2._dp*DPI*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(TEST_BOX, dp)
            expected = cmplx(cos(argument), sin(argument), kind=dp)
            ratio = cmplx(shifted(h, k), kind=dp)/cmplx(base(h, k), kind=dp)
            call assert_true(abs(ratio - expected) < 3.e-5_dp, &
                &'Fourier shift phase has the wrong sign or native-pixel scale')
        end do
        call workspace%kill
    end subroutine test_shift_phase_sign

    ! Independently evaluate 1-NCC and verify that positive particle-amplitude
    ! scaling changes neither the objective nor its five pose derivatives.
    subroutine test_ncc_objective_formula()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        real, allocatable :: volume(:, :, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex :: prediction(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: cross_sum, particle, model
        real(dp) :: truth_rotation(3, 3), candidate_rotation(3, 3)
        real(dp) :: truth_shift(2), candidate_shift(2), scales(3)
        real(dp) :: objective, oracle, baseline_objective
        real(dp) :: gradient(5), baseline_gradient(5)
        real(dp) :: particle_power, prediction_power
        integer :: h, i, k

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        truth_rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift = [0.31_dp, -0.24_dp]
        candidate_rotation = real(euler2m([20., 36.2, 28.7]), dp)
        candidate_shift = [-0.08_dp, 0.06_dp]
        call workspace%predict_unweighted(truth_rotation, truth_shift, observed)
        call workspace%predict_unweighted(candidate_rotation, candidate_shift, prediction)
        scales = [0.5_dp, 1._dp, 2._dp]

        do i = 1, size(scales)
            call prepare_unweighted_particle(workspace, real(scales(i), kind=kind(observed))*observed, &
                &data, [2, 4])
            call workspace%prepared_objective_gradient(candidate_rotation, candidate_shift, data, &
                &objective, gradient, POSE_CONT_OBJECTIVE_CART_NCC)

            cross_sum = cmplx(0._dp, 0._dp, kind=dp)
            particle_power = 0._dp
            prediction_power = 0._dp
            do k = -TEST_BOX/2, TEST_BOX/2
                do h = -TEST_BOX/2, TEST_BOX/2
                    if (h*h + k*k < 2**2 .or. h*h + k*k > 4**2) cycle
                    particle = cmplx(scales(i)*observed(h, k), kind=dp)
                    model = cmplx(prediction(h, k), kind=dp)
                    cross_sum = cross_sum + conjg(particle)*model
                    particle_power = particle_power + real(conjg(particle)*particle, dp)
                    prediction_power = prediction_power + real(conjg(model)*model, dp)
                end do
            end do
            oracle = 1._dp - real(cross_sum, dp)/sqrt(particle_power*prediction_power)
            call assert_true(abs(objective - oracle) < NCC_FORMULA_TOL, &
                &'Cartesian NCC backend disagrees with the independent formula')
            if (i == 1) then
                baseline_objective = objective
                baseline_gradient = gradient
            else
                call assert_true(abs(objective - baseline_objective) < 1.e-12_dp .and. &
                    &maxval(abs(gradient - baseline_gradient)) < 1.e-10_dp, &
                    &'Cartesian NCC changed under positive particle-amplitude scaling')
            end if
        end do
        call workspace%kill()
    end subroutine test_ncc_objective_formula

    subroutine test_five_parameter_gradient()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        real, allocatable :: volume(:, :, :)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3, 3), truth_rotation(3, 3), minus_rotation(3, 3), plus_rotation(3, 3)
        real(dp) :: shift(2), truth_shift(2), minus_shift(2), plus_shift(2)
        real(dp) :: objective, objective_minus, objective_plus, gradient(5), errors(5), basis(3), step
        integer :: axis, objective_kind

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = real(euler2m([17., 31., 23.]), dp)
        truth_rotation = right_increment_rotation(rotation, [0.01_dp, -0.008_dp, 0.006_dp])
        shift = [0.21_dp, -0.16_dp]
        truth_shift = [0.25_dp, -0.18_dp]
        call workspace%predict_unweighted(truth_rotation, truth_shift, observed)
        call prepare_unweighted_particle(workspace, observed, data, [2, 4])
        do objective_kind = POSE_CONT_OBJECTIVE_CART_NCC, POSE_CONT_OBJECTIVE_CART_EUCLID
            call workspace%prepared_objective_gradient(rotation, shift, data, objective, gradient, &
                &objective_kind)
            do axis = 1, 5
                basis = 0._dp
                minus_rotation = rotation
                plus_rotation = rotation
                minus_shift = shift
                plus_shift = shift
                if (axis <= 3) then
                    step = 1.e-4_dp
                    basis(axis) = 1._dp
                    minus_rotation = right_increment_rotation(rotation, -step*basis)
                    plus_rotation = right_increment_rotation(rotation, step*basis)
                else
                    step = 1.e-3_dp
                    minus_shift(axis - 3) = minus_shift(axis - 3) - step
                    plus_shift(axis - 3) = plus_shift(axis - 3) + step
                end if
                call workspace%prepared_objective_gradient(minus_rotation, minus_shift, data, &
                    &objective_minus, gradient, objective_kind)
                call workspace%prepared_objective_gradient(plus_rotation, plus_shift, data, &
                    &objective_plus, gradient, objective_kind)
                call workspace%prepared_objective_gradient(rotation, shift, data, objective, gradient, &
                    &objective_kind)
                errors(axis) = abs((objective_plus - objective_minus)/(2._dp*step) - gradient(axis))/ &
                    &max(1.e-5_dp, abs(gradient(axis)))
            end do
            call assert_true(all(errors < GRADIENT_TOL), &
                &'a five-parameter objective-gradient component failed centred differences')
        end do
        call workspace%kill
    end subroutine test_five_parameter_gradient

    !> Compare the Cartesian gather with SIMPLE's established projector kernel
    !! at identical rotated 3-D coordinates from the same physical reference.
    subroutine test_matched_projector_boundary()
        type(cartesian_pose_refiner) :: workspace
        type(image) :: volume_image
        type(projector) :: pftc_projector
        real, allocatable :: volume(:, :, :)
        complex :: cartesian(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex :: pftc(-TEST_BOX/2:TEST_BOX/2, -TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: cross_sum
        real(dp) :: rotation(3, 3), cartesian_power, pftc_power
        real(dp) :: correlation, gain, relative_l2, residual_power
        real(sp) :: loc(3)
        integer :: h, k, radius_squared

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        call volume_image%new([TEST_BOX, TEST_BOX, TEST_BOX], 1.)
        call volume_image%set_rmat(volume, .false.)
        call pftc_projector%new([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, &
            &OSMPL_PAD_FAC*TEST_BOX], 1.)
        call volume_image%pad_fft(pftc_projector)
        call pftc_projector%expand_cmat()

        rotation = real(euler2m([17., 31., 23.]), dp)
        call workspace%predict_unweighted(rotation, [0._dp, 0._dp], cartesian)
        pftc = cmplx(0., 0.)
        cross_sum = cmplx(0._dp, 0._dp, kind=dp)
        cartesian_power = 0._dp
        pftc_power = 0._dp
        residual_power = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                radius_squared = h*h + k*k
                if (radius_squared < 4 .or. radius_squared > (TEST_BOX/2 - 1)**2) cycle
                loc = real(matmul(real([h, k, 0], dp), rotation), sp)
                pftc(h, k) = pftc_projector%interp_fcomp_oversamp(loc)
                cross_sum = cross_sum + conjg(cmplx(pftc(h, k), kind=dp))* &
                    &cmplx(cartesian(h, k), kind=dp)
                pftc_power = pftc_power + abs(cmplx(pftc(h, k), kind=dp))**2
                cartesian_power = cartesian_power + abs(cmplx(cartesian(h, k), kind=dp))**2
                residual_power = residual_power + &
                    &abs(cmplx(cartesian(h, k) - pftc(h, k), kind=dp))**2
            end do
        end do
        correlation = real(cross_sum, dp)/sqrt(pftc_power*cartesian_power)
        gain = real(cross_sum, dp)/pftc_power
        relative_l2 = sqrt(residual_power/pftc_power)
        write (*, '(a,3(1x,es12.4))') &
            &'POSE_CONT_MATCHED_PROJECTOR', correlation, gain, relative_l2
        call assert_true(correlation >= 1._dp - 1.e-6_dp .and. &
            &abs(gain - 1._dp) <= 1.e-5_dp .and. relative_l2 <= 2.e-5_dp, &
            &'PFTC projector kernel and Cartesian gather disagree at a matched boundary')

        call workspace%kill()
        call pftc_projector%kill_expanded()
        call pftc_projector%kill()
        call volume_image%kill()
    end subroutine test_matched_projector_boundary

    subroutine test_rotation_increment()
        real(dp) :: rotation(3, 3), updated(3, 3), identity(3, 3)
        real(dp) :: determinant, input_determinant, input_orthogonality, updated_orthogonality

        rotation = real(euler2m([19., 37., 28.]), dp)
        updated = right_increment_rotation(rotation, [0.013_dp, -0.017_dp, 0.011_dp])
        identity = identity_rotation()
        input_orthogonality = sqrt(sum((matmul(transpose(rotation), rotation) - identity)**2))
        input_determinant = determinant3(rotation)
        updated_orthogonality = sqrt(sum((matmul(transpose(updated), updated) - identity)**2))
        determinant = determinant3(updated)
        call assert_true(updated_orthogonality <= input_orthogonality + ORTHOGONAL_TOL, &
            &'right rotation increment increased the input orthogonality error')
        call assert_true(abs(determinant - 1._dp) <= abs(input_determinant - 1._dp) + ORTHOGONAL_TOL, &
            &'right rotation increment increased the input determinant error')
    end subroutine test_rotation_increment

    pure function determinant3(matrix) result(determinant)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp) :: determinant
        determinant = matrix(1, 1)*(matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
            &matrix(1, 2)*(matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            &matrix(1, 3)*(matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant3

end module pose_cont_refinement_numerics_test
