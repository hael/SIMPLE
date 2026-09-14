module pose_cont_refinement_numerics_test
use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
use pose_cont_refinement_test_helpers, only: assert_true, build_test_volume, &
    &identity_rotation, prepare_unweighted_particle, TEST_BOX
use simple_defs, only: dp, DPI, KBALPHA, KBWINSZ, OSMPL_PAD_FAC
use simple_core_module_api, only: euler2m
use simple_cartesian_pose_refiner, only: cartesian_pose_refiner, cartesian_pose_data, &
    &right_increment_rotation
use simple_ctf, only: ctf
use simple_gridding, only: kb_stencil_centered_crop_inv_envelope_1d
use simple_kbinterpol, only: kbinterpol
use simple_type_defs, only: ctfparams, ctfvars, CTFFLAG_FLIP, CTFFLAG_NO, CTFFLAG_YES
implicit none
private
public :: run_pose_cont_numerics

real(dp), parameter :: GRADIENT_TOL = 3.e-2_dp
real(dp), parameter :: ORTHOGONAL_TOL = 2.e-12_dp

contains

    subroutine run_pose_cont_numerics()
        call test_prepared_particle_contract()
        call test_shift_phase_sign()
        call test_five_parameter_gradient()
        call test_rotation_increment()
        write(*,'(a)') 'POSE_CONT_REFINEMENT_NUMERICS: PASS'
    end subroutine run_pose_cont_numerics

    subroutine test_prepared_particle_contract()
        type(cartesian_pose_refiner) :: workspace, corrected_workspace
        type(cartesian_pose_data) :: data
        type(ctfparams) :: ctfparms
        type(ctf) :: tfun
        type(ctfvars) :: ctfvals
        type(kbinterpol) :: kbwin
        real, allocatable :: volume(:,:,:), corrected_volume(:,:,:), inv1d(:)
        real, allocatable :: sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real :: sigma2(0:TEST_BOX/2), short_sigma2(0:TEST_BOX/2-2), angle, cval, v
        complex :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex :: corrected_prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex :: raw_observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), shift(2), objective, gradient(5)
        integer :: effective_range(2), h, i, j, k, shell

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = real(euler2m([17.,31.,23.]),dp)
        shift = [0.23_dp,-0.17_dp]
        call workspace%predict_unweighted(rotation,shift,prediction)

        call prepare_unweighted_particle(workspace,prediction,data,[2,TEST_BOX/2-1])
        call assert_true(data%is_valid(),'valid unweighted particle was rejected')
        call assert_true(all(data%get_shell_range() == [2,TEST_BOX/2-1]), &
            &'prepared particle retained the wrong shell range')
        call workspace%prepared_objective_gradient(rotation,shift,data,objective,gradient)
        call assert_true(objective < 1.e-10_dp .and. maxval(abs(gradient)) < 1.e-8_dp, &
            &'exact unweighted particle has nonzero objective or gradient')
        call workspace%prepared_sigma_contribution(rotation,shift,data,sigma_contrib,ref_pow,ptcl_pow,v)
        call assert_true(maxval(abs(sigma_contrib)) < 1.e-8 .and. abs(v) < 1.e-8, &
            &'exact particle has nonzero unwhitened residual accounting')

        sigma2 = [(1.+0.1*real(shell),shell=0,TEST_BOX/2)]
        ctfparms%smpd = 1.
        ctfparms%kv = 300.
        ctfparms%cs = 2.7
        ctfparms%fraca = 0.1
        ctfparms%dfx = 1.4
        ctfparms%dfy = 1.65
        ctfparms%angast = 23.
        ctfparms%phshift = 0.37
        ctfparms%ctfflag = CTFFLAG_YES
        tfun = ctf(ctfparms%smpd,ctfparms%kv,ctfparms%cs,ctfparms%fraca)
        call tfun%init(ctfparms%dfx,ctfparms%dfy,ctfparms%angast)
        ctfvals = tfun%get_ctfvars(ctfparms%phshift)
        raw_observed = cmplx(0.,0.)
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if( h*h+k*k > (TEST_BOX/2)**2 ) cycle
                angle = 0.
                if( h /= 0 .or. k /= 0 ) angle = atan2(real(k),real(h))
                cval = tfun%eval_canonical(real(h*h+k*k)/real(TEST_BOX*TEST_BOX), &
                    &angle,ctfvals%phshift)
                raw_observed(h,k) = cval*prediction(h,k)
            enddo
        enddo
        call workspace%prepare_particle(raw_observed,ctfparms,sigma2,[2,TEST_BOX/2-1],data)
        call workspace%prepared_objective_gradient(rotation,shift,data,objective,gradient)
        call assert_true(objective < 1.e-9_dp .and. maxval(abs(gradient)) < 1.e-7_dp, &
            &'CTF and shell whitening changed an exact particle match')

        ctfparms%ctfflag = CTFFLAG_FLIP
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if( h*h+k*k > (TEST_BOX/2)**2 ) cycle
                angle = 0.
                if( h /= 0 .or. k /= 0 ) angle = atan2(real(k),real(h))
                cval = abs(tfun%eval_canonical(real(h*h+k*k)/real(TEST_BOX*TEST_BOX), &
                    &angle,ctfvals%phshift))
                raw_observed(h,k) = cval*prediction(h,k)
            enddo
        enddo
        call workspace%prepare_particle(raw_observed,ctfparms,sigma2,[2,TEST_BOX/2-1],data)
        call workspace%prepared_objective_gradient(rotation,shift,data,objective,gradient)
        call assert_true(objective < 1.e-9_dp .and. maxval(abs(gradient)) < 1.e-7_dp, &
            &'phase-flipped CTF changed an exact particle match')

        short_sigma2 = 1.
        call workspace%prepare_particle(prediction,ctfparms,short_sigma2,[2,TEST_BOX/2],data)
        effective_range = data%get_shell_range()
        call assert_true(data%is_valid() .and. effective_range(2) == TEST_BOX/2-2, &
            &'short noise spectrum did not cap the active shell range')

        sigma2(TEST_BOX/2) = -1.
        call workspace%prepare_particle(prediction,ctfparms,sigma2,[2,TEST_BOX/2],data)
        call assert_true(.not. data%is_valid(),'invalid active noise variance was accepted')
        sigma2 = 1.
        sigma2(1) = ieee_value(0.,ieee_quiet_nan)
        ctfparms%ctfflag = CTFFLAG_NO
        call workspace%prepare_particle(prediction,ctfparms,sigma2,[2,TEST_BOX/2],data)
        call assert_true(data%is_valid(),'invalid variance outside the active range rejected the particle')

        kbwin = kbinterpol(KBWINSZ,KBALPHA)
        call kb_stencil_centered_crop_inv_envelope_1d(kbwin,OSMPL_PAD_FAC*TEST_BOX,TEST_BOX,inv1d)
        allocate(corrected_volume,source=volume)
        do k = 1, TEST_BOX
            do j = 1, TEST_BOX
                do i = 1, TEST_BOX
                    corrected_volume(i,j,k) = corrected_volume(i,j,k)*inv1d(i)*inv1d(j)*inv1d(k)
                enddo
            enddo
        enddo
        call workspace%new_inverse_envelope_reference(volume)
        call corrected_workspace%new_physical_reference(corrected_volume)
        call workspace%predict_unweighted(identity_rotation(),[0._dp,0._dp],prediction)
        call corrected_workspace%predict_unweighted(identity_rotation(),[0._dp,0._dp],corrected_prediction)
        call assert_true(maxval(abs(prediction-corrected_prediction)) < 1.e-6, &
            &'inverse-envelope constructor did not apply the correction exactly once')
        call workspace%kill
        call corrected_workspace%kill
    end subroutine test_prepared_particle_contract

    subroutine test_shift_phase_sign()
        type(cartesian_pose_refiner) :: workspace
        real, allocatable :: volume(:,:,:)
        complex :: base(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex :: shifted(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: expected, ratio
        real(dp) :: rotation(3,3), shift(2), argument
        integer, parameter :: modes(2,3) = reshape([1,0,0,1,2,-1],[2,3])
        integer :: h, imode, k

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = identity_rotation()
        shift = [0.35_dp,-0.27_dp]
        call workspace%predict_unweighted(rotation,[0._dp,0._dp],base)
        call workspace%predict_unweighted(rotation,shift,shifted)
        do imode = 1, size(modes,2)
            h = modes(1,imode)
            k = modes(2,imode)
            call assert_true(abs(base(h,k)) > 1.e-6,'shift-sign test selected an absent mode')
            argument = 2._dp*DPI*(real(h,dp)*shift(1)+real(k,dp)*shift(2))/real(TEST_BOX,dp)
            expected = cmplx(cos(argument),sin(argument),kind=dp)
            ratio = cmplx(shifted(h,k),kind=dp)/cmplx(base(h,k),kind=dp)
            call assert_true(abs(ratio-expected) < 3.e-5_dp, &
                &'Fourier shift phase has the wrong sign or native-pixel scale')
        enddo
        call workspace%kill
    end subroutine test_shift_phase_sign

    subroutine test_five_parameter_gradient()
        type(cartesian_pose_refiner) :: workspace
        type(cartesian_pose_data) :: data
        real, allocatable :: volume(:,:,:)
        complex :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), truth_rotation(3,3), minus_rotation(3,3), plus_rotation(3,3)
        real(dp) :: shift(2), truth_shift(2), minus_shift(2), plus_shift(2)
        real(dp) :: objective, objective_minus, objective_plus, gradient(5), errors(5), basis(3), step
        integer :: axis, stencil_switches

        call build_test_volume(volume)
        call workspace%new_physical_reference(volume)
        rotation = real(euler2m([17.,31.,23.]),dp)
        truth_rotation = right_increment_rotation(rotation,[0.01_dp,-0.008_dp,0.006_dp])
        shift = [0.21_dp,-0.16_dp]
        truth_shift = [0.25_dp,-0.18_dp]
        call workspace%predict_unweighted(truth_rotation,truth_shift,observed)
        call prepare_unweighted_particle(workspace,observed,data,[2,4])
        call workspace%prepared_objective_gradient(rotation,shift,data,objective,gradient)
        do axis = 1, 5
            basis = 0._dp
            minus_rotation = rotation
            plus_rotation = rotation
            minus_shift = shift
            plus_shift = shift
            if( axis <= 3 )then
                step = 1.e-4_dp
                basis(axis) = 1._dp
                minus_rotation = right_increment_rotation(rotation,-step*basis)
                plus_rotation = right_increment_rotation(rotation,step*basis)
                stencil_switches = workspace%count_stencil_switches(rotation,minus_rotation,[2,4])+ &
                    &workspace%count_stencil_switches(rotation,plus_rotation,[2,4])
                call assert_true(stencil_switches == 0, &
                    &'rotation finite difference crossed an interpolation-stencil boundary')
            else
                step = 1.e-3_dp
                minus_shift(axis-3) = minus_shift(axis-3)-step
                plus_shift(axis-3) = plus_shift(axis-3)+step
            endif
            call workspace%prepared_objective_gradient(minus_rotation,minus_shift,data, &
                &objective_minus,gradient)
            call workspace%prepared_objective_gradient(plus_rotation,plus_shift,data, &
                &objective_plus,gradient)
            call workspace%prepared_objective_gradient(rotation,shift,data,objective,gradient)
            errors(axis) = abs((objective_plus-objective_minus)/(2._dp*step)-gradient(axis))/ &
                &max(1.e-5_dp,abs(gradient(axis)))
        enddo
        call assert_true(all(errors < GRADIENT_TOL), &
            &'a five-parameter objective-gradient component failed centred differences')
        call workspace%kill
    end subroutine test_five_parameter_gradient

    subroutine test_rotation_increment()
        real(dp) :: rotation(3,3), updated(3,3), identity(3,3)
        real(dp) :: determinant, input_determinant, input_orthogonality, updated_orthogonality

        rotation = real(euler2m([19.,37.,28.]),dp)
        updated = right_increment_rotation(rotation,[0.013_dp,-0.017_dp,0.011_dp])
        identity = identity_rotation()
        input_orthogonality = sqrt(sum((matmul(transpose(rotation),rotation)-identity)**2))
        input_determinant = determinant3(rotation)
        updated_orthogonality = sqrt(sum((matmul(transpose(updated),updated)-identity)**2))
        determinant = determinant3(updated)
        call assert_true(updated_orthogonality <= input_orthogonality+ORTHOGONAL_TOL, &
            &'right rotation increment increased the input orthogonality error')
        call assert_true(abs(determinant-1._dp) <= abs(input_determinant-1._dp)+ORTHOGONAL_TOL, &
            &'right rotation increment increased the input determinant error')
    end subroutine test_rotation_increment

    pure function determinant3(matrix) result(determinant)
        real(dp), intent(in) :: matrix(3,3)
        real(dp) :: determinant
        determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3)-matrix(2,3)*matrix(3,2))- &
            &matrix(1,2)*(matrix(2,1)*matrix(3,3)-matrix(2,3)*matrix(3,1))+ &
            &matrix(1,3)*(matrix(2,1)*matrix(3,2)-matrix(2,2)*matrix(3,1))
    end function determinant3

end module pose_cont_refinement_numerics_test
