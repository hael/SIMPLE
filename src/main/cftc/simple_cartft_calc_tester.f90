!@descr: unit tests for the Cartesian Fourier calculator of continuous pose refinement (simple_cartft_calc)
! Particle preparation, objectives and gradients, the gather against the PFTC projector, reference
! volumes and their files, batch preparation, the sigma owner and its polish fallback, the
! per-shell sigma2 against the polar calculator (noise-only, model power, and the group sigma2 of
! noisy particles) and the canonical sigma2 round trip. Test ids (E*, N*) are those of the
! pose_cont refactoring plan, section 10.
module simple_cartft_calc_tester
use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value, ieee_is_finite
use, intrinsic :: iso_fortran_env, only: int64, real32
use simple_defs,            only: dp, sp, DPI, OSMPL_PAD_FAC
use simple_core_module_api, only: euler2m, string, file_exists, del_file
use simple_cartft_calc,     only: cartft_calc, CART_REFVOLS_FORMAT_VERSION
use simple_builder,         only: builder
use simple_matcher_refvol_utils, only: read_reprojection_model, remove_ref_section_files
use simple_matcher_ptcl_batch, only: prep_cart_batch, prep_sigmas_objfun
use simple_euclid_sigma2,   only: euclid_sigma2
!$ use omp_lib,             only: omp_get_max_threads, omp_set_num_threads
use simple_refine3D_fnames, only: refine3D_cart_refvols_fname
use simple_cartft_pose_opt, only: cartft_pose_opt, right_increment_rotation, CARTFT_ACCEPTED, CARTFT_BOUND_REJECTED
use simple_cmdline,         only: cmdline
use simple_ctf,             only: ctf
use simple_image,           only: image
use simple_memoize_ft_maps, only: memoize_ft_maps, forget_ft_maps
use simple_parameters,      only: parameters
use simple_polarft_calc,    only: polarft_calc, vol_pad2ref_pfts_opt
use simple_projector,       only: projector
use simple_oris,            only: oris
use simple_sp_project,      only: sp_project
use simple_sigma2_state_file, only: sigma2_state_header, sigma2_state_init_header, sigma2_state_create_candidate, &
    &sigma2_state_write_local_range, sigma2_state_read_particles, sigma2_state_read_groups, sigma2_state_read_header, &
    &SIGMA2_GROUP_GLOBAL, SIGMA2_PROV_RESIDUAL
use simple_sigma2_state,    only: sigma2_state_candidate_path, sigma2_state_range_path, sigma2_state_prepare_update, &
    &sigma2_state_merge_local_ranges, sigma2_state_reduce_groups, sigma2_state_commit
use simple_matcher_2Dprep,  only: prepimg4align_cart
use simple_type_defs,       only: ctfparams, ctfvars, CTFFLAG_FLIP, CTFFLAG_NO, CTFFLAG_YES, OBJFUN_CC, OBJFUN_EUCLID
use simple_test_utils
implicit none
private
public :: run_all_cartft_calc_tests

integer,  parameter :: TEST_BOX        = 24
integer,  parameter :: SMALL_BOX       = 16
integer(int64), parameter :: FIXTURE_DIGEST = 3301_int64
real(dp), parameter :: GRADIENT_TOL    = 3.e-2_dp
! single-precision samples summed in double precision: a few single-precision ulps
real(dp), parameter :: FORMULA_TOL     = 20._dp*real(epsilon(1.), dp)
real,     parameter :: LINEARITY_TOL   = 1.e-5

contains

    subroutine run_all_cartft_calc_tests()
        write(*,'(A)') '**** running all Cartesian calculator tests ****'
        write(*,'(A)') 'test_particle_preparation_contract'
        call test_particle_preparation_contract()
        write(*,'(A)') 'test_sigma_shell_contract'
        call test_sigma_shell_contract()
        write(*,'(A)') 'test_shift_phase_sign'
        call test_shift_phase_sign()
        write(*,'(A)') 'test_cc_objective_formula'
        call test_cc_objective_formula()
        write(*,'(A)') 'test_euclid_objective_formula'
        call test_euclid_objective_formula()
        write(*,'(A)') 'test_five_parameter_gradient'
        call test_five_parameter_gradient()
        write(*,'(A)') 'test_matched_projector_boundary'
        call test_matched_projector_boundary()
        write(*,'(A)') 'test_reference_lifecycle'
        call test_reference_lifecycle()
        write(*,'(A)') 'test_even_odd_references'
        call test_even_odd_references()
        write(*,'(A)') 'test_sigma_endpoint_contracts'
        call test_sigma_endpoint_contracts()
        write(*,'(A)') 'test_polar_grid_pose_agreement'
        call test_polar_grid_pose_agreement()
        write(*,'(A)') 'test_observation_preparation'
        call test_observation_preparation()
        write(*,'(A)') 'test_reference_volume_file'
        call test_reference_volume_file()
        write(*,'(A)') 'test_reference_volume_handoff'
        call test_reference_volume_handoff()
        write(*,'(A)') 'test_observation_established_path'
        call test_observation_established_path()
        write(*,'(A)') 'test_batch_preparation'
        call test_batch_preparation()
        write(*,'(A)') 'test_concurrent_evaluation'
        call test_concurrent_evaluation()
        write(*,'(A)') 'test_sigma_owner'
        call test_sigma_owner()
        write(*,'(A)') 'test_polar_sigma_agreement'
        call test_polar_sigma_agreement()
        write(*,'(A)') 'test_polish_sigma_fallback'
        call test_polish_sigma_fallback()
        write(*,'(A)') 'test_polar_sigma_signal'
        call test_polar_sigma_signal()
        write(*,'(A)') 'test_polar_sigma_group'
        call test_polar_sigma_group()
        write(*,'(A)') 'test_canonical_sigma_round_trip'
        call test_canonical_sigma_round_trip()
    end subroutine run_all_cartft_calc_tests

    ! FIXTURES

    !> Four Gaussian blobs in a TEST_BOX cube (the refiner tester's volume), or in a cube of
    !! box when given.
    subroutine build_test_volume( volume, box )
        real, allocatable, intent(out) :: volume(:,:,:)
        integer, optional, intent(in)  :: box
        real, parameter :: centres(3,4) = reshape([-5., -3., 2., 4., 5., -3., 0., -6., -5., 3., -2., 6.], [3,4])
        real, parameter :: sigmas(4)     = [2., 2.5, 1.8, 2.2]
        real, parameter :: amplitudes(4) = [1., 0.8, 0.6, 0.5]
        real    :: centre, dx, dy, dz
        integer :: blob, i, j, k, n
        n = TEST_BOX
        if( present(box) ) n = box
        allocate(volume(n,n,n), source=0.)
        centre = real(n)/2. + 0.5
        do k = 1, n
            do j = 1, n
                do i = 1, n
                    do blob = 1, 4
                        dx = real(i) - centre - centres(1,blob)
                        dy = real(j) - centre - centres(2,blob)
                        dz = real(k) - centre - centres(3,blob)
                        volume(i,j,k) = volume(i,j,k) + amplitudes(blob)*exp(-(dx*dx + dy*dy + dz*dz)/(2.*sigmas(blob)**2))
                    end do
                end do
            end do
        end do
    end subroutine build_test_volume

    !> A smooth SMALL_BOX volume (the adapter tester's volume).
    subroutine build_small_volume( volume )
        real, allocatable, intent(out) :: volume(:,:,:)
        integer :: i, j, k
        allocate(volume(SMALL_BOX,SMALL_BOX,SMALL_BOX))
        do k = 1, SMALL_BOX
            do j = 1, SMALL_BOX
                do i = 1, SMALL_BOX
                    volume(i,j,k) = sin(0.11*real(i) + 0.17*real(j) - 0.07*real(k)) + 0.02*real(i*j - k) + 0.003*real(i*k)
                end do
            end do
        end do
    end subroutine build_small_volume

    !> A one-state calculator with the same volume in both halves and one particle slot.
    subroutine new_test_calc( calc, volume )
        type(cartft_calc), intent(inout) :: calc
        real,              intent(in)    :: volume(:,:,:)
        call calc%new(1, size(volume,1), 1)
        call calc%set_ref(1, .true.,  volume)
        call calc%set_ref(1, .false., volume)
    end subroutine new_test_calc

    !> Slot 1 for OBJFUN_EUCLID from an observation without CTF and with unit sigma2.
    subroutine set_unweighted_particle( calc, observed, kfromto )
        type(cartft_calc), intent(inout) :: calc
        complex,           intent(in)    :: observed(:,:)
        integer, optional, intent(in)    :: kfromto(2)
        type(ctfparams) :: no_ctf
        real, allocatable :: sigma2(:)
        integer :: active_range(2), box
        box = calc%get_box()
        allocate(sigma2(0:box/2), source=1.)
        no_ctf%ctfflag = CTFFLAG_NO
        active_range   = [2, box/2]
        if( present(kfromto) ) active_range = kfromto
        call calc%set_ptcl(1, observed, no_ctf, sigma2, active_range)
    end subroutine set_unweighted_particle

    !> sum |X|^2/sigma2 over the pixels whose shell nint(|(h,k)|) lies in kfromto (unit sigma2
    !! when absent): the normalizer of the euclid loss.
    function whitened_power( observed, kfromto, sigma2 ) result( power )
        complex,        intent(in) :: observed(-TEST_BOX/2:,-TEST_BOX/2:)
        integer,        intent(in) :: kfromto(2)
        real, optional, intent(in) :: sigma2(0:)
        real(dp) :: power
        integer  :: h, k, shell
        power = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                if( present(sigma2) )then
                    power = power + abs(cmplx(observed(h,k), kind=dp))**2/real(sigma2(shell), dp)
                else
                    power = power + abs(cmplx(observed(h,k), kind=dp))**2
                endif
            end do
        end do
    end function whitened_power

    pure function identity_rotation() result( rotation )
        real(dp) :: rotation(3,3)
        rotation      = 0._dp
        rotation(1,1) = 1._dp
        rotation(2,2) = 1._dp
        rotation(3,3) = 1._dp
    end function identity_rotation

    ! TESTS

    ! E4
    subroutine test_particle_preparation_contract()
        type(cartft_calc) :: calc
        type(ctfparams)   :: ctfparms
        type(ctf)         :: tfun
        type(ctfvars)     :: ctfvals
        real, allocatable :: volume(:,:,:), sigma_contrib(:), ref_pow(:), ptcl_pow(:)
        real     :: sigma2(0:TEST_BOX/2), angle, cval, v
        complex  :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex  :: raw_observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), shift(2), objective, gradient(5), power
        integer  :: h, k, shell
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        rotation = real(euler2m([17., 31., 23.]), dp)
        shift    = [0.23_dp, -0.17_dp]
        call calc%predict(1, .true., rotation, shift, prediction)
        call set_unweighted_particle(calc, prediction, [2, TEST_BOX/2 - 1])
        call assert_true(calc%ptcl_is_valid(1), 'valid unweighted particle was rejected')
        call assert_true(all(calc%get_ptcl_kfromto(1) == [2, TEST_BOX/2 - 1]), 'prepared particle retained the wrong shell range')
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
        ! the thresholds of the unnormalized half squared residual this check was written for,
        ! carried into the normalized loss: L = 2E/P and grad L = 2 grad E/P, P the whitened power
        power = whitened_power(prediction, [2, TEST_BOX/2 - 1])
        call assert_true(objective*power/2._dp < 1.e-10_dp .and. maxval(abs(gradient))*power/2._dp < 1.e-8_dp, &
            &'exact unweighted particle has nonzero objective or gradient')
        call calc%sigma_contribution(1, .true., 1, rotation, shift, sigma_contrib, ref_pow, ptcl_pow, v)
        call assert_true(maxval(abs(sigma_contrib)) < 1.e-8 .and. abs(v) < 1.e-8, &
            &'exact particle has nonzero unwhitened residual accounting')
        sigma2           = [(1. + 0.1*real(shell), shell = 0, TEST_BOX/2)]
        ctfparms%smpd    = 1.
        ctfparms%kv      = 300.
        ctfparms%cs      = 2.7
        ctfparms%fraca   = 0.1
        ctfparms%dfx     = 1.4
        ctfparms%dfy     = 1.65
        ctfparms%angast  = 23.
        ctfparms%phshift = 0.37
        ctfparms%ctfflag = CTFFLAG_YES
        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        ctfvals = tfun%get_ctfvars(ctfparms%phshift)
        raw_observed = cmplx(0., 0.)
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if( h*h + k*k > (TEST_BOX/2)**2 ) cycle
                angle = 0.
                if( h /= 0 .or. k /= 0 ) angle = atan2(real(k), real(h))
                ! CTFFLAG_YES: the observation is phase-flipped before the mask (O4), so an exact
                ! observation carries abs(CTF)
                cval = abs(tfun%eval_canonical(real(h*h + k*k)/real(TEST_BOX*TEST_BOX), angle, ctfvals%phshift))
                raw_observed(h,k) = cval*prediction(h,k)
            end do
        end do
        call calc%set_ptcl(1, raw_observed, ctfparms, sigma2, [2, TEST_BOX/2 - 1])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
        power = whitened_power(raw_observed, [2, TEST_BOX/2 - 1], sigma2)
        call assert_true(objective*power/2._dp < 1.e-9_dp .and. maxval(abs(gradient))*power/2._dp < 1.e-7_dp, &
            &'CTF and shell whitening changed an exact phase-flipped particle match')
        ! the same exact match under the correlation objective, prepared without sigma2
        call calc%set_ptcl(1, raw_observed, ctfparms, [2, TEST_BOX/2 - 1])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
        call assert_true(abs(objective) < 1.e-6_dp .and. maxval(abs(gradient)) < 1.e-5_dp, &
            &'CTF changed an exact phase-flipped particle match under cc')
        ctfparms%ctfflag = CTFFLAG_FLIP
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                if( h*h + k*k > (TEST_BOX/2)**2 ) cycle
                angle = 0.
                if( h /= 0 .or. k /= 0 ) angle = atan2(real(k), real(h))
                cval = abs(tfun%eval_canonical(real(h*h + k*k)/real(TEST_BOX*TEST_BOX), angle, ctfvals%phshift))
                raw_observed(h,k) = cval*prediction(h,k)
            end do
        end do
        call calc%set_ptcl(1, raw_observed, ctfparms, sigma2, [2, TEST_BOX/2 - 1])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
        power = whitened_power(raw_observed, [2, TEST_BOX/2 - 1], sigma2)
        call assert_true(objective*power/2._dp < 1.e-9_dp .and. maxval(abs(gradient))*power/2._dp < 1.e-7_dp, &
            &'phase-flipped CTF changed an exact particle match')
        call calc%kill
    end subroutine test_particle_preparation_contract

    ! E5
    subroutine test_sigma_shell_contract()
        type(cartft_calc) :: calc
        type(ctfparams)   :: ctfparms
        real, allocatable :: volume(:,:,:)
        real     :: sigma2(0:TEST_BOX/2), short_sigma2(0:TEST_BOX/2 - 2)
        complex  :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), shift(2)
        integer  :: effective_range(2)
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        rotation = real(euler2m([17., 31., 23.]), dp)
        shift    = [0.23_dp, -0.17_dp]
        call calc%predict(1, .true., rotation, shift, prediction)
        ctfparms%ctfflag = CTFFLAG_NO
        short_sigma2 = 1.
        call calc%set_ptcl(1, prediction, ctfparms, short_sigma2, [2, TEST_BOX/2])
        effective_range = calc%get_ptcl_kfromto(1)
        call assert_true(calc%ptcl_is_valid(1) .and. effective_range(2) == TEST_BOX/2 - 2, &
            &'short noise spectrum did not cap the active shell range')
        sigma2 = 1.
        sigma2(TEST_BOX/2) = -1.
        call calc%set_ptcl(1, prediction, ctfparms, sigma2, [2, TEST_BOX/2])
        call assert_true(.not. calc%ptcl_is_valid(1), 'invalid active noise variance was accepted')
        sigma2 = 1.
        sigma2(1) = ieee_value(0., ieee_quiet_nan)
        call calc%set_ptcl(1, prediction, ctfparms, sigma2, [2, TEST_BOX/2])
        call assert_true(calc%ptcl_is_valid(1), 'invalid variance outside the active range rejected the particle')
        call calc%kill
    end subroutine test_sigma_shell_contract

    ! E7
    subroutine test_shift_phase_sign()
        type(cartft_calc) :: calc
        real, allocatable :: volume(:,:,:)
        complex     :: base(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex     :: shifted(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: expected, ratio
        real(dp)    :: rotation(3,3), shift(2), argument
        integer, parameter :: modes(2,3) = reshape([1, 0, 0, 1, 2, -1], [2,3])
        integer     :: h, imode, k
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        rotation = identity_rotation()
        shift    = [0.35_dp, -0.27_dp]
        call calc%predict(1, .true., rotation, [0._dp, 0._dp], base)
        call calc%predict(1, .true., rotation, shift, shifted)
        do imode = 1, size(modes,2)
            h = modes(1,imode)
            k = modes(2,imode)
            call assert_true(abs(base(h,k)) > 1.e-6, 'shift-sign test selected an absent mode')
            argument = 2._dp*DPI*(real(h, dp)*shift(1) + real(k, dp)*shift(2))/real(TEST_BOX, dp)
            expected = cmplx(cos(argument), sin(argument), kind=dp)
            ratio    = cmplx(shifted(h,k), kind=dp)/cmplx(base(h,k), kind=dp)
            call assert_true(abs(ratio - expected) < 3.e-5_dp, 'Fourier shift phase has the wrong sign or native-pixel scale')
        end do
        call calc%kill
    end subroutine test_shift_phase_sign

    ! N3 (rewrites E8): 1 - cc with a uniform weight per Cartesian pixel over the pixels whose
    ! shell nint(|(h,k)|) lies in the range, against a closed form over the test arrays;
    ! invariant to the particle gain; the slot is prepared and evaluated without any sigma2
    ! (the cc preparation takes none, C5), and its score is cc
    subroutine test_cc_objective_formula()
        type(cartft_calc) :: calc
        type(ctfparams)   :: no_ctf
        real, allocatable :: volume(:,:,:)
        complex     :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex     :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: cross_sum, particle, model
        real(dp)    :: truth_rotation(3,3), candidate_rotation(3,3)
        real(dp)    :: truth_shift(2), candidate_shift(2), scales(3)
        real(dp)    :: objective, oracle, baseline_objective
        real(dp)    :: gradient(5), baseline_gradient(5)
        real(dp)    :: particle_power, prediction_power
        integer     :: h, i, k, shell
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        no_ctf%ctfflag     = CTFFLAG_NO
        truth_rotation     = real(euler2m([19., 37., 28.]), dp)
        truth_shift        = [0.31_dp, -0.24_dp]
        candidate_rotation = real(euler2m([20., 36.2, 28.7]), dp)
        candidate_shift    = [-0.08_dp, 0.06_dp]
        call calc%predict(1, .true., truth_rotation, truth_shift, observed)
        call calc%predict(1, .true., candidate_rotation, candidate_shift, prediction)
        scales = [0.5_dp, 1._dp, 2._dp]
        do i = 1, size(scales)
            call calc%set_ptcl(1, real(scales(i), kind=kind(observed))*observed, no_ctf, [2, 4])
            call assert_int(OBJFUN_CC, calc%get_ptcl_objfun(1), 'a slot prepared without sigma2 is not a cc slot')
            call calc%objective_gradient(1, .true., 1, candidate_rotation, candidate_shift, objective, gradient)
            cross_sum        = cmplx(0._dp, 0._dp, kind=dp)
            particle_power   = 0._dp
            prediction_power = 0._dp
            do k = -TEST_BOX/2, TEST_BOX/2
                do h = -TEST_BOX/2, TEST_BOX/2
                    shell = nint(sqrt(real(h*h + k*k)))
                    if( shell < 2 .or. shell > 4 ) cycle
                    particle = cmplx(scales(i)*observed(h,k), kind=dp)
                    model    = cmplx(prediction(h,k), kind=dp)
                    cross_sum        = cross_sum + conjg(particle)*model
                    particle_power   = particle_power + real(conjg(particle)*particle, dp)
                    prediction_power = prediction_power + real(conjg(model)*model, dp)
                end do
            end do
            oracle = 1._dp - real(cross_sum, dp)/sqrt(particle_power*prediction_power)
            call assert_true(abs(objective - oracle) < FORMULA_TOL, 'Cartesian cc disagrees with the independent formula')
            call assert_true(abs(calc%score(1, objective) - (1._dp - oracle)) < FORMULA_TOL, 'cc score is not the correlation')
            if( i == 1 )then
                baseline_objective = objective
                baseline_gradient  = gradient
            else
                call assert_true(abs(objective - baseline_objective) < 1.e-12_dp .and. &
                    &maxval(abs(gradient - baseline_gradient)) < 1.e-10_dp, &
                    &'Cartesian cc changed under positive particle-amplitude scaling')
            endif
        end do
        call calc%kill
    end subroutine test_cc_objective_formula

    ! N4: the euclid loss L = sum |X - M|^2/sigma2 / sum |X|^2/sigma2 over the shells of the
    ! range, against a brute-force sum with a shell-dependent sigma2; the score exp(-L); L is
    ! invariant to a uniform scaling of sigma2 and changes under a shell-dependent scaling
    subroutine test_euclid_objective_formula()
        type(cartft_calc) :: calc
        type(ctfparams)   :: no_ctf
        real, allocatable :: volume(:,:,:)
        complex     :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex     :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: particle, model
        real        :: sigma2(0:TEST_BOX/2)
        real(dp)    :: truth_rotation(3,3), candidate_rotation(3,3), truth_shift(2), candidate_shift(2)
        real(dp)    :: objective, scaled_objective, reshaped_objective, oracle, gradient(5)
        real(dp)    :: residual_sum, particle_sum
        integer     :: h, k, shell
        integer, parameter :: kfromto(2) = [2, 6]
        ! single-precision whitened samples: a uniform rescaling changes the rounding only
        real(dp), parameter :: SCALING_TOL = 1.e-5_dp
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        no_ctf%ctfflag     = CTFFLAG_NO
        truth_rotation     = real(euler2m([19., 37., 28.]), dp)
        truth_shift        = [0.31_dp, -0.24_dp]
        candidate_rotation = real(euler2m([20., 36.2, 28.7]), dp)
        candidate_shift    = [-0.08_dp, 0.06_dp]
        call calc%predict(1, .true., truth_rotation, truth_shift, observed)
        call calc%predict(1, .true., candidate_rotation, candidate_shift, prediction)
        sigma2 = [(0.5 + 0.1*real(shell), shell = 0, TEST_BOX/2)]
        call calc%set_ptcl(1, observed, no_ctf, sigma2, kfromto)
        call assert_int(OBJFUN_EUCLID, calc%get_ptcl_objfun(1), 'a slot prepared with sigma2 is not a euclid slot')
        call calc%objective_gradient(1, .true., 1, candidate_rotation, candidate_shift, objective, gradient)
        residual_sum = 0._dp
        particle_sum = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell < kfromto(1) .or. shell > kfromto(2) ) cycle
                particle = cmplx(observed(h,k), kind=dp)
                model    = cmplx(prediction(h,k), kind=dp)
                residual_sum = residual_sum + abs(particle - model)**2/real(sigma2(shell), dp)
                particle_sum = particle_sum + abs(particle)**2/real(sigma2(shell), dp)
            end do
        end do
        oracle = residual_sum/particle_sum
        call assert_true(oracle > 1.e-3_dp, 'euclid fixture has no residual')
        call assert_true(abs(objective - oracle) < FORMULA_TOL*max(1._dp, oracle), &
            &'Cartesian euclid loss disagrees with the brute-force sum')
        call assert_true(abs(calc%score(1, objective) - exp(-oracle)) < FORMULA_TOL, 'euclid score is not exp(-L)')
        call calc%set_ptcl(1, observed, no_ctf, 3.7*sigma2, kfromto)
        call calc%objective_gradient(1, .true., 1, candidate_rotation, candidate_shift, scaled_objective, gradient)
        call assert_true(abs(scaled_objective - objective) < SCALING_TOL*objective, &
            &'euclid loss changed under a uniform scaling of sigma2')
        call calc%set_ptcl(1, observed, no_ctf, [(sigma2(shell)*(1. + 0.5*real(shell)), shell = 0, TEST_BOX/2)], kfromto)
        call calc%objective_gradient(1, .true., 1, candidate_rotation, candidate_shift, reshaped_objective, gradient)
        call assert_true(abs(reshaped_objective - objective) > 1.e-3_dp*objective, &
            &'euclid loss did not change under a shell-dependent scaling of sigma2')
        call calc%kill
    end subroutine test_euclid_objective_formula

    ! N5 (extends E9): the analytic five-parameter gradients of both objectives against centred
    ! differences, without CTF, with CTF (CTFFLAG_YES, the flipped observation) and with a
    ! phase-flipped stack (CTFFLAG_FLIP)
    subroutine test_five_parameter_gradient()
        type(cartft_calc) :: calc
        type(ctfparams)   :: ctfparms
        type(ctf)         :: tfun
        type(ctfvars)     :: ctfvals
        real, allocatable :: volume(:,:,:)
        complex  :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex  :: prediction(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real     :: sigma2(0:TEST_BOX/2), angle
        real(dp) :: rotation(3,3), truth_rotation(3,3), minus_rotation(3,3), plus_rotation(3,3)
        real(dp) :: shift(2), truth_shift(2), minus_shift(2), plus_shift(2)
        real(dp) :: objective, objective_minus, objective_plus, gradient(5), errors(5), basis(3), step
        integer  :: axis, iobj, ictf, h, k, shell
        integer, parameter :: objfuns(2) = [OBJFUN_CC, OBJFUN_EUCLID], ctfflags(3) = [CTFFLAG_NO, CTFFLAG_YES, CTFFLAG_FLIP]
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        rotation       = real(euler2m([17., 31., 23.]), dp)
        truth_rotation = right_increment_rotation(rotation, [0.01_dp, -0.008_dp, 0.006_dp])
        shift          = [0.21_dp, -0.16_dp]
        truth_shift    = [0.25_dp, -0.18_dp]
        sigma2         = [(1. + 0.1*real(shell), shell = 0, TEST_BOX/2)]
        call calc%predict(1, .true., truth_rotation, truth_shift, prediction)
        ctfparms%smpd    = 1.
        ctfparms%kv      = 300.
        ctfparms%cs      = 2.7
        ctfparms%fraca   = 0.1
        ctfparms%dfx     = 1.4
        ctfparms%dfy     = 1.65
        ctfparms%angast  = 23.
        ctfparms%phshift = 0.
        tfun = ctf(ctfparms%smpd, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
        call tfun%init(ctfparms%dfx, ctfparms%dfy, ctfparms%angast)
        ctfvals = tfun%get_ctfvars(ctfparms%phshift)
        do ictf = 1, size(ctfflags)
            ctfparms%ctfflag = ctfflags(ictf)
            observed = prediction
            if( ctfparms%ctfflag /= CTFFLAG_NO )then
                do k = -TEST_BOX/2, TEST_BOX/2
                    do h = -TEST_BOX/2, TEST_BOX/2
                        angle = 0.
                        if( h /= 0 .or. k /= 0 ) angle = atan2(real(k), real(h))
                        observed(h,k) = abs(tfun%eval_canonical(real(h*h + k*k)/real(TEST_BOX*TEST_BOX), angle, &
                            &ctfvals%phshift))*prediction(h,k)
                    end do
                end do
            endif
            do iobj = 1, size(objfuns)
                if( objfuns(iobj) == OBJFUN_CC )then
                    call calc%set_ptcl(1, observed, ctfparms, [2, 6])
                else
                    call calc%set_ptcl(1, observed, ctfparms, sigma2, [2, 6])
                endif
                call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
                do axis = 1, 5
                    basis          = 0._dp
                    minus_rotation = rotation
                    plus_rotation  = rotation
                    minus_shift    = shift
                    plus_shift     = shift
                    if( axis <= 3 )then
                        step = 1.e-4_dp
                        basis(axis)    = 1._dp
                        minus_rotation = right_increment_rotation(rotation, -step*basis)
                        plus_rotation  = right_increment_rotation(rotation,  step*basis)
                    else
                        step = 1.e-3_dp
                        minus_shift(axis-3) = minus_shift(axis-3) - step
                        plus_shift(axis-3)  = plus_shift(axis-3)  + step
                    endif
                    call calc%objective_gradient(1, .true., 1, minus_rotation, minus_shift, objective_minus, gradient)
                    call calc%objective_gradient(1, .true., 1, plus_rotation,  plus_shift,  objective_plus,  gradient)
                    call calc%objective_gradient(1, .true., 1, rotation, shift, objective, gradient)
                    errors(axis) = abs((objective_plus - objective_minus)/(2._dp*step) - gradient(axis))/&
                        &max(1.e-5_dp, abs(gradient(axis)))
                end do
                call assert_true(all(errors < GRADIENT_TOL), 'a five-parameter objective-gradient component failed centred differences')
            end do
        end do
        call calc%kill
    end subroutine test_five_parameter_gradient

    ! E10: the Cartesian gather against SIMPLE's projector kernel at identical rotated 3D
    ! coordinates from the same physical reference (permanent guard for C10)
    subroutine test_matched_projector_boundary()
        type(cartft_calc) :: calc
        type(image)       :: volume_image
        type(projector)   :: pftc_projector
        real, allocatable :: volume(:,:,:)
        complex     :: cartesian(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex     :: pftc(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex(dp) :: cross_sum
        real(dp)    :: rotation(3,3), cartesian_power, pftc_power
        real(dp)    :: correlation, gain, relative_l2, residual_power
        real(sp)    :: loc(3)
        integer     :: h, k, radius_squared
        call build_test_volume(volume)
        call new_test_calc(calc, volume)
        call volume_image%new([TEST_BOX, TEST_BOX, TEST_BOX], 1.)
        call volume_image%set_rmat(volume, .false.)
        call pftc_projector%new([OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX, OSMPL_PAD_FAC*TEST_BOX], 1.)
        call volume_image%pad_fft(pftc_projector)
        call pftc_projector%expand_cmat()
        rotation = real(euler2m([17., 31., 23.]), dp)
        call calc%predict(1, .true., rotation, [0._dp, 0._dp], cartesian)
        pftc            = cmplx(0., 0.)
        cross_sum       = cmplx(0._dp, 0._dp, kind=dp)
        cartesian_power = 0._dp
        pftc_power      = 0._dp
        residual_power  = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                radius_squared = h*h + k*k
                if( radius_squared < 4 .or. radius_squared > (TEST_BOX/2 - 1)**2 ) cycle
                loc       = real(matmul(real([h, k, 0], dp), rotation), sp)
                pftc(h,k) = pftc_projector%interp_fcomp_oversamp(loc)
                cross_sum = cross_sum + conjg(cmplx(pftc(h,k), kind=dp))*cmplx(cartesian(h,k), kind=dp)
                pftc_power      = pftc_power + abs(cmplx(pftc(h,k), kind=dp))**2
                cartesian_power = cartesian_power + abs(cmplx(cartesian(h,k), kind=dp))**2
                residual_power  = residual_power + abs(cmplx(cartesian(h,k) - pftc(h,k), kind=dp))**2
            end do
        end do
        correlation = real(cross_sum, dp)/sqrt(pftc_power*cartesian_power)
        gain        = real(cross_sum, dp)/pftc_power
        relative_l2 = sqrt(residual_power/pftc_power)
        write(*,'(a,3(1x,es12.4))') 'CARTFT_MATCHED_PROJECTOR', correlation, gain, relative_l2
        call assert_true(correlation >= 1._dp - 1.e-6_dp .and. abs(gain - 1._dp) <= 1.e-5_dp .and. relative_l2 <= 2.e-5_dp, &
            &'PFTC projector kernel and Cartesian gather disagree at a matched boundary')
        call calc%kill
        call pftc_projector%kill_expanded()
        call pftc_projector%kill()
        call volume_image%kill()
    end subroutine test_matched_projector_boundary

    ! E19 (lifecycle part): references of both halves, invalid states, rebuild and kill
    subroutine test_reference_lifecycle()
        type(cartft_calc) :: calc
        real, allocatable :: even_volume(:,:,:), odd_volume(:,:,:)
        call build_small_volume(even_volume)
        odd_volume = -0.5*even_volume
        call calc%new(1, SMALL_BOX, 1)
        call assert_true(.not. calc%ref_exists(1, .true.) .and. .not. calc%ref_exists(1, .false.), &
            &'a new calculator reported a reference before set_ref')
        call calc%set_ref(1, .true.,  even_volume)
        call calc%set_ref(1, .false., odd_volume)
        call assert_true(calc%ref_exists(1, .true.) .and. calc%ref_exists(1, .false.), &
            &'calculator did not hold both half-set references')
        call assert_true(.not. calc%ref_exists(0, .true.) .and. .not. calc%ref_exists(2, .false.), &
            &'calculator accepted an invalid state')
        ! rebuilding replaces the previous references
        call calc%new(1, SMALL_BOX, 1)
        call assert_true(.not. calc%ref_exists(1, .true.), 'rebuilt calculator kept a previous reference')
        call calc%set_ref(1, .true.,  even_volume)
        call calc%set_ref(1, .false., odd_volume)
        call assert_true(calc%ref_exists(1, .true.) .and. calc%ref_exists(1, .false.), 'calculator could not be rebuilt')
        call calc%kill
        call assert_true(.not. calc%ref_exists(1, .true.) .and. .not. calc%ref_exists(1, .false.), &
            &'killed calculator remained ready')
        call assert_true(.not. calc%ptcl_is_valid(1), 'killed calculator kept a particle slot')
        call calc%kill
    end subroutine test_reference_lifecycle

    ! N1: new/kill repeated; even and odd references from different volumes give different
    ! predictions for the same pose. Expected value: the odd volume is -0.5 times the even one,
    ! so by linearity of the transform and the gather its prediction is -0.5 times the even one.
    subroutine test_even_odd_references()
        type(cartft_calc) :: calc
        real, allocatable :: even_volume(:,:,:), odd_volume(:,:,:)
        complex  :: even_pred(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex  :: odd_pred(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), shift(2)
        real     :: scale
        integer  :: icycle
        call build_test_volume(even_volume)
        odd_volume = -0.5*even_volume
        rotation   = real(euler2m([17., 31., 23.]), dp)
        shift      = [0.23_dp, -0.17_dp]
        do icycle = 1, 3
            call calc%new(2, TEST_BOX, 2)
            call calc%set_ref(2, .true.,  even_volume)
            call calc%set_ref(2, .false., odd_volume)
            call assert_true(calc%ref_exists(2, .true.) .and. calc%ref_exists(2, .false.) .and. &
                &.not. calc%ref_exists(1, .true.) .and. .not. calc%ref_exists(1, .false.), &
                &'state references were not held as set')
            call calc%predict(2, .true.,  rotation, shift, even_pred)
            call calc%predict(2, .false., rotation, shift, odd_pred)
            scale = maxval(abs(even_pred))
            call assert_true(scale > 1.e-3, 'even prediction is empty')
            call assert_true(maxval(abs(odd_pred - even_pred)) > 0.5*scale, &
                &'even and odd references gave the same prediction')
            call assert_true(maxval(abs(odd_pred + 0.5*even_pred)) <= LINEARITY_TOL*scale, &
                &'odd prediction is not the gather of the odd volume')
            call calc%kill
            call assert_true(.not. calc%ref_exists(2, .true.), 'killed calculator remained ready')
        end do
    end subroutine test_even_odd_references

    ! E18: the sigma contribution of an accepted transaction differs from the seed's; that of a
    ! rolled-back transaction reproduces the seed's exactly (a euclid slot: only euclid has sigma2)
    subroutine test_sigma_endpoint_contracts()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        real, allocatable :: volume(:,:,:), terminal_sigma(:), seed_sigma(:), rollback_sigma(:), ref_pow(:), ptcl_pow(:)
        complex  :: observed(-SMALL_BOX/2:SMALL_BOX/2,-SMALL_BOX/2:SMALL_BOX/2)
        real(dp) :: truth_rotation(3,3), truth_shift(2), seed_rotation(3,3), seed_shift(2), rotation(3,3), shift(2)
        real     :: v
        call build_small_volume(volume)
        call new_test_calc(calc, volume)
        truth_rotation = real(euler2m([19., 37., 28.]), dp)
        truth_shift    = [0.31_dp, -0.24_dp]
        call calc%predict(1, .true., truth_rotation, truth_shift, observed)
        call set_unweighted_particle(calc, observed, [2, SMALL_BOX/2 - 1])
        seed_rotation = real(euler2m([20., 36.2, 28.7]), dp)
        seed_shift    = [-0.08_dp, 0.06_dp]
        rotation = seed_rotation
        shift    = seed_shift
        call opt%new(SMALL_BOX, SMALL_BOX, 5., 15.)
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call assert_int(CARTFT_ACCEPTED, opt%get_status(), 'sigma endpoint fixture did not produce an accepted transaction')
        call calc%sigma_contribution(1, .true., 1, rotation, shift, terminal_sigma, ref_pow, ptcl_pow, v)
        call calc%sigma_contribution(1, .true., 1, seed_rotation, seed_shift, seed_sigma, ref_pow, ptcl_pow, v)
        call assert_int(size(seed_sigma), size(terminal_sigma), 'accepted terminal sigma has the wrong shell count')
        call assert_true(all(terminal_sigma >= 0.), 'accepted terminal sigma contains a negative contribution')
        call assert_true(any(abs(terminal_sigma - seed_sigma) > 10.*epsilon(1.)), &
            &'accepted and seed poses produced the same sigma contribution')
        rotation = truth_rotation
        shift    = seed_shift
        call opt%new(SMALL_BOX, SMALL_BOX, 1.e-30, 15.)
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call assert_int(CARTFT_BOUND_REJECTED, opt%get_status(), 'sigma rollback fixture did not produce a rejected transaction')
        call calc%sigma_contribution(1, .true., 1, rotation, shift, rollback_sigma, ref_pow, ptcl_pow, v)
        call calc%sigma_contribution(1, .true., 1, truth_rotation, seed_shift, seed_sigma, ref_pow, ptcl_pow, v)
        call assert_int(size(seed_sigma), size(rollback_sigma), 'rollback sigma has the wrong shell count')
        call assert_true(all(rollback_sigma == seed_sigma), 'rollback pose did not reproduce the seed-pose sigma contribution')
        call opt%kill
        call calc%kill
    end subroutine test_sigma_endpoint_contracts

    ! N6: for a noise-free particle placed on a grid pose (an in-plane rotation of the polar grid
    ! at one projection direction), the Cartesian objective evaluated over the grid poses peaks
    ! where the polar objective peaks, for both objectives (gen_objfun_vals; for cc also
    ! gen_corr_for_rot_8, the ring weighting of the continuous polar stages). Reference and
    ! particle are the same central sections on both sides; the polar in-plane index irot maps
    ! to e3 = 360 - rot(irot), as assign_ori stores it. Expected value: the placed index.
    subroutine test_polar_grid_pose_agreement()
        ! parameters%new requires a box above 26
        integer, parameter :: POLAR_BOX = 32, PLACED_INDEX = 7, KFROMTO(2) = [2, 10]
        real,    parameter :: E1 = 37., E2 = 63.
        type(cartft_calc)  :: calc
        type(polarft_calc) :: pftc
        type(parameters), target :: p
        type(cmdline)      :: cline
        type(ctfparams)    :: no_ctf
        type(image)        :: ref_img, ptcl_img
        real, allocatable  :: volume(:,:,:), vals(:), sigma2_noise(:,:)
        complex  :: ref_plane(-POLAR_BOX/2:POLAR_BOX/2,-POLAR_BOX/2:POLAR_BOX/2)
        complex  :: ptcl_plane(-POLAR_BOX/2:POLAR_BOX/2,-POLAR_BOX/2:POLAR_BOX/2)
        real     :: sigma2(0:POLAR_BOX/2)
        real(dp) :: objective, gradient(5)
        real(dp), allocatable :: cart_objectives(:), polar_corrs(:)
        integer  :: iobj, irot, nrots, pdim(3)
        integer, parameter :: objfuns(2) = [OBJFUN_CC, OBJFUN_EUCLID]
        character(len=6), parameter :: objfun_names(2) = ['cc    ', 'euclid']
        call build_test_volume(volume, POLAR_BOX)
        call new_test_calc(calc, volume)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = 1.
        do iobj = 1, size(objfuns)
            call cline%set('box',     real(POLAR_BOX))
            call cline%set('smpd',    1.0)
            call cline%set('mskdiam', 20.0)
            call cline%set('nptcls',  1.0)
            call cline%set('nthr',    1.0)
            call cline%set('ctf',     'no')
            call cline%set('objfun',  trim(objfun_names(iobj)))
            call p%new(cline, silent=.true.)
            call pftc%new(p, 1, [1,1], KFROMTO)
            nrots = pftc%get_nrots()
            pdim  = pftc%get_pdim_srch()
            ! the reference at e3 = 0, the particle at the grid rotation PLACED_INDEX
            call calc%predict(1, .true., real(euler2m([E1, E2, 0.]), dp), [0._dp, 0._dp], ref_plane)
            call calc%predict(1, .true., grid_rotation(PLACED_INDEX), [0._dp, 0._dp], ptcl_plane)
            call plane_to_image(ref_plane, ref_img)
            call plane_to_image(ptcl_plane, ptcl_img)
            call ref_img%memoize4polarize(pdim)
            call pftc%polarize_ref_pft(ref_img, 1, iseven=.true., pdim=pdim, oversamp=.false.)
            call pftc%polarize_ptcl_pft(ptcl_img, 1, pdim=pdim, oversamp=.false.)
            call pftc%set_eo(1, .true.)
            call pftc%memoize_refs
            call pftc%memoize_ptcls
            if( objfuns(iobj) == OBJFUN_EUCLID )then
                allocate(sigma2_noise(KFROMTO(1):KFROMTO(2),1), source=1.)
                call pftc%assign_sigma2_noise(sigma2_noise)
                call pftc%memoize_sqsum_ptcl(1)
                call calc%set_ptcl(1, ptcl_plane, no_ctf, sigma2, KFROMTO)
            else
                call calc%set_ptcl(1, ptcl_plane, no_ctf, KFROMTO)
            endif
            allocate(vals(nrots), cart_objectives(nrots))
            call pftc%gen_objfun_vals(1, 1, [0., 0.], vals)
            do irot = 1, nrots
                call calc%objective_gradient(1, .true., 1, grid_rotation(irot), [0._dp, 0._dp], objective, gradient)
                cart_objectives(irot) = objective
            end do
            call assert_int(PLACED_INDEX, maxloc(vals, dim=1), 'polar objective does not peak at the placed grid pose')
            call assert_int(maxloc(vals, dim=1), minloc(cart_objectives, dim=1), &
                &'Cartesian and polar '//trim(objfun_names(iobj))//' select different grid poses')
            if( objfuns(iobj) == OBJFUN_CC )then
                allocate(polar_corrs(nrots))
                do irot = 1, nrots
                    polar_corrs(irot) = pftc%gen_corr_for_rot_8(1, 1, irot)
                end do
                call assert_int(minloc(cart_objectives, dim=1), maxloc(polar_corrs, dim=1), &
                    &'Cartesian cc and the ring-weighted polar cc select different grid poses')
                deallocate(polar_corrs)
            endif
            deallocate(vals, cart_objectives)
            if( allocated(sigma2_noise) ) deallocate(sigma2_noise)
            call pftc%kill
            call ref_img%kill
            call ptcl_img%kill
            call cline%kill
        end do
        call calc%kill

      contains

        function grid_rotation( irot ) result( rotation )
            integer, intent(in) :: irot
            real(dp) :: rotation(3,3)
            rotation = real(euler2m([E1, E2, 360. - pftc%get_rot(irot)]), dp)
        end function grid_rotation

    end subroutine test_polar_grid_pose_agreement

    !> A Fourier-space image whose expand_ft is the Hermitian full-disk plane.
    subroutine plane_to_image( plane, img )
        complex,     intent(in)    :: plane(:,:)
        type(image), intent(inout) :: img
        integer :: h, k, box
        box = size(plane,1) - 1
        call img%new([box, box, 1], 1.0)
        call img%zero_and_flag_ft()
        do k = -box/2, box/2 - 1
            do h = 0, box/2
                call img%set_fcomp([h, k, 0], img%comp_addr_phys(h, k, 0), plane(h + box/2 + 1, k + box/2 + 1))
            end do
        end do
    end subroutine plane_to_image

    ! N28: one particle with a stored shift and CTFFLAG_YES. The prepared observation equals
    ! prepimg4align's image before its padding (O4): the same fused noise normalization, clip,
    ! shift by minus the stored shift and phase flip, then the mask; ifft_mask_pad_fft, the step
    ! prepimg4align ends with, masks identically, its padded image cropping to that image. With
    ! O5 (a) the observation is that image times taper(i)*taper(j), and the taper equals the
    ! response of polarize_oversamp's normalized stencil at an on-grid sample: the polar samples
    ! of the zero-padded masked image at 0 and 90 degrees equal the tapered observation.
    subroutine test_observation_preparation()
        integer, parameter :: NATIVE_BOX = 2*SMALL_BOX, PFTSZ = 8, KFROMTO(2) = [2, SMALL_BOX/2 - 1]
        real,    parameter :: CROP_SMPD = 1.5, MSKRAD = 6., NATIVE_SHIFT(2) = [1.5, -2.25]
        type(cartft_calc) :: calc
        type(image)       :: raw, oracle, work, oracle_work, padded_work, masked, padded
        type(ctfparams)   :: ctf_native, ctf_crop, ctf_out
        type(ctf)         :: tfun
        logical, allocatable :: noise_mask(:,:,:)
        complex, allocatable :: observed(:,:), untapered(:,:), tapered(:,:)
        real,    allocatable :: taper(:), masked_rmat(:,:,:), padded_rmat(:,:,:), tapered_rmat(:,:,:), ones(:)
        complex  :: pft(PFTSZ,KFROMTO(1):KFROMTO(2))
        real     :: crop_shift(2), scale, shvec(2)
        integer  :: i, j, k, offset
        call calc%new(1, SMALL_BOX, 1)
        taper = calc%get_ptcl_taper()
        call assert_int(SMALL_BOX, size(taper), 'taper does not span the box')
        call assert_true(abs(taper(SMALL_BOX/2 + 1) - 1.) < 1.e-6 .and. all(taper > 0.) .and. all(taper <= 1. + 1.e-6) .and. &
            &all(abs(taper(2:SMALL_BOX) - taper(SMALL_BOX:2:-1)) < 1.e-6), 'taper is not unit-centred, positive and symmetric')
        ! a smooth native particle with a stored shift and CTF
        call raw%new([NATIVE_BOX, NATIVE_BOX, 1], CROP_SMPD*real(SMALL_BOX)/real(NATIVE_BOX))
        do j = 1, NATIVE_BOX
            do i = 1, NATIVE_BOX
                call raw%set_rmat_at(i, j, 1, sin(0.21*real(i)) + cos(0.16*real(j)) + 0.01*real(i*j) + &
                    &2.*exp(-real((i-14)**2 + (j-19)**2)/18.))
            end do
        end do
        call oracle%copy(raw)
        call masked%disc([NATIVE_BOX, NATIVE_BOX, 1], raw%get_smpd(), 0.4*real(NATIVE_BOX), noise_mask)
        call masked%kill
        ctf_native%smpd    = raw%get_smpd()
        ctf_native%kv      = 300.
        ctf_native%cs      = 2.7
        ctf_native%fraca   = 0.1
        ctf_native%dfx     = 1.2
        ctf_native%dfy     = 1.35
        ctf_native%angast  = 31.
        ctf_native%phshift = 0.
        ctf_native%ctfflag = CTFFLAG_YES
        crop_shift = NATIVE_SHIFT*real(SMALL_BOX)/real(NATIVE_BOX)
        call work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        call oracle_work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        call padded_work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        call padded%new([OSMPL_PAD_FAC*SMALL_BOX, OSMPL_PAD_FAC*SMALL_BOX, 1], CROP_SMPD)
        call work%memoize_mask_coords()
        call oracle_work%memoize_mask_coords()
        call padded_work%memoize_mask_coords()
        ! prepimg4align's steps: shvec = -stored shift on the cropped grid, CTF at smpd_crop;
        ! the fused phase flip reads the memoized Fourier maps of the cropped box
        call memoize_ft_maps([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        shvec    = -crop_shift
        ctf_crop = ctf_native
        ctf_crop%smpd = CROP_SMPD
        tfun = ctf(ctf_crop%smpd, ctf_crop%kv, ctf_crop%cs, ctf_crop%fraca)
        call oracle%norm_noise_fft_clip_shift_ctf_flip(noise_mask, oracle_work, shvec, tfun, ctf_crop)
        call padded_work%copy(oracle_work)
        call padded_work%ifft_mask_pad_fft(MSKRAD, padded)
        call oracle_work%ifft_mask_fft(MSKRAD)
        untapered = oracle_work%expand_ft()
        ! the padded image of prepimg4align crops to the unpadded masked image
        call oracle_work%ifft()
        call padded%ifft()
        masked_rmat = oracle_work%get_rmat()
        padded_rmat = padded%get_rmat()
        offset = (OSMPL_PAD_FAC*SMALL_BOX - SMALL_BOX)/2
        scale  = maxval(abs(masked_rmat))
        call assert_true(maxval(abs(padded_rmat(offset+1:offset+SMALL_BOX,offset+1:offset+SMALL_BOX,1:1) - masked_rmat)) &
            &<= 1.e-5*scale, 'prepimg4align padded image does not crop to the masked image')
        ! the Cartesian observation without and with the taper
        allocate(ones(SMALL_BOX), source=1.)
        call prepimg4align_cart(raw, noise_mask, work, MSKRAD, CROP_SMPD, crop_shift, ones, ctf_native, &
            &observed, ctf_out)
        scale = maxval(abs(untapered))
        call assert_true(maxval(abs(observed - untapered)) <= 2.e-5*scale, &
            &'observation differs from prepimg4align before its padding')
        call assert_true(abs(ctf_out%smpd - CROP_SMPD) <= epsilon(CROP_SMPD) .and. ctf_out%ctfflag == CTFFLAG_YES, &
            &'observation did not carry the cropped CTF parameters')
        ! the fused preparation normalizes its input in place: rebuild the raw particle
        call raw%new([NATIVE_BOX, NATIVE_BOX, 1], CROP_SMPD*real(SMALL_BOX)/real(NATIVE_BOX))
        do j = 1, NATIVE_BOX
            do i = 1, NATIVE_BOX
                call raw%set_rmat_at(i, j, 1, sin(0.21*real(i)) + cos(0.16*real(j)) + 0.01*real(i*j) + &
                    &2.*exp(-real((i-14)**2 + (j-19)**2)/18.))
            end do
        end do
        call prepimg4align_cart(raw, noise_mask, work, MSKRAD, CROP_SMPD, crop_shift, taper, ctf_native, &
            &observed, ctf_out)
        tapered_rmat = masked_rmat
        do j = 1, SMALL_BOX
            do i = 1, SMALL_BOX
                tapered_rmat(i,j,1) = masked_rmat(i,j,1)*taper(i)*taper(j)
            end do
        end do
        call oracle_work%set_rmat(tapered_rmat, .false.)
        call oracle_work%fft()
        tapered = oracle_work%expand_ft()
        call assert_true(maxval(abs(observed - tapered)) <= 2.e-5*scale, &
            &'observation is not the masked image times the stencil taper')
        call forget_ft_maps
        ! the polar on-grid samples of the zero-padded, untapered masked image (polarize_oversamp)
        call oracle_work%set_rmat(masked_rmat, .false.)
        call padded%new([OSMPL_PAD_FAC*SMALL_BOX, OSMPL_PAD_FAC*SMALL_BOX, 1], CROP_SMPD)
        call oracle_work%pad(padded, backgr=0.)
        call padded%fft()
        call padded%memoize4polarize_oversamp([PFTSZ, KFROMTO(1), KFROMTO(2)])
        call padded%polarize_oversamp(pft)
        do k = KFROMTO(1), KFROMTO(2)
            ! angle index 1 samples (0,-k), index PFTSZ/2+1 samples (k,0), both on the native lattice
            call assert_true(abs(pft(1,k) - observed(0,-k)) <= 1.e-3*scale .and. &
                &abs(pft(PFTSZ/2 + 1,k) - observed(k,0)) <= 1.e-3*scale, &
                &'the taper is not the response of the polar stencil at an on-grid sample')
        end do
        call calc%kill
        call raw%kill
        call oracle%kill
        call work%kill
        call oracle_work%kill
        call padded_work%kill
        call padded%kill
    end subroutine test_observation_preparation

    !> Two states of distinct volumes: state 1 the blobs (even) and -0.5 times them (odd),
    !! state 2 the blobs mirrored along x (even) and 1.5 times those (odd).
    subroutine build_two_state_volumes( volumes )
        real, allocatable, intent(out) :: volumes(:,:,:,:,:) !< (box,box,box,half,state)
        real, allocatable :: volume(:,:,:)
        call build_test_volume(volume)
        allocate(volumes(TEST_BOX,TEST_BOX,TEST_BOX,2,2))
        volumes(:,:,:,1,1) = volume
        volumes(:,:,:,2,1) = -0.5*volume
        volumes(:,:,:,1,2) = volume(TEST_BOX:1:-1,:,:)
        volumes(:,:,:,2,2) = 1.5*volume(TEST_BOX:1:-1,:,:)
    end subroutine build_two_state_volumes

    !> True when every state and half of the two calculators gives the same prediction bit for bit.
    logical function same_predictions( calc_a, calc_b, nstates )
        type(cartft_calc), intent(in) :: calc_a, calc_b
        integer,           intent(in) :: nstates
        complex  :: pred_a(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex  :: pred_b(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3)
        integer  :: state, ihalf
        same_predictions = .true.
        rotation = real(euler2m([17., 31., 23.]), dp)
        do state = 1, nstates
            do ihalf = 1, 2
                call calc_a%predict(state, ihalf == 1, rotation, [0.23_dp, -0.17_dp], pred_a)
                call calc_b%predict(state, ihalf == 1, rotation, [0.23_dp, -0.17_dp], pred_b)
                same_predictions = same_predictions .and. all(pred_a == pred_b)
            end do
        end do
    end function same_predictions

    ! N8 (file): the reference-volume files of section 6.6 written from staged volumes and read
    ! back give predictions identical, bit for bit, to references built from the same volumes,
    ! for every state and half; the header carries format version, box, state count, kfrom, kto
    ! and sampling, and the reader returns that band; the compatibility function is false for
    ! a wrong version, box, state count, sampling, an out-of-range band, and even and odd
    ! headers that differ. Expected values: the written arrays; hand-built headers.
    subroutine test_reference_volume_file()
        character(len=*), parameter :: FEVEN = 'tmp_cartft_calc_tester_refvols_even.bin'
        character(len=*), parameter :: FODD  = 'tmp_cartft_calc_tester_refvols_odd.bin'
        real,    parameter :: SMPD = 1.3
        integer, parameter :: KFROMTO(2) = [2, 9]
        type(cartft_calc) :: writer, reader, direct
        type(string)      :: fname_even, fname_odd
        real, allocatable :: volumes(:,:,:,:,:)
        character(len=32) :: field
        integer :: header(5), good(5), kfromto_read(2), state, ihalf, funit
        real    :: smpd_read
        call build_two_state_volumes(volumes)
        fname_even = FEVEN
        fname_odd  = FODD
        call writer%new(2, TEST_BOX, 1)
        call direct%new(2, TEST_BOX, 1)
        do state = 1, 2
            do ihalf = 1, 2
                call writer%set_refvol(state, ihalf == 1, volumes(:,:,:,ihalf,state))
                call direct%set_ref(state, ihalf == 1, volumes(:,:,:,ihalf,state))
            end do
        end do
        call writer%write(fname_even, fname_odd, KFROMTO, SMPD)
        call writer%kill
        ! the header as written
        open(newunit=funit, file=FEVEN, access='stream', action='read', status='old')
        read(funit, pos=1) header, smpd_read
        close(funit)
        call assert_true(all(header == [CART_REFVOLS_FORMAT_VERSION, TEST_BOX, 2, KFROMTO(1), KFROMTO(2)]) .and. &
            &smpd_read == SMPD, 'reference-volume header does not carry version, box, nstates, kfrom, kto and sampling')
        call reader%read(fname_even, fname_odd, 2, TEST_BOX, 1, SMPD, kfromto_read)
        call assert_true(all(kfromto_read == KFROMTO), 'reader did not return the band limit of the header')
        call assert_true(same_predictions(reader, direct, 2), 'read references differ from references of the same volumes')
        ! the compatibility rule, on hand-built headers against the reader's box and state count
        good = [CART_REFVOLS_FORMAT_VERSION, TEST_BOX, 2, 2, 9]
        call assert_true(reader%cart_refvols_header_compatible(good, SMPD, good, SMPD, SMPD, field) .and. field == '', &
            &'a compatible header was rejected')
        header = good
        header(1) = CART_REFVOLS_FORMAT_VERSION + 1
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), 'wrong version accepted')
        call assert_char('format_version', trim(field), 'wrong version not named')
        header = good
        header(2) = TEST_BOX + 2
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), 'wrong box accepted')
        header = good
        header(3) = 1
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), 'wrong nstates accepted')
        call assert_false(reader%cart_refvols_header_compatible(good, 1.31, good, 1.31, SMPD, field), 'wrong sampling accepted')
        call assert_char('smpd_crop', trim(field), 'wrong sampling not named')
        header = good
        header(4) = 0
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), 'kfrom below 1 accepted')
        header = good
        header(4:5) = [5, 4]
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), 'kto below kfrom accepted')
        header = good
        header(5) = TEST_BOX/2 + 1
        call assert_false(reader%cart_refvols_header_compatible(header, SMPD, header, SMPD, SMPD, field), &
            &'kto beyond fdim(box)-1 accepted')
        call assert_char('kfrom/kto', trim(field), 'out-of-range band not named')
        header = good
        header(5) = 8
        call assert_false(reader%cart_refvols_header_compatible(good, SMPD, header, SMPD, SMPD, field), &
            &'different even and odd headers accepted')
        call assert_false(reader%cart_refvols_header_compatible(good, SMPD, good, 1.31, SMPD, field), &
            &'different even and odd sampling accepted')
        call reader%kill
        call direct%kill
        call del_file(FEVEN)
        call del_file(FODD)
        call assert_true(.not. file_exists(FEVEN) .and. .not. file_exists(FODD), 'test reference-volume files survived')
    end subroutine test_reference_volume_file

    ! N8 (reader) and N9: through the matcher's reader of a Cartesian pass
    ! (read_reprojection_model with l_cart_refine), the files named by
    ! refine3D_cart_refvols_fname build build%cftc with predictions identical to references
    ! of the same volumes, and the reader adopts kfromto and lp of the header; then
    ! remove_ref_section_files removes both files and no file of the test survives.
    subroutine test_reference_volume_handoff()
        real,    parameter :: SMPD = 1.3
        integer, parameter :: KFROMTO(2) = [3, 10]
        type(cartft_calc) :: writer, direct
        type(builder)     :: build
        type(parameters)  :: params
        type(string)      :: fname_even, fname_odd
        real, allocatable :: volumes(:,:,:,:,:)
        integer :: state, ihalf
        call build_two_state_volumes(volumes)
        fname_even = refine3D_cart_refvols_fname('even')
        fname_odd  = refine3D_cart_refvols_fname('odd')
        call writer%new(2, TEST_BOX, 1)
        call direct%new(2, TEST_BOX, 1)
        do state = 1, 2
            do ihalf = 1, 2
                call writer%set_refvol(state, ihalf == 1, volumes(:,:,:,ihalf,state))
                call direct%set_ref(state, ihalf == 1, volumes(:,:,:,ihalf,state))
            end do
        end do
        call writer%write(fname_even, fname_odd, KFROMTO, SMPD)
        call writer%kill
        params%l_cart_refine = .true.
        params%nstates   = 2
        params%nthr      = 1
        params%box       = TEST_BOX
        params%box_crop  = TEST_BOX
        params%smpd      = SMPD
        params%smpd_crop = SMPD
        params%kfromto   = [2, 11]
        call read_reprojection_model(params, build, 1)
        call assert_true(all(params%kfromto == KFROMTO), 'the reader did not adopt kfromto of the header')
        call assert_real(real(TEST_BOX)*SMPD/real(KFROMTO(2)), params%lp, 1.e-4, 'the reader did not adopt lp of the header')
        call assert_true(same_predictions(build%cftc, direct, 2), 'build%cftc differs from references of the same volumes')
        call build%cftc%kill
        call direct%kill
        call remove_ref_section_files
        call assert_true(.not. file_exists(fname_even) .and. .not. file_exists(fname_odd), &
            &'remove_ref_section_files left a Cartesian reference-volume file')
    end subroutine test_reference_volume_handoff

    !> A smooth raw particle of box n with a Gaussian blob, and optional added noise.
    subroutine build_raw_particle( img, n, smpd, seed_offset )
        type(image), intent(inout) :: img
        integer,     intent(in)    :: n
        real,        intent(in)    :: smpd
        integer,     intent(in)    :: seed_offset
        integer :: i, j
        call img%new([n, n, 1], smpd)
        do j = 1, n
            do i = 1, n
                call img%set_rmat_at(i, j, 1, sin(0.21*real(i+seed_offset)) + cos(0.16*real(j)) + 0.01*real(i*j) + &
                    &2.*exp(-real((i-n/2+seed_offset)**2 + (j-n/2-2)**2)/18.))
            end do
        end do
    end subroutine build_raw_particle

    ! E13 (observation part, moved from the adapter tester in Phase 6): without the taper (a
    ! unit taper) and at a zero stored shift the Cartesian observation of prepimg4align_cart is
    ! the established particle path (norm_noise_fft_clip_shift, ifft_mask_fft) on the full
    ! redundant disk, with the cropped sampling in its CTF parameters
    subroutine test_observation_established_path()
        integer, parameter :: NATIVE_BOX = 32
        real,    parameter :: CROP_SMPD = 1.5
        type(image)        :: raw, oracle, work, oracle_work
        type(ctfparams)    :: input_ctf, output_ctf
        logical, allocatable :: noise_mask(:,:,:)
        complex, allocatable :: observed(:,:), expected(:,:)
        real,    allocatable :: unit_taper(:)
        call build_raw_particle(raw, NATIVE_BOX, CROP_SMPD*real(SMALL_BOX)/real(NATIVE_BOX), 0)
        call oracle%copy(raw)
        call work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        call oracle_work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
        call work%memoize_mask_coords()
        allocate(noise_mask(NATIVE_BOX, NATIVE_BOX, 1), source=.false.)
        input_ctf%smpd    = raw%get_smpd()
        input_ctf%ctfflag = CTFFLAG_NO
        allocate(unit_taper(SMALL_BOX), source=1.)
        call prepimg4align_cart(raw, noise_mask, work, 6., CROP_SMPD, [0., 0.], unit_taper, input_ctf, observed, output_ctf)
        call oracle%norm_noise_fft_clip_shift(noise_mask, oracle_work, [0., 0.])
        call oracle_work%ifft_mask_fft(6.)
        expected = oracle_work%expand_ft()
        call assert_true(maxval(abs(observed - expected)) <= 2.e-5, 'observation disagrees with the established particle path')
        call assert_true(all(lbound(observed) == [-SMALL_BOX/2, -SMALL_BOX/2]) .and. &
            &all(ubound(observed) == [SMALL_BOX/2, SMALL_BOX/2]), 'observation is not the full redundant disk')
        call assert_true(abs(output_ctf%smpd - CROP_SMPD) <= epsilon(CROP_SMPD), 'observation did not carry the cropped sampling')
        call raw%kill
        call oracle%kill
        call work%kill
        call oracle_work%kill
    end subroutine test_observation_established_path

    !> The inputs of a small CTF batch of NB particles: raw images of box 2*SMALL_BOX, stored
    !! shifts on the cropped grid, CTF parameters and a per-particle sigma2.
    subroutine build_batch_fixture( nb, raws, crop_shifts, ctfs, sigma2, noise_mask )
        integer,                      intent(in)  :: nb
        type(image),     allocatable, intent(out) :: raws(:)
        real,            allocatable, intent(out) :: crop_shifts(:,:), sigma2(:,:)
        type(ctfparams), allocatable, intent(out) :: ctfs(:)
        logical,         allocatable, intent(out) :: noise_mask(:,:,:)
        type(image) :: masker
        integer :: i, shell
        allocate(raws(nb), crop_shifts(2,nb), ctfs(nb), sigma2(0:SMALL_BOX/2,nb))
        do i = 1, nb
            call build_raw_particle(raws(i), 2*SMALL_BOX, 0.75, i)
            crop_shifts(:,i) = [0.4*real(i) - 1., 0.7 - 0.3*real(i)]
            ctfs(i)%smpd    = 0.75
            ctfs(i)%kv      = 300.
            ctfs(i)%cs      = 2.7
            ctfs(i)%fraca   = 0.1
            ctfs(i)%dfx     = 1.0 + 0.1*real(i)
            ctfs(i)%dfy     = 1.1 + 0.1*real(i)
            ctfs(i)%angast  = 10.*real(i)
            ctfs(i)%phshift = 0.
            ctfs(i)%ctfflag = CTFFLAG_YES
            sigma2(:,i) = [(1. + 0.05*real(shell*i), shell = 0, SMALL_BOX/2)]
        end do
        call masker%disc([2*SMALL_BOX, 2*SMALL_BOX, 1], 0.75, 0.4*real(2*SMALL_BOX), noise_mask)
        call masker%kill
    end subroutine build_batch_fixture

    !> Objectives and gradients of every slot at one fixed pose, for comparing preparations.
    subroutine slot_signatures( calc, nb, values )
        type(cartft_calc), intent(in)  :: calc
        integer,           intent(in)  :: nb
        real(dp),          intent(out) :: values(6,nb)
        real(dp) :: gradient(5)
        integer  :: i
        do i = 1, nb
            call calc%objective_gradient(1, .true., i, real(euler2m([20., 36.2, 28.7]), dp), [0.1_dp, -0.2_dp], &
                &values(1,i), gradient)
            values(2:6,i) = gradient
        end do
    end subroutine slot_signatures

    ! N11 and N12: the batch preparation of production (prep_cart_batch) fills slot i exactly as
    ! the single-particle path called in the test (prepimg4align_cart on a copy, then set_ptcl),
    ! under both objectives; with a team of three threads it gives the serial result bit for
    ! bit. Expected values: that path, and the serial run.
    subroutine test_batch_preparation()
        integer, parameter :: NB = 5, KFROMTO(2) = [2, SMALL_BOX/2 - 1]
        real,    parameter :: CROP_SMPD = 1.5, MSKRAD = 6.
        type(cartft_calc) :: batch_calc, single_calc, team_calc
        type(image),     allocatable :: raws(:)
        type(image)      :: raw_copy, work
        type(ctfparams), allocatable :: ctfs(:)
        type(ctfparams)  :: ctf_out
        logical,         allocatable :: noise_mask(:,:,:)
        real,            allocatable :: crop_shifts(:,:), sigma2(:,:), volume(:,:,:), taper(:)
        complex,         allocatable :: observed(:,:)
        real(dp) :: batch_values(6,NB), single_values(6,NB), team_values(6,NB)
        integer  :: i, iobj, nthr_saved
        call build_small_volume(volume)
        call build_batch_fixture(NB, raws, crop_shifts, ctfs, sigma2, noise_mask)
        nthr_saved = 1
        !$ nthr_saved = omp_get_max_threads()
        do iobj = 1, 2
            call new_batch_calc(batch_calc)
            call new_batch_calc(single_calc)
            call new_batch_calc(team_calc)
            taper = single_calc%get_ptcl_taper()
            !$ call omp_set_num_threads(1)
            if( iobj == 1 )then
                call prep_cart_batch(batch_calc, NB, raws, noise_mask, SMALL_BOX, CROP_SMPD, MSKRAD, crop_shifts, ctfs, KFROMTO)
            else
                call prep_cart_batch(batch_calc, NB, raws, noise_mask, SMALL_BOX, CROP_SMPD, MSKRAD, crop_shifts, ctfs, KFROMTO, &
                    &sigma2)
            endif
            !$ call omp_set_num_threads(3)
            if( iobj == 1 )then
                call prep_cart_batch(team_calc, NB, raws, noise_mask, SMALL_BOX, CROP_SMPD, MSKRAD, crop_shifts, ctfs, KFROMTO)
            else
                call prep_cart_batch(team_calc, NB, raws, noise_mask, SMALL_BOX, CROP_SMPD, MSKRAD, crop_shifts, ctfs, KFROMTO, &
                    &sigma2)
            endif
            !$ call omp_set_num_threads(nthr_saved)
            ! the single-particle path, serially
            call work%new([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
            call work%memoize_mask_coords()
            call memoize_ft_maps([SMALL_BOX, SMALL_BOX, 1], CROP_SMPD)
            do i = 1, NB
                call raw_copy%copy(raws(i))
                call prepimg4align_cart(raw_copy, noise_mask, work, MSKRAD, CROP_SMPD, crop_shifts(:,i), taper, ctfs(i), &
                    &observed, ctf_out)
                if( iobj == 1 )then
                    call single_calc%set_ptcl(i, observed, ctf_out, KFROMTO)
                else
                    call single_calc%set_ptcl(i, observed, ctf_out, sigma2(:,i), KFROMTO)
                endif
            end do
            call forget_ft_maps
            call slot_signatures(batch_calc,  NB, batch_values)
            call slot_signatures(single_calc, NB, single_values)
            call slot_signatures(team_calc,   NB, team_values)
            call assert_true(all(batch_values == single_values), &
                &'batch preparation differs from the single-particle preparation (N11)')
            call assert_true(all(team_values == batch_values), 'threaded batch preparation differs from the serial one (N12)')
            call assert_true(batch_calc%get_ptcl_objfun(1) == merge(OBJFUN_CC, OBJFUN_EUCLID, iobj == 1), &
                &'batch preparation chose the wrong objective')
            call batch_calc%kill
            call single_calc%kill
            call team_calc%kill
            call work%kill
            call raw_copy%kill
        end do
        do i = 1, NB
            call raws(i)%kill
        end do

      contains

        subroutine new_batch_calc( calc )
            type(cartft_calc), intent(inout) :: calc
            call calc%new(1, SMALL_BOX, 1)
            call calc%set_ref(1, .true.,  volume)
            call calc%set_ref(1, .false., volume)
            call calc%new_ptcls(NB)
        end subroutine new_batch_calc

    end subroutine test_batch_preparation

    ! N26: objective, gradient and residual evaluated concurrently on different particles by a
    ! team of three threads equal the serial results bit for bit (the calculator is intent(in)
    ! in every evaluation). Expected value: the serial results.
    subroutine test_concurrent_evaluation()
        integer, parameter :: NP = 6, KFROMTO(2) = [2, TEST_BOX/2 - 1]
        type(cartft_calc) :: calc
        type(ctfparams)   :: no_ctf
        real, allocatable :: volume(:,:,:), sigma_serial(:,:), sigma_team(:,:), contrib(:), ref_pow(:), ptcl_pow(:)
        complex  :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real     :: sigma2(0:TEST_BOX/2), v
        real(dp) :: serial(6,NP), team(6,NP), gradient(5), rotation(3,3)
        integer  :: i, shell
        call build_test_volume(volume)
        call calc%new(1, TEST_BOX, NP)
        call calc%set_ref(1, .true., volume)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = [(1. + 0.1*real(shell), shell = 0, TEST_BOX/2)]
        do i = 1, NP
            call calc%predict(1, .true., real(euler2m([10.*real(i), 30., 5.*real(i)]), dp), [0.2_dp, -0.1_dp], observed)
            if( mod(i,2) == 0 )then
                call calc%set_ptcl(i, observed, no_ctf, KFROMTO)
            else
                call calc%set_ptcl(i, observed, no_ctf, sigma2, KFROMTO)
            endif
        end do
        allocate(sigma_serial(KFROMTO(1):KFROMTO(2),NP), sigma_team(KFROMTO(1):KFROMTO(2),NP), source=0.)
        do i = 1, NP
            rotation = real(euler2m([10.*real(i) + 1., 31., 5.*real(i)]), dp)
            call calc%objective_gradient(1, .true., i, rotation, [0.1_dp, 0._dp], serial(1,i), gradient)
            serial(2:6,i) = gradient
            if( mod(i,2) /= 0 )then
                call calc%sigma_contribution(1, .true., i, rotation, [0.1_dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
                sigma_serial(:,i) = contrib
            endif
        end do
        !$omp parallel do num_threads(3) default(shared) private(i,rotation,gradient,contrib,ref_pow,ptcl_pow,v) schedule(static,1)
        do i = 1, NP
            rotation = real(euler2m([10.*real(i) + 1., 31., 5.*real(i)]), dp)
            call calc%objective_gradient(1, .true., i, rotation, [0.1_dp, 0._dp], team(1,i), gradient)
            team(2:6,i) = gradient
            if( mod(i,2) /= 0 )then
                call calc%sigma_contribution(1, .true., i, rotation, [0.1_dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
                sigma_team(:,i) = contrib
            endif
        end do
        !$omp end parallel do
        call assert_true(all(team == serial), 'concurrent objective and gradient differ from the serial results')
        call assert_true(all(sigma_team == sigma_serial), 'concurrent residuals differ from the serial results')
        call calc%kill
    end subroutine test_concurrent_evaluation

    ! N13: under euclid the sigma owner's Cartesian calc_sigma2 at a committed pose stores the
    ! per-shell residual (mean squared residual over two per sample, a brute-force sum in the
    ! test); an exact match stores zero. A Cartesian cc pass allocates, reads and writes no
    ! sigma2: prep_sigmas_objfun leaves the sigma owner unbuilt (C5).
    subroutine test_sigma_owner()
        integer, parameter :: KFROMTO(2) = [2, 6]
        type(cartft_calc)    :: calc
        type(euclid_sigma2)  :: esig
        type(parameters), target :: params
        type(builder)        :: build
        type(ctfparams)      :: no_ctf
        type(string)         :: fname
        real, allocatable    :: volume(:,:,:), stored(:)
        complex  :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        complex  :: model(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real     :: sigma2(0:TEST_BOX/2)
        real(dp) :: truth(3,3), committed(3,3), oracle(KFROMTO(1):KFROMTO(2)), counts(KFROMTO(1):KFROMTO(2))
        integer  :: h, k, shell
        call build_test_volume(volume)
        call calc%new(1, TEST_BOX, 1)
        call calc%set_ref(1, .true., volume)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = 1.
        truth     = real(euler2m([19., 37., 28.]), dp)
        committed = real(euler2m([20., 36.2, 28.7]), dp)
        call calc%predict(1, .true., truth, [0.3_dp, -0.2_dp], observed)
        call calc%set_ptcl(1, observed, no_ctf, sigma2, KFROMTO)
        params%fromp    = 7
        params%top      = 7
        params%kfromto  = KFROMTO
        params%box      = TEST_BOX
        fname = 'tmp_cartft_calc_tester_sigma2.bin'
        call esig%new(params, fname, TEST_BOX)
        call esig%allocate_ptcls
        call esig%calc_sigma2(calc, 1, 7, 1, .true., committed, [0.1_dp, 0._dp])
        stored = esig%get_sigma2_part(7)
        call calc%predict(1, .true., committed, [0.1_dp, 0._dp], model)
        oracle = 0._dp
        counts = 0._dp
        do k = -TEST_BOX/2, TEST_BOX/2
            do h = -TEST_BOX/2, TEST_BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell < KFROMTO(1) .or. shell > KFROMTO(2) ) cycle
                oracle(shell) = oracle(shell) + abs(cmplx(observed(h,k) - model(h,k), kind=dp))**2
                counts(shell) = counts(shell) + 1._dp
            end do
        end do
        oracle = oracle/(2._dp*counts)
        call assert_true(all(abs(real(stored, dp) - oracle) <= 1.e-5_dp*maxval(oracle)), &
            &'Cartesian calc_sigma2 is not the per-shell residual at the committed pose')
        call esig%calc_sigma2(calc, 1, 7, 1, .true., truth, [0.3_dp, -0.2_dp])
        stored = esig%get_sigma2_part(7)
        call assert_true(all(abs(stored) <= 1.e-6*real(maxval(oracle))), 'an exact match stored a nonzero sigma2')
        call esig%kill
        ! a Cartesian cc pass: no sigma2 state at all
        params%l_cart_refine = .true.
        params%cc_objfun     = OBJFUN_CC
        params%cc_emit_sigma = 'no'
        call prep_sigmas_objfun(params, build)
        call assert_false(allocated(build%esig%sigma2_noise), 'a cc pass allocated sigma2')
        call calc%kill
    end subroutine test_sigma_owner

    ! N14 (with the particle-side taper of O5 (a)): for a noise-only particle at a grid pose
    ! (a zero reference, so the residual is the particle), the per-shell sigma2 from the
    ! Cartesian residual (prepimg4align_cart with the taper, cartft_calc%sigma_contribution)
    ! agrees with the polar one (prepimg4align's steps with the padding, polarize_oversamp,
    ! gen_sigma_contrib) within the standard error of a shell; without the taper their ratio
    ! is the taper's noise-power factor in every shell. Criterion declared before the first run
    ! (journal_phase_6.md): mean_s |d_s|/SE_s <= 1, max_s <= 3, SE_s = sqrt(2/n_s).
    subroutine test_polar_sigma_agreement()
        ! the box, sampling, mask and band of N6, so the polar calculator reuses its FFT plans
        ! (the fast tier's budget)
        integer, parameter :: BOX = 32, KFROMTO(2) = [2, 10]
        real,    parameter :: SMPD = 1.0, MSKRAD = 10.
        type(cartft_calc)  :: calc
        type(polarft_calc) :: pftc
        type(parameters), target :: p
        type(cmdline)      :: cline
        type(image)        :: noise, polar_work, polar_pad, cart_work, zero_pad, masked
        type(ctfparams)    :: no_ctf, ctf_out
        logical, allocatable :: noise_mask(:,:,:)
        complex, allocatable :: observed(:,:)
        real,    allocatable :: zero_volume(:,:,:), taper(:), unit_taper(:), contrib(:), ref_pow(:), ptcl_pow(:)
        real,    allocatable :: mrmat(:,:,:)
        real     :: polar_sigma(KFROMTO(1):KFROMTO(2)), cart_tapered(KFROMTO(1):KFROMTO(2)), cart_plain(KFROMTO(1):KFROMTO(2))
        real     :: sigma2(0:BOX/2), v
        real(dp) :: se(KFROMTO(1):KFROMTO(2)), d(KFROMTO(1):KFROMTO(2)), factor, num, den
        integer  :: pdim(3), h, k, shell, i, j, n_s(KFROMTO(1):KFROMTO(2))
        ! the polar calculator from parameters
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', 2.*MSKRAD*SMPD)
        call cline%set('nptcls',  1.0)
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline, silent=.true.)
        call set_fixed_seed(20261001)
        call pftc%new(p, 1, [1,1], KFROMTO)
        pdim = pftc%get_pdim_interp()
        ! one noise-only particle
        call noise%new([BOX, BOX, 1], SMPD)
        call noise%gauran(0., 1.)
        call noise%memoize_mask_coords
        call masked%disc([BOX, BOX, 1], SMPD, MSKRAD, noise_mask)
        call masked%kill
        no_ctf%smpd    = SMPD
        no_ctf%ctfflag = CTFFLAG_NO
        ! polar: prepimg4align's steps (noise normalization, mask, pad by 2), polarize_oversamp
        call polar_work%new([BOX, BOX, 1], SMPD)
        call polar_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, 1], SMPD)
        call zero_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, 1], SMPD)
        call masked%copy(noise)
        call masked%norm_noise_fft_clip_shift(noise_mask, polar_work, [0., 0.])
        call polar_work%ifft_mask_pad_fft(MSKRAD, polar_pad)
        call polar_pad%memoize4polarize_oversamp(pdim)
        call pftc%polarize_ptcl_pft(polar_pad, 1, pdim=pdim, oversamp=.true.)
        call zero_pad%zero_and_flag_ft()
        call pftc%polarize_ref_pft(zero_pad, 1, iseven=.true., pdim=pdim, oversamp=.true.)
        call pftc%set_eo(1, .true.)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call pftc%gen_sigma_contrib(1, 1, [0., 0.], 1, sigma_contrib=polar_sigma)
        ! Cartesian: the same particle with and without the taper, against a zero reference
        allocate(zero_volume(BOX,BOX,BOX), source=0.)
        call calc%new(1, BOX, 1)
        call calc%set_ref(1, .true., zero_volume)
        taper = calc%get_ptcl_taper()
        allocate(unit_taper(BOX), source=1.)
        sigma2 = 1.
        call cart_work%new([BOX, BOX, 1], SMPD)
        call masked%copy(noise)
        call prepimg4align_cart(masked, noise_mask, cart_work, MSKRAD, SMPD, [0., 0.], taper, no_ctf, observed, ctf_out)
        call calc%set_ptcl(1, observed, ctf_out, sigma2, KFROMTO)
        call calc%sigma_contribution(1, .true., 1, real(euler2m([0., 0., 0.]), dp), [0._dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
        cart_tapered = contrib
        call masked%copy(noise)
        call prepimg4align_cart(masked, noise_mask, cart_work, MSKRAD, SMPD, [0., 0.], unit_taper, no_ctf, observed, ctf_out)
        call calc%set_ptcl(1, observed, ctf_out, sigma2, KFROMTO)
        call calc%sigma_contribution(1, .true., 1, real(euler2m([0., 0., 0.]), dp), [0._dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
        cart_plain = contrib
        ! the taper's noise-power factor of this realization (Parseval over the masked image)
        call cart_work%ifft()
        mrmat = cart_work%get_rmat()
        num = 0._dp
        den = 0._dp
        do j = 1, BOX
            do i = 1, BOX
                num = num + (real(mrmat(i,j,1), dp)*taper(i)*taper(j))**2
                den = den + real(mrmat(i,j,1), dp)**2
            end do
        end do
        factor = num/den
        ! standard error of a shell mean of n_s pixels, n_s/2 independent (Hermitian pairs)
        n_s = 0
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell >= KFROMTO(1) .and. shell <= KFROMTO(2) ) n_s(shell) = n_s(shell) + 1
            end do
        end do
        se = sqrt(2._dp/real(n_s, dp))
        d  = (real(cart_tapered, dp) - real(polar_sigma, dp))/real(cart_tapered, dp)
        write(*,'(a,2(1x,f8.4),a,f8.4)') 'CARTFT_N14 mean/max |d|/SE with taper:', sum(abs(d)/se)/real(size(d)), &
            &maxval(abs(d)/se), '; taper noise-power factor', factor
        call assert_true(sum(abs(d)/se)/real(size(d)) <= 1._dp .and. maxval(abs(d)/se) <= 3._dp, &
            &'Cartesian and polar per-shell sigma2 disagree beyond the standard error (N14)')
        d = real(polar_sigma, dp)/real(cart_plain, dp) - factor
        write(*,'(a,2(1x,f8.4))') 'CARTFT_N14 untapered polar/Cartesian ratio min/max:', &
            &minval(real(polar_sigma, dp)/real(cart_plain, dp)), maxval(real(polar_sigma, dp)/real(cart_plain, dp))
        call assert_true(all(abs(d)/factor <= 3._dp*se), &
            &'without the taper the polar/Cartesian sigma2 ratio is not the taper noise-power factor (N14)')
        call pftc%kill
        call calc%kill
        call noise%kill
        call masked%kill
        call polar_work%kill
        call polar_pad%kill
        call zero_pad%kill
        call cart_work%kill
        call cline%kill
    end subroutine test_polar_sigma_agreement

    ! N30: a polish pass starts its sigma2 rows from the range the discrete pass wrote into the open
    ! transaction (prep_sigmas_objfun with l_cont_polish), overlays the residual of a particle whose
    ! Cartesian slot is valid and keeps the discrete row of one whose slot is invalid; the candidate
    ! the transaction commits holds no row of the committed generation.
    subroutine test_polish_sigma_fallback()
        integer, parameter :: NP = 3, KFROMTO(2) = [2, 6], NSHELL = TEST_BOX/2
        real,    parameter :: COMMITTED_VALUE = 5.
        character(len=*), parameter :: COMMITTED = 'tmp_cartft_polish_sigma2_state.bin'
        type(cartft_calc)         :: calc
        type(parameters), target  :: params
        type(builder),    target  :: build
        type(ctfparams)           :: no_ctf
        type(string)              :: candidate, range_path
        real, allocatable         :: volume(:,:,:), polished(:)
        real(real32), allocatable :: rows(:,:)
        real(real32) :: discrete(NSHELL,NP)
        complex      :: observed(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real         :: sigma2(0:TEST_BOX/2)
        real(dp)     :: truth(3,3), committed_pose(3,3)
        integer      :: i, status
        logical      :: outside_band
        character(len=128) :: message
        call new_sigma2_fixture(COMMITTED, NP, TEST_BOX, COMMITTED_VALUE, params, build)
        ! the transaction the discrete pass opened and the residuals it wrote for part 1
        candidate = sigma2_state_candidate_path(COMMITTED, 2_int64)
        range_path = sigma2_state_range_path(COMMITTED, 2_int64, 1, 1)
        call sigma2_state_prepare_update(COMMITTED, candidate%to_char(), status, message)
        call assert_int(0, status, 'open the sigma2 transaction: '//trim(message))
        do i = 1, NP
            discrete(:,i) = real(7 + i, real32)
        end do
        call sigma2_state_write_local_range(range_path%to_char(), 2_int64, FIXTURE_DIGEST, 1, discrete, 1, NSHELL, status, message)
        call assert_int(0, status, 'write the discrete range: '//trim(message))
        ! the polish: slot 1 prepared for particle 2; slot 2, never prepared, for particle 3
        params%kfromto       = KFROMTO
        params%l_cont_polish = .true.
        call prep_sigmas_objfun(params, build)
        call build_test_volume(volume)
        call calc%new(1, TEST_BOX, 2)
        call calc%set_ref(1, .true., volume)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2         = 1.
        truth          = real(euler2m([19., 37., 28.]), dp)
        committed_pose = real(euler2m([20., 36.2, 28.7]), dp)
        call calc%predict(1, .true., truth, [0.3_dp, -0.2_dp], observed)
        call calc%set_ptcl(1, observed, no_ctf, sigma2, KFROMTO)
        call build%esig%calc_sigma2(calc, 1, 2, 1, .true., committed_pose, [0.1_dp, 0._dp])
        call build%esig%calc_sigma2(calc, 2, 3, 1, .true., committed_pose, [0._dp, 0._dp])
        polished = build%esig%get_sigma2_part(2)
        call build%esig%write_sigma2
        call sigma2_state_merge_local_ranges(candidate%to_char(), [range_path], [(.true., i = 1, NP)], status, message)
        call assert_int(0, status, 'merge the polish range: '//trim(message))
        call sigma2_state_read_particles(candidate%to_char(), 1, NP, rows, status, message)
        call assert_int(0, status, 'read the candidate rows: '//trim(message))
        call assert_true(all(rows(:,1) == discrete(:,1)), 'the polish changed the row of a particle it did not evaluate')
        call assert_true(all(abs(rows(KFROMTO(1):KFROMTO(2),2) - polished) <= 1.e-6*maxval(polished)), &
            &'the polish row is not the Cartesian residual inside the band')
        outside_band = all(rows(:KFROMTO(1)-1,2) == discrete(:KFROMTO(1)-1,2)) .and. &
            &all(rows(KFROMTO(2)+1:,2) == discrete(KFROMTO(2)+1:,2))
        call assert_true(outside_band, 'the polish row outside the band is not the discrete residual')
        call assert_true(all(rows(:,3) == discrete(:,3)), 'an invalid Cartesian slot did not keep the discrete residual')
        call assert_true(all(rows /= real(COMMITTED_VALUE, real32)), 'the candidate holds a row of the committed generation')
        call build%esig%kill
        call calc%kill
        call build%spproj%kill
        call remove_sigma2_files(COMMITTED, 2)
    end subroutine test_polish_sigma_fallback

    ! N31: model amplitude and CTF handling, polar against Cartesian. One volume gives both
    ! references through their production projections; a noise-free particle with an astigmatic CTF
    ! goes through both production preparations; both evaluate at the same grid pose. Gated: the
    ! per-shell power of the CTF-modulated model, ring at radius k against pixels with nint(r) = k
    ! (from a model of the fixture over 300 random poses: band total within 3.6 %, shells within
    ! 0.80-1.25): band total within 5 %, every shell within [0.75, 1.3]. The residual of a noise-free
    ! particle is all structure, which the two samplings need not agree on; it is printed only (N33).
    subroutine test_polar_sigma_signal()
        integer,  parameter :: BOX = 48, KFROMTO(2) = [3, 10], PLACED_INDEX = 9
        real,     parameter :: SMPD = 2.0, MSKRAD = 14., NOISE_RADIUS = 22., E1 = 37., E2 = 63.
        real,     parameter :: PSI_OFFSET = 1.4, TRUE_SHIFT(2) = [0.43, -0.71], DFX = 0.15, DFY = 0.18, ANGAST = 31.
        real(dp), parameter :: TOTAL_TOL = 0.05_dp, SHELL_RANGE(2) = [0.75_dp, 1.3_dp]
        type(cartft_calc)  :: calc
        type(polarft_calc) :: pftc
        type(parameters), target :: p
        type(sp_project), target :: spproj
        type(cmdline)      :: cline
        type(oris)         :: eulspace
        type(projector)    :: vol_pad
        type(image)        :: vol_img, ptcl, noise, masked, polar_work, polar_pad, cart_work
        type(ctf)          :: tfun
        type(ctfparams)    :: ctfparms, ctf_out
        logical, allocatable :: noise_mask(:,:,:)
        complex, allocatable :: observed(:,:)
        real,    allocatable :: volume(:,:,:), taper(:), contrib(:), ref_pow(:), ptcl_pow(:), prmat(:,:,:), nrmat(:,:,:)
        complex  :: plane(-BOX/2:BOX/2,-BOX/2:BOX/2)
        real     :: polar_sigma(KFROMTO(1):KFROMTO(2)), sigma2(0:BOX/2), v, step, centre, e3
        real     :: polar_ref_pow(KFROMTO(1):KFROMTO(2)), polar_ptcl_pow(KFROMTO(1):KFROMTO(2))
        real(dp) :: ratio(KFROMTO(1):KFROMTO(2)), n_s(KFROMTO(1):KFROMTO(2)), truth(3,3), evaluated(3,3), total
        real     :: noise_amp
        integer  :: pdim(3), h, k, i, j, shell, n_noise, n_mask
        ! the polar calculator; its CTF matrices come from a one-particle project below
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', 2.*MSKRAD*SMPD)
        call cline%set('nptcls',  1.0)
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline, silent=.true.)
        p%nspace = 1
        call set_fixed_seed(20261002)
        call pftc%new(p, 1, [1,1], KFROMTO)
        call pftc%set_with_ctf(.true.)
        pdim      = pftc%get_pdim_interp()
        step      = 360./real(pftc%get_nrots())
        e3        = 360. - pftc%get_rot(PLACED_INDEX)
        evaluated = real(euler2m([E1, E2, e3]), dp)
        truth     = real(euler2m([E1, E2, e3 - PSI_OFFSET*step]), dp)
        ! both references from one volume, each through its production projection
        call build_test_volume(volume, BOX)
        call vol_img%new([BOX, BOX, BOX], SMPD)
        call vol_img%set_rmat(volume, .false.)
        call vol_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX], SMPD)
        call vol_img%pad_fft(vol_pad)
        call vol_pad%expand_cmat()
        call eulspace%new(1, is_ptcl=.false.)
        call eulspace%set_euler(1, [E1, E2, 0.])
        call vol_pad2ref_pfts_opt(pftc, vol_pad, eulspace, 1, .true.)
        call calc%new(1, BOX, 1)
        call calc%set_ref(1, .true., volume)
        ! the particle: the projection at the true pose times the CTF, and noise beyond the soft mask
        ! scaled to unit deviation over the normalization's background (r > MSKRAD), so the residual
        ! is the pose mismatch rather than a gain
        call calc%predict(1, .true., truth, real(TRUE_SHIFT, dp), plane)
        call ptcl%new([BOX, BOX, 1], SMPD)
        call ptcl%zero_and_flag_ft()
        do k = -BOX/2, BOX/2 - 1
            do h = 0, BOX/2
                call ptcl%set_fcomp([h, k, 0], ptcl%comp_addr_phys(h, k, 0), plane(h, k))
            end do
        end do
        call memoize_ft_maps([BOX, BOX, 1], SMPD)
        ctfparms%smpd    = SMPD
        ctfparms%kv      = 300.
        ctfparms%cs      = 2.7
        ctfparms%fraca   = 0.1
        ctfparms%dfx     = DFX
        ctfparms%dfy     = DFY
        ctfparms%angast  = ANGAST
        ctfparms%phshift = 0.
        ctfparms%ctfflag = CTFFLAG_YES
        tfun = ctf(SMPD, ctfparms%kv, ctfparms%cs, ctfparms%fraca)
        call ptcl%apply_ctf(tfun, 'ctf', ctfparms)
        call ptcl%ifft()
        call noise%new([BOX, BOX, 1], SMPD)
        call noise%gauran(0., 1.)
        prmat   = ptcl%get_rmat()
        nrmat   = noise%get_rmat()
        centre  = real(BOX/2 + 1)
        n_noise = 0
        n_mask  = 0
        do j = 1, BOX
            do i = 1, BOX
                if( sqrt((real(i) - centre)**2 + (real(j) - centre)**2) > MSKRAD       ) n_mask  = n_mask  + 1
                if( sqrt((real(i) - centre)**2 + (real(j) - centre)**2) > NOISE_RADIUS ) n_noise = n_noise + 1
            end do
        end do
        noise_amp = sqrt(real(n_mask)/real(n_noise))
        do j = 1, BOX
            do i = 1, BOX
                if( sqrt((real(i) - centre)**2 + (real(j) - centre)**2) > NOISE_RADIUS ) &
                    &prmat(i,j,1) = prmat(i,j,1) + noise_amp*nrmat(i,j,1)
            end do
        end do
        call ptcl%set_rmat(prmat, .false.)
        call ptcl%memoize_mask_coords
        call masked%disc([BOX, BOX, 1], SMPD, MSKRAD, noise_mask)
        call masked%kill
        ! polar: prepimg4align's steps (noise normalization, phase flip, mask, pad by 2), the CTF
        ! matrices of a one-particle project, the residual at grid rotation PLACED_INDEX
        call polar_work%new([BOX, BOX, 1], SMPD)
        call polar_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, 1], SMPD)
        call polar_work%memoize_mask_coords
        call masked%copy(ptcl)
        call masked%norm_noise_fft_clip_shift_ctf_flip(noise_mask, polar_work, [0., 0.], tfun, ctfparms)
        call polar_work%ifft_mask_pad_fft(MSKRAD, polar_pad)
        call polar_pad%memoize4polarize_oversamp(pdim)
        call pftc%polarize_ptcl_pft(polar_pad, 1, pdim=pdim, oversamp=.true.)
        call spproj%os_stk%new(1, is_ptcl=.false.)
        call spproj%os_stk%set(1, 'stk',   'n31_particles.mrcs')
        call spproj%os_stk%set(1, 'fromp', 1)
        call spproj%os_stk%set(1, 'top',   1)
        call spproj%os_stk%set(1, 'nptcls_stk', 1)
        call spproj%os_stk%set(1, 'box',   BOX)
        call spproj%os_stk%set(1, 'smpd',  SMPD)
        call spproj%os_stk%set(1, 'ctf',   'yes')
        call spproj%os_stk%set(1, 'kv',    ctfparms%kv)
        call spproj%os_stk%set(1, 'cs',    ctfparms%cs)
        call spproj%os_stk%set(1, 'fraca', ctfparms%fraca)
        call spproj%os_ptcl3D%new(1, is_ptcl=.true.)
        call spproj%os_ptcl3D%set(1, 'stkind', 1)
        call spproj%os_ptcl3D%set(1, 'indstk', 1)
        call spproj%os_ptcl3D%set(1, 'dfx',    DFX)
        call spproj%os_ptcl3D%set(1, 'dfy',    DFY)
        call spproj%os_ptcl3D%set(1, 'angast', ANGAST)
        call spproj%os_ptcl3D%set_state(1, 1)
        call pftc%create_polar_absctfmats(spproj, 'ptcl3D')
        call pftc%set_eo(1, .true.)
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call pftc%gen_sigma_contrib(1, 1, [0., 0.], PLACED_INDEX, sigma_contrib=polar_sigma, &
            &ref_pow=polar_ref_pow, ptcl_pow=polar_ptcl_pow)
        ! Cartesian: prepimg4align_cart with the taper, the residual at the same pose
        taper = calc%get_ptcl_taper()
        call cart_work%new([BOX, BOX, 1], SMPD)
        call cart_work%memoize_mask_coords
        call masked%copy(ptcl)
        call prepimg4align_cart(masked, noise_mask, cart_work, MSKRAD, SMPD, [0., 0.], taper, ctfparms, observed, ctf_out)
        sigma2 = 1.
        call calc%set_ptcl(1, observed, ctf_out, sigma2, KFROMTO)
        call calc%sigma_contribution(1, .true., 1, evaluated, [0._dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
        call forget_ft_maps
        call assert_true(all(ref_pow > 0.) .and. all(contrib > 0.), 'N31 fixture: model or residual vanishes in a shell')
        ! the Cartesian estimator's own shell weights: its pixel counts
        n_s = 0._dp
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell >= KFROMTO(1) .and. shell <= KFROMTO(2) ) n_s(shell) = n_s(shell) + 1._dp
            end do
        end do
        ratio = real(polar_ref_pow, dp)/real(ref_pow, dp)
        total = sum(n_s*real(polar_ref_pow, dp))/sum(n_s*real(ref_pow, dp))
        write(*,'(a,*(1x,f7.4))') 'CARTFT_N31 polar/Cartesian model power per shell:', ratio
        write(*,'(a,f8.4)') 'CARTFT_N31 model power band total ratio', total
        write(*,'(a,*(1x,f7.4))') 'CARTFT_N31 (diagnostic) polar/Cartesian observation power per shell:', polar_ptcl_pow/ptcl_pow
        write(*,'(a,*(1x,f7.4))') 'CARTFT_N31 (diagnostic) polar/Cartesian residual per shell:', polar_sigma/contrib
        call assert_true(abs(total - 1._dp) <= TOTAL_TOL, &
            &'polar and Cartesian model power differ over the band beyond the sampling difference (N31)')
        call assert_true(all(ratio >= SHELL_RANGE(1) .and. ratio <= SHELL_RANGE(2)), &
            &'polar and Cartesian model power differ in a shell beyond the sampling difference (N31)')
        call pftc%kill
        call calc%kill
        call vol_pad%kill_expanded()
        call vol_pad%kill()
        call vol_img%kill
        call ptcl%kill
        call noise%kill
        call masked%kill
        call polar_work%kill
        call polar_pad%kill
        call cart_work%kill
        call eulspace%kill
        call spproj%kill
        call cline%kill
    end subroutine test_polar_sigma_signal

    ! N33: the group sigma2 that mixed polar/Cartesian generations reduce agrees between the two
    ! representations under realistic conditions: NG particles at in-plane angles off the polar grid
    ! (true shift zero, the shift a search would find) and astigmatic CTFs, white noise over the
    ! whole image at a signal-to-noise ratio of SNR, a mask clear of the particle; both production
    ! paths, both evaluated at the nearest-below grid pose. Criterion of N14 for the group mean:
    ! mean_s |d_s|/SE_s <= 1 and max_s <= 3, d_s relative, SE_s = sqrt(2/(n_s*NG)). Both sides see
    ! the same noise, so d_s is mostly systematic and the gate bounds it at about SE_s (4-9 %).
    subroutine test_polar_sigma_group()
        integer, parameter :: BOX = 48, KFROMTO(2) = [3, 10], PLACED_INDEX = 9, NG = 16
        real,    parameter :: SMPD = 2.0, MSKRAD = 18., SNR = 0.1
        type(cartft_calc)  :: calc
        type(polarft_calc) :: pftc
        type(parameters), target :: p
        type(sp_project), target :: spproj
        type(cmdline)      :: cline
        type(oris)         :: eulspace
        type(projector)    :: vol_pad
        type(image)        :: vol_img, ptcl, noise, masked, polar_work, polar_pad, cart_work
        type(ctf)          :: tfun
        type(ctfparams)    :: ctfparms(NG), ctf_out
        logical, allocatable :: noise_mask(:,:,:)
        complex, allocatable :: observed(:,:)
        real,    allocatable :: volume(:,:,:), taper(:), contrib(:), ref_pow(:), ptcl_pow(:), prmat(:,:,:), nrmat(:,:,:)
        complex  :: plane(-BOX/2:BOX/2,-BOX/2:BOX/2)
        real     :: polar_sigma(KFROMTO(1):KFROMTO(2)), sigma2(0:BOX/2), v, step, e3, gain, mean_sig, var_sig
        real     :: e1(NG), e2(NG), psi_off(NG)
        real(dp) :: polar_sum(KFROMTO(1):KFROMTO(2)), cart_sum(KFROMTO(1):KFROMTO(2)), n_s(KFROMTO(1):KFROMTO(2))
        real(dp) :: d(KFROMTO(1):KFROMTO(2)), se(KFROMTO(1):KFROMTO(2)), truth(3,3), evaluated(3,3)
        integer  :: pdim(3), h, k, i, shell
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    SMPD)
        call cline%set('mskdiam', 2.*MSKRAD*SMPD)
        call cline%set('nptcls',  real(NG))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline, silent=.true.)
        p%nspace = NG
        call set_fixed_seed(20261003)
        call pftc%new(p, NG, [1,NG], KFROMTO)
        call pftc%set_with_ctf(.true.)
        pdim = pftc%get_pdim_interp()
        step = 360./real(pftc%get_nrots())
        e3   = 360. - pftc%get_rot(PLACED_INDEX)
        ! the group: deterministic, well spread poses, in-plane offsets and defoci
        do i = 1, NG
            e1(i)       = 360.*modulo(0.6180340*real(i), 1.)
            e2(i)       = 20. + 140.*modulo(0.4142136*real(i), 1.)
            psi_off(i)  = 0.3 + 1.2*modulo(0.7320508*real(i), 1.)
            ctfparms(i)%smpd    = SMPD
            ctfparms(i)%kv      = 300.
            ctfparms(i)%cs      = 2.7
            ctfparms(i)%fraca   = 0.1
            ctfparms(i)%dfx     = 0.15 + 0.10*modulo(0.3141593*real(i), 1.)
            ctfparms(i)%dfy     = 1.2*ctfparms(i)%dfx
            ctfparms(i)%angast  = 180.*modulo(0.1234567*real(i), 1.)
            ctfparms(i)%phshift = 0.
            ctfparms(i)%ctfflag = CTFFLAG_YES
        end do
        tfun = ctf(SMPD, 300., 2.7, 0.1)
        call memoize_ft_maps([BOX, BOX, 1], SMPD)
        ! the signal scaled once (references and particles alike) to variance SNR over the box
        call build_test_volume(volume, BOX)
        call calc%new(1, BOX, NG)
        call calc%set_ref(1, .true., volume)
        call make_signal(1)
        prmat    = ptcl%get_rmat()
        mean_sig = sum(prmat)/real(BOX*BOX)
        var_sig  = sum((prmat - mean_sig)**2)/real(BOX*BOX)
        gain     = sqrt(SNR/var_sig)
        volume   = gain*volume
        call forget_ft_maps
        call calc%new(1, BOX, NG)
        call calc%set_ref(1, .true., volume)
        call vol_img%new([BOX, BOX, BOX], SMPD)
        call vol_img%set_rmat(volume, .false.)
        call vol_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX], SMPD)
        call vol_img%pad_fft(vol_pad)
        call vol_pad%expand_cmat()
        call eulspace%new(NG, is_ptcl=.false.)
        do i = 1, NG
            call eulspace%set_euler(i, [e1(i), e2(i), 0.])
        end do
        call vol_pad2ref_pfts_opt(pftc, vol_pad, eulspace, 1, .true.)
        ! every particle through both production preparations
        call memoize_ft_maps([BOX, BOX, 1], SMPD)
        call masked%disc([BOX, BOX, 1], SMPD, MSKRAD, noise_mask)
        call masked%kill
        call noise%new([BOX, BOX, 1], SMPD)
        call polar_work%new([BOX, BOX, 1], SMPD)
        call polar_pad%new([OSMPL_PAD_FAC*BOX, OSMPL_PAD_FAC*BOX, 1], SMPD)
        call cart_work%new([BOX, BOX, 1], SMPD)
        call polar_work%memoize_mask_coords
        call cart_work%memoize_mask_coords
        call polar_pad%memoize4polarize_oversamp(pdim)
        taper  = calc%get_ptcl_taper()
        sigma2 = 1.
        cart_sum = 0._dp
        do i = 1, NG
            call make_signal(i)
            call noise%gauran(0., 1.)
            prmat = ptcl%get_rmat()
            nrmat = noise%get_rmat()
            call ptcl%set_rmat(prmat + nrmat, .false.)
            call ptcl%memoize_mask_coords
            call masked%copy(ptcl)
            call masked%norm_noise_fft_clip_shift_ctf_flip(noise_mask, polar_work, [0., 0.], tfun, ctfparms(i))
            call polar_work%ifft_mask_pad_fft(MSKRAD, polar_pad)
            call pftc%polarize_ptcl_pft(polar_pad, i, pdim=pdim, oversamp=.true.)
            call masked%copy(ptcl)
            call prepimg4align_cart(masked, noise_mask, cart_work, MSKRAD, SMPD, [0., 0.], taper, ctfparms(i), observed, ctf_out)
            call calc%set_ptcl(i, observed, ctf_out, sigma2, KFROMTO)
            evaluated = real(euler2m([e1(i), e2(i), e3]), dp)
            call calc%sigma_contribution(1, .true., i, evaluated, [0._dp, 0._dp], contrib, ref_pow, ptcl_pow, v)
            cart_sum = cart_sum + real(contrib, dp)
        end do
        call forget_ft_maps
        ! the polar residuals, with the CTF matrices of a project of the group
        call spproj%os_stk%new(1, is_ptcl=.false.)
        call spproj%os_stk%set(1, 'stk',   'n33_particles.mrcs')
        call spproj%os_stk%set(1, 'fromp', 1)
        call spproj%os_stk%set(1, 'top',   NG)
        call spproj%os_stk%set(1, 'nptcls_stk', NG)
        call spproj%os_stk%set(1, 'box',   BOX)
        call spproj%os_stk%set(1, 'smpd',  SMPD)
        call spproj%os_stk%set(1, 'ctf',   'yes')
        call spproj%os_stk%set(1, 'kv',    300.)
        call spproj%os_stk%set(1, 'cs',    2.7)
        call spproj%os_stk%set(1, 'fraca', 0.1)
        call spproj%os_ptcl3D%new(NG, is_ptcl=.true.)
        do i = 1, NG
            call spproj%os_ptcl3D%set(i, 'stkind', 1)
            call spproj%os_ptcl3D%set(i, 'indstk', i)
            call spproj%os_ptcl3D%set(i, 'dfx',    ctfparms(i)%dfx)
            call spproj%os_ptcl3D%set(i, 'dfy',    ctfparms(i)%dfy)
            call spproj%os_ptcl3D%set(i, 'angast', ctfparms(i)%angast)
            call spproj%os_ptcl3D%set_state(i, 1)
            call pftc%set_eo(i, .true.)
        end do
        call pftc%create_polar_absctfmats(spproj, 'ptcl3D')
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        polar_sum = 0._dp
        do i = 1, NG
            call pftc%gen_sigma_contrib(i, i, [0., 0.], PLACED_INDEX, sigma_contrib=polar_sigma)
            polar_sum = polar_sum + real(polar_sigma, dp)
        end do
        call assert_true(all(cart_sum > 0._dp) .and. all(polar_sum > 0._dp), 'N33 fixture: a group residual vanishes in a shell')
        ! the group means against the standard error of a group mean
        n_s = 0._dp
        do k = -BOX/2, BOX/2
            do h = -BOX/2, BOX/2
                shell = nint(sqrt(real(h*h + k*k)))
                if( shell >= KFROMTO(1) .and. shell <= KFROMTO(2) ) n_s(shell) = n_s(shell) + 1._dp
            end do
        end do
        se = sqrt(2._dp/(n_s*real(NG, dp)))
        d  = (polar_sum - cart_sum)/cart_sum
        write(*,'(a,*(1x,f7.4))') 'CARTFT_N33 group polar/Cartesian sigma2 per shell:', polar_sum/cart_sum
        write(*,'(a,2(1x,f8.4))') 'CARTFT_N33 mean/max |d|/SE:', sum(abs(d)/se)/real(size(d)), maxval(abs(d)/se)
        call assert_true(sum(abs(d)/se)/real(size(d)) <= 1._dp .and. maxval(abs(d)/se) <= 3._dp, &
            &'polar and Cartesian group sigma2 disagree beyond the standard error (N33)')
        call pftc%kill
        call calc%kill
        call vol_pad%kill_expanded()
        call vol_pad%kill()
        call vol_img%kill
        call ptcl%kill
        call noise%kill
        call masked%kill
        call polar_work%kill
        call polar_pad%kill
        call cart_work%kill
        call eulspace%kill
        call spproj%kill
        call cline%kill

      contains

        !> ptcl: the signal of particle i, its projection at the true pose times its CTF, real space
        subroutine make_signal( i )
            integer, intent(in) :: i
            truth = real(euler2m([e1(i), e2(i), e3 - psi_off(i)*step]), dp)
            call calc%predict(1, .true., truth, [0._dp, 0._dp], plane)
            call ptcl%new([BOX, BOX, 1], SMPD)
            call ptcl%zero_and_flag_ft()
            do k = -BOX/2, BOX/2 - 1
                do h = 0, BOX/2
                    call ptcl%set_fcomp([h, k, 0], ptcl%comp_addr_phys(h, k, 0), plane(h, k))
                end do
            end do
            call ptcl%apply_ctf(tfun, 'ctf', ctfparms(i))
            call ptcl%ifft()
        end subroutine make_signal

    end subroutine test_polar_sigma_group

    ! N32: the canonical round trip of a generation a Cartesian pass writes and a polar pass reads:
    ! write (calc_sigma2, write_sigma2), merge, reduce, commit, then reload through the polar sigma
    ! owner (read_part, read_groups). The exact-zero rule: a record with zero shells (an observation
    ! empty there against a zero reference) stays in the particle rows, the reducer excludes it from
    ! the groups, and its particle reloads the positive group of its half.
    subroutine test_canonical_sigma_round_trip()
        integer, parameter :: BOX = 32, NP = 4, KFROMTO(2) = [2, 10], NSHELL = BOX/2, ZERO_FROM = 4
        real,    parameter :: SEED = 3.
        character(len=*), parameter :: COMMITTED = 'tmp_cartft_round_trip_sigma2_state.bin'
        type(cartft_calc)         :: calc
        type(parameters), target  :: p
        type(builder),    target  :: build
        type(cmdline)             :: cline
        type(ctfparams)           :: no_ctf
        type(sigma2_state_header) :: header
        type(string)              :: candidate, range_path
        real, allocatable         :: volume(:,:,:), zero_volume(:,:,:), band(:)
        real(real32), allocatable :: groups(:,:,:)
        complex  :: observed(-BOX/2:BOX/2,-BOX/2:BOX/2)
        real     :: sigma2(0:BOX/2), written(NSHELL,NP), mean_even(NSHELL)
        real(dp) :: truth(3,3), committed_pose(3,3)
        logical  :: active(NP)
        integer  :: eo(NP), group_ids(NP), i, h, k, status
        character(len=128) :: message
        call cline%set('box',     real(BOX))
        call cline%set('smpd',    1.0)
        call cline%set('mskdiam', 20.0)
        call cline%set('nptcls',  real(NP))
        call cline%set('nthr',    1.0)
        call cline%set('ctf',     'no')
        call cline%set('objfun',  'euclid')
        call p%new(cline, silent=.true.)
        call new_sigma2_fixture(COMMITTED, NP, BOX, SEED, p, build)
        p%kfromto = KFROMTO
        candidate = sigma2_state_candidate_path(COMMITTED, 2_int64)
        range_path = sigma2_state_range_path(COMMITTED, 2_int64, 1, 1)
        call sigma2_state_prepare_update(COMMITTED, candidate%to_char(), status, message)
        call assert_int(0, status, 'open the sigma2 transaction: '//trim(message))
        ! write: a Cartesian euclid pass, three particles at poses off their truth, the fourth with
        ! an observation empty from shell ZERO_FROM on against the zero reference of state 2
        call prep_sigmas_objfun(p, build)
        call build_test_volume(volume, BOX)
        allocate(zero_volume(BOX,BOX,BOX), source=0.)
        call calc%new(2, BOX, NP)
        call calc%set_ref(1, .true., volume)
        call calc%set_ref(2, .true., zero_volume)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2 = 1.
        do i = 1, NP
            truth          = real(euler2m([19. + 5.*real(i), 37., 28.]), dp)
            committed_pose = real(euler2m([20. + 5.*real(i), 36.2, 28.7]), dp)
            call calc%predict(1, .true., truth, [0.3_dp, -0.2_dp], observed)
            if( i == NP )then
                do k = -BOX/2, BOX/2
                    do h = -BOX/2, BOX/2
                        if( nint(sqrt(real(h*h + k*k))) >= ZERO_FROM ) observed(h,k) = cmplx(0., 0.)
                    end do
                end do
            endif
            call calc%set_ptcl(i, observed, no_ctf, sigma2, KFROMTO)
            call build%esig%calc_sigma2(calc, i, i, merge(2, 1, i == NP), .true., committed_pose, [0.1_dp, 0._dp])
            written(:,i) = SEED
            written(KFROMTO(1):KFROMTO(2),i) = build%esig%get_sigma2_part(i)
        end do
        call assert_true(all(written(ZERO_FROM:KFROMTO(2),NP) == 0.) .and. all(written(KFROMTO(1):ZERO_FROM-1,NP) > 0.), &
            &'N32 fixture: the fourth record does not have exactly zero shells')
        call build%esig%write_sigma2
        ! reduce and commit
        call sigma2_state_merge_local_ranges(candidate%to_char(), [range_path], [(.true., i = 1, NP)], status, message)
        call assert_int(0, status, 'merge the Cartesian range: '//trim(message))
        active    = .true.
        eo        = [(modulo(i-1, 2), i = 1, NP)]
        group_ids = 1
        call sigma2_state_reduce_groups(candidate%to_char(), active, eo, group_ids, status, message)
        call assert_int(0, status, 'a record with zero shells is skipped, not fatal: '//trim(message))
        call sigma2_state_commit(candidate%to_char(), COMMITTED, active, eo, group_ids, status, message)
        call assert_int(0, status, 'commit the Cartesian generation: '//trim(message))
        call sigma2_state_read_header(COMMITTED, header, status, message)
        call assert_true(status == 0 .and. header%generation == 2_int64, 'the Cartesian generation was not committed')
        call sigma2_state_read_groups(COMMITTED, groups, status, message)
        call assert_int(0, status, 'read the committed groups: '//trim(message))
        mean_even = 0.5*(written(:,1) + written(:,3))
        call assert_true(all(abs(groups(:,1,1) - mean_even) <= 1.e-6*mean_even), &
            &'the even group is not the mean of its records')
        call assert_true(all(abs(groups(:,2,1) - written(:,2)) <= 1.e-6*written(:,2)), &
            &'the odd group is not the mean of its valid records alone')
        ! reload: the polar sigma owner of the next pass
        call build%esig%kill
        p%l_cart_refine = .false.
        call build%pftc%new(p, 1, [1, NP], KFROMTO)
        call prep_sigmas_objfun(p, build)
        do i = 1, NP
            band = build%esig%get_sigma2_part(i)
            call assert_true(all(band == written(KFROMTO(1):KFROMTO(2),i)), 'a particle record did not reload as written')
            call assert_true(all(build%esig%sigma2_noise(:,i) == groups(:,eo(i)+1,1)), &
                &'a particle did not reload the group of its half')
            call assert_true(all(ieee_is_finite(build%esig%sigma2_noise(:,i))) .and. all(build%esig%sigma2_noise(:,i) > 0.), &
                &'a reloaded alignment sigma2 is not finite and positive')
        end do
        call build%esig%kill
        call build%pftc%kill
        call calc%kill
        call build%spproj%kill
        call cline%kill
        call remove_sigma2_files(COMMITTED, 2)
    end subroutine test_canonical_sigma_round_trip

    ! FIXTURES OF THE CANONICAL SIGMA2 STATE

    !> A committed canonical sigma2 state at path (generation 1, one global group, shells 1..box/2 of
    !! np particles, every value seed) and a builder whose project registers it, np particles in
    !! state 1 with alternating halves; params gets what the sigma owner of a Cartesian euclid pass reads.
    subroutine new_sigma2_fixture( path, np, box, seed, params, build )
        character(len=*),       intent(in)    :: path
        integer,                intent(in)    :: np, box
        real,                   intent(in)    :: seed
        type(parameters),       intent(inout) :: params
        type(builder),  target, intent(inout) :: build
        type(sigma2_state_header) :: header
        type(string) :: candidate, range_path
        real(real32) :: spectra(box/2,np)
        logical      :: active(np)
        integer      :: eo(np), groups(np), i, status
        character(len=128) :: message
        call remove_sigma2_files(path, 3)
        candidate = sigma2_state_candidate_path(path, 1_int64)
        range_path = sigma2_state_range_path(path, 1_int64, 1, 1)
        call sigma2_state_init_header(header, 1, box/2, np, box, 1.0, 1, SIGMA2_GROUP_GLOBAL, 1_int64, &
            &FIXTURE_DIGEST, SIGMA2_PROV_RESIDUAL)
        call sigma2_state_create_candidate(candidate%to_char(), header, status, message)
        call assert_int(0, status, 'sigma2 fixture candidate: '//trim(message))
        spectra = real(seed, real32)
        call sigma2_state_write_local_range(range_path%to_char(), 1_int64, FIXTURE_DIGEST, 1, spectra, 1, box/2, status, message)
        call assert_int(0, status, 'sigma2 fixture range: '//trim(message))
        call sigma2_state_merge_local_ranges(candidate%to_char(), [range_path], [(.true., i = 1, np)], status, message)
        call assert_int(0, status, 'sigma2 fixture merge: '//trim(message))
        active = .true.
        eo     = [(modulo(i-1, 2), i = 1, np)]
        groups = 1
        call sigma2_state_reduce_groups(candidate%to_char(), active, eo, groups, status, message)
        call assert_int(0, status, 'sigma2 fixture reduction: '//trim(message))
        call sigma2_state_commit(candidate%to_char(), path, active, eo, groups, status, message)
        call assert_int(0, status, 'sigma2 fixture commit: '//trim(message))
        call del_file(range_path)
        call build%spproj%set_sigma2_state_path(string(path))
        call build%spproj%os_ptcl3D%new(np, is_ptcl=.true.)
        do i = 1, np
            call build%spproj%os_ptcl3D%set_state(i, 1)
            call build%spproj%os_ptcl3D%set(i, 'eo', eo(i))
        end do
        build%spproj_field => build%spproj%os_ptcl3D
        params%cc_objfun     = OBJFUN_EUCLID
        params%cc_emit_sigma = 'no'
        params%l_cart_refine = .true.
        params%l_cont_polish = .false.
        params%l_sigma_glob  = .true.
        params%fromp         = 1
        params%top           = np
        params%part          = 1
        params%numlen        = 1
        params%box           = box
    end subroutine new_sigma2_fixture

    !> The committed state at path and the candidates and part-1 ranges of generations 1..ngen.
    subroutine remove_sigma2_files( path, ngen )
        character(len=*), intent(in) :: path
        integer,          intent(in) :: ngen
        integer(int64) :: g
        call del_file(path)
        do g = 1_int64, int(ngen, int64)
            call del_file(sigma2_state_candidate_path(path, g))
            call del_file(sigma2_state_range_path(path, g, 1, 1))
        end do
    end subroutine remove_sigma2_files

end module simple_cartft_calc_tester
