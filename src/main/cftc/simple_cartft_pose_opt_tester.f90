!@descr: unit tests for the continuous pose optimizer (simple_cartft_pose_opt)
! The rotation increment (E11), the stage solvers (E12), invalid inputs (E14), the NCC solver
! (E15), both routes' stage records (E16: shift then joint on request, the joint stage alone by
! default), rollback (E17), the ori entry point and its scaling (E13, E19, N2), recovery under
! each objective (N7) and the bounds trs and athres_cont (O6). Ids: pose_cont plan, section 10.
module simple_cartft_pose_opt_tester
use simple_defs,            only: dp, DPI
use simple_defs_ori,        only: I_E1, I_E2, I_E3, I_X, I_Y, N_PTCL_ORIPARAMS
use simple_core_module_api, only: euler2m
use simple_ori,             only: ori
use simple_cartft_calc,     only: cartft_calc
use simple_cartft_pose_opt, only: cartft_pose_opt, right_increment_rotation, CARTFT_STAGE_SHIFT, CARTFT_STAGE_JOINT, &
    &CARTFT_NOT_ATTEMPTED, CARTFT_INVALID_PREPARATION, CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT, &
    &CARTFT_NO_RELIABLE_UPDATE, CARTFT_BOUND_REJECTED
use simple_type_defs,       only: ctfparams, CTFFLAG_NO, OBJFUN_CC, OBJFUN_EUCLID
use simple_test_utils
implicit none
private
public :: run_all_cartft_pose_opt_tests

integer,  parameter :: TEST_BOX       = 24
integer,  parameter :: SMALL_BOX      = 16
real(dp), parameter :: ORTHOGONAL_TOL = 2.e-12_dp
real(dp), parameter :: ROTATION_TOL   = 8.e-4_dp
real(dp), parameter :: SHIFT_TOL      = 2.e-3_dp
! the seed rotations come from single-precision Euler angles, orthonormal to ~1e-7; the
! trace form of the geodesic distance used as the expected value then carries an error of
! about 1e-7*3/sin(angle), 3e-6 rad at the 1.8 degree motion of the fixture
real(dp), parameter :: MOTION_TOL     = 1.e-4_dp

type :: fixture_record
    real(dp) :: truth_rotation(3,3), truth_shift(2), seed_rotation(3,3), seed_shift(2)
end type fixture_record

contains

    subroutine run_all_cartft_pose_opt_tests()
        write(*,'(A)') '**** running all pose optimizer tests ****'
        write(*,'(A)') 'test_rotation_increment'
        call test_rotation_increment()
        write(*,'(A)') 'test_shift_solver'
        call test_shift_solver()
        write(*,'(A)') 'test_joint_solver'
        call test_joint_solver()
        write(*,'(A)') 'test_tiny_accepted_reduction_stop'
        call test_tiny_accepted_reduction_stop()
        write(*,'(A)') 'test_ncc_solver'
        call test_ncc_solver()
        write(*,'(A)') 'test_invalid_and_unobservable_inputs'
        call test_invalid_and_unobservable_inputs()
        write(*,'(A)') 'test_route_transaction_contracts'
        call test_route_transaction_contracts()
        write(*,'(A)') 'test_rollback_transaction_contracts'
        call test_rollback_transaction_contracts()
        write(*,'(A)') 'test_ori_round_trip'
        call test_ori_round_trip()
        write(*,'(A)') 'test_objective_recovery'
        call test_objective_recovery()
        write(*,'(A)') 'test_parameter_bounds'
        call test_parameter_bounds()
    end subroutine run_all_cartft_pose_opt_tests

    ! FIXTURES

    !> Four Gaussian blobs in a TEST_BOX cube (the refiner tester's volume).
    subroutine build_test_volume( volume )
        real, allocatable, intent(out) :: volume(:,:,:)
        real, parameter :: centres(3,4) = reshape([-5., -3., 2., 4., 5., -3., 0., -6., -5., 3., -2., 6.], [3,4])
        real, parameter :: sigmas(4)     = [2., 2.5, 1.8, 2.2]
        real, parameter :: amplitudes(4) = [1., 0.8, 0.6, 0.5]
        real    :: centre, dx, dy, dz
        integer :: blob, i, j, k
        allocate(volume(TEST_BOX,TEST_BOX,TEST_BOX), source=0.)
        centre = real(TEST_BOX)/2. + 0.5
        do k = 1, TEST_BOX
            do j = 1, TEST_BOX
                do i = 1, TEST_BOX
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

    !> A one-state calculator with the volume in both halves and particle slot 1 holding the
    !! noise-free prediction at the truth pose (no CTF, shells kfromto) for objfun (unit sigma2
    !! under OBJFUN_EUCLID); the seed is the perturbed pose the refiner and adapter testers
    !! start from.
    subroutine new_fixture( calc, volume, kfromto, gain, objfun, fix )
        type(cartft_calc),    intent(inout) :: calc
        real,                 intent(in)    :: volume(:,:,:)
        integer,              intent(in)    :: kfromto(2)
        real,                 intent(in)    :: gain
        integer,              intent(in)    :: objfun
        type(fixture_record), intent(out)   :: fix
        type(ctfparams)      :: no_ctf
        complex, allocatable :: observed(:,:)
        real,    allocatable :: sigma2(:)
        integer :: box
        box = size(volume,1)
        call calc%new(1, box, 1)
        call calc%set_ref(1, .true.,  volume)
        call calc%set_ref(1, .false., volume)
        fix%truth_rotation = real(euler2m([19., 37., 28.]), dp)
        fix%truth_shift    = [0.31_dp, -0.24_dp]
        fix%seed_rotation  = real(euler2m([20., 36.2, 28.7]), dp)
        fix%seed_shift     = [-0.08_dp, 0.06_dp]
        allocate(observed(-box/2:box/2,-box/2:box/2))
        call calc%predict(1, .true., fix%truth_rotation, fix%truth_shift, observed)
        allocate(sigma2(0:box/2), source=1.)
        no_ctf%ctfflag = CTFFLAG_NO
        if( objfun == OBJFUN_CC )then
            call calc%set_ptcl(1, gain*observed, no_ctf, kfromto)
        else
            call calc%set_ptcl(1, gain*observed, no_ctf, sigma2, kfromto)
        endif
    end subroutine new_fixture

    !> Particle slot 1 holding the noise-free truth prediction centred on the stored shift
    !! seed_shift (cropped pixels), as the prepared observation is (O4): the model shift that
    !! matches it is the increment truth_shift - seed_shift.
    subroutine set_centred_particle( calc, fix, seed_shift, kfromto, objfun )
        type(cartft_calc),    intent(inout) :: calc
        type(fixture_record), intent(in)    :: fix
        real(dp),             intent(in)    :: seed_shift(2)
        integer,              intent(in)    :: kfromto(2), objfun
        type(ctfparams)      :: no_ctf
        complex, allocatable :: observed(:,:)
        real,    allocatable :: sigma2(:)
        integer :: box
        box = calc%get_box()
        allocate(observed(-box/2:box/2,-box/2:box/2))
        call calc%predict(1, .true., fix%truth_rotation, fix%truth_shift - seed_shift, observed)
        no_ctf%ctfflag = CTFFLAG_NO
        if( objfun == OBJFUN_CC )then
            call calc%set_ptcl(1, observed, no_ctf, kfromto)
        else
            allocate(sigma2(0:box/2), source=1.)
            call calc%set_ptcl(1, observed, no_ctf, sigma2, kfromto)
        endif
    end subroutine set_centred_particle

    pure function rotation_distance( left, right ) result( distance )
        real(dp), intent(in) :: left(3,3), right(3,3)
        real(dp) :: distance, cosine
        cosine   = 0.5_dp*(sum(left*right) - 1._dp)
        distance = acos(max(-1._dp, min(1._dp, cosine)))
    end function rotation_distance

    pure function identity_rotation() result( rotation )
        real(dp) :: rotation(3,3)
        rotation      = 0._dp
        rotation(1,1) = 1._dp
        rotation(2,2) = 1._dp
        rotation(3,3) = 1._dp
    end function identity_rotation

    pure function determinant3( matrix ) result( determinant )
        real(dp), intent(in) :: matrix(3,3)
        real(dp) :: determinant
        determinant = matrix(1,1)*(matrix(2,2)*matrix(3,3) - matrix(2,3)*matrix(3,2)) - &
            &matrix(1,2)*(matrix(2,1)*matrix(3,3) - matrix(2,3)*matrix(3,1)) + &
            &matrix(1,3)*(matrix(2,1)*matrix(3,2) - matrix(2,2)*matrix(3,1))
    end function determinant3

    subroutine get_stage_counts( opt, istage, status, niter, nattempted, naccepted, nbound_hits, &
        &max_rotation_step, max_shift_step, objective_before, objective_after )
        type(cartft_pose_opt), intent(in)  :: opt
        integer,               intent(in)  :: istage
        integer,               intent(out) :: status, niter, nattempted, naccepted, nbound_hits
        real(dp),              intent(out) :: max_rotation_step, max_shift_step, objective_before, objective_after
        call opt%get_stage(istage, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, objective_before, objective_after)
    end subroutine get_stage_counts

    ! TESTS

    ! E11
    subroutine test_rotation_increment()
        real(dp), parameter :: angle = 0.02_dp
        real(dp) :: rotation(3,3), updated(3,3), identity(3,3)
        real(dp) :: axis_rotation(3,3), expected(3,3)
        real(dp) :: determinant, input_determinant, input_orthogonality, updated_orthogonality
        rotation = real(euler2m([19., 37., 28.]), dp)
        updated  = right_increment_rotation(rotation, [angle, 0._dp, 0._dp])
        identity = identity_rotation()
        axis_rotation      = identity
        axis_rotation(2,2) = cos(angle)
        axis_rotation(2,3) = -sin(angle)
        axis_rotation(3,2) = sin(angle)
        axis_rotation(3,3) = cos(angle)
        expected = matmul(rotation, axis_rotation)
        input_orthogonality   = sqrt(sum((matmul(transpose(rotation), rotation) - identity)**2))
        input_determinant     = determinant3(rotation)
        updated_orthogonality = sqrt(sum((matmul(transpose(updated), updated) - identity)**2))
        determinant           = determinant3(updated)
        call assert_true(updated_orthogonality <= input_orthogonality + ORTHOGONAL_TOL, &
            &'right rotation increment increased the input orthogonality error')
        call assert_true(abs(determinant - 1._dp) <= abs(input_determinant - 1._dp) + ORTHOGONAL_TOL, &
            &'right rotation increment increased the input determinant error')
        call assert_true(maxval(abs(updated - expected)) <= ORTHOGONAL_TOL, &
            &'rotation increment does not use the declared right-handed local-axis convention')
    end subroutine test_rotation_increment

    ! E12 (shift-only stage)
    subroutine test_shift_solver()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: shift(2), max_rotation_step, max_shift_step, objective_before, objective_after
        real(dp) :: gradient_exact(5)
        integer  :: status, niter, nattempted, naccepted, nbound_hits
        call build_test_volume(volume)
        call new_fixture(calc, volume, [2, TEST_BOX/2], 1., OBJFUN_EUCLID, fix)
        call opt%new(TEST_BOX, TEST_BOX, 5., 15., maxits=20, shift_step_bound=1._dp)
        shift = [-0.15_dp, 0.12_dp]
        call opt%refine_shift(calc, 1, .true., 1, fix%truth_rotation, shift)
        call get_stage_counts(opt, CARTFT_STAGE_SHIFT, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, objective_before, objective_after)
        call assert_true(status == CARTFT_ACCEPTED .and. sqrt(sum((shift - fix%truth_shift)**2)) < SHIFT_TOL, &
            &'shift-only LM did not recover a known shift')
        call assert_true(naccepted > 0 .and. max_shift_step <= 1._dp + epsilon(1._dp), &
            &'shift-only LM violated its accepted-step contract')
        call assert_true(niter < 20, 'shift-only LM exhausted its iteration limit after converging')
        shift = fix%truth_shift
        call calc%objective_gradient(1, .true., 1, fix%truth_rotation, shift, objective_before, gradient_exact)
        call assert_true(objective_before <= real(epsilon(1.), dp)**2, &
            &'exact shift fixture exceeds the single-precision residual floor')
        call opt%refine_shift(calc, 1, .true., 1, fix%truth_rotation, shift)
        call get_stage_counts(opt, CARTFT_STAGE_SHIFT, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, objective_before, objective_after)
        call assert_true(status == CARTFT_NO_IMPROVEMENT, 'exact shift did not report no improvement')
        call assert_true(all(shift == fix%truth_shift), 'exact shift was not retained')
        call assert_true(nattempted == 0 .and. naccepted == 0, 'exact shift attempted a roundoff-level update')
        call opt%kill
        call calc%kill
    end subroutine test_shift_solver

    ! E12 (joint stage): recovery, exact pose, masks, cumulative guard
    subroutine test_joint_solver()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), frozen_rotation(3,3), shift(2), frozen_shift(2)
        real(dp) :: objective_before, objective_after, gradient(5)
        real(dp) :: max_rotation_step, max_shift_step, stage_before, stage_after
        integer  :: status, niter, nattempted, naccepted, nbound_hits
        call build_test_volume(volume)
        call new_fixture(calc, volume, [2, TEST_BOX/2], 1., OBJFUN_EUCLID, fix)
        call opt%new(TEST_BOX, TEST_BOX, 5., 15., maxits=20, rotation_scale=0.10_dp, shift_step_bound=1._dp)
        rotation = fix%seed_rotation
        shift    = fix%seed_shift
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift)
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, stage_before, stage_after)
        call assert_true(status == CARTFT_ACCEPTED .and. objective_after < objective_before, &
            &'joint LM did not accept an objective-reducing pose')
        call assert_true(rotation_distance(rotation, fix%truth_rotation) < ROTATION_TOL .and. &
            &sqrt(sum((shift - fix%truth_shift)**2)) < SHIFT_TOL, 'joint LM did not recover the known five-parameter pose')
        call assert_true(max_rotation_step <= 0.10_dp + epsilon(1._dp) .and. max_shift_step <= 1._dp + epsilon(1._dp), &
            &'joint LM exceeded a configured proposal bound')
        call assert_true(niter < 20, 'joint LM exhausted its iteration limit after converging')
        rotation = fix%truth_rotation
        shift    = fix%truth_shift
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift)
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_JOINT) == CARTFT_NO_IMPROVEMENT .and. &
            &all(rotation == fix%truth_rotation) .and. all(shift == fix%truth_shift), 'exact joint pose was not retained')
        ! shifts only
        rotation        = fix%seed_rotation
        shift           = fix%seed_shift
        frozen_rotation = rotation
        frozen_shift    = shift
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift, active=[.false., .false., .false., .true., .true.])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_JOINT) == CARTFT_ACCEPTED .and. &
            &objective_after < objective_before .and. sqrt(sum((shift - frozen_shift)**2)) > 1.e-10_dp .and. &
            &all(rotation == frozen_rotation), 'shift-only joint solve did not improve only active shifts')
        ! rotations only
        rotation        = fix%seed_rotation
        shift           = fix%seed_shift
        frozen_rotation = rotation
        frozen_shift    = shift
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift, active=[.true., .true., .true., .false., .false.])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_JOINT) == CARTFT_ACCEPTED .and. &
            &objective_after < objective_before .and. rotation_distance(rotation, frozen_rotation) > 1.e-10_dp .and. &
            &all(shift == frozen_shift), 'rotation-only joint solve did not improve only active rotations')
        ! cumulative guard anchored at the input pose
        call opt%new(TEST_BOX, TEST_BOX, 1.e-12, real(1.e-12_dp*180._dp/DPI), maxits=20, rotation_scale=0.10_dp, &
            &shift_step_bound=1._dp)
        rotation        = fix%seed_rotation
        shift           = fix%seed_shift
        frozen_rotation = rotation
        frozen_shift    = shift
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift, anchor_rotmat=frozen_rotation, anchor_shift=frozen_shift)
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, stage_before, stage_after)
        call assert_true(nbound_hits > 0, 'cumulative guard test did not exercise an out-of-bound proposal')
        call assert_true(status == CARTFT_BOUND_REJECTED .and. niter == 8, 'cumulative guard did not stop after eight rejected proposals')
        call assert_true(all(rotation == frozen_rotation) .and. all(shift == frozen_shift) .and. &
            &objective_after == objective_before, 'cumulative-bound rejection changed the complete input pose')
        call opt%kill
        call calc%kill
    end subroutine test_joint_solver

    ! E12 (stopping): a material rotation residual with only a minuscule shift step permitted,
    ! so the first accepted reduction is positive but below the relative-reduction threshold
    subroutine test_tiny_accepted_reduction_stop()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), frozen_rotation(3,3), shift(2), initial_shift(2)
        real(dp) :: objective_before, objective_after, relative_reduction, gradient(5)
        real(dp) :: max_rotation_step, max_shift_step, stage_before, stage_after
        integer  :: status, niter, nattempted, naccepted, nbound_hits
        call build_test_volume(volume)
        call new_fixture(calc, volume, [2, TEST_BOX/2], 1., OBJFUN_EUCLID, fix)
        rotation        = fix%seed_rotation
        shift           = fix%seed_shift
        frozen_rotation = rotation
        initial_shift   = shift
        call opt%new(TEST_BOX, TEST_BOX, 5., 15., maxits=20, rotation_scale=0.10_dp, shift_step_bound=1.e-8_dp)
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift, active=[.false., .false., .false., .true., .true.])
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, stage_before, stage_after)
        relative_reduction = (objective_before - objective_after)/objective_before
        call assert_true(status == CARTFT_ACCEPTED .and. niter == 1, 'LM did not stop after the first tiny accepted reduction')
        call assert_true(relative_reduction > 0._dp .and. relative_reduction < 1.e-6_dp, &
            &'tiny accepted reduction fixture does not exercise the production threshold')
        call assert_true(all(rotation == frozen_rotation) .and. any(shift /= initial_shift), &
            &'tiny accepted reduction changed inactive coordinates or lost its endpoint')
        call opt%kill
        call calc%kill
    end subroutine test_tiny_accepted_reduction_stop

    ! E15
    subroutine test_ncc_solver()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), shift(2), objective_before, objective_after, gradient(5)
        call build_test_volume(volume)
        call new_fixture(calc, volume, [2, TEST_BOX/2], 2., OBJFUN_CC, fix)
        call opt%new(TEST_BOX, TEST_BOX, 5., 15., maxits=20, rotation_scale=0.10_dp, shift_step_bound=1._dp)
        rotation = fix%seed_rotation
        shift    = fix%seed_shift
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_before, gradient)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift)
        call calc%objective_gradient(1, .true., 1, rotation, shift, objective_after, gradient)
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_JOINT) == CARTFT_ACCEPTED .and. objective_after < objective_before, &
            &'Cartesian NCC LM did not accept a correlation-improving pose')
        call assert_true(rotation_distance(rotation, fix%truth_rotation) < ROTATION_TOL .and. &
            &sqrt(sum((shift - fix%truth_shift)**2)) < SHIFT_TOL, 'Cartesian NCC LM did not recover the known gain-scaled pose')
        call opt%kill
        call calc%kill
    end subroutine test_ncc_solver

    ! E14
    subroutine test_invalid_and_unobservable_inputs()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(ctfparams)       :: no_ctf
        real     :: sigma2(0:TEST_BOX/2), zero_volume(TEST_BOX,TEST_BOX,TEST_BOX)
        complex  :: zero_plane(-TEST_BOX/2:TEST_BOX/2,-TEST_BOX/2:TEST_BOX/2)
        real(dp) :: rotation(3,3), original_rotation(3,3), shift(2), original_shift(2)
        zero_volume    = 0.
        zero_plane     = cmplx(0., 0.)
        no_ctf%ctfflag = CTFFLAG_NO
        sigma2    = 1.
        sigma2(2) = -1.
        call calc%new(1, TEST_BOX, 1)
        call calc%set_ref(1, .true., zero_volume)
        call calc%set_ptcl(1, zero_plane, no_ctf, sigma2, [2, TEST_BOX/2])
        call assert_true(.not. calc%ptcl_is_valid(1), 'invalid sigma data was accepted')
        sigma2 = 1.
        ! both objectives normalize by the particle power, so an observation without power is an
        ! invalid preparation; the unobservable case is a zero reference against a particle
        call calc%set_ptcl(1, zero_plane, no_ctf, sigma2, [2, TEST_BOX/2])
        call assert_true(.not. calc%ptcl_is_valid(1), 'an observation without power was accepted')
        call calc%set_ptcl(1, zero_plane + cmplx(1., -0.5), no_ctf, sigma2, [2, TEST_BOX/2])
        rotation = real(euler2m([19., 37., 28.]), dp)
        shift    = [0.2_dp, -0.1_dp]
        original_rotation = rotation
        original_shift    = shift
        call opt%new(TEST_BOX, TEST_BOX, 5., 15., maxits=10, rotation_scale=0.10_dp, shift_step_bound=1._dp)
        call opt%refine_joint(calc, 1, .true., 1, rotation, shift)
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_JOINT) == CARTFT_NO_RELIABLE_UPDATE .and. &
            &all(rotation == original_rotation) .and. all(shift == original_shift), 'unobservable joint solve changed the input pose')
        call opt%kill
        call calc%kill
    end subroutine test_invalid_and_unobservable_inputs

    ! E16: the shift-then-joint route (shift_first) completes both stages and the stage records
    ! chain the objectives of the transaction; the default (production) route runs the joint
    ! stage alone from the seed and commits an improving pose
    subroutine test_route_transaction_contracts()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), shift(2), objective_before, objective_after, rot_motion, shift_motion
        real(dp) :: sh_rot_step, sh_shift_step, sh_before, sh_after, jt_rot_step, jt_shift_step, jt_before, jt_after
        integer  :: sh_status, sh_niter, sh_att, sh_acc, sh_hits, jt_status, jt_niter, jt_att, jt_acc, jt_hits
        call build_small_volume(volume)
        call new_fixture(calc, volume, [2, SMALL_BOX/2 - 1], 1., OBJFUN_CC, fix)
        call assert_true(calc%ptcl_is_valid(1) .and. all(calc%get_ptcl_kfromto(1) == [2, SMALL_BOX/2 - 1]), &
            &'transaction fixture did not prepare a valid particle on the requested shells')
        ! shift stage, then joint stage
        call opt%new(SMALL_BOX, SMALL_BOX, 5., 15., shift_first=.true.)
        rotation = fix%seed_rotation
        shift    = fix%seed_shift
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call opt%get_objectives(objective_before, objective_after)
        call opt%get_motion(rot_motion, shift_motion)
        call get_stage_counts(opt, CARTFT_STAGE_SHIFT, sh_status, sh_niter, sh_att, sh_acc, sh_hits, &
            &sh_rot_step, sh_shift_step, sh_before, sh_after)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, jt_status, jt_niter, jt_att, jt_acc, jt_hits, &
            &jt_rot_step, jt_shift_step, jt_before, jt_after)
        call assert_true(opt%get_status() == CARTFT_ACCEPTED .and. objective_after < objective_before, &
            &'transaction did not commit an improving pose')
        call assert_true(any(sh_status == [CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT]) .and. &
            &any(jt_status == [CARTFT_ACCEPTED, CARTFT_NO_IMPROVEMENT]), 'transaction did not complete both LM stages')
        call assert_true(sh_before == objective_before .and. jt_before == sh_after .and. jt_after == objective_after, &
            &'stage records do not chain the transaction objectives')
        call assert_true(sh_att > 0 .and. jt_att > 0 .and. sh_rot_step == 0._dp, 'stage records do not account for both stages')
        call assert_true(abs(rot_motion - rotation_distance(rotation, fix%seed_rotation)) < MOTION_TOL .and. &
            &abs(shift_motion - sqrt(sum((shift - fix%seed_shift)**2))) < 1.e-12_dp, &
            &'reported motion is not the seed-to-result motion')
        ! the production route: the joint stage alone, from the seed
        call opt%new(SMALL_BOX, SMALL_BOX, 5., 15.)
        rotation = fix%seed_rotation
        shift    = fix%seed_shift
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call opt%get_objectives(objective_before, objective_after)
        call get_stage_counts(opt, CARTFT_STAGE_SHIFT, sh_status, sh_niter, sh_att, sh_acc, sh_hits, &
            &sh_rot_step, sh_shift_step, sh_before, sh_after)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, jt_status, jt_niter, jt_att, jt_acc, jt_hits, &
            &jt_rot_step, jt_shift_step, jt_before, jt_after)
        call assert_true(sh_status == CARTFT_NOT_ATTEMPTED .and. sh_att == 0, 'the default route ran the shift-only stage')
        call assert_true(opt%get_status() == CARTFT_ACCEPTED .and. objective_after < objective_before .and. &
            &jt_before == objective_before .and. jt_after == objective_after .and. jt_shift_step > 0._dp, &
            &'the joint stage alone did not commit an improving pose with shifts from the seed')
        call opt%kill
        call calc%kill
    end subroutine test_route_transaction_contracts

    ! E17: bound rejection, finite no-improvement and invalid preparation return the seed bit for bit
    subroutine test_rollback_transaction_contracts()
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), shift(2), max_rotation_step, max_shift_step, stage_before, stage_after
        integer  :: status, niter, nattempted, naccepted, sh_hits, jt_hits
        call build_small_volume(volume)
        call new_fixture(calc, volume, [2, SMALL_BOX/2 - 1], 1., OBJFUN_CC, fix)
        call opt%new(SMALL_BOX, SMALL_BOX, 1.e-30, 15.)
        rotation = fix%truth_rotation
        shift    = fix%seed_shift
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call get_stage_counts(opt, CARTFT_STAGE_SHIFT, status, niter, nattempted, naccepted, sh_hits, &
            &max_rotation_step, max_shift_step, stage_before, stage_after)
        call get_stage_counts(opt, CARTFT_STAGE_JOINT, status, niter, nattempted, naccepted, jt_hits, &
            &max_rotation_step, max_shift_step, stage_before, stage_after)
        call assert_true(opt%get_status() == CARTFT_BOUND_REJECTED .and. sh_hits + jt_hits > 0, &
            &'transaction did not report a cumulative shift-bound rejection')
        call assert_true(all(rotation == fix%truth_rotation) .and. all(shift == fix%seed_shift), &
            &'bound-rejected transaction did not preserve the input pose')
        call opt%new(SMALL_BOX, SMALL_BOX, 5., 15.)
        rotation = fix%truth_rotation
        shift    = fix%truth_shift
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call assert_true(opt%get_status() == CARTFT_NO_IMPROVEMENT, 'exact-pose transaction did not report finite no-improvement')
        call assert_true(all(rotation == fix%truth_rotation) .and. all(shift == fix%truth_shift), &
            &'non-improving transaction changed the input pose')
        ! a slot that was never prepared
        call calc%new(1, SMALL_BOX, 1)
        call calc%set_ref(1, .true., volume)
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call assert_true(opt%get_status() == CARTFT_INVALID_PREPARATION .and. &
            &all(rotation == fix%truth_rotation) .and. all(shift == fix%truth_shift), &
            &'invalid particle preparation changed the input pose')
        call opt%kill
        call calc%kill
    end subroutine test_rollback_transaction_contracts

    ! E19 (round trip), E13 (scaling), N2: the ori entry point on a cropped calculator
    ! (native box twice the cropped box). A perturbed seed in native pixels is refined to the
    ! truth expressed in native pixels; only Euler angles and shift change; a rejected solve
    ! leaves the record bit-identical.
    subroutine test_ori_round_trip()
        integer, parameter :: NATIVE_BOX = 2*SMALL_BOX
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        type(ori)             :: o
        real, allocatable :: volume(:,:,:)
        real     :: prec_before(N_PTCL_ORIPARAMS), prec_after(N_PTCL_ORIPARAMS), native_scale
        logical  :: pose_slot(N_PTCL_ORIPARAMS)
        integer  :: eo
        call build_small_volume(volume)
        call new_fixture(calc, volume, [2, SMALL_BOX/2 - 1], 1., OBJFUN_CC, fix)
        native_scale = real(NATIVE_BOX)/real(SMALL_BOX)
        pose_slot = .false.
        pose_slot([I_E1, I_E2, I_E3, I_X, I_Y]) = .true.
        do eo = 0, 1
            ! accepted solve from a perturbed seed, half eo (both halves hold the same volume);
            ! the slot holds the observation centred on the stored shift (O4)
            call make_seed_ori(o, [20., 36.2, 28.7], native_scale*real(fix%seed_shift), eo)
            call set_centred_particle(calc, fix, real(o%get_2Dshift(), dp)/native_scale, [2, SMALL_BOX/2 - 1], OBJFUN_CC)
            call o%ori2prec(prec_before)
            call opt%new(NATIVE_BOX, SMALL_BOX, 5., 15.)
            call opt%refine(calc, 1, o)
            call o%ori2prec(prec_after)
            call assert_int(CARTFT_ACCEPTED, opt%get_status(), 'ori solve from a perturbed seed was not accepted')
            call assert_true(rotation_distance(real(o%get_mat(), dp), fix%truth_rotation) < ROTATION_TOL, &
                &'ori solve did not recover the rotation')
            call assert_true(sqrt(sum((real(o%get_2Dshift(), dp) - native_scale*fix%truth_shift)**2)) < native_scale*SHIFT_TOL, &
                &'ori solve did not return the shift in native pixels')
            call assert_true(all(prec_after == prec_before .or. pose_slot), 'ori solve changed a non-pose field')
            call assert_true(o%get_state() == 1 .and. o%get_eo() == eo .and. o%get_proj() == 7, &
                &'ori solve changed state, half or projection')
            ! rejected solve: the record is unchanged bit for bit
            call make_seed_ori(o, [20., 36.2, 28.7], native_scale*real(fix%seed_shift), eo)
            call o%ori2prec(prec_before)
            call opt%new(NATIVE_BOX, SMALL_BOX, 1.e-30, 15.)
            call opt%refine(calc, 1, o)
            call o%ori2prec(prec_after)
            call assert_int(CARTFT_BOUND_REJECTED, opt%get_status(), 'bounded ori solve was not rejected')
            call assert_true(all(prec_after == prec_before), 'rejected ori solve changed the record')
        end do
        call o%kill
        call opt%kill
        call calc%kill
    end subroutine test_ori_round_trip

    ! N7: recovery of an injected rotation and shift under each objective through the ori entry
    ! point of production: native box twice the cropped box, the observation centred on the
    ! stored shift, the solve from a zero increment committing stored shift plus increment.
    ! Expected value: the injected (truth) pose in native pixels.
    subroutine test_objective_recovery()
        integer, parameter :: NATIVE_BOX = 2*TEST_BOX
        integer, parameter :: objfuns(2) = [OBJFUN_CC, OBJFUN_EUCLID]
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        type(ori)             :: o
        real, allocatable :: volume(:,:,:)
        real(dp) :: objective_before, objective_after
        real     :: native_scale
        integer  :: iobj
        call build_test_volume(volume)
        native_scale = real(NATIVE_BOX)/real(TEST_BOX)
        do iobj = 1, size(objfuns)
            call new_fixture(calc, volume, [2, TEST_BOX/2], 1., objfuns(iobj), fix)
            call make_seed_ori(o, [20., 36.2, 28.7], native_scale*real(fix%seed_shift), 1)
            call set_centred_particle(calc, fix, real(o%get_2Dshift(), dp)/native_scale, [2, TEST_BOX/2], objfuns(iobj))
            call assert_int(objfuns(iobj), calc%get_ptcl_objfun(1), 'recovery fixture prepared the wrong objective')
            call opt%new(NATIVE_BOX, TEST_BOX, 5., 15.)
            call opt%refine(calc, 1, o)
            call opt%get_objectives(objective_before, objective_after)
            call assert_int(CARTFT_ACCEPTED, opt%get_status(), 'injected pose was not recovered: solve not accepted')
            call assert_true(objective_after < objective_before, 'recovery did not lower the objective')
            call assert_true(rotation_distance(real(o%get_mat(), dp), fix%truth_rotation) < ROTATION_TOL, &
                &'injected rotation was not recovered')
            call assert_true(sqrt(sum((real(o%get_2Dshift(), dp) - native_scale*fix%truth_shift)**2)) < native_scale*SHIFT_TOL, &
                &'injected shift was not recovered as stored shift plus increment')
            call calc%kill
        end do
        call o%kill
        call opt%kill
    end subroutine test_objective_recovery

    ! O6: the total rotation bound is athres_cont (degrees) and the total shift bound trs (native
    ! pixels); trs = 0 freezes the shifts as l_doshift=.false. does in the polar branch.
    ! Expected values: the seed shift unchanged under trs = 0 while the rotation improves; an
    ! accepted motion within the rotation bound; the 1.8 degree seed error rejected under a
    ! 0.5 degree bound only if no improving pose exists inside it.
    subroutine test_parameter_bounds()
        integer, parameter :: NATIVE_BOX = 2*SMALL_BOX
        type(cartft_calc)     :: calc
        type(cartft_pose_opt) :: opt
        type(fixture_record)  :: fix
        type(ori)             :: o
        real, allocatable :: volume(:,:,:)
        real(dp) :: rotation(3,3), shift(2), rot_motion, shift_motion
        real     :: seed_native(2), native_scale
        call build_small_volume(volume)
        call new_fixture(calc, volume, [2, SMALL_BOX/2 - 1], 1., OBJFUN_CC, fix)
        native_scale = real(NATIVE_BOX)/real(SMALL_BOX)
        ! trs = 0: rotation-only refinement, the stored shift is kept bit for bit
        call make_seed_ori(o, [20., 36.2, 28.7], native_scale*real(fix%seed_shift), 0)
        call set_centred_particle(calc, fix, real(o%get_2Dshift(), dp)/native_scale, [2, SMALL_BOX/2 - 1], OBJFUN_CC)
        seed_native = o%get_2Dshift()
        call opt%new(NATIVE_BOX, SMALL_BOX, 0., 15.)
        call opt%refine(calc, 1, o)
        call opt%get_motion(rot_motion, shift_motion)
        call assert_int(CARTFT_ACCEPTED, opt%get_status(), 'rotation-only solve under trs=0 was not accepted')
        call assert_true(all(o%get_2Dshift() == seed_native) .and. shift_motion == 0._dp .and. rot_motion > 0._dp, &
            &'trs=0 did not freeze the shift while refining the rotation')
        call assert_true(opt_stage_status(opt, CARTFT_STAGE_SHIFT) == CARTFT_NOT_ATTEMPTED, &
            &'trs=0 ran the shift-only stage')
        ! the rotation bound limits the committed rotation from the seed
        rotation = fix%seed_rotation
        shift    = fix%seed_shift
        call new_fixture(calc, volume, [2, SMALL_BOX/2 - 1], 1., OBJFUN_CC, fix)
        call opt%new(SMALL_BOX, SMALL_BOX, 5., 0.5)
        call opt%refine_pose(calc, 1, .true., 1, rotation, shift)
        call opt%get_motion(rot_motion, shift_motion)
        call assert_true(any(opt%get_status() == [CARTFT_ACCEPTED, CARTFT_BOUND_REJECTED]), &
            &'rotation-bounded solve neither committed inside the bound nor rejected')
        call assert_true(rot_motion <= 0.5_dp*DPI/180._dp + 1.e-9_dp, 'committed rotation exceeds the rotation bound')
        call assert_true(rotation_distance(rotation, fix%truth_rotation) > 1._dp*DPI/180._dp, &
            &'rotation-bounded solve reached a truth outside its bound')
        call o%kill
        call opt%kill
        call calc%kill
    end subroutine test_parameter_bounds

    subroutine make_seed_ori( o, euler, native_shift, eo )
        type(ori), intent(inout) :: o
        real,      intent(in)    :: euler(3), native_shift(2)
        integer,   intent(in)    :: eo
        call o%new(.true.)
        call o%set_euler(euler)
        call o%set_shift(native_shift)
        call o%set_state(1)
        call o%set('eo',    real(eo))
        call o%set('proj',  7.)
        call o%set('corr',  0.42)
        call o%set('inpl',  13.)
        call o%set('class', 3.)
        call o%set('w',     0.8)
        call o%set('dist',  2.5)
        call o%set('frac',  50.)
    end subroutine make_seed_ori

    integer function opt_stage_status( opt, istage ) result( status )
        type(cartft_pose_opt), intent(in) :: opt
        integer,               intent(in) :: istage
        integer  :: niter, nattempted, naccepted, nbound_hits
        real(dp) :: max_rotation_step, max_shift_step, objective_before, objective_after
        call opt%get_stage(istage, status, niter, nattempted, naccepted, nbound_hits, &
            &max_rotation_step, max_shift_step, objective_before, objective_after)
    end function opt_stage_status

end module simple_cartft_pose_opt_tester
