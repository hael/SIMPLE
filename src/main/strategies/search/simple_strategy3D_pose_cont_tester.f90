!@descr: unit tests for standalone pose_cont strategy policy (simple_strategy3D_pose_cont)
module simple_strategy3D_pose_cont_tester
use simple_cmdline, only: cmdline
use simple_ori, only: ori
use simple_parameters, only: parameters
use simple_strategy3D_pose_cont, only: pose_cont_seed_is_valid, pose_cont_sigma_is_enabled
use simple_test_utils
use simple_type_defs, only: OBJFUN_CC, OBJFUN_EUCLID
implicit none
private
public :: run_all_strategy3D_pose_cont_tests

contains

    subroutine run_all_strategy3D_pose_cont_tests()
        write (*, '(A)') '**** running all pose_cont strategy tests ****'
        write (*, '(A)') 'test_outer_objective_sigma_policy'
        call test_outer_objective_sigma_policy()
        write (*, '(A)') 'test_seed_contract'
        call test_seed_contract()
    end subroutine run_all_strategy3D_pose_cont_tests

    subroutine test_outer_objective_sigma_policy()
        type(cmdline) :: cline
        type(parameters) :: params

        call cline%set('objfun', 'euclid')
        call cline%set('pose_cont', 'yes')
        call params%new(cline, silent=.true.)
        call assert_int(OBJFUN_EUCLID, params%cc_objfun, &
            &'Euclidean pose-cont changed the outer objective')
        call assert_char('no', trim(params%cc_emit_sigma), &
            &'Euclidean pose-cont unexpectedly changed cc_emit_sigma')
        call assert_true(pose_cont_sigma_is_enabled(params%cc_objfun, params%cc_emit_sigma), &
            &'Euclidean pose-cont disabled sigma emission')
        call cline%kill()

        call cline%set('objfun', 'cc')
        call cline%set('pose_cont', 'yes')
        call params%new(cline, silent=.true.)
        call assert_int(OBJFUN_CC, params%cc_objfun, &
            &'CC post-matcher pose-cont changed the outer objective')
        call assert_char('yes', trim(params%cc_emit_sigma), &
            &'CC post-matcher pose-cont did not enable cc_emit_sigma')
        call assert_true(pose_cont_sigma_is_enabled(params%cc_objfun, params%cc_emit_sigma), &
            &'CC post-matcher pose-cont disabled sigma emission')
        call cline%kill()

        call cline%set('objfun', 'cc')
        call cline%set('refine', 'pose_cont')
        call params%new(cline, silent=.true.)
        call assert_int(OBJFUN_CC, params%cc_objfun, &
            &'CC standalone pose-cont changed the outer objective')
        call assert_char('yes', trim(params%cc_emit_sigma), &
            &'CC standalone pose-cont did not enable cc_emit_sigma')
        call assert_true(pose_cont_sigma_is_enabled(params%cc_objfun, params%cc_emit_sigma), &
            &'CC standalone pose-cont disabled sigma emission')
        call cline%kill()

        call cline%set('objfun', 'cc')
        call cline%set('pose_cont', 'no')
        call params%new(cline, silent=.true.)
        call assert_int(OBJFUN_CC, params%cc_objfun, &
            &'ordinary CC changed the outer objective')
        call assert_char('no', trim(params%cc_emit_sigma), &
            &'ordinary CC unexpectedly enabled cc_emit_sigma')
        call assert_false(pose_cont_sigma_is_enabled(params%cc_objfun, params%cc_emit_sigma), &
            &'ordinary CC unexpectedly enabled pose-cont sigma emission')
        call cline%kill()
    end subroutine test_outer_objective_sigma_policy

    ! Identity is a valid initialized pose; explicit state/half metadata, not
    ! nonzero Euler coordinates, defines readiness for the standalone class.
    subroutine test_seed_contract()
        type(ori) :: seed

        call seed%set_euler([0., 0., 0.])
        call seed%set_shift([1.25, -0.75])
        call seed%set('state', 1.)
        call seed%set('eo', 0.)
        call seed%set('proj', 1.)
        call assert_true(pose_cont_seed_is_valid(seed), &
            &'standalone pose strategy rejected a valid identity seed')
        call seed%set('proj', 0.)
        call assert_false(pose_cont_seed_is_valid(seed), &
            &'standalone pose strategy accepted a missing projection seed')
        call seed%kill()
    end subroutine test_seed_contract

end module simple_strategy3D_pose_cont_tester
