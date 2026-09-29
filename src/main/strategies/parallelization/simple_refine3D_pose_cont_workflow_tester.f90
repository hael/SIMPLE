!@descr: unit tests for refine3D pose_cont commander and workflow stage policy
module simple_refine3D_pose_cont_workflow_tester
use simple_commander_base, only: commander_base
use simple_commander_refine3D_pose_cont, only: commander_refine3D_pose_cont
use simple_cmdline, only: cmdline
use simple_fileio, only: append2basename
use simple_parameters, only: parameters
use simple_refine3D_pose_cont_workflow, only: refine3D_pose_cont_workflow
use simple_sp_project, only: sp_project
use simple_string, only: string
use simple_syslib, only: del_file, file_exists
use simple_test_utils
implicit none
private
public :: run_all_refine3D_pose_cont_workflow_tests

integer, parameter :: MAX_RECORDED_STAGES = 3
character(len=*), parameter :: PROJECT_FILE = 'tmp_pose_cont_workflow.simple'
character(len=*), parameter :: SIGMA_FILE = 'tmp_pose_cont_workflow_sigma.bin'

type, extends(commander_base) :: recording_refine3D_executor
    integer :: ncalls = 0
    character(len=16) :: refine(MAX_RECORDED_STAGES) = ''
    character(len=3) :: pose_cont(MAX_RECORDED_STAGES) = ''
    character(len=3) :: inpl_cont(MAX_RECORDED_STAGES) = ''
    character(len=16) :: route(MAX_RECORDED_STAGES) = ''
    integer :: nsample(MAX_RECORDED_STAGES) = 0
    integer :: startit(MAX_RECORDED_STAGES) = 0
contains
    procedure :: execute => record_refine3D_stage
end type recording_refine3D_executor

contains

    subroutine run_all_refine3D_pose_cont_workflow_tests()
        write (*, '(A)') '**** running all refine3D pose_cont workflow tests ****'
        write (*, '(A)') 'test_registration_stage_policy'
        call test_registration_stage_policy()
        write (*, '(A)') 'test_post_matcher_stage_policy'
        call test_post_matcher_stage_policy()
        write (*, '(A)') 'test_standalone_final_stage_policy'
        call test_standalone_final_stage_policy()
    end subroutine run_all_refine3D_pose_cont_workflow_tests

    subroutine test_registration_stage_policy()
        type(commander_refine3D_pose_cont) :: commander
        type(cmdline) :: cline
        type(string) :: value

        call cline%set('refine', 'greedy')
        call cline%set('pose_cont', 'no')
        call cline%set('inpl_cont', 'yes')
        call cline%set('pose_cont_route', 'shift_then_joint')
        call commander%configure_registration_stage(cline)
        value = cline%get_carg('refine')
        call assert_char('greedy', value%to_char(), &
            &'registration changed the global greedy search')
        value = cline%get_carg('pose_cont')
        call assert_char('yes', value%to_char(), &
            &'registration did not enable pose_cont polishing')
        value = cline%get_carg('inpl_cont')
        call assert_char('no', value%to_char(), &
            &'registration did not disable inpl_cont')
        value = cline%get_carg('pose_cont_route')
        call assert_char('joint', value%to_char(), &
            &'registration did not select the joint route')
        call cline%kill()
    end subroutine test_registration_stage_policy

    subroutine test_post_matcher_stage_policy()
        type(refine3D_pose_cont_workflow) :: workflow
        type(recording_refine3D_executor) :: executor
        type(parameters) :: params
        type(cmdline) :: cline

        params%maxits = 4
        call cline%set('refine', 'prob_neigh')
        call workflow%new('post_matcher')
        call workflow%execute_main(params, cline, executor)
        call assert_int(1, executor%ncalls, 'post_matcher executes one main refinement call')
        call assert_char('prob_neigh', executor%refine(1), &
            &'post_matcher changed the probabilistic matcher')
        call assert_char('yes', executor%pose_cont(1), &
            &'post_matcher did not enable pose_cont polishing')
        call assert_char('no', executor%inpl_cont(1), &
            &'post_matcher did not disable inpl_cont')
        call assert_char('joint', executor%route(1), &
            &'post_matcher did not select the joint route')
        call workflow%kill()
        call cline%kill()
    end subroutine test_post_matcher_stage_policy

    subroutine test_standalone_final_stage_policy()
        type(refine3D_pose_cont_workflow) :: workflow
        type(recording_refine3D_executor) :: executor
        type(parameters) :: params
        type(cmdline) :: cline
        type(string) :: checkpoint_project, checkpoint_sigma

        call create_standalone_fixture()
        params%maxits = 3
        call cline%set('refine', 'prob_neigh')
        call cline%set('projfile', PROJECT_FILE)
        call workflow%new('standalone_final')
        call workflow%execute_main(params, cline, executor)
        call assert_int(1, executor%ncalls, 'standalone_final executes one probabilistic main call')
        call assert_char('prob_neigh', executor%refine(1), &
            &'standalone_final changed the probabilistic main matcher')
        call assert_char('no', executor%pose_cont(1), &
            &'standalone_final enabled pose_cont inside its main loop')
        call assert_char('no', executor%inpl_cont(1), &
            &'standalone_final main loop did not disable inpl_cont')
        call assert_char('joint', executor%route(1), &
            &'standalone_final main loop did not retain the joint route')

        call workflow%execute_pre_final(cline, executor)
        call assert_int(2, executor%ncalls, &
            &'standalone_final did not execute exactly one terminal call')
        call assert_char('pose_cont', executor%refine(2), &
            &'standalone_final terminal call did not select refine=pose_cont')
        call assert_char('no', executor%pose_cont(2), &
            &'standalone terminal call recursively enabled post-matcher pose_cont')
        call assert_char('no', executor%inpl_cont(2), &
            &'standalone terminal call did not disable inpl_cont')
        call assert_char('joint', executor%route(2), &
            &'standalone terminal call did not select the joint route')
        call assert_int(3, executor%nsample(2), &
            &'standalone terminal call did not select every active particle')
        call assert_int(4, executor%startit(2), &
            &'standalone terminal call did not immediately follow the main loop')

        checkpoint_project = append2basename(string(PROJECT_FILE), '_pre_pose_cont')
        checkpoint_sigma = append2basename(string(SIGMA_FILE), '_pre_pose_cont')
        call assert_true(file_exists(checkpoint_project), &
            &'standalone_final did not preserve the preterminal project')
        call assert_true(file_exists(checkpoint_sigma), &
            &'standalone_final did not preserve the preterminal sigma state')
        call workflow%kill()
        call cline%kill()
        call remove_standalone_fixture(checkpoint_project, checkpoint_sigma)
    end subroutine test_standalone_final_stage_policy

    subroutine record_refine3D_stage(self, cline)
        class(recording_refine3D_executor), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        type(string) :: value
        integer :: call_index, endit

        self%ncalls = self%ncalls + 1
        call_index = self%ncalls
        if (call_index > MAX_RECORDED_STAGES) return
        if (cline%defined('refine')) then
            value = cline%get_carg('refine')
            self%refine(call_index) = value%to_char()
        end if
        if (cline%defined('pose_cont')) then
            value = cline%get_carg('pose_cont')
            self%pose_cont(call_index) = value%to_char()
        end if
        if (cline%defined('inpl_cont')) then
            value = cline%get_carg('inpl_cont')
            self%inpl_cont(call_index) = value%to_char()
        end if
        if (cline%defined('pose_cont_route')) then
            value = cline%get_carg('pose_cont_route')
            self%route(call_index) = value%to_char()
        end if
        if (cline%defined('nsample')) self%nsample(call_index) = cline%get_iarg('nsample')
        if (cline%defined('startit')) self%startit(call_index) = cline%get_iarg('startit')
        if (cline%defined('startit')) then
            endit = cline%get_iarg('startit')
        else
            endit = cline%get_iarg('maxits')
        end if
        call cline%set('endit', endit)
    end subroutine record_refine3D_stage

    subroutine create_standalone_fixture()
        type(sp_project) :: project
        type(string) :: checkpoint_project, checkpoint_sigma
        integer :: unit, ios

        checkpoint_project = append2basename(string(PROJECT_FILE), '_pre_pose_cont')
        checkpoint_sigma = append2basename(string(SIGMA_FILE), '_pre_pose_cont')
        call remove_standalone_fixture(checkpoint_project, checkpoint_sigma)
        call project%os_ptcl2D%new(4, is_ptcl=.true.)
        call project%os_ptcl3D%new(4, is_ptcl=.true.)
        call project%os_ptcl2D%set(1, 'state', 1.)
        call project%os_ptcl2D%set(2, 'state', 1.)
        call project%os_ptcl2D%set(3, 'state', 0.)
        call project%os_ptcl2D%set(4, 'state', 2.)
        call project%os_ptcl3D%set(1, 'state', 1.)
        call project%os_ptcl3D%set(2, 'state', 1.)
        call project%os_ptcl3D%set(3, 'state', 0.)
        call project%os_ptcl3D%set(4, 'state', 2.)
        call project%update_projinfo(string(PROJECT_FILE))
        open (newunit=unit, file=SIGMA_FILE, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'standalone workflow sigma fixture opens')
        if (ios == 0) then
            write (unit, '(A)') 'pose-cont workflow sigma fixture'
            close (unit)
        end if
        call project%set_sigma2_state_path(string(SIGMA_FILE))
        call project%write(string(PROJECT_FILE))
        call project%kill()
    end subroutine create_standalone_fixture

    subroutine remove_standalone_fixture(checkpoint_project, checkpoint_sigma)
        type(string), intent(in) :: checkpoint_project, checkpoint_sigma
        call del_file(PROJECT_FILE)
        call del_file(SIGMA_FILE)
        call del_file(checkpoint_project)
        call del_file(checkpoint_sigma)
    end subroutine remove_standalone_fixture

end module simple_refine3D_pose_cont_workflow_tester
