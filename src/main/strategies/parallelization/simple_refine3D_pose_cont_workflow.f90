!@descr: stage lifecycle for the automated refine3D Cartesian-pose workflow
module simple_refine3D_pose_cont_workflow
use simple_commander_base, only: commander_base
use simple_cmdline, only: cmdline
use simple_defs, only: logfhandle
use simple_error, only: simple_exception
use simple_fileio, only: append2basename, simple_copy_file
use simple_parameters, only: parameters
use simple_sp_project, only: sp_project
use simple_string, only: string
use simple_syslib, only: file_exists
implicit none
private
#include "simple_local_flags.inc"

public :: refine3D_pose_cont_workflow

type :: refine3D_pose_cont_workflow
    character(len=16), private :: mode = 'off'
    integer, private :: next_iter = 1
contains
    procedure :: new
    procedure :: kill
    procedure :: execute_main
    procedure :: execute_pre_final
    procedure, private :: execute_main_stage
    procedure, private :: execute_single_iteration
    procedure, private :: prepare_single_iteration
    procedure, private :: commit_stage_iteration
    procedure, private :: write_preterminal_checkpoint
end type refine3D_pose_cont_workflow

contains

    subroutine new(self, mode)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        character(len=*), intent(in) :: mode
        call self%kill
        select case (trim(mode))
        case ('post_matcher', 'standalone_final')
            self%mode = trim(mode)
        case default
            THROW_HARD('unsupported refine3D pose-cont workflow mode')
        end select
    end subroutine new

    subroutine kill(self)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        self%mode = 'off'
        self%next_iter = 1
    end subroutine kill

    subroutine execute_main(self, params, parent_cline, refine3D_executor)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        type(parameters), intent(in) :: params
        class(cmdline), intent(inout) :: parent_cline
        class(commander_base), intent(inout) :: refine3D_executor

        select case (trim(self%mode))
        case ('post_matcher')
            call self%execute_main_stage(params, parent_cline, refine3D_executor, .true.)
        case ('standalone_final')
            call self%execute_main_stage(params, parent_cline, refine3D_executor, .false.)
        case default
            THROW_HARD('refine3D pose-cont workflow is not initialized')
        end select
    end subroutine execute_main

    subroutine execute_pre_final(self, parent_cline, refine3D_executor)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        class(cmdline), intent(inout) :: parent_cline
        class(commander_base), intent(inout) :: refine3D_executor

        select case (trim(self%mode))
        case ('post_matcher')
            return
        case ('standalone_final')
            call self%write_preterminal_checkpoint(parent_cline)
            call self%execute_single_iteration(parent_cline, refine3D_executor, &
                &'pose_cont', 'terminal')
        case default
            THROW_HARD('refine3D pose-cont workflow is not initialized')
        end select
    end subroutine execute_pre_final

    subroutine execute_main_stage(self, params, parent_cline, refine3D_executor, pose_cont)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        type(parameters), intent(in) :: params
        class(cmdline), intent(inout) :: parent_cline
        class(commander_base), intent(inout) :: refine3D_executor
        logical, intent(in) :: pose_cont
        type(cmdline) :: stage_cline

        stage_cline = parent_cline
        call stage_cline%set('inpl_cont', 'no')
        call stage_cline%set('pose_cont_route', 'joint')
        if (pose_cont) then
            call stage_cline%set('pose_cont', 'yes')
        else
            call stage_cline%set('pose_cont', 'no')
        end if
        if (trim(self%mode) == 'standalone_final') call stage_cline%set('keepvol', 'yes')
        call stage_cline%set('maxits', params%maxits)
        call refine3D_executor%execute(stage_cline)
        call self%commit_stage_iteration(parent_cline, stage_cline)
        call stage_cline%kill
    end subroutine execute_main_stage

    subroutine execute_single_iteration(self, parent_cline, refine3D_executor, refine_mode, stage_label)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        class(cmdline), intent(inout) :: parent_cline
        class(commander_base), intent(inout) :: refine3D_executor
        character(len=*), intent(in) :: refine_mode, stage_label
        type(cmdline) :: stage_cline

        write (logfhandle, '(A,I0,A,A)') '>>> REFINE3D_POSE_CONT CARTESIAN STAGE ITERATION ', &
            &self%next_iter, ': ', trim(stage_label)
        stage_cline = parent_cline
        call self%prepare_single_iteration(stage_cline, refine_mode, self%next_iter)
        call refine3D_executor%execute(stage_cline)
        call self%commit_stage_iteration(parent_cline, stage_cline)
        call stage_cline%kill
    end subroutine execute_single_iteration

    subroutine prepare_single_iteration(self, stage_cline, refine_mode, iter)
        class(refine3D_pose_cont_workflow), intent(in) :: self
        class(cmdline), intent(inout) :: stage_cline
        character(len=*), intent(in) :: refine_mode
        integer, intent(in) :: iter
        type(sp_project) :: sampling_project
        integer :: nactive

        if (trim(self%mode) == 'off') THROW_HARD('pose-cont stage preparation requires an initialized workflow')
        call stage_cline%set('prg', 'refine3D')
        call stage_cline%set('refine', refine_mode)
        call stage_cline%set('inpl_cont', 'no')
        call stage_cline%set('pose_cont', 'no')
        call stage_cline%set('pose_cont_route', 'joint')
        call stage_cline%set('trail_rec', 'no')
        call stage_cline%set('maxits', 1)
        call stage_cline%set('minits', 1)
        call stage_cline%set('startit', iter)
        call stage_cline%set('which_iter', iter)
        call stage_cline%set('extr_iter', iter)
        call stage_cline%delete('update_frac')
        call sampling_project%read(stage_cline%get_carg('projfile'))
        nactive = sampling_project%count_state_gt_zero()
        call sampling_project%kill
        if (nactive < 1) THROW_HARD('standalone pose-cont stage has no active particles')
        call stage_cline%set('nsample', nactive)
        call stage_cline%delete('endit')
        call stage_cline%delete('continue')
    end subroutine prepare_single_iteration

    subroutine commit_stage_iteration(self, parent_cline, stage_cline)
        class(refine3D_pose_cont_workflow), intent(inout) :: self
        class(cmdline), intent(inout) :: parent_cline
        class(cmdline), intent(in) :: stage_cline
        integer :: last_iter

        if (.not. stage_cline%defined('endit')) &
            &THROW_HARD('refine3D pose-cont child stage did not report its final iteration')
        last_iter = stage_cline%get_iarg('endit')
        call parent_cline%set('endit', last_iter)
        self%next_iter = last_iter + 1
    end subroutine commit_stage_iteration

    !> Preserve a recoverable project and canonical sigma state immediately
    !! before the terminal Cartesian pass. keepvol=yes retains the preceding
    !! probabilistic stage's iterative references.
    subroutine write_preterminal_checkpoint(self, parent_cline)
        class(refine3D_pose_cont_workflow), intent(in) :: self
        class(cmdline), intent(in) :: parent_cline
        type(sp_project) :: checkpoint_project
        type(string) :: project_path, checkpoint_path, sigma_path, sigma_checkpoint
        logical :: sigma_found

        if (trim(self%mode) /= 'standalone_final') &
            &THROW_HARD('preterminal checkpoint is only valid for standalone_final')
        if (.not. parent_cline%defined('projfile')) THROW_HARD('standalone_final checkpoint requires projfile')
        project_path = parent_cline%get_carg('projfile')
        checkpoint_path = append2basename(project_path, '_pre_pose_cont')
        call checkpoint_project%read(project_path)
        call checkpoint_project%get_sigma2_state_path(sigma_path, sigma_found)
        if (.not. sigma_found) THROW_HARD('standalone_final checkpoint has no registered canonical sigma state')
        if (.not. file_exists(sigma_path)) THROW_HARD('standalone_final canonical sigma state does not exist')
        sigma_checkpoint = append2basename(sigma_path, '_pre_pose_cont')
        call simple_copy_file(sigma_path, sigma_checkpoint)
        call checkpoint_project%set_sigma2_state_path(sigma_checkpoint)
        call checkpoint_project%write(checkpoint_path)
        call checkpoint_project%kill
        write (logfhandle, '(A,A)') '>>> REFINE3D_POSE_CONT PRETERMINAL PROJECT: ', checkpoint_path%to_char()
        write (logfhandle, '(A,A)') '>>> REFINE3D_POSE_CONT PRETERMINAL SIGMA2: ', sigma_checkpoint%to_char()
    end subroutine write_preterminal_checkpoint

end module simple_refine3D_pose_cont_workflow
