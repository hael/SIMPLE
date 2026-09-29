!@descr: refine3D_auto variant that owns the public Cartesian pose-cont workflow
module simple_commander_refine3D_pose_cont
use simple_commander_base, only: commander_base
use simple_commanders_refine3D, only: commander_refine3D_auto
use simple_cmdline, only: cmdline
use simple_error, only: simple_exception
use simple_parameters, only: parameters
use simple_refine3D_pose_cont_workflow, only: refine3D_pose_cont_workflow
implicit none
private

#include "simple_local_flags.inc"

public :: commander_refine3D_pose_cont

type, extends(commander_refine3D_auto) :: commander_refine3D_pose_cont
    type(refine3D_pose_cont_workflow), private :: workflow
contains
    procedure :: configure_workflow_defaults => configure_pose_cont_defaults
    procedure :: configure_registration_stage => configure_pose_cont_registration_stage
    procedure :: execute_main_stage => execute_pose_cont_main_stage
    procedure :: before_final_reconstruction => before_pose_cont_final_reconstruction
end type commander_refine3D_pose_cont

contains

    subroutine configure_pose_cont_defaults(self, cline, minits_default)
        class(commander_refine3D_pose_cont), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        integer, intent(in) :: minits_default
        integer :: maxits_user

        if (.not. cline%defined('objfun')) call cline%set('objfun', 'euclid')
        call cline%set('inpl_cont', 'no')
        call cline%set('pose_cont', 'no')
        call cline%set('pose_cont_route', 'joint')
        call cline%set('oritype', 'ptcl3D')
        call cline%set('projrec', 'no')
        if (cline%defined('maxits')) then
            maxits_user = cline%get_iarg('maxits')
            if (maxits_user < 1) THROW_HARD('maxits must be >= 1 for REFINE3D_POSE_CONT')
            if (cline%defined('minits')) then
                call cline%set('minits', min(maxits_user, max(1, cline%get_iarg('minits'))))
            else
                call cline%set('minits', min(maxits_user, minits_default))
            end if
        else if (cline%defined('minits')) then
            call cline%set('minits', max(minits_default, cline%get_iarg('minits')))
        else
            call cline%set('minits', minits_default)
        end if
    end subroutine configure_pose_cont_defaults

    subroutine configure_pose_cont_registration_stage(self, cline)
        class(commander_refine3D_pose_cont), intent(inout) :: self
        class(cmdline), intent(inout) :: cline

        call cline%set('inpl_cont', 'no')
        call cline%set('pose_cont', 'yes')
        call cline%set('pose_cont_route', 'joint')
    end subroutine configure_pose_cont_registration_stage

    subroutine execute_pose_cont_main_stage(self, params, cline, refine3D_executor)
        class(commander_refine3D_pose_cont), intent(inout) :: self
        type(parameters), intent(in) :: params
        class(cmdline), intent(inout) :: cline
        class(commander_base), intent(inout) :: refine3D_executor

        call self%workflow%new(params%pose_cont_mode)
        call self%workflow%execute_main(params, cline, refine3D_executor)
    end subroutine execute_pose_cont_main_stage

    subroutine before_pose_cont_final_reconstruction(self, cline, refine3D_executor)
        class(commander_refine3D_pose_cont), intent(inout) :: self
        class(cmdline), intent(inout) :: cline
        class(commander_base), intent(inout) :: refine3D_executor

        call self%workflow%execute_pre_final(cline, refine3D_executor)
        call self%workflow%kill
    end subroutine before_pose_cont_final_reconstruction

end module simple_commander_refine3D_pose_cont
