!@descr: the pose student's thin commanders: learn_pose_student (the teacher's rows and the training) and predict_poses (the mode, the cone and the confidence of unseen particles)
module simple_commanders_pose_student
use simple_commander_base,        only: commander_base
use simple_cmdline,               only: cmdline
use simple_strategy_pose_student, only: run_learn_pose_student, run_predict_poses
use simple_pose_student,          only: NSPACE_STUDENT
implicit none
private
public :: commander_learn_pose_student, commander_predict_poses
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_learn_pose_student
  contains
    procedure :: execute => exec_learn_pose_student
end type commander_learn_pose_student

type, extends(commander_base) :: commander_predict_poses
  contains
    procedure :: execute => exec_predict_poses
end type commander_predict_poses

contains

    !> The teacher's rows of the project's selected particles at the student grid (nspace) and band (lp) against the
    !! references vol1, with their resolvability and their features, and the student trained on them, written as the
    !! model fname with its weights beside it (simple_private_exec). The dense table's mode and objective are set
    !! here: the student has no key of its own
    subroutine exec_learn_pose_student( self, cline )
        class(commander_learn_pose_student), intent(inout) :: self
        class(cmdline),                      intent(inout) :: cline
        real, parameter :: TRS_DEFAULT = 5.  ! the shift polish's box (pixels)
        if( .not. cline%defined('mkdir')  ) call cline%set('mkdir',  'no')
        if( .not. cline%defined('nspace') ) call cline%set('nspace', NSPACE_STUDENT)
        if( .not. cline%defined('trs')    ) call cline%set('trs',    TRS_DEFAULT)
        if( .not. cline%defined('fname')  ) call cline%set('fname',  'pose_student.txt')
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'euclid')
        call cline%set('refine',  'prob')
        call cline%set('center',  'no')
        call run_learn_pose_student(cline)
    end subroutine exec_learn_pose_student

    !> The student's prediction for the project's selected particles under the model infile: per particle the mode
    !! direction, the cone holding 99% of the posterior and the confidence, to the text table beside fname; the
    !! project is not written (simple_private_exec)
    subroutine exec_predict_poses( self, cline )
        class(commander_predict_poses), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        if( .not. cline%defined('mkdir') ) call cline%set('mkdir', 'no')
        if( .not. cline%defined('fname') ) call cline%set('fname', 'pose_student.txt')
        call cline%set('oritype', 'ptcl3D')
        call cline%set('objfun',  'euclid')
        call run_predict_poses(cline)
    end subroutine exec_predict_poses

end module simple_commanders_pose_student
