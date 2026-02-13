module simple_ui_api_other
use simple_ui_api_modules
implicit none

type(ui_program), target :: match_stacks
type(ui_program), target :: mkdir_
type(ui_program), target :: normalize_
type(ui_program), target :: print_ui_stream
type(ui_program), target :: split_
type(ui_program), target :: split_stack

contains

    subroutine register_simple_ui_other( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call prgtab%set_ui_program('match_stacks',    match_stacks)
        call prgtab%set_ui_program('mkdir_',          mkdir_)
        call prgtab%set_ui_program('normalize_',      normalize_)
        call prgtab%set_ui_program('print_ui_stream', print_ui_stream)
        call prgtab%set_ui_program('split_',          split_)
        call prgtab%set_ui_program('split_stack',     split_stack)
    end subroutine register_simple_ui_other

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_match_stacks
        ! PROGRAM SPECIFICATION
        call match_stacks%new(&
        &'match_stack',&                                              ! name
        &'Analysis of class averages with k-medoids',&                ! descr_short
        &'is a program for analyzing class averages with k-medoids',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.false.)                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call match_stacks%add_input(UI_IMG, stk,  required_override=.true.)
        call match_stacks%add_input(UI_IMG, stk2, required_override=.true.)
        ! parameter input/output
         ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call match_stacks%add_input(UI_SRCH, 'clust_crit', 'multi', 'Clustering criterion', 'Clustering criterion(sig|sig_clust|cc|res|hybrid){hybrid}',&
        &'(sig|sig_clust|cc|res|hybrid){hybrid}', .false., 'cc')
        ! filter controls
        call match_stacks%add_input(UI_FILT, hp, required_override=.true.)
        call match_stacks%add_input(UI_FILT, lp, required_override=.true.)
        ! mask controls
        call match_stacks%add_input(UI_MASK, mskdiam)
        ! computer controls
        call match_stacks%add_input(UI_COMP, nthr)
    end subroutine new_match_stacks

    subroutine new_mkdir_
        ! PROGRAM SPECIFICATION
        call mkdir_%new(&
        &'mkdir',&                                                       ! name
        &'Make directory',&                                              ! descr_short
        &'is a program for making an automatically numbered directory',& ! descr_long
        &'simple_exec',&                                                 ! executable
        &.false.)                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call mkdir_%add_input(UI_PARM, 'dir', 'dir', 'Name of directory to create', 'Name of directory name to create(e.g. X_dir/)', 'e.g. X_dir/', .true., ' ')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_mkdir_

    subroutine new_normalize
        ! PROGRAM SPECIFICATION
        call normalize_%new(&
        &'normalize',&                         ! name
        &'Normalize volume/stack',&            ! descr_short
        &'is a program for normalization of MRC or SPIDER stacks and volumes',&
        &'simple_exec',&                       ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call normalize_%add_input(UI_PARM, smpd)
        call normalize_%add_input(UI_PARM, 'norm',       'binary', 'Normalize',       'Statistical normalization: avg=zero, var=1(yes|no){no}',    '(yes|no){no}', .false., 'no')
        call normalize_%add_input(UI_PARM, 'noise_norm', 'binary', 'Noise normalize', 'Statistical normalization based on background(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! alternative inputs
        call normalize_%add_input(UI_ALT, 'stk',  'file', 'Stack to normalize',  'Stack of images to normalize', 'e.g. imgs.mrc', .false., '')
        call normalize_%add_input(UI_ALT, 'vol1', 'file', 'Volume to normalize', 'Volume to normalize',          'e.g. vol.mrc',  .false., '')
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call normalize_%add_input(UI_MASK, mskdiam)
        ! computer controls
        call normalize_%add_input(UI_COMP, nthr)
    end subroutine new_normalize

    subroutine new_split_
        call split_%new(&
        &'split',&                                   ! name
        &'Split stack into substacks',&              ! descr_short
        &'is a program for splitting a stack into evenly partitioned substacks',& ! descr_long
        &'simple_exec',&                             ! executable
        &.false.)                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call split_%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call split_%add_input(UI_PARM, smpd)
        ! computer controls
        call split_%add_input(UI_COMP, nparts)
    end subroutine new_split_


    subroutine new_split_stack
        ! PROGRAM SPECIFICATION
        call split_stack%new(&
        &'split_stack',&                                              ! name
        &'split stack in project',&                                   ! descr_short
        &'is a program for splitting a stack into nparts substacks',& ! descr_long
        &'simple_exec',&                                              ! executable
        &.true.)                                                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call split_stack%add_input(UI_PARM, 'nparts', 'num', 'Number of parts balanced splitting of the stack', '# parts', '# parts', .true., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        ! <empty>
    end subroutine new_split_stack

end module simple_ui_api_other
