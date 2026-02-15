!@descr: "single_ui_trajectory" UI api (concrete implementation)
module single_ui_trajectory
use simple_ui_modules
implicit none

type(ui_program), target :: extract_substk
type(ui_program), target :: graphene_subtr
type(ui_program), target :: import_trajectory
type(ui_program), target :: trajectory_denoise
type(ui_program), target :: trajectory_make_projavgs
type(ui_program), target :: trajectory_reconstruct3D
type(ui_program), target :: trajectory_swap_stack

contains

    subroutine construct_single_trajectory_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_extract_substk(prgtab)
        call new_graphene_subtr(prgtab)
        call new_import_trajectory(prgtab)
        call new_trajectory_denoise(prgtab)
        call new_trajectory_make_projavgs(prgtab)
        call new_trajectory_reconstruct3D(prgtab)
        call new_trajectory_swap_stack(prgtab)
    end subroutine construct_single_trajectory_programs

    subroutine print_single_trajectory_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('TRAJECTORY:', C_UNDERLINED)
        write(logfhandle,'(A)') extract_substk%name%to_char()
        write(logfhandle,'(A)') graphene_subtr%name%to_char()
        write(logfhandle,'(A)') import_trajectory%name%to_char()
        write(logfhandle,'(A)') trajectory_denoise%name%to_char()
        write(logfhandle,'(A)') trajectory_make_projavgs%name%to_char()
        write(logfhandle,'(A)') trajectory_reconstruct3D%name%to_char()
        write(logfhandle,'(A)') trajectory_swap_stack%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_trajectory_programs

    subroutine new_extract_substk( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call extract_substk%new(&
        &'extract_substk',&                                                                                            ! name
        &'extraction of a substack segment of time-series of metallic nanoparticles',&                                 ! descr_short
        &'is a shared-memory workflow for extraction of a substack segment of time-series of metallic nanoparticles',& ! descr_long
        &'single_exec',&                                                                                               ! executable
        &.true., gui_advanced=.false.)                                                                                 ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call extract_substk%add_input(UI_PARM, projfile)
        call extract_substk%add_input(UI_PARM, 'fromp', 'num', 'From index', 'Start index for stack copy', 'start index', .false., 1.0)
        call extract_substk%add_input(UI_PARM, 'top',   'num', 'To index', 'Stop index for stack copy', 'stop index', .false., 1.0)
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
        call add_ui_program('extract_substk', extract_substk, prgtab)
    end subroutine new_extract_substk

    subroutine new_graphene_subtr( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call graphene_subtr%new(&
        &'graphene_subtr',&                                  ! name
        &'Removes graphene Fourier peaks in time-series',&   ! descr_short
        &'Removes graphene Fourier peaks in time-series',&   ! descr_long
        &'single_exec',&                                     ! executable
        &.false., gui_advanced=.false.)                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call graphene_subtr%add_input(UI_IMG, stk_traj)
        call graphene_subtr%add_input(UI_IMG, stk_backgr)
        call graphene_subtr%add_input(UI_IMG, outstk)
        ! parameter input/output
        call graphene_subtr%add_input(UI_PARM, smpd)
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call graphene_subtr%add_input(UI_COMP, nthr)
        call add_ui_program('graphene_subtr', graphene_subtr, prgtab)
    end subroutine new_graphene_subtr

    subroutine new_import_trajectory( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call import_trajectory%new(&
        &'import_trajectory',&                    ! name
        &'Imports time-series particles stack',&            ! descr_short
        &'is a workflow for importing time-series data',&   ! descr_long
        &'single_exec',&                                    ! executable
        &.true., gui_advanced=.false.)                      ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call import_trajectory%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        call import_trajectory%add_input(UI_PARM, smpd)
        call import_trajectory%add_input(UI_PARM, deftab)
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
        call add_ui_program('import_trajectory', import_trajectory, prgtab)
    end subroutine new_import_trajectory

    subroutine new_trajectory_denoise( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_denoise%new(&
        &'trajectory_denoise',&                                       ! name
        &'kPCA-based denoising',&                                     ! descr_short
        &'is a program for kPCA-based denoising of an image stack',&  ! descr_long
        &'single_exec',&                                              ! executable
        &.false., gui_advanced=.false.)                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_denoise%add_input(UI_IMG, 'stk',  'file', 'Stack to denoise',  'Stack of images to denoise', 'e.g. stk.mrcs', .true., '')
        call trajectory_denoise%add_input(UI_IMG, outstk)
        ! parameter input/output
        call trajectory_denoise%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call trajectory_denoise%add_input(UI_FILT, 'neigs', 'num', 'Number of eigencomponents, corresponding to the number of classes in the stack', 'Number of eigencomponents, corresponding to the number of classes in the stack', '# eigenvecs', .false., 500.0)
        ! mask controls
        ! <empty>
        ! computer controls
        call trajectory_denoise%add_input(UI_COMP, nthr)
        call add_ui_program('trajectory_denoise', trajectory_denoise, prgtab)
    end subroutine new_trajectory_denoise

    subroutine new_trajectory_make_projavgs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_make_projavgs%new(&
        &'trajectory_make_projavgs',&                                                    ! name
        &'Align & average the first few frames of the time-series',&                     ! descr_short
        &'is a program for aligning & averaging the first few frames of the time-series&
        & to accomplish SNR enhancement for particle identification',&                   ! descr_long
        &'single_exec',&                                                                 ! executable
        &.true., gui_advanced=.false.)                                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call trajectory_make_projavgs%add_input(UI_PARM, nspace)
        call trajectory_make_projavgs%add_input(UI_PARM, 'athres', 'num', 'Angular threshold (degrees)', 'Angular threshold (degrees)', 'in degrees{10}', .false., 10.)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        call trajectory_make_projavgs%add_input(UI_MASK, mskdiam)
        ! computer controls
        call trajectory_make_projavgs%add_input(UI_COMP, nthr)
        call add_ui_program('trajectory_make_projavgs', trajectory_make_projavgs, prgtab)
    end subroutine new_trajectory_make_projavgs

    subroutine new_trajectory_reconstruct3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_reconstruct3D%new(&
        &'trajectory_reconstruct3D',&                                    ! name
        &'Time windowed 3D reconstruction from oriented particles',&     ! descr_long
        &'Time windowed 3D reconstruction from oriented particles',&
        &'single_exec',&                                                 ! executable
        &.true., gui_advanced=.false.)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call trajectory_reconstruct3D%add_input(UI_PARM, 'stepsz',  'num', 'Time window size (# frames){500}', 'Time window size (# frames) for windowed 3D rec{500}', 'give # frames',  .false., 500.)
        call trajectory_reconstruct3D%add_input(UI_PARM, 'fromp', 'num', 'From particle index', 'Start index for 3D reconstruction', 'start index', .false., 1.0)
        call trajectory_reconstruct3D%add_input(UI_PARM, 'top',   'num', 'To particle index', 'Stop index for 3D reconstruction', 'stop index', .false., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        call trajectory_reconstruct3D%add_input(UI_SRCH, pgrp)
        ! filter controls
        ! <empty>
        ! mask controls
        call trajectory_reconstruct3D%add_input(UI_MASK, mskdiam)
        call trajectory_reconstruct3D%add_input(UI_MASK, mskfile)
        ! computer controls
        call trajectory_reconstruct3D%add_input(UI_COMP, nparts, required_override=.false.)
        call trajectory_reconstruct3D%add_input(UI_COMP, nthr)
        call add_ui_program('trajectory_reconstruct3D', trajectory_reconstruct3D, prgtab)
    end subroutine new_trajectory_reconstruct3D

    subroutine new_trajectory_swap_stack( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call trajectory_swap_stack%new(&
        &'trajectory_swap_stack',&                                        ! name
        &'Substitutes stack into an existing project',&                   ! descr_short
        &'is a program for substituting stack into an existing project',& ! descr_long
        &'single_exec',&                                                  ! executable
        &.true., gui_advanced=.false.)                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call trajectory_swap_stack%add_input(UI_IMG, stk, required_override=.true.)
        ! parameter input/output
        ! <empty>
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
        call add_ui_program('trajectory_swap_stack', trajectory_swap_stack, prgtab)
    end subroutine new_trajectory_swap_stack

end module single_ui_trajectory
