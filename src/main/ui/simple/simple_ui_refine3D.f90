!@descr: module defining the user interfaces for 3D refinement programs in the simple_exec suite
module simple_ui_refine3D
use simple_ui_modules
implicit none

type(ui_program), target :: automask
type(ui_program), target :: postprocess
type(ui_program), target :: reconstruct3D
type(ui_program), target :: refine3D
type(ui_program), target :: refine3D_auto

contains

    subroutine construct_refine3D_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_automask(prgtab)
        call new_postprocess(prgtab)
        call new_reconstruct3D(prgtab)
        call new_refine3D(prgtab)
        call new_refine3D_auto(prgtab)
    end subroutine construct_refine3D_programs

    subroutine print_refine3D_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('REFINE 3D WORKFLOWS:', C_UNDERLINED)
        write(logfhandle,'(A)') automask%name%to_char()
        write(logfhandle,'(A)') postprocess%name%to_char()
        write(logfhandle,'(A)') reconstruct3D%name%to_char()
        write(logfhandle,'(A)') refine3D%name%to_char()
        write(logfhandle,'(A)') refine3D_auto%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_refine3D_programs

    subroutine new_automask( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call automask%new(&
        &'automask',&                                    ! name
        &'envelope masking',&                            ! descr_short
        &'is a program for automated envelope masking',& ! descr_long
        &'simple_exec',&                                 ! executable
        &.false.)                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call automask%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call automask%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call automask%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call automask%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 8.)
        ! mask controls
        call automask%add_input(UI_MASK, 'binwidth', 'num', 'Envelope binary layers width',&
        &'Binary layers grown for molecular envelope in pixels{1}', 'Molecular envelope binary layers width in pixels{1}', .false., 1.)
        call automask%add_input(UI_MASK, 'thres', 'num', 'Volume threshold',&
        &'Volume threshold for envelope mask generation', 'Volume threshold, give 0 if unknown', .false., 0.)
        call automask%add_input(UI_MASK, 'edge', 'num', 'Envelope mask soft edge',&
        &'Cosine edge size for softening molecular envelope in pixels{6}', '# pixels cosine edge{6}', .false., 6.)
        call automask%add_input(UI_MASK, automsk)
        ! computer controls
        call automask%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('automask', automask, prgtab)
    end subroutine new_automask

    subroutine new_postprocess( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call postprocess%new(&
        &'postprocess',&                                                      ! name
        &'Post-processing of volume',&                                        ! descr_short
        &'is a program for map post-processing. Use program volops to estimate the B-factor with the Guinier plot',& ! descr_long
        &'simple_exec',&                                                      ! executable
        &.true.)                                                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call postprocess%add_input(UI_PARM, 'state', 'num', 'State to postprocess', 'State to postprocess{1}', 'Input state{1}', .false., 1.0)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        call postprocess%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit for map filtering', 'Low-pass limit for map filtering', 'low-pass limit in Angstroms', .false., 20.)
        call postprocess%add_input(UI_FILT, bfac)
        call postprocess%add_input(UI_FILT, mirr)
        ! mask controls
        call postprocess%add_input(UI_MASK, mskdiam)
        ! computer controls
        call postprocess%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('postprocess', postprocess, prgtab)
    end subroutine new_postprocess

    subroutine new_reconstruct3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call reconstruct3D%new(&
        &'reconstruct3D',&                                               ! name
        &'3D reconstruction from oriented particles',&                   ! descr_short
        &'is a distributed workflow for reconstructing volumes from MRC and SPIDER stacks,&
        & given input orientations and state assignments. The algorithm is based on direct Fourier inversion&
        & with a Kaiser-Bessel (KB) interpolation kernel',&
        &'simple_exec',&                                                 ! executable
        &.true.)                                                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call reconstruct3D%add_input(UI_SRCH, pgrp)
        call reconstruct3D%add_input(UI_SRCH, frac)
        ! filter controls
        call reconstruct3D%add_input(UI_FILT, envfsc)
        call reconstruct3D%add_input(UI_FILT, wiener)
        ! mask controls
        call reconstruct3D%add_input(UI_MASK, mskdiam)
        call reconstruct3D%add_input(UI_MASK, mskfile)
        ! computer controls
        call reconstruct3D%add_input(UI_COMP, nparts, required_override=.false.)
        call reconstruct3D%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('reconstruct3D', reconstruct3D, prgtab)
    end subroutine new_reconstruct3D

    subroutine new_refine3D( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D%new(&
        &'refine3D',&                                                                               ! name
        &'3D refinement',&                                                                          ! descr_short
        &'is a distributed workflow for 3D refinement based on probabilistic projection matching',& ! descr_long
        &'simple_exec',&                                                                            ! executable
        &.true.,&                                                                                   ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "search,filter,mask,compute")                     ! GUI
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D%add_input(UI_IMG, 'vol1', 'file', 'Reference volume', 'Reference volume for creating polar 2D central &
        & sections for particle image matching', 'input volume e.g. vol.mrc', .false., 'vol1.mrc')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D%add_input(UI_SRCH, nspace, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, trs, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes', gui_submenu="search")
        call refine3D%add_input(UI_SRCH, maxits, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, update_frac, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, frac, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, pgrp, gui_submenu="search", gui_advanced=.false.)
        call refine3D%add_input(UI_SRCH, nstates, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, objfun, gui_submenu="search")
        call refine3D%add_input(UI_SRCH, 'refine', 'multi', 'Refinement mode', 'Refinement mode(shc_smpl|neigh|prob|prob_state){shc}', '(snhc|shc|neigh|shc_neigh){shc}',&
        &.false., 'shc', gui_submenu="search")
        call refine3D%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}', '(yes|no){no}', .false.,&
        &'no', gui_submenu="search")
        call refine3D%add_input(UI_SRCH, sigma_est, gui_submenu="search")
        ! filter controls
        call refine3D%add_input(UI_FILT, hp, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{30}', .false., 30., gui_submenu="filter")
        call refine3D%add_input(UI_FILT, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms', .false., 20., gui_submenu="filter")
        call refine3D%add_input(UI_FILT, 'lpstop', 'num', 'Low-pass limit for frequency limited refinement', 'Low-pass limit used to limit the resolution &
        &to avoid possible overfitting', 'low-pass limit in Angstroms', .false., 1.0, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, lplim_crit, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, lp_backgr, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, envfsc, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 12., gui_submenu="filter")
        call refine3D%add_input(UI_FILT, wiener, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, ml_reg, gui_submenu="filter")
        call refine3D%add_input(UI_FILT, combine_eo, gui_submenu="filter")
        ! mask controls
        call refine3D%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        call refine3D%add_input(UI_MASK, mskfile, gui_submenu="mask")
        call refine3D%add_input(UI_MASK, automsk, gui_submenu="mask")
        ! computer controls
        call refine3D%add_input(UI_COMP, nparts, required_override=.false., gui_submenu="compute", gui_advanced=.false.)
        call refine3D%add_input(UI_COMP, nthr,                              gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('refine3D', refine3D, prgtab)
    end subroutine new_refine3D

    subroutine new_refine3D_auto( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D_auto%new(&
        &'refine3D_auto',&                                                                          ! name
        &'automated single-state 3D refinement',&                                                   ! descr_short
        &'is an automated workflow for single-state 3D refinement based on probabilistic projection matching',& ! descr_long
        &'simple_exec',&                                                                            ! executable
        &.true.,&                                                                                   ! requires sp_project
        &gui_advanced=.false., gui_submenu_list = "search,filter,mask,compute")                     ! GUI
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D_auto%add_input(UI_SRCH, maxits,      required_override=.true., gui_submenu="search")
        call refine3D_auto%add_input(UI_SRCH, update_frac, required_override=.true., gui_submenu="search")
        call refine3D_auto%add_input(UI_SRCH, pgrp,                                  gui_submenu="search", gui_advanced=.false.)
        call refine3D_auto%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}', '(yes|no){no}', .false.,&
        &'no', gui_submenu="search")
        ! filter controls
        call refine3D_auto%add_input(UI_FILT, 'amsklp', 'num', 'Low-pass limit for envelope mask generation',&
        & 'Low-pass limit for envelope mask generation in Angstroms', 'low-pass limit in Angstroms', .false., 12., gui_submenu="filter")
        call refine3D_auto%add_input(UI_FILT, combine_eo, gui_submenu="filter")
        call refine3D_auto%add_input(UI_FILT, 'res_target', 'num', 'Resolution target (in A)',&
        & 'Resolutiuon target in Angstroms', 'Resolution target in Angstroms', .false., 3., gui_submenu="filter")
        ! mask controls
        call refine3D_auto%add_input(UI_MASK, mskdiam, gui_submenu="mask", gui_advanced=.false.)
        ! computer controls
        call refine3D_auto%add_input(UI_COMP, nparts, gui_submenu="compute", gui_advanced=.false.)
        call refine3D_auto%add_input(UI_COMP, nthr, gui_submenu="compute", gui_advanced=.false.)
        ! add to ui_hash
        call add_ui_program('refine3D_auto', refine3D_auto, prgtab)
    end subroutine new_refine3D_auto

end module simple_ui_refine3D
