!@descr: "single_ui_api_nano3D" UI api (concrete implementation)
module single_ui_api_nano3D
use simple_ui_api_modules
implicit none

type(ui_program), target :: autorefine3D_nano
type(ui_program), target :: refine3D_nano

contains

    subroutine print_single_nano3D_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('NANO 3D:', C_UNDERLINED)
        write(logfhandle,'(A)') autorefine3D_nano%name%to_char()
        write(logfhandle,'(A)') refine3D_nano%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_single_nano3D_programs

    subroutine new_autorefine3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call autorefine3D_nano%new(&
        &'autorefine3D_nano',&                                                            ! name
        &'auto 3D refinement of metallic nanoparticles',&                                 ! descr_short
        &'is a distributed workflow for automated 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! descr_long
        &'single_exec',&                                                                  ! executable
        &.true., gui_advanced=.false.)                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call autorefine3D_nano%add_input(UI_IMG, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .true., '')
        ! parameter input/output
        call autorefine3D_nano%add_input(UI_PARM, smpd)
        call autorefine3D_nano%add_input(UI_PARM, element)
        ! alternative inputs
        ! <empty>
        ! search controls
        call autorefine3D_nano%add_input(UI_SRCH, nspace)
        call autorefine3D_nano%add_input(UI_SRCH, trs)
        call autorefine3D_nano%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call autorefine3D_nano%add_input(UI_SRCH, 'maxits', 'num', 'Max iterations', 'Maximum number of iterations', 'Max # iterations{5}', .false., 5.)
        call autorefine3D_nano%add_input(UI_SRCH, pgrp)
        call autorefine3D_nano%add_input(UI_SRCH, nrestarts)
        ! filter controls
        call autorefine3D_nano%add_input(UI_FILT, hp)
        call autorefine3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5.)
        call autorefine3D_nano%add_input(UI_FILT, 'lp', 'num', 'Initial low-pass limit', 'Initial low-pass limit', 'low-pass limit in Angstroms{1.5}', .true., 1.5)
        ! mask controls
        call autorefine3D_nano%add_input(UI_MASK, mskdiam)
        ! computer controls
        call autorefine3D_nano%add_input(UI_COMP, nthr)
        call autorefine3D_nano%add_input(UI_COMP, script)
        call add_ui_program('autorefine3D_nano', autorefine3D_nano, prgtab)
    end subroutine new_autorefine3D_nano

    subroutine new_refine3D_nano( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call refine3D_nano%new(&
        &'refine3D_nano',&                                                                                                    ! name
        &'3D refinement of metallic nanoparticles',&                                                                          ! descr_short
        &'is a distributed workflow for 3D refinement of metallic nanoparticles based on probabilistic projection matching',& ! descr_long
        &'single_exec',&                                                                                                      ! executable
        &.true., gui_advanced=.false.)                                                                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call refine3D_nano%add_input(UI_IMG, 'vol1', 'file', 'FCC reference volume', 'FCC lattice reference volume for creating polar 2D central &
        & sections for nanoparticle image matching', 'input volume e.g. vol.mrc', .false., '')
        call refine3D_nano%add_input(UI_IMG, 'vol_odd',  'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .false., '')
        call refine3D_nano%add_input(UI_IMG, 'vol_even', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .false., '')
        ! parameter input/output
        ! <empty>
        ! alternative inputs
        ! <empty>
        ! search controls
        call refine3D_nano%add_input(UI_SRCH, nspace)
        call refine3D_nano%add_input(UI_SRCH, trs)
        call refine3D_nano%add_input(UI_SRCH, 'center', 'binary', 'Center reference volume(s)', 'Center reference volume(s) by their &
        &center of gravity and map shifts back to the particles(yes|no){yes}', '(yes|no){yes}', .false., 'yes')
        call refine3D_nano%add_input(UI_SRCH, maxits)
        call refine3D_nano%add_input(UI_SRCH, update_frac)
        call refine3D_nano%add_input(UI_SRCH, frac)
        call refine3D_nano%add_input(UI_SRCH, pgrp)
        call refine3D_nano%add_input(UI_SRCH, 'continue', 'binary', 'Continue previous refinement', 'Continue previous refinement(yes|no){no}', '(yes|no){no}', .false., 'no')
        ! filter controls
        call refine3D_nano%add_input(UI_FILT, hp)
        call refine3D_nano%add_input(UI_FILT, 'cenlp', 'num', 'Centering low-pass limit', 'Limit for low-pass filter used in binarisation &
        &prior to determination of the center of gravity of the reference volume(s) and centering', 'centering low-pass limit in &
        &Angstroms{5}', .false., 5.)
        call refine3D_nano%add_input(UI_FILT, 'lp', 'num', 'Static low-pass limit', 'Static low-pass limit', 'low-pass limit in Angstroms{1.0}', .false., 1.)
        ! mask controls
        call refine3D_nano%add_input(UI_MASK, mskdiam)
        call refine3D_nano%add_input(UI_MASK, mskfile)
        ! computer controls
        call refine3D_nano%add_input(UI_COMP, nparts, required_override=.false.)
        call refine3D_nano%add_input(UI_COMP, nthr)
        call add_ui_program('refine3D_nano', refine3D_nano, prgtab)
    end subroutine new_refine3D_nano

end module single_ui_api_nano3D
