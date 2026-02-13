!@descr: "resolution" UI api (concrete implementation)
module simple_ui_api_resolution
use simple_ui_api_modules
implicit none

type(ui_program), target :: fsc
type(ui_program), target :: clin_fsc

contains

    subroutine register_simple_ui_resolution(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call add_ui_program('fsc',      fsc,      prgtab)
        call add_ui_program('clin_fsc', clin_fsc, prgtab)
    end subroutine register_simple_ui_resolution

! ============================================================
! Constructors moved from simple_user_interface.f90
! ============================================================

    subroutine new_clin_fsc
        ! PROGRAM SPECIFICATION
        call clin_fsc%new(&
        &'clin_fsc', &                                                          ! name
        &'Calculate FSC between the two input volumes using common lines',&     ! descr_short
        &'is a program for calculating the FSC between the two input volumes',& ! descr_long
        &'simple_exec',&                                                        ! executable
        &.false.)                                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call clin_fsc%add_input(UI_PARM, projfile)
        call clin_fsc%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        call clin_fsc%add_input(UI_SRCH, nspace)
        call clin_fsc%add_input(UI_SRCH, pgrp)
        ! filter controls
        call clin_fsc%add_input(UI_FILT, lp, required_override=.true.)
        ! mask controls
        call clin_fsc%add_input(UI_MASK, mskdiam)
        call clin_fsc%add_input(UI_MASK, mskfile)
        ! computer controls
        call clin_fsc%add_input(UI_COMP, nthr)
    end subroutine new_clin_fsc


    subroutine new_fsc
        ! PROGRAM SPECIFICATION
        call fsc%new(&
        &'fsc', &                                                               ! name
        &'Calculate FSC between the two input volumes',&                        ! descr_short
        &'is a program for calculating the FSC between the two input volumes',& ! descr_long
        &'simple_exec',&                                                        ! executable
        &.false.)                                                               ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call fsc%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call fsc%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        ! parameter input/output
        call fsc%add_input(UI_PARM, smpd)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        hp%required = .false.
        lp%required = .false.
        call fsc%add_input(UI_FILT, hp)
        call fsc%add_input(UI_FILT, lp)
        ! mask controls
        call fsc%add_input(UI_MASK, mskdiam)
        call fsc%add_input(UI_MASK, mskfile)
        ! computer controls
        call fsc%add_input(UI_COMP, nthr)
    end subroutine new_fsc


end module simple_ui_api_resolution
