!@descr: module defining the user interfaces for resolution estimation programs in the simple_exec suite
module simple_ui_res
use simple_ui_modules
implicit none

type(ui_program), target :: fsc
type(ui_program), target :: fsc_area_score

contains

    subroutine construct_res_programs(prgtab)
        class(ui_hash), intent(inout) :: prgtab
        call new_fsc(prgtab)
        call new_fsc_area_score(prgtab)
    end subroutine construct_res_programs

    subroutine print_res_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('RESOLUTION ESTIMATION:', C_UNDERLINED)
        write(logfhandle,'(A)') fsc%name%to_char()
        write(logfhandle,'(A)') fsc_area_score%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_res_programs

    subroutine new_fsc( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! computer controls
        call fsc%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('fsc', fsc, prgtab)
    end subroutine new_fsc

    subroutine new_fsc_area_score( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call fsc_area_score%new(&
        &'fsc_area_score', &
        &'Calculate a conical FSC area ratio from two half maps',&
        &'is a program for calculating a CryoSPARC-like conical FSC area ratio from two half maps',&
        &'simple_exec',&
        &.false.)
        call fsc_area_score%add_input(UI_IMG, 'vol1', 'file', 'Odd volume',  'Odd volume',  'vol1.mrc file', .true., '')
        call fsc_area_score%add_input(UI_IMG, 'vol2', 'file', 'Even volume', 'Even volume', 'vol2.mrc file', .true., '')
        call fsc_area_score%add_input(UI_PARM, smpd)
        call fsc_area_score%add_input(UI_PARM, 'nspace', 'num', 'Number of cone axes', &
            &'Number of Fibonacci-sampled cone axes used for conical FSC area scoring{256}', '# cone axes{256}', .false., 256.)
        call fsc_area_score%add_input(UI_PARM, 'athres', 'num', 'Cone half-angle', &
            &'Cone half-angle in degrees for conical FSC area scoring{20}', 'degrees{20}', .false., 20.)
        call fsc_area_score%add_input(UI_FILT, lplim_crit)
        call fsc_area_score%add_input(UI_MASK, automsk)
        call fsc_area_score%add_input(UI_MASK, mskdiam)
        call fsc_area_score%add_input(UI_PARM, 'fbody', 'string', 'Output file body', &
            &'File body for fsc_area_score output tables', 'file body{fsc_area_score}', .false., 'fsc_area_score')
        call fsc_area_score%add_input(UI_COMP, nthr)
        call add_ui_program('fsc_area_score', fsc_area_score, prgtab)
    end subroutine new_fsc_area_score

end module simple_ui_res
