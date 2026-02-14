!@descr: "validation" UI api (concrete implementation)
module simple_ui_api_validation
use simple_ui_api_modules
implicit none

type(ui_program), target :: model_validation
type(ui_program), target :: mini_stream
type(ui_program), target :: check_refpick

contains

    subroutine print_validation_programs(logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('VALIDATION:', C_UNDERLINED)
        write(logfhandle,'(A)') model_validation%name%to_char()
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') check_refpick%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_validation_programs

    subroutine new_check_refpick( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call check_refpick%new(&
        &'check_refpick',&                                           ! name
        &'validation of reference-based picking',&                   ! descr_short
        &'is a program for validation of reference-based picking',&  ! descr_long
        &'simple_exec',&                                             ! executable
        &.false.)                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call check_refpick%add_input(UI_IMG, 'filetab', 'file', 'List of files', 'List of files (*.mrcs) to process', 'e.g. mics.txt', .true., '')
        call check_refpick%add_input(UI_IMG, pickrefs, required_override=.true.)
        ! parameter input/output
        call check_refpick%add_input(UI_PARM, smpd,    required_override=.true.)
        call check_refpick%add_input(UI_PARM, pcontrast)
        call check_refpick%add_input(UI_PARM, kv,      required_override=.true.)
        call check_refpick%add_input(UI_PARM, cs,      required_override=.true.)
        call check_refpick%add_input(UI_PARM, fraca)
        ! alternative inputs
        ! <empty>
        ! search controls
        call check_refpick%add_input(UI_SRCH, 'nptcls_per_cls','num',   'Number of particles per class', 'Number of particles per class{200}', '# particles per class{200}', .false., 200.)
        call check_refpick%add_input(UI_SRCH, pick_roi)
        call check_refpick%add_input(UI_SRCH, particle_density)
        call check_refpick%add_input(UI_SRCH, nboxes_max)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call check_refpick%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('check_refpick', check_refpick, prgtab)
    end subroutine new_check_refpick

    subroutine new_mini_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mini_stream%new(&
        &'mini_stream',&                                ! name
        &'standalone mini_stream for a quick look',&    ! descr_short
        &'is a program for doing a standalone mini_stream for a quick look',&  ! descr_long
        &'simple_exec',&                                ! executable
        &.false.)                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call mini_stream%add_input(UI_IMG, 'filetab',    'file', 'List of files', 'List of files (*.mrcs) to process', 'e.g. mics.txt', .true., '')
        ! parameter input/output
        call mini_stream%add_input(UI_PARM, smpd, required_override=.true.)
        call mini_stream%add_input(UI_PARM, pcontrast)
        call mini_stream%add_input(UI_PARM, kv,   required_override=.true.)
        call mini_stream%add_input(UI_PARM, cs,   required_override=.true.)
        call mini_stream%add_input(UI_PARM, fraca)
        call mini_stream%add_input(UI_PARM, moldiam_max)
        ! alternative inputs
        ! <empty>
        ! search controls
        call mini_stream%add_input(UI_SRCH, 'nptcls_per_cls','num',   'Number of particles per class', 'Number of particles per class{200}', '# particles per class{200}', .false., 200.)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call mini_stream%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('mini_stream', mini_stream, prgtab)
    end subroutine new_mini_stream

    subroutine new_model_validation( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call model_validation%new(&
        &'model_validation', &                                                                        ! name
        &'Validation of atomic model',&                                                               ! descr_short
        &'is a program to validate the PDB atomic model given a 3D experimental density map in MRC',& ! descr long
        &'simple_exec',&                                                                              ! executable
        &.false.)                                                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call model_validation%add_input(UI_IMG, 'vol1', 'file', 'Experimental volume',  'Experimental volume',  'vol.mrc file', .true., '')
        call model_validation%add_input(UI_IMG, 'pdbfile', 'file', 'PDB input coordinates file', 'Input coordinates file in PDB format', 'PDB file e.g. molecule.pdb', .true., 'molecule.pdb')
        ! parameter input/output
        call model_validation%add_input(UI_PARM, smpd)
        call model_validation%add_input(UI_PARM, smpd_target)
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('model_validation', model_validation, prgtab)
    end subroutine new_model_validation

end module simple_ui_api_validation
