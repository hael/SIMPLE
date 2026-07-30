!@descr: module defining the user interfaces for validation programs in the simple_exec suite
module simple_ui_validate
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('validate', 'Validation', 150)
type(ui_program), target :: check_refpick
type(ui_program), target :: mini_stream
type(ui_program), target :: model_validate

contains

    subroutine construct_validate_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_check_refpick(prgtab)
        call new_mini_stream(prgtab)
        call new_model_validate(prgtab)
    end subroutine construct_validate_programs

    subroutine print_validate_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('VALIDATION:', C_UNDERLINED)
        write(logfhandle,'(A)') check_refpick%name%to_char()
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') model_validate%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_validate_programs

    subroutine new_check_refpick( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call check_refpick%new(&
        &'check_refpick',&                                           ! name
        &'validation of reference-based picking',&                   ! summary
        &'is a program for validation of reference-based picking',&  ! help
        &'simple_exec',&                                             ! executable
        &.false.)                                                    ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call check_refpick%add_input(UI_FILE, 'filetab', 'file', 'List of files', 'List of files (*.mrcs) to process', 'e.g. mics.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call check_refpick%add_input(UI_IMG, pickrefs, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call check_refpick%add_input(UI_PARM, smpd,    required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call check_refpick%add_input(UI_PARM, pcontrast, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_PARM, kv,      required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call check_refpick%add_input(UI_PARM, cs,      required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call check_refpick%add_input(UI_PARM, fraca, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_PARM, fit_phshift, &
        &visibility=UI_VIS_DEVELOPER)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call check_refpick%add_input(UI_SRCH, 'nptcls_per_cls','num',   'Number of particles per class', 'Number of particles per class{200}', '# particles per class{200}', .false., 200., &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, pick_roi, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, particle_density, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, nboxes_max, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, phshift_min, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, phshift_max, &
        &visibility=UI_VIS_DEVELOPER)
        call check_refpick%add_input(UI_SRCH, phshift_step, &
        &visibility=UI_VIS_DEVELOPER)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call check_refpick%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('check_refpick', check_refpick, prgtab, UI_CATEGORY)
    end subroutine new_check_refpick

    subroutine new_mini_stream( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call mini_stream%new(&
        &'mini_stream',&                                ! name
        &'standalone mini_stream for a quick look',&    ! summary
        &'is a program for doing a standalone mini_stream for a quick look',&  ! help
        &'simple_exec',&                                ! executable
        &.false., &
        &visibility=UI_VIS_ADVANCED)                                       ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call mini_stream%add_input(UI_FILE, 'filetab',    'file', 'List of files', 'List of files (*.mrcs) to process', 'e.g. mics.txt', .true., '', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call mini_stream%add_input(UI_PARM, smpd, required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call mini_stream%add_input(UI_PARM, pcontrast, &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_PARM, kv,   required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call mini_stream%add_input(UI_PARM, cs,   required_override=.true., &
        &visibility=UI_VIS_STANDARD)
        call mini_stream%add_input(UI_PARM, fraca, &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_PARM, moldiam_max, &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_PARM, fit_phshift, &
        &visibility=UI_VIS_ADVANCED)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        call mini_stream%add_input(UI_SRCH, 'nptcls_per_cls','num',   'Number of particles per class', 'Number of particles per class{200}', '# particles per class{200}', .false., 200., &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_SRCH, phshift_min, &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_SRCH, phshift_max, &
        &visibility=UI_VIS_ADVANCED)
        call mini_stream%add_input(UI_SRCH, phshift_step, &
        &visibility=UI_VIS_ADVANCED)
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call mini_stream%add_input(UI_COMP, nthr, &
        &visibility=UI_VIS_STANDARD)
        ! add to ui_hash
        call add_ui_program('mini_stream', mini_stream, prgtab, UI_CATEGORY)
    end subroutine new_mini_stream

    subroutine new_model_validate( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        ! PROGRAM SPECIFICATION
        call model_validate%new(&
        &'model_validate', &                                                                          ! name
        &'Validate an atomic model against an experimental density map',& ! summary
        &'is a program to validate the PDB atomic model given a 3D experimental density map in MRC',& ! descr long
        &'simple_exec',&                                                                              ! executable
        &.false.)                                                                                     ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call model_validate%add_input(UI_IMG, 'vol1', 'file', 'Experimental volume',  'Experimental volume',  'vol.mrc file', .true., '', &
        &visibility=UI_VIS_STANDARD)
        call model_validate%add_input(UI_FILE, 'pdbfile', 'file', 'PDB input coordinates file', 'Input coordinates file in PDB format', 'PDB file e.g. molecule.pdb', .true., 'molecule.pdb', &
        &visibility=UI_VIS_STANDARD)
        ! parameter input/output
        call model_validate%add_input(UI_PARM, smpd, &
        &visibility=UI_VIS_STANDARD)
        call model_validate%add_input(UI_PARM, smpd_target, &
        &visibility=UI_VIS_STANDARD)
        ! <no additional inputs>
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! mask controls
        ! computer controls
        ! add to ui_hash
        call add_ui_program('model_validate', model_validate, prgtab, UI_CATEGORY)
    end subroutine new_model_validate

end module simple_ui_validate
