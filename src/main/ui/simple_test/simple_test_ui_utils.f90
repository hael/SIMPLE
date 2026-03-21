!@descr: module defining the user interfaces for utils  programs in the simple_test_exec suite
module simple_test_ui_utils
use simple_ui_modules
implicit none

type(ui_program), target :: ansi_colors
type(ui_program), target :: binoris_test
type(ui_program), target :: binoris_io_test
type(ui_program), target :: cif2mrc
type(ui_program), target :: cif2pdb
type(ui_program), target :: cmdline
type(ui_program), target :: install
type(ui_program), target :: nice
type(ui_program), target :: pdb2mrc
type(ui_program), target :: serialize
type(ui_program), target :: stringmatch

contains

    subroutine construct_test_utils_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_ansi_colors(tsttab)
        call new_binoris_test(tsttab)
        call new_binoris_io_test(tsttab)
        call new_cif2mrc(tsttab)
        call new_cif2pdb(tsttab)
        call new_cmdline(tsttab)
        call new_install(tsttab)
        call new_nice(tsttab)
        call new_pdb2mrc(tsttab)
        call new_serialize(tsttab)
        call new_stringmatch(tsttab)
    end subroutine construct_test_utils_programs

    subroutine print_test_utils_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('UTILS:', C_UNDERLINED)
        write(logfhandle,'(A)') ansi_colors%name%to_char()
        write(logfhandle,'(A)') binoris_test%name%to_char()
        write(logfhandle,'(A)') binoris_io_test%name%to_char()
        write(logfhandle,'(A)') cif2mrc%name%to_char()
        write(logfhandle,'(A)') cif2pdb%name%to_char()
        write(logfhandle,'(A)') cmdline%name%to_char()
        write(logfhandle,'(A)') install%name%to_char()
        write(logfhandle,'(A)') nice%name%to_char()
        write(logfhandle,'(A)') pdb2mrc%name%to_char()
        write(logfhandle,'(A)') serialize%name%to_char()
        write(logfhandle,'(A)') stringmatch%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_utils_programs

    subroutine new_ansi_colors( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call ansi_colors%new(&
        &'ansi_colors',&                       ! name
        &'ansi_colors ',&                      ! descr_short
        &'is a test program for ansi colors',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call ansi_colors%add_input(UI_IO, )
        ! parameter input/output
        !call ansi_colors%add_input(UI_IMG, )
        ! alternative inputs
        !call ansi_colors%add_input(UI_PARM, )
        ! search controls
        !call ansi_colors%add_input(UI_SRCH, )
        ! filter controls
        !call ansi_colors%add_input(UI_FILT, )
        ! mask controls
        !call ansi_colors%add_input(UI_MASK, )
        ! computer controls
        !call ansi_colors%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('ansi_colors', ansi_colors, tsttab)
    end subroutine new_ansi_colors

    subroutine new_binoris_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call binoris_test%new(&
        &'binoris_test',&                      ! name
        &'binoris_test ',&                     ! descr_short
        &'is a test program for binoris',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call binoris_test%add_input(UI_IO, )
        ! parameter input/output
        !call binoris_test%add_input(UI_IMG, )
        ! alternative inputs
        !call binoris_test%add_input(UI_PARM, )
        ! search controls
        !call binoris_test%add_input(UI_SRCH, )
        ! filter controls
        !call binoris_test%add_input(UI_FILT, )
        ! mask controls
        !call binoris_test%add_input(UI_MASK, )
        ! computer controls
        !call binoris_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('binoris_test', binoris_test, tsttab)
    end subroutine new_binoris_test

    subroutine new_binoris_io_test( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call binoris_io_test%new(&
        &'binoris_io_test',&                   ! name
        &'binoris_io_test ',&                  ! descr_short
        &'is a test program for binoris input/output',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call binoris_io_test%add_input(UI_IO, )
        ! parameter input/output
        !call binoris_io_test%add_input(UI_IMG, )
        ! alternative inputs
        !call binoris_io_test%add_input(UI_PARM, )
        ! search controls
        !call binoris_io_test%add_input(UI_SRCH, )
        ! filter controls
        !call binoris_io_test%add_input(UI_FILT, )
        ! mask controls
        !call binoris_io_test%add_input(UI_MASK, )
        ! computer controls
        !call binoris_io_test%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('binoris_io_test', binoris_io_test, tsttab)
    end subroutine new_binoris_io_test

    subroutine new_cif2mrc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call cif2mrc%new(&
        &'cif2mrc',&                         ! name
        &'cif2mrc',&                         ! descr_short
        &'is a test program for PDBx/mmCIF to MRC',&
        &'simple_test_exec',&                ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call cmdline%add_input(UI_IO, )
        ! parameter input/output
        !call cmdline%add_input(UI_IMG, )
        ! alternative inputs
        !call cmdline%add_input(UI_PARM, )
        ! search controls
        !call cmdline%add_input(UI_SRCH, )
        ! filter controls
        !call cmdline%add_input(UI_FILT, )
        ! mask controls
        !call cmdline%add_input(UI_MASK, )
        ! computer controls
        !call cmdline%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('cif2mrc', cif2mrc, tsttab)
    end subroutine new_cif2mrc

    subroutine new_cif2pdb( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call cif2pdb%new(&
        &'cif2pdb',&                         ! name
        &'test for cif2pdb',&                        ! descr_short
        &'is a test program for PDBx/mmCIF to PDB convertion',&
        &'simple_test_exec',&                ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call cmdline%add_input(UI_IO, )
        ! parameter input/output
        !call cmdline%add_input(UI_IMG, )
        ! alternative inputs
        !call cmdline%add_input(UI_PARM, )
        ! search controls
        !call cmdline%add_input(UI_SRCH, )
        ! filter controls
        !call cmdline%add_input(UI_FILT, )
        ! mask controls
        !call cmdline%add_input(UI_MASK, )
        ! computer controls
        !call cmdline%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('cif2pdb', cif2pdb, tsttab)
    end subroutine new_cif2pdb

    subroutine new_cmdline( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call cmdline%new(&
        &'cmdline',&                         ! name
        &'cmdline ',&                        ! descr_short
        &'is a test program for cmdline',&
        &'simple_test_exec',&                ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call cmdline%add_input(UI_IO, )
        ! parameter input/output
        !call cmdline%add_input(UI_IMG, )
        ! alternative inputs
        !call cmdline%add_input(UI_PARM, )
        ! search controls
        !call cmdline%add_input(UI_SRCH, )
        ! filter controls
        !call cmdline%add_input(UI_FILT, )
        ! mask controls
        !call cmdline%add_input(UI_MASK, )
        ! computer controls
        !call cmdline%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('cmdline', cmdline, tsttab)
    end subroutine new_cmdline

    subroutine new_install( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call install%new(&
        &'install',&                         ! name
        &'install ',&                        ! descr_short
        &'is a test program for install',&
        &'simple_test_exec',&                ! executable
        &.false.)                            ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call install%add_input(UI_IO, )
        ! parameter input/output
        !call install%add_input(UI_IMG, )
        ! alternative inputs
        !call install%add_input(UI_PARM, )
        ! search controls
        !call install%add_input(UI_SRCH, )
        ! filter controls
        !call install%add_input(UI_FILT, )
        ! mask controls
        !call install%add_input(UI_MASK, )
        ! computer controls
        !call install%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('install', install, tsttab)
    end subroutine new_install

    subroutine new_nice( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call nice%new(&
        &'nice',&                         ! name
        &'nice ',&                        ! descr_short
        &'is a test program for NICE',&
        &'simple_test_exec',&             ! executable
        &.false.)                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call nice%add_input(UI_IO, )
        ! parameter input/output
        !call nice%add_input(UI_IMG, )
        ! alternative inputs
        !call nice%add_input(UI_PARM, )
        ! search controls
        !call nice%add_input(UI_SRCH, )
        ! filter controls
        !call nice%add_input(UI_FILT, )
        ! mask controls
        !call nice%add_input(UI_MASK, )
        ! computer controls
        !call nice%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('nice', nice, tsttab)
    end subroutine new_nice

    subroutine new_pdb2mrc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call pdb2mrc%new(&
        &'pdb2mrc',&                         ! name
        &'pdb2mrc ',&                     ! descr_short
        &'is a test program for pdb2mrc',&
        &'simple_test_exec',&             ! executable
        &.false.)                         ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call pdb2mrc%add_input(UI_IO, )
        ! parameter input/output
        !call pdb2mrc%add_input(UI_IMG, )
        ! alternative inputs
        !call pdb2mrc%add_input(UI_PARM, )
        ! search controls
        !call pdb2mrc%add_input(UI_SRCH, )
        ! filter controls
        !call pdb2mrc%add_input(UI_FILT, )
        ! mask controls
        !call pdb2mrc%add_input(UI_MASK, )
        ! computer controls
        !call pdb2mrc%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('pdb2mrc', pdb2mrc, tsttab)
    end subroutine new_pdb2mrc

    subroutine new_serialize( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call serialize%new(&
        &'serialize',&                         ! name
        &'serialize ',&                        ! descr_short
        &'is a test program for serialize',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call serialize%add_input(UI_IO, )
        ! parameter input/output
        !call serialize%add_input(UI_IMG, )
        ! alternative inputs
        !call serialize%add_input(UI_PARM, )
        ! search controls
        !call serialize%add_input(UI_SRCH, )
        ! filter controls
        !call serialize%add_input(UI_FILT, )
        ! mask controls
        !call serialize%add_input(UI_MASK, )
        ! computer controls
        !call serialize%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('serialize', serialize, tsttab)
    end subroutine new_serialize

    subroutine new_stringmatch( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call stringmatch%new(&
        &'stringmatch',&                       ! name
        &'stringmatch ',&                      ! descr_short
        &'is a test program for stringmatch',&
        &'simple_test_exec',&                  ! executable
        &.false.)                              ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call stringmatch%add_input(UI_IO, )
        ! parameter input/output
        !call stringmatch%add_input(UI_IMG, )
        ! alternative inputs
        !call stringmatch%add_input(UI_PARM, )
        ! search controls
        !call stringmatch%add_input(UI_SRCH, )
        ! filter controls
        !call stringmatch%add_input(UI_FILT, )
        ! mask controls
        !call stringmatch%add_input(UI_MASK, )
        ! computer controls
        !call stringmatch%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('stringmatch', stringmatch, tsttab)
    end subroutine new_stringmatch

end module simple_test_ui_utils
