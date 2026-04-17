!@descr: module defining the user interfaces for miscellaneous programs in the simple_exec suite
module simple_ui_other
use simple_ui_modules
implicit none

type(ui_program), target :: cif2pdb
type(ui_program), target :: fractionate_movies
type(ui_program), target :: split_
type(ui_program), target :: split_stack

contains

    subroutine construct_other_programs( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call new_cif2pdb(prgtab)
        call new_fractionate_movies(prgtab)
        call new_split_(prgtab)
        call new_split_stack(prgtab)
    end subroutine construct_other_programs

    subroutine print_other_programs( logfhandle )
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('OTHER UTILITIES:', C_UNDERLINED)
        write(logfhandle,'(A)') cif2pdb%name%to_char()
        write(logfhandle,'(A)') fractionate_movies%name%to_char()
        write(logfhandle,'(A)') split_%name%to_char()
        write(logfhandle,'(A)') split_stack%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_other_programs

    subroutine new_cif2pdb( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call cif2pdb%new(&
        &'cif2pdb',&                                       ! name
        &'convert PDBx/mmCIF to PDB',&                     ! descr_short
        &'is a program for converting PDBx/mmCIF to PDB',& ! descr_long
        &'simple_exec',&                                   ! executable
        &.false.)                                          ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call cif2pdb%add_input(UI_IMG, 'ciffile', 'file', 'PDBx/mmCIF input coordinates file', 'Input coordinates file in PDBx/mmCIF format', 'PDBx/mmCIF file e.g. molecule.cif', .true., 'molecule.cif')
        ! parameter input/output
        ! computer controls
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>               
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('cif2pdb', cif2pdb, prgtab)
    end subroutine new_cif2pdb

    subroutine new_fractionate_movies( prgtab )
        class(ui_hash), intent(inout) :: prgtab
        call fractionate_movies%new(&
        &'fractionate_movies', &
        &'Re-generate micrographs from selected movie frames',&
        &'is a distributed program for re-generating micrographs from a subset of movie frames',&
        &'simple_exec',&
        &.true.)
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        ! <empty>
        ! parameter input/output
        call fractionate_movies%add_input(UI_PARM, 'fromf', 'num', 'Starting frame', &
        & 'Starting movie frame for micrograph re-generation', 'frame index{1}', .false., 1.0)
        call fractionate_movies%add_input(UI_PARM, 'tof', 'num', 'Final frame', &
        & 'Final movie frame for micrograph re-generation(0=all)', 'frame index{0}', .false., 0.0)
        call fractionate_movies%add_input(UI_PARM, 'mcconvention', 'str', 'Movie alignment convention', &
        & 'Movie alignment and naming convention(simple|unblur|relion|motioncorr|cryosparc|cs){simple}', &
        & '(simple|unblur|relion|motioncorr|cryosparc|cs){simple}', .false., 'simple')
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! computer controls
        call fractionate_movies%add_input(UI_COMP, nparts)
        call fractionate_movies%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('fractionate_movies', fractionate_movies, prgtab)
    end subroutine new_fractionate_movies

    subroutine new_split_( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! alternative inputs
        ! <empty>
        ! search controls
        ! <empty>               
        ! filter controls
        ! <empty>
        ! mask controls
        ! <empty>
        ! add to ui_hash
        call add_ui_program('split', split_, prgtab)
    end subroutine new_split_

    subroutine new_split_stack( prgtab )
        class(ui_hash), intent(inout) :: prgtab
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
        ! add to ui_hash
        call add_ui_program('split_stack', split_stack, prgtab)
    end subroutine new_split_stack

end module simple_ui_other
