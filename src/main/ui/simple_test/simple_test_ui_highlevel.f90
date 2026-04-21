!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(ui_program), target :: mini_stream
type(ui_program), target :: simulated_workflow
type(ui_program), target :: simulate_particles
type(ui_program), target :: reproject
type(ui_program), target :: subproject_distr
type(ui_program), target :: movie_ppca_subproject_distr

contains

    subroutine construct_test_highlevel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mini_stream(tsttab)
        call new_simulate_particles(tsttab)
        call new_simulated_workflow(tsttab)
        call new_reproject(tsttab)
        call new_subproject_distr(tsttab)
        call new_movie_ppca_subproject_distr(tsttab)
    end subroutine construct_test_highlevel_programs

    subroutine print_test_highlevel_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('HIGH-LEVEL:', C_UNDERLINED)
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') simulate_particles%name%to_char()
        write(logfhandle,'(A)') simulated_workflow%name%to_char()
        write(logfhandle,'(A)') reproject%name%to_char()
        write(logfhandle,'(A)') subproject_distr%name%to_char()
        write(logfhandle,'(A)') movie_ppca_subproject_distr%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_highlevel_programs

    subroutine new_mini_stream( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mini_stream%new(&
        &'mini_stream',&                       ! name
        &'mini_stream',&                       ! descr_short
        &'is a test program for mini_stream',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call mini_stream%add_input(UI_IO, )
        ! parameter input/output
        !call mini_stream%add_input(UI_IMG, )
        ! alternative inputs
        !call mini_stream%add_input(UI_PARM, )
        ! search controls
        !call mini_stream%add_input(UI_SRCH, )
        ! filter controls
        !call mini_stream%add_input(UI_FILT, )
        ! mask controls
        !call mini_stream%add_input(UI_MASK, )
        ! computer controls
        !call mini_stream%add_input(UI_COMP, )
        ! add to ui_hash
        call add_ui_program('mini_stream', mini_stream, tsttab)
    end subroutine new_mini_stream

    subroutine new_reproject( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call reproject%new(&
        &'reproject',&                              ! name
        &'reproject',&                              ! descr_short
        &'is a test program for reproject',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('reproject', reproject, tsttab)
    end subroutine new_reproject

    subroutine new_simulate_particles( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulate_particles%new(&
        &'simulate_particles',&                     ! name
        &'simulate_particles',&                     ! descr_short
        &'is a test program for simulate_particles',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('simulate_particles', simulate_particles, tsttab)
    end subroutine new_simulate_particles

    subroutine new_simulated_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulated_workflow%new(&
        &'simulated_workflow',&                          ! name
        &'self-contained simulated workflow test',&      ! descr_short
        &'is a self-contained simulated workflow test',&
        &'simple_test_exec',&                            ! executable
        &.false.)                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call simulated_workflow%add_input(UI_IO, )
        ! parameter input/output
        !call simulated_workflow%add_input(UI_IMG, )
        ! alternative inputs
        !call simulated_workflow%add_input(UI_PARM, )
        ! search controls
        !call simulated_workflow%add_input(UI_SRCH, )
        ! filter controls
        !call simulated_workflow%add_input(UI_FILT, )
        ! mask controls
        !call simulated_workflow%add_input(UI_MASK, )
        ! computer controls
        !call simulated_workflow%add_input(UI_COMP, nthr)
        ! add to ui_hash
        call add_ui_program('simulated_workflow', simulated_workflow, tsttab)
    end subroutine new_simulated_workflow

    subroutine new_subproject_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call subproject_distr%new(&
        &'subproject_distr',&                                      ! name
        &'test subproject split, parallel exec & merge',&          ! descr_short
        &'Integration test: split project into subprojects, run in parallel via generate_scripts_subprojects, merge back',&
        &'simple_test_exec',&                                      ! executable
        &.true.)                                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('subproject_distr', subproject_distr, tsttab)
    end subroutine new_subproject_distr

    subroutine new_movie_ppca_subproject_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call movie_ppca_subproject_distr%new(&
        &'movie_ppca_subproject_distr',&                                  ! name
        &'split movie chunks + ppca denoise in parallel',&                ! descr_short
        &'Integration test: split filetab movies into equal chunks, build chunk stacks, denoise each chunk with '&
        &'ppca_denoise in parallel subprojects, merge outputs',&
        &'simple_test_exec',&                                             ! executable
        &.false.)   
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call movie_ppca_subproject_distr%add_input(UI_IMG, 'filetab', 'file', 'List of individual movie frame files', 'List of frame files (*.mrcs) to import', 'e.g. movie_frames.txt', .true., '')
        call movie_ppca_subproject_distr%add_input(UI_PARM, 'smpd', 'real', 'SMPD', 'SMPD parameter', 'e.g. 1.3', .true., '')
        call movie_ppca_subproject_distr%add_input(UI_COMP, 'nparts', 'integer', 'Number of parts', 'Number of parts to split the movie frames into', 'e.g. 4', .true., '')
        call add_ui_program('movie_ppca_subproject_distr', movie_ppca_subproject_distr, tsttab)
    end subroutine new_movie_ppca_subproject_distr

end module simple_test_ui_highlevel
