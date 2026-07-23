!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(ui_program), target :: mini_stream
type(ui_program), target :: simulated_workflow
type(ui_program), target :: simulate_particles
type(ui_program), target :: reproject
type(ui_program), target :: subproject_distr
type(ui_program), target :: ptcls_ppca_subproject_distr
type(ui_program), target :: flex_preimage_identity
type(ui_program), target :: flex_preimage_basis_ab
type(ui_program), target :: pcg_recon_ctf_free

contains

    subroutine construct_test_highlevel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mini_stream(tsttab)
        call new_simulate_particles(tsttab)
        call new_simulated_workflow(tsttab)
        call new_reproject(tsttab)
        call new_subproject_distr(tsttab)
        call new_ptcls_ppca_subproject_distr(tsttab)
        call new_flex_preimage_identity(tsttab)
        call new_flex_preimage_basis_ab(tsttab)
        call new_pcg_recon_ctf_free(tsttab)
    end subroutine construct_test_highlevel_programs

    subroutine print_test_highlevel_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('HIGH-LEVEL:', C_UNDERLINED)
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') simulate_particles%name%to_char()
        write(logfhandle,'(A)') simulated_workflow%name%to_char()
        write(logfhandle,'(A)') reproject%name%to_char()
        write(logfhandle,'(A)') subproject_distr%name%to_char()
        write(logfhandle,'(A)') ptcls_ppca_subproject_distr%name%to_char()
        write(logfhandle,'(A)') flex_preimage_identity%name%to_char()
        write(logfhandle,'(A)') flex_preimage_basis_ab%name%to_char()
        write(logfhandle,'(A)') pcg_recon_ctf_free%name%to_char()
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
        call simulated_workflow%add_input(UI_PARM, 'system', 'multi', 'Embedded molecular system', &
            &'Embedded coordinates used to generate the simulated data(6vxx|1jxy)', &
            &'(6vxx|1jxy)', .true., '')
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

    subroutine new_ptcls_ppca_subproject_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call ptcls_ppca_subproject_distr%new(&
        &'ptcls_ppca_subproject_distr',&                                  ! name
        &'split particle chunks + ppca denoise in parallel',&                ! descr_short
        &'Integration test: split filetab particles into equal chunks, build chunk stacks, denoise each chunk with '&
        &'ppca_denoise in parallel subprojects, merge outputs',&
        &'simple_test_exec',&                                             ! executable
        &.false.)   
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ptcls_ppca_subproject_distr%add_input(UI_IMG, 'filetab', 'file', 'List of individual particle files', 'List of particle files (*.mrcs) to import', 'e.g. particle_frames.txt', .true., '')
        call ptcls_ppca_subproject_distr%add_input(UI_PARM, 'smpd', 'real', 'SMPD', 'SMPD parameter', 'e.g. 1.3', .true., '')
        call ptcls_ppca_subproject_distr%add_input(UI_COMP, 'nparts', 'integer', 'Number of parts', 'Number of parts to split the particle files into', 'e.g. 4', .true., '')
        call add_ui_program('ptcls_ppca_subproject_distr', ptcls_ppca_subproject_distr, tsttab)
    end subroutine new_ptcls_ppca_subproject_distr

    subroutine new_flex_preimage_identity( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call flex_preimage_identity%new(&
        &'flex_preimage_identity',&
        &'flex pre-image reconstruction identity diagnostic',&
        &'Compares flex observation/operator planes with refine3D planes, then checks that constant z=1 reproduces reconstruct3D from the registered project',&
        &'simple_test_exec',&
        &.true.)
        call flex_preimage_identity%add_input(UI_IMG, 'vol1', 'file', &
            'Mean volume', 'Fixed mean volume used for the residual identity reconstruction', &
            'input volume e.g. rec_final_state01.mrc', .true., '')
        call flex_preimage_identity%add_input(UI_PARM, projfile, required_override=.true.)
        call flex_preimage_identity%add_input(UI_PARM, nspace, required_override=.true.)
        call flex_preimage_identity%add_input(UI_COMP, nthr, required_override=.false.)
        call add_ui_program('flex_preimage_identity', flex_preimage_identity, tsttab)
    end subroutine new_flex_preimage_identity

    subroutine new_flex_preimage_basis_ab( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call flex_preimage_basis_ab%new(&
        &'flex_preimage_basis_ab',&
        &'flex pre-image raw/canonical latent-basis A/B diagnostic',&
        &'Reconstructs the same Nyström representatives from raw graph eigenfunctions and an inverse-target whitened basis',&
        &'simple_test_exec',&
        &.true.)
        call flex_preimage_basis_ab%add_input(UI_IMG, 'vol1', 'file', &
            'Mean volume', 'Fixed mean volume used for both residual reconstructions', &
            'input volume e.g. rec_final_state01.mrc', .true., '')
        call flex_preimage_basis_ab%add_input(UI_PARM, projfile, required_override=.true.)
        call flex_preimage_basis_ab%add_input(UI_PARM, nspace, required_override=.true.)
        call flex_preimage_basis_ab%add_input(UI_FILT, 'neigs', 'num', 'Diffusion eigenpair scan limit', &
            'Positive number of graph eigenpairs supplied to the A/B test', '# eigenpairs', .true., 12.0)
        call flex_preimage_basis_ab%add_input(UI_FILT, 'npreimages', 'num', 'Representative states', &
            'Number of k-medoids reconstructed by both A/B arms', '# states', .false., 8.0)
        call flex_preimage_basis_ab%add_input(UI_FILT, 'icm', 'binary', 'ICM rank selection', &
            'Set no to test every eigenpair returned by neigs', 'yes/no', .false., 'yes')
        call flex_preimage_basis_ab%add_input(UI_COMP, nthr, required_override=.false.)
        call add_ui_program('flex_preimage_basis_ab', flex_preimage_basis_ab, tsttab)
    end subroutine new_flex_preimage_basis_ab

    subroutine new_pcg_recon_ctf_free( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call pcg_recon_ctf_free%new(&
        &'pcg_recon_ctf_free',&                      ! name
        &'CTF-free PCG reconstruction operator validation',&
        &'Milestone 0 of the CTF/sigma-weighted PCG reconstruction design note: in-memory '&
        &'adjoint dot-product test, normal-operator test, and no-CTF/no-noise synthetic '&
        &'recovery test against a deterministic phantom. Self-contained, no project required.',&
        &'simple_test_exec',&                       ! executable
        &.false.)                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('pcg_recon_ctf_free', pcg_recon_ctf_free, tsttab)
    end subroutine new_pcg_recon_ctf_free

end module simple_test_ui_highlevel
