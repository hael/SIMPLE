!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('highlevel', 'High-Level', 40)
type(ui_program), target :: mini_stream
type(ui_program), target :: simulated_workflow
type(ui_program), target :: simulate_particles
type(ui_program), target :: reproject
type(ui_program), target :: subproject_distr
type(ui_program), target :: ptcls_ppca_subproject_distr
type(ui_program), target :: pcg_recon
type(ui_program), target :: pcg_frac_update
type(ui_program), target :: pcg_priors
type(ui_program), target :: rec3D_backends

contains

    subroutine construct_test_highlevel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mini_stream(tsttab)
        call new_simulate_particles(tsttab)
        call new_simulated_workflow(tsttab)
        call new_reproject(tsttab)
        call new_subproject_distr(tsttab)
        call new_ptcls_ppca_subproject_distr(tsttab)
        call new_pcg_recon(tsttab)
        call new_pcg_frac_update(tsttab)
        call new_pcg_priors(tsttab)
        call new_rec3D_backends(tsttab)
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
        write(logfhandle,'(A)') pcg_recon%name%to_char()
        write(logfhandle,'(A)') pcg_frac_update%name%to_char()
        write(logfhandle,'(A)') pcg_priors%name%to_char()
        write(logfhandle,'(A)') rec3D_backends%name%to_char()
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
        ! <no additional inputs>
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
        call add_ui_program('mini_stream', mini_stream, tsttab, UI_CATEGORY)
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
        call add_ui_program('reproject', reproject, tsttab, UI_CATEGORY)
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
        call add_ui_program('simulate_particles', simulate_particles, tsttab, UI_CATEGORY)
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
            &'Embedded coordinates used to generate the simulated data(6vxx|1jxy)','', .true., '', &
        &choices=ui_choices([character(len=4) :: '6vxx', '1jxy']))
        call simulated_workflow%add_input(UI_PARM, 'picker', 'multi', 'Picker under test', &
            &'Particle picker used by the simulated workflow(segdiam|new){segdiam}','', .false., 'segdiam', &
        &choices=ui_choices([character(len=7) :: 'segdiam', 'new']))
        call simulated_workflow%add_input(UI_PARM, 'refpick_backend', 'multi', 'Reference-picker backend', &
            &'Reference-picker implementation used when picker=new(legacy|optimized|compare){legacy}','', .false., 'legacy', &
        &choices=ui_choices([character(len=9) :: 'legacy', 'optimized', 'compare']))
        ! <no additional inputs>
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
        call add_ui_program('simulated_workflow', simulated_workflow, tsttab, UI_CATEGORY)
    end subroutine new_simulated_workflow

    subroutine new_subproject_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call subproject_distr%new(&
        &'subproject_distr',&                                      ! name
        &'test subproject split, parallel exec & merge',&          ! descr_short
        &'Integration test: split project into subprojects, run in parallel via generate_scripts_subprojects, merge back',&
        &'simple_test_exec',&                                      ! executable
        &.false.)                                                  ! requires sp_project
        ! add to ui_hash
        call add_ui_program('subproject_distr', subproject_distr, tsttab, UI_CATEGORY)
    end subroutine new_subproject_distr

    subroutine new_ptcls_ppca_subproject_distr( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call ptcls_ppca_subproject_distr%new(&
        &'ptcls_ppca_subproject_distr',&                                  ! name
        &'split particle chunks + ppca denoise in parallel',&                ! descr_short
        &'Integration test: split filetab particles into equal chunks, build chunk stacks, denoise each chunk with '//&
        &'ppca_denoise in parallel subprojects, merge outputs',&
        &'simple_test_exec',&                                             ! executable
        &.false.)   
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        call ptcls_ppca_subproject_distr%add_input(UI_FILE, 'filetab', 'file', 'List of individual particle files', 'List of particle files (*.mrcs) to import', 'e.g. particle_frames.txt', .true., '')
        call ptcls_ppca_subproject_distr%add_input(UI_PARM, 'smpd', 'real', 'SMPD', 'SMPD parameter', 'e.g. 1.3', .true., '')
        call ptcls_ppca_subproject_distr%add_input(UI_COMP, 'nparts', 'integer', 'Number of parts', 'Number of parts to split the particle files into', 'e.g. 4', .true., '')
        call add_ui_program('ptcls_ppca_subproject_distr', ptcls_ppca_subproject_distr, tsttab, UI_CATEGORY)
    end subroutine new_ptcls_ppca_subproject_distr

    subroutine new_pcg_recon( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call pcg_recon%new(&
        &'pcg_recon',&                               ! name
        &'PCG reconstruction operator and solver validation',&
        &'Single in-memory gate for the PCG reconstruction path, against a deterministic '//&
        &'asymmetric phantom and no project. Nine fail-fast stages in increasing cost: adjoint '//&
        &'dot-product identity with a trivial transfer and then with a real astigmatic CTF, shift '//&
        &'and sigma2; normal-operator symmetry and positive-definiteness; synthetic recovery; '//&
        &'kernelized-vs-matrix-free equivalence reported for all voxels and for an interior '//&
        &'region; kernel invariants (shift-invariant, CTF-dependent) and the sampling-density '//&
        &'preconditioner; streaming batch accumulation reproducing the monolithic solve; and '//&
        &'deapodization against ENVELOPE-FREE observations, the one stage that does not commit '//&
        &'an inverse crime; and C2 coordinate replication equivalence.',&
        &'simple_test_exec',&                        ! executable
        &.false.)                                    ! requires sp_project
        ! add to ui_hash
        call add_ui_program('pcg_recon', pcg_recon, tsttab, UI_CATEGORY)
    end subroutine new_pcg_recon

    subroutine new_pcg_frac_update( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call pcg_frac_update%new(&
        &'pcg_frac_update',&
        &'Project-backed PCG fractional-update validation',&
        &'Uses reconstruct3D inputs and internally partitions every populated state/half into two '//&
        &'deterministic complementary subsets. Validates raw (B,D) additivity, u/f continuation '//&
        &'weighting, full-mass ensemble preservation, and continuation artifact replay.',&
        &'simple_test_exec',&
        &.true.)
        call pcg_frac_update%add_input(UI_PARM, 'box_crop', 'num', 'Reconstruction box', &
        &'Even Fourier-cropped reconstruction box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0)
        call pcg_frac_update%add_input(UI_PARM, 'projrec', 'binary', 'Projection-direction reconstruction', &
        &'Accepted for reconstruct3D command compatibility; PCG validation requires no', '', .false., 'no', &
        &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call pcg_frac_update%add_input(UI_SRCH, trs)
        call pcg_frac_update%add_input(UI_SRCH, pgrp)
        call pcg_frac_update%add_input(UI_SRCH, objfun)
        call pcg_frac_update%add_input(UI_SRCH, ptcl_src)
        call pcg_frac_update%add_input(UI_FILT, ml_reg)
        call pcg_frac_update%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map', &
        &'Accepted for reconstruct3D command compatibility; no production maps are written by this test', &
        &'', .false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call pcg_frac_update%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Iterations used for the continuation replay comparison', 'iterations{2}', .false., 2.)
        call pcg_frac_update%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Tolerance used for the continuation replay comparison', 'tolerance{0}', .false., 0.0)
        call pcg_frac_update%add_input(UI_MASK, mskdiam)
        call pcg_frac_update%add_input(UI_COMP, nparts, required_override=.false.)
        call pcg_frac_update%add_input(UI_COMP, nthr)
        call add_ui_program('pcg_frac_update', pcg_frac_update, tsttab, UI_CATEGORY)
    end subroutine new_pcg_frac_update

    subroutine new_pcg_priors( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call pcg_priors%new(&
        &'pcg_priors',&
        &'PCG prior-operator unit gate',&
        &'In-memory, project-free unit gate for the PCG prior operators (pcg_priors.md). Currently '//&
        &'validates the graded solvent-flatness precision Q_s: normalization contract (unit mean '//&
        &'diagonal on the effective support), adjoint identity, positive semidefiniteness, constant '//&
        &'null space, zero action at zero solvent confidence, graded-edge continuity, the '//&
        &'finite-difference gradient of the penalty, and composition with the masked normal operator.',&
        &'simple_test_exec',&
        &.false.)
        call add_ui_program('pcg_priors', pcg_priors, tsttab, UI_CATEGORY)
    end subroutine new_pcg_priors

    subroutine new_rec3D_backends( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call rec3D_backends%new(&
        &'rec3D_backends',&
        &'Same-inputs gridding vs PCG reconstruction comparison',&
        &'Reconstructs one fixed set of particles/orientations/sigma2 with reconstruct3D using the gridding '//&
        &'and the PCG backend, in a numbered execution directory named after the run settings (run it inside '//&
        &'a refine3D run directory; the sigma2 group files and any NU evidence envelope are symlinked in; '//&
        &'mkdir=no runs in the current directory instead), and compares the merged maps: per-shell amplitude '//&
        &'ratio and FSC between '//&
        &'backends, and the radial real-space profile ratio that exposes a deapodization mismatch. Writes '//&
        &'recvol_stateXX_gridding.mrc and recvol_stateXX_pcg.mrc. Measurement only, no thresholds '//&
        &'(doc/implementation_notes/drop_legacy_box_division.md, plan step 2).',&
        &'simple_test_exec',&
        &.true.)
        call rec3D_backends%add_input(UI_PARM, 'box_crop', 'num', 'Reconstruction box', &
        &'Even Fourier-cropped reconstruction box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0)
        call rec3D_backends%add_input(UI_SRCH, trs)
        call rec3D_backends%add_input(UI_SRCH, pgrp)
        call rec3D_backends%add_input(UI_SRCH, objfun)
        call rec3D_backends%add_input(UI_FILT, ml_reg)
        call rec3D_backends%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'PCG iterations for the comparison', 'iterations{2}', .false., 2.)
        call rec3D_backends%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Tolerance for the comparison', 'tolerance{0}', .false., 0.0)
        call rec3D_backends%add_input(UI_FILT, 'pcg_solvent_lambda_rel', 'num', 'PCG solvent-flatness prior strength', &
        &'Solvent-flatness prior strength relative to the PCG data scale, applied on the pcg leg; requires '//&
        &'ml_reg (objfun=euclid) and the NU evidence envelope; 0 = off', 'strength{0}', .false., 0.0)
        call rec3D_backends%add_input(UI_IMG, 'vol1', 'file', 'Ground-truth volume', &
        &'Known volume the particles were simulated from; enables the radial recon/truth table', &
        &'e.g. truth.mrc', .false., '')
        call rec3D_backends%add_input(UI_FILT, 'lp', 'num', 'Low-pass limit for the truth comparison', &
        &'Applied identically to the truth and both reconstructions before the radial comparison', &
        &'Angstroms{0}', .false., 0.0)
        call rec3D_backends%add_input(UI_MASK, mskdiam)
        call rec3D_backends%add_input(UI_COMP, nthr)
        call add_ui_program('rec3D_backends', rec3D_backends, tsttab, UI_CATEGORY)
    end subroutine new_rec3D_backends

end module simple_test_ui_highlevel
