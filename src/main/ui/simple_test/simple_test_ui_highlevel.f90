!@descr: module defining the user interfaces for highlevel test programs in the simple_test_exec suite
module simple_test_ui_highlevel
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('highlevel', 'High-level tests', 40)
type(ui_program), target :: mini_stream
type(ui_program), target :: simulated_workflow
type(ui_program), target :: simulate_particles
type(ui_program), target :: solve3D_addon
type(ui_program), target :: solve3D_addon_snapshots
type(ui_program), target :: cont_refine3D_1jxy
type(ui_program), target :: cont_refine3D_bench
type(ui_program), target :: polar_matching_bench
type(ui_program), target :: multistate_assignment
type(ui_program), target :: multistate_workflows
type(ui_program), target :: state_fixture
type(ui_program), target :: pcg_recon
type(ui_program), target :: pcg_frac_update
type(ui_program), target :: rec3D_backends
type(ui_program), target :: single_atoms_stats
type(ui_program), target :: species_discovery
type(ui_program), target :: single_workflow
type(ui_program), target :: stream_preproc
type(ui_program), target :: flex_pca_blobs
type(ui_program), target :: write_state_weights_labels
type(ui_program), target :: write_state_weights_mixed
type(ui_program), target :: nightly
type(ui_program), target :: profile
type(ui_program), target :: identity

contains

    subroutine construct_test_highlevel_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_mini_stream(tsttab)
        call new_simulate_particles(tsttab)
        call new_solve3D_addon(tsttab)
        call new_solve3D_addon_snapshots(tsttab)
        call new_cont_refine3D_1jxy(tsttab)
        call new_cont_refine3D_bench(tsttab)
        call new_polar_matching_bench(tsttab)
        call new_multistate_assignment(tsttab)
        call new_multistate_workflows(tsttab)
        call new_state_fixture(tsttab)
        call new_simulated_workflow(tsttab)
        call new_pcg_recon(tsttab)
        call new_pcg_frac_update(tsttab)
        call new_rec3D_backends(tsttab)
        call new_single_atoms_stats(tsttab)
        call new_species_discovery(tsttab)
        call new_single_workflow(tsttab)
        call new_stream_preproc(tsttab)
        call new_flex_pca_blobs(tsttab)
        call new_write_state_weights(tsttab)
        call new_nightly(tsttab)
        call new_profile(tsttab)
        call new_identity(tsttab)
    end subroutine construct_test_highlevel_programs

    subroutine print_test_highlevel_programs( logfhandle)
        integer, intent(in) :: logfhandle
        write(logfhandle,'(A)') format_str('High-level tests:', C_UNDERLINED)
        write(logfhandle,'(A)') mini_stream%name%to_char()
        write(logfhandle,'(A)') simulate_particles%name%to_char()
        write(logfhandle,'(A)') solve3D_addon%name%to_char()
        write(logfhandle,'(A)') solve3D_addon_snapshots%name%to_char()
        write(logfhandle,'(A)') cont_refine3D_1jxy%name%to_char()
        write(logfhandle,'(A)') cont_refine3D_bench%name%to_char()
        write(logfhandle,'(A)') polar_matching_bench%name%to_char()
        write(logfhandle,'(A)') multistate_assignment%name%to_char()
        write(logfhandle,'(A)') multistate_workflows%name%to_char()
        write(logfhandle,'(A)') state_fixture%name%to_char()
        write(logfhandle,'(A)') simulated_workflow%name%to_char()
        write(logfhandle,'(A)') pcg_recon%name%to_char()
        write(logfhandle,'(A)') pcg_frac_update%name%to_char()
        write(logfhandle,'(A)') rec3D_backends%name%to_char()
        write(logfhandle,'(A)') single_atoms_stats%name%to_char()
        write(logfhandle,'(A)') species_discovery%name%to_char()
        write(logfhandle,'(A)') single_workflow%name%to_char()
        write(logfhandle,'(A)') stream_preproc%name%to_char()
        write(logfhandle,'(A)') flex_pca_blobs%name%to_char()
        write(logfhandle,'(A)') write_state_weights_labels%name%to_char()
        write(logfhandle,'(A)') write_state_weights_mixed%name%to_char()
        write(logfhandle,'(A)') nightly%name%to_char()
        write(logfhandle,'(A)') profile%name%to_char()
        write(logfhandle,'(A)') identity%name%to_char()
        write(logfhandle,'(A)') ''
    end subroutine print_test_highlevel_programs

    subroutine new_mini_stream( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call mini_stream%new(&
        &'mini_stream',&                       ! name
        &'validate mini-stream with the embedded 6VXX and 1JXY systems',&
        &'runs molecular 6VXX and 1JXY suites and validates CTF recovery, picking, extraction, '//&
        &'solve2D class averages, and shape ranking',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Mini Stream Validation')                                   ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        call mini_stream%add_input(UI_PARM, 'suite', 'str', 'Run one molecular suite', &
            &'Choose 6vxx or 1jxy for a focused run; use list to print the accepted names', &
            &'', .false., '')
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

    subroutine new_simulate_particles( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulate_particles%new(&
        &'simulate_particles',&                     ! name
        &'quantitatively validate particle simulation from embedded 6VXX',&
        &'checks volume and stack metadata, every image variance, orientation diversity, shifts, '//&
        &'and CTF parameter bounds against deterministic inputs; a mixture of three states with given '//&
        &'populations (exact counts, each image from its own state''s volume); and the multi-state test '//&
        &'fixture in its families and modes, its cache and the point groups of its truths',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Particle Simulation Validation')                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('simulate_particles', simulate_particles, tsttab, UI_CATEGORY)
    end subroutine new_simulate_particles

    subroutine new_solve3D_addon( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call solve3D_addon%new(&
        &'solve3D_addon',&                          ! name
        &'validate solve3D_addon end to end on simulated particles',&
        &'solve3D on a seeded subset of simulated particles of a symmetry-broken 6VXX map, then '//&
        &'solve3D_addon on all of them, checked against the simulation truth (poses, map) and the base run',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='solve3D Add-On Validation')                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('solve3D_addon', solve3D_addon, tsttab, UI_CATEGORY)
    end subroutine new_solve3D_addon

    subroutine new_solve3D_addon_snapshots( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call solve3D_addon_snapshots%new(&
        &'solve3D_addon_snapshots',&
        &'generate cumulative projects for solve3D_addon stream tests',&
        &'splits all input particle images into append-only physical stacks and writes snapshot1.simple through snapshotN.simple',&
        &'simple_test_exec',&
        &.true., display_name='solve3D Add-On Snapshots')
        call solve3D_addon_snapshots%add_input(UI_PARM, 'nptcls_base', 'num', 'Base particle count', &
            &'Particles in snapshot1 and the parent solve3D run', 'particles', .true., 0.)
        call solve3D_addon_snapshots%add_input(UI_PARM, 'nsnapshots', 'num', 'Snapshot count', &
            &'Number of cumulative projects, including the base project', 'projects', .true., 0.)
        call add_ui_program('solve3D_addon_snapshots', solve3D_addon_snapshots, tsttab, UI_CATEGORY)
    end subroutine new_solve3D_addon_snapshots

    subroutine new_cont_refine3D_1jxy( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call cont_refine3D_1jxy%new(&
        &'cont_refine3D_1jxy',&                     ! name
        &'validate continuous Cartesian pose refinement end to end on simulated 1JYX particles',&
        &'refine3D refine=cont (euclid, cc, distributed) and refine3D_auto pose_cont=only on 1 000 simulated 1JYX '//&
        &'particles from poses perturbed by 15 degrees and 2 pixels, and as the continuation of a polar refine3D_auto, '//&
        &'checked against the simulation truth',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Continuous Refinement Validation')                                   ! requires sp_project
        ! add to ui_hash
        call add_ui_program('cont_refine3D_1jxy', cont_refine3D_1jxy, tsttab, UI_CATEGORY)
    end subroutine new_cont_refine3D_1jxy

    subroutine new_cont_refine3D_bench( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call cont_refine3D_bench%new(&
        &'cont_refine3D_bench',&                    ! name
        &'benchmark continuous Cartesian pose refinement on simulated 1JYX particles',&
        &'One refine3D refine=cont pass per regime and objective on simulated 1JYX particles at box 256, '//&
        &'against the reconstruction at the truth poses: '//&
        &'polish (poses 2 degrees and 1 A from the truth, lp 3 A) and capture (15 degrees and 2.6 A, lp 8 A). '//&
        &'A measurement with no thresholds and no CTest entry: run directories, BENCH/ and the pose errors '//&
        &'(metrics.tsv) stay in ./cont_refine3D_bench '//&
        &'(doc/implementation_notes/completed/cont_lm_acceleration.md, section 6.3)',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Continuous Refinement Benchmark')                                ! requires sp_project
        call cont_refine3D_bench%add_input(UI_PARM, 'suite', 'str', 'Regime', &
            &'Run one regime: polish or capture; list prints them. Default runs both', '', .false., '')
        call cont_refine3D_bench%add_input(UI_PARM, 'objfun', 'str', 'Objective', &
            &'Run one objective: euclid or cc. Default runs both', '', .false., '')
        call cont_refine3D_bench%add_input(UI_PARM, 'nptcls', 'num', 'Particles', &
            &'Number of simulated particles', '# particles{10000}', .false., 10000.)
        call cont_refine3D_bench%add_input(UI_PARM, 'snr', 'num', 'SNR', &
            &'Signal-to-noise ratio of the simulated particles; cryo-EM particles are 0.02 or better', &
            &'signal-to-noise ratio{0.05}', .false., 0.05)
        call cont_refine3D_bench%add_input(UI_COMP, nparts, required_override=.false.)
        call cont_refine3D_bench%add_input(UI_COMP, nthr, required_override=.false.)
        ! add to ui_hash
        call add_ui_program('cont_refine3D_bench', cont_refine3D_bench, tsttab, UI_CATEGORY)
    end subroutine new_cont_refine3D_bench

    subroutine new_polar_matching_bench( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call polar_matching_bench%new(&
        &'polar_matching_bench',&                   ! name
        &'benchmark the layouts of the dense probability table''s euclid scoring',&
        &'Synthetic references and CTF-weighted particles at production polar sizes (box 180 at 2.055 A, '//&
        &'crop 160, mask 270 A, lp 6 A), scored as the dense probability table scores them, per layout, '//&
        &'case (1000 and 2500 directions; 1, 2 and 4 states) and team (1, 5 and nthr threads). '//&
        &'A measurement with no thresholds and no CTest entry: rows are appended to ./polar_matching_bench.tsv '//&
        &'(./polar_matching_check.tsv for suite=check); run with OPENBLAS_NUM_THREADS=1 '//&
        &'(doc/implementation_notes/planned/polar_matching_layout_bench.md)',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Polar Matching Benchmark')                                       ! requires sp_project
        call polar_matching_bench%add_input(UI_PARM, 'suite', 'str', 'Layout', &
            &'Run one layout (per_pair, replica, blocked, states_inner, matrix, matrix_t, per_pair_cc), the identity '//&
            &'checks of the layouts (check) or of the production kernels (identity, 32 particles unless nptcls is '//&
            &'given, one-thread results written to ./polar_kernel_identity_*.bin); list prints them. '//&
            &'Default times all euclid layouts', '', .false., '')
        call polar_matching_bench%add_input(UI_PARM, 'nspace', 'num', 'Directions', &
            &'Run one number of projection directions. Default runs 1000 and 2500', '', .false., 0.)
        call polar_matching_bench%add_input(UI_PARM, 'nstates', 'num', 'States', &
            &'Run one number of states. Default runs 1, 2 and 4', '', .false., 0.)
        call polar_matching_bench%add_input(UI_PARM, 'nptcls', 'num', 'Particles', &
            &'Particles per measurement (at most 200 at one thread)', '# particles{2000}', .false., 2000.)
        call polar_matching_bench%add_input(UI_COMP, nthr, required_override=.false.)
        ! add to ui_hash
        call add_ui_program('polar_matching_bench', polar_matching_bench, tsttab, UI_CATEGORY)
    end subroutine new_polar_matching_bench

    subroutine new_multistate_assignment( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call multistate_assignment%new(&
        &'multistate_assignment',&                  ! name
        &'validate multi-state refine3D against the known states of a simulated blob phantom',&
        &'refine3D nstates=3 with refine=prob and with refine=prob_neigh on 1 500 simulated particles of three '//&
        &'blob phantoms (populations 50, 30 and 20 %), from the truths low-passed to 20 A and random initial '//&
        &'states, the state labels, maps and poses checked against the simulation truth, with shuffled-label and '//&
        &'duplicated-maps controls in every refinement and an identical-states control on request',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Multi-State Assignment Validation')                            ! requires sp_project
        call multistate_assignment%add_input(UI_PARM, 'suite', 'str', 'Sub-suite', &
            &'Run one case: prob, prob_neigh or identical (the identical-states control); list prints them. '//&
            &'Default runs prob and prob_neigh', '', .false., '')
        call multistate_assignment%add_input(UI_COMP, nthr, required_override=.false.)
        ! add to ui_hash
        call add_ui_program('multistate_assignment', multistate_assignment, tsttab, UI_CATEGORY)
    end subroutine new_multistate_assignment

    subroutine new_multistate_workflows( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call multistate_workflows%new(&
        &'multistate_workflows',&                   ! name
        &'validate multi-state solve3D and refine3D_states against known beta-galactosidase states',&
        &'On 3 000 simulated particles of the tetramer, the tetramer less a chain and a dimer of '//&
        &'beta-galactosidase (populations 50, 30 and 20 %): solve3D_cavgs then solve3D cavg_ini_ext=yes on the '//&
        &'class averages of an oracle 2D classification, solve3D from the particles, and refine3D_states '//&
        &'from a single-state solution; the state labels, the maps (docked in both hands) and the poses '//&
        &'checked against the simulation truth, with shuffled-label and duplicated-maps controls; the solve3D '//&
        &'route from the particles is held to floors, the other two are robustness profiles',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Multi-State Workflow Validation')                              ! requires sp_project
        call multistate_workflows%add_input(UI_PARM, 'suite', 'str', 'Sub-suite', &
            &'Run one workflow: cavgs, shc or states; <workflow>_identical runs it on identical states (a control), '//&
            &'hard_<workflow> on the conformational states (a robustness profile); list prints them. '//&
            &'Default (the CTest entry) runs shc, the route held to floors', '', .false., '')
        call multistate_workflows%add_input(UI_COMP, nthr, required_override=.false.)
        ! add to ui_hash
        call add_ui_program('multistate_workflows', multistate_workflows, tsttab, UI_CATEGORY)
    end subroutine new_multistate_workflows

    subroutine new_state_fixture( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call state_fixture%new(&
        &'state_fixture',&                          ! name
        &'validate the multi-state test fixture: its modes, its cache and the point groups of its truths',&
        &'The blob family in its search, known-pose and exact-view modes, each built twice (the second from the '//&
        &'cache), with exact state counts and the project each mode promises; the compositional and '//&
        &'conformational beta-galactosidase states, each truth holding its declared point group and not the next '//&
        &'larger one',&
        &'simple_test_exec',&                       ! executable
        &.false., display_name='Multi-State Fixture Validation')                               ! requires sp_project
        call state_fixture%add_input(UI_COMP, nthr, required_override=.false.)
        ! add to ui_hash
        call add_ui_program('state_fixture', state_fixture, tsttab, UI_CATEGORY)
    end subroutine new_state_fixture

    subroutine new_simulated_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        ! PROGRAM SPECIFICATION
        call simulated_workflow%new(&
        &'simulated_workflow',&                          ! name
        &'validate simulated workflows for embedded 6VXX and 1JXY systems',&
        &'runs one molecular suite through preprocessing and initial-model reconstruction, then validates the final volume',&
        &'simple_test_exec',&                            ! executable
        &.false., display_name='Simulated Workflow Validation')                                        ! requires sp_project
        ! INPUT PARAMETER SPECIFICATIONS
        ! image input/output
        !call simulated_workflow%add_input(UI_IO, )
        ! parameter input/output
        call simulated_workflow%add_input(UI_PARM, 'suite', 'multi', 'Molecular suite', &
            &'The simulated system: 6vxx (SARS-CoV-2 spike, c3) or 1jxy (beta-galactosidase, d2); '//&
            &'list prints the accepted names(6vxx|1jxy|list)', '', .true., '', &
            &choices=ui_choices([character(len=4) :: '6vxx', '1jxy', 'list']))
        call simulated_workflow%add_input(UI_PARM, 'picker', 'multi', 'Picker under test', &
            &'Particle picker used by the simulated workflow(segdiam|new){segdiam}','', .false., 'segdiam', &
        &choices=ui_choices([character(len=7) :: 'segdiam', 'new']), preserve_default=.true.)
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
        &.false., display_name='PCG Reconstruction Validation')                                    ! requires sp_project
        ! add to ui_hash
        call add_ui_program('pcg_recon', pcg_recon, tsttab, UI_CATEGORY)
    end subroutine new_pcg_recon

    subroutine new_pcg_frac_update( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call pcg_frac_update%new(&
        &'pcg_frac_update',&
        &'Project-backed PCG fractional-update validation',&
        &'Uses reconstruct3D inputs (projfile=, or without it a project of known inputs the multi-state '//&
        &'test fixture builds: 600 blob particles at their true poses, objfun=cc) and internally partitions '//&
        &'every populated state/half into two '//&
        &'deterministic complementary subsets. Validates raw (B,D) additivity, u/f continuation '//&
        &'weighting, full-mass ensemble preservation, and continuation artifact replay.',&
        &'simple_test_exec',&
        &.false., display_name='PCG Fractional Update Validation')
        call pcg_frac_update%add_input(UI_PARM, 'box_crop', 'num', 'Reconstruction box', &
        &'Even Fourier-cropped reconstruction box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0)
        call pcg_frac_update%add_input(UI_SRCH, trs)
        call pcg_frac_update%add_input(UI_SRCH, pgrp, required_override=.false.)
        call pcg_frac_update%add_input(UI_SRCH, objfun)
        call pcg_frac_update%add_input(UI_FILT, ml_reg)
        call pcg_frac_update%add_input(UI_FILT, 'postprocess', 'binary', 'Postprocess final map', &
        &'Accepted for reconstruct3D command compatibility; no production maps are written by this test', &
        &'', .false., 'no', choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call pcg_frac_update%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'Iterations used for the continuation replay comparison', 'iterations{2}', .false., 2.)
        call pcg_frac_update%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Tolerance used for the continuation replay comparison', 'tolerance{0}', .false., 0.0)
        call pcg_frac_update%add_input(UI_MASK, mskdiam, required_override=.false.)
        call pcg_frac_update%add_input(UI_COMP, nparts, required_override=.false.)
        call pcg_frac_update%add_input(UI_COMP, nthr)
        call add_ui_program('pcg_frac_update', pcg_frac_update, tsttab, UI_CATEGORY)
    end subroutine new_pcg_frac_update


    subroutine new_rec3D_backends( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call rec3D_backends%new(&
        &'rec3D_backends',&
        &'Same-inputs gridding vs PCG reconstruction comparison',&
        &'Minimal invocation is projfile+pgrp+mskdiam+nthr: defaults to objfun=euclid ml_reg=yes '//&
        &'maxits_pcg=5 rtol=1e-3. '//&
        &'Reconstructs one fixed set of particles/orientations/sigma2 with reconstruct3D using the gridding '//&
        &'and the PCG backend, in a numbered execution directory named after the run settings (run it inside '//&
        &'a refine3D run directory; the sigma2 store is the project''s registered one; '//&
        &'mkdir=no runs in the current directory instead), and compares the merged maps: per-shell amplitude '//&
        &'ratio and FSC between '//&
        &'backends, and the radial real-space profile ratio that exposes a deapodization mismatch. Writes '//&
        &'recvol_stateXX_gridding.mrc and recvol_stateXX_pcg.mrc. Measurement only, no thresholds '//&
        &'(doc/implementation_notes/completed/drop_legacy_box_division.md, plan step 2).',&
        &'simple_test_exec',&
        &.true., display_name='Reconstruction Backend Comparison')
        call rec3D_backends%add_input(UI_PARM, 'box_crop', 'num', 'Reconstruction box', &
        &'Even Fourier-cropped reconstruction box; native project geometry remains authoritative', &
        &'pixels{native box}', .false., 0.0)
        call rec3D_backends%add_input(UI_SRCH, trs)
        call rec3D_backends%add_input(UI_SRCH, pgrp)
        call rec3D_backends%add_input(UI_SRCH, objfun)
        call rec3D_backends%add_input(UI_FILT, ml_reg)
        call rec3D_backends%add_input(UI_FILT, 'maxits_pcg', 'num', 'PCG maximum iterations', &
        &'PCG iterations for the comparison', 'iterations{5}', .false., 5.)
        call rec3D_backends%add_input(UI_FILT, 'rtol', 'num', 'PCG relative residual tolerance', &
        &'Tolerance for the comparison', 'tolerance{1e-3}', .false., 1.e-3)
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

    subroutine new_single_atoms_stats( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call single_atoms_stats%new(&
        &'single_atoms_stats',&
        &'SINGLE nanoparticle atom statistics',&
        &'validates atom detection and statistics against a simulated nanoparticle',&
        &'simple_test_exec',&
        &.false., display_name='SINGLE Atom Statistics Test')
        call single_atoms_stats%add_input(UI_PARM, 'smpd', 'num', 'Sampling distance', &
        &'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms{0.358}', .false., 0.358)
        call single_atoms_stats%add_input(UI_FILT, 'element', 'str', 'Atom element name: Au, Pt etc.', &
        &'Atom element name: Au, Pt etc.', 'atom composition{Pt}', .false., 'Pt')
        call add_ui_program('single_atoms_stats', single_atoms_stats, tsttab, UI_CATEGORY)
    end subroutine new_single_atoms_stats

    subroutine new_species_discovery( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call species_discovery%new(&
        &'species_discovery',&
        &'species discovery in atom detection',&
        &'validates the species list and the residual recovery and species call of detect_atoms on Pt, Pt/Ni and Pt/Al nanoparticles rendered with element kernels and surface B factors',&
        &'simple_test_exec',&
        &.false., display_name='Species Discovery Test')
        call species_discovery%add_input(UI_COMP, nthr)
        call add_ui_program('species_discovery', species_discovery, tsttab, UI_CATEGORY)
    end subroutine new_species_discovery

    subroutine new_single_workflow( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call single_workflow%new(&
        &'single_workflow',&
        &'validate SINGLE reconstruction for FCC Pt and wurtzite CdSe',&
        &'runs focused Pt and CdSeW nanoparticle suites and validates each final map against its known truth',&
        &'simple_test_exec',&
        &.false., display_name='SINGLE Workflow Validation')
        call single_workflow%add_input(UI_PARM, 'suite', 'multi', 'Nanoparticle suite', &
        &'The simulated nanoparticle: fcc (platinum) or wurtzite (CdSe); list prints the accepted '//&
        &'names(fcc|wurtzite|list)', '', .true., '', &
        &choices=ui_choices([character(len=8) :: 'fcc', 'wurtzite', 'list']))
        call single_workflow%add_input(UI_PARM, 'smpd', 'num', 'Sampling distance', &
        &'Distance between neighbouring pixels in Angstroms', 'pixel size in Angstroms{0.358}', .false., 0.358)
        call add_ui_program('single_workflow', single_workflow, tsttab, UI_CATEGORY)
    end subroutine new_single_workflow

    subroutine new_stream_preproc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call stream_preproc%new('stream_preproc', 'stream-preproc', &
            &'generates five synthetic movies and validates streaming motion-correction and CTF outputs', &
            &'simple_test_exec', .false., display_name='Stream Preprocessing Validation')
        call add_ui_program('stream_preproc', stream_preproc, tsttab, UI_CATEGORY)
    end subroutine new_stream_preproc

    subroutine new_flex_pca_blobs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call flex_pca_blobs%new(&
        &'flex_pca_blobs',&
        &'FLEX PCA shared/distributed workflow validation',&
        &'runs one deterministic two-state phantom through shared memory and two real local workers, '//&
        &'then checks each result against truth and records cross-mode agreement; then refines each '//&
        &'truth-matched state with refine3D_auto state=X, with the frozen FLEX weights and with the hard labels; '//&
        &'then runs cls_expansion and refine3D_states_cavgs on three beta-galactosidase states in exact views '//&
        &'and checks the class-average '//&
        &'states against truth',&
        &'simple_test_exec',&
        &.false., display_name='FLEX PCA Workflow Validation')
        call flex_pca_blobs%add_input(UI_COMP, nthr)
        call add_ui_program('flex_pca_blobs', flex_pca_blobs, tsttab, UI_CATEGORY)
    end subroutine new_flex_pca_blobs

    !> manual test tools (not registered with CTest): state weight sets for weighted reconstruct3D runs
    subroutine new_write_state_weights( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call write_state_weights_labels%new(&
        &'write_state_weights_labels',&
        &'state weight set equal to the hard labels',&
        &'publishes a 0/1 state weight set equal to the project''s hard state labels, for weighted '//&
        &'reconstruct3D comparisons against hard runs',&
        &'simple_test_exec',&
        &.true., display_name='Write Hard-Label State Weights')
        call write_state_weights_labels%add_input(UI_SRCH, nstates)
        call add_ui_program('write_state_weights_labels', write_state_weights_labels, tsttab, UI_CATEGORY)
        call write_state_weights_mixed%new(&
        &'write_state_weights_mixed',&
        &'fractional state weight set from the hard labels',&
        &'publishes a fractional state weight set: 0.7 on each particle''s labelled state and 0.3 spread over '//&
        &'the other states, for fractional reconstruct3D tests',&
        &'simple_test_exec',&
        &.true., display_name='Write Mixed State Weights')
        call write_state_weights_mixed%add_input(UI_SRCH, nstates)
        call add_ui_program('write_state_weights_mixed', write_state_weights_mixed, tsttab, UI_CATEGORY)
    end subroutine new_write_state_weights

    subroutine new_nightly( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call nightly%new(&
        &'nightly',&
        &'the whole high-level suite, packed by CPU reservation',&
        &'runs every correctness entry of the test registry once, each in its own directory under dir= as a child '//&
        &'process in its own session with its deadline, packed so that the peak CPU reservations of the running '//&
        &'entries stay within ncunits= cores; writes nightly_report.tsv (one row per entry and gate metric, with '//&
        &'outcome, wall time and peak CPU use) and fails if an entry or a gate failed',&
        &'simple_test_exec',&
        &.false., display_name='Nightly Test Suite')
        call nightly%add_input(UI_PARM, 'ncunits', 'num', 'Core budget', &
            &'Cores the running entries may reserve together{the processors of the machine}', '# cores', .false., 0.)
        call nightly%add_input(UI_PARM, 'dir', 'dir', 'Output directory', &
            &'Directory of the run{nightly_<date>_<time>}', 'directory', .false., '')
        call add_ui_program('nightly', nightly, tsttab, UI_CATEGORY)
    end subroutine new_nightly

    subroutine new_profile( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call profile%new(&
        &'profile',&
        &'robustness profile of high-level test entries over seeds and thread counts',&
        &'runs a registry entry (or all) at seed offsets 0 to nseeds-1 and each thread count of nthr=, as children '//&
        &'of the test dispatcher packed within ncunits= cores, and writes per entry and thread count the '//&
        &'distribution of every gate metric (runs, median, robust spread, worst, failure rate against its floor); '//&
        &'with dir1= and dir2= it compares two profile directories instead',&
        &'simple_test_exec',&
        &.false., display_name='Test Robustness Profile')
        call profile%add_input(UI_PARM, 'target', 'str', 'Registry entry', &
            &'The registry entries to profile, comma-separated, or all', 'entry', .false., '')
        call profile%add_input(UI_PARM, 'suite', 'str', 'Sub-suite', &
            &'Replaces the suite= of the entry', 'sub-suite', .false., '')
        call profile%add_input(UI_PARM, 'nseeds', 'num', 'Seed offsets', &
            &'Number of seed offsets, 0 to nseeds-1 (0 draws what CTest draws)', '# seeds', .false., 20.)
        call profile%add_input(UI_PARM, 'nthr', 'str', 'Thread counts', &
            &'Comma-separated OMP_NUM_THREADS settings{the entry''s CTest setting}', 'e.g. 1,8', .false., '')
        call profile%add_input(UI_PARM, 'ncunits', 'num', 'Core budget', &
            &'Cores the running executions may reserve together{the processors of the machine}', '# cores', .false., 0.)
        call profile%add_input(UI_PARM, 'dir', 'dir', 'Output directory', &
            &'Directory of the profile{profile_<date>_<time>}', 'directory', .false., '')
        call profile%add_input(UI_PARM, 'keep_outputs', 'binary', 'Keep the outputs', &
            &'Keep every execution''s outputs; no deletes them as it ends (run.log and metrics.tsv stay), and the disk '//&
            &'budget counts only the running executions', '', .false., 'yes', &
            &choices=ui_choices([character(len=3) :: 'yes', 'no']))
        call profile%add_input(UI_PARM, 'dir1', 'dir', 'First profile', &
            &'A profile directory to compare with dir2=', 'directory', .false., '')
        call profile%add_input(UI_PARM, 'dir2', 'dir', 'Second profile', &
            &'A profile directory to compare with dir1=', 'directory', .false., '')
        call add_ui_program('profile', profile, tsttab, UI_CATEGORY)
    end subroutine new_profile

    subroutine new_identity( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call identity%new(&
        &'identity',&
        &'byte identity of two run directories of a high-level entry',&
        &'compares the outputs the test registry declares for the entry in two run directories byte for byte, '//&
        &'ignoring the declared run-specific fields (the run identifier, path digests and checksum of '//&
        &'solve3D_manifest.txt and the paths of its artifacts, the layout digest of sigma2_state.bin; project '//&
        &'files by the records of their particle and class segments); writes identity_report.txt and fails on '//&
        &'any difference',&
        &'simple_test_exec',&
        &.false., display_name='Test Run Identity')
        call identity%add_input(UI_PARM, 'target', 'str', 'Registry entry', &
            &'The registry entry whose declared outputs are compared', 'entry', .true., '')
        call identity%add_input(UI_PARM, 'dir1', 'dir', 'First run directory', &
            &'The run directory of the first execution', 'directory', .true., '')
        call identity%add_input(UI_PARM, 'dir2', 'dir', 'Second run directory', &
            &'The run directory of the second execution', 'directory', .true., '')
        call add_ui_program('identity', identity, tsttab, UI_CATEGORY)
    end subroutine new_identity

end module simple_test_ui_highlevel
