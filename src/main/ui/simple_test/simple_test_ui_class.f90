!@descr: user interfaces of the unit-test suites: the fast gate (unit_<area>), its umbrella (units), the library suites (lib_<area>) and the platform-tier forked-process suite
module simple_test_ui_class
use simple_ui_modules
implicit none

type(category_descriptor), parameter :: UI_CATEGORY = category_descriptor('class', 'Unit tests', 10)
type(ui_program), target :: units
type(ui_program), target :: unit_core
type(ui_program), target :: unit_ori
type(ui_program), target :: unit_image
type(ui_program), target :: unit_numerics
type(ui_program), target :: unit_project
type(ui_program), target :: unit_ui
type(ui_program), target :: unit_ipc
type(ui_program), target :: unit_reconstruction
type(ui_program), target :: lib_reconstruction
type(ui_program), target :: unit_pftc_align2D3D
type(ui_program), target :: unit_cart_align3D
type(ui_program), target :: lib_cart_align3D
type(ui_program), target :: unit_heterogeneity
type(ui_program), target :: unit_parallel
type(ui_program), target :: unit_single
type(ui_program), target :: lib_single
type(ui_program), target :: lib_stream
type(ui_program), target :: lib_heterogeneity
type(ui_program), target :: flex_gpu
type(ui_program), target :: openmp_offload
type(ui_program), target :: forked_process

contains

    subroutine construct_test_class_programs( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call new_units(tsttab)
        call new_unit_core(tsttab)
        call new_unit_ori(tsttab)
        call new_unit_image(tsttab)
        call new_unit_numerics(tsttab)
        call new_unit_project(tsttab)
        call new_unit_ui(tsttab)
        call new_unit_ipc(tsttab)
        call new_unit_reconstruction(tsttab)
        call new_lib_reconstruction(tsttab)
        call new_unit_pftc_align2D3D(tsttab)
        call new_unit_cart_align3D(tsttab)
        call new_lib_cart_align3D(tsttab)
        call new_unit_heterogeneity(tsttab)
        call new_lib_heterogeneity(tsttab)
        call new_unit_parallel(tsttab)
        call new_unit_single(tsttab)
        call new_lib_single(tsttab)
        call new_lib_stream(tsttab)
        call new_flex_gpu(tsttab)
        call new_openmp_offload(tsttab)
        call new_forked_process(tsttab)
    end subroutine construct_test_class_programs

    subroutine new_units( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call units%new(&
        &'units',&
        &'every fast-gate unit suite in sequence',&
        &'runs all unit_<area> suites in one process; a developer convenience, the build gate runs the area suites separately',&
        &'simple_test_exec',&
        &.false.)
        call add_ui_program('units', units, tsttab, UI_CATEGORY)
    end subroutine new_units

    subroutine new_unit_core( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_core%new(&
        &'unit_core',&
        &'unit tests: core containers, strings, file I/O and the command line',&
        &'is the fast-gate unit suite for core containers, strings, file I/O and the command line',&
        &'simple_test_exec',&
        &.false.)
        call unit_core%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (string, syslib, fileio, character_hash, hash, value_reference_hash, linked_list, record_list, command_line)', '', .false., '')
        call add_ui_program('unit_core', unit_core, tsttab, UI_CATEGORY)
    end subroutine new_unit_core

    subroutine new_unit_ori( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_ori%new(&
        &'unit_ori',&
        &'unit tests: orientations',&
        &'is the fast-gate unit suite for orientations',&
        &'simple_test_exec',&
        &.false.)
        call unit_ori%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (orientation, orientation_collection, orientation_data, euler_shift)', '', .false., '')
        call add_ui_program('unit_ori', unit_ori, tsttab, UI_CATEGORY)
    end subroutine new_unit_ori

    subroutine new_unit_image( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_image%new(&
        &'unit_image',&
        &'unit tests: images, Fourier transforms and B-spline smoothing',&
        &'is the fast-gate unit suite for images, Fourier transforms and B-spline smoothing',&
        &'simple_test_exec',&
        &.false.)
        call unit_image%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (image, image_header, fourier_iterator, b_spline_smoother_2d, b_spline_smoother_3d, masks, binary_image, segmentation, trailing_reconstruction_blend, ctf, image_serialisation)', '', .false., '')
        call add_ui_program('unit_image', unit_image, tsttab, UI_CATEGORY)
    end subroutine new_unit_image

    subroutine new_unit_numerics( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_numerics%new(&
        &'unit_numerics',&
        &'unit tests: numerics: variance, random draws, fitting, clustering',&
        &'is the fast-gate unit suite for numerics: variance, random draws, fitting, clustering',&
        &'simple_test_exec',&
        &.false.)
        call unit_numerics%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (online_variance, random_draws, affinity_propagation, hierarchical_clustering, cavg_quality_relations, diffusion_map_graphs)', '', .false., '')
        call add_ui_program('unit_numerics', unit_numerics, tsttab, UI_CATEGORY)
    end subroutine new_unit_numerics

    subroutine new_unit_project( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_project%new(&
        &'unit_project',&
        &'unit tests: projects, STAR files, class compatibility, sieving, motion gain',&
        &'is the fast-gate unit suite for projects, STAR files, class compatibility, sieving, motion gain',&
        &'simple_test_exec',&
        &.false.)
        call unit_project%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (star_file, project_merge, class_compatibility, particle_sieve, 2d_search_space_map_i/o, motion_gain)', '', .false., '')
        call add_ui_program('unit_project', unit_project, tsttab, UI_CATEGORY)
    end subroutine new_unit_project

    subroutine new_unit_ui( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_ui%new(&
        &'unit_ui',&
        &'unit tests: the UI registry and GUI metadata',&
        &'is the fast-gate unit suite for the UI registry and GUI metadata',&
        &'simple_test_exec',&
        &.false.)
        call unit_ui%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (ui_json, gui_metadata, gui_assembler, ui_hash, ui_visibility)', '', .false., '')
        call add_ui_program('unit_ui', unit_ui, tsttab, UI_CATEGORY)
    end subroutine new_unit_ui

    subroutine new_unit_ipc( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_ipc%new(&
        &'unit_ipc',&
        &'unit tests: localhost IPC: sockets, HTTP POST and persistent-worker messaging',&
        &'is the fast-gate unit suite for localhost IPC: sockets, HTTP POST and persistent-worker messaging',&
        &'simple_test_exec',&
        &.false.)
        call unit_ipc%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (ipc_tcp_socket, http_post, persistent_worker_message, persistent_worker_server)', '', .false., '')
        call add_ui_program('unit_ipc', unit_ipc, tsttab, UI_CATEGORY)
    end subroutine new_unit_ipc

    subroutine new_unit_reconstruction( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_reconstruction%new(&
        &'unit_reconstruction',&
        &'unit tests: 3D reconstruction backends and observation noise',&
        &'is the fast-gate unit suite for 3D reconstruction: the rec3D backend selector and the Gaussian observation-noise contracts',&
        &'simple_test_exec',&
        &.false.)
        call unit_reconstruction%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (rec3d_backend, observation_noise, class_average_accumulator)', '', .false., '')
        call add_ui_program('unit_reconstruction', unit_reconstruction, tsttab, UI_CATEGORY)
    end subroutine new_unit_reconstruction

    subroutine new_lib_reconstruction( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call lib_reconstruction%new(&
        &'lib_reconstruction',&
        &'library tests: half-set PCG reconstruction against gridding',&
        &'is the nightly library suite for 3D reconstruction: independent half-set PCG solves, lambda sweep and FSC against gridding',&
        &'simple_test_exec',&
        &.false.)
        call lib_reconstruction%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (pcg_half_set)', '', .false., '')
        call add_ui_program('lib_reconstruction', lib_reconstruction, tsttab, UI_CATEGORY)
    end subroutine new_lib_reconstruction

    subroutine new_unit_pftc_align2D3D( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_pftc_align2D3D%new(&
        &'unit_pftc_align2D3D',&
        &'unit tests: registration on the polar Fourier transform (2D and 3D)',&
        &'is the fast-gate unit suite for polar Fourier registration: the continuous in-plane evaluators and joint route, and the refine3D in-plane search state and policy',&
        &'simple_test_exec',&
        &.false.)
        call unit_pftc_align2D3D%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (polar_correlation, continuous_in_plane, refine3d_in_plane_state, 2d_probability_table_i/o, sigma2_state, class_average_registration)', '', .false., '')
        call add_ui_program('unit_pftc_align2D3D', unit_pftc_align2D3D, tsttab, UI_CATEGORY)
    end subroutine new_unit_pftc_align2D3D

    subroutine new_unit_cart_align3D( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_cart_align3D%new(&
        &'unit_cart_align3D',&
        &'unit tests: Cartesian (continuous) 3D registration',&
        &'is the fast-gate unit suite for Cartesian 3D registration: the neutral Cartesian Fourier layer, the five-parameter pose refiner and its refine3D adapter',&
        &'simple_test_exec',&
        &.false.)
        call unit_cart_align3D%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (cartesian_fourier, pose_refiner, pose_adapter)', '', .false., '')
        call add_ui_program('unit_cart_align3D', unit_cart_align3D, tsttab, UI_CATEGORY)
    end subroutine new_unit_cart_align3D

    subroutine new_lib_cart_align3D( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call lib_cart_align3D%new(&
        &'lib_cart_align3D',&
        &'library tests: Cartesian pose refinement on simulated 1JYX particles',&
        &'is the nightly library suite for Cartesian 3D registration: 5000 simulated 1JYX particles refined from perturbed poses and reconstructed',&
        &'simple_test_exec',&
        &.false.)
        call lib_cart_align3D%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (pose_1jyx_recovery)', '', .false., '')
        call add_ui_program('lib_cart_align3D', lib_cart_align3D, tsttab, UI_CATEGORY)
    end subroutine new_lib_cart_align3D

    subroutine new_unit_single( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_single%new(&
        &'unit_single',&
        &'unit tests: SINGLE (nanoparticles, atomic models)',&
        &'is the fast-gate unit suite for SINGLE: the atoms module and the C-alpha candidate search on a synthetic three-residue map',&
        &'simple_test_exec',&
        &.false.)
        call unit_single%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (atoms, c_alpha_finder)', '', .false., '')
        call add_ui_program('unit_single', unit_single, tsttab, UI_CATEGORY)
    end subroutine new_unit_single

    subroutine new_lib_single( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call lib_single%new(&
        &'lib_single',&
        &'library tests: SINGLE pipelines',&
        &'is the nightly library suite for SINGLE: the Pt nanoparticle atoms pipeline (simulate, detect, statistics) and the C-alpha benchmark on the built-in 6VXX and 1JYX models',&
        &'simple_test_exec',&
        &.false.)
        call lib_single%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (nanoparticle_atoms, c_alpha_molecules, pdb2mrc)', '', .false., '')
        call add_ui_program('lib_single', lib_single, tsttab, UI_CATEGORY)
    end subroutine new_lib_single

    subroutine new_lib_stream( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call lib_stream%new(&
        &'lib_stream',&
        &'library tests: in-process stream stages',&
        &'is the nightly library suite for the stream stages that run in-process: optics assignment on two beam-shift clusters, picking-reference generation and reference picking with extraction on a synthetic micrograph',&
        &'simple_test_exec',&
        &.false.)
        call lib_stream%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (optics_assignment, picking_references, pick_and_extract)', '', .false., '')
        call add_ui_program('lib_stream', lib_stream, tsttab, UI_CATEGORY)
    end subroutine new_lib_stream

    subroutine new_unit_parallel( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_parallel%new(&
        &'unit_parallel',&
        &'unit tests: distributed execution',&
        &'is the fast-gate unit suite for distributed execution: the job controller on the local backend (scripts only, nothing submitted) and the installation-path policy of the queue-system environment',&
        &'simple_test_exec',&
        &.false.)
        call unit_parallel%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (qsys_control, qsys_environment)', '', .false., '')
        call add_ui_program('unit_parallel', unit_parallel, tsttab, UI_CATEGORY)
    end subroutine new_unit_parallel

    subroutine new_unit_heterogeneity( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_heterogeneity%new(&
        &'unit_heterogeneity',&
        &'unit tests: heterogeneity analysis (flex_pca)',&
        &'is the fast-gate unit suite for flex_pca: latent model, state weights, deconvolution of 4000 particles and the PCG M-step operator at box 32 with the baseline solve',&
        &'simple_test_exec',&
        &.false.)
        call unit_heterogeneity%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (flex_pca, flex_pcg_operator)', '', .false., '')
        call add_ui_program('unit_heterogeneity', unit_heterogeneity, tsttab, UI_CATEGORY)
    end subroutine new_unit_heterogeneity

    subroutine new_lib_heterogeneity( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call lib_heterogeneity%new(&
        &'lib_heterogeneity',&
        &'library tests: flex_pca deconvolution and PCG operator',&
        &'is the nightly library suite for flex_pca: deconvolution of 20000 particles at realistic noise, the PCG M-step operator at box 64 against the exact Gram and the PCG solve sweep at box 32',&
        &'simple_test_exec',&
        &.false.)
        call lib_heterogeneity%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (flex_pca_deconvolution_20k, flex_pcg_operator_64, flex_pcg_solve_sweep)', '', .false., '')
        call add_ui_program('lib_heterogeneity', lib_heterogeneity, tsttab, UI_CATEGORY)
    end subroutine new_lib_heterogeneity

    subroutine new_flex_gpu( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call flex_gpu%new(&
        &'flex_gpu',&
        &'CUDA-C flex kernels against the CPU batch path',&
        &'compares the CUDA-C flex insertion, coupled, banked, psample and E-step kernels with the CPU path; needs a USE_FLEX_CUDA build and a device, platform label',&
        &'simple_test_exec',&
        &.false.)
        call add_ui_program('flex_gpu', flex_gpu, tsttab, UI_CATEGORY)
    end subroutine new_flex_gpu

    subroutine new_openmp_offload( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call openmp_offload%new(&
        &'openmp_offload',&
        &'OpenMP target offload on a device',&
        &'checks OpenMP target offload, data persistence and asynchronous execution, cuFFT against FFTW, cuBLAS and the KB device forms on device device= with nthr= host threads; needs a USE_OPENMP_OFFLOAD build and a device, platform label',&
        &'simple_test_exec',&
        &.false.)
        call add_ui_program('openmp_offload', openmp_offload, tsttab, UI_CATEGORY)
    end subroutine new_openmp_offload

    subroutine new_forked_process( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call forked_process%new(&
        &'forked_process',&
        &'unit tests of forked child processes',&
        &'exercises real child processes with clock-based polling (the forked-process lifecycle and the stream heartbeat); not part of the build gate, run under the platform label',&
        &'simple_test_exec',&
        &.false.)
        call forked_process%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this suite to run alone (forked_process, stream_heartbeat)', '', .false., '')
        call add_ui_program('forked_process', forked_process, tsttab, UI_CATEGORY)
    end subroutine new_forked_process

end module simple_test_ui_class
