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
type(ui_program), target :: unit_pftc_registration2D3D
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
        call new_unit_pftc_registration2D3D(tsttab)
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
            &'One sub-suite of this area to run alone (image, image_header, fourier_iterator, fourier_shift_search, b_spline_smoother_2d, b_spline_smoother_3d)', '', .false., '')
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
            &'One sub-suite of this area to run alone (online_variance, multinomial_random_draw, straight_line_fit, affinity_propagation, hierarchical_clustering)', '', .false., '')
        call add_ui_program('unit_numerics', unit_numerics, tsttab, UI_CATEGORY)
    end subroutine new_unit_numerics

    subroutine new_unit_project( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_project%new(&
        &'unit_project',&
        &'unit tests: projects, STAR files, class compatibility, sieving, motion gain, atoms',&
        &'is the fast-gate unit suite for projects, STAR files, class compatibility, sieving, motion gain, atoms',&
        &'simple_test_exec',&
        &.false.)
        call unit_project%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (star_file, project_merge, class_compatibility, particle_sieve, 2d_search_space_map_i/o, motion_gain, atoms)', '', .false., '')
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
            &'One sub-suite of this area to run alone (ui_json, gui_metadata, gui_assembler)', '', .false., '')
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
            &'One sub-suite of this area to run alone (rec3d_backend, observation_noise)', '', .false., '')
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

    subroutine new_unit_pftc_registration2D3D( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call unit_pftc_registration2D3D%new(&
        &'unit_pftc_registration2D3D',&
        &'unit tests: registration on the polar Fourier transform (2D and 3D)',&
        &'is the fast-gate unit suite for polar Fourier registration: the continuous in-plane evaluators and joint route, and the refine3D in-plane search state and policy',&
        &'simple_test_exec',&
        &.false.)
        call unit_pftc_registration2D3D%add_input(UI_PARM, 'suite', 'str', 'Run one sub-suite', &
            &'One sub-suite of this area to run alone (continuous_in_plane, refine3d_in_plane_state)', '', .false., '')
        call add_ui_program('unit_pftc_registration2D3D', unit_pftc_registration2D3D, tsttab, UI_CATEGORY)
    end subroutine new_unit_pftc_registration2D3D

    subroutine new_forked_process( tsttab )
        class(ui_hash), intent(inout) :: tsttab
        call forked_process%new(&
        &'forked_process',&
        &'unit tests of forked child processes',&
        &'exercises real child processes with clock-based polling; not part of the build gate, run under the platform label',&
        &'simple_test_exec',&
        &.false.)
        call add_ui_program('forked_process', forked_process, tsttab, UI_CATEGORY)
    end subroutine new_forked_process

end module simple_test_ui_class
