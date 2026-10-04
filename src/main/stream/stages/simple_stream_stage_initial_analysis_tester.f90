!@descr: unit tests for the steps of the stream initial-analysis stage (simple_stream_stage_initial_analysis)
! Each test assembles a stage in a fresh fixture directory from init_params and init_gui only:
! no queue (init_queue is not called, and no step that submits a job is run), no waits, and a
! settle time of -1 so the upstream watcher takes projects written in the same second. The command
! line names no registered program (so params%new neither requires a project nor makes an output
! directory) and carries qsys_name=local for the computing environment of the stage's project.
! Picking, extraction, solve2D/3D, the sieve and the reprojection of the references need real
! micrographs and jobs; they are left to the high-level stream tests.
module simple_stream_stage_initial_analysis_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                    only: MRC_EXT, STK_EXT, JPG_EXT, METADATA_EXT
use simple_defs_stream,                   only: DIR_STREAM, DIR_STREAM_COMPLETED, STREAM_NMOVS_SET, OPENING2D_PICKREFS
use simple_string,                        only: string
use simple_string_utils,                  only: int2str, int2str_pad
use simple_fileio,                        only: file_exists, simple_getcwd, swap_suffix
use simple_syslib,                        only: simple_mkdir
use simple_imghead,                       only: find_ldim_nptcls
use simple_cmdline,                       only: cmdline
use simple_image,                         only: image
use simple_sp_project,                    only: sp_project
use simple_qsys_async_job,                only: ASYNC_JOB_IDLE
use simple_gui_metadata_utils,            only: max_metadata_size
use simple_gui_metadata_types,            only: GUI_METADATA_STREAM_UPDATE_TYPE, GUI_METADATA_STREAM_INITIAL_PICKING_TYPE,&
                                               &GUI_METADATA_STREAM_OPENING2D_TYPE, GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE
use simple_gui_metadata_stream_update,    only: gui_metadata_stream_update
use simple_gui_metadata_stream_picking,   only: gui_metadata_stream_picking
use simple_gui_metadata_cavg2D,           only: gui_metadata_cavg2D
use simple_stream_pipe,                   only: stream_pipe
use simple_stream_stage_initial_analysis, only: stream_stage_initial_analysis
implicit none
private
public :: run_all_stream_stage_initial_analysis_tests

character(len=*), parameter :: TEST_PROJFILE  = 'test_initial_analysis.simple'
character(len=*), parameter :: UPSTREAM       = 'preproc'   ! the preprocessing stage directory (dir_target)
real,             parameter :: SMPD           = 1.3
integer,          parameter :: NPTCLS_PER_MIC = 10
integer,          parameter :: BOX_CAVG       = 16
real,             parameter :: GRADIENT       = 0.001  ! per pixel along x, so no image is constant
real,             parameter :: PIXEL_TOL      = 0.01   ! pixel (1,1) is the class value + GRADIENT
integer,          parameter :: ALL_ACCEPTED(STREAM_NMOVS_SET) = 1
! cycle 1 steps, as in the stage
integer,          parameter :: INIT_SETUP = 0, INIT_PICK = 1

contains

    subroutine run_all_stream_stage_initial_analysis_tests()
        write(*,'(A)') '**** running all stream initial analysis stage tests ****'
        call test_init_params()
        call test_paths()
        call test_attach_upstream_waits()
        call test_import_projects()
        call test_rebuild_init_mics()
        call test_cycle1_setup_waits()
        call test_iterate_passes()
        call test_gui_selection_ends_stage()
        call test_published_pickrefs_are_final()
        call test_status_messages()
        call test_balance_classes()
        call test_find_final_solve3D_dir()
        call test_finished()
    end subroutine run_all_stream_stage_initial_analysis_tests

    !> the stage's project is made with a computing environment, and the stage directory is the
    !! working directory
    subroutine test_init_params()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(sp_project)                    :: proj
        type(string)                        :: cwd_saved, root, cwd, qsys
        integer                             :: nfail0
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call assert_false(allocated(stage%params), 'a fresh stage leaves parameters unallocated')
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'parameter initialization allocates owned state')
        call assert_int(1, stage%params%nthr, 'owned parameters retain the parsed thread count')
        call simple_getcwd(cwd)
        call assert_true(file_exists(string(TEST_PROJFILE)), 'the stage project is written')
        if( file_exists(string(TEST_PROJFILE)) )then
            call proj%read_segment('compenv', string(TEST_PROJFILE))
            qsys = proj%compenv%get_str(1, 'qsys_name')
            call assert_char('local', qsys%to_char(), 'the project carries the computing environment')
            call proj%kill
        endif
        call assert_char(cwd%to_char(), stage%cwd%to_char(), 'the stage directory is the working directory')
        call assert_int(ASYNC_JOB_IDLE, stage%job%status(), 'no job is started')
        call assert_int(1, stage%icycle, 'the stage starts in cycle 1')
        call stage%kill
        call assert_false(allocated(stage%params), 'stage cleanup releases parameters')
        call stage%kill
        call assert_false(allocated(stage%params), 'parameter cleanup is idempotent')
        call stage%init_params(cline)
        call assert_true(allocated(stage%params), 'partial initialization recreates owned parameters')
        call stage%kill
        call assert_false(allocated(stage%params), 'partial-stage cleanup releases parameters')
        call stage%spproj%kill
        call stage%cwd%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    subroutine test_paths()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root, cwd, fname, expected
        integer                             :: nfail0
        write(*,'(A)') 'test_paths'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_paths', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call simple_getcwd(cwd)
        fname    = stage%cycle_projfile(1)
        expected = cwd//'/'//DIR_STREAM//'init/init'//METADATA_EXT
        call assert_char(expected%to_char(), fname%to_char(), 'the cycle 1 project')
        fname    = stage%cycle_projfile(2)
        expected = cwd//'/'//DIR_STREAM//'all/all'//METADATA_EXT
        call assert_char(expected%to_char(), fname%to_char(), 'the cycle 2 project')
        fname    = stage%all_projfile(3)
        expected = cwd//'/'//DIR_STREAM//'all/00003'//METADATA_EXT
        call assert_char(expected%to_char(), fname%to_char(), 'the local copy of the third project of the "all" set')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_paths

    !> the stage waits, logging once, until preprocessing has made its completed-projects folder
    subroutine test_attach_upstream_waits()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root
        integer                             :: nfail0
        write(*,'(A)') 'test_attach_upstream_waits'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_attach', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no completed-projects folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call make_upstream()
        call stage%attach_upstream()
        call assert_true(stage%l_attached,       'attached once the folder exists')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_attach_upstream_waits

    !> one record per accepted micrograph of each newly completed project, each project once
    subroutine test_import_projects()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root
        integer                             :: nfail0
        write(*,'(A)') 'test_import_projects'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_import', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, [1,0,1,1,0])
        call write_completed_project(2, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_projects()
        call assert_int(8, stage%project_list%size(), 'one record per accepted micrograph')
        call assert_int(8, stage%n_mics_imported,     'the accepted micrographs are counted')
        call assert_int(8 * NPTCLS_PER_MIC, stage%n_ptcls_imported, 'their particles are counted')
        call stage%import_projects()
        call assert_int(8, stage%project_list%size(), 'a project is imported once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_projects

    !> the cycle 1 micrographs are the accepted micrographs of every imported project; rebuilding
    !! replaces them rather than appending
    subroutine test_rebuild_init_mics()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root
        integer                             :: nfail0
        write(*,'(A)') 'test_rebuild_init_mics'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_rebuild', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, [1,0,1,1,0])
        call write_completed_project(2, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_projects()
        call stage%rebuild_init_mics()
        call assert_int(8, stage%spproj%os_mic%get_noris(),         'the accepted micrographs of both projects')
        call assert_int(8, stage%spproj%os_mic%count_state_gt_zero(), 'only accepted micrographs')
        call stage%rebuild_init_mics()
        call assert_int(8, stage%spproj%os_mic%get_noris(),         'rebuilding replaces the segment')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_rebuild_init_mics

    !> cycle 1 sets up its project and reports, then waits for enough micrographs without
    !! starting a job or sending anything more
    subroutine test_cycle1_setup_waits()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(stream_pipe)                   :: reader
        character(len=:), allocatable       :: buffer
        type(string)                        :: cwd_saved, root
        integer(c_int)                      :: fds(2)
        integer                             :: nfail0
        write(*,'(A)') 'test_cycle1_setup_waits'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_cycle1', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%run_cycle1()
        call assert_int(INIT_PICK, stage%step1, 'set up; waiting for micrographs')
        call assert_true(file_exists(stage%cycle_projfile(1)), 'the cycle 1 project is written')
        call assert_int(ASYNC_JOB_IDLE, stage%job%status(), 'no extraction is started')
        call assert_int(1, stage%icycle, 'still cycle 1')
        call assert_true(reader%receive(buffer), 'a picking status is sent')
        call assert_int(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, meta_type_of(buffer), 'it is a picking status')
        call assert_true(reader%receive(buffer), 'a 2D status is sent')
        call assert_int(GUI_METADATA_STREAM_OPENING2D_TYPE, meta_type_of(buffer), 'it is a 2D status')
        call assert_false(reader%receive(buffer), 'nothing more is sent')
        call stage%run_cycle1()
        call assert_int(INIT_PICK, stage%step1, 'too few micrographs: still waiting')
        call assert_false(reader%receive(buffer), 'waiting sends nothing')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_cycle1_setup_waits

    !> the public loop: a pass while preprocessing has no output yet, then a pass that attaches,
    !! imports and sets cycle 1 up
    subroutine test_iterate_passes()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root
        integer                             :: nfail0
        write(*,'(A)') 'test_iterate_passes'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_iterate', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no upstream output, not attached')
        call assert_int(INIT_SETUP, stage%step1, 'pass 1: cycle 1 not set up')
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call stage%iterate()
        call assert_true(stage%l_attached, 'pass 2: attached')
        call assert_int(STREAM_NMOVS_SET, stage%n_mics_imported, 'pass 2: the micrographs are imported')
        call assert_int(INIT_PICK, stage%step1, 'pass 2: cycle 1 set up and waiting')
        call assert_false(stage%finished(), 'pass 2: not finished')
        call assert_false(allocated(stage%sieve), 'waiting for extraction does not allocate the sieve')
        allocate(stage%sieve)
        call stage%kill
        call assert_false(allocated(stage%sieve), 'stage cleanup releases an inactive sieve')
        call stage%kill ! idempotence
        call assert_false(allocated(stage%sieve), 'repeated cleanup leaves the sieve unallocated')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_passes

    !> a GUI selection of references for cycle 1 publishes the selected class averages as the file
    !! reference picking waits for, sends them back and ends the stage; a selection without a
    !! cycle, one for a cycle without class averages, one naming no class, and indices out of
    !! range are ignored
    subroutine test_gui_selection_ends_stage()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(stream_pipe)                   :: writer, reader
        type(gui_metadata_stream_update)    :: update
        type(gui_metadata_cavg2D)           :: cavg
        type(image)                         :: img
        character(len=:), allocatable       :: buffer
        type(string)                        :: cwd_saved, root, cwd, refs
        integer(c_int)                      :: gui_in(2), gui_out(2)
        integer                             :: nfail0, ldim(3), nrefs, imsg
        write(*,'(A)') 'test_gui_selection_ends_stage'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_selection', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call simple_getcwd(cwd)
        call simple_mkdir(cwd//'/'//DIR_STREAM)
        call simple_mkdir(cwd//'/'//DIR_STREAM//'init')
        call write_cavgs_project(stage%cycle_projfile(1), string('cavgs.mrc'), [1,1,1,1], [1,1,1,1])
        call open_loopback(gui_in)
        call open_loopback(gui_out)
        call stage%init_gui(int(gui_in(1)), int(gui_out(2)))
        call writer%new(-1, int(gui_in(2)), max_metadata_size(), 'test GUI')
        call reader%new(int(gui_out(1)), -1, max_metadata_size(), 'test GUI')
        ! a selection without a cycle
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_pickrefs_selection([1,2])
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call stage%apply_gui_updates()
        call assert_false(stage%l_done,           'a selection without a cycle is ignored')
        call assert_false(reader%receive(buffer), 'and sends nothing')
        ! cycle 2, which has no class averages yet
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_pickrefs_cycle(2)
        call update%set_pickrefs_selection([1,2])
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call stage%apply_gui_updates()
        call assert_false(stage%l_done,                                   'a selection from a cycle without class averages is ignored')
        call assert_false(file_exists(string(OPENING2D_PICKREFS)),        'and publishes nothing')
        call assert_false(reader%receive(buffer),                         'and sends nothing')
        ! no index is a class of cycle 1
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_pickrefs_cycle(1)
        call update%set_pickrefs_selection([9])
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call stage%apply_gui_updates()
        call assert_false(stage%l_done,                                   'a selection naming no class is ignored')
        call assert_false(file_exists(string(OPENING2D_PICKREFS)),        'and publishes nothing')
        call assert_false(reader%receive(buffer),                         'and sends nothing')
        ! classes 2 and 4 of cycle 1, and an index out of range
        call update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call update%set_pickrefs_cycle(1)
        call update%set_pickrefs_selection([2,4,9])
        call update%serialise(buffer)
        call writer%send(buffer)
        call update%kill
        call stage%apply_gui_updates()
        call assert_true(stage%l_done,   'the selection ends the stage')
        call assert_true(stage%finished(), 'the stage is finished')
        refs = OPENING2D_PICKREFS
        call assert_true(file_exists(refs), 'the selected class averages are published as the picking references')
        call assert_true(file_exists(swap_suffix(OPENING2D_PICKREFS, JPG_EXT, STK_EXT)), 'with their sprite sheet')
        call assert_false(file_exists(string('pickrefs_selection'//STK_EXT)), 'and no stack is left under another name')
        if( file_exists(refs) )then
            call find_ldim_nptcls(refs, ldim, nrefs)
            call assert_int(2, nrefs, 'two references: the out-of-range index is ignored')
            if( nrefs == 2 )then
                call img%new([BOX_CAVG,BOX_CAVG,1], SMPD, wthreads=.false.)
                call img%read(refs, 1)
                call assert_real(2., img%get_rmat_at(1,1,1), PIXEL_TOL, 'the first reference is class 2')
                call img%read(refs, 2)
                call assert_real(4., img%get_rmat_at(1,1,1), PIXEL_TOL, 'the second reference is class 4')
                call img%kill
            endif
        endif
        do imsg = 1,2
            call assert_true(reader%receive(buffer), 'reference '//int2str(imsg)//' is sent to the GUI')
            if( .not. allocated(buffer) ) cycle
            call assert_int(GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE, meta_type_of(buffer), 'as a picking reference')
            if( meta_type_of(buffer) /= GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE ) cycle
            cavg = transfer(buffer, cavg)
            call assert_int(imsg, cavg%get_idx(),   'indexed in the written stack')
            call assert_int(2,    cavg%get_i_max(), 'of two')
        enddo
        call assert_false(reader%receive(buffer), 'one message per reference')
        call writer%kill
        call reader%kill
        call stage%kill
        call close_loopback(gui_in)
        call close_loopback(gui_out)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_gui_selection_ends_stage

    !> published picking references are final: a restarted stage sends them to the GUI again and is
    !! finished at once, and nothing publishes over them
    subroutine test_published_pickrefs_are_final()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(stream_pipe)                   :: reader
        type(gui_metadata_cavg2D)           :: cavg
        character(len=:), allocatable       :: buffer
        type(string)                        :: cwd_saved, root, refs, later
        integer(c_int)                      :: gui_out(2)
        integer                             :: nfail0, ldim(3), nrefs, imsg
        logical                             :: l_published
        write(*,'(A)') 'test_published_pickrefs_are_final'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_pickrefs_final', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(gui_out)
        call stage%init_gui(-1, int(gui_out(2)))
        call reader%new(int(gui_out(1)), -1, max_metadata_size(), 'test GUI')
        call stage%restore_pickrefs()
        call assert_false(stage%l_done,           'nothing published: the stage runs its plan')
        call assert_false(reader%receive(buffer), 'and sends nothing')
        ! references an earlier run published
        refs = OPENING2D_PICKREFS
        call write_class_stack(refs, 2, 0.)
        call stage%restore_pickrefs()
        call assert_true(stage%finished(), 'published references: a restarted stage is finished at once')
        do imsg = 1,2
            call assert_true(reader%receive(buffer), 'reference '//int2str(imsg)//' is sent to the GUI again')
            if( .not. allocated(buffer) ) cycle
            call assert_int(GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE, meta_type_of(buffer), 'as a picking reference')
            if( meta_type_of(buffer) /= GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE ) cycle
            cavg = transfer(buffer, cavg)
            call assert_int(2, cavg%get_i_max(), 'of two')
        enddo
        call assert_false(reader%receive(buffer), 'one message per reference')
        ! a later result, from the 3D route or a selection, is not published over them
        later = 'later_refs'//STK_EXT
        call write_class_stack(later, 3, 0.)
        call stage%publish_pickrefs(later, 'THE TEST', l_published)
        call assert_false(l_published,            'references are published once per run')
        call assert_true(file_exists(later),      'the later stack is left where it is')
        call assert_false(reader%receive(buffer), 'and nothing is sent')
        call find_ldim_nptcls(refs, ldim, nrefs)
        call assert_int(2, nrefs,                 'the published references are unchanged')
        call reader%kill
        call stage%kill
        call close_loopback(gui_out)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_published_pickrefs_are_final

    !> one message per status call, carrying the stage's counts
    subroutine test_status_messages()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(stream_pipe)                   :: reader
        type(gui_metadata_stream_picking)   :: picking
        character(len=:), allocatable       :: buffer
        type(string)                        :: cwd_saved, root, stage_name
        integer(c_int)                      :: fds(2)
        integer                             :: nfail0, nimported, naccepted, nrejected, nptcls, nppm, tlast, box
        logical                             :: l_assigned
        write(*,'(A)') 'test_status_messages'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        stage%box = 64
        call stage%send_picking_status(string('picking particles'))
        call assert_true(reader%receive(buffer), 'a picking status is sent')
        if( allocated(buffer) )then
            call assert_int(GUI_METADATA_STREAM_INITIAL_PICKING_TYPE, meta_type_of(buffer), 'it is a picking status')
            if( meta_type_of(buffer) == GUI_METADATA_STREAM_INITIAL_PICKING_TYPE )then
                picking    = transfer(buffer, picking)
                l_assigned = picking%get(stage_name, nimported, naccepted, nrejected, nptcls, nppm, tlast, box)
                call assert_char('picking particles', stage_name%to_char(), 'the stage name')
                call assert_int(stage%nmics_target, nimported, 'the micrograph target of cycle 1')
                call assert_int(0,  naccepted, 'no micrograph accepted yet')
                call assert_int(64, box,       'the box')
            endif
        endif
        call stage%send_opening2D_status(string('classifying particles'), 64, 1)
        call assert_true(reader%receive(buffer), 'a 2D status is sent')
        if( allocated(buffer) ) call assert_int(GUI_METADATA_STREAM_OPENING2D_TYPE, meta_type_of(buffer), 'it is a 2D status')
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_status_messages

    !> the selected classes are replicated in proportion to their populations up to 501 rows, the
    !! remainder going to the largest fractions; the rejected class is left out, and the even/odd
    !! stacks, the cls3D segment, the output entry and the project follow
    subroutine test_balance_classes()
        integer, parameter :: STATES(4) = [1, 0, 1, 1]
        integer, parameter :: POPS(4)   = [1, 5, 3, 6]
        type(stream_stage_initial_analysis) :: stage
        type(sp_project)                    :: proj, written
        type(image)                         :: img
        type(string)                        :: cwd_saved, root, cwd, projfile, balanced, stk
        integer                             :: nfail0, ldim(3), n, ncls
        real                                :: smpd_here
        write(*,'(A)') 'test_balance_classes'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_balance', cwd_saved, root)
        call simple_getcwd(cwd)
        projfile = cwd//'/all'//METADATA_EXT
        call write_cavgs_project(projfile, string('cavgs_iter010.mrc'), STATES, POPS, l_evenodd=.true.)
        call proj%read(projfile)
        call stage%balance_classes(proj, projfile, string('balance_classes/all'))
        ! 498 extra rows over populations 1, 3 and 6: 49.8, 149.4 and 298.8 -> 49, 149 and 298, and
        ! the two left over to the two largest fractions, the 0.8s, over the 0.4 (populations whose
        ! deciding fractions tie on paper would depend on single-precision rounding)
        call assert_int(501, proj%os_cls2D%get_noris(),           'balanced to 501 rows')
        call assert_int(51,  count_class(proj, 1), 'class 1 (population 1): 1 + 49 + 1 rows')
        call assert_int(0,   count_class(proj, 2), 'the rejected class is left out')
        call assert_int(150, count_class(proj, 3), 'class 3 (population 3): 1 + 149 rows')
        call assert_int(300, count_class(proj, 4), 'class 4 (population 6): 1 + 298 + 1 rows')
        call assert_int(501, proj%os_cls2D%count_state_gt_zero(), 'every row is selected')
        call assert_int(501, proj%os_cls3D%get_noris(),           'the cls3D segment follows')
        balanced = 'balance_classes/all/cavgs_balanced'//MRC_EXT
        call assert_true(file_exists(balanced), 'the balanced stack is written')
        if( file_exists(balanced) )then
            call find_ldim_nptcls(balanced, ldim, n)
            call assert_int(501, n, 'the balanced stack has 501 images')
            if( n == 501 )then
                call img%new([BOX_CAVG,BOX_CAVG,1], SMPD, wthreads=.false.)
                call img%read(balanced, 1)
                call assert_real(1., img%get_rmat_at(1,1,1), PIXEL_TOL, 'class 1 comes first')
                call img%read(balanced, 52)
                call assert_real(3., img%get_rmat_at(1,1,1), PIXEL_TOL, 'then class 3')
                call img%read(balanced, 202)
                call assert_real(4., img%get_rmat_at(1,1,1), PIXEL_TOL, 'then class 4')
                call img%read(balanced, 501)
                call assert_real(4., img%get_rmat_at(1,1,1), PIXEL_TOL, 'class 4 last')
                call img%kill
            endif
        endif
        call check_stack_size(string('balance_classes/all/cavgs_balanced_even'//MRC_EXT), 501, 'the even stack follows')
        call check_stack_size(string('balance_classes/all/cavgs_balanced_odd'//MRC_EXT),  501, 'the odd stack follows')
        call proj%get_cavgs_stk(stk, ncls, smpd_here, fail=.false.)
        call assert_int(501, ncls, 'the class-average output entry counts 501')
        call assert_true(stk%has_substr('cavgs_balanced'//MRC_EXT), 'and names the balanced stack')
        call written%read(projfile)
        call assert_int(501, written%os_cls2D%get_noris(), 'the project is written')
        call written%kill
        call proj%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_balance_classes

    !> the result of solve3D_cavgs is in the highest-numbered '<n>_solve3D_cavgs' directory
    subroutine test_find_final_solve3D_dir()
        type(stream_stage_initial_analysis) :: stage
        type(string)                        :: cwd_saved, root, final_dir
        integer                             :: nfail0
        write(*,'(A)') 'test_find_final_solve3D_dir'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_final_dir', cwd_saved, root)
        call stage%find_final_solve3D_cavgs_dir(final_dir)
        call assert_int(0, final_dir%strlen(), 'no restart directory: empty')
        call simple_mkdir('1_solve3D_cavgs')
        call simple_mkdir('3_solve3D_cavgs')
        call simple_mkdir('10_solve3D_cavgs')
        call simple_mkdir('x_solve3D_cavgs')       ! no number
        call simple_mkdir('solve3D_cavgs')         ! no number
        call simple_mkdir('12_solve3D_cavgs_old')  ! the suffix does not end the name
        call stage%find_final_solve3D_cavgs_dir(final_dir)
        call assert_char('10_solve3D_cavgs', final_dir%to_char(), 'the highest number, compared as a number')
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_find_final_solve3D_dir

    subroutine test_finished()
        type(stream_stage_initial_analysis) :: stage
        type(cmdline)                       :: cline
        type(string)                        :: cwd_saved, root
        integer                             :: nfail0
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('ia_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call assert_false(allocated(stage%sieve), 'a fresh stage does not allocate the sieve')
        allocate(stage%sieve)
        stage%l_sieve_active = .true.
        stage%l_done = .true.
        call assert_true(stage%finished(), 'finished once the references are written')
        call stage%kill
        call assert_false(allocated(stage%sieve), 'stage cleanup releases an active sieve')
        call assert_false(allocated(stage%params), 'active-sieve cleanup also releases parameters')
        call assert_false(stage%l_sieve_active, 'stage cleanup resets sieve activity')
        call stage%kill
        call assert_false(allocated(stage%sieve), 'active-sieve cleanup remains idempotent')
        call assert_false(allocated(stage%params), 'repeated active cleanup leaves parameters unallocated')
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'a cleaned stage can reinitialize parameters')
        call assert_int(100, stage%params%nptcls_per_cls, 'reinitialization preserves parameter parsing')
        call assert_false(stage%finished(), 'reinitialization starts unfinished')
        call assert_false(allocated(stage%sieve), 'reinitialization keeps the sieve lazy')
        call stage%kill
        call assert_false(allocated(stage%params), 'reinitialized parameters are released')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    ! the stage's command line, under a program name outside every UI table, with a local queue
    ! system for the computing environment create_stream_project writes into the stage's project
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',            'stream_initial_analysis_stage_tester')
        call cline%set('mkdir',          'no')
        call cline%set('projfile',       TEST_PROJFILE)
        call cline%set('outdir',         '')
        call cline%set('nthr',           1)
        call cline%set('dir_target',     UPSTREAM)
        call cline%set('qsys_name',      'local')
        call cline%set('nptcls_per_cls', 100)
    end subroutine set_test_cline

    ! a stage from init_params and init_gui (no queue, no pipe), with no waits and a settle time
    ! that takes files written in the same second
    subroutine make_test_stage( stage, cline )
        type(stream_stage_initial_analysis), intent(inout) :: stage
        type(cmdline),                       intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s = -1
        stage%wait_s   = 0
        stage%l_exists = .true.
    end subroutine make_test_stage

    ! the folders preprocessing makes, which the stage waits for
    subroutine make_upstream()
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
    end subroutine make_upstream

    ! completed preprocessing project number @p id: five micrographs with the given states
    subroutine write_completed_project( id, states )
        integer, intent(in) :: id, states(STREAM_NMOVS_SET)
        type(sp_project) :: proj
        type(string)     :: fname, cwd, stem
        integer          :: imic
        fname = UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        call simple_getcwd(cwd)
        call proj%os_mic%new(STREAM_NMOVS_SET, is_ptcl=.false.)
        do imic = 1,STREAM_NMOVS_SET
            stem = cwd//'/'//UPSTREAM//'/mic_'//int2str(id)//'_'//int2str(imic)
            call proj%os_mic%set_state(imic, states(imic))
            call proj%os_mic%set(imic, 'intg',      stem//'_intg.mrc')
            call proj%os_mic%set(imic, 'imgkind',   'mic')
            call proj%os_mic%set(imic, 'importind', real((id - 1) * STREAM_NMOVS_SET + imic))
            call proj%os_mic%set(imic, 'smpd',      SMPD)
            call proj%os_mic%set(imic, 'nptcls',    NPTCLS_PER_MIC)
        enddo
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end subroutine write_completed_project

    ! a project with a class-average stack (image i is i plus a small gradient), and optionally its
    ! even (10+i) and odd (20+i) stacks, and a cls2D segment with the given states and populations
    subroutine write_cavgs_project( projfile, stk, states, pops, l_evenodd )
        class(string),     intent(in) :: projfile, stk
        integer,           intent(in) :: states(:), pops(:)
        logical, optional, intent(in) :: l_evenodd
        type(sp_project) :: proj
        integer          :: icls, ncls
        ncls = size(states)
        call write_class_stack(stk, ncls, 0.)
        if( present(l_evenodd) )then
            if( l_evenodd )then
                call write_class_stack(swap_ext(stk, '_even'//MRC_EXT), ncls, 10.)
                call write_class_stack(swap_ext(stk, '_odd'//MRC_EXT),  ncls, 20.)
            endif
        endif
        call proj%add_cavgs2os_out(stk, SMPD, imgkind='cavg')
        call proj%os_cls2D%new(ncls, is_ptcl=.false.)
        do icls = 1,ncls
            call proj%os_cls2D%set(icls, 'class', icls)
            call proj%os_cls2D%set_state(icls, states(icls))
            call proj%os_cls2D%set(icls, 'pop',   pops(icls))
        enddo
        call proj%update_projinfo(projfile)
        call proj%write(projfile)
        call proj%kill
    end subroutine write_cavgs_project

    ! @p n images of BOX_CAVG^2; image i is offset + i plus GRADIENT per pixel along x. Not
    ! constant: an image of one value makes the JPEG writer's (x - lo)/(hi - lo) a 0/0.
    subroutine write_class_stack( fname, n, offset )
        class(string), intent(in) :: fname
        integer,       intent(in) :: n
        real,          intent(in) :: offset
        type(image) :: img
        real        :: rmat(BOX_CAVG,BOX_CAVG,1)
        integer     :: i, ix
        call img%new([BOX_CAVG,BOX_CAVG,1], SMPD, wthreads=.false.)
        do i = 1,n
            do ix = 1,BOX_CAVG
                rmat(ix,:,1) = offset + real(i) + GRADIENT * real(ix)
            enddo
            call img%set_rmat(rmat, .false.)
            call img%write(fname, i, del_if_exists=(i == 1))
        enddo
        call img%kill
    end subroutine write_class_stack

    ! @p fname with its '.mrc' extension replaced by @p ext
    function swap_ext( fname, ext ) result( swapped )
        class(string),    intent(in) :: fname
        character(len=*), intent(in) :: ext
        type(string) :: swapped
        character(len=:), allocatable :: f
        f       = fname%to_char()
        swapped = f(1:len(f) - len(MRC_EXT))//ext
    end function swap_ext

    integer function count_class( proj, icls )
        type(sp_project), intent(inout) :: proj
        integer,          intent(in)    :: icls
        integer :: i
        count_class = 0
        do i = 1,proj%os_cls2D%get_noris()
            if( proj%os_cls2D%get_int(i, 'class') == icls ) count_class = count_class + 1
        enddo
    end function count_class

    subroutine check_stack_size( fname, nexpected, msg )
        class(string),    intent(in) :: fname
        integer,          intent(in) :: nexpected
        character(len=*), intent(in) :: msg
        integer :: ldim(3), n
        call assert_true(file_exists(fname), msg//': written')
        if( .not. file_exists(fname) ) return
        call find_ldim_nptcls(fname, ldim, n)
        call assert_int(nexpected, n, msg)
    end subroutine check_stack_size

    ! the GUI metadata type tag at the start of a serialised message
    integer function meta_type_of( buffer )
        character(len=*), intent(in) :: buffer
        meta_type_of = transfer(buffer, meta_type_of)
    end function meta_type_of

    ! a pipe with a non-blocking read end
    subroutine open_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc, flags
        fds   = -1
        rc    = c_pipe(fds)
        call assert_int(0, int(rc), 'pipe() succeeds')
        flags = c_fcntl(fds(1), F_GETFL, 0_c_int)
        rc    = c_fcntl(fds(1), F_SETFL, ior(flags, O_NONBLOCK))
        call assert_int(0, int(rc), 'the read end is made non-blocking')
    end subroutine open_loopback

    subroutine close_loopback( fds )
        integer(c_int), intent(inout) :: fds(2)
        integer(c_int) :: rc
        rc = c_close(fds(1))
        call assert_int(0, int(rc), 'the read end closes')
        rc = c_close(fds(2))
        call assert_int(0, int(rc), 'the write end closes')
        fds = -1
    end subroutine close_loopback

end module simple_stream_stage_initial_analysis_tester
