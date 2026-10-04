!@descr: unit tests for the steps of the stream reference-picking stage (simple_stream_stage_refpick)
! Each test assembles a stage in a fresh fixture directory from init_params, init_job_dirs,
! build_worker_cline and init_gui: no queue (init_queue is not called, and no step that submits a
! job or runs make_pickrefs is run), no waits, and a settle time of -1. The command line names no
! registered program and carries qsys_name=local for the stage project's computing environment.
! The upstream is a preprocessing directory whose completed projects hold five micrographs.
! Submission, make_pickrefs and the pick_extract jobs are left to the high-level stream tests.
module simple_stream_stage_refpick_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                  only: TERM_STREAM, METADATA_EXT
use simple_defs_stream,                 only: DIR_STREAM, DIR_STREAM_COMPLETED, STREAM_NMOVS_SET
use simple_string,                      only: string
use simple_string_utils,                only: int2str, int2str_pad
use simple_fileio,                      only: file_exists, del_file, simple_getcwd, simple_touch
use simple_syslib,                      only: dir_exists, simple_mkdir
use simple_cmdline,                     only: cmdline
use simple_sp_project,                  only: sp_project
use simple_gui_metadata_utils,          only: max_metadata_size
use simple_gui_metadata_types,          only: GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE
use simple_gui_metadata_stream_picking, only: gui_metadata_stream_picking
use simple_stream_pipe,                 only: stream_pipe
use simple_optics_maps,                 only: publish_optics_map
use simple_stream_stage_refpick,        only: stream_stage_refpick
implicit none
private
public :: run_all_stream_stage_refpick_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_refpick.simple'
character(len=*), parameter :: UPSTREAM      = 'preproc'               ! the preprocessing stage directory (dir_target)
character(len=*), parameter :: PICKREFS_IN   = 'refs/selected_refs.mrc' ! the input picking references
real,             parameter :: SMPD          = 1.3
integer,          parameter :: ALL_ACCEPTED(STREAM_NMOVS_SET) = 1

contains

    subroutine run_all_stream_stage_refpick_tests()
        write(*,'(A)') '**** running all stream reference picking stage tests ****'
        call test_init_params()
        call test_attach_upstream()
        call test_pickrefs_available()
        call test_create_set_project()
        call test_import_finished_sets()
        call test_write_project()
        call test_optics_map_applied()
        call test_without_optics_map()
        call test_restart_history()
        call test_restart_clear()
        call test_send_status()
        call test_iterate_waits()
        call test_finished()
    end subroutine run_all_stream_stage_refpick_tests

    subroutine test_init_params()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root, cwd
        integer                    :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_init_params'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_init_params', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_true(allocated(stage%params), 'initialization allocates the owned parameters')
        call simple_getcwd(cwd)
        call assert_char('stream', trim(stage%params%split_mode), 'split_mode is stream after init_params')
        call assert_char(cwd%to_char(), stage%cwd%to_char(),     'the stage directory is the working directory')
        call assert_int(0, stage%spproj%os_mic%get_noris(),      'the project starts without micrographs')
        call assert_false(stage%l_restart,                       'no output directory given: not a restart')
        call assert_true(dir_exists(string(DIR_STREAM)),          'the job folder is made')
        call assert_true(dir_exists(string(DIR_STREAM_COMPLETED)),'the completed folder is made')
        call stage%kill
        call assert_false(allocated(stage%params), 'cleanup releases the owned parameters')
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_init_params

    !> the stage waits for preprocessing's completed folder; once attached, a restart's upstream
    !! projects are in the watcher history
    subroutine test_attach_upstream()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root, upstream1
        integer                    :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_attach_upstream'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_attach', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call assert_false(stage%l_attached,      'no completed folder: not attached')
        call assert_true(stage%l_waiting_logged, 'the wait is logged')
        call make_upstream()
        upstream1 = write_upstream_project(1, ALL_ACCEPTED)
        allocate(stage%restored_sources(2))
        stage%restored_sources(1) = upstream1
        stage%restored_sources(2) = ''
        call stage%attach_upstream()
        call assert_true(stage%l_attached, 'attached once the folder exists')
        call assert_true(stage%project_buff%is_past(upstream1), 'a restart''s upstream project is in the history')
        call assert_false(allocated(stage%restored_sources),     'the restart''s list is used once')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_attach_upstream

    !> the stage waits for the input picking references
    subroutine test_pickrefs_available()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root
        integer                    :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_pickrefs_available'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_pickrefs', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%pickrefs_available(), 'no picking references yet')
        call simple_mkdir('refs')
        call simple_touch(string(PICKREFS_IN))
        call assert_true(stage%pickrefs_available(), 'the picking references are found')
        call assert_true(stage%pickrefs_available(), 'and stay found')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_pickrefs_available

    !> an upstream project becomes a job set of its accepted micrographs, without the picking
    !! preprocessing outputs; the pixel size is taken from it; a project with nothing accepted
    !! makes no set
    subroutine test_create_set_project()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(sp_project)           :: set_proj
        type(string)               :: cwd_saved, root, upstream1, upstream2, job_dir, val
        integer                    :: nfail0, nselected
        allocate(stage)
        write(*,'(A)') 'test_create_set_project'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_set_project', cwd_saved, root)
        call make_upstream()
        upstream1 = write_upstream_project(1, [1,0,1,1,1])
        upstream2 = write_upstream_project(2, [0,0,0,0,0])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        job_dir = stage%sets%get_job_dir()
        call stage%create_set_project(upstream1, nselected)
        call assert_int(4, nselected, 'the four accepted micrographs are selected')
        call assert_int(1, stage%sets%get_counter(), 'one set is written')
        call assert_int(4, stage%n_mics_submitted,   'they count as submitted')
        call assert_real(SMPD, stage%params%smpd, 1.e-6, 'the pixel size comes from the micrographs')
        val = stage%cline_exec%get_carg('projfile')
        call assert_char('00001.simple', val%to_char(), 'the worker command line names the set')
        call assert_int(4, stage%cline_exec%get_iarg('top'), 'over its four micrographs')
        call set_proj%read_segment('mic', job_dir//'/00001.simple')
        call assert_int(4, set_proj%os_mic%get_noris(), 'the set holds the accepted micrographs')
        if( set_proj%os_mic%get_noris() == 4 )then
            call assert_false(set_proj%os_mic%isthere(1, 'mic_den'), 'without the picking preprocessing outputs')
            call assert_true(set_proj%os_mic%isthere(1, 'intg'),     'with the rest of their fields')
        endif
        call set_proj%kill
        call stage%create_set_project(upstream2, nselected)
        call assert_int(0, nselected, 'nothing accepted: nothing selected')
        call assert_int(1, stage%sets%get_counter(), 'and no set written')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_create_set_project

    !> finished sets with micrographs move to the completed folder and are imported, unchanged; an
    !! empty one stays
    subroutine test_import_finished_sets()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(sp_project)           :: set_proj, done_proj
        type(string)               :: cwd_saved, root, job_dir, completed_dir, done(2)
        integer                    :: nfail0, n_imported
        allocate(stage)
        write(*,'(A)') 'test_import_finished_sets'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_import', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        job_dir       = stage%sets%get_job_dir()
        completed_dir = stage%sets%get_completed_dir()
        call make_extracted_set(set_proj, [2,3,1], 0.)
        call stage%sets%write_set(set_proj, stage%cline_exec, 3)
        call set_proj%kill
        call make_extracted_set(set_proj, [integer ::], 0.)
        call stage%sets%write_set(set_proj, stage%cline_exec, 1)
        call set_proj%kill
        done(1) = job_dir//'/00001.simple'
        done(2) = job_dir//'/00002.simple'
        call stage%import_finished_sets(done, n_imported)
        call assert_int(3, n_imported,                          'the three micrographs of set 1 are imported')
        call assert_int(3, stage%spproj%os_mic%get_noris(),     'the project holds them')
        call assert_int(6, stage%nptcls_glob,                   'their particles are counted')
        call assert_int(1, size(stage%set_projects),            'set 1 is remembered for the project write')
        call assert_true(file_exists(completed_dir//'/00001.simple'), 'set 1 moves to the completed folder')
        call assert_true(file_exists(done(2)),                         'the empty set 2 stays')
        if( file_exists(completed_dir//'/00001.simple') )then
            call done_proj%read_segment('out', completed_dir//'/00001.simple')
            call assert_int(1, done_proj%os_out%get_noris(), 'its output segment is left as the job wrote it')
            call done_proj%kill
        endif
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_finished_sets

    !> the project gets one stack per micrograph of every imported set, in import order, with the
    !! particle ranges renumbered and each particle pointing at its stack
    subroutine test_write_project()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(sp_project)           :: set_proj, written
        type(string)               :: cwd_saved, root, sets(2)
        integer, parameter         :: STKINDS(9) = [1,1,1,2,2,3,3,3,3]
        integer                    :: nfail0, n_imported, iptcl
        logical                    :: l_ok
        allocate(stage)
        write(*,'(A)') 'test_write_project'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_write_project', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        sets(1) = stage%sets%get_completed_dir()//'/00001.simple'
        sets(2) = stage%sets%get_completed_dir()//'/00002.simple'
        call make_extracted_set(set_proj, [3,2], 0.)    ! particles tagged 1..5
        call set_proj%update_projinfo(sets(1))
        call set_proj%write(sets(1))
        call set_proj%kill
        call make_extracted_set(set_proj, [4], 5.)      ! particles tagged 6..9
        call set_proj%update_projinfo(sets(2))
        call set_proj%write(sets(2))
        call set_proj%kill
        call stage%import_sets(sets, n_imported)
        call stage%write_project()
        call written%read(stage%params%projfile)
        call assert_int(3, written%os_mic%get_noris(),  'the three micrographs')
        call assert_int(3, written%os_stk%get_noris(),  'one stack per micrograph')
        if( written%os_stk%get_noris() == 3 )then
            call assert_int(1, written%os_stk%get_fromp(1), 'stack 1 starts at particle 1')
            call assert_int(3, written%os_stk%get_top(1),   'and ends at 3')
            call assert_int(4, written%os_stk%get_fromp(2), 'stack 2 continues at 4')
            call assert_int(5, written%os_stk%get_top(2),   'and ends at 5')
            call assert_int(6, written%os_stk%get_fromp(3), 'stack 3, from the second set, continues at 6')
            call assert_int(9, written%os_stk%get_top(3),   'and ends at 9')
        endif
        call assert_int(9, written%os_ptcl2D%get_noris(), 'the nine particles')
        call assert_int(9, written%os_ptcl3D%get_noris(), 'in the 3D segment as well')
        if( written%os_ptcl2D%get_noris() == 9 )then
            l_ok = .true.
            do iptcl = 1,9
                l_ok = l_ok .and. nint(written%os_ptcl2D%get(iptcl, 'dfx')) == iptcl
                l_ok = l_ok .and. written%os_ptcl2D%get_int(iptcl, 'stkind') == STKINDS(iptcl)
            enddo
            call assert_true(l_ok, 'the particles keep their order and point at their stack')
        endif
        call written%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_write_project

    !> with a map published by optics assignment, the micrographs, stacks and particles take its
    !! groups by import index, in the STAR file and in the project
    subroutine test_optics_map_applied()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(sp_project)           :: set_proj, written
        type(string)               :: cwd_saved, root, sets(1)
        integer                    :: nfail0, n_imported
        allocate(stage)
        write(*,'(A)') 'test_optics_map_applied'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_optics_map', cwd_saved, root)
        call simple_mkdir('optics')
        call publish_test_optics_map(string('optics'))
        call set_test_cline(cline)
        call cline%set('optics_dir', 'optics')
        call make_test_stage(stage, cline)
        sets(1) = stage%sets%get_completed_dir()//'/00001.simple'
        call make_extracted_set(set_proj, [2,1,3], 0., l_importind=.true.)
        call set_proj%update_projinfo(sets(1))
        call set_proj%write(sets(1))
        call set_proj%kill
        call stage%import_sets(sets, n_imported)
        call stage%write_mic_star()
        call assert_true(file_exists(string('micrographs.star')), 'the STAR file is written')
        call assert_int(2, stage%spproj%os_optics%get_noris(), 'the map''s two optics groups')
        if( stage%spproj%os_mic%get_noris() == 3 )then
            call assert_int(1, stage%spproj%os_mic%get_int(1, 'ogid'), 'micrograph 1 (import index 1) is in group 1')
            call assert_int(2, stage%spproj%os_mic%get_int(2, 'ogid'), 'micrograph 2 is in group 2')
            call assert_int(2, stage%spproj%os_mic%get_int(3, 'ogid'), 'micrograph 3 is in group 2')
        endif
        call stage%write_project()
        call written%read(stage%params%projfile)
        call assert_int(2, written%os_optics%get_noris(), 'the project holds the optics groups')
        if( written%os_stk%get_noris() == 3 .and. written%os_ptcl2D%get_noris() == 6 )then
            call assert_int(1, written%os_stk%get_int(1, 'ogid'),    'stack 1 takes its micrograph''s group')
            call assert_int(2, written%os_stk%get_int(3, 'ogid'),    'stack 3 takes its micrograph''s group')
            call assert_int(1, written%os_ptcl2D%get_int(1, 'ogid'), 'a particle of micrograph 1 is in group 1')
            call assert_int(2, written%os_ptcl2D%get_int(6, 'ogid'), 'a particle of micrograph 3 is in group 2')
            call assert_int(2, written%os_ptcl3D%get_int(6, 'ogid'), 'in the 3D segment as well')
        else
            call assert_true(.false., 'the project holds the three stacks and six particles')
        endif
        call written%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_optics_map_applied

    !> before optics assignment has published a map, every micrograph is in one optics group
    subroutine test_without_optics_map()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(sp_project)           :: set_proj
        type(string)               :: cwd_saved, root, sets(1)
        integer                    :: nfail0, n_imported, imic
        allocate(stage)
        write(*,'(A)') 'test_without_optics_map'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_no_optics_map', cwd_saved, root)
        call simple_mkdir('optics')
        call set_test_cline(cline)
        call cline%set('optics_dir', 'optics')
        call make_test_stage(stage, cline)
        sets(1) = stage%sets%get_completed_dir()//'/00001.simple'
        call make_extracted_set(set_proj, [2,1,3], 0., l_importind=.true.)
        call set_proj%update_projinfo(sets(1))
        call set_proj%write(sets(1))
        call set_proj%kill
        call stage%import_sets(sets, n_imported)
        call stage%write_mic_star()
        call assert_true(file_exists(string('micrographs.star')), 'the STAR file is written')
        call assert_int(1, stage%spproj%os_optics%get_noris(), 'one optics group')
        call assert_true(all([(stage%spproj%os_mic%get_int(imic, 'ogid') == 1, imic=1,3)]), 'every micrograph is in it')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_without_optics_map

    !> restart: the completed sets are imported again and their upstream projects are put in the
    !! watcher history; an unfinished set is dropped and its upstream project is picked again
    subroutine test_restart_history()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root, upstream1, upstream2, job_dir, done
        integer                    :: nfail0, nselected
        allocate(stage)
        write(*,'(A)') 'test_restart_history'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_restart', cwd_saved, root)
        call make_upstream()
        upstream1 = write_upstream_project(1, [1,1,1,0,0])
        upstream2 = write_upstream_project(2, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        job_dir = stage%sets%get_job_dir()
        call stage%create_set_project(upstream1, nselected)
        call stage%create_set_project(upstream2, nselected)
        call stage%sets%complete(job_dir//'/00001.simple', done) ! set 2 is left unfinished
        call stage%resume_previous_run()
        call assert_int(3, stage%spproj%os_mic%get_noris(), 'the completed set''s micrographs are imported again')
        call assert_int(1, size(stage%set_projects),        'the completed set is remembered')
        call assert_int(1, stage%sets%get_counter(),        'numbering continues after the completed set')
        call assert_false(file_exists(job_dir//'/00002.simple'), 'the unfinished set is dropped')
        call stage%attach_upstream()
        call assert_true(stage%project_buff%is_past(upstream1),  'the completed set''s upstream project is in the history')
        call assert_false(stage%project_buff%is_past(upstream2), 'the unfinished set''s upstream project is not')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_history

    !> a restart with clear=yes discards the completed sets and starts the numbering again
    subroutine test_restart_clear()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root, upstream1, job_dir, completed_dir, done
        integer                    :: nfail0, nselected
        allocate(stage)
        write(*,'(A)') 'test_restart_clear'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_restart_clear', cwd_saved, root)
        call make_upstream()
        upstream1 = write_upstream_project(1, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        job_dir       = stage%sets%get_job_dir()
        completed_dir = stage%sets%get_completed_dir()
        call stage%create_set_project(upstream1, nselected)
        call stage%sets%complete(job_dir//'/00001.simple', done)
        stage%params%clear = 'yes'
        call stage%resume_previous_run()
        call assert_int(0, stage%spproj%os_mic%get_noris(),     'nothing is imported')
        call assert_int(0, stage%sets%get_counter(),            'numbering starts again')
        call assert_false(file_exists(completed_dir//'/00001.simple'), 'the completed set is discarded')
        call assert_true(dir_exists(completed_dir),             'the completed folder is ready')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_clear

    subroutine test_send_status()
        class(stream_stage_refpick), allocatable        :: stage
        type(cmdline)                     :: cline
        type(stream_pipe)                 :: reader
        type(gui_metadata_stream_picking) :: picking
        character(len=:), allocatable     :: buffer
        type(string)                      :: cwd_saved, root, stage_name
        integer(c_int)                    :: fds(2)
        integer                           :: nfail0, nimported, naccepted, nrejected, nptcls, nppm, tlast, box, meta_type
        logical                           :: l_assigned
        allocate(stage)
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%init_gui(-1, int(fds(2)))
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%spproj%os_mic%new(3, is_ptcl=.false.)
        stage%n_mics_submitted = 7
        stage%nptcls_glob      = 42
        stage%params%box       = 128
        call stage%send_status(string('picking and extracting micrographs'))
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE, meta_type, 'it is a reference-picking status')
            if( meta_type == GUI_METADATA_STREAM_REFERENCE_PICKING_TYPE )then
                picking    = transfer(buffer, picking)
                l_assigned = picking%get(stage_name, nimported, naccepted, nrejected, nptcls, nppm, tlast, box)
                call assert_int(7,   nimported, 'micrographs submitted to picking')
                call assert_int(3,   naccepted, 'micrographs imported')
                call assert_int(42,  nptcls,    'particles extracted')
                call assert_int(128, box,       'the extraction box')
            endif
        endif
        call assert_false(reader%receive(buffer), 'one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> the public loop waits for the upstream folder, then for the picking references, submitting
    !! nothing meanwhile
    subroutine test_iterate_waits()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root, upstream1
        integer                    :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_iterate_waits'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_iterate', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no upstream output, not attached')
        call make_upstream()
        upstream1 = write_upstream_project(1, ALL_ACCEPTED)
        call stage%iterate()
        call assert_true(stage%l_attached,        'pass 2: attached')
        call assert_false(stage%l_pickrefs_found, 'pass 2: waiting for the picking references')
        call assert_int(0, stage%sets%get_counter(), 'pass 2: nothing submitted')
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_waits

    subroutine test_finished()
        class(stream_stage_refpick), allocatable :: stage
        type(cmdline)              :: cline
        type(string)               :: cwd_saved, root
        integer                    :: nfail0
        allocate(stage)
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('rp_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call simple_touch(TERM_STREAM)
        call assert_true(stage%finished(), 'the termination file finishes the stage')
        call del_file(TERM_STREAM)
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    ! the stage's command line, as its commander normalises it, under a program name outside
    ! every UI table, with a local queue system for the stage project's computing environment
    subroutine set_test_cline( cline )
        type(cmdline), intent(inout) :: cline
        call cline%set('prg',        'stream_refpick_stage_tester')
        call cline%set('mkdir',      'no')
        call cline%set('projfile',   TEST_PROJFILE)
        call cline%set('outdir',     '')
        call cline%set('nthr',       1)
        call cline%set('nparts',     2)
        call cline%set('numlen',     5)
        call cline%set('stream',     'yes')
        call cline%set('dir_target', UPSTREAM)
        call cline%set('qsys_name',  'local')
        call cline%set('pickrefs',   PICKREFS_IN)
    end subroutine set_test_cline

    ! a stage without a queue or a pipe, with no waits and a settle time that takes files written
    ! in the same second
    subroutine make_test_stage( stage, cline )
        class(stream_stage_refpick), intent(inout) :: stage
        type(cmdline),              intent(inout) :: cline
        call stage%init_params(cline)
        call stage%init_job_dirs()
        call stage%build_worker_cline(cline)
        call stage%init_gui(-1, -1)
        stage%settle_s = -1
        stage%wait_s   = 0
        stage%l_exists = .true.
    end subroutine make_test_stage

    ! the folders preprocessing makes
    subroutine make_upstream()
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
    end subroutine make_upstream

    ! completed preprocessing project number @p id: five micrographs with the given states, a
    ! picking-preprocessing output, and the pixel size; returns its absolute path
    function write_upstream_project( id, states ) result( fname )
        integer, intent(in) :: id, states(STREAM_NMOVS_SET)
        type(string)     :: fname
        type(sp_project) :: proj
        type(string)     :: cwd, stem
        integer          :: imic
        call simple_getcwd(cwd)
        fname = cwd//'/'//UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT
        call proj%os_mic%new(STREAM_NMOVS_SET, is_ptcl=.false.)
        do imic = 1,STREAM_NMOVS_SET
            stem = cwd//'/'//UPSTREAM//'/mic_'//int2str(id)//'_'//int2str(imic)
            call proj%os_mic%set_state(imic, states(imic))
            call proj%os_mic%set(imic, 'intg',    stem//'_intg.mrc')
            call proj%os_mic%set(imic, 'mic_den', stem//'_den.mrc')
            call proj%os_mic%set(imic, 'imgkind', 'mic')
            call proj%os_mic%set(imic, 'smpd',    SMPD)
        enddo
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end function write_upstream_project

    ! the project a pick_extract job leaves: one micrograph and one stack per entry of
    ! @p nptcls_mic, the particles tagged (dfx) tag_offset+1, +2, ..., and one output entry;
    ! with @p l_importind the micrographs have import indices 1, 2, ... and CTF constants
    subroutine make_extracted_set( proj, nptcls_mic, tag_offset, l_importind )
        type(sp_project),  intent(inout) :: proj
        integer,           intent(in)    :: nptcls_mic(:)
        real,              intent(in)    :: tag_offset
        logical, optional, intent(in)    :: l_importind
        integer :: nmics, imic, iptcl, fromp
        nmics = size(nptcls_mic)
        call proj%os_mic%new(nmics, is_ptcl=.false.)
        call proj%os_stk%new(nmics, is_ptcl=.false.)
        call proj%os_ptcl2D%new(sum(nptcls_mic), is_ptcl=.true.)
        fromp = 1
        do imic = 1,nmics
            call proj%os_mic%set_state(imic, 1)
            call proj%os_mic%set(imic, 'nptcls', nptcls_mic(imic))
            if( present(l_importind) )then
                if( l_importind )then
                    call proj%os_mic%set(imic, 'importind', imic)
                    call proj%os_mic%set(imic, 'smpd',      SMPD)
                    call proj%os_mic%set(imic, 'cs',        2.7)
                    call proj%os_mic%set(imic, 'kv',        300.)
                    call proj%os_mic%set(imic, 'fraca',     0.1)
                endif
            endif
            call proj%os_stk%set(imic, 'fromp',  fromp)
            call proj%os_stk%set(imic, 'top',    fromp + nptcls_mic(imic) - 1)
            do iptcl = fromp,fromp + nptcls_mic(imic) - 1
                call proj%os_ptcl2D%set(iptcl, 'dfx', tag_offset + real(iptcl))
                call proj%os_ptcl2D%set_stkind(iptcl, imic)
            enddo
            fromp = fromp + nptcls_mic(imic)
        enddo
        call proj%os_out%new(1, is_ptcl=.false.)
        call proj%os_out%set(1, 'imgkind', 'boxes')
    end subroutine make_extracted_set

    ! map 1 in @p dir, as optics assignment publishes it: import indices 1, 2, 3 in optics groups
    ! 1, 2, 2, and the two groups' rows
    subroutine publish_test_optics_map( dir )
        class(string), intent(in) :: dir
        integer, parameter :: OGIDS(3) = [1, 2, 2]
        type(sp_project) :: proj
        integer          :: imic, igroup
        call proj%os_mic%new(3, is_ptcl=.false.)
        do imic = 1,3
            call proj%os_mic%set(imic, 'importind', imic)
            call proj%os_mic%set(imic, 'ogid',      OGIDS(imic))
        enddo
        call proj%os_optics%new(2, is_ptcl=.false.)
        do igroup = 1,2
            call proj%os_optics%set(igroup, 'ogid',   igroup)
            call proj%os_optics%set(igroup, 'ogname', 'opticsgroup'//int2str(igroup))
            call proj%os_optics%set(igroup, 'state',  1)
            call proj%os_optics%set(igroup, 'pop',    count(OGIDS == igroup))
            call proj%os_optics%set(igroup, 'smpd',   SMPD)
            call proj%os_optics%set(igroup, 'cs',     2.7)
            call proj%os_optics%set(igroup, 'kv',     300.)
            call proj%os_optics%set(igroup, 'fraca',  0.1)
        enddo
        call publish_optics_map(proj, dir, 1, 5)
        call proj%kill
    end subroutine publish_test_optics_map

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

end module simple_stream_stage_refpick_tester
