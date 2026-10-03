!@descr: unit tests for the steps of the stream optics-assignment stage (simple_stream_stage_optics)
! Each test builds a stage with its own new() in a fresh fixture directory, from a command line
! that names no registered program (so params%new neither requires a project nor makes an output
! directory), with no waits and a settle time of -1 so the upstream watcher takes projects written
! in the same second. The upstream is a preprocessing directory with spprojs/ and
! spprojs_completed/; its projects hold five micrographs in two beam-shift clusters, near (0,0)
! twice and near (5,5) three times, as in the stream optics test of simple_stream_tester.
module simple_stream_stage_optics_tester
use, intrinsic :: iso_c_binding, only: c_int
use unix,                        only: c_pipe, c_close, c_fcntl, F_GETFL, F_SETFL, O_NONBLOCK
use simple_test_utils
use simple_defs_fname,                            only: TERM_STREAM, METADATA_EXT, OPTICS_MAP_PREFIX, TXT_EXT
use simple_defs_stream,                           only: DIR_STREAM, DIR_STREAM_COMPLETED, STREAM_NMOVS_SET
use simple_string,                                only: string
use simple_string_utils,                          only: int2str, int2str_pad
use simple_fileio,                                only: file_exists, del_file, simple_getcwd, simple_touch
use simple_syslib,                                only: simple_mkdir
use simple_cmdline,                               only: cmdline
use simple_sp_project,                            only: sp_project
use simple_gui_metadata_utils,                    only: max_metadata_size
use simple_gui_metadata_types,                    only: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE,&
                                                       &GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE
use simple_gui_metadata_optics_group,             only: gui_metadata_optics_group
use simple_gui_metadata_stream_optics_assignment, only: gui_metadata_stream_optics_assignment
use simple_stream_pipe,                           only: stream_pipe
use simple_stream_stage_optics,                   only: stream_stage_optics
implicit none
private
public :: run_all_stream_stage_optics_tests

character(len=*), parameter :: TEST_PROJFILE = 'test_optics_stage.simple'
character(len=*), parameter :: UPSTREAM      = 'preproc'   ! the preprocessing stage directory (dir_target)
real,             parameter :: SMPD = 1.3, CS = 2.7, KV = 300.0, FRACA = 0.1
real,             parameter :: TILT_THRES = 0.5
real,             parameter :: SHIFT_X(STREAM_NMOVS_SET)    = [0.0, 0.1, 5.0, 5.1, 4.9]
real,             parameter :: SHIFT_Y(STREAM_NMOVS_SET)    = [0.0,-0.1, 5.0, 4.9, 5.1]
integer,          parameter :: TILT_GROUP(STREAM_NMOVS_SET) = [1, 1, 2, 2, 2]
integer,          parameter :: ALL_ACCEPTED(STREAM_NMOVS_SET) = 1

contains

    subroutine run_all_stream_stage_optics_tests()
        write(*,'(A)') '**** running all stream optics assignment stage tests ****'
        call test_new_starts_empty()
        call test_restart_removes_termination_file()
        call test_attach_upstream_waits()
        call test_import_new_projects()
        call test_assign_and_publish()
        call test_beamtilt_from_command_line()
        call test_map_ids_continue_after_restart()
        call test_iterate_passes()
        call test_send_status()
        call test_send_group_shifts()
        call test_finished()
    end subroutine run_all_stream_stage_optics_tests

    !> a restart imports every upstream project again: the micrograph segment starts empty even
    !! when the stage's project holds micrographs; no map yet, and maps live in the stage directory
    subroutine test_new_starts_empty()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root, cwd
        integer                   :: nfail0
        write(*,'(A)') 'test_new_starts_empty'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_new', cwd_saved, root)
        call write_mic_project(string(TEST_PROJFILE), 1, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call simple_getcwd(cwd)
        call assert_int(0, stage%spproj%os_mic%get_noris(), 'the micrograph segment starts empty')
        call assert_int(0, stage%map_id,                    'no optics map in the stage directory')
        call assert_char(cwd%to_char(), stage%map_dir%to_char(), 'maps are published in the stage directory')
        call assert_false(stage%l_attached,                 'not attached to the upstream before the first pass')
        call assert_false(stage%finished(),                 'a fresh stage is not finished')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_new_starts_empty

    !> a restart (the output directory exists) removes a termination file left by the previous
    !! run, which would otherwise stop the stage at once. A unit test cannot enter a run
    !! directory (mkdir acts only for programs in the UI tables), so this checks the removal in
    !! the directory finished() reads, not the move into the stage directory itself.
    subroutine test_restart_removes_termination_file()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        write(*,'(A)') 'test_restart_removes_termination_file'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_restart_term', cwd_saved, root)
        call simple_mkdir('previous_run')
        call simple_touch(TERM_STREAM)
        call set_test_cline(cline)
        call cline%set('outdir', 'previous_run')
        call make_test_stage(stage, cline)
        call assert_false(file_exists(TERM_STREAM), 'the termination file is removed')
        call assert_false(stage%finished(),         'the restarted stage is not finished')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restart_removes_termination_file

    !> the stage waits, logging once, until preprocessing has made its completed-projects folder
    subroutine test_attach_upstream_waits()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        write(*,'(A)') 'test_attach_upstream_waits'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_attach', cwd_saved, root)
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

    !> every micrograph of a newly completed project is imported, rejected ones included; a
    !! project is imported once, and projects completed later are picked up by a later pass
    subroutine test_import_new_projects()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0, nimported
        write(*,'(A)') 'test_import_new_projects'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_import', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call write_completed_project(2, [1,0,1,0,1])
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call assert_int(2*STREAM_NMOVS_SET, nimported,                  'both projects are imported')
        call assert_int(2*STREAM_NMOVS_SET, stage%spproj%os_mic%get_noris(), 'the segment holds all their micrographs')
        call assert_int(2*STREAM_NMOVS_SET - 2, stage%spproj%os_mic%count_state_gt_zero(),&
            &'rejected micrographs are imported and stay rejected')
        call stage%import_new_projects(nimported)
        call assert_int(0, nimported, 'a project is imported once')
        call write_completed_project(3, ALL_ACCEPTED)
        call stage%import_new_projects(nimported)
        call assert_int(STREAM_NMOVS_SET, nimported, 'a project completed later is imported by a later pass')
        call assert_int(3*STREAM_NMOVS_SET, stage%spproj%os_mic%get_noris(), 'the segment grows')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_import_new_projects

    !> grouping the two shift clusters: two optics groups, the STAR files, the project with the
    !! groups, and the first optics map
    subroutine test_assign_and_publish()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(sp_project)          :: proj
        type(string)              :: cwd_saved, root
        integer                   :: nfail0, nimported
        write(*,'(A)') 'test_assign_and_publish'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_assign', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call stage%assign_and_publish()
        call check_two_groups(stage%spproj)
        call assert_true(file_exists(string('optics.star')),      'optics.star is written')
        call assert_true(file_exists(string('micrographs.star')), 'micrographs.star is written')
        call assert_int(1, stage%map_id, 'the first map is number 1')
        call assert_true(file_exists(map_file(stage%map_dir, 1)), 'the first optics map is published')
        call assert_true(file_exists(stage%map_dir//'/'//OPTICS_MAP_PREFIX//'1'//METADATA_EXT), 'with its optics segment')
        call assert_false(file_exists(map_file(stage%map_dir, 1)//'.tmp'),                    'and no table left under a temporary name')
        call assert_false(file_exists(stage%map_dir//'/'//OPTICS_MAP_PREFIX//'1.tmp'),         'nor a segment')
        call proj%read_segment('optics', stage%params%projfile)
        call assert_int(2, proj%os_optics%get_noris(), 'the project holds the two groups')
        call proj%read_segment('mic', stage%params%projfile)
        call assert_int(STREAM_NMOVS_SET, proj%os_mic%get_noris(), 'the project holds the micrographs')
        call proj%kill
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_assign_and_publish

    !> beamtilt=yes on the command line groups by tilt group first: five micrographs with the
    !! same beam-image shift form one optics group per tilt group, and one group without it
    subroutine test_beamtilt_from_command_line()
        real, parameter :: NO_SHIFT(STREAM_NMOVS_SET) = 0.
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0, nimported
        write(*,'(A)') 'test_beamtilt_from_command_line'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_beamtilt', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED, NO_SHIFT, NO_SHIFT)
        ! with beam tilt
        call set_test_cline(cline)
        call cline%set('beamtilt', 'yes')
        call make_test_stage(stage, cline)
        call assert_char('yes', trim(stage%params%beamtilt), 'beamtilt reaches the parameters')
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call stage%assign_and_publish()
        call assert_int(2, stage%spproj%os_optics%get_noris(), 'with beam tilt: one group per tilt group')
        if( stage%spproj%os_optics%get_noris() == 2 )then
            call assert_true(stage%spproj%os_mic%get_int(1, 'ogid') == stage%spproj%os_mic%get_int(2, 'ogid'),&
                &'with beam tilt: tilt group 1 is one optics group')
            call assert_true(stage%spproj%os_mic%get_int(1, 'ogid') /= stage%spproj%os_mic%get_int(3, 'ogid'),&
                &'with beam tilt: the tilt groups are different optics groups')
        endif
        call stage%kill
        call cline%kill
        ! without (the default)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_char('no', trim(stage%params%beamtilt), 'beamtilt defaults to no')
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call stage%assign_and_publish()
        call assert_int(1, stage%spproj%os_optics%get_noris(), 'without beam tilt: one group for the one shift')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_beamtilt_from_command_line

    !> a restarted stage continues the map ids from the newest map in its directory, so readers
    !! never take an older map for the newest; the older maps are kept
    subroutine test_map_ids_continue_after_restart()
        type(stream_stage_optics) :: stage, restarted
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0, nimported
        write(*,'(A)') 'test_map_ids_continue_after_restart'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_restart', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call stage%assign_and_publish()
        call stage%assign_and_publish()
        call assert_int(2, stage%map_id, 'two maps published')
        call stage%kill
        call make_test_stage(restarted, cline)
        call assert_int(2, restarted%map_id, 'the restarted stage starts from the newest map')
        call assert_int(0, restarted%spproj%os_mic%get_noris(), 'the restarted stage starts with no micrographs')
        call restarted%attach_upstream()
        call restarted%import_new_projects(nimported)
        call assert_int(STREAM_NMOVS_SET, nimported, 'the completed project is imported again')
        call restarted%assign_and_publish()
        call assert_int(3, restarted%map_id, 'the next map is number 3')
        call assert_true(file_exists(map_file(restarted%map_dir, 3)), 'map 3 is published')
        call assert_true(file_exists(map_file(restarted%map_dir, 1)), 'older maps are kept')
        call restarted%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_map_ids_continue_after_restart

    !> the public loop: a pass while preprocessing has no output yet, a pass that attaches, imports
    !! and publishes, an idle pass, then finished once nmics micrographs are in; finalize writes
    !! the project
    subroutine test_iterate_passes()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(sp_project)          :: proj
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        write(*,'(A)') 'test_iterate_passes'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_iterate', cwd_saved, root)
        call simple_mkdir(UPSTREAM)
        call set_test_cline(cline, nmics=STREAM_NMOVS_SET)
        call make_test_stage(stage, cline)
        call stage%iterate()
        call assert_false(stage%l_attached, 'pass 1: no upstream output, not attached')
        call assert_false(stage%finished(), 'pass 1: not finished')
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call stage%iterate()
        call assert_true(stage%l_attached, 'pass 2: attached')
        call assert_int(STREAM_NMOVS_SET, stage%spproj%os_mic%get_noris(), 'pass 2: the micrographs are imported')
        call assert_int(2, stage%spproj%os_optics%get_noris(), 'pass 2: two optics groups')
        call assert_int(1, stage%map_id, 'pass 2: the first map is published')
        call assert_true(stage%finished(), 'pass 2: nmics reached, finished')
        call stage%iterate()
        call assert_int(1, stage%map_id, 'pass 3: nothing new, no new map')
        call assert_int(STREAM_NMOVS_SET, stage%spproj%os_mic%get_noris(), 'pass 3: nothing new imported')
        call stage%finalize()
        call proj%read(stage%params%projfile)
        call assert_int(STREAM_NMOVS_SET, proj%os_mic%get_noris(),  'finalize: the project holds the micrographs')
        call assert_int(2,                proj%os_optics%get_noris(), 'finalize: the project holds the groups')
        call proj%kill
        call stage%kill
        call stage%kill ! idempotence
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_iterate_passes

    !> one status message per call, counting accepted micrographs and optics groups
    subroutine test_send_status()
        type(stream_stage_optics)                   :: stage
        type(cmdline)                               :: cline
        type(stream_pipe)                           :: reader
        type(gui_metadata_stream_optics_assignment) :: status
        character(len=:), allocatable               :: buffer
        type(string)                                :: cwd_saved, root, stage_name
        integer(c_int)                              :: fds(2)
        integer                                     :: nfail0, imic, meta_type, nassigned, ngroups, nimported, tlast
        logical                                     :: l_assigned
        write(*,'(A)') 'test_send_status'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_status', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call open_loopback(fds)
        call stage%pipe%new(-1, int(fds(2)), max_metadata_size(), 'optics stage under test')
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%spproj%os_mic%new(3, is_ptcl=.false.)
        do imic = 1,3
            call stage%spproj%os_mic%set_state(imic, merge(0, 1, imic == 2))
        enddo
        call stage%spproj%os_optics%new(2, is_ptcl=.false.)
        call stage%send_status()
        call assert_true(reader%receive(buffer), 'a status message is sent')
        if( allocated(buffer) )then
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE, meta_type, 'it is an optics-assignment status')
            if( meta_type == GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE )then
                status     = transfer(buffer, status)
                l_assigned = status%get(stage_name, nassigned, ngroups, tlast, nimported)
                call assert_int(2, nassigned, 'micrographs assigned: the accepted ones')
                call assert_int(2, ngroups,   'optics groups assigned')
                call assert_int(2, nimported, 'micrographs imported: the accepted ones')
            endif
        endif
        call assert_false(reader%receive(buffer), 'exactly one message per call')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_status

    !> one message per optics group, each with the shifts of that group's micrographs
    subroutine test_send_group_shifts()
        type(stream_stage_optics)       :: stage
        type(cmdline)                   :: cline
        type(stream_pipe)               :: reader
        type(gui_metadata_optics_group) :: group
        character(len=:), allocatable   :: buffer
        real,             allocatable   :: xs(:), ys(:)
        type(string)                    :: cwd_saved, root
        integer(c_int)                  :: fds(2)
        integer                         :: nfail0, nimported, iframe, meta_type, i, i_max, n_shifts, npoints(2)
        logical                         :: l_assigned
        write(*,'(A)') 'test_send_group_shifts'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_shifts', cwd_saved, root)
        call make_upstream()
        call write_completed_project(1, ALL_ACCEPTED)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call stage%attach_upstream()
        call stage%import_new_projects(nimported)
        call stage%assign_and_publish()
        call open_loopback(fds)
        call stage%pipe%new(-1, int(fds(2)), max_metadata_size(), 'optics stage under test')
        call reader%new(int(fds(1)), -1, max_metadata_size(), 'test reader')
        call stage%send_group_shifts()
        npoints = 0
        do iframe = 1,2
            call assert_true(reader%receive(buffer), 'a message for group '//int2str(iframe))
            if( .not. allocated(buffer) ) cycle
            meta_type = transfer(buffer, meta_type)
            call assert_int(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE, meta_type, 'it is an optics-group message')
            if( meta_type /= GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE ) cycle
            group      = transfer(buffer, group)
            l_assigned = group%get(i, i_max, xs, ys, n_shifts)
            call assert_int(2, i_max, 'each message counts the two groups')
            if( i >= 1 .and. i <= 2 ) npoints(i) = n_shifts
        enddo
        call assert_false(reader%receive(buffer), 'one message per group, no more')
        call assert_int(STREAM_NMOVS_SET, sum(npoints),  'every micrograph shift is sent once')
        call assert_int(2, minval(npoints), 'the (0,0) group sends two shifts')
        call assert_int(3, maxval(npoints), 'the (5,5) group sends three shifts')
        call reader%kill
        call stage%kill
        call close_loopback(fds)
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_send_group_shifts

    !> the stage finishes on the stream's termination file or once enough micrographs are imported
    subroutine test_finished()
        type(stream_stage_optics) :: stage
        type(cmdline)             :: cline
        type(string)              :: cwd_saved, root
        integer                   :: nfail0
        write(*,'(A)') 'test_finished'
        nfail0 = tests_failed
        call enter_fixture('optics_stage_finished', cwd_saved, root)
        call set_test_cline(cline)
        call make_test_stage(stage, cline)
        call assert_false(stage%finished(), 'a fresh stage is not finished')
        call simple_touch(TERM_STREAM)
        call assert_true(stage%finished(), 'the termination file finishes the stage')
        call del_file(TERM_STREAM)
        stage%l_nmics_reached = .true.
        call assert_true(stage%finished(), 'reaching the requested micrographs finishes the stage')
        call stage%kill
        call cline%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_finished

    ! ---- fixtures ------------------------------------------------------------

    ! the stage's command line, as its commander normalises it, under a program name outside
    ! every UI table, so params%new neither requires a project nor makes an output directory
    subroutine set_test_cline( cline, nmics )
        type(cmdline),     intent(inout) :: cline
        integer, optional, intent(in)    :: nmics
        call cline%set('prg',        'stream_optics_stage_tester')
        call cline%set('mkdir',      'no')
        call cline%set('projfile',   TEST_PROJFILE)
        call cline%set('outdir',     '')
        call cline%set('nthr',       1)
        call cline%set('dir_target', UPSTREAM)
        call cline%set('tilt_thres', TILT_THRES)
        if( present(nmics) ) call cline%set('nmics', nmics)
    end subroutine set_test_cline

    ! a stage from its own new(), with no waits and a settle time that takes files written in
    ! the same second; the project file is made first, so new() does not build one from the
    ! environment
    subroutine make_test_stage( stage, cline )
        type(stream_stage_optics), intent(inout) :: stage
        type(cmdline),             intent(inout) :: cline
        type(sp_project) :: proj
        if( .not. file_exists(string(TEST_PROJFILE)) )then
            call proj%update_projinfo(string(TEST_PROJFILE))
            call proj%write(string(TEST_PROJFILE))
            call proj%kill
        endif
        call stage%new(cline)
        stage%settle_s = -1
        stage%wait_s   = 0
    end subroutine make_test_stage

    ! the folders preprocessing makes, which the stage waits for
    subroutine make_upstream()
        call simple_mkdir(UPSTREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM)
        call simple_mkdir(UPSTREAM//'/'//DIR_STREAM_COMPLETED)
    end subroutine make_upstream

    ! completed preprocessing project number @p id, as preprocessing moves it into spprojs_completed
    subroutine write_completed_project( id, states, shiftx, shifty )
        integer,        intent(in) :: id, states(STREAM_NMOVS_SET)
        real, optional, intent(in) :: shiftx(STREAM_NMOVS_SET), shifty(STREAM_NMOVS_SET)
        call write_mic_project(string(UPSTREAM//'/'//DIR_STREAM_COMPLETED//int2str_pad(id, 5)//METADATA_EXT), id, states,&
            &shiftx, shifty)
    end subroutine write_completed_project

    ! a project of five preprocessed micrographs, by default in the two shift clusters
    subroutine write_mic_project( fname, id, states, shiftx, shifty )
        class(string),  intent(in) :: fname
        integer,        intent(in) :: id, states(STREAM_NMOVS_SET)
        real, optional, intent(in) :: shiftx(STREAM_NMOVS_SET), shifty(STREAM_NMOVS_SET)
        type(sp_project) :: proj
        type(string)     :: cwd, stem
        real             :: xs(STREAM_NMOVS_SET), ys(STREAM_NMOVS_SET)
        integer          :: imic
        xs = SHIFT_X
        ys = SHIFT_Y
        if( present(shiftx) ) xs = shiftx
        if( present(shifty) ) ys = shifty
        call simple_getcwd(cwd)
        call proj%os_mic%new(STREAM_NMOVS_SET, is_ptcl=.false.)
        do imic = 1,STREAM_NMOVS_SET
            stem = cwd//'/'//UPSTREAM//'/mic_'//int2str(id)//'_'//int2str(imic)
            call proj%os_mic%set_state(imic, states(imic))
            call proj%os_mic%set(imic, 'movie',     stem//'.mrc')
            call proj%os_mic%set(imic, 'intg',      stem//'_intg.mrc')
            call proj%os_mic%set(imic, 'imgkind',   'mic')
            call proj%os_mic%set(imic, 'importind', real((id - 1) * STREAM_NMOVS_SET + imic))
            call proj%os_mic%set(imic, 'smpd',      SMPD)
            call proj%os_mic%set(imic, 'cs',        CS)
            call proj%os_mic%set(imic, 'kv',        KV)
            call proj%os_mic%set(imic, 'fraca',     FRACA)
            call proj%os_mic%set(imic, 'dfx',       1.5)
            call proj%os_mic%set(imic, 'dfy',       1.6)
            call proj%os_mic%set(imic, 'angast',    0.0)
            call proj%os_mic%set(imic, 'phshift',   0.0)
            call proj%os_mic%set(imic, 'ctfres',    4.0)
            call proj%os_mic%set(imic, 'shiftx',    xs(imic))
            call proj%os_mic%set(imic, 'shifty',    ys(imic))
            call proj%os_mic%set(imic, 'tiltgrp',   real(TILT_GROUP(imic)))
        enddo
        call proj%update_projinfo(fname)
        call proj%write(fname)
        call proj%kill
    end subroutine write_mic_project

    ! the expectations of the stream optics test: micrographs 1-2 and 3-5 in two different
    ! groups, of two and three micrographs
    subroutine check_two_groups( spproj )
        type(sp_project), intent(inout) :: spproj
        integer :: group_a, group_b
        logical :: l_groups
        call assert_int(2, spproj%os_optics%get_noris(), 'two optics groups')
        if( spproj%os_optics%get_noris() /= 2 .or. spproj%os_mic%get_noris() /= STREAM_NMOVS_SET ) return
        group_a  = spproj%os_mic%get_int(1, 'ogid')
        group_b  = spproj%os_mic%get_int(3, 'ogid')
        l_groups = group_a /= group_b .and. all([group_a, group_b] >= 1) .and. all([group_a, group_b] <= 2)
        call assert_true(l_groups, 'the two shift clusters are in different groups')
        if( .not. l_groups ) return
        call assert_int(group_a, spproj%os_mic%get_int(2, 'ogid'), 'micrograph 2 is in the (0,0) group')
        call assert_int(group_b, spproj%os_mic%get_int(4, 'ogid'), 'micrograph 4 is in the (5,5) group')
        call assert_int(group_b, spproj%os_mic%get_int(5, 'ogid'), 'micrograph 5 is in the (5,5) group')
        call assert_int(2, nint(spproj%os_optics%get(group_a, 'pop')), 'the (0,0) group holds two micrographs')
        call assert_int(3, nint(spproj%os_optics%get(group_b, 'pop')), 'the (5,5) group holds three micrographs')
    end subroutine check_two_groups

    function map_file( dir, id ) result( fname )
        class(string), intent(in) :: dir
        integer,       intent(in) :: id
        type(string) :: fname
        fname = dir//'/'//OPTICS_MAP_PREFIX//int2str(id)//TXT_EXT
    end function map_file

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

end module simple_stream_stage_optics_tester
