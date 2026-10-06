!@descr: unit tests for the pieces of the stream master: the GUI's commands, the metadata store, the stages' pipes
! The GUI's answer is parsed from JSON text, the store is fed serialised metadata objects, and a
! stage's pipes are opened and used from both sides in this process. run_stream_heartbeat_tests
! reads live forked processes into the GUI heartbeat's records, so it runs in the forked_process
! entry. Forking the stages and the HTTP link are left to the high-level stream tests.
module simple_stream_master_tester
use simple_test_utils
use simple_string,                     only: string
use simple_cmdline,                    only: cmdline
use simple_commander_base,             only: commander_base
use simple_gui_metadata_api,           only: gui_metadata_cavg2D, gui_metadata_vol3D, gui_metadata_stream_pool2D,&
                                            &sprite_sheet_pos, GUI_METADATA_STREAM_POOL2D_TYPE,&
                                            &GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                            &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE
use simple_gui_metadata_utils,         only: max_metadata_size
use simple_stream_pipe,                only: stream_pipe
use simple_stream_master_stage_ids,    only: NSTAGES, STAGE_PREPROCESS, STAGE_POOL2D, STAGE_SOLVE3D, stage_gui_key,&
                                            &stage_job_name, master_fds, stage_fds
use simple_stream_master_stage,        only: stream_master_stage, fork_gui_status
use simple_gui_assembler,              only: gui_stage_status, GUI_STAGE_STATUS_RUNNING, GUI_STAGE_STATUS_FINISHED
use simple_forked_process,             only: forked_process, FORK_POLL_TIME
use unix,                              only: c_usleep
use simple_stream_master_meta_store,   only: stream_master_meta_store
use simple_stream_master_gui_commands, only: stream_master_gui_commands
use simple_gui_metadata_stream_update, only: MAX_PICKREFS_SELECTION, MAX_SNAPSHOT2D_SELECTION, MAX_SNAPSHOT2D_FNAME_LEN
use simple_string_utils,               only: int2str
implicit none
private
public :: run_all_stream_master_tests, run_stream_heartbeat_tests

!> a stage's commander for the tests, which never fork the stage
type, extends(commander_base) :: noop_commander
contains
    procedure :: execute => execute_noop
end type noop_commander

contains

    subroutine run_all_stream_master_tests()
        write(*,'(A)') '**** running all stream master tests ****'
        call test_stage_names()
        call test_gui_commands()
        call test_gui_commands_fresh_each_answer()
        call test_gui_commands_invalid()
        call test_gui_commands_oversized()
        call test_gui_commands_snapshot_name()
        call test_store_status()
        call test_store_list()
        call test_store_volume()
        call test_store_clear_stage()
        call test_stage_pipes()
        call test_update_dedupe()
        call test_skipped_stays_skipped()
    end subroutine run_all_stream_master_tests

    subroutine test_stage_names()
        write(*,'(A)') 'test_stage_names'
        call assert_char('preprocess',            stage_gui_key(STAGE_PREPROCESS), 'the GUI key of preprocessing')
        call assert_char('pool2D',                stage_gui_key(STAGE_POOL2D),     'of pool 2D')
        call assert_char('solve3D_multistate', stage_gui_key(STAGE_SOLVE3D), 'of multistate 3D')
        call assert_char('classification_2D',     stage_job_name(STAGE_POOL2D),    'pool 2D''s folder')
        call assert_int(7, NSTAGES, 'seven stages')
    end subroutine test_stage_names

    !> a stop of one stage, a restart of another, and updates
    subroutine test_gui_commands()
        type(stream_master_gui_commands) :: commands
        integer, allocatable      :: selection(:)
        type(string)              :: fname
        integer                   :: snapshot_id, iteration
        write(*,'(A)') 'test_gui_commands'
        call assert_true(commands%parse('{"terminate":false,"terminate_pool2D":true,"restart_preprocess":true,'//&
            &'"ctfresthreshold":5.5,"mskdiam2D":180.0,"snapshot2D":{"id":3,"iteration":12,"selection":[1,4],'//&
            &'"filename":"snap.simple"}}'), 'the answer is parsed')
        call assert_false(commands%l_terminate_all,                'the stream is not stopped')
        call assert_true(commands%l_terminate(STAGE_POOL2D),       'pool 2D is stopped')
        call assert_int(1, count(commands%l_terminate),            'and no other stage')
        call assert_true(commands%l_restart(STAGE_PREPROCESS),     'preprocessing is restarted')
        call assert_int(1, count(commands%l_restart),              'and no other stage')
        call assert_true(commands%update%assigned(),               'there are updates')
        call assert_real(5.5, commands%update%get_ctfres_update(),    1.e-5, 'the CTF resolution threshold')
        call assert_real(180., commands%update%get_mskdiam2D_update(), 1.e-5, 'the 2D mask diameter')
        call assert_true(commands%update%has_snapshot2D_update(),  'a snapshot request')
        call commands%update%get_snapshot2D_update(snapshot_id, iteration, selection, fname)
        call assert_int(3,  snapshot_id, 'its id')
        call assert_int(12, iteration,   'its iteration')
        call assert_true(all(selection == [1, 4]),        'its classes')
        call assert_char('snap.simple', fname%to_char(),  'its file')
        call commands%kill()
    end subroutine test_gui_commands

    !> an answer carries only its own updates: what an earlier one brought is not sent again
    subroutine test_gui_commands_fresh_each_answer()
        type(stream_master_gui_commands) :: commands
        write(*,'(A)') 'test_gui_commands_fresh_each_answer'
        call assert_true(commands%parse('{"mskdiam2D":180.0,"snapshot2D":{"id":3,"iteration":12,"selection":[1],'//&
            &'"filename":"snap.simple"}}'), 'the first answer is parsed')
        call assert_true(commands%parse('{"terminate":true}'), 'the second answer is parsed')
        call assert_true(commands%l_terminate_all,                 'it stops the stream')
        call assert_false(commands%update%assigned(),              'it has no updates')
        call assert_false(commands%update%has_snapshot2D_update(), 'the earlier snapshot request is gone')
        call assert_real(0., commands%update%get_mskdiam2D_update(), 1.e-6, 'and the earlier mask diameter')
        call commands%kill()
    end subroutine test_gui_commands_fresh_each_answer

    subroutine test_gui_commands_invalid()
        type(stream_master_gui_commands) :: commands
        write(*,'(A)') 'test_gui_commands_invalid'
        call assert_false(commands%parse('not json'),  'an answer that is not JSON is refused')
        call assert_false(commands%l_terminate_all,     'and asks nothing')
        call assert_false(any(commands%l_restart),      'no restart')
        call commands%kill()
    end subroutine test_gui_commands_invalid

    !> selections larger than the update holds are dropped whole, and the master does not stop;
    !! the rest of the answer is applied
    subroutine test_gui_commands_oversized()
        type(stream_master_gui_commands) :: commands
        write(*,'(A)') 'test_gui_commands_oversized'
        call assert_true(commands%parse('{"mskdiam2D":180.0,"pickrefs_cycle":1,"pickrefs_selection":'//&
            &int_list(MAX_PICKREFS_SELECTION + 1)//',"snapshot2D":{"id":3,"iteration":12,"selection":'//&
            &int_list(MAX_SNAPSHOT2D_SELECTION + 1)//',"filename":"snap.simple"}}'), 'the answer is parsed')
        call assert_int(0, commands%update%get_pickrefs_selection_length(), 'the oversized reference selection is dropped')
        call assert_false(commands%update%has_snapshot2D_update(),          'and so is the oversized snapshot')
        call assert_real(180., commands%update%get_mskdiam2D_update(), 1.e-5, 'the other updates are kept')
        call assert_true(commands%parse('{"pickrefs_cycle":1,"pickrefs_selection":'//&
            &int_list(MAX_PICKREFS_SELECTION)//'}'), 'a full selection is parsed')
        call assert_int(MAX_PICKREFS_SELECTION, commands%update%get_pickrefs_selection_length(), 'and kept')
        call commands%kill()
    end subroutine test_gui_commands_oversized

    !> a snapshot whose name is not a bare *.simple file name is dropped (p06 makes a folder of
    !! it), and the rest of the answer is applied
    subroutine test_gui_commands_snapshot_name()
        type(stream_master_gui_commands) :: commands
        character(len=:), allocatable    :: long_name
        write(*,'(A)') 'test_gui_commands_snapshot_name'
        call assert_false(snapshot_kept(commands, '../snap.simple'), 'a name with a directory is dropped')
        call assert_real(180., commands%update%get_mskdiam2D_update(), 1.e-5, 'the other updates are kept')
        call assert_false(snapshot_kept(commands, 'snap.txt'),       'a name without .simple is dropped')
        call assert_false(snapshot_kept(commands, '.simple'),        'a name of .simple alone is dropped')
        long_name = repeat('s', MAX_SNAPSHOT2D_FNAME_LEN)//'.simple'
        call assert_false(snapshot_kept(commands, long_name),        'an overlong name is dropped')
        long_name = repeat('s', MAX_SNAPSHOT2D_FNAME_LEN - len('.simple'))//'.simple'
        call assert_true(snapshot_kept(commands, long_name),         'a name of the full length is kept')
        call assert_true(snapshot_kept(commands, 'snapshot_3.simple'), 'and the names NICE sends')
        call commands%kill()

    contains

        logical function snapshot_kept( cmds, fname )
            type(stream_master_gui_commands), intent(inout) :: cmds
            character(len=*),                 intent(in)    :: fname
            call assert_true(cmds%parse('{"mskdiam2D":180.0,"snapshot2D":{"id":3,"iteration":12,'//&
                &'"selection":[1,4],"filename":"'//fname//'"}}'), 'the answer is parsed')
            snapshot_kept = cmds%update%has_snapshot2D_update()
        end function snapshot_kept

    end subroutine test_gui_commands_snapshot_name

    !> a stage's status replaces the previous one
    subroutine test_store_status()
        type(stream_master_meta_store)          :: store
        type(gui_metadata_stream_pool2D) :: status
        character(len=:), allocatable    :: buffer
        type(string)                     :: stage
        integer :: iter, nimported, naccepted, nrejected, tlast, mskdiam
        real    :: mskscale, res
        logical :: l_user_input, l_assigned
        write(*,'(A)') 'test_store_status'
        call store%new()
        call status%new(GUI_METADATA_STREAM_POOL2D_TYPE)
        call status%set(stage=string('first'), iteration=1, particles_imported=10, particles_accepted=8,&
            &particles_rejected=2, mskdiam=150, mskscale=1., resolution=9.)
        call status%serialise(buffer)
        call store%store(buffer)
        call status%set(stage=string('second'), iteration=2, particles_imported=20, particles_accepted=18,&
            &particles_rejected=2, mskdiam=150, mskscale=1., resolution=8.)
        call status%serialise(buffer)
        call store%store(buffer)
        l_assigned = store%pool2D%get(stage, iter, nimported, naccepted, nrejected, tlast, l_user_input, mskdiam,&
            &mskscale, res)
        call assert_char('second', stage%to_char(), 'the latest status is kept')
        call assert_int(2,  iter,      'with its iteration')
        call assert_int(20, nimported, 'and its counts')
        ! a frame with the right tag and the wrong length (a desynchronised pipe) is dropped
        call status%set(stage=string('third'), iteration=3, particles_imported=30, particles_accepted=28,&
            &particles_rejected=2, mskdiam=150, mskscale=1., resolution=7.)
        call status%serialise(buffer)
        call store%store(buffer(1:len(buffer)-4))
        l_assigned = store%pool2D%get(stage, iter, nimported, naccepted, nrejected, tlast, l_user_input, mskdiam,&
            &mskscale, res)
        call assert_char('second', stage%to_char(), 'a frame shorter than its type is dropped')
        call status%kill()
        call store%kill()
    end subroutine test_store_status

    !> list items go to slot i of a list of i_max, remade when i_max changes; an item out of range
    !! is dropped
    subroutine test_store_list()
        type(stream_master_meta_store) :: store
        write(*,'(A)') 'test_store_list'
        call store%new()
        call store%store(cavg_message(2, 3, 7))
        call assert_true(allocated(store%pool2D_cavgs), 'the list is made on the first item')
        if( allocated(store%pool2D_cavgs) )then
            call assert_int(3, size(store%pool2D_cavgs), 'with i_max slots')
            call assert_int(7, cavg_idx(store%pool2D_cavgs(2)), 'the item is in slot i')
        endif
        call store%store(cavg_message(4, 3, 9))
        call assert_int(3, size(store%pool2D_cavgs), 'an item out of range is dropped')
        call store%store(cavg_message(5, 5, 11))
        call assert_int(5, size(store%pool2D_cavgs), 'a new i_max remakes the list')
        call assert_int(11, cavg_idx(store%pool2D_cavgs(5)), 'with the new item in place')
        call store%kill()
        call assert_false(allocated(store%pool2D_cavgs), 'kill releases the lists')
    end subroutine test_store_list

    !> a volume sent with reprojection tiles arrives without them (they are not serialised)
    subroutine test_store_volume()
        type(stream_master_meta_store)       :: store
        type(gui_metadata_vol3D)      :: vol
        type(gui_metadata_cavg2D)     :: tiles(3)
        character(len=:), allocatable :: buffer
        integer :: i
        write(*,'(A)') 'test_store_volume'
        call store%new()
        call vol%new(GUI_METADATA_VOL3D_TYPE)
        call vol%set(string('reprojs.jpg'), string('vol.mrc'), string(''), string(''), string(''), 2, 64, 1.5, 2, 3)
        do i = 1,3
            call tiles(i)%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE)
        enddo
        call vol%set_reprojtiles(tiles)
        call vol%serialise(buffer)
        call store%store(buffer)
        call assert_true(allocated(store%solve3D_vols), 'the volume list is made')
        if( allocated(store%solve3D_vols) )then
            call assert_int(3, size(store%solve3D_vols),           'one slot per state')
            call assert_int(2, store%solve3D_vols(2)%get_state(), 'the volume is in its slot')
        endif
        call vol%kill()
        call store%kill()
    end subroutine test_store_volume

    !> a restarted stage's lists are dropped, and only that stage's
    subroutine test_store_clear_stage()
        type(stream_master_meta_store) :: store
        type(gui_metadata_vol3D)       :: vol
        character(len=:), allocatable  :: buffer
        write(*,'(A)') 'test_store_clear_stage'
        call store%new()
        call vol%new(GUI_METADATA_VOL3D_TYPE)
        call vol%set(string('reprojs.jpg'), string('vol.mrc'), string(''), string(''), string(''), 1, 64, 1.5, 1, 2)
        call vol%serialise(buffer)
        call store%store(buffer)
        call store%store(cavg_message(1, 2, 7))
        call assert_true(allocated(store%solve3D_vols), 'a volume is stored')
        call store%clear_stage(STAGE_SOLVE3D)
        call assert_false(allocated(store%solve3D_vols), 'restarting multistate 3D drops its volumes')
        call assert_true(allocated(store%pool2D_cavgs),  'and leaves the pool''s class averages')
        call store%clear_stage(STAGE_POOL2D)
        call assert_false(allocated(store%pool2D_cavgs), 'restarting pool 2D drops them')
        call vol%kill()
        call store%kill()
    end subroutine test_store_clear_stage

    !> a stage keeps the commander it is given; its messages reach the master's reader, the master's
    !! updates reach the stage, and a restart's discard empties both pipes
    subroutine test_stage_pipes()
        class(stream_master_stage), allocatable :: proc
        type(stream_pipe)             :: stage_side
        class(cmdline), allocatable :: cline
        type(noop_commander)          :: commander
        character(len=:), allocatable :: buffer
        integer :: fd_write, fd_read, fd_master_read, fd_master_write
        write(*,'(A)') 'test_stage_pipes'
        allocate(proc)
        allocate(cline)
        call proc%new(STAGE_POOL2D, commander, cline, .true., max_metadata_size())
        call assert_true(allocated(proc%fork%commander), 'the stage keeps its commander')
        call stage_fds(STAGE_POOL2D, fd_write, fd_read)
        call assert_true(fd_write >= 0 .and. fd_read >= 0, 'the stage''s pipes are open')
        call stage_side%new(fd_read, fd_write, max_metadata_size(), 'stage side')
        call stage_side%send('status')
        call assert_true(proc%reader%receive(buffer), 'the master reads what the stage sends')
        if( allocated(buffer) ) call assert_char('status', buffer, 'intact')
        call proc%writer%send('update')
        call assert_true(stage_side%receive(buffer), 'the stage reads what the master sends')
        if( allocated(buffer) ) call assert_char('update', buffer, 'intact')
        call stage_side%send('stale message')
        call proc%writer%send('stale update')
        call proc%discard_pipes(max_metadata_size())
        call assert_false(proc%reader%receive(buffer), 'after the discard the master reads nothing')
        call assert_false(stage_side%receive(buffer),  'nor does the next process of the stage')
        call stage_side%kill()
        call proc%kill()
        call master_fds(STAGE_POOL2D, fd_master_read, fd_master_write)
        call assert_int(-1, fd_master_read,  'kill closes the pipes')
        call assert_int(-1, fd_master_write, 'both of them')
        call assert_false(allocated(proc%fork%commander), 'and releases the commander')
        call cline%kill()
    end subroutine test_stage_pipes

    !> an update goes to a stage only when it differs from the last one the stage was sent; a stop
    !! request is remembered (no more updates until the next start), and kill forgets both
    subroutine test_update_dedupe()
        class(stream_master_stage), allocatable :: proc
        class(cmdline), allocatable :: cline
        type(noop_commander)                    :: commander
        write(*,'(A)') 'test_update_dedupe'
        allocate(proc)
        allocate(cline)
        call proc%new(STAGE_POOL2D, commander, cline, .true., max_metadata_size())
        call assert_true(proc%is_new_update('thresholds'), 'the first update is new')
        proc%last_update = 'thresholds'
        call assert_false(proc%is_new_update('thresholds'), 'a repeated update is not')
        call assert_true(proc%is_new_update('threshold2'), 'a changed update is')
        call assert_true(proc%is_new_update('thresholds '), 'so is a longer one, whatever it ends with')
        call assert_false(proc%l_stop_requested, 'a new stage has no stop request')
        call proc%request_stop()
        call assert_true(proc%l_stop_requested, 'a stop request is remembered')
        call proc%kill()
        call assert_false(allocated(proc%last_update), 'kill forgets the last update')
        call assert_false(proc%l_stop_requested, 'and the stop request')
        deallocate(proc)
        call cline%kill()
    end subroutine test_update_dedupe

    !> a skipped stage is never started: start() leaves it skipped and forks nothing
    subroutine test_skipped_stays_skipped()
        class(stream_master_stage), allocatable :: proc
        class(cmdline), allocatable :: cline
        type(noop_commander)                    :: commander
        write(*,'(A)') 'test_skipped_stays_skipped'
        allocate(proc)
        allocate(cline)
        call proc%new(STAGE_PREPROCESS, commander, cline, .true., max_metadata_size())
        call assert_false(proc%is_skipped(), 'a new stage is not skipped')
        call proc%skip()
        call assert_true(proc%is_skipped(),  'a skipped stage is skipped')
        call assert_false(proc%is_running(), 'and not running')
        call proc%start()
        call assert_true(proc%is_skipped(),  'starting a skipped stage leaves it skipped')
        call assert_false(proc%is_running(), 'and forks nothing')
        call proc%kill()
        deallocate(proc)
        call cline%kill()
    end subroutine test_skipped_stays_skipped

    ! ---- the GUI heartbeat over live processes ---------------------------------

    !> Reads live forked processes into the GUI heartbeat's records. They fork real children, so
    !! they run in the forked_process platform entry, not with the other master tests.
    subroutine run_stream_heartbeat_tests()
        write(*,'(A)') '**** running all stream heartbeat tests ****'
        call test_fork_gui_status()
    end subroutine run_stream_heartbeat_tests

    !> a running child reports running, with its pid and start time and no stop time; once
    !! stopped, finished, with its stop time
    subroutine test_fork_gui_status()
        type(forked_process)   :: fork
        type(gui_stage_status) :: stage_status
        integer                :: rc
        write(*,'(A)') 'test_fork_gui_status'
#if defined(_WIN32)
        write(*,'(A)') 'skipped: forked processes are unavailable on Windows'
#else
        call fork%start(name=string('TEST_HEARTBEAT_STAGE'))
        rc = c_usleep(FORK_POLL_TIME * 5)
        stage_status = fork_gui_status(fork)
        call assert_int(GUI_STAGE_STATUS_RUNNING, stage_status%status, 'a running child reports running')
        call assert_true(stage_status%pid > 0,                         'with its pid')
        call assert_true(stage_status%starttime > 0,                   'and its start time')
        call assert_int(0, stage_status%stoptime,                      'and no stop time')
        call fork%terminate()
        call fork%await_final_status()
        stage_status = fork_gui_status(fork)
        call assert_int(GUI_STAGE_STATUS_FINISHED, stage_status%status, 'a stopped child reports finished')
        call assert_true(stage_status%stoptime > 0,                     'with its stop time')
#endif
    end subroutine test_fork_gui_status

    ! ---- fixtures ------------------------------------------------------------

    subroutine execute_noop( self, cline )
        class(noop_commander), intent(inout) :: self
        class(cmdline),        intent(inout) :: cline
    end subroutine execute_noop

    ! a serialised pool 2D class average, item i of i_max, of class idx
    function cavg_message( i, i_max, idx ) result( buffer )
        integer, intent(in) :: i, i_max, idx
        character(len=:), allocatable :: buffer
        type(gui_metadata_cavg2D) :: cavg
        call cavg%new(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
        call cavg%set(path=string('cavgs.jpg'), mrcpath=string('cavgs.mrc'), idx=idx,&
            &sprite=sprite_sheet_pos(x=0., y=0., h=100, w=100), i=i, i_max=i_max)
        call cavg%serialise(buffer)
        call cavg%kill()
    end function cavg_message

    integer function cavg_idx( cavg )
        type(gui_metadata_cavg2D), intent(in) :: cavg
        type(string)           :: path, mrcpath
        type(sprite_sheet_pos) :: sprite
        logical                :: l_assigned
        l_assigned = cavg%get(path, mrcpath, cavg_idx, sprite)
    end function cavg_idx

    ! the JSON array [1,2,...,n]
    function int_list( n ) result( list )
        integer, intent(in) :: n
        character(len=:), allocatable :: list
        integer :: i
        list = '['
        do i = 1,n
            list = list//int2str(i)
            if( i < n ) list = list//','
        enddo
        list = list//']'
    end function int_list

end module simple_stream_master_tester
