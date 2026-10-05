!@descr: unit tests for the numbered job sets of the stream stages (simple_stream_job_sets)
! Naming, writing, completing and restoring sets, in a fresh fixture directory each. Queueing,
! scheduling and collecting need a queue and are covered by the stages' high-level tests.
module simple_stream_job_sets_tester
use simple_test_utils
use simple_defs_fname,      only: STDERROUT_DIR
use simple_string,          only: string
use simple_fileio,          only: file_exists, simple_getcwd
use simple_syslib,          only: dir_exists
use simple_cmdline,         only: cmdline
use simple_sp_project,      only: sp_project
use simple_stream_job_sets, only: stream_job_sets
implicit none
private
public :: run_all_stream_job_sets_tests

character(len=*), parameter :: JOB_FOLDER       = 'spprojs'
character(len=*), parameter :: COMPLETED_FOLDER = 'spprojs_completed'

contains

    subroutine run_all_stream_job_sets_tests()
        write(*,'(A)') '**** running all stream job set tests ****'
        call test_new_makes_folders()
        call test_write_set()
        call test_complete()
        call test_restore()
        call test_restore_past_unfinished()
    end subroutine run_all_stream_job_sets_tests

    subroutine test_new_makes_folders()
        type(stream_job_sets) :: sets
        type(string)          :: cwd_saved, root, cwd, dir
        integer               :: nfail0
        write(*,'(A)') 'test_new_makes_folders'
        nfail0 = tests_failed
        call enter_fixture('job_sets_new', cwd_saved, root)
        call simple_getcwd(cwd)
        call sets%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        call assert_true(dir_exists(string(JOB_FOLDER)),                     'the job folder is made')
        call assert_true(dir_exists(string(JOB_FOLDER//'/'//STDERROUT_DIR)), 'with the folder for the job output')
        call assert_true(dir_exists(string(COMPLETED_FOLDER)),               'the completed folder is made')
        dir = sets%get_job_dir()
        call assert_char(cwd%to_char()//'/'//JOB_FOLDER, dir%to_char(), 'the job folder is kept as an absolute path')
        dir = sets%get_completed_dir()
        call assert_char(cwd%to_char()//'/'//COMPLETED_FOLDER, dir%to_char(), 'so is the completed folder')
        call assert_int(0, sets%get_counter(), 'no set yet')
        call sets%kill
        call sets%kill ! idempotence
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_new_makes_folders

    !> sets are numbered and named in turn, written in the job folder with their project info and
    !! origin, and the worker command line points at the latest one
    subroutine test_write_set()
        type(stream_job_sets) :: sets
        type(sp_project)      :: set_proj, written
        type(cmdline)         :: cline_worker
        type(string)          :: cwd_saved, root, val, job_dir
        integer               :: nfail0
        write(*,'(A)') 'test_write_set'
        nfail0 = tests_failed
        call enter_fixture('job_sets_write', cwd_saved, root)
        call sets%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        job_dir = sets%get_job_dir()
        call set_proj%os_mic%new(3, is_ptcl=.false.)
        call sets%write_set(set_proj, cline_worker, 3, source=string('/upstream/00042.simple'))
        call assert_int(1, sets%get_counter(), 'the first set is number 1')
        call assert_true(file_exists(job_dir//'/00001.simple'), 'it is written in the job folder as 00001.simple')
        call sets%write_set(set_proj, cline_worker, 2)
        call assert_int(2, sets%get_counter(), 'the second set is number 2')
        call assert_true(file_exists(job_dir//'/00002.simple'), 'and written as 00002.simple')
        val = cline_worker%get_carg('projfile')
        call assert_char('00002.simple', val%to_char(), 'the worker command line names the latest set')
        val = cline_worker%get_carg('projname')
        call assert_char('00002', val%to_char(), 'and its project name')
        call assert_int(1, cline_worker%get_iarg('fromp'), 'from its first item')
        call assert_int(2, cline_worker%get_iarg('top'),   'to its last')
        call written%read_segment('projinfo', job_dir//'/00001.simple')
        ! the project file name, not projname: the orientation reader stores a numeric-looking value
        ! such as the name 00001 as a real, which get_str does not return
        val = written%projinfo%get_str(1, 'projfile')
        call assert_char('00001.simple', val%to_char(), 'the set records its project file')
        val = written%projinfo%get_str(1, 'cwd')
        call assert_char(job_dir%to_char(), val%to_char(), 'and the job folder as its directory')
        call written%read_segment('mic', job_dir//'/00001.simple')
        call assert_int(3, written%os_mic%get_noris(), 'its micrographs are written')
        call written%kill
        call set_proj%kill
        call cline_worker%kill
        call sets%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_write_set

    subroutine test_complete()
        type(stream_job_sets) :: sets
        type(sp_project)      :: set_proj
        type(cmdline)         :: cline_worker
        type(string)          :: cwd_saved, root, job_dir, completed_dir, completed
        integer               :: nfail0
        write(*,'(A)') 'test_complete'
        nfail0 = tests_failed
        call enter_fixture('job_sets_complete', cwd_saved, root)
        call sets%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        job_dir       = sets%get_job_dir()
        completed_dir = sets%get_completed_dir()
        call sets%write_set(set_proj, cline_worker, 1)
        call sets%complete(job_dir//'/00001.simple', completed)
        call assert_char(completed_dir%to_char()//'/00001.simple', completed%to_char(), 'the new path is in the completed folder')
        call assert_true(file_exists(completed),                     'the set is there')
        call assert_false(file_exists(job_dir//'/00001.simple'),     'and no longer in the job folder')
        call set_proj%kill
        call cline_worker%kill
        call sets%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_complete

    !> restart: numbering continues after the highest completed set; each set's origin comes back
    !! ('' when none was recorded); the unfinished sets are dropped, and the job folder is ready
    subroutine test_restore()
        type(stream_job_sets)     :: sets, restarted
        type(sp_project)          :: set_proj
        type(cmdline)             :: cline_worker
        type(string), allocatable :: completed(:), sources(:)
        type(string)              :: cwd_saved, root, job_dir, done
        integer                   :: nfail0, i
        logical                   :: l_found
        write(*,'(A)') 'test_restore'
        nfail0 = tests_failed
        call enter_fixture('job_sets_restore', cwd_saved, root)
        call sets%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        job_dir = sets%get_job_dir()
        ! sets 1 to 7; 3 (with an origin) and 7 (without) complete, 4 left unfinished
        do i = 1,7
            if( i == 3 )then
                call sets%write_set(set_proj, cline_worker, 1, source=string('/upstream/00011.simple'))
            else
                call sets%write_set(set_proj, cline_worker, 1)
            endif
        enddo
        call sets%complete(job_dir//'/00003.simple', done)
        call sets%complete(job_dir//'/00007.simple', done)
        call sets%kill
        call restarted%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        call restarted%restore(completed, sources)
        call assert_int(2, size(completed), 'the two completed sets come back')
        call assert_int(7, restarted%get_counter(), 'numbering continues after the highest completed set')
        l_found = .false.
        do i = 1,size(completed)
            if( completed(i)%has_substr('00003.simple') )then
                call assert_char('/upstream/00011.simple', sources(i)%to_char(), 'set 3 returns its origin')
                l_found = .true.
            else
                call assert_int(0, sources(i)%strlen(), 'set 7 recorded no origin')
            endif
        enddo
        call assert_true(l_found, 'set 3 is among the completed sets')
        call assert_false(file_exists(job_dir//'/00004.simple'), 'the unfinished sets are dropped')
        call assert_true(file_exists(job_dir//'_unfinished1/00004.simple'), 'with their folder set aside, which a job may still use')
        call assert_true(dir_exists(job_dir//'/'//STDERROUT_DIR), 'the job folder is ready for new sets')
        call restarted%write_set(set_proj, cline_worker, 1)
        call assert_true(file_exists(job_dir//'/00008.simple'), 'the next set is number 8')
        call set_proj%kill
        call cline_worker%kill
        call restarted%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restore

    !> restart: numbering continues past an unfinished set numbered above the highest completed one
    subroutine test_restore_past_unfinished()
        type(stream_job_sets)     :: sets, restarted
        type(sp_project)          :: set_proj
        type(cmdline)             :: cline_worker
        type(string), allocatable :: completed(:)
        type(string)              :: cwd_saved, root, job_dir, done
        integer                   :: nfail0, i
        write(*,'(A)') 'test_restore_past_unfinished'
        nfail0 = tests_failed
        call enter_fixture('job_sets_restore_unfinished', cwd_saved, root)
        call sets%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        job_dir = sets%get_job_dir()
        ! sets 1 to 5; 3 complete, 5 still queued or running when the stage stopped
        do i = 1,5
            call sets%write_set(set_proj, cline_worker, 1)
        enddo
        call sets%complete(job_dir//'/00003.simple', done)
        call sets%kill
        call restarted%new(string(JOB_FOLDER), string(COMPLETED_FOLDER), 5)
        call restarted%restore(completed)
        call assert_int(1, size(completed), 'the completed set comes back')
        call assert_int(5, restarted%get_counter(), 'numbering continues after the highest set, completed or not')
        call restarted%write_set(set_proj, cline_worker, 1)
        call assert_true(file_exists(job_dir//'/00006.simple'), 'so a new set never takes an unfinished set''s number')
        call set_proj%kill
        call cline_worker%kill
        call restarted%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_restore_past_unfinished

end module simple_stream_job_sets_tester
