!@descr: unit tests for the history of the stream watcher (simple_stream_watcher)
! The history holds the basenames of the files a watcher has reported; is_past finds them by binary
! search in a sorted index. Each test works in a fresh fixture directory, with files touched there
! (add2history only takes files that exist).
module simple_stream_watcher_tester
use simple_test_utils
use simple_string,         only: string
use simple_string_utils,   only: int2str, int2str_pad
use simple_fileio,         only: simple_getcwd, simple_touch
use simple_stream_watcher, only: stream_watcher
implicit none
private
public :: run_all_stream_watcher_tests

contains

    subroutine run_all_stream_watcher_tests()
        write(*,'(A)') '**** running all stream watcher tests ****'
        call test_history_lookup()
        call test_history_growth()
        call test_rate_after_restored_history()
    end subroutine run_all_stream_watcher_tests

    !> the history a restart restores before the first watch is not counted as movies detected
    !! since: the rate counts only what is added after the first watch
    subroutine test_rate_after_restored_history()
        integer, parameter   :: NRESTORED = 100
        type(stream_watcher) :: watcher
        type(string), allocatable :: movies(:)
        type(string)         :: cwd_saved, root, cwd
        integer              :: nfail0, i, nmovies
        write(*,'(A)') 'test_rate_after_restored_history'
        nfail0 = tests_failed
        call enter_fixture('watcher_rate', cwd_saved, root)
        call simple_getcwd(cwd)
        watcher = stream_watcher(-1, cwd)
        do i = 1,NRESTORED
            call simple_touch(string('restored_'//int2str_pad(i, 3)//'.mrc'))
            call watcher%add2history(cwd//'/restored_'//int2str_pad(i, 3)//'.mrc')
        enddo
        call watcher%watch(nmovies, movies)
        call assert_int(0, nmovies, 'the restored movies are not new')
        do i = 1,2
            call simple_touch(string('new_'//int2str(i)//'.mrc'))
            call watcher%add2history(cwd//'/new_'//int2str(i)//'.mrc')
        enddo
        call sleep(1)
        call watcher%watch(nmovies, movies)
        call assert_true(watcher%rate <= 7200, 'two movies a second or less, not the restored hundred')
        call assert_true(watcher%rate > 0,     'the new ones count')
        call watcher%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_rate_after_restored_history

    !> files added in any order are found, by basename, whatever their directory; a file not added
    !! is not; a file added twice is kept once; clearing empties the history
    subroutine test_history_lookup()
        character(len=*), parameter :: NAMES(5) = [character(len=12) :: 'movie_c.mrc', 'movie_a.mrc',&
            &'movie_e.mrc', 'movie_b.mrc', 'movie_d.mrc']
        type(stream_watcher) :: watcher
        type(string)         :: cwd_saved, root, cwd
        integer              :: nfail0, i
        write(*,'(A)') 'test_history_lookup'
        nfail0 = tests_failed
        call enter_fixture('watcher_history', cwd_saved, root)
        call simple_getcwd(cwd)
        watcher = stream_watcher(-1, cwd)
        do i = 1,size(NAMES)
            call simple_touch(string(trim(NAMES(i))))
            call watcher%add2history(cwd//'/'//trim(NAMES(i)))
        enddo
        call assert_int(5, watcher%n_history, 'five files in the history')
        do i = 1,size(NAMES)
            call assert_true(watcher%is_past(string(trim(NAMES(i)))), trim(NAMES(i))//' is found')
        enddo
        call assert_true(watcher%is_past(string('/elsewhere/movie_d.mrc')), 'by basename, whatever the directory')
        call assert_false(watcher%is_past(string('movie_f.mrc')),          'a file never added is not found')
        call assert_false(watcher%is_past(string('movie_.mrc')),           'nor a name sorting before all')
        call assert_false(watcher%is_past(string('movie_z.mrc')),          'nor one sorting after all')
        call watcher%add2history(cwd//'/movie_b.mrc')
        call assert_int(5, watcher%n_history, 'a file added twice is kept once')
        call watcher%clear_history()
        call assert_int(0, watcher%n_history,                      'clearing empties the history')
        call assert_false(watcher%is_past(string('movie_a.mrc')), 'and nothing is found')
        call watcher%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_history_lookup

    !> a history past its first capacity grows and keeps every file findable
    subroutine test_history_growth()
        integer, parameter   :: NFILES = 1500
        type(stream_watcher) :: watcher
        type(string)         :: cwd_saved, root, cwd, fname
        integer              :: nfail0, i, j, nfound
        write(*,'(A)') 'test_history_growth'
        nfail0 = tests_failed
        call enter_fixture('watcher_growth', cwd_saved, root)
        call simple_getcwd(cwd)
        watcher = stream_watcher(-1, cwd)
        ! added in an order unrelated to their names' order
        do i = 1,NFILES
            j     = mod(i * 7919, NFILES) + 1
            fname = 'mic_'//int2str_pad(j, 5)//'.mrc'
            call simple_touch(fname)
            call watcher%add2history(cwd//'/'//fname%to_char())
        enddo
        call assert_int(NFILES, watcher%n_history, int2str(NFILES)//' files in the history')
        nfound = 0
        do i = 1,NFILES
            fname = 'mic_'//int2str_pad(i, 5)//'.mrc'
            if( watcher%is_past(fname) ) nfound = nfound + 1
        enddo
        call assert_int(NFILES, nfound, 'every one is found')
        call assert_false(watcher%is_past(string('mic_'//int2str_pad(NFILES + 1, 5)//'.mrc')), 'and no other')
        call watcher%kill
        call leave_fixture(cwd_saved, root, nfail0)
    end subroutine test_history_growth

end module simple_stream_watcher_tester
