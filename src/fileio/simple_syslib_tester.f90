module simple_syslib_tester
use simple_syslib
use simple_string
use simple_defs
use simple_test_utils
implicit none
private
public :: run_all_syslib_tests
type(string) :: test_root    ! base test directory
type(string) :: original_cwd ! cwd to restore

contains

    subroutine run_all_syslib_tests()
        write(*,'(A)') '**** running all syslib tests ****'
        call setup_test_env()
        call test_simple_touch_and_file_exists()
        call test_simple_mkdir_dir_exists_rmdir()
        call test_simple_rename_and_del_file()
        call test_syslib_symlink()
        call test_simple_file_stat()
        call test_is_io_and_is_open()
        call test_is_file_open()
        call test_simple_getcwd_and_chdir()
        call test_simple_list_dirs()
        call test_simple_list_files()
        call test_simple_list_files_regexp()
        call test_simple_abspath()
        call test_getenv_wrapper()
        call test_exec_cmdline()
        call test_get_process_id()
        call test_syslib_c2fortran_string()
        call test_find_next_int_dir_prefix()
        call test_print_slurm_env_smoke()
        call teardown_test_env()
        ! call report_summary()
    end subroutine run_all_syslib_tests

    !---------------- test environment helpers ----------------

    subroutine setup_test_env()
        type(string) :: cwd
        integer      :: pid, istat
        write(*,'(A)') 'setup_test_env'
        call simple_getcwd(cwd)
        original_cwd = cwd
        pid       = get_process_id()
        test_root = cwd // '/' // 'simple_syslib_test_' // string(pid)
        if (.not. dir_exists(test_root)) then
            call simple_mkdir(test_root, verbose=.false.)
            call assert_true(dir_exists(test_root), 'setup: test_root directory created')
        end if
        call simple_chdir(test_root, status=istat)
        call assert_int(0, istat, 'setup: chdir to test_root')
    end subroutine setup_test_env

    subroutine teardown_test_env()
        integer :: istat
        write(*,'(A)') 'teardown_test_env'
        ! go back to original CWD
        call simple_chdir(original_cwd, status=istat)
        call assert_int(0, istat, 'teardown: chdir back to original CWD')
        ! remove test_root directory and its contents as best we can
        ! (tests should have cleaned up most things already)
        call simple_rmdir(test_root, status=istat)
    end subroutine teardown_test_env

    !---------------- simple_touch / file_exists / del_file ----------------

    subroutine test_simple_touch_and_file_exists()
        type(string) :: fname
        integer :: istat
        write(*,'(A)') 'test_simple_touch_and_file_exists'
        fname = 'touch_test.txt'
        call assert_true(.not. file_exists(fname), 'touch: file does not exist at start')
        call simple_touch(fname, status=istat)
        call assert_int(0, istat, 'simple_touch status')
        call assert_true(file_exists(fname), 'touch: file exists after simple_touch')
        call del_file(fname)
        call assert_true(.not. file_exists(fname), 'del_file removes file')
    end subroutine test_simple_touch_and_file_exists

    !---------------- simple_mkdir / dir_exists / simple_rmdir ----------------

    subroutine test_simple_mkdir_dir_exists_rmdir()
        type(string) :: dname
        integer :: istat
        write(*,'(A)') 'test_simple_mkdir_dir_exists_rmdir'
        dname = 'test_dir'
        call assert_true(.not. dir_exists(dname), 'mkdir: dir does not exist at start')
        call simple_mkdir(dname, verbose=.false.)
        call assert_true(dir_exists(dname), 'mkdir: directory exists after simple_mkdir')
        call simple_rmdir(dname, status=istat)
        call assert_int(0, istat, 'rmdir: status')
        call assert_true(.not. dir_exists(dname), 'rmdir: directory removed')
    end subroutine test_simple_mkdir_dir_exists_rmdir

    !---------------- simple_rename + del_file ----------------

    subroutine test_simple_rename_and_del_file()
        type(string) :: src, dst
        integer :: istat
        write(*,'(A)') 'test_simple_rename_and_del_file'
        src = 'rename_src.txt'
        dst = 'rename_dst.txt'
        call simple_touch(src, status=istat)
        call assert_int(0, istat, 'rename: touch src')
        call assert_true(file_exists(src), 'rename: src exists')
        call assert_true(.not. file_exists(dst), 'rename: dst does not exist initially')
        call simple_rename(src, dst, overwrite=.true.)
        call assert_true(.not. file_exists(src), 'rename: src gone after rename')
        call assert_true(file_exists(dst), 'rename: dst exists after rename')
        call del_file(dst)
        call assert_true(.not. file_exists(dst), 'rename: dst deleted with del_file')
    end subroutine test_simple_rename_and_del_file

    !---------------- syslib_symlink ----------------

    subroutine test_syslib_symlink()
        type(string) :: target, link
        integer      :: istat
        write(*,'(A)') 'test_syslib_symlink'
        target = 'symlink_target.txt'
        link   = 'symlink_link.txt'
        call simple_touch(target, status=istat)
        call assert_int(0, istat, 'symlink: touch target')
        call assert_true(file_exists(target), 'symlink: target exists')
        call syslib_symlink(target, link, status=istat)
        call assert_int(0, istat, 'symlink: syslib_symlink status')
        call assert_true(file_exists(link), 'symlink: link appears in filesystem')
        call del_file(link)
        call del_file(target)
    end subroutine test_syslib_symlink

    !---------------- simple_file_stat ----------------

    subroutine test_simple_file_stat()
        type(string) :: fname
        integer, allocatable :: buf(:)
        integer :: istat
        write(*,'(A)') 'test_simple_file_stat'
        fname = 'stat_test.txt'
        call simple_touch(fname, status=istat)
        call assert_int(0, istat, 'stat: touch file')
        call simple_file_stat(fname, istat, buf)
        call assert_int(0, istat, 'stat: simple_file_stat status')
        call assert_true(allocated(buf), 'stat: buffer allocated')
        call assert_int(13, size(buf), 'stat: buffer size 13')
        if (allocated(buf)) deallocate(buf)
        call del_file(fname)
    end subroutine test_simple_file_stat

    !---------------- is_io / is_open ----------------

    subroutine test_is_io_and_is_open()
        integer :: u, ios
        write(*,'(A)') 'test_is_io_and_is_open'
        call assert_true(is_io(ERROR_UNIT),  'is_io: ERROR_UNIT')
        call assert_true(is_io(OUTPUT_UNIT), 'is_io: OUTPUT_UNIT')
        call assert_true(is_io(INPUT_UNIT),  'is_io: INPUT_UNIT')
        ! scratch file
        open(newunit=u, status='scratch', iostat=ios)
        call assert_int(0, ios, 'is_open: open scratch')
        call assert_true(is_open(u), 'is_open: scratch file is open')
        close(u, iostat=ios)
        call assert_int(0, ios, 'is_open: close scratch')
        call assert_true(.not. is_open(u), 'is_open: closed unit reports not open')
    end subroutine test_is_io_and_is_open

    !---------------- is_file_open ----------------

    subroutine test_is_file_open()
        type(string) :: fname
        integer :: u, ios
        write(*,'(A)') 'test_is_file_open'
        fname = 'open_test.txt'
        open(newunit=u, file=fname%to_char(), status='replace', iostat=ios)
        call assert_int(0, ios, 'is_file_open: open file')
        call assert_true(is_file_open(fname), 'is_file_open: file is reported open')
        close(u, iostat=ios)
        call assert_int(0, ios, 'is_file_open: close file')
        call assert_true(.not. is_file_open(fname), 'is_file_open: file not open after close')
        call del_file(fname)
    end subroutine test_is_file_open

    !---------------- simple_getcwd / simple_chdir ----------------

    subroutine test_simple_getcwd_and_chdir()
        type(string) :: cwd_now, cwd_saved
        type(string) :: subdir
        integer      :: istat
        write(*,'(A)') 'test_simple_getcwd_and_chdir'
        call simple_getcwd(cwd_saved)
        call assert_true(.not. cwd_saved%is_blank(), 'getcwd: not blank')
        subdir = 'cwd_test_dir'
        call simple_mkdir(subdir, verbose=.false.)
        call simple_chdir(subdir, status=istat)
        call assert_int(0, istat, 'chdir: change into subdir')
        call simple_getcwd(cwd_now)
        call assert_true(cwd_now%has_substr(subdir), 'chdir: cwd contains subdir name')
        call simple_chdir(cwd_saved, status=istat)
        call assert_int(0, istat, 'chdir: back to saved cwd')
        call simple_rmdir(subdir, status=istat)
    end subroutine test_simple_getcwd_and_chdir

    !---------------- simple_list_dirs ----------------

    subroutine test_simple_list_dirs()
        type(string) :: d1, d2
        type(string), allocatable :: dirs(:)
        integer :: i, stat
        write(*,'(A)') 'test_simple_list_dirs'
        d1 = 'ldir_1'
        d2 = 'ldir_2'
        call simple_mkdir(d1, verbose=.false.)
        call simple_mkdir(d2, verbose=.false.)
        dirs = simple_list_dirs('.', status=stat)
        call assert_int(0, stat, 'list_dirs: status')
        ! just check that both names appear at least once
        call assert_true(any([ (dirs(i) == d1, i=1,size(dirs)) ]), 'list_dirs: contains d1')
        call assert_true(any([ (dirs(i) == d2, i=1,size(dirs)) ]), 'list_dirs: contains d2')
        if (allocated(dirs)) deallocate(dirs)
        call simple_rmdir(d1, status=stat)
        call simple_rmdir(d2, status=stat)
    end subroutine test_simple_list_dirs

    !---------------- simple_list_files ----------------

    subroutine test_simple_list_files()
        type(string), allocatable :: list(:)
        type(string) :: f1, f2
        integer      :: istat
        write(*,'(A)') 'test_simple_list_files'
        f1 = 'file_a.txt'
        f2 = 'file_b.txt'
        call simple_touch(f1, status=istat)
        call simple_touch(f2, status=istat)
        call simple_list_files('file_*.txt', list)
        call assert_true(allocated(list), 'list_files: list allocated')
        call assert_true(any(list == f1), 'list_files: contains f1')
        call assert_true(any(list == f2), 'list_files: contains f2')
        if (allocated(list)) deallocate(list)
        call del_file(f1)
        call del_file(f2)
    end subroutine test_simple_list_files

    !---------------- simple_list_files_regexp ----------------

    subroutine test_simple_list_files_regexp()
        type(string), allocatable :: list(:)
        type(string) :: d, f1, f2
        integer      :: istat
        write(*,'(A)') 'test_simple_list_files_regexp'
        d  = 'regdir'
        call simple_mkdir(d, verbose=.false.)
        f1 = d//'/'//'data_001.log'
        f2 = d//'/'//'note.txt'
        call simple_touch(f1, status=istat)
        call simple_touch(f2, status=istat)
        call simple_list_files_regexp(d, 'data_.*\.log', list)
        call assert_true(allocated(list), 'list_files_regexp: list allocated')
        call assert_true(any(list == f1), 'list_files_regexp: contains matching file')
        call assert_true(.not. any(list == f2), 'list_files_regexp: does not contain non-matching file')
        if (allocated(list)) deallocate(list)
        call del_file(f1)
        call del_file(f2)
        call simple_rmdir(d, status=istat)
    end subroutine test_simple_list_files_regexp

    !---------------- simple_abspath ----------------

    subroutine test_simple_abspath()
        type(string) :: fname, absname
        integer :: istat
        write(*,'(A)') 'test_simple_abspath'
        fname = 'abs_test.txt'
        call simple_touch(fname, status=istat)
        call assert_int(0, istat, 'abspath: touch file')
        absname = simple_abspath(fname, status=istat, check_exists=.true.)
        call assert_int(0, istat, 'abspath: status')
        ! crude check: absolute path should start with '/'
        call assert_true(absname%has_substr('/'), 'abspath: contains "/"')
        call del_file(fname)
    end subroutine test_simple_abspath

    !---------------- simple_getenv wrappers ----------------

    subroutine test_getenv_wrapper()
        integer      :: status
        type(string) :: s
        write(*,'(A)') 'test_getenv_wrappers'
        ! choose a variable that almost always exists, e.g. PATH
        s = simple_getenv('PATH', status)
        call assert_true(status == 0 .or. status == 1, 'getenv_1: status 0 or 1')
        if (status == 0) then
            call assert_true(s%strlen_trim() > 0, 'getenv_1: non-empty PATH')
        end if
    end subroutine test_getenv_wrapper

    !---------------- exec_cmdline ----------------

    subroutine test_exec_cmdline()
        integer :: exitstat
        write(*,'(A)') 'test_exec_cmdline'
        ! benign command that should succeed on most systems
        call exec_cmdline('true', waitflag=.true., suppress_errors=.false., exitstat=exitstat)
        call assert_int(0, exitstat, 'exec_cmdline: true command exit status')
        ! Another command that writes a file
        call exec_cmdline('echo hello > exec_test.txt', waitflag=.true., suppress_errors=.true., exitstat=exitstat)
        call assert_true(file_exists('exec_test.txt'), 'exec_cmdline: created file')
        call del_file('exec_test.txt')
    end subroutine test_exec_cmdline

    !---------------- get_process_id ----------------

    subroutine test_get_process_id()
        integer(4) :: pid
        write(*,'(A)') 'test_get_process_id'
        pid = get_process_id()
        call assert_true(pid > 0, 'get_process_id: positive pid')
    end subroutine test_get_process_id

    !---------------- syslib_c2fortran_string ----------------

    subroutine test_syslib_c2fortran_string()
        character(len=32) :: s
        integer           :: l
        write(*,'(A)') 'test_syslib_c2fortran_string'
        s = 'hello'//char(0)//'garbage'
        call syslib_c2fortran_string(s, len=l)
        call assert_int(5, l, 'c2fortran_string: length up to null terminator')
        call assert_char('hello', s(1:l), 'c2fortran_string: content up to null terminator')
    end subroutine test_syslib_c2fortran_string

    !---------------- find_next_int_dir_prefix ----------------

    subroutine test_find_next_int_dir_prefix()
        type(string) :: base, d1, d2, last
        integer :: next
        write(*,'(A)') 'test_find_next_int_dir_prefix'
        base = '.'
        d1 = '1run'
        d2 = '3run'
        call simple_mkdir(d1, verbose=.false.)
        call simple_mkdir(d2, verbose=.false.)
        next = find_next_int_dir_prefix(base, last_prev_dir=last)
        call assert_int(4, next, 'find_next_int_dir_prefix: returns max+1 (3+1=4)')
        call assert_true(last%to_char() == d2%to_char() .or. last%to_char() == d1%to_char(), &
            'find_next_int_dir_prefix: last_prev_dir is one of the integer dirs')
        call simple_rmdir(d1)
        call simple_rmdir(d2)
    end subroutine test_find_next_int_dir_prefix

    !---------------- print_slurm_env (smoke test) ----------------

    subroutine test_print_slurm_env_smoke()
        write(*,'(A)') 'test_print_slurm_env_smoke'
        call print_slurm_env()
        ! just ensure no crash
        call assert_true(.true., 'print_slurm_env: smoke test')
    end subroutine test_print_slurm_env_smoke

end module simple_syslib_tester
