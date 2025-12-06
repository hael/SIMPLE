module simple_fileio_tester
use simple_defs
use simple_string
use simple_string_utils
use simple_syslib
use simple_fileio
use simple_test_utils
implicit none
private
public :: run_all_fileio_tests

contains

    subroutine run_all_fileio_tests()
        write(*,'(A)') '**** running all simple_fileio tests ****'
        call test_fname2ext_and_basename()
        call test_add2fbody_rm_swap_getfbody()
        call test_make_dir_and_file_names()
        call test_filepath_overloads()
        call test_fname2iter()
        call test_fname2format()
        call test_get_fpath_and_stemname()
        call test_get_relative_path()
        call test_fopen_and_fclose_basic()
        call test_nlines_and_filelength()
        call test_arr2file_sp_dp_roundtrip()
        call test_arr2txtfile_roundtrip()
        call test_rmat_lmat_roundtrip()
        call test_read_write_filetable()
        call test_write_singlelineoftext_and_read_exit_code()
        call test_simple_copy_file()
        ! --- NEW: OS / syslib-related tests ---
        call test_simple_mkdir_dir_exists_chdir_getcwd_rmdir()
        call test_simple_touch_rename_abspath()
        call test_simple_list_files_and_regexp()
        ! call report_summary()
    end subroutine run_all_fileio_tests

    !------------------------------------------------------------
    !  Name / extension helpers
    !------------------------------------------------------------

    subroutine test_fname2ext_and_basename()
        type(string) :: s, ext, base
        write(*,'(A)') 'test_fname2ext_and_basename'
        s    = string('/tmp/myfile.dat')
        ext  = fname2ext(s)
        base = basename(s)
        call assert_char('dat',  ext%to_char(),        'fname2ext: .dat')
        call assert_char('myfile.dat', base%to_char(), 'basename: myfile.dat')
        s    = string('noext')
        ext  = fname2ext(s)
        call assert_true(ext .eq. '',                  'fname2ext: empty string')
        call s%kill
        call ext%kill
        call base%kill
    end subroutine test_fname2ext_and_basename

    subroutine test_add2fbody_rm_swap_getfbody()
        type(string) :: fname, suf, ins, res, rmstr
        write(*,'(A)') 'test_add2fbody_rm_swap_getfbody'
        fname = string('path/body_old.ext')
        suf   = string('.ext')
        ins   = string('_new')
        res = add2fbody(fname, suf, ins)
        call assert_char('path/body_old_new.ext', res%to_char(), 'add2fbody: insert before suffix')
        rmstr = string('_new')
        res = rm_from_fbody(res, suffix=suf, str=rmstr)
        call assert_char('path/body_old.ext', res%to_char(),     'rm_from_fbody: body without insert')
        res = swap_suffix(res, suffix=string('.dat'), old_suffix=suf)
        call assert_char('path/body_old.dat', res%to_char(),     'swap_suffix: replace .ext by .dat')
        res = get_fbody(res, suffix=string('dat'))
        call assert_char('path/body_old', res%to_char(),         'get_fbody: body before .ext')
        call fname%kill
        call suf%kill
        call ins%kill
        call res%kill
        call rmstr%kill
    end subroutine test_add2fbody_rm_swap_getfbody

    subroutine test_make_dir_and_file_names()
        type(string), allocatable :: dirs(:), files(:)
        integer :: i
        write(*,'(A)') 'test_make_dir_and_file_names'
        dirs = make_dirnames('run_', 3)
        call assert_int(3, size(dirs), 'make_dirnames: size=3')
        call assert_char('run_1', dirs(1)%to_char(), 'make_dirnames: first')
        call assert_char('run_3', dirs(3)%to_char(), 'make_dirnames: last')
        files = make_filenames('img_', 3, '.dat')
        call assert_int(3, size(files), 'make_filenames: size=3')
        ! default numlen depends on len(int2str(n)), so for n=3 we expect 1 digit
        call assert_char('img_1.dat', files(1)%to_char(), 'make_filenames: img_1.dat')
        call assert_char('img_3.dat', files(3)%to_char(), 'make_filenames: img_3.dat')
        do i=1,size(dirs)
            call dirs(i)%kill
        end do
        do i=1,size(files)
            call files(i)%kill
        end do
        deallocate(dirs, files)
    end subroutine test_make_dir_and_file_names

    subroutine test_filepath_overloads()
        type(string) :: s1, s2, s3, fname
        write(*,'(A)') 'test_filepath_overloads'
        s1 = string('/tmp')
        s2 = string('subdir')
        s3 = string('file.txt')
        fname = filepath(s1, s2)
        call assert_char('/tmp'//PATH_SEPARATOR//'subdir', fname%to_char(), 'filepath_1: 2 args')
        fname = filepath(s1, s2, s3)
        call assert_char('/tmp'//PATH_SEPARATOR//'subdir'//PATH_SEPARATOR//'file.txt', fname%to_char(), 'filepath_1: 3 args')
        fname = filepath(s1, 'name')
        call assert_char('/tmp'//PATH_SEPARATOR//'name', fname%to_char(), 'filepath_2: string,char')
        fname = filepath('/root', 'x')
        call assert_char('/root'//PATH_SEPARATOR//'x', fname%to_char(), 'filepath_3: char,char')
        fname = filepath('/data', s2)
        call assert_char('/data'//PATH_SEPARATOR//'subdir', fname%to_char(), 'filepath_4: char,string')
        call s1%kill
        call s2%kill
        call s3%kill
        call fname%kill
    end subroutine test_filepath_overloads

    subroutine test_fname2iter()
        type(string) :: s
        integer :: it
        write(*,'(A)') 'test_fname2iter'
        s = string('run_iter0100_extra.mrc')
        it = fname2iter(s)
        call assert_int(100, it, 'fname2iter: 0100 -> 100')
        s = string('noiter_here.mrc')
        it = fname2iter(s)
        call assert_int(0, it, 'fname2iter: no iter -> 0')
        call s%kill
    end subroutine test_fname2iter

    subroutine test_fname2format()
        type(string) :: s
        character(1) :: fmt
        write(*,'(A)') 'test_fname2format'
        s = string('image.img')
        fmt = fname2format(s)
        call assert_char('I', fmt, 'fname2format: img -> I')
        s = string('data.mrc')
        fmt = fname2format(s)
        call assert_char('M', fmt, 'fname2format: mrc -> M')
        s = string('list.txt')
        fmt = fname2format(s)
        call assert_char('T', fmt, 'fname2format: txt -> T')
        s = string('unknown.xyz')
        fmt = fname2format(s)
        call assert_char('N', fmt, 'fname2format: unknown -> N')
        call s%kill
    end subroutine test_fname2format

    subroutine test_get_fpath_and_stemname()
        type(string) :: fname, path, stem
        write(*,'(A)') 'test_get_fpath_and_stemname'
        fname = string('/tmp/subdir/file.txt')
        path = get_fpath(fname)
        call assert_char('/tmp/subdir'//PATH_SEPARATOR, path%to_char(), 'get_fpath: path including separator')
        stem = stemname(fname)
        call assert_char('/tmp/subdir', stem%to_char(), 'stemname: path without last element')
        fname = string('justfile')
        path = get_fpath(fname)
        call assert_char(PATH_HERE, path%to_char(), 'get_fpath: no separator -> PATH_HERE')
        call fname%kill
        call path%kill
        call stem%kill
    end subroutine test_get_fpath_and_stemname

    subroutine test_get_relative_path()
        type(string) :: full, root, rel
        integer :: trimlen
        write(*,'(A)') 'test_get_relative_path'
        full = string('/home/user/project/data/file.dat')
        root = string('project/')
        trimlen = 0
        rel  = get_relative_path(full, root, trimlength=trimlen)
        call assert_true(trimlen > 0, 'get_relative_path: trimlength > 0')
        call assert_char('project/data/file.dat', rel%to_char(), 'get_relative_path: substring from root')
        ! Case with no match -> path unchanged
        full = string('/other/path/file.dat')
        root = string('/home/user')
        trimlen = 0
        rel = get_relative_path(full, root)
        call assert_char('/other/path/file.dat', rel%to_char(), 'get_relative_path: no match -> unchanged')
    end subroutine test_get_relative_path

    !------------------------------------------------------------
    !  File IO helpers (Fortran-only)
    !------------------------------------------------------------

    subroutine test_fopen_and_fclose_basic()
        type(string) :: fname
        integer :: funit, ios
        write(*,'(A)') 'test_fopen_and_fclose_basic'
        fname = string('tmp_fileio_fopen_basic.dat')
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'fopen: basic replace/write iostat=0')
        call assert_true(is_funit_open(funit), 'fopen: unit is open')
        write(funit, '(A)') 'hello'
        call fclose(funit)
        call assert_true(.not. is_funit_open(funit), 'fclose: unit closed')
        call del_file(fname)
        call fname%kill
    end subroutine test_fopen_and_fclose_basic

    subroutine test_nlines_and_filelength()
        type(string) :: fname
        integer :: funit, ios
        integer :: n, lenbytes
        write(*,'(A)') 'test_nlines_and_filelength'
        fname = string('tmp_fileio_nlines.dat')
        call fopen(funit, fname, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'nlines: open for write')
        write(funit,'(A)') 'line1'
        write(funit,'(A)') '# comment'
        write(funit,'(A)') 'line2'
        call fclose(funit)
        n = nlines(fname)
        ! nlines reads character(1) and applies str_is_comment; here only 'line1' and 'line2' count
        call assert_int(2, n, 'nlines: counts non-comment lines')
        lenbytes = filelength(fname)
        call assert_true(lenbytes > 0, 'filelength: positive')
        call del_file(fname)
        call fname%kill
    end subroutine test_nlines_and_filelength

    subroutine test_arr2file_sp_dp_roundtrip()
        type(string) :: fsp, fdp
        real, allocatable     :: arr_sp(:), arr_sp_back(:)
        real(dp), allocatable :: arr_dp(:), arr_dp_back(:)
        integer :: i
        write(*,'(A)') 'test_arr2file_sp_dp_roundtrip'
        allocate(arr_sp(3))
        arr_sp = [1.0, 2.0, 3.5]
        allocate(arr_dp(3))
        arr_dp = [1.0_dp, 2.0_dp, 3.5_dp]
        fsp = string('tmp_arr2file_sp.bin')
        fdp = string('tmp_arr2file_dp.bin')
        call arr2file(arr_sp, fsp)
        call arr2file(arr_dp, fdp)
        arr_sp_back = file2rarr(fsp)
        arr_dp_back = file2drarr(fdp)
        call assert_int(3, size(arr_sp_back), 'file2rarr: size')
        call assert_int(3, size(arr_dp_back), 'file2drarr: size')
        do i=1,3
            call assert_true(abs(arr_sp(i)-arr_sp_back(i)) < 1.0e-6, 'file2rarr: element match')
            call assert_true(abs(arr_dp(i)-arr_dp_back(i)) < 1.0e-12_dp, 'file2drarr: element match')
        end do
        call del_file(fsp)
        call del_file(fdp)
        call fsp%kill
        call fdp%kill
        deallocate(arr_sp, arr_sp_back, arr_dp, arr_dp_back)
    end subroutine test_arr2file_sp_dp_roundtrip

    subroutine test_arr2txtfile_roundtrip()
        type(string) :: fsp, fint
        real,    allocatable :: arr_sp(:), arr_sp_back(:)
        integer, allocatable :: arr_i(:), arr_i_back(:)
        integer :: i, funit, ios
        write(*,'(A)') 'test_arr2txtfile_roundtrip'
        allocate(arr_sp(3))
        arr_sp = [1.0, 2.5, -3.0]
        allocate(arr_i(3))
        arr_i = [1, 2, 3]
        fsp  = string('tmp_arr2txtfile_sp.txt')
        fint = string('tmp_arr2txtfile_int.txt')
        call arr2txtfile(arr_sp,  fsp)
        call arr2txtfile(arr_i,  fint)
        ! read back reals
        call fopen(funit, fsp, status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'arr2txtfile sp: open for read')
        allocate(arr_sp_back(3))
        do i=1,3
            read(funit,*,iostat=ios) arr_sp_back(i)
            call assert_int(0, ios, 'arr2txtfile sp: read iostat=0')
        end do
        call fclose(funit)
        ! read back ints
        call fopen(funit, fint, status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'arr2txtfile int: open for read')
        allocate(arr_i_back(3))
        do i=1,3
            read(funit,*,iostat=ios) arr_i_back(i)
            call assert_int(0, ios, 'arr2txtfile int: read iostat=0')
        end do
        call fclose(funit)
        do i=1,3
            call assert_true(abs(arr_sp(i)-arr_sp_back(i)) < 1.0e-6, 'arr2txtfile sp: element match')
            call assert_int(arr_i(i), arr_i_back(i), 'arr2txtfile int: element match')
        end do
        call del_file(fsp)
        call del_file(fint)
        call fsp%kill
        call fint%kill
        deallocate(arr_sp, arr_sp_back, arr_i, arr_i_back)
    end subroutine test_arr2txtfile_roundtrip

    subroutine test_rmat_lmat_roundtrip()
        type(string) :: fR, fL
        real,    allocatable :: rmat(:,:), rback(:,:)
        logical, allocatable :: lmat(:,:), lback(:,:)
        integer :: i,j
        write(*,'(A)') 'test_rmat_lmat_roundtrip'
        allocate(rmat(2,3))
        rmat = reshape([1.0,2.0,3.0,4.0,5.0,6.0], shape(rmat))
        allocate(lmat(2,2))
        lmat = reshape([.true.,.false.,.false.,.true.], shape(lmat))
        fR = string('tmp_rmat2file.bin')
        fL = string('tmp_lmat2file.bin')
        call rmat2file(rmat, fR)
        call lmat2file(lmat, fL)
        call file2rmat(fR, rback)
        call file2lmat(fL, lback)
        call assert_int(2, size(rback,1), 'file2rmat: dim1')
        call assert_int(3, size(rback,2), 'file2rmat: dim2')
        do j=1,3
            do i=1,2
                call assert_true(abs(rmat(i,j)-rback(i,j)) < 1.0e-6, 'file2rmat: element match')
            end do
        end do
        call assert_int(2, size(lback,1), 'file2lmat: dim1')
        call assert_int(2, size(lback,2), 'file2lmat: dim2')
        do j=1,2
            do i=1,2
                call assert_true(lmat(i,j) .eqv. lback(i,j), 'file2lmat: element match')
            end do
        end do
        call del_file(fR)
        call del_file(fL)
        call fR%kill
        call fL%kill
        deallocate(rmat, rback, lmat, lback)
    end subroutine test_rmat_lmat_roundtrip

    subroutine test_read_write_filetable()
        type(string) :: table, a, b, c, a_read, c_read
        type(string), allocatable :: names_out(:), names_in(:)
        integer :: i
        write(*,'(A)') 'test_read_write_filetable'
        table = string('tmp_filetable.txt')
        a = string('fileA.mrc')
        b = string('fileB.mrc')
        c = string('fileC.mrc')
        call simple_touch(a)
        call simple_touch(b)
        call simple_touch(c)
        allocate(names_out(3))
        names_out = [a, b, c]
        call write_filetable(table, names_out)
        call read_filetable(table, names_in)
        call assert_int(3, size(names_in), 'filetable: size=3')
        a_read = basename(names_in(1))
        c_read = basename(names_in(3))
        call assert_char(a%to_char(), a_read%to_char(), 'filetable: first')
        call assert_char(c%to_char(), c_read%to_char(), 'filetable: last')
        do i=1,size(names_out)
            call names_out(i)%kill
        end do
        do i=1,size(names_in)
            call names_in(i)%kill
        end do
        deallocate(names_out, names_in)
        call del_file(table)
        call table%kill
        call a%kill
        call b%kill
        call c%kill
        call a_read%kill
        call c_read%kill
    end subroutine test_read_write_filetable

    subroutine test_write_singlelineoftext_and_read_exit_code()
        type(string) :: fname_text, fname_exit
        type(string) :: line
        integer :: funit, ios, exit_code
        logical :: err
        write(*,'(A)') 'test_write_singlelineoftext_and_read_exit_code'
        fname_text = string('tmp_singleline.txt')
        fname_exit = string('tmp_exitcode.txt')
        line = string('hello world')
        call write_singlelineoftext(fname_text, line)
        call fopen(funit, fname_text, status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'write_singlelineoftext: reopen')
        call line%readline(funit, ios)
        call assert_int(0, ios, 'write_singlelineoftext: read')
        call fclose(funit)
        call assert_char('hello world', line%to_char(), 'write_singlelineoftext: content')
        ! exit code file: one integer in i32 format
        call fopen(funit, fname_exit, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'read_exit_code: open for write')
        write(funit,'(i32)') 42
        call fclose(funit)
        call read_exit_code(fname_exit, exit_code, err)
        call assert_true(.not. err, 'read_exit_code: err=.false.')
        call assert_int(42, exit_code, 'read_exit_code: value')
        call del_file(fname_text)
        call del_file(fname_exit)
        call fname_text%kill
        call fname_exit%kill
        call line%kill
    end subroutine test_write_singlelineoftext_and_read_exit_code

    subroutine test_simple_copy_file()
        type(string) :: src, dst
        integer :: funit, ios
        character(len=64) :: text
        character(len=64) :: text2
        write(*,'(A)') 'test_simple_copy_file'
        src = string('tmp_copy_src.dat')
        dst = string('tmp_copy_dst.dat')
        call fopen(funit, src, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'simple_copy_file: open src')
        write(text,'(A)') 'copy_me'
        write(funit,'(A)') trim(text)
        call fclose(funit)
        call simple_copy_file(src, dst)
        call fopen(funit, dst, status='old', action='read', iostat=ios)
        call assert_int(0, ios, 'simple_copy_file: open dst')
        read(funit,'(A)',iostat=ios) text2
        call assert_int(0, ios, 'simple_copy_file: read dst')
        call fclose(funit)
        call assert_char('copy_me', trim(text2), 'simple_copy_file: content identical')
        call del_file(src)
        call del_file(dst)
        call src%kill
        call dst%kill
    end subroutine test_simple_copy_file

    subroutine test_simple_mkdir_dir_exists_chdir_getcwd_rmdir()
        type(string) :: testdir, cwd_before, cwd_after
        integer      :: status
        type(string) :: base_after
        write(*,'(A)') 'test_simple_mkdir_dir_exists_chdir_getcwd_rmdir'
        ! Pick a simple scratch subdirectory name
        testdir = string('tmp_syslib_dir_test')
        ! Get current working directory
        call simple_getcwd(cwd_before)
        call assert_true(cwd_before%strlen_trim() > 0, 'simple_getcwd: cwd_before non-empty')
        ! Create directory
        call simple_mkdir(testdir)
        call assert_true(dir_exists(testdir), 'simple_mkdir + dir_exists: directory exists')
        ! Change into new directory
        call simple_chdir(testdir, status)
        call assert_int(0, status, 'simple_chdir: status==0')
        ! Check we’re in that directory via simple_getcwd + basename
        call simple_getcwd(cwd_after)
        base_after = basename(cwd_after)
        call assert_char(testdir%to_char(), base_after%to_char(), 'simple_chdir: basename(cwd_after)==testdir')
        ! Change back to original cwd
        call simple_chdir(cwd_before, status)
        call assert_int(0, status, 'simple_chdir back: status==0')
        ! Remove directory and confirm it’s gone
        call simple_rmdir(testdir, status)
        call assert_int(0, status, 'simple_rmdir: status==0 (or 0 when no error)')
        call assert_true(.not. dir_exists(testdir), 'simple_rmdir + dir_exists: directory removed')
        call testdir%kill
        call cwd_before%kill
        call cwd_after%kill
        call base_after%kill
    end subroutine test_simple_mkdir_dir_exists_chdir_getcwd_rmdir

    subroutine test_simple_touch_rename_abspath()
        type(string) :: f1, f2, abs2, base2
        integer :: status_here
        write(*,'(A)') 'test_simple_touch_rename_abspath'
        f1 = string('tmp_sys_touch_1.dat')
        f2 = string('tmp_sys_touch_2.dat')
        ! Touch first file
        call simple_touch(f1, status_here)
        call assert_int(0, status_here, 'simple_touch: status==0')
        call assert_true(file_exists(f1), 'simple_touch + file_exists: f1 exists')
        ! Rename to f2 (and overwrite if needed)
        call simple_rename(f1, f2, overwrite=.true.)
        call assert_true(.not. file_exists(f1), 'simple_rename: f1 no longer exists')
        call assert_true(file_exists(f2),       'simple_rename: f2 exists')
        ! Absolute path of f2 should have basename f2
        abs2 = simple_abspath(f2, status_here, check_exists=.true.)
        call assert_int(0, status_here, 'simple_abspath: status==0')
        base2 = basename(abs2)
        call assert_char(f2%to_char(), base2%to_char(), 'simple_abspath + basename: basename(abs2)==f2')
        call del_file(f2)
        call f1%kill
        call f2%kill
        call abs2%kill
        call base2%kill
    end subroutine test_simple_touch_rename_abspath

    subroutine test_simple_list_files_and_regexp()
        type(string) :: f1, f2, f3, dir
        type(string), allocatable :: list(:)
        integer :: funit, ios, i
        logical :: found1, found2, found3
        write(*,'(A)') 'test_simple_list_files_and_regexp'
        dir = string('.')
        f1  = string('a1.txt')
        f2  = string('a2.txt')
        f3  = string('b1.log')
        ! Create three small files
        call fopen(funit, f1, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'list_files: open a1')
        write(funit,'(A)') 'a1'
        call fclose(funit)
        call fopen(funit, f2, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'list_files: open a2')
        write(funit,'(A)') 'a2'
        call fclose(funit)
        call fopen(funit, f3, status='replace', action='write', iostat=ios)
        call assert_int(0, ios, 'list_files: open b1')
        write(funit,'(A)') 'b1'
        call fclose(funit)
        ! ---- simple_list_files with pattern '*.txt' ----
        call simple_list_files('*.txt', list)
        call assert_true(allocated(list), 'simple_list_files: list allocated')
        found1 = .false.
        found2 = .false.
        found3 = .false.
        do i=1,size(list)
            if (list(i)%to_char() == f1%to_char()) found1 = .true.
            if (list(i)%to_char() == f2%to_char()) found2 = .true.
            if (list(i)%to_char() == f3%to_char()) found3 = .true.
        end do
        call assert_true(found1, 'simple_list_files: finds a1.txt')
        call assert_true(found2, 'simple_list_files: finds a2.txt')
        call assert_true(.not. found3, 'simple_list_files: does not include b1.log')
        if (allocated(list)) then
            do i=1,size(list)
                call list(i)%kill
            end do
            deallocate(list)
        end if
        ! ---- simple_list_files_regexp on "." with regexp for a*.txt ----
        call simple_list_files_regexp(dir, '^a[0-9]+\.txt$', list)
        ! The implementation prefixes with dir%to_char()//'/' so we expect './a1.txt', './a2.txt'
        found1 = .false.
        found2 = .false.
        found3 = .false.
        do i=1,size(list)
            if (list(i)%has_substr('a1.txt')) found1 = .true.
            if (list(i)%has_substr('a2.txt')) found2 = .true.
            if (list(i)%has_substr('b1.log')) found3 = .true.
        end do
        call assert_true(found1, 'simple_list_files_regexp: finds a1.txt')
        call assert_true(found2, 'simple_list_files_regexp: finds a2.txt')
        call assert_true(.not. found3, 'simple_list_files_regexp: does not include b1.log')
        if (allocated(list)) then
            do i=1,size(list)
                call list(i)%kill
            end do
            deallocate(list)
        end if
        ! Cleanup
        call del_file(f1)
        call del_file(f2)
        call del_file(f3)
        call f1%kill
        call f2%kill
        call f3%kill
        call dir%kill
    end subroutine test_simple_list_files_and_regexp

end module simple_fileio_tester
