program simple_test_syslib
include 'simple_lib.f08'
implicit none

#include "simple_local_flags.inc"
#include "simple_timer.h"
real :: hbwsize, te
integer :: policy, io_stat, n,i, pid
character(len=:),allocatable :: res(:), aname
character(len=STDLEN) :: simple_path_str, cur_working_dir
integer(8) :: version(3), t1, t2
integer(8) :: vmRSS, vmHWM, vmPeak, vmSize
logical :: passed, dir_e
character(len=8)      :: datestr
character(len=STDLEN) :: folder,  oldCWDfolder, curDir
character(len=30)     :: fname
integer(kind=8) :: valueRSS,valuePeak,valueSize,valueHWM
global_verbose=.true.
call seed_rnd
call date_and_time(date=datestr)
folder = trim('./SIMPLE_TEST_SYSLIB_'//datestr)

call simple_mkdir( trim(folder) , status=io_stat)
if(io_stat/=0) call simple_stop("mkdir failed")
print *," Changing directory to ", folder
call simple_chdir( trim(folder),  oldCWDfolder)
call simple_getcwd(curDir)
print *," Current working directory ", curDir
print *," Previous working directory ", oldCWDfolder

passed= .true.
global_debug=.true.
global_verbose=.true.
print *, '>>>'
print *, '>>> Preprocessor defined variables '
print *, '>>>'
#if defined (PGI)
print *, '  PGI  COMPILER identified. '
#elif defined (INTEL)
print *, '  Intel  COMPILER identified. '
#elif defined (GNU)
print *, '  GNU COMPILER identified. '
#else
print *, '  NO COMPILER identified. '
#endif

if(FORTRAN_COMPILER == FC_GNU_COMPILER)then
    print *, '  Fortran Compiler: Gfortran '
else if(FORTRAN_COMPILER == FC_PGI_COMPILER) then
    print *, '  Platform OS: PGI. '
else if(FORTRAN_COMPILER == FC_INTEL_COMPILER) then
    print *, '  Platform OS: Intel. '
endif

#if _OpenACC
print *, '  OpenACC identified. '
#endif
#if _OpenMP
 print *, '  OpenMP identified. '
#endif
if(OS_PLATFORM == OS_LINUX)then
    print *, '  Platform OS: Linux. '
else if(OS_PLATFORM == OS_MACOSX) then
    print *, '  Platform OS: MacOSX. '
endif


print *, '>>>'
print *, '>>> Syslib function print_compiler_info '
print *, '>>>'
call print_compiler_info()

print *, '>>>'
print *, '>>> cmake variables defined in SimpleGitVersion.h: '
print *, '>>>'
print *, ''
print *, '    SIMPLE build type: ',  trim(adjustl(BUILD_TYPE))
print *, '    SIMPLE build description: ',  trim(adjustl(build_descr))
print *, '    SIMPLE VERSION: ', trim(adjustl(SIMPLE_LIB_VERSION))
print *, '    SIMPLE_PATH  (installation directory): ', trim(adjustl(SIMPLE_PATH))
print *, '    SIMPLE source  directory: ', trim(adjustl(SIMPLE_SOURCE_PATH))
print *, '    SIMPLE build directory: ', trim(adjustl(SIMPLE_BUILD_PATH))
print *, '    COMPILER: ', trim(adjustl(FC_COMPILER))
write(*,'(A,I0,A,I0,A,I0)') '     COMPILER VERSION: ', FC_COMPILER_VERSION(1), &
    '.',FC_COMPILER_VERSION(2), '.',FC_COMPILER_VERSION(3)
print *, ''

print *, '>>>'
print *, '>>> Syslib functions: Fortran Name / C Name(if used)'
print *, '>>>'
print *, '>>> Syslib function Test 1: ENVIRONMENT FUNCTIONS '
print *, '>>> Syslib function Test 1a: simple_isenv '
if( simple_isenv('SIMPLE_PATH') )then
    print *, '     simple_isenv found SIMPLE_PATH'
else
    print *, '     simple_isenv failed to find SIMPLE_PATH'
    call simple_stop("simple_isenv failed",__FILENAME__,__LINE__)
end if

print *, '>>> Syslib function Test 1b: simple_getenv/getenv '

io_stat = simple_getenv('SIMPLE_PATH',simple_path_str)
if (io_stat /= 0) call simple_stop("simple_getenv failed",__FILENAME__,__LINE__)
if(len_trim(simple_path_str)==0) then
    call simple_stop("simple_getenv   SIMPLE_PATH str len zero")
else
    print *, '     simple_getenv found SIMPLE_PATH: ', trim(simple_path_str)
end if

print *, '>>> Syslib function Test 1c: simple_getcwd/getcwd '

call simple_getcwd (cur_working_dir)
print *, '      CWD:', cur_working_dir

print *, '>>> Syslib function Test 1d: simple_full_path / canonicalize_file_name'
aname = simple_full_path("./")
print *,'      absolute path "./": ', trim(aname)
if(allocated(aname)) deallocate(aname)

print *, '>>>'
print *, '>>> Syslib function Test 2: FILE OPERATIONS FUNCTIONS '
print *, '>>>'
print *, '>>> Syslib function Test 2a: simple_touch '

call simple_touch(trim('test_syslib.file'))
!call exec_cmdline("if [ ! -f test_syslib.file ];then  echo ' simple_touch FAILED!!!!';else echo 'simple_touch WORKED';fi " )
if  (.not. file_exists('test_syslib.file'))then
    call simple_stop("simple_touch failed ",__FILENAME__,__LINE__)
else
    print *, '      simple_touch success'
endif

print *, '>>> Syslib function Test 2b: simple_rename '

io_stat = simple_rename(trim('test_syslib.file'), trim('test_syslib.file2'))
!call exec_cmdline("if [ -f test_syslib.file ] && [ ! -f test_syslib.file2 ];&
!    &then  echo ' simple_rename FAILED!!!!';else echo 'simple_rename WORKED';fi " )
if (io_stat /= 0 .or. (.not. file_exists('test_syslib.file2')).or. (file_exists('test_syslib.file'))) then
    call simple_stop("simple_rename failed ",__FILENAME__,__LINE__)
else
    print *, '      simple_rename success'
endif
call exec_cmdline("rm  -f test_syslib.file* " )


print *, '>>> Syslib function Test 2c: simple_copy_file '
call exec_cmdline("rm -rf testcopy*.mrc SIMPLE_TEST_FILEIO_*; simple_test_fileio 2>/dev/null >/dev/null")
TBLOCK()
call simple_copy_file(trim('SIMPLE_TEST_FILEIO_'//datestr//'/cubes.mrc'), trim('testcopy.mrc'), status=io_stat)
TSTOP()
if (io_stat /= 0 .or. (.not. file_exists('testcopy.mrc'))) then
    call simple_stop("simple_copy_file failed ",__FILENAME__,__LINE__)
else
    call exec_cmdline("diff -q testcopy.mrc SIMPLE_TEST_FILEIO_"//datestr//"/cubes.mrc" )
endif
TBLOCK()
call syslib_copy_file(trim('SIMPLE_TEST_FILEIO_'//datestr//'/cubes.mrc'), trim('testcopy1.mrc'), status=io_stat)
TSTOP()
if (io_stat /= 0 .or. (.not. file_exists('testcopy1.mrc'))) then
    call simple_stop("simple_copy_file failed ",__FILENAME__,__LINE__)
else
    call exec_cmdline("diff -q testcopy1.mrc SIMPLE_TEST_FILEIO_"//datestr//"/cubes.mrc" )
endif

TBLOCK()
call syslib_copy_file2(trim('SIMPLE_TEST_FILEIO_'//datestr//'/cubes.mrc'), trim('testcopy2.mrc'), status=io_stat)
TSTOP()
if (io_stat /= 0 .or. (.not. file_exists('testcopy2.mrc'))) then
    call simple_stop("simple_copy_file failed ",__FILENAME__,__LINE__)
else
    call exec_cmdline("diff -q testcopy2.mrc SIMPLE_TEST_FILEIO_"//datestr//"/cubes.mrc" )
endif
TBLOCK()
call syslib_copy_file_stream(trim('SIMPLE_TEST_FILEIO_'//datestr//'/cubes.mrc'), trim('testcopy3.mrc'), status=io_stat)
TSTOP()
if (io_stat /= 0 .or. (.not. file_exists('testcopy3.mrc'))) then
    call simple_stop("simple_copy_file failed ",__FILENAME__,__LINE__)
else
    call exec_cmdline("diff -q testcopy3.mrc SIMPLE_TEST_FILEIO_"//datestr//"/cubes.mrc" )
endif
! TBLOCK()
! call syslib_copy_file_direct(trim('SIMPLE_TEST_FILEIO_'//datestr//'/cubes.mrc'), trim('testcopy4.mrc'), status=io_stat)
! TSTOP()
! !call exec_cmdline("if [ ! -f testcopy.mrc ];then  &
! !    &echo ' simple_copy_file FAILED!!!!' exit 1; &
! !    &else  echo 'simple_copy_file WORKED'; fi")
! call exec_cmdline("diff -q testcopy4.mrc SIMPLE_TEST_FILEIO_"//datestr//"/cubes.mrc" )
! if (io_stat /= 0 .or. (.not. file_exists('testcopy4.mrc'))) then
!     call simple_stop("simple_copy_file failed ",__FILENAME__,__LINE__)
! else
!     print *, '    simple_copy_file success'
! endif
call exec_cmdline("rm  -f testcopy*.mrc; rm -rf SIMPLE_TEST_FILEIO* " )


print *, '>>>'
print *, '>>> Syslib function Test 3: DIRECTORY FUNCTIONS '
print *, '>>>'
print *, '>>> Syslib function Test 3a: isdir '
print *, '     isdir(CWD) ', isdir(trim(adjustl(cur_working_dir)), len_trim(adjustl(cur_working_dir)))
print *, '>>> Syslib function Test 3b: dir_exists '
dir_e = dir_exists(cur_working_dir)
print *, '     dir_exists (expecting true)', dir_e
call exec_cmdline("rm -rf test_syslib")
dir_e = dir_exists('test_syslib')
print *, '     dir_exists (expecting false) ', dir_e


print *, '>>> Syslib function Test 3c: simple_mkdir / mkdir'
call simple_mkdir(trim('test_syslib'), ignore=.false., status=io_stat)
if (io_stat /= 0 .or. (.not. dir_exists('test_syslib'))) then
    call simple_stop("simple_mkdir failed ",__FILENAME__,__LINE__)
else
    print *, '      simple_mkdir success ✓'
endif


print *, '>>> Syslib function Test 3d: simple_rmdir '

call simple_rmdir('test_syslib', status=io_stat)
if (io_stat /= 0 .or. (dir_exists('test_syslib'))) then
    call simple_stop("simple_rmdir failed ",__FILENAME__,__LINE__)
else
    print *, '    simple_rmdir test_syslib removal success'
endif
call exec_cmdline("[ -d test_syslib ] && echo ' simple_rmdir FAILED!!!!' " )


print *, '>>> Syslib function Test 3e: simple_mkdir / mkdir  test_syslib/temp1/temp2'
call simple_mkdir('test_syslib/temp1/temp2', ignore=.false., status=io_stat)
if (io_stat /= 0 .or. (.not. dir_exists('test_syslib/temp1/temp2'))) then
    call simple_stop("simple_mkdir failed ",__FILENAME__,__LINE__)
else
    print *, '    simple_mkdir test_syslib/temp1/temp2 success'
endif

print *, '>>> Syslib function Test 3f: simple_rmdir Recursive'

call simple_rmdir('test_syslib/temp1/temp2', status=io_stat)
if (io_stat /= 0 .or. (dir_exists('test_syslib/temp1/temp2'))) then
    call simple_stop("simple_rmdir failed ",__FILENAME__,__LINE__)
else
    print *, '    simple_rmdir test_syslib/temp1/temp2 removal success'
endif
call simple_rmdir('./test_syslib', status=io_stat)
if (io_stat /= 0 .or. (dir_exists('test_syslib'))) then
    call simple_stop("simple_rmdir non-empty failed ",__FILENAME__,__LINE__)
else
    print *, '    simple_rmdir test_syslib non-empty removal success'
endif

call exec_cmdline("mkdir -p dir1/dir1.\{1,2\}; &
    &mkdir -p dir1/dir1.2/dir1.2.1; &
    &mkdir -p dir2/dir2.1; &
    &mkdir -p dir2/dir2.2/dir2.2.\{1,2\}; &
    &mkdir -p dir3/dir3.1; mkdir -p dir\{4,5\}; &
    &touch dir1/dir1.1/file.txt; &
    &touch dir1/dir1.2/file.txt; &
    &touch dir2/dir2.2/file.c; &
    &touch dir2/dir2.2/dir2.2.2/file.f90; &
    &touch dir3/file.js; &
    &touch dir3/dir3.1/file.c; touch file1.txt file2.txt")

print *, '>>>'
print *, '>>> Syslib function Test 4: FILE LISTING FUNCITONS '
print *, '>>>'
print *, '>>> Syslib function Test 4a: simple_list_files / get_file_list -- Root directory'
res = simple_list_files ("/", status=io_stat)
if (io_stat /= 0 .or. .not.allocated(res) .or. (file_exists('__simple_filelist__'))) then
    print *, '      simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files '/' failed ",__FILENAME__,__LINE__)
    passed=.false.
else
    print *, '      simple_list_files: args("/")                  :success'
    if(debug)then
        print *,"    file list num files:  size ", size(res)
        print *,"    files: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)


print *, '>>> Syslib function Test 4b: simple_list_files / get_file_list -- Current directory Empty '
res = simple_list_files ("", status=io_stat)
if  (io_stat /= 0 .or. .not.allocated(res) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files '' failed ",__FILENAME__,__LINE__)
    passed=.false.
else
    print *, '      simple_list_files: args("")                   :success'
    if(debug)then
        print *,"    file list num files: ", n, " size ", size(res)
        print *,"    files: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)


print *, '>>> Syslib function Test 4c: simple_list_files / get_file_list -- Current directory with slash '
call exec_cmdline("simple_test_fileio 2> /dev/null > /dev/null ")
res = simple_list_files ("./", status=io_stat)
if (io_stat /= 0 .or. .not.allocated(res) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files './' failed ",__FILENAME__,__LINE__)
    passed=.false.
else
    print *, '      simple_list_files: args("./")                 :success'

    if(debug)then
        print *,"    num files: ", n, " size ", size(res)
        print *,"    files: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)

print *, '>>> Syslib function Test 4d: simple_list_files / glob_file_list -- single glob (this will emulate "ls *" ) '
call exec_cmdline("simple_test_fileio 2> /dev/null > /dev/null")
res = simple_list_files (glob="SIMPLE_TEST_FILEIO*/*.txt", status=io_stat)
if (io_stat /= 0 .or. .not.allocated(res) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_list_files or glob_file_list failed '
    call simple_stop("simple_list_files glob='*' failed ",__FILENAME__,__LINE__)
    passed=.false.
else
    print *, '      simple_list_files: args(glob="*")             :success'

    if(debug)then
        print *,"    num files: ", size(res)
        print *,"    files: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)


print *, '>>> Syslib function Test 4e: simple_glob_list_tofile / glob_file_list -- advanced glob  '
call del_file('tmp_file_list.txt')
io_stat = simple_glob_list_tofile(glob="SIMPLE_TEST_FILEIO*/*.txt", outfile='tmp_file_list.txt')
if  (io_stat /= 0 .or. (.not.file_exists('tmp_file_list.txt')))then
    print *, '<<<  simple_glob_list_tofile or glob_file_list failed '
    call simple_stop("simple_glob_list_tofile failed",__FILENAME__,__LINE__)
    passed=.false.
else
    print *, '      simple_glob_list_tofile: args("SIMPLE_TEST_FILEIO*/*.txt", outfile="tmp_file_list.txt")    :success'

    if(debug)then
        print *,"    num files: ", size(res)
        print *,"    files: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
        call exec_cmdline('cat tmp_file_list.txt')
    endif
    call del_file('tmp_file_list.txt')
endif
if(allocated(res))deallocate(res)


print *, '>>> Syslib function Test 4f:   simple_list_dirs / list_dirs '
res = simple_list_dirs(".", status=io_stat)
if (io_stat /= 0 .or. .not.allocated(res) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<   simple_list_dirs/get_file_list "." failed '
    call simple_stop('simple_list_dirs failed ',__FILENAME__,__LINE__)
    passed=.false.
else
    if(debug)then
        print *,"    num dirs: ", size(res)
        print *,"    dirs: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)


print *, '>>>'
print *, '>>> Syslib function Test 5: DELETE FUNCTIONS '
print *, '>>>'
print *, '>>> Syslib function Test 5a: simple_del_files / glob_file_list -- (emulate: rm -f SIMPLE_TEST_FILEIO*/tmp*.txt)'
call exec_cmdline("simple_test_fileio 2> /dev/null > /dev/null ")
n = simple_del_files (glob="SIMPLE_TEST_FILEIO*/tmp*.txt", dlist=res, status=io_stat)
if  (io_stat /= 0 .or. (.not.allocated(res)) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_del_files / glob_file_list rm -f SIMPLE_TEST_FILEIO*/tmp*.txt FAILED'
    call simple_stop('simple_del_files  FAILED',__FILENAME__,__LINE__)
    passed=.false.
else
    if(debug)then
        print *,"   num of deleted files: ", n, " size ", size(res)
        print *,"   deleted files: "
        do i=1, size(res)
            if(file_exists(trim(res(i)))) write(*,'(a,": File still exists")')  trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)



print *, '>>> Syslib function Test 5b: simple_rm_force / glob_file_list -- (emulate: rm -rf SIMPLE_TEST_FILEIO*)'

call exec_cmdline("simple_test_fileio 2> /dev/null > /dev/null")
n = simple_rm_force (glob="SIMPLE_TEST_FILEIO*", dlist=res, status=io_stat)
if  (io_stat /= 0 .or. (.not.allocated(res)) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_rm_force or glob_rm_rf failed '
    call simple_stop( ' simple_rm_force failed ',__FILENAME__,__LINE__)
    passed=.false.
else
    if(debug)then
        print *,"   num of deleted files and dirs: ", n, " size ", size(res)
        print *,"   deleted files and dirs: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)


print *, '>>> Syslib function Test 5c: simple_rm_force / glob_rm_all -- (emulate: rm -rf ./dir*)'

n = simple_rm_force (glob="./dir*", dlist=res, status=io_stat)
if  (io_stat /= 0 .or. (.not.allocated(res)) .or. (file_exists('__simple_filelist__'))) then
    print *, '<<<  simple_rm_force or glob_rm_rf failed '
    call simple_stop( ' simple_rm_force failed ',__FILENAME__,__LINE__)
    passed=.false.
else
    if(debug)then
        print *,"   num of deleted files and dirs: ", n, " size ", size(res)
        print *,"   deleted files and dirs: "
        do i=1, size(res)
            write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
            !       print *, res(i)
        enddo
    endif
endif
if(allocated(res))deallocate(res)

print *, '>>>'
print *, '>>> Syslib function Test 6: Memory check system calls '
print *, '>>>'
print *, '>>> Syslib function Test 6a: sysinfo '
call simple_sysinfo_usage(valueRSS,valuePeak,valueSize,valueHWM)
print *," simple_sysinfo_usage :"
print *," Total usable main memory size (bytes):", valuePeak
print *," Amount of shared memory          ", valueSize
print *," Memory used by buffers           ", valueRSS
print *," High water mark (shared buffers) ", valueHWM

print *, '>>> Syslib function Test 6b: simple_mem_usage '
call simple_mem_usage(vmRSS, vmPeak, vmSize, vmHWM)
print *, '>>> Memory Usage:'
print *, '>>>    vmRSS (Resident set size, i.e. current RAM used) ', vmRSS, 'kB'
print *, '>>>    vmHWM (High water mark RSS, max RAM by process ) ', vmHWM, 'kB'
print *, '>>>    vmSize (total shared Memory of process)          ', vmSize, 'kB'
print *, '>>>    vmPeak (peak vmSize)                             ', vmPeak, 'kB'
print *, '>>> Syslib function Test 6c: simple_dump_mem_usage: memory dumped to file '
call simple_dump_mem_usage()


print *, '>>>'
print *, '>>> Syslib function Test 7: Proccess Management '
print *, '>>>'
print *, '>>> Syslib function Test 7a: exec_subprocess / subprocess'
call exec_subprocess('echo "" && echo "Beginning exec_subprocess call" && '//&
    &'echo "Hostname: $(hostname)" && '//&
    &'echo "Operating System: $(uname -a)" && '//&
    &'echo "Who is Logged in: $(who -u)" &&'//&
    &'echo "Current Directory: $(pwd)" &&'//&
    &' ls -al && '//&
    &'echo "Ending exec_subprocess call"  ', pid)
print *, '   subprocess returned PID ', pid
print *, '>>>  Syslib function Test 7b: wait_pid'
io_stat = wait_pid(pid)
print *, '      wait_pid status', io_stat

! print *, '>>> Syslib function print_fftw_version  '
! call print_fftw_version()

#if defined(INTEL)
! print *,"  Is  High Bandwidth Memory  available? ", hbw_availability()
! hbwsize= get_hbw_size()
! print *,"  What is the High Bandwidth Memory size? ", hbwsize
! print *," Check Fast memory policy.."
! policy = fastmem_policy()
#endif

print *, '>>>'
print *, '>>> Syslib tests completed successfully'
print *, '>>>'

! contains
!     function fftw3_lib_version() result(result_str)
!         use iso_c_binding
!         implicit none
!         character(len=STDLEN) :: result_str
!         character,pointer,dimension(:) :: fftw_version_str
!         integer :: i,sz
!         CALL C_F_POINTER(get_fftw_version(), fftw_version_str, [ 255 ])
!         sz=len_trim(fft_version_str)
!         do i=1,sz
!             result_str(i:i+1)=fftw_version_str(i)
!             endif
!         write(*,'(a,a)') 'FFTW library version ', trim(fftw_version_str)
!     end

!   use iso_c_binding, only: C_CHAR, C_NULL_CHAR
!   interface
!     subroutine print_c(string) bind(C, name="print_C")
!       use iso_c_binding, only: c_char
!       character(kind=c_char) :: string(*)
!     end subroutine print_c
!   end interface
!   call print_c(C_CHAR_"Hello World"//C_NULL_CHAR)

end program simple_test_syslib
