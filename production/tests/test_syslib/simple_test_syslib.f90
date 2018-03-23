
program simple_test_syslib
include 'simple_lib.f08'
implicit none

! interface
! subroutine get_fftw_version(L, index, k) bind(c, name="fftw_version")
!   use, intrinsic :: iso_c_binding
!   type(c_ptr), value :: L
!   integer(kind=c_int), value :: index
!   character(kind=c_char), dimension(*) :: k
! end subroutine get_fftw_version
! interface

real :: hbwsize
integer :: policy, io_stat, n,i
character(len=:),allocatable :: res(:), aname
character(len=STDLEN) :: simple_path_str, cur_working_dir
integer(8) :: version(3)
integer(8) :: vmRSS, vmHWM, vmPeak, vmSize
logical :: passed
passed= .true.

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
#if _OpenACC
    print *, '  OpenACC identified. '
#endif
#if _OpenMP
    print *, '  OpenMP identified. '
#endif

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
write(*,'(A,I0,A,I0,A,I0)') '    COMPILER VERSION: ', FC_COMPILER_VERSION(1), &
    '.',FC_COMPILER_VERSION(2), '.',FC_COMPILER_VERSION(3)
print *, ''

print *, '>>>'
print *, '>>> Syslib functions: Fortran Name / C Name(if used)'
print *, '>>>'

print *, '>>>'
print *, '>>> Syslib function Test 1: simple_isenv '
print *, '>>>'
if( simple_isenv('SIMPLE_PATH') )then
    print *, '     simple_isenv found SIMPLE_PATH'
else
    print *, '     simple_isenv failed to find SIMPLE_PATH'
    call simple_stop("simple_isenv failed")
end if
print *, '>>>'
print *, '>>> Syslib function Test 1: simple_getenv/getenv '
print *, '>>>'
io_stat = simple_getenv('SIMPLE_PATH',simple_path_str)
if (io_stat /= 0) call simple_stop("simple_getenv failed")
if(len_trim(simple_path_str)==0) then
    print*,'SIMPLE_PATH not defined in shell env'
else
     print *, '     simple_getenv found SIMPLE_PATH: ', trim(simple_path_str)
end if


print *, '>>>'
print *, '>>> Syslib function Test 2: simple_getcwd/getcwd '
print *, '>>>'
call simple_getcwd (cur_working_dir)
print *, '     simple_getcwd returned ', cur_working_dir

print *, '>>>'
print *, '>>> Syslib function Test 3: simple_mkdir / mkdir'
print *, '>>>'
call simple_mkdir('test_syslib', ignore=.false., status=io_stat)
if (io_stat /= 0) then
    call simple_stop("simple_mkdir failed ")
else
    print *, '    simple_mkdir test_syslib success'
endif
call simple_mkdir('test_syslib/temp1/temp2', ignore=.false., status=io_stat)
if (io_stat /= 0) then
    call simple_stop("simple_mkdir failed ")
else
    print *, '    simple_mkdir test_syslib/temp1/temp2 success'
endif

print *, '>>>'
print *, '>>> Syslib function Test 4: simple_rmdir '
print *, '>>>'
call simple_rmdir('test_syslib', status=io_stat)
if (io_stat /= 0) then
    call simple_stop("simple_rmdir failed ")
else
    print *, '    simple_rmdir test_syslib success'
endif

!! New interface: get_file_list, get_subdir_list, subprocess, glob_file_list, show_dir_content_recursive
!! New OS calls:  simple_list_dirs, simple_list_files, simple_rmdir, simple_del_files, exec_subprocess

print *, '>>>'
print *, '>>> Syslib function Test 5: simple_list_files / get_file_list -- Root directory'
print *, '>>>'
if (io_stat == 0 .and. allocated(res))then
    print *,"    file list num files:  size ", size(res)
    print *,"    files: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files '/' failed ")
    passed=.false.
endif
if(allocated(res))deallocate(res)


print *, '>>>'
print *, '>>> Syslib function Test 5: simple_list_files / get_file_list -- Current directory Empty '
print *, '>>>'
res = simple_list_files ("", status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"    file list num files: ", n, " size ", size(res)
    print *,"    files: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files '' failed ")
    passed=.false.
endif
if(allocated(res))deallocate(res)
print *, '>>>'
print *, '>>> Syslib function Test 5: simple_list_files / get_file_list -- Current directory with slash '
print *, '>>>'
res = simple_list_files ("./", status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"    num files: ", n, " size ", size(res)
    print *,"    files: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_list_files or get_file_list failed '
    call simple_stop("simple_list_files './' failed ")
    passed=.false.
endif
if(allocated(res))deallocate(res)
print *, '>>>'
print *, '>>> Syslib function Test 5: simple_list_files / glob_file_list -- single glob (this will emulate "ls *" ) '
print *, '>>>'
res = simple_list_files (glob="*", status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"    num files: ", size(res)
    print *,"    files: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_list_files or glob_file_list failed '
    call simple_stop("simple_list_files glob='*' failed ")
    passed=.false.
endif
if(allocated(res))deallocate(res)


print *, '>>>'
print *, '>>> Syslib function Test 5: simple_glob_list_tofile / glob_file_list -- advanced glob  '
print *, '>>>'
io_stat = simple_glob_list_tofile(glob=simple_path_str//"/scripts/*.pl", outfile='perl_scripts_list.tmp')
if (io_stat == 0 .and. allocated(res))then
    print *,"    num files: ", size(res)
    print *,"    files: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
    call del_file('perl_scripts_list.tmp')
else
    print *, '<<<  simple_glob_list_tofile or glob_file_list failed '
    call simple_stop("simple_glob_list_tofile failed")
    passed=.false.
endif
if(allocated(res))deallocate(res)


print *, '>>>'
print *, '>>> Syslib function Test 6:   simple_list_dirs / list_dirs '
print *, '>>>'
res = simple_list_dirs(".", status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"    num dirs: ", size(res)
    print *,"    dirs: "
      do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<   simple_list_dirs/get_file_list "." failed '
    call simple_stop('simple_list_dirs failed ')
    passed=.false.
endif
if(allocated(res))deallocate(res)

call exec_cmdline("simple_test_fileio 2>&1 > /dev/null")
print *, '>>>'
print *, '>>> Syslib function Test 6: simple_del_files / glob_file_list -- (emulate: rm -f SIMPLE_TEST_FILEIO*/tmp*.txt)'
print *, '>>>'
n = simple_del_files (glob="SIMPLE_TEST_FILEIO*/tmp*.txt", dlist=res, status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"   num of deleted files: ", n, " size ", size(res)
    print *,"   deleted files: "
    do i=1, size(res)
        if(file_exists(trim(res(i)))) write(*,'(a,": File still exists")')  trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_del_files / glob_file_list rm -f SIMPLE_TEST_FILEIO*/tmp*.txt FAILED'
    call simple_stop('simple_del_files  FAILED')
    passed=.false.
endif
if(allocated(res))deallocate(res)

call exec_cmdline("simple_test_fileio 2>&1 > /dev/null")
print *, '>>>'
print *, '>>> Syslib function Test 6: simple_rm_force / glob_file_list -- (emulate: rm -rf SIMPLE_TEST_FILEIO*)'
print *, '>>>'
n = simple_rm_force (glob="SIMPLE_TEST_FILEIO*", dlist=res, status=io_stat)
if (io_stat == 0 .and. allocated(res))then
    print *,"   num of deleted files and dirs: ", n, " size ", size(res)
    print *,"   deleted files and dirs: "
    do i=1, size(res)
        write(*,'(i30,":",a)') len_trim(res(i)), trim(adjustl(res(i)))
 !       print *, res(i)
    enddo
else
    print *, '<<<  simple_rm_force or glob_rm_rf failed '
    call simple_stop( ' simple_rm_force failed ')
    passed=.false.
endif
if(allocated(res))deallocate(res)


print *, '>>>'
print *, '>>> Syslib function Test 7: simple_sleep / sleep'
print *, '>>>'
call simple_sleep(0)
print *, '>>>'
print *, '>>> Syslib function Test 8: simple_mem_usage '
print *, '>>>'
call simple_mem_usage(vmRSS, vmPeak, vmSize, vmHWM)
print *, '>>> Memory Usage:'
print *, '>>>    vmRSS (Resident set size, i.e. current RAM used) ', vmRSS, 'kB'
print *, '>>>    vmHWM (High water mark RSS, max RAM by process ) ', vmHWM, 'kB'
print *, '>>>    vmSize (total shared Memory of process)          ', vmSize, 'kB'
print *, '>>>    vmPeak (peak vmSize)                             ', vmPeak, 'kB'
print *, '>>> Syslib function simple_dump_mem_usage: memory dumped to file '
call simple_dump_mem_usage()

print *, '>>>'
print *, '>>> Syslib function Test 9: simple_full_path / canonicalize_file_name'
print *, '>>>'
aname = simple_full_path("./")
print *," result of simple_absolute_path: ", trim(aname)
if(allocated(aname)) deallocate(aname)



! print *, '>>> Syslib function print_fftw_version  '
! call print_fftw_version()

#if defined(INTEL)
! print *,"  Is  High Bandwidth Memory  available? ", hbw_availability()
! hbwsize= get_hbw_size()
! print *,"  What is the High Bandwidth Memory size? ", hbwsize
! print *," Check Fast memory policy.."
! policy = fastmem_policy()

#endif

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

end program
