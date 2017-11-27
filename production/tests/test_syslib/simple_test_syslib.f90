
program test_syslib
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
integer :: policy, io_stat
character(len=STDLEN) :: simple_path_str
integer(8) :: version(3)


#if defined (PGI)
    print *, '  PGI  COMPILER identified. '
#elif defined (INTEL)
    print *, '  Intel  COMPILER identified. '
#elif defined (GNU)
    print *, '  GNU COMPILER identified. '
#else
    print *, ' NO COMPILER identified. '
#endif

print *, '>>> cmake variables defined in SimpleGitVersion.h:'
print *, '    SIMPLE build type: ',  trim(adjustl(BUILD_TYPE))
print *, '    SIMPLE build description: ',  trim(adjustl(build_descr))
print *, '    SIMPLE VERSION: ', trim(adjustl(SIMPLE_LIB_VERSION))
print *, '    SIMPLE_PATH  (installation directory): ', trim(adjustl(SIMPLE_PATH))
print *, '    SIMPLE source  directory: ', trim(adjustl(SIMPLE_SOURCE_PATH))
print *, '    SIMPLE build directory: ', trim(adjustl(SIMPLE_BUILD_PATH))
print *, '    COMPILER: ', trim(adjustl(FC_COMPILER))
write(*,'(A,I0,A,I0,A,I0)') '    COMPILER VERSION: ', FC_COMPILER_VERSION(1), &
    '.',FC_COMPILER_VERSION(2), '.',FC_COMPILER_VERSION(3)


print *, '>>> Syslib functions: '
print *, '>>> Syslib function simple_getenv '
io_stat = simple_getenv('SIMPLE_PATH',simple_path_str)
if (io_stat /= 0) call simple_error_check()
if(len_trim(simple_path_str)==0) then
    print*,'SIMPLE_PATH not defined in shell env'
end if

print *, '>>> Syslib function simple_sleep '
call simple_sleep(1)
print *, '>>> Syslib function print_compiler_info '
call print_compiler_info()

! print *, '>>> Syslib function print_fftw_version '
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

end program test_syslib
