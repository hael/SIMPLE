
program test_syslib
#include "simple_lib.f08"


implicit none
real :: hbwsize
integer :: policy
character(len=STDLEN), allocatable :: simple_path_str
integer(8) :: version(3)

#if defined (PGI)
    print *, '  PGI  COMPILER found. '
#elif defined (INTEL)
    print *, '  Intel  COMPILER found. '
#elif defined (GNU)
    print *, '  GNU COMPILER found. '
#else
    print *, ' NO COMPILER found. '
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
simple_path_str = simple_getenv('SIMPLE_PATH')
if(.not. allocated(simple_path_str)) stop 'SIMPLE_PATH not defined in shell env'
print *, '>>> Syslib function simple_sleep '
call simple_sleep(1) 
print *, '>>> Syslib function print_compiler_info '
call print_compiler_info()

#if defined(INTEL)
! print *,"  Is  High Bandwidth Memory  available? ", hbw_availability()
! hbwsize= get_hbw_size()
! print *,"  What is the High Bandwidth Memory size? ", hbwsize
! print *," Check Fast memory policy.."
! policy = fastmem_policy()

#endif
end program test_syslib
