program test_syslib
use simple_defs
use simple_syslib
implicit none
real :: hbwsize
integer :: policy
character(len=STDLEN), allocatable :: simple_path_str
simple_path_str = simple_getenv('SIMPLE_PATH')
if(.not. allocated(simple_path_str)) stop 'SIMPLE_PATH not defined in shell env'

call simple_sleep(1)

call print_compiler_info()

#if defined(INTEL)
print *,"  Is  High Bandwidth Memory  available? ", hbw_availability()
hbwsize= get_hbw_size()
print *,"  What is the High Bandwidth Memory size? ", hbwsize
print *," Check Fast memory policy.."
policy = fastmem_policy()

#endif
end program test_syslib
