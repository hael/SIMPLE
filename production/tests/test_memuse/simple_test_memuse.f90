program simple_test_memuse
include 'simple_lib.f08'
implicit none

integer(kind=8) :: valueRSS, valuePeak, valueSize, valueHWM
real, allocatable  :: a(:)
integer, parameter :: SZ = 1000

allocate(a(SZ))

print *, 'sizeof(a) ', sizeof(a)

! call simple_sysinfo_usage(valueRSS,valuePeak,valueSize,valueHWM)
call simple_mem_usage(valueRSS,valuePeak,valueSize,valueHWM)

print *, 'valueRSS',  valueRSS
print *, 'valuePeak', valuePeak
print *, 'valueSize', valueSize
print *, 'valueHWM',  valueHWM

 call simple_dump_mem_usage

end program simple_test_memuse
