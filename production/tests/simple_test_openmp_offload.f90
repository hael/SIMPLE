program simple_test_openmp_offload
use simple_core_module_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
implicit none
integer, parameter :: N = 10000
type(cmdline)      :: cline
type(parameters)   :: p
integer(dp) :: t
real(dp)    :: rt_cpu, rt_gpu, rt_cpus
real, allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
integer :: i,j,k, num_devices,nteams,nthreads, device_id,ndevices
integer :: default_device, host_device, initial_device
logical :: is_host
if( command_argument_count() < 2 )then
    write(logfhandle,'(a)') 'ERROR! Usage: simple_test_openmp_offload device=x nthr=y'
else
  call cline%parse_oldschool
endif
call cline%checkvar('nthr',    1)
call cline%checkvar('device',  2)
call cline%check
call p%new(cline)

! do we have at least on device?
ndevices = omp_get_num_devices()
print *, 'Number of devices:                            ', ndevices
if( ndevices < 1 ) then
    stop 'Fatal error: No device identified'
endif
if( (p%device > ndevices) .or. (p%device < 0) ) then
    stop 'Fatal error: Invalid device ID provided'
endif

! Host device sanity checks
! omp_get_device_num() on the host equals omp_get_initial_device();
! use omp_get_initial_device() to avoid missing Fortran binding in libgomp.so
initial_device = omp_get_initial_device()
host_device    = initial_device
print *, 'Host device ID from omp_get_initial_device:   ', initial_device
if( initial_device /= host_device ) then
    stop 'Fatal error: incompatible device IDs'
endif
if( host_device == p%device ) then
    stop 'Fatal error: params%device must be the ID of the GPU, not the CPU host'
endif

! #ifdef should not be required, temporary
#ifdef USE_OPENMP_OFFLOAD
! Hooking up device
default_device = omp_get_default_device()
print *, 'Default device ID from omp_get_default_device:', default_device
print *, 'User provided device ID:                      ', p%device
call omp_set_default_device(p%device)

! pragma functionality test
!$omp target map(tofrom:nteams,nthreads,is_host,device_id)
is_host   = omp_is_initial_device()
nteams    = omp_get_num_teams()
nthreads  = omp_get_num_threads()
device_id = omp_get_device_num() 
!$omp end target 
if( is_host ) then
    stop 'Fatal error: offloading failed, code is running on host'
endif
print *, 'Device ID from omp_get_device_num:            ', device_id
if( device_id /= p%device ) then
    stop 'Fatal error: offloading failed, could not set the correct device'
endif
print *, 'Setup test successful'

! numerical test

allocate(a(N,N), b(N,N), c(N,N), d(N,N), source=0.)
a(:,1) = (/(real(10*i), i=1,N)/) / real(N)
do i = 2,N
    a(:,i) = a(:,1)
enddo
b = a / 2.

! single cpu
t = tic()
d =     log(1. + 2.*sqrt(a) + exp(-b))**2.0
d = d + log(1. + 2.*sqrt(a) + exp(-b))**2.0
d = d - log(1. + 2.*sqrt(a) + exp(-b))**2.0
d = 3. * exp(-cos(d * pi))
rt_cpu = toc(t)

! openmp cpu
t = tic()
!$omp parallel private(i,j) default(shared) proc_bind(close) 
!$omp do collapse(2) schedule(static)
do j = 1,N
  do i = 1,N
        c(i,j) = log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end do
!$omp do collapse(2) schedule(static)
do j = 1,N
  do i = 1,N
        c(i,j) = c(i,j) + log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end do
!$omp do collapse(2) schedule(static)
do j = 1,N
  do i = 1,N
        c(i,j) = c(i,j) - log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end do
!$omp do collapse(2) schedule(static)
do j = 1,N
    do i = 1,N
        c(i,j) = 3. * exp(-cos(c(i,j) * pi))
    enddo
enddo
!$omp end do
!$omp end parallel
rt_cpus = toc(t)

! openmp offload
t = tic()
!$omp target teams map(to:a,b) map(tofrom:c)
!$omp loop
do j = 1,N
  do i = 1,N
        c(i,j) = log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end loop
!$omp loop
do j = 1,N
  do i = 1,N
        c(i,j) = c(i,j) + log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end loop
!$omp loop
do j = 1,N
  do i = 1,N
        c(i,j) = c(i,j) - log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
    enddo
enddo
!$omp end loop
!$omp loop
do j = 1,N
    do i = 1,N
        c(i,j) = 3. * exp(-cos(c(i,j) * pi))
    enddo
enddo
!$omp end loop
!$omp end target teams
rt_gpu = toc(t)
print *, 'Numerical test c(1,1) = ', c(1,1), ' c(N,N) = ', c(N,N)
print *, 'Numerical test d(1,1) = ', d(1,1), ' d(N,N) = ', d(N,N)
if( any(abs(c-d) > 1e-4) ) then
    stop 'Fatal error: numerical test failed, results do not match'
endif
print *, 'Numerical test successful'
print *,rt_cpu,  ' seconds on 1 CPU'
print *,rt_cpus, ' seconds on CPUs'
print *,rt_gpu,  ' seconds on GPU'
#endif
end program simple_test_openmp_offload

