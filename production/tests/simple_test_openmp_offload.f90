! OMP_TARGET_OFFLOAD=MANDATORY simple_test_openmp_offload nthr= device=
program simple_test_openmp_offload
use simple_core_module_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_gpu_utils
use omp_lib
implicit none
type(cmdline)      :: cline
type(parameters)   :: p
integer(dp) :: t
real(dp)    :: rt_cpu, rt_gpu, rt_cpus
integer :: num_devices,nteams,nthreads, device_id,ndevices
integer :: default_device, host_device
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
call set_offload_device(p%device)

#ifdef USE_OPENMP_OFFLOAD
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

call numerical_test
call persistence_test
call async_test

contains

    subroutine persistence_test
        integer,  parameter :: n = 1000000
        real(dp), dimension(:), allocatable :: a, b
        integer :: i,j,k
        print *, 'Starting persistence test'
        allocate(a(n), b(n))
        a = 1.0d0
        b = 2.0d0
        !$omp target enter data map(to: a, b)
        !$omp target teams distribute parallel do map(alloc: a, b)
        do i = 1, n
            a(i) = a(i) + b(i)
        end do
        !$omp end target teams distribute parallel do
        ! a & b should still be on device, a set to 3
        ! next setting b to 6
        !$omp target teams distribute parallel do map(alloc: a, b)
        do i = 1, n
            b(i) = a(i) * 2.0d0
        end do
        !$omp end target teams distribute parallel do
        !$omp target exit data map(from: a, b) 
        print *, '  Final result (b(1)) : ', b(1)
        print *, '  Final result (b(N)) : ', b(N)
        if( abs(b(1) - 6.0d0) > 1e-10 ) then
            stop 'Fatal error: persistence test failed, data did not persist on device'
        endif
        if( abs(b(N) - 6.0d0) > 1e-10 ) then
            stop 'Fatal error: persistence test failed, data did not persist on device'
        endif
        deallocate(a, b)
        print *, 'Persistence test successful'
    end subroutine persistence_test

    subroutine numerical_test
        integer, parameter :: N = 10000
        real, allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
        integer :: i,j
        print *, 'Starting numerical test and benchmark'
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
        !$omp loop collapse(2)
        do j = 1,N
        do i = 1,N
                c(i,j) = log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
            enddo
        enddo
        !$omp end loop
        !$omp loop collapse(2)
        do j = 1,N
        do i = 1,N
                c(i,j) = c(i,j) + log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
            enddo
        enddo
        !$omp end loop
        !$omp loop collapse(2)
        do j = 1,N
        do i = 1,N
                c(i,j) = c(i,j) - log(1. + 2.*sqrt(a(i,j)) + exp(-b(i,j)))**2.0
            enddo
        enddo
        !$omp end loop
        !$omp loop collapse(2)
        do j = 1,N
            do i = 1,N
                c(i,j) = 3. * exp(-cos(c(i,j) * pi))
            enddo
        enddo
        !$omp end loop
        !$omp end target teams
        rt_gpu = toc(t)
        print *, '  Numerical test c(1,1) = ', c(1,1), ' c(N,N) = ', c(N,N)
        print *, '  Numerical test d(1,1) = ', d(1,1), ' d(N,N) = ', d(N,N)
        if( any(abs(c-d) > 1e-4) ) then
            stop 'Fatal error: numerical test failed, results do not match'
        endif
        print *, 'Numerical test successful'
        print *,rt_cpu,  ' seconds on 1 CPU'
        print *,rt_cpus, ' seconds on CPUs'
        print *,rt_gpu,  ' seconds on GPU'
    end subroutine numerical_test

    subroutine async_test
        integer,    parameter :: n = 500000
        integer,    parameter :: repeats = 10
        real(dp), allocatable :: a(:), b(:), c(:)
        real(dp) :: tcpu, tgpu, tserial, tasync
        integer  :: i,j
        allocate(a(n),b(n),c(n))
        a(:) = (/(real(i,dp)/real(n,dp), i = 1, n)/)
        c(:) = a

        ! gpu only
        t = tic()
        !$omp target teams map(to: a) private(i,j)
        !$omp loop
        do j = 1, repeats
        do i = 1, n
            a(i) = sqrt(log(1. +(2.0*exp(-a(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        !$omp end loop
        !$omp end target teams
        tgpu = toc(t)
        print *,'  GPU execution: ', tgpu, ' seconds'

        ! cpu only
        b(:) = c
        t = tic()
        do j = 1, 5*repeats
        do i = 1, n
            b(i) = sqrt(log(1. +(2.0*exp(-b(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        tcpu = toc(t)
        print *,'  CPU execution: ', tcpu, ' seconds'

        ! asynchronous
        a(:) = c; b(:) = c
        t = tic()
        !$omp parallel sections num_threads(2) private(i,j)
        !$omp section
        !$omp target teams map(to: a)
        !$omp loop
        do j = 1,repeats
        do i = 1, n
            a(i) = sqrt(log(1. +(2.0*exp(-a(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        !$omp end loop
        !$omp end target teams
        !$omp section
        do j = 1, 5*repeats
        do i = 1, n
            b(i) = sqrt(log(1. +(2.0*exp(-b(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        !$omp end parallel sections
        tasync = toc(t)
        print *,'  Asynchronous execution: ', tasync, ' seconds'

        ! serial
        a(:) = c; b(:) = c
        t = tic()
        !$omp target teams map(to: a) private(i,j)
        !$omp loop
        do j = 1, repeats
        do i = 1, n
            a(i) = sqrt(log(1. +(2.0*exp(-a(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        !$omp end loop
        !$omp end target teams
        do j = 1, 5*repeats
        do i = 1, n
            b(i) = sqrt(log(1. +(2.0*exp(-b(i)-0.001))**2)*3.) - 0.1
        end do
        end do
        tserial = toc(t)
        print *,'  Serial execution: ', tserial, ' seconds'
        deallocate(a,b,c)
        if( tserial > 0.9*(tcpu+tgpu) .and. tasync>0.9*(max(tcpu,tgpu)) ) then
            print *,'Asynchronous execution test successful'
        else
            stop 'Fatal error: asynchronous execution should not be slower than serial execution'
        endif
    end subroutine async_test

#endif
end program simple_test_openmp_offload

