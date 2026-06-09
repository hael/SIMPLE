! OMP_TARGET_OFFLOAD=MANDATORY simple_test_openmp_offload nthr= device=
program simple_test_openmp_offload
use, intrinsic :: iso_c_binding, only: c_float, c_int, c_ptr, c_loc, c_f_pointer, c_associated, c_null_ptr
use simple_core_module_api
use simple_parameters, only: parameters
use simple_cmdline,    only: cmdline
use simple_gpu_utils
use simple_fftw3
use omp_lib
use omp_lib_kinds
implicit none

type dummytype
    integer           :: i
    real, allocatable :: a(:)
end type dummytype

type(cmdline)      :: cline
type(parameters)   :: p
integer(dp) :: t
real(dp)    :: rt_cpu, rt_gpu, rt_cpus
integer :: nteams, nthreads, device_id
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

print *, 'From OpenMP Number of teams  : ', nteams
print *, 'From OpenMP Number of threads: ', nthreads
call print_gpu_specs(p%device)
call numerical_test
call persistence_test
call persistence_test2
call async_test
call test_fftw_vs_cufft_1d_in_place
call test_fftw_vs_cufft_2d_in_place
call test_fftw_vs_cufft_2d_in_place_roundtrip_timing
call test_fftw_vs_cufft_3d_in_place_roundtrip_timing
call test_fftw_vs_cufft_1d_many
call test_fftw_vs_cufft_2d_many
call test_fftw_vs_cufft_1d_many_c2r
call test_fftw_vs_cufft_movie
call test_cublas_sgemm
call test_kb
call test_pointers
call banner('All tests passed successfully')
contains

    subroutine persistence_test
        integer,  parameter :: n = 1000000
        real(dp), dimension(:), allocatable :: a, b
        integer :: i
        call banner('TEST: Memory persistence test')
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

    subroutine persistence_test2
        integer,  parameter :: n = 256
        real(sp), dimension(:), allocatable :: b
        type(dummytype) :: dt
        integer :: i
        call banner('TEST: Memory persistence test of type bound allocatable')
        ! dummy type init
        dt%i = 42
        allocate(dt%a(n),source=1.)
        allocate(b(n))
        b = 2.
        !$omp target enter data map(to: dt%a, b)
        !$omp target teams distribute parallel do map(alloc: dt%a, b)
        do i = 1, n
            dt%a(i) = dt%a(i) + b(i)
        end do
        !$omp end target teams distribute parallel do
        ! a & b should still be on device, a set to 3
        ! next setting b to 6
        !$omp target teams distribute parallel do map(alloc: dt%a, b)
        do i = 1, n
            b(i) = dt%a(i) * 2.0d0
        end do
        !$omp end target teams distribute parallel do
        !$omp target exit data map(from: dt%a, b) 
        print *, '  Final result (b(1)) : ', b(1)
        print *, '  Final result (b(N)) : ', b(N)
        if( abs(b(1) - 6.0d0) > 1e-10 ) then
            stop 'Fatal error: persistence test2 failed, data did not persist on device'
        endif
        if( abs(b(N) - 6.0d0) > 1e-10 ) then
            stop 'Fatal error: persistence test2 failed, data did not persist on device'
        endif
        deallocate(dt%a, b)
        print *, 'Persistence test2 successful'
    end subroutine persistence_test2

    subroutine numerical_test
        integer, parameter :: N = 10000
        real, allocatable :: a(:,:), b(:,:), c(:,:), d(:,:)
        integer :: i,j
        call banner('TEST: Numerical test and benchmark')
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
        call banner('TEST: Asynchronous execution test')
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

    subroutine banner(msg)
        character(len=*), intent(in) :: msg
        print *
        print *, '============================================================'
        print *, trim(msg)
        print *, '============================================================'
    end subroutine banner

    subroutine report_fftw(label, err, tol_local, passed)
        character(len=*), intent(in) :: label
        real(sp), intent(in) :: err, tol_local
        logical, intent(out) :: passed
        passed = (err <= tol_local)
        if (passed) then
            print *, 'PASS: ', trim(label), '  max error = ', err, '  tol = ', tol_local
        else
            print *, 'FAIL: ', trim(label), '  max error = ', err, '  tol = ', tol_local
            stop 0
        end if
    end subroutine report_fftw

    function relerr_c1d(a, b) result(err)
        complex(sp), intent(in) :: a(:), b(:)
        real(sp) :: err, denom
        denom = max(1.0_sp, maxval(abs(b)))
        err = maxval(abs(a - b)) / denom
    end function relerr_c1d

    function relerr_c2d(a, b) result(err)
        complex(sp), intent(in) :: a(:,:), b(:,:)
        real(sp) :: err, denom
        denom = max(1.0_sp, maxval(abs(b)))
        err = maxval(abs(a - b)) / denom
    end function relerr_c2d

    subroutine print_cpair(label, a, b, nshow)
        character(len=*), intent(in) :: label
        complex(sp), intent(in) :: a(:), b(:)
        integer, intent(in) :: nshow
        integer :: k
        print *, trim(label)
        print *, '  index        CUFFT R2C output(real,imag)   FFTW R2C output(real,imag)'
        do k = 1, min(nshow, size(a), size(b))
          print '(A,I0,A,2(ES14.6,1X),A,2(ES14.6,1X))', '  [', k, ']  ', &
            real(a(k),sp), aimag(a(k)), '    ', real(b(k),sp), aimag(b(k))
        end do
    end subroutine print_cpair

    subroutine init_dummy_1d(x, n, offset)
        real(sp), intent(out) :: x(:)
        integer,  intent(in)  :: n
        real(sp), intent(in), optional :: offset
        real(sp) :: off
        integer  :: i
        off = 0.0_sp
        if (present(offset)) off = offset
        do i = 1, n
          x(i) = sin(2.0_sp*pi*real(i-1,sp)/real(n,sp)) + &
                 0.25_sp*cos(4.0_sp*pi*real(i-1,sp)/real(n,sp)) + off
        end do
    end subroutine init_dummy_1d

    subroutine init_dummy_2d(x, nx, ny, offset)
        real(sp), intent(out) :: x(:,:)
        integer,  intent(in)  :: nx, ny
        real(sp), intent(in), optional :: offset
        real(sp) :: off
        integer  :: i, j
        off = 0.0_sp
        if (present(offset)) off = offset
        do j = 1, ny
          do i = 1, nx
            x(i,j) = sin(2.0_sp*pi*real(i-1,sp)/real(nx,sp)) + &
                     0.5_sp*cos(2.0_sp*pi*real(j-1,sp)/real(ny,sp)) + off
          end do
        end do
    end subroutine init_dummy_2d

    subroutine init_dummy_3d(x, nx, ny, nz)
        real(sp), intent(out) :: x(:,:,:)
        integer,  intent(in)  :: nx, ny, nz
        integer  :: i, j, k
        do k = 1, nz
          do j = 1, ny
            do i = 1, nx
              x(i,j,k) = sin(2.0_sp*pi*real(i-1,sp)/real(nx,sp)) + &
                         0.5_sp*cos(2.0_sp*pi*real(j-1,sp)/real(ny,sp)) + &
                         0.25_sp*sin(2.0_sp*pi*real(k-1,sp)/real(nz,sp))
            end do
          end do
        end do
    end subroutine init_dummy_3d

    subroutine test_fftw_vs_cufft_1d_in_place
        integer, parameter :: n1   = 16
        integer, parameter :: n1c  = n1/2 + 1
        integer, parameter :: n1pad = 2*n1c
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        integer(c_int) :: plan_cufft
        type(c_ptr)    :: plan_fftw, p_r, p_c
        integer(c_int) :: ierr
        real(sp), allocatable, target :: buf1_cufft(:), buf1_fftw(:)
        complex(sp), pointer :: buf1_cufft_c(:), buf1_fftw_c(:)
        real(sp) :: maxerr
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 1D IN-PLACE R2C')
        allocate(buf1_cufft(n1pad), buf1_fftw(n1pad))
        call c_f_pointer(c_loc(buf1_cufft(1)), buf1_cufft_c, [n1c])
        call c_f_pointer(c_loc(buf1_fftw(1)),  buf1_fftw_c,  [n1c])
        buf1_cufft = 0.0_sp
        buf1_fftw  = 0.0_sp
        call init_dummy_1d(buf1_cufft(1:n1), n1)
        buf1_fftw(1:n1) = buf1_cufft(1:n1)
        ! Plans
        ierr = cufftPlan1d(plan_cufft, int(n1,c_int), CUFFT_R2C, 1_c_int)
        ! print *, 'cufftPlan1d 1D IP comparison R2C ierr = ', ierr, ' handle = ', plan_cufft
        if (ierr /= CUFFT_SUCCESS) stop 58
        plan_fftw = fftwf_plan_dft_r2c_1d(int(n1,c_int), buf1_fftw, buf1_fftw_c, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw)) stop 59

        call fftwf_execute_dft_r2c(plan_fftw, buf1_fftw, buf1_fftw_c)
        !$omp target data map(tofrom:buf1_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf1_cufft(1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf1_cufft(1)), 0)
          ierr = cufftExecR2C(plan_cufft, p_r, p_c)
          if (ierr /= CUFFT_SUCCESS) stop 60
        !$omp end target data

        call print_cpair('1D in-place R2C output coefficients, all bins:', &
          buf1_cufft_c, buf1_fftw_c, n1c)
        maxerr = relerr_c1d(buf1_cufft_c, buf1_fftw_c)
        print *, 'FFTW/CUFFT 1D in-place relative R2C output error = ', maxerr
        call report_fftw('FFTW VS CUFFT 1D IN-PLACE R2C', maxerr, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= CUFFT_SUCCESS) stop 62
        call fftwf_destroy_plan(plan_fftw)
        deallocate(buf1_cufft, buf1_fftw)
    end subroutine test_fftw_vs_cufft_1d_in_place

    subroutine test_fftw_vs_cufft_2d_in_place
        integer, parameter :: nx     = 4096
        integer, parameter :: ny     = 4096
        integer, parameter :: n2cx   = nx/2 + 1
        integer, parameter :: n2padx = 2*n2cx
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        integer(c_int) :: plan_cufft
        type(c_ptr)    :: plan_fftw, p_r, p_c
        integer(c_int) :: ierr
        real(sp), allocatable, target :: buf2_cufft(:,:), buf2_fftw(:,:)
        complex(sp), pointer :: buf2_cufft_c(:,:), buf2_fftw_c(:,:)
        real(sp) :: maxerr
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 2D IN-PLACE R2C')
        allocate(buf2_cufft(n2padx,ny), buf2_fftw(n2padx,ny))
        call c_f_pointer(c_loc(buf2_cufft(1,1)), buf2_cufft_c, [n2cx,ny])
        call c_f_pointer(c_loc(buf2_fftw(1,1)),  buf2_fftw_c,  [n2cx,ny])
        buf2_cufft = 0.0_sp
        buf2_fftw  = 0.0_sp
        call init_dummy_2d(buf2_cufft(1:nx,1:ny), nx, ny)
        buf2_fftw(1:nx,1:ny) = buf2_cufft(1:nx,1:ny)
        ! plan
        ierr = cufftPlan2d(plan_cufft, int(ny,c_int), int(nx,c_int), CUFFT_R2C)
        ! print *, 'cufftPlan2d 2D IP comparison R2C ierr = ', ierr, ' handle = ', plan_cufft
        if (ierr /= CUFFT_SUCCESS) stop 63
        plan_fftw = fftwf_plan_dft_r2c_2d(int(ny,c_int), int(nx,c_int), &
          buf2_fftw, buf2_fftw_c, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw)) stop 64
        ! execute
        call fftwf_execute_dft_r2c(plan_fftw, buf2_fftw, buf2_fftw_c)
        !$omp target data map(tofrom:buf2_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1)), 0)
          ierr = cufftExecR2C(plan_cufft, p_r, p_c)
          if (ierr /= CUFFT_SUCCESS) stop 65
        !$omp end target data
        ! report
        call print_cpair('2D in-place R2C output coefficients, row 1, first four bins:', &
          buf2_cufft_c(:,1), buf2_fftw_c(:,1), 4)
        maxerr = relerr_c2d(buf2_cufft_c, buf2_fftw_c)
        print *, 'FFTW/CUFFT 2D in-place relative R2C output error = ', maxerr
        call report_fftw('FFTW VS CUFFT 2D IN-PLACE R2C', maxerr, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= CUFFT_SUCCESS) stop 67
        call fftwf_destroy_plan(plan_fftw)
        deallocate(buf2_cufft, buf2_fftw)
    end subroutine test_fftw_vs_cufft_2d_in_place

    subroutine test_fftw_vs_cufft_2d_in_place_roundtrip_timing
        integer, parameter :: nx     = 4096
        integer, parameter :: ny     = 4096
        integer, parameter :: n2cx   = nx/2 + 1
        integer, parameter :: n2padx = 2*n2cx
        integer, parameter :: repeats = 5
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        real(sp), parameter :: scale = 1.0_sp / real(nx*ny, sp)
        integer(c_int) :: plan_cufft_fwd, plan_cufft_inv
        type(c_ptr)    :: plan_fftw_fwd, plan_fftw_inv, p_r, p_c
        integer(c_int) :: ierr
        integer        :: i, j, rep
        real(sp), allocatable, target :: buf2_cufft(:,:), buf2_fftw(:,:)
        real(sp), allocatable :: xref(:,:)
        complex(sp), pointer :: buf2_fftw_c(:,:)
        real(sp) :: maxerr, maxdiff
        real(dp) :: fftw_time, cufft_time
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 2D IN-PLACE R2C->C2R TIMING')
        allocate(buf2_cufft(n2padx,ny), buf2_fftw(n2padx,ny), xref(nx,ny))
        call c_f_pointer(c_loc(buf2_fftw(1,1)), buf2_fftw_c, [n2cx,ny])
        call init_dummy_2d(xref, nx, ny)
        buf2_cufft = 0.0_sp
        buf2_fftw  = 0.0_sp
        buf2_cufft(1:nx,1:ny) = xref
        buf2_fftw(1:nx,1:ny)  = xref
        ! Plans
        ierr = cufftPlan2d(plan_cufft_fwd, int(ny,c_int), int(nx,c_int), CUFFT_R2C)
        if (ierr /= CUFFT_SUCCESS) stop 84
        ierr = cufftPlan2d(plan_cufft_inv, int(ny,c_int), int(nx,c_int), CUFFT_C2R)
        if (ierr /= CUFFT_SUCCESS) stop 85
        plan_fftw_fwd = fftwf_plan_dft_r2c_2d(int(ny,c_int), int(nx,c_int), &
          buf2_fftw, buf2_fftw_c, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_fwd)) stop 86
        plan_fftw_inv = fftwf_plan_dft_c2r_2d(int(ny,c_int), int(nx,c_int), &
          buf2_fftw_c, buf2_fftw, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_inv)) stop 87
        ! execute
        t = tic()
        do rep = 1, repeats
          call fftwf_execute_dft_r2c(plan_fftw_fwd, buf2_fftw, buf2_fftw_c)
          call fftwf_execute_dft_c2r(plan_fftw_inv, buf2_fftw_c, buf2_fftw)
          buf2_fftw(1:nx,1:ny) = buf2_fftw(1:nx,1:ny) * scale
        end do
        fftw_time = toc(t)
        !$omp target data map(tofrom:buf2_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1)), 0)
          ierr = cudaDeviceSynchronize()
          if (ierr /= CUDA_SUCCESS) stop 88
          t = tic()
          do rep = 1, repeats
            ierr = cufftExecR2C(plan_cufft_fwd, p_r, p_c)
            if (ierr /= CUFFT_SUCCESS) stop 89
            ierr = cufftExecC2R(plan_cufft_inv, p_c, p_r)
            if (ierr /= CUFFT_SUCCESS) stop 90
            ierr = cudaDeviceSynchronize()
            if (ierr /= CUDA_SUCCESS) stop 91
            !$omp target teams distribute parallel do collapse(2)
            do j = 1, ny
              do i = 1, nx
                buf2_cufft(i,j) = buf2_cufft(i,j) * scale
              end do
            end do
            !$omp end target teams distribute parallel do
          end do
          cufft_time = toc(t)
        !$omp end target data
        ! report
        maxerr  = maxval(abs(buf2_cufft(1:nx,1:ny) - xref))
        maxdiff = maxval(abs(buf2_cufft(1:nx,1:ny) - buf2_fftw(1:nx,1:ny)))
        print *, 'FFTW 2D in-place R2C->C2R total time     = ', fftw_time, ' seconds'
        print *, 'CUFFT 2D in-place R2C->C2R total time    = ', cufft_time, ' seconds'
        print *, 'FFTW 2D in-place R2C->C2R avg time       = ', fftw_time / real(repeats,dp), ' seconds'
        print *, 'CUFFT 2D in-place R2C->C2R avg time      = ', cufft_time / real(repeats,dp), ' seconds'
        print *, 'CUFFT 2D in-place round-trip max error   = ', maxerr
        print *, 'FFTW/CUFFT 2D in-place round-trip diff   = ', maxdiff
        call report_fftw('FFTW VS CUFFT 2D IN-PLACE R2C->C2R', maxdiff, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft_fwd)
        if (ierr /= CUFFT_SUCCESS) stop 93
        ierr = cufftDestroy(plan_cufft_inv)
        if (ierr /= CUFFT_SUCCESS) stop 94
        call fftwf_destroy_plan(plan_fftw_fwd)
        call fftwf_destroy_plan(plan_fftw_inv)
        deallocate(buf2_cufft, buf2_fftw, xref)
    end subroutine test_fftw_vs_cufft_2d_in_place_roundtrip_timing

    subroutine test_fftw_vs_cufft_3d_in_place_roundtrip_timing
        integer, parameter :: nx     = 256
        integer, parameter :: ny     = 256
        integer, parameter :: nz     = 256
        integer, parameter :: n3cx   = nx/2 + 1
        integer, parameter :: n3padx = 2*n3cx
        integer, parameter :: repeats = 3
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        real(sp), parameter :: scale = 1.0_sp / real(nx*ny*nz, sp)
        integer(c_int) :: plan_cufft_fwd, plan_cufft_inv
        type(c_ptr)    :: plan_fftw_fwd, plan_fftw_inv, p_r, p_c
        integer(c_int) :: ierr
        integer        :: i, j, k, rep
        real(sp), allocatable, target :: buf3_cufft(:,:,:), buf3_fftw(:,:,:)
        real(sp), allocatable :: xref(:,:,:)
        complex(sp), pointer :: buf3_fftw_c(:,:,:)
        real(sp) :: maxerr, maxdiff
        real(dp) :: fftw_time, cufft_time
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 3D IN-PLACE R2C->C2R TIMING')
        allocate(buf3_cufft(n3padx,ny,nz), buf3_fftw(n3padx,ny,nz), xref(nx,ny,nz))
        call c_f_pointer(c_loc(buf3_fftw(1,1,1)), buf3_fftw_c, [n3cx,ny,nz])
        call init_dummy_3d(xref, nx, ny, nz)
        buf3_cufft = 0.0_sp
        buf3_fftw  = 0.0_sp
        buf3_cufft(1:nx,1:ny,1:nz) = xref
        buf3_fftw(1:nx,1:ny,1:nz)  = xref
        ! Planning
        ierr = cufftPlan3d(plan_cufft_fwd, int(nz,c_int), int(ny,c_int), int(nx,c_int), CUFFT_R2C)
        if (ierr /= CUFFT_SUCCESS) stop 95
        ierr = cufftPlan3d(plan_cufft_inv, int(nz,c_int), int(ny,c_int), int(nx,c_int), CUFFT_C2R)
        if (ierr /= CUFFT_SUCCESS) stop 96
        ierr = fftwf_init_threads()
        call fftwf_plan_with_nthreads(p%nthr)
        plan_fftw_fwd = fftwf_plan_dft_r2c_3d(int(nz,c_int), int(ny,c_int), int(nx,c_int), &
          buf3_fftw, buf3_fftw_c, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_fwd)) stop 97
        plan_fftw_inv = fftwf_plan_dft_c2r_3d(int(nz,c_int), int(ny,c_int), int(nx,c_int), &
          buf3_fftw_c, buf3_fftw, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_inv)) stop 98
        ! Execution
        t = tic()
        do rep = 1, repeats
          call fftwf_execute_dft_r2c(plan_fftw_fwd, buf3_fftw, buf3_fftw_c)
          call fftwf_execute_dft_c2r(plan_fftw_inv, buf3_fftw_c, buf3_fftw)
          buf3_fftw(1:nx,1:ny,1:nz) = buf3_fftw(1:nx,1:ny,1:nz) * scale
        end do
        fftw_time = toc(t)
        !$omp target data map(tofrom:buf3_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf3_cufft(1,1,1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf3_cufft(1,1,1)), 0)
          ierr = cudaDeviceSynchronize()
          if (ierr /= CUDA_SUCCESS) stop 99
          t = tic()
          do rep = 1, repeats
            ierr = cufftExecR2C(plan_cufft_fwd, p_r, p_c)
            if (ierr /= CUFFT_SUCCESS) stop 100
            ierr = cufftExecC2R(plan_cufft_inv, p_c, p_r)
            if (ierr /= CUFFT_SUCCESS) stop 101
            ierr = cudaDeviceSynchronize()
            if (ierr /= CUDA_SUCCESS) stop 102
            !$omp target teams distribute parallel do collapse(3)
            do k = 1, nz
              do j = 1, ny
                do i = 1, nx
                  buf3_cufft(i,j,k) = buf3_cufft(i,j,k) * scale
                end do
              end do
            end do
            !$omp end target teams distribute parallel do
          end do
          cufft_time = toc(t)
        !$omp end target data
        maxerr  = maxval(abs(buf3_cufft(1:nx,1:ny,1:nz) - xref))
        maxdiff = maxval(abs(buf3_cufft(1:nx,1:ny,1:nz) - buf3_fftw(1:nx,1:ny,1:nz)))
        print *, 'FFTW 3D in-place R2C->C2R total time     = ', fftw_time, ' seconds'
        print *, 'CUFFT 3D in-place R2C->C2R total time    = ', cufft_time, ' seconds'
        print *, 'FFTW 3D in-place R2C->C2R avg time       = ', fftw_time / real(repeats,dp), ' seconds'
        print *, 'CUFFT 3D in-place R2C->C2R avg time      = ', cufft_time / real(repeats,dp), ' seconds'
        print *, 'CUFFT 3D in-place round-trip max error   = ', maxerr
        print *, 'FFTW/CUFFT 3D in-place round-trip diff   = ', maxdiff
        call report_fftw('FFTW VS CUFFT 3D IN-PLACE R2C->C2R', maxdiff, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft_fwd)
        if (ierr /= CUFFT_SUCCESS) stop 104
        ierr = cufftDestroy(plan_cufft_inv)
        if (ierr /= CUFFT_SUCCESS) stop 105
        call fftwf_destroy_plan(plan_fftw_fwd)
        call fftwf_destroy_plan(plan_fftw_inv)
        deallocate(buf3_cufft, buf3_fftw, xref)
        call fftwf_plan_with_nthreads(1)
    end subroutine test_fftw_vs_cufft_3d_in_place_roundtrip_timing

    subroutine test_fftw_vs_cufft_1d_many
        integer, parameter :: n1 = 16
        integer, parameter :: n1c = n1/2 + 1
        integer, parameter :: n1pad = 2*n1c
        integer, parameter :: nbatch = 3
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        integer(c_int) :: plan_cufft
        type(c_ptr)    :: plan_fftw, p_r, p_c
        integer(c_int) :: ierr
        integer(c_int) :: n(1), inembed(1), onembed(1)
        integer(c_int) :: rank, howmany, istride, ostride, idist, odist
        integer :: batch_idx
        real(sp), allocatable, target :: buf1_cufft(:,:), buf1_fftw(:,:)
        complex(sp), pointer :: buf1_cufft_c(:,:), buf1_fftw_c(:,:)
        real(sp) :: maxerr
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 1D MANY IN-PLACE R2C')
        allocate(buf1_cufft(n1pad,nbatch), buf1_fftw(n1pad,nbatch))
        call c_f_pointer(c_loc(buf1_cufft(1,1)), buf1_cufft_c, [n1c,nbatch])
        call c_f_pointer(c_loc(buf1_fftw(1,1)),  buf1_fftw_c,  [n1c,nbatch])
        buf1_cufft = 0.0_sp
        buf1_fftw  = 0.0_sp
        do batch_idx = 1, nbatch
            call init_dummy_1d(buf1_cufft(1:n1,batch_idx), n1, 0.1_sp*real(batch_idx-1,sp))
        end do
        buf1_fftw(1:n1,:) = buf1_cufft(1:n1,:)
        ! Plan
        rank       = 1_c_int
        n(1)       = int(n1, c_int)
        howmany    = int(nbatch, c_int)
        inembed(1) = int(n1, c_int)
        onembed(1) = int(n1c, c_int)
        istride    = 1_c_int
        ostride    = 1_c_int
        idist      = int(n1pad, c_int)
        odist      = int(n1c, c_int)
        ierr = cufftPlanMany(plan_cufft, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C, howmany)
        ! print *, 'cufftPlanMany 1D comparison R2C ierr = ', ierr, ' handle = ', plan_cufft
        if (ierr /= CUFFT_SUCCESS) stop 68
        plan_fftw = fftwf_plan_many_dft_r2c(rank, n, howmany, buf1_fftw, inembed, istride, idist, buf1_fftw_c, onembed, ostride, odist, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw)) stop 69
        ! execute
        call fftwf_execute_dft_r2c(plan_fftw, buf1_fftw, buf1_fftw_c)
        !$omp target data map(tofrom:buf1_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf1_cufft(1,1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf1_cufft(1,1)), 0)
          ierr = cufftExecR2C(plan_cufft, p_r, p_c)
          if (ierr /= CUFFT_SUCCESS) stop 70
        !$omp end target data
        ! report
        do batch_idx = 1, nbatch
          call print_cpair('1D many R2C output coefficients, batch '//trim(int2str(batch_idx))//':', &
            buf1_cufft_c(:,batch_idx), buf1_fftw_c(:,batch_idx), n1c)
        end do
        maxerr = relerr_c2d(buf1_cufft_c, buf1_fftw_c)
        print *, 'FFTW/CUFFT 1D many relative R2C output error = ', maxerr
        call report_fftw('FFTW VS CUFFT 1D MANY R2C', maxerr, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= CUFFT_SUCCESS) stop 72
        call fftwf_destroy_plan(plan_fftw)
        deallocate(buf1_cufft, buf1_fftw)
    end subroutine test_fftw_vs_cufft_1d_many

    subroutine test_fftw_vs_cufft_2d_many
        integer, parameter :: nx = 32
        integer, parameter :: ny = 24
        integer, parameter :: n2cx = nx/2 + 1
        integer, parameter :: n2padx = 2*n2cx
        integer, parameter :: nbatch = 2
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        integer(c_int) :: plan_cufft
        type(c_ptr)    :: plan_fftw, p_r, p_c
        integer(c_int) :: ierr
        integer(c_int) :: n(2), inembed(2), onembed(2)
        integer(c_int) :: rank, howmany, istride, ostride, idist, odist
        integer :: batch_idx
        real(sp), allocatable, target :: buf2_cufft(:,:,:), buf2_fftw(:,:,:)
        complex(sp), pointer :: buf2_cufft_c(:,:,:), buf2_fftw_c(:,:,:)
        real(sp) :: maxerr
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 2D MANY IN-PLACE R2C')
        allocate(buf2_cufft(n2padx,ny,nbatch), buf2_fftw(n2padx,ny,nbatch))
        call c_f_pointer(c_loc(buf2_cufft(1,1,1)), buf2_cufft_c, [n2cx,ny,nbatch])
        call c_f_pointer(c_loc(buf2_fftw(1,1,1)),  buf2_fftw_c,  [n2cx,ny,nbatch])
        buf2_cufft = 0.0_sp
        buf2_fftw  = 0.0_sp
        do batch_idx = 1, nbatch
            call init_dummy_2d(buf2_cufft(1:nx,1:ny,batch_idx), nx, ny, 0.05_sp*real(batch_idx-1,sp))
        end do
        buf2_fftw(1:nx,1:ny,:) = buf2_cufft(1:nx,1:ny,:)
        ! Plans
        rank       = 2_c_int
        n          = [int(ny, c_int), int(nx, c_int)]
        howmany    = int(nbatch, c_int)
        inembed    = [int(ny, c_int), int(nx, c_int)]
        onembed    = [int(ny, c_int), int(n2cx, c_int)]
        istride    = 1_c_int
        ostride    = 1_c_int
        idist      = int(n2padx*ny, c_int)
        odist      = int(n2cx*ny, c_int)
        ierr = cufftPlanMany(plan_cufft, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_R2C_MANY, howmany)
        if (ierr /= CUFFT_SUCCESS) stop 73
        plan_fftw = fftwf_plan_many_dft_r2c(rank, n, howmany, buf2_fftw, inembed, istride, idist, buf2_fftw_c, onembed, ostride, odist, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw)) stop 74
        call fftwf_execute_dft_r2c(plan_fftw, buf2_fftw, buf2_fftw_c)
        !$omp target data map(tofrom:buf2_cufft)
            p_r = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1,1)), 0)
            p_c = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1,1)), 0)
            ierr = cufftExecR2C(plan_cufft, p_r, p_c)
            if (ierr /= CUFFT_SUCCESS) stop 75
        !$omp end target data
        do batch_idx = 1, nbatch
            call print_cpair('2D many R2C output coefficients, batch '//trim(int2str(batch_idx))//' row 1, first four bins:', &
            buf2_cufft_c(:,1,batch_idx), buf2_fftw_c(:,1,batch_idx), 4)
        end do
        maxerr = maxval(abs(buf2_cufft_c - buf2_fftw_c)) / max(1.0_sp, maxval(abs(buf2_fftw_c)))
        print *, 'FFTW/CUFFT 2D many relative R2C output error = ', maxerr
        call report_fftw('FFTW VS CUFFT 2D MANY R2C', maxerr, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft)
        if (ierr /= CUFFT_SUCCESS) stop 77
        call fftwf_destroy_plan(plan_fftw)
        deallocate(buf2_cufft, buf2_fftw)
    end subroutine test_fftw_vs_cufft_2d_many

    subroutine test_fftw_vs_cufft_1d_many_c2r
        integer, parameter :: n1 = 16
        integer, parameter :: n1c = n1/2 + 1
        integer, parameter :: n1pad = 2*n1c
        integer, parameter :: nbatch = 3
        real(sp), parameter :: tol_cmp = 1.0e-4_sp
        integer(c_int) :: plan_cufft_inv
        type(c_ptr)    :: plan_fftw_fwd, plan_fftw_inv, p_r, p_c
        integer(c_int) :: ierr
        integer(c_int) :: n(1), inembed(1), onembed(1)
        integer(c_int) :: rank, howmany, istride, ostride, idist, odist
        integer :: batch_idx
        real(sp), allocatable, target :: buf1_cufft(:,:), buf1_fftw(:,:)
        complex(sp), pointer :: buf1_cufft_c(:,:), buf1_fftw_c(:,:)
        real(sp) :: maxerr
        logical  :: passed
        call banner('TEST: FFTW VS CUFFT 1D MANY IN-PLACE C2R')
        allocate(buf1_cufft(n1pad,nbatch), buf1_fftw(n1pad,nbatch))
        call c_f_pointer(c_loc(buf1_cufft(1,1)), buf1_cufft_c, [n1c,nbatch])
        call c_f_pointer(c_loc(buf1_fftw(1,1)),  buf1_fftw_c,  [n1c,nbatch])
        buf1_fftw = 0.0_sp
        do batch_idx = 1, nbatch
            call init_dummy_1d(buf1_fftw(1:n1,batch_idx), n1, 0.1_sp*real(batch_idx-1,sp))
        end do
        ! Plan
        rank       = 1_c_int
        n(1)       = int(n1, c_int)
        howmany    = int(nbatch, c_int)
        istride    = 1_c_int
        ostride    = 1_c_int
        ! Forward many R2C on FFTW side to generate physically consistent C2R inputs.
        inembed(1) = int(n1, c_int)
        onembed(1) = int(n1c, c_int)
        idist      = int(n1pad, c_int)
        odist      = int(n1c, c_int)
        plan_fftw_fwd = fftwf_plan_many_dft_r2c(rank, n, howmany, buf1_fftw, inembed, istride, idist, buf1_fftw_c, onembed, ostride, odist, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_fwd)) stop 78
        call fftwf_execute_dft_r2c(plan_fftw_fwd, buf1_fftw, buf1_fftw_c)
        call fftwf_destroy_plan(plan_fftw_fwd)
        ! Copy generated complex values the CUFFT in-place buffer.
        buf1_cufft = 0.0_sp
        buf1_cufft_c = buf1_fftw_c
        ! Build many C2R plans (in-place layout).
        inembed(1) = int(n1c, c_int)
        onembed(1) = int(n1, c_int)
        idist      = int(n1c, c_int)
        odist      = int(n1pad, c_int)
        plan_fftw_inv = fftwf_plan_many_dft_c2r(rank, n, howmany, buf1_fftw_c, inembed, istride, idist, buf1_fftw, onembed, ostride, odist, FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_inv)) stop 79

        ierr = cufftPlanMany(plan_cufft_inv, rank, n, inembed, istride, idist, onembed, ostride, odist, CUFFT_C2R, howmany)
        if (ierr /= CUFFT_SUCCESS) stop 80

        call fftwf_execute_dft_c2r(plan_fftw_inv, buf1_fftw_c, buf1_fftw)
        !$omp target data map(tofrom:buf1_cufft)
          p_r = omp_get_mapped_ptr(c_loc(buf1_cufft(1,1)), 0)
          p_c = omp_get_mapped_ptr(c_loc(buf1_cufft(1,1)), 0)
          ierr = cufftExecC2R(plan_cufft_inv, p_c, p_r)
          if (ierr /= CUFFT_SUCCESS) stop 81
        !$omp end target data

        maxerr = maxval(abs(buf1_cufft(1:n1,:) - buf1_fftw(1:n1,:))) / max(1.0_sp, maxval(abs(buf1_fftw(1:n1,:))))
        print *, 'FFTW/CUFFT 1D many relative C2R output error = ', maxerr
        call report_fftw('FFTW VS CUFFT 1D MANY C2R', maxerr, tol_cmp, passed)
        ierr = cufftDestroy(plan_cufft_inv)
        if (ierr /= CUFFT_SUCCESS) stop 83
        call fftwf_destroy_plan(plan_fftw_inv)
        deallocate(buf1_cufft, buf1_fftw)
    end subroutine test_fftw_vs_cufft_1d_many_c2r

    subroutine test_fftw_vs_cufft_movie
        integer,            parameter :: nx      = 4096
        integer,            parameter :: ny      = 4096
        integer,            parameter :: nxsc    = 2048
        integer,            parameter :: nysc    = 2048
        integer,            parameter :: nframes = 50
        integer,            parameter :: n2cx    = nx/2 + 1
        integer,            parameter :: n2padx  = 2*n2cx
        real(sp),           parameter :: tol_cmp = 1.0e-4_sp
        real(sp), allocatable, target :: buf2_cufft(:,:,:), buf2_fftw(:,:,:)
        complex(sp),          pointer :: buf2_cufft_c(:,:,:), buf2_fftw_c(:,:,:)
        type(c_ptr)    :: plan_fftw_r2c, p_r, p_c
        integer(c_int) :: plan_cu_r2c
        integer(c_int) :: ierr
        integer        :: i
        call banner('BENCH: FFTW VS CUFFT 2D R2C 50 FRAMES')
        allocate(buf2_cufft(n2padx,ny,nframes), buf2_fftw(n2padx,ny,nframes))
        call c_f_pointer(c_loc(buf2_cufft(1,1,1)), buf2_cufft_c, [n2cx,ny,nframes])
        call c_f_pointer(c_loc(buf2_fftw(1,1,1)),  buf2_fftw_c,  [n2cx,ny,nframes])
        buf2_cufft = 0.0_sp
        buf2_fftw  = 0.0_sp
        do i = 1, nframes
            call init_dummy_2d(buf2_cufft(1:nx,1:ny,i), nx, ny, 0.05_sp*real(i-1,sp))
            buf2_cufft(1:nx,1:ny,i) = buf2_cufft(1:nx,1:ny,i) + real(i-1) / real(nframes-1, sp)
        end do
        buf2_fftw(1:nx,1:ny,:) = buf2_cufft(1:nx,1:ny,:)
        ! Plans
        ierr = cufftPlan2d(plan_cu_r2c, ny, nx, CUFFT_R2C)
        if (ierr /= CUFFT_SUCCESS) stop 73
        ierr = fftwf_init_threads()
        call fftwf_plan_with_nthreads(1)
        plan_fftw_r2c = fftwf_plan_dft_r2c_2d(int(ny,c_int), int(nx,c_int), buf2_fftw(:,:,1), buf2_fftw_c(:,:,1), FFTW_ESTIMATE)
        if (.not. c_associated(plan_fftw_r2c)) stop 74

        t = tic()
        !$omp parallel do private(i)
        do i = 1,nframes
            call fftwf_execute_dft_r2c(plan_fftw_r2c, buf2_fftw(:,:,i), buf2_fftw_c(:,:,i))
        end do
        !$omp end parallel do
        rt_cpu = toc(t)
        t = tic()
        !$omp target data map(to:buf2_cufft)
        do i = 1,nframes
            p_r = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1,i)), 0)
            p_c = omp_get_mapped_ptr(c_loc(buf2_cufft(1,1,i)), 0)
            ierr = cufftExecR2C(plan_cu_r2c, p_r, p_c)
            if (ierr /= CUFFT_SUCCESS) stop 75
        enddo
        !$omp end target data
        rt_gpu = toc(t)
        print *, '  FFTW 2D many R2C total time     = ', rt_cpu, ' seconds'
        print *, '  CUFFT 2D many R2C total time    = ', rt_gpu, ' seconds'
        ierr = cufftDestroy(plan_cu_r2c)
        if (ierr /= CUFFT_SUCCESS) stop 77
        call fftwf_destroy_plan(plan_fftw_r2c)
        deallocate(buf2_cufft, buf2_fftw)
    end subroutine test_fftw_vs_cufft_movie

    subroutine test_cublas_sgemm
        integer, parameter :: m = 4
        integer, parameter :: n = 5
        integer, parameter :: k = 3
        real(c_float), allocatable, target :: a(:,:), b(:,:), c(:,:), cref(:,:)
        real(c_float) :: alpha, beta, maxerr
        type(c_ptr) :: handle, d_a, d_b, d_c
        integer(c_int) :: ierr
        integer :: i, j, id
        call banner('TEST: CUBLAS SGEMM')
        id = p%device
        allocate(a(m,k), b(k,n), c(m,n), cref(m,n))
        do j = 1, k
            do i = 1, m
                a(i,j) = real(10*i + j, c_float) / 7.0_c_float
            end do
        end do
        do j = 1, n
            do i = 1, k
                b(i,j) = real(i - 2*j, c_float) / 5.0_c_float
            end do
        end do
        c     = 0.0_c_float
        cref  = matmul(a, b)
        alpha = 1.0_c_float
        beta  = 0.0_c_float
        handle = c_null_ptr
        ierr = cublasCreate(handle)
        if( ierr /= CUBLAS_STATUS_SUCCESS ) stop 'Fatal error: cublasCreate failed'
        if( .not. c_associated(handle) ) stop 'Fatal error: cublasCreate returned null handle'
        !$omp target data map(to:a,b,id) map(from:c)
        d_a = omp_get_mapped_ptr(c_loc(a(1,1)), id)
        d_b = omp_get_mapped_ptr(c_loc(b(1,1)), id)
        d_c = omp_get_mapped_ptr(c_loc(c(1,1)), id)
        if( .not. c_associated(d_a) ) stop 'Fatal error: CUBLAS SGEMM A is not mapped'
        if( .not. c_associated(d_b) ) stop 'Fatal error: CUBLAS SGEMM B is not mapped'
        if( .not. c_associated(d_c) ) stop 'Fatal error: CUBLAS SGEMM C is not mapped'
        ierr = cublasSgemm(handle, CUBLAS_OP_N, CUBLAS_OP_N, int(m,c_int), int(n,c_int), int(k,c_int),&
            &alpha, d_a, int(m,c_int), d_b, int(k,c_int), beta, d_c, int(m,c_int))
        if( ierr /= CUBLAS_STATUS_SUCCESS ) stop 'Fatal error: cublasSgemm failed'
        ierr = cudaDeviceSynchronize()
        if( ierr /= CUDA_SUCCESS ) stop 'Fatal error: cudaDeviceSynchronize failed after cublasSgemm'
        !$omp end target data
        maxerr = maxval(abs(c - cref))
        print *, 'CUBLAS SGEMM max error = ', maxerr
        if( maxerr > 1.0e-5_c_float ) stop 'Fatal error: CUBLAS SGEMM accuracy test failed'
        ierr = cublasDestroy(handle)
        if( ierr /= CUBLAS_STATUS_SUCCESS ) stop 'Fatal error: cublasDestroy failed'
        deallocate(a, b, c, cref)
        print *, 'PASS: CUBLAS SGEMM test successful'
    end subroutine test_cublas_sgemm

    subroutine test_kb
      use simple_kbinterpol, only: apod_device
      integer, parameter :: n = 15
      type(kbinterpol) :: kb
      real,allocatable :: x(:),y(:)
      real    :: z, diff
      integer :: i
      call banner('TEST: KB INTERPOLATION')
      allocate(x(n),y(n),source=0.0)
      kb = kbinterpol(KBWINSZ,KBALPHA)
      do i = 1,n
          x(i) = real(i-1)/real(n-1)
      end do
      !$omp target teams map(to:x) map(to:kb) map(from:y)
      !$omp loop
      do i = 1, n
          y(i) = apod_device(kb, x(i))
      enddo
      !$omp end loop
      !$omp end target teams
      do i = 1, n
          z = kb%apod(x(i))
          write(*,'(A,I0,A,F9.6,A,F9.6,A,F9.6)')'  ',i,'  x = ', x(i),'  gpu apod(x) = ', y(i),'  cpu apod(x) = ', z
          diff = abs(y(i) - z) 
      end do
      diff = diff / real(n)
      if( diff > 1e-6 ) then
          stop 'Fatal error: KB interpolation test failed, results do not match'
      endif
      print *, 'PASS: KB interpolation test successful, average difference = ', diff
    end subroutine test_kb

    subroutine test_pointers
        integer,  parameter :: nx      = 64
        integer,  parameter :: ny      = nx
        integer,  parameter :: nz      = nx
        real,    allocatable :: reven(:,:,:), rodd(:,:,:)
        complex, allocatable :: ceven(:,:,:), codd(:,:,:)
        complex :: cavg_even, cavg_odd
        real    :: avg_even, avg_odd
        logical :: flag(nx), passed
        integer :: i
        call banner('TEST: POINTER')
        allocate(reven(nx,ny,nz), rodd(nx,ny,nz), source=0.)
        allocate(ceven(nx,ny,nz), codd(nx,ny,nz), source=CMPLX_ZERO)
        flag(1:nx:2) = .true.
        flag(2:nx:2) = .false.
        call pointer_kernel(nx,ny,nz, reven,rodd, ceven,codd, flag)
        passed = .true.
        do i = 1, nx
            avg_even  = sum(reven(i,:,:)) / real(ny*nz)
            avg_odd   = sum(rodd(i,:,:))  / real(ny*nz)
            cavg_even = sum(ceven(i,:,:)) / real(ny*nz)
            cavg_odd  = sum(codd(i,:,:))  / real(ny*nz)
            if( flag(i) ) then
                if( .not.is_equal(avg_even, real(i)) ) then
                    print *, 'FAIL: pointer test failed for avg_even at i = ', i, ' value = ', avg_even
                    passed = .false.
                endif
                if( .not.is_equal(avg_odd, 0.) ) then
                    print *, 'FAIL: pointer test failed for avg_odd at i = ', i, ' value = ', avg_odd
                    passed = .false.
                endif
                if( .not.is_equal(imag(cavg_even), real(i)) ) then
                    print *, 'FAIL: pointer test failed for cavg_even at i = ', i, ' value = ', cavg_even
                    passed = .false.
                endif
                if( .not.is_equal(imag(cavg_odd), 0.) ) then
                    print *, 'FAIL: pointer test failed for cavg_odd at i = ', i, ' value = ', cavg_odd
                    passed = .false.
                endif
            else
                if( .not.is_equal(avg_even, 0.) ) then
                    print *, 'FAIL: pointer test failed for avg_even at i = ', i, ' value = ', avg_even
                    passed = .false.
                endif
                if( .not.is_equal(avg_odd, real(i)) ) then
                    print *, 'FAIL: pointer test failed for avg_odd at i = ', i, ' value = ', avg_odd
                    passed = .false.
                endif
                if( .not.is_equal(imag(cavg_even), 0.) ) then
                    print *, 'FAIL: pointer test failed for cavg_even at i = ', i, ' value = ', cavg_even
                    passed = .false.
                endif
                if( .not.is_equal(imag(cavg_odd), real(i)) ) then
                    print *, 'FAIL: pointer test failed for cavg_odd at i = ', i, ' value = ', cavg_odd
                    passed = .false.
                endif
            endif
            if( .not.is_equal(real(cavg_even), 0.) ) then
                print *, 'FAIL: pointer test failed for cavg_even at i = ', i, ' value = ', cavg_even
                passed = .false.
            endif
            if( .not.is_equal(real(cavg_odd), 0.) ) then
                print *, 'FAIL: pointer test failed for cavg_odd at i = ', i, ' value = ', cavg_odd
                passed = .false.
            endif
        enddo
        if( passed ) then
            print *, 'PASS: pointer test successful'
        else
            stop 'FAIL: pointer test failed'
        endif
    end subroutine test_pointers

    subroutine pointer_kernel(nx,ny,nz, re, ro, ce, co, flag)
        integer,         intent(in)    :: nx, ny, nz
        real,    target, intent(inout) :: re(nx,ny,nz), ro(nx,ny,nz)
        complex, target, intent(inout) :: ce(nx,ny,nz), co(nx,ny,nz)
        logical,         intent(in)    :: flag(nx)
        real,    pointer :: p_r(:,:,:)
        complex, pointer :: p_c(:,:,:)
        integer :: i,j,k
        !$omp target enter data map(to: re, ro, ce, co)
        !$omp target teams distribute parallel do&
        !$omp& map(to: flag) default(shared) private(i,j,k,p_r,p_c)
        do i = 1,nx
            ! pointer to e/o
            if( flag(i) ) then
                p_r => re
                p_c => ce
            else
                p_r => ro
                p_c => co
            endif
            ! dummy
            do j = 1, ny
                do k = 1, nz
                    p_r(i,j,k) = real(i)
                    p_c(i,j,k) = cmplx(0., real(i))
                end do
            end do
        enddo
        !$omp end target teams distribute parallel do
        !$omp target exit data map(from: re, ro, ce, co)
    end subroutine pointer_kernel

#endif
end program simple_test_openmp_offload
