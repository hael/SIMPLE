!@descr: for all parallel tests
module simple_commanders_test_parallel
use simple_commanders_api
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_coarrays
  contains
    procedure :: execute      => exec_test_coarrays
end type commander_test_coarrays

type, extends(commander_base) :: commander_test_openacc
  contains
    procedure :: execute      => exec_test_openacc
end type commander_test_openacc

type, extends(commander_base) :: commander_test_openmp
  contains
    procedure :: execute      => exec_test_openmp
end type commander_test_openmp

type, extends(commander_base) :: commander_test_simd
  contains
    procedure :: execute      => exec_test_simd
end type commander_test_simd

contains

subroutine exec_test_coarrays( self, cline )
    ! use simple_stream_api   
    ! use simple_image,   only: image
    ! use simple_imgfile, only: imgfile
     class(commander_test_coarrays),    intent(inout) :: self
     class(cmdline),                    intent(inout) :: cline
    ! integer     :: ldim(3), i, j
    ! real        :: smpd
    ! type(image) :: img
    ! integer     :: my_rank, num_procs
    ! integer     :: a[*]
    ! ! Define a coarray variable
    ! integer, codimension[*] :: counter
    ! ! Get the rank (ID) and number of processors
    ! my_rank   = this_image()
    ! num_procs = num_images()
    ! if( this_image() == 1 )then
    !     ldim = [120,120,1]
    !     call img%new(ldim, smpd)
    !     call img%square(20)
    !     call img%write(string('squares_mrc.mrc'),1)
    ! endif
    ! ! Write out a message from each rank
    ! print *, "Hello from processor    ", my_rank,"out of",num_procs
    ! ! Synchronize to ensure all 'write' executions are done
    ! sync all
    ! ! Increment the counter on the first rank
    ! if( my_rank == 1 ) counter[1] = counter[1] + 1
    ! ! Again, ensure all processes are synchronized
    ! sync all
    ! ! Print the incremented value from rank 1
    ! if( my_rank == 1 )then
    !     print *, "The counter value is now", counter[1]
    ! endif
    ! a = 0
    ! if( this_image() == 1 ) then
    !     a = 1
    !     print *, 'Image ', this_image(), 'has a value', a
    !     print *, 'Image ', this_image(), 'sending new value to image 2.'
    !     a[2] = 2 * a
    ! endif
    ! sync all
    ! if( this_image() == 2 ) then
    !     a = 1
    !     print *, 'Image ', this_image(), 'now has a value', a
    !     print *, 'Image ', this_image(), 'sending new value to image 1.'
    !     a[1] = 2 * a
    ! endif
    ! sync all
    ! print *, 'Image ', this_image(), &
    !   'sees that image 1 now has a value ', a[1]
    ! sync all
    ! call sleep(10)
    ! sync all
    call simple_end('**** SIMPLE_TEST_COARRAYS_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_coarrays

subroutine exec_test_openacc( self, cline )
    class(commander_test_openacc),    intent(inout) :: self
    class(cmdline),                   intent(inout) :: cline
!    integer           :: i, n
!    real, allocatable :: x(:), y(:)
!    real              :: start, finish
!    n = 1000000000
!    allocate(x(n), y(n))
!    !$acc kernels
!    x(:) = 1
!    y(:) = 2
!    !$acc end kernels
!    call cpu_time(start)
!    call acc_saxpy(x, y, n, 0.5)
!    call cpu_time(finish)
!    print '("OpenACC saxpy time = ",f6.3," seconds.")', (finish - start)
!    call cpu_time(start)
!    call par_saxpy(x, y, n, 0.5)
!    call cpu_time(finish)
!    print '("Parallel saxpy time = ",f6.3," seconds.")', (finish - start)
!    call cpu_time(start)
!    call seq_saxpy(x, y, n, 0.5)
!    call cpu_time(finish)
!    print '("Sequential saxpy time = ",f6.3," seconds.")', (finish - start)
!    deallocate(x, y)
!    call simple_end('**** SIMPLE_TEST_OPENACC_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_openacc

subroutine exec_test_openmp( self, cline )
    use simple_image
    class(commander_test_openmp),    intent(inout) :: self
    class(cmdline),                  intent(inout) :: cline
    integer :: vals(100), i,j, counts(10), correct_counts(10)
    ! reference, no openmp
    vals    = 0
    counts  = 0
    do i=1,100
        vals(i) = nint(real(i)/real(10)) + 1
    enddo
    do j=1,10
        do i = 1,100
            if(vals(i)==j)then
                counts(j) = counts(j)+1
            endif
        enddo
    enddo
    correct_counts = counts
    ! safe
    vals    = 0
    counts  = 0
    !$omp parallel default(shared) proc_bind(close) private(i,j)
    !$omp do schedule(static)
    do i=1,100
        vals(i) = nint(real(i)/real(10)) + 1
    enddo
    !$omp end do
    !$omp do schedule(static)
    do j=1,10
        do i = 1,100
            if(vals(i)==j)then
                counts(j) = counts(j)+1
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel
    if(all(counts==correct_counts))then
        print *,'passed scenario one'
    else
        print *,'failed scenario one'
    endif
    ! unsafe, nowait
    vals    = 0
    counts  = 0
    !$omp parallel default(shared) proc_bind(close) private(i,j)
    !$omp do schedule(static)
    do i=1,100
        vals(i) = nint(real(i)/real(10)) + 1
    enddo
    !$omp end do nowait
    !$omp do schedule(static)
    do j=1,10
        do i = 1,100
            if(vals(i)==j)then
                counts(j) = counts(j)+1
            endif
        enddo
    enddo
    !$omp end do
    !$omp end parallel

    if(all(counts==correct_counts))then
        print *,'passed scenario two'
    else
        print *,'failed scenario two'
    endif
    call simple_end('**** SIMPLE_TEST_OPENMP_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_openmp

subroutine exec_test_simd( self, cline )
    use simple_timer
    class(commander_test_simd),    intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    integer, parameter      :: N=1000000000
    real                    :: a(N), b(N), c(N), t1, t2
    integer(timer_int_kind) :: t_loop, t_loop_simd
    real(timer_int_kind)    :: rt_loop, rt_loop_simd
    integer :: i
    a = 0.
    b = 0.
    c = 0.
    ! t_loop = tic()
    ! do i=1,N
    !     a(i) = b(i) + c(i)
    ! end do
    ! rt_loop = toc(t_loop)
    ! print *, 'time(loop): ', rt_loop
    ! t_loop_simd = tic()
    ! !$omp simd
    ! do i=1,N
    !     a(i) = b(i) + c(i)
    ! end do
    ! !$omp end simd
    ! rt_loop_simd = toc(t_loop_simd)
    ! print *, 'time(loop_simd): ', rt_loop_simd
    ! print *, 'speedup with simd: ', rt_loop / rt_loop_simd
    t_loop = tic()
    do i=1,N
        t1 = func1(b(i), c(i))
        t2 = func2(b(i), c(i))
        a(i) = t1 + t2
    end do
    rt_loop = toc(t_loop)
    print *, 'time(loop): ', rt_loop

    t_loop_simd = tic()
    !$omp simd private(t1,t2)
    do i=1,N
        t1 = func1(b(i), c(i))
        t2 = func2(b(i), c(i))
        a(i) = t1 + t2
    end do
    !$omp end simd
    rt_loop_simd = toc(t_loop_simd)
    print *, 'time(loop_simd): ', rt_loop_simd
    print *, 'speedup with simd: ', rt_loop / rt_loop_simd
    call simple_end('**** SIMPLE_TEST_SIMD_WORKFLOW NORMAL STOP ****')

    contains

        real function func1( a, b )
            real, intent(in) :: a, b
            func1 = a * b
        end function func1

        real function func2( a, b )
            real, intent(in) :: a, b
            func2 = (a + b)**2.0
        end function func2

end subroutine exec_test_simd

end module simple_commanders_test_parallel
