! serial version
! Max error at iteration         3372  was    9.99832153E-03
! Total time was    9.27200031      seconds

! naive parallel version
! Max error at iteration         3372  was    9.99832153E-03
! Total time was    8.74800014      seconds
program simple_test_laplace_solver
use simple_defs
use openacc
implicit none
integer,  parameter :: width=1000, height=1000
real(sp), parameter :: temp_tol=0.01
integer             :: i, j, iter
real(sp) :: worst_dt, start_time, stop_time
real(sp), dimension(0:height+1,0:width+1) :: temp, temp_prev
call cpu_time(start_time)
print *, 'start time ', start_time
call initialize
iter     = 1
worst_dt = 100.
!$acc data copy(temp_prev), create(temp)
do while( worst_dt > temp_tol )
    !$acc kernels
    do j=1,width
        do i=1,height
            temp(i,j) = 0.25 * (temp_prev(i+1,j) + temp_prev(i-1,j) + temp_prev(i,j+1) + temp_prev(i,j-1))
        end do
    end do
    !$acc end kernels
    worst_dt = 0.
    !$acc kernels
    do j=1,width
        do i=1,height
            worst_dt = max(abs(temp(i,j) - temp_prev(i,j)), worst_dt)
            temp_prev(i,j) = temp(i,j)
        end do
    end do
    !$acc end kernels
    if( mod(iter,100).eq.0 )then
        !$acc update host(temp)
        call track_progress
    endif
    iter = iter + 1
end do
!$acc end data
call cpu_time(stop_time)
print *, 'Max error at iteration ', iter - 1, ' was ', worst_dt
print *, 'Total time was ', stop_time - start_time, ' seconds'

contains

  subroutine initialize
      integer :: i, j
      temp_prev = 0.
      do i=0,height+1
         temp_prev(i,0) = 0.
         temp_prev(i,width+1) = (100./height) * i
      enddo
      do j=0,width+1
         temp_prev(0,j) = 0.
         temp_prev(height+1,j) = ((100.)/width) * j
      enddo
  end subroutine initialize

  subroutine track_progress
      integer :: i
      print *, '---------- Iteration number: ', iter, ' ---------------'
      do i=5,0,-1
         write (*,'("(",i4,",",i4,"):",f6.2,"  ")',advance='no') height-i,width-i,temp(height-i,width-i)
      enddo
      print *
  end subroutine track_progress

end program simple_test_laplace_solver
