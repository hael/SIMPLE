!
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 10th of October 2014.
!
! Name:
! simple_sorting - Various utilities and sorting routine for both CPU and GPU
!                  for other modules.
!
! Description:
! simple_sorting provides sorting routine for arrays and data sets on both
! CPU and GPU. When on GPU depending if implemented in CUDA or OpenCL a #ifdef
! is set for now. To be used in other modules for sorting large datasets and
! data arrays.
!*******************************************************************************
!
module simple_sorting

  use simple_defs

  implicit none

#ifdef CUDA
!TODO: insert the CUDA sorting variables and macros 
!#else MPI
!TODO: insert the MPI directives and layout and macros 
#endif

  interface
     !TODO: insert the subroutines and functions and procedure to interface
  end interface

contains
!*******************************************************************************
! DESCRIPTION
! subroutine to ShellSort an array of fixed size on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine shellSort_cpu(shell_sorted, timers, ntim)
    implicit none
    !global variables
    integer                        :: ntim            !# of named timers 
    integer                        :: shell_sorted(0:*) !Index of sorted timers
    character(len=*)               :: timers(0:*)       !timers_name array
    
    !local variables
    integer                        :: increment       !Incrementer for the sorting procedure
    !counters
    integer                        :: i,j,k
    !start of the execution commands

    !Sort the timers by the shellsort method.
    !Worst case O(n^2), best case O(n*log(n)) to speed up access.
    
    increment = 1
    do
       increment = 3*increment+1
       if (increment > ntim+1) exit
    end do

    shell_sorted(ntim) = ntim

    do
       increment = increment/3
       do i = increment, ntim
          k = shell_sorted(i)
          j = i
          do
             if (timers(shell_sorted(j-increment)) <= timers(k)) exit
             shell_sorted(j) = shell_sorted(j-increment)
             j = j-increment
             if (j < increment) exit
          end do
          shell_sorted(j) = k
       end do
       if (increment <= 1) exit
    end do

    return
  end subroutine shellSort_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to quick sort an array of fixed size on CPU.
!
!*******************************************************************************
! SOURCE
  subroutine quickSort_cpu(quick_sorted, a, low, high)
    implicit none
    !global variables
    integer                        :: low, high
    double precision               :: a(*)            !Input array
    double precision               :: quick_sorted(*) !Quick sorted array
    !local variables

    !start of the execution commands
    !TODO: implement the quick sort algorithm

    return
  end subroutine quickSort_cpu
!*******************************************************************************
! DESCRIPTION
! subroutine to greet and good bye the sorter on CPU.
!
!*******************************************************************************
! SOURCE
  !Helloer of the sorting module
  subroutine hello_sorting_cpu()
    implicit none
    write(*,*)"Hello sorting on CPU"
    return
  end subroutine hello_sorting_cpu
  !goodbyer of the sorting module
  subroutine bye_sorting_cpu()
    implicit none
    write(*,*)"Good bye sorting on CPU"
    return
  end subroutine bye_sorting_cpu

end module simple_sorting
