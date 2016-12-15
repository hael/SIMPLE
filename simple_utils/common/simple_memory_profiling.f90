! 
!*******************************************************************************
!     Author: Frederic D.R. Bonnet date: 10th of July 2016.
!
! Name:
! simple_memory_profiling - Module to profile memory
!
!*******************************************************************************
!
module simple_memory_profiling
  use simple_defs
  implicit none

  public

  !public data structure for the memstat
  type, public :: memstat
     character(len=36) :: routine,array
     integer(kind=8) :: memory,peak
  end type memstat
  !public data structure for the state of the memory
  type, public :: memory_state
     integer :: memalloc     !Number of recorded allocations
     integer :: memdealloc   !Number of recorded deallocations
     type(memstat) :: memloc !Memory state in local routine
     type(memstat) :: memtot !Memory global state in profiler
  end type memory_state

  real(sp) :: memorylimit = 0.e0 !Limit of the memory allowed, in Gb

contains
  !TODO: insert the handlers for the profiling here

  !Set a memory limit
  subroutine file_set_memory_limit(limit)
    real, intent(in) :: limit
    memorylimit = limit
  end subroutine file_set_memory_limit

end module simple_memory_profiling

  
