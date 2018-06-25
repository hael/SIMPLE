! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module simple_qsort
implicit none
private
public :: qsort, quicksort_m
private :: partition

contains

recursive subroutine qsort(A)
  real, intent(in out), dimension(:) :: A
  integer(8) :: iq

  if(size(A) > 1) then
     call partition(A, iq)
     call qsort(A(:iq-1))
     call qsort(A(iq:))
  endif
end subroutine

subroutine partition(A, marker)
  real, intent(in out), dimension(:) :: A
  integer(8), intent(out) :: marker
  integer(8) :: i, j
  real :: temp
  real :: x      ! pivot point
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do
end subroutine partition




recursive subroutine quicksort_m(a, lo0, hi)
    real, intent(inout) :: a(:)
    integer(8), intent(in) :: lo0, hi
    integer(8) :: i, j, lo
    real :: pivot, t
    lo = lo0
    i = lo
    j = hi
    do while (i < hi)
        pivot = a((lo+hi)/2)
        do while (i <= j)
            do while (a(i) < pivot)
                i = i + 1
            end do
            do while (a(j) > pivot)
                j = j - 1
            end do
            if (i <= j) then
                t = a(i)
                a(i) = a(j)
                a(j) = t
                i = i + 1
                j = j - 1
            end if
        end do
        if (lo < j) call quicksort_m(a, lo, j)
        lo = i
        j = hi
    end do
end subroutine quicksort_m
end module simple_qsort
