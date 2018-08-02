
module simple_sort
    implicit none
    private

    !! Quicksort
    public :: qsortf
    private :: partition
    !! Quicksort microbench 
    public :: quicksort_m_sp, quicksort_m_dp, quicksort_m_int
    !! Bitonic sort
    public :: bitSort, bitonic_sort   
    !! Radix sort
    public :: radixsort  
    logical, parameter :: noassert = .true.
contains

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

    recursive subroutine qsortf(A)
        real, intent(in out), dimension(:) :: A
        integer(8) :: iq

        if(size(A) > 1) then
            call partition(A, iq)
            call qsortf(A(:iq-1))
            call qsortf(A(iq:))
        endif
    end subroutine qsortf

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


    recursive subroutine quicksort_m_dp(a, lo0, hi)
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
            if (lo < j) call quicksort_m_dp(a, lo, j)
            lo = i
            j = hi
        end do
    end subroutine quicksort_m_dp

    recursive subroutine quicksort_m_sp(a, lo0, hi)
        real, intent(inout) :: a(:)
        integer(4), intent(in) :: lo0, hi
        integer(4) :: i, j, lo
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
            if (lo < j) call quicksort_m_sp(a, lo, j)
            lo = i
            j = hi
        end do
    end subroutine quicksort_m_sp

    recursive subroutine quicksort_m_int(a, lo0, hi)
        integer, intent(inout) :: a(:)
        integer(4), intent(in) :: lo0, hi
        integer(4) :: i, j, lo
        integer :: pivot, t
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
            if (lo < j) call quicksort_m_int(a, lo, j)
            lo = i
            j = hi
        end do
    end subroutine quicksort_m_int



    pure subroutine assertion(cond)
        implicit none
        logical, intent(in) :: cond
        real, volatile :: r 

        if ( noassert ) return   !<=== Arjen M special note: compiler will detect .true. and eliminate assertion function
        r = 1.0
        if (.not. cond) r = r / 0.0
    end subroutine assertion

    pure subroutine bitonic_kernel(up, a, p,  q) 
        real,    intent(inout):: a(:)
        integer, intent(in)  :: p,q
        integer d, i,id
        real tmp
        logical, intent(in):: up

        d=ishft(1,-(p-q))

        do i=1, size(a)
            ! up = iand(ishft((i-1) , p) , 2) == 0
            id=ior(i,d)
            if ((iand(i-1, d) == 0) .and. ((a(i) > a(id)) .eqv. up)) then
                tmp = a(i); a(i) = a(id); a(id) = tmp;
            endif
        enddo
    end subroutine bitonic_kernel

    recursive subroutine bitSort(up, aN,  a) 
        real,    intent(inout):: a(:)
        integer, intent(in)   :: aN
        logical, intent(inout):: up
        integer i,j,logn
        logn=int(log10(real(size(a)))/log10(2.0))
        call assertion((aN == size(a)) .and. (aN == ishft(1,-logn)))
        !        assert a.length == 1 << logn;
        do  i=1,logn
            do j=1, i
                call bitonic_kernel(up, a, i, j)
            end do
        end do
    end subroutine bitSort

    subroutine bitonic_sort(up, x,p1,p2)
        real,intent(inout) :: x(:)
        logical, intent(inout):: up
        integer, optional :: p1,p2
        logical s1,s2
        if(.not.present(p1))then
            p1=1
            p2=size(x)
        endif
        if (p2-p1 > 1)then
            s1=.true.;s2=.false.
            call bitonic_sort(s1, x, p1, p1-1+(p2-p1)/2)
            call bitonic_sort(s2, x, p1+(p2-p1)/2, p2)
            call bitonic_merge(up, x, p1, p2)
        end if
    end subroutine bitonic_sort
    subroutine bitonic_merge(up, x,p1, p2)
        logical, intent(inout) :: up
        real, intent(inout):: x(:)
        integer, intent(in) :: p1,p2
        ! assume input x is bitonic, and sorted list is returned 
        if (p2-p1 > 1) then
            call bitonic_compare(up, x,p1,p2)
            call bitonic_merge(up, x, p1, p1+(p2-p1)/2)
            call bitonic_merge(up, x, p1+(p2-p1)/2,p2)
        end if
    end subroutine bitonic_merge
    subroutine  bitonic_compare(up, x,p1,p2)
        logical, intent(in):: up
        real, intent(inout) :: x(:)
        integer, intent(in) :: p1,p2
        real::tmp
        integer i, dist
        dist = p1 + (p2-p1)/2
        do i=p1,dist  
            if ((x(i) > x(i + dist)) .eqv. up)then
                tmp=x(i);
                x(i)=x(i + dist) 
                x(i + dist)=x(i)
            end if
        end do
    end subroutine bitonic_compare


    !! radix sort -- only for integers
    subroutine radixsort (ix, iw, n)      
        implicit none
        integer n
        integer :: ix(n), iw(n)
        integer :: i                        ! count bits
        integer :: ilim                     ! bits in an integer
        integer :: j                        ! count array elements
        integer :: p1old, p0old, p1, p0     ! indices to ones and zeros
        integer :: swap
        logical odd                         ! even or odd bit position

        !      if (n < 2) return      ! validate
        !
        ilim = bit_size(i)    !get the fixed number of bits
        ! alternate between putting data into iw and into ix
        p1 = n+1
        p0 = n                ! read from 1, n on first pass thru
        odd = .false.
        do i = 0, ilim-2
            p1old = p1
            p0old = p0         ! save the value from previous bit
            p1 = n+1
            p0 = 0                 ! start a fresh count for next bit

            if (odd) then
                do j = 1, p0old, +1             ! copy data from the zeros
                    if ( btest(iw(j), i) ) then
                        p1 = p1 - 1
                        ix(p1) = iw(j)
                    else
                        p0 = p0 + 1
                        ix(p0) = iw(j)
                    end if
                end do
                do j = n, p1old, -1             ! copy data from the ones
                    if ( btest(iw(j), i) ) then
                        p1 = p1 - 1
                        ix(p1) = iw(j)
                    else
                        p0 = p0 + 1
                        ix(p0) = iw(j)
                    end if
                end do

            else 
                do j = 1, p0old, +1             ! copy data from the zeros
                    if ( btest(ix(j), i) ) then
                        p1 = p1 - 1
                        iw(p1) = ix(j)
                    else
                        p0 = p0 + 1
                        iw(p0) = ix(j)
                    end if
                end do
                do j = n, p1old, -1            ! copy data from the ones
                    if ( btest(ix(j), i) ) then
                        p1 = p1 - 1
                        iw(p1) = ix(j)
                    else
                        p0 = p0 + 1
                        iw(p0) = ix(j)
                    end if
                end do
            end if  ! even or odd i

            odd = .not. odd
        end do  ! next i
        !        the sign bit
        p1old = p1
        p0old = p0
        p1 = n+1
        p0 = 0 

        !          if sign bit is set, send to the zero end
        do j = 1, p0old, +1
            if ( btest(iw(j), ilim-1) ) then 
                p0 = p0 + 1
                ix(p0) = iw(j)
            else
                p1 = p1 - 1
                ix(p1) = iw(j)
            end if
        end do
        do j = n, p1old, -1
            if ( btest(iw(j), ilim-1) ) then
                p0 = p0 + 1
                ix(p0) = iw(j)
            else
                p1 = p1 - 1
                ix(p1) = iw(j)
            end if
        end do
        !       reverse the order of the greater value partition
        p1old = p1
        do j = n, (p1old+n)/2+1, -1
            swap = ix(j)
            ix(j) = ix(p1)
            ix(p1) = swap
            p1 = p1 + 1
        end do
        return
    end subroutine radixsort

!! https://rosettacode.org/wiki/Sorting_algorithms/Merge_sort
subroutine Merge(A,NA,B,NB,C,NC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   integer, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   integer, intent(in)     :: B(NB)
   integer, intent(in out) :: C(NC)
 
   integer :: I,J,K
 
   I = 1; J = 1; K = 1;
   do while(I <= NA .and. J <= NB)
      if (A(I) <= B(J)) then
         C(K) = A(I)
         I = I+1
      else
         C(K) = B(J)
         J = J+1
      endif
      K = K + 1
   enddo
   do while (I <= NA)
      C(K) = A(I)
      I = I + 1
      K = K + 1
   enddo
   return
 
end subroutine merge
 
recursive subroutine MergeSort(A,N,T)
 
   integer, intent(in) :: N
   integer, dimension(N), intent(in out) :: A
   integer, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort(A,NA,T)
   call MergeSort(A(NA+1),NB,T)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge(T,NA,A(NA+1),NB,A,N)
   endif
   return
 
end subroutine MergeSort

subroutine Merge_r4(A,NA,B,NB,C,NC)
 
   integer, intent(in) :: NA,NB,NC         ! Normal usage: NA+NB = NC
   real, intent(in out) :: A(NA)        ! B overlays C(NA+1:NC)
   real, intent(in)     :: B(NB)
   real, intent(in out) :: C(NC)
 
   integer :: i,j,k
 
   i = 1; j = 1; k = 1;
   do while(i <= NA .and. j <= NB)
      if (A(i) <= B(j)) then
         C(k) = A(i)
         i = i+1
      else
         C(k) = B(j)
         j = j+1
      endif
      k = k + 1
   enddo
   do while (i <= NA)
      C(k) = A(i)
      i = i + 1
      k = k + 1
   enddo
   return
 
end subroutine merge_r4
 
recursive subroutine MergeSort_r4(A,N,T)
 
   integer, intent(in) :: N
   real, dimension(N), intent(in out) :: A
   real, dimension((N+1)/2), intent (out) :: T
 
   integer :: NA,NB,V
 
   if (N < 2) return
   if (N == 2) then
      if (A(1) > A(2)) then
         V = A(1)
         A(1) = A(2)
         A(2) = V
      endif
      return
   endif      
   NA=(N+1)/2
   NB=N-NA
 
   call MergeSort_r4(A,NA,T)
   call MergeSort_r4(A(NA+1),NB,T)
 
   if (A(NA) > A(NA+1)) then
      T(1:NA)=A(1:NA)
      call Merge_r4(T,NA,A(NA+1),NB,A,N)
   endif
   return
 
end subroutine MergeSort_r4


end module simple_sort
