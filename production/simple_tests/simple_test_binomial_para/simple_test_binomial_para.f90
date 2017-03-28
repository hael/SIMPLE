program test_binomial_para
use simple_strings, only: int2str
use simple_defs
implicit none

integer, parameter                 :: NPARTS=10 
character(len=STDLEN), allocatable :: labels(:)
integer :: i
allocate(labels(NPARTS))
do i=1,NPARTS
    labels(i) = int2str(i)
end do
call fibonacci_partition(labels)

contains

    recursive subroutine fibonacci_partition( labels )
        character(len=STDLEN), allocatable :: labels(:)
        character(len=STDLEN), allocatable :: reduced_labels(:)
        integer :: n, nred, cnt
        n = size(labels)
        if( n == 1 ) return
        if( mod(n,2) == 0. )then
            nred = n/2
        else
            nred = (n-1)/2 + 1
        endif
        allocate(reduced_labels(nred))
        cnt = 0
        do i=1,n,2
            cnt = cnt+1
            if( i+1 > n )then
                reduced_labels(cnt) = trim(labels(i))
            else
                reduced_labels(cnt) = trim(labels(i))//trim(labels(i+1))
            endif
            write(*,'(a,a)',advance='no') trim(reduced_labels(cnt)), ' '
        end do
        write(*,*) ''
        ! labels now become reduced labels
        deallocate(labels)
        allocate(labels(nred), source=reduced_labels)
        deallocate(reduced_labels)
        ! recursion
        call fibonacci_partition(labels)
    end subroutine fibonacci_partition

end program test_binomial_para
