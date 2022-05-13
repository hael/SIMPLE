! bayesian optimization (Pelikan's style)
module simple_opt_bayesian
include 'simple_lib.f08'
implicit none

contains
    ! simple translation of rb's compute_count_for_edges
    subroutine compute_count_for_edges(pop, indexes, counts)
        integer, intent(in)    :: indexes(:)
        integer, intent(in)    :: pop(:,:)
        integer, intent(inout) :: counts(:)
        integer :: k, l, index, val
        counts = 0
        do k = 1, size(pop, 1)
            index = 1
            do l = 1, size(indexes)
                val = indexes(size(indexes) - l + 1) + 1 
                if (pop(k, val) == 1) then
                    index = index + 2**(l-1)
                endif
            enddo
            counts(index) = counts(index) + 1
        enddo
    end subroutine compute_count_for_edges

    function fact(n) result(ret)
        integer, intent(in) :: n
        integer :: ret
        ret = 1
        if( n > 1 ) ret = n*fact(n-1)
    end function fact

    function k2equation(node, candidates, pop) result(total)
        integer, intent(in)  :: node
        integer, intent(in)  :: candidates(:)
        integer, intent(in)  :: pop(:,:)
        integer, allocatable :: indexes(:)
        integer, allocatable :: counts(:)
        real    :: total, rs
        integer :: k, a1, a2
        allocate(indexes(size(candidates)+1), counts(int(2**(size(candidates)+1))), source=0)
        indexes(1)                    = node
        indexes(2:size(candidates)+1) = candidates
        call compute_count_for_edges(pop, indexes, counts)
        write(*, *) "counts = ", counts
        total = -1.
        do k = 0, int(size(counts)/2)-1
            a1 = counts(k*2 + 1)
            a2 = counts(k*2 + 2)
            rs = 1.0*fact(a1)*fact(a2)/fact((a1+a2)+1)
            if( total < 0. )then
                total = rs
            else
                total = total*rs
            endif
        enddo
    end function k2equation
end module simple_opt_bayesian
