! bayesian optimization (Pelikan's style)
module simple_opt_bayesian
include 'simple_lib.f08'
implicit none

contains
    ! simple translation of rb's compute_count_for_edges
    subroutine compute_count_for_edges(str_len, pop, pop_size, indexes, ind_size, counts)
        integer, intent(in)    :: str_len, pop_size, ind_size
        integer, intent(in)    :: indexes(:)
        integer, intent(in)    :: pop(:,:)
        integer, intent(inout) :: counts(:)

        integer :: k, l, index, val

        counts = 0
        do k = 1, size(pop, 1)
            index = 1
            do l = 1, size(indexes)
                val = indexes(size(indexes) - l + 1) + 1
                write(*, *) 'temp = ', pop(k, val) 
                if (pop(k, val) == 1) then
                    index = index + 2**(l-1)
                endif
            enddo
            counts(index) = counts(index) + 1
        enddo
    end subroutine compute_count_for_edges
end module simple_opt_bayesian
