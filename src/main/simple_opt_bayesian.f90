! bayesian optimization (Pelikan's style)
module simple_opt_bayesian
use simple_stack
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
                val = indexes(size(indexes) - l + 1)
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

    function path_exist(i, j, graph) result(val)
        integer, intent(in)  :: i
        integer, intent(in)  :: j
        real   , intent(in)  :: graph(:,:)
        logical              :: val
        type(stack)          :: stk, visited
        integer              :: k, l
        call stk%new
        call visited%new
        call stk%push(i)
        do while (.not. stk%is_empty())
            if( stk%contains(j) )then
                val = .true.
                return
            endif
            k = stk%pop()
            if( visited%contains(k) ) cycle
            call visited%push(k)
            do l = 1,size(graph,2)
                if( graph(k,l) > 0. .and. .not. visited%contains(l) )then
                    call stk%push(l)
                endif
            enddo
        enddo
        val = .false.
    end function path_exist

    function can_add_edge(i, j, graph) result(val)
        integer, intent(in) :: i
        integer, intent(in) :: j
        real   , intent(in) :: graph(:,:)
        logical :: val
        val = (.not. ( graph(i,j) > 0. .or. path_exist(j,i,graph)))
    end function can_add_edge

    subroutine get_viable_parents(node, graph, viable)
        integer    , intent(in)    :: node
        real       , intent(in)    :: graph(:,:)
        type(stack), intent(inout) :: viable
        integer :: k
        do k = 1, size(graph,1)
            if( (.not. node == k) .and. can_add_edge(node, k, graph) ) call viable%push(k)
        enddo
    end subroutine get_viable_parents

    subroutine compute_gains(node, graph, population, gains, max_in)
        integer,           intent(in)    :: node
        real   ,           intent(in)    :: graph(:,:)
        integer,           intent(in)    :: population(:,:)
        real   ,           intent(inout) :: gains(:)
        integer, optional, intent(in)    :: max_in
        integer              :: max, k, l, k_in_cnt, node_in_cnt
        integer, allocatable :: node_in(:)
        type(stack) :: viable
        if( present(max_in) )then
            max = max_in
        else
            max = 2
        endif
        call viable%new()
        call get_viable_parents(node, graph, viable)
        gains = -1.
        ! get the nodes pointing to node
        node_in_cnt = 0
        do l = 1, size(graph,1)
            if( graph(l,node) > 0. ) node_in_cnt = node_in_cnt + 1
        enddo
        allocate(node_in(node_in_cnt + 1), source=0)
        node_in_cnt = 1
        do l = 1, size(graph,1)
            if( graph(l,node) > 0. )then
                node_in(node_in_cnt) = l
                node_in_cnt          = node_in_cnt + 1
            endif
        enddo
        do k = 1, size(gains)
            ! compute the nodes pointing to k
            k_in_cnt = 0
            do l = 1, size(graph,1)
                if( graph(l,k) > 0.) k_in_cnt = k_in_cnt + 1
            enddo
            if( k_in_cnt < max .and. viable%contains(k) )then
                node_in(size(node_in)) = k
                gains(k)               = k2equation(node, node_in, population)
            endif
        enddo
    end subroutine compute_gains

    subroutine construct_network(population, prob_size, graph, max_edges_in)
        integer,           intent(in)    :: population(:,:)
        integer,           intent(in)    :: prob_size
        real   ,           intent(inout) :: graph(:,:)
        integer, optional, intent(in)    :: max_edges_in
        real    :: gains(prob_size), max
        integer :: k, l, m, from, to, max_edges
        max_edges = 3*size(population, 1)
        if( .not. present(max_edges_in) ) max_edges = max_edges_in
        graph = 0.
        gains = 0.
        do k = 1, max_edges
            max  = -1
            from = -1
            to   = -1
            do l = 1, size(graph, 1)
                call compute_gains(l, graph, population, gains)
                do m = 1, size(gains)
                    if( gains(m) > max )then
                        from = l
                        to   = m
                        max  = gains(m)
                    endif
                enddo
            enddo
            if( max <= 0.0 ) return
            graph(from, to) = 1.
        enddo
    end subroutine construct_network
end module simple_opt_bayesian
