! bayesian optimization (Pelikan's style)
module simple_opt_bayesian
use simple_stack
include 'simple_lib.f08'
implicit none

type bit_cost
    integer, allocatable :: bitstr(:)
    integer              :: cost
end type bit_cost

abstract interface
    function objective_func(bitstring) result(val)
        integer, intent(in)  :: bitstring(:)
        real                 :: val
    end function objective_func
end interface

interface
subroutine qsort_C(array, elem_count, elem_size, compare) bind(C,name="qsort")     
    use iso_c_binding, only: c_ptr, c_size_t, c_funptr
    implicit none
    type(c_ptr)      , value :: array       ! C-pointer to the first entry of the array
    integer(c_size_t), value :: elem_count  ! Number of elements in the array
    integer(c_size_t), value :: elem_size   ! Size of each element, according to c_sizeof()
    type(c_funptr)   , value :: compare     ! c_funptr to the user-provided comparison function
end subroutine qsort_C
end interface

contains
    ! simple translation of rb's compute_count_for_edges
    subroutine compute_count_for_edges(obj_func, population, indexes, counts)
        procedure(objective_func), pointer       :: obj_func
        integer    ,               intent(inout) :: population(:,:)
        type(stack),               intent(in)    :: indexes
        integer    ,               intent(inout) :: counts(:)
        integer :: k, l, index, val
        real    :: cur_cost
        counts = 0
        do k = 1, size(population, 1)
            index = 1
            do l = 1, indexes%size_of()
                val = indexes%get_at(indexes%size_of() - l + 1)
                cur_cost = obj_func(population(k, :))
                population(k, val) = 1 - population(k, val)
                if( obj_func(population(k, :)) < cur_cost ) index = index + 2**(l-1)
                population(k, val) = 1 - population(k, val)
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

    function k2equation(obj_func, node, candidates, pop) result(total)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(in)    :: node
        integer,                   intent(in)    :: candidates(:)
        integer,                   intent(inout) :: pop(:,:)
        type(stack)          :: indexes
        integer, allocatable :: counts(:)
        real                 :: total, rs
        integer              :: k, a1, a2
        allocate(counts(int(2**(size(candidates)+1))), source=0)
        call indexes%new()
        call indexes%push(node)
        call indexes%push(candidates)
        call compute_count_for_edges(obj_func, pop, indexes, counts)
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

    subroutine compute_gains(obj_func, node, graph, population, gains, max_in)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(in)    :: node
        real   ,                   intent(in)    :: graph(:,:)
        integer,                   intent(inout) :: population(:,:)
        real   ,                   intent(inout) :: gains(:)
        integer,         optional, intent(in)    :: max_in
        integer              :: max, k, l, k_in_cnt, node_in_cnt
        integer, allocatable :: node_in(:)
        type(stack)          :: viable
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
                gains(k)               = k2equation(obj_func, node, node_in, population)
            endif
        enddo
    end subroutine compute_gains

    subroutine construct_network(obj_func, population, prob_size, graph, max_edges_in)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(inout) :: population(:,:)
        integer,                   intent(in)    :: prob_size
        real   ,                   intent(inout) :: graph(:,:)
        integer,         optional, intent(in)    :: max_edges_in
        real    :: gains(prob_size), max
        integer :: k, l, m, from, to, max_edges
        max_edges = 3*size(population, 1)
        if( present(max_edges_in) ) max_edges = max_edges_in
        graph = 0.
        gains = 0.
        do k = 1, max_edges
            max  = -1
            from = -1
            to   = -1
            do l = 1, size(graph, 1)
                call compute_gains(obj_func, l, graph, population, gains)
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

    subroutine topological_ordering(graph, ordered)
        real   , intent(in)    :: graph(:,:)
        integer, intent(inout) :: ordered(:)
        integer, allocatable   :: in_count(:)
        integer     :: prob_size, k, l, cur_item
        type(stack) :: stk
        prob_size = size(graph, 1)
        allocate(in_count(prob_size), source=0)
        do k = 1, prob_size
            in_count(k) = 0
            do l = 1, prob_size
                if( graph(l, k) > 0. ) in_count(k) = in_count(k) + 1
            enddo
        enddo
        call stk%new()
        do k = 1, prob_size
            if( in_count(k) == 0 ) then
                call stk%push(k)
            endif
        enddo
        do k = 1, prob_size
            cur_item = stk%pop()
            do l = 1, prob_size
                if( graph(cur_item, l) > 0. )then
                    in_count(l) = in_count(l) - 1
                    if( in_count(l) <= 0 ) call stk%push(l)
                endif
            enddo
            ordered(k) = cur_item
        enddo
    end subroutine topological_ordering

    function marginal_probability(node, population) result(val)
        integer, intent(in) :: node
        integer, intent(in) :: population(:,:)
        real    :: val
        integer :: k
        val = 0.
        do k = 1, size(population, 1)
            val = val + population(k, node)
        enddo
        val = val/size(population, 1)
    end function marginal_probability

    function count_in(node, graph) result(val)
        integer, intent(in) :: node
        real   , intent(in) :: graph(:,:)
        integer             :: val, k
        val = 0
        do k = 1, size(graph, 1)
            if( graph(k, node) > 0. ) val = val + 1
        enddo
    end function count_in

    function calculate_probability(obj_func, node, bitstring, graph, ordered, population) result(val)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(in)    :: node
        integer,                   intent(inout) :: bitstring(:)
        real   ,                   intent(in)    :: graph(:,:)
        integer,                   intent(in)    :: ordered(:)
        integer,                   intent(inout) :: population(:,:)
        real                 :: val, cur_cost
        type(stack)          :: indexes
        integer              :: k, index, i1, i2, cnt
        integer, allocatable :: counts(:)
        if( count_in(node, graph) == 0 )then
            val = marginal_probability(node, population)
            return
        endif
        allocate(counts(2**size(bitstring)), source=0)
        call indexes%new()
        call indexes%push(node)
        ! add the 'in' to node
        do k = 1, size(graph, 1)
            if( graph(k, node) > 0. ) call indexes%push(k)
        enddo
        call compute_count_for_edges(obj_func, population, indexes, counts)
        index = 0
        cnt   = 0
        do k = 1, size(graph, 1)
            if( graph(ordered(k), node) > 0. )then
                cur_cost = obj_func(bitstring)
                bitstring(ordered(k)) = 1 - bitstring(ordered(k))
                if( obj_func(bitstring) < cur_cost ) index = index + 2**cnt
                bitstring(ordered(k)) = 1 - bitstring(ordered(k))
                cnt = cnt + 1
            endif
        enddo
        i1  = index + 1*2**(count_in(node, graph)) + 1
        i2  = index + 0*2**(count_in(node, graph)) + 1
        if( counts(i1) == 0 .and. counts(i2) == 0 )then
            val = 0.0
        else
            val = 1.0*counts(i1)/(counts(i1) + counts(i2))*1.0
        endif
    end function calculate_probability

    subroutine probabilistic_logic_sample(obj_func, graph, ordered, population, bitstring, seed)
        procedure(objective_func), pointer       :: obj_func
        real   ,                   intent(in)    :: graph(:,:)
        integer,                   intent(in)    :: ordered(:)
        integer,                   intent(inout) :: population(:,:)
        integer,                   intent(inout) :: bitstring(:)
        integer,         optional, intent(in)    :: seed
        integer :: k, ordered_ind
        real    :: rand_val, prob_val
        bitstring = 0
        if( present(seed) ) call srand(seed)
        do k = 1, size(graph, 1)
            ordered_ind = ordered(k)
            rand_val    = rand()
            prob_val    = calculate_probability(obj_func, ordered_ind, bitstring, graph, ordered, population) 
            if( rand_val < prob_val ) bitstring(ordered_ind) = 1
        enddo
    end subroutine probabilistic_logic_sample

    subroutine sample_from_network(obj_func, population, graph, num_samples, samples, seed)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(inout) :: population(:,:)
        real   ,                   intent(in)    :: graph(:,:)
        integer,                   intent(in)    :: num_samples
        integer,                   intent(inout) :: samples(:, :)
        integer,         optional, intent(in)    :: seed
        integer, allocatable :: ordered(:), bitstring(:)
        integer              :: k
        allocate(ordered(  size(graph, 1)), source=0)
        allocate(bitstring(size(graph, 1)), source=-1)
        call topological_ordering(graph, ordered)
        if( present(seed) )then
            call srand(seed)
        else
            call srand(time())
        endif
        do k = 1, num_samples
            call probabilistic_logic_sample(obj_func, graph, ordered, population, bitstring)
            samples(k, :) = bitstring
        enddo
    end subroutine sample_from_network

    subroutine bayesian_search(obj_func, num_bits, max_iter, pop_size, select_size, num_child, best, seed)
        procedure(objective_func), pointer       :: obj_func
        integer,                   intent(in)    :: num_bits
        integer,                   intent(in)    :: max_iter
        integer,                   intent(in)    :: pop_size
        integer,                   intent(in)    :: select_size
        integer,                   intent(in)    :: num_child
        integer,                   intent(inout) :: best(:)
        integer,         optional, intent(in)    :: seed
        integer, allocatable :: children(:, :), selected(:, :)
        integer              :: k, l, m, best_cost, ind, cnt
        real   , allocatable :: network(:, :)
        integer, parameter   :: BITS_IN_BYTE = 8
        logical              :: converged
        type(bit_cost), allocatable, target :: pop_cost(:)
        allocate(selected(select_size, num_bits), children(num_child, num_bits), source=0)
        allocate(network(num_bits, num_bits), source=0.)
        allocate(pop_cost(pop_size-select_size+num_child))
        do k = 1, pop_size-select_size+num_child
            allocate(pop_cost(k)%bitstr(num_bits), source=0)
        enddo
        if( present(seed) )then
            call srand(seed)
        else
            call srand(time())
        endif
        do k = 1, pop_size
            pop_cost(k)%bitstr = 0
            do l = 1, num_bits
                if( rand() < 0.5 ) pop_cost(k)%bitstr(l) = 1
            enddo
            pop_cost(k)%cost = obj_func(pop_cost(k)%bitstr)
        enddo
        do k = pop_size+1,pop_size-select_size+num_child
            pop_cost(k)%bitstr = 0
            pop_cost(k)%cost   = 0
        enddo
        ! In-place sort of the pop_cost
        call qsort_C( c_loc(pop_cost(1)), & 
                      size(pop_cost, kind=c_size_t), &
                      int(storage_size(pop_cost(1))/BITS_IN_BYTE, kind=c_size_t), &
                      c_funloc(cost_compare) )
        best      = pop_cost(1)%bitstr
        best_cost = pop_cost(1)%cost
        write(*, *) 'current cost = ', best_cost
        do k = 1, max_iter
            do l = 1, select_size
                selected(l, :) = pop_cost(l)%bitstr
            enddo
            call construct_network(obj_func, selected, num_bits, network)
            call sample_from_network(obj_func, selected, network, num_child, children)
            do l = 1, num_child
                ind = pop_size-select_size+l
                pop_cost(ind)%bitstr = children(l, :)
                pop_cost(ind)%cost   = obj_func(pop_cost(ind)%bitstr)
            enddo
            ! In-place sort of the pop_cost
            call qsort_C( c_loc(pop_cost(1)), & 
                            size(pop_cost, kind=c_size_t), &
                            int(storage_size(pop_cost(1))/BITS_IN_BYTE, kind=c_size_t), &
                            c_funloc(cost_compare) )
            if( pop_cost(1)%cost >= best_cost )then
                best      = pop_cost(1)%bitstr
                best_cost = pop_cost(1)%cost
            endif
            write(*, *) 'current cost = ', best_cost
            converged = .false.
            cnt = 0
            do l = 1, size(pop_cost)
                if( .not. (all(pop_cost(l)%bitstr .eq. pop_cost(1)%bitstr)) ) cnt = cnt + 1
            enddo
            if( cnt == 0 ) converged = .true.
            if( converged .or. best_cost == num_bits ) return
        enddo
    end subroutine bayesian_search

    function cost_compare(i1ptr, i2ptr) result(sgn) bind(C)
        type(c_ptr), value, intent(in) :: i1ptr, i2ptr
        type(bit_cost),     pointer    :: i1, i2
        integer(c_int) :: sgn
        call c_f_pointer(i1ptr, i1)
        call c_f_pointer(i2ptr, i2)
        ! The user defines what 'less than', 'equal', 'greater than' means by setting 'sgn'
        if(i1%cost <  i2%cost) sgn =  1_c_int
        if(i1%cost == i2%cost) sgn =  0_c_int
        if(i1%cost >  i2%cost) sgn = -1_c_int
    end function

    subroutine test_qsort_C
        type(bit_cost), target :: array(5)
        integer, parameter     :: BITS_IN_BYTE = 8
        integer :: k
        integer :: cost(5) = [4, 3, 5, 1, 2]
        ! To be sorted
        do k = 1, 5
            allocate(array(k)%bitstr(4), source=0)
            array(k)%bitstr = [1, 2, 3, k]
            array(k)%cost   = cost(k)
        enddo
        ! In-place sort of the array
        call qsort_C( c_loc(array(1)), & 
                      size(array, kind=c_size_t), &
                      int(storage_size(array(1))/BITS_IN_BYTE, kind=c_size_t), &
                      c_funloc(cost_compare) )
        do k = 1, 5
            write(*, *) array(k)%cost, array(k)%bitstr
        enddo
    end subroutine

end module simple_opt_bayesian
