program simple_test_ori_bayesian
    include 'simple_lib.f08'
    use simple_opt_bayesian
    use simple_stack
    implicit none
    procedure(objective_func), pointer :: obj_func => null()
    integer, allocatable :: best(:), goal(:)
    integer              :: num_bits, max_iter, pop_size, select_size, num_child, k
    num_bits    = 20
    max_iter    = 100
    pop_size    = 100
    select_size = 25
    num_child   = 35
    allocate(best(num_bits), goal(num_bits), source=0)
    call srand(time())
    do k = 1, num_bits
        if( rand() > 0.2 ) goal(k) = 1
    enddo
    print *, 'GOAL = ', goal
    obj_func => objective
    call bayesian_search(obj_func, num_bits, max_iter, pop_size, select_size, num_child, best)
    if( all(best .eq. goal) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    print *, 'result = ', best
contains
    function objective(bitstring) result(val)
        integer, intent(in) :: bitstring(:)
        real    :: val
        integer :: k
        val = 0
        do k = 1, size(goal)
            if( bitstring(k) == goal(k) ) val = val + 1
        enddo
    end function objective
end program simple_test_ori_bayesian