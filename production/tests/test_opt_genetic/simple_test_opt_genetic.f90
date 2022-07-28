program simple_test_opt_genetic
    include 'simple_lib.f08'
    use simple_opt_genetic
    implicit none
    procedure(objective_func), pointer :: obj_func => null()
    integer, allocatable :: best(:), goal(:)
    integer              :: num_bits, max_iter, pop_size, k
    real                 :: best_val, cross_rate, mut_rate
    num_bits    = 20
    max_iter    = 100
    pop_size    = 100
    cross_rate  = 0.9
    mut_rate    = 1./(num_bits)
    allocate(best(num_bits), goal(num_bits), source=0)
    call srand(time())
    do k = 1, num_bits
        if( rand() > 0.2 ) goal(k) = 1
    enddo
    print *, 'GOAL = ', goal
    obj_func => objective
    call genetic_opt(obj_func, num_bits, max_iter, pop_size, cross_rate, mut_rate, best, best_val)
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
end program simple_test_opt_genetic