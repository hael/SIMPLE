program simple_test_bayesian
    include 'simple_lib.f08'
    use simple_opt_bayesian
    implicit none
    
    ! unit tests for 'compute_count_for_edges'
    integer :: pop(3, 4), indexes(3), counts(2**3), candidates(1)
    real    :: total
    indexes  = (/ 0, 2, 3 /)
    pop(1,:) = (/ 0, 0, 1, 1 /)
    pop(2,:) = (/ 0, 0, 1, 1 /)
    pop(3,:) = (/ 0, 0, 0, 1 /)
    counts   = 0
    call compute_count_for_edges(pop, indexes, counts)
    write(*, *) 'counts = ', counts
    candidates = (/ 1 /)
    total = k2equation(0, candidates, pop)
    write(*, *) 'fact(3) = ', fact(3)
    write(*, *) 'total = ', total
end program simple_test_bayesian