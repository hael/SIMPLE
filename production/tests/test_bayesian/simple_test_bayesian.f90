program simple_test_bayesian
    include 'simple_lib.f08'
    use simple_opt_bayesian
    implicit none
    
    ! unit tests for 'compute_count_for_edges'
    integer :: pop(3, 4), indexes(3), counts(2**3)
    indexes  = (/ 0, 2, 3 /)
    pop(1,:) = (/ 1, 0, 1, 1 /)
    pop(2,:) = (/ 1, 0, 1, 0 /)
    pop(3,:) = (/ 1, 0, 0, 0 /)
    counts   = 0
    call compute_count_for_edges(3, pop, 2, indexes, 3, counts)

    write(*, *) 'counts = ', counts
end program simple_test_bayesian