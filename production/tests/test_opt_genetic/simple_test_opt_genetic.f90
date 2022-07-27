program simple_test_opt_genetic
    include 'simple_lib.f08'
    use simple_opt_genetic
    implicit none
    integer, allocatable :: best(:)
    integer              :: num_bits, max_iter, pop_size
    real                 :: best_val, cross_rate, mut_rate
    num_bits    = 20
    max_iter    = 100
    pop_size    = 100
    cross_rate  = 0.9
    mut_rate    = 1./(num_bits)
    allocate(best(num_bits), source=0)
    call genetic_opt(num_bits, max_iter, pop_size, cross_rate, mut_rate, best, best_val)
end program simple_test_opt_genetic