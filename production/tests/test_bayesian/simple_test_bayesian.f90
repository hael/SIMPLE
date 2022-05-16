program simple_test_bayesian
    include 'simple_lib.f08'
    use simple_opt_bayesian
    use simple_stack
    implicit none
    integer     :: population(3, 4), indexes(3), counts(2**3), candidates(1)
    real        :: total
    real        :: graph(4, 4)
    real        :: gains(4)
    ! unit tests for 'compute_count_for_edges'
    write(*, *) 'Testing counting for each edge:'
    indexes         = (/ 0, 2, 3 /)
    population(1,:) = (/ 0, 0, 1, 1 /)
    population(2,:) = (/ 0, 0, 1, 1 /)
    population(3,:) = (/ 0, 0, 0, 1 /)
    counts   = 0
    call compute_count_for_edges(population, indexes, counts)
    write(*, *) 'counts = ', counts
    write(*, *) 'Testing Dirichlet metrics:'
    candidates = (/ 1 /)
    total = k2equation(0, candidates, population)
    write(*, *) 'fact(3) = ', fact(3)
    write(*, *) 'total = ', total
    write(*, *) 'Testing path_exist of a graph:'
    graph(1,:) = (/ 0, 1, 0, 0 /)
    graph(2,:) = (/ 0, 0, 1, 0 /)
    graph(3,:) = (/ 0, 0, 0, 0 /)
    graph(4,:) = (/ 0, 0, 0, 0 /)
    write(*, *) path_exist(1, 3, graph), path_exist(1, 4, graph)
    write(*, *) 'Testing can_add_edge of a graph:'
    write(*, *) can_add_edge(3, 1, graph), can_add_edge(3, 4, graph), can_add_edge(1, 3, graph), can_add_edge(2, 3, graph)
    write(*, *) 'Testing compute_gains of a graph:'
    population(1,:) = (/ 1, 1, 1, 1 /)
    population(2,:) = (/ 0, 0, 1, 0 /)
    population(3,:) = (/ 0, 0, 1, 0 /)
    graph(1,:) = (/ 0., 1., 1., 1. /)
    graph(2,:) = (/ 0., 0., 1., 0. /)
    graph(3,:) = (/ 0., 0., 0., 0. /)
    graph(4,:) = (/ 0., 0., 0., 0. /)
    gains      = 0.
    call compute_gains(4, graph, population, gains) 
    write(*, *) 'gains = ', gains
    write(*, *) 'Testing constructing network:'
    population(1,:) = (/ 0, 0, 1, 0 /)
    population(2,:) = (/ 0, 0, 1, 0 /)
    population(3,:) = (/ 0, 0, 0, 0 /)
    graph = 0.
    call construct_network(population, 4, graph, size(population,1))
    write(*, *) graph
end program simple_test_bayesian