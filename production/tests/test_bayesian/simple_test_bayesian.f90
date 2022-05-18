program simple_test_bayesian
    include 'simple_lib.f08'
    use simple_opt_bayesian
    use simple_stack
    implicit none
    integer     :: population(3, 4), counts(2**3), candidates(1), population15(15, 4)
    real        :: res
    real        :: graph(4, 4)
    real        :: gains(4)
    integer     :: ordered(4), node, bitstring(4)
    type(stack) :: indexes
    integer     :: samples(5,4), k
    ! unit tests for 'compute_count_for_edges'
    write(*, *) 'Testing counting for each edge:'
    call indexes%new()
    call indexes%push((/ 1, 2 /))
    population(1,:) = (/ 1, 0, 1, 1 /)
    population(2,:) = (/ 1, 1, 0, 0 /)
    population(3,:) = (/ 1, 0, 0, 1 /)
    counts   = 0
    call compute_count_for_edges(population, indexes, counts)
    write(*, *) 'counts = ', counts
    if( all(counts .eq. (/ 0, 0, 2, 1, 0, 0, 0, 0 /)) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing Dirichlet metrics
    write(*, *) 'Testing Dirichlet metrics:'
    candidates = (/ 1 /)
    node       = 0
    res = k2equation(node, candidates, population)
    write(*, *) 'total = ', res
    if( res == 0.25 )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing path_exist
    write(*, *) 'Testing path_exist of a graph:'
    graph(1,:) = (/ 0, 1, 0, 0 /)
    graph(2,:) = (/ 0, 0, 1, 0 /)
    graph(3,:) = (/ 0, 0, 0, 0 /)
    graph(4,:) = (/ 0, 0, 0, 0 /)
    write(*, *) path_exist(1, 3, graph), path_exist(1, 4, graph)
    if( path_exist(1, 3, graph) .eqv. .true. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    if( path_exist(1, 4, graph) .eqv. .false. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing can_add_edge
    write(*, *) 'Testing can_add_edge of a graph:'
    write(*, *) can_add_edge(3, 1, graph), can_add_edge(3, 4, graph), can_add_edge(1, 3, graph), can_add_edge(2, 3, graph)
    if( can_add_edge(3, 1, graph) .eqv. .false. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    if( can_add_edge(3, 4, graph) .eqv. .true. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    if( can_add_edge(1, 3, graph) .eqv. .true. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    if( can_add_edge(2, 3, graph) .eqv. .false. )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing compute_gains
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
    if( all(gains .eq. (/ -1., 1./6, -1., -1. /)) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing constructing network
    write(*, *) 'Testing constructing network:'
    population(1,:) = (/ 0, 0, 1, 0 /)
    population(2,:) = (/ 0, 0, 1, 0 /)
    population(3,:) = (/ 0, 0, 0, 0 /)
    graph = 0.
    call construct_network(population, 4, graph, size(population,1))
    write(*, *) graph
    if( all(pack(graph, .true.) .eq. (/ 0., 0., 1., 0., 1., 0., 1., 0., 0., 0., 0., 0., 1., 1., 0., 0. /)) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing topological ordering of the graph
    write(*, *) 'Testing topological ordering of the graph:'
    graph(1,:) = (/ 0., 1., 0., 0. /)
    graph(2,:) = (/ 0., 0., 0., 0. /)
    graph(3,:) = (/ 1., 1., 0., 1. /)
    graph(4,:) = (/ 1., 0., 0., 0. /)
    ordered    = 0
    call topological_ordering(graph, ordered)
    write(*, *) ordered
    if( all(ordered .eq. (/ 3, 4, 1, 2 /)) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED!'
    endif
    ! Testing marginal_probability
    write(*, *) 'Testing marginal probability:'
    population(1,:) = (/ 1, 1, 1, 1 /)
    population(2,:) = (/ 0, 1, 1, 0 /)
    population(3,:) = (/ 0, 0, 1, 0 /)
    node            = 1
    res             = marginal_probability(node, population)
    if( abs(res - 1./3.) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    node            = 2
    res             = marginal_probability(node, population)
    if( abs(res - 2./3.) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    node            = 3
    res             = marginal_probability(node, population)
    if( abs(res - 3./3.) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    ! testing calculate_probability
    write(*, *) 'Testing probability calculation:'
    population(1,:) = (/ 1, 0, 1, 1 /)
    population(2,:) = (/ 1, 0, 1, 0 /)
    population(3,:) = (/ 0, 0, 0, 1 /)
    graph(1,:)      = (/ 0., 1., 1., 1. /)
    graph(2,:)      = (/ 0., 0., 1., 0. /)
    graph(3,:)      = (/ 0., 0., 0., 0. /)
    graph(4,:)      = (/ 0., 1., 0., 0. /)
    node            = 1
    bitstring       = (/ -1, -1, -1, -1 /)
    res = calculate_probability(node, bitstring, graph, population)
    if( abs(res - 2./3.) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    population(1,:) = (/ 1, 0, 1, 1 /)
    population(2,:) = (/ 1, 0, 1, 0 /)
    population(3,:) = (/ 0, 0, 0, 1 /)
    graph(1,:)      = (/ 0., 1., 1., 1. /)
    graph(2,:)      = (/ 0., 0., 1., 0. /)
    graph(3,:)      = (/ 0., 0., 0., 0. /)
    graph(4,:)      = (/ 0., 1., 0., 0. /)
    node            = 4
    bitstring       = (/ 0, -1, -1, -1 /)
    res = calculate_probability(node, bitstring, graph, population)
    if( abs(res - 1.) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    population15( 1,:) = (/ 1, 0, 1, 1 /)
    population15( 2,:) = (/ 1, 0, 1, 1 /)
    population15( 3,:) = (/ 1, 0, 1, 1 /)
    population15( 4,:) = (/ 1, 0, 1, 1 /)
    population15( 5,:) = (/ 1, 0, 1, 1 /)
    population15( 6,:) = (/ 1, 0, 1, 1 /)
    population15( 7,:) = (/ 1, 0, 1, 0 /)
    population15( 8,:) = (/ 1, 0, 1, 0 /)
    population15( 9,:) = (/ 1, 0, 1, 0 /)
    population15(10,:) = (/ 1, 0, 1, 0 /)
    population15(11,:) = (/ 1, 0, 1, 0 /)
    population15(12,:) = (/ 1, 0, 1, 0 /)
    population15(13,:) = (/ 1, 0, 1, 0 /)
    population15(14,:) = (/ 1, 0, 1, 0 /)
    population15(15,:) = (/ 1, 0, 1, 0 /)
    graph(1,:)      = (/ 0., 1., 1., 0. /)
    graph(2,:)      = (/ 0., 0., 1., 0. /)
    graph(3,:)      = (/ 0., 0., 0., 0. /)
    graph(4,:)      = (/ 1., 1., 0., 0. /)
    node            = 4
    bitstring       = (/ -1, -1, -1, -1 /)
    res = calculate_probability(node, bitstring, graph, population15)
    if( abs(res - .4) < epsilon(res) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', res
    endif
    ! Testing probabilistic_logic_sample
    write(*, *) 'Testing probabilistic_logic_sample:'
    population(1,:) = (/ 0, 1, 1, 0 /)
    population(2,:) = (/ 0, 0, 1, 1 /)
    population(3,:) = (/ 1, 0, 0, 0 /)
    graph(1,:)      = (/ 0., 1., 1., 0. /)
    graph(2,:)      = (/ 0., 0., 1., 1. /)
    graph(3,:)      = (/ 0., 0., 0., 1. /)
    graph(4,:)      = (/ 0., 0., 0., 0. /)
    bitstring       = (/ -1, -1, -1, -1 /)
    ordered         = (/ 1, 2, 3, 4 /)
    call probabilistic_logic_sample(graph, ordered, population, bitstring, 86465)
    if( all(bitstring .eq. (/ 0, 1, 1, 0 /)) )then
        write(*, *) 'PASSED!'
    else
        write(*, *) 'FAILED! result = ', bitstring
    endif
    ! Testing sample_from_network
    write(*, *) 'Testing sample_from_network:'
    population(1,:) = (/ 1, 1, 1, 1 /)
    population(2,:) = (/ 1, 1, 0, 0 /)
    population(3,:) = (/ 0, 1, 0, 0 /)
    graph(1,:)      = (/ 0., 1., 1., 0. /)
    graph(2,:)      = (/ 0., 0., 0., 1. /)
    graph(3,:)      = (/ 0., 1., 0., 1. /)
    graph(4,:)      = (/ 0., 0., 0., 0. /)
    call sample_from_network(population, graph, 5, samples)
    do k = 1, size(samples,1)
        write(*, *) samples(k, :)
    enddo
    write(*, *) 'BARELY PASSED! These are random!'
end program simple_test_bayesian