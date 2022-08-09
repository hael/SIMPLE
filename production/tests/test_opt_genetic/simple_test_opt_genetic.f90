program simple_test_opt_genetic
    include 'simple_lib.f08'
    use simple_opt_genetic
    implicit none
    procedure(objective_func), pointer :: obj_func => null()
    integer, allocatable :: best(:), goal(:)
    integer              :: nptcls, max_iter, pop_size, k, nrefs
    real                 :: best_val, cross_rate, mut_rate, val, bounds(2)
    call srand(time())
    ! testing genetic algo with continuous function f(x) = 2 - (x-4)^2
    bounds      = [-5., 15.]
    nrefs       = 2            ! only need two bits here (0 and 1)
    nptcls      = 10
    max_iter    = 20
    pop_size    = 30
    cross_rate  = 0.9
    mut_rate    = 1./nptcls
    obj_func    => objective_cont
    allocate(best(nptcls), goal(nptcls), source=0)
    call genetic_opt(obj_func, nptcls, nrefs, max_iter, pop_size, cross_rate, mut_rate, best, best_val,&
                    &tol=1e-3)
    print *, 'result = ', to_int(best)
    ! testing genetic algo with discrete cost function
    nrefs = 5
    do k = 1, nptcls
        goal(k) = floor(rand()*nrefs)
    enddo
    print *, 'GOAL = ', goal
    obj_func => objective
    call genetic_opt(obj_func, nptcls, nrefs, max_iter, pop_size, cross_rate, mut_rate, best, best_val)
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

    function to_int(bitstring) result(val)
        integer, intent(in) :: bitstring(:)
        real    :: val
        integer :: k
        val = 0.
        do k = 1, size(bitstring)
            if( bitstring(k) == 1 ) val = val + 2**(k-1)    
        enddo
        val = bounds(1) + val*(bounds(2) - bounds(1))/(2**size(bitstring)-1.)
    end function to_int

    function objective_cont(bitstring) result(val)
        integer, intent(in) :: bitstring(:)
        real    :: val
        integer :: k
        val = 2. - (to_int(bitstring) - (-4.))**2
    end function objective_cont
end program simple_test_opt_genetic