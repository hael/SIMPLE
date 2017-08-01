program simple_test_bforce_opt
    use simple_bforce_opt, only: bforce_opt
    use simple_opt_spec,   only: opt_spec
    use simple_testfuns
    implicit none
    type(opt_spec)              :: spec
    type(bforce_opt)            :: bforce
    real                        :: lowest_cost, limits(2,2), range(2), gmin
    procedure(testfun), pointer :: pfun
    pfun = get_testfun(7, 2, gmin, range)
    limits(1,1) = range(1)
    limits(1,2) = range(2)
    limits(2,1) = range(1)
    limits(2,2) = range(2)
    call spec%specify('bforce', 2, limits=limits, stepsz=[0.01,0.01])
    call spec%set_costfun(pfun)
    call bforce%new(spec)
    call bforce%minimize(spec, lowest_cost)
    print *, 'euclid:', sqrt(sum(spec%x**2.0))
    write(*,*) 'diff glob: ', abs(lowest_cost)
end program simple_test_bforce_opt
