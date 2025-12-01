program simple_test_lbfgsb
    include 'simple_lib.f08'
    use simple_optimizer,   only: optimizer
    use simple_opt_factory, only: opt_factory
    use simple_opt_spec,    only: opt_spec
    implicit none
    class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
    integer,          parameter :: NDIM=1, NRESTARTS=1
    type(opt_factory) :: ofac                           ! the optimization factory object
    type(opt_spec)    :: spec                           ! the optimizer specification object
    character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
    real              :: lims(NDIM,2), lowest_cost
    str_opts  = 'lbfgsb'
    lims(1,1) = -5.
    lims(1,2) =  5.
    call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfct)                                      ! set pointer to costfun
    call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    spec%x = 0.5
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
    write(*, *) lowest_cost, spec%x
    call opt_ptr%kill
    deallocate(opt_ptr)

  contains

    function costfct( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        r = (x(1) - 1.)**2
    end function

    subroutine gradfct( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        grad(1) = 2. * (x(1) - 1.)
    end subroutine

end program simple_test_lbfgsb
