program simple_test_lbfgsb_cosine
    include 'simple_lib.f08'
    use simple_optimizer,            only: optimizer
    use simple_opt_factory,          only: opt_factory
    use simple_opt_spec,             only: opt_spec
    implicit none
    class(optimizer), pointer   :: opt_ptr=>null()      ! the generic optimizer object
    integer,          parameter :: NDIM=2, NRESTARTS=10
    type(opt_factory) :: ofac                           ! the optimization factory object
    type(opt_spec)    :: spec                           ! the optimizer specification object
    character(len=8)  :: str_opts                       ! string descriptors for the NOPTS optimizers
    real              :: lims(NDIM,2), lowest_cost, y_norm(2), y(2)
    str_opts  = 'lbfgsb'
    lims(:,1) = -100.
    lims(:,2) =  100.
    y         = [.75, .25]
    y_norm    = y / sqrt(sum(y**2))
    call spec%specify(str_opts, NDIM, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfct)                                      ! set pointer to costfun
    call spec%set_gcostfun(gradfct)                                     ! set pointer to gradient of costfun
    spec%opt_callback => callback_fct
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    spec%x = [0.5, 0.75]
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function
    spec%x = spec%x * sqrt(sum(y**2)) / sqrt(sum(spec%x**2))
    write(*, *) lowest_cost, spec%x, y
    call opt_ptr%kill
    deallocate(opt_ptr)

  contains

    function costfct( fun_self, x, d ) result( r )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(in)    :: x(d)
        real                    :: r
        real :: x_norm(d)
        x_norm = x / sqrt(sum(x**2))
        r      = acos(sum(x_norm * y_norm)) * 180. / PI
        print *, 'cost = ', r
    end function

    subroutine gradfct( fun_self, x, grad, d )
        class(*), intent(inout) :: fun_self
        integer,  intent(in)    :: d
        real,     intent(inout) :: x(d)
        real,     intent(out)   :: grad(d)
        real :: abs_x, x_norm(d)
        abs_x  = sqrt(sum(x**2))
        x_norm = x / abs_x
        grad   = - 180. / PI * y_norm * (abs_x**2 - x**2) / abs_x**3 / sqrt(1. - sum(x_norm * y_norm)**2)
        grad   = grad * sqrt(sum(y**2)) / sqrt(sum(grad**2))
    end subroutine

    subroutine callback_fct( self )
        class(*), intent(inout) :: self
        spec%x = spec%x * sqrt(sum(y**2)) / sqrt(sum(spec%x**2))
    end subroutine callback_fct

end program simple_test_lbfgsb_cosine