program simple_test_math
    include 'simple_lib.f08'
    use simple_math,                 only: butterworth
    use simple_testfuns_butterworth, only: testfun_butterworth, ButterworthCost
    use simple_optimizer,   only: optimizer
    use simple_opt_factory, only: opt_factory
    use simple_opt_spec,    only: opt_spec
    implicit none
    
    integer :: n = 8
    real    :: s = 1., fc = 1.
    real    :: val

    integer, parameter :: ndim = 1, NRESTARTS = 1
    real    :: lims(ndim,2), lowest_cost
    procedure(testfun_butterworth), pointer :: costfun_ptr      !< pointer 2 test function
    type(opt_factory)         :: ofac                 ! the optimization factory object
    type(opt_spec)            :: spec                 ! the optimizer specification object
    class(optimizer), pointer :: opt_ptr=>null()      ! the generic optimizer object
    character(len=8)          :: str_opts             ! string descriptors for the NOPTS optimizers
    class(*), pointer         :: fun_self => null()
    
    ! first Butterworth function test
    val = butterworth(s, n, fc)
    if (abs(val - 1./sqrt(2.)) <= 10*epsilon(val)) then
        write(*, *) 'Test (s = 1, fc = 1, n = 8) passed!'
    else
        write(*, *) epsilon(val), abs(val - 1./sqrt(2.)), 'value of Butterworth poly at freq = 1 (cut-off frequency = 1) should be 1/sqrt(2)!'
    end if

    ! second Butterworth function test
    fc  = 2.
    s   = 2.
    val = butterworth(s, n, fc)
    if (abs(val - 1./sqrt(2.)) <= 10*epsilon(val)) then
        write(*, *) 'Test (s = 2, fc = 2, n = 8) passed!'
    else
        write(*, *) epsilon(val), abs(val - 1./sqrt(2.)), 'value of Butterworth poly at freq = 2 (cut-off frequency = 2) should be 1/sqrt(2)!'
    end if

    ! Test the optimizer
    costfun_ptr => ButterworthCost
    str_opts  = 'de'
    lims(1,1) = -5.
    lims(1,2) =  10.
    call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfun_ptr)                                  ! set pointer to costfun
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function

    write(*, *) lowest_cost

    call opt_ptr%kill
    deallocate(opt_ptr)
end program simple_test_math
    