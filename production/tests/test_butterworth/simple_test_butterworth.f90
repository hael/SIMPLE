program simple_test_butterworth
    include 'simple_lib.f08'
    use simple_butterworth
    use simple_optimizer,   only: optimizer
    use simple_opt_factory, only: opt_factory
    use simple_opt_spec,    only: opt_spec
    implicit none
    
    integer :: n = 8
    real    :: s = 1., fc = 1.
    real    :: val(2)

    integer, parameter :: ndim = 1, NRESTARTS = 1
    
    procedure(fun_butterworth),     pointer   :: costfun_ptr      !< pointer 2 cost function
    class(optimizer),               pointer   :: opt_ptr=>null()  ! the generic optimizer object
    character(len=*),               parameter :: target_img_name  = 'TargetImg.mrc'
    character(len=*),               parameter :: obj_img_name     = 'NoisyObj1.mrc'
    integer,                        parameter :: box    = 202
    real,                           parameter :: smpd   = 1.275
    type(opt_factory)   :: ofac                           ! the optimization factory object
    type(opt_spec)      :: spec                           ! the optimizer specification object
    character(len=8)    :: str_opts                       ! string descriptors for the NOPTS optimizers
    real                :: lims(ndim,2), lowest_cost
   
    call even_img%new([box,box,box], smpd)
    call even_img%read(target_img_name)

    call odd_img%new([box,box,box], smpd)
    call odd_img%read(obj_img_name)

    call ker_odd_img%new([box,box,box], smpd)
    call ker_even_img%new([box,box,box], smpd)
    call ker_der_odd_img%new([box,box,box], smpd)
    call ker_der_even_img%new([box,box,box], smpd)
    
    write(*, *) 'Butterworth unit tests:'
    ! first Butterworth function test
    val = butterworth(s, n, fc)
    if (abs(val(1) - 1./sqrt(2.)) <= 100*epsilon(val(1))) then
        write(*, *) 'Test (s = 1, fc = 1, n = 8) passed!'
    else
        write(*, *) epsilon(val(1)), abs(val(1) - 1./sqrt(2.)), 'value of Butterworth poly at freq = 1 (cut-off frequency = 1) should be 1/sqrt(2)!'
    end if

    ! second Butterworth function test
    fc  = 2.
    s   = 2.
    val = butterworth(s, n, fc)
    if (abs(val(1) - 1./sqrt(2.)) <= 100*epsilon(val(1))) then
        write(*, *) 'Test (s = 2, fc = 2, n = 8) passed!'
    else
        write(*, *) epsilon(val(1)), abs(val(1) - 1./sqrt(2.)), 'value of Butterworth poly at freq = 2 (cut-off frequency = 2) should be 1/sqrt(2)!'
    end if
    
    ! Test the optimizer
    write(*, *)
    write(*, *) 'Simple cut-off frequency optimization test:'
    costfun_ptr  => butterworth_cost
    str_opts  = 'lbfgsb'
    lims(1,1) =  1.
    lims(1,2) =  50.
    call spec%specify(str_opts, ndim, limits=lims, nrestarts=NRESTARTS) ! make optimizer spec
    call spec%set_costfun(costfun_ptr)                                  ! set pointer to costfun
    call spec%set_gcostfun(butterworth_gcost)                           ! set pointer to gradient of costfun         
    call ofac%new(spec, opt_ptr)                                        ! generate optimizer object with the factory
    spec%x    = 1.                                                      ! set initial guess
    call opt_ptr%minimize(spec, opt_ptr, lowest_cost)                   ! minimize the test function

    write(*, *) 'cost = ', lowest_cost, '; x = ', spec%x

    call opt_ptr%kill
    deallocate(opt_ptr)
end program simple_test_butterworth