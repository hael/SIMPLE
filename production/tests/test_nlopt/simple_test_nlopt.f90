program simple_test_nlopt
    include 'simple_lib.f08'
    use nlopt_wrap,                 only : nlopt_opt, nlopt_func, create, destroy
    use nlopt_enum,                 only : NLOPT_SUCCESS, algorithm_from_string
    use simple_cartft_corrcalc,     only: cartft_corrcalc
    use simple_eval_cartftcc,       only: eval_cartftcc
    use simple_cmdline,             only: cmdline
    use simple_builder,             only: builder
    use simple_parameters,          only: parameters
    use simple_cftcc_shsrch_grad,   only: cftcc_shsrch_grad
    use simple_projector_hlev
    use simple_timer
    use simple_oris
    use simple_image
    implicit none
    integer, parameter :: wp = kind(0.0d0)
    type :: constraint_data
        real(wp) :: d(2)
    end type
    type(nlopt_opt) :: opt
    real(wp)        :: lb(2), ub(2), x(2), minf
    integer         :: stat
    type(constraint_data), target    :: d1, d2
    real(wp),              parameter :: xtol = 1.0e-8_wp
    type(parameters)         :: p
    type(cartft_corrcalc)    :: cftcc
    type(cmdline)            :: cline
    type(builder)            :: b
    integer,     parameter   :: NSPACE=1    ! set to 1 for fast test
    real                     :: corrs(NSPACE), corrs2(NSPACE), grad(2, NSPACE), lims(2,2), cxy(3)
    type(image), allocatable :: imgs(:)
    real,        allocatable :: pshifts(:,:)
    integer                  :: iref, iptcl, loc(1), iter
    type(eval_cartftcc)      :: evalcc
    type(cftcc_shsrch_grad)  :: grad_carshsrch_obj
    ! setting up for shift search
    if( command_argument_count() < 3 )then
        write(logfhandle,'(a)') 'simple_test_nlopt lp=xx smpd=yy nthr=zz vol1=vol1.mrc opt=oo'
        stop
    endif
    call cline%parse_oldschool
    call cline%checkvar('lp',   1)
    call cline%checkvar('smpd', 2)
    call cline%checkvar('nthr', 3)
    call cline%checkvar('vol1', 4)
    call cline%checkvar('opt',  5)
    call cline%set('ctf', 'no')
    call cline%set('match_filt', 'no')
    call cline%check
    call p%new(cline)
    p%kfromto(1) = 3
    p%kfromto(2) = calc_fourier_index(p%lp, p%box, p%smpd)
    p%kstop      = p%kfromto(2)
    call b%build_general_tbox(p, cline)
    ! prep projections
    call b%vol%read(p%vols(1))
    call b%eulspace%new(NSPACE, is_ptcl=.false.)
    call b%eulspace%spiral
    p%nptcls = NSPACE
    imgs     = reproject(b%vol, b%eulspace)
    call b%vol%fft
    call b%vol%expand_cmat(KBALPHA) ! necessary for re-projection
    ! prep correlator
    call cftcc%new(p%nptcls, [1, p%nptcls], .false.)
    do iref = 1, p%nptcls
        call imgs(iref)%fft
        call cftcc%set_ref(iref, imgs(iref), .true.)
    end do
    do iptcl = 1,p%nptcls
        call cftcc%set_ptcl(iptcl,imgs(iptcl))
    end do
    ! prep evaluator
    call evalcc%new(b%vol, b%vol, NSPACE)
    allocate(pshifts(p%nptcls, 2), source=0.)
    do iref = 1,p%nptcls
        ! mapping ran3 (0 to 1) to [-5, 5]
        pshifts(iref, 1) = floor(ran3()*10.99) - 5
        pshifts(iref, 2) = floor(ran3()*10.99) - 5
        call evalcc%set_ori(iref, b%eulspace%get_euler(iref), pshifts(iref, :))
    end do
    ! nlopt-f default unit test
    call create(opt, algorithm_from_string(trim(p%opt)), 2)
    call opt%get_lower_bounds(lb)
    lb(2) = 0.0_wp
    call opt%set_lower_bounds(lb)
    d1%d = [+2.0_wp, +0.0_wp]
    d2%d = [-1.0_wp, +1.0_wp]
    associate(f   => nlopt_func(nloptf_myfunc), &
            & fc1 => nlopt_func(nloptf_myconstraint, d1), &
            & fc2 => nlopt_func(nloptf_myconstraint, d2))
        call opt%set_min_objective(f)
        call opt%add_inequality_constraint(fc1, 1.0e-8_wp)
        call opt%add_inequality_constraint(fc2, 1.0e-8_wp)
        call opt%set_xtol_rel(xtol)
        x = [1.234_wp, 5.678_wp]
        call opt%optimize(x, minf, stat)
    end associate
    if (stat < NLOPT_SUCCESS) then
        write(*, '(a)') "NLopt failed!"
        stop 1
    endif
    print *, '---NLOPT-F default test---'
    write(*, '(a, *(1x, g0))') "Found minimum at", x
    write(*, '(a, *(1x, g0))') "Minimum value is", minf
    call destroy(opt)
    ! a unit test with no constraint and with derivative-free optimizer
    call create(opt, algorithm_from_string(trim(p%opt)), 2)
    lb(1) = -100.0_wp
    lb(2) = -100.0_wp
    call opt%set_lower_bounds(lb)
    ub(1) = 100.0_wp
    ub(2) = 100.0_wp
    call opt%set_upper_bounds(ub)
    associate(f => nlopt_func(xy_mixed_func))
        call opt%set_min_objective(f)
        call opt%set_xtol_rel(xtol)
        x = [1.234_wp, 5.678_wp]
        call opt%optimize(x, minf, stat)
    end associate
    if (stat < NLOPT_SUCCESS) then
        write(*, '(a)') "NLopt failed!"
        stop 1
    endif
    print *, '---NLOpt derivative-free unit test---'
    write(*, '(a, *(1x, g0))') "Found minimum at", x
    write(*, '(a, *(1x, g0))') "Minimum value is", minf
    call destroy(opt)
    ! testing with basic NLOpt's gradient optimizer
    iptcl = 1
    iref  = 1
    print *, '---Shift search test using NLOpt optimizer---'
    print *, 'initial shift = ', pshifts(iref, :)
    call create(opt, algorithm_from_string(trim(p%opt)), 2)
    lb(1) = -6.0_wp
    lb(2) = -6.0_wp
    call opt%set_lower_bounds(lb)
    ub(1) = 6.0_wp
    ub(2) = 6.0_wp
    call opt%set_upper_bounds(ub)
    associate(f => nlopt_func(shift_grad_func))
        call opt%set_min_objective(f)
        call opt%set_xtol_rel(xtol)
        x = pshifts(iref, :)
        call opt%optimize(x, minf, stat)
    end associate
    if (stat < NLOPT_SUCCESS) then
        write(*, '(a)') "NLopt failed!"
        stop 1
    endif
    write(*, '(a, *(1x, g0))') "Found minimum at", x
    write(*, '(a, *(1x, g0))') "Minimum value is", minf
    call destroy(opt)
    ! testing with lbfgsb
    print *, '---Shift search test using SIMPLE lbfgsb optimizer---'
    lims(:,1) = -5.
    lims(:,2) =  5.
    call grad_carshsrch_obj%new(lims)
    call grad_carshsrch_obj%set_indices(1, 1)
    cxy = grad_carshsrch_obj%minimize()
    write(*, '(a, *(1x, g0))') "Found minimum at", cxy(2:)
    write(*, '(a, *(1x, g0))') "Minimum value is", cxy(1)

contains
    function nloptf_myfunc(x, gradient, func_data) result(f)
        real(wp), intent(in)              :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in),    optional :: func_data
        real(wp) :: f
        if (present(gradient)) then
            gradient(1) = 0.0_wp
            gradient(2) = 0.5_wp / sqrt(x(2))
        endif
        f = sqrt(x(2))
    end function nloptf_myfunc
    
    function nloptf_myconstraint(x, gradient, func_data) result(f)
        real(wp), intent(in)              :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in),    optional :: func_data
        real(wp) :: f
        select type(func_data)
            type is(constraint_data)
            associate(a => func_data%d(1), b => func_data%d(2))
                if (present(gradient)) then
                    gradient(1) = 3.0_wp * a * (a*x(1) + b)**2
                    gradient(2) = -1.0_wp
                endif
                f = (a*x(1) + b)**3 - x(2)
            end associate
        end select
    end function nloptf_myconstraint

    function xy_mixed_func(x, gradient, func_data) result(f)
        real(wp), intent(in)              :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in),    optional :: func_data
        real(wp) :: f
        if (present(gradient)) then
            gradient(1) = 2*(x(1) - 2.)
            gradient(2) = 3*(x(2) - 4.)**3
        endif
        f = (x(1)-2.)**2 + (x(2)-4.)**4
    end function xy_mixed_func

    function shift_grad_func(x, gradient, func_data) result(f)
        real(wp), intent(in)              :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in),    optional :: func_data
        real(wp) :: f
        call evalcc%set_ori(iref, b%eulspace%get_euler(iref), real(x))
        call evalcc%project_and_correlate(iptcl, corrs, grad)
        if (present(gradient)) then
            gradient(1) = - grad(1,1)
            gradient(2) = - grad(2,1)
        endif
        f = - corrs(1)
    end function shift_grad_func
end program simple_test_nlopt