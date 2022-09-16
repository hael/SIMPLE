program simple_test_nlopt
    use nlopt_wrap, only : nlopt_opt, nlopt_func, create, destroy
    use nlopt_enum, only : NLOPT_SUCCESS, algorithm_from_string
    implicit none
    integer, parameter :: wp = kind(0.0d0)
    type :: constraint_data
        real(wp) :: d(2)
    end type
    
    type(nlopt_opt) :: opt
    real(wp) :: lb(2), x(2), minf
    integer :: stat
    type(constraint_data), target :: d1, d2
    real(wp), parameter :: xtol = 1.0e-4_wp
  
    call create(opt, algorithm_from_string("LD_MMA"), 2)
  
    call opt%get_lower_bounds(lb)
  
    lb(2) = 0.0_wp
    call opt%set_lower_bounds(lb)
  
    d1%d = [+2.0_wp, +0.0_wp]
    d2%d = [-1.0_wp, +1.0_wp]
    associate(&
            & f => nlopt_func(myfunc), &
            & fc1 => nlopt_func(myconstraint, d1), &
            & fc2 => nlopt_func(myconstraint, d2))
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
  
    write(*, '(a, *(1x, g0))') "Found minimum at", x
    write(*, '(a, *(1x, g0))') "Minimum value is", minf
  
    call destroy(opt)

contains
    function myfunc(x, gradient, func_data) result(f)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in), optional :: func_data
        real(wp) :: f
    
        if (present(gradient)) then
        gradient(1) = 0.0_wp
        gradient(2) = 0.5_wp / sqrt(x(2))
        endif
        f = sqrt(x(2))
    end function myfunc
    
    function myconstraint(x, gradient, func_data) result(f)
        real(wp), intent(in) :: x(:)
        real(wp), intent(inout), optional :: gradient(:)
        class(*), intent(in), optional :: func_data
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
    end function myconstraint

end program simple_test_nlopt