program simple_test_gsl
    use, intrinsic :: iso_c_binding
    use fgsl
implicit none
    integer(fgsl_int) :: iter, status
    real(fgsl_double), target :: par(5)
    real(fgsl_double), target :: x(2)
    type(fgsl_vector) :: xvec, fvec
    type(fgsl_multimin_fdfminimizer) :: mmin_fdfmin
    type(fgsl_multimin_function_fdf) :: mmin_fdf   
    character(kind=fgsl_char,len=fgsl_strmax) :: name
    type(c_ptr) :: ptr
    real(fgsl_double), parameter :: eps5 = 1.0d-5
    integer(fgsl_int), parameter :: itmax_root = 1000
    real(fgsl_double), pointer :: xptr(:)
    write (*,*) "minimizing the function f(x,y) = 10 * (x-1)^2 + 20 * (y-2)^2 + 30"
    write (*,*) "using Fletcher-Reeves conjugate gradient algorithm"
    x = [5.0_fgsl_double, 7.0_fgsl_double]
    xvec = fgsl_vector_init(1.0_fgsl_double)
    status = fgsl_vector_align(x, 2_fgsl_size_t, xvec, 2_fgsl_size_t, 0_fgsl_size_t, 1_fgsl_size_t)
    mmin_fdfmin = fgsl_multimin_fdfminimizer_alloc(fgsl_multimin_fdfminimizer_conjugate_fr, 2_fgsl_size_t)   
    iter = 0
    par = [1.0_fgsl_double, 2.0_fgsl_double, 10.0_fgsl_double, 20.0_fgsl_double, 30.0_fgsl_double]
    ptr = c_loc(par)
    mmin_fdf = fgsl_multimin_function_fdf_init(my_f,my_df,my_fdf,2_fgsl_size_t,ptr)
    status = fgsl_multimin_fdfminimizer_set(mmin_fdfmin, mmin_fdf, xvec, 0.01_fgsl_double, eps5)    
    call fgsl_vector_free(xvec)
    xvec = fgsl_multimin_fdfminimizer_x(mmin_fdfmin)     
    status = fgsl_vector_align(xptr, xvec)               
    iter = 0
    write (*, *) "       step                    x                         y                      fval                     error"
    do
        iter = iter + 1
        status = fgsl_multimin_fdfminimizer_iterate(mmin_fdfmin)
        if (status /= fgsl_success .or. iter > itmax_root) then
            exit
        end if
        write (*, *) iter, xptr(1), xptr(2), fgsl_multimin_fdfminimizer_minimum(mmin_fdfmin), &
            & abs(fgsl_multimin_fdfminimizer_minimum(mmin_fdfmin) - 30.0)
        fvec = fgsl_multimin_fdfminimizer_gradient(mmin_fdfmin)
        status = fgsl_multimin_test_gradient(fvec, eps5)
        if (status == fgsl_success) then 
            exit
        end if
    end do    
contains
    function my_f(v, params) bind(c)
        type(c_ptr), value :: v
        type(c_ptr), value :: params
        real(fgsl_double) :: my_f       
        real(fgsl_double), pointer :: p(:), pv(:)
        type(fgsl_vector) :: vec
        integer(fgsl_int) :: status
        real(fgsl_double) :: x, y
        call fgsl_obj_c_ptr(vec, v)
        status = fgsl_vector_align(pv, vec)
        call c_f_pointer(params, p, [ 2 ])        
        x = pv(1)
        y = pv(2)
        my_f = p(3) * (x - p(1)) * (x - p(1)) + p(4) * (y - p(2)) * (y - p(2)) + p(5)
    end function my_f
    subroutine my_df(v, params, df)
        type(c_ptr), value :: v, params, df
        type(fgsl_vector) :: vec, df_vec
        real(fgsl_double), pointer :: p(:), pv(:), pdf(:)
        integer(fgsl_int) :: status
        real(fgsl_double) :: x, y
        call fgsl_obj_c_ptr(vec, v)
        status = fgsl_vector_align(pv, vec)
        call fgsl_obj_c_ptr(df_vec, df)
        status = fgsl_vector_align(pdf, df_vec)
        call c_f_pointer(params, p, [2])
        x = pv(1)
        y = pv(2)
        pdf(1) = 2.0 * p(3) * ( x - p(1) )
        pdf(2) = 2.0 * p(4) * ( y - p(2) )       
    end subroutine my_df
    subroutine my_fdf (v, params, f, df)
        type(c_ptr), value :: v, params, df
        real(fgsl_double) :: f       
        f = my_f(v, params)
        call my_df(v, params, df)
    end subroutine my_fdf   
end program simple_test_gsl
