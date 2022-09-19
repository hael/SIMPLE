program simple_test_autodiff
    use ADLib_NumParameters_m
    use ADdnSVM_m
    implicit none
    type(dnS_t) :: r, th, x
    real        :: d1(2), d0(2)
    r  = Variable( Val=2.0_Rkind, nVar=2, iVar=1, nderiv=1 )
    th = Variable( Val=1.0_Rkind, nVar=2, iVar=2, nderiv=1 )
    x  = test_ad_func(r, th)
    d1 = get_d1(x)
    d0 = get_d0(x)
    write(out_unitp,*) '[x]:', d0
    write(out_unitp,*) '[dx/dr, dx/dth]:', d1
    call set_dnS(r, 0.5_Rkind)
    x  = test_ad_func(r, th)
    d1 = get_d1(x)
    write(out_unitp,*) '[dx/dr, dx/dth]:', d1
contains
    function test_ad_func(r, th) result(f)
        type(dnS_t), intent(in) :: r, th
        type(dnS_t)             :: f
        integer                 :: i
        f = 0._Rkind
        do i = 1, 2
            if( i == 1 )then
                f = f + r*cos(i*th)
            else
                f = f - r*cos(i*th)
            endif
        enddo
    end function test_ad_func
end program simple_test_autodiff