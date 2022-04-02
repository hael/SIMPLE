
module simple_testfuns_butterworth
    use simple_defs
    implicit none
    
    !>  \brief  defines the test function interface
    abstract interface
        function testfun_butterworth( fun_self, vec, D ) result( cost )
            class(*), intent(inout) :: fun_self
            integer,  intent(in)    :: D
            real,     intent(in)    :: vec(D)
            real                    :: cost
        end function
    end interface
    
    contains
        function ButterworthCost( fun_self, x, d ) result( r )
            class(*), intent(inout) :: fun_self        
            integer, intent(in) :: d
            real, intent(in)    :: x(d)
            integer :: i
            real :: r
            r = -2.
            do i=1,d
                r = r + (x(i)-3.)**2
            end do
        end function
end module
    