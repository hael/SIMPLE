!@descr: online moments estimation
! Welford's update: the running mean and the sum of squared deviations from it, in double precision,
! so a large offset in single-precision input costs no accuracy.
module simple_online_var
use simple_defs
implicit none

public :: online_var
private

type :: online_var
    private
    real(dp) :: cnt  = 0.d0 !< nr of observations
    real(dp) :: mean = 0.d0 !< running mean
    real(dp) :: var  = 0.d0 !< sum of squared deviations from the running mean
  contains
    procedure :: add
    procedure :: get_var
    procedure :: get_mean
end type online_var

contains

    !>  \brief  updates the mean and variance
    subroutine add( self, x )
        class(online_var), intent(inout) :: self
        real(sp),          intent(in)    :: x !< new input
        real(dp) :: delta
        self%cnt  = self%cnt + 1.d0
        delta     = dble(x) - self%mean
        self%mean = self%mean + delta / self%cnt
        self%var  = self%var + delta * (dble(x) - self%mean)
    end subroutine add

    !>  \brief  the mean (0 without observations)
    function get_mean( self ) result( mean )
        class(online_var), intent(in) :: self
        real(sp) :: mean
        mean = real(self%mean, kind=sp)
    end function get_mean

    !>  \brief  the sample variance, n-1 in the denominator (0 below two observations)
    function get_var( self ) result( var )
        class(online_var), intent(in) :: self
        real(sp) :: var
        if( self%cnt > 1.d0 )then
            var = real(self%var / (self%cnt - 1.d0), kind=sp)
        else
            var = 0.
        endif
    end function get_var

end module simple_online_var
