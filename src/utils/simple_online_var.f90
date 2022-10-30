! online moments estimation
module simple_online_var
use simple_defs
implicit none

public :: online_var, test_online_var
private

type :: online_var
    private
    real(dp) :: sumw = 0.d0 !< sum of weights
    real(dp) :: mean = 0.d0 !< moving average (real part)
    real(dp) :: var  = 0.d0 !< variance (real part)
    real(dp) :: cnt  = 0.d0 !< nr of observations
  contains
    procedure, private :: add_1, add_2
    generic            :: add => add_1, add_2
    procedure          :: reset_mean
    procedure          :: serialize
    procedure          :: unserialize
    procedure          :: get_var
    procedure          :: get_mean
end type online_var

contains

    !>  \brief  updates the mean and variance
    subroutine add_1( self, x )
        class(online_var), intent(inout) :: self
        real(sp),          intent(in)    :: x !< new input
        real(dp) :: temp, delta, R
        self%sumw = self%sumw + 1.d0
        delta     = dble(x) - self%mean
        self%mean = self%mean + delta / self%sumw
        self%var  = self%var + delta * (dble(x) - self%mean)
        self%cnt  = self%cnt + 1.d0
    end subroutine add_1

    !>  \brief  updates the mean and variance (weighted)
    subroutine add_2( self, x, w )
        class(online_var), intent(inout) :: self
        real(sp),          intent(in)    :: x !< new input
        real(sp),          intent(in)    :: w !< new input weight
        real(dp) :: temp, delta, R, ww
        ww        = dble(w)
        temp      = ww + self%sumw
        delta     = dble(x) - self%mean
        R         = delta * ww / temp
        self%mean = self%mean + R
        self%var  = self%var + self%sumw * delta * R
        self%sumw = temp
        self%cnt  = self%cnt + 1.d0
    end subroutine add_2

    !>  \brief  for re-setting the mean (needed for cyclic variables)
    subroutine reset_mean( self, mean )
        class(online_var), intent(inout) :: self
        real(sp), intent(in)             :: mean !< new mean
        self%mean = dble(mean)
    end subroutine reset_mean

    !>  \brief  4 serialization of the object
    function serialize( self ) result( arr )
        class(online_var), intent(in) :: self
        real(dp) :: arr(4)
        arr(1) = self%sumw
        arr(2) = self%mean
        arr(3) = self%var
        arr(4) = self%cnt
    end function serialize

    !>  \brief  4 serialization of the object
    subroutine unserialize( self, arr )
        class(online_var), intent(inout) :: self
        real(dp), intent(in) :: arr(4)
        self%sumw = arr(1)
        self%mean = arr(2)
        self%var  = arr(3)
        self%cnt  = arr(4)
    end subroutine unserialize

    !>  \brief  2 get the real mean
    function get_mean( self ) result( mean )
        class(online_var), intent(in) :: self
        real(sp) :: mean
        mean = real(self%mean, kind=sp)
    end function get_mean

    !>  \brief  2 get the real variance
    function get_var( self ) result( var )
        class(online_var), intent(inout) :: self
        real(dp) :: var_n
        real(sp) :: var
        if( self%cnt > 1.d0 )then
            var_n = self%var / self%sumw
            var   = real(var_n * self%cnt / (self%cnt - 1.d0), kind=sp)
        else
            var = 0.
        endif
    end function get_var

    !>  \brief  is the unit test associated with this class
    subroutine test_online_var
        use simple_stat, only: moment
        use simple_rnd,  only: gasdev
        real    :: samples(10000), ave, sdev, var, mv(2)
        integer :: i
        logical :: err
        type(online_var) :: owv, owv2
        write(logfhandle,'(a)') '**info(simple_online_var_unit_test): testing all functionality'
        do i=1,10000
            samples(i) = gasdev(5., 2.)
        end do
        call moment(samples, ave, sdev, var, err )
        write(logfhandle,*) 'classical estimation, ave/sdev:', ave, sdev
        do i=1,10000
            call owv%add(samples(i))
        end do
        owv2 = owv
        write(logfhandle,*) 'online estimation, ave/sdev:',  owv2%get_mean(), sqrt(owv2%get_var())
        if( abs(owv2%get_mean() - ave) < 0.005 .and. abs(sqrt(owv2%get_var()) - sdev) < 0.005 )then
            write(logfhandle,'(a)') 'SIMPLE_ONLINE_var_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
        endif
    end subroutine test_online_var

end module simple_online_var
