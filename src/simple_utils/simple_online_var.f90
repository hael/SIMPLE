! online moments estimation
module simple_online_var
use simple_defs ! singleton
implicit none

public :: online_var, test_online_var
private

type :: online_var
    private
    real(dp) :: sumw=0. !< sum of weights
    real(dp) :: mean=0. !< moving average (real part)
    real(dp) :: var=0.  !< variance (real part)
    real(dp) :: cnt=0.  !< nr of observations
  contains
    procedure :: add
    procedure :: reset_mean
    procedure :: finalize
    procedure :: serialize
    procedure :: unserialize
    procedure :: get_var
    procedure :: get_mean
    procedure :: get_mean_var
end type online_var

contains

    !>  \brief  updates the mean and variance
    subroutine add( self, x, w )
        class(online_var), intent(inout) :: self
        real(sp), intent(in)             :: x !< new input
        real(sp), intent(in), optional   :: w !< new input weight
        real(dp) :: temp, delta, R, ww
        ww = 1.
        if( present(w) )then
            ww        = dble(w)
            temp      = ww+self%sumw
            delta     = dble(x)-self%mean
            R         = delta*ww/temp
            self%mean = self%mean+R
            self%var  = self%var+self%sumw*delta*R
            self%sumw = temp
        else
            self%sumw = self%sumw+1.d0
            delta     = dble(x)-self%mean
            self%mean = self%mean+delta/self%sumw
            self%var  = self%var+delta*(dble(x)-self%mean)
        endif
        self%cnt = self%cnt+1.d0
    end subroutine add
    
    !>  \brief  for re-setting the mean (needed for cyclic variables)
    subroutine reset_mean( self, mean ) 
        class(online_var), intent(inout) :: self
        real(sp), intent(in)             :: mean !< new mean
        self%mean = dble(mean)
    end subroutine reset_mean
    
    !>  \brief  finalizes the variance
    subroutine finalize( self )
        class(online_var), intent(inout) :: self
        real(dp) :: var_n
        if( self%cnt > 1.d0 )then
            var_n = self%var/self%sumw
            self%var = var_n*(self%cnt)/(self%cnt-1.d0)
        else
            self%var = 0.d0
        endif
    end subroutine finalize
    
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
    
    !>  \brief  2 get the real variance
    function get_var( self ) result( var )
        class(online_var), intent(inout) :: self
        real(dp) :: var_n
        real(sp) :: var
        if( self%cnt > 1.d0 )then
            var_n = self%var/self%sumw
            var = real(var_n*(self%cnt)/(self%cnt-1.d0))
        else
            var = 0.
        endif
    end function get_var
    
    !>  \brief  2 get the real mean
    function get_mean( self ) result( mean )
        class(online_var), intent(in) :: self
        real(sp) :: mean
        mean = real(self%mean)
    end function get_mean
    
    !>  \brief  get mean and var
    function get_mean_var( self ) result( mv )
        class(online_var), intent(in) :: self
        real(sp) :: mv(2)
        mv(1) = real(self%mean)
        mv(2) = real(self%var)
    end function get_mean_var
    
    !>  \brief  is the unit test associated with this class
    subroutine test_online_var
        use simple_stat, only: moment
        use simple_rnd,  only: gasdev
        real    :: samples(10000), ave, sdev, var, mv(2)
        integer :: i
        logical :: err
        type(online_var) :: owv, owv2
        write(*,'(a)') '**info(simple_online_var_unit_test): testing all functionality'
        do i=1,10000
            samples(i) = gasdev(5., 2.)
        end do
        call moment(samples, ave, sdev, var, err )
        write(*,*) 'classical estimation, ave/sdev:', ave, sdev
        do i=1,10000
            call owv%add(samples(i))
        end do
        call owv%finalize
        owv2 = owv
        mv = owv2%get_mean_var()
        write(*,*) 'online estimation, ave/sdev:', mv(1), sqrt(mv(2))
        if( abs(mv(1)-ave) < 0.005 .and. abs(sqrt(mv(2))-sdev) < 0.005 )then
            write(*,'(a)') 'SIMPLE_ONLINE_var_UNIT_TEST COMPLETED SUCCESSFULLY ;-)'
        endif
    end subroutine test_online_var

end module simple_online_var
