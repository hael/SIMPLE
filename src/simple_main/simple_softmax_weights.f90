module simple_softmax_weights
implicit none

public :: softmax_weights, test_softmax_weights
private

type softmax_weights
    private
    integer           :: n
    real              :: t=1. ! temperature factor, hardens the distribution close to 0
    real, allocatable :: w(:)
  contains
    procedure, private :: get_1
    procedure, private :: get_2
    generic            :: get => get_1, get_2
    procedure          :: plot_weights
    procedure          :: kill
end type

interface softmax_weights
    module procedure constructor
end interface

contains
    
    function constructor( n, t ) result( self )
        use simple_jiffys, only: alloc_err
        integer,        intent(in) :: n
        real, optional, intent(in) :: t
        integer                    :: alloc_stat, i
        real                       :: wsum, a
        type(softmax_weights)      :: self
        call self%kill
        if( present(t) )then
            if( t > 1. .or. t < 0.01 ) stop 'temperature out of range; constructor; simple_softmax_weights'
            self%t = t
        endif
        self%n = n
        allocate( self%w(n), stat=alloc_stat )
        call alloc_err( 'constructor; softmax_weights', alloc_stat )
        wsum = 0.
        do i=0,n-1
            a           = real(n-i)/real(n)
            self%w(i+1) = exp(a/self%t)
            wsum        = wsum+self%w(i+1) ! sum 4 normalization
        end do
        do i=1,n
            self%w(i) = self%w(i)/wsum   ! normalized exponential weight
        end do
    end function
    
    function get_1( self, i ) result( w )
        class(softmax_weights), intent(in) :: self
        integer, intent(in) :: i
        real :: w
        if( i < 0 .or. i > self%n ) stop 'weight index out of range; get_1; simple_softmax_weights'
        w = self%w(i)
    end function
    
    function get_2( self ) result( w )
        use simple_jiffys, only: alloc_err
        class(softmax_weights), intent(in) :: self
        real, allocatable                 :: w(:)
        integer                           :: alloc_stat
        allocate( w(size(self%w)), source=self%w, stat=alloc_stat )
        call alloc_err( 'get_2; softmax_weights', alloc_stat )
    end function
    
    subroutine plot_weights( self )
        use gnufor2, only: plot
        class(softmax_weights), intent(in) :: self
        real    :: x(self%n)
        integer :: i
        do i=1,self%n
            x(i) = real(i)
        end do
        call plot(x,self%w)
    end subroutine
    
    subroutine kill( self )
        class(softmax_weights), intent(inout) :: self
        if( allocated(self%w) ) deallocate(self%w)
    end subroutine
    
    subroutine test_softmax_weights
        type(softmax_weights) :: lew
        lew = softmax_weights(500,1.)
        call lew%plot_weights
        lew = softmax_weights(500,0.8)
        call lew%plot_weights
        lew = softmax_weights(500,0.6)
        call lew%plot_weights
        lew = softmax_weights(500,0.4)
        call lew%plot_weights
        lew = softmax_weights(500,0.2)
        call lew%plot_weights
        lew = softmax_weights(500,0.1)
        call lew%plot_weights
        lew = softmax_weights(500,0.08)
        call lew%plot_weights
        lew = softmax_weights(500,0.06)
        call lew%plot_weights
        lew = softmax_weights(500,0.04)
        call lew%plot_weights
        lew = softmax_weights(500,0.02)
        call lew%plot_weights
    end subroutine

end module simple_softmax_weights
