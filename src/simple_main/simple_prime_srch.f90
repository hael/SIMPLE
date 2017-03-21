module simple_prime_srch
use simple_defs
use simple_params, only: params
use simple_jiffys, only: alloc_err
implicit none

type prime_srch
    private
    class(params), pointer :: pp => null()     !< pointer to parameters
    integer                :: nrots = 0        !< number of in-plane rotations in polar representation    
    real, allocatable      :: angtab(:)        !< list of in-plane angles
    logical                :: exists = .false. !< 2 indicate existence
  contains
    ! CONSTRUCTOR 
    procedure :: new
    ! GETTERS
    procedure :: roind
    procedure :: rot
    ! DESTRUCTOR
    procedure :: kill
end type prime_srch

interface prime_srch
    module procedure constructor
end interface prime_srch

contains
    
    !>  \brief  is a constructor
    function constructor( p ) result( self )
        class(params), target, intent(in) :: p !< parameters
        type(prime_srch) :: self
        call self%new(p)
    end function constructor
    
    !>  \brief  is a constructor
    subroutine new( self, p )
        use simple_math, only: round2even, rad2deg
        class(prime_srch),     intent(inout) :: self  !< instance
        class(params), target, intent(in)    :: p     !< parameters
        integer :: i, alloc_stat
        real    :: dang
        call self%kill
        self%pp    => p
        self%nrots = round2even(twopi*real(p%ring2))
        ! generate the array of in-plane angles
        allocate(self%angtab(self%nrots), stat=alloc_stat)
        call alloc_err('In: new; simple_prime_srch', alloc_stat)
        dang = twopi/real(self%nrots)
        forall( i=1:self%nrots ) self%angtab(i) = rad2deg(real(i-1)*dang)
        self%exists = .true.
    end subroutine new
    
    !>  \brief calculates the in-plane rotational index for the rot in-plane angle
    function roind( self, rot ) result( ind )
        class(prime_srch), intent(in) :: self
        real,              intent(in) :: rot
        real    :: dists(self%nrots)
        integer :: ind, loc(1)
        dists = abs(self%angtab-rot)
        where(dists>180.)dists = 360.-dists
        loc = minloc(dists)
        ind = loc(1)
    end function roind
    
    !>  \brief calculates the in-plane angle of in-plane rotational index rot
    function rot( self, ind ) result( r )
        class(prime_srch), intent(in) :: self
        integer,           intent(in) :: ind
        real :: r
        r = self%angtab(ind)
        if( r == 360. ) r = 0.
    end function rot
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(prime_srch), intent(inout) :: self !< instance
        if( self%exists )then
            deallocate(self%angtab)
            self%pp     => null()
            self%nrots  =  0
            self%exists = .false.
        endif
    end subroutine kill

end module simple_prime_srch
