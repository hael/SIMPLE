!>  \brief  SIMPLE spatial_median class
module simple_spatial_median 
use simple_bfgs_opt, only: bfgs_opt
use simple_opt_spec, only: opt_spec
implicit none

public :: spatial_median, test_spatial_median
private

type :: spatial_median
    private
    real, pointer  :: points(:,:)=>null() !< point cloud
    type(bfgs_opt) :: opt                 !< optimizer
    type(opt_spec) :: ospec               !< optimizer specification
    integer        :: n=0, ndim=0         !< nr of points, dim of point
    logical        :: existence=.false.   !< to indicate existence
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! CALCULATORS
    procedure :: find_median
    procedure :: closest2median
    ! DESTRUCTOR
    procedure :: kill
end type

interface spatial_median
    module procedure constructor
end interface

real, parameter :: TOL=1e-10
real, pointer   :: ptr_points(:,:) => null()
integer         :: n_glob

contains

    !>  \brief  constructor
    function constructor( points ) result( self )
        type(spatial_median) :: self
        real, intent(in), target :: points(:,:)
        call self%new(points)
    end function constructor
    
    !>  \brief  constructor
    subroutine new( self, points )
        class(spatial_median), intent(inout) :: self
        real, intent(in), target :: points(:,:)
        self%n    = size(points,1) ! nr of points
        self%ndim = size(points,2) ! point dimension
        self%points => points      ! set pointer to point cloud
        call self%ospec%specify('bfgs', self%ndim, ftol=TOL, gtol=TOL)
        call self%opt%new(self%ospec)
    end subroutine new
    
    !>  \brief  does the hard work
    function find_median( self ) result( median )
        class(spatial_median), intent(inout) :: self
        real, allocatable :: median(:)
        real :: dist
        ! set globals
        ptr_points => self%points
        n_glob     = self%n
        ! specify 
        call self%ospec%set_costfun(dist2points)
        call self%ospec%set_gcostfun(gdist2points)
        ! initialize with point closest to median
        self%ospec%x = self%points(self%closest2median(),:)
        call self%opt%minimize(self%ospec, dist)
        allocate(median(self%ndim), source=self%ospec%x)
    end function find_median

    !>  \brief  returns the point closest to the spatial median
    !!          can be used as a starting point 4 optimization
    function closest2median( self ) result( this )
        use simple_math, only: euclid
        class(spatial_median), intent(in) :: self
        integer :: this, i, j, loc(1)
        real    :: dists(self%n)
        do i=1,self%n
            dists(i) = 0.
            do j=1,self%n
                if( i /= j )then
                    dists(i) = dists(i)+euclid(self%points(i,:),self%points(j,:))
                endif
            end do
        end do
        loc = minloc(dists)
        this = loc(1)
    end function closest2median
  
    !>  \brief  destructor
    subroutine kill( self )
        class(spatial_median), intent(inout) :: self
        if( self%existence )then
            self%points => null()
            call self%opt%kill
            call self%ospec%kill
            self%n = 0
            self%ndim = 0
            self%existence = .false.
        endif
    end subroutine kill

    ! PRIVATES

    !>  \brief  distance between vec and the point cloud
    function dist2points( vec, D ) result( dist ) 
        use simple_math, only: euclid
        integer, intent(in) :: D
        real,    intent(in) :: vec(D)
        real    :: dist
        integer :: i
        dist = 0.
        do i=1,n_glob
            dist = dist+euclid(vec,ptr_points(i,:))
        end do
    end function dist2points
    
    !>  \brief  gradient of the above function 
    function gdist2points( vec, D ) result( grad )
        use simple_math, only: euclid
        integer, intent(in)    :: D
        real,    intent(inout) :: vec(D)
        real    :: grad(D), dist
        integer :: i
        grad = 0.
        do i=1,n_glob
            dist = euclid(vec,ptr_points(i,:))
            if( dist > TOL ) grad(:) = grad(:)+(vec(:)-ptr_points(i,:))/dist
        end do
    end function gdist2points
    
    ! UNIT TEST
    
    !>  \brief  unit test
    subroutine test_spatial_median
        use simple_math, only: euclid
        use simple_rnd,  only: mnorm_smp
        real    :: Imat(5,5), mean(5), dist
        real    :: points(100,5), centroid(5), median(5)
        integer :: i,j
        type(spatial_median) :: sm
        do j=1,100
            ! generate 100 points with mean vec 1,2,3,4,5 and unit variance
            Imat=0.; do i=1,5 ; Imat(i,i)=1. ; mean(i)=real(i) ; end do
            do i=1,100
                points(i,:) = mnorm_smp(Imat, 5, mean)
            end do
            ! calculate the centroid of the point cloud
            do i=1,5
                centroid(i) = sum(points(:,i))/100.
            end do
            ! calculate the spatial median of the point cloud
            sm = spatial_median(points)
            median = sm%find_median()
            dist = euclid(median, [1.,2.,3.,4.,5.]) 
            if( dist < 0.5 )then
                ! alles ok
            else
                write(*,*) 'dist:', dist
                stop 'test_spatial_median failed!'
            endif
        end do
        write(*,'(a)') 'SIMPLE_SPATIAL_MEDIAN UNIT TEST COMPLETED ;-)'
    end subroutine test_spatial_median
  
end module simple_spatial_median
