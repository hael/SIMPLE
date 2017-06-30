!==Class simple_centre_clust
!
!> simple_centre_clust is the SIMPLE class for centre-based clustering of feature vectors.
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2012-02-02.
!
!==Changes are documented below
!
module simple_centre_clust
use simple_oris,     only: oris
use simple_jiffys,   only: alloc_err
use simple_ran_tabu, only: ran_tabu
implicit none

public :: centre_clust
private

type centre_clust
    private
    type(ran_tabu)       :: rt                             !< object for random number generation
    integer, allocatable :: pops(:)                        !< class populations
    real, allocatable    :: sums(:,:), dists(:), avgs(:,:) !< stuff for calculation
    integer, allocatable :: srch_order(:)                  !< random search order 
    real, pointer        :: vecs(:,:)=>null()              !< data vectors
    integer              :: N                              !< nr of ptcls
    integer              :: D                              !< dimension of data vec
    integer              :: ncls                           !< nr of classes
    class(oris), pointer :: o_ptr                          !< ponter to orientation data struct
    logical              :: existence=.false.              !< indicates existence 
  contains
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: get_avgs
    procedure :: set_avgs
    ! PUBLIC
    procedure :: srch_greedy
    procedure :: srch_shc
    ! PRIVATE
    procedure, private :: calc_avgs
    procedure, private :: subtr_ptcl
    procedure, private :: add_ptcl
    procedure, private :: calc_dists
    procedure, private :: sq_dist
    procedure, private :: cost
    procedure :: kill
end type centre_clust

interface centre_clust
    module procedure constructor
end interface centre_clust

contains

    !>  \brief  is a constructor
    function constructor( vecs, o, N, D, ncls ) result( self )
        use simple_oris,  only: oris
        class(oris), intent(in), target :: o          !< for storing cluster info
        integer,     intent(in)         :: N, D, ncls !< params
        real,        intent(in), target :: vecs(N,D)  !< data vectors
        type(centre_clust)              :: self
        call self%new( vecs, o, N, D, ncls ) 
    end function constructor

    !>  \brief  is a constructor
    subroutine new( self, vecs, o, N, D, ncls )
        use simple_oris,  only: oris
        class(centre_clust), intent(inout)      :: self       !< object
        integer,             intent(in)         :: N, D, ncls !< params
        real,                intent(in), target :: vecs(N,D)  !< data vectors
        class(oris),         intent(in), target :: o          !< for storing cluster info
        integer :: alloc_stat
        real    :: x
        call self%kill
        self%vecs  => vecs
        self%o_ptr => o
        self%N     = N
        self%D     = D
        self%ncls  = ncls
        ! construct composites
        allocate( self%pops(ncls), self%sums(ncls,D), self%srch_order(ncls),&
        self%dists(ncls), self%avgs(ncls,D), stat=alloc_stat )
        call alloc_err( "new; simple_centre_clust", alloc_stat )
        self%pops       = 0 
        self%sums       = 0.
        self%srch_order = 0
        self%dists      = huge(x)
        self%avgs       = 0.
        self%rt = ran_tabu(ncls)
        self%existence = .true.
    end subroutine new
    
    ! GETTERS/SETTERS
    
    !>  \brief  is for gettign the averages
    function get_avgs( self ) result( avgs )
        class(centre_clust), intent(in) :: self
        real, allocatable :: avgs(:,:)
        integer :: alloc_stat
        allocate(avgs(size(self%avgs,1),size(self%avgs,2)), stat=alloc_stat)
        call alloc_err('get_avgs; simple_centre_clust', alloc_stat)
        avgs = self%avgs
    end function get_avgs
    
    !>  \brief  is for setting the averages
    subroutine set_avgs( self, avgs )
        class(centre_clust), intent(inout) :: self
        real,                intent(in)    :: avgs(self%ncls,self%D)
        self%avgs = avgs
    end subroutine set_avgs
    
    ! PUBLIC
    
    !>  \brief  does the k-means clustering
    subroutine srch_greedy( self, maxits )
        class(centre_clust), intent(inout) :: self
        integer,             intent(in)    :: maxits
        integer :: i, it, cls, loc(1)
        real    :: adist, adist_prev, x
        write(*,'(A)') '>>> K-MEANS CLUSTERING'
        it = 1
        adist = huge(x)
        call self%calc_avgs 
        do
            adist_prev = adist
            adist = self%cost()
            ! if( it == 1 .or. mod(it,5) == 0 )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F7.3)") 'Iteration:', it, 'Cost:', adist
            ! endif
            do i=1,self%N
                cls = nint(self%o_ptr%get(i, 'class'))
                call self%subtr_ptcl(i,cls)
                call self%calc_dists(i)
                loc = minloc(self%dists)
                call self%o_ptr%set(i, 'class', real(loc(1)))  
                call self%add_ptcl(i,loc(1))
            end do
            if( abs(adist_prev-adist) < 0.0001 .or. it == maxits ) exit
            it = it+1
        end do
    end subroutine srch_greedy
    
    !>  \brief  does the shc clustering
    subroutine srch_shc( self, maxits )
        use simple_opt_subs, only: shc_selector
        class(centre_clust), intent(inout) :: self
        integer,             intent(in)    :: maxits
        integer                            :: i, it, cls, loc(1)
        real                               :: adist, adist_prev, x
        write(*,'(A)') '>>> SHC CLUSTERING'
        it = 1
        adist = huge(x)
        call self%calc_avgs 
        do
            adist_prev = adist
            adist = self%cost()
            ! if( it == 1 .or. mod(it,5) == 0 )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F7.3)") 'Iteration:', it, 'Cost:', adist
            ! endif
            do i=1,self%N
                cls = nint(self%o_ptr%get(i, 'class'))
                call self%subtr_ptcl(i,cls)
                call self%calc_dists(i)
                loc(1) = shc_selector(self%dists, cls)
                call self%o_ptr%set(i, 'class', real(loc(1)))  
                call self%add_ptcl(i,loc(1))
            end do
            if( abs(adist_prev-adist) < 0.0001 .or. it == maxits ) exit
            it = it+1
        end do
    end subroutine srch_shc

    ! PRIVATE
    
    !>  \brief  is for calculating averages given the clustering solution
    subroutine calc_avgs( self )
        class(centre_clust), intent(inout) :: self
        integer :: k, i, cls
        self%sums = 0.
        self%pops = 0
        do i=1,self%N
            cls = nint(self%o_ptr%get(i, 'class'))
            if( cls /= 0 )then
                self%sums(cls,:) = self%sums(cls,:)+self%vecs(i,:)
                self%pops(cls) = self%pops(cls)+1
            endif
        end do
        self%avgs = 0.
        do k=1,self%ncls
            if( self%pops(k) > 1 )then
                self%avgs(k,:) = self%sums(k,:)/real(self%pops(k))
            else
                self%avgs(k,:) = 0.
            endif
        end do
    end subroutine calc_avgs
  
    !>  \brief  is for subtracting the contribution from ptcl, assumes that ptcl is read
    subroutine subtr_ptcl( self, i, cls )
        class(centre_clust), intent(inout) :: self
        integer,             intent(in)    :: i, cls
        if( cls /= 0 )then
            self%sums(cls,:) = self%sums(cls,:)-self%vecs(i,:)
            self%pops(cls)   = self%pops(cls)-1
            self%avgs(cls,:) = self%sums(cls,:)/real(self%pops(cls))
        endif
    end subroutine subtr_ptcl
    
    !>  \brief  is for adding the contribution from ptcl i, assumes that ptcl is read
    subroutine add_ptcl( self, i, cls )
        class(centre_clust), intent(inout) :: self
        integer,             intent(in)    :: i, cls
        self%sums(cls,:) = self%sums(cls,:)+self%vecs(i,:)
        self%pops(cls)   = self%pops(cls)+1
        self%avgs(cls,:) = self%sums(cls,:)/real(self%pops(cls))
    end subroutine add_ptcl
    
    !>  \brief  is for calculating all distances, assumes that  ptcl is read
    subroutine calc_dists( self, i )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(centre_clust), intent(inout) :: self
        integer,             intent(in)    :: i
        integer :: k
        real :: x
        !$omp parallel do schedule(static) default(shared) private(k,x) proc_bind(close)
        do k=1,self%ncls
            if( self%pops(k) > 0 )then
                self%dists(k) = self%sq_dist(i,k)
            else
                self%dists(k) = huge(x) 
            endif
        end do
        !$omp end parallel do
    end subroutine calc_dists
    
    !>  \brief  is for calculating the square distance between data vec and average
    function sq_dist( self, i, k ) result( dist )
        use simple_math, only: euclid
        class(centre_clust), intent(in) :: self
        integer,             intent(in) :: i
        integer,             intent(in) :: k
        real :: dist
        dist = euclid(self%vecs(i,:),self%avgs(k,:))**2.
    end function sq_dist

    !>  \brief  is the cost function being minimized
    function cost( self ) result( adist )
        class(centre_clust), intent(inout) :: self
        real    :: adist
        integer :: i, cnt, cls
        adist = 0.
        cnt = 0 
        do i=1,self%N
            cls = nint(self%o_ptr%get(i, 'class'))  
            if( cls /= 0  )then
                if( self%pops(cls) > 1 )then 
                    adist = adist+self%sq_dist(i,cls)
                    cnt = cnt+1
                endif
            endif
        end do
        adist = adist/real(cnt)
    end function cost
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(centre_clust), intent(inout) :: self
        if( self%existence )then
            deallocate(self%pops, self%sums, self%srch_order, self%dists, self%avgs)
            call self%rt%kill 
            self%existence = .false.
        endif
    end subroutine kill

end module simple_centre_clust
