!> simple_kmeans is the SIMPLE class for k-means refinement of the clustering solution in o. 
! The code is distributed with the hope that it will be useful, but _WITHOUT_ _ANY_ _WARRANTY_. 
! Redistribution or modification is regulated by the GNU General Public License. 
! *Author:* Hans Elmlund, 2012-02-02.
!
!==Changes are documented below
!
module simple_kmeans
use simple_oris,   only: oris
use simple_jiffys, only: alloc_err
implicit none

public :: kmeans, test_kmeans
private

type kmeans
    private
    integer, allocatable  :: pops(:)                        !< class populations
    real, allocatable     :: sums(:,:), dists(:), avgs(:,:) !< stuff for calculation
    real, pointer         :: vecs(:,:)=>null()              !< data vectors
    integer               :: N                              !< nr of ptcls
    integer               :: D                              !< dimension of data vec
    integer               :: ncls                           !< nr of classes
    class(oris), pointer  :: o_ptr                          !< ponter to orientation data struct
    logical               :: existence=.false.              !< indicates existence 
  contains
    procedure :: new
    ! GETTERS/SETTERS
    procedure :: get_avgs
    procedure :: set_avgs
    ! PUBLIC
    procedure :: refine
    ! PRIVATE
    procedure, private :: calc_avgs
    procedure, private :: assign_class
    procedure, private :: subtr_ptcl
    procedure, private :: add_ptcl
    procedure, private :: calc_dists
    procedure, private :: sq_dist
    procedure, private :: cost
    procedure :: kill
end type

interface kmeans
    module procedure constructor
end interface

contains

    !>  \brief  is a constructor
    function constructor( vecs, o, N, D, ncls ) result( self )
        use simple_oris,  only: oris
        class(oris), intent(in), target :: o          !< for storing cluster info
        integer, intent(in)             :: N, D, ncls !< params
        real, intent(in), target        :: vecs(N,D)  !< data vectors
        type(kmeans)                    :: self
        call self%new( vecs, o, N, D, ncls ) 
    end function

    !>  \brief  is a constructor
    subroutine new( self, vecs, o, N, D, ncls )
        use simple_oris,  only: oris
        class(kmeans), intent(inout)    :: self       !< object
        integer, intent(in)             :: N, D, ncls !< params
        real, intent(in), target        :: vecs(N,D)  !< data vectors
        class(oris), intent(in), target :: o          !< for storing cluster info
        integer                         :: alloc_stat
        real                            :: x
        call self%kill
        self%vecs  => vecs
        self%o_ptr => o
        self%N     = N
        self%D     = D
        self%ncls  = ncls
        allocate( self%pops(ncls), self%sums(ncls,D),&
        self%dists(ncls), self%avgs(ncls,D), stat=alloc_stat )
        call alloc_err( "new; simple_kmeans", alloc_stat )
        self%pops  = 0
        self%sums  = 0.
        self%dists = huge(x)
        self%avgs  = 0.
        self%existence = .true.
    end subroutine
    
    ! GETTERS/SETTERS
    
    !>  \brief  is for gettign the averages
    function get_avgs( self ) result( avgs )
        class(kmeans), intent(in) :: self
        real, allocatable :: avgs(:,:)
        integer :: alloc_stat
        allocate(avgs(size(self%avgs,1),size(self%avgs,2)), stat=alloc_stat)
        call alloc_err('get_avgs; simple_kmeans', alloc_stat)
        avgs = self%avgs
    end function
    
    !>  \brief  is for setting the averages
    subroutine set_avgs( self, avgs )
        class(kmeans), intent(inout) :: self
        real, intent(in)             :: avgs(self%ncls,self%D)
        self%avgs = avgs
    end subroutine
    
    ! PUBLIC
    
    !>  \brief  this one does it all
    subroutine refine( self, maxits )
        class(kmeans), intent(inout) :: self
        integer, intent(in)          :: maxits
        integer                      :: i, it
        real                         :: adist, adist_prev, x
        write(*,'(A)') '>>> K-MEANS REFINEMENT'
        it = 1
        adist = huge(x)
        call self%calc_avgs 
        do
            adist_prev = adist
            adist = self%cost()
            if( it == 1 .or. mod(it,5) == 0 )then
                write(*,"(1X,A,1X,I3,1X,A,1X,F7.3)") 'Iteration:', it, 'Cost:', adist
            endif
            do i=1,self%N
                call self%assign_class(i)
            end do
            if( abs(adist_prev-adist) < 0.0001 .or. it == maxits ) exit
            it = it+1
        end do
    end subroutine

    ! PRIVATE
    
    !>  \brief  is for calculating averages given the clustering solution, assumes that data stack is open
    subroutine calc_avgs( self )
        class(kmeans), intent(inout) :: self
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
    end subroutine
    
    !>  \brief  is for assigning class to ptcl, assumes that data stack is open
    subroutine assign_class( self, i )
        class(kmeans), intent(inout) :: self
        integer, intent(in) :: i
        integer :: loc(1), cls
        cls = nint(self%o_ptr%get(i, 'class'))
        call self%subtr_ptcl(i,cls)
        call self%calc_dists(i)
        loc = minloc(self%dists)
        call self%o_ptr%set(i, 'class', real(loc(1)))  
        call self%add_ptcl(i,loc(1))
    end subroutine
  
    !>  \brief  is for subtracting the contribution from ptcl, assumes that ptcl is read
    subroutine subtr_ptcl( self, i, cls )
        class(kmeans), intent(inout) :: self
        integer, intent(in)          :: i
        integer, intent(in)          :: cls
        if( cls /= 0 )then
            self%sums(cls,:) = self%sums(cls,:)-self%vecs(i,:)
            self%pops(cls)   = self%pops(cls)-1
            self%avgs(cls,:) = self%sums(cls,:)/real(self%pops(cls))
        endif
    end subroutine
    
    !>  \brief  is for adding the contribution from ptcl i, assumes that ptcl is read
    subroutine add_ptcl( self, i, cls )
        class(kmeans), intent(inout) :: self
        integer, intent(in)          :: i
        integer, intent(in)          :: cls
        self%sums(cls,:) = self%sums(cls,:)+self%vecs(i,:)
        self%pops(cls)   = self%pops(cls)+1
        self%avgs(cls,:) = self%sums(cls,:)/real(self%pops(cls))
    end subroutine
    
    !>  \brief  is for claculating all distances, assumes that  ptcl is read
    subroutine calc_dists( self, i )
        class(kmeans), intent(inout) :: self
        integer, intent(in) :: i
        integer :: k
        real :: x
        do k=1,self%ncls
            if( self%pops(k) > 0 )then
                self%dists(k) = self%sq_dist(i,k)
            else
                self%dists(k) = huge(x) 
            endif
        end do
    end subroutine
    
    !>  \brief  is for calculating the square distance between data vec and average, assumes that ptcl is read
    function sq_dist( self, i, k ) result( dist )
        use simple_math, only: euclid
        class(kmeans), intent(in) :: self
        integer, intent(in)       :: i
        integer, intent(in)       :: k
        real :: dist
        dist = euclid(self%vecs(i,:),self%avgs(k,:))**2.
    end function

    !>  \brief  is the cost function being minimized
    function cost( self ) result( adist )
        class(kmeans), intent(inout) :: self
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
    end function
    
    ! TEST KMEANS
    
    !>  \brief  is the kmeans unit test
    subroutine test_kmeans
        !$ use omp_lib
        !$ use omp_lib_kinds
        use simple_pair_dtab, only: pair_dtab
        use simple_math
        use simple_hac
        real :: datavecs(900,5)
        type(pair_dtab) :: pd
        type(oris)      :: o
        type(hac)       :: hacls
        type(kmeans)    :: kmcls
        integer :: i, ia, ib
        write(*,'(a)') '**info(simple_kmeans_unit_test): testing all functionality'
        call o%new(900)
        pd = pair_dtab(100)
        call hacls%new(900, o, 100)
        ! make data
        do i=1,300
            datavecs(i,:) = 1.
        end do
        do i=301,600
            datavecs(i,:) = 3.
        end do
        do i=601,900
            datavecs(i,:) = 5.
        end do
        !$omp parallel default(shared) private(ib) proc_bind(close)
        do ia=1,100-1
            !$omp do schedule(static) 
            do ib=ia+1,100
                call pd%set_pair_d(ia, ib, euclid(datavecs(hacls%get_node(ia),:),datavecs(hacls%get_node(ib),:)))
            end do
            !$omp end do nowait
        end do
        !$omp end parallel
        call hacls%cluster(pd, 'pdfile.bin', 3, 1 )
        call o%write('test_kmeans_hcl.txt')
        call kmcls%new(datavecs, o, 900, 5, 3)
        call kmcls%refine(100)
        call o%write( 'test_kmeans_kmeans.txt' )
        call kmcls%kill
        call hacls%kill
        call o%kill
        call pd%kill
        write(*,'(a)') 'SIMPLE_KMEANS_UNIT_TEST COMPLETED WITHOUT TERMINAL BUGS ;-)'
        write(*,'(a)') 'PLEASE, INSPECT THE RESULTS'
    end subroutine
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(kmeans), intent(inout) :: self
        if( self%existence )then
            deallocate(self%pops, self%sums, self%dists, self%avgs )
            self%existence = .false.
        endif
    end subroutine

end module simple_kmeans
