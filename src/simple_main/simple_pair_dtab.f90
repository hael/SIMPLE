!==Class simple_pair_dtab

module simple_pair_dtab
use simple_jiffys, only: alloc_err
implicit none

public :: pair_dtab
private

type :: pair_d
    real, pointer :: d=>null()
end type

type pair_dtab
    private
    integer                   :: N=0, minloc(2)=0
    type(pair_d), allocatable :: pd(:,:)
    real, allocatable         :: minvals(:), row(:)
    integer, allocatable      :: minpos(:)
    real                      :: minval
    logical                   :: existence=.false.
  contains
    ! CONSTRUCTOR
    procedure :: new
    ! I/O
    procedure :: write
    procedure :: read
    ! GETTERS/SETTERS
    procedure :: get_nobjs
    procedure :: set_pair_d
    procedure :: get_dmin_pair
    procedure :: get_dmin
    procedure :: pair_d_exist
    procedure :: get_pair_d
    ! MODIFIERS
    procedure :: build
    procedure :: clear
    procedure :: alloc_pair_d
    procedure :: merge_pair_d
    procedure :: remove_pair_d
    procedure :: linear_ordering
    ! PRIVATE STUFF
    procedure :: update_nn
    procedure :: update_best
    ! DESTRUCTOR
    procedure :: kill
end type

interface pair_dtab
    module procedure constructor
end interface 

contains

    !>  \brief  constructs a new pair weight table
    function constructor( Nobjs ) result( self )
        integer, intent(in) :: Nobjs
        type(pair_dtab)     :: self
        call self%new( Nobjs )
    end function

    !>  \brief  constructs a new pair weight table
    subroutine new( self, Nobjs )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: Nobjs
        integer                         :: alloc_stat
        real                            :: x
        call self%kill
        allocate(self%pd(Nobjs,Nobjs), self%minvals(Nobjs),&
        self%minpos(Nobjs), self%row(Nobjs), stat=alloc_stat)
        call alloc_err('In: make_pair_dtab, module: simple_pair_dtab.f90', alloc_stat)
        self%minvals   = huge(x)
        self%minval    = huge(x)
        self%N         = Nobjs
        self%minpos    = 0
        self%row       = 0.
        self%existence = .true.
    end subroutine
    
    !>  \brief  is for flushing a table to disk
    subroutine write( self, fname )
        class(pair_dtab) :: self
        character(len=*) :: fname
        integer          :: wtl, file_stat, i, j
        real             :: x
        ! Check rec-length and open file
        inquire(iolength=wtl) self%row
        open(unit=12, file=fname, status='replace', iostat=file_stat,&
        access='direct', action='write', form='unformatted', recl=wtl )
        if( file_stat /= 0 )then ! cannot open file
            write(*,*) 'Cannot open file: ', fname
            write(*,*) 'In: flush_pair_dtab, module: simple_pair_dtab.f90'
            stop
        endif
        do i=1,self%N
            do j=1,self%N
                if( i == j )then
                    self%row(j) = huge(x)
                else
                    self%row(j) = self%get_pair_d(i, j)
                endif
            end do
            write(12,rec=i) self%row
        end do
        close(12)
    end subroutine
    
    !>  \brief  is for recovering a flushed table 
    subroutine read( self, fname )
        class(pair_dtab) :: self
        character(len=*) :: fname
        integer          :: wtl, file_stat, i, j
        ! Check rec-length and open file
        inquire(iolength=wtl) self%row
        open(unit=12, file=fname, status='old', iostat=file_stat,&
        access='direct', action='read', form='unformatted', recl=wtl )
        if( file_stat /= 0 )then ! cannot open file
            write(*,*) 'Cannot open file: ', fname
            write(*,*) 'In: recover_pair_dtab, module: simple_pair_dtab.f90'
            stop
        endif
        do i=1,self%N
            read(12,rec=i) self%row
            do j=1,self%N
                if( i /= j )then
                    call set_pair_d(self,i,j,self%row(j))
                endif
            end do
        end do
        close(12)
    end subroutine
    
    ! GETTERS/SETTERS
    
    !>  \brief  is for getting the number of objects
    function get_nobjs( self ) result( n )
        class(pair_dtab), intent(in) :: self
        integer                      :: n
        n = self%n
    end function
    
    !>  \brief  is for setting a pair distance
    subroutine set_pair_d( self, i, j, d )
        class(pair_dtab)    :: self
        integer, intent(in) :: i, j
        real, intent(in)    :: d
        integer             :: alloc_stat
        if( j == i )then
            return
        else
            if( .not. associated(self%pd(i,j)%d) )then
                allocate( self%pd(i,j)%d, stat=alloc_stat ) ! make space for the new distance
                call alloc_err('In: set_pair_d, module: simple_pair_dtab.f90', alloc_stat)
            endif
            self%pd(i,j)%d = d
            self%pd(j,i)%d => self%pd(i,j)%d ! Introduction of mirror symmetry in the lookup
            ! update nearest neighs
            if( d < self%minvals(i) )then
                self%minvals(i) = d
                self%minpos(i)  = j
            endif    
            if( d < self%minvals(j) )then
                self%minvals(j) = d
                self%minpos(j)  = i
            endif
            ! update best  
            if( d < self%minval )then
                self%minval    = d
                self%minloc(1) = i
                self%minloc(2) = j
            endif
        endif     
    end subroutine
    
    !>  \brief  is for getting the minimum dist
    subroutine get_dmin_pair( self, k, l, d )
        class(pair_dtab), intent(in) :: self
        integer, intent(out)         :: k, l
        real, intent(out)            :: d
        k = self%minloc(1)
        l = self%minloc(2)
        d = self%minval
    end subroutine
    
    !>  \brief  is for getting the minimum dist
    subroutine get_dmin( self, k, l, d )
        class(pair_dtab), intent(in) :: self
        integer, intent(in)          :: k
        integer, intent(out)         :: l
        real, intent(out)            :: d
        integer                      :: i
        real                         :: x
        l = 0
        d = huge(x)
        do i=1,self%N
            if( associated(self%pd(k,i)%d) )then
                 if( self%pd(k,i)%d < d )then
                    d = self%pd(k,i)%d
                    l = i
                 endif
            endif
        end do
    end subroutine
    
    !>  \brief  checks existence of a pair weight
    function pair_d_exist( self, i, j ) result(bool)
        class(pair_dtab) , intent(in) :: self
        integer, intent(in)           :: i, j
        logical                       :: bool
        if( j == i )then
            bool = .false.
        else
            bool = associated(self%pd(i,j)%d)
        endif
    end function
    
    !>  \brief  is for getting a pair weight
    function get_pair_d( self, i, j ) result( d )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: i, j
        real                            :: d, x
        d = huge(x)
        if( pair_d_exist(self,i,j) )then
            d = self%pd(i,j)%d
        else
            write(*,'(a)') 'WARNING: distance:', i, j, 'does not exist!'
        endif
    end function
    
    ! MODIFIERS
    
    !>  \brief  builds a table of distances given self and a distance function
    subroutine build( self, distfun )
        !$ use omp_lib
        !$ use omp_lib_kinds
        class(pair_dtab), intent(inout) :: self
        interface
            function distfun( i, j ) result( d )
                integer, intent(in) :: i, j
                real :: d
            end function
        end interface
        integer :: ia, ib
        real    :: dist
        ! first we allocate the table so that we immediately see what the memory usage is
        write(*,'(A)') '>>> ALLOCATING DISTANCE PD TABLE'
        !$omp parallel default(shared) private(ib) proc_bind(close)
        do ia=1,self%N-1
            !$omp do schedule(static) 
            do ib=ia+1,self%N
                call self%alloc_pair_d(ia,ib)    
            end do
            !$omp end do nowait
        end do
        ! then we fill it
        !$omp single
        write(*,'(A)') '>>> BUILDING DISTANCE PD TABLE'
        !$omp end single nowait
        do ia=1,self%N-1
            !$omp do schedule(static) 
            do ib=ia+1,self%N
                dist = distfun(ia,ib)
                call self%set_pair_d(ia, ib, distfun(ia,ib))
            end do
            !$omp end do nowait
        end do
        !$omp end parallel
    end subroutine
    
    !>  \brief  is for clearing the table
    subroutine clear( self )
        class(pair_dtab), intent(inout) :: self
        integer :: i, j
        if( self%existence )then
            do i=1,self%N-1
                do j=i+1,self%N
                    if( associated(self%pd(i,j)%d) )then
                        deallocate(self%pd(i,j)%d)
                        nullify(self%pd(i,j)%d, self%pd(j,i)%d) ! mirror symmetric nullification
                    endif
                end do
            end do
            self%minpos = 0
            self%minloc = 0
        endif
    end subroutine
    
    !>  \brief  is for allocating a pair dist
    subroutine alloc_pair_d( self, i, j )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: i, j
        integer                         :: alloc_stat
        if( j == i )then
            return
        else
            allocate( self%pd(i,j)%d, stat=alloc_stat ) ! make space for the new distance
            call alloc_err('In: alloc_pair_d, module: simple_pair_dtab.f90', alloc_stat)
        endif
    end subroutine

    !>  \brief  is for merging two pairs using the group average method
    subroutine merge_pair_d( self, i, j )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: i, j
        integer                         :: k
        real                            :: x
        ! merge
        do k=1,self%N
            if( k /= i .and. k /= j ) then
                if( associated(self%pd(i,k)%d) .and. associated(self%pd(j,k)%d) )then
                    self%pd(i,k)%d = 0.5*(self%pd(i,k)%d+self%pd(j,k)%d)
                    self%pd(k,i)%d = self%pd(i,k)%d
                endif
            endif
            if( pair_d_exist(self,j,k) )then
                deallocate( self%pd(j,k)%d )
                nullify( self%pd(j,k)%d, self%pd(k,j)%d ) ! mirror symmetric nullification
            endif
        end do
        ! update nearest neighs
        call update_nn( self, i )
        ! remove j from nearest neighs
        self%minvals(j) = huge(x)
        ! translate the positions in neigh
        do k=1,self%N
            if(self%minpos(k) == j) self%minpos(k) = i
        end do
        ! find new best 
        call update_best( self )
    end subroutine
    
    !>  \brief  is for removing a pair w
    subroutine remove_pair_d( self, i, j )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: i, j
        if( pair_d_exist(self,i,j) )then
            deallocate( self%pd(i,j)%d )
            nullify( self%pd(i,j)%d, self%pd(j,i)%d ) ! mirror symmetric nullification
        endif
        call update_nn( self, i )
        call update_nn( self, j )
        call update_best( self )
    end subroutine
    
    !>  \brief  is for providing a linear ordering according to the pairwise distances
    function linear_ordering( self ) result( order_out )
        class(pair_dtab), intent(inout) :: self
        integer, allocatable :: order(:), order_out(:)
        integer :: i, j, n, k, l, inds(2), oldj, oldk, cnt, alloc_stat
        real    :: d, dj, dk
        n = self%get_nobjs()
        allocate(order(-n:n), order_out(n), stat=alloc_stat)
        call alloc_err('order; simple_paird_tab', alloc_stat)
        order = 0
        call self%get_dmin_pair(k, l, d)
        order(1) = k
        order(2) = l
        inds(1)  = 1
        inds(2)  = 2
        call self%remove_pair_d(k, l)
        do i=3,n
            oldj = order(inds(1))
            oldk = order(inds(2))
            call self%get_dmin(oldj, j, dj)
            call self%get_dmin(oldk, k, dk)
            if( dj <= dk )then
                inds(1) = inds(1)-1
                order(inds(1)) = j
                call self%remove_pair_d(oldj, j)
            else
                inds(2) = inds(2)+1
                order(inds(2)) = k
                call self%remove_pair_d(oldk, k)
            endif
        end do
        cnt = 0
        do i=-n,n
            if( order(i) > 0 )then
                cnt = cnt+1
                order_out(cnt) = order(i)
            endif
        end do
        deallocate(order)
    end function
    
    ! PRIVATE STUFF
    
    !>  \brief  is for updating nearest neighs
    subroutine update_nn( self, i )
        class(pair_dtab), intent(inout) :: self
        integer, intent(in)             :: i
        integer                         :: k, loc(1)
        real                            :: x
        do k=1,self%N
            if( associated(self%pd(i,k)%d) )then
                self%row(k) = self%pd(i,k)%d
            else
                self%row(k) = huge(x)
            endif 
        end do
        loc = minloc(self%row)
        self%minpos(i) = loc(1)
        if( associated(self%pd(i,self%minpos(i))%d) )then
            self%minvals(i) = self%pd(i,self%minpos(i))%d
        else
            self%minvals(i) = huge(x)
        endif
    end subroutine
    
    !>  \brief  is for finding new best
    subroutine update_best( self )
        class(pair_dtab), intent(inout) :: self
        integer                         :: loc(1)
        loc            = minloc(self%minvals)
        self%minval    = self%minvals(loc(1))
        self%minloc(1) = loc(1)
        self%minloc(2) = self%minpos(self%minloc(1))
    end subroutine
    
    ! DESTRUCTOR
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(pair_dtab), intent(inout) :: self
        if( self%existence )then
            deallocate(self%pd, self%minvals, self%minpos, self%row)
            self%existence = .false.
        endif
    end subroutine

end module simple_pair_dtab
