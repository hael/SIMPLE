!==Class simple_sll
!
! simple_sll is a runtime polymorphic singly linked list class. 
! The code is distributed with the hope that it will be useful, but 
! _WITHOUT_ _ANY_ _WARRANTY_. Redistribution or modification is
! regulated by the GNU General Public License. 
! *Author:* Dominika Elmlund, 2009-05-25.
! 
!==Changes are documented below
!
!* incorporated in the _SIMPLE_ library, HE 2009-06-25
!
module simple_sll
use simple_defs
use simple_arr, only: arr
implicit none

public :: sll
private

type sll_node
! contains the array _content_ and a pointer _next_ to the nextcoming node
    type(arr) :: content
    type(sll_node), pointer :: next=>null()
end type sll_node

type sll
! contains the list. In this implementation I have separated the head of the 
! list from the rest of the list (Donald Knuth-style) 
    private
    integer                 :: list_size=0
    type(sll_node), pointer :: head=>null()
  contains
    procedure :: new
    procedure :: add
    procedure :: get
    procedure :: set 
    procedure :: del
    procedure, private :: assign
    generic :: assignment(=) => assign
    procedure :: append
    procedure :: display
    procedure :: size
    procedure :: kill
end type sll

interface sll
    module procedure constructor
end interface sll

contains

    !>  \brief  is a constructor that allocates the head of the list 
    !! and nullifies its pointer to the nextcoming node
    function constructor() result(self)
        type(sll) :: self
        call self%new
    end function constructor

    !>  \brief  is a constructor that allocates the head of the list 
    !! and nullifies its pointer to the nextcoming node
    subroutine new(self)
        class(sll), intent(inout) :: self
        if( associated(self%head) ) call self%kill
        allocate(self%head)     ! allocate memory for the object
        nullify(self%head%next) ! start with an empty list
    end subroutine new
    
    !>  \brief  does polymorphic addition of a node to the end of 
    !! the singly linked list and updates the list size
    subroutine add( self, iarr, rarr )
        class(sll),        intent(inout) :: self
        integer, optional, intent(in)    :: iarr(:)
        real, optional,    intent(in)    :: rarr(:)
        type(sll_node), pointer          :: prev, curr
        ! initialization, begin at the 0:th position
        prev => self%head 
        curr => prev%next
        do while( associated(curr) )! find location to insert new node
          prev => curr
          curr => curr%next
        end do
        allocate( curr ) ! insert it at the end of the list
        if( present(iarr) ) curr%content = iarr
        if( present(rarr) ) curr%content = rarr
        self%list_size = self%list_size+1  
        ! Redirect the old last
        prev%next => curr
    end subroutine add
    
    !>  \brief  is a polymorphic getter
    subroutine get( self, pos, iarr, rarr )
        class(sll),                     intent(in)  :: self
        integer,                        intent(in)  :: pos
        integer, allocatable, optional, intent(out) :: iarr(:)
        real,    allocatable, optional, intent(out) :: rarr(:)
        type(sll_node), pointer :: curr
        integer                 :: counter
        if ( pos < 1 .or. pos > self%list_size ) then
            write(*,*) 'Variable pos is out of range!'
            write(*,*) 'get; simple_sll'
            stop
        endif
        curr => self%head%next
        counter = 0
        do ! find location of node
            counter = counter+1
            if( counter == pos ) exit
            curr => curr%next
        end do
        if( present(iarr) )then
            if( allocated(iarr) ) deallocate(iarr)
            iarr = curr%content%iget()
        endif
        if( present(rarr) )then
            if( allocated(rarr) ) deallocate(rarr)
            rarr = curr%content%rget()
        endif
    end subroutine get
    
    !>  \brief  is a polymorphic setter
    subroutine set( self, pos, iarr, rarr )
        class(sll),        intent(in) :: self
        integer,           intent(in) :: pos
        integer, optional, intent(in) :: iarr(:)
        real,    optional, intent(in) :: rarr(:)
        type(sll_node), pointer       :: curr
        integer                       :: counter
        if ( pos < 1 .or. pos > self%list_size ) then
          write(*,*) 'Variable pos is out of range!'
          write(*,*) 'In: set_sll_node, module: simple_sll.f90'
          stop
        endif
        curr => self%head%next
        counter = 0
        do ! find location of node
            counter = counter+1
            if( counter == pos ) exit
            curr => curr%next
        end do
        if( present(iarr) ) curr%content = iarr
        if( present(rarr) ) curr%content = rarr
    end subroutine set
    
    !>  \brief  deallocates a sll node and redirects the pointer to next node
    subroutine del( self, pos )
        class(sll), intent(inout) :: self
        integer,    intent(in)    :: pos
        type(sll_node), pointer   :: prev, curr
        integer                   :: counter
        if ( pos < 1 .or. pos > self%list_size ) then
            write(*,*) 'Variable pos is out of range!'
            write(*,*) 'get; simple_sll'
            stop
        endif
        ! initialization, begin at the 0:th position
        prev => self%head
        curr => prev%next
        if( .not. associated(curr) ) then ! end of list
            self%list_size = 0
            return
        endif     
        counter = 0
        do ! find node to delete
          counter = counter+1
          if( pos == counter ) then
            exit
          else ! move to the next node of the list
            prev => curr
            curr => curr%next
          endif
        end do
        ! delete the node
        if( associated( curr%next ) ) then       
          prev%next => curr%next          ! redirect pointer
          call curr%content%kill          ! free space for the content    
          deallocate( curr )              ! free space for node
          nullify( curr )
        else
          call curr%content%kill          ! free space for the list object
          deallocate( curr )              ! free space for node
          nullify( curr, prev%next )
        endif
        self%list_size = self%list_size-1 ! update the list size
    end subroutine del
    
    !>  \brief  clones a sll
    subroutine assign( self, self_in )
        class(sll), intent(in)    :: self_in
        class(sll), intent(inout) :: self
        call kill( self )
        call self%new
        self%list_size = self_in%list_size
        ! make resulting list a replica of self_in
        self%head%next => self_in%head%next
    end subroutine assign
   
    !>  \brief is joining two singly linked lists together. The first node of 
    !! the second list (_self2_) is joined with the last node of the first list
    !! (_self1_). The input lists do not exist anymore after this operation.
    function append( self1, self2 ) result( self )
        class(sll), intent(inout)  :: self1
        class(sll), intent(inout)  :: self2
        type(sll)                  :: self
        type(sll_node), pointer    :: prev, curr
        ! Make resulting list
        call self%new
        self%list_size = self1%list_size + self2%list_size
        ! make resulting list a replica of self1
        self%head%next => self1%head%next
        ! remove list 1
        deallocate( self1%head )
        nullify( self1%head )
        self1%list_size = 0
        ! do the choka choka
        prev => self%head
        curr => prev%next
        do while( associated(curr) )
          prev => curr
          curr => curr%next
        end do
        prev%next => self2%head%next
        ! remove list 2
        deallocate( self2%head )
        nullify( self2%head )
        self2%list_size = 0
    end function append
  
    !>  \brief  is for printing
    subroutine display( self )
        class(sll), intent(in)  :: self
        type(sll_node), pointer :: curr
        integer                 :: pos
        ! initialization, begin at the 1:th position 
        curr => self%head%next
        pos  = 0 ! with pos set to 0 
        do while( associated( curr ) )
          pos = pos+1
          write( *,* ) 'DATA IN NODE: ', pos
          call curr%content%display
          curr => curr%next
        end do
    end subroutine display
    
    !>  \brief  returns the size of the list
    pure function size( self ) result( list_size )
        class(sll), intent(in) :: self
        integer                :: list_size
        list_size = self%list_size   
    end function size
    
    !>  \brief  is a destructor
    subroutine kill( self )
        class(sll), intent(inout) :: self
        integer                   :: i
        if( self%list_size >= 0 )then
            do i=1,self%list_size
                call self%del(1)
            end do
        endif
        if( associated(self%head) )then
            deallocate(self%head)
            nullify(self%head)
        endif
    end subroutine kill

end module simple_sll
