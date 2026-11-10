!@descr: indexed tree topology for top-down splits and bottom-up joins
module simple_indexed_tree
use simple_error, only: simple_exception
implicit none
private
public :: indexed_tree
#include "simple_local_flags.inc"

! IDs are stable, one-based array indices until new/kill; zero means no link.
! Multiple roots are allowed while building a hierarchy. Payloads belong to callers.
type :: indexed_tree
    private
    integer :: nnodes = 0
    integer, allocatable :: parent(:), first_child(:), last_child(:), next_sibling(:)
  contains
    procedure :: new
    procedure :: kill
    procedure :: add_root
    procedure :: add_child
    procedure :: join_roots
    procedure :: get_nnodes
    procedure :: get_parent
    procedure :: get_children
    procedure :: get_roots
    procedure :: get_leaves
    procedure, private :: check_node
    procedure, private :: reserve
    procedure, private :: link_child
end type indexed_tree

contains

    ! Start with nroots isolated nodes, numbered 1:nroots; zero makes an empty tree.
    subroutine new( self, nroots )
        class(indexed_tree), intent(inout) :: self
        integer,             intent(in)    :: nroots
        if( nroots < 0 ) THROW_HARD('indexed_tree: negative root count; new')
        call self%kill
        call self%reserve(nroots)
        self%nnodes = nroots
    end subroutine new

    subroutine kill( self )
        class(indexed_tree), intent(inout) :: self
        if( allocated(self%parent)       ) deallocate(self%parent)
        if( allocated(self%first_child)  ) deallocate(self%first_child)
        if( allocated(self%last_child)   ) deallocate(self%last_child)
        if( allocated(self%next_sibling) ) deallocate(self%next_sibling)
        self%nnodes = 0
    end subroutine kill

    integer function add_root( self ) result( inode )
        class(indexed_tree), intent(inout) :: self
        if( self%nnodes == huge(self%nnodes) ) THROW_HARD('indexed_tree: node index overflow; add_root')
        inode = self%nnodes + 1
        call self%reserve(inode)
        self%nnodes = inode
    end function add_root

    ! Append a new leaf to an existing node, preserving child insertion order.
    integer function add_child( self, iparent ) result( inode )
        class(indexed_tree), intent(inout) :: self
        integer,             intent(in)    :: iparent
        call self%check_node(iparent)
        inode = self%add_root()
        call self%link_child(iparent, inode)
    end function add_child

    ! Add a parent above at least two distinct roots, in the supplied order.
    ! Requiring roots keeps existing subtrees intact and prevents cycles.
    integer function join_roots( self, roots ) result( inode )
        class(indexed_tree), intent(inout) :: self
        integer,             intent(in)    :: roots(:)
        integer :: i
        if( size(roots) < 2 ) THROW_HARD('indexed_tree: need at least two roots; join_roots')
        do i = 1, size(roots)
            call self%check_node(roots(i))
            if( self%parent(roots(i)) /= 0 ) THROW_HARD('indexed_tree: node already has a parent; join_roots')
            if( any(roots(:i-1) == roots(i)) ) THROW_HARD('indexed_tree: duplicate root; join_roots')
        end do
        ! Validate the whole list before changing the topology.
        inode = self%add_root()
        do i = 1, size(roots)
            call self%link_child(inode, roots(i))
        end do
    end function join_roots

    pure integer function get_nnodes( self ) result( n )
        class(indexed_tree), intent(in) :: self
        n = self%nnodes
    end function get_nnodes

    integer function get_parent( self, inode ) result( iparent )
        class(indexed_tree), intent(in) :: self
        integer,             intent(in) :: inode
        call self%check_node(inode)
        iparent = self%parent(inode)
    end function get_parent

    ! Return a copy in insertion order; a leaf has an allocated, empty result.
    function get_children( self, inode ) result( children )
        class(indexed_tree), intent(in) :: self
        integer,             intent(in) :: inode
        integer, allocatable :: children(:)
        integer :: child, n, i
        call self%check_node(inode)
        n     = 0
        child = self%first_child(inode)
        do while( child /= 0 )
            n     = n + 1
            child = self%next_sibling(child)
        end do
        allocate(children(n))
        child = self%first_child(inode)
        do i = 1, n
            children(i) = child
            child       = self%next_sibling(child)
        end do
    end function get_children

    ! Root and leaf lists use increasing node IDs, independent of child order.
    function get_roots( self ) result( roots )
        class(indexed_tree), intent(in) :: self
        integer, allocatable :: roots(:)
        integer :: i
        if( self%nnodes == 0 )then
            allocate(roots(0))
            return
        endif
        roots = pack([(i, i=1,self%nnodes)], self%parent(:self%nnodes) == 0)
    end function get_roots

    function get_leaves( self ) result( leaves )
        class(indexed_tree), intent(in) :: self
        integer, allocatable :: leaves(:)
        integer :: i
        if( self%nnodes == 0 )then
            allocate(leaves(0))
            return
        endif
        leaves = pack([(i, i=1,self%nnodes)], self%first_child(:self%nnodes) == 0)
    end function get_leaves

    subroutine check_node( self, inode )
        class(indexed_tree), intent(in) :: self
        integer,             intent(in) :: inode
        if( inode < 1 .or. inode > self%nnodes ) THROW_HARD('indexed_tree: node index out of range')
    end subroutine check_node

    ! Called only with a validated parent and an unattached root.
    subroutine link_child( self, iparent, child )
        class(indexed_tree), intent(inout) :: self
        integer,             intent(in)    :: iparent, child
        integer :: last
        last = self%last_child(iparent)
        if( last == 0 )then
            self%first_child(iparent) = child
        else
            self%next_sibling(last) = child
        endif
        self%last_child(iparent) = child
        self%parent(child)       = iparent
    end subroutine link_child

    subroutine reserve( self, needed )
        class(indexed_tree), intent(inout) :: self
        integer,             intent(in)    :: needed
        integer :: capacity, newcap
        capacity = 0
        if( allocated(self%parent) ) capacity = size(self%parent)
        if( needed <= capacity ) return
        newcap = max(16, needed)
        if( capacity <= shiftr(huge(capacity), 1) ) newcap = max(newcap, 2 * capacity)
        call grow(self%parent)
        call grow(self%first_child)
        call grow(self%last_child)
        call grow(self%next_sibling)
      contains
        subroutine grow( values )
            integer, allocatable, intent(inout) :: values(:)
            integer, allocatable :: tmp(:)
            allocate(tmp(newcap), source=0)
            if( self%nnodes > 0 ) tmp(:self%nnodes) = values(:self%nnodes)
            call move_alloc(tmp, values)
        end subroutine grow
    end subroutine reserve

end module simple_indexed_tree
