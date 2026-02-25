!@descr: simple binary tree structure for hierarchical clustering results
module simple_binary_tree
use simple_error, only: simple_exception
implicit none
private
public :: bt_node, binary_tree

#include "simple_local_flags.inc"

type :: bt_node
    integer :: left_idx   = 0
    integer :: right_idx  = 0
    integer :: parent_idx = 0
    integer :: level      = 0
    integer :: ref_idx    = 0   ! GLOBAL ref id (medoid for internal nodes)
    integer :: node_idx   = 0   ! 1..n_nodes (2*nrefs-1)
end type bt_node

type :: binary_tree
    private
    type(bt_node), allocatable :: nodes(:)
    integer :: root_idx = 0
    integer :: nrefs    = 0
    logical :: exists   = .false.
contains
    procedure :: kill
    procedure :: get_n_nodes
    procedure :: get_nrefs
    procedure :: get_root_node
    procedure :: get_node
    procedure :: is_leaf
    procedure :: get_height
    procedure, private :: height_recursive
    procedure :: build_from_hclust
end type binary_tree

contains

    subroutine kill(self)
        class(binary_tree), intent(inout) :: self
        integer :: k
        if (.not. self%exists) return
        if (allocated(self%nodes)) then
            do k = 1, size(self%nodes)
                self%nodes(k)%left_idx   = 0
                self%nodes(k)%right_idx  = 0
                self%nodes(k)%parent_idx = 0
                self%nodes(k)%level      = 0
                self%nodes(k)%ref_idx    = 0
                self%nodes(k)%node_idx   = 0
            end do
            deallocate(self%nodes)
        end if

        self%root_idx = 0
        self%nrefs    = 0
        self%exists   = .false.
    end subroutine kill

    pure integer function get_n_nodes(self) result(n)
        class(binary_tree), intent(in) :: self
        if (allocated(self%nodes)) then
            n = size(self%nodes)
        else
            n = 0
        end if
    end function get_n_nodes

    pure integer function get_nrefs(self) result(n)
        class(binary_tree), intent(in) :: self
        n = self%nrefs
    end function get_nrefs

    pure function get_root_node(self) result(root_node)
        class(binary_tree), intent(in) :: self
        type(bt_node) :: root_node
        root_node%left_idx   = 0
        root_node%right_idx  = 0
        root_node%parent_idx = 0
        root_node%level      = 0
        root_node%ref_idx    = 0
        root_node%node_idx   = 0
        if (.not. allocated(self%nodes)) return
        if (self%root_idx < 1 .or. self%root_idx > size(self%nodes)) return
        root_node = self%nodes(self%root_idx)
    end function get_root_node

    pure function get_node(self, inode) result(node)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: inode
        type(bt_node) :: node
        node%left_idx   = 0
        node%right_idx  = 0
        node%parent_idx = 0
        node%level      = 0
        node%ref_idx    = 0
        node%node_idx   = 0
        if (.not. allocated(self%nodes)) return
        if (inode < 1 .or. inode > size(self%nodes)) return
        node = self%nodes(inode)
    end function get_node

    pure logical function is_leaf(self, inode)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: inode
        is_leaf = .false.
        if (.not. allocated(self%nodes)) return
        if (inode < 1 .or. inode > size(self%nodes)) return
        is_leaf = (self%nodes(inode)%left_idx == 0 .and. self%nodes(inode)%right_idx == 0)
    end function is_leaf

    pure integer function get_height(self) result(h)
        class(binary_tree), intent(in) :: self
        h = 0
        if (.not. self%exists) return
        if (.not. allocated(self%nodes)) return
        if (self%root_idx == 0) return
        h = height_recursive(self, self%root_idx)
    end function get_height

    pure recursive integer function height_recursive(self, idx) result(h)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: idx
        integer :: hl, hr, l, r
        if (idx == 0) then
            h = 0
            return
        end if
        l  = self%nodes(idx)%left_idx
        r  = self%nodes(idx)%right_idx
        hl = height_recursive(self, l)
        hr = height_recursive(self, r)
        h  = 1 + max(hl, hr)
    end function height_recursive

    ! --------------------------
    ! Build from merge matrix
    ! Computes medoids using temporary membership lists (not stored in bt_node).
    ! --------------------------
    subroutine build_from_hclust(self, merge_mat, refs, dist_mat)
        class(binary_tree), intent(inout) :: self
        integer,            intent(in)    :: merge_mat(:,:)   ! (2, nrefs-1)
        integer,            intent(in)    :: refs(:)          ! length nrefs (GLOBAL ids)
        real,               intent(in)    :: dist_mat(:,:)    ! (nrefs,nrefs) LOCAL to refs ordering
        type :: int_list
            integer, allocatable :: v(:)    ! LOCAL indices (1..nrefs)
        end type int_list
        type(int_list), allocatable :: members(:)
        integer,        allocatable :: tmp(:)
        integer :: nrefs, n_total, k, s, l, r, p, m, i, j, best_local
        real    :: best_sum, sum
        call self%kill()
        nrefs = size(refs)
        if (nrefs < 1) then
            THROW_HARD("build_from_hclust: empty refs")
        end if
        if (size(merge_mat,1) /= 2 .or. size(merge_mat,2) /= nrefs-1) then
            THROW_HARD("build_from_hclust: merge_mat must be (2, nrefs-1)")
        end if
        if (size(dist_mat,1) /= nrefs .or. size(dist_mat,2) /= nrefs) then
            THROW_HARD("build_from_hclust: dist_mat must be (nrefs,nrefs)")
        end if
        self%nrefs = nrefs
        n_total    = 2*nrefs - 1
        self%root_idx = n_total
        allocate(self%nodes(n_total), members(n_total))
        ! init nodes + empty member lists
        do k = 1, n_total
            self%nodes(k)%node_idx   = k
            self%nodes(k)%left_idx   = 0
            self%nodes(k)%right_idx  = 0
            self%nodes(k)%parent_idx = 0
            self%nodes(k)%level      = 0
            self%nodes(k)%ref_idx    = 0
        end do
        ! leaves: node k corresponds to local ref k
        do k = 1, nrefs
            self%nodes(k)%ref_idx = refs(k)   ! GLOBAL id
            allocate(members(k)%v(1))
            members(k)%v = [k]
        end do
        ! internal nodes: wire + compute medoid for each merged cluster
        m = nrefs
        do s = 1, nrefs - 1
            l = merge_mat(1, s)
            r = merge_mat(2, s)
            p = m + s
            self%nodes(p)%left_idx   = l
            self%nodes(p)%right_idx  = r
            self%nodes(l)%parent_idx = p
            self%nodes(r)%parent_idx = p
            self%nodes(p)%level      = s
            ! union membership
            allocate(tmp(size(members(l)%v) + size(members(r)%v)))
            tmp = [members(l)%v, members(r)%v]
            call move_alloc(tmp, members(p)%v)
            ! medoid from dist_mat over members(p)%v
            best_local = -1
            best_sum   = huge(1.0)
            do i = 1, size(members(p)%v)
                sum = 0.0
                do j = 1, size(members(p)%v)
                    sum = sum + dist_mat(members(p)%v(i), members(p)%v(j))
                end do
                if (sum < best_sum) then
                    best_sum   = sum
                    best_local = members(p)%v(i)
                end if
            end do
            if (best_local < 1 .or. best_local > nrefs) then
                THROW_HARD("build_from_hclust: computed medoid index out of range")
            end if
            self%nodes(p)%ref_idx = refs(best_local)
        end do
        ! free temporary membership lists
        do k = 1, n_total
            if (allocated(members(k)%v)) deallocate(members(k)%v)
        end do
        deallocate(members)
        self%exists = .true.
    end subroutine build_from_hclust

end module simple_binary_tree