!@descr: simple binary tree structure for hierarchical clustering results
module simple_binary_tree
implicit none
private
public :: bt_node, binary_tree

type :: bt_node
    integer :: left_idx   = 0
    integer :: right_idx  = 0
    integer :: parent_idx = 0
    integer, allocatable :: subset(:)   ! LOCAL indices (1..nref)
    integer :: level    = 0
    integer :: ref_idx  = 0             ! GLOBAL ref id (medoid for internal nodes)
    integer :: node_idx = 0             ! 1..n_nodes
end type bt_node

abstract interface
    subroutine node_visitor(node)
        import :: bt_node
        type(bt_node), intent(in) :: node
    end subroutine node_visitor
end interface

type :: binary_tree
    private
    type(bt_node), allocatable :: nodes(:)
    integer,      allocatable :: refs(:)     ! LOCAL->GLOBAL mapping (length nref)
    integer :: root_idx = 0
    integer :: nref     = 0
    logical :: exists   = .false.
contains
    procedure          :: kill
    procedure          :: n_nodes
    procedure          :: get_root_idx
    procedure          :: get_node
    procedure          :: is_leaf
    procedure          :: get_subset_size
    procedure          :: get_node_ref
    procedure          :: get_children_idx
    procedure          :: get_children_ref
    procedure          :: get_nref
    procedure          :: local_to_global_ref
    procedure          :: get_medoid
    procedure          :: get_height
    procedure, private :: height_recursive
    ! build from hierarchical clustering merge matrix + original reference ids
    procedure          :: build_from_hclust
    ! traversals
    procedure          :: traverse_preorder
    procedure          :: traverse_inorder
    procedure          :: traverse_postorder
    procedure          :: traverse_levelorder
end type binary_tree

contains

    subroutine kill(self)
        class(binary_tree), intent(inout) :: self
        integer :: k
        if (.not. self%exists) return
        if (allocated(self%nodes)) then
            do k = 1, size(self%nodes)
                if (allocated(self%nodes(k)%subset)) deallocate(self%nodes(k)%subset)
                self%nodes(k)%left_idx   = 0
                self%nodes(k)%right_idx  = 0
                self%nodes(k)%parent_idx = 0
                self%nodes(k)%level      = 0
                self%nodes(k)%ref_idx    = 0
                self%nodes(k)%node_idx   = 0
            end do
            deallocate(self%nodes)
        end if
        if (allocated(self%refs)) deallocate(self%refs)
        self%root_idx = 0
        self%nref     = 0
        self%exists   = .false.
    end subroutine kill

    pure integer function n_nodes(self)
        class(binary_tree), intent(in) :: self
        if (allocated(self%nodes)) then
            n_nodes = size(self%nodes)
        else
            n_nodes = 0
        end if
    end function n_nodes

    pure integer function get_root_idx(self)
        class(binary_tree), intent(in) :: self
        get_root_idx = self%root_idx
    end function get_root_idx

    pure integer function get_nref(self)
        class(binary_tree), intent(in) :: self
        get_nref = self%nref
    end function get_nref

    pure function get_node(self, inode) result(node)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: inode
        type(bt_node) :: node
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

    pure integer function get_subset_size(self, inode)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: inode
        get_subset_size = 0
        if (.not. allocated(self%nodes)) return
        if (inode < 1 .or. inode > size(self%nodes)) return
        if (.not. allocated(self%nodes(inode)%subset)) return
        get_subset_size = size(self%nodes(inode)%subset)
    end function get_subset_size

    pure integer function get_node_ref(self, inode) result(ref)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: inode
        ref = 0
        if (.not. allocated(self%nodes)) return
        if (inode < 1 .or. inode > size(self%nodes)) return
        ref = self%nodes(inode)%ref_idx
    end function get_node_ref

    pure subroutine get_children_idx(self, inode, left_idx, right_idx)
        class(binary_tree), intent(in)  :: self
        integer,            intent(in)  :: inode
        integer,            intent(out) :: left_idx, right_idx
        left_idx  = 0
        right_idx = 0
        if (.not. allocated(self%nodes)) return
        if (inode < 1 .or. inode > size(self%nodes)) return
        left_idx  = self%nodes(inode)%left_idx
        right_idx = self%nodes(inode)%right_idx
    end subroutine get_children_idx

    pure subroutine get_children_ref(self, inode, left_ref, right_ref)
        class(binary_tree), intent(in)  :: self
        integer,            intent(in)  :: inode
        integer,            intent(out) :: left_ref, right_ref
        integer :: lidx, ridx
        left_ref  = 0
        right_ref = 0
        call self%get_children_idx(inode, lidx, ridx)
        if (lidx /= 0) left_ref  = self%nodes(lidx)%ref_idx
        if (ridx /= 0) right_ref = self%nodes(ridx)%ref_idx
    end subroutine get_children_ref

    pure integer function local_to_global_ref(self, local_k) result(global_ref)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: local_k
        global_ref = 0
        if (.not. allocated(self%refs)) return
        if (local_k < 1 .or. local_k > size(self%refs)) return
        global_ref = self%refs(local_k)
    end function local_to_global_ref

    pure integer function get_medoid(self) result(ref)
        class(binary_tree), intent(in) :: self
        ref = 0
        if (.not. self%exists) return
        if (self%root_idx == 0) return
        ref = self%nodes(self%root_idx)%ref_idx
    end function get_medoid

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
        l = self%nodes(idx)%left_idx
        r = self%nodes(idx)%right_idx
        hl = height_recursive(self, l)
        hr = height_recursive(self, r)
        h = 1 + max(hl, hr)   ! "levels" counted as nodes along path
    end function height_recursive

    ! --------------------------
    ! Build from merge matrix
    ! --------------------------
    subroutine build_from_hclust(self, merge_mat, refs, dist_mat)
        class(binary_tree), intent(inout) :: self
        integer,            intent(in)    :: merge_mat(:,:)   ! (2, nref-1)
        integer,            intent(in)    :: refs(:)          ! length nref (GLOBAL ids)
        real,               intent(in)    :: dist_mat(:,:)    ! (nref,nref) LOCAL to refs ordering
        integer :: nref, n_total, k, s, l, r, p, m, i, j
        integer :: best_local
        integer, allocatable :: tmp(:)
        real :: best_sum, sum
        call self%kill()
        nref    = size(refs)
        self%nref = nref
        n_total = 2*nref - 1
        allocate(self%refs(nref))
        self%refs = refs
        allocate(self%nodes(n_total))
        ! init all nodes
        do k = 1, n_total
            self%nodes(k)%node_idx   = k
            self%nodes(k)%left_idx   = 0
            self%nodes(k)%right_idx  = 0
            self%nodes(k)%parent_idx = 0
            self%nodes(k)%level      = 0
            self%nodes(k)%ref_idx    = 0
            if (allocated(self%nodes(k)%subset)) deallocate(self%nodes(k)%subset)
        end do
        ! leaves 1..nref
        do k = 1, nref
            self%nodes(k)%ref_idx = refs(k)     ! GLOBAL ref id
            allocate(self%nodes(k)%subset(1))
            self%nodes(k)%subset = [k]          ! LOCAL index (1..nref)
        end do
        ! internal nodes nref+1 .. n_total
        m = nref
        do s = 1, nref - 1
            l = merge_mat(1, s)
            r = merge_mat(2, s)
            p = m + s
            self%nodes(p)%left_idx   = l
            self%nodes(p)%right_idx  = r
            self%nodes(l)%parent_idx = p
            self%nodes(r)%parent_idx = p
            self%nodes(p)%level      = s
            ! union subset (still LOCAL indices)
            allocate(tmp(size(self%nodes(l)%subset) + size(self%nodes(r)%subset)))
            tmp = [self%nodes(l)%subset, self%nodes(r)%subset]
            call move_alloc(tmp, self%nodes(p)%subset)
            ! compute medoid using LOCAL dist_mat (indexed by LOCAL indices)
            best_local = -1
            best_sum   = huge(1.0)
            do i = 1, size(self%nodes(p)%subset)
                sum = 0.0
                do j = 1, size(self%nodes(p)%subset)
                    sum = sum + dist_mat(self%nodes(p)%subset(i), self%nodes(p)%subset(j))
                end do
                if (sum < best_sum) then
                    best_sum   = sum
                    best_local = self%nodes(p)%subset(i)
                end if
            end do
            self%nodes(p)%ref_idx = refs(best_local)   ! store GLOBAL id
        end do
        self%root_idx = n_total
        self%exists   = .true.
    end subroutine build_from_hclust

    ! --------------------------
   ! Traversals (visitor-based)
   ! --------------------------
    subroutine traverse_preorder(self, visit)
        class(binary_tree), intent(in) :: self
        procedure(node_visitor) :: visit
        call preorder_rec(self, self%root_idx, visit)
    end subroutine traverse_preorder

    recursive subroutine preorder_rec(self, idx, visit)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: idx
        procedure(node_visitor) :: visit
        if (idx == 0) return
        call visit(self%nodes(idx))
        call preorder_rec(self, self%nodes(idx)%left_idx,  visit)
        call preorder_rec(self, self%nodes(idx)%right_idx, visit)
    end subroutine preorder_rec

    subroutine traverse_inorder(self, visit)
        class(binary_tree), intent(in) :: self
        procedure(node_visitor) :: visit
        call inorder_rec(self, self%root_idx, visit)
    end subroutine traverse_inorder

    recursive subroutine inorder_rec(self, idx, visit)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: idx
        procedure(node_visitor) :: visit
        if (idx == 0) return
        call inorder_rec(self, self%nodes(idx)%left_idx, visit)
        call visit(self%nodes(idx))
        call inorder_rec(self, self%nodes(idx)%right_idx, visit)
    end subroutine inorder_rec

    subroutine traverse_postorder(self, visit)
        class(binary_tree), intent(in) :: self
        procedure(node_visitor) :: visit
        call postorder_rec(self, self%root_idx, visit)
    end subroutine traverse_postorder

    recursive subroutine postorder_rec(self, idx, visit)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: idx
        procedure(node_visitor) :: visit
        if (idx == 0) return
        call postorder_rec(self, self%nodes(idx)%left_idx,  visit)
        call postorder_rec(self, self%nodes(idx)%right_idx, visit)
        call visit(self%nodes(idx))
    end subroutine postorder_rec

    subroutine traverse_levelorder(self, visit)
        class(binary_tree), intent(in) :: self
        procedure(node_visitor) :: visit
        integer, allocatable :: q(:)
        integer :: head, tail, cur, l, r
        if (self%root_idx == 0) return
        allocate(q(self%n_nodes()))
        head = 1; tail = 1
        q(1) = self%root_idx
        do while (head <= tail)
            cur = q(head); head = head + 1
            call visit(self%nodes(cur))
            l = self%nodes(cur)%left_idx
            r = self%nodes(cur)%right_idx
            if (l /= 0) then
                tail = tail + 1; q(tail) = l
            end if
            if (r /= 0) then
                tail = tail + 1; q(tail) = r
            end if
        end do
        deallocate(q)
    end subroutine traverse_levelorder

    pure integer function find_node_by_ref(self, ref_idx) result(found)
        class(binary_tree), intent(in) :: self
        integer,            intent(in) :: ref_idx
        integer :: k
        found = 0
        if (.not. allocated(self%nodes)) return
        do k = 1, size(self%nodes)
            if (self%nodes(k)%ref_idx == ref_idx) then
                found = k
                return
            end if
        end do
    end function find_node_by_ref

end module simple_binary_tree
