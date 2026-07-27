!@descr: deterministic flat-array k-d tree and exact Euclidean k-nearest-neighbor queries
module simple_kd_tree
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
implicit none
private
#include "simple_local_flags.inc"

public :: kd_tree, knn_table

type :: knn_table
    integer :: n = 0
    integer :: k = 0
    integer, allocatable :: neighbor(:,:)
    real,    allocatable :: distance2(:,:)
contains
    procedure :: kill => kill_knn_table
end type knn_table

type :: kd_tree
    private
    integer :: ndim      = 0
    integer :: npoints   = 0
    integer :: nnodes    = 0
    integer :: root      = 0
    integer :: leaf_size = 32
    integer, allocatable :: permutation(:)
    integer, allocatable :: left(:), right(:)
    integer, allocatable :: first(:), count(:)
    integer, allocatable :: split_dim(:)
    real,    allocatable :: split_value(:)
contains
    procedure :: build       => build_kd_tree
    procedure :: query       => query_kd_tree
    procedure :: query_all   => query_all_kd_tree
    procedure :: get_nnodes  => get_kd_tree_nnodes
    procedure :: get_height  => get_kd_tree_height
    procedure :: kill        => kill_kd_tree
end type kd_tree

contains

    subroutine kill_knn_table( self )
        class(knn_table), intent(inout) :: self
        if( allocated(self%neighbor)  ) deallocate(self%neighbor)
        if( allocated(self%distance2) ) deallocate(self%distance2)
        self%n = 0
        self%k = 0
    end subroutine kill_knn_table

    subroutine build_kd_tree( self, features, leaf_size )
        class(kd_tree), intent(inout) :: self
        real,           intent(in)    :: features(:,:)
        integer, optional, intent(in) :: leaf_size
        integer :: i, initial_capacity
        call self%kill()
        self%ndim    = size(features,1)
        self%npoints = size(features,2)
        if( self%ndim < 1 .or. self%npoints < 1 ) THROW_HARD('kd_tree requires a nonempty feature table')
        if( .not.all(ieee_is_finite(features)) ) THROW_HARD('kd_tree feature table contains nonfinite values')
        self%leaf_size = 32
        if( present(leaf_size) ) self%leaf_size = max(1,leaf_size)
        allocate(self%permutation(self%npoints))
        self%permutation = [(i,i=1,self%npoints)]
        initial_capacity = max(16,4*((self%npoints+self%leaf_size-1)/self%leaf_size)+1)
        call reserve_nodes(self,initial_capacity)
        call build_node(self,features,1,self%npoints,self%root)
    end subroutine build_kd_tree

    recursive subroutine build_node( self, features, lo, hi, inode )
        class(kd_tree), intent(inout) :: self
        real,           intent(in)    :: features(:,:)
        integer,        intent(in)    :: lo,hi
        integer,        intent(out)   :: inode
        integer :: axis,mid,left_node,right_node
        inode = new_node(self)
        self%first(inode) = lo
        self%count(inode) = hi-lo+1
        if( self%count(inode) <= self%leaf_size ) return
        axis = widest_dimension(self,features,lo,hi)
        mid  = (lo+hi)/2
        call select_kth(self,features,axis,lo,hi,mid)
        self%split_dim(inode)   = axis
        self%split_value(inode) = features(axis,self%permutation(mid))
        call build_node(self,features,lo,mid,left_node)
        call build_node(self,features,mid+1,hi,right_node)
        self%left(inode)  = left_node
        self%right(inode) = right_node
    end subroutine build_node

    integer function widest_dimension( self, features, lo, hi ) result(axis)
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:)
        integer,        intent(in) :: lo,hi
        real(dp) :: xmin,xmax,span,best_span
        integer :: q,p
        axis = 1
        best_span = -1.0_dp
        do q=1,self%ndim
            xmin = real(features(q,self%permutation(lo)),dp)
            xmax = xmin
            do p=lo+1,hi
                xmin = min(xmin,real(features(q,self%permutation(p)),dp))
                xmax = max(xmax,real(features(q,self%permutation(p)),dp))
            end do
            span = xmax-xmin
            if( span > best_span )then
                best_span = span
                axis = q
            endif
        end do
    end function widest_dimension

    subroutine select_kth( self, features, axis, lo_in, hi_in, kth )
        class(kd_tree), intent(inout) :: self
        real,           intent(in)    :: features(:,:)
        integer,        intent(in)    :: axis,lo_in,hi_in,kth
        integer :: lo,hi,i,store,pivot_pos,pivot_id
        lo = lo_in
        hi = hi_in
        do while( lo < hi )
            pivot_pos = median_of_three_position(self,features,axis,lo,(lo+hi)/2,hi)
            pivot_id  = self%permutation(pivot_pos)
            call swap_int(self%permutation(pivot_pos),self%permutation(hi))
            store = lo
            do i=lo,hi-1
                if( point_less(features,axis,self%permutation(i),pivot_id) )then
                    call swap_int(self%permutation(store),self%permutation(i))
                    store = store+1
                endif
            end do
            call swap_int(self%permutation(store),self%permutation(hi))
            if( store == kth )then
                return
            else if( kth < store )then
                hi = store-1
            else
                lo = store+1
            endif
        end do
    end subroutine select_kth

    integer function median_of_three_position( self, features, axis, p1, p2, p3 ) result(pos)
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:)
        integer,        intent(in) :: axis,p1,p2,p3
        integer :: a,b,c
        a = self%permutation(p1)
        b = self%permutation(p2)
        c = self%permutation(p3)
        if( point_less(features,axis,a,b) )then
            if( point_less(features,axis,b,c) )then
                pos = p2
            else if( point_less(features,axis,a,c) )then
                pos = p3
            else
                pos = p1
            endif
        else
            if( point_less(features,axis,a,c) )then
                pos = p1
            else if( point_less(features,axis,b,c) )then
                pos = p3
            else
                pos = p2
            endif
        endif
    end function median_of_three_position

    pure logical function point_less( features, axis, i, j ) result(is_less)
        real,    intent(in) :: features(:,:)
        integer, intent(in) :: axis,i,j
        if( features(axis,i) < features(axis,j) )then
            is_less = .true.
        else if( features(axis,i) > features(axis,j) )then
            is_less = .false.
        else
            is_less = i < j
        endif
    end function point_less

    subroutine query_kd_tree( self, features, query, k, neighbors, distances2, exclude )
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:),query(:)
        integer,        intent(in) :: k
        integer, allocatable, intent(out) :: neighbors(:)
        real,    allocatable, intent(out) :: distances2(:)
        integer, optional, intent(in) :: exclude
        integer, allocatable :: heap_ids(:)
        real(dp), allocatable :: heap_d2(:)
        integer :: exclude_here,heap_size,visited
        call validate_query(self,features,query,k)
        exclude_here = 0
        if( present(exclude) ) exclude_here = exclude
        if( exclude_here < 0 .or. exclude_here > self%npoints ) THROW_HARD('kd_tree exclude index outside feature table')
        if( k > self%npoints-merge(1,0,exclude_here>0) ) THROW_HARD('kd_tree query requests too many neighbors')
        allocate(heap_ids(k),heap_d2(k),neighbors(k),distances2(k))
        heap_size = 0
        visited   = 0
        call query_node(self,features,query,exclude_here,self%root,k,heap_ids,heap_d2,heap_size,visited)
        if( heap_size /= k ) THROW_HARD('kd_tree query did not find the requested number of neighbors')
        call sort_neighbor_pairs(heap_ids,heap_d2,k)
        neighbors  = heap_ids
        distances2 = real(heap_d2)
        deallocate(heap_ids,heap_d2)
    end subroutine query_kd_tree

    subroutine query_all_kd_tree( self, features, k, table, mean_nodes_visited )
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:)
        integer,        intent(in) :: k
        type(knn_table), intent(inout) :: table
        real, optional, intent(out) :: mean_nodes_visited
        integer, allocatable :: heap_ids(:)
        real(dp), allocatable :: heap_d2(:)
        integer :: i,heap_size,visited
        integer(kind=8) :: total_visited
        if( any(shape(features)/=[self%ndim,self%npoints]) ) THROW_HARD('kd_tree query feature shape mismatch')
        if( k < 1 .or. k >= self%npoints ) THROW_HARD('kd_tree all-point query requires 1 <= k < npoints')
        call table%kill()
        table%n = self%npoints
        table%k = k
        allocate(table%neighbor(k,self%npoints),source=0)
        allocate(table%distance2(k,self%npoints),source=0.)
        total_visited = 0_8
        !$omp parallel default(shared) private(i,heap_ids,heap_d2,heap_size,visited) reduction(+:total_visited)
        allocate(heap_ids(k),heap_d2(k))
        !$omp do schedule(static)
        do i=1,self%npoints
            heap_size = 0
            visited   = 0
            call query_node(self,features,features(:,i),i,self%root,k,heap_ids,heap_d2,heap_size,visited)
            if( heap_size /= k ) THROW_HARD('kd_tree all-point query returned too few neighbors')
            call sort_neighbor_pairs(heap_ids,heap_d2,k)
            table%neighbor(:,i)  = heap_ids
            table%distance2(:,i) = real(heap_d2)
            total_visited = total_visited+int(visited,kind=8)
        end do
        !$omp end do
        deallocate(heap_ids,heap_d2)
        !$omp end parallel
        if( present(mean_nodes_visited) )then
            mean_nodes_visited = real(total_visited,sp)/real(max(1,self%npoints),sp)
        endif
    end subroutine query_all_kd_tree

    recursive subroutine query_node( self, features, query, exclude, inode, k, heap_ids, heap_d2, heap_size, visited )
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:),query(:)
        integer,        intent(in) :: exclude,inode,k
        integer,        intent(inout) :: heap_ids(:),heap_size,visited
        real(dp),       intent(inout) :: heap_d2(:)
        integer :: pos,point_id,near_node,far_node,q
        real(dp) :: d2,delta,lower_bound
        if( inode == 0 ) return
        visited = visited+1
        if( self%left(inode) == 0 .and. self%right(inode) == 0 )then
            do pos=self%first(inode),self%first(inode)+self%count(inode)-1
                point_id = self%permutation(pos)
                if( point_id == exclude ) cycle
                d2 = 0.0_dp
                do q=1,self%ndim
                    delta = real(query(q),dp)-real(features(q,point_id),dp)
                    d2 = d2+delta*delta
                end do
                call heap_consider(heap_ids,heap_d2,heap_size,k,point_id,d2)
            end do
            return
        endif
        delta = real(query(self%split_dim(inode)),dp)-real(self%split_value(inode),dp)
        if( delta <= 0.0_dp )then
            near_node = self%left(inode)
            far_node  = self%right(inode)
        else
            near_node = self%right(inode)
            far_node  = self%left(inode)
        endif
        call query_node(self,features,query,exclude,near_node,k,heap_ids,heap_d2,heap_size,visited)
        lower_bound = delta*delta
        if( heap_size < k .or. lower_bound <= heap_d2(1) )then
            call query_node(self,features,query,exclude,far_node,k,heap_ids,heap_d2,heap_size,visited)
        endif
    end subroutine query_node

    subroutine heap_consider( ids, d2s, heap_size, capacity, candidate_id, candidate_d2 )
        integer,  intent(inout) :: ids(:),heap_size
        real(dp), intent(inout) :: d2s(:)
        integer,  intent(in)    :: capacity,candidate_id
        real(dp), intent(in)    :: candidate_d2
        integer :: pos,parent
        if( heap_size < capacity )then
            heap_size = heap_size+1
            ids(heap_size) = candidate_id
            d2s(heap_size) = candidate_d2
            pos = heap_size
            do while( pos > 1 )
                parent = pos/2
                if( .not.pair_worse(d2s(pos),ids(pos),d2s(parent),ids(parent)) ) exit
                call swap_pair(ids,d2s,pos,parent)
                pos = parent
            end do
        else if( pair_better(candidate_d2,candidate_id,d2s(1),ids(1)) )then
            ids(1) = candidate_id
            d2s(1) = candidate_d2
            call heap_sift_down(ids,d2s,1,capacity)
        endif
    end subroutine heap_consider

    subroutine heap_sift_down( ids, d2s, root_in, heap_size )
        integer,  intent(inout) :: ids(:)
        real(dp), intent(inout) :: d2s(:)
        integer,  intent(in)    :: root_in,heap_size
        integer :: root,child,worst
        root = root_in
        do
            child = 2*root
            if( child > heap_size ) exit
            worst = child
            if( child+1 <= heap_size )then
                if( pair_worse(d2s(child+1),ids(child+1),d2s(child),ids(child)) ) worst = child+1
            endif
            if( .not.pair_worse(d2s(worst),ids(worst),d2s(root),ids(root)) ) exit
            call swap_pair(ids,d2s,root,worst)
            root = worst
        end do
    end subroutine heap_sift_down

    subroutine sort_neighbor_pairs( ids, d2s, n )
        integer,  intent(inout) :: ids(:)
        real(dp), intent(inout) :: d2s(:)
        integer,  intent(in)    :: n
        integer :: i,j,id
        real(dp) :: d2
        do i=2,n
            id = ids(i)
            d2 = d2s(i)
            j = i-1
            do while( j >= 1 )
                if( .not.pair_better(d2,id,d2s(j),ids(j)) ) exit
                ids(j+1) = ids(j)
                d2s(j+1) = d2s(j)
                j = j-1
            end do
            ids(j+1) = id
            d2s(j+1) = d2
        end do
    end subroutine sort_neighbor_pairs

    pure logical function pair_better( d2_a, id_a, d2_b, id_b ) result(better)
        real(dp), intent(in) :: d2_a,d2_b
        integer,  intent(in) :: id_a,id_b
        better = d2_a < d2_b .or. (d2_a == d2_b .and. id_a < id_b)
    end function pair_better

    pure logical function pair_worse( d2_a, id_a, d2_b, id_b ) result(worse)
        real(dp), intent(in) :: d2_a,d2_b
        integer,  intent(in) :: id_a,id_b
        worse = d2_a > d2_b .or. (d2_a == d2_b .and. id_a > id_b)
    end function pair_worse

    subroutine swap_pair( ids, d2s, i, j )
        integer,  intent(inout) :: ids(:)
        real(dp), intent(inout) :: d2s(:)
        integer,  intent(in)    :: i,j
        integer :: id
        real(dp) :: d2
        id = ids(i); ids(i) = ids(j); ids(j) = id
        d2 = d2s(i); d2s(i) = d2s(j); d2s(j) = d2
    end subroutine swap_pair

    pure subroutine validate_query( self, features, query, k )
        class(kd_tree), intent(in) :: self
        real,           intent(in) :: features(:,:),query(:)
        integer,        intent(in) :: k
        if( self%root == 0 ) error stop 'kd_tree query on an unbuilt tree'
        if( any(shape(features)/=[self%ndim,self%npoints]) ) error stop 'kd_tree query feature shape mismatch'
        if( size(query) /= self%ndim ) error stop 'kd_tree query dimension mismatch'
        if( k < 1 ) error stop 'kd_tree query requires k >= 1'
    end subroutine validate_query

    integer function new_node( self ) result(inode)
        class(kd_tree), intent(inout) :: self
        if( self%nnodes == size(self%left) ) call grow_nodes(self)
        self%nnodes = self%nnodes+1
        inode = self%nnodes
        self%left(inode)        = 0
        self%right(inode)       = 0
        self%first(inode)       = 0
        self%count(inode)       = 0
        self%split_dim(inode)   = 0
        self%split_value(inode) = 0.
    end function new_node

    subroutine reserve_nodes( self, capacity )
        class(kd_tree), intent(inout) :: self
        integer,        intent(in)    :: capacity
        allocate(self%left(capacity),self%right(capacity),self%first(capacity),self%count(capacity), &
            &self%split_dim(capacity),source=0)
        allocate(self%split_value(capacity),source=0.)
    end subroutine reserve_nodes

    subroutine grow_nodes( self )
        class(kd_tree), intent(inout) :: self
        integer, allocatable :: left(:),right(:),first(:),count(:),split_dim(:)
        real,    allocatable :: split_value(:)
        integer :: old_capacity,new_capacity
        old_capacity = size(self%left)
        new_capacity = max(old_capacity+1,2*old_capacity)
        allocate(left(new_capacity),right(new_capacity),first(new_capacity),count(new_capacity), &
            &split_dim(new_capacity),source=0)
        allocate(split_value(new_capacity),source=0.)
        left(:old_capacity)        = self%left
        right(:old_capacity)       = self%right
        first(:old_capacity)       = self%first
        count(:old_capacity)       = self%count
        split_dim(:old_capacity)   = self%split_dim
        split_value(:old_capacity) = self%split_value
        call move_alloc(left,self%left)
        call move_alloc(right,self%right)
        call move_alloc(first,self%first)
        call move_alloc(count,self%count)
        call move_alloc(split_dim,self%split_dim)
        call move_alloc(split_value,self%split_value)
    end subroutine grow_nodes

    pure integer function get_kd_tree_nnodes( self ) result(n)
        class(kd_tree), intent(in) :: self
        n = self%nnodes
    end function get_kd_tree_nnodes

    pure integer function get_kd_tree_height( self ) result(height)
        class(kd_tree), intent(in) :: self
        if( self%root == 0 )then
            height = 0
        else
            height = node_height(self,self%root)
        endif
    end function get_kd_tree_height

    recursive pure integer function node_height( self, inode ) result(height)
        class(kd_tree), intent(in) :: self
        integer,        intent(in) :: inode
        if( inode == 0 )then
            height = 0
        else
            height = 1+max(node_height(self,self%left(inode)),node_height(self,self%right(inode)))
        endif
    end function node_height

    subroutine kill_kd_tree( self )
        class(kd_tree), intent(inout) :: self
        if( allocated(self%permutation) ) deallocate(self%permutation)
        if( allocated(self%left) ) deallocate(self%left)
        if( allocated(self%right) ) deallocate(self%right)
        if( allocated(self%first) ) deallocate(self%first)
        if( allocated(self%count) ) deallocate(self%count)
        if( allocated(self%split_dim) ) deallocate(self%split_dim)
        if( allocated(self%split_value) ) deallocate(self%split_value)
        self%ndim      = 0
        self%npoints   = 0
        self%nnodes    = 0
        self%root      = 0
        self%leaf_size = 32
    end subroutine kill_kd_tree

    pure subroutine swap_int( a, b )
        integer, intent(inout) :: a,b
        integer :: tmp
        tmp = a
        a = b
        b = tmp
    end subroutine swap_int

end module simple_kd_tree
