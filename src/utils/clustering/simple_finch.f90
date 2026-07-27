!@descr: FINCH clustering hierarchy using deterministic exact first-neighbor search
module simple_finch
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_kd_tree, only: kd_tree,knn_table
implicit none
private
#include "simple_local_flags.inc"
integer, parameter :: FINCH_MAX_WARD_SEEDS = 1024

public :: finch_hierarchy
public :: fit_finch
public :: finch_partition_from_first_neighbors
public :: finch_representatives
public :: select_finch_level
public :: refine_finch_level

type :: ward_edge_heap
    private
    integer :: n = 0
    real(dp), allocatable :: cost(:)
    integer, allocatable :: node_a(:),node_b(:)
contains
    procedure :: push => ward_heap_push
    procedure :: pop  => ward_heap_pop
    procedure :: kill => kill_ward_heap
end type ward_edge_heap

type :: finch_hierarchy
    private
    integer :: npoints = 0
    integer :: nlevels = 0
    integer :: maps_used = 0
    integer, allocatable :: nclusters(:)
    integer(kind=8), allocatable :: map_offsets(:)
    integer, allocatable :: parent_maps(:)
contains
    procedure :: get_npoints          => finch_get_npoints
    procedure :: get_nlevels          => finch_get_nlevels
    procedure :: get_nclusters        => finch_get_nclusters
    procedure :: get_labels           => finch_get_labels
    procedure :: get_finest_labels    => finch_get_finest_labels
    procedure :: get_coarsest_labels  => finch_get_coarsest_labels
    procedure :: get_stored_map_count => finch_get_stored_map_count
    procedure :: kill                 => kill_finch_hierarchy
end type finch_hierarchy

contains

    subroutine fit_finch( features, hierarchy, leaf_size )
        real,                  intent(in)    :: features(:,:)
        type(finch_hierarchy), intent(inout) :: hierarchy
        integer, optional,     intent(in)    :: leaf_size
        type(kd_tree) :: tree
        type(knn_table) :: neighbors
        real, allocatable :: current_features(:,:),next_features(:,:)
        integer, allocatable :: first_neighbor(:),level_labels(:),current_weights(:),next_weights(:)
        integer :: n,leaf,current_n,next_n,max_levels
        n = size(features,2)
        if( size(features,1) < 1 .or. n < 1 ) THROW_HARD('FINCH requires a nonempty feature table')
        if( .not.all(ieee_is_finite(features)) ) THROW_HARD('FINCH feature table contains nonfinite values')
        call hierarchy%kill()
        hierarchy%npoints = n
        max_levels = 2
        if( n > 1 ) max_levels = 2+ceiling(log(real(n,dp))/log(2.0_dp))
        allocate(hierarchy%nclusters(max_levels),source=0)
        allocate(hierarchy%map_offsets(max_levels+1),source=0_8)
        allocate(hierarchy%parent_maps(max(1,2*n)),source=0)
        hierarchy%map_offsets(1) = 1_8
        if( n == 1 )then
            call append_level(hierarchy,[1],1)
            call trim_hierarchy(hierarchy)
            return
        endif
        leaf = 32
        if( present(leaf_size) ) leaf = max(1,leaf_size)
        call tree%build(features,leaf)
        call tree%query_all(features,1,neighbors)
        allocate(first_neighbor(n),source=neighbors%neighbor(1,:))
        call finch_partition_from_first_neighbors(first_neighbor,level_labels,next_n)
        call append_level(hierarchy,level_labels,next_n)
        call tree%kill()
        call neighbors%kill()
        deallocate(first_neighbor)
        if( next_n <= 1 )then
            deallocate(level_labels)
            call trim_hierarchy(hierarchy)
            return
        endif
        call weighted_cluster_centroids(features,level_labels,next_n,next_features,next_weights)
        call move_alloc(next_features,current_features)
        call move_alloc(next_weights,current_weights)
        current_n = next_n
        deallocate(level_labels)
        do while( current_n > 1 )
            call tree%build(current_features,leaf)
            call tree%query_all(current_features,1,neighbors)
            allocate(first_neighbor(current_n),source=neighbors%neighbor(1,:))
            call finch_partition_from_first_neighbors(first_neighbor,level_labels,next_n)
            call tree%kill()
            call neighbors%kill()
            deallocate(first_neighbor)
            if( next_n == 1 )then
                deallocate(level_labels)
                exit
            endif
            if( next_n >= current_n ) THROW_HARD('FINCH hierarchy did not reduce the cluster count')
            call append_level(hierarchy,level_labels,next_n)
            call weighted_cluster_centroids(current_features,level_labels,next_n,next_features,next_weights,current_weights)
            call move_alloc(next_features,current_features)
            call move_alloc(next_weights,current_weights)
            current_n = next_n
            deallocate(level_labels)
        end do
        if( allocated(current_features) ) deallocate(current_features)
        if( allocated(current_weights)  ) deallocate(current_weights)
        call trim_hierarchy(hierarchy)
    end subroutine fit_finch

    subroutine finch_partition_from_first_neighbors( first_neighbor, labels, nclusters )
        integer,              intent(in)  :: first_neighbor(:)
        integer, allocatable, intent(out) :: labels(:)
        integer,              intent(out) :: nclusters
        integer, allocatable :: parent(:),tree_size(:),root_label(:)
        integer :: n,i,root
        n = size(first_neighbor)
        if( n < 1 ) THROW_HARD('FINCH first-neighbor table is empty')
        allocate(labels(n),source=0)
        if( n == 1 )then
            if( first_neighbor(1) /= 1 .and. first_neighbor(1) /= 0 ) &
                &THROW_HARD('FINCH singleton first-neighbor index is invalid')
            labels(1) = 1
            nclusters = 1
            return
        endif
        if( any(first_neighbor < 1) .or. any(first_neighbor > n) ) THROW_HARD('FINCH first-neighbor index outside table')
        do i=1,n
            if( first_neighbor(i) == i ) THROW_HARD('FINCH first-neighbor table contains a self edge')
        end do
        allocate(parent(n),tree_size(n),root_label(n),source=0)
        parent = [(i,i=1,n)]
        tree_size = 1
        do i=1,n
            call union_components(parent,tree_size,i,first_neighbor(i))
        end do
        nclusters = 0
        do i=1,n
            root = find_root(parent,i)
            if( root_label(root) == 0 )then
                nclusters = nclusters+1
                root_label(root) = nclusters
            endif
            labels(i) = root_label(root)
        end do
        deallocate(parent,tree_size,root_label)
    end subroutine finch_partition_from_first_neighbors

    subroutine finch_representatives( features, labels, representatives )
        real,                 intent(in)  :: features(:,:)
        integer,              intent(in)  :: labels(:)
        integer, allocatable, intent(out) :: representatives(:)
        real, allocatable :: centroids(:,:)
        integer, allocatable :: populations(:)
        real(dp), allocatable :: best_d2(:)
        real(dp) :: d2,delta
        integer :: n,nclusters,i,q,c
        n = size(features,2)
        if( size(labels) /= n .or. n < 1 .or. size(features,1) < 1 ) THROW_HARD('invalid FINCH representative inputs')
        if( minval(labels) /= 1 ) THROW_HARD('FINCH representative labels must start at one')
        nclusters = maxval(labels)
        if( any(labels > nclusters) ) THROW_HARD('invalid FINCH representative labels')
        call weighted_cluster_centroids(features,labels,nclusters,centroids,populations)
        if( any(populations < 1) ) THROW_HARD('FINCH representative labels are not contiguous')
        allocate(representatives(nclusters),source=0)
        allocate(best_d2(nclusters),source=huge(1.0_dp))
        do i=1,n
            c = labels(i)
            d2 = 0.0_dp
            do q=1,size(features,1)
                delta = real(features(q,i),dp)-real(centroids(q,c),dp)
                d2 = d2+delta*delta
            end do
            if( d2 < best_d2(c) .or. (d2 == best_d2(c) .and. &
                &(representatives(c) == 0 .or. i < representatives(c))) )then
                best_d2(c) = d2
                representatives(c) = i
            endif
        end do
        deallocate(centroids,populations,best_d2)
    end subroutine finch_representatives

    subroutine select_finch_level( features, hierarchy, selected_level, selection_scores )
        real,                  intent(in)  :: features(:,:)
        type(finch_hierarchy), intent(in)  :: hierarchy
        integer,               intent(out) :: selected_level
        real(dp), allocatable, optional, intent(out) :: selection_scores(:)
        real(dp), allocatable :: scores(:),log_variances(:)
        real, allocatable :: current_centroids(:,:),next_centroids(:,:)
        integer, allocatable :: labels(:),current_counts(:),next_counts(:)
        real(dp) :: sse,merge_sse,total_norm2,variance_floor,nobs,cluster_ratio,delta,max_score
        integer :: n,d,level,current_k,next_k,i,q,start,parent
        n = size(features,2)
        d = size(features,1)
        if( n /= hierarchy%npoints .or. n < 1 .or. d < 1 ) &
            &THROW_HARD('invalid FINCH automatic-selection inputs')
        if( hierarchy%nlevels < 1 ) THROW_HARD('cannot select a level from an empty FINCH hierarchy')
        if( .not.all(ieee_is_finite(features)) ) THROW_HARD('FINCH selector feature table contains nonfinite values')
        allocate(scores(hierarchy%nlevels),log_variances(hierarchy%nlevels))
        nobs = real(n,dp)*real(d,dp)
        start = int(hierarchy%map_offsets(1))
        allocate(labels(n),source=hierarchy%parent_maps(start:start+n-1))
        current_k = hierarchy%nclusters(1)
        call weighted_cluster_centroids(features,labels,current_k,current_centroids,current_counts)
        sse = 0.0_dp
        total_norm2 = 0.0_dp
        do i=1,n
            do q=1,d
                delta = real(features(q,i),dp)-real(current_centroids(q,labels(i)),dp)
                sse = sse+delta*delta
                total_norm2 = total_norm2+real(features(q,i),dp)**2
            end do
        end do
        variance_floor = epsilon(1.0_dp)*max(1.0_dp,total_norm2/max(1.0_dp,nobs))
        log_variances(1) = log(max(sse,nobs*variance_floor)/nobs)
        deallocate(labels)
        do level=2,hierarchy%nlevels
            next_k = hierarchy%nclusters(level)
            start = int(hierarchy%map_offsets(level))
            call weighted_cluster_centroids(current_centroids, &
                &hierarchy%parent_maps(start:start+current_k-1),next_k, &
                &next_centroids,next_counts,current_counts)
            merge_sse = 0.0_dp
            do i=1,current_k
                parent = hierarchy%parent_maps(start+i-1)
                do q=1,d
                    delta = real(current_centroids(q,i),dp)-real(next_centroids(q,parent),dp)
                    merge_sse = merge_sse+real(current_counts(i),dp)*delta*delta
                end do
            end do
            sse = sse+merge_sse
            log_variances(level) = log(max(sse,nobs*variance_floor)/nobs)
            call move_alloc(next_centroids,current_centroids)
            call move_alloc(next_counts,current_counts)
            current_k = next_k
        end do
        selected_level = 1
        if( hierarchy%nlevels == 1 )then
            scores(1) = 0.0_dp
        else
            ! Select the finer side of the dominant hierarchy distortion knee.
            ! Normalizing by log(K_l/K_l+1) makes uneven FINCH level sizes
            ! comparable.  A broad knee can span adjacent transitions, so use
            ! its coarsest level whose score is within 90 percent of the peak.
            do level=1,hierarchy%nlevels-1
                cluster_ratio = real(hierarchy%nclusters(level),dp) &
                    &/real(hierarchy%nclusters(level+1),dp)
                scores(level) = (log_variances(level+1)-log_variances(level)) &
                    &/log(cluster_ratio)
            end do
            scores(hierarchy%nlevels) = 0.0_dp
            max_score = maxval(scores(:hierarchy%nlevels-1))
            do level=1,hierarchy%nlevels-1
                if( scores(level) >= 0.9_dp*max_score ) selected_level = level
            end do
        endif
        if( present(selection_scores) )then
            call move_alloc(scores,selection_scores)
        else
            deallocate(scores)
        endif
        deallocate(current_centroids,current_counts,log_variances)
    end subroutine select_finch_level

    subroutine refine_finch_level( features, hierarchy, initial_level, labels, nclusters, &
        &merge_costs, gap_scores, ward_penalty, validity_scores )
        real,                  intent(in)  :: features(:,:)
        type(finch_hierarchy), intent(in)  :: hierarchy
        integer,               intent(in)  :: initial_level
        integer, allocatable,  intent(out) :: labels(:)
        integer,               intent(out) :: nclusters
        real(dp), allocatable, optional, intent(out) :: merge_costs(:),gap_scores(:)
        real(dp), optional,    intent(out) :: ward_penalty
        real(dp), allocatable, optional, intent(out) :: validity_scores(:)
        type(ward_edge_heap) :: heap
        real, allocatable :: centroids(:,:)
        real(dp), allocatable :: costs(:),normalized_costs(:),gaps(:),validity(:),node_sse(:)
        integer, allocatable :: fine_labels(:),node_count(:),node_rep(:)
        integer, allocatable :: merge_rep_a(:),merge_rep_b(:),component_parent(:),root_label(:)
        logical, allocatable :: active(:)
        real(dp) :: cost,cost_floor,delta,penalty_value,best_validity,validity_tolerance
        integer :: n,d,nlevels,fine_level,khi,nmerges,max_nodes
        integer :: i,q,a,b,new_node,merge_index,selected_merges,root
        n = size(features,2)
        d = size(features,1)
        nlevels = hierarchy%nlevels
        if( n /= hierarchy%npoints .or. n < 1 .or. d < 1 ) THROW_HARD('invalid FINCH refinement inputs')
        if( initial_level < 1 .or. initial_level > nlevels ) THROW_HARD('FINCH refinement level outside range')
        if( .not.all(ieee_is_finite(features)) ) THROW_HARD('FINCH refinement feature table contains nonfinite values')
        if( nlevels == 1 )then
            call hierarchy%get_labels(1,labels)
            nclusters = hierarchy%nclusters(1)
            if( present(merge_costs) ) allocate(merge_costs(0))
            if( present(gap_scores)  ) allocate(gap_scores(0))
            if( present(ward_penalty) ) ward_penalty = 0.0_dp
            if( present(validity_scores) ) allocate(validity_scores(1),source=0.0_dp)
            return
        endif
        fine_level = max(1,initial_level-2)
        do while( hierarchy%nclusters(fine_level) > FINCH_MAX_WARD_SEEDS .and. fine_level < nlevels )
            fine_level = fine_level+1
        end do
        call hierarchy%get_labels(fine_level,fine_labels)
        khi = hierarchy%nclusters(fine_level)
        nmerges = khi-1
        if( nmerges < 1 ) THROW_HARD('FINCH refinement requires at least two seed clusters')
        max_nodes = khi+nmerges
        allocate(centroids(d,max_nodes),source=0.)
        allocate(node_sse(max_nodes),source=0.0_dp)
        allocate(node_count(max_nodes),node_rep(max_nodes),source=0)
        allocate(active(max_nodes),source=.false.)
        call initialize_refinement_nodes(features,fine_labels,khi,centroids(:,:khi),node_count(:khi))
        do i=1,n
            do q=1,d
                delta = real(features(q,i),dp)-real(centroids(q,fine_labels(i)),dp)
                node_sse(fine_labels(i)) = node_sse(fine_labels(i))+delta*delta
            end do
        end do
        do i=1,khi
            node_rep(i) = i
            active(i) = .true.
        end do
        do a=1,khi-1
            do b=a+1,khi
                call heap%push(ward_merge_cost(centroids,node_count,a,b),a,b)
            end do
        end do
        allocate(costs(nmerges),normalized_costs(nmerges),validity(nmerges+1))
        allocate(merge_rep_a(nmerges),merge_rep_b(nmerges))
        validity(1) = ward_validity_score(centroids,node_count,node_sse,active,khi)
        new_node = khi
        do merge_index=1,nmerges
            do
                call heap%pop(cost,a,b)
                if( active(a) .and. active(b) ) exit
            end do
            new_node = new_node+1
            costs(merge_index) = cost
            merge_rep_a(merge_index) = node_rep(a)
            merge_rep_b(merge_index) = node_rep(b)
            node_count(new_node) = node_count(a)+node_count(b)
            node_sse(new_node) = node_sse(a)+node_sse(b)+cost
            normalized_costs(merge_index) = cost/real(node_count(new_node),dp)
            node_rep(new_node) = min(node_rep(a),node_rep(b))
            do q=1,d
                centroids(q,new_node) = real( &
                    &(real(node_count(a),dp)*real(centroids(q,a),dp) &
                    &+real(node_count(b),dp)*real(centroids(q,b),dp)) &
                    &/real(node_count(new_node),dp),sp)
            end do
            active(a) = .false.
            active(b) = .false.
            active(new_node) = .true.
            do i=1,new_node-1
                if( .not.active(i) ) cycle
                call heap%push(ward_merge_cost(centroids,node_count,i,new_node),i,new_node)
            end do
            validity(merge_index+1) = ward_validity_score(centroids,node_count,node_sse,active,new_node)
        end do
        call heap%kill()
        allocate(gaps(max(0,nmerges-1)),source=0.0_dp)
        selected_merges = max(0,khi-hierarchy%nclusters(initial_level))
        if( nmerges > 1 )then
            cost_floor = epsilon(1.0_dp)*max(1.0_dp,maxval(normalized_costs))
            do i=1,nmerges-1
                gaps(i) = log((max(0.0_dp,normalized_costs(i+1))+cost_floor) &
                    &/(max(0.0_dp,normalized_costs(i))+cost_floor))
            end do
        endif
        best_validity = maxval(validity(:nmerges))
        validity_tolerance = 1.0e-6_dp*max(1.0_dp,abs(best_validity))
        do i=0,nmerges-1
            if( validity(i+1) >= best_validity-validity_tolerance ) selected_merges = i
        end do
        allocate(component_parent(khi),root_label(khi),source=0)
        component_parent = [(i,i=1,khi)]
        root_label = 1
        do i=1,selected_merges
            call union_components(component_parent,root_label,merge_rep_a(i),merge_rep_b(i))
        end do
        root_label = 0
        allocate(labels(n),source=0)
        nclusters = 0
        do i=1,n
            root = find_root(component_parent,fine_labels(i))
            if( root_label(root) == 0 )then
                nclusters = nclusters+1
                root_label(root) = nclusters
            endif
            labels(i) = root_label(root)
        end do
        if( nclusters /= khi-selected_merges ) THROW_HARD('FINCH refinement cluster count mismatch')
        penalty_value = 0.0_dp
        if( selected_merges >= 1 .and. selected_merges < nmerges )then
            penalty_value = sqrt(max(0.0_dp,costs(selected_merges)) &
                &*max(0.0_dp,costs(selected_merges+1)))
        endif
        if( present(ward_penalty) ) ward_penalty = penalty_value
        if( present(merge_costs) )then
            call move_alloc(costs,merge_costs)
        else
            deallocate(costs)
        endif
        if( present(gap_scores) )then
            call move_alloc(gaps,gap_scores)
        else
            deallocate(gaps)
        endif
        if( present(validity_scores) )then
            call move_alloc(validity,validity_scores)
        else
            deallocate(validity)
        endif
        deallocate(fine_labels,centroids,node_count,node_rep,active)
        deallocate(normalized_costs,node_sse,merge_rep_a,merge_rep_b,component_parent,root_label)
    end subroutine refine_finch_level

    subroutine initialize_refinement_nodes( features, labels, nclusters, centroids, populations )
        real,    intent(in)  :: features(:,:)
        integer, intent(in)  :: labels(:),nclusters
        real,    intent(out) :: centroids(:,:)
        integer, intent(out) :: populations(:)
        real, allocatable :: computed_centroids(:,:)
        integer, allocatable :: computed_populations(:)
        call weighted_cluster_centroids(features,labels,nclusters,computed_centroids,computed_populations)
        centroids = computed_centroids
        populations = computed_populations
        deallocate(computed_centroids,computed_populations)
    end subroutine initialize_refinement_nodes

    real(dp) function ward_merge_cost( centroids, populations, a, b ) result(cost)
        real,    intent(in) :: centroids(:,:)
        integer, intent(in) :: populations(:),a,b
        real(dp) :: delta,distance2
        integer :: q
        distance2 = 0.0_dp
        do q=1,size(centroids,1)
            delta = real(centroids(q,a),dp)-real(centroids(q,b),dp)
            distance2 = distance2+delta*delta
        end do
        cost = real(populations(a),dp)*real(populations(b),dp) &
            &/real(populations(a)+populations(b),dp)*distance2
    end function ward_merge_cost

    real(dp) function ward_validity_score( centroids, populations, node_sse, active, nused ) result(score)
        real,    intent(in) :: centroids(:,:)
        integer, intent(in) :: populations(:),nused
        real(dp), intent(in) :: node_sse(:)
        logical, intent(in) :: active(:)
        type(kd_tree) :: tree
        type(knn_table) :: neighbors
        real, allocatable :: active_centroids(:,:)
        integer, allocatable :: active_nodes(:)
        real(dp) :: radius,separation,denom,weighted_score
        integer(kind=8) :: total_population
        integer :: i,j,k
        k = count(active(:nused))
        if( k < 2 )then
            score = -huge(1.0_dp)
            return
        endif
        allocate(active_centroids(size(centroids,1),k),active_nodes(k))
        j = 0
        do i=1,nused
            if( .not.active(i) ) cycle
            j = j+1
            active_centroids(:,j) = centroids(:,i)
            active_nodes(j) = i
        end do
        call tree%build(active_centroids)
        call tree%query_all(active_centroids,1,neighbors)
        weighted_score = 0.0_dp
        total_population = 0_8
        do j=1,k
            i = active_nodes(j)
            radius = sqrt(max(0.0_dp,node_sse(i))/real(populations(i),dp))
            separation = sqrt(max(0.0_dp,real(neighbors%distance2(1,j),dp)))
            denom = max(radius,separation,epsilon(1.0_dp))
            weighted_score = weighted_score+real(populations(i),dp)*(separation-radius)/denom
            total_population = total_population+int(populations(i),kind=8)
        end do
        score = weighted_score/real(max(1_8,total_population),dp)
        call tree%kill()
        call neighbors%kill()
        deallocate(active_centroids,active_nodes)
    end function ward_validity_score

    subroutine ward_heap_push( self, cost, node_a, node_b )
        class(ward_edge_heap), intent(inout) :: self
        real(dp),              intent(in)    :: cost
        integer,               intent(in)    :: node_a,node_b
        integer :: pos,parent,a,b
        a = min(node_a,node_b)
        b = max(node_a,node_b)
        if( .not.allocated(self%cost) ) call reserve_ward_heap(self,16)
        if( self%n == size(self%cost) ) call grow_ward_heap(self)
        self%n = self%n+1
        pos = self%n
        self%cost(pos) = cost
        self%node_a(pos) = a
        self%node_b(pos) = b
        do while( pos > 1 )
            parent = pos/2
            if( .not.ward_edge_better(self%cost(pos),self%node_a(pos),self%node_b(pos), &
                &self%cost(parent),self%node_a(parent),self%node_b(parent)) ) exit
            call swap_ward_edges(self,pos,parent)
            pos = parent
        end do
    end subroutine ward_heap_push

    subroutine ward_heap_pop( self, cost, node_a, node_b )
        class(ward_edge_heap), intent(inout) :: self
        real(dp),              intent(out)   :: cost
        integer,               intent(out)   :: node_a,node_b
        integer :: pos,child,best
        if( self%n < 1 ) THROW_HARD('FINCH Ward refinement heap is empty')
        cost = self%cost(1)
        node_a = self%node_a(1)
        node_b = self%node_b(1)
        self%cost(1) = self%cost(self%n)
        self%node_a(1) = self%node_a(self%n)
        self%node_b(1) = self%node_b(self%n)
        self%n = self%n-1
        pos = 1
        do
            child = 2*pos
            if( child > self%n ) exit
            best = child
            if( child+1 <= self%n )then
                if( ward_edge_better(self%cost(child+1),self%node_a(child+1),self%node_b(child+1), &
                    &self%cost(child),self%node_a(child),self%node_b(child)) ) best = child+1
            endif
            if( .not.ward_edge_better(self%cost(best),self%node_a(best),self%node_b(best), &
                &self%cost(pos),self%node_a(pos),self%node_b(pos)) ) exit
            call swap_ward_edges(self,pos,best)
            pos = best
        end do
    end subroutine ward_heap_pop

    pure logical function ward_edge_better( cost_a, a1, a2, cost_b, b1, b2 ) result(better)
        real(dp), intent(in) :: cost_a,cost_b
        integer,  intent(in) :: a1,a2,b1,b2
        better = cost_a < cost_b .or. (cost_a == cost_b .and. &
            &(a1 < b1 .or. (a1 == b1 .and. a2 < b2)))
    end function ward_edge_better

    subroutine reserve_ward_heap( self, capacity )
        class(ward_edge_heap), intent(inout) :: self
        integer,               intent(in)    :: capacity
        allocate(self%cost(capacity),source=0.0_dp)
        allocate(self%node_a(capacity),self%node_b(capacity),source=0)
    end subroutine reserve_ward_heap

    subroutine grow_ward_heap( self )
        class(ward_edge_heap), intent(inout) :: self
        real(dp), allocatable :: cost(:)
        integer, allocatable :: node_a(:),node_b(:)
        integer :: old_capacity,new_capacity
        old_capacity = size(self%cost)
        new_capacity = max(old_capacity+1,2*old_capacity)
        allocate(cost(new_capacity),source=0.0_dp)
        allocate(node_a(new_capacity),node_b(new_capacity),source=0)
        cost(:old_capacity) = self%cost
        node_a(:old_capacity) = self%node_a
        node_b(:old_capacity) = self%node_b
        call move_alloc(cost,self%cost)
        call move_alloc(node_a,self%node_a)
        call move_alloc(node_b,self%node_b)
    end subroutine grow_ward_heap

    subroutine swap_ward_edges( self, i, j )
        class(ward_edge_heap), intent(inout) :: self
        integer,               intent(in)    :: i,j
        real(dp) :: cost
        integer :: node
        cost = self%cost(i); self%cost(i) = self%cost(j); self%cost(j) = cost
        node = self%node_a(i); self%node_a(i) = self%node_a(j); self%node_a(j) = node
        node = self%node_b(i); self%node_b(i) = self%node_b(j); self%node_b(j) = node
    end subroutine swap_ward_edges

    subroutine kill_ward_heap( self )
        class(ward_edge_heap), intent(inout) :: self
        if( allocated(self%cost) ) deallocate(self%cost)
        if( allocated(self%node_a) ) deallocate(self%node_a)
        if( allocated(self%node_b) ) deallocate(self%node_b)
        self%n = 0
    end subroutine kill_ward_heap

    subroutine weighted_cluster_centroids( features, labels, nclusters, centroids, populations, weights )
        real,                 intent(in)  :: features(:,:)
        integer,              intent(in)  :: labels(:),nclusters
        real, allocatable,    intent(out) :: centroids(:,:)
        integer, allocatable, intent(out) :: populations(:)
        integer, optional,    intent(in)  :: weights(:)
        integer, allocatable :: member_counts(:),offsets(:),cursor(:),members(:)
        real(dp) :: accum
        integer :: n,i,c,pos,q,w
        n = size(features,2)
        if( size(labels) /= n .or. nclusters < 1 ) THROW_HARD('invalid FINCH centroid dimensions')
        if( any(labels < 1) .or. any(labels > nclusters) ) THROW_HARD('FINCH centroid label outside range')
        if( present(weights) )then
            if( size(weights) /= n .or. any(weights < 1) ) THROW_HARD('invalid FINCH centroid weights')
        endif
        allocate(member_counts(nclusters),offsets(nclusters+1),cursor(nclusters),source=0)
        do i=1,n
            member_counts(labels(i)) = member_counts(labels(i))+1
        end do
        if( any(member_counts < 1) ) THROW_HARD('FINCH centroid labels are not contiguous')
        offsets(1) = 1
        do c=1,nclusters
            offsets(c+1) = offsets(c)+member_counts(c)
        end do
        cursor = offsets(:nclusters)
        allocate(members(n))
        do i=1,n
            c = labels(i)
            members(cursor(c)) = i
            cursor(c) = cursor(c)+1
        end do
        allocate(centroids(size(features,1),nclusters),source=0.)
        allocate(populations(nclusters),source=0)
        !$omp parallel do default(shared) private(c,pos,i,q,w,accum) schedule(static)
        do c=1,nclusters
            do pos=offsets(c),offsets(c+1)-1
                i = members(pos)
                w = 1
                if( present(weights) ) w = weights(i)
                populations(c) = populations(c)+w
            end do
            do q=1,size(features,1)
                accum = 0.0_dp
                do pos=offsets(c),offsets(c+1)-1
                    i = members(pos)
                    w = 1
                    if( present(weights) ) w = weights(i)
                    accum = accum+real(w,dp)*real(features(q,i),dp)
                end do
                centroids(q,c) = real(accum/real(populations(c),dp),sp)
            end do
        end do
        !$omp end parallel do
        deallocate(member_counts,offsets,cursor,members)
    end subroutine weighted_cluster_centroids

    recursive integer function find_root( parent, node ) result(root)
        integer, intent(inout) :: parent(:)
        integer, intent(in)    :: node
        if( parent(node) == node )then
            root = node
        else
            parent(node) = find_root(parent,parent(node))
            root = parent(node)
        endif
    end function find_root

    subroutine union_components( parent, tree_size, a, b )
        integer, intent(inout) :: parent(:),tree_size(:)
        integer, intent(in)    :: a,b
        integer :: root_a,root_b
        root_a = find_root(parent,a)
        root_b = find_root(parent,b)
        if( root_a == root_b ) return
        if( tree_size(root_a) < tree_size(root_b) .or. &
            &(tree_size(root_a) == tree_size(root_b) .and. root_b < root_a) )then
            parent(root_a) = root_b
            tree_size(root_b) = tree_size(root_b)+tree_size(root_a)
        else
            parent(root_b) = root_a
            tree_size(root_a) = tree_size(root_a)+tree_size(root_b)
        endif
    end subroutine union_components

    subroutine append_level( hierarchy, labels, nclusters )
        type(finch_hierarchy), intent(inout) :: hierarchy
        integer,               intent(in)    :: labels(:),nclusters
        integer :: level,start,finish
        level = hierarchy%nlevels+1
        if( level > size(hierarchy%nclusters) ) THROW_HARD('FINCH hierarchy level capacity exceeded')
        start  = hierarchy%maps_used+1
        finish = hierarchy%maps_used+size(labels)
        if( finish > size(hierarchy%parent_maps) ) THROW_HARD('FINCH hierarchy map capacity exceeded')
        hierarchy%parent_maps(start:finish) = labels
        hierarchy%maps_used = finish
        hierarchy%nlevels = level
        hierarchy%nclusters(level) = nclusters
        hierarchy%map_offsets(level)   = int(start,kind=8)
        hierarchy%map_offsets(level+1) = int(finish+1,kind=8)
    end subroutine append_level

    subroutine trim_hierarchy( hierarchy )
        type(finch_hierarchy), intent(inout) :: hierarchy
        integer, allocatable :: nclusters(:),parent_maps(:)
        integer(kind=8), allocatable :: offsets(:)
        allocate(nclusters(hierarchy%nlevels),source=hierarchy%nclusters(:hierarchy%nlevels))
        allocate(offsets(hierarchy%nlevels+1),source=hierarchy%map_offsets(:hierarchy%nlevels+1))
        allocate(parent_maps(hierarchy%maps_used),source=hierarchy%parent_maps(:hierarchy%maps_used))
        call move_alloc(nclusters,hierarchy%nclusters)
        call move_alloc(offsets,hierarchy%map_offsets)
        call move_alloc(parent_maps,hierarchy%parent_maps)
    end subroutine trim_hierarchy

    subroutine finch_get_labels( self, level, labels )
        class(finch_hierarchy), intent(in) :: self
        integer,                intent(in) :: level
        integer, allocatable,   intent(out) :: labels(:)
        integer :: i,l,start
        if( level < 1 .or. level > self%nlevels ) THROW_HARD('FINCH hierarchy level outside range')
        allocate(labels(self%npoints))
        start = int(self%map_offsets(1))
        labels = self%parent_maps(start:start+self%npoints-1)
        do l=2,level
            start = int(self%map_offsets(l))
            do i=1,self%npoints
                labels(i) = self%parent_maps(start+labels(i)-1)
            end do
        end do
    end subroutine finch_get_labels

    subroutine finch_get_finest_labels( self, labels )
        class(finch_hierarchy), intent(in) :: self
        integer, allocatable, intent(out) :: labels(:)
        call self%get_labels(1,labels)
    end subroutine finch_get_finest_labels

    subroutine finch_get_coarsest_labels( self, labels )
        class(finch_hierarchy), intent(in) :: self
        integer, allocatable, intent(out) :: labels(:)
        call self%get_labels(self%nlevels,labels)
    end subroutine finch_get_coarsest_labels

    pure integer function finch_get_npoints( self ) result(n)
        class(finch_hierarchy), intent(in) :: self
        n = self%npoints
    end function finch_get_npoints

    pure integer function finch_get_nlevels( self ) result(n)
        class(finch_hierarchy), intent(in) :: self
        n = self%nlevels
    end function finch_get_nlevels

    pure integer function finch_get_nclusters( self, level ) result(n)
        class(finch_hierarchy), intent(in) :: self
        integer, intent(in) :: level
        n = 0
        if( level >= 1 .and. level <= self%nlevels ) n = self%nclusters(level)
    end function finch_get_nclusters

    pure integer function finch_get_stored_map_count( self ) result(n)
        class(finch_hierarchy), intent(in) :: self
        n = self%maps_used
    end function finch_get_stored_map_count

    subroutine kill_finch_hierarchy( self )
        class(finch_hierarchy), intent(inout) :: self
        if( allocated(self%nclusters)  ) deallocate(self%nclusters)
        if( allocated(self%map_offsets) ) deallocate(self%map_offsets)
        if( allocated(self%parent_maps) ) deallocate(self%parent_maps)
        self%npoints = 0
        self%nlevels = 0
        self%maps_used = 0
    end subroutine kill_finch_hierarchy

end module simple_finch
