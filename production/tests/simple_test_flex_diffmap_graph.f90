!@descr: validates angularly gated registered-residual kNN construction
program simple_test_flex_diffmap_graph
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, build_gated_euclidean_knn_graph, &
    &find_gated_euclidean_neighbors_rows, build_gated_euclidean_graph_from_neighbors, build_euclidean_knn_graph, &
    &projection_occupancy_weights
use simple_diffusion_maps, only: embed_graph
implicit none
type(diffmap_graph) :: graph,graph_parts,graph_many,graph_balanced
real :: features(2,6), features_many(2,24), dirs(3,2), cmean,angle
integer, allocatable :: nbrs1(:,:),nbrs2(:,:),nbrs(:,:),nc1(:),nc2(:),nc(:)
real, allocatable :: d2s1(:,:),d2s2(:,:),d2s(:,:),view_weights(:)
real, allocatable :: coords(:,:),raw_coords(:,:),eigvals(:),eigenfunctions(:,:),nystrom_coords(:,:)
real :: weighted_degree(6),perron(6),lhs
integer :: proj(6),uneven_proj(6),cmin,cmax,i,p,j,noccupied,occ_min,occ_max

features(:,1) = [0.0,0.0]
features(:,2) = [0.1,0.0]
features(:,3) = [0.0,0.1]
features(:,4) = [9.9,10.0]
features(:,5) = [10.0,9.9]
features(:,6) = [10.0,10.0]
proj = [1,1,1,2,2,2]
dirs(:,1) = [0.,0.,1.]
dirs(:,2) = [1.,0.,0.]

call build_gated_euclidean_knn_graph(features,proj,dirs,2,2,graph,cmin,cmax,cmean)
if( graph%n /= 6 ) stop 'gated graph particle count mismatch'
if( graph%k_nn /= 2 ) stop 'gated graph k_nn mismatch'
if( cmin /= 2 .or. cmax /= 2 .or. abs(cmean-2.) > 1.e-6 ) stop 'gated graph candidate cap mismatch'
do i=1,6
    do p=graph%rowptr(i),graph%rowptr(i+1)-1
        if( proj(graph%colind(p)) /= proj(i) ) stop 'angular gate admitted a distant projection bin'
    end do
end do
call find_gated_euclidean_neighbors_rows(features,proj,dirs,2,2,[1,2,3],nbrs1,d2s1,nc1)
call find_gated_euclidean_neighbors_rows(features,proj,dirs,2,2,[4,5,6],nbrs2,d2s2,nc2)
allocate(nbrs(2,6),d2s(2,6),nc(6))
nbrs(:,:3)=nbrs1; nbrs(:,4:)=nbrs2
d2s(:,:3)=d2s1; d2s(:,4:)=d2s2
nc(:3)=nc1; nc(4:)=nc2
call build_gated_euclidean_graph_from_neighbors(6,nbrs,d2s,nc,graph_parts)
if( any(graph_parts%rowptr/=graph%rowptr) ) stop 'distributed graph row pointers differ'
if( any(graph_parts%colind/=graph%colind) ) stop 'distributed graph neighbors differ'
if( maxval(abs(graph_parts%w-graph%w))>1.e-6 ) stop 'distributed graph weights differ'

uneven_proj=[1,1,1,2,2,3]
call projection_occupancy_weights(uneven_proj,4,view_weights,noccupied,occ_min,occ_max)
if( noccupied/=3 .or. occ_min/=1 .or. occ_max/=3 ) stop 'view occupancy diagnostics mismatch'
if( abs(sum(view_weights)-6.)>1.e-6 ) stop 'view weights do not have unit mean'
if( abs(sum(view_weights(1:3))-2.)>1.e-6 .or. abs(sum(view_weights(4:5))-2.)>1.e-6 .or. &
    &abs(view_weights(6)-2.)>1.e-6 ) stop 'occupied projection bins do not have equal mass'
call build_gated_euclidean_knn_graph(features,proj,dirs,2,2,graph_balanced,cmin,cmax,cmean, &
    &sample_weights=view_weights)
weighted_degree=0.
do i=1,graph_balanced%n
    do p=graph_balanced%rowptr(i),graph_balanced%rowptr(i+1)-1
        j=graph_balanced%colind(p)
        weighted_degree(i)=weighted_degree(i)+graph_balanced%w(p)*view_weights(j)
    end do
end do
perron=sqrt(view_weights*weighted_degree)
do i=1,graph_balanced%n
    lhs=0.
    do p=graph_balanced%rowptr(i),graph_balanced%rowptr(i+1)-1
        j=graph_balanced%colind(p)
        lhs=lhs+graph_balanced%wnorm(p)*perron(j)
    end do
    if( abs(lhs-perron(i))>2.e-5 ) stop 'view-balanced symmetric operator has invalid Perron vector'
end do
deallocate(view_weights)
call graph_balanced%kill()

call embed_graph(graph,2,coords,eigvals,raw_coords,eigenfunctions,nystrom_coords)
if( any(shape(eigenfunctions)/=[2,6]) .or. any(shape(nystrom_coords)/=[2,6]) ) &
    &stop 'diffusion spectral output shape mismatch'
if( maxval(abs(eigenfunctions-nystrom_coords))>1.e-4 ) stop 'training-node Nystrom coefficients differ from eigenfunctions'
deallocate(coords,raw_coords,eigvals,eigenfunctions,nystrom_coords)
do i=1,size(features_many,2)
    angle=2.*acos(-1.)*real(i-1)/real(size(features_many,2))
    features_many(:,i)=[cos(angle),sin(angle)]
end do
call build_euclidean_knn_graph(features_many,6,graph_many)
call embed_graph(graph_many,21,coords,eigvals)
if( size(eigvals)/=21 .or. any(shape(coords)/=[21,24]) ) stop 'diffusion scan was capped below requested rank'
deallocate(coords,eigvals)
call graph_many%kill()
deallocate(nbrs1,nbrs2,nbrs,d2s1,d2s2,d2s,nc1,nc2,nc)
call graph_parts%kill()
call graph%kill()
call simple_end('**** SIMPLE_FLEX_DIFFMAP_GRAPH TEST NORMAL STOP ****')

end program simple_test_flex_diffmap_graph
