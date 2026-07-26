!@descr: validates angularly gated registered-residual kNN construction
program simple_test_flex_diffmap_graph
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, build_gated_euclidean_knn_graph, &
    &find_gated_euclidean_neighbors_rows, build_gated_euclidean_graph_from_neighbors, build_euclidean_knn_graph
use simple_diffusion_maps, only: embed_graph
implicit none
type(diffmap_graph) :: graph,graph_parts,graph_many,graph_density
real :: features(2,6),features_many(2,24),features_density(1,32),dirs(3,2),cmean,angle
integer, allocatable :: nbrs1(:,:),nbrs2(:,:),nbrs(:,:),nc1(:),nc2(:),nc(:)
real, allocatable :: d2s1(:,:),d2s2(:,:),d2s(:,:)
real, allocatable :: coords(:,:),raw_coords(:,:),eigvals(:),eigenfunctions(:,:),nystrom_coords(:,:),stationary(:)
real :: phi0(32),markov_lhs,max_markov_residual,max_raw_residual
integer :: proj(6),cmin,cmax,i,p,k,j

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
if( graph%ncomponents()/=2 ) stop 'angularly gated graph component count mismatch'

! A connected, nonuniformly sampled interval exercises the distinction between
! symmetric eigenvectors and the right eigenfunctions of the Markov operator.
do i=1,24
    features_density(1,i)=0.4*real(i-1)/24.
end do
do i=25,32
    features_density(1,i)=0.4+0.6*real(i-25)/7.
end do
call build_euclidean_knn_graph(features_density,6,graph_density,alpha=1.0)
if( graph_density%ncomponents()/=1 ) stop 'density-normalized regression graph is disconnected'
call embed_graph(graph_density,4,coords,eigvals,raw_coords,eigenfunctions,nystrom_coords,stationary)
if( any(shape(eigenfunctions)/=[4,32]) .or. any(shape(nystrom_coords)/=[4,32]) ) &
    &stop 'diffusion spectral output shape mismatch'
if( maxval(abs(eigenfunctions-nystrom_coords))>1.e-4 ) stop 'training-node Nystrom coefficients differ from eigenfunctions'
if( abs(sum(stationary)-1.)>1.e-5 .or. minval(stationary)<=0. ) stop 'invalid diffusion stationary measure'
phi0=sqrt(stationary)
max_markov_residual=0.
max_raw_residual=0.
do k=1,size(eigenfunctions,1)
    do i=1,graph_density%n
        markov_lhs=0.
        do p=graph_density%rowptr(i),graph_density%rowptr(i+1)-1
            j=graph_density%colind(p)
            markov_lhs=markov_lhs+graph_density%wnorm(p)*phi0(j)*eigenfunctions(k,j)
        end do
        markov_lhs=markov_lhs/phi0(i)
        max_markov_residual=max(max_markov_residual,abs(markov_lhs-eigvals(k)*eigenfunctions(k,i)))
        max_raw_residual=max(max_raw_residual,abs(raw_coords(k,i)-eigvals(k)*eigenfunctions(k,i)))
    end do
end do
if( max_markov_residual>5.e-4 ) stop 'diffusion coordinates are not right Markov eigenfunctions'
if( max_raw_residual>5.e-5 ) stop 'raw diffusion coordinates do not equal lambda times right eigenfunctions'
deallocate(coords,raw_coords,eigvals,eigenfunctions,nystrom_coords,stationary)
call graph_density%kill()
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
