!@descr: validates angularly gated registered-residual kNN construction
program simple_test_flex_diffmap_graph
use simple_core_module_api
use simple_diff_map_graphs, only: diffmap_graph, build_gated_euclidean_knn_graph
implicit none
type(diffmap_graph) :: graph
real :: features(2,6), dirs(3,2), cmean
integer :: proj(6), cmin, cmax, i, p

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
call graph%kill()
call simple_end('**** SIMPLE_FLEX_DIFFMAP_GRAPH TEST NORMAL STOP ****')

end program simple_test_flex_diffmap_graph
