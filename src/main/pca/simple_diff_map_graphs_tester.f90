!@descr: unit tests for the angularly gated kNN diffusion-map graph engine (simple_diff_map_graphs, simple_diffusion_maps)
! Six particles in two tight feature clusters and two projection bins: the gated graph keeps its
! neighbours inside the projection bin, a graph assembled from row blocks equals the whole one,
! view-occupancy weights give every occupied bin equal mass with unit mean, the view-balanced
! symmetric operator has the Perron vector sqrt(w * weighted degree), the Nystrom coefficients of
! the training nodes equal the eigenfunctions, and the spectral scan is not capped below the
! requested rank (24 points on a circle, rank 21).
module simple_diff_map_graphs_tester
use simple_diff_map_graphs, only: diffmap_graph, build_gated_euclidean_knn_graph, &
    &find_gated_euclidean_neighbors_rows, build_gated_euclidean_graph_from_neighbors, build_euclidean_knn_graph, &
    &projection_occupancy_weights
use simple_diffusion_maps,  only: embed_graph
use simple_test_utils
implicit none
private
public :: run_all_diff_map_graphs_tests

contains

    subroutine run_all_diff_map_graphs_tests()
        write(*,'(A)') '**** running all diffusion-map graph tests ****'
        call test_gated_graph()
        call test_view_balanced_operator()
        call test_spectral_embedding()
    end subroutine run_all_diff_map_graphs_tests

    !> two clusters of three particles, cluster = projection bin
    subroutine two_clusters( features, proj, dirs )
        real,    intent(out) :: features(2,6), dirs(3,2)
        integer, intent(out) :: proj(6)
        features(:,1) = [0.0,0.0]
        features(:,2) = [0.1,0.0]
        features(:,3) = [0.0,0.1]
        features(:,4) = [9.9,10.0]
        features(:,5) = [10.0,9.9]
        features(:,6) = [10.0,10.0]
        proj = [1,1,1,2,2,2]
        dirs(:,1) = [0.,0.,1.]
        dirs(:,2) = [1.,0.,0.]
    end subroutine two_clusters

    !> the angular gate and the block-row assembly used by distributed builds
    subroutine test_gated_graph()
        type(diffmap_graph) :: graph, graph_parts
        real    :: features(2,6), dirs(3,2), cmean
        integer :: proj(6), cmin, cmax, i, p, nforeign
        integer, allocatable :: nbrs1(:,:), nbrs2(:,:), nbrs(:,:), nc1(:), nc2(:), nc(:)
        real,    allocatable :: d2s1(:,:), d2s2(:,:), d2s(:,:)
        write(*,'(A)') 'test_gated_graph'
        call two_clusters(features, proj, dirs)
        call build_gated_euclidean_knn_graph(features,proj,dirs,2,2,graph,cmin,cmax,cmean)
        call assert_int(6, graph%n,    'gated graph: particle count')
        call assert_int(2, graph%k_nn, 'gated graph: k_nn')
        call assert_true(cmin == 2 .and. cmax == 2 .and. abs(cmean-2.) <= 1.e-6, 'gated graph: candidate cap of 2 per particle')
        nforeign = 0
        do i = 1,6
            do p = graph%rowptr(i),graph%rowptr(i+1)-1
                if( proj(graph%colind(p)) /= proj(i) ) nforeign = nforeign + 1
            end do
        end do
        call assert_int(0, nforeign, 'the angular gate admits no neighbour from a distant projection bin')
        call find_gated_euclidean_neighbors_rows(features,proj,dirs,2,2,[1,2,3],nbrs1,d2s1,nc1)
        call find_gated_euclidean_neighbors_rows(features,proj,dirs,2,2,[4,5,6],nbrs2,d2s2,nc2)
        allocate(nbrs(2,6), d2s(2,6), nc(6))
        nbrs(:,:3) = nbrs1; nbrs(:,4:) = nbrs2
        d2s(:,:3)  = d2s1;  d2s(:,4:)  = d2s2
        nc(:3)     = nc1;   nc(4:)     = nc2
        call build_gated_euclidean_graph_from_neighbors(6,nbrs,d2s,nc,graph_parts)
        call assert_true(all(graph_parts%rowptr == graph%rowptr), 'graph from row blocks: row pointers equal the whole build')
        call assert_true(all(graph_parts%colind == graph%colind), 'graph from row blocks: neighbours equal the whole build')
        call assert_true(maxval(abs(graph_parts%w - graph%w)) <= 1.e-6, 'graph from row blocks: weights equal the whole build')
        call graph_parts%kill()
        call graph%kill()
        deallocate(nbrs1, nbrs2, nbrs, d2s1, d2s2, d2s, nc1, nc2, nc)
    end subroutine test_gated_graph

    !> occupancy weights and the Perron vector of the view-balanced symmetric operator
    subroutine test_view_balanced_operator()
        type(diffmap_graph) :: graph_balanced
        real    :: features(2,6), dirs(3,2), cmean, weighted_degree(6), perron(6), lhs, maxdev
        real, allocatable :: view_weights(:)
        integer :: proj(6), uneven_proj(6), cmin, cmax, i, p, j, noccupied, occ_min, occ_max
        write(*,'(A)') 'test_view_balanced_operator'
        call two_clusters(features, proj, dirs)
        uneven_proj = [1,1,1,2,2,3]
        call projection_occupancy_weights(uneven_proj,4,view_weights,noccupied,occ_min,occ_max)
        call assert_true(noccupied == 3 .and. occ_min == 1 .and. occ_max == 3, 'view occupancy: three bins, one to three members')
        call assert_real(6., sum(view_weights), 1.e-6, 'view weights have unit mean')
        call assert_true(abs(sum(view_weights(1:3))-2.) <= 1.e-6 .and. abs(sum(view_weights(4:5))-2.) <= 1.e-6 .and. &
            &abs(view_weights(6)-2.) <= 1.e-6, 'every occupied projection bin has equal mass')
        call build_gated_euclidean_knn_graph(features,proj,dirs,2,2,graph_balanced,cmin,cmax,cmean, &
            &sample_weights=view_weights)
        weighted_degree = 0.
        do i = 1,graph_balanced%n
            do p = graph_balanced%rowptr(i),graph_balanced%rowptr(i+1)-1
                j = graph_balanced%colind(p)
                weighted_degree(i) = weighted_degree(i) + graph_balanced%w(p)*view_weights(j)
            end do
        end do
        perron = sqrt(view_weights*weighted_degree)
        maxdev = 0.
        do i = 1,graph_balanced%n
            lhs = 0.
            do p = graph_balanced%rowptr(i),graph_balanced%rowptr(i+1)-1
                j   = graph_balanced%colind(p)
                lhs = lhs + graph_balanced%wnorm(p)*perron(j)
            end do
            maxdev = max(maxdev, abs(lhs-perron(i)))
        end do
        call assert_true(maxdev <= 2.e-5, 'the view-balanced symmetric operator has sqrt(w * degree) as Perron vector')
        deallocate(view_weights)
        call graph_balanced%kill()
    end subroutine test_view_balanced_operator

    !> Nystrom coefficients of the training nodes and the rank of the spectral scan
    subroutine test_spectral_embedding()
        type(diffmap_graph) :: graph, graph_many
        real    :: features(2,6), dirs(3,2), cmean, features_many(2,24), angle
        real, allocatable :: coords(:,:), raw_coords(:,:), eigvals(:), eigenfunctions(:,:), nystrom_coords(:,:)
        integer :: proj(6), cmin, cmax, i
        write(*,'(A)') 'test_spectral_embedding'
        call two_clusters(features, proj, dirs)
        call build_gated_euclidean_knn_graph(features,proj,dirs,2,2,graph,cmin,cmax,cmean)
        call embed_graph(graph,2,coords,eigvals,raw_coords,eigenfunctions,nystrom_coords)
        call assert_true(all(shape(eigenfunctions) == [2,6]) .and. all(shape(nystrom_coords) == [2,6]), &
            &'diffusion spectral output has rank x particles shape')
        if( all(shape(eigenfunctions) == shape(nystrom_coords)) )then
            call assert_true(maxval(abs(eigenfunctions-nystrom_coords)) <= 1.e-4, &
                &'Nystrom coefficients of the training nodes equal the eigenfunctions')
        endif
        deallocate(coords, raw_coords, eigvals, eigenfunctions, nystrom_coords)
        call graph%kill()
        do i = 1,size(features_many,2)
            angle = 2.*acos(-1.)*real(i-1)/real(size(features_many,2))
            features_many(:,i) = [cos(angle),sin(angle)]
        end do
        call build_euclidean_knn_graph(features_many,6,graph_many)
        call embed_graph(graph_many,21,coords,eigvals)
        call assert_true(size(eigvals) == 21 .and. all(shape(coords) == [21,24]), &
            &'the diffusion scan returns the requested rank 21, not a capped one')
        deallocate(coords, eigvals)
        call graph_many%kill()
    end subroutine test_spectral_embedding

end module simple_diff_map_graphs_tester
