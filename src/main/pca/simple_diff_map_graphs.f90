!@descr: graph construction helpers for diffusion-map
module simple_diff_map_graphs
use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
use simple_core_module_api
use simple_builder,          only: builder
use simple_cluster2D_strategy, only: cluster2D_inmem_strategy
use simple_cmdline,          only: cmdline
use simple_default_clines,   only: set_cluster2D_defaults
use simple_image,            only: image
use simple_imgarr_utils,     only: copy_imgarr, dealloc_imgarr, extract_imgarr, read_stk_into_imgarr, write_imgarr
use simple_ori,              only: ori
use simple_parameters,       only: parameters
use simple_polarft_calc,     only: polarft_calc
use simple_sp_project,       only: sp_project
implicit none
#include "simple_local_flags.inc"

private
public :: build_so3_split_affinity

contains

    subroutine build_so3_split_affinity( params, spproj, pinds, imgs, aff, theta, shift_x, shift_y )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        type(image),          intent(inout) :: imgs(:)
        real, allocatable,    intent(out)   :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        character(len=STDLEN) :: graph_mode
        graph_mode = lowercase(trim(params%so3_graph))
        select case(trim(graph_mode))
            case('projection_registration')
                call build_projection_registered_affinity(params, spproj, pinds, imgs, aff, theta, shift_x, shift_y)
            case('cluster2d')
                call build_cluster2d_affinity(params, spproj, pinds, imgs, aff, theta, shift_x, shift_y)
            case DEFAULT
                THROW_HARD('so3_graph must be cluster2d or projection_registration')
        end select
    end subroutine build_so3_split_affinity

    subroutine build_projection_registered_affinity( params, spproj, pinds, imgs, aff, theta, shift_x, shift_y )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        type(image),          intent(inout) :: imgs(:)
        real, allocatable,    intent(out)   :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        type(parameters), target :: params_pft
        type(polarft_calc) :: pftc
        real, allocatable  :: gauge_angle(:), gauge_img_dist(:), gauge_shift_x(:), gauge_shift_y(:)
        real               :: corr, img_dist, angle, sh(2)
        integer            :: n, i, center_loc
        n = size(pinds)
        allocate(gauge_angle(n), gauge_img_dist(n), gauge_shift_x(n), gauge_shift_y(n), source=0.)
        center_loc = choose_orientation_medoid(spproj, pinds)
        call prepare_so3_pft(params, imgs, params_pft, pftc)
        do i = 1,n
            call align_pair(params, pftc, center_loc, i, corr, img_dist, angle, sh)
            gauge_img_dist(i) = img_dist
            gauge_angle(i)    = angle
            gauge_shift_x(i)  = sh(1)
            gauge_shift_y(i)  = sh(2)
        end do
        call build_gauge_knn_graph(params, spproj, pinds, gauge_angle, gauge_shift_x, gauge_shift_y, gauge_img_dist, &
            aff, theta, shift_x, shift_y)
        call pftc%kill
        if( allocated(gauge_angle)    ) deallocate(gauge_angle)
        if( allocated(gauge_img_dist) ) deallocate(gauge_img_dist)
        if( allocated(gauge_shift_x)  ) deallocate(gauge_shift_x)
        if( allocated(gauge_shift_y)  ) deallocate(gauge_shift_y)
    end subroutine build_projection_registered_affinity

    subroutine build_cluster2d_affinity( params, spproj, pinds, imgs, aff, theta, shift_x, shift_y )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        type(image),          intent(inout) :: imgs(:)
        real, allocatable,    intent(out)   :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        type(image), allocatable :: region_imgs(:), class_avgs(:), rep_imgs(:)
        real, allocatable        :: gauge_angle(:), gauge_img_dist(:), gauge_shift_x(:), gauge_shift_y(:)
        real, allocatable        :: all_gauge_angle(:), all_gauge_img_dist(:), all_gauge_shift_x(:), all_gauge_shift_y(:)
        real, allocatable        :: class_angle(:,:), class_img_dist(:,:), class_shift_x(:,:), class_shift_y(:,:)
        real, allocatable        :: local_aff(:,:), local_theta(:,:), local_shift_x(:,:), local_shift_y(:,:)
        integer, allocatable     :: region_labels(:), region_inds(:), region_pinds(:), labels(:)
        integer, allocatable     :: region_anchor(:)
        integer     :: n, nregions, iregion, nr, nclasses, first_anchor_local
        n = size(pinds)
        if( n <= 1 )then
            allocate(aff(n,n), theta(n,n), shift_x(n,n), shift_y(n,n), source=0.)
            return
        endif
        call build_projection_region_labels(params, spproj, pinds, region_labels, nregions)
        allocate(aff(n,n), theta(n,n), shift_x(n,n), shift_y(n,n), source=0.)
        allocate(all_gauge_angle(n), all_gauge_img_dist(n), all_gauge_shift_x(n), all_gauge_shift_y(n), source=0.)
        allocate(region_anchor(nregions), source=0)
        allocate(rep_imgs(nregions))
        do iregion = 1,nregions
            call inds_for_label(region_labels, iregion, region_inds)
            nr = size(region_inds)
            if( nr == 0 ) cycle
            allocate(region_pinds(nr))
            region_pinds = pinds(region_inds)
            if( nr <= 2 )then
                allocate(labels(nr))
                allocate(gauge_angle(nr), gauge_img_dist(nr), gauge_shift_x(nr), gauge_shift_y(nr), source=0.)
                labels = 1
                region_imgs = extract_imgarr(imgs, region_inds)
                call build_class_averages(region_imgs, labels, 1, class_avgs)
            else
                nclasses    = local_cluster2d_nclasses(params, nr)
                region_imgs = extract_imgarr(imgs, region_inds)
                call run_local_cluster2d(params, region_pinds, region_imgs, nclasses, labels, class_avgs, &
                    gauge_angle, gauge_shift_x, gauge_shift_y, gauge_img_dist)
                nclasses = maxval(labels)
                if( .not. allocated(class_avgs) .or. size(class_avgs) /= nclasses )then
                    if( allocated(class_avgs) ) call dealloc_imgarr(class_avgs)
                    call build_class_averages(region_imgs, labels, nclasses, class_avgs)
                endif
            endif
            call build_class_gauge_links(params, class_avgs, class_angle, class_shift_x, class_shift_y, class_img_dist)
            call build_gauge_knn_graph(params, spproj, region_pinds, gauge_angle, gauge_shift_x, gauge_shift_y, &
                gauge_img_dist, local_aff, local_theta, local_shift_x, local_shift_y, &
                labels, class_angle, class_shift_x, class_shift_y, class_img_dist)
            first_anchor_local = choose_orientation_medoid(spproj, region_pinds)
            call merge_local_affinity(region_inds, local_aff, local_theta, local_shift_x, local_shift_y, aff, theta, shift_x, shift_y)
            all_gauge_angle(region_inds)    = gauge_angle
            all_gauge_img_dist(region_inds) = gauge_img_dist
            all_gauge_shift_x(region_inds)  = gauge_shift_x
            all_gauge_shift_y(region_inds)  = gauge_shift_y
            region_anchor(iregion) = region_inds(first_anchor_local)
            call rep_imgs(iregion)%copy(class_avgs(labels(first_anchor_local)))
            if( allocated(gauge_angle)    ) deallocate(gauge_angle)
            if( allocated(gauge_img_dist) ) deallocate(gauge_img_dist)
            if( allocated(gauge_shift_x)  ) deallocate(gauge_shift_x)
            if( allocated(gauge_shift_y)  ) deallocate(gauge_shift_y)
            if( allocated(class_angle)    ) deallocate(class_angle)
            if( allocated(class_shift_x)  ) deallocate(class_shift_x)
            if( allocated(class_shift_y)  ) deallocate(class_shift_y)
            if( allocated(class_img_dist) ) deallocate(class_img_dist)
            if( allocated(labels)         ) deallocate(labels)
            if( allocated(class_avgs)     ) call dealloc_imgarr(class_avgs)
            if( allocated(region_imgs)    ) call dealloc_imgarr(region_imgs)
            if( allocated(local_aff)      ) deallocate(local_aff)
            if( allocated(local_theta)    ) deallocate(local_theta)
            if( allocated(local_shift_x)  ) deallocate(local_shift_x)
            if( allocated(local_shift_y)  ) deallocate(local_shift_y)
            if( allocated(region_pinds)   ) deallocate(region_pinds)
            if( allocated(region_inds)    ) deallocate(region_inds)
        end do
        call connect_region_anchors(params, spproj, pinds, region_anchor, rep_imgs, all_gauge_angle, all_gauge_shift_x, &
            all_gauge_shift_y, all_gauge_img_dist, aff, theta, shift_x, shift_y)
        if( allocated(region_labels)      ) deallocate(region_labels)
        if( allocated(all_gauge_angle)    ) deallocate(all_gauge_angle)
        if( allocated(all_gauge_img_dist) ) deallocate(all_gauge_img_dist)
        if( allocated(all_gauge_shift_x)  ) deallocate(all_gauge_shift_x)
        if( allocated(all_gauge_shift_y)  ) deallocate(all_gauge_shift_y)
        if( allocated(region_anchor)      ) deallocate(region_anchor)
        if( allocated(rep_imgs)           ) call dealloc_imgarr(rep_imgs)
    end subroutine build_cluster2d_affinity

    subroutine build_projection_region_labels( params, spproj, pinds, labels, nregions )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        integer, allocatable, intent(out)   :: labels(:)
        integer,              intent(out)   :: nregions
        integer, allocatable                :: raw_labels(:)
        call assign_projection_coordinate_regions(spproj, pinds, raw_labels)
        call connectedize_region_labels(params, spproj, pinds, raw_labels, labels, nregions)
        if( allocated(raw_labels) ) deallocate(raw_labels)
    end subroutine build_projection_region_labels

    subroutine assign_projection_coordinate_regions( spproj, pinds, labels )
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        integer, allocatable, intent(out)   :: labels(:)
        integer, allocatable :: center_refs(:)
        type(ori) :: ptcl_ori, ref_ori
        real      :: dist, best_dist
        integer   :: n, nrefs, ncenters, iref, icen, i, best_center
        n = size(pinds)
        allocate(labels(n), source=1)
        if( n <= 1 ) return
        nrefs = spproj%os_cls3D%get_noris()
        if( nrefs < 1 ) THROW_HARD('cls3D segment required for SO3 projection coordinates')
        if( .not. spproj%os_cls3D%isthere('e1') .or. .not. spproj%os_cls3D%isthere('e2') )then
            THROW_HARD('cls3D e1/e2 coordinates required for SO3 projection regions')
        endif
        allocate(center_refs(nrefs), source=0)
        ncenters = 0
        do iref = 1,nrefs
            if( spproj%os_cls3D%isthere(iref, 'state') )then
                if( spproj%os_cls3D%get_state(iref) <= 0 ) cycle
            endif
            ncenters = ncenters + 1
            center_refs(ncenters) = iref
        end do
        if( ncenters < 1 ) THROW_HARD('no active cls3D projection coordinates for SO3 projection regions')
        do i = 1,n
            call spproj%os_ptcl3D%get_ori(pinds(i), ptcl_ori)
            call ptcl_ori%e3set(0.)
            best_center = 1
            best_dist   = huge(1.)
            do icen = 1,ncenters
                call spproj%os_cls3D%get_ori(center_refs(icen), ref_ori)
                call ref_ori%e3set(0.)
                dist = ptcl_ori .euldist. ref_ori
                if( dist < best_dist )then
                    best_dist   = dist
                    best_center = icen
                endif
                call ref_ori%kill
            end do
            labels(i) = best_center
            call ptcl_ori%kill
        end do
        if( allocated(center_refs) ) deallocate(center_refs)
    end subroutine assign_projection_coordinate_regions

    subroutine connectedize_region_labels( params, spproj, pinds, raw_labels, labels, nregions )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:), raw_labels(:)
        integer, allocatable, intent(out)   :: labels(:)
        integer,              intent(out)   :: nregions
        logical, allocatable :: edge_mask(:,:)
        call build_projection_knn_graph(spproj, max(1, params%k_nn), pinds, edge_mask)
        call split_labels_by_connectivity(raw_labels, edge_mask, labels, nregions)
        if( allocated(edge_mask) ) deallocate(edge_mask)
    end subroutine connectedize_region_labels

    subroutine split_labels_by_connectivity( raw_labels, edge_mask, labels, nregions )
        integer,              intent(in)  :: raw_labels(:)
        logical,              intent(in)  :: edge_mask(:,:)
        integer, allocatable, intent(out) :: labels(:)
        integer,              intent(out) :: nregions
        logical, allocatable :: visited(:)
        integer, allocatable :: stack(:)
        integer :: n, seed, v, u, top, raw_label
        n = size(raw_labels)
        allocate(labels(n), stack(n))
        allocate(visited(n))
        labels = 0
        stack  = 0
        visited = .false.
        nregions = 0
        do seed = 1,n
            if( visited(seed) ) cycle
            nregions  = nregions + 1
            raw_label = raw_labels(seed)
            top = 1
            stack(top) = seed
            visited(seed) = .true.
            labels(seed)  = nregions
            do while( top > 0 )
                v = stack(top)
                top = top - 1
                do u = 1,n
                    if( visited(u) ) cycle
                    if( raw_labels(u) /= raw_label ) cycle
                    if( .not. edge_mask(v,u) ) cycle
                    visited(u) = .true.
                    labels(u)  = nregions
                    top = top + 1
                    stack(top) = u
                end do
            end do
        end do
        if( allocated(visited) ) deallocate(visited)
        if( allocated(stack)   ) deallocate(stack)
    end subroutine split_labels_by_connectivity

    subroutine inds_for_label( labels, label, inds )
        integer,              intent(in)  :: labels(:), label
        integer, allocatable, intent(out) :: inds(:)
        integer :: i, cnt
        allocate(inds(count(labels == label)))
        cnt = 0
        do i = 1,size(labels)
            if( labels(i) == label )then
                cnt = cnt + 1
                inds(cnt) = i
            endif
        end do
    end subroutine inds_for_label

    integer function local_cluster2d_nclasses( params, n ) result(nclasses)
        type(parameters), intent(in) :: params
        integer,          intent(in) :: n
        if( n <= 2 )then
            nclasses = 1
        else
            nclasses = min(max(2, default_local_split_nclasses(params)), n - 1)
        endif
    end function local_cluster2d_nclasses

    integer function default_local_split_nclasses( params ) result(nclasses)
        type(parameters), intent(in) :: params
        if( params%ncls > 1 .and. params%ncls <= max(params%nsubcls_max, params%nsubcls_min) )then
            nclasses = params%ncls
        else
            nclasses = max(2, params%nsubcls_min)
        endif
    end function default_local_split_nclasses

    subroutine merge_local_affinity( inds, local_aff, local_theta, local_shift_x, local_shift_y, aff, theta, shift_x, shift_y )
        integer, intent(in)    :: inds(:)
        real,    intent(in)    :: local_aff(:,:), local_theta(:,:), local_shift_x(:,:), local_shift_y(:,:)
        real,    intent(inout) :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        integer :: i, j, gi, gj
        do i = 1,size(inds)
            gi = inds(i)
            do j = 1,size(inds)
                if( local_aff(i,j) <= 0. ) cycle
                gj = inds(j)
                aff(gi,gj)   = max(aff(gi,gj), local_aff(i,j))
                theta(gi,gj) = local_theta(i,j)
                shift_x(gi,gj) = local_shift_x(i,j)
                shift_y(gi,gj) = local_shift_y(i,j)
            end do
        end do
    end subroutine merge_local_affinity

    subroutine connect_region_anchors( params, spproj, pinds, region_anchor, rep_imgs, gauge_angle, gauge_shift_x, gauge_shift_y, &
            gauge_img_dist, aff, theta, shift_x, shift_y )
        type(parameters), intent(in)    :: params
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:), region_anchor(:)
        type(image),      intent(inout) :: rep_imgs(:)
        real,             intent(in)    :: gauge_angle(:), gauge_shift_x(:), gauge_shift_y(:), gauge_img_dist(:)
        real,             intent(inout) :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        type(polarft_calc) :: pftc
        type(parameters), target :: params_pft
        integer, allocatable :: edge_i(:), edge_j(:)
        real, allocatable :: edge_img_dist(:), edge_angle(:), edge_shift_x(:), edge_shift_y(:)
        real :: corr, img_dist, angle, sh(2), eps_img, w
        integer :: nregions, nedges, e, gi, gj
        nregions = size(region_anchor)
        if( nregions <= 1 ) return
        call build_region_mst(spproj, pinds, region_anchor, edge_i, edge_j, nedges)
        if( nedges <= 0 )then
            if( allocated(edge_i) ) deallocate(edge_i)
            if( allocated(edge_j) ) deallocate(edge_j)
            return
        endif
        allocate(edge_img_dist(nedges), edge_angle(nedges), edge_shift_x(nedges), edge_shift_y(nedges), source=0.)
        call prepare_so3_pft(params, rep_imgs, params_pft, pftc)
        eps_img = 0.
        do e = 1,nedges
            call align_pair(params, pftc, edge_i(e), edge_j(e), corr, img_dist, angle, sh)
            gi = region_anchor(edge_i(e))
            gj = region_anchor(edge_j(e))
            edge_img_dist(e) = max(img_dist, 0.) + 0.5 * (max(gauge_img_dist(gi), 0.) + max(gauge_img_dist(gj), 0.))
            edge_angle(e)    = wrap_angle_pm180(angle + gauge_angle(gj) - gauge_angle(gi))
            edge_shift_x(e)  = sh(1) + gauge_shift_x(gj) - gauge_shift_x(gi)
            edge_shift_y(e)  = sh(2) + gauge_shift_y(gj) - gauge_shift_y(gi)
            eps_img = eps_img + edge_img_dist(e)
        end do
        eps_img = max(eps_img / real(max(nedges,1)), 1.e-6)
        do e = 1,nedges
            gi = region_anchor(edge_i(e))
            gj = region_anchor(edge_j(e))
            w = exp(-edge_img_dist(e) / eps_img)
            aff(gi,gj)   = max(aff(gi,gj), w)
            aff(gj,gi)   = max(aff(gj,gi), w)
            theta(gi,gj) = deg2rad(edge_angle(e))
            theta(gj,gi) = -theta(gi,gj)
            shift_x(gi,gj) = edge_shift_x(e)
            shift_x(gj,gi) = -edge_shift_x(e)
            shift_y(gi,gj) = edge_shift_y(e)
            shift_y(gj,gi) = -edge_shift_y(e)
        end do
        call pftc%kill
        if( allocated(edge_i)        ) deallocate(edge_i)
        if( allocated(edge_j)        ) deallocate(edge_j)
        if( allocated(edge_img_dist) ) deallocate(edge_img_dist)
        if( allocated(edge_angle)    ) deallocate(edge_angle)
        if( allocated(edge_shift_x)  ) deallocate(edge_shift_x)
        if( allocated(edge_shift_y)  ) deallocate(edge_shift_y)
    end subroutine connect_region_anchors

    subroutine build_region_mst( spproj, pinds, region_anchor, edge_i, edge_j, nedges )
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:), region_anchor(:)
        integer, allocatable, intent(out)   :: edge_i(:), edge_j(:)
        integer,              intent(out)   :: nedges
        logical, allocatable :: in_tree(:)
        real, allocatable :: best_dist(:)
        integer, allocatable :: best_parent(:)
        real :: dist, min_dist
        integer :: nregions, i, j, next_region
        nregions = size(region_anchor)
        allocate(edge_i(max(0,nregions-1)), edge_j(max(0,nregions-1)))
        edge_i = 0
        edge_j = 0
        nedges = 0
        if( nregions <= 1 ) return
        allocate(in_tree(nregions), best_dist(nregions), best_parent(nregions))
        in_tree = .false.
        best_dist = huge(best_dist)
        best_parent = 0
        in_tree(1) = .true.
        do j = 2,nregions
            best_dist(j)   = projection_dir_dist(spproj, pinds(region_anchor(1)), pinds(region_anchor(j)))
            best_parent(j) = 1
        end do
        do while( count(in_tree) < nregions )
            min_dist = huge(min_dist)
            next_region = 0
            do i = 1,nregions
                if( in_tree(i) ) cycle
                if( best_dist(i) < min_dist )then
                    min_dist = best_dist(i)
                    next_region = i
                endif
            end do
            if( next_region == 0 ) exit
            nedges = nedges + 1
            edge_i(nedges) = best_parent(next_region)
            edge_j(nedges) = next_region
            in_tree(next_region) = .true.
            do j = 1,nregions
                if( in_tree(j) ) cycle
                dist = projection_dir_dist(spproj, pinds(region_anchor(next_region)), pinds(region_anchor(j)))
                if( dist < best_dist(j) )then
                    best_dist(j)   = dist
                    best_parent(j) = next_region
                endif
            end do
        end do
        if( allocated(in_tree)     ) deallocate(in_tree)
        if( allocated(best_dist)   ) deallocate(best_dist)
        if( allocated(best_parent) ) deallocate(best_parent)
    end subroutine build_region_mst

    subroutine build_projection_knn_graph( spproj, k_nn, pinds, edge_mask )
        type(sp_project), intent(inout)          :: spproj
        integer,          intent(in)             :: k_nn, pinds(:)
        logical, allocatable, intent(out)        :: edge_mask(:,:)
        integer, allocatable :: nbrs(:,:)
        real, allocatable    :: dists(:,:)
        real :: dist
        integer :: n, i, j, m, k_used
        n = size(pinds)
        allocate(edge_mask(n,n), source=.false.)
        if( n <= 1 ) return
        k_used = min(max(1, k_nn), n-1)
        allocate(nbrs(k_used,n), source=0)
        allocate(dists(k_used,n), source=huge(1.))
        do i = 1,n
            do j = 1,n
                if( i == j ) cycle
                dist = projection_dir_dist(spproj, pinds(i), pinds(j))
                call insert_neighbor(j, dist, nbrs(:,i), dists(:,i))
            end do
        end do
        do i = 1,n
            do m = 1,k_used
                j = nbrs(m,i)
                if( j < 1 ) cycle
                edge_mask(i,j) = .true.
                edge_mask(j,i) = .true.
            end do
        end do
        if( allocated(nbrs) ) deallocate(nbrs)
        if( allocated(dists) ) deallocate(dists)
    end subroutine build_projection_knn_graph

    real function projection_dir_dist( spproj, pind_i, pind_j ) result(dist)
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pind_i, pind_j
        type(ori) :: oi, oj
        call spproj%os_ptcl3D%get_ori(pind_i, oi)
        call spproj%os_ptcl3D%get_ori(pind_j, oj)
        call oi%e3set(0.)
        call oj%e3set(0.)
        dist = oi .euldist. oj
        call oi%kill
        call oj%kill
    end function projection_dir_dist

    subroutine run_local_cluster2d( params, pinds, imgs, nclasses_req, labels, class_avgs, gauge_angle, gauge_shift_x, &
            gauge_shift_y, gauge_img_dist )
        type(parameters),     intent(in)    :: params
        integer,              intent(in)    :: pinds(:), nclasses_req
        type(image),          intent(inout) :: imgs(:)
        integer, allocatable, intent(out)   :: labels(:)
        type(image), allocatable, intent(out) :: class_avgs(:)
        real, allocatable,    intent(out)   :: gauge_angle(:), gauge_shift_x(:), gauge_shift_y(:), gauge_img_dist(:)
        type(cluster2D_inmem_strategy) :: c2d_strategy
        type(parameters) :: c2d_params
        type(builder)    :: c2d_build
        type(cmdline)    :: c2d_cline
        type(sp_project) :: local_proj, result_proj
        type(oris)       :: local_os
        type(ori)        :: o
        type(ctfparams)  :: ctfvars
        type(string) :: cwd_before, cwd_now, workdir, stk_fname, projfile, cavg_fname, dot
        integer :: n, i, iter, ldim(3), istat, next_dir, nclasses, nlocal_maxits
        logical :: converged
        real :: smpd, mskdiam_local, corr
        n    = size(pinds)
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        call simple_getcwd(cwd_before)
        dot      = '.'
        next_dir = find_next_int_dir_prefix(dot)
        workdir  = int2str(next_dir)//'_so3_local_cluster2d_'//int2str(pinds(1))
        call simple_mkdir(workdir, verbose=.false.)
        call simple_chdir(workdir, status=istat)
        if( istat /= 0 ) THROW_HARD('failed to enter local SO3 cluster2D directory')
        call simple_getcwd(cwd_now)
        CWD_GLOB = cwd_now%to_char()

        stk_fname = 'particles.mrcs'
        projfile  = 'local_cluster2d.simple'
        call write_imgarr(imgs, stk_fname)
        call local_proj%update_projinfo(projfile)
        call local_os%new(n, is_ptcl=.true.)
        do i = 1,n
            call o%new(is_ptcl=.true.)
            call o%set('state', 1.)
            call o%set('w',     1.)
            call o%set('x',     0.)
            call o%set('y',     0.)
            call o%e3set(0.)
            call local_os%set_ori(i, o)
            call o%kill
        end do
        ctfvars%ctfflag = CTFFLAG_NO
        ctfvars%smpd    = smpd
        ctfvars%kv      = max(params%kv, 300.)
        ctfvars%cs      = max(params%cs, 2.7)
        ctfvars%fraca   = max(params%fraca, 0.1)
        call local_proj%add_single_stk(stk_fname, ctfvars, local_os)
        call local_proj%write(projfile)

        mskdiam_local = params%mskdiam
        if( mskdiam_local <= 0. ) mskdiam_local = real(ldim(1)) * smpd
        nlocal_maxits = max(1, params%so3_local_cluster2d_maxits)
        call c2d_cline%set('prg',           'cluster2D')
        call c2d_cline%set('projfile',      projfile%to_char())
        call c2d_cline%set('ncls',          nclasses_req)
        call c2d_cline%set('mskdiam',       mskdiam_local)
        call c2d_cline%set('smpd',          smpd)
        call c2d_cline%set('box',           ldim(1))
        call c2d_cline%set('ctf',           'no')
        call c2d_cline%set('objfun',        'cc')
        call c2d_cline%set('ml_reg',        'no')
        call c2d_cline%set('maxits',        nlocal_maxits)
        call c2d_cline%set('startit',       1)
        call c2d_cline%set('nthr',          max(1, params%nthr))
        call c2d_cline%set('trs',           params%trs)
        call c2d_cline%set('mkdir',         'no')
        call c2d_cline%set('restore_cavgs', 'yes')
        call set_cluster2D_defaults(c2d_cline)
        call c2d_strategy%initialize(c2d_params, c2d_build, c2d_cline)
        c2d_params%which_iter = c2d_params%startit - 1
        if( c2d_cline%defined('extr_iter') )then
            c2d_params%extr_iter = c2d_params%extr_iter - 1
        else
            c2d_params%extr_iter = c2d_params%startit - 1
        endif
        do iter = 1,nlocal_maxits
            c2d_params%which_iter = c2d_params%which_iter + 1
            c2d_params%extr_iter  = c2d_params%extr_iter  + 1
            call c2d_strategy%execute_iteration(c2d_params, c2d_build, c2d_cline, converged)
            call c2d_strategy%finalize_iteration(c2d_params, c2d_build)
            if( converged ) exit
        end do
        call c2d_strategy%finalize_run(c2d_params, c2d_build, c2d_cline)
        call c2d_strategy%cleanup(c2d_params)
        call result_proj%read(projfile)
        allocate(labels(n), source=0)
        allocate(gauge_angle(n), gauge_shift_x(n), gauge_shift_y(n), gauge_img_dist(n), source=0.)
        do i = 1,n
            labels(i) = max(1, result_proj%os_ptcl2D%get_int(i, 'class'))
            gauge_angle(i) = result_proj%os_ptcl2D%e3get(i)
            if( result_proj%os_ptcl2D%isthere('x') ) gauge_shift_x(i) = result_proj%os_ptcl2D%get(i, 'x')
            if( result_proj%os_ptcl2D%isthere('y') ) gauge_shift_y(i) = result_proj%os_ptcl2D%get(i, 'y')
            corr = 0.
            if( result_proj%os_ptcl2D%isthere('corr') ) corr = result_proj%os_ptcl2D%get(i, 'corr')
            corr = min(1., max(-1., corr))
            gauge_img_dist(i) = max(0., 1. - corr)
        end do
        nclasses = maxval(labels)
        cavg_fname = CAVGS_ITER_FBODY//int2str_pad(c2d_params%which_iter,3)//MRC_EXT
        if( file_exists(cavg_fname) ) class_avgs = read_stk_into_imgarr(cavg_fname)
        if( allocated(class_avgs) .and. size(class_avgs) /= nclasses ) call dealloc_imgarr(class_avgs)
        call result_proj%kill
        call local_proj%kill
        call local_os%kill
        call c2d_build%kill_general_tbox()
        call c2d_build%kill_strategy2D_tbox()
        call c2d_cline%kill()
        call simple_chdir(cwd_before, status=istat)
        if( istat /= 0 ) THROW_HARD('failed to restore cwd after local SO3 cluster2D')
        call simple_getcwd(cwd_now)
        CWD_GLOB = cwd_now%to_char()
    end subroutine run_local_cluster2d

    subroutine build_class_averages( imgs, labels, nclasses, class_avgs )
        type(image),          intent(inout) :: imgs(:)
        integer,              intent(in)    :: labels(:), nclasses
        type(image), allocatable, intent(out) :: class_avgs(:)
        integer, allocatable :: class_pops(:)
        integer :: ldim(3), i, icls
        real :: smpd
        ldim = imgs(1)%get_ldim()
        smpd = imgs(1)%get_smpd()
        allocate(class_avgs(nclasses))
        allocate(class_pops(nclasses), source=0)
        do icls = 1,nclasses
            call class_avgs(icls)%new(ldim, smpd)
            call class_avgs(icls)%zero_and_unflag_ft
        end do
        do i = 1,size(imgs)
            icls = labels(i)
            class_pops(icls) = class_pops(icls) + 1
            call class_avgs(icls)%add(imgs(i))
        end do
        do icls = 1,nclasses
            if( class_pops(icls) > 0 )then
                call class_avgs(icls)%div(real(class_pops(icls)))
            else
                call class_avgs(icls)%copy(imgs(1))
            endif
        end do
        if( allocated(class_pops) ) deallocate(class_pops)
    end subroutine build_class_averages

    subroutine build_class_gauge_links( params, class_avgs, class_angle, class_shift_x, class_shift_y, class_img_dist )
        type(parameters),  intent(in)    :: params
        type(image),       intent(inout) :: class_avgs(:)
        real, allocatable, intent(out)   :: class_angle(:,:), class_shift_x(:,:), class_shift_y(:,:), class_img_dist(:,:)
        type(polarft_calc) :: pftc
        type(parameters), target :: params_pft
        real :: corr, img_dist, angle, sh(2)
        integer :: nclasses, icls, jcls
        nclasses = size(class_avgs)
        allocate(class_angle(nclasses,nclasses), class_shift_x(nclasses,nclasses), class_shift_y(nclasses,nclasses), &
            class_img_dist(nclasses,nclasses), source=0.)
        if( nclasses <= 1 ) return
        call prepare_so3_pft(params, class_avgs, params_pft, pftc)
        do icls = 1,nclasses-1
            do jcls = icls+1,nclasses
                call align_pair(params, pftc, icls, jcls, corr, img_dist, angle, sh)
                class_angle(icls,jcls)    = angle
                class_angle(jcls,icls)    = -angle
                class_shift_x(icls,jcls)  = sh(1)
                class_shift_x(jcls,icls)  = -sh(1)
                class_shift_y(icls,jcls)  = sh(2)
                class_shift_y(jcls,icls)  = -sh(2)
                class_img_dist(icls,jcls) = img_dist
                class_img_dist(jcls,icls) = img_dist
            end do
        end do
        call pftc%kill
    end subroutine build_class_gauge_links

    subroutine build_gauge_knn_graph( params, spproj, pinds, gauge_angle, gauge_shift_x, gauge_shift_y, gauge_img_dist, &
            aff, theta, shift_x, shift_y, class_labels, class_angle, class_shift_x, class_shift_y, class_img_dist )
        type(parameters),     intent(in)    :: params
        type(sp_project),     intent(inout) :: spproj
        integer,              intent(in)    :: pinds(:)
        real,                 intent(in)    :: gauge_angle(:), gauge_shift_x(:), gauge_shift_y(:), gauge_img_dist(:)
        real, allocatable,    intent(out)   :: aff(:,:), theta(:,:), shift_x(:,:), shift_y(:,:)
        integer, optional,    intent(in)    :: class_labels(:)
        real,    optional,    intent(in)    :: class_angle(:,:), class_shift_x(:,:), class_shift_y(:,:), class_img_dist(:,:)
        logical, allocatable :: edge_mask(:,:)
        real, allocatable :: so3_mat(:,:)
        real :: edge_img_dist, edge_angle, edge_shift_x, edge_shift_y, eps_img, eps_so3, w
        integer :: n, i, j, ci, cj, nedges
        logical :: use_cluster_gauge
        n = size(pinds)
        use_cluster_gauge = present(class_labels) .and. present(class_angle) .and. present(class_shift_x) .and. &
            present(class_shift_y) .and. present(class_img_dist)
        call build_orientation_knn_graph(spproj, params%k_nn, pinds, edge_mask, so3_mat)
        eps_img = 0.
        eps_so3 = 0.
        nedges  = 0
        do i = 1,n-1
            do j = i+1,n
                if( .not. edge_mask(i,j) ) cycle
                call gauge_edge_terms(i, j, edge_img_dist, edge_angle, edge_shift_x, edge_shift_y)
                eps_img = eps_img + edge_img_dist
                eps_so3 = eps_so3 + so3_mat(i,j)**2
                nedges  = nedges + 1
            end do
        end do
        eps_img = max(eps_img / real(max(nedges,1)), 1.e-6)
        eps_so3 = max(eps_so3 / real(max(nedges,1)), 1.e-6)
        allocate(aff(n,n), theta(n,n), shift_x(n,n), shift_y(n,n), source=0.)
        do i = 1,n-1
            do j = i+1,n
                if( .not. edge_mask(i,j) ) cycle
                call gauge_edge_terms(i, j, edge_img_dist, edge_angle, edge_shift_x, edge_shift_y)
                w = exp(-edge_img_dist / eps_img) * exp(-(so3_mat(i,j)**2) / eps_so3)
                aff(i,j)   = w
                aff(j,i)   = w
                theta(i,j) = deg2rad(edge_angle)
                theta(j,i) = -theta(i,j)
                shift_x(i,j) = edge_shift_x
                shift_x(j,i) = -edge_shift_x
                shift_y(i,j) = edge_shift_y
                shift_y(j,i) = -edge_shift_y
            end do
        end do
        if( allocated(edge_mask) ) deallocate(edge_mask)
        if( allocated(so3_mat)   ) deallocate(so3_mat)
        contains
            subroutine gauge_edge_terms( ii, jj, edge_img_dist_out, edge_angle_out, edge_shift_x_out, edge_shift_y_out )
                integer, intent(in) :: ii, jj
                real, intent(out)   :: edge_img_dist_out, edge_angle_out, edge_shift_x_out, edge_shift_y_out
                edge_img_dist_out = 0.5 * (max(gauge_img_dist(ii), 0.) + max(gauge_img_dist(jj), 0.))
                edge_angle_out    = wrap_angle_pm180(gauge_angle(jj) - gauge_angle(ii))
                edge_shift_x_out  = gauge_shift_x(jj) - gauge_shift_x(ii)
                edge_shift_y_out  = gauge_shift_y(jj) - gauge_shift_y(ii)
                if( use_cluster_gauge )then
                    ci = class_labels(ii)
                    cj = class_labels(jj)
                    if( ci /= cj )then
                        edge_img_dist_out = edge_img_dist_out + max(class_img_dist(ci,cj), 0.)
                        edge_angle_out    = wrap_angle_pm180(class_angle(ci,cj) + gauge_angle(jj) - gauge_angle(ii))
                        edge_shift_x_out  = class_shift_x(ci,cj) + gauge_shift_x(jj) - gauge_shift_x(ii)
                        edge_shift_y_out  = class_shift_y(ci,cj) + gauge_shift_y(jj) - gauge_shift_y(ii)
                    endif
                endif
            end subroutine gauge_edge_terms
    end subroutine build_gauge_knn_graph

    subroutine build_orientation_knn_graph( spproj, k_nn, pinds, edge_mask, so3_mat )
        type(sp_project), intent(inout)          :: spproj
        integer,          intent(in)             :: k_nn, pinds(:)
        logical, allocatable, intent(out)        :: edge_mask(:,:)
        real,    allocatable, intent(out)        :: so3_mat(:,:)
        integer, allocatable :: nbrs(:,:)
        real, allocatable    :: so3_dists(:,:)
        integer :: n, i, j, m, k_used
        n = size(pinds)
        allocate(edge_mask(n,n), source=.false.)
        allocate(so3_mat(n,n), source=0.)
        k_used = min(max(1, k_nn), n-1)
        call find_so3_neighbors(spproj, pinds, k_used, nbrs, so3_dists)
        do i = 1,n
            do m = 1,k_used
                j = nbrs(m,i)
                if( j < 1 ) cycle
                edge_mask(i,j) = .true.
                edge_mask(j,i) = .true.
                so3_mat(i,j) = so3_dists(m,i)
                so3_mat(j,i) = so3_dists(m,i)
            end do
        end do
        if( allocated(nbrs) ) deallocate(nbrs)
        if( allocated(so3_dists) ) deallocate(so3_dists)
    end subroutine build_orientation_knn_graph

    subroutine find_so3_neighbors( spproj, pinds, k_used, nbrs, so3_dists )
        type(sp_project), intent(inout)      :: spproj
        integer, intent(in)                  :: pinds(:), k_used
        integer, allocatable, intent(out)    :: nbrs(:,:)
        real,    allocatable, intent(out)    :: so3_dists(:,:)
        type(ori) :: oi, oj
        real :: dist
        integer :: i, j
        allocate(nbrs(k_used,size(pinds)), source=0)
        allocate(so3_dists(k_used,size(pinds)), source=0.)
        so3_dists = huge(so3_dists)
        do i = 1,size(pinds)
            call spproj%os_ptcl3D%get_ori(pinds(i), oi)
            do j = 1,size(pinds)
                if( i == j ) cycle
                call spproj%os_ptcl3D%get_ori(pinds(j), oj)
                dist = oi .euldist. oj
                call insert_neighbor(j, dist, nbrs(:,i), so3_dists(:,i))
                call oj%kill
            end do
            call oi%kill
        end do
    end subroutine find_so3_neighbors

    subroutine insert_neighbor( ind, dist, inds, dists )
        integer, intent(in)    :: ind
        real,    intent(in)    :: dist
        integer, intent(inout) :: inds(:)
        real,    intent(inout) :: dists(:)
        integer :: pos, k
        pos = 0
        do k = 1,size(dists)
            if( dist < dists(k) )then
                pos = k
                exit
            endif
        end do
        if( pos == 0 ) return
        do k = size(dists),pos+1,-1
            dists(k) = dists(k-1)
            inds(k)  = inds(k-1)
        end do
        dists(pos) = dist
        inds(pos)  = ind
    end subroutine insert_neighbor

    subroutine prepare_so3_pft( params, imgs, params_pft, pftc )
        type(parameters),  intent(in)    :: params
        type(image),       intent(inout) :: imgs(:)
        type(parameters), target, intent(out) :: params_pft
        type(polarft_calc), intent(inout) :: pftc
        type(image), allocatable :: imgs_work(:)
        integer :: ldim(3), kfromto(2), pdim_srch(3), i
        real    :: smpd
        params_pft             = params
        params_pft%objfun      = 'cc'
        params_pft%cc_objfun   = OBJFUN_CC
        params_pft%ctf         = 'no'
        ldim                   = imgs(1)%get_ldim()
        smpd                   = imgs(1)%get_smpd()
        params_pft%ldim        = ldim
        params_pft%box         = ldim(1)
        params_pft%box_crop    = ldim(1)
        params_pft%smpd        = smpd
        params_pft%smpd_crop   = smpd
        if( params_pft%pftsz < 1 ) params_pft%pftsz = magic_pftsz((real(ldim(1)) - COSMSKHALFWIDTH) / 2.)
        kfromto(1) = max(2, calc_fourier_index(params_pft%hp, params_pft%box, params_pft%smpd))
        kfromto(2) =        calc_fourier_index(params_pft%lp, params_pft%box, params_pft%smpd)
        if( kfromto(2) < kfromto(1) ) THROW_HARD('invalid Fourier limits for SO3 graph registration')
        call pftc%new(params_pft, size(imgs), [1,size(imgs)], kfromto)
        pdim_srch = pftc%get_pdim_srch()
        imgs_work = copy_imgarr(imgs)
        call imgs_work(1)%memoize4polarize(pdim_srch)
        do i = 1,size(imgs_work)
            call imgs_work(i)%fft()
            call pftc%polarize_ref_pft( imgs_work(i), i, iseven=.true., pdim=pdim_srch, oversamp=.false.)
            call pftc%polarize_ptcl_pft(imgs_work(i), i,               pdim=pdim_srch, oversamp=.false.)
            call imgs_work(i)%ifft()
        end do
        call pftc%memoize_refs
        call pftc%memoize_ptcls
        call dealloc_imgarr(imgs_work)
    end subroutine prepare_so3_pft

    subroutine align_pair( params, pftc, iref, iptcl_loc, corr, img_dist, angle, sh )
        type(parameters),  intent(in)    :: params
        type(polarft_calc), intent(inout) :: pftc
        integer,           intent(in)    :: iref, iptcl_loc
        real,              intent(out)   :: corr, img_dist, angle, sh(2)
        real, allocatable :: vals(:)
        real :: best_corr, trial_shift(2)
        integer :: nrots, loc, best_loc, ix, iy, trs_i
        nrots = pftc%get_nrots()
        allocate(vals(nrots))
        best_corr = -huge(best_corr)
        best_loc  = 1
        sh        = 0.
        trs_i     = nint(params%trs)
        do ix = -trs_i, trs_i
            do iy = -trs_i, trs_i
                trial_shift = [real(ix), real(iy)]
                call pftc%gen_objfun_vals(iref, iptcl_loc, trial_shift, vals)
                loc = maxloc(vals, dim=1)
                if( vals(loc) > best_corr )then
                    best_corr = vals(loc)
                    best_loc  = loc
                    sh        = trial_shift
                endif
            end do
        end do
        corr = best_corr
        if( .not. ieee_is_finite(corr) ) corr = 0.
        corr     = min(1., max(-1., corr))
        img_dist = max(0., 1. - corr)
        angle    = pftc%get_rot(best_loc)
        deallocate(vals)
    end subroutine align_pair

    integer function choose_orientation_medoid( spproj, pinds ) result(center_ind)
        type(sp_project), intent(inout) :: spproj
        integer,          intent(in)    :: pinds(:)
        type(ori) :: oi, oj
        real      :: dist, dsum, best_sum
        integer   :: i, j, n
        n = size(pinds)
        center_ind = 1
        if( n <= 1 ) return
        best_sum = huge(best_sum)
        do i = 1,n
            call spproj%os_ptcl3D%get_ori(pinds(i), oi)
            dsum = 0.
            do j = 1,n
                if( i == j ) cycle
                call spproj%os_ptcl3D%get_ori(pinds(j), oj)
                dist = oi .euldist. oj
                dsum = dsum + dist
                call oj%kill
            end do
            if( dsum < best_sum )then
                best_sum   = dsum
                center_ind = i
            endif
            call oi%kill
        end do
    end function choose_orientation_medoid

    real function wrap_angle_pm180( angle_deg ) result(wrapped)
        real, intent(in) :: angle_deg
        wrapped = modulo(angle_deg + 180., 360.) - 180.
    end function wrap_angle_pm180

end module simple_diff_map_graphs
