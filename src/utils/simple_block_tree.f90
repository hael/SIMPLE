module simple_block_tree
use simple_core_module_api
use simple_srchspace_map, only: srchspace_map
use simple_multi_dendro,  only: multi_dendro
use simple_binary_tree,   only: bt_node
implicit none

public :: gen_eulspace_block_tree
public :: srch_eul_bl_tree_greedy
public :: srch_eul_bl_tree_stoch
public :: srch_eul_bl_tree_prob
private

integer, parameter :: LINK_SINGLE   = 1
integer, parameter :: LINK_COMPLETE = 2
integer, parameter :: LINK_AVERAGE  = 3
integer, parameter :: LEAF_BUDGET   = 128
integer, parameter :: MAX_LEVELS    = 64

contains

    function gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym) result(block_tree)
        class(oris), intent(in)    :: eulspace, eulspace_sub
        class(sym),  intent(inout) :: pgrpsym
        type(ori)                  :: o, o_sub, oi, oj, osym
        type(srchspace_map)        :: mapper
        type(multi_dendro)         :: block_tree
        integer, allocatable       :: labels(:), refs(:)
        real,    allocatable       :: distmat(:,:), sub_distmat(:,:)
        integer                    :: i, j, ntrees, itree, nrefs, nspace, nspace_sub
        real                       :: inplrotdist, dtmp
        nspace     = eulspace%get_noris()
        nspace_sub = eulspace_sub%get_noris()
        allocate(distmat(nspace_sub, nspace))
        !$omp parallel do default(shared) proc_bind(close) private(i,j,o,o_sub,osym,inplrotdist,dtmp) schedule(static)
        do i = 1, nspace_sub
            call eulspace_sub%get_ori(i, o_sub)
            do j = 1, nspace
                call eulspace%get_ori(j, o)
                call pgrpsym%sym_dists(o_sub, o, osym, dtmp, inplrotdist)
                distmat(i,j) = dtmp
            end do
        end do
        !$omp end parallel do
        call mapper%new(nspace, nspace_sub, distmat)
        labels = mapper%get_full2sub_map()
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,i,j,oi,oj,osym,inplrotdist,dtmp) schedule(static)
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            write(*,'(a,1x,i0,a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            allocate(sub_distmat(nrefs, nrefs), source=0.0)
            if (nrefs > 1) then
                do i = 1, nrefs - 1
                    call eulspace%get_ori(refs(i), oi)        
                    do j = i + 1, nrefs
                        call eulspace%get_ori(refs(j), oj)
                        call pgrpsym%sym_dists(oi, oj, osym, dtmp, inplrotdist)
                        sub_distmat(i,j) = dtmp
                        sub_distmat(j,i) = dtmp
                    end do
                end do
            end if
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(sub_distmat)
            if (allocated(refs)) deallocate(refs)
        end do
        !$omp end parallel do
        deallocate(distmat)
        if (allocated(labels)) deallocate(labels)
    end function gen_eulspace_block_tree

    subroutine srch_eul_bl_tree_greedy( otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        type(ori)     :: o, osym
        integer       :: inode, level, lidx, ridx, local_k, j
        real          :: dist_left, dist_right, inplrotdist, dist
        type(bt_node) :: node_L, node_R, node_cur
        inode = block_tree%get_root_node_idx(itree)
        do level = 1, MAX_LEVELS
            if (inode == 0) exit
            if (block_tree%is_leaf(itree, inode)) exit
            if (block_tree%get_subset_size(itree, inode) <= LEAF_BUDGET) exit
            call block_tree%get_children_idx(itree, inode, lidx, ridx)
            if (lidx == 0 .and. ridx == 0) exit
            dist_left  = huge(1.0)
            dist_right = huge(1.0)
            if (lidx /= 0) then
                node_L = block_tree%get_node_by_idx(itree, lidx)
                call eulspace%get_ori(node_L%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_left, inplrotdist)
            end if
            if (ridx /= 0) then
                node_R = block_tree%get_node_by_idx(itree, ridx)
                call eulspace%get_ori(node_R%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_right, inplrotdist)
            end if
            if (dist_left <= dist_right) then
                inode = lidx
            else
                inode = ridx
            end if
        end do
        ! ---- (final) brute-force scan the selected cluster ----
        dist_min = huge(1.0)
        best_ref = 0
        if (inode /= 0) then
            node_cur = block_tree%get_node_by_idx(itree, inode)
            do j = 1, size(node_cur%subset)
                local_k  = node_cur%subset(j)                              ! LOCAL index
                best_ref = block_tree%local_to_global_ref(itree, local_k)  ! GLOBAL ref id
                if (best_ref == 0) cycle
                call eulspace%get_ori(best_ref, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist, inplrotdist)
                if (dist < dist_min) then
                    dist_min = dist
                end if
            end do
        end if
    end subroutine srch_eul_bl_tree_greedy

    subroutine srch_eul_bl_tree_stoch( otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        type(ori)     :: o, osym
        integer       :: inode, level, lidx, ridx, local_k, j, max_levs, rnd_levs
        real          :: dist_left, dist_right, inplrotdist, dist
        type(bt_node) :: node_L, node_R, node_cur
        inode = block_tree%get_root_node_idx(itree)
        max_levs = block_tree%get_tree_height(itree)
        if( max_levs > 2)then
            rnd_levs = irnd_uni(max_levs)
            rnd_levs = max(2, rnd_levs) ! ensure at least 2 levels of descent
        else
            rnd_levs = max_levs
        end if
        do level = 1, rnd_levs
            if (inode == 0) exit
            if (block_tree%is_leaf(itree, inode)) exit
            if (block_tree%get_subset_size(itree, inode) <= LEAF_BUDGET) exit
            call block_tree%get_children_idx(itree, inode, lidx, ridx)
            if (lidx == 0 .and. ridx == 0) exit
            dist_left  = huge(1.0)
            dist_right = huge(1.0)
            if (lidx /= 0) then
                node_L = block_tree%get_node_by_idx(itree, lidx)
                call eulspace%get_ori(node_L%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_left, inplrotdist)
            end if
            if (ridx /= 0) then
                node_R = block_tree%get_node_by_idx(itree, ridx)
                call eulspace%get_ori(node_R%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_right, inplrotdist)
            end if
            if (dist_left <= dist_right) then
                inode = lidx
            else
                inode = ridx
            end if
        end do
        ! ---- (final) brute-force scan the selected cluster ----
        dist_min = huge(1.0)
        best_ref = 0
        if (inode /= 0) then
            node_cur = block_tree%get_node_by_idx(itree, inode)
            do j = 1, size(node_cur%subset)
                local_k  = node_cur%subset(j)                              ! LOCAL index
                best_ref = block_tree%local_to_global_ref(itree, local_k)  ! GLOBAL ref id
                if (best_ref == 0) cycle
                call eulspace%get_ori(best_ref, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist, inplrotdist)
                if (dist < dist_min) then
                    dist_min = dist
                end if
            end do
        end if
    end subroutine srch_eul_bl_tree_stoch

    subroutine srch_eul_bl_tree_prob( otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        type(ori)     :: o, osym
        integer       :: inode, level, lidx, ridx, local_k, j
        real          :: dist_left, dist_right, inplrotdist, dist
        type(bt_node) :: node_L, node_R, node_cur
        inode = block_tree%get_root_node_idx(itree)
        do level = 1, MAX_LEVELS
            if (inode == 0) exit
            if (block_tree%is_leaf(itree, inode)) exit
            if (block_tree%get_subset_size(itree, inode) <= LEAF_BUDGET) exit
            call block_tree%get_children_idx(itree, inode, lidx, ridx)
            if (lidx == 0 .and. ridx == 0) exit
            dist_left  = huge(1.0)
            dist_right = huge(1.0)
            if (lidx /= 0) then
                node_L = block_tree%get_node_by_idx(itree, lidx)
                call eulspace%get_ori(node_L%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_left, inplrotdist)
            end if
            if (ridx /= 0) then
                node_R = block_tree%get_node_by_idx(itree, ridx)
                call eulspace%get_ori(node_R%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_right, inplrotdist)
            end if
            ! ---- probabilistic choice 
            if (lidx /= 0 .and. ridx /= 0) then
                if (sample_two(dist_left, dist_right) == 1) then
                    inode = lidx
                else
                    inode = ridx
                end if
            else if (lidx /= 0) then
                inode = lidx
            else
                inode = ridx
            end if
        end do
        ! ---- (final) brute-force scan the selected cluster ----
        dist_min = huge(1.0)
        best_ref = 0
        if (inode /= 0) then
            node_cur = block_tree%get_node_by_idx(itree, inode)
            do j = 1, size(node_cur%subset)
                local_k  = node_cur%subset(j)                              ! LOCAL index
                best_ref = block_tree%local_to_global_ref(itree, local_k)  ! GLOBAL ref id
                if (best_ref == 0) cycle
                call eulspace%get_ori(best_ref, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist, inplrotdist)
                if (dist < dist_min) dist_min = dist
            end do
        end if
        contains
            function sample_two(p1, p2) result(which)
                real, intent(in) :: p1, p2
                integer          :: which
                real             :: r, psum
                psum = p1 + p2
                if (psum <= 0.0) then
                    which = 1
                    return
                end if
                r = ran3()
                if (r < p1 / psum) then
                    which = 1
                else
                    which = 2
                end if
            end function sample_two
    end subroutine srch_eul_bl_tree_prob
end module simple_block_tree