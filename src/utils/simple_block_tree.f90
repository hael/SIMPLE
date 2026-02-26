module simple_block_tree
use simple_core_module_api
use simple_srchspace_map,    only: srchspace_map
use simple_multi_dendro,     only: multi_dendro
use simple_binary_tree,      only: bt_node
use simple_strategy2D_utils, only: calc_cluster_cavgs_dmat
use simple_image,            only: image
use simple_parameters,       only: parameters
implicit none

public :: gen_eulspace_block_tree
public :: gen_eulspace_block_tree2D
public :: srch_eul_bl_tree_exhaustive
public :: srch_eul_bl_tree
public :: srch_eul_bl_tree_prob
private
#include "simple_local_flags.inc"

integer, parameter :: LINK_SINGLE   = 1
integer, parameter :: LINK_COMPLETE = 2
integer, parameter :: LINK_AVERAGE  = 3
logical, parameter :: DEBUG = .true.

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
        !$omp parallel do default(shared) proc_bind(close) private(i,j,o,o_sub,osym,dtmp,inplrotdist) schedule(static)
        do i = 1, nspace_sub
            call eulspace_sub%get_ori(i, o_sub)
            do j = 1, nspace
                call eulspace%get_ori(j, o)
                call pgrpsym%sym_dists(o_sub, o, osym, dtmp, inplrotdist)
                distmat(i,j) = dtmp
            end do
        end do
        !$omp end parallel do
        call normalize_minmax(distmat)
        call mapper%new(nspace, nspace_sub, distmat)
        labels = mapper%get_full2sub_map()
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,i,j,oi,oj,osym,inplrotdist,dtmp) schedule(static)
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            if( nrefs == 0 ) then
                write(*,'(a,1x,i0)') 'TREE ', itree, ': EMPTY, SKIPPING...'
                cycle
            end if
            write(*,'(a,1x,i0,a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            allocate(sub_distmat(nrefs, nrefs), source=0.0)  
            do i = 1, nrefs - 1
                call eulspace%get_ori(refs(i), oi)        
                do j = i + 1, nrefs
                    call eulspace%get_ori(refs(j), oj)
                    call pgrpsym%sym_dists(oi, oj, osym, dtmp, inplrotdist)
                    sub_distmat(i,j) = dtmp
                    sub_distmat(j,i) = dtmp
                end do
            end do
            call normalize_minmax(sub_distmat)
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(sub_distmat)
            deallocate(refs)
        end do
        !$omp end parallel do
        deallocate(distmat)
        if (allocated(labels)) deallocate(labels)
        if( DEBUG) print *, 'Finished building block tree.'
    end function gen_eulspace_block_tree

    function gen_eulspace_block_tree2D(eulspace, eulspace_sub, pgrpsym, refimgs, params) result(block_tree)
        class(oris), intent(in)        :: eulspace, eulspace_sub
        class(image), intent(in)       :: refimgs(:)
        class(parameters), intent(in)  :: params
        class(sym),  intent(inout)     :: pgrpsym
        type(ori)                  :: o, o_sub, oi, oj, osym
        type(srchspace_map)        :: mapper
        type(multi_dendro)         :: block_tree
        type(image), allocatable   :: sub_imgs(:)     
        integer, allocatable       :: labels(:), refs(:)
        real,    allocatable       :: distmat(:,:), sub_distmat(:,:)
        integer                    :: i, j, ntrees, itree, nrefs, nspace, nspace_sub
        real                       :: inplrotdist, dtmp, oa_minmax(2)
        nspace     = eulspace%get_noris()
        nspace_sub = eulspace_sub%get_noris()
        allocate(distmat(nspace_sub, nspace))
        !$omp parallel do default(shared) proc_bind(close) private(i,j,o,o_sub,osym,dtmp,inplrotdist) schedule(static)
        do i = 1, nspace_sub
            call eulspace_sub%get_ori(i, o_sub)
            do j = 1, nspace
                call eulspace%get_ori(j, o)
                call pgrpsym%sym_dists(o_sub, o, osym, dtmp, inplrotdist)
                distmat(i,j) = dtmp
            end do
        end do
        !$omp end parallel do
        call normalize_minmax(distmat)
        call mapper%new(nspace, nspace_sub, distmat)
        labels = mapper%get_full2sub_map()
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,sub_imgs,i,j,oi,oj,osym,inplrotdist,dtmp) schedule(static)
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            if( nrefs == 0 ) then
                write(*,'(a,1x,i0)') 'TREE ', itree, ': EMPTY, SKIPPING...'
                cycle
            end if
            write(*,'(a,1x,i0,a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            allocate(sub_distmat(nrefs, nrefs), source=0.0)  
            allocate(sub_imgs(nrefs))
            ! need images corresponding to refs in sub tree
            sub_imgs = refimgs(refs)
            oa_minmax  = [0.,1.]
            sub_distmat = calc_cluster_cavgs_dmat(params, sub_imgs, oa_minmax, 'cc')
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(sub_distmat, refs, sub_imgs)
        end do
        !$omp end parallel do
        deallocate(distmat)
        if (allocated(labels)) deallocate(labels)
        if( DEBUG) print *, 'Finished building block tree.'
    end function gen_eulspace_block_tree2D

    ! Exhaustive search over *leaf nodes only* reachable from the root.
    ! Intended for testing/validation against greedy descent.
    subroutine srch_eul_bl_tree_exhaustive(otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        type(ori)     :: o, osym
        real          :: inplrotdist, dist
        type(bt_node) :: node_root, node_cur
        integer, allocatable :: stack(:)
        integer :: top, inode
        integer :: ref, nfun, pop
        ! Initialize outputs
        best_ref = 0
        dist_min = huge(1.0)
        node_root = block_tree%get_root_node(itree)
        if (node_root%node_idx == 0) then
            THROW_HARD('srch_eul_bl_tree_exhaustive_leaves: empty tree / invalid itree')
        end if
        allocate(stack(64)) ! whatever, chat wanted 64
        top = 1
        stack(top) = node_root%node_idx
        nfun = 0
        pop = block_tree%get_tree_pop(itree)
        do while (top > 0)
            inode = stack(top)
            top = top - 1
            if (inode == 0) cycle
            node_cur = block_tree%get_node(itree, inode)
            ! Only evaluate leaves
            if (block_tree%is_leaf(itree, inode)) then
                ref = node_cur%ref_idx
                if (ref == 0) then
                    THROW_HARD('srch_eul_bl_tree_exhaustive_leaves: leaf has ref_idx=0')
                end if
                call eulspace%get_ori(ref, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist, inplrotdist)
                nfun = nfun + 1
                if (dist < dist_min) then
                    dist_min = dist
                    best_ref = ref
                end if
            else
                ! Push children
                if (node_cur%left_idx /= 0)  call push_stack(stack, top, node_cur%left_idx)
                if (node_cur%right_idx /= 0) call push_stack(stack, top, node_cur%right_idx)
            end if
        end do
        if (allocated(stack)) deallocate(stack)

        contains

            subroutine push_stack(stk, top, val)
                integer, allocatable, intent(inout) :: stk(:)
                integer,              intent(inout) :: top
                integer,              intent(in)    :: val
                integer, allocatable :: tmp(:)
                integer :: n
                if (val == 0) return
                n = size(stk)
                if (top + 1 > n) then
                    allocate(tmp(2*n))
                    tmp(1:n) = stk
                    tmp(n+1:) = 0
                    call move_alloc(tmp, stk)
                end if
                top = top + 1
                stk(top) = val
            end subroutine push_stack

    end subroutine srch_eul_bl_tree_exhaustive

    subroutine srch_eul_bl_tree( otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min, l_greedy)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        logical,            intent(in)    :: l_greedy
        type(ori)     :: o, osym
        integer       :: inode, level, local_k, j, max_levs, endit
        real          :: dist_left, dist_right, inplrotdist, dist, dist_subspace
        type(bt_node) :: node_L, node_R, node_cur, node_root
        dist_subspace = dist_min ! true on input
        node_root = block_tree%get_root_node(itree)
        if( node_root%ref_idx == 0 ) THROW_HARD('root node has no ref_idx')
        call eulspace%get_ori(node_root%ref_idx, o)
        call pgrpsym%sym_dists(otrial, o, osym, dist_min, inplrotdist)
        inode    = node_root%node_idx
        max_levs = block_tree%get_tree_height(itree)
        endit    = max_levs
        if( max_levs > 2 .and. .not. l_greedy )then
            endit = irnd_uni(max_levs)
            endit = max(2, endit) ! ensure at least 2 levels of descent
        end if
        level = 0
        do
            level = level + 1
            if (inode == 0) exit
            if (block_tree%is_leaf(itree, inode)) exit
            node_cur = block_tree%get_node(itree, inode)
            if (node_cur%left_idx == 0 .and. node_cur%right_idx == 0) exit
            dist_left  = huge(1.0)
            dist_right = huge(1.0)
            if (node_cur%left_idx /= 0) then
                node_L = block_tree%get_node(itree, node_cur%left_idx)
                call eulspace%get_ori(node_L%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_left, inplrotdist)
            end if
            if (node_cur%right_idx /= 0) then
                node_R = block_tree%get_node(itree, node_cur%right_idx)
                call eulspace%get_ori(node_R%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_right, inplrotdist)
            end if
            ! travel down the closer child
            if (dist_left <= dist_right) then
                inode = node_cur%left_idx
            else
                inode = node_cur%right_idx
            end if
            ! best-first search
            if( any([dist_left, dist_right] < dist_min) ) then
                dist_min = min(dist_left, dist_right)
                best_ref = merge(node_L%ref_idx, node_R%ref_idx, dist_left <= dist_right)
            end if
            if (.not. l_greedy .and. level > endit) exit
        end do
         if( dist_min >= dist_subspace)then
            ! fallback to exhaustive descent if we ended up in a bad leaf
            call srch_eul_bl_tree_exhaustive(otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        end if
    end subroutine srch_eul_bl_tree

    subroutine srch_eul_bl_tree_prob( otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        class(ori),         intent(in)    :: otrial
        class(oris),        intent(in)    :: eulspace
        class(sym),         intent(inout) :: pgrpsym
        type(multi_dendro), intent(in)    :: block_tree
        integer,            intent(in)    :: itree
        integer,            intent(out)   :: best_ref
        real,               intent(out)   :: dist_min
        type(ori)     :: o, osym
        integer       :: inode, level, local_k, j
        real          :: dist_left, dist_right, inplrotdist, dist, p_left, p_right, dist_subspace
        type(bt_node) :: node_L, node_R, node_cur, node_root
        dist_subspace = dist_min ! true on input
        node_root = block_tree%get_root_node(itree)
        if( node_root%ref_idx == 0 ) THROW_HARD('root node has no ref_idx')
        call eulspace%get_ori(node_root%ref_idx, o)
        call pgrpsym%sym_dists(otrial, o, osym, dist_min, inplrotdist)
        inode = node_root%node_idx
        do
            if (inode == 0) exit
            if (block_tree%is_leaf(itree, inode)) exit
            node_cur = block_tree%get_node(itree, inode)
            if( node_cur%left_idx == 0 .and. node_cur%right_idx == 0 ) exit
            dist_left  = huge(1.0)
            dist_right = huge(1.0)
            if (node_cur%left_idx /= 0) then
                node_L = block_tree%get_node(itree, node_cur%left_idx)
                call eulspace%get_ori(node_L%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_left, inplrotdist)
            end if
            if (node_cur%right_idx /= 0) then
                node_R = block_tree%get_node(itree, node_cur%right_idx)
                call eulspace%get_ori(node_R%ref_idx, o)
                call pgrpsym%sym_dists(otrial, o, osym, dist_right, inplrotdist)
            end if
            p_left  = exp(-dist_left)  ! assumes normalized distances [0,1)
            p_right = exp(-dist_right) ! assumes normalized distances [0,1)
            ! travel down the most likely child
            if( sample_two(p_left, p_right) == 1 )then
                inode = node_cur%left_idx
            else
                inode = node_cur%right_idx
            end if
            ! best-first search
            if( any([dist_left, dist_right] < dist_min) ) then
                dist_min = min(dist_left, dist_right)
                best_ref = merge(node_L%ref_idx, node_R%ref_idx, dist_left <= dist_right)
            end if
        end do
        if( dist_min >= dist_subspace)then
            ! fallback to exhaustive descent if we ended up in a bad leaf due to randomness
            call srch_eul_bl_tree_exhaustive(otrial, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
        end if
            
        contains

            function sample_two(p1, p2) result(which)
                real, intent(in) :: p1, p2
                integer          :: which
                real             :: r, psum
                psum = p1 + p2
                if (psum <= 0.0) then
                    ! policy: choose uniformly between 1 and 2 if both zero
                    which = merge(1,2,ran3() < 0.5)
                    return
                end if
                r = ran3() ! assumed in [0,1)
                which = merge(1,2,r < p1/psum)
            end function sample_two

    end subroutine srch_eul_bl_tree_prob

end module simple_block_tree