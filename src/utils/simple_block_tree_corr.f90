module simple_block_tree_corr
use simple_core_module_api
use simple_multi_dendro, only: multi_dendro
use simple_image,        only: image
use simple_parameters,   only: parameters
implicit none

public :: gen_eulspace_block_tree_corr
public :: gen_corr_block_tree_aff_prop
public :: gen_corr_block_tree_hac
private
#include "simple_local_flags.inc"

logical, parameter :: DEBUG = .false.

contains

    function gen_corr_block_tree_hac(refimgs, params, ntrees_target) result(block_tree)
        use simple_corrmat, only: calc_inpl_invariant_cc_nomirr
        use simple_hclust,  only: hclust
        class(image),      intent(inout) :: refimgs(:)
        class(parameters), intent(in)    :: params
        integer,           intent(in)    :: ntrees_target
        type(multi_dendro)   :: block_tree
        real,    allocatable :: distmat(:,:), sub_distmat(:,:)
        real,    allocatable :: height(:)
        integer, allocatable :: labels(:), refs(:)
        integer, allocatable :: merge_mat(:,:)
        type(hclust)         :: hc
        integer              :: nspace, ntrees, itree, i, j, nrefs
        nspace = size(refimgs)
        allocate(distmat(nspace, nspace), source = 0.0)
        ! Calculate full distance matrix
        distmat = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, refimgs)
        distmat = 1.0 - distmat
        call normalize_minmax(distmat)
        ntrees = min(max(1, ntrees_target), nspace)
        allocate(labels(nspace), source=1)
        if (nspace > 1) then
            allocate(merge_mat(2, nspace-1), height(nspace-1))
            ! Full HAC then cut to requested number of splits.
            call hc%new(nspace, distmat, LINK_AVERAGE)
            call hc%cluster(merge_mat, height, labels, ntrees)
            call hc%kill()
        else
            ntrees = 1
        end if
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        if(DEBUG) write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        ! Build trees
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,i,j) schedule(static)
        do itree = 1, ntrees
            refs = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            if(nrefs == 0) then
                if(DEBUG) write(*,'(a,1x,i0)') 'TREE ', itree, ': EMPTY, SKIPPING...'
                cycle
            end if
            allocate(sub_distmat(nrefs, nrefs), source = 0.0)
            do i = 1, nrefs
                do j = 1, nrefs
                    if(i /= j) sub_distmat(i,j) = distmat(refs(i), refs(j))
                end do
            end do
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            if(DEBUG) write(*,'(a,1x,i0,a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            deallocate(refs, sub_distmat)
        end do
        !$omp end parallel do
        if(allocated(merge_mat)) deallocate(merge_mat)
        if(allocated(height))    deallocate(height)
        deallocate(distmat, labels)
        if(DEBUG) print *, 'Finished building even-split trees.'
    end function gen_corr_block_tree_hac

     function gen_eulspace_block_tree_corr(eulspace, eulspace_sub, pgrpsym, refimgs, params) result(block_tree)
        use simple_corrmat, only: calc_inpl_invariant_cc_nomirr
        use simple_eulspace_neigh_map, only: eulspace_neigh_map
        class(oris),       intent(in)    :: eulspace, eulspace_sub
        class(image),      intent(in)    :: refimgs(:)
        class(parameters), intent(in)    :: params
        class(sym),        intent(inout) :: pgrpsym
        type(eulspace_neigh_map)   :: neigh_map
        type(multi_dendro)         :: block_tree
        type(image), allocatable   :: sub_imgs(:)     
        integer, allocatable       :: labels(:), refs(:)
        real,    allocatable       :: sub_distmat(:,:)
        integer                    :: ntrees, itree, nrefs
        call neigh_map%new(eulspace, eulspace_sub, pgrpsym)
        labels = neigh_map%get_full2sub_map()
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        if(DEBUG) write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            if( nrefs == 0 ) then
                if(DEBUG) write(*,'(a,1x,i0)') 'TREE ', itree, ': EMPTY, SKIPPING...'
                cycle
            end if
            if(DEBUG) write(*,'(a,1x,i0,a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            allocate(sub_distmat(nrefs, nrefs), source=0.0)  
            allocate(sub_imgs(nrefs))
            ! getting images corresponding to each tree
            sub_imgs = refimgs(refs)
            sub_distmat = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, sub_imgs)
            sub_distmat = 1 - sub_distmat
            call normalize_minmax(sub_distmat)
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(sub_distmat, refs, sub_imgs)
        end do
        call neigh_map%kill
        if (allocated(labels)) deallocate(labels)
        if(DEBUG) print *, 'Finished building block tree.'
    end function gen_eulspace_block_tree_corr

    function gen_corr_block_tree_aff_prop(refimgs, params) result(block_tree)
        use simple_aff_prop, only: aff_prop 
        use simple_corrmat,  only: calc_inpl_invariant_cc_nomirr
        class(image),      intent(inout) :: refimgs(:)
        class(parameters), intent(in)    :: params
        type(multi_dendro)   :: block_tree
        type(aff_prop)       :: affprop
        real,    allocatable :: distmat(:,:), sub_distmat(:,:)
        integer, allocatable :: labels(:), centers(:), refs(:)
        integer              :: nspace, nspace_sub, ntrees, itree, nrefs, i, j 
        real                 :: simsum
        ! affinty propagation exemplars are coarse, members are fine 
        nspace  = size(refimgs)
        allocate(distmat(nspace, nspace))
        distmat = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, refimgs)
        call normalize_minmax(distmat)
        ! normalize 
        call affprop%new(nspace, distmat)
        call affprop%propagate(centers, labels, simsum)
        nspace_sub  = size(centers) 
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,i,j) schedule(static)
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = size(refs)
            if( nrefs == 0 ) then
                if(DEBUG) write(*,'(a,1x,i0)') 'TREE ', itree, ': EMPTY, SKIPPING...'
                cycle
            end if
            allocate(sub_distmat(nrefs,nrefs), source = 0.)
            do i = 1, nrefs 
                do j = 1, nrefs 
                    if(i /= j) sub_distmat(i,j) = distmat(refs(i), refs(j))
                end do 
            end do 
            sub_distmat = 1 - sub_distmat
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(refs,sub_distmat)
        end do
        !$omp end parallel do
    end function gen_corr_block_tree_aff_prop

end module simple_block_tree_corr
