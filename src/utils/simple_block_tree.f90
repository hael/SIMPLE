module simple_block_tree
use simple_core_module_api
use simple_srchspace_map, only: srchspace_map
use simple_multi_dendro,  only: multi_dendro
implicit none

public :: gen_eulspace_block_tree
private

integer, parameter :: LINK_SINGLE   = 1
integer, parameter :: LINK_COMPLETE = 2
integer, parameter :: LINK_AVERAGE  = 3

contains

    function gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym ) result(block_tree)
        class(oris), intent(in)    :: eulspace, eulspace_sub
        class(sym),  intent(inout) :: pgrpsym
        type(ori)                  :: o, o_sub, oi, oj, osym
        type(srchspace_map)        :: mapper
        type(multi_dendro)         :: block_tree
        integer, allocatable       :: labels(:), refs(:)
        real,    allocatable       :: distmat(:,:), sub_distmat(:,:)
        integer                    :: i, j, ntrees, itree, nrefs, nspace, nspace_sub
        real                       :: inplrotdist
        nspace     = eulspace%get_noris()
        nspace_sub = eulspace_sub%get_noris()
        allocate(distmat(nspace_sub, nspace))
        !$omp parallel do default(shared) proc_bind(close) private(i,j,o,o_sub,osym,inplrotdist) schedule(static)
        do i = 1,nspace_sub
            call eulspace_sub%get_ori(i, o_sub)
            do j = 1,nspace
                call eulspace%get_ori(j, o)
                distmat(i,j) = o_sub.euldist.o
                call pgrpsym%sym_dists(o_sub, o, osym, distmat(i,j), inplrotdist)
            end do
        end do
        !$omp end parallel do 
        call mapper%new(nspace, nspace_sub, distmat)
        labels = mapper%get_full2sub_map()
        call block_tree%new(labels)
        ntrees = block_tree%get_n_trees()
        write(*,'(a,1x,i0)') 'NUMBER OF TREES :', ntrees
        !$omp parallel do default(shared) proc_bind(close) private(itree,refs,nrefs,sub_distmat,i,j,oi,oj,osym,inplrotdist) schedule(static)
        do itree = 1, ntrees
            refs  = block_tree%get_tree_refs(itree)
            nrefs = block_tree%get_tree_pop(itree)
            write(*,'(a,1x,i0)') 'TREE ', itree, ': NUMBER OF REFS :', nrefs
            allocate(sub_distmat(nrefs,nrefs), source=0.0)
            do i = 1, nrefs - 1
                call eulspace%get_ori(i, oi)
                do j = i + 1, nrefs
                    call eulspace%get_ori(j, oj)
                    call pgrpsym%sym_dists(oi, oj, osym, sub_distmat(i,j), inplrotdist)
                    sub_distmat(j,i) = sub_distmat(i,j)
                enddo
            enddo
            call block_tree%build_tree_from_subdistmat(itree, refs, sub_distmat, LINK_AVERAGE)
            deallocate(sub_distmat)
        end do
        !$omp end parallel do
        deallocate(distmat)
    end function gen_eulspace_block_tree

end module simple_block_tree
