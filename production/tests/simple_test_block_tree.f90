!@descr: test for refine=neigh modes
program simple_test_block_tree
use simple_core_module_api
use simple_srchspace_map, only: srchspace_map
use simple_multi_dendro,  only: multi_dendro
implicit none

character(len=*), parameter :: PGRP       = 'c1'
integer,          parameter :: NSPACE     = 5000
integer,          parameter :: NSPACE_SUB = 300

logical              :: lclosest(NSPACE)
integer              :: i, j
type(ori)            :: o, o_sub, o_close
type(srchspace_map)  :: mapper
type(multi_dendro)   :: block_tree
integer, allocatable :: labels(:)
real, allocatable    :: distmat(:,:)
type(oris)           :: eulspace, eulspace_sub
call eulspace    %new(NSPACE,     is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call eulspace%spiral
call eulspace_sub%spiral
allocate(distmat(NSPACE_SUB, NSPACE))
do i = 1,NSPACE_SUB
    call eulspace_sub%get_ori(i, o_sub)
    do j = 1,NSPACE
        call eulspace%get_ori(j, o)
        distmat(i,j) = o_sub.euldist.o
     end do
end do
call mapper%new(NSPACE, NSPACE_SUB, distmat)




end program simple_test_block_tree
