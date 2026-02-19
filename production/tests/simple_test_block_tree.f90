!@descr: test for block-tree search space decomposition and mapping
program simple_test_block_tree
use simple_core_module_api
use simple_multi_dendro,  only: multi_dendro
use simple_block_tree,    only: gen_eulspace_block_tree
implicit none

character(len=*), parameter :: PGRP       = 'c1'
integer,          parameter :: NSPACE     = 5000
integer,          parameter :: NSPACE_SUB = 500
integer,          parameter :: NSAMPLE    = 1000

type(sym)            :: pgrpsym 
type(multi_dendro)   :: block_tree
type(oris)           :: eulspace, eulspace_sub

integer   :: i, j, ind_min, irnd, imed, left_ref_idx, right_ref_idx
type(ori) :: osmp, o, osym
real      :: inplrotdist, dist, dist_min, dist_left, dist_right, dist_subspace

call pgrpsym%new(PGRP)
call eulspace    %new(NSPACE,     is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)

block_tree = gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym)

do i = 1, NSAMPLE
    irnd = irnd_uni(NSPACE)
    call eulspace%get_ori(irnd, osmp)
    ! find closest sub-space point
    dist_min = huge(1.0)
    do j = 1, NSPACE_SUB
        call eulspace_sub%get_ori(j, o)
        call pgrpsym%sym_dists(osmp, o, osym, dist, inplrotdist)

        ! print *, 'SAMPLE ', irnd, ': sub-space point ', j, ' dist = ', dist

        if( dist < dist_min) then
            dist_min = dist
            ind_min  = j
        end if
    end do
    dist_subspace = dist_min

    ! print *, 'to here: closest sub-space point = ', ind_min, ' dist = ', dist_min
    imed = block_tree%get_medoid(ind_min)

    ! print *, 'SAMPLE ', irnd, ': closest sub-space point = ', ind_min, ' medoid = ', imed, ' dist = ', dist_min

    call eulspace%get_ori(imed, o)
    call pgrpsym%sym_dists(osmp, o, osym, dist_min, inplrotdist)
    do
        call block_tree%get_left_right_idxs(imed, left_ref_idx, right_ref_idx)
        dist_left  = huge(1.0)
        if( left_ref_idx /= 0 )then
            call eulspace%get_ori(left_ref_idx, o)
            call pgrpsym%sym_dists(osmp, o, osym, dist_left, inplrotdist)
        endif
        dist_right = huge(1.0)
        if( right_ref_idx /= 0 )then
            call eulspace%get_ori(right_ref_idx, o)
            call pgrpsym%sym_dists(osmp, o, osym, dist_right, inplrotdist)
        endif
        if( any([dist_left, dist_right] < dist_min) ) then
             if( dist_left < dist_right) then
                dist_min = dist_left
                imed     = left_ref_idx
             else
                dist_min = dist_right
                imed     = right_ref_idx
             end if
        else
            exit
        end if
    end do
    print *, 'SAMPLE ', irnd, ': closest sub-space point = ', ind_min, ' medoid = ', imed, ' dist = ', dist_min, ' dist_subspace = ', dist_subspace
end do


end program simple_test_block_tree
