!@descr: test for block-tree search space decomposition and mapping
program simple_test_block_tree
use simple_core_module_api
use simple_multi_dendro,  only: multi_dendro
use simple_block_tree     ! use all in there
use simple_timer
implicit none

character(len=*), parameter :: PGRP       = 'c1'
integer,          parameter :: NSPACE     = 5000
integer,          parameter :: NSPACE_SUB = 300
integer,          parameter :: NSAMPLE    = 1000

type(sym)          :: pgrpsym
type(multi_dendro) :: block_tree
type(oris)         :: eulspace, eulspace_sub
integer            :: i, j, ind_min, irnd
integer            :: itree, best_ref
type(ori)          :: osmp, o, osym
real               :: inplrotdist, dist, dist_min, dist_subspace, dists(NSAMPLE)
real               :: dist_min_tot, dist_subspace_tot, speedup, dist_min_exhaustive_tot
integer(timer_int_kind) :: total_time
real(timer_int_kind)    :: elapsed_tot = 0.0, exhaustive_time = 0.0

call pgrpsym%new(PGRP)
call eulspace    %new(NSPACE,     is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)
! this transformation is necessary to make the subspace a bona fide subset of the full space
call eulspace_sub%replace_with_closest(eulspace)

print *, 'Building block tree...'

block_tree = gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym)
total_time = tic()
dist_min_tot = 0.0
dist_subspace_tot = 0.0
do i = 1, NSAMPLE
    irnd = irnd_uni(NSPACE)
    call eulspace%get_ori(irnd, osmp)
    ! ---- (1) coarse: find closest sub-space point ----
    dist_min = huge(1.0)
    ind_min  = 1
    do j = 1, NSPACE_SUB
        call eulspace_sub%get_ori(j, o)
        call pgrpsym%sym_dists(osmp, o, osym, dist, inplrotdist)
        if (dist < dist_min) then
            dist_min = dist
            ind_min  = j
        end if
    end do
    dist_subspace = dist_min
    itree = ind_min
    ! call srch_eul_bl_tree_exhaustive(osmp, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
    !  call srch_eul_bl_tree(osmp, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min, l_greedy=.false.)
    call srch_eul_bl_tree_prob(osmp, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
    dists(i)          = dist_min
    dist_min_tot      = dist_min_tot + dist_min
    dist_subspace_tot = dist_subspace_tot + dist_subspace
    print *, 'SAMPLE ', irnd, ': itree=', itree, ' dist=', dist_min, ' dist_subspace=', dist_subspace
end do
elapsed_tot = toc(total_time)
total_time  = tic()
dist_min_exhaustive_tot = 0.0
do i = 1, NSAMPLE
    irnd = irnd_uni(NSPACE)
    call eulspace%get_ori(irnd, osmp)
    ! ---- (1) coarse: find closest sub-space point ----
    dist_min = huge(1.0)
    ind_min  = 1
    do j = 1, NSPACE
        call eulspace%get_ori(j, o)
        call pgrpsym%sym_dists(osmp, o, osym, dist, inplrotdist)
        if (dist < dist_min) then
            dist_min = dist
            ind_min  = j
        end if
    end do
    dist_min_exhaustive_tot = dist_min_exhaustive_tot + dist_min
end do
exhaustive_time = toc(total_time)
speedup = exhaustive_time / elapsed_tot
print *, 'AVERAGE DISTANCE TO CLOSEST SUB-SPACE POINT = ', dist_subspace_tot / NSAMPLE
print *, 'AVERAGE SEARCH TIME (ms)                    = ', elapsed_tot / NSAMPLE * 1e3
print *, 'AVERAGE EXHAUSTIVE SEARCH TIME (ms)         = ', exhaustive_time / NSAMPLE * 1e3
print *, 'SPEEDUP OVER EXHAUSTIVE SEARCH              = ', speedup
print *, 'AVERAGE EXHAUSTIVE DISTANCE                 = ', dist_min_exhaustive_tot / NSAMPLE
print *, '*****FINAL AVERAGE DISTANCE*****            = ', dist_min_tot / NSAMPLE
print *, '*****FINAL MEDIAN  DISTANCE*****            = ', median(dists)
end program simple_test_block_tree