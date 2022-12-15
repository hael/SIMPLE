program simple_test_neigh
include 'simple_lib.f08'
implicit none
character(len=*), parameter :: PGRP       = 'c2'
integer,          parameter :: NSPACE     = 20000
integer,          parameter :: NSPACE_SUB = 500
logical            :: lnns(NSPACE), lclosest(NSPACE)
real               :: nnn(NSPACE_SUB), angres, tot_nnn, euldists(NSPACE_SUB), inplrotdist
integer            :: i, closest, subspace_inds(NSPACE_SUB)
type(sym)          :: pgrpsym 
type(ori)          :: o, o_close, osym
type(stats_struct) :: nnn_stats, euldist_stats
type(oris)         :: eulspace, eulspace_sub
call pgrpsym%new(PGRP)
call eulspace    %new(NSPACE,     is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)
lclosest = .false.
write(*,'(a)') '>>> OBTAINING STATISTICS ABOUT THE SUBSPACE'
do i = 1,NSPACE_SUB
    call progress_gfortran(i, NSPACE_SUB)
    call eulspace_sub%get_ori(i, o)
    closest           = pgrpsym%find_closest_proj(eulspace, o)
    call eulspace%get_ori(closest, o_close)
    call pgrpsym%sym_dists(o, o_close, osym, euldists(i), inplrotdist)
    lclosest(closest) = .true.
    subspace_inds(i)  = closest
end do
call calc_stats(euldists, euldist_stats)
write(*,'(a,1x,f15.6)') '# UNIQUE PROJECTIONS IDENTIFIED :', real(count(lclosest))
write(*,'(a,1x,f15.6)') '% UNIQUE PROJECTIONS IDENTIFIED :', (real(count(lclosest)) / real(NSPACE_SUB)) * 100.
write(*,'(a)') '>>> DISTANCE STATISTICS FOR CORRESPONDING PROJECTION DIRECTIONS'
write(*,'(a,1x,f15.6)') 'AVERAGE EULDIST (IN DEGREES) :', euldist_stats%avg
write(*,'(a,1x,f15.6)') 'SDEV    EULDIST (IN DEGREES) :', euldist_stats%sdev
write(*,'(a,1x,f15.6)') 'MINIMUM EULDIST (IN DEGREES) :', euldist_stats%minv
write(*,'(a,1x,f15.6)') 'MAXIMUM EULDIST (IN DEGREES) :', euldist_stats%maxv
write(*,'(a)') '>>> OBTAINING STATISTICS ABOUT THE NEAREST NEIGHBORS'
angres = eulspace_sub%find_angres()
write(*,'(a,1x,f15.6)') 'ANGULAR RESOLUTION OF SUBSPACE (IN DEGREES) :', angres
tot_nnn = 0.
do i = 1,NSPACE_SUB
    call progress_gfortran(i, NSPACE_SUB)
    call eulspace_sub%get_ori(i, o)
    call pgrpsym%nearest_proj_neighbors(eulspace, o, angres, lnns)
    nnn(i) = real(count(lnns))
    tot_nnn = tot_nnn + nnn(i)
end do
call calc_stats(nnn, nnn_stats)
write(*,'(a,1x,f15.6)') 'TOTAL   # NN :', tot_nnn
write(*,'(a,1x,f15.6)') 'AVERAGE # NN :', nnn_stats%avg
write(*,'(a,1x,f15.6)') 'SDEV    # NN :', nnn_stats%sdev
write(*,'(a,1x,f15.6)') 'MINIMUM # NN :', nnn_stats%minv
write(*,'(a,1x,f15.6)') 'MAXIMUM # NN :', nnn_stats%maxv
end program simple_test_neigh
