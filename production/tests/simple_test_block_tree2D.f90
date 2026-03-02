program simple_test_block_tree2D
! Test to see if can find arbitrary reprojection in tree structure
use simple_core_module_api
use simple_commanders_api
use simple_multi_dendro,            only: multi_dendro
use simple_molecule_data,           only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
use simple_atoms,                   only: atoms
use simple_simple_volinterp,        only: reproject
use simple_parameters,              only: parameters
use simple_image,                   only: image                 
use simple_corrmat,                 only: calc_inpl_invariant_cc_nomirr
use simple_timer
use simple_block_tree
implicit none 
type(parameters)                :: p1
type(oris)                      :: eulspace, eulspace_sub
type(ori)                       :: osmp, o, osym
type(cmdline)                   :: cline
type(image)                     :: vol
character(len=:), allocatable   :: cmd  
real, allocatable               :: cc_mat(:,:)
real                            :: oa_minmax (2) = [0.,1.]
type(image), allocatable        :: proj_arr(:)
character(len=*), parameter     :: PGRP       = 'c1'
type(string)                    :: pdb_file, vol_file
type(molecule_data)             :: mol
type(atoms)                     :: molecule
type(parameters)                :: params
type(sym)                       :: pgrpsym
type(multi_dendro)              :: block_tree
integer(timer_int_kind)         :: build_time, search_time, full_cc, cc_exhaust
real(timer_int_kind)            :: build_time_rt, search_time_rt, full_cc_rt, cc_exhaust_rt
integer,          parameter     :: NSPACE     = 500
integer,          parameter     :: NSPACE_SUB = 30
integer,          parameter     :: NSAMPLE    = 1000
integer     :: i, j, best_ref, itree, ind_min, irnd
real        :: inplrotdist, dist, dist_min, dist_subspace, dists(NSAMPLE)
real        :: dist_min_tot, dist_subspace_tot, speedup, dist_min_exhaustive_tot

! SPIRAL REPROJECTIONS
params%smpd   = 3.
params%lp     = 3.
params%hp     = 20.
params%trs    = 10.
params%ctf    = 'no'
params%objfun = 'cc'
mol = sars_cov2_spkgp_6vxx()
pdb_file = '6VXX.pdb'
vol_file = '6VXX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, params%smpd, mol=mol)
call find_ldim_nptcls(vol_file, params%ldim, params%nptcls)
params%box = params%ldim(1)
params%mskdiam = params%smpd * params%ldim(1) / 2.
params%msk = params%box / 2.
call params%new(cline)
call vol%new(params%ldim, params%smpd)
call vol%read(vol_file)
call pgrpsym%new(PGRP)
call eulspace%new(NSPACE, is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)
call eulspace_sub%replace_with_closest(eulspace)
proj_arr   = reproject(vol, eulspace)
build_time = tic()
block_tree = gen_eulspace_block_tree2d(eulspace, eulspace_sub, pgrpsym, proj_arr, params)
! block_tree = gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym)
build_time_rt = toc(build_time)
print *, 'build_time', build_time_rt
! SEARCHING
dist_min_tot = 0.0
dist_subspace_tot = 0.0
search_time = tic()
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
    call srch_eul_bl_tree(osmp, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min, l_greedy=.false.)
    ! call srch_eul_bl_tree_prob(osmp, eulspace, pgrpsym, block_tree, itree, best_ref, dist_min)
    dists(i)          = dist_min
    dist_min_tot      = dist_min_tot + dist_min
    dist_subspace_tot = dist_subspace_tot + dist_subspace
    print *, 'SAMPLE ', irnd, ': itree=', itree, ' dist=', dist_min, ' dist_subspace=', dist_subspace
end do
print *, 'AVG DIST TO COARSE', dist_subspace_tot / NSAMPLE 
print *, 'AVG DIST TO FINE',   dist_min_tot / NSAMPLE
search_time_rt = toc(search_time)
print *, 'AVG SEARCH TIME', search_time_rt / NSAMPLE
! exahustive search over full corrmat 
full_cc = tic()
cc_mat = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, proj_arr)
full_cc_rt = toc(full_cc)
print *, 'FULL CORRMAT CALCULATION TIME', full_cc_rt
! cc_exhaust = tic()
! do i = 1, NSAMPLE
!     irnd = irnd_uni(NSPACE)
!     call eulspace%get_ori(irnd, osmp)
!     ! ---- (1) coarse: find closest sub-space point ----
!     dist_min = huge(1.0)
!     ind_min  = 1
!     do j = 1, NSPACE
!         call eulspace%get_ori(j, o)
!         call pgrpsym%sym_dists(osmp, o, osym, dist, inplrotdist)
!         if (dist < dist_min) then
!             dist_min = dist
!             ind_min  = j
!         end if
!     end do
!     dist_min_exhaustive_tot = dist_min_exhaustive_tot + dist_min
! end do
! cc_exhaust_rt = toc(cc_exhaust)
! print *, ' >>> AVG TREE SEARCH TIME:', search_time_rt / NSAMPLE 
! print *, '>>> AVG EXAHUSTIVE TIME:', cc_exhaust_rt / NSAMPLE 
! print *, '>>> AVG ERROR:'
! exhausive search over angular distances 
end program simple_test_block_tree2D
    