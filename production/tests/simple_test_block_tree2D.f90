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
use simple_block_tree
use simple_corrmat,                 only: calc_inpl_invariant_cc_nomirr
use simple_strategy2D_utils,        only: calc_cluster_cavgs_dmat
implicit none 
type(parameters)                :: p1
type(oris)                      :: eulspace, eulspace_sub
type(cmdline)                   :: cline
integer                         ::  NSPACE = 50, NSPACE_SUB = 10
type(image)                     :: vol
character(len=:), allocatable   :: cmd  
real, allocatable               :: dist_mat_cc(:,:)
real                            :: oa_minmax (2) = [0.,1.]
type(image), allocatable        :: proj_arr(:)
character(len=*), parameter     :: PGRP       = 'c1'
type(string)                    :: pdb_file, vol_file
type(molecule_data)             :: mol
type(atoms)                     :: molecule
type(parameters)                :: params
type(sym)                       :: pgrpsym
type(multi_dendro)              :: block_tree
! SPIRAL REPROJECTIONS
params%smpd   = 3.
params%lp     = 3.
params%trs    = 5.
params%ctf    = 'no'
params%objfun = 'cc'
mol = sars_cov2_spkgp_6vxx()
pdb_file = '6VXX.pdb'
vol_file = '6VXX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, params%smpd, mol=mol)
call find_ldim_nptcls(vol_file, params%ldim, params%nptcls)
params%box = params%ldim(1)
params%mskdiam = params%smpd * params%ldim(1) / 2.
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
! block_tree = gen_eulspace_block_tree2d(eulspace, eulspace_sub, pgrpsym, proj_arr, params)
! block_tree = gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym)
! allocate(dist_mat_fm(dist_mat_cc(NSPACE,NSPACE))
! dist_mat_cc = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, proj_arr)
dist_mat_cc = calc_cluster_cavgs_dmat(params, proj_arr, oa_minmax, 'cc')

end program simple_test_block_tree2D
    