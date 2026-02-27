program simple_test_tree_srch
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
implicit none 
type(parameters)                :: p1
type(oris)                      :: eulspace, eulspace_sub
type(cmdline)                   :: cline
integer                         ::  NSPACE = 5000, NSPACE_SUB = 300
type(image)                     :: vol
character(len=:), allocatable   :: cmd  
real, allocatable               :: dist_mat_cc(:,:)
type(image), allocatable        :: proj_arr(:)
character(len=*), parameter     :: PGRP       = 'c1'
type(string)                    :: pdb_file, vol_file
type(molecule_data)             :: mol
type(atoms)                     :: molecule
type(parameters)                :: params
type(sym)                       :: pgrpsym
type(multi_dendro)              :: block_tree
! SPIRAL REPROJECTIONS
call params%new(cline)
params%smpd = 3.
params%hp   = 20.
params%lp   = 3.
params%trs  = 5.
params%ctf  = 'no'
mol = sars_cov2_spkgp_6vxx()
pdb_file = '6VXX.pdb'
vol_file = '6VXX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, params%smpd, mol=mol)
call find_ldim_nptcls(vol_file, params%ldim, params%nptcls)
params%box = params%ldim(1)
params%mskdiam = params%smpd * params%ldim(1) / 2.
! params%msk = params%box / 2.
params%msk_crop = round2even(params%mskdiam / params%smpd / 2.)
print *, params%box, params%msk_crop, cline%defined('pftsz')
! need to check how msk_crop is set
call vol%new(params%ldim, params%smpd)
call vol%read(vol_file)
call pgrpsym%new(PGRP)
call eulspace%new(NSPACE, is_ptcl=.false.)
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)
call eulspace_sub%replace_with_closest(eulspace)
proj_arr   = reproject(vol, eulspace, NSPACE)
! block_tree = gen_eulspace_block_tree2d(eulspace, eulspace_sub, pgrpsym, proj_arr, params)
! block_tree = gen_eulspace_block_tree(eulspace, eulspace_sub, pgrpsym)
! allocate(dist_mat_fm(dist_mat_cc(NSPACE,NSPACE))
! dist_mat_cc = calc_inpl_invariant_cc_nomirr(params, params%hp, params%lp, params%trs, proj_arr)

end program simple_test_tree_srch
    