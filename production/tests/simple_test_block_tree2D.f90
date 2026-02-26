program simple_test_block_tree2D
use simple_core_module_api
use simple_commanders_api
use simple_multi_dendro,            only: multi_dendro
use simple_molecule_data,           only: molecule_data, betagal_1jyx, sars_cov2_spkgp_6vxx
use simple_atoms,                   only: atoms
use simple_simple_volinterp,        only: reproject
use simple_parameters,              only: parameters
use simple_image,                   
use simple_block_tree
implicit none
type(string)                        :: pdb_file, vol_file
type(molecule_data)                 :: mol
type(atoms)                         :: molecule
type(oris)                          :: eulspace, eulspace_sub
type(image)                         :: vol 
type(image), allocatable            :: refimgs(:), refimgs_sub(:)
type(sym)                           :: pgrpsym
type(multi_dendro)                  :: block_tree
type(parameters)                    :: params
real, parameter                     :: smpd = 5.0
integer                             :: nspace = 5000, nspace_sub = 300
character(len=*), parameter         :: PGRP       = 'c1'
! populate relevant params for pftcc calculation
params%ctf = 'no'
params%trs = 5.0
params%hp  = 20. 
params%lp  = 3. 
params%smpd = smpd
! spiral reprojections of test volume
mol = sars_cov2_spkgp_6vxx()
pdb_file = '6VXX.pdb'
vol_file = '6VXX.mrc'
call molecule%pdb2mrc( pdb_file, vol_file, params%smpd, mol=mol )
call find_ldim_nptcls(vol_file, params%ldim, params%nptcls)
call vol%new(params%ldim, params%smpd)
call vol%read(vol_file)
call pgrpsym%new(PGRP)
call eulspace    %new(NSPACE,     is_ptcl=.false.)
call eulspace%spiral()
call eulspace_sub%new(NSPACE_SUB, is_ptcl=.false.)
call pgrpsym%build_refspiral(eulspace)
call pgrpsym%build_refspiral(eulspace_sub)
call eulspace_sub%replace_with_closest(eulspace)
! need to reproject based on eulspace spiral
! allocate(refimgs(nspace))
refimgs = reproject(vol, eulspace)
! block_tree = gen_eulspace_block_tree2D(eulspace, eulspace_sub, pgrpsym, refimgs, params)
end program simple_test_block_tree2D