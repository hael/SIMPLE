! Standalone end-to-end test for a simulated single-particle workflow.
program simple_test_1jyx_abinitio
use simple_core_module_api
use simple_atoms,                    only: atoms
use simple_molecule_data,            only: molecule_data, betagal_1jyx
use simple_cmdline,                  only: cmdline
use simple_commanders_sim,           only: commander_simulate_particles
use simple_commanders_project_core,  only: commander_new_project
use simple_commanders_project_ptcl,  only: commander_import_particles
use simple_commanders_abinitio2D,    only: commander_abinitio2D
use simple_commanders_abinitio,      only: commander_abinitio3D
use simple_procimgstk,               only: add_noise_imgfile
use simple_imghead,                  only: find_ldim_nptcls
use simple_sp_project,               only: sp_project
use simple_ui,                       only: make_ui
implicit none
#include "simple_local_flags.inc"

character(len=*), parameter :: WORKDIR    = 'test_1jyx_abinitio'
character(len=*), parameter :: PROJNAME   = 'onejyx_abinitio'
character(len=*), parameter :: PROJFILE   = PROJNAME//'.simple'
character(len=*), parameter :: VOL_FILE   = '1JYX.mrc'
character(len=*), parameter :: CLEAN_FILE = '1JYX_particles_clean.mrcs'
character(len=*), parameter :: NOISY_FILE = '1JYX_particles_noisy.mrcs'
character(len=*), parameter :: ORI_FILE   = '1JYX_particles_oris.txt'

real    :: smpd, mskdiam, noise_snr
integer :: nptcls, ncls, nthr
type(cmdline) :: user_args, cline_sim, cline_new_project, cline_import
type(cmdline) :: cline_abinitio2D, cline_abinitio3D
type(commander_simulate_particles) :: xsimulate_particles
type(commander_new_project)        :: xnew_project
type(commander_import_particles)   :: ximport_particles
type(commander_abinitio2D)         :: xabinitio2D
type(commander_abinitio3D)         :: xabinitio3D
type(molecule_data)                :: mol
type(atoms)                        :: molecule
type(sp_project)                   :: spproj
type(string) :: original_cwd, workflow_root, project_root, project_path
type(string) :: volume_path, clean_path, noisy_path, orientations_path
integer      :: ldim(3), nimgs, status

! Defaults keep the test reasonably small while leaving enough particles for
! both classifications.  They can be overridden with old-style key=value args.
smpd      = 1.3
mskdiam   = 180.0
noise_snr = 0.10
nptcls    = 200
ncls      = 4
nthr      = 4
if( command_argument_count() > 0 )then
    call user_args%parse_oldschool
    if( user_args%defined('smpd')     ) smpd      = user_args%get_rarg('smpd')
    if( user_args%defined('mskdiam')  ) mskdiam   = user_args%get_rarg('mskdiam')
    if( user_args%defined('snr')      ) noise_snr = user_args%get_rarg('snr')
    if( user_args%defined('nptcls')   ) nptcls    = user_args%get_iarg('nptcls')
    if( user_args%defined('ncls')     ) ncls      = user_args%get_iarg('ncls')
    if( user_args%defined('nthr')     ) nthr      = user_args%get_iarg('nthr')
endif
if( smpd <= 0.0 )      THROW_HARD('smpd must be positive')
if( mskdiam <= 0.0 )   THROW_HARD('mskdiam must be positive')
if( noise_snr <= 0.0 ) THROW_HARD('snr must be positive')
if( nptcls < 8 )       THROW_HARD('nptcls must be at least 8')
if( ncls < 1 .or. ncls >= nptcls ) THROW_HARD('ncls must be in [1,nptcls)')
if( nthr < 1 )         THROW_HARD('nthr must be positive')

! Direct commander calls need the normal SIMPLE UI metadata for mkdir=yes.
call make_ui
call simple_getcwd(original_cwd)
if( file_exists(WORKDIR) )then
    call simple_rmdir(WORKDIR, status)
    if( status /= 0 ) THROW_HARD('could not reset '//WORKDIR)
endif
call simple_mkdir(WORKDIR)
call simple_chdir(WORKDIR, status)
if( status /= 0 ) THROW_HARD('could not enter '//WORKDIR)
call simple_getcwd(workflow_root)

write(logfhandle,'(a)') '>>> Step 1/6: generate a volume from embedded 1JYX beta-galactosidase coordinates'
mol = betagal_1jyx()
call molecule%pdb2mrc(volfile=string(VOL_FILE), smpd=smpd, mol=mol, center_pdb=.true.)
call molecule%kill()
volume_path = simple_abspath(string(VOL_FILE))
call find_ldim_nptcls(volume_path, ldim, nimgs)
if( any(ldim < 1) .or. ldim(3) <= 1 ) THROW_HARD('1JYX volume has invalid dimensions')

write(logfhandle,'(a)') '>>> Step 2/6: simulate a clean stack at random orientations'
call cline_sim%set('prg',        'simulate_particles')
call cline_sim%set('mkdir',                       'no')
call cline_sim%set('vol1',                volume_path)
call cline_sim%set('outstk',               CLEAN_FILE)
call cline_sim%set('outfile',                 ORI_FILE)
call cline_sim%set('smpd',                       smpd)
call cline_sim%set('mskdiam',                 mskdiam)
call cline_sim%set('nptcls',                   nptcls)
call cline_sim%set('nthr',                       nthr)
call cline_sim%set('pgrp',                       'c1')
call cline_sim%set('ctf',                         'no')
call cline_sim%set('snr',                         10.0) ! simimg adds no noise for SNR >= 5
call cline_sim%set('bfac',                          0.0)
call cline_sim%set('sherr',                        0.0)
call xsimulate_particles%execute(cline_sim)
call cline_sim%kill()
clean_path        = simple_abspath(string(CLEAN_FILE))
orientations_path = simple_abspath(string(ORI_FILE))
call find_ldim_nptcls(clean_path, ldim, nimgs)
if( nimgs /= nptcls ) THROW_HARD('clean simulated stack has the wrong particle count')

write(logfhandle,'(a,f7.3)') '>>> Step 3/6: add Gaussian noise to the stack; SNR = ', noise_snr
call add_noise_imgfile(clean_path, string(NOISY_FILE), noise_snr, smpd)
noisy_path = simple_abspath(string(NOISY_FILE))
call find_ldim_nptcls(noisy_path, ldim, nimgs)
if( nimgs /= nptcls ) THROW_HARD('noisy simulated stack has the wrong particle count')

write(logfhandle,'(a)') '>>> Step 4/6: create a SIMPLE project and import the noisy particles'
call cline_new_project%set('prg',          'new_project')
call cline_new_project%set('projname',          PROJNAME)
call cline_new_project%set('qsys_name',          'local')
call xnew_project%execute(cline_new_project)
call cline_new_project%kill()
call simple_getcwd(project_root)
project_path = simple_abspath(string(PROJFILE))
call cline_import%set('prg',             'import_particles')
call cline_import%set('mkdir',                           'no')
call cline_import%set('projfile',               project_path)
call cline_import%set('stk',                       noisy_path)
call cline_import%set('oritab',             orientations_path)
call cline_import%set('smpd',                            smpd)
call cline_import%set('ctf',                              'no')
call ximport_particles%execute(cline_import)
call cline_import%kill()

write(logfhandle,'(a)') '>>> Step 5/6: run ab initio 2D classification'
call cline_abinitio2D%set('prg',          'abinitio2D')
call cline_abinitio2D%set('projfile',     project_path)
call cline_abinitio2D%set('mkdir',               'yes')
call cline_abinitio2D%set('mskdiam',           mskdiam)
call cline_abinitio2D%set('ncls',                 ncls)
call cline_abinitio2D%set('nstages',                 1)
call cline_abinitio2D%set('nits_per_stage',          1)
call cline_abinitio2D%set('nthr',                 nthr)
call xabinitio2D%execute(cline_abinitio2D)
call cline_abinitio2D%kill()
call capture_stage_project('abinitio2D')
call spproj%read(project_path)
if( spproj%os_cls2D%get_noris() < 1 ) THROW_HARD('abinitio2D produced no class averages')
call spproj%kill()

write(logfhandle,'(a)') '>>> Step 6/6: run ab initio 3D reconstruction'
call cline_abinitio3D%set('prg',            'abinitio3D')
call cline_abinitio3D%set('projfile',       project_path)
call cline_abinitio3D%set('mkdir',                 'yes')
call cline_abinitio3D%set('pgrp',                   'c1')
call cline_abinitio3D%set('pgrp_start',             'c1')
call cline_abinitio3D%set('mskdiam',             mskdiam)
call cline_abinitio3D%set('nstates',                   1)
call cline_abinitio3D%set('multivol_mode',      'single')
call cline_abinitio3D%set('filt_mode',            'none')
call cline_abinitio3D%set('automsk',                'no')
call cline_abinitio3D%set('nstages',                   1)
call cline_abinitio3D%set('nsample',              nptcls)
call cline_abinitio3D%set('nparts',                    1)
call cline_abinitio3D%set('nthr',                   nthr)
call xabinitio3D%execute(cline_abinitio3D)
call cline_abinitio3D%kill()
call capture_stage_project('abinitio3D')
call spproj%read(project_path)
if( spproj%get_nptcls() /= nptcls ) THROW_HARD('final project has the wrong particle count')
if( spproj%os_ptcl3D%get_noris() /= nptcls ) THROW_HARD('abinitio3D produced incomplete particle orientations')
call spproj%kill()

call simple_chdir(original_cwd, status)
if( status /= 0 ) THROW_HARD('could not restore the original working directory')
write(logfhandle,'(a)') '>>> Results: '//workflow_root%to_char()
call simple_end('**** SIMPLE_TEST_1JYX_ABINITIO NORMAL STOP ****')

contains

    subroutine capture_stage_project( stage )
        character(len=*), intent(in) :: stage
        if( file_exists(PROJFILE) ) project_path = simple_abspath(string(PROJFILE))
        call simple_chdir(project_root, status)
        if( status /= 0 ) THROW_HARD('could not leave the '//trim(stage)//' directory')
    end subroutine capture_stage_project

end program simple_test_1jyx_abinitio
