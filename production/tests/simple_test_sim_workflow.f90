program simple_test_sim_workflow
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_atoms,        only: commander_pdb2mrc
use simple_commanders_volops,       only: commander_reproject
use simple_commanders_pick,         only: commander_make_pickrefs, commander_pick
use simple_commanders_sim,          only: commander_simulate_particles, commander_simulate_movie
use simple_commanders_preprocess,   only: commander_ctf_estimate, commander_motion_correct
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cluster2D,    only: commander_cluster2D
implicit none
character(len=*), parameter         :: filetab_file='filetab.txt'
type(cmdline)                       :: cline_pdb2mrc, cline_projection, cline_make_pick_refs, cline_sim_mov, cline_ctf_est, cline_mot_corr
type(cmdline)                       :: cline_new_project, cline_import_movies, cline_segpick, cline_refpick, cline_sim_ptcls, cline_abinitio2D
type(cmdline)                       :: cline_cluster2D 
type(parameters)                    :: params
type(commander_new_project)         :: xnew_project
type(commander_pdb2mrc)             :: xpdb2mrc
type(commander_reproject)           :: xreproject
type(commander_make_pickrefs)       :: xmakepickrefs
type(commander_simulate_movie)      :: xsimov
type(commander_motion_correct)      :: xmotcorr
type(commander_ctf_estimate)        :: xctf_estimate
type(commander_import_movies)       :: ximport_movies
type(commander_pick)                :: xsegpick, xrefpick
type(commander_simulate_particles)  :: xsim_ptcls      
type(commander_abinitio2D)          :: xabinitio2D
type(commander_cluster2D)           :: xcluster2D
integer                             :: rc
type(string)                        :: cmd, projfile
logical                             :: mrc_exists
if( command_argument_count() .ne. 0 )then
   write(logfhandle,'(a)') 'ERROR! Usage: simple_test_sim_workflow'
   call exit(-1)
endif
! Download pdb & generate volume and reprojections
inquire(file="1JYX.mrc", exist=mrc_exists)
if( .not. mrc_exists )then
   write(*, *) 'Downloading the example dataset...'
   cmd = 'curl -s -o 1JYX.pdb https://files.rcsb.org/download/1JYX.pdb'
   call execute_command_line(cmd%to_char(), exitstat=rc)
   print *, 'Converting .pdb to .mrc...'
   call cline_pdb2mrc%set('smpd',            1.)
   call cline_pdb2mrc%set('pdbfile', '1JYX.pdb')
   call cline_pdb2mrc%checkvar('smpd',        1)
   call cline_pdb2mrc%checkvar('pdbfile',     2)
   call cline_pdb2mrc%check()
   call xpdb2mrc%execute_safe(cline_pdb2mrc)
   call cline_pdb2mrc%kill()
   cmd = 'rm 1JYX.pdb'
   call execute_command_line(cmd%to_char(), exitstat=rc)
endif

! Generate reprojections
print *, 'Projecting 1JYX.mrc...'
call cline_projection%set('vol1', '1JYX.mrc')
call cline_projection%set('smpd',        1.3)
call cline_projection%set('pgrp',       'c1')
call cline_projection%set('mskdiam',    180.)
call cline_projection%set('nspace',      10.)
call cline_projection%set('nthr',        16.)
call xreproject%execute_safe(cline_projection)
call cline_projection%kill()

! Make pick references
call cline_make_pick_refs%set('pickrefs', 'reprojs.mrcs')
call cline_make_pick_refs%set('nthr',                16.)
call xmakepickrefs%execute_safe(cline_make_pick_refs)
call cline_make_pick_refs%kill()

! Simulate movie
call cline_sim_mov%set('stk', 'reprojs.mrcs')
call cline_sim_mov%set('xdim',          4096)
call cline_sim_mov%set('ydim',          4096)
call cline_sim_mov%set('nthr',           16.)
call xsimov%execute_safe(cline_sim_mov)
call cline_sim_mov%kill()

! Project creation
call cline_new_project%set('projname',  params%projname)
call xnew_project%execute_safe(cline_new_project)
call cline_new_project%kill()
projfile = params%projname//'.simple'

! movie import
call cline_import_movies%set('prg',                'import_movies')
call cline_import_movies%set('mkdir',                        'yes')
call cline_import_movies%set('cs',                       params%cs)
call cline_import_movies%set('fraca',                 params%fraca)
call cline_import_movies%set('kv',                       params%kv)
call cline_import_movies%set('smpd',                   params%smpd)
call cline_import_movies%set('filetab',               filetab_file)
call cline_import_movies%set('ctf',                          'yes')
call ximport_movies%execute_safe(cline_import_movies)
call cline_import_movies%kill()

! Motion correction - algorithm iso
call cline_mot_corr%set('stk', 'simulated_movies.mrcs')
call cline_mot_corr%set('algorithm',    'iso')
call cline_mot_corr%set('nparts',           4)
call cline_mot_corr%set('nthr',            16)
call xmotcorr%execute_safe(cline_mot_corr)
call cline_mot_corr%kill()
! Motion correction - algorithm patch
call cline_mot_corr%set('stk',      'simulated_movies.mrcs')
call cline_mot_corr%set('algorithm',                'patch')
call cline_mot_corr%set('nparts',                         4)
call cline_mot_corr%set('nthr',                          16)
call xmotcorr%execute_safe(cline_mot_corr)
call cline_mot_corr%kill()

! CTF estimate - patch yes
call cline_ctf_est%set('ctfpatch', 'yes')
call cline_ctf_est%set('nparts',       4)
call cline_ctf_est%set('nthr',        16)
call xctf_estimate%execute_safe(cline_ctf_est)
call cline_ctf_est%kill()
! CTF estimate - patch yes
call cline_ctf_est%set('prg',     'ctf_estimate')
call cline_ctf_est%set('ctfpatch',          'no')
call cline_ctf_est%set('nparts',               4)
call cline_ctf_est%set('nthr',                16)
call xctf_estimate%execute_safe(cline_ctf_est)
call cline_ctf_est%kill()

! Segmentation-based picking
call cline_segpick%set('prg',       'pick')
call cline_segpick%set('picker', 'segdiam')
call cline_segpick%set('nparts',         4)
call cline_segpick%set('nthr',          16)      
call xsegpick%execute_safe(cline_segpick)
call cline_segpick%kill()  

! Reference-based picking
 call cline_refpick%set('prg',               'pick')
 call cline_refpick%set('pickrefs', 'pickrefs.mrcs')
 call cline_refpick%set('nparts',                 4)
 call cline_refpick%set('nthr',                  16)      
 call xrefpick%execute_safe(cline_refpick)
 call cline_refpick%kill()  

! 2D analysis
! generate simulate particles
call cline_sim_ptcls%set('prg', 'simulate_particles')
call cline_sim_ptcls%set('vol1',          '1JXY.mrc')
call cline_sim_ptcls%set('ctf',                'yes')
call cline_sim_ptcls%set('nptcls',               500)
call cline_sim_ptcls%set('smpd',                 1.3)
call cline_sim_ptcls%set('snr',                  0.5)
call cline_sim_ptcls%set('pgrp',                'c1')
call cline_sim_ptcls%set('mskdiam',              180)
call cline_sim_ptcls%set('nthr',                  16)
call xsim_ptcls%execute_safe(cline_sim_ptcls)
call cline_sim_ptcls%kill()
! abinitio2D
 call cline_abinitio2D%set('prg', 'abinitio2D')
 call cline_abinitio2D%set('mskdiam',      180)
 call cline_abinitio2D%set('ncls',         100)
 call cline_abinitio2D%set('nthr',          16)
 call xabinitio2D%execute_safe(cline_abinitio2D)
 call cline_abinitio2D%kill()

! cluster2D
call cline_cluster2D%set('prg', 'cluster2D')
call cline_cluster2D%set('mskdiam',     180)
call cline_cluster2D%set('ncls',        100)
call cline_cluster2D%set('nthr',         16)
call xcluster2D%execute_safe(cline_abinitio2D)
call cline_cluster2D%kill()

end program simple_test_sim_workflow 
