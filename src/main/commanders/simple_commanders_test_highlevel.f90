!@descr: for all highlevel tests
module simple_commanders_test_highlevel
use simple_commanders_api
use simple_stream_api
use simple_commanders_project_core, only: commander_new_project, commander_selection
use simple_commanders_project_mov,  only: commander_import_movies
use simple_commanders_atoms,        only: commander_pdb2mrc
use simple_commanders_volops,       only: commander_reproject
use simple_commanders_pick,         only: commander_make_pickrefs, commander_pick
use simple_commanders_sim,          only: commander_simulate_particles, commander_simulate_movie
use simple_commanders_preprocess,   only: commander_ctf_estimate, commander_motion_correct, commander_preprocess_distr
use simple_commanders_abinitio2D,   only: commander_abinitio2D
use simple_commanders_cluster2D,    only: commander_cluster2D
use simple_micproc,                 only: sample_filetab
use simple_commanders_validate,     only: commander_mini_stream
implicit none
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_test_mini_stream
  contains
    procedure :: execute      => exec_test_mini_stream
end type commander_test_mini_stream

type, extends(commander_base) :: commander_test_simulated_workflow
  contains
    procedure :: execute      => exec_test_simulated_workflow
end type commander_test_simulated_workflow

contains

subroutine exec_test_mini_stream( self, cline )
    class(commander_test_mini_stream),  intent(inout) :: self
    class(cmdline),                     intent(inout) :: cline
    real,         parameter          :: CTFRES_THRES = 8.0, ICE_THRES = 1.0, OVERSHOOT = 1.2
    type(string), allocatable        :: dataset_cmds(:)
    type(string), allocatable        :: micstab(:), filetab(:), movfnames(:)
    type(string)                     :: output_dir, imgkind
    integer,      allocatable        :: orimap(:)
    type(cmdline)                    :: cline_dataset, cline_new_project, cline_import_movies, cline_preprocess
    type(cmdline)                    :: cline_select, cline_mini_stream
    type(parameters)                 :: params
    type(commander_new_project)      :: xnew_project
    type(commander_preprocess_distr) :: xpreprocess
    type(commander_import_movies)    :: ximport_movies
    type(commander_mini_stream)      :: xmini_stream
    type(commander_selection)        :: xsel
    type(sp_project)                 :: spproj
    type(stream_watcher)             :: movie_buff
    integer                          :: i, ndata_sets, n_nonzero, nmovf
    type(string)                     :: abspath, projfile
    character(len=*), parameter      :: filetab_file='filetab.txt'
    ! Parsing
    if( command_argument_count() < 1 )then
        write(logfhandle,'(a)') 'ERROR! Usage: simple_test_mini_stream fname=filetab.txt'
        call exit(-1)
    else 
        call cline%parse_oldschool
    endif
    call cline%checkvar('fname',        1)
    call cline%check()
    call params%new(cline)
    call read_filetable(params%fname, dataset_cmds)
    ndata_sets=size(dataset_cmds)
    ! projname=system_name smpd=1.3 cs=2.7 kv=300 fraca=0.1 total_dose=53 dir_movies=/usr/local/data/movies gainref=gainref.mrc nparts=4 nthr=16 moldiam_max=200 nram=100
    call simple_getcwd(abspath)
    output_dir=abspath
    do i = 1, ndata_sets
        call cline_dataset%read(dataset_cmds(i)%to_char())
        call cline_dataset%checkvar('projname',        1)
        call cline_dataset%checkvar('smpd',            2)
        call cline_dataset%checkvar('cs',              3)
        call cline_dataset%checkvar('kv',              4)
        call cline_dataset%checkvar('fraca',           5)
        call cline_dataset%checkvar('total_dose',      6)
        call cline_dataset%checkvar('dir_movies',      7)
        call cline_dataset%checkvar('gainref',         8)
        call cline_dataset%checkvar('nparts',          9)
        call cline_dataset%checkvar('nthr',           10)
        call cline_dataset%checkvar('moldiam_max',    11)
        call cline_dataset%checkvar('nran',           12)
        call cline_dataset%check()
        call params%new(cline_dataset)
        call cline_dataset%kill()
        ! project creation
        call cline_new_project%set('projname',  params%projname)
        call xnew_project%execute_safe(cline_new_project)
        call cline_new_project%kill()
        projfile = params%projname//'.simple'
        ! create filetab with a subset of overshoot randomly selected movies
        movie_buff = stream_watcher(1,params%dir_movies)
        call movie_buff%watch(nmovf, movfnames)
        call movie_buff%kill
        filetab = sample_filetab(movfnames, ceiling(real(params%nran)*OVERSHOOT))
        call write_filetable(string(filetab_file), filetab)
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
        ! check either movies or micrographs
        call spproj%read(projfile)
        imgkind = spproj%get_mic_kind(1)
        if( imgkind.eq.'intg' )then
            ! nothing to do
        else
            ! preprocess
            call simple_chdir(output_dir//'/'//params%projname%to_char())
            call cline_preprocess%set('prg',                  'preprocess')
            call cline_preprocess%set('mkdir',                       'yes')
            call cline_preprocess%set('gainref',            params%gainref)
            call cline_preprocess%set('total_dose',      params%total_dose)
            call cline_preprocess%set('dfmin',               DFMIN_DEFAULT)
            call cline_preprocess%set('dfmax',               DFMAX_DEFAULT)
            call cline_preprocess%set('hp',                            30.)
            call cline_preprocess%set('lp',                             2.)
            call cline_preprocess%set('mcpatch',                      'no')
            call cline_preprocess%set('nparts',              params%nparts)
            call cline_preprocess%set('nthr',                  params%nthr)
            call cline_preprocess%check()
            call xpreprocess%execute_safe(cline_preprocess)
            call cline_preprocess%kill()
        endif
        ! reject based on CTF resolution and ice score
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call cline_select%delete('nran')
        call cline_select%set('prg',                           'selection')
        call cline_select%set('mkdir',                               'yes')
        call cline_select%set('oritype',                             'mic')
        call cline_select%set('ctfresthreshold',              CTFRES_THRES)
        call cline_select%set('icefracthreshold',                ICE_THRES)
        call xsel%execute_safe(cline_select)
        call cline_select%kill()
        ! state=0/1 should now be in project file on disk
        ! re-run for random selection
        call spproj%read(projfile)
        n_nonzero = spproj%get_n_insegment_state('mic', 1)
        if( n_nonzero > params%nran )then
            ! make random selection
            call simple_chdir(output_dir//'/'//params%projname%to_char())
            call cline_select%delete('ctfresthreshold')
            call cline_select%delete('icefracthreshold')
            call cline_select%set('prg',                       'selection')
            call cline_select%set('mkdir',                           'yes') 
            call cline_select%set('oritype',                         'mic')
            call cline_select%set('nran',                      params%nran)
            call xsel%execute_safe(cline_select)
            call cline_select%kill()
        endif
        call spproj%read(projfile)
        call spproj%get_mics_table(micstab, orimap)
        call simple_chdir(output_dir//'/'//params%projname%to_char())
        call write_filetable(string('intgs.txt'),micstab)
        ! mini stream 
        call cline_mini_stream%set('prg',                    'mini_stream')
        call cline_mini_stream%set('mkdir',                          'yes')
        call cline_mini_stream%set('filetab',                  'intgs.txt')
        call cline_mini_stream%set('smpd',                     params%smpd)
        call cline_mini_stream%set('fraca',                   params%fraca)
        call cline_mini_stream%set('kv',                         params%kv)
        call cline_mini_stream%set('cs',                         params%cs)
        call cline_mini_stream%set('moldiam_max',       params%moldiam_max)
        call cline_mini_stream%set('nparts',                 params%nparts)
        call cline_mini_stream%set('nthr',                     params%nthr)
        call xmini_stream%execute_safe(cline_mini_stream)
        call cline_dataset%kill()
        call cline_new_project%kill()
        call cline_import_movies%kill()
        call cline_preprocess%kill()
        call cline_mini_stream%kill()
        call cline_dataset%kill()
        call simple_chdir(output_dir)
    enddo
    call simple_end('**** SIMPLE_TEST_MINI_STREAM_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_mini_stream

subroutine exec_test_simulated_workflow( self, cline )
    class(commander_test_simulated_workflow), intent(inout) :: self
    class(cmdline),                           intent(inout) :: cline
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
    call simple_end('**** SIMPLE_TEST_SIMULATED_WORKFLOW_WORKFLOW NORMAL STOP ****')
end subroutine exec_test_simulated_workflow

end module simple_commanders_test_highlevel
