! concrete commander: stream processing routines
module simple_commander_stream_wflows
include 'simple_lib.f08'
use simple_cmdline,        only: cmdline
use simple_commander_base, only: commander_base
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_qsys_funs,      only: qsys_cleanup
use simple_parameters,     only: parameters, params_glob
implicit none

public :: preprocess_stream_commander
public :: cluster2D_stream_distr_commander
public :: pick_extract_stream_distr_commander
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: preprocess_stream_commander
  contains
    procedure :: execute      => exec_preprocess_stream
end type preprocess_stream_commander
type, extends(commander_base) :: cluster2D_stream_distr_commander
  contains
    procedure :: execute      => exec_cluster2D_stream_distr
end type cluster2D_stream_distr_commander
type, extends(commander_base) :: pick_extract_stream_distr_commander
  contains
    procedure :: execute      => exec_pick_extract_stream_distr
end type pick_extract_stream_distr_commander

contains

    subroutine exec_preprocess_stream( self, cline )
        use simple_moviewatcher,         only: moviewatcher
        use simple_qsys_funs,            only: qsys_cleanup
        use simple_commander_preprocess, only: preprocess_commander
        class(preprocess_stream_commander), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)                       :: params
        integer,                   parameter   :: SHORTTIME = 60   ! folder watched every minute
        integer,                   parameter   :: LONGTIME  = 600  ! time lag after which a movie is processed
        class(cmdline),            allocatable :: completed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj
        character(len=LONGSTRLEN), allocatable :: movies(:), prev_movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker, imgkind
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
        character(len=LONGSTRLEN)              :: movie
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
        integer                                :: nptcls, nptcls_prev, nmovs, nmovs_prev, cnt, i, n_os_out
        logical                                :: l_pick
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call cline%set('numlen', real(5))
        call cline%set('stream','yes')
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        ! read in movies
        call spproj%read( params%projfile )
        ! picking
        l_pick = .false.
        if( cline%defined('refs') .or. cline%defined('vol1') ) l_pick = .true.
        ! check for previously processed movies
        call spproj%get_movies_table(prev_movies)
        ! output directories
        output_dir = PATH_HERE
        output_dir_ctf_estimate   = filepath(trim(output_dir), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(output_dir), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        if( l_pick )then
            output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
            output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
            call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
            call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        endif
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! prepares picking references
        if( l_pick )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call cline_make_pickrefs%set('stream','no')
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', trim(PICKREFS)//params%ext)
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME)
        call spproj%get_movies_table(prev_movies)
        call movie_buff%add_to_history(prev_movies)
        call spproj%get_mics_table(prev_movies)
        call movie_buff%add_to_history(prev_movies)
        ! start watching
        prev_stacksz = 0
        nmovies      = 0
        iter         = 0
        do
            if( file_exists(trim(TERM_STREAM)) )then
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> PREPROCES STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(SHORTTIME)
            enddo
            iter = iter + 1
            call movie_buff%watch( nmovies, movies )
            ! append movies to processing stack
            if( nmovies > 0 )then
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    call create_individual_project
                    call qenv%qscripts%add_to_streaming( cline )
                enddo
            endif
            ! stream scheduling
            call qenv%qscripts%schedule_streaming( qenv%qdescr )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I5)')'>>> MOVIES/MICROGRAPHS TO PROCESS: ', stacksz
            endif
            ! completed jobs update the current project
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                ! append new processed movies to project
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                nptcls_prev = spproj%get_nptcls()
                nmovs_prev  = spproj%os_mic%get_noris()
                do icline=1,size(completed_jobs_clines)
                    stream_spprojfile = completed_jobs_clines(icline)%get_carg('projfile')
                    call stream_spproj%read( stream_spprojfile )
                    call spproj%append_project(stream_spproj, 'mic')
                    if( l_pick )then
                        call spproj%append_project(stream_spproj, 'stk')
                        ! transfer ptcl2D box coordinates
                        do i=1,stream_spproj%os_ptcl2D%get_noris()
                            cnt = spproj%os_ptcl2D%get_noris()-stream_spproj%os_ptcl2D%get_noris()+i
                            ! box coordinates
                            if(stream_spproj%os_ptcl2D%isthere(i,'xpos'))then
                                call spproj%os_ptcl2D%set(cnt,'xpos',stream_spproj%os_ptcl2D%get(i,'xpos'))
                            endif
                            if(stream_spproj%os_ptcl2D%isthere(i,'ypos'))then
                                call spproj%os_ptcl2D%set(cnt,'ypos',stream_spproj%os_ptcl2D%get(i,'ypos'))
                            endif
                        enddo
                    endif
                    call stream_spproj%kill()
                    deallocate(stream_spprojfile)
                enddo
                nptcls = spproj%get_nptcls()
                nmovs  = spproj%os_mic%get_noris()
                ! write
                call spproj%write
                ! update for 2d streaming
                call update_projects_list
                deallocate(completed_jobs_clines)
            endif
            ! wait
            call simple_sleep(SHORTTIME)
        end do
        ! termination
        call spproj%write
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            subroutine update_projects_list
                type(cmdline) :: cline_mov
                character(len=:),          allocatable :: fname, abs_fname
                character(len=LONGSTRLEN), allocatable :: old_fnames(:), fnames(:)
                integer :: i, n_spprojs, n_old
                n_spprojs = size(completed_jobs_clines)
                if( n_spprojs == 0 )return
                if( file_exists(STREAM_SPPROJFILES) )then
                    ! append
                    call read_filetable(STREAM_SPPROJFILES, old_fnames)
                    n_old = size(old_fnames)
                    allocate(fnames(n_spprojs+n_old))
                    fnames(1:n_old) = old_fnames(:)
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
                        fnames(n_old+i) = trim(abs_fname)
                        deallocate(abs_fname)
                    enddo
                else
                    ! first write
                    allocate(fnames(n_spprojs))
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 2')
                        fnames(i) = trim(abs_fname)
                        deallocate(abs_fname)
                    enddo
                endif
                call write_filetable(STREAM_SPPROJFILES, fnames)
            end subroutine update_projects_list

            subroutine create_individual_project
                type(sp_project)              :: spproj_here
                type(cmdline)                 :: cline_here
                type(ctfparams)               :: ctfvars
                character(len=STDLEN)         :: ext, movie_here
                character(len=LONGSTRLEN)     :: projname, projfile
                movie_here = basename(trim(movie))
                ext        = fname2ext(trim(movie_here))
                projname   = 'preprocess_'//trim(get_fbody(trim(movie_here), trim(ext)))
                projfile   = trim(projname)//trim(METADATA_EXT)
                call cline_here%set('projname', trim(projname))
                call cline_here%set('projfile', trim(projfile))
                call spproj_here%update_projinfo(cline_here)
                spproj_here%compenv  = spproj%compenv
                spproj_here%jobproc  = spproj%jobproc
                ctfvars%ctfflag      = CTFFLAG_YES
                ctfvars%smpd         = params%smpd
                ctfvars%cs           = params%cs
                ctfvars%kv           = params%kv
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = params%phaseplate.eq.'yes'
                call spproj_here%add_single_movie(trim(movie), ctfvars)
                call spproj_here%write
                call spproj_here%kill
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
            end subroutine create_individual_project

    end subroutine exec_preprocess_stream

    subroutine exec_cluster2D_stream_distr( self, cline )
        use simple_projection_frcs,        only: projection_frcs
        use simple_commander_distr_wflows, only: cluster2D_distr_commander, make_cavgs_distr_commander
        use simple_commander_cluster2D,    only: rank_cavgs_commander
        use simple_image,                  only: image
        use simple_oris,                   only: oris
        use simple_ori,                    only: ori
        class(cluster2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        !integer,               parameter   :: CCRES_NPTCLS_LIM    = 10000 ! # of ptcls required to turn on objfun=ccres
        integer,               parameter   :: WAIT_WATCHER        = 30    ! seconds prior to new stack detection
        integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
        integer,               parameter   :: MINBOXSZ            = 72   ! minimum boxsize for scaling
        real                               :: SMPD_TARGET         = 5.    ! target sampling distance
        ! dev settings
        ! integer,               parameter   :: CCRES_NPTCLS_LIM    = 10000 ! # of ptcls required to turn on objfun=ccres
        ! integer,               parameter   :: WAIT_WATCHER        = 30    ! seconds prior to new stack detection
        ! integer,               parameter   :: ORIGPROJ_WRITEFREQ  = 60    ! 10mins, Frequency at which the original project file should be updated
        character(len=STDLEN), parameter   :: MICS_SELECTION_FILE = 'stream2D_selection.txt'
        character(len=STDLEN), parameter   :: PROJFILE_BUFFER     = 'buffer.simple'
        character(len=STDLEN), parameter   :: PROJFILE_POOL       = 'pool.simple'
        character(len=STDLEN), parameter   :: PROJFILE2D          = 'cluster2D.simple'
        character(len=STDLEN), parameter   :: SCALE_DIR           = './scaled_stks/'
        logical,               parameter   :: debug_here = .false.
        type(parameters)                   :: params
        type(make_cavgs_distr_commander)   :: xmake_cavgs
        type(rank_cavgs_commander)         :: xrank_cavgs
        type(cmdline)                      :: cline_cluster2D, cline_cluster2D_buffer
        type(cmdline)                      :: cline_make_cavgs, cline_rank_cavgs
        type(sp_project)                   :: orig_proj, stream_proj, buffer_proj, pool_proj
        type(ctfparams)                    :: ctfvars
        type(oris)                         :: os_stk
        type(ori)                          :: o_stk
        character(LONGSTRLEN), allocatable :: spproj_list(:), stk_list(:)
        character(len=:),      allocatable :: spproj_list_fname, orig_projfile, stk
        character(len=STDLEN)              :: str_iter, refs_glob, refs_glob_ranked
        character(len=LONGSTRLEN)          :: buffer_dir
        real    :: orig_smpd, msk, scale_factor, orig_msk, smpd
        integer :: iter, orig_box, box, nptcls_glob, iproj, ncls_glob, n_transfers
        integer :: nptcls_glob_prev, n_spprojs, orig_nparts, last_injection, nparts
        integer :: origproj_time, max_ncls, nptcls_per_buffer, buffer_ptcls_range(2), pool_iter
        logical :: do_autoscale, l_maxed, buffer_exists, do_wait
        ! seed the random number generator
        call seed_rnd
        ! set oritype
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'ptcl2D')
        call cline%set('stream','yes') ! only for parameters determination
        call params%new(cline)
        ! sanity
        if( .not.file_exists(params%projfile) )then
            THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        endif
        orig_projfile = trim(params%projfile)
        if( .not.file_exists(params%dir_target) )then
            THROW_HARD('folder: '//trim(params%dir_target)//' does not exist!')
        endif
        call cline%set('stream','no') ! was only for parameters determination
        call cline%set('mkdir','no')
        ! init
        do_autoscale      = params%autoscale.eq.'yes'
        orig_nparts       = params%nparts
        max_ncls          = floor(real(params%ncls)/real(params%ncls_start))*params%ncls_start ! effective maximum # of classes
        nptcls_per_buffer = params%nptcls_per_cls*params%ncls_start         ! # of particles in each buffer
        buffer_exists     = .false.                                         ! whether the buffer exists
        l_maxed           = .false.                                         ! whether all chunks have been merged
        do_wait           = .true.
        spproj_list_fname = filepath(trim(params%dir_target), trim(STREAM_SPPROJFILES))
        ncls_glob         = 0
        n_transfers       = 0
        ! for microscopes that don't work too good, automatically turned off after 2 hours
        if(.not.cline%defined('time_inactive'))params%time_inactive = 2*3600
        ! init command-lines
        call cline%delete('lp') ! gold-standard refinement
        cline_cluster2D         = cline
        cline_cluster2D_buffer  = cline
        cline_make_cavgs        = cline
        ! buffer classification
        ! down-scaling for fast execution, aggressive stochastic optimisation, no match filter, bi-linear interpolation,
        ! no incremental learning, objective function default is standard cross-correlation (cc), no centering
        call cline_cluster2D_buffer%set('prg',       'cluster2D')
        call cline_cluster2D_buffer%set('projfile',  trim(PROJFILE_BUFFER))
        call cline_cluster2D_buffer%set('projname',  trim(get_fbody(trim(PROJFILE_BUFFER),trim('simple'))))
        call cline_cluster2D_buffer%set('objfun',    'cc')
        call cline_cluster2D_buffer%set('wfun',      'bilinear')
        call cline_cluster2D_buffer%set('center',    'no')
        call cline_cluster2D_buffer%set('match_filt','no')
        call cline_cluster2D_buffer%set('autoscale', 'no')
        call cline_cluster2D_buffer%set('stream',    'yes') ! more stringent convergence
        call cline_cluster2D_buffer%set('refine',    'snhc')
        call cline_cluster2D_buffer%delete('update_frac')
        ! pool classification
        ! down-scaling for fast execution, stochastic optimisation, optional match filter, bi-linear interpolation,
        ! no incremental learning, objective function is standard cross-correlation (cc)
        call cline_cluster2D%set('prg',       'cluster2D')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('extr_iter', 100.)
        call cline_cluster2D%set('trs',       MINSHIFT)
        call cline_cluster2D%set('wfun',      'bilinear')
        call cline_cluster2D%set('projfile',  trim(PROJFILE_POOL))
        call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE_POOL),trim('simple'))))
        call cline_cluster2D%set('objfun',    'cc')
        if( cline%defined('objfun') )then
            if( cline%get_carg('objfun').eq.'ccres' )then
                call cline_cluster2D%set('objfun', 'ccres')
                call cline_cluster2D%set('bfac',   1000.)
            endif
        endif
        if( .not.cline%defined('match_filt') ) call cline_cluster2D%set('match_filt','no')
        call cline_cluster2D%delete('update_frac')
        ! scaling
        call cline_make_cavgs%set('prg',    'make_cavgs')
        call cline_make_cavgs%set('wfun',   'bilinear')
        call cline_make_cavgs%delete('autoscale')
        call cline_make_cavgs%delete('remap_cls')
        ! WAIT FOR FIRST STACKS
        nptcls_glob = 0
        do
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_mics
                    write(logfhandle,'(A,I8,A,A)')'>>> # OF PARTICLES: ', nptcls_glob, ' : ',cast_time_char(simple_gettime())
                    call flush(6)
                    if( nptcls_glob > nptcls_per_buffer )then
                        exit ! Enough particles to initiate cluster2D
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! transfer projects info, rename & rewrite
        call orig_proj%read(params%projfile)
        pool_proj%projinfo = orig_proj%projinfo
        pool_proj%compenv  = orig_proj%compenv
        if( orig_proj%jobproc%get_noris()>0 ) pool_proj%jobproc = orig_proj%jobproc
        call pool_proj%projinfo%delete_entry('projname')
        call pool_proj%projinfo%delete_entry('projfile')
        call pool_proj%update_projinfo(cline_cluster2D)
        ! getting general parameters from the first sp_project
        call stream_proj%read(trim(spproj_list(1)))
        orig_box  = stream_proj%get_box()
        orig_smpd = stream_proj%get_smpd()
        orig_msk  = params%msk
        call stream_proj%kill
        params%smpd_targets2D(1) = max(orig_smpd, params%lp*LP2SMPDFAC)
        if( do_autoscale )then
            if( orig_box < MINBOXSZ )then
                do_autoscale = .false.
            else
                call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                if( box < MINBOXSZ )then
                    SMPD_TARGET = orig_smpd * real(orig_box) / real(MINBOXSZ)
                    call autoscale(orig_box, orig_smpd, SMPD_TARGET, box, smpd, scale_factor)
                endif
                if( box == orig_box ) do_autoscale = .false.
            endif
        endif
        if( do_autoscale )then
            msk = orig_msk * scale_factor
        else
            smpd = orig_smpd
            box  = orig_box
            msk  = orig_msk
            scale_factor = 1.
        endif
        ! FIRST IMPORT
        allocate(stk_list(n_spprojs))
        call os_stk%new(n_spprojs)
        do iproj=1,n_spprojs
            call stream_proj%read(spproj_list(iproj))
            o_stk           = stream_proj%os_stk%get_ori(1)
            ctfvars         = stream_proj%get_ctfparams('ptcl2D', 1)
            ctfvars%smpd    = smpd
            stk             = stream_proj%get_stkname(1)
            stk_list(iproj) = trim(stk)
            call o_stk%set_ctfvars(ctfvars)
            call os_stk%set_ori(iproj, o_stk)
        enddo
        ! updates & scales stacks
        call scale_stks( stk_list ) ! must come first as names updated
        call pool_proj%add_stktab(stk_list, os_stk)
        nptcls_glob = pool_proj%get_nptcls()
        if( file_exists(MICS_SELECTION_FILE) )then
            call flag_selection
        else
            do iproj=1,pool_proj%os_mic%get_noris()
                call pool_proj%os_mic%set(iproj,'state',1.)
            enddo
            do iproj=1,pool_proj%os_stk%get_noris()
                call pool_proj%os_stk%set(iproj,'state',1.)
            enddo
        endif
        call pool_proj%write
        ! generates buffer project
        buffer_ptcls_range = 0
        call gen_buffer_from_pool
        ! MAIN LOOP
        last_injection = simple_gettime()
        origproj_time  = last_injection
        pool_iter      = 1
        do iter = 1,9999
            str_iter  = int2str_pad(iter,3)
            pool_iter = min(999,pool_iter)
            write(logfhandle,'(A,I3)')'>>> WAIT CYCLE ',iter
            if( is_timeout(simple_gettime()) )exit
            if( buffer_exists )then
                ! ten iterations of the buffer
                call classify_buffer
                ! particles selection
                call reject_from_buffer
                ! book keeping
                call transfer_buffer_to_pool
                buffer_exists = .false.
                ! cleanup
                call buffer_proj%kill
                call simple_rmdir(buffer_dir)
            endif
            ! one iteration of the whole dataset when not converged,
            ! or after transfer from buffer or each 5 iterations
            if( .not.cline_cluster2D%defined('converged') .or. mod(iter,5)==0 )then
                call classify_pool
                pool_iter = pool_iter+1
                ! rejection each 5 iterations
                if( mod(iter,5)==0 )then
                    call reject_from_pool
                    call pool_proj%write
                endif
            endif
            ! termination and/or pause
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> CLUSTER2D STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(WAIT_WATCHER)
            enddo
            if( file_exists(TERM_STREAM) )then
                write(logfhandle,'(A,A)')'>>> TERMINATING CLUSTER2D STREAM ',cast_time_char(simple_gettime())
                exit
            endif
            if( do_wait ) call simple_sleep(WAIT_WATCHER)
            ! detect new project files
            call append_new_mics
            ! optionally builds new buffer
            call gen_buffer_from_pool
            ! update original project
            if( simple_gettime()-origproj_time > ORIGPROJ_WRITEFREQ )then
                call update_orig_proj
                origproj_time = simple_gettime()
            endif
        enddo
        call qsys_cleanup
        ! updates original project
        call update_orig_proj
        ! class averages at original sampling
        if ( do_autoscale )then
            nptcls_glob = orig_proj%get_nptcls()
            nparts      = calc_nparts(orig_proj, nptcls_glob)
            call orig_proj%kill
            call cline_make_cavgs%set('nparts', real(nparts))
            call cline_make_cavgs%set('ncls',   real(ncls_glob))
            call cline_make_cavgs%set('refs',   refs_glob)
            call xmake_cavgs%execute(cline_make_cavgs)
            call orig_proj%read(orig_projfile)
        endif
        call orig_proj%add_cavgs2os_out(refs_glob, orig_smpd)
        call orig_proj%add_frcs2os_out(trim(FRCS_FILE), 'frc2D')
        call orig_proj%write_segment_inside('out',orig_projfile)
        call orig_proj%kill()
        ! ranking
        refs_glob_ranked = add2fbody(refs_glob,params%ext,'_ranked')
        call cline_rank_cavgs%set('projfile', orig_projfile)
        call cline_rank_cavgs%set('stk',      refs_glob)
        call cline_rank_cavgs%set('outstk',   trim(refs_glob_ranked))
        call xrank_cavgs%execute(cline_rank_cavgs)
        ! cleanup
        call qsys_cleanup
        call simple_rmdir(SCALE_DIR)
        call del_file(PROJFILE_POOL)
        call del_file(PROJFILE_BUFFER)
        call del_file(PROJFILE2D)
        ! end gracefully
        call simple_end('**** SIMPLE_CLUSTER2D_STREAM NORMAL STOP ****')
        contains

            subroutine reject_from_buffer
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                character(STDLEN)    :: refs_buffer
                real                 :: ave, sdev, maxv, minv
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, endit, ncls_here, cnt
                if( debug_here ) print *,'in reject from_buffer'; call flush(6)
                ncls_rejected   = 0
                nptcls_rejected = 0
                ncls_here       = buffer_proj%os_cls2D%get_noris()
                boxmatch        = find_boxmatch(box, msk)
                endit           = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer     = trim(buffer_dir)//'/cavgs_iter'//int2str_pad(endit,3)//params%ext
                allocate(cls_mask(ncls_here), source=.true.)
                if( debug_here )then
                    ! variance
                    call img%new([box,box,1],smpd)
                    do icls=1,ncls_here
                        call img%read(refs_buffer,icls)
                        call img%stats(ave, sdev, maxv, minv)
                        call buffer_proj%os_cls2D%set(icls,'sdev',sdev)
                    enddo
                endif
                ! resolution and correlation
                call buffer_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,params%ndev)
                if( debug_here ) call buffer_proj%os_cls2D%write('buffer_'//trim(str_iter)//'.txt')
                if( any(cls_mask) )then
                    ncls_rejected = count(.not.cls_mask)
                    do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                        if( buffer_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(buffer_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call buffer_proj%os_ptcl2D%set(iptcl,'state',0.)
                    enddo
                    do icls=1,ncls_here
                        if( .not.cls_mask(icls) )then
                            call buffer_proj%os_cls2D%set(icls,'pop',0.)
                            call buffer_proj%os_cls2D%set(icls,'corr',-1.)
                        endif
                    enddo
                    call img%new([box,box,1],smpd)
                    cnt = 0
                    do icls=1,ncls_here
                        if( cls_mask(icls) ) cycle
                        cnt = cnt+1
                        !if( debug_here )then
                            call img%read(refs_buffer,icls)
                            call img%write('rejected_'//int2str(pool_iter)//'.mrc',cnt)
                        !endif
                        img = 0.
                        call img%write(refs_buffer,icls)
                    enddo
                    call img%read(refs_buffer, params%ncls_start)
                    call img%write(refs_buffer,params%ncls_start)
                    call img%kill
                    deallocate(cls_mask)
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM BUFFER: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                endif
                if( debug_here ) print *,'end reject from_buffer'; call flush(6)
            end subroutine reject_from_buffer

            subroutine reject_from_pool
                type(image)          :: img
                logical, allocatable :: cls_mask(:)
                real                 :: ave, sdev, minv, maxv, ndev_here
                integer              :: nptcls_rejected, ncls_rejected, iptcl
                integer              :: boxmatch, icls, cnt
                if( debug_here ) print *,'in reject from_pool'; call flush(6)
                ncls_rejected   = 0
                nptcls_rejected = 0
                boxmatch        = find_boxmatch(box, msk)
                allocate(cls_mask(ncls_glob), source=.true.)
                if( debug_here )then
                    call img%new([box,box,1],smpd)
                    do icls=1,ncls_glob
                        if( pool_proj%os_cls2D%get_state(icls)==0 .or. pool_proj%os_cls2D%get(icls,'pop')<1. ) cycle
                        call img%read(refs_glob,icls)
                        call img%stats(ave, sdev, maxv, minv)
                        call pool_proj%os_cls2D%set(icls,'sdev',sdev)
                    enddo
                endif
                if( debug_here )call pool_proj%os_cls2D%write('pool_'//int2str(pool_iter)//'.txt')
                ! correlation & resolution
                ndev_here = 1.5*params%ndev ! less stringent rejection
                call pool_proj%os_cls2D%find_best_classes(boxmatch,smpd,params%lpthresh,cls_mask,ndev_here)
                if( debug_here )call pool_proj%os_cls2D%write('pool_aftersel_'//int2str(pool_iter)//'.txt')
                if( .not.all(cls_mask) )then
                    ncls_rejected = count(.not.cls_mask)
                    do iptcl=1,pool_proj%os_ptcl2D%get_noris()
                        if( pool_proj%os_ptcl2D%get_state(iptcl) == 0 )cycle
                        icls = nint(pool_proj%os_ptcl2D%get(iptcl,'class'))
                        if( cls_mask(icls) ) cycle
                        nptcls_rejected = nptcls_rejected+1
                        call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                    enddo
                    if( nptcls_rejected > 0 )then
                        do icls=1,ncls_glob
                            if( .not.cls_mask(icls) )then
                                call pool_proj%os_cls2D%set(icls,'pop',0.)
                                call pool_proj%os_cls2D%set(icls,'corr',-1.)
                            endif
                        enddo
                        cnt = 0
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_mask(icls) ) cycle
                            cnt = cnt+1
                            !if( debug_here )then
                                call img%read(refs_glob,icls)
                                call img%write('rejected_pool_'//int2str(pool_iter)//'.mrc',cnt)
                            !endif
                            img = 0.
                            call img%write(refs_glob,icls)
                        enddo
                        call img%read(refs_glob, ncls_glob)
                        call img%write(refs_glob, ncls_glob)
                        call img%kill
                        deallocate(cls_mask)
                        write(logfhandle,'(A,I4,A,I6,A)')'>>> REJECTED FROM POOL: ',nptcls_rejected,' PARTICLES IN ',ncls_rejected,' CLUSTERS'
                    endif
                else
                    write(logfhandle,'(A,I4,A,I6,A)')'>>> NO PARTICLES FLAGGED FOR REJECTION FROM POOL'
                endif
                if( debug_here ) print *,'end reject from_pool'; call flush(6)
            end subroutine reject_from_pool

            subroutine append_new_mics
                integer :: n_spprojs_prev, n_new_spprojs, cnt, iptcl, n_new_ptcls
                if( debug_here ) print *,'in append_new_mics'; call flush(6)
                if( .not.is_file_open(spproj_list_fname) )then
                    do_wait        = .false.
                    n_spprojs_prev = n_spprojs
                    n_spprojs      = nlines(spproj_list_fname)
                    n_new_spprojs  = n_spprojs - n_spprojs_prev
                    if( n_new_spprojs > 0 )then
                        ! fetch new stacks
                        n_new_ptcls = 0
                        if(allocated(spproj_list))deallocate(spproj_list)
                        if(allocated(stk_list))   deallocate(stk_list)
                        allocate(stk_list(n_new_spprojs))
                        call read_filetable(spproj_list_fname, spproj_list)
                        cnt = 0
                        do iproj=n_spprojs_prev+1,n_spprojs
                            cnt = cnt + 1
                            call stream_proj%read(spproj_list(iproj))
                            n_new_ptcls   = n_new_ptcls + stream_proj%get_nptcls()
                            stk           = stream_proj%get_stkname(1)
                            stk_list(cnt) = trim(stk)
                        enddo
                        if( n_new_ptcls > 0 )then
                            call scale_stks( stk_list )
                            ! update project with new images
                            nptcls_glob = pool_proj%get_nptcls() + n_new_ptcls
                            cnt = 0
                            do iproj=n_spprojs_prev+1,n_spprojs
                                cnt = cnt + 1
                                call stream_proj%read(spproj_list(iproj))
                                ctfvars      = stream_proj%get_ctfparams('ptcl2D', 1)
                                ctfvars%smpd = smpd
                                call pool_proj%add_stk(stk_list(cnt), ctfvars)
                            enddo
                            do iptcl=nptcls_glob-n_new_ptcls+1,nptcls_glob
                                call pool_proj%os_ptcl2D%set(iptcl,'state',0.) ! deactivate by default
                            enddo
                            write(logfhandle,'(A,I8,A,A)')'>>> # OF PARTICLES: ', nptcls_glob, ' ; ',cast_time_char(simple_gettime())
                            last_injection = simple_gettime()
                            call cline_cluster2D%delete('converged') ! reactivates pool classification
                        endif
                    endif
                    ! updates pool with selection
                    if( file_exists(MICS_SELECTION_FILE) ) call flag_selection
                else
                    do_wait = .true.
                endif
                if( debug_here ) print *,'end append_new_mics'; call flush(6)
            end subroutine append_new_mics

            !>  flags mics and stacks segements with selection
            subroutine flag_selection
                character(len=:), allocatable :: mic_name, mic_name_from_proj
                type(oris)                    :: mics_sel
                integer                       :: nmics, imic, iproj
                logical                       :: included
                nmics = nlines(MICS_SELECTION_FILE)
                call mics_sel%new(nmics)
                call mics_sel%read(MICS_SELECTION_FILE)
                do iproj=1,pool_proj%os_mic%get_noris()
                    if(.not.pool_proj%os_mic%isthere('intg'))cycle
                    call pool_proj%os_mic%getter(iproj,'intg',mic_name_from_proj)
                    ! check whether corresponding mic is selected
                    included = .true.
                    do imic=1,nmics
                        call mics_sel%getter(imic,'intg',mic_name)
                        if( trim(mic_name).eq.trim(mic_name_from_proj) )then
                            included  = mics_sel%get_state(imic) == 1
                            exit
                        endif
                    enddo
                    if( included )then
                        call pool_proj%os_mic%set(iproj,'state',1.)
                        call pool_proj%os_stk%set(iproj,'state',1.)
                    else
                        write(logfhandle,'(A,A)')'>>> DESELECTING MICROGRAPH: ',trim(mic_name)
                        call pool_proj%os_mic%set(iproj,'state',0.)
                        call pool_proj%os_stk%set(iproj,'state',0.)
                    endif
                enddo
                call mics_sel%kill
            end subroutine flag_selection

            integer function calc_nparts(spproj, nptcls_in)
                type(sp_project), intent(in) :: spproj
                integer,          intent(in) :: nptcls_in
                integer :: i, tailsz
                ! adjust number of parts
                tailsz = 0
                do i=nptcls_in,1,-1
                    if( spproj%os_ptcl2D%get_state(i) == 0 )then
                        tailsz = tailsz + 1
                    else
                        exit
                    endif
                enddo
                do calc_nparts=orig_nparts,1,-1
                    if(real(nptcls_in)/real(calc_nparts) > real(tailsz) )exit
                enddo
            end function calc_nparts

            !>  runs iterations of cluster2D for the buffer in a separate folder
            subroutine classify_buffer
                use simple_commander_hlev_wflows, only: cluster2D_autoscale_commander
                type(cluster2D_distr_commander) :: xcluster2D_distr
                integer :: nptcls, nparts
                if( debug_here ) print *,'in classify_buffer'; call flush(6)
                write(logfhandle,'(A)')'>>> 2D CLASSIFICATION OF NEW BUFFER'
                ! directory structure
                call chdir('buffer2D')
                ! init
                nptcls = buffer_proj%get_nptcls()
                nparts = calc_nparts(buffer_proj, nptcls)
                call buffer_proj%kill
                ! cluster2d execution
                params_glob%projfile = trim('./'//trim(PROJFILE_BUFFER))
                call cline_cluster2D_buffer%set('startit',   1.)
                call cline_cluster2D_buffer%set('maxits',    12.) ! guaranties 3 iterations with withdrawal
                call cline_cluster2D_buffer%set('extr_iter', real(MAX_EXTRLIM2D-2))
                call cline_cluster2D_buffer%set('ncls',      real(params%ncls_start))
                call cline_cluster2D_buffer%set('nparts',    real(nparts))
                call cline_cluster2D_buffer%set('box',       real(box))
                call cline_cluster2D_buffer%set('msk',       real(box/2)-3.)
                call cline_cluster2D_buffer%delete('trs')
                call cline_cluster2D_buffer%delete('endit')
                call cline_cluster2D_buffer%delete('converged')
                params_glob%nptcls = nptcls
                call xcluster2D_distr%execute(cline_cluster2D_buffer)
                call getcwd(buffer_dir)
                call rename(PROJFILE_BUFFER, '../'//trim(PROJFILE_BUFFER))
                call chdir('..')
                call buffer_proj%read(PROJFILE_BUFFER)
                params_glob%nparts   = orig_nparts
                params_glob%projfile = trim(orig_projfile)
                if( debug_here ) print *,'end classify_buffer'; call flush(6)
            end subroutine classify_buffer

            !>  runs one iteration of cluster2D for the merged buffers in the cwd
            subroutine classify_pool
                type(cluster2D_distr_commander)  :: xcluster2D_distr
                type(sp_project)                 :: spproj2D
                character(len=:), allocatable    :: prev_refs
                integer :: iptcl, nptcls, nptcls_sel, istk, nstks, nptcls_glob
                ! transfer from pool
                call cline_cluster2D%set('projfile',  trim(PROJFILE2D))
                call cline_cluster2D%set('projname',  trim(get_fbody(trim(PROJFILE2D),trim('simple'))))
                spproj2D%projinfo = pool_proj%projinfo
                spproj2D%compenv  = pool_proj%compenv
                if( pool_proj%jobproc%get_noris()>0 ) spproj2D%jobproc = pool_proj%jobproc
                call spproj2D%projinfo%delete_entry('projname')
                call spproj2D%projinfo%delete_entry('projfile')
                call spproj2D%update_projinfo(cline_cluster2D)
                nptcls_glob = pool_proj%get_nptcls()
                do iptcl=nptcls_glob,1,-1
                    if( pool_proj%os_ptcl2D%get_state(iptcl) == 1 )exit
                enddo
                nstks  = nint(pool_proj%os_ptcl2D%get(iptcl, 'stkind'))
                nptcls = nint(pool_proj%os_stk%get(nstks, 'top'))
                call spproj2D%os_stk%new(nstks)
                call spproj2D%os_ptcl2D%new(nptcls)
                do iptcl=1,nptcls
                    call spproj2D%os_ptcl2D%set_ori(iptcl, pool_proj%os_ptcl2D%get_ori(iptcl))
                enddo
                do istk=1,nstks
                    call spproj2D%os_stk%set_ori(istk, pool_proj%os_stk%get_ori(istk))
                enddo
                spproj2D%os_cls2D = pool_proj%os_cls2D
                spproj2D%os_out   = pool_proj%os_out
                call spproj2D%write(PROJFILE2D)
                nptcls_sel = spproj2D%os_ptcl2D%get_noris(consider_state=.true.)
                write(logfhandle,'(A,I8,A,I4,A)')'>>> 2D CLASSIFICATION OF POOL: ',nptcls_sel,' PARTICLES IN ',ncls_glob,' CLUSTERS'
                ! cluster2d execution
                params_glob%projfile = trim(PROJFILE_POOL)
                call cline_cluster2D%set('startit', real(pool_iter))
                call cline_cluster2D%set('maxits',  real(pool_iter))
                call cline_cluster2D%set('ncls',    real(ncls_glob))
                call cline_cluster2D%set('box',     real(box))
                call cline_cluster2D%set('msk',     real(msk))
                call cline_cluster2D%set('refs',    trim(refs_glob))
                call cline_cluster2D%set('frcs',    trim(FRCS_FILE))
                call cline_cluster2D%delete('endit')
                call cline_cluster2D%delete('converged')
                params_glob%nptcls = nptcls
                call xcluster2D_distr%execute(cline_cluster2D)
                params_glob%projfile = trim(orig_projfile)
                refs_glob            = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter,3))//trim(params%ext)
                params_glob%nptcls   = nptcls_glob
                ! transfer back to pool
                call spproj2D%read(PROJFILE2D)
                do iptcl=1,nptcls
                    call pool_proj%os_ptcl2D%set_ori(iptcl, spproj2D%os_ptcl2D%get_ori(iptcl))
                enddo
                do istk=1,nstks
                    call pool_proj%os_stk%set_ori(istk, spproj2D%os_stk%get_ori(istk))
                enddo
                pool_proj%os_cls2D = spproj2D%os_cls2D
                pool_proj%os_out   = spproj2D%os_out
                ! removes previous references
                if(pool_iter>1)then
                    prev_refs = trim(CAVGS_ITER_FBODY)//trim(int2str_pad(pool_iter-1,3))//trim(params%ext)
                    call del_file(prev_refs)
                    call del_file(add2fbody(trim(prev_refs),params%ext,'_even'))
                    call del_file(add2fbody(trim(prev_refs),params%ext,'_odd'))
                endif
                ! cleanup
                call spproj2D%kill
            end subroutine classify_pool

            !>  append the classified buffer to the pool
            subroutine transfer_buffer_to_pool
                type(projection_frcs)         :: frcs_glob, frcs_buffer, frcs_prev
                type(image)                   :: img
                type(ori)                     :: o
                character(len=:), allocatable :: refs_buffer, stkout, stkin
                integer,          allocatable :: cls_pop(:), cls_buffer_pop(:), pinds(:)
                integer                       :: endit, iptcl, ind, state, ncls_here, icls, i, cnt, stat
                real                          :: stkind
                n_transfers = n_transfers+1
                write(logfhandle,'(A,I4)')'>>> TRANSFER BUFFER PARTICLES CLASSIFICATION TO POOL #',n_transfers
                ! max # of classes reached ?
                l_maxed = ncls_glob >= max_ncls
                ! updates # of classes
                if( .not.l_maxed ) ncls_here = ncls_glob+params%ncls_start
                ! transfer class parameters to pool
                if( l_maxed )then
                    cls_pop = nint(pool_proj%os_cls2D%get_all('pop'))
                else
                    if( ncls_glob == 0 )then
                        pool_proj%os_cls2D = buffer_proj%os_cls2D
                    else
                        call pool_proj%os_cls2D%reallocate(ncls_here)
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            o   = buffer_proj%os_cls2D%get_ori(icls)
                            call o%set('class',real(ind))
                            call pool_proj%os_cls2D%set_ori(ind,o)
                        enddo
                    endif
                endif
                ! transfer particles parameters to pool
                do iptcl=1,buffer_proj%os_ptcl2D%get_noris()
                    o = buffer_proj%os_ptcl2D%get_ori(iptcl)
                    ! greedy search
                    call o%set('updatecnt', 1.)
                    ! updates class
                    if( l_maxed )then
                        icls = irnd_uni(ncls_glob)
                        do while(cls_pop(icls)==0)
                            icls = irnd_uni(ncls_glob)
                        enddo
                        call o%set('class',real(icls))
                    else
                        icls = nint(o%get('class'))
                        call o%set('class',real(ncls_glob+icls))
                    endif
                    ind    = buffer_ptcls_range(1)+iptcl-1
                    stkind = pool_proj%os_ptcl2D%get(ind,'stkind')
                    call o%set('stkind',stkind) ! preserves stack index
                    call pool_proj%os_ptcl2D%set_ori(ind,o)
                enddo
                do iptcl=ind+1,pool_proj%os_ptcl2D%get_noris()
                    call pool_proj%os_ptcl2D%set(iptcl,'state',0.)
                enddo
                ! transfer references to pool
                endit       = nint(cline_cluster2D_buffer%get_rarg('endit'))
                refs_buffer = trim(buffer_dir)//'/cavgs_iter'//int2str_pad(endit,3)//params%ext
                if( .not.l_maxed )then
                    if( ncls_glob == 0 )then
                        ! first time
                        refs_glob = 'start_cavgs'//params%ext
                        call simple_copy_file(refs_buffer,refs_glob)
                        stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                        stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                        call simple_copy_file(stkin,stkout)
                        stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                        stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call simple_copy_file(stkin,stkout)
                        call simple_copy_file(trim(buffer_dir)//'/'//trim(FRCS_FILE),FRCS_FILE)
                    else
                        ! class averages &
                        call img%new([box,box,1],smpd)
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call img%read(refs_buffer, icls)
                            call img%write(refs_glob, ind)
                        enddo
                        stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                        stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call img%read(stkin, icls)
                            call img%write(stkout, ind)
                        enddo
                        stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                        stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                        do icls=1,params%ncls_start
                            ind = ncls_glob+icls
                            call img%read(stkin, icls)
                            call img%write(stkout, ind)
                        enddo
                        ! FRCs
                        state = 1
                        call frcs_prev%new(ncls_glob, box, smpd, state)
                        call frcs_buffer%new(params%ncls_start, box, smpd, state)
                        call frcs_glob%new(ncls_here, box, smpd, state)
                        call frcs_prev%read(FRCS_FILE)
                        call frcs_buffer%read(trim(buffer_dir)//'/'//trim(FRCS_FILE))
                        do icls=1,ncls_glob
                            call frcs_glob%set_frc(icls,frcs_prev%get_frc(icls, box, state), state)
                        enddo
                        ind = 0
                        do icls=ncls_glob+1,ncls_here
                            ind = ind + 1
                            call frcs_glob%set_frc(icls,frcs_buffer%get_frc(ind, box, state), state)
                        enddo
                        call frcs_glob%write(FRCS_FILE)
                        ! cleanup
                        call frcs_prev%kill
                        call frcs_glob%kill
                        call frcs_buffer%kill
                        call img%kill
                    endif
                    ! global parameters
                    ncls_glob = ncls_here
                endif
                ! remapping
                if( l_maxed .and. trim(params%remap_cls).eq.'yes' )then
                    if( any(cls_pop==0) )then
                        state = 1
                        cls_buffer_pop = nint(buffer_proj%os_cls2D%get_all('pop'))
                        call frcs_buffer%new(params%ncls_start, box, smpd, state)
                        call frcs_glob%new(ncls_glob, box, smpd, state)
                        call frcs_glob%read(FRCS_FILE)
                        call frcs_buffer%read(trim(buffer_dir)//'/'//trim(FRCS_FILE))
                        call img%new([box,box,1],smpd)
                        do icls=1,ncls_glob
                            if( cls_pop(icls)>0 )cycle
                            cnt = 0
                            ind = irnd_uni(params%ncls_start)
                            do while( cls_buffer_pop(ind)==0 )
                                ind = irnd_uni(params%ncls_start)
                                 ! pathological cases and insufficient # of classes
                                cnt = cnt + 1
                                if( cnt>params%ncls_start )exit
                            enddo
                            if( cnt>params%ncls_start )cycle
                            cls_buffer_pop(ind) = 0 ! excludes from being picked again
                            ! classes
                            call img%read(refs_buffer, ind)
                            call img%write(refs_glob, icls)
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_even')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                            call img%read(stkin, ind)
                            call img%write(stkout, icls)
                            stkin  = add2fbody(trim(refs_buffer),params%ext,'_odd')
                            stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                            call img%read(stkin, ind)
                            call img%write(stkout, icls)
                            o = buffer_proj%os_cls2D%get_ori(ind)
                            call o%set('class',real(icls))
                            call pool_proj%os_cls2D%set_ori(icls,o)
                            ! frcs
                            call frcs_glob%set_frc(icls,frcs_buffer%get_frc(ind, box, state), state)
                            ! assignments
                            call buffer_proj%os_ptcl2D%get_pinds(ind,'class',pinds,consider_w=.false.)
                            pinds = pinds + buffer_ptcls_range(1)
                            do i=1,size(pinds)
                                iptcl = pinds(i)
                                call pool_proj%os_ptcl2D%set(iptcl,'class',real(icls))
                            enddo
                        enddo
                        call img%read(refs_glob, ncls_glob)
                        call img%write(refs_glob, ncls_glob)
                        stkout = add2fbody(trim(refs_glob),params%ext,'_even')
                        call img%read(stkout, ncls_glob)
                        call img%write(stkout, ncls_glob)
                        stkout = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call img%read(stkout, ncls_glob)
                        call img%write(stkout, ncls_glob)
                        call frcs_glob%write(FRCS_FILE)
                        call frcs_glob%kill
                        call frcs_buffer%kill
                        call img%kill
                    endif
                endif
                ! preserve buffer
                stat = rename(refs_buffer,'cavgs_buffer'//int2str_pad(n_transfers,4)//params%ext)
                ! cleanup
                call o%kill
                if(allocated(cls_pop))deallocate(cls_pop)
                if(allocated(cls_buffer_pop))deallocate(cls_buffer_pop)
                if( debug_here )print *,'end transfer_buffer_to_pool'; call flush(6)
            end subroutine transfer_buffer_to_pool

            subroutine read_mics
                character(len=:), allocatable :: mic_name, mic_name_from_proj
                type(oris)                    :: mics_sel
                integer                       :: nptcls, nmics, imic
                logical                       :: included, do_selection
                if( debug_here )print *,'in read_mics'; call flush(6)
                do_selection = file_exists(MICS_SELECTION_FILE)
                call read_filetable(spproj_list_fname, spproj_list)
                if( do_selection )then
                    nmics = nlines(MICS_SELECTION_FILE)
                    call mics_sel%new(nmics)
                    call mics_sel%read(MICS_SELECTION_FILE)
                endif
                ! determine number of particles
                nptcls_glob_prev = nptcls_glob
                nptcls_glob      = 0
                if( allocated(spproj_list) )then
                    n_spprojs = size(spproj_list)
                    do iproj = 1,n_spprojs
                        call stream_proj%read(spproj_list(iproj))
                        nptcls   = stream_proj%get_nptcls()
                        included = .true. ! included by default
                        if( do_selection )then
                            call stream_proj%os_mic%getter(1,'intg',mic_name_from_proj)
                            ! check whether corresponding mic is selected
                            do imic=1,nmics
                                call mics_sel%getter(imic,'intg',mic_name)
                                if( trim(mic_name).eq.trim(mic_name_from_proj) )then
                                    included  = mics_sel%get_state(imic) == 1
                                    exit
                                endif
                            enddo
                        endif
                        if( included ) nptcls_glob = nptcls_glob + nptcls
                    enddo
                    call stream_proj%kill
                else
                    n_spprojs = 0
                endif
                if(do_selection)call mics_sel%kill
                if( debug_here )print *,'end read_mics'; call flush(6)
            end subroutine read_mics

            !> builds and writes new buffer project from the pool
            subroutine gen_buffer_from_pool
                type(ctfparams) :: ctfparms
                character(len=:), allocatable :: stkname
                integer :: iptcl, stkind, ind_in_stk, imic, micind, nptcls_here
                integer :: fromp, top, istk, state, nmics_here, nptcls_tot
                if( debug_here )print *,'in gen_buffer_from_pool'; call flush(6)
                buffer_exists = .false.
                call buffer_proj%kill
                ! to account for directory structure
                call simple_rmdir('buffer2D') ! clean previous round
                call simple_mkdir('buffer2D')
                call simple_chdir('buffer2D')
                ! determines whether there are enough particles for a buffer
                if( buffer_ptcls_range(2) > 0 )then
                    call pool_proj%map_ptcl_ind2stk_ind('ptcl2D', buffer_ptcls_range(2), stkind, ind_in_stk)
                else
                    stkind = 0
                endif
                micind = stkind + 1
                nptcls_here = 0
                do istk=micind,pool_proj%os_stk%get_noris()
                    if(.not.pool_proj%os_stk%isthere(istk,'state'))then
                        write(logfhandle,*)'error: missing state flag; gen_buffer_from_pool'
                        stop
                    endif
                    state = pool_proj%os_stk%get_state(istk)
                    if( state == 0 )cycle
                    fromp       = nint(pool_proj%os_stk%get(istk,'fromp'))
                    top         = nint(pool_proj%os_stk%get(istk,'top'))
                    nptcls_here = nptcls_here + top-fromp+1
                enddo
                if( nptcls_here < nptcls_per_buffer )then
                    call simple_chdir('..')
                    return
                endif
                ! builds buffer if enough particles
                buffer_proj%projinfo = pool_proj%projinfo
                buffer_proj%compenv  = pool_proj%compenv
                if( pool_proj%jobproc%get_noris()>0 ) buffer_proj%jobproc = pool_proj%jobproc
                call buffer_proj%projinfo%delete_entry('projname')
                call buffer_proj%projinfo%delete_entry('projfile')
                call buffer_proj%update_projinfo(cline_cluster2D_buffer)
                istk        = 0
                nptcls_here = 0
                nmics_here  = 0
                nptcls_tot  = 0
                do imic=micind,pool_proj%os_stk%get_noris()
                    state = pool_proj%os_stk%get_state(imic)
                    istk  = istk+1
                    call pool_proj%os_stk%getter(imic,'stk',stkname)
                    ! to account for directory structure
                    if( stkname(1:1) /= '/') stkname = '../'//trim(stkname)
                    ! import
                    ctfparms = pool_proj%os_stk%get_ctfvars(imic)
                    call buffer_proj%add_stk(stkname, ctfparms)
                    fromp      = nint(buffer_proj%os_stk%get(istk,'fromp'))
                    top        = nint(buffer_proj%os_stk%get(istk,'top'))
                    call buffer_proj%os_stk%set(istk,'state',real(state))
                    do iptcl=fromp,top
                        call buffer_proj%os_ptcl2D%set(iptcl,'state',real(state))
                    enddo
                    if( state == 0 ) cycle
                    nmics_here  = nmics_here+1
                    nptcls_here = nptcls_here + top-fromp+1
                    if( nptcls_here > nptcls_per_buffer )exit
                enddo
                buffer_ptcls_range(1) = buffer_ptcls_range(2)+1
                buffer_ptcls_range(2) = buffer_ptcls_range(1)+top-1
                call buffer_proj%write
                call simple_chdir('..')
                buffer_exists = .true.
                write(logfhandle,'(A,I4,A,I6,A)')'>>> BUILT NEW BUFFER WITH ', nmics_here, ' MICROGRAPHS, ',nptcls_here, ' PARTICLES'
            end subroutine gen_buffer_from_pool

            logical function is_timeout( time_now )
                integer, intent(in) :: time_now
                is_timeout = .false.
                if(time_now-last_injection > params%time_inactive)then
                    write(logfhandle,'(A,A)')'>>> TIME LIMIT WITHOUT NEW IMAGES REACHED: ',cast_time_char(time_now)
                    is_timeout = .true.
                    if( .not.cline_cluster2D%defined('converged') )is_timeout = .false.
                else if(time_now-last_injection > 3600)then
                    write(logfhandle,'(A,A)')'>>> OVER ONE HOUR WITHOUT NEW PARTICLES: ',cast_time_char(time_now)
                    call flush(6)
                endif
                return
            end function is_timeout

            subroutine update_orig_proj
                type(ori)                          :: o
                type(oris)                         :: o_stks
                type(ctfparams)                    :: ctfparms
                character(LONGSTRLEN), allocatable :: stk_list(:), mic_list(:)
                real                               :: rstate
                integer                            :: nprev, n2append, cnt
                logical                            :: has_mics
                write(logfhandle,'(A,A)')'>>> UPDATING PROJECT AT: ',cast_time_char(simple_gettime())
                call orig_proj%read(orig_projfile)
                nprev    = orig_proj%os_stk%get_noris()
                n2append = n_spprojs-nprev
                has_mics = .false.
                if( n2append > 0 )then
                    allocate(stk_list(n2append),mic_list(n2append))
                    call o_stks%new(n2append)
                    cnt = 0
                    do iproj=nprev+1,n_spprojs
                        cnt = cnt + 1
                        ! stk
                        call stream_proj%read(spproj_list(iproj))
                        o             = stream_proj%os_stk%get_ori(1)
                        stk_list(cnt) = trim(stream_proj%get_stkname(1))
                        ctfparms      = stream_proj%get_ctfparams('ptcl2D', 1)
                        call o%set_ctfvars(ctfparms)
                        call o_stks%set_ori(cnt,o)
                        ! mic
                        if( stream_proj%get_nintgs() == 1 )then
                            has_mics      = .true.
                            ctfparms      = stream_proj%get_micparams(1)
                            mic_list(cnt) = trim(stream_proj%os_mic%get_static(1,'intg'))
                        endif
                    enddo
                    call orig_proj%add_stktab(stk_list, o_stks)
                    if( has_mics )then
                        call orig_proj%add_movies(mic_list, ctfparms)
                    endif
                    ! updates 2D, os_out not updated as not correct scale
                    do iproj=nprev+1,n_spprojs
                        rstate = pool_proj%os_stk%get(iproj,'state')
                        call orig_proj%os_stk%set(iproj,'state', rstate)
                        if( has_mics )call orig_proj%os_mic%set(iproj,'state', rstate)
                    enddo
                    orig_proj%os_cls2D  = pool_proj%os_cls2D
                    orig_proj%os_ptcl2D = pool_proj%os_ptcl2D
                    if( do_autoscale )call orig_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
                    call orig_proj%write
                endif
                ! cleanup
                call o%kill
                call o_stks%kill
                if(allocated(mic_list))deallocate(mic_list)
                if(allocated(stk_list))deallocate(stk_list)
                if( debug_here )print *,'end update_orig_proj'; call flush(6)
            end subroutine update_orig_proj

            subroutine scale_stks( stk_fnames )
                character(len=*), allocatable, intent(inout) :: stk_fnames(:)
                character(len=*), parameter :: SCALE_FILETAB = 'stkscale.txt'
                character(len=:), allocatable :: fname
                type(qsys_env) :: qenv
                type(cmdline)  :: cline_scale
                integer        :: istk
                if( .not.do_autoscale )return
                if( .not.allocated(stk_fnames) )return
                call simple_mkdir(SCALE_DIR, errmsg= "commander_stream_wflows:: cluster2D_stream_distr scale_stks")
                call qenv%new(params%nparts)
                call cline_scale%set('prg',        'scale')
                call cline_scale%set('smpd',       orig_smpd)
                call cline_scale%set('box',        real(orig_box))
                call cline_scale%set('newbox',     real(box))
                call cline_scale%set('filetab',    trim(SCALE_FILETAB))
                call cline_scale%set('nthr',       real(params%nthr))
                call cline_scale%set('dir_target', trim(SCALE_DIR))
                call cline_scale%set('stream',     'yes')
                call write_filetable(trim(SCALE_FILETAB), stk_fnames)
                call qenv%exec_simple_prg_in_queue(cline_scale, 'JOB_FINISHED_1')
                call qsys_cleanup
                do istk=1,size(stk_fnames)
                    fname            = add2fbody(stk_fnames(istk), params%ext, SCALE_SUFFIX)
                    stk_fnames(istk) = filepath(trim(SCALE_DIR), basename(fname))
                enddo
                call del_file(SCALE_FILETAB)
            end subroutine scale_stks

    end subroutine exec_cluster2D_stream_distr

    subroutine exec_pick_extract_stream_distr( self, cline )
        class(pick_extract_stream_distr_commander), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        integer,          parameter   :: WAIT_WATCHER        = 10    ! seconds prior to new stack detection
        !integer,          parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
        type(parameters)                    :: params
        type(cmdline)                       :: cline_pick_extract, cline_make_pickrefs
        type(sp_project)                    :: orig_proj, stream_proj
        type(qsys_env)                      :: qenv
        class(cmdline),         allocatable :: completed_jobs_clines(:)
        character(LONGSTRLEN),  allocatable :: spproj_list(:)
        character(len=:),       allocatable :: spproj_list_fname, output_dir, output_dir_picker, output_dir_extract, projfname
        integer :: iter, origproj_time, tnow, iproj, icline, nptcls, prev_stacksz, stacksz
        integer :: last_injection, n_spprojs, n_spprojs_prev, n_newspprojs, nmics
        ! output command line executed
        write(logfhandle,'(a)') '>>> COMMAND LINE EXECUTED'
        write(logfhandle,*) trim(cmdline_glob)
        ! set oritype & defaults
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call cline%set('stream','yes')
        call cline%set('numlen', 5.)
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        if( .not.file_exists(params%projfile) )THROW_HARD('project file: '//trim(params%projfile)//' does not exist!')
        spproj_list_fname = filepath(trim(params%dir_target),trim(STREAM_SPPROJFILES))
        ! read project info
        call orig_proj%read(params%projfile)
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! output directories
        output_dir = PATH_HERE
        output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
        output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
        call simple_mkdir(trim(output_dir),errmsg="commander_stream_wflows :: exec_pick_extract_stream_distr;  ")
        call simple_mkdir(trim(output_dir_picker),errmsg="commander_stream_wflows :: exec_pick_extract_stream_distr;  ")
        call simple_mkdir(trim(output_dir_extract),errmsg="commander_stream_wflows :: exec_pick_extract_stream_distr;  ")
        ! init command-lines
        cline_pick_extract  = cline
        cline_make_pickrefs = cline
        call cline_pick_extract%set('prg', 'pick_extract')
        call cline_pick_extract%delete('projname')
        ! prepares picking references
        call cline_make_pickrefs%set('prg','make_pickrefs')
        call cline_make_pickrefs%set('stream','no')
        call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
        call cline_pick_extract%set('refs', trim(PICKREFS)//params%ext)
        call cline_pick_extract%delete('vol1')
        ! wait for the first stacks
        last_injection  = simple_gettime()
        origproj_time   = last_injection
        prev_stacksz    = 0
        n_spprojs       = 0
        n_spprojs_prev  = 0
        do iter = 1,999999
            tnow = simple_gettime()
            if(tnow-last_injection > params%time_inactive .and. stacksz==0)then
                write(logfhandle,*)'>>> TIME LIMIT WITHOUT NEW MICROGRAPHS REACHED'
                exit
            endif
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_filetable(spproj_list_fname, spproj_list)
                    if( allocated(spproj_list) )n_spprojs = size(spproj_list)
                    n_newspprojs = n_spprojs - n_spprojs_prev
                    if( n_newspprojs > 0 )then
                        ! copy projects and add to processing stack
                        do iproj = n_spprojs_prev+1, n_spprojs
                            call stream_proj%read(spproj_list(iproj))
                            projfname  = filepath(PATH_HERE, basename(spproj_list(iproj)))
                            call stream_proj%write(projfname)
                            call cline_pick_extract%set('projfile', trim(projfname))
                            call qenv%qscripts%add_to_streaming( cline_pick_extract )
                            call stream_proj%kill
                            deallocate(projfname)
                        enddo
                        n_spprojs_prev = n_spprojs
                    endif
                    ! streaming scheduling
                    call qenv%qscripts%schedule_streaming( qenv%qdescr )
                    stacksz = qenv%qscripts%get_stacksz()
                    if( stacksz .ne. prev_stacksz )then
                        prev_stacksz = stacksz
                        write(logfhandle,'(A,I5)')'>>> MICROGRAPHS TO PROCESS: ', stacksz
                    endif
                    ! completed jobs update the current project
                    if( qenv%qscripts%get_done_stacksz() > 0 )then
                        call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                        do icline=1,size(completed_jobs_clines)
                            projfname = completed_jobs_clines(icline)%get_carg('projfile')
                            call stream_proj%read( projfname )
                            call orig_proj%append_project(stream_proj, 'mic')
                            call orig_proj%append_project(stream_proj, 'stk')
                            call stream_proj%kill()
                            deallocate(projfname)
                        enddo
                        nptcls = orig_proj%get_nptcls()
                        nmics  = orig_proj%get_nintgs()
                        write(logfhandle,'(A,I8)')'>>> NEW MICROGRAPHS COUNT: ', nmics
                        write(logfhandle,'(A,I8)')'>>> NEW PARTICLES   COUNT: ', nptcls
                        call orig_proj%write
                        deallocate(completed_jobs_clines)
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_PICK_EXTRACT_STREAM NORMAL STOP ****')
    end subroutine exec_pick_extract_stream_distr


end module simple_commander_stream_wflows
