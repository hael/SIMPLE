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
        integer,                   parameter   :: LONGTIME  = 600
        class(cmdline),            allocatable :: completed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj
        character(len=LONGSTRLEN), allocatable :: movies(:), prev_movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract, stream_spprojfile
        character(len=LONGSTRLEN)              :: movie
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, icline
        integer                                :: nptcls, nptcls_prev, nmovs, nmovs_prev
        logical                                :: l_pick
        if( .not. cline%defined('oritype') ) call cline%set('oritype', 'mic')
        call cline%set('numlen', real(5))
        call cline%set('stream','yes')
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        l_pick = cline%defined('refs')
        ! read in movies
        call spproj%read( params%projfile )
        ! check for previously processed movies
        call spproj%get_movies_table(prev_movies)
        ! output directories
        output_dir = PATH_HERE
        output_dir_ctf_estimate   = filepath(trim(output_dir), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(output_dir), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        if( l_pick )then
            output_dir_picker  = filepath(trim(output_dir), trim(DIR_PICKER))
            output_dir_extract = filepath(trim(output_dir), trim(DIR_EXTRACT))
            call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
            call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        endif
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true. )
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME, prev_movies)
        ! start watching
        prev_stacksz = 0
        nmovies      = 0
        iter         = 0
        do
            if( file_exists(trim(TERM_STREAM)) )then
                write(*,'(A)')'>>> TERMINATING PREPROCESSING STREAM'
                exit
            endif
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
                write(*,'(A,I5)')'>>> MOVIES TO PROCESS: ', stacksz
            endif
            ! completed jobs update the current project
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                ! append new processed movies to project
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                nptcls_prev = spproj%get_nptcls()
                nmovs_prev  = spproj%get_nmovies()
                do icline=1,size(completed_jobs_clines)
                    stream_spprojfile = completed_jobs_clines(icline)%get_carg('projfile')
                    call stream_spproj%read( stream_spprojfile )
                    call spproj%append_project(stream_spproj, 'mic')
                    if( l_pick )call spproj%append_project(stream_spproj, 'stk')
                    call stream_spproj%kill()
                    deallocate(stream_spprojfile)
                enddo
                nptcls = spproj%get_nptcls()
                nmovs  = spproj%get_nmovies()
                ! write
                if( nmovs == 0 )then
                    ! first write
                    call spproj%write
                else
                    ! write inside
                    call spproj%write_segment_inside('mic', fromto=[nmovs_prev+1, nmovs])
                    if( l_pick )then
                        call spproj%write_segment_inside('stk',    fromto=[nmovs_prev+1, nmovs])
                        call spproj%write_segment_inside('ptcl2D', fromto=[nptcls_prev+1, nptcls])
                        call spproj%write_segment_inside('ptcl3D', fromto=[nptcls_prev+1, nptcls])
                    endif
                endif
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
                        call simple_abspath(fname, abs_fname, errmsg='preprocess_stream :: update_projects_list 1')
                        fnames(n_old+i) = trim(abs_fname)
                        deallocate(abs_fname)
                    enddo
                else
                    ! first write
                    allocate(fnames(n_spprojs))
                    do i=1,n_spprojs
                        cline_mov = completed_jobs_clines(i)
                        fname     = trim(cline_mov%get_carg('projfile'))
                        call simple_abspath(fname, abs_fname, errmsg='preprocess_stream :: update_projects_list 2')
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
                call spproj_here%write()
                call spproj_here%kill()
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
            end subroutine create_individual_project

    end subroutine exec_preprocess_stream

    subroutine exec_cluster2D_stream_distr( self, cline )
        use simple_commander_distr_wflows, only: cluster2D_distr_commander, make_cavgs_distr_commander
        use simple_image,                  only: image
        use simple_oris,                   only: oris
        use simple_ori,                    only: ori
        class(cluster2D_stream_distr_commander), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        character(len=:), allocatable :: WORK_PROJFILE
        integer,          parameter   :: SHIFTSRCH_PTCLSLIM  = 2000  ! # of ptcls required to turn on shift search
        integer,          parameter   :: SHIFTSRCH_ITERLIM   = 10     ! # of iterations prior to turn on shift search
        integer,          parameter   :: CCRES_NPTCLS_LIM    = 10000 ! # of ptcls required to turn on objfun=ccres
        integer,          parameter   :: WAIT_WATCHER        = 60    ! seconds prior to new stack detection
        integer,          parameter   :: MAXNCLS             = 1000  ! maximum # of classes
        integer,          parameter   :: ORIGPROJ_WRITEFREQ  = 900   ! 15mins, Frequency at which the original project file should be updated
        type(parameters)                    :: params
        type(cluster2D_distr_commander)     :: xcluster2D_distr
        type(make_cavgs_distr_commander)    :: xmake_cavgs
        type(cmdline)                       :: cline_cluster2D, cline_make_cavgs
        type(sp_project)                    :: orig_proj, work_proj, stream_proj
        type(ctfparams)                     :: ctfvars
        type(oris)                          :: os_stk
        type(ori)                           :: o_stk
        character(LONGSTRLEN),  allocatable :: spproj_list(:), stk_list(:)
        character(len=:),       allocatable :: spproj_list_fname, stk
        character(len=STDLEN)               :: str_iter, refs_glob
        real    :: orig_smpd, msk, scale_factor, orig_msk, smpd
        integer :: iter, origproj_time, orig_box, box, nptcls_glob, ncls_glob, tnow, iproj, n_new_spprojs
        integer :: nptcls_glob_prev, ncls_glob_prev, last_injection, n_spprojs, n_spprojs_prev, n_new_ptcls
        logical :: do_autoscale, work_proj_has_changed, wait_for_new_ptcls, l_ccres
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
        if( .not.file_exists(params%dir_target) )then
            THROW_HARD('folder: '//trim(params%dir_target)//' does not exist!')
        endif
        call cline%set('stream','no') ! was only for parameters determination
        call cline%set('mkdir','no')
        ! init
        do_autoscale = params%autoscale.eq.'yes'
        allocate(WORK_PROJFILE, source='cluster2D_stream_tmproj.simple')
        l_ccres = .false.
        ! filepath creates spproj_list_fname with checks
        spproj_list_fname = filepath(trim(params%dir_target), trim(STREAM_SPPROJFILES))
        ! for microscopes that don't work too good
        if(.not.cline%defined('time_inactive'))params%time_inactive = 24*3600
        ! init command-lines
        cline_cluster2D  = cline
        cline_make_cavgs = cline
        call cline_cluster2D%set('prg',       'cluster2D')
        call cline_cluster2D%set('autoscale', 'no')
        call cline_cluster2D%set('extr_iter', 100.) ! no extremal randomization
        call cline_cluster2D%delete('frac_update')
        call cline_cluster2D%delete('projname')
        call cline_cluster2D%set('projfile', trim(WORK_PROJFILE))
        call cline_cluster2D%set('projname', trim(get_fbody(trim(WORK_PROJFILE),trim('simple'))))
        call cline_make_cavgs%set('prg',    'make_cavgs')
        call cline_make_cavgs%delete('autoscale')
        ! wait for the first stacks
        nptcls_glob = 0
        do
            if( file_exists(spproj_list_fname) )then
                if( .not.is_file_open(spproj_list_fname) )then
                    call read_filetable(spproj_list_fname, spproj_list)
                    nptcls_glob_prev = nptcls_glob
                    nptcls_glob      = 0
                    if( allocated(spproj_list) )then
                        ! determine number of particles
                        n_spprojs = size(spproj_list)
                        do iproj = 1,n_spprojs
                            call stream_proj%read(spproj_list(iproj))
                            nptcls_glob = nptcls_glob + stream_proj%get_nptcls()
                            call stream_proj%kill
                        enddo
                    endif
                    write(*,'(A,I8,A,A)')'>>> PARTICLES COUNT: ', nptcls_glob, ' : ',cast_time_char(simple_gettime())
                    call flush(6)
                    if( nptcls_glob > params%ncls_start * params%nptcls_per_cls )then
                        exit ! Enough particles to initiate cluster2D
                    endif
                endif
            endif
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! transfer project info, rename & rewrite
        call orig_proj%read(params%projfile)
        work_proj%projinfo = orig_proj%projinfo
        work_proj%compenv  = orig_proj%compenv
        if( orig_proj%jobproc%get_noris()>0 ) work_proj%jobproc = orig_proj%jobproc
        call work_proj%projinfo%delete_entry('projname')
        call work_proj%projinfo%delete_entry('projfile')
        call work_proj%update_projinfo(cline_cluster2D) ! name change
        call work_proj%write                            ! & write
        ! getting general parameters from the first sp_project
        call stream_proj%read(trim(spproj_list(1)))
        orig_box  = stream_proj%get_box()
        orig_smpd = stream_proj%get_smpd()
        orig_msk  = params%msk
        !params%smpd_targets2D(1) = max(orig_smpd, params%lp*LP2SMPDFAC*0.8)
        ! such that 2*msk=58 & boxmatch ~64
        params%smpd_targets2D(1) = max(orig_smpd, orig_smpd*orig_msk/29.)
        if( do_autoscale )then
            call autoscale(orig_box, orig_smpd, params%smpd_targets2D(1), box, smpd, scale_factor)
            if( box == orig_box ) do_autoscale = .false.
        endif
        if( do_autoscale )then
            msk = orig_msk * scale_factor
        else
            smpd = orig_smpd
            box  = orig_box
            msk  = orig_msk
            scale_factor = 1.
        endif
        call cline_cluster2D%set('projfile', trim(WORK_PROJFILE))
        call cline_cluster2D%set('box',      real(box))
        call cline_cluster2D%set('msk',      real(msk))
        call stream_proj%kill
        ! prep for new stacks
        allocate(stk_list(n_spprojs))
        call os_stk%new(n_spprojs) ! for original project stktab import
        do iproj=1,n_spprojs
            call stream_proj%read(spproj_list(iproj))
            o_stk           = stream_proj%os_stk%get_ori(1)
            ctfvars         = stream_proj%get_ctfparams('ptcl2D', 1)
            stk             = stream_proj%get_stkname(1)
            stk_list(iproj) = trim(stk)
            call o_stk%set_ctfvars(ctfvars)
            call os_stk%set_ori(iproj, o_stk)
        enddo
        ! updates original project
        call write_filetable('stktab.txt', stk_list)
        call orig_proj%add_stktab('stktab.txt', os_stk)
        call orig_proj%write
        call orig_proj%kill
        ! scale & updates stack list
        call scale_stks( stk_list )
        ! updates current project
        call write_filetable('stktab.txt', stk_list)
        call os_stk%set_all2single('smpd', smpd)
        call work_proj%add_stktab('stktab.txt', os_stk)
        call work_proj%write
        call del_file('stktab.txt')
        ! MAIN LOOP
        nptcls_glob    = work_proj%get_nptcls()
        ncls_glob      = nint(real(nptcls_glob) / real(params%nptcls_per_cls))
        last_injection = simple_gettime()
        origproj_time  = last_injection
        do iter = 1, 999
            str_iter  = int2str_pad(iter,3)
            ! time handling
            if( is_timeout(simple_gettime()) )exit
            ! CLUSTER2D
            call cline_cluster2D%delete('endit')
            call cline_cluster2D%set('startit', real(iter))
            call cline_cluster2D%set('maxits',  real(iter))
            call cline_cluster2D%set('ncls',    real(ncls_glob))
            call cline_cluster2D%set('nparts',  real(min(ncls_glob,params%nparts)))
            call cline_cluster2D%set('stream',  'yes') ! has to be updated at each iteration
            if( iter > 1 ) call cline_cluster2D%set('refs', trim(refs_glob))
            ! shift limit
            if( nptcls_glob>SHIFTSRCH_PTCLSLIM .and. iter>SHIFTSRCH_ITERLIM )call cline_cluster2D%set('trs', MINSHIFT)
            ! objective function
            if( .not.l_ccres )then
                if( nptcls_glob > CCRES_NPTCLS_LIM .and. iter > SHIFTSRCH_ITERLIM )then
                    l_ccres = .true.
                    call cline_cluster2D%set('objfun', 'ccres')
                    write(*,'(A)')'>>> SWITCHED TO CCRES OBJECTIVE FUNCTION'
                endif
            endif
            ! execute
            params_glob%nptcls = nptcls_glob
            call xcluster2D_distr%execute(cline_cluster2D)
            wait_for_new_ptcls    = .false.
            if( cline_cluster2D%defined('converged') )then
                wait_for_new_ptcls = (cline_cluster2D%get_carg('converged').eq.'yes') .and. is_even(iter)
                call cline_cluster2D%delete('converged')
            endif
            ! update
            call work_proj%kill
            call work_proj%read(trim(WORK_PROJFILE))
            ! current references file name
            refs_glob = trim(CAVGS_ITER_FBODY)//trim(str_iter)//trim(params%ext)
            ! remap zero-population classes
            work_proj_has_changed = .false.
            call remap_empty_classes
            if( work_proj_has_changed )call work_proj%write_segment_inside('ptcl2D',trim(WORK_PROJFILE))
            ! termination and/or pause
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                write(*,'(A,A)')'>>> CLUSTER2D STREAM PAUSED ',cast_time_char(simple_gettime())
                call simple_sleep(WAIT_WATCHER)
            enddo
            if( file_exists(trim(TERM_STREAM)) )then
                write(*,'(A,A)')'>>> TERMINATING CLUSTER2D STREAM ',cast_time_char(simple_gettime())
                exit
            endif
            ! wait for new images
            tnow = simple_gettime()
            if( wait_for_new_ptcls )write(*,'(A,A)')'>>> WAITING FOR NEW PARTICLES ',cast_time_char(tnow)
            do while( wait_for_new_ptcls )
                if( is_timeout(tnow) )exit
                call simple_sleep(WAIT_WATCHER)
                if( is_file_open(spproj_list_fname) )cycle
                if( nlines(spproj_list_fname) > n_spprojs )wait_for_new_ptcls = .false.
            enddo
            if( is_timeout(tnow) )exit
            ! handles whether new individual project files have appeared
            if( .not.is_file_open(spproj_list_fname) )then
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
                    do iproj=n_spprojs_prev+1,n_spprojs
                        call stream_proj%read(spproj_list(iproj))
                        n_new_ptcls = n_new_ptcls + stream_proj%get_nptcls()
                        stk = stream_proj%get_stkname(1)
                        stk_list(iproj-n_spprojs_prev) = trim(stk)
                    enddo
                    if( n_new_ptcls > 0 )then
                        ! scale new stacks
                        call scale_stks( stk_list )
                        ! updates counters
                        nptcls_glob_prev = nptcls_glob
                        nptcls_glob      = nptcls_glob + n_new_ptcls
                        ncls_glob_prev   = ncls_glob
                        ncls_glob        = min(MAXNCLS, nint(real(nptcls_glob)/real(params%nptcls_per_cls)))
                        ! remap in case of increased number of classes using not updated project!
                        call map_new_ptcls
                        ! then only update project with new images
                        do iproj=n_spprojs_prev+1,n_spprojs
                            call stream_proj%read(spproj_list(iproj))
                            ctfvars = stream_proj%get_ctfparams('ptcl2D', 1)
                            if( do_autoscale )ctfvars%smpd = ctfvars%smpd / scale_factor
                            call work_proj%add_stk(stk_list(iproj-n_spprojs_prev), ctfvars)
                        enddo
                        call work_proj%write  !!!
                        ! logging
                        write(*,'(A,I8,A,A)')'>>> NEW PARTICLES COUNT: ', nptcls_glob, ' ; ',cast_time_char(simple_gettime())
                        last_injection = simple_gettime()
                    endif
                endif
            endif
            ! update original project
            if( simple_gettime()-origproj_time > ORIGPROJ_WRITEFREQ )then
                write(*,'(A,A)')'>>> UPDATING PROJECT FILE ', trim(params%projfile)
                call update_orig_proj
                origproj_time = simple_gettime()
            endif
            ! wait
            call simple_sleep(WAIT_WATCHER)
        enddo
        ! cleanup
        call qsys_cleanup
        ! updates original project
        n_spprojs_prev = n_spprojs
        call update_orig_proj
        ! class averages at original sampling
        if ( do_autoscale )then
            call cline_make_cavgs%set('ncls', real(ncls_glob))
            call cline_make_cavgs%set('refs', refs_glob)
            call xmake_cavgs%execute(cline_make_cavgs) ! need be distributed
        endif
        call orig_proj%add_cavgs2os_out(refs_glob, orig_smpd)
        call orig_proj%write_segment_inside('out',fromto=[1,1])
        ! cleanup
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_DISTR_CLUSTER2D_STREAM NORMAL STOP ****')

        contains

            logical function is_timeout( time_now )
                integer, intent(in) :: time_now
                is_timeout = .false.
                if(time_now-last_injection > params%time_inactive)then
                    write(*,*)'>>> TIME LIMIT WITHOUT NEW IMAGES REACHED: ',cast_time_char(time_now)
                    is_timeout = .true.
                else if(time_now-last_injection > 3600)then
                    write(*,*)'>>> OVER ONE HOUR WITHOUT NEW PARTICLES: ',cast_time_char(time_now)
                    call flush(6)
                endif
                return
            end function is_timeout

            subroutine update_orig_proj
                ! assumes work_proj is to correct dimension
                call orig_proj%read(params%projfile)
                do iproj=n_spprojs_prev+1,n_spprojs
                    call stream_proj%read(spproj_list(iproj))
                    stk     = stream_proj%get_stkname(1)
                    ctfvars = stream_proj%get_ctfparams('ptcl2D', 1)
                    call orig_proj%add_stk(stk, ctfvars)
                enddo
                ! updates 2D & wipes 3D segment
                if( do_autoscale )call work_proj%os_ptcl2D%mul_shifts( 1./scale_factor )
                orig_proj%os_ptcl2D = work_proj%os_ptcl2D
                orig_proj%os_ptcl3D = work_proj%os_ptcl3D ! to wipe the 3d segment
                call orig_proj%add_cavgs2os_out(refs_glob, orig_smpd)
                call orig_proj%write()
            end subroutine update_orig_proj

            subroutine scale_stks( stk_fnames )
                character(len=*), allocatable, intent(inout) :: stk_fnames(:)
                character(len=*), parameter :: SCALE_FILETAB = 'stkscale.txt'
                character(len=*), parameter :: SCALE_DIR     = './scaled_stks/'
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

            !>  empty classes re-mapping, commits the project to disk
            subroutine remap_empty_classes
                integer, allocatable  :: fromtocls(:,:)
                type(image)           :: img_cavg
                character(len=STDLEN) :: stk
                real                  :: smpd
                integer               :: icls
                call work_proj%os_ptcl2D%fill_empty_classes(ncls_glob, fromtocls)
                if( allocated(fromtocls) )then
                    ! updates document later
                    work_proj_has_changed = .true.
                    smpd = work_proj%get_smpd()
                    ! updates classes
                    call img_cavg%new([box, box,1], smpd)
                    do icls = 1, size(fromtocls, dim=1)
                        ! cavg
                        call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                        call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                        ! even & odd
                        stk = add2fbody(trim(refs_glob),params%ext,'_even')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                        stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                    enddo
                    ! stack size preservation
                    call img_cavg%read(trim(refs_glob), ncls_glob)
                    call img_cavg%write(trim(refs_glob), ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_even')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    ! cleanup
                    call img_cavg%kill
                    deallocate(fromtocls)
                endif
            end subroutine remap_empty_classes

            subroutine map_new_ptcls
                use simple_projection_frcs, only: projection_frcs
                type(projection_frcs) :: frcs_prev, frcs
                type(image)           :: img_cavg
                integer, allocatable  :: fromtocls(:,:)
                real                  :: smpd
                integer               :: icls, state
                character(len=STDLEN) :: stk
                state = 1
                if( ncls_glob.eq.ncls_glob_prev )return
                ! updates references & FRCs
                write(*,'(A,I8)')'>>> NEW CLASSES   COUNT: ', ncls_glob
                call work_proj%os_ptcl2D%fill_empty_classes(ncls_glob, fromtocls)
                if( allocated(fromtocls) )then
                    ! references
                    smpd = work_proj%get_smpd()
                    call img_cavg%new([box,box,1], smpd)
                    do icls = 1, size(fromtocls, dim=1)
                        ! cavg
                        call img_cavg%read(trim(refs_glob), fromtocls(icls, 1))
                        call img_cavg%write(trim(refs_glob), fromtocls(icls, 2))
                        ! even & odd
                        stk = add2fbody(trim(refs_glob),params%ext,'_even')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                        stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                        call img_cavg%read(stk, fromtocls(icls, 1))
                        call img_cavg%write(stk, fromtocls(icls, 2))
                    enddo
                    ! stack size preservation
                    call img_cavg%read(trim(refs_glob), ncls_glob)
                    call img_cavg%write(trim(refs_glob), ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_even')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    stk = add2fbody(trim(refs_glob),params%ext,'_odd')
                    call img_cavg%read(stk, ncls_glob)
                    call img_cavg%write(stk, ncls_glob)
                    ! FRCs
                    call frcs_prev%new(ncls_glob_prev, box, smpd, state)
                    call frcs%new(ncls_glob, box, smpd, state)
                    call frcs_prev%read(FRCS_FILE)
                    do icls = 1, ncls_glob_prev
                        call frcs%set_frc(icls,&
                        &frcs_prev%get_frc(icls, box, state), state)
                    enddo
                    do icls = 1, size(fromtocls, dim=1)
                        call frcs%set_frc( fromtocls(icls,2),&
                        &frcs%get_frc(fromtocls(icls,1), box, state), state)
                    enddo
                    call frcs%write(FRCS_FILE)
                    call frcs%kill
                    call frcs_prev%kill
                    deallocate(fromtocls)
                    call img_cavg%kill
                endif
            end subroutine map_new_ptcls

    end subroutine exec_cluster2D_stream_distr

    subroutine exec_pick_extract_stream_distr( self, cline )
        class(pick_extract_stream_distr_commander), intent(inout) :: self
        class(cmdline),                             intent(inout) :: cline
        integer,          parameter   :: WAIT_WATCHER        = 10    ! seconds prior to new stack detection
        !integer,          parameter   :: ORIGPROJ_WRITEFREQ  = 600   ! 10mins, Frequency at which the original project file should be updated
        type(parameters)                    :: params
        type(cmdline)                       :: cline_pick_extract
        type(sp_project)                    :: orig_proj, stream_proj
        type(qsys_env)                      :: qenv
        class(cmdline),         allocatable :: completed_jobs_clines(:)
        character(LONGSTRLEN),  allocatable :: spproj_list(:)
        character(len=:),       allocatable :: spproj_list_fname, output_dir, output_dir_picker, output_dir_extract, projfname
        integer :: iter, origproj_time, tnow, iproj, icline, nptcls, prev_stacksz, stacksz
        integer :: last_injection, n_spprojs, n_spprojs_prev, n_newspprojs, nmics
        ! output command line executed
        write(*,'(a)') '>>> COMMAND LINE EXECUTED'
        write(*,*) trim(cmdline_glob)
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
        cline_pick_extract = cline
        call cline_pick_extract%set('prg', 'pick_extract')
        call cline_pick_extract%delete('projname')
        ! wait for the first stacks
        last_injection  = simple_gettime()
        origproj_time   = last_injection
        prev_stacksz    = 0
        n_spprojs       = 0
        n_spprojs_prev  = 0
        do iter = 1,999999
            tnow = simple_gettime()
            if(tnow-last_injection > params%time_inactive)then
                write(*,*)'>>> TIME LIMIT WITHOUT NEW MICROGRAPHS REACHED'
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
                        write(*,'(A,I5)')'>>> MICROGRAPHS TO PROCESS: ', stacksz
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
                        write(*,'(A,I8)')'>>> NEW MICROGRAPHS COUNT: ', nmics
                        write(*,'(A,I8)')'>>> NEW PARTICLES   COUNT: ', nptcls
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
