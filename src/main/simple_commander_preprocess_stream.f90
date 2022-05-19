! concrete commander: pre-processing routines
module simple_commander_preprocess_stream
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,                        only: cmdline
use simple_parameters,                     only: parameters, params_glob
use simple_commander_base,                 only: commander_base
use simple_image,                          only: image
use simple_ori,                            only: ori
use simple_oris,                           only: oris
use simple_sp_project,                     only: sp_project
use simple_qsys_env,                       only: qsys_env
use simple_stack_io,                       only: stack_io
use simple_starproject,                    only: starproject
use simple_commander_cluster2D_stream_dev
use simple_qsys_funs
use simple_commander_preprocess
implicit none

public :: preprocess_commander_stream_dev

private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: dir_preprocess = trim(PATH_HERE)//'spprojs/'

type, extends(commander_base) :: preprocess_commander_stream_dev
  contains
    procedure :: execute      => exec_preprocess_stream_dev
end type preprocess_commander_stream_dev

contains

    subroutine exec_preprocess_stream_dev( self, cline )
        use simple_moviewatcher, only: moviewatcher
        use simple_timer
        class(preprocess_commander_stream_dev), intent(inout) :: self
        class(cmdline),                         intent(inout) :: cline
        type(parameters)                       :: params
        integer,                   parameter   :: WAITTIME        = 5   ! folder watched every 5 seconds
        integer,                   parameter   :: LONGTIME        = 300  ! time lag after which a movie is processed
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .false.
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        class(cmdline),            allocatable :: completed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj
        type(starproject)                      :: starproj
        character(len=LONGSTRLEN), allocatable :: movies(:)
        character(len=LONGSTRLEN), allocatable, target :: completed_fnames(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract
        character(len=LONGSTRLEN)              :: movie
        real                                   :: pickref_scale
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection
        integer                                :: cnt, n_imported, nptcls_glob
        logical                                :: l_pick, l_movies_left, l_haschanged, l_cluster2d
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        if( .not. cline%defined('oritype')         ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('mkdir')           ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('walltime')        ) call cline%set('walltime',  29.0*60.0) ! 29 minutes
        ! motion correction
        if( .not. cline%defined('trs')             ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')         ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')          ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')            ) call cline%set('bfac',             50.)
        call cline%set('groupframes',     'no')
        if( .not. cline%defined('mcconvention')    ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_fraction')    ) call cline%set('eer_fraction',     20.)
        if( .not. cline%defined('eer_upsampling')  ) call cline%set('eer_upsampling',    1.)
        if( .not. cline%defined('algorithm')       ) call cline%set('algorithm',    'patch')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')         ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate') ) call cline%set('hp_ctf_estimate',  30.)
        if( .not. cline%defined('lp_ctf_estimate') ) call cline%set('lp_ctf_estimate',   5.)
        if( .not. cline%defined('dfmin')           ) call cline%set('dfmin',            0.3)
        if( .not. cline%defined('dfmax')           ) call cline%set('dfmax',            5.0)
        if( .not. cline%defined('ctfpatch')        ) call cline%set('ctfpatch',       'yes')
        ! picking
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('ndev')            ) call cline%set('ndev',              2.)
        if( .not. cline%defined('thres')           ) call cline%set('thres',            24.)
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( cline%defined('refs') .and. cline%defined('vol1') )then
            THROW_HARD('REFS and VOL1 cannot be both provided!')
        endif
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        ! 2D classification
        if( .not. cline%defined('lpthresh')  ) call cline%set('lpthresh',      30.0)
        if( .not. cline%defined('ndev2D')    ) call cline%set('ndev2D',         1.5)
        if( .not. cline%defined('wiener')    ) call cline%set('wiener',   'partial')
        if( .not. cline%defined('autoscale') ) call cline%set('autoscale',    'yes')
        ! master parameters
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        if( cline%defined('dir_prev') .and. .not.file_exists(params%dir_prev) )then
            THROW_HARD('Directory '//trim(params%dir_prev)//' does not exist!')
        endif
        ! master project file
        call spproj%read( params%projfile )
        call spproj%update_projinfo(cline)
        if( spproj%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        call spproj%write(micspproj_fname)
        ! picking
        l_pick = .false.
        if( cline%defined('refs') .or. cline%defined('vol1') ) l_pick = .true.
        ! output directories
        output_dir = trim(PATH_HERE)//trim(dir_preprocess)
        call simple_mkdir(output_dir)
        output_dir_ctf_estimate   = filepath(trim(PATH_HERE), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(PATH_HERE), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        if( l_pick )then
            output_dir_picker  = filepath(trim(PATH_HERE), trim(DIR_PICKER))
            output_dir_extract = filepath(trim(PATH_HERE), trim(DIR_EXTRACT))
            call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
            call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        endif
        call cline%set('dir','../')
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! prepares picking references
        if( l_pick )then
            cline_make_pickrefs = cline
            call cline_make_pickrefs%set('prg','make_pickrefs')
            call cline_make_pickrefs%set('stream','no')
            call cline_make_pickrefs%delete('ncls')
            call cline_make_pickrefs%delete('mskdiam')
            if( cline_make_pickrefs%defined('eer_upsampling') )then
                pickref_scale = real(params%eer_upsampling) * params%scale
                call cline_make_pickrefs%set('scale',pickref_scale)
            endif
            call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
            call cline%set('refs', '../'//trim(PICKREFS)//trim(params%ext))
            call cline%delete('vol1')
            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
        endif
        ! prep for 2D classification
        l_cluster2D = .false.
        if( l_pick )then
            call init_cluster2D_stream(cline, spproj, trim(PICKREFS)//trim(params%ext), micspproj_fname, l_cluster2D)
        endif
        call cline%delete('ncls')
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME)
        ! import previous runs
        nptcls_glob = 0
        call import_prev_streams
        ! start watching
        last_injection = simple_gettime()
        prev_stacksz  = 0
        nmovies       = 0
        iter          = 0
        n_imported    = 0
        l_movies_left = .false.
        l_haschanged  = .false.
        do
            ! termination & pausing
            if( file_exists(trim(TERM_STREAM)) )then
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            do while( file_exists(trim(PAUSE_STREAM)) )
                if( file_exists(trim(TERM_STREAM)) ) exit
                call write_singlelineoftext(PAUSE_STREAM, 'PAUSED')
                write(logfhandle,'(A,A)')'>>> PREPROCESS STREAM PAUSED ',cast_time_char(simple_gettime())
                call sleep(WAITTIME)
            enddo
            iter = iter + 1
            call movie_buff%watch( nmovies, movies )
            ! append movies to processing stack
            if( nmovies > 0 )then
                cnt = 0
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    if( movie_buff%is_past(movie) )cycle
                    call create_individual_project
                    call qenv%qscripts%add_to_streaming( cline )
                    call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
                    call movie_buff%add2history( movies(imovie) )
                    cnt = cnt+1
                    if( cnt == min(params%nparts,nmovies) ) exit
                enddo
                l_movies_left = cnt .ne. nmovies
            else
                l_movies_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I5)')'>>> MOVIES TO PROCESS:                ', stacksz
            endif
            ! fetch completed jobs list & updates of cluster2D_stream
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                call update_projects_list( completed_fnames, n_imported )
                do cnt = 1,size(completed_jobs_clines)
                    call completed_jobs_clines(cnt)%kill
                enddo
                deallocate(completed_jobs_clines)
            else
                n_imported = 0 ! newly imported
            endif
            ! project update
            if( n_imported > 0 )then
                n_imported = spproj%os_mic%get_noris()
                write(logfhandle,'(A,I5)')'>>> # MOVIES PROCESSED & IMPORTED:    ',n_imported
                if( l_pick ) write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:         ',nptcls_glob
                ! write project
                call spproj%write_segment_inside('mic',micspproj_fname)
                ! write starfile snapshot
                if (spproj%os_mic%get_noris() > 0) then
                    if( file_exists("micrographs.star") ) call del_file("micrographs.star")
                    call starproj%assign_optics(cline, spproj)
                    call starproj%export_mics(cline, spproj)
                end if
                last_injection = simple_gettime()
                l_haschanged   = .true.
                n_imported     = spproj%os_mic%get_noris()
            else
                ! wait & write snapshot
                if( l_cluster2D )then
                    ! cf 2D section
                else
                    if( .not.l_movies_left )then
                        if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                            ! write project when inactive...
                            call write_project
                            l_haschanged = .false.
                        else
                            ! ...or wait
                            call sleep(WAITTIME)
                        endif
                    endif
                endif
            endif
            ! 2D classification section
            if( l_cluster2D )then
                call update_chunks
                call update_pool_status
                call update_pool
                call reject_from_pool
                call write_project_stream2D(.true.,.false.)
                call import_chunks_into_pool
                call classify_pool
                call update_projects_mask(completed_fnames)
                call classify_new_chunks(completed_fnames)
            endif
        end do
        ! termination
        if( l_cluster2D )then
            call terminate_stream2D
        else
            call write_project
        endif
        call spproj%kill
        ! cleanup
        call qsys_cleanup
        call del_file(micspproj_fname)
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            subroutine write_project()
                logical, allocatable :: stk_mask(:)
                integer, allocatable :: states(:)
                integer              :: iproj,nptcls,istk,fromp,top,i,iptcl,nstks,n,nmics
                write(logfhandle,'(A)')'>>> PROJECT UPDATE'
                nmics = spproj%os_mic%get_noris()
                call spproj%write_segment_inside('mic', params%projfile)
                if( l_pick )then
                    if( DEBUG_HERE ) t0 = tic()
                    ! stacks
                    allocate(stk_mask(nmics))
                    allocate(states(nmics))
                    do iproj = 1,nmics
                        stk_mask(iproj) = nint(spproj%os_mic%get(iproj,'nptcls')) > 0
                        states(iproj)   = spproj%os_mic%get_state(iproj)
                    enddo
                    nstks = count(stk_mask)
                    call spproj%os_stk%new(nstks, is_ptcl=.false.)
                    nptcls = 0
                    istk   = 0
                    fromp  = 0
                    top    = 0
                    do iproj = 1,nmics
                        if( .not.stk_mask(iproj) ) cycle
                        istk = istk+1
                        call stream_spproj%read_segment('stk', completed_fnames(iproj))
                        call stream_spproj%os_stk%set(1, 'state', real(states(iproj)))
                        n      = nint(stream_spproj%os_stk%get(1,'nptcls'))
                        fromp  = nptcls + 1
                        nptcls = nptcls + n
                        top    = nptcls
                        call spproj%os_stk%transfer_ori(istk,stream_spproj%os_stk,1)
                        call spproj%os_stk%set(istk, 'fromp',real(fromp))
                        call spproj%os_stk%set(istk, 'top',  real(top))
                    enddo
                    call spproj%write_segment_inside('stk', params%projfile)
                    call spproj%os_stk%reset
                    call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                    call spproj%os_ptcl3D%new(nptcls, is_ptcl=.true.)
                    ! particles 2D
                    istk   = 0
                    iptcl  = 0
                    do iproj = 1,nmics
                        if( .not.stk_mask(iproj) ) cycle
                        istk = istk+1
                        call stream_spproj%read_segment('ptcl2D', completed_fnames(iproj))
                        nptcls = stream_spproj%os_ptcl2D%get_noris()
                        do i = 1,nptcls
                            iptcl = iptcl + 1
                            call spproj%os_ptcl2D%transfer_ori(iptcl,stream_spproj%os_ptcl2D,i)
                            call spproj%os_ptcl2D%set(iptcl, 'stkind', real(istk))
                            call spproj%os_ptcl2D%set(iptcl, 'state',  real(states(iproj)))
                        enddo
                        call stream_spproj%kill
                    enddo
                    call spproj%write_segment_inside('ptcl2D', params%projfile)
                    call spproj%os_ptcl2D%reset
                    ! particles 3D
                    istk   = 0
                    iptcl  = 0
                    do iproj = 1,nmics
                        if( .not.stk_mask(iproj) ) cycle
                        istk = istk+1
                        call stream_spproj%read_segment('ptcl3D', completed_fnames(iproj))
                        nptcls = stream_spproj%os_ptcl3D%get_noris()
                        do i = 1,nptcls
                            iptcl = iptcl + 1
                            call spproj%os_ptcl3D%transfer_ori(iptcl,stream_spproj%os_ptcl3D,i)
                            call spproj%os_ptcl3D%set(iptcl, 'stkind', real(istk))
                            call spproj%os_ptcl3D%set(iptcl, 'state',  real(states(iproj)))
                        enddo
                        call stream_spproj%kill
                    enddo
                    call spproj%write_segment_inside('ptcl3D', params%projfile)
                    write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:         ',spproj%os_ptcl2D%get_noris()
                endif
                call spproj%write_non_data_segments(params%projfile)
                ! cleanup, we preserve os_mic
                call spproj%os_stk%kill
                call spproj%os_ptcl2D%kill
                call spproj%os_ptcl3D%kill
                ! benchmark
                if( DEBUG_HERE )then
                    rt_write = toc(t0)
                    print *,'rt_write  : ', rt_write; call flush(6)
                endif
            end subroutine write_project

            ! returns list of completed jobs + updates for cluster2D_stream
            subroutine update_projects_list( completedfnames, nimported )
                character(len=LONGSTRLEN), allocatable, intent(inout) :: completedfnames(:)
                integer,                                intent(inout) :: nimported
                type(sp_project)                       :: streamspproj
                character(len=:),          allocatable :: fname, abs_fname
                character(len=LONGSTRLEN), allocatable :: old_fnames(:)
                logical, allocatable :: spproj_mask(:)
                integer :: i, n_spprojs, n_old, nptcls_here, state, j, n2import, nprev_imports, n_completed, nptcls
                n_completed = 0
                nimported   = 0
                ! projects to import
                n_spprojs = size(completed_jobs_clines)
                if( n_spprojs == 0 )return
                allocate(spproj_mask(n_spprojs),source=.true.)
                ! previously completed projects
                n_old = 0 ! on first import
                if( allocated(completedfnames) ) n_old = size(completed_fnames)
                if( l_pick )then
                    do i = 1,n_spprojs
                        ! flags zero-picked mics that will not be imported
                        fname = trim(output_dir)//trim(completed_jobs_clines(i)%get_carg('projfile'))
                        call check_nptcls(fname,nptcls_here,state)
                        spproj_mask(i) = (nptcls_here > 0) .and. (state > 0)
                    enddo
                endif
                n2import = count(spproj_mask)
                if( l_pick )then 
                    if( n2import /= n_spprojs )then
                        write(logfhandle,'(A,I3,A)')'>>> NO PARTICLES FOUND IN ',n_spprojs-n2import,' MICROGRAPHS'
                    endif
                endif
                if( n2import > 0 )then
                    n_completed = n_old + n2import
                    nimported   = n2import
                    nprev_imports = spproj%os_mic%get_noris() ! should be equal to n_old ?
                    if( nprev_imports == 0 )then
                        call spproj%os_mic%new(n2import, is_ptcl=.false.) ! first time
                        allocate(completedfnames(n_completed))
                    else
                        call spproj%os_mic%reallocate(n_completed)
                        old_fnames = completed_fnames(:)
                        deallocate(completed_fnames)
                        allocate(completedfnames(n_completed))
                        if( n_old > 0 )then
                            completedfnames(1:n_old) = old_fnames(:)
                        endif
                        deallocate(old_fnames)
                    endif
                    j = 0
                    nptcls_here = 0
                    do i=1,n_spprojs
                        if( .not.spproj_mask(i) ) cycle
                        j = j+1
                        fname     = trim(output_dir)//trim(completed_jobs_clines(i)%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
                        completedfnames(n_old+j) = trim(abs_fname)
                        call streamspproj%read_segment('mic', abs_fname)
                        if( l_pick )then
                            nptcls = nint(streamspproj%os_mic%get(1,'nptcls'))
                            if( nptcls > 0 )then
                                nptcls_here = nptcls_here + nptcls
                                call streamspproj%read_segment('stk', abs_fname)
                                call update_ori_path(streamspproj%os_stk, 1, 'stk')
                                call streamspproj%write_segment_inside('stk', abs_fname)
                            endif
                        endif
                        call spproj%os_mic%transfer_ori(n_old+j, streamspproj%os_mic, 1)
                        call streamspproj%kill
                        ! update paths such that relative paths are with respect to root folder
                        call update_ori_path(spproj%os_mic, n_old+j, 'starfile')
                        call update_ori_path(spproj%os_mic, n_old+j, 'intg')
                        call update_ori_path(spproj%os_mic, n_old+j, 'forctf')
                        call update_ori_path(spproj%os_mic, n_old+j, 'thumb')
                        call update_ori_path(spproj%os_mic, n_old+j, 'mceps')
                        call update_ori_path(spproj%os_mic, n_old+j, 'ctfeps')
                        call update_ori_path(spproj%os_mic, n_old+j, 'ctfjpg')
                        if( l_pick ) call update_ori_path(spproj%os_mic, n_old+j, 'boxfile')
                    enddo
                    if( l_pick ) nptcls_glob = nptcls_glob + nptcls_here
                else
                    nimported  = 0
                    return
                endif
            end subroutine update_projects_list

            subroutine update_ori_path(os, i, key)
                class(oris),      intent(inout) :: os
                integer,          intent(in)    :: i
                character(len=*), intent(in)    :: key
                character(len=:), allocatable :: fname
                character(len=LONGSTRLEN)     :: newfname
                if( os%isthere(i,key) )then
                    call os%getter(i,key,fname)
                    if( len_trim(fname) > 3 )then
                        if( fname(1:3) == '../' ) fname = trim(fname(4:))
                    endif
                    call make_relativepath(CWD_GLOB, fname, newfname)
                    call os%set(i, key, newfname)
                endif
            end subroutine update_ori_path

            subroutine check_nptcls( fname, nptcls, state )
                character(len=*), intent(in)  :: fname
                integer,          intent(out) :: nptcls, state
                type(sp_project) :: spproj_here
                call spproj_here%read_segment('mic',fname)
                nptcls = nint(spproj_here%os_mic%get(1,'nptcls'))
                state  = nint(spproj_here%os_mic%get(1,'state'))
                call spproj_here%kill
            end subroutine check_nptcls

            subroutine create_individual_project
                type(sp_project)              :: spproj_here
                type(cmdline)                 :: cline_here
                type(ctfparams)               :: ctfvars
                character(len=STDLEN)         :: ext, movie_here
                character(len=LONGSTRLEN)     :: projname, projfile
                character(len=XLONGSTRLEN)    :: cwd, cwd_old
                cwd_old = trim(cwd_glob)
                call chdir(output_dir)
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                movie_here = basename(trim(movie))
                ext        = fname2ext(trim(movie_here))
                projname   = trim(PREPROCESS_PREFIX)//trim(get_fbody(trim(movie_here), trim(ext)))
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
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
                call chdir(cwd_old)
                cwd_glob = trim(cwd_old)
                call spproj_here%kill
                call cline_here%kill
            end subroutine create_individual_project

            !>  import previous run to the current project based on past single project files
            ! TO DOUBLE-CHECK
            subroutine import_prev_streams
                use simple_ori, only: ori
                type(sp_project) :: streamspproj
                type(ori)        :: o, o_stk
                character(len=LONGSTRLEN), allocatable :: sp_files(:)
                character(len=:), allocatable :: mic, mov
                logical,          allocatable :: spproj_mask(:)
                integer :: iproj,nprojs,icnt,nptcls
                logical :: err
                if( .not.cline%defined('dir_prev') ) return
                err = .false.
                call simple_list_files_regexp(params%dir_prev,'^'//trim(PREPROCESS_PREFIX)//'.*\.simple$',sp_files)
                nprojs = size(sp_files)
                if( nprojs < 1 ) return
                allocate(spproj_mask(nprojs),source=.false.)
                nptcls = 0
                do iproj = 1,nprojs
                    call streamspproj%read_segment('mic', sp_files(iproj) )
                    if( streamspproj%os_mic%get_noris() /= 1 )then
                        THROW_WARN('Ignoring previous project'//trim(sp_files(iproj)))
                        cycle
                    endif
                    if( .not. streamspproj%os_mic%isthere(1,'intg') )cycle
                    if( l_pick )then
                        if( streamspproj%os_mic%get(1,'nptcls') < 0.5 )cycle
                    endif
                    spproj_mask(iproj) = .true.
                enddo
                if( count(spproj_mask) == 0 )then
                    nptcls_glob = 0
                    return
                endif
                icnt = 0
                do iproj = 1,nprojs
                    if( .not.spproj_mask(iproj) )cycle
                    call streamspproj%read_segment('mic',sp_files(iproj))
                    call streamspproj%os_mic%get_ori(1, o)
                    ! import mic segment
                    call movefile2folder('intg',        output_dir_motion_correct, o, err)
                    call movefile2folder('forctf',      output_dir_motion_correct, o, err)
                    call movefile2folder('thumb',       output_dir_motion_correct, o, err)
                    call movefile2folder('mc_starfile', output_dir_motion_correct, o, err)
                    call movefile2folder('mceps',       output_dir_motion_correct, o, err)
                    call movefile2folder('ctfjpg',      output_dir_ctf_estimate,   o, err)
                    call movefile2folder('ctfdoc',      output_dir_ctf_estimate,   o, err)
                    call movefile2folder('ctfeps',      output_dir_ctf_estimate,   o, err)
                    if( l_pick )then
                        ! import mic & updates stk segment
                        call movefile2folder('boxfile', output_dir_picker, o, err)
                        nptcls = nptcls + nint(o%get('nptcls'))
                        call streamspproj%os_mic%set_ori(1, o)
                        if( .not.err )then
                            call streamspproj%read_segment('stk', sp_files(iproj))
                            if( streamspproj%os_stk%get_noris() == 1 )then
                                call streamspproj%os_stk%get_ori(1, o_stk)
                                call movefile2folder('stk', output_dir_extract, o_stk, err)
                                call streamspproj%os_stk%set_ori(1, o_stk)
                                call streamspproj%read_segment('ptcl2D', sp_files(iproj))
                                call streamspproj%read_segment('ptcl3D', sp_files(iproj))
                            endif
                        endif
                    else
                        ! import mic segment
                        call streamspproj%os_mic%set_ori(1, o)
                    endif
                    ! add to history
                    call o%getter('movie', mov)
                    call o%getter('intg', mic)
                    call movie_buff%add2history(mov)
                    call movie_buff%add2history(mic)
                    ! write updated individual project file
                    call streamspproj%write(basename(sp_files(iproj)))
                    ! count
                    icnt = icnt + 1
                enddo
                if( icnt > 0 )then
                    ! updating STREAM_SPPROJFILES for Cluster2D_stream
                    allocate(completed_jobs_clines(icnt))
                    icnt = 0
                    do iproj = 1,nprojs
                        if(spproj_mask(iproj))then
                            icnt = icnt+1
                            call completed_jobs_clines(icnt)%set('projfile',basename(sp_files(iproj)))
                            call completed_jobs_clines(icnt)%printline
                        endif
                    enddo
                    call update_projects_list(completed_fnames, n_imported)
                    deallocate(completed_jobs_clines)
                endif
                call o%kill
                call o_stk%kill
                call streamspproj%kill
                write(*,'(A,I3)')'>>> IMPORTED PREVIOUS PROCESSED MOVIES: ', icnt
            end subroutine import_prev_streams

            subroutine movefile2folder(key, folder, o, err)
                character(len=*), intent(in)    :: key, folder
                class(ori),       intent(inout) :: o
                logical,          intent(out)   :: err
                character(len=:), allocatable :: src
                character(len=LONGSTRLEN) :: dest,reldest
                integer :: iostat
                err = .false.
                if( .not.o%isthere(key) )then
                    err = .true.
                    return
                endif
                call o%getter(key,src)
                if( .not.file_exists(src) )then
                    err = .true.
                    return
                endif
                dest   = trim(folder)//'/'//basename(src)
                iostat = rename(src,dest)
                if( iostat /= 0 )then
                    THROW_WARN('Ignoring '//trim(src))
                    return
                endif
                iostat = rename(src,reldest)
                call make_relativepath(CWD_GLOB,dest,reldest)
                call o%set(key,reldest)
            end subroutine movefile2folder

    end subroutine exec_preprocess_stream_dev

end module simple_commander_preprocess_stream
