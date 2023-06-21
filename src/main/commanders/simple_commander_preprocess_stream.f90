! concrete commander: pre-processing routines
module simple_commander_preprocess_stream
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,                        only: cmdline
use simple_parameters,                     only: parameters, params_glob
use simple_commander_base,                 only: commander_base
use simple_image,                          only: image
use simple_sp_project,                     only: sp_project
use simple_qsys_env,                       only: qsys_env
use simple_stack_io,                       only: stack_io
use simple_starproject,                    only: starproject
use simple_commander_cluster2D_stream_dev
use simple_qsys_funs
use simple_commander_preprocess
use simple_progress
use FoX_dom
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
        integer,                   parameter   :: WAITTIME        = 3    ! movie folder watched every WAITTIME seconds
        integer,                   parameter   :: LONGTIME        = 300  ! time lag after which a movie is processed
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .false.
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj, tmp_spproj
        type(starproject)                      :: starproj
        character(len=LONGSTRLEN), allocatable :: movies(:)
        character(len=LONGSTRLEN), allocatable :: completed_fnames(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_picker
        character(len=:),          allocatable :: output_dir_motion_correct, output_dir_extract
        character(len=LONGSTRLEN)              :: movie
        real                                   :: pickref_scale
        integer                                :: nchunks_imported_glob, nchunks_imported, box_extract
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection, iproj
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, ncls_in, nmic_star
        logical                                :: l_pick, l_movies_left, l_haschanged, l_cluster2d, l_nchunks_maxed, l_whether2D
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        if( .not. cline%defined('oritype')          ) call cline%set('oritype',        'mic')
        if( .not. cline%defined('mkdir')            ) call cline%set('mkdir',          'yes')
        if( .not. cline%defined('walltime')         ) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! motion correction
        if( .not. cline%defined('trs')              ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')             ) call cline%set('bfac',             50.)
        call cline%set('groupframes',     'no')
        if( .not. cline%defined('mcconvention')     ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')   ) call cline%set('eer_upsampling',    1.)
        if( .not. cline%defined('algorithm')        ) call cline%set('algorithm',    'patch')
        if( .not. cline%defined('mcpatch')          ) call cline%set('mcpatch',        'yes')
        if( .not. cline%defined('mcpatch_thres')    ) call cline%set('mcpatch_thres','  yes')
        if( .not. cline%defined('tilt_thres')       ) call cline%set('tilt_thres',      0.05)
        ! ctf estimation
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512.)
        if( .not. cline%defined('hp_ctf_estimate')  ) call cline%set('hp_ctf_estimate',  HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate')  ) call cline%set('lp_ctf_estimate',  LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',            DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',            DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',         'yes')
        if( .not. cline%defined('ctfresthreshold')  ) call cline%set('ctfresthreshold',  CTFRES_THRESHOLD)
        if( .not. cline%defined('icefracthreshold') ) call cline%set('icefracthreshold', ICEFRAC_THRESHOLD)
        ! picking
        if( .not. cline%defined('picker')          ) call cline%set('picker',         'old')
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          20.)
        if( .not. cline%defined('ndev')            ) call cline%set('ndev',              2.)
        if( .not. cline%defined('thres')           ) call cline%set('thres',            24.)
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('extractfrommov')  ) call cline%set('extractfrommov',  'no')
        ! 2D classification
        if( .not. cline%defined('lpthres')     ) call cline%set('lpthres',       30.0)
        if( .not. cline%defined('ndev2D')      ) call cline%set('ndev2D',         1.5)
        if( .not. cline%defined('wiener')      ) call cline%set('wiener',   'partial')
        if( .not. cline%defined('autoscale')   ) call cline%set('autoscale',    'yes')
        if( .not. cline%defined('nonuniform')  ) call cline%set('nonuniform',    'no')
        if( .not. cline%defined('nparts_chunk')) call cline%set('nparts_chunk',   1.0)
        if( .not. cline%defined('nchunks')     ) call cline%set('nchunks',        2.0)
        if( .not. cline%defined('prune')       ) call cline%set('prune',         'no')
        if( .not. cline%defined('reject_cls')  ) call cline%set('reject_cls',   'yes')
        if( .not. cline%defined('objfun')      ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',        'no')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ncls_in = 0
        if( cline%defined('ncls') )then
            ! to circumvent parameters class stringency, restored after params%new
            ncls_in = nint(cline%get_rarg('ncls'))
            call cline%delete('ncls')
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        if( ncls_in > 0 ) call cline%set('ncls', real(ncls_in))
        if( cline%defined('dir_prev') .and. .not.file_exists(params%dir_prev) )then
            THROW_HARD('Directory '//trim(params%dir_prev)//' does not exist!')
        endif
        if( .not.cline%defined('nparts_pool') )then
            params_glob%nparts_pool = params_glob%nparts
            call cline%set('nparts_pool', real(params_glob%nparts_pool))
        endif
        if( .not. cline%defined('maxnchunks') .or. params_glob%maxnchunks < 1 )then
            params_glob%maxnchunks = huge(params_glob%maxnchunks)
        endif
        call cline%delete('maxnchunks')
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj%read( params%projfile )
        call spproj%update_projinfo(cline)
        if( spproj%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        call spproj%write(micspproj_fname)
        ! picking
        l_pick = .false.
        if( cline%defined('picker') )then
            select case(trim(params%picker))
            case('old')
                if(.not.cline%defined('pickrefs')) THROW_HARD('PICKREFS required for picker=old')
            case('new')
                if(cline%defined('pickrefs'))then
                    if( .not. cline%defined('mskdiam') )then
                        THROW_HARD('New picker requires mask diameter (in A) in conjunction with pickrefs')
                    endif
                else
                    if( .not.cline%defined('moldiam') )then
                        THROW_HARD('MOLDIAM required for picker=new')
                    endif
                endif
            end select
            l_pick = .true.
        endif
        ! is 2D classification going to be performed & required parameters present
        l_whether2D = .false.
        if( l_pick ) call check_params_for_cluster2D(cline, l_whether2D)
        ! output directories
        output_dir = trim(PATH_HERE)//trim(dir_preprocess)
        call simple_mkdir(output_dir)
        call simple_mkdir(trim(output_dir)//trim(STDERROUT_DIR))
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
        if( l_pick .and. cline%defined('pickrefs') )then
            if( trim(params%picker).eq.'old' )then
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
                call cline%set('pickrefs', '../'//trim(PICKREFS)//trim(params%ext))
                write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
                call qsys_cleanup
            endif
        endif
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME)
        ! import previous runs
        nptcls_glob = 0
        call import_prev_streams
        ! start watching
        nchunks_imported_glob = 0
        last_injection        = simple_gettime()
        prev_stacksz          = 0
        nmovies               = 0
        iter                  = 0
        n_imported            = 0
        n_failed_jobs         = 0
        n_added               = 0
        nmic_star             = 0
        l_movies_left         = .false.
        l_haschanged          = .false.
        l_nchunks_maxed       = .false.
        l_cluster2D           = .false.
        do
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            iter = iter + 1
            ! detection of new movies
            if( l_nchunks_maxed )then
                call movie_buff%kill
                nmovies = 0
            else
                call movie_buff%watch( nmovies, movies )
            endif
            ! append movies to processing stack
            if( nmovies > 0 )then
                cnt = 0
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    if( movie_buff%is_past(movie) )cycle
                    call create_individual_project(movie)
                    call qenv%qscripts%add_to_streaming( cline )
                    call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
                    call movie_buff%add2history( movies(imovie) )
                    cnt = cnt+1
                    n_added = n_added+1
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
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz
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
            ! failed jobs
            if( qenv%qscripts%get_failed_stacksz() > 0 )then
                call qenv%qscripts%get_stream_fail_stack( failed_jobs_clines, n_fail_iter )
                if( n_fail_iter > 0 )then
                    n_failed_jobs = n_failed_jobs + n_fail_iter
                    do cnt = 1,n_fail_iter
                        call failed_jobs_clines(cnt)%kill
                    enddo
                    deallocate(failed_jobs_clines)
                endif
            endif
            ! project update
            if( n_imported > 0 )then
                n_imported = spproj%os_mic%get_noris()
                write(logfhandle,'(A,I5)')                         '>>> # MOVIES PROCESSED & IMPORTED:     ',n_imported
                if( l_pick ) write(logfhandle,'(A,I8)')            '>>> # PARTICLES EXTRACTED:          ',nptcls_glob
                write(logfhandle,'(A,I3,A1,I3)')                   '>>> # OF COMPUTING UNITS IN USE/TOTAL:   ',qenv%get_navail_computing_units(),'/',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I5)') '>>> # FAILED JOBS                :     ',n_failed_jobs
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                call spproj%write_segment_inside('mic',micspproj_fname)
                last_injection = simple_gettime()
                l_haschanged   = .true.
                n_imported     = spproj%os_mic%get_noris()
                ! always write micrographs snapshot if less than 1000 mics, else every 100
                if( n_imported < 1000 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_migrographs_starfile
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_migrographs_starfile
                    nmic_star = n_imported
                endif
                ! init cluster2D
                if( l_whether2D .and.(.not.l_cluster2D) )then
                    do iproj = 1,size(completed_fnames)
                        call tmp_spproj%kill
                        call tmp_spproj%read_segment('stk',completed_fnames(iproj))
                        if( .not.tmp_spproj%os_stk%isthere(iproj,'box') ) cycle
                        box_extract = nint(tmp_spproj%os_stk%get(iproj,'box')) ! getting particle size from first project
                        call init_cluster2D_stream(cline, spproj, box_extract, micspproj_fname, l_cluster2D)
                        call cline%delete('job_memory_per_task2D')
                        call cline%delete('qsys_partition2D')
                        call cline%delete('ncls')
                        exit
                    enddo
                    call tmp_spproj%kill
                endif
            else
                ! wait & write snapshot
                if( l_cluster2D )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        call update_user_params(cline)
                        call write_migrographs_starfile
                        l_haschanged = .false.
                    endif
                else
                    if( .not.l_movies_left )then
                        if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                            ! write project when inactive...
                            call write_project
                            call update_user_params(cline)
                            call write_migrographs_starfile
                            l_haschanged = .false.
                        else
                            ! ...or wait
                            call sleep(WAITTIME)
                        endif
                    endif
                endif
            endif
            ! read beamtilts if not 2D
            if(.not. l_cluster2D .and. cline%defined('dir_meta')) call read_xml_beamtilts()
            ! 2D classification section
            if( l_cluster2D )then
                call update_user_params(cline)
                call update_chunks
                call update_pool_status
                call update_pool
                call update_user_params(cline)
                call reject_from_pool
                call read_pool_xml_beamtilts()
                call assign_pool_optics(cline, propagate = .false.)
                ! call reject_from_pool_user
                if( .not.l_nchunks_maxed )then
                    call write_project_stream2D(.true.)
                    call import_chunks_into_pool(.false., nchunks_imported)
                    nchunks_imported_glob = nchunks_imported_glob + nchunks_imported
                    l_nchunks_maxed       = nchunks_imported_glob >= params_glob%maxnchunks
                    call classify_pool
                    call update_projects_mask(completed_fnames)
                    call classify_new_chunks(completed_fnames)
                else
                    ! # of chunks is above desired threshold
                    if( is_pool_available() ) exit
                endif
                
                call sleep(WAITTIME)
            endif
        end do
        ! termination
        if( l_cluster2D )then
            call read_pool_xml_beamtilts()
            call assign_pool_optics(cline, propagate = .true.)
            call terminate_stream2D(.true.)
        else
            call write_project
        endif
        call update_user_params(cline)
        call write_migrographs_starfile
        ! cleanup
        call spproj%kill
        call qsys_cleanup
        call del_file(micspproj_fname)
        ! end gracefully
        call simple_end('**** SIMPLE_PREPROCESS_STREAM NORMAL STOP ****')
        contains

            !>  write starfile snapshot
            subroutine write_migrographs_starfile
                integer(timer_int_kind)      :: ms0
                real(timer_int_kind)         :: ms_assign, ms_export
                
                if (spproj%os_mic%get_noris() > 0) then
                    if( .not. l_cluster2D ) then
                        if( DEBUG_HERE ) ms0 = tic()
                        call starproj%assign_optics(cline, spproj)
                        if( DEBUG_HERE )then
                            ms_assign = toc(ms0)
                            print *,'ms_assign  : ', ms_assign; call flush(6)
                        endif
                    end if
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%export_mics(cline, spproj)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                    if(allocated(starproj%tiltinfo)) deallocate(starproj%tiltinfo)
                end if
                
            end subroutine write_migrographs_starfile

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
                        call stream_spproj%os_stk%set_state(1, states(iproj))
                        n      = nint(stream_spproj%os_stk%get(1,'nptcls'))
                        fromp  = nptcls + 1
                        nptcls = nptcls + n
                        top    = nptcls
                        call spproj%os_stk%transfer_ori(istk,stream_spproj%os_stk,1)
                        call spproj%os_stk%set(istk, 'fromp',real(fromp))
                        call spproj%os_stk%set(istk, 'top',  real(top))
                    enddo
                    call spproj%write_segment_inside('stk', params%projfile)
                    call spproj%os_stk%kill
                    ! particles
                    call spproj%os_ptcl2D%new(nptcls, is_ptcl=.true.)
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
                            call spproj%os_ptcl2D%set_stkind(iptcl, istk)
                            call spproj%os_ptcl2D%set_state(iptcl, states(iproj))
                        enddo
                        call stream_spproj%kill
                    enddo
                    write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:          ',spproj%os_ptcl2D%get_noris()
                    call spproj%write_segment_inside('ptcl2D', params%projfile)
                    spproj%os_ptcl3D = spproj%os_ptcl2D
                    call spproj%os_ptcl2D%kill
                    call spproj%os_ptcl3D%delete_2Dclustering
                    call spproj%write_segment_inside('ptcl3D', params%projfile)
                    call spproj%os_ptcl3D%kill
                endif
                call spproj%write_non_data_segments(params%projfile)
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
                do i = 1,n_spprojs
                    ! flags zero-picked mics that will not be imported
                    fname = trim(output_dir)//trim(completed_jobs_clines(i)%get_carg('projfile'))
                    call check_nptcls(fname, nptcls_here, state, l_pick)
                    if( l_pick )then
                        spproj_mask(i) = (nptcls_here > 0) .and. (state > 0)
                    else
                        spproj_mask(i) = state > 0
                    endif
                enddo
                n2import      = count(spproj_mask)
                n_failed_jobs = n_failed_jobs + (n_spprojs-n2import)
                if( l_pick .and. n2import /= n_spprojs )then
                    write(logfhandle,'(A,I3,A)')'>>> NO PARTICLES FOUND IN ',n_spprojs-n2import,' MICROGRAPHS'
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
                            if( nptcls > 0 ) nptcls_here = nptcls_here + nptcls
                        endif
                        call spproj%os_mic%transfer_ori(n_old+j, streamspproj%os_mic, 1)
                        call streamspproj%kill
                        ! update paths such that relative paths are with respect to root folder
                        call update_path(spproj%os_mic, n_old+j, 'mc_starfile')
                        call update_path(spproj%os_mic, n_old+j, 'intg')
                        call update_path(spproj%os_mic, n_old+j, 'forctf')
                        call update_path(spproj%os_mic, n_old+j, 'thumb')
                        call update_path(spproj%os_mic, n_old+j, 'mceps')
                        call update_path(spproj%os_mic, n_old+j, 'ctfdoc')
                        call update_path(spproj%os_mic, n_old+j, 'ctfjpg')
                        if( l_pick ) call update_path(spproj%os_mic, n_old+j, 'boxfile')
                    enddo
                    if( l_pick ) nptcls_glob = nptcls_glob + nptcls_here
                else
                    nimported  = 0
                    return
                endif
                ! call write_filetable(STREAM_SPPROJFILES, completedfnames)
            end subroutine update_projects_list

            subroutine create_individual_project( movie )
                character(len=*), intent(in) :: movie
                type(sp_project)             :: spproj_here
                type(cmdline)                :: cline_here
                type(ctfparams)              :: ctfvars
                character(len=STDLEN)        :: ext, movie_here
                character(len=LONGSTRLEN)    :: projname, projfile,xmlfile,xmldir
                character(len=XLONGSTRLEN)   :: cwd, cwd_old
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
                if(cline%defined('dir_meta')) then
                    xmlfile = basename(trim(movie))
                    xmldir = cline%get_carg('dir_meta')
                    if(index(xmlfile, '_fractions') > 0) then
                        xmlfile = xmlfile(:index(xmlfile, '_fractions') - 1)
                    end if
                    if(index(xmlfile, '_EER') > 0) then
                        xmlfile = xmlfile(:index(xmlfile, '_EER') - 1)
                    end if
                    xmlfile = trim(adjustl(xmldir)) // '/' // trim(adjustl(xmlfile)) // '.xml'
                    call spproj_here%os_mic%set(1, "meta", trim(adjustl(xmlfile)))
                    call spproj_here%os_mic%set(1, "tiltx", 0.0)
                    call spproj_here%os_mic%set(1, "tilty", 0.0)
                end if
                call spproj_here%write
                call cline%set('projname', trim(projname))
                call cline%set('projfile', trim(projfile))
                call chdir(cwd_old)
                cwd_glob = trim(cwd_old)
                call spproj_here%kill
                call cline_here%kill
            end subroutine create_individual_project

            !>  import previous run to the current project based on past single project files
            subroutine import_prev_streams
                type(sp_project) :: streamspproj
                type(ori)        :: o, o_stk
                character(len=LONGSTRLEN), allocatable :: sp_files(:)
                character(len=:), allocatable :: mic, mov, dir
                logical,          allocatable :: spproj_mask(:)
                integer :: iproj,nprojs,icnt,nptcls
                logical :: err
                if( .not.cline%defined('dir_prev') ) return
                err = .false.
                dir = filepath(params%dir_prev, 'spprojs/')
                call simple_list_files_regexp(dir,'^'//trim(PREPROCESS_PREFIX)//'.*\.simple$',sp_files)
                if( .not.allocated(sp_files) )then
                    write(logfhandle,'(A)') '>>> Could not find previously processed movies'
                    return
                endif
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
                    call movefile2folder('intg',        dir, output_dir_motion_correct, o, err)
                    call movefile2folder('forctf',      dir, output_dir_motion_correct, o, err)
                    call movefile2folder('thumb',       dir, output_dir_motion_correct, o, err)
                    call movefile2folder('mc_starfile', dir, output_dir_motion_correct, o, err)
                    call movefile2folder('mceps',       dir, output_dir_motion_correct, o, err)
                    call movefile2folder('ctfjpg',      dir, output_dir_ctf_estimate,   o, err)
                    call movefile2folder('ctfdoc',      dir, output_dir_ctf_estimate,   o, err)
                    if( l_pick )then
                        ! import mic & updates stk segment
                        call movefile2folder('boxfile', dir, output_dir_picker, o, err)
                        nptcls = nptcls + nint(o%get('nptcls'))
                        call streamspproj%os_mic%set_ori(1, o)
                        if( .not.err )then
                            call streamspproj%read_segment('stk', sp_files(iproj))
                            if( streamspproj%os_stk%get_noris() == 1 )then
                                call streamspproj%os_stk%get_ori(1, o_stk)
                                call movefile2folder('stk', dir, output_dir_extract, o_stk, err)
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
                    call streamspproj%write(trim(dir_preprocess)//basename(sp_files(iproj)))
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

            subroutine movefile2folder(key, input_dir, folder, o, err)
                character(len=*), intent(in)    :: key, input_dir, folder
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
                src = filepath(input_dir, src)
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
                iostat = rename(src,dest)
                ! files will be used/imported from the spprojs folder
                call make_relativepath(trim(CWD_GLOB)//'/spprojs/',dest,reldest)
                call o%set(key,reldest)
            end subroutine movefile2folder
            
            subroutine read_xml_beamtilts()
                type(Node), pointer :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
                integer :: i
                do i = 1, spproj%os_mic%get_noris()
                    if ( spproj%os_mic%get(i, "tiltx") == 0.0 .and. spproj%os_mic%get(i, "tilty") == 0.0) then
                        if(file_exists(spproj%os_mic%get_static(i, "meta"))) then
                            write(logfhandle, *) "stream reading " // trim(adjustl(spproj%os_mic%get_static(i,"meta")))
                            xmldoc => parseFile(trim(adjustl(spproj%os_mic%get_static(i,"meta"))))
                            beamtiltnode => item(getElementsByTagname(xmldoc, "BeamShift"),0)
                            beamtiltnodex => item(getElementsByTagname(beamtiltnode, "a:_x"), 0)
                            beamtiltnodey => item(getElementsByTagname(beamtiltnode, "a:_y"), 0)
                            call spproj%os_mic%set(i, "tiltx", str2real(getTextContent(beamtiltnodex)))
                            call spproj%os_mic%set(i, "tilty", str2real(getTextContent(beamtiltnodey)))
                            call destroy(xmldoc)
                        endif
                    endif
                end do
            end subroutine read_xml_beamtilts

    end subroutine exec_preprocess_stream_dev

    ! Common utility functions

    subroutine check_nptcls( fname, nptcls, state, ptcl_check )
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: nptcls, state
        logical,          intent(in)  :: ptcl_check
        type(sp_project) :: spproj_here
        integer :: nmics, nstks
        state  = 0
        call spproj_here%read_data_info(fname, nmics, nstks, nptcls)
        if( ptcl_check )then
            if( nmics /= 1 )then
                THROW_WARN('No micrograph for: '//trim(fname)//'. Skipping')
                nptcls = 0
            else if( nptcls < 1 )then
                THROW_WARN('No particles extracted for: '//trim(fname)//'. Skipping')
                nptcls = 0
            else
                call spproj_here%read_segment('mic',fname)
                nptcls = nint(spproj_here%os_mic%get(1,'nptcls'))
                state  = spproj_here%os_mic%get_state(1)
                call spproj_here%kill
            endif
        else
            if( nmics /= 1 )then
                THROW_WARN('No micrograph for: '//trim(fname)//'. Skipping')
            else
                call spproj_here%read_segment('mic',fname)
                state  = spproj_here%os_mic%get_state(1)
                call spproj_here%kill
            endif
        endif
    end subroutine check_nptcls

end module simple_commander_preprocess_stream
