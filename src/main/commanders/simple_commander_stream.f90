! concrete commander: streaming pre-processing routines
module simple_commander_stream
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,        only: cmdline
use simple_parameters,     only: parameters, params_glob
use simple_commander_base, only: commander_base
use simple_sp_project,     only: sp_project
use simple_qsys_env,       only: qsys_env
use simple_starproject,    only: starproject
use simple_guistats,       only: guistats
use simple_qsys_funs
use simple_commander_preprocess
use simple_progress
implicit none

public :: commander_stream_preprocess
public :: commander_multipick_cluster2D
public :: commander_pick_extract_cluster2D

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_stream_preprocess
  contains
    procedure :: execute => exec_stream_preprocess
end type commander_stream_preprocess

type, extends(commander_base) :: commander_multipick_cluster2D
  contains
    procedure :: execute => exec_multipick_cluster2D
end type commander_multipick_cluster2D

type, extends(commander_base) :: commander_pick_extract_cluster2D
  contains
    procedure :: execute => exec_pick_extract_cluster2D
end type commander_pick_extract_cluster2D

! module constants
character(len=STDLEN), parameter :: DIR_STREAM           = trim(PATH_HERE)//'spprojs/'           ! location for projects to be processed
character(len=STDLEN), parameter :: DIR_STREAM_COMPLETED = trim(PATH_HERE)//'spprojs_completed/' ! location for projects processed
character(len=STDLEN), parameter :: USER_PARAMS     = 'stream_user_params.txt'                   ! really necessary here? - left in for now
integer,               parameter :: NMOVS_SET       = 5                                          ! number of movies processed at once (>1)
integer,               parameter :: LONGTIME        = 120                                        ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME        = 3    ! movie folder watched every WAITTIME seconds

! module variables
integer :: movies_set_counter = 0

contains

    subroutine exec_stream_preprocess( self, cline )
        use simple_moviewatcher, only: moviewatcher
        class(commander_stream_preprocess), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .false.
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(cmdline)                          :: cline_exec
        type(qsys_env)                         :: qenv
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj_glob    ! global project
        type(starproject)                      :: starproj
        character(len=LONGSTRLEN), allocatable :: movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_motion_correct
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection
        integer                                :: cnt, n_imported, n_added, n_failed_jobs, n_fail_iter, nmic_star, iset
        logical                                :: l_movies_left, l_haschanged
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('reject_mics', 'no')
        if( .not. cline%defined('walltime')         ) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! motion correction
        call cline%set('groupframes', 'no')
        if( .not. cline%defined('trs')              ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')             ) call cline%set('bfac',             50.)
        if( .not. cline%defined('mcconvention')     ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')   ) call cline%set('eer_upsampling',    1.)
        if( .not. cline%defined('algorithm')        ) call cline%set('algorithm',    'patch')
        if( .not. cline%defined('mcpatch')          ) call cline%set('mcpatch',        'yes')
        if( .not. cline%defined('mcpatch_thres')    ) call cline%set('mcpatch_thres','  yes')
        if( .not. cline%defined('tilt_thres')       ) call cline%set('tilt_thres',      0.05)
        if( .not. cline%defined('beamtilt')         ) call cline%set('beamtilt',        'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',         512.)
        if( .not. cline%defined('hp_ctf_estimate')  ) call cline%set('hp_ctf_estimate', HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate')  ) call cline%set('lp_ctf_estimate', LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',           DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',           DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',        'yes')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            ! restart
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME, params%dir_movies)
        ! restart
        movies_set_counter = 0  ! global number of movies set
        nmic_star          = 0
        if( cline%defined('dir_exec') )then
            call import_previous_projects
            nmic_star = spproj_glob%os_mic%get_noris()
            call write_mic_star_and_field(write_field=.true.)
        endif
        ! output directories
        call simple_mkdir(trim(PATH_HERE)//trim(DIR_STREAM_COMPLETED))
        output_dir = trim(PATH_HERE)//trim(DIR_STREAM)
        call simple_mkdir(output_dir)
        call simple_mkdir(trim(output_dir)//trim(STDERROUT_DIR))
        output_dir_ctf_estimate   = filepath(trim(PATH_HERE), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(PATH_HERE), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream :: exec_preprocess_stream;  ")
        call cline%set('dir','../')
        ! initialise progress monitor
        call progressfile_init()
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! Infinite loop
        last_injection = simple_gettime()
        prev_stacksz   = 0
        iter           = 0
        n_imported     = 0
        n_failed_jobs  = 0
        n_added        = 0
        l_movies_left  = .false.
        l_haschanged   = .false.
        cline_exec     = cline
        call cline_exec%set('fromp',1)
        call cline_exec%set('top',  NMOVS_SET)
        ! guistats init
        call gui_stats%init
        do
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            iter = iter + 1
            ! detection of new movies
            call movie_buff%watch( nmovies, movies, max_nmovies=params%nparts*NMOVS_SET )
            ! append movies to processing stack
            if( nmovies >= NMOVS_SET )then
                cnt = 0
                do iset = 1,nmovies,NMOVS_SET
                    call create_movies_set_project(movies(iset:iset+NMOVS_SET-1))
                    call qenv%qscripts%add_to_streaming( cline_exec )
                    do imovie = 1,NMOVS_SET
                        call movie_buff%add2history( movies(iset+imovie-1) )
                        cnt     = cnt     + 1
                        n_added = n_added + 1 ! global number of movie sets
                    enddo
                    if( cnt == min(params%nparts*NMOVS_SET,nmovies) ) exit
                enddo
                write(logfhandle,'(A,I4,A,A)')'>>> ',cnt,' NEW MOVIES ADDED; ', cast_time_char(simple_gettime())
                l_movies_left = cnt .ne. nmovies
            else
                l_movies_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz*NMOVS_SET
                ! guistats
                call gui_stats%set('micrographs', 'movies', int2str(spproj_glob%os_mic%get_noris()) // '/' // int2str(stacksz*NMOVS_SET + spproj_glob%os_mic%get_noris()), primary=.true.)
            endif
            ! fetch completed jobs list
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                call update_projects_list( n_imported )
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
                n_imported = spproj_glob%os_mic%get_noris()
                write(logfhandle,'(A,I8)')                         '>>> # MOVIES PROCESSED & IMPORTED       : ',n_imported
                write(logfhandle,'(A,I3,A2,I3)')                   '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                call gui_stats%set('micrographs', 'movies',  int2str(n_imported) // '/' // int2str(stacksz + spproj_glob%os_mic%get_noris()), primary=.true.)
                call gui_stats%set('micrographs', 'compute', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
                if( n_failed_jobs > 0 ) call gui_stats%set('micrographs', 'rejected', n_failed_jobs, primary=.true.)
                if(spproj_glob%os_mic%isthere("ctfres")) then
                    call gui_stats%set('micrographs', 'avg_ctf_res', spproj_glob%os_mic%get_avg("ctfres"), primary=.true.)
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                last_injection = simple_gettime()
                ! guistats
                call gui_stats%set_now('micrographs', 'last_new_movie')
                if(spproj_glob%os_mic%isthere('thumb')) then
                    call gui_stats%set('micrographs', 'latest_micrograph', trim(adjustl(CWD_GLOB))//'/'//&
                        &trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(),'thumb'))), thumbnail=.true.)
                end if
                l_haschanged = .true.
                n_imported   = spproj_glob%os_mic%get_noris()
                ! always write micrographs snapshot if less than 1000 mics, else every 100
                if( n_imported < 1000 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_mic_star_and_field
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_mic_star_and_field
                    nmic_star = n_imported
                endif
            else
                ! wait & write snapshot
                if( .not.l_movies_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive...
                        call update_user_params(cline)
                        call write_mic_star_and_field
                        l_haschanged = .false.
                    else
                        ! ...or wait
                        call sleep(WAITTIME)
                    endif
                endif
            endif
            ! read beamtilts
            if( cline%defined('dir_meta')) call read_xml_beamtilts(spproj_glob)
            ! guistats
            call gui_stats%write_json
        end do
        ! termination
        call update_user_params(cline)
        call write_mic_star_and_field(write_field=.true.)
        ! final stats
        call gui_stats%hide('micrographs', 'compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_PREPROC NORMAL STOP ****')
        contains

            subroutine write_mic_star_and_field( write_field )
                logical, optional, intent(in) :: write_field
                logical :: l_wfield
                l_wfield = .false.
                if( present(write_field) ) l_wfield = write_field
                call write_migrographs_starfile
                if( l_wfield )then
                    call spproj_glob%write_segment_inside('mic', params%projfile)
                    call spproj_glob%write_non_data_segments(params%projfile)
                endif
            end subroutine write_mic_star_and_field

            !>  write starfile snapshot
            subroutine write_migrographs_starfile
                integer(timer_int_kind)      :: ms0
                real(timer_int_kind)         :: ms_assign, ms_export
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%assign_optics(cline, spproj_glob)
                    if( DEBUG_HERE )then
                        ms_assign = toc(ms0)
                        print *,'ms_assign  : ', ms_assign; call flush(6)
                    endif
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%export_mics(cline, spproj_glob)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                    if(allocated(starproj%tiltinfo)) deallocate(starproj%tiltinfo)
                end if
            end subroutine write_migrographs_starfile

            ! returns list of completed jobs
            subroutine update_projects_list( nimported )
                integer,                   intent(out) :: nimported
                type(sp_project), allocatable          :: streamspprojs(:)
                character(len=LONGSTRLEN), allocatable :: completed_fnames(:)
                character(len=:),          allocatable :: fname, abs_fname
                logical,                   allocatable :: mics_mask(:)
                integer :: i, n_spprojs, n_old, j, n2import, n_completed, iproj, nmics, imic, cnt
                n_completed = 0
                nimported   = 0
                n_spprojs = size(completed_jobs_clines) ! projects to import
                if( n_spprojs == 0 )return
                n_old = spproj_glob%os_mic%get_noris()       ! previously processed mmovies
                nmics = NMOVS_SET * n_spprojs           ! incoming number of processed movies
                allocate(streamspprojs(n_spprojs), completed_fnames(n_spprojs), mics_mask(nmics))
                ! read all
                do iproj = 1,n_spprojs
                    fname     = trim(output_dir)//trim(completed_jobs_clines(iproj)%get_carg('projfile'))
                    abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
                    completed_fnames(iproj) = trim(abs_fname)
                    call streamspprojs(iproj)%read_segment('mic', abs_fname)
                    cnt = 0
                    do imic = (iproj-1)*NMOVS_SET+1, iproj*NMOVS_SET
                        cnt = cnt + 1
                        mics_mask(imic) = streamspprojs(iproj)%os_mic%get_state(cnt) == 1
                    enddo
                enddo
                n2import      = count(mics_mask)
                n_failed_jobs = n_failed_jobs + (nmics-n2import)
                if( n2import > 0 )then
                    ! reallocate global project
                    n_completed = n_old + n2import
                    nimported   = n2import
                    if( n_old == 0 )then
                        ! first time
                        call spproj_glob%os_mic%new(n2import, is_ptcl=.false.)
                    else
                        call spproj_glob%os_mic%reallocate(n_completed)
                    endif
                    ! actual import
                    imic = 0
                    j    = n_old
                    do iproj = 1,n_spprojs
                        do i = 1,NMOVS_SET
                            imic = imic+1
                            if( mics_mask(imic) )then
                                j   = j + 1
                                ! transfer info
                                call spproj_glob%os_mic%transfer_ori(j, streamspprojs(iproj)%os_mic, i)
                                ! update paths such that relative paths are with respect to root folder
                                call update_path(spproj_glob%os_mic, j, 'mc_starfile')
                                call update_path(spproj_glob%os_mic, j, 'intg')
                                call update_path(spproj_glob%os_mic, j, 'forctf')
                                call update_path(spproj_glob%os_mic, j, 'thumb')
                                call update_path(spproj_glob%os_mic, j, 'mceps')
                                call update_path(spproj_glob%os_mic, j, 'ctfdoc')
                                call update_path(spproj_glob%os_mic, j, 'ctfjpg')
                            endif
                        enddo
                        call streamspprojs(iproj)%kill
                    enddo
                endif
                ! finally we move the completed projects to appropriate directory
                do iproj = 1,n_spprojs
                    imic = (iproj-1)*NMOVS_SET+1
                    if( any(mics_mask(imic:imic+NMOVS_SET-1)) )then
                        fname = trim(DIR_STREAM_COMPLETED)//trim(basename(completed_fnames(iproj)))
                        call simple_rename(completed_fnames(iproj), fname)
                    endif
                enddo
                ! cleanup
                do iproj = 1,n_spprojs
                    call completed_jobs_clines(iproj)%kill
                enddo
                deallocate(completed_jobs_clines)
                deallocate(streamspprojs,mics_mask,completed_fnames)
            end subroutine update_projects_list

            subroutine create_movies_set_project( movie_names )
                character(len=*), intent(in) :: movie_names(NMOVS_SET)
                type(sp_project)             :: spproj_here
                type(ctfparams)              :: ctfvars
                character(len=LONGSTRLEN)    :: projname, projfile,xmlfile,xmldir
                character(len=XLONGSTRLEN)   :: cwd, cwd_old
                integer :: imovie
                cwd_old = trim(cwd_glob)
                call chdir(output_dir)
                call simple_getcwd(cwd)
                cwd_glob = trim(cwd)
                ! movies set
                movies_set_counter = movies_set_counter + 1
                projname   = int2str_pad(movies_set_counter,params%numlen)
                projfile   = trim(projname)//trim(METADATA_EXT)
                call cline_exec%set('projname', trim(projname))
                call cline_exec%set('projfile', trim(projfile))
                call spproj_here%update_projinfo(cline_exec)
                spproj_here%compenv  = spproj_glob%compenv
                spproj_here%jobproc  = spproj_glob%jobproc
                ! movies parameters
                ctfvars%ctfflag      = CTFFLAG_YES
                ctfvars%smpd         = params%smpd
                ctfvars%cs           = params%cs
                ctfvars%kv           = params%kv
                ctfvars%fraca        = params%fraca
                ctfvars%l_phaseplate = params%phaseplate.eq.'yes'
                call spproj_here%add_movies(movie_names(1:NMOVS_SET), ctfvars, verbose = .false.)
                if(cline%defined('dir_meta')) then
                    xmldir = cline%get_carg('dir_meta')
                    do imovie = 1,NMOVS_SET
                        xmlfile = basename(trim(movie_names(imovie)))
                        if(index(xmlfile, '_fractions') > 0) xmlfile = xmlfile(:index(xmlfile, '_fractions') - 1)
                        if(index(xmlfile, '_EER') > 0)       xmlfile = xmlfile(:index(xmlfile, '_EER') - 1)
                        xmlfile = trim(adjustl(xmldir))//'/'//trim(adjustl(xmlfile))//'.xml'
                        call spproj_here%os_mic%set(imovie, "meta", trim(adjustl(xmlfile)))
                        call spproj_here%os_mic%set(imovie, "tiltx", 0.0)
                        call spproj_here%os_mic%set(imovie, "tilty", 0.0)
                    enddo
                end if
                call spproj_here%write
                call chdir(cwd_old)
                cwd_glob = trim(cwd_old)
                call spproj_here%kill
            end subroutine create_movies_set_project

            !>  import previous movies and updates global project & variables
            subroutine import_previous_projects
                type(sp_project),          allocatable :: spprojs(:)
                character(len=LONGSTRLEN), allocatable :: completed_fnames(:)
                character(len=:),          allocatable :: fname
                logical,                   allocatable :: mics_mask(:)
                character(len=LONGSTRLEN)              :: moviename
                integer :: n_spprojs, iproj, nmics, imic, jmic, cnt, iostat,id
                ! previously completed projects
                call simple_list_files_regexp(DIR_STREAM_COMPLETED, '\.simple$', completed_fnames)
                if( .not.allocated(completed_fnames) )then
                    return ! nothing was previously completed
                endif
                n_spprojs = size(completed_fnames)
                ! import into global project
                allocate(spprojs(n_spprojs), mics_mask(n_spprojs*NMOVS_SET))
                jmic = 0
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%read_segment('mic', completed_fnames(iproj))
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        jmic = jmic + 1
                        mics_mask(jmic) = spprojs(iproj)%os_mic%get_state(imic) == 1
                    enddo
                enddo
                nmics = count(mics_mask)
                if( nmics ==0 )then
                    ! nothing to import
                    do iproj = 1,n_spprojs
                        call spprojs(iproj)%kill
                        fname = trim(DIR_STREAM_COMPLETED)//trim(completed_fnames(iproj))
                        call del_file(fname)
                    enddo
                    deallocate(spprojs,mics_mask)
                    return
                endif
                call spproj_glob%os_mic%new(nmics, is_ptcl=.false.)
                jmic = 0
                cnt  = 0
                do iproj = 1,n_spprojs
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        cnt = cnt + 1
                        if( mics_mask(cnt) )then
                            jmic = jmic + 1
                            call spproj_glob%os_mic%transfer_ori(jmic, spprojs(iproj)%os_mic, imic)
                        endif
                    enddo
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
                ! update global movie set counter
                movies_set_counter = 0
                do iproj = 1,n_spprojs
                    fname = basename_safe(completed_fnames(iproj))
                    fname = trim(get_fbody(trim(fname),trim(METADATA_EXT),separator=.false.))
                    call str2int(fname, iostat, id)
                    if( iostat==0 ) movies_set_counter = max(movies_set_counter, id)
                enddo
                ! add previous movies to history
                do imic = 1,spproj_glob%os_mic%get_noris()
                    moviename = spproj_glob%os_mic%get_static(imic,'movie')
                    call movie_buff%add2history(moviename)
                enddo
                ! tidy files
                call simple_rmdir(DIR_STREAM)
                write(logfhandle,'(A,I6,A)')'>>> IMPORTED ',nmics,' PREVIOUSLY PROCESSED MOVIES'
            end subroutine import_previous_projects

    end subroutine exec_stream_preprocess

    subroutine exec_multipick_cluster2D( self, cline )
        use simple_commander_preprocess, only: pick_commander_distr
        use simple_commander_preprocess, only: extract_commander_distr
        use simple_commander_preprocess, only: make_pickrefs_commander
        use simple_commander_cluster2D,  only: cleanup2D_commander_hlev
        class(commander_multipick_cluster2D), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        integer, parameter :: NPTCLS_PER_CLS = 600
        type(pick_commander_distr)     :: xpick
        type(extract_commander_distr)  :: xextract
        type(make_pickrefs_commander)  :: xmake_pickrefs
        type(cleanup2D_commander_hlev) :: xcleanup2D
        type(parameters) :: params
        type(sp_project) :: spproj_glob
        type(cmdline)    :: cline_multipick, cline_pick, cline_extract
        type(cmdline)    :: cline_pickrefs, cline_cleanup2D
        type(ran_tabu)   :: rt
        character(len=:), allocatable :: cavgs
        integer,          allocatable :: states(:), states_backup(:), vec(:)
        real    :: smpd
        integer :: nmics, ncls, cnt, nmics_sel, imic
        if( .not. cline%defined('mkdir')       ) call cline%set('mkdir',       'yes')
        if( .not. cline%defined('pcontrast')   ) call cline%set('pcontrast', 'black')
        if( .not. cline%defined('thres')       ) call cline%set('thres',         24.)
        if( .not. cline%defined('pick_roi')    ) call cline%set('pick_roi',     'no')
        if( .not. cline%defined('backgr_subtr')) call cline%set('backgr_subtr', 'no')
        if( .not. cline%defined('nran')        ) call cline%set('nran',          50.)
        if( .not. cline%defined('walltime')    ) call cline%set('walltime',29.0*60.0) ! 29 minutes
        call cline%set('picker',  'new')
        call cline%set('oritype', 'mic')
        if( .not.cline%defined('moldiam') .and. .not.cline%defined('multi_moldiams') )then
            THROW_HARD('MOLDIAM or MULTI_MOLDIAMS must be defined!')
        endif
        call params%new(cline)
        ! sanity check
        call spproj_glob%read(params%projfile)
        nmics = spproj_glob%get_nintgs()
        if( nmics == 0 ) THROW_HARD('No micrograph to process! exec_pick_distr')
        call spproj_glob%update_projinfo(cline)
        call spproj_glob%write_segment_inside('projinfo')
        ! set mkdir to no (to avoid nested directory structure)
        call cline%set('mkdir', 'no')
        if( cline%defined('nparts') )then
            params%nparts = min(params%nparts, nmics)
            call cline%set('nparts', params%nparts)
        endif
        params%numlen = len(int2str(params%nparts))
        call cline%set('numlen', real(params%numlen))
        ! more sanity checks
        if( trim(params%pick_roi).eq.'yes' )then
            params%backgr_subtr = 'yes'
            call cline%set('backgr_subtr', params%backgr_subtr)
        endif
        ! Micrographs rejection and random selection
        nmics         = spproj_glob%os_mic%get_noris()
        states        = nint(spproj_glob%os_mic%get_all('state'))
        states_backup = states
        if( cline%defined('ctfresthreshold') )then
            do imic = 1, nmics
                if( states(imic) == 1 )then
                    if(spproj_glob%os_mic%get(imic,'ctfres') > params%ctfresthreshold) states(imic) = 0
                endif
            enddo
        endif
        if( cline%defined('icefracthreshold') )then
            do imic = 1, nmics
                if( states(imic) == 1 )then
                    if(spproj_glob%os_mic%get(imic,'icefrac') > params%icefracthreshold) states(imic) = 0
                endif
            enddo
        endif
        nmics_sel     = count(states == 1)
        if( nmics_sel < params%nran )then
            THROW_HARD('Insufficient number of micrographs!')
        endif
        rt = ran_tabu(nmics_sel)
        allocate(vec(nmics_sel),source=0)
        vec(1:params%nran) = 1
        call rt%shuffle(vec)
        cnt = 0
        do imic = 1,nmics
            if( states(imic)==1 )then
                cnt = cnt+1
                call spproj_glob%os_mic%set(imic,'state',real(vec(cnt)))
            endif
        enddo
        call rt%kill
        deallocate(vec)
        call spproj_glob%write_segment_inside('mic',params%projfile)
        ! multi pick
        write(logfhandle,'(A)')'>>> PERFORMING MULTI-DIAMETER PICKING'
        cline_multipick = cline
        call cline_multipick%set('prg', 'pick')
        call xpick%execute(cline_multipick)
        call qsys_cleanup
        ! diameter selection
        ! interaction with gui here
        ! ... TBD: moldiam and any other relevant option
        params%moldiam = 138.5 !!!!!!!!!
        ! single pick
        write(logfhandle,'(A)')'>>> PERFORMING SINGLE DIAMETER PICKING'
        call cline%delete('nmoldiams')
        call cline%delete('moldiam_max')
        call cline%delete('multi_moldiam')
        call cline%set('moldiam', params%moldiam)
        cline_pick = cline
        call cline_pick%set('prg', 'pick')
        call xpick%execute(cline_pick)
        call qsys_cleanup
        ! extraction
        write(logfhandle,'(A)')'>>> PERFORMING EXTRACTION'
        cline_extract = cline
        call cline_extract%set('prg', 'extract')
        call xextract%execute(cline_extract)
        ! 2D classification
        write(logfhandle,'(A)')'>>> PERFORMING 2D CLASSIFICATION'
        params%mskdiam = 1.1*params%moldiam
        call spproj_glob%read_segment('stk', params%projfile)
        params%nptcls = spproj_glob%get_nptcls()
        params%ncls   = min(30,ceiling(real(params%nptcls)/real(NPTCLS_PER_CLS)))
        call spproj_glob%os_ptcl2D%kill
        write(logfhandle,'(A,I6)')'>>> # OF PARTICLES: ', params%nptcls
        write(logfhandle,'(A,I6)')'>>> # OF CLASSES  : ', params%ncls
        cline_cleanup2D = cline
        call cline_cleanup2D%set('prg',     'cleanup2D')
        call cline_cleanup2D%set('oritype', 'ptcl2D')
        call cline_cleanup2D%set('mskdiam', params%mskdiam)
        call cline_cleanup2D%set('ncls',    params%ncls)
        call xcleanup2D%execute(cline_cleanup2D)
        call qsys_cleanup
        ! restores states
        call spproj_glob%os_mic%set_all('state',real(states_backup))
        call spproj_glob%write_segment_inside('mic',params%projfile)
        ! class averages selection here
        ! ...
        ! picking references
        write(logfhandle,'(A)')'>>> GENERATING PICKING REFERENCES'
        call spproj_glob%read_segment('out', params%projfile)
        call spproj_glob%get_cavgs_stk(cavgs, ncls, smpd, imgkind='cavg')
        cline_pickrefs = cline
        call cline_pickrefs%set('prg',  'make_pickrefs')
        call cline_pickrefs%set('pickrefs', cavgs)
        call xmake_pickrefs%execute_shmem(cline_pickrefs)
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        call simple_end('**** SIMPLE_MULTIPICK_CLUSTER2D NORMAL STOP ****')
    end subroutine exec_multipick_cluster2D

    subroutine exec_pick_extract_cluster2D( self, cline )
        use simple_moviewatcher, only: moviewatcher
        use simple_stream_chunk, only: micproj_record
        use simple_commander_cluster2D_stream_dev
        use simple_timer
        class(commander_pick_extract_cluster2D), intent(inout) :: self
        class(cmdline),                          intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .false.
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(micproj_record),      allocatable :: micproj_records(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs, cline_pick_extract
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob, stream_spproj, tmp_spproj
        type(starproject)                      :: starproj
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=:),          allocatable :: output_dir, output_dir_extract, output_dir_picker
        character(len=LONGSTRLEN)              :: cwd_job
        integer                                :: nmics_sel, nmics_rej, nmics_rejected_glob
        integer                                :: nchunks_imported_glob, nchunks_imported, box_extract
        integer                                :: nprojects, stacksz, prev_stacksz, iter, last_injection, iproj
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, ncls_in, nmic_star
        logical                                :: l_templates_provided, l_projects_left, l_haschanged, l_cluster2d, l_nchunks_maxed, l_whether2D
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        call cline%set('oritype',   'mic')
        call cline%set('mkdir',     'yes')
        call cline%set('autoscale', 'yes')
        if( .not.cline%defined('dir_target') )then
            THROW_HARD('DIR_TARGET must be defined!')
        endif
        if( .not. cline%defined('walltime')) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! micrograph selection
        if( .not. cline%defined('reject_mics')     ) call cline%set('reject_mics',      'yes')
        if( .not. cline%defined('ctfresthreshold') ) call cline%set('ctfresthreshold',  CTFRES_THRESHOLD_STREAM)
        if( .not. cline%defined('icefracthreshold')) call cline%set('icefracthreshold', ICEFRAC_THRESHOLD_STREAM)
        ! picking
        if( .not. cline%defined('picker')      ) call cline%set('picker',         'new')
        if( .not. cline%defined('lp_pick')     ) call cline%set('lp_pick',         PICK_LP_DEFAULT)
        if( .not. cline%defined('ndev')        ) call cline%set('ndev',              2.)
        if( .not. cline%defined('thres')       ) call cline%set('thres',            24.)
        if( .not. cline%defined('pick_roi')    ) call cline%set('pick_roi',        'no')
        if( .not. cline%defined('backgr_subtr')) call cline%set('backgr_subtr',    'no')
        ! extraction
        if( .not. cline%defined('pcontrast')     ) call cline%set('pcontrast',    'black')
        if( .not. cline%defined('extractfrommov')) call cline%set('extractfrommov',  'no')
        ! 2D classification
        call cline%set('wiener', 'full')
        if( .not. cline%defined('lpthres')     ) call cline%set('lpthres',       30.0)
        if( .not. cline%defined('ndev2D')      ) call cline%set('ndev2D', CLS_REJECT_STD)
        if( .not. cline%defined('nonuniform')  ) call cline%set('nonuniform',    'no')
        if( .not. cline%defined('nparts_chunk')) call cline%set('nparts_chunk',   1.0)
        if( .not. cline%defined('nchunks')     ) call cline%set('nchunks',        2.0)
        if( .not. cline%defined('prune')       ) call cline%set('prune',         'no')
        if( .not. cline%defined('reject_cls')  ) call cline%set('reject_cls',   'yes')
        if( .not. cline%defined('objfun')      ) call cline%set('objfun',    'euclid')
        if( .not. cline%defined('ml_reg')      ) call cline%set('ml_reg',        'no')
        if( .not. cline%defined('rnd_cls_init')) call cline%set('rnd_cls_init',  'no')
        if( .not. cline%defined('remove_chunks'))call cline%set('remove_chunks','yes')
        if( .not. cline%defined('kweight_chunk'))call cline%set('kweight_chunk','default')
        if( .not. cline%defined('kweight_pool') )call cline%set('kweight_pool', 'default')
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
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        cline_pick_extract = cline
        call cline_pick_extract%set('prg', 'pick_extract')
        call cline_pick_extract%set('dir','../')
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
        if( trim(params%pick_roi).eq.'yes' )then
            params%backgr_subtr = 'yes'
            call cline%set('backgr_subtr', params%backgr_subtr)
        endif
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('pick_extract_cluster2D must start from an empty project (eg from root project folder)')
        ! picking
        l_templates_provided = cline%defined('pickrefs')
        if( cline%defined('picker') )then
            select case(trim(params%picker))
            case('old')
                if( .not.l_templates_provided ) THROW_HARD('PICKREFS required for picker=old')
            case('new')
                if( l_templates_provided )then
                    if( .not. cline%defined('mskdiam') )then
                        THROW_HARD('New picker requires mask diameter (in A) in conjunction with pickrefs')
                    endif
                else if( .not.cline%defined('moldiam') )then
                    THROW_HARD('MOLDIAM required for picker=new reference-free picking')
                endif
            case DEFAULT
                THROW_HARD('Unsupported picker')
            end select
        endif
        ! output directories
        output_dir = trim(PATH_HERE)//trim(DIR_STREAM)
        call simple_mkdir(output_dir)
        call simple_mkdir(trim(output_dir)//trim(STDERROUT_DIR))
        call simple_mkdir(trim(PATH_HERE)//trim(DIR_STREAM_COMPLETED))
        output_dir_picker  = filepath(trim(PATH_HERE), trim(DIR_PICKER))
        output_dir_extract = filepath(trim(PATH_HERE), trim(DIR_EXTRACT))
        call simple_mkdir(output_dir_picker,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_extract,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! prepares picking references
        if( l_templates_provided )then
            if( trim(params%picker).eq.'old' )then
                cline_make_pickrefs = cline
                call cline_make_pickrefs%set('prg','make_pickrefs')
                call cline_make_pickrefs%set('stream','no')
                call cline_make_pickrefs%delete('ncls')
                call cline_make_pickrefs%delete('mskdiam')
                call qenv%exec_simple_prg_in_queue(cline_make_pickrefs, 'MAKE_PICKREFS_FINISHED')
                call cline%set('pickrefs', '../'//trim(PICKREFS_FBODY)//trim(params%ext))
                write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
                call qsys_cleanup
            endif
        endif
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! import previous runs
        nptcls_glob = 0             ! global number of particles
        ! call import_prev_streams
        ! 2D classification
        l_whether2D = .false.
        l_cluster2D = .false.
        call check_params_for_cluster2D_dev(cline, l_whether2D)
        ! Infinitie loop
        nchunks_imported_glob = 0
        last_injection        = simple_gettime()
        prev_stacksz          = 0
        nprojects             = 0
        iter                  = 0
        n_imported            = 0   ! global number of imported processed micrographs
        n_failed_jobs         = 0
        n_added               = 0   ! global number of micrographs added to preocessing stack
        nmic_star             = 0
        nmics_rejected_glob   = 0   ! global number of micrographs rejected
        l_projects_left       = .false.
        l_haschanged          = .false.
        l_nchunks_maxed       = .false.
        ! guistats init
        call gui_stats%init
        do
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            iter = iter + 1
            ! detection of new projects
            if( l_nchunks_maxed )then
                call project_buff%kill
                nprojects = 0
            else
                call project_buff%watch( nprojects, projects, max_nmovies=params%nparts )
            endif
            ! append projects to processing stack
            if( nprojects > 0 )then
                cnt = 0
                do iproj = 1, nprojects
                    call create_individual_project(projects(iproj), nmics_sel, nmics_rej)
                    call project_buff%add2history(projects(iproj))
                    if( nmics_sel > 0 ) call qenv%qscripts%add_to_streaming(cline_pick_extract)
                    cnt     = cnt     + nmics_sel
                    n_added = n_added + nmics_sel
                    nmics_rejected_glob = nmics_rejected_glob + nmics_rej
                    if( cnt == min(params%nparts*NMOVS_SET, nprojects) ) exit
                enddo
                write(logfhandle,'(A,I4,A,A)')'>>> ',cnt,' NEW MICROGRAPHS ADDED; ', cast_time_char(simple_gettime())
                l_projects_left = cnt .ne. nprojects
            else
                l_projects_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz
                ! guistats
                call gui_stats%set('micrographs', 'movies', int2str(spproj_glob%os_mic%get_noris()) // '/' // int2str(stacksz + spproj_glob%os_mic%get_noris()), primary=.true.)
            endif
            ! fetch completed jobs list & updates of cluster2D_stream
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                call update_projects_list( micproj_records, n_imported )
                call completed_jobs_clines(:)%kill
                deallocate(completed_jobs_clines)
            else
                n_imported = 0 ! newly imported
            endif
            ! failed jobs
            if( qenv%qscripts%get_failed_stacksz() > 0 )then
                call qenv%qscripts%get_stream_fail_stack( failed_jobs_clines, n_fail_iter )
                if( n_fail_iter > 0 )then
                    n_failed_jobs = n_failed_jobs + n_fail_iter
                    call failed_jobs_clines(:)%kill
                    deallocate(failed_jobs_clines)
                endif
            endif
            ! project update
            if( n_imported > 0 )then
                n_imported = spproj_glob%os_mic%get_noris()
                write(logfhandle,'(A,I8)')       '>>> # MOVIES PROCESSED & IMPORTED       : ',n_imported
                write(logfhandle,'(A,I8)')       '>>> # PARTICLES EXTRACTED               : ',nptcls_glob
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                call gui_stats%set('micrographs', 'movies',  int2str(n_imported) // '/' // int2str(stacksz + spproj_glob%os_mic%get_noris()), primary=.true.)
                call gui_stats%set('micrographs', 'compute', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
                call gui_stats%set('micrographs', 'ptcls', nptcls_glob, primary=.true.)
                if( n_failed_jobs > 0 ) call gui_stats%set('micrographs', 'rejected', n_failed_jobs, primary=.true.)
                if(spproj_glob%os_mic%isthere("ctfres")) then
                    call gui_stats%set('micrographs', 'avg_ctf_res', spproj_glob%os_mic%get_avg("ctfres"), primary=.true.)
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                ! call spproj_glob%write_segment_inside('mic',micspproj_fname) !??
                last_injection = simple_gettime()
                ! guistats
                call gui_stats%set_now('micrographs', 'last_new_movie')
                if(spproj_glob%os_mic%isthere('thumb')) then
                    call gui_stats%set('micrographs', 'latest_micrograph', trim(adjustl(cwd_job)) // '/' // trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'thumb'))), thumbnail=.true.)
                end if
                l_haschanged = .true.
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
                    call tmp_spproj%read_segment('stk',micproj_records(1)%projname)
                    box_extract = nint(tmp_spproj%os_stk%get(micproj_records(1)%micind,'box'))
                    call init_cluster2D_stream_dev(cline, spproj_glob, box_extract, micspproj_fname, l_cluster2D)
                    call cline%delete('job_memory_per_task2D')
                    call cline%delete('qsys_partition2D')
                    call cline%delete('ncls')
                    call tmp_spproj%kill
                endif
            else
                ! wait & write snapshot
                if( l_cluster2D )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        call update_user_params_dev(cline)
                        call write_migrographs_starfile
                        l_haschanged = .false.
                    endif
                else
                    if( .not.l_projects_left )then
                        if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                            ! write project when inactive
                            call write_project
                            call update_user_params_dev(cline)
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
            if(.not. l_cluster2D .and. cline%defined('dir_meta')) call read_xml_beamtilts(spproj_glob) !!! to update
            ! 2D classification section
            if( l_cluster2D )then
                call update_user_params_dev(cline)
                call update_chunks_dev
                call update_pool_status_dev
                call update_pool_dev
                call update_user_params_dev(cline)
                call reject_from_pool_dev
                call read_pool_xml_beamtilts_dev()
                call assign_pool_optics_dev(cline, propagate = .false.)
                call reject_from_pool_user_dev
                if( .not.l_nchunks_maxed )then
                    call write_project_stream2D_dev(.true.)
                    call import_chunks_into_pool_dev(.false., nchunks_imported)
                    nchunks_imported_glob = nchunks_imported_glob + nchunks_imported
                    l_nchunks_maxed       = nchunks_imported_glob >= params_glob%maxnchunks
                    call classify_pool_dev
                    call update_projects_mask_dev(micproj_records)    ! merge these
                    call classify_new_chunks_dev(micproj_records)     ! two ?
                else
                    ! # of chunks is above desired threshold
                    if( is_pool_available_dev() ) exit
                endif
                call sleep(WAITTIME)
            endif
            ! guistats
            if(file_exists(POOLSTATS_FILE)) call gui_stats%merge(POOLSTATS_FILE)
            call gui_stats%write_json
            call sleep(WAITTIME)
        end do
        ! termination
        if( l_cluster2D )then
            call read_pool_xml_beamtilts_dev()
            call assign_pool_optics_dev(cline, propagate = .true.)
            call terminate_stream2D_dev(.true.)
            if( get_pool_iter_dev() == 0 )then
                ! iteration one was never completed so imported particles & micrographs need be written
                call write_project
            endif
        else
            call write_project
        endif
        call update_user_params(cline)
        call write_migrographs_starfile
        ! final stats
        if(file_exists(POOLSTATS_FILE)) call gui_stats%merge(POOLSTATS_FILE, delete = .true.)
        call gui_stats%hide('micrographs', 'compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT_CLUSTER2D NORMAL STOP ****')
        contains

            !>  write starfile snapshot
            subroutine write_migrographs_starfile
                integer(timer_int_kind)      :: ms0
                real(timer_int_kind)         :: ms_assign, ms_export
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( .not. l_cluster2D ) then
                        if( DEBUG_HERE ) ms0 = tic()
                        call starproj%assign_optics(cline, spproj_glob)
                        if( DEBUG_HERE )then
                            ms_assign = toc(ms0)
                            print *,'ms_assign  : ', ms_assign; call flush(6)
                        endif
                    end if
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%export_mics(cline, spproj_glob)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                    if(allocated(starproj%tiltinfo)) deallocate(starproj%tiltinfo)
                end if
            end subroutine write_migrographs_starfile

            !> For writing the project when the 2D classification has not occured
            subroutine write_project()
                logical, allocatable :: stk_mask(:)
                integer, allocatable :: states(:), fromps(:)
                integer              :: nptcls,istk,fromp,top,i,iptcl,nstks,n,nmics,imic,micind
                character(len=:), allocatable :: prev_projname
                write(logfhandle,'(A)')'>>> PROJECT UPDATE'
                nmics = spproj_glob%os_mic%get_noris()
                call spproj_glob%write_segment_inside('mic', params%projfile)
                if( DEBUG_HERE ) t0 = tic()
                ! stacks
                stk_mask = spproj_glob%os_mic%get_all('nptcls') > 0.5
                states   = nint(spproj_glob%os_mic%get_all('state'))
                nstks    = count(stk_mask)
                allocate(fromps(nmics), source=0)
                call spproj_glob%os_stk%new(nstks, is_ptcl=.false.)
                nptcls = 0
                istk   = 0
                fromp  = 0
                top    = 0
                prev_projname = ''
                do imic = 1,nmics
                    if( .not.stk_mask(imic) ) cycle
                    istk = istk+1
                    if( trim(micproj_records(imic)%projname) /= prev_projname )then
                        call stream_spproj%kill
                        call stream_spproj%read_segment('stk', micproj_records(imic)%projname)
                        prev_projname = trim(micproj_records(imic)%projname)
                    endif
                    micind = micproj_records(imic)%micind
                    call stream_spproj%os_stk%set_state(micind, states(imic))
                    fromps(imic) = nint(stream_spproj%os_stk%get(micind,'fromp'))
                    n            = nint(stream_spproj%os_stk%get(micind,'nptcls'))
                    fromp        = nptcls + 1
                    nptcls       = nptcls + n
                    top          = nptcls
                    call spproj_glob%os_stk%transfer_ori(istk,stream_spproj%os_stk, micind)
                    call spproj_glob%os_stk%set(istk, 'fromp',real(fromp))
                    call spproj_glob%os_stk%set(istk, 'top',  real(top))
                enddo
                call spproj_glob%write_segment_inside('stk', params%projfile)
                call spproj_glob%os_stk%kill
                call stream_spproj%kill
                ! particles
                call spproj_glob%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                istk  = 0
                iptcl = 0
                prev_projname = ''
                do imic = 1,nmics
                    if( .not.stk_mask(imic) ) cycle
                    istk = istk+1
                    if( trim(micproj_records(imic)%projname) /= prev_projname )then
                        call stream_spproj%kill
                        call stream_spproj%read_segment('ptcl2D', micproj_records(imic)%projname)
                        prev_projname = trim(micproj_records(imic)%projname)
                    endif
                    nptcls = micproj_records(imic)%nptcls
                    fromp  = fromps(imic)
                    top    = fromp + nptcls - 1
                    do i = fromp,top
                        iptcl = iptcl + 1
                        call spproj_glob%os_ptcl2D%transfer_ori(iptcl,stream_spproj%os_ptcl2D,i)
                        call spproj_glob%os_ptcl2D%set_stkind(iptcl, istk)
                        call spproj_glob%os_ptcl2D%set_state(iptcl, states(imic))
                    enddo
                enddo
                call stream_spproj%kill
                write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:          ',spproj_glob%os_ptcl2D%get_noris()
                call spproj_glob%write_segment_inside('ptcl2D', params%projfile)
                spproj_glob%os_ptcl3D = spproj_glob%os_ptcl2D
                call spproj_glob%os_ptcl2D%kill
                call spproj_glob%os_ptcl3D%delete_2Dclustering
                call spproj_glob%write_segment_inside('ptcl3D', params%projfile)
                call spproj_glob%os_ptcl3D%kill
                call spproj_glob%write_non_data_segments(params%projfile)
                ! benchmark
                if( DEBUG_HERE )then
                    rt_write = toc(t0)
                    print *,'rt_write  : ', rt_write; call flush(6)
                endif
            end subroutine write_project

            ! updates global project, returns list of processed micrographs
            ! & updates for cluster2D_stream
            subroutine update_projects_list( records, nimported )
                type(micproj_record), allocatable, intent(inout) :: records(:)
                integer,                           intent(inout) :: nimported
                type(sp_project),     allocatable :: spprojs(:)
                type(micproj_record), allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname
                integer :: n_spprojs, n_old, j, nprev_imports, n_completed, nptcls, nmics, imic
                n_completed = 0
                nimported   = 0
                ! previously imported
                n_old = 0 ! on first import
                if( allocated(records) ) n_old = size(records)
                ! projects to import
                n_spprojs = size(completed_jobs_clines)
                if( n_spprojs == 0 )return
                allocate(spprojs(n_spprojs))
                ! because pick_extract purges state=0 and nptcls=0 mics,
                ! all mics can be assumed associated with particles
                nmics = 0
                do iproj = 1,n_spprojs
                    fname = trim(output_dir)//trim(completed_jobs_clines(iproj)%get_carg('projfile'))
                    call spprojs(iproj)%read_segment('mic', fname)
                    nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                enddo
                if( nmics == 0 )then
                    ! nothing to import
                else
                    ! import micrographs
                    n_completed   = n_old + nmics
                    nimported     = nmics
                    nprev_imports = spproj_glob%os_mic%get_noris()
                    ! reallocate global project
                    if( nprev_imports == 0 )then
                        call spproj_glob%os_mic%new(nmics, is_ptcl=.false.) ! first import
                        allocate(micproj_records(nmics))
                    else
                        call spproj_glob%os_mic%reallocate(n_completed)
                        old_records = micproj_records(:)
                        deallocate(micproj_records)
                        allocate(micproj_records(n_completed))
                        if( n_old > 0 ) micproj_records(1:n_old) = old_records(:)
                        deallocate(old_records)
                    endif
                    ! transfer micrographs parameters
                    j = n_old
                    do iproj = 1,n_spprojs
                        do imic = 1,spprojs(iproj)%os_mic%get_noris()
                            j = j + 1
                            nptcls     = nint(spprojs(iproj)%os_mic%get(imic,'nptcls'))
                            nptcls_glob = nptcls_glob + nptcls ! global update
                            fname      = trim(output_dir)//trim(completed_jobs_clines(iproj)%get_carg('projfile'))
                            abs_fname  = simple_abspath(fname, errmsg='pick_extract_cluster2D :: update_projects_list 1')
                            micproj_records(j)%projname = trim(abs_fname)
                            micproj_records(j)%micind   = imic
                            micproj_records(j)%nptcls   = nptcls
                            call spproj_glob%os_mic%transfer_ori(j, spprojs(iproj)%os_mic, imic)
                        enddo
                    enddo
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_projects_list

            subroutine create_individual_project( project_fname, nselected, nrejected )
                character(len=*), intent(in)  :: project_fname
                integer,          intent(out) :: nselected, nrejected
                type(sp_project)              :: tmp_proj, spproj_here
                integer,         allocatable  :: states(:)
                character(len=STDLEN)         :: proj_fname, projname, projfile
                character(len=LONGSTRLEN)     :: path
                integer :: imic, nmics, cnt
                nselected = 0
                nrejected = 0
                call tmp_proj%read_segment('mic', project_fname)
                states    = nint(tmp_proj%os_mic%get_all('state'))
                nmics     = count(states==1)
                nselected = nmics
                nrejected = tmp_proj%os_mic%get_noris() - nselected
                if( nmics == 0 )then
                    call tmp_proj%kill
                    return ! nothing to add to queue
                endif
                ! micrograph rejection
                if( trim(params%reject_mics).eq.'yes' )then
                    do imic = 1,tmp_proj%os_mic%get_noris()
                        if( states(imic) == 0 ) cycle
                        if( tmp_proj%os_mic%isthere(imic, 'ctfres') )then
                            if( tmp_proj%os_mic%get(imic,'ctfres') > (params%ctfresthreshold-0.001) ) states(imic) = 0
                        end if
                        if( states(imic) == 0 ) cycle
                        if( tmp_proj%os_mic%isthere(imic, 'icefrac') )then
                            if( tmp_proj%os_mic%get(imic,'icefrac') > (params%icefracthreshold-0.001) ) states(imic) = 0
                        end if
                    enddo
                    nmics     = count(states==1)
                    nselected = nmics
                    nrejected = tmp_proj%os_mic%get_noris() - nselected
                    if( nmics == 0 )then
                        call tmp_proj%kill
                        return ! nothing to add to queue
                    endif
                endif
                ! as per update_projinfo
                path       = trim(cwd_glob)//'/'//trim(output_dir)
                proj_fname = basename(project_fname)
                projname   = trim(get_fbody(trim(proj_fname), trim(METADATA_EXT), separator=.false.))
                projfile   = trim(projname)//trim(METADATA_EXT)
                call spproj_here%projinfo%new(1, is_ptcl=.false.)
                call spproj_here%projinfo%set(1,'projname', trim(projname))
                call spproj_here%projinfo%set(1,'projfile', trim(projfile))
                call spproj_here%projinfo%set(1,'cwd',      trim(path))
                ! from current global project
                spproj_here%compenv   = spproj_glob%compenv
                spproj_here%jobproc   = spproj_glob%jobproc
                spproj_here%os_optics = spproj_glob%os_optics !??
                call spproj_here%os_mic%new(nmics,is_ptcl=.false.)
                cnt = 0
                do imic = 1,tmp_proj%os_mic%get_noris()
                    if( states(imic) == 0 ) cycle
                    cnt = cnt+1
                    call update_path(tmp_proj%os_mic, imic, 'mc_starfile')
                    call update_path(tmp_proj%os_mic, imic, 'intg')
                    call update_path(tmp_proj%os_mic, imic, 'forctf')
                    call update_path(tmp_proj%os_mic, imic, 'thumb')
                    call update_path(tmp_proj%os_mic, imic, 'mceps')
                    call update_path(tmp_proj%os_mic, imic, 'ctfdoc')
                    call update_path(tmp_proj%os_mic, imic, 'ctfjpg')
                    call spproj_here%os_mic%transfer_ori(cnt, tmp_proj%os_mic, imic)
                enddo
                nselected = cnt
                nrejected = tmp_proj%os_mic%get_noris() - nselected
                ! update for execution
                call cline_pick_extract%set('projname', trim(projname))
                call cline_pick_extract%set('projfile', trim(projfile))
                call cline_pick_extract%set('fromp',    1)
                call cline_pick_extract%set('top',      nselected)
                call spproj_here%write(trim(path)//'/'//trim(projfile))
                call spproj_here%kill
                call tmp_proj%kill
            end subroutine create_individual_project

            !>  import previous run to the current project based on past single project files
            subroutine import_prev_streams
                ! type(sp_project) :: streamspproj
                ! type(ori)        :: o, o_stk
                ! character(len=LONGSTRLEN), allocatable :: sp_files(:)
                ! character(len=:), allocatable :: mic, mov, dir
                ! logical,          allocatable :: spproj_mask(:)
                ! integer :: iproj,nprojs,icnt,nptcls
                ! logical :: err
                ! if( .not.cline%defined('dir_prev') ) return
                ! err = .false.
                ! dir = filepath(params%dir_prev, 'spprojs/')
                ! call simple_list_files_regexp(dir,'^'//trim(PREPROCESS_PREFIX)//'.*\.simple$',sp_files)
                ! if( .not.allocated(sp_files) )then
                !     write(logfhandle,'(A)') '>>> Could not find previously processed movies'
                !     return
                ! endif
                ! nprojs = size(sp_files)
                ! if( nprojs < 1 ) return
                ! allocate(spproj_mask(nprojs),source=.false.)
                ! nptcls = 0
                ! do iproj = 1,nprojs
                !     call streamspproj%read_segment('mic', sp_files(iproj) )
                !     if( streamspproj%os_mic%get_noris() /= 1 )then
                !         THROW_WARN('Ignoring previous project'//trim(sp_files(iproj)))
                !         cycle
                !     endif
                !     if( .not. streamspproj%os_mic%isthere(1,'intg') )cycle
                !     if( l_pick )then
                !         if( streamspproj%os_mic%get(1,'nptcls') < 0.5 )cycle
                !     endif
                !     spproj_mask(iproj) = .true.
                ! enddo
                ! if( count(spproj_mask) == 0 )then
                !     nptcls_glob = 0
                !     return
                ! endif
                ! icnt = 0
                ! do iproj = 1,nprojs
                !     if( .not.spproj_mask(iproj) )cycle
                !     call streamspproj%read_segment('mic',sp_files(iproj))
                !     call streamspproj%os_mic%get_ori(1, o)
                !     ! import mic segment
                !     call movefile2folder('intg',        dir, output_dir_motion_correct, o, err)
                !     call movefile2folder('forctf',      dir, output_dir_motion_correct, o, err)
                !     call movefile2folder('thumb',       dir, output_dir_motion_correct, o, err)
                !     call movefile2folder('mc_starfile', dir, output_dir_motion_correct, o, err)
                !     call movefile2folder('mceps',       dir, output_dir_motion_correct, o, err)
                !     call movefile2folder('ctfjpg',      dir, output_dir_ctf_estimate,   o, err)
                !     call movefile2folder('ctfdoc',      dir, output_dir_ctf_estimate,   o, err)
                !     if( l_pick )then
                !         ! import mic & updates stk segment
                !         call movefile2folder('boxfile', dir, output_dir_picker, o, err)
                !         nptcls = nptcls + nint(o%get('nptcls'))
                !         call streamspproj%os_mic%set_ori(1, o)
                !         if( .not.err )then
                !             call streamspproj%read_segment('stk', sp_files(iproj))
                !             if( streamspproj%os_stk%get_noris() == 1 )then
                !                 call streamspproj%os_stk%get_ori(1, o_stk)
                !                 call movefile2folder('stk', dir, output_dir_extract, o_stk, err)
                !                 call streamspproj%os_stk%set_ori(1, o_stk)
                !                 call streamspproj%read_segment('ptcl2D', sp_files(iproj))
                !                 call streamspproj%read_segment('ptcl3D', sp_files(iproj))
                !             endif
                !         endif
                !     else
                !         ! import mic segment
                !         call streamspproj%os_mic%set_ori(1, o)
                !     endif
                !     ! add to history
                !     call o%getter('movie', mov)
                !     call o%getter('intg', mic)
                !     call movie_buff%add2history(mov)
                !     call movie_buff%add2history(mic)
                !     ! write updated individual project file
                !     call streamspproj%write(trim(dir_preprocess)//basename(sp_files(iproj)))
                !     ! count
                !     icnt = icnt + 1
                ! enddo
                ! if( icnt > 0 )then
                !     ! updating STREAM_SPPROJFILES for Cluster2D_stream
                !     allocate(completed_jobs_clines(icnt))
                !     icnt = 0
                !     do iproj = 1,nprojs
                !         if(spproj_mask(iproj))then
                !             icnt = icnt+1
                !             call completed_jobs_clines(icnt)%set('projfile',basename(sp_files(iproj)))
                !         endif
                !     enddo
                !     call update_projects_list(completed_fnames, n_imported)
                !     deallocate(completed_jobs_clines)
                ! endif
                ! call o%kill
                ! call o_stk%kill
                ! call streamspproj%kill
                ! write(*,'(A,I3)')'>>> IMPORTED PREVIOUS PROCESSED MOVIES: ', icnt
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

            subroutine update_path(os, i, key)
                class(oris),      intent(inout) :: os
                integer,          intent(in)    :: i
                character(len=*), intent(in)    :: key
                character(len=:), allocatable :: fname
                character(len=LONGSTRLEN)     :: newfname
                if( os%isthere(i,key) )then
                    call os%getter(i,key,fname)
                    if( len_trim(fname) > 3 )then
                        if( fname(1:3) == '../' ) fname = trim(params%dir_target)//'/'//trim(fname(4:))
                    endif
                    call make_relativepath(trim(CWD_GLOB)//'/'//trim(DIR_STREAM), fname, newfname)
                    call os%set(i, key, newfname)
                endif
            end subroutine update_path

    end subroutine exec_pick_extract_cluster2D

    ! PRIVATE UTILITIES

    subroutine read_xml_beamtilts( spproj )
        use FoX_dom
        type(sp_project), intent(inout) :: spproj
        type(Node), pointer :: xmldoc, beamtiltnode, beamtiltnodex, beamtiltnodey
        integer :: i
        do i = 1, spproj%os_mic%get_noris()
            if ( is_equal(spproj%os_mic%get(i, "tiltx"), 0.) .and. is_equal(spproj%os_mic%get(i, "tilty"), 0.0)) then
                if(file_exists(spproj%os_mic%get_static(i, "meta"))) then
                    write(logfhandle, *) "stream reading " // trim(adjustl(spproj%os_mic%get_static(i,"meta")))
                    xmldoc => parseFile(trim(adjustl(spproj%os_mic%get_static(i,"meta"))))
                    beamtiltnode  => item(getElementsByTagname(xmldoc, "BeamShift"),0)
                    beamtiltnodex => item(getElementsByTagname(beamtiltnode, "a:_x"), 0)
                    beamtiltnodey => item(getElementsByTagname(beamtiltnode, "a:_y"), 0)
                    call spproj%os_mic%set(i, "tiltx", str2real(getTextContent(beamtiltnodex)))
                    call spproj%os_mic%set(i, "tilty", str2real(getTextContent(beamtiltnodey)))
                    call destroy(xmldoc)
                endif
            endif
        end do
    end subroutine read_xml_beamtilts

    !> updates current parameters with user input
    subroutine update_user_params( cline_here )
        type(cmdline), intent(inout) :: cline_here
        type(oris) :: os
        real       :: tilt_thres
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS) )then
            if( os%isthere(1,'tilt_thres') ) then
                tilt_thres = os%get(1,'tilt_thres')
                if( abs(tilt_thres-params_glob%tilt_thres) > 0.001) then
                     if(tilt_thres < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO LOW: ',tilt_thres
                     else if(tilt_thres > 1) then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO HIGH: ',tilt_thres
                     else
                         params_glob%tilt_thres = tilt_thres
                         call cline_here%set('tilt_thres', params_glob%tilt_thres)
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES UPDATED TO: ',tilt_thres
                     endif
                endif
            endif
            call del_file(USER_PARAMS)
        endif
        call os%kill
    end subroutine update_user_params

    ! Common utility functions

    subroutine update_path(os, i, key)
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
    end subroutine update_path

end module simple_commander_stream
