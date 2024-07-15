! concrete commander: streaming pre-processing routines
module simple_commander_stream
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters, params_glob
use simple_commander_base,     only: commander_base
use simple_sp_project,         only: sp_project
use simple_qsys_env,           only: qsys_env
use simple_starproject_stream, only: starproject_stream
use simple_guistats,           only: guistats
use simple_moviewatcher,       only: moviewatcher
use simple_stream_chunk,       only: micproj_record
use simple_commander_cluster2D_stream_dev
use simple_qsys_funs
use simple_commander_preprocess
use simple_progress
use simple_timer
implicit none

public :: commander_stream_preprocess
public :: commander_stream_gen_picking_refs
public :: commander_stream_pick_extract
public :: commander_stream_assign_optics
public :: commander_stream_cluster2D

private
#include "simple_local_flags.inc"

type, extends(commander_base) :: commander_stream_preprocess
  contains
    procedure :: execute => exec_stream_preprocess
end type commander_stream_preprocess

type, extends(commander_base) :: commander_stream_pick_extract
  contains
    procedure :: execute => exec_stream_pick_extract
end type commander_stream_pick_extract

type, extends(commander_base) :: commander_stream_gen_picking_refs
  contains
    procedure :: execute => exec_stream_gen_picking_refs
end type commander_stream_gen_picking_refs

type, extends(commander_base) :: commander_stream_assign_optics
  contains
    procedure :: execute => exec_stream_assign_optics
end type commander_stream_assign_optics

type, extends(commander_base) :: commander_stream_cluster2D
  contains
    procedure :: execute => exec_stream_cluster2D
end type commander_stream_cluster2D

! module constants
character(len=STDLEN), parameter :: DIR_STREAM           = trim(PATH_HERE)//'spprojs/'           ! location for projects to be processed
character(len=STDLEN), parameter :: DIR_STREAM_COMPLETED = trim(PATH_HERE)//'spprojs_completed/' ! location for projects processed
character(len=STDLEN), parameter :: USER_PARAMS     = 'stream_user_params.txt'                   
integer,               parameter :: NMOVS_SET       = 5                                          ! number of movies processed at once (>1)
integer,               parameter :: LONGTIME        = 60                                         ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME        = 10                                         ! movie folder watched every WAITTIME seconds
integer,               parameter :: SHORTWAIT       = 2                                          ! movie folder watched every SHORTTIME seconds in shmem

contains

    subroutine exec_stream_preprocess( self, cline )
        use simple_motion_correct, only: flip_gain
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
        type(starproject_stream)               :: starproj_stream
        character(len=LONGSTRLEN), allocatable :: movies(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_motion_correct
        character(len=STDLEN)                  :: preproc_nthr_env, preproc_part_env
        integer                                :: movies_set_counter, import_counter
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection, nsets, i,j
        integer                                :: cnt, n_imported, n_added, n_failed_jobs, n_fail_iter, nmic_star, iset, envlen
        logical                                :: l_movies_left, l_haschanged, pause_import
        real                                   :: avg_tmp, preproc_nthr
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('reject_mics', 'no')
        call cline%set('groupframes', 'no')
        if( .not. cline%defined('walltime')         ) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! motion correction
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
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NTHR, preproc_nthr_env, envlen)
        if(envlen > 0) then
            read(preproc_nthr_env,*) preproc_nthr
            call cline%set('nthr', preproc_nthr)
        end if
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params%split_mode = 'stream'
        params%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME, params%dir_movies)
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('movies',      'movies_imported',      int2str(0), primary=.true.)
        call gui_stats%set('movies',      'movies_processed',     int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs',          int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs_rejected', int2str(0), primary=.true.)
        call gui_stats%set('compute',     'compute_in_use',       int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! restart
        movies_set_counter = 0  ! global number of movies set
        import_counter     = 0  ! global import id
        nmic_star          = 0
        if( cline%defined('dir_exec') )then
            call del_file(TERM_STREAM)
            call cline%delete('dir_exec')
            call import_previous_projects
            nmic_star = spproj_glob%os_mic%get_noris()
            call write_mic_star_and_field(write_field=.true.)
            ! guistats
            call gui_stats%set('movies',      'movies_imported',      int2commastr(nmic_star),              primary=.true.)
            call gui_stats%set('movies',      'movies_processed',     int2commastr(nmic_star) // ' (100%)', primary=.true.)
            call gui_stats%set('micrographs', 'micrographs',          int2commastr(nmic_star),              primary=.true.)
            if(spproj_glob%os_mic%isthere("ctfres")) then
                avg_tmp = spproj_glob%os_mic%get_avg("ctfres")
                if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 7.0) then
                    call gui_stats%set('micrographs', 'avg_ctf_resolution', avg_tmp, primary=.true., alert=.true., alerttext='average CTF resolution &
                        &lower than expected for high resolution structure determination', notify=.false.)
                else
                    call gui_stats%set('micrographs', 'avg_ctf_resolution', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                end if
            end if
            if(spproj_glob%os_mic%isthere("icefrac")) then
                avg_tmp = spproj_glob%os_mic%get_avg("icefrac")
                if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 1.0) then
                    call gui_stats%set('micrographs', 'avg_ice_score', avg_tmp, primary=.true., alert=.true., alerttext='average ice score &
                        &greater than expected for high resolution structure determination', notify=.false.)
                else
                    call gui_stats%set('micrographs', 'avg_ice_score', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                end if
            end if
            if(spproj_glob%os_mic%isthere("astig")) then
                avg_tmp = spproj_glob%os_mic%get_avg("astig")
                if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 0.1) then
                    call gui_stats%set('micrographs', 'avg_astigmatism', avg_tmp, primary=.true., alert=.true., alerttext='average astigmatism &
                        &greater than expected for high resolution structure determination', notify=.false.)
                else
                    call gui_stats%set('micrographs', 'avg_astigmatism', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                end if
            end if
            if(spproj_glob%os_mic%isthere('thumb')) then
                call gui_stats%set('latest', '', trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(),'thumb'))), thumbnail=.true.)
            end if
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
        call get_environment_variable(SIMPLE_STREAM_PREPROC_PARTITION, preproc_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(1,stream=.true.,qsys_partition=trim(preproc_part_env))
        else
            call qenv%new(1,stream=.true.)
        end if
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
                nsets = floor(real(nmovies) / real(NMOVS_SET))
                cnt   = 0
                do iset = 1,nsets
                    i = (iset-1)*NMOVS_SET+1
                    j = iset*NMOVS_SET
                    call create_movies_set_project(movies(i:j))
                    call qenv%qscripts%add_to_streaming( cline_exec )
                    do imovie = i,j
                        call movie_buff%add2history( movies(imovie) )
                        cnt     = cnt     + 1
                        n_added = n_added + 1 ! global number of movie sets
                    enddo
                    if( cnt == min(params%nparts*NMOVS_SET,nmovies) ) exit
                enddo
                write(logfhandle,'(A,I4,A,A)')'>>> ',cnt,' NEW MOVIES ADDED; ', cast_time_char(simple_gettime())
                l_movies_left = cnt .ne. nmovies
                ! guistats
                call gui_stats%set('movies', 'movies_imported', int2commastr(movie_buff%n_history), primary=.true.)
                call gui_stats%set_now('movies', 'last_movie_imported')
            else
                l_movies_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
            ! guistats
            call gui_stats%set('compute', 'compute_in_use', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz*NMOVS_SET
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
                call gui_stats%set('movies',      'movies_processed', int2commastr(n_imported) // ' (' // int2str(ceiling(100.0 * real(n_imported) / real(movie_buff%n_history))) // '%)', primary=.true.)
                call gui_stats%set('micrographs', 'micrographs',      int2commastr(n_imported), primary=.true.)
                if( n_failed_jobs > 0 ) call gui_stats%set('micrographs', 'micrographs_rejected', n_failed_jobs, primary=.true.)
                if(spproj_glob%os_mic%isthere("ctfres")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("ctfres")
                    if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 7.0) then
                        call gui_stats%set('micrographs', 'avg_ctf_resolution', avg_tmp, primary=.true., alert=.true., alerttext='average CTF resolution &
                            &lower than expected for high resolution structure determination', notify=.false.)
                    else
                        call gui_stats%set('micrographs', 'avg_ctf_resolution', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                    end if
                end if
                if(spproj_glob%os_mic%isthere("icefrac")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("icefrac")
                    if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 1.0) then
                        call gui_stats%set('micrographs', 'avg_ice_score', avg_tmp, primary=.true., alert=.true., alerttext='average ice score &
                            &greater than expected for high resolution structure determination', notify=.false.)
                    else
                        call gui_stats%set('micrographs', 'avg_ice_score', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                    end if
                end if
                if(spproj_glob%os_mic%isthere("astig")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("astig")
                    if(spproj_glob%os_mic%get_noris() > 50 .and. avg_tmp > 0.1) then
                        call gui_stats%set('micrographs', 'avg_astigmatism', avg_tmp, primary=.true., alert=.true., alerttext='average astigmatism &
                            &greater than expected for high resolution structure determination', notify=.false.)
                    else
                        call gui_stats%set('micrographs', 'avg_astigmatism', avg_tmp, primary=.true., alert=.false., notify=.true., notifytext='tick')
                    end if
                end if
                if(spproj_glob%os_mic%isthere('thumb')) then
                    call gui_stats%set('latest', '', trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(),'thumb'))), thumbnail=.true.)
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                last_injection = simple_gettime()
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
            ! guistats
            call gui_stats%write_json
        end do
        ! termination
        call update_user_params(cline)
        call write_mic_star_and_field(write_field=.true., copy_optics=.true.)
        ! final stats
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_PREPROC NORMAL STOP ****')
        contains

            subroutine write_mic_star_and_field( write_field, copy_optics )
                logical, optional, intent(in) :: write_field, copy_optics
                logical :: l_wfield, l_copy_optics
                l_wfield      = .false.
                l_copy_optics = .false.
                if( present(write_field) ) l_wfield      = write_field
                if( present(copy_optics) ) l_copy_optics = copy_optics
                if(l_copy_optics) then
                    call copy_micrographs_optics()
                    call write_migrographs_starfile(optics_set = .true.)
                else
                    call write_migrographs_starfile
                end if
                if( l_wfield )then
                    call spproj_glob%write_segment_inside('mic', params%projfile)
                    call spproj_glob%write_non_data_segments(params%projfile)
                endif
            end subroutine write_mic_star_and_field

            subroutine copy_micrographs_optics
                integer(timer_int_kind) ::ms0
                real(timer_int_kind)    :: ms_copy_optics
                type(sp_project)        :: spproj_optics
                if( params%projfile_optics .ne. '' .and. file_exists('../' // trim(params%projfile_optics)) ) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call spproj_optics%read('../' // trim(params%projfile_optics))
                    call starproj_stream%copy_optics(spproj_glob, spproj_optics)
                    call spproj_optics%kill()
                    if( DEBUG_HERE )then
                        ms_copy_optics = toc(ms0)
                        print *,'ms_copy_optics  : ', ms_copy_optics; call flush(6)
                    endif
                end if
            end subroutine copy_micrographs_optics

            !>  write starfile snapshot
            subroutine write_migrographs_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ms0
                real(timer_int_kind)          :: ms_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                end if
            end subroutine write_migrographs_starfile

            ! returns list of completed jobs
            subroutine update_projects_list( nimported )
                integer,                   intent(out) :: nimported
                type(sp_project),          allocatable :: streamspprojs(:)
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
                    call streamspprojs(iproj)%read_segment('mic', completed_fnames(iproj))
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
                                ! From now on all MC/CTF metadata use absolute path
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'mc_starfile')
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'intg')
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'thumb')
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'mceps')
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'ctfdoc')
                                call update_relative_path_to_absolute(streamspprojs(iproj)%os_mic, i, 'ctfjpg')
                                ! transfer info
                                call spproj_glob%os_mic%transfer_ori(j, streamspprojs(iproj)%os_mic, i)
                            endif
                        enddo
                        call streamspprojs(iproj)%write_segment_inside('mic', completed_fnames(iproj))
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
                call completed_jobs_clines(:)%kill
                deallocate(completed_jobs_clines,streamspprojs,mics_mask,completed_fnames)
            end subroutine update_projects_list

            subroutine update_relative_path_to_absolute(os, i, key)
                class(oris),      intent(inout) :: os
                integer,          intent(in)    :: i
                character(len=*), intent(in)    :: key
                character(len=:), allocatable :: fname
                character(len=LONGSTRLEN)     :: newfname
                if( os%isthere(i,key) )then
                    call os%getter(i,key,fname)
                    if( fname(1:1) == '/' )then
                        ! already absolute path
                        call os%set(i,key,fname)
                    else
                        ! is relative to ./spprojs
                        newfname = simple_abspath(fname(4:len_trim(fname)))
                        call os%set(i,key,newfname)
                    endif
                endif
            end subroutine update_relative_path_to_absolute

            subroutine create_movies_set_project( movie_names )
                character(len=LONGSTRLEN), intent(in) :: movie_names(NMOVS_SET)
                type(sp_project)             :: spproj_here
                type(ctfparams)              :: ctfvars
                character(len=LONGSTRLEN)    :: projname, projfile, xmlfile, xmldir
                character(len=XLONGSTRLEN)   :: cwd, cwd_old
                integer :: imov
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
                do imov = 1,NMOVS_SET
                    import_counter = import_counter + 1
                    call spproj_here%os_mic%set(imov, "importind", real(import_counter))
                    call spproj_here%os_mic%set(imov, "tiltgrp",   0.0)
                    call spproj_here%os_mic%set(imov, "shiftx",    0.0)
                    call spproj_here%os_mic%set(imov, "shifty",    0.0)
                    call spproj_here%os_mic%set(imov, "flsht",     0.0)
                    if(cline%defined('dir_meta')) then
                        xmldir = cline%get_carg('dir_meta')
                        xmlfile = basename(trim(movie_names(imov)))
                        if(index(xmlfile, '_fractions') > 0) xmlfile = xmlfile(:index(xmlfile, '_fractions') - 1)
                        if(index(xmlfile, '_EER') > 0)       xmlfile = xmlfile(:index(xmlfile, '_EER') - 1)
                        xmlfile = trim(adjustl(xmldir))//'/'//trim(adjustl(xmlfile))//'.xml'
                        call spproj_here%os_mic%set(imov, "meta", trim(adjustl(xmlfile)))
                    end if
                enddo
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
                ! update import id counter
                import_counter = spproj_glob%os_mic%get_noris()
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

    subroutine exec_stream_pick_extract( self, cline )
        use simple_histogram,    only: histogram
        class(commander_stream_pick_extract), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(make_pickrefs_commander)          :: xmake_pickrefs
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time triggers writing of project file
        logical,                   parameter   :: DEBUG_HERE      = .true.
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(micproj_record),      allocatable :: micproj_records(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_make_pickrefs, cline_pick_extract
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob, stream_spproj
        type(starproject_stream)               :: starproj_stream
        type(histogram)                        :: histogram_moldiams
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=:),          allocatable :: odir, odir_extract, odir_picker, odir_completed
        character(len=:),          allocatable :: odir_picker_init, odir_picker_init_completed
        character(len=LONGSTRLEN)              :: cwd_job, latest_boxfile
        character(len=STDLEN)                  :: pick_nthr_env, pick_part_env
        real,                      allocatable :: moldiams(:)
        real                                   :: pick_nthr
        integer                                :: nmics_sel, nmics_rej, nmics_rejected_glob, pick_extract_set_counter
        integer                                :: nmics, nprojects, stacksz, prev_stacksz, iter, last_injection, iproj, envlen
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, nmic_star
        logical                                :: l_templates_provided, l_projects_left, l_haschanged, l_multipick, l_extract, l_once
        logical                                :: l_multipick_init, l_multipick_refine, pause_import
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        call cline%set('picker',  'new')
        if( .not. cline%defined('outdir')          ) call cline%set('outdir',           '')
        if( .not. cline%defined('walltime')        ) call cline%set('walltime',         29.0*60.0) ! 29 minutes
        ! micrograph selection
        if( .not. cline%defined('reject_mics')     ) call cline%set('reject_mics',      'yes')
        if( .not. cline%defined('ctfresthreshold') ) call cline%set('ctfresthreshold',  CTFRES_THRESHOLD_STREAM)
        if( .not. cline%defined('icefracthreshold')) call cline%set('icefracthreshold', ICEFRAC_THRESHOLD_STREAM)
        if( .not. cline%defined('astigthreshold'  )) call cline%set('astigthreshold',   ASTIG_THRESHOLD_STREAM)
        ! picking
        if( .not. cline%defined('lp_pick')         ) call cline%set('lp_pick',          PICK_LP_DEFAULT)
        if( .not. cline%defined('pick_roi')        ) call cline%set('pick_roi',         'no')
        if( .not. cline%defined('backgr_subtr')    ) call cline%set('backgr_subtr',     'no')
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',        'black')
        if( .not. cline%defined('extractfrommov')  ) call cline%set('extractfrommov',   'no')
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_PICK_NTHR, pick_nthr_env, envlen)
        if(envlen > 0) then
            read(pick_nthr_env,*) pick_nthr
            call cline%set('nthr', pick_nthr)
        end if
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params%split_mode = 'stream'
        params%ncunits    = params%nparts
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        ! picking
        l_multipick = cline%defined('nmoldiams')
        l_multipick_init   = .false.
        l_multipick_refine = .false.
        if( l_multipick )then
            l_extract            = .false.
            l_templates_provided = .false.
            l_multipick_refine   = .false.
            l_multipick_init     = cline%defined('ninit')
            write(logfhandle,'(A)')'>>> PERFORMING MULTI-DIAMETER PICKING'
            moldiams = equispaced_vals(params%moldiam, params%moldiam_max, params%nmoldiams)
            call histogram_moldiams%new(moldiams)
            deallocate(moldiams)
            ! remove existing files (restart)
            if(file_exists("micrographs.star")) call del_file("micrographs.star")
            if(file_exists("micrographs_init.star")) call del_file("micrographs_init.star")
            if(file_exists("pick.star")) call del_file("pick.star")
            if(file_exists("pick_init.star")) call del_file("pick_init.star")
        else
            l_extract            = .true.
            l_templates_provided = cline%defined('pickrefs')
            if( l_templates_provided )then
                if( .not.file_exists(params%pickrefs) ) THROW_HARD('Could not find: '//trim(params%pickrefs))
                write(logfhandle,'(A)')'>>> PERFORMING REFERENCE-BASED PICKING'
                if( cline%defined('moldiam') )then
                    call cline%delete('moldiam')
                    write(logfhandle,'(A)')'>>> MOLDIAM IGNORED'
                endif
            else if( .not.cline%defined('moldiam') )then
                THROW_HARD('MOLDIAM required for picker=new reference-free picking')
                write(logfhandle,'(A)')'>>> PERFORMING SINGLE DIAMETER PICKING'
            endif
        endif
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('micrographs', 'micrographs_imported', int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs_rejected', int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs_picked',   int2str(0), primary=.true.)
        if( l_multipick_init ) then
            call gui_stats%set('current_search', 'type', 'initial global search')
            call gui_stats%set('current_search', 'range', int2str(floor(params%moldiam)) // 'Å - ' // int2str(floor(params%moldiam_max)) // 'Å')
            call gui_stats%set('current_search', 'status', 'running')
        endif
        call gui_stats%set('compute',     'compute_in_use',       int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! directories structure & restart
        odir                       = trim(DIR_STREAM)
        odir_completed             = trim(DIR_STREAM_COMPLETED)
        odir_picker                = trim(PATH_HERE)//trim(DIR_PICKER)
        odir_picker_init           = trim(odir_picker)//'init/'
        odir_picker_init_completed = trim(odir_completed)//'init/'
        odir_extract               = trim(PATH_HERE)//trim(DIR_EXTRACT)
        pick_extract_set_counter = 0    ! global counter of projects to be processed
        nptcls_glob              = 0    ! global number of particles
        nmics_rejected_glob      = 0    ! global number of micrographs rejected
        nmic_star                = 0
        if( cline%defined('dir_exec') )then
            call del_file(TERM_STREAM)
            call cline%delete('dir_exec')
            call simple_rmdir(odir)
            if( l_multipick )then
                ! removes directory structure
                call simple_rmdir(odir_completed)
                call simple_rmdir(odir_picker)
                call simple_rmdir(odir_picker_init)
                call simple_rmdir(odir_picker_init_completed)
                call simple_rmdir(odir_extract)
            else
                ! import previous run and updates stats for gui
                call import_previous_mics( micproj_records )
                if( allocated(micproj_records) )then
                    nptcls_glob = sum(micproj_records(:)%nptcls)
                    nmic_star   = spproj_glob%os_mic%get_noris()
                    call gui_stats%set('micrographs', 'micrographs_imported', int2commastr(nmic_star),              primary=.true.)
                    call gui_stats%set('micrographs', 'micrographs_picked',   int2commastr(nmic_star) // ' (100%)', primary=.true.)
                    if( spproj_glob%os_mic%isthere("nptcls") ) then
                        call gui_stats%set('micrographs', 'avg_number_picks', ceiling(spproj_glob%os_mic%get_avg("nptcls")), primary=.true.)
                    end if
                    call gui_stats%set('particles', 'total_extracted_particles', nptcls_glob, primary=.true.)
                    if(spproj_glob%os_mic%isthere('intg') .and. spproj_glob%os_mic%isthere('boxfile')) then
                        latest_boxfile = trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'boxfile'))
                        if(file_exists(trim(latest_boxfile))) call gui_stats%set('latest', '', trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'intg')), thumbnail=.true., boxfile=trim(latest_boxfile))
                    end if
                endif
            endif
        endif
        ! make directories structure
        call simple_mkdir(odir)
        call simple_mkdir(trim(odir)//trim(STDERROUT_DIR))
        call simple_mkdir(odir_completed)
        call simple_mkdir(odir_picker)
        if( l_multipick_init ) then
            call simple_mkdir(odir_picker_init)
            call simple_mkdir(odir_picker_init_completed)
        endif
        if( l_extract ) call simple_mkdir(odir_extract)
        ! initialise progress monitor
        call progressfile_init()
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PICK_PARTITION, pick_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(1,stream=.true.,qsys_partition=trim(pick_part_env))
        else
            call qenv%new(1,stream=.true.)
        end if
        ! command line for execution
        cline_pick_extract = cline
        call cline_pick_extract%set('prg','pick_extract')
        if( l_multipick_init ) then
            call cline_pick_extract%set('dir', filepath(PATH_PARENT,odir_picker_init))
        else
            call cline_pick_extract%set('dir', PATH_PARENT)
        endif
        if( l_extract )then
            call cline_pick_extract%set('extract','yes')
        else
            call cline_pick_extract%set('extract','no')
        endif
        ! ugly single use flag for backwards compatibility, will need to go
        call cline_pick_extract%set('newstream','yes')
        ! Infinite loop
        last_injection        = simple_gettime()
        prev_stacksz          = 0
        iter                  = 0
        n_imported            = 0   ! global number of imported processed micrographs
        n_failed_jobs         = 0
        n_added               = 0   ! global number of micrographs added to processing stack
        l_projects_left       = .false.
        l_haschanged          = .false.
        l_once                = .true.
        pause_import          = .false.
        do
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING STREAM PICK_EXTRACT'
                exit
            endif
            iter = iter + 1
            ! switch to diameter refinement
            if( l_multipick_refine .and. params_glob%updated .eq. 'yes' .and. params_glob%moldiam_refine .gt. 0.0) then
                write(logfhandle,'(A,I3)') '>>> REFINING MOLECULAR DIAMETER        : ', int(params_glob%moldiam_refine)
                !! necessary??
                odir_picker = filepath(trim(PATH_HERE), trim(DIR_PICKER))
                call simple_mkdir(odir_picker, errmsg="commander_stream :: exec_stream_pick_extract;  ")
                !!
                params_glob%updated = 'no'
                params%moldiam      = params_glob%moldiam_refine - 50.0
                params%moldiam_max  = params_glob%moldiam_refine + 50.0
                params%nmoldiams    = 11.0
                call cline_pick_extract%set('moldiam',     params%moldiam)
                call cline_pick_extract%set('moldiam_max', params%moldiam_max)
                call cline_pick_extract%set('nmoldiams',   params%nmoldiams)
                call cline_pick_extract%set('dir','../')
                ! undo already processed micrographs
                call spproj_glob%os_mic%kill()
                call project_buff%clear_history()
                if(allocated(micproj_records)) deallocate(micproj_records)
                ! reset histogram
                moldiams = equispaced_vals(params%moldiam, params%moldiam_max, params%nmoldiams)
                call histogram_moldiams%kill
                call histogram_moldiams%new(moldiams)
                deallocate(moldiams)
                call gui_stats%set('current_search', 'type', 'refinement')
                call gui_stats%set('current_search', 'range', int2str(floor(params%moldiam)) // 'Å - ' // int2str(floor(params%moldiam_max)) // 'Å')
                call gui_stats%set('current_search', 'estimated_diameter', '')
                call gui_stats%set('current_search', 'status', 'running')
                call gui_stats%set('latest', '', '')
                call gui_stats%write_json
                pause_import = .false.
            endif
            ! pause import after ninit mics for global search
            if(l_multipick_init .and. .not. pause_import .and. n_added >= params%ninit) then
                write(logfhandle,'(A,A,A)') '>>> NEW MICROGRAPH IMPORT PAUSED AFTER ', int2str(params%ninit), ' MICROGRAPHS WHILE INITIAL SEARCH IS PERFORMED';
                call gui_stats%set('micrographs', 'micrographs_imported', int2str(project_buff%n_history * NMOVS_SET) // "(paused)", primary=.true.)
                pause_import = .true.
            endif
            ! detection of new projects
            call project_buff%watch( nprojects, projects, max_nmovies=params%nparts )
            ! append projects to processing stack
            if( .not. pause_import .and. nprojects > 0 )then
                cnt   = 0
                nmics = 0
                do iproj = 1, nprojects
                    call create_individual_project(projects(iproj), nmics_sel, nmics_rej)
                    if( l_once )then
                        ! prepares picking references here as we did not have pixel size before
                        if( l_templates_provided )then
                            cline_make_pickrefs = cline
                            call cline_make_pickrefs%set('prg',   'make_pickrefs')
                            call cline_make_pickrefs%set('stream','no')
                            call cline_make_pickrefs%set('smpd',  params%smpd)
                            call xmake_pickrefs%execute_shmem(cline_make_pickrefs)
                            call cline_pick_extract%set('pickrefs', '../'//trim(PICKREFS_FBODY)//trim(params%ext))
                            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
                            call qsys_cleanup
                        endif
                        l_once = .false.
                    endif
                    call project_buff%add2history(projects(iproj))
                    if( nmics_sel > 0 )then
                        call qenv%qscripts%add_to_streaming(cline_pick_extract)
                        call qenv%qscripts%schedule_streaming( qenv%qdescr, path=odir )
                        cnt   = cnt   + 1
                        nmics = nmics + nmics_sel
                    endif
                    n_added             = n_added + nmics_sel
                    nmics_rejected_glob = nmics_rejected_glob + nmics_rej
                    if( cnt == min(params%nparts, nprojects) ) exit
                enddo
                if( nmics > 0 )then
                    write(logfhandle,'(A,I4,A,A)')'>>> ',nmics,' NEW MICROGRAPHS ADDED; ',cast_time_char(simple_gettime())
                endif
                l_projects_left = cnt .ne. nprojects
                ! guistats
                call gui_stats%set('micrographs', 'micrographs_imported', int2commastr(project_buff%n_history * NMOVS_SET), primary=.true.)
                call gui_stats%set_now('micrographs', 'last_micrograph_imported')
            else
                l_projects_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=odir )
            stacksz = qenv%qscripts%get_stacksz()
            ! guistats
            call gui_stats%set('compute', 'compute_in_use', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz                          ! # of projects
                stacksz      = qenv%qscripts%get_stack_range()  ! # of micrographs
                write(logfhandle,'(A,I6)')'>>> MICROGRAPHS TO PROCESS:                 ', stacksz
            endif
            ! fetch completed jobs list & updates
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
                write(logfhandle,'(A,I8)')                '>>> # MICROGRAPHS PROCESSED & IMPORTED  : ',n_imported
                if( l_extract ) write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED               : ',nptcls_glob
                if( l_multipick )then
                    call histogram_moldiams%plot('moldiams')
                    if( l_multipick_init )  call starproj_stream%stream_export_pick_diameters(params%outdir, histogram_moldiams, filename="pick_init.star")
                    if( l_multipick_refine) call starproj_stream%stream_export_pick_diameters(params%outdir, histogram_moldiams)
                    write(logfhandle,'(A,F8.2)') '>>> ESTIMATED MOLECULAR DIAMETER        : ',histogram_moldiams%mean()
                    call gui_stats%set('current_search', 'estimated_diameter', int2str(floor(histogram_moldiams%mean())) // 'Å')
                endif
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                call gui_stats%set('micrographs', 'micrographs_picked', int2commastr(n_imported) // ' (' // int2str(ceiling(100.0 * real(n_imported) / real(project_buff%n_history * NMOVS_SET))) // '%)', primary=.true.)
                call gui_stats%set('micrographs', 'micrographs_rejected', int2commastr(n_failed_jobs + nmics_rejected_glob) // ' (' // int2str(floor(100.0 * real(n_failed_jobs + nmics_rejected_glob) / real(project_buff%n_history * NMOVS_SET))) // '%)', primary=.true.)
                if(spproj_glob%os_mic%isthere("nptcls") .and. .not. l_multipick) then
                    call gui_stats%set('micrographs', 'avg_number_picks', ceiling(spproj_glob%os_mic%get_avg("nptcls")), primary=.true.)
                end if
                if(.not. l_multipick) call gui_stats%set('particles', 'total_extracted_particles', int2commastr(nptcls_glob), primary=.true.)
                if(spproj_glob%os_mic%isthere('intg') .and. spproj_glob%os_mic%isthere('boxfile')) then
                    latest_boxfile = trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'boxfile'))
                    if(file_exists(trim(latest_boxfile))) call gui_stats%set('latest', '', trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'intg')), thumbnail=.true., boxfile=trim(latest_boxfile))
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                last_injection = simple_gettime()
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
            else
                ! write snapshot
                if( .not.l_projects_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive
                        call write_project
                        call update_user_params(cline)
                        call write_migrographs_starfile
                        l_haschanged = .false.
                    endif
                endif
            endif
            if(l_multipick_init .and. pause_import .and. n_imported + n_failed_jobs >= params%ninit) then
                write(logfhandle,'(A)')'>>> WAITING FOR USER INPUT REFINEMENT DIAMETER'
                call gui_stats%set('current_search', 'status', 'waiting for refinement diameter')
                l_multipick_init   = .false.
                l_multipick_refine = .true.
            endif
            call update_user_params(cline)
            ! guistats
            call gui_stats%write_json
            call sleep(WAITTIME)
        end do
        ! termination
        call write_project
        call update_user_params(cline)
        call copy_micrographs_optics
        call write_migrographs_starfile(optics_set=.true.)
        call write_particles_starfile(optics_set=.true.)
        ! kill ptcls now starfiles written
        call spproj_glob%os_stk%kill
        call spproj_glob%os_ptcl2D%kill
        ! final stats
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT NORMAL STOP ****')
        contains

            subroutine copy_micrographs_optics
                integer(timer_int_kind) ::ms0
                real(timer_int_kind)    :: ms_copy_optics
                type(sp_project)        :: spproj_optics
                if( params%projfile_optics .ne. '' .and. file_exists('../' // trim(params%projfile_optics)) ) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call spproj_optics%read('../' // trim(params%projfile_optics))
                    call starproj_stream%copy_optics(spproj_glob, spproj_optics)
                    call spproj_optics%kill()
                    call spproj_glob%write
                    if( DEBUG_HERE )then
                        ms_copy_optics = toc(ms0)
                        print *,'ms_copy_optics  : ', ms_copy_optics; call flush(6)
                    endif
                end if
            end subroutine copy_micrographs_optics

            !>  write starfile snapshot
            subroutine write_migrographs_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ms0
                real(timer_int_kind)          :: ms_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( DEBUG_HERE ) ms0 = tic()
                    if(l_multipick_init) then
                        call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set, filename="micrographs_init.star")
                    else
                        call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set)
                    endif
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                end if
            end subroutine write_migrographs_starfile

            subroutine write_particles_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ptcl0
                real(timer_int_kind)          :: ptcl_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_ptcl2D%get_noris() > 0) then
                    if( DEBUG_HERE ) ptcl0 = tic()
                    call starproj_stream%stream_export_particles_2D(spproj_glob, params%outdir, optics_set=l_optics_set, verbose=.true.)
                    if( DEBUG_HERE )then
                        ptcl_export = toc(ptcl0)
                        print *,'ptcl_export  : ', ptcl_export; call flush(6)
                    endif
                end if
            end subroutine write_particles_starfile

            subroutine write_project()
                integer, allocatable :: fromps(:)
                integer              :: nptcls,fromp,top,i,iptcl,nmics,imic,micind
                character(len=:), allocatable :: prev_projname
                write(logfhandle,'(A)')'>>> PROJECT UPDATE'
                if( DEBUG_HERE ) t0 = tic()
                ! micrographs
                nmics = spproj_glob%os_mic%get_noris()
                call spproj_glob%write_segment_inside('mic', params%projfile)
                if( l_extract )then
                    ! stacks
                    allocate(fromps(nmics), source=0)
                    call spproj_glob%os_stk%new(nmics, is_ptcl=.false.)
                    nptcls        = 0
                    fromp         = 0
                    top           = 0
                    prev_projname = ''
                    do imic = 1,nmics
                        if( trim(micproj_records(imic)%projname) /= trim(prev_projname) )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('stk', micproj_records(imic)%projname)
                            prev_projname = trim(micproj_records(imic)%projname)
                        endif
                        micind = micproj_records(imic)%micind
                        fromps(imic) = nint(stream_spproj%os_stk%get(micind,'fromp')) ! fromp from individual project
                        fromp        = nptcls + 1
                        nptcls       = nptcls + micproj_records(imic)%nptcls
                        top          = nptcls
                        call spproj_glob%os_stk%transfer_ori(imic, stream_spproj%os_stk, micind)
                        call spproj_glob%os_stk%set(imic, 'fromp',real(fromp))
                        call spproj_glob%os_stk%set(imic, 'top',  real(top))
                    enddo
                    call spproj_glob%write_segment_inside('stk', params%projfile)
                    ! particles
                    call spproj_glob%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                    iptcl         = 0
                    prev_projname = ''
                    do imic = 1,nmics
                        if( trim(micproj_records(imic)%projname) /= prev_projname )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('ptcl2D', micproj_records(imic)%projname)
                            prev_projname = trim(micproj_records(imic)%projname)
                        endif
                        fromp = fromps(imic)
                        top   = fromp + micproj_records(imic)%nptcls - 1
                        do i = fromp,top
                            iptcl = iptcl + 1
                            call spproj_glob%os_ptcl2D%transfer_ori(iptcl, stream_spproj%os_ptcl2D, i)
                            call spproj_glob%os_ptcl2D%set_stkind(iptcl, imic)
                        enddo
                    enddo
                    call stream_spproj%kill
                    write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:          ',spproj_glob%os_ptcl2D%get_noris()
                    call spproj_glob%write_segment_inside('ptcl2D', params%projfile)
                    spproj_glob%os_ptcl3D = spproj_glob%os_ptcl2D
                    call spproj_glob%os_ptcl3D%delete_2Dclustering
                    call spproj_glob%write_segment_inside('ptcl3D', params%projfile)
                    call spproj_glob%os_ptcl3D%kill
                endif
                call spproj_glob%write_non_data_segments(params%projfile)
                ! benchmark
                if( DEBUG_HERE )then
                    rt_write = toc(t0)
                    print *,'rt_write  : ', rt_write; call flush(6)
                endif
            end subroutine write_project

            ! updates global project, returns records of processed micrographs
            subroutine update_projects_list( records, nimported )
                type(micproj_record), allocatable, intent(inout) :: records(:)
                integer,                           intent(inout) :: nimported
                type(sp_project),     allocatable :: spprojs(:)
                type(micproj_record), allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname, new_fname
                type(sp_project)                  :: tmpproj
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
                    fname = filepath(odir, completed_jobs_clines(iproj)%get_carg('projfile'))
                    call tmpproj%read_segment('projinfo', fname)
                    if(l_multipick_refine .and. tmpproj%projinfo%isthere("init") .and. tmpproj%projinfo%get(1, "init") .eq. 1.0) then
                        if(file_exists(trim(fname))) call del_file(trim(fname))
                        call spprojs(iproj)%os_mic%new(0, .false.)
                    else
                        call spprojs(iproj)%read_segment('mic', fname)
                       ! call spprojs(iproj)%read_segment('projinfo', fname)
                        nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                    endif
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
                        call move_alloc(micproj_records, old_records)
                        allocate(micproj_records(n_completed))
                        if( n_old > 0 ) micproj_records(1:n_old) = old_records(:)
                        deallocate(old_records)
                    endif
                    ! update records and global project
                    j = n_old
                    do iproj = 1,n_spprojs
                        if( spprojs(iproj)%os_mic%get_noris() == 0 ) cycle
                        ! move project to appropriate directory
                        fname = filepath(odir, completed_jobs_clines(iproj)%get_carg('projfile'))
                        if(tmpproj%projinfo%isthere("init") .and. tmpproj%projinfo%get(1, "init") .eq. 1.0) then
                            new_fname = filepath(odir_picker_init_completed, completed_jobs_clines(iproj)%get_carg('projfile'))
                        else
                            new_fname = filepath(odir_completed, completed_jobs_clines(iproj)%get_carg('projfile'))
                        endif    
                        call simple_rename(fname, new_fname)
                        abs_fname = simple_abspath(new_fname, errmsg='stream pick_extract :: update_projects_list 1')
                        ! records & project
                        do imic = 1,spprojs(iproj)%os_mic%get_noris()
                            j = j + 1
                            micproj_records(j)%projname = trim(abs_fname)
                            micproj_records(j)%micind   = imic
                            if( l_multipick )then
                                micproj_records(j)%nptcls = 0
                            else
                                nptcls      = nint(spprojs(iproj)%os_mic%get(imic,'nptcls'))
                                nptcls_glob = nptcls_glob + nptcls ! global update
                                micproj_records(j)%nptcls = nptcls
                            endif
                            call spproj_glob%os_mic%transfer_ori(j, spprojs(iproj)%os_mic, imic)
                        enddo
                    enddo
                    if( l_multipick )then
                        ! update molecular diameters
                        do j = n_old+1,spproj_glob%os_mic%get_noris()
                            if( spproj_glob%os_mic%isthere(j, 'moldiam') )then
                                call histogram_moldiams%update(spproj_glob%os_mic%get(j, 'moldiam'))
                            endif
                        enddo
                    endif
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
                call tmpproj%kill
            end subroutine update_projects_list

            ! prepares project for processing and performs micrograph selection
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
                allocate(states(tmp_proj%os_mic%get_noris()), source=nint(tmp_proj%os_mic%get_all('state')))
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
                        if( states(imic) == 0 ) cycle
                        if( tmp_proj%os_mic%isthere(imic, 'astig') )then
                            if( tmp_proj%os_mic%get(imic,'astig') > (params%astigthreshold-0.001) ) states(imic) = 0
                        end if
                    enddo
                    nmics     = count(states==1)
                    nselected = nmics
                    nrejected = tmp_proj%os_mic%get_noris() - nselected
                    if( nselected == 0 )then
                        call tmp_proj%kill
                        return ! nothing to add to queue
                    endif
                endif
                ! as per update_projinfo
                path       = trim(cwd_glob)//'/'//trim(odir)
                proj_fname = basename(project_fname)
                projname   = trim(get_fbody(trim(proj_fname), trim(METADATA_EXT), separator=.false.))
                projfile   = trim(projname)//trim(METADATA_EXT)
                call spproj_here%projinfo%new(1, is_ptcl=.false.)
                call spproj_here%projinfo%set(1,'projname', trim(projname))
                call spproj_here%projinfo%set(1,'projfile', trim(projfile))
                call spproj_here%projinfo%set(1,'cwd',      trim(path))
                if(l_multipick_init) call spproj_here%projinfo%set(1,'init', 1.0)
                ! from current global project
                spproj_here%compenv = spproj_glob%compenv
                spproj_here%jobproc = spproj_glob%jobproc
                ! import micrographs & updates path to files
                call spproj_here%os_mic%new(nmics,is_ptcl=.false.)
                cnt = 0
                do imic = 1,tmp_proj%os_mic%get_noris()
                    if( states(imic) == 0 ) cycle
                    cnt = cnt+1
                    call spproj_here%os_mic%transfer_ori(cnt, tmp_proj%os_mic, imic)
                enddo
                nselected = cnt
                nrejected = tmp_proj%os_mic%get_noris() - nselected
                ! Gather pixel size once
                if( l_once ) params%smpd = spproj_here%os_mic%get(1,'smpd')
                ! update for execution
                pick_extract_set_counter = pick_extract_set_counter + 1
                projname = int2str_pad(pick_extract_set_counter,params%numlen)
                projfile = trim(projname)//trim(METADATA_EXT)
                call cline_pick_extract%set('projname', trim(projname))
                call cline_pick_extract%set('projfile', trim(projfile))
                call cline_pick_extract%set('fromp',    1)
                call cline_pick_extract%set('top',      nselected)
                call spproj_here%write(trim(path)//'/'//trim(projfile))
                call spproj_here%kill
                call tmp_proj%kill
            end subroutine create_individual_project

            !>  import previous run to the current project and reselect micrographs
            subroutine import_previous_mics( records )
                type(micproj_record),      allocatable, intent(inout) :: records(:)
                type(sp_project),          allocatable :: spprojs(:)
                character(len=LONGSTRLEN), allocatable :: completed_fnames(:)
                character(len=:),          allocatable :: fname
                logical,                   allocatable :: mics_mask(:)
                integer :: n_spprojs, iproj, nmics, imic, jmic, iostat,id, nsel_mics, irec
                pick_extract_set_counter = 0
                ! previously completed projects
                call simple_list_files_regexp(odir_completed, '\.simple$', completed_fnames)
                if( .not.allocated(completed_fnames) )then
                    return ! nothing was previously completed
                endif
                n_spprojs = size(completed_fnames)
                allocate(spprojs(n_spprojs))
                ! read projects micrographs
                nmics = 0
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%read_segment('mic', completed_fnames(iproj))
                    nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                enddo
                ! selection, because pick_extract purges state=0 and nptcls=0 mics,
                ! all mics can be assumed associated with particles
                allocate(mics_mask(nmics))
                jmic = 0
                do iproj = 1,n_spprojs
                    do imic = 1, spprojs(iproj)%os_mic%get_noris()
                        jmic            = jmic+1
                        mics_mask(jmic) = .true.
                        if( spprojs(iproj)%os_mic%isthere(imic, 'ctfres') )then
                            if( spprojs(iproj)%os_mic%get(imic,'ctfres') > (params%ctfresthreshold-0.001) ) mics_mask(jmic) = .false.
                        end if
                        if( .not.mics_mask(jmic) ) cycle
                        if( spprojs(iproj)%os_mic%isthere(imic, 'icefrac') )then
                            if( spprojs(iproj)%os_mic%get(imic,'icefrac') > (params%icefracthreshold-0.001) ) mics_mask(jmic) = .false.
                        end if
                        if( .not.mics_mask(jmic) ) cycle
                        if( spprojs(iproj)%os_mic%isthere(imic, 'astig') )then
                            if( spprojs(iproj)%os_mic%get(imic,'astig') > (params%astigthreshold-0.001) ) mics_mask(jmic) = .false.
                        end if
                    enddo
                enddo
                nsel_mics = count(mics_mask)
                if( nsel_mics == 0 )then
                    ! nothing to import
                    do iproj = 1,n_spprojs
                        call spprojs(iproj)%kill
                        fname = filepath(odir_completed, completed_fnames(iproj))
                        call del_file(fname)
                    enddo
                    deallocate(spprojs)
                    nmics_rejected_glob = nmics
                    return
                endif
                ! updates global records & project
                allocate(records(nsel_mics))
                call spproj_glob%os_mic%new(nsel_mics, is_ptcl=.false.)
                irec = 0
                jmic = 0
                do iproj = 1,n_spprojs
                    do imic = 1, spprojs(iproj)%os_mic%get_noris()
                        jmic = jmic+1
                        if( mics_mask(jmic) )then
                            irec = irec + 1
                            records(irec)%projname = trim(completed_fnames(iproj))
                            records(irec)%micind   = imic
                            records(irec)%nptcls   = nint(spprojs(iproj)%os_mic%get(imic, 'nptcls'))
                            call spproj_glob%os_mic%transfer_ori(irec, spprojs(iproj)%os_mic, imic)
                        endif
                    enddo
                    call spprojs(iproj)%kill
                enddo
                ! update global set counter
                do iproj = 1,n_spprojs
                    fname = basename_safe(completed_fnames(iproj))
                    fname = trim(get_fbody(trim(fname),trim(METADATA_EXT),separator=.false.))
                    call str2int(fname, iostat, id)
                    if( iostat==0 ) pick_extract_set_counter = max(pick_extract_set_counter, id)
                enddo
                nmics_rejected_glob = nmics - nsel_mics
                ! add previous projects to history
                do iproj = 1,n_spprojs
                    call project_buff%add2history(completed_fnames(iproj))
                enddo
                write(logfhandle,'(A,I6,A)')'>>> IMPORTED ',nsel_mics,' PREVIOUSLY PROCESSED MICROGRAPHS'
            end subroutine import_previous_mics

    end subroutine exec_stream_pick_extract

    subroutine exec_stream_gen_picking_refs( self, cline )
        class(commander_stream_gen_picking_refs), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        integer,                   parameter   :: MAXPOP_DEFAULT  = 200000 ! # of particles after wich picking is stopped
        logical,                   parameter   :: DEBUG_HERE      = .FALSE.
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(micproj_record),      allocatable :: micproj_records(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_extract, cline_pick_extract
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob, stream_spproj
        type(starproject_stream)               :: starproj_stream
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=:),          allocatable :: output_dir, output_dir_extract, output_dir_picker
        character(len=LONGSTRLEN)              :: cwd_job, latest_boxfile
        character(len=STDLEN)                  :: refgen_nthr_env, refgen_part_env
        real                                   :: refgen_nthr
        integer                                :: extract_set_counter    ! Internal counter of projects to be processed
        integer                                :: nmics_sel, nmics_rej, nmics_rejected_glob
        integer                                :: nmics, nprojects, stacksz, prev_stacksz, iter, iproj, envlen
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, nmic_star
        logical                                :: l_haschanged, l_once, l_pick_extract
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('picker',       'new')
        call cline%set('reject_cls',   'no')
        call cline%set('autoscale',    'yes')
        call cline%set('kweight_chunk','default')
        call cline%set('kweight_pool', 'default')
        call cline%set('prune',        'no')
        call cline%set('wiener',       'full')
        if( .not. cline%defined('outdir')          ) call cline%set('outdir',           '')
        if( .not. cline%defined('walltime')        ) call cline%set('walltime',         29*60) ! 29 minutes
        ! micrograph selection
        if( .not. cline%defined('reject_mics')     ) call cline%set('reject_mics',      'yes')
        if( .not. cline%defined('ctfresthreshold') ) call cline%set('ctfresthreshold',  CTFRES_THRESHOLD_STREAM)
        if( .not. cline%defined('icefracthreshold')) call cline%set('icefracthreshold', ICEFRAC_THRESHOLD_STREAM)
        if( .not. cline%defined('astigthreshold'  )) call cline%set('astigthreshold',   ASTIG_THRESHOLD_STREAM)
        ! extraction
        if( .not. cline%defined('backgr_subtr')    ) call cline%set('backgr_subtr',     'no')
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',        'black')
        if( .not. cline%defined('extractfrommov')  ) call cline%set('extractfrommov',   'no')
        ! 2D classification
        if( .not. cline%defined('ncls')            ) call cline%set('ncls',           30)
        if( .not. cline%defined('nptcls_per_cls')  ) call cline%set('nptcls_per_cls', 500)
        if( .not. cline%defined('ml_reg')          ) call cline%set('ml_reg',         'no')
        if( .not. cline%defined('refine')          ) call cline%set('refine',         'snhc_smpl')
        if( .not. cline%defined('maxpop')          ) call cline%set('maxpop',         MAXPOP_DEFAULT)
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_REFGEN_NTHR, refgen_nthr_env, envlen)
        if(envlen > 0) then
            read(refgen_nthr_env,*) refgen_nthr
            call cline%set('nthr', refgen_nthr)
            call cline%set('nthr2D', refgen_nthr)
        end if
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params%split_mode = 'stream'
        params%ncunits    = params%nparts
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        call cline%delete('maxpop')
        if( params%nptcls_per_cls*params%ncls > params%maxpop )then
            THROW_HARD('Incorrect inputs: MAXPOP < NPTCLS_PER_CLS x NCLS')
        endif
        ! determine whether we are picking/extracting or extracting only
        if( .not.dir_exists(trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED)) )then
            THROW_HARD('Invalid DIR_TARGET 1')
        endif
        l_pick_extract = .false.
        if( dir_exists(trim(params%dir_target)//'/'//trim(DIR_MOTION_CORRECT)) )then
            l_pick_extract = .true.
        else
            if( .not.dir_exists(trim(params%dir_target)//'/'//trim(DIR_PICKER)) )then
                THROW_HARD('Invalid DIR_TARGET 2')
            endif
        endif
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        ! force local if nparts < 2
        if(params%nparts < 2) then
            call spproj_glob%compenv%set(1, 'qsys_name', 'local')
            call spproj_glob%write_non_data_segments(trim(params%projfile))
        endif
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! 2D parameters
        params%nchunks     = 0
        params%ncls_start  = params%ncls
        params%ml_reg_pool = params%ml_reg
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('particles', 'particles_extracted',         0,            primary=.true.)
        call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        call gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! directories
        output_dir         = trim(PATH_HERE)//trim(DIR_STREAM)
        output_dir_picker  = filepath(trim(PATH_HERE), trim(DIR_PICKER))
        output_dir_extract = filepath(trim(PATH_HERE), trim(DIR_EXTRACT))
        ! restart
        extract_set_counter = 0
        nptcls_glob         = 0     ! global number of particles
        nmics_rejected_glob = 0     ! global number of micrographs rejected
        nmic_star           = 0
        if( cline%defined('dir_exec') )then
            ! no restart here, the folder is wiped
            call cline%delete('dir_exec')
            call cleanup_root_folder(all=.true.)
            call del_file(micspproj_fname)
            call simple_rmdir(output_dir)
            call simple_rmdir(STDERROUT_DIR)
            call simple_rmdir(trim(output_dir)//trim(STDERROUT_DIR))
            call simple_rmdir(trim(PATH_HERE)//trim(DIR_STREAM_COMPLETED))
            call simple_rmdir(output_dir_picker)
            call simple_rmdir(output_dir_extract)
        endif
        ! output directories
        call simple_mkdir(output_dir)
        call simple_mkdir(trim(output_dir)//trim(STDERROUT_DIR))
        call simple_mkdir(trim(PATH_HERE)//trim(DIR_STREAM_COMPLETED))
        call simple_mkdir(output_dir_picker, errmsg="commander_stream :: exec_stream_pick_extract;  ")
        call simple_mkdir(output_dir_extract,errmsg="commander_stream :: exec_stream_pick_extract;  ")
        ! initialise progress monitor
        call progressfile_init()
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_REFGEN_PARTITION, refgen_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(1,stream=.true.,qsys_partition=trim(refgen_part_env))
        else
            call qenv%new(1,stream=.true.)
        end if
        ! command line for execution
        if( l_pick_extract )then
            cline_pick_extract = cline
            call cline_pick_extract%set('prg','pick_extract')
            call cline_pick_extract%set('dir','../')
            call cline_pick_extract%set('extract','yes')
            if( cline%defined('box_extract') )then
                params%box = params%box_extract
                call cline_pick_extract%delete('box_extract')
                call cline_pick_extract%set('box', params%box)
            else
                params%box = 0
            endif
            ! ugly single use flag for backwards compatibility, will need to go
            call cline_pick_extract%set('newstream','yes')
        else
            cline_extract = cline
            call cline_extract%set('prg','extract')
            call cline_extract%set('dir','../'//trim(trim(DIR_EXTRACT)))
            call cline_extract%set('extract','yes')
            if( cline%defined('box_extract') )then
                params%box = params%box_extract
                call cline_extract%delete('box_extract')
                call cline_extract%set('box', params%box)
            else
                params%box = 0
            endif
            ! ugly single use flag for backwards compatibility, will need to go
            call cline_extract%set('newstream','yes')
        endif
        ! Infinite loop
        prev_stacksz          = 0
        iter                  = 0
        n_imported            = 0   ! global number of imported processed micrographs
        n_failed_jobs         = 0
        n_added               = 0   ! global number of micrographs added to processing stack
        l_haschanged          = .false.
        l_once                = .true.
        do
            if( file_exists(trim(TERM_STREAM)) .or. file_exists(STREAM_REJECT_CLS))then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING STREAM GEN_PICKING_REFERENCES'
                exit
            endif
            iter = iter + 1
            ! detection of new projects
            if( nptcls_glob > params%maxpop .and. project_buff%does_exist() )then
                ! limit reached
                call project_buff%kill
                write(logfhandle,'(A)')'>>> PICKING/EXTRACTION OF PARTICLES WILL STOP FROM NOW'
            endif
            call project_buff%watch( nprojects, projects, max_nmovies=params%nparts )
            ! append projects to processing stack
            if( nprojects > 0 )then
                cnt   = 0
                nmics = 0
                do iproj = 1, nprojects
                    call create_individual_project(projects(iproj), nmics_sel, nmics_rej)
                    call project_buff%add2history(projects(iproj))
                    if( nmics_sel > 0 )then
                        if( l_pick_extract )then
                            call qenv%qscripts%add_to_streaming(cline_pick_extract)
                        else
                            call qenv%qscripts%add_to_streaming(cline_extract)
                        endif
                        call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
                        cnt   = cnt   + 1
                        nmics = nmics + nmics_sel
                    endif
                    n_added             = n_added + nmics_sel
                    nmics_rejected_glob = nmics_rejected_glob + nmics_rej
                    if( cnt == min(params%nparts, nprojects) ) exit
                enddo
                if( nmics > 0 )then
                    write(logfhandle,'(A,I4,A,A)')'>>> ',nmics,' NEW MICROGRAPHS ADDED; ',cast_time_char(simple_gettime())
                endif
                ! guistats 
                call gui_stats%set_now('particles', 'last_particles_imported')
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
            ! guistats
            call gui_stats%set('compute', 'compute_in_use', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz                          ! # of projects
                stacksz      = qenv%qscripts%get_stack_range()  ! # of micrographs
                write(logfhandle,'(A,I6)')'>>> MICROGRAPHS TO PROCESS:                 ', stacksz
            endif
            ! fetch completed jobs list & updates
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                call update_records_with_project( micproj_records, n_imported )
                call completed_jobs_clines(:)%kill
                deallocate(completed_jobs_clines)
            else
                n_imported = 0
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
            ! ! project update
            if( n_imported > 0 )then
                n_imported = spproj_glob%os_mic%get_noris()
                write(logfhandle,'(A,I8)')       '>>> # MICROGRAPHS PROCESSED & IMPORTED  : ',n_imported
                write(logfhandle,'(A,I8)')       '>>> # PARTICLES EXTRACTED               : ',nptcls_glob
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                call gui_stats%set('particles', 'particles_extracted', int2commastr(nptcls_glob), primary=.true.)
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
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
            endif
            ! 2D section
            if( l_once .and. (nptcls_glob > params%nptcls_per_cls*params%ncls) )then
                if( .not.cline%defined('mskdiam') ) params%mskdiam = 0.85 * real(params%box) * params%smpd
                call init_cluster2D_stream_dev( cline, spproj_glob, params%box, micspproj_fname, reference_generation=.true. )
                l_once = .false.
            endif
            call update_pool_status_dev
            call update_pool_dev
            call import_records_into_pool( micproj_records )
            call classify_pool_dev
            if(  params%nparts > 1 ) then
                call sleep(WAITTIME)
            else
                call sleep(SHORTWAIT)
            endif
            ! guistats
            if(file_exists(POOLSTATS_FILE)) then
                call gui_stats%merge(POOLSTATS_FILE)
            end if
            call gui_stats%write_json
        end do
        ! termination
        if( l_once )then
            ! nothing to write
        else
            call terminate_stream2D_dev
        endif
        if(file_exists(STREAM_REJECT_CLS)) call write_pool_cls_selected_user_dev
        call gui_stats%delete('latest', '')
        call gui_stats%set('selected references', '', trim(cwd_glob) // '/' // STREAM_SELECTED_REFS // trim(JPG_EXT), thumbnail = .true.)
        !call write_project
        !call update_user_params(cline)
        ! call copy_micrographs_optics
        ! call write_migrographs_starfile(optics_set=.true.)
        ! call write_particles_starfile(optics_set=.true.)
        ! final stats
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_GENERATE_PICKING_REFERENCES NORMAL STOP ****')
        contains

            subroutine copy_micrographs_optics
                integer(timer_int_kind) ::ms0
                real(timer_int_kind)    :: ms_copy_optics
                type(sp_project)        :: spproj_optics
                if( params%projfile_optics .ne. '' .and. file_exists('../' // trim(params%projfile_optics)) ) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call spproj_optics%read('../' // trim(params%projfile_optics))
                    call starproj_stream%copy_optics(spproj_glob, spproj_optics)
                    call spproj_optics%kill()
                    call spproj_glob%write
                    if( DEBUG_HERE )then
                        ms_copy_optics = toc(ms0)
                        print *,'ms_copy_optics  : ', ms_copy_optics; call flush(6)
                    endif
                end if
            end subroutine copy_micrographs_optics

            !>  write starfile snapshot
            subroutine write_migrographs_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ms0
                real(timer_int_kind)          :: ms_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                end if
            end subroutine write_migrographs_starfile

            subroutine write_particles_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ptcl0
                real(timer_int_kind)          :: ptcl_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_ptcl2D%get_noris() > 0) then
                    if( DEBUG_HERE ) ptcl0 = tic()
                    call starproj_stream%stream_export_particles_2D(spproj_glob, params%outdir, optics_set=l_optics_set, verbose=.true.)
                    if( DEBUG_HERE )then
                        ptcl_export = toc(ptcl0)
                        print *,'ptcl_export  : ', ptcl_export; call flush(6)
                    endif
                end if
            end subroutine write_particles_starfile

            ! updates global project, returns records of processed micrographs
            subroutine update_records_with_project( records, nimported )
                type(micproj_record), allocatable, intent(inout) :: records(:)
                integer,                           intent(inout) :: nimported
                type(sp_project),     allocatable :: spprojs(:)
                type(micproj_record), allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname, new_fname
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
                ! because extract purges state=0 and nptcls=0 mics,
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
                        call move_alloc(micproj_records, old_records)
                        allocate(micproj_records(n_completed))
                        if( n_old > 0 ) micproj_records(1:n_old) = old_records(:)
                        deallocate(old_records)
                    endif
                    ! update records and global project
                    j = n_old
                    do iproj = 1,n_spprojs
                        if( spprojs(iproj)%os_mic%get_noris() == 0 ) cycle
                        ! move project to appropriate directory
                        fname     = trim(output_dir)//trim(completed_jobs_clines(iproj)%get_carg('projfile'))
                        new_fname = trim(DIR_STREAM_COMPLETED)//trim(completed_jobs_clines(iproj)%get_carg('projfile'))
                        call simple_rename(fname, new_fname)
                        abs_fname = simple_abspath(new_fname, errmsg='stream gen_picking_refs :: update_projects_list 1')
                        ! records & project
                        do imic = 1,spprojs(iproj)%os_mic%get_noris()
                            j = j + 1
                            micproj_records(j)%projname = trim(abs_fname)
                            micproj_records(j)%micind   = imic
                            nptcls      = nint(spprojs(iproj)%os_mic%get(imic,'nptcls'))
                            nptcls_glob = nptcls_glob + nptcls ! global update
                            micproj_records(j)%nptcls     = nptcls
                            micproj_records(j)%nptcls_sel = nptcls
                            call spproj_glob%os_mic%transfer_ori(j, spprojs(iproj)%os_mic, imic)
                        enddo
                    enddo
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_records_with_project

            ! prepares project for processing and performs micrograph selection
            subroutine create_individual_project( project_fname, nselected, nrejected )
                use simple_nrtxtfile, only: nrtxtfile
                character(len=*), intent(in)  :: project_fname
                integer,          intent(out) :: nselected, nrejected
                type(sp_project)              :: tmp_proj, spproj_here
                type(nrtxtfile)               :: boxfile_multi, boxfile
                real,            allocatable  :: entries(:,:), dists(:)
                integer,         allocatable  :: states(:)
                logical,         allocatable  :: mask(:)
                character(len=STDLEN)         :: proj_fname, projname, projfile
                character(len=LONGSTRLEN)     :: path, boxfname
                real    :: selected_moldiam
                integer :: imic, nmics, cnt, l, nl, nb, i,j,ind, vals(5)
                nselected = 0
                nrejected = 0
                call tmp_proj%read_segment('mic', project_fname)
                allocate(states(tmp_proj%os_mic%get_noris()), source=nint(tmp_proj%os_mic%get_all('state')))
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
                        if( states(imic) == 0 ) cycle
                        if( tmp_proj%os_mic%isthere(imic, 'astig') )then
                            if( tmp_proj%os_mic%get(imic,'astig') > (params%astigthreshold-0.001) ) states(imic) = 0
                        end if
                    enddo
                    nmics     = count(states==1)
                    nselected = nmics
                    nrejected = tmp_proj%os_mic%get_noris() - nselected
                    if( nselected == 0 )then
                        call tmp_proj%kill
                        return ! nothing to add to queue
                    endif
                endif
                if( l_pick_extract )then
                    do imic = 1,tmp_proj%os_mic%get_noris()
                        if( states(imic) == 1 )then
                            params%smpd = tmp_proj%os_mic%get(imic, 'smpd')
                            exit
                        endif
                    enddo
                else
                    ! select diameter, multi moldiam boxfile format: pick_gau%report_multipick
                    if( .not.tmp_proj%os_mic%isthere('boxfile') ) return
                    do imic = 1,tmp_proj%os_mic%get_noris()
                        if( states(imic) /= 1 ) cycle
                        boxfname = trim(tmp_proj%os_mic%get_static(imic,'boxfile'))
                        if( .not. file_exists(boxfname) )then
                            THROW_WARN('FILE NOT FOUND: '//trim(boxfname))
                            cycle
                        endif
                        call boxfile_multi%new(boxfname, 1)
                        if( boxfile_multi%get_nrecs_per_line() /= 6 )then
                            THROW_WARN('INVALID FORMAT FOR: '//trim(boxfname))
                            states(imic) = 0
                            call boxfile_multi%kill
                            cycle
                        endif
                        nl = boxfile_multi%get_ndatalines()
                        if(allocated(entries)) deallocate(entries,dists)
                        allocate(entries(6,nl),dists(nl),source=0.)
                        do l = 1,nl
                            call boxfile_multi%readNextDataLine(entries(:,l))
                        enddo
                        call boxfile_multi%kill
                        dists = entries(4,:) - params%moldiam
                        if( any(dists > 0.) )then
                            selected_moldiam = minval(entries(4,:), mask=dists>-0.001)
                        else
                            selected_moldiam = maxval(entries(4,:))
                        endif
                        allocate(mask(nl),source=.false.)
                        !Cyril check
                        if(params%box == 0) then
                            mask = abs(entries(4,:)-selected_moldiam) < 0.001
                        else
                            mask = abs(entries(3,:)-params%box) < 0.001
                        endif
                        nb   = count(mask)
                        if( nb < 1 )then
                            states(imic) = 0
                        else
                            ! update to global extraction box size
                            if( .not.cline%defined('box_extract') )then
                                ind = 0
                                do j = 1,nl
                                    if( mask(j) )then
                                        ind = j
                                        exit
                                    endif
                                enddo
                                if( ind == 0 ) THROW_HARD('BOX SIZE ERROR')
                                if( params%box==0 )then
                                    params%box  = nint(entries(3,ind)) ! first time
                                else
                                    if( nint(entries(3,ind)) /= params%box )then
                                        THROW_HARD('INCONSISTENT EXTRACTION BOX SIZES '//trim(boxfname)//' '//int2str(params%box)//' '//int2str(nint(entries(3,ind))))
                                    endif
                                endif
                            endif
                            params%smpd = tmp_proj%os_mic%get(imic, 'smpd')
                            ! write new boxfile
                            path = trim(cwd_glob)//'/'//trim(output_dir_picker)//'/'//basename(boxfname)
                            call boxfile%new(path, 2, wanted_recs_per_line=5)
                            vals(5) = -3
                            do i = 1,nl
                                if( mask(i) )then
                                    vals(1:2) = nint(entries(1:2,i))
                                    vals(3:4) = nint(entries(3,i))
                                    call boxfile%write(vals)
                                endif
                            enddo
                            call boxfile%kill
                            call tmp_proj%os_mic%set(imic,'boxfile',path)
                        endif
                        deallocate(mask)
                    enddo
                endif
                if( all(states==0) )then
                    call tmp_proj%kill
                    return
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
                spproj_here%compenv = spproj_glob%compenv
                spproj_here%jobproc = spproj_glob%jobproc
                ! import
                call spproj_here%os_mic%new(nmics,is_ptcl=.false.)
                cnt = 0
                do imic = 1,tmp_proj%os_mic%get_noris()
                    if( states(imic) == 0 ) cycle
                    cnt = cnt+1
                    call spproj_here%os_mic%transfer_ori(cnt, tmp_proj%os_mic, imic)
                enddo
                nselected = cnt
                nrejected = tmp_proj%os_mic%get_noris() - nselected
                ! update for execution
                extract_set_counter = extract_set_counter + 1
                projname = int2str_pad(extract_set_counter,params%numlen)
                projfile = trim(projname)//trim(METADATA_EXT)
                if( l_pick_extract )then
                    call cline_pick_extract%set('projname', trim(projname))
                    call cline_pick_extract%set('projfile', trim(projfile))
                    call cline_pick_extract%set('fromp',    1)
                    call cline_pick_extract%set('top',      nselected)
                else
                    call cline_extract%set('projname', trim(projname))
                    call cline_extract%set('projfile', trim(projfile))
                    call cline_extract%set('fromp',    1)
                    call cline_extract%set('top',      nselected)
                endif
                call spproj_here%write(trim(path)//'/'//trim(projfile))
                call spproj_here%kill
                call tmp_proj%kill
            end subroutine create_individual_project

    end subroutine exec_stream_gen_picking_refs

    subroutine exec_stream_assign_optics( self, cline )
        class(commander_stream_assign_optics), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj, spproj_part
        type(starproject_stream)               :: starproj_stream
        character(len=LONGSTRLEN), allocatable :: projects(:)
        integer                                :: nprojects, iproj, iori, new_oris, nimported
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('dir_target') ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('outdir')     ) call cline%set('outdir', '')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exist: '//trim(cline%get_carg('dir_exec')))
            endif
            call del_file(TERM_STREAM)
        endif
        ! master parameters
        call params%new(cline)
        ! master project file
        call spproj%read( params%projfile )
        call spproj%update_projinfo(cline)
        if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(0, .false.)
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! initialise progress monitor
        call progressfile_init()
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('micrographs', 'micrographs_imported',  int2str(0), primary=.true.)
        call gui_stats%set('groups',      'optics_group_assigned', int2str(0), primary=.true.)
        ! Infinite loop
        nimported = 0
        do
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING STREAM ASSIGN OPTICS'
                exit
            endif
            ! detection of new projects
            call project_buff%watch( nprojects, projects, max_nmovies=50 )
            ! append projects to processing stack
            if( nprojects > 0 )then
                nimported = spproj%os_mic%get_noris()
                if(nimported > 0) then
                    new_oris  =  nimported + nprojects * NMOVS_SET
                    call spproj%os_mic%reallocate(new_oris)
                else
                    new_oris = nprojects * NMOVS_SET
                    call spproj%os_mic%new(new_oris, .false.)
                end if
                do iproj = 1, nprojects
                    call project_buff%add2history(projects(iproj))
                    call spproj_part%read(trim(projects(iproj)))
                    do iori = 1, NMOVS_SET
                        nimported = nimported + 1
                        call spproj%os_mic%transfer_ori(nimported, spproj_part%os_mic, iori)
                    end do
                    call spproj_part%kill()
                enddo
                write(logfhandle,'(A,I4,A,A)')'>>> ' , nprojects * NMOVS_SET, ' NEW MICROGRAPHS IMPORTED; ',cast_time_char(simple_gettime())
                call starproj_stream%stream_export_optics(spproj, params%outdir)
                call starproj_stream%stream_export_micrographs(spproj, params%outdir, optics_set=.true.)
                ! guistats
                call gui_stats%set('micrographs', 'micrographs_imported', int2str(0), primary=.true.)
                call gui_stats%set('groups', 'optics_group_assigned', spproj%os_optics%get_noris(), primary=.true.)
            else
                call sleep(WAITTIME) ! may want to increase as 3s default
            endif
            call update_user_params(cline)
            if(params%updated .eq. 'yes') then
                call starproj_stream%stream_export_optics(spproj, params%outdir)
                ! guistats
                call gui_stats%set('groups', 'optics_group_assigned', spproj%os_optics%get_noris(), primary=.true.)
                params%updated = 'no'
            end if
            call gui_stats%write_json
        end do
        if(allocated(projects)) deallocate(projects)
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj%write
        call spproj%kill
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS NORMAL STOP ****')
    end subroutine exec_stream_assign_optics

    subroutine exec_stream_cluster2D( self, cline )
        class(commander_stream_cluster2D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        type(oris)                             :: moldiamori
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .true.
        type(micproj_record),      allocatable :: micproj_records(:)
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=LONGSTRLEN)              :: cwd_job
        integer                                :: nmics_rejected_glob
        integer                                :: nchunks_imported_glob, nchunks_imported, nprojects, iter
        integer                                :: n_imported, n_added, nptcls_glob, n_failed_jobs, ncls_in, nmic_star
        integer                                :: pool_iter_last_chunk_imported, pool_iter_max_chunk_imported
        logical                                :: l_nchunks_maxed, l_pause, l_params_updated
        real                                   :: nptcls_pool, moldiam
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('kweight_chunk','default')
        call cline%set('kweight_pool', 'default')
        call cline%set('prune',        'no')
        call cline%set('wiener',       'full')
        if( .not. cline%defined('dir_target')   ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('walltime')     ) call cline%set('walltime',     29*60) ! 29 minutes
        if( .not. cline%defined('lpthres')      ) call cline%set('lpthres',      RES_THRESHOLD_STREAM)
        if( .not. cline%defined('ndev')         ) call cline%set('ndev',         CLS_REJECT_STD)
        if( .not. cline%defined('reject_cls')   ) call cline%set('reject_cls',   'yes')
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',       'euclid')
        if( .not. cline%defined('ml_reg')       ) call cline%set('ml_reg',       'no')
        if( .not. cline%defined('tau')          ) call cline%set('tau',          5)
        if( .not. cline%defined('rnd_cls_init') ) call cline%set('rnd_cls_init', 'no')
        if( .not. cline%defined('remove_chunks')) call cline%set('remove_chunks','yes')
        if( .not. cline%defined('refine')       ) call cline%set('refine',       'snhc_smpl')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ncls_in = 0
        if( cline%defined('ncls') )then
            ! to circumvent parameters class stringency, restored after params%new
            ncls_in = nint(cline%get_rarg('ncls'))
            call cline%delete('ncls')
        endif        
        ! master parameters
        call cline%set('numlen', 5)
        call cline%set('stream','yes')
        call params%new(cline)
        params%nthr2D       = params%nthr
        params%ml_reg_chunk = trim(params%ml_reg)
        params%ml_reg_pool  = trim(params%ml_reg)
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        if( ncls_in > 0 ) call cline%set('ncls', real(ncls_in))
        ! limit to # of chunks
        if( .not. cline%defined('maxnchunks') .or. params%maxnchunks < 1 )then
            params%maxnchunks = huge(params%maxnchunks)
        endif
        call cline%delete('maxnchunks')
        ! restart
        if( cline%defined('dir_exec') )then
            call cline%delete('dir_exec')
            call del_file(micspproj_fname)
            call cleanup_root_folder
        endif
         ! mskdiam
        if( .not. cline%defined('mskdiam') ) then
            if( .not. file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) THROW_HARD('either mskdiam must be given or '// trim(STREAM_MOLDIAM) // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            if( .not. moldiamori%isthere(1, "moldiam") ) THROW_HARD( 'moldiam missing from ' // trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            moldiam = moldiamori%get(1, "moldiam")
            call moldiamori%kill
            call cline%set('mskdiam', moldiam * 1.2)
            params%mskdiam = moldiam * 1.2
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
        endif
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        call spproj_glob%update_projinfo(cline)
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true.)
        ! Infinite loop
        nptcls_glob           = 0   ! global number of particles
        nchunks_imported_glob = 0   ! global number of completed chunks
        nprojects             = 0
        iter                  = 0
        n_imported            = 0   ! global number of imported processed micrographs
        n_failed_jobs         = 0
        nmic_star             = 0
        nmics_rejected_glob   = 0   ! global number of micrographs rejected
        l_nchunks_maxed       = .false.
        l_pause               = .false.
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('particles', 'particles_imported',          0,            primary=.true.)
        call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        call gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
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
                ! watch & update global records
                call project_buff%watch(nprojects, projects, max_nmovies=10*params%nparts)
            endif
            ! update global records
            if( nprojects > 0 )then
                call update_records_with_project(projects, n_imported )
                call project_buff%add2history(projects)
            endif
            ! project update
            if( nprojects > 0 )then
                n_imported = size(micproj_records)
                write(logfhandle,'(A,I6,I8)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ',n_imported, nptcls_glob
                ! guistats
                call gui_stats%set('particles', 'particles_imported', int2commastr(nptcls_glob), primary=.true.)
                call gui_stats%set_now('particles', 'last_particles_imported')
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! remove this?
                if( n_imported < 1000 )then
                    call update_user_params(cline)
                else if( n_imported > nmic_star + 100 )then
                    call update_user_params(cline)
                    nmic_star = n_imported
                endif
            endif
            ! 2D classification section
            call update_user_params_dev(cline, l_params_updated)
            if( l_params_updated ) l_pause = .false.
            call update_chunks_dev
            if( l_pause )then
                call update_user_params_dev(cline, l_params_updated)
                if( l_params_updated ) l_pause = .false.
            else
                call update_pool_status_dev
                call update_pool_dev
                call update_user_params_dev(cline, l_params_updated)
                call reject_from_pool_dev
            endif
            call reject_from_pool_user_dev
            if( l_nchunks_maxed )then
                ! # of chunks is above desired number
                if( is_pool_available_dev() .and. (get_pool_iter()>=pool_iter_max_chunk_imported+10) ) exit
                call classify_pool_dev
            else
                call import_chunks_into_pool_dev( nchunks_imported )
                if( nchunks_imported > 0 )then
                    nchunks_imported_glob         = nchunks_imported_glob + nchunks_imported
                    pool_iter_last_chunk_imported = get_pool_iter()
                    l_pause = .false.
                endif
                if( nchunks_imported_glob >= params%maxnchunks )then
                    if( .not.l_nchunks_maxed ) pool_iter_max_chunk_imported = get_pool_iter()
                    l_nchunks_maxed = .true.
                endif
                if( l_pause )then
                    ! skipping pool classification
                else
                    if( get_pool_iter() >= pool_iter_last_chunk_imported+10 )then
                        ! pause pool classification & rejection in absence of new chunks, resumes
                        ! when new chunks are added or classification parameters have been updated
                        l_pause = is_pool_available_dev()
                    endif
                endif
                if( .not.l_pause ) call classify_pool_dev
                call classify_new_chunks_dev(micproj_records)
            endif
            call sleep(WAITTIME)
            ! guistats
            if(file_exists(POOLSTATS_FILE)) then
                call gui_stats%merge(POOLSTATS_FILE)
                call gui_stats%get('particles', 'particles_processed', nptcls_pool)
                if(nptcls_pool > 0.0) call gui_stats%set('particles', 'particles_processed', int2commastr(floor(nptcls_pool)) // ' (' // int2str(ceiling(100.0 * real(nptcls_pool) / real(nptcls_glob))) // '%)')
            end if
            call gui_stats%write_json
        end do
        ! termination
        call terminate_stream2D_dev
        call update_user_params(cline) !!??
        ! final stats
        if(file_exists(POOLSTATS_FILE)) call gui_stats%merge(POOLSTATS_FILE, delete = .true.)
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_CLUSTER2D NORMAL STOP ****')
        contains

            ! updates global records
            subroutine update_records_with_project( projectnames, n_imported )
                character(len=LONGSTRLEN), allocatable, intent(in)  :: projectnames(:)
                integer,                                intent(out) :: n_imported
                type(sp_project),     allocatable :: spprojs(:)
                type(micproj_record), allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname
                integer :: iproj, n_spprojs, n_old, irec, n_completed, nptcls, nmics, imic, n_ptcls, first
                n_imported = 0
                n_ptcls    = 0
                if( .not.allocated(projectnames) ) return
                n_spprojs  = size(projectnames)
                if( n_spprojs == 0 )return
                n_old = 0 ! on first import
                if( allocated(micproj_records) ) n_old = size(micproj_records)
                allocate(spprojs(n_spprojs))
                ! because pick_extract purges state=0 and nptcls=0 mics,
                ! all mics can be assumed associated with particles
                nmics = 0
                first = 0
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%read_segment('mic', trim(projectnames(iproj)))
                    nmics = nmics + spprojs(iproj)%os_mic%get_noris()
                    if( (first == 0) .and. (nmics > 0) ) first = iproj
                enddo
                if( nmics == 0 ) return
                ! Updates global parameters once and init 2D
                if( n_old == 0 )then
                    params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                    call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                    params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                    call init_cluster2D_stream_dev(cline, spproj_glob, params%box, micspproj_fname)
                    call cline%delete('ncls')
                endif
                ! import micrographs
                n_completed = n_old + nmics
                n_imported  = nmics
                ! reallocate records
                if( n_old == 0 )then
                    allocate(micproj_records(nmics))
                else
                    call move_alloc(micproj_records, old_records)
                    allocate(micproj_records(n_completed))
                    micproj_records(1:n_old) = old_records(:)
                    deallocate(old_records)
                endif
                ! update global records and some global variables
                irec = n_old
                do iproj = 1,n_spprojs
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        irec      = irec + 1
                        nptcls    = nint(spprojs(iproj)%os_mic%get(imic,'nptcls'))
                        n_ptcls   = n_ptcls + nptcls ! global update
                        fname     = trim(projectnames(iproj))
                        abs_fname = simple_abspath(fname, errmsg='stream_cluster2D :: update_projects_list 1')
                        micproj_records(irec)%projname   = trim(abs_fname)
                        micproj_records(irec)%micind     = imic
                        micproj_records(irec)%nptcls     = nptcls
                        micproj_records(irec)%nptcls_sel = nptcls
                        micproj_records(irec)%included   = .false.
                    enddo
                enddo
                nptcls_glob = nptcls_glob + n_ptcls ! global update
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_records_with_project

    end subroutine exec_stream_cluster2D

    ! PRIVATE UTILITIES

    !> updates current parameters with user input
    subroutine update_user_params( cline_here )
        type(cmdline), intent(inout) :: cline_here
        type(oris) :: os
        real       :: tilt_thres, beamtilt, astigthreshold, ctfresthreshold, icefracthreshold
        real       :: moldiam_refine
        call os%new(1, is_ptcl=.false.)
        if( file_exists(USER_PARAMS) )then
            call os%read(USER_PARAMS)
            if( os%isthere(1,'tilt_thres') ) then
                tilt_thres = os%get(1,'tilt_thres')
                if( abs(tilt_thres-params_glob%tilt_thres) > 0.001) then
                     if(tilt_thres < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO LOW: ',tilt_thres
                     else if(tilt_thres > 1) then
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES TOO HIGH: ',tilt_thres
                     else
                         params_glob%tilt_thres = tilt_thres
                         params_glob%updated    = 'yes'
                         call cline_here%set('tilt_thres', params_glob%tilt_thres)
                         write(logfhandle,'(A,F8.2)')'>>> OPTICS TILT_THRES UPDATED TO: ',tilt_thres
                     endif
                endif
            endif
            if( os%isthere(1,'beamtilt') ) then
                beamtilt = os%get(1,'beamtilt')
                if( beamtilt .eq. 1.0 ) then
                    params_glob%beamtilt = 'yes'
                    params_glob%updated  = 'yes'
                    call cline_here%set('beamtilt', params_glob%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO USE BEAMTILT'
                else if( beamtilt .eq. 0.0 ) then
                    params_glob%beamtilt = 'no'
                    params_glob%updated  = 'yes'
                    call cline_here%set('beamtilt', params_glob%beamtilt)
                    write(logfhandle,'(A)')'>>> OPTICS ASSIGNMENT UDPATED TO IGNORE BEAMTILT'   
                else
                    write(logfhandle,'(A,F8.2)')'>>> OPTICS UPDATE INVALID BEAMTILT VALUE: ',beamtilt
                endif
            endif
            if( os%isthere(1,'astigthreshold') ) then
                astigthreshold = os%get(1,'astigthreshold')
                if( abs(astigthreshold-params_glob%astigthreshold) > 0.001) then
                     if(astigthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD TOO LOW: ',astigthreshold
                     else if(astigthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM TOO HIGH: ',astigthreshold
                     else
                         params_glob%astigthreshold = astigthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('astigthreshold', params_glob%astigthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ',astigthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'ctfresthreshold') ) then
                ctfresthreshold = os%get(1,'ctfresthreshold')
                if( abs(ctfresthreshold-params_glob%ctfresthreshold) > 0.001) then
                     if(ctfresthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD TOO LOW: ',ctfresthreshold
                     else if(ctfresthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION TOO HIGH: ',ctfresthreshold
                     else
                         params_glob%ctfresthreshold = ctfresthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('ctfresthreshold', params_glob%ctfresthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ',ctfresthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'icefracthreshold') ) then
                icefracthreshold = os%get(1,'icefracthreshold')
                if( abs(icefracthreshold-params_glob%icefracthreshold) > 0.001) then
                     if(icefracthreshold < 0.01)then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD TOO LOW: ',icefracthreshold
                     else if(icefracthreshold > 100.0) then
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION TOO HIGH: ',icefracthreshold
                     else
                         params_glob%icefracthreshold = icefracthreshold
                         params_glob%updated    = 'yes'
                         call cline_here%set('icefracthreshold', params_glob%icefracthreshold)
                         write(logfhandle,'(A,F8.2)')'>>> ICE FRACTION THRESHOLD UPDATED TO: ',icefracthreshold
                     endif
                endif
            endif
            if( os%isthere(1,'moldiam_refine') ) then
                moldiam_refine = os%get(1,'moldiam_refine')
                if( abs(moldiam_refine-params_glob%moldiam_refine) > 0.001) then
                     if(moldiam_refine < 10)then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO LOW: ' , moldiam_refine
                     else if(moldiam_refine > 1000) then
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER TOO HIGH: ', moldiam_refine
                     else
                         params_glob%moldiam_refine = moldiam_refine
                         params_glob%updated        = 'yes'
                         write(logfhandle,'(A,F8.2)')'>>> MULTIPICK REFINE DIAMETER UPDATED TO: ', moldiam_refine
                     endif
                endif
            endif
            call del_file(USER_PARAMS)
        endif
        call os%kill
    end subroutine update_user_params

end module simple_commander_stream
