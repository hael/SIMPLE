! concrete commander: streaming pre-processing routines
module simple_commander_stream
include 'simple_lib.f08'
use simple_binoris_io
use simple_cmdline,            only: cmdline
use simple_parameters,         only: parameters
use simple_commander_base,     only: commander_base
use simple_sp_project,         only: sp_project
use simple_qsys_env,           only: qsys_env
use simple_starproject_stream, only: starproject_stream
use simple_guistats,           only: guistats
use simple_moviewatcher,       only: moviewatcher
use simple_stream_utils
use simple_stream_communicator
use simple_commander_cluster2D_stream
use simple_qsys_funs
use simple_commander_preprocess
use simple_progress
use simple_nrtxtfile
use simple_imgproc
use simple_timer
use simple_nice
implicit none

real, parameter, dimension(21)  :: astig_bins     = [0.0, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(19)  :: ctfres_bins    = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 11.0, 12.0, 13.0, 14.0, 15.0, 16.0, 17.0, 18.0, 19.0, 20.0]
real, parameter, dimension(21)  :: icescore_bins  = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5, 1.6, 1.7, 1.8, 1.9, 2.0]

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
integer,               parameter :: NMOVS_SET      = 5                                           ! number of movies processed at once (>1)
integer,               parameter :: NMOVS_SET_TIFF = 3                                           ! number of TIFF movies processed at once (>1)
integer,               parameter :: LONGTIME       = 60                                          ! time lag after which a movie/project is processed
integer,               parameter :: WAITTIME       = 10                                          ! movie folder watched every WAITTIME seconds
integer,               parameter :: SHORTWAIT      = 2                                           ! movie folder watched every SHORTTIME seconds in shmem

contains

    subroutine exec_stream_preprocess( self, cline )
        use simple_motion_correct_utils, only: flip_gain
        class(commander_stream_preprocess), intent(inout) :: self
        class(cmdline),                     intent(inout) :: cline
        type(parameters)                                  :: params
        type(stream_http_communicator)                    :: http_communicator
        type(json_value),                   pointer       :: latest_micrographs, histograms
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
        character(len=STDLEN)                  :: preproc_nthr_env, preproc_part_env, preproc_nparts_env
        integer                                :: movies_set_counter, import_counter
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection, nsets, i, j, i_thumb, i_max
        integer                                :: cnt, n_imported, n_added, n_failed_jobs, n_fail_iter, nmic_star, iset, envlen
        logical                                :: l_movies_left, l_haschanged
        logical(LK)                            :: found 
        real                                   :: avg_tmp, stat_dfx_threshold, stat_dfy_threshold
        real                                   :: stat_astig_threshold, stat_icefrac_threshold, stat_ctfres_threshold
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('reject_mics', 'no')
        if( .not. cline%defined('walltime')         ) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! motion correction
        if( .not. cline%defined('trs')              ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart',           8.)
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',            5.)
        if( .not. cline%defined('bfac')             ) call cline%set('bfac',             50.)
        if( .not. cline%defined('mcconvention')     ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')   ) call cline%set('eer_upsampling',     1)
        if( .not. cline%defined('algorithm')        ) call cline%set('algorithm',    'patch')
        if( .not. cline%defined('mcpatch')          ) call cline%set('mcpatch',        'yes')
        if( .not. cline%defined('mcpatch_thres')    ) call cline%set('mcpatch_thres','  yes')
        if( .not. cline%defined('tilt_thres')       ) call cline%set('tilt_thres',      0.05)
        if( .not. cline%defined('beamtilt')         ) call cline%set('beamtilt',        'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp_ctf_estimate')  ) call cline%set('hp_ctf_estimate', HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate')  ) call cline%set('lp_ctf_estimate', LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',           DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',           DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',        'yes')
        if( .not. cline%defined('ctfresthreshold') )  call cline%set('ctfresthreshold',  CTFRES_THRESHOLD_STREAM)
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NTHR, preproc_nthr_env, envlen)
        if( envlen > 0)  call cline%set('nthr', str2int(preproc_nthr_env))
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NPARTS, preproc_nparts_env, envlen)
        if( envlen > 0 ) call cline%set('nparts', str2int(preproc_nparts_env))
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'preprocess')
            call cline%set('projfile', 'preprocess.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params%split_mode = 'stream'
        params%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "preprocessing")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME, params%dir_movies)
        ! restart
        movies_set_counter = 0  ! global number of movies set
        import_counter     = 0  ! global import id
        nmic_star          = 0
        if( cline%defined('dir_exec') )then
            call del_file(TERM_STREAM)
            call cline%delete('dir_exec')
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage", "importing previously processed data", found)
            call http_communicator%send_jobstats()
            call import_previous_projects
            nmic_star = spproj_glob%os_mic%get_noris()
            call write_mic_star_and_field(write_field=.true.)
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "movies_imported",  spproj_glob%os_mic%get_noris(), found)
            call http_communicator%json%update(http_communicator%job_json, "movies_processed", spproj_glob%os_mic%get_noris(), found)
            call http_communicator%json%update(http_communicator%job_json, "movies_rejected",  0,                              found)
            call communicator_clear_histograms()
            if(spproj_glob%os_mic%isthere("ctfres")) then
                avg_tmp = spproj_glob%os_mic%get_avg("ctfres")
                call http_communicator%json%update(http_communicator%job_json, "average_ctf_res", dble(avg_tmp), found)
                call communicator_add_histogram("ctfres")
                !nice_communicator%view_micrographs%cutoff_ctf_res = params%ctfresthreshold
            end if
            if(spproj_glob%os_mic%isthere("icefrac")) then
                avg_tmp = spproj_glob%os_mic%get_avg("icefrac")
                call http_communicator%json%update(http_communicator%job_json, "average_ice_score", dble(avg_tmp), found)
                call communicator_add_histogram("icefrac")
               ! nice_communicator%view_micrographs%cutoff_ice_score = params%icefracthreshold
            end if
            if(spproj_glob%os_mic%isthere("astig")) then
                avg_tmp = spproj_glob%os_mic%get_avg("astig")
                call http_communicator%json%update(http_communicator%job_json, "average_astigmatism", dble(avg_tmp), found)
                call communicator_add_histogram("astig")
               ! nice_communicator%view_micrographs%cutoff_astigmatism = params%astigthreshold
            end if
            if(spproj_glob%os_mic%isthere('thumb')) then
              !  call gui_stats%set('latest', '', trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(),'thumb'))), thumbnail=.true.)
            end if
            call http_communicator%send_jobstats()
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
        stat_dfx_threshold     = 0.
        stat_dfy_threshold     = 0.
        stat_astig_threshold   = 0.
        stat_icefrac_threshold = 0.
        stat_ctfres_threshold  = 0.
        call cline_exec%set('fromp',1)
        call cline_exec%set('top',  NMOVS_SET)
        do
            if( file_exists(trim(TERM_STREAM)) .or. http_communicator%exit )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
             !   call http_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_PREPROC USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage",               "finding and processing new movies", found)    
            call http_communicator%json%update(http_communicator%job_json, "cutoff_ctf_res",      dble(params%ctfresthreshold),        found)
            call http_communicator%json%update(http_communicator%job_json, "cutoff_ice_score",    dble(params%icefracthreshold),       found)
            call http_communicator%json%update(http_communicator%job_json, "cutoff_astigmatism",  dble(params%astigthreshold),         found)
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
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "movies_imported",      movie_buff%n_history, found)
                call http_communicator%json%update(http_communicator%job_json, "last_movie_imported",  stream_datestr(),     found)
            else
                l_movies_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
            stacksz = qenv%qscripts%get_stacksz()
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
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "movies_processed", n_imported, found)
                if( n_failed_jobs > 0 )  call http_communicator%json%update(http_communicator%job_json, "movies_rejected",  n_failed_jobs, found)
                call communicator_clear_histograms()
                if(spproj_glob%os_mic%isthere("ctfres")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("ctfres")
                    call http_communicator%json%update(http_communicator%job_json, "average_ctf_res", dble(avg_tmp), found)
                    call communicator_add_histogram("ctfres")
                end if
                if(spproj_glob%os_mic%isthere("icefrac")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("icefrac")
                    call http_communicator%json%update(http_communicator%job_json, "average_ice_score", dble(avg_tmp), found)
                    call communicator_add_histogram("icefrac")
                end if
                if(spproj_glob%os_mic%isthere("astig")) then
                    avg_tmp = spproj_glob%os_mic%get_avg("astig")
                    call http_communicator%json%update(http_communicator%job_json, "average_astigmatism", dble(avg_tmp), found)
                    call communicator_add_histogram("astig")
                end if
                if(spproj_glob%os_mic%isthere('thumb')) then
                    i_max = 10
                    if(spproj_glob%os_mic%get_noris() < 10) i_max = spproj_glob%os_mic%get_noris()
                    ! create an empty latest_micrographs json array
                    call http_communicator%json%remove(latest_micrographs, destroy=.true.)
                    call http_communicator%json%create_array(latest_micrographs, "latest_micrographs")
                    call http_communicator%json%add(http_communicator%job_json, latest_micrographs)
                    do i_thumb=0, i_max - 1
                        call communicator_add_micrograph(trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris() - i_thumb, "thumb"))))
                    end do
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                last_injection = simple_gettime()
                l_haschanged = .true.
                n_imported   = spproj_glob%os_mic%get_noris()
                ! sliding window test
                if(n_imported > 100) then
                    if(stat_dfx_threshold == 0) then
                        call set_stat_thresholds()
                    else
                        call test_stat_thresholds()
                    endif
                endif
                ! always write micrographs snapshot if less than 1000 mics, else every 100
                if( n_imported < 1000 .and. l_haschanged )then
  !!                 ! call update_user_params(cline, nice_communicator%update_arguments)
                    call write_mic_star_and_field
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
  !!                !  call update_user_params(cline, nice_communicator%update_arguments)
                    call write_mic_star_and_field
                    nmic_star = n_imported
                endif
            else
                ! wait & write snapshot
                if( .not.l_movies_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive...
   !!                 !    call update_user_params(cline, nice_communicator%update_arguments)
                        call write_mic_star_and_field
                        l_haschanged = .false.
                    else
                        ! ...or wait
                        call sleep(WAITTIME)
                    endif
                endif
            endif
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        call http_communicator%json%update(http_communicator%job_json, "stage", "terminating", found) 
        call http_communicator%send_jobstats()
        call update_user_params(cline)
        call write_mic_star_and_field(write_field=.true., copy_optics=.true.)
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_PREPROC NORMAL STOP ****')
        contains

            subroutine mics_window_stats(key, fromto, avg, sdev)
                character(len=*), intent(in)    :: key
                integer,          intent(in)    :: fromto(2)
                real,             intent(inout) :: avg, sdev
                real,             allocatable   :: arr(:)
                arr = spproj_glob%os_mic%get_all(key, fromto)
                call avg_sdev(arr, avg, sdev)
                if(allocated(arr)) deallocate(arr)
            end subroutine mics_window_stats

            subroutine set_stat_thresholds()
                real    :: avg_window, sdev_window
                integer :: window(2)
                window(1) = 1
                window(2) = 100
                call mics_window_stats("astig", window, avg_window, sdev_window)
                stat_astig_threshold = avg_window + sdev_window
                call mics_window_stats("ctfres", window, avg_window, sdev_window)
                stat_ctfres_threshold = avg_window + sdev_window
                call mics_window_stats("icefrac", window, avg_window, sdev_window)
                stat_icefrac_threshold = avg_window + sdev_window
                call mics_window_stats("dfx", window, avg_window, sdev_window)
                stat_dfx_threshold = avg_window + sdev_window
                call mics_window_stats("dfy", window, avg_window, sdev_window)
                stat_dfy_threshold = avg_window + sdev_window
            end subroutine set_stat_thresholds

            subroutine test_stat_thresholds()
                real    :: avg_window, sdev_window, avg_glob
                integer :: window(2)
                window(1) = spproj_glob%os_mic%get_noris() - 100
                window(2) = spproj_glob%os_mic%get_noris()
                call mics_window_stats("astig", window, avg_window, sdev_window)
                if(avg_window > stat_astig_threshold) then
                    avg_glob = spproj_glob%os_mic%get_avg("astig")
                 !   call gui_stats%set('micrographs', 'avg_astigmatism', avg_glob, primary=.true., alert=.true., alerttext='average astigmatism &
                  !          &has increased by more that 1 sigma within the past 100 micrographs. Please check collection', notify=.false.)
                endif
                call mics_window_stats("ctfres", window, avg_window, sdev_window)
                if(avg_window > stat_ctfres_threshold) then
                    avg_glob = spproj_glob%os_mic%get_avg("ctfres")
                  !  call gui_stats%set('micrographs', 'avg_ctf_resolution', avg_glob, primary=.true., alert=.true., alerttext='average ctf resolution &
                   !         &has decreased by more that 1 sigma within the past 100 micrographs. Please check collection', notify=.false.)
                endif
                call mics_window_stats("icefrac", window, avg_window, sdev_window)
                if(avg_window > stat_icefrac_threshold) then
                    avg_glob = spproj_glob%os_mic%get_avg("icefrac")
                  !  call gui_stats%set('micrographs', 'avg_ice_score', avg_glob, primary=.true., alert=.true., alerttext='average ice score &
                   !         &has increased by more that 1 sigma within the past 100 micrographs. Please check collection', notify=.false.)
                endif
                call mics_window_stats("dfx", window, avg_window, sdev_window)
                if(avg_window > stat_dfx_threshold) then
                    avg_glob = spproj_glob%os_mic%get_avg("dfx")
                  !  call gui_stats%set('micrographs', 'avg_defocus_x', avg_glob, primary=.true., alert=.true., alerttext='average defocus in x &
                   !         &has increased by more that 1 sigma within the past 100 micrographs. Please check collection', notify=.false.)
                else
                !    call gui_stats%delete('micrographs', 'avg_defocus_x')
                endif
                call mics_window_stats("dfy", window, avg_window, sdev_window)
                if(avg_window > stat_dfy_threshold) then
                    avg_glob = spproj_glob%os_mic%get_avg("dfy")
                  !  call gui_stats%set('micrographs', 'avg_defocus_y', avg_glob, primary=.true., alert=.true., alerttext='average defocus in y &
                  !          &has increased by more that 1 sigma within the past 100 micrographs. Please check collection', notify=.false.)
                else
                 !   call gui_stats%delete('micrographs', 'avg_defocus_y')
                endif
            end subroutine test_stat_thresholds

            subroutine write_mic_star_and_field( write_field, copy_optics )
                logical, optional, intent(in) :: write_field, copy_optics
                logical :: l_wfield, l_copy_optics
                l_wfield      = .false.
                l_copy_optics = .false.
                if( present(write_field) ) l_wfield      = write_field
                if( present(copy_optics) ) l_copy_optics = copy_optics
                if(l_copy_optics) then
                    call starproj_stream%copy_micrographs_optics(spproj_glob, verbose=DEBUG_HERE)
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=.true.)
                else
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir)
                end if
                if( l_wfield )then
                    call spproj_glob%write_segment_inside('mic', params%projfile)
                    call spproj_glob%write_non_data_segments(params%projfile)
                endif
            end subroutine write_mic_star_and_field

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
                                j = j + 1
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
                    id    = str2int(fname)
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

            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",               "initialising")
                call http_communicator%json%add(http_communicator%job_json, "movies_imported",     0)
                call http_communicator%json%add(http_communicator%job_json, "movies_processed",    0)
                call http_communicator%json%add(http_communicator%job_json, "movies_rejected",     0)
                call http_communicator%json%add(http_communicator%job_json, "average_ctf_res",     dble(0.0))
                call http_communicator%json%add(http_communicator%job_json, "average_ice_score",   dble(0.0))
                call http_communicator%json%add(http_communicator%job_json, "average_astigmatism", dble(0.0))
                call http_communicator%json%add(http_communicator%job_json, "cutoff_ctf_res",      dble(params%ctfresthreshold))
                call http_communicator%json%add(http_communicator%job_json, "cutoff_ice_score",    dble(params%icefracthreshold))
                call http_communicator%json%add(http_communicator%job_json, "cutoff_astigmatism",  dble(params%astigthreshold))
                call http_communicator%json%add(http_communicator%job_json, "last_movie_imported", "")
                call http_communicator%json%create_array(latest_micrographs, "latest_micrographs")
                call http_communicator%json%add(http_communicator%job_json, latest_micrographs)
                call http_communicator%json%create_object(histograms,        "histograms")
                call http_communicator%json%add(http_communicator%job_json, histograms)
            end subroutine communicator_init

            subroutine communicator_add_micrograph(path)
                character(*),     intent(in) :: path
                type(json_value), pointer    :: micrograph
                call http_communicator%json%create_object(micrograph, "")
                call http_communicator%json%add(micrograph, "path", path)
                call http_communicator%json%add(latest_micrographs, micrograph)
            end subroutine communicator_add_micrograph

            subroutine communicator_clear_histograms()
                type(json_value), pointer    :: micrograph
                call http_communicator%json%remove(histograms, destroy=.true.)
                call http_communicator%json%create_object(histograms, "histograms")
                call http_communicator%json%add(http_communicator%job_json, histograms)
            end subroutine communicator_clear_histograms

            subroutine communicator_add_histogram(key)
                character(*),     intent(in) :: key
                type(json_value), pointer    :: new_histogram, histogram_labels, histogram_data
                type(histogram)              :: key_histogram
                real,          allocatable   :: histogram_rvec(:)
                integer                      :: ikey, ibin
                if(key == "ctfres") then
                    allocate(histogram_rvec(size(ctfres_bins)))
                    histogram_rvec = ctfres_bins
                else if(key == "astig") then
                    allocate(histogram_rvec(size(astig_bins)))
                    histogram_rvec = astig_bins
                else if(key == "icefrac") then
                    allocate(histogram_rvec(size(icescore_bins)))
                    histogram_rvec = icescore_bins
                else
                    return
                endif
                call key_histogram%new(histogram_rvec)
                call key_histogram%zero()
                do ikey = 1, spproj_glob%os_mic%get_noris()
                    call key_histogram%update(spproj_glob%os_mic%get(ikey, key))
                enddo
                call http_communicator%json%create_object(new_histogram, key)
                call http_communicator%json%create_array(histogram_labels, "labels")
                call http_communicator%json%add(new_histogram, histogram_labels)
                call http_communicator%json%create_array(histogram_data,   "data")
                call http_communicator%json%add(new_histogram, histogram_data)
                do ibin=1, size(histogram_rvec)
                    call http_communicator%json%add(histogram_labels, "", dble(histogram_rvec(ibin)))
                    call http_communicator%json%add(histogram_data,   "", int(key_histogram%get(ibin)))
                end do
                call http_communicator%json%add(histograms, new_histogram)
                call key_histogram%kill()
                deallocate(histogram_rvec)
            end subroutine communicator_add_histogram


    end subroutine exec_stream_preprocess

    subroutine exec_stream_pick_extract( self, cline )
        use simple_histogram,    only: histogram
        use simple_image
        class(commander_stream_pick_extract), intent(inout) :: self
        class(cmdline),                       intent(inout) :: cline
        type(make_pickrefs_commander)          :: xmake_pickrefs
        type(parameters)                       :: params
        type(simple_nice_communicator)         :: nice_communicator
        type(guistats)                         :: gui_stats
        integer,                   parameter   :: INACTIVE_TIME = 900  ! inactive time triggers writing of project file
        logical,                   parameter   :: DEBUG_HERE    = .false.
        logical,                   allocatable :: jobs_done(:), jobs_submitted(:)
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(projrecord),          allocatable :: projrecords(:), projrecords_main(:)
        type(qsys_env),            pointer     :: qenv
        type(qsys_env),            target      :: qenv_main, qenv_interactive
        type(cmdline)                          :: cline_make_pickrefs, cline_pick_extract
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob, stream_spproj, spproj_tmp, interactive_spproj
        type(starproject_stream)               :: starproj_stream
        type(stream_http_communicator)         :: http_communicator
        type(histogram)                        :: histogram_moldiams
        type(image)                            :: mic_thumb
        type(json_value),          pointer     :: latest_picked_micrographs, latest_extracted_particles, picking_templates, picking_diameters
        type(nrtxtfile)                        :: boxsize_file
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=:),          allocatable :: odir, odir_extract, odir_picker, odir_completed
        character(len=:),          allocatable :: odir_interactive, odir_interactive_picker, odir_interactive_completed
        character(len=LONGSTRLEN)              :: cwd_job, latest_boxfile
        character(len=STDLEN)                  :: pick_nthr_env, pick_part_env
        real,                      allocatable :: moldiams(:), saved_boxsize(:)
        real                                   :: jpg_scale, pickrefs_thumbnail_scale, rnd
        integer,                   allocatable :: complete_search_diameters(:), active_search_diameters(:), refined_search_diameters(:)
        integer                                :: nmics_sel, nmics_rej, nmics_rejected_glob, pick_extract_set_counter, i_max, i_thumb, i
        integer                                :: nmics, nprojects, stacksz, prev_stacksz, iter, last_injection, iproj, envlen, imic
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, nmic_star, thumbid_offset
        integer                                :: n_pickrefs, thumbcount, xtile, ytile, xtiles, ytiles
        integer                                :: ldim_mic(3)
        logical                                :: l_templates_provided, l_projects_left, l_haschanged, l_multipick, l_extract, l_once
        logical                                :: pause_import, l_interactive, interactive_waiting, found
        integer(timer_int_kind) :: t0
        real(timer_int_kind)    :: rt_write
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        call cline%set('picker',  'new')
        if( .not. cline%defined('outdir')          ) call cline%set('outdir',           '')
        if( .not. cline%defined('walltime')        ) call cline%set('walltime',         29*60) ! 29 minutes
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
        if(envlen > 0)  call cline%set('nthr', str2int(pick_nthr_env))
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
         ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'pick_extract')
            call cline%set('projfile', 'pick_extract.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params%split_mode = 'stream'
        params%ncunits    = params%nparts
        call simple_getcwd(cwd_job)
        call cline%set('mkdir', 'no')
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, "")
       ! call nice_communicator%cycle()
        ! picking
        l_multipick   = cline%defined('nmoldiams')
        l_interactive = params%interactive == 'yes'
        ! http communicator init
        if(l_interactive) then
            call http_communicator%create(params%niceprocid, params%niceserver, "initial_picking")
        else
            call http_communicator%create(params%niceprocid, params%niceserver, "pick_extract")   
        endif
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
        endif
        if( l_multipick )then
            interactive_waiting  = .false.
            l_extract            = .false.
            l_templates_provided = .false.
            write(logfhandle,'(A)')'>>> PERFORMING MULTI-DIAMETER PICKING'
            moldiams = equispaced_vals(params%moldiam, params%moldiam_max, params%nmoldiams)
            call histogram_moldiams%new(moldiams)
            deallocate(moldiams)
            ! store diameters
            if(.not. allocated(active_search_diameters)) allocate(active_search_diameters(0))
            do i=1, histogram_moldiams%get_nbins()
                active_search_diameters = [active_search_diameters, int(histogram_moldiams%get_x(i))]
            end do
            ! remove existing files (restart)
            if(file_exists("micrographs.star"))      call del_file("micrographs.star")
            if(file_exists("micrographs_init.star")) call del_file("micrographs_init.star")
            if(file_exists("pick.star"))             call del_file("pick.star")
            if(file_exists("pick_init.star"))        call del_file("pick_init.star")
        else
            l_extract            = .true.
            l_templates_provided = cline%defined('pickrefs')
            if( l_templates_provided )then
                if( .not.file_exists(params%pickrefs) ) then
                    if(params%clear .eq. "yes") then
                        write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 120 MINUTES FOR '//trim(params%pickrefs)
                        do i=1, 720
                            if(file_exists(trim(params%pickrefs))) exit
                            call sleep(10)
                        end do
                    else
                        THROW_HARD('Could not find: '//trim(params%pickrefs))
                    end if
                endif
                if(file_exists(swap_suffix(trim(params%pickrefs), TXT_EXT, STK_EXT))) then
                    call boxsize_file%new(swap_suffix(trim(params%pickrefs), TXT_EXT, STK_EXT), 1, 1)
                    allocate(saved_boxsize(1))
                    call boxsize_file%readNextDataLine(saved_boxsize)
                    call cline%set('box', int(saved_boxsize(1)))
                    write(logfhandle,'(A)')'>>> USING SAVED BOXSIZE OF ' // int2str(int(saved_boxsize(1)))
                    deallocate(saved_boxsize)
                endif
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
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
        endif
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('micrographs', 'micrographs_imported', int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs_rejected', int2str(0), primary=.true.)
        call gui_stats%set('micrographs', 'micrographs_picked',   int2str(0), primary=.true.)
        if( l_multipick ) then
            call gui_stats%set('current_search', 'type', 'initial global search')
            call gui_stats%set('current_search', 'range', int2str(floor(params%moldiam)) // 'Å - ' // int2str(floor(params%moldiam_max)) // 'Å')
            call gui_stats%set('current_search', 'status', 'running')
        endif
        call gui_stats%set('compute', 'compute_in_use', int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! nice
       ! call nice_communicator%update_pick(micrographs_imported=0, micrographs_picked=0, micrographs_rejected=0)
       ! call nice_communicator%cycle()
        ! directories structure & restart
        odir                       = trim(DIR_STREAM)
        odir_completed             = trim(DIR_STREAM_COMPLETED)
        odir_picker                = trim(PATH_HERE)//trim(DIR_PICKER)
        odir_interactive           = 'interactive/'//trim(DIR_STREAM)
        odir_interactive_picker    = 'interactive/'//trim(DIR_PICKER)
        odir_interactive_completed = 'interactive/'//trim(DIR_STREAM_COMPLETED)
        odir_extract               = trim(PATH_HERE)//trim(DIR_EXTRACT)
        pick_extract_set_counter = 0    ! global counter of projects to be processed
        nptcls_glob              = 0    ! global number of particles
        nmics_rejected_glob      = 0    ! global number of micrographs rejected
        nmic_star                = 0
        thumbid_offset           = 0
        thumbcount               = 0
        if( cline%defined('dir_exec') )then
            call del_file(TERM_STREAM)
            call cline%delete('dir_exec')
            call simple_rmdir(odir)
            ! nice stats
            nice_communicator%stat_root%stage = "importing previously processed data"
            call nice_communicator%cycle()
            if( l_multipick .or. params%clear .eq. "yes")then
                ! removes directory structure
                call simple_rmdir(odir_completed)
                call simple_rmdir(odir_picker)
                call simple_rmdir(odir_extract)
                call simple_rmdir(odir_interactive)
                call simple_rmdir(odir_interactive_picker)
                call simple_rmdir(odir_interactive_completed)
            else
                ! import previous run and updates stats for gui
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "stage", "importing previously processed data", found)
                call http_communicator%send_jobstats()
                call import_previous_mics( projrecords )
                if( allocated(projrecords) )then
                    nptcls_glob = sum(projrecords(:)%nptcls)
                    nmic_star   = spproj_glob%os_mic%get_noris()
                    ! nice
                  !  call nice_communicator%update_pick(micrographs_imported=spproj_glob%os_mic%get_noris(), micrographs_picked=spproj_glob%os_mic%get_noris())
                   ! call nice_communicator%cycle()
                    ! http stats
                    call http_communicator%json%update(http_communicator%job_json, "micrographs_imported",  spproj_glob%os_mic%get_noris(), found)
                    call http_communicator%json%update(http_communicator%job_json, "micrographs_processed", spproj_glob%os_mic%get_noris(), found)
                    call http_communicator%json%update(http_communicator%job_json, "movies_rejected",       0,                              found)

                endif
            endif
        endif
        ! make directories structure
        call simple_mkdir(odir)
        call simple_mkdir(trim(odir)//trim(STDERROUT_DIR))
        call simple_mkdir(odir_completed)
        call simple_mkdir(odir_picker)
        if( l_multipick ) then
            call simple_mkdir(odir_interactive)
            call simple_mkdir(odir_interactive_picker)
            call simple_mkdir(odir_interactive_completed)
        endif
        if( l_extract ) call simple_mkdir(odir_extract)
        ! initialise progress monitor
        call progressfile_init()
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PICK_PARTITION, pick_part_env, envlen)
        if(envlen > 0) then
            call qenv_main%new(1,stream=.true.,qsys_partition=trim(pick_part_env))
            call qenv_interactive%new(1,stream=.true.,qsys_partition=trim(pick_part_env))
        else
            call qenv_main%new(1,stream=.true.)
            call qenv_interactive%new(1,stream=.true.)
        end if
        if(l_interactive) then
            qenv => qenv_interactive
        else
            qenv => qenv_main
        endif
        ! command line for execution
        cline_pick_extract = cline
        call cline_pick_extract%set('prg','pick_extract')
        if( l_multipick .and. l_interactive ) then
            odir           = odir_interactive
            odir_completed = odir_interactive_completed
            odir_picker    = odir_interactive_picker
            call cline_pick_extract%set('dir', PATH_PARENT)
        else
            call cline_pick_extract%set('dir', PATH_PARENT)
        endif
        if( l_extract )then
            call cline_pick_extract%set('extract','yes')
        else
            call cline_pick_extract%set('extract','no')
        endif
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
            if( file_exists(trim(TERM_STREAM)) .or. http_communicator%exit) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
              !  call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! go back to interactive multipick
            if(l_multipick .and. .not. l_interactive .and. params%interactive == 'yes') then
                call qenv%qscripts%clear_stack()
                call qenv%qscripts%get_jobs_status(jobs_done, jobs_submitted)
                if(count(jobs_done) < params%nparts) then
                    if(.not. pause_import) write(logfhandle, *) ">>>  WAITING FOR ALL PARTS TO COMPLETE"
                    pause_import = .true.
                    ! nice stats
                    nice_communicator%stat_root%stage = "waiting for running parts to complete"
                    nice_communicator%stat_root%user_input = .false.
                    call nice_communicator%cycle()
                else
                    write(logfhandle, *) ">>>  RE-ENTERING INTERACTIVE MODE"
                    if(allocated(projrecords)) then
                        do iproj=1, size(projrecords)
                            if(file_exists(trim(projrecords(iproj)%projname))) then
                                call simple_rename(trim(projrecords(iproj)%projname), filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                ! update for execution
                                call spproj_tmp%read(filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                do imic=1, spproj_tmp%os_mic%get_noris()
                                    if(file_exists(trim(spproj_tmp%os_mic%get_static(imic, 'boxfile')))) call del_file(trim(spproj_tmp%os_mic%get_static(imic, 'boxfile')))
                                end do
                              
                                call spproj_tmp%kill
                            end if
                        end do
                    end if
                    odir           = odir_interactive
                    odir_completed = odir_interactive_completed
                    odir_picker    = odir_interactive_picker
                    l_interactive  = .true.
                    params%updated        = 'no'
                    params%moldiam_refine = 0.0
                !    params%moldiam        = moldiam_interactive
                !    params%nmoldiams      = nmoldiams_interactive
                    pause_import   = .true.
                    interactive_waiting = .true.
                    call kill_projrecords(projrecords_main)
                    allocate(projrecords_main(size(projrecords)))
                    projrecords_main(:) = projrecords(:)
                    call kill_projrecords(projrecords)
                    call spproj_glob%os_mic%kill()
                    call spproj_glob%os_mic%copy(interactive_spproj%os_mic, is_ptcl=.false.)
                    call import_previous_mics( projrecords )
                    if( allocated(projrecords) )then
                        nptcls_glob = sum(projrecords(:)%nptcls)
                        nmic_star   = spproj_glob%os_mic%get_noris()
                    endif
                    qenv => qenv_interactive
                    n_imported = 0
                    n_failed_jobs = size(projrecords) - n_imported
                    thumbid_offset = 100
                    if(spproj_glob%os_mic%isthere('intg') .and. spproj_glob%os_mic%isthere('xdim') .and. spproj_glob%os_mic%isthere('ydim') \
                    .and. spproj_glob%os_mic%isthere('smpd') .and. spproj_glob%os_mic%isthere('boxfile')) then
                        i_max = 20
                        if(spproj_glob%os_mic%get_noris() < i_max) i_max = spproj_glob%os_mic%get_noris()
                        do i_thumb=0, i_max - 1
                            ! Generate jpegs 
                            i = spproj_glob%os_mic%get_noris() - i_thumb
                            ldim_mic(1) = int(spproj_glob%os_mic%get(i, 'xdim'))
                            ldim_mic(2) = int(spproj_glob%os_mic%get(i, 'ydim'))
                            ldim_mic(3) = 1
                            call mic_thumb%new(ldim_mic, spproj_glob%os_mic%get(i, 'smpd'))
                            call mic_thumb%read(trim(spproj_glob%os_mic%get_static(i, 'intg')))
                            ! resize
                            jpg_scale = 1000.0 / ldim_mic(2)
                            call mic_thumb%fft
                            call mic_thumb%clip_inplace([2 * (nint(0.5 * jpg_scale * ldim_mic(1))), 2 * (nint(0.5 * jpg_scale * ldim_mic(2))), 1]) ! ensures new img has even dims
                            call mic_thumb%ifft
                            call mic_thumb%write_jpg(trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT)
                            call mic_thumb%kill
                            if(i_thumb == 0) then
                                call nice_communicator%update_pick(carousel=.true., clear_carousel=.true., thumbnail=trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT,&
                                thumbnail_id=i + thumbid_offset, thumbnail_static_id=i_thumb + 1, boxfile=trim(spproj_glob%os_mic%get_static(i, 'boxfile')), scale=jpg_scale)  
                            else
                                call nice_communicator%update_pick(carousel=.true., thumbnail=trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT,&
                                thumbnail_id=i + thumbid_offset, thumbnail_static_id=i_thumb + 1, boxfile=trim(spproj_glob%os_mic%get_static(i, 'boxfile')), scale=jpg_scale)
                            end if
                        end do
                    end if

                end if
            end if
            if(l_interactive .and. l_multipick) then
                ! pause import if more that ninit mics imported
                if(n_added >= params%ninit .and. .not. pause_import ) then
                    write(logfhandle,'(A,A,A)') '>>> NEW MICROGRAPH IMPORT PAUSED AFTER ', int2str(params%ninit), ' MICROGRAPHS WHILE INITIAL PICKING IS PERFORMED';
                    ! http
                    call http_communicator%json%update(http_communicator%job_json, "stage", "initial search on " // int2str(project_buff%n_history * NMOVS_SET) // " micrographs", found)
                    pause_import = .true.
                end if
                ! restart pick if moldiam_refine updated in append mode
                if(interactive_waiting .and. params%updated .eq. 'yes' .and. (params%moldiam_refine .gt. 0.0 .or. params%moldiam_refine .lt. 0.0 )) then
                    if(params%moldiam_refine .gt. 0.0 ) then
                        write(logfhandle,'(A,I3)') '>>> REFINING MOLECULAR DIAMETER        : ', int(params%moldiam_refine)
                        call http_communicator%json%update(http_communicator%job_json, "stage", "refining pick diameter", found)
                    else if (params%moldiam_refine .eq. -1.0 ) then
                        write(logfhandle,'(A)') '>>> INCREASING SEARCH RANGE'
                        call http_communicator%json%update(http_communicator%job_json, "stage", "increasing search range", found)
                    else if (params%moldiam_refine .eq. -2.0 ) then
                        write(logfhandle,'(A)') '>>> DECREASING SEARCH RANGE'
                        call http_communicator%json%update(http_communicator%job_json, "stage", "decreasing search range", found)
                    end if
                     ! http
                    call http_communicator%json%update(http_communicator%job_json, "user_input", .false., found)
                    ! update params
                    if(params%moldiam_refine .eq. -2.0) then
                        !decrease search range
                        params%moldiam_max = minval(complete_search_diameters) - 20
                        params%moldiam     = 20
                        params%nmoldiams   = 4
                    else if(params%moldiam_refine .eq. -1.0) then
                        !increase search range
                        params%moldiam     = maxval(complete_search_diameters) + 100
                        params%moldiam_max = params%moldiam + 100
                        params%nmoldiams   = 2
                    else
                        params%moldiam     = params%moldiam_refine - 40.0
                        params%moldiam_max = params%moldiam_refine + 40.0
                        params%nmoldiams   = 5
                    end if
                    !reset
                    params%updated        = 'no'
                    params%moldiam_refine = 0.0
                    ! update cline for picker
                    call cline_pick_extract%set('moldiam',     params%moldiam)
                    call cline_pick_extract%set('moldiam_max', params%moldiam_max)
                    call cline_pick_extract%set('nmoldiams',   params%nmoldiams)
                    !call cline_pick_extract%set('dir','../')
                    call cline_pick_extract%set('append','yes')
                    ! undo already processed micrographs
                    call spproj_glob%os_mic%kill()
                    if(allocated(projrecords)) then
                        do iproj=1, size(projrecords)
                            if(file_exists(trim(projrecords(iproj)%projname))) then
                                call simple_rename(trim(projrecords(iproj)%projname), filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                ! update for execution
                                call spproj_tmp%read(filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                call cline_pick_extract%set('projname', trim(spproj_tmp%projinfo%get_static(1, 'projname')))
                                call cline_pick_extract%set('projfile', basename(trim(projrecords(iproj)%projname)))
                                call cline_pick_extract%set('fromp',    1)
                                call cline_pick_extract%set('top',      spproj_tmp%os_mic%get_noris())
                                call spproj_tmp%kill
                                call qenv%qscripts%add_to_streaming(cline_pick_extract)
                                call qenv%qscripts%schedule_streaming( qenv%qdescr, path=odir )
                            end if
                        end do
                        deallocate(projrecords)
                    end if
                    n_imported = 0
                    ! reset histogram
                    moldiams = equispaced_vals(params%moldiam, params%moldiam_max, params%nmoldiams)
                    call histogram_moldiams%kill
                    call histogram_moldiams%new(moldiams)
                    deallocate(moldiams)
                    ! store search diameters
                    if(allocated(active_search_diameters))  deallocate(active_search_diameters)
                    if(allocated(refined_search_diameters)) deallocate(refined_search_diameters)
                    allocate(active_search_diameters(0))
                    do i=1, histogram_moldiams%get_nbins()
                        active_search_diameters = [active_search_diameters, int(histogram_moldiams%get_x(i))]
                    end do
                    ! keep import paused
                    pause_import = .true.
                    interactive_waiting = .false.
                    thumbid_offset = 100
                end if
                ! exit interactive mode
                if(params%updated .eq. 'yes' .and. params%interactive == 'no') then
                    ! undo already processed micrographs
                    write(logfhandle, '(A)') ">>> EXITING INTERACTIVE MODE"
                    write(logfhandle, '(A)') ">>> PICKING USING :"
                    call interactive_spproj%os_mic%copy(spproj_glob%os_mic, is_ptcl=.false.)
                    call spproj_glob%os_mic%kill()
                    call project_buff%clear_history()
                    qenv => qenv_main
                    !reset
                    params%updated        = 'no'
                    n_imported = 0
                    ! reset paths
                    odir            = trim(DIR_STREAM)
                    odir_completed  = trim(DIR_STREAM_COMPLETED)
                    odir_picker     = trim(PATH_HERE)//trim(DIR_PICKER)
                    l_interactive = .false.
                    pause_import  = .false.
                    interactive_waiting = .false.
                    params%nmoldiams = 1.0
                    if(params%moldiam .gt. 0.0)then
                        write(logfhandle, '(A,F12.2)') "        GAUSSIAN DIAMETER : ", params%moldiam
                        call cline_pick_extract%set('moldiam',    params%moldiam)
                    else
                        call cline_pick_extract%set('moldiam',    0.0)
                    end if
                    if(params%ring == 'yes')then
                        write(logfhandle, '(A,F12.2)') "        RING DIAMETER : ", params%moldiam_ring
                        call cline_pick_extract%set('moldiam_ring',    params%moldiam_ring)
                    else
                        call cline_pick_extract%delete('moldiam_ring')
                    end if
                    call cline_pick_extract%set('nmoldiams',  params%nmoldiams)
                    call cline_pick_extract%set('ring',       params%ring)
                    call cline_pick_extract%delete('moldiam_max')
                    nice_communicator%stat_root%stage = "picking with user selected parameters"
                    nice_communicator%stat_root%user_input = .false.
                    call nice_communicator%update_pick(gaussian_diameter=int(params%moldiam_refine))
                    ! http
                    call http_communicator%json%update(http_communicator%job_json, "stage", "picking with selected parameters", found)
                    call http_communicator%json%update(http_communicator%job_json, "user_input", .false., found)
                    call http_communicator%json%remove(latest_picked_micrographs, destroy=.true.)
                    call http_communicator%json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                    call http_communicator%json%add(http_communicator%job_json, latest_picked_micrographs)
                    call kill_projrecords(projrecords)
                    if(allocated(projrecords_main)) then
                        allocate(projrecords(size(projrecords_main)))
                        projrecords(:) = projrecords_main(:)
                    end if
                    call qenv%qscripts%clear_stack()
                    if(allocated(projrecords)) then
                        do iproj=1, size(projrecords)
                            if(file_exists(trim(projrecords(iproj)%projname))) then
                                call simple_rename(trim(projrecords(iproj)%projname), filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                ! update for execution
                                call spproj_tmp%read(filepath(odir, basename(trim(projrecords(iproj)%projname))))
                                do imic=1, spproj_tmp%os_mic%get_noris()
                                    if(file_exists(trim(spproj_tmp%os_mic%get_static(imic, 'boxfile')))) call del_file(trim(spproj_tmp%os_mic%get_static(imic, 'boxfile')))
                                end do
                                call spproj_tmp%kill
                            end if
                        end do
                        call kill_projrecords(projrecords)
                    end if
                    thumbid_offset = 0
                end if
            end if
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
                            ! nice stats
                            nice_communicator%stat_root%stage = "preparing picking references"
                            call nice_communicator%cycle()
                            cline_make_pickrefs = cline
                            call cline_make_pickrefs%set('prg',   'make_pickrefs')
                            call cline_make_pickrefs%set('stream','no')
                            call cline_make_pickrefs%set('smpd',  params%smpd) ! is the required output sampling distance
                            call xmake_pickrefs%execute_safe(cline_make_pickrefs)
                            call cline_pick_extract%set('pickrefs', '../'//trim(PICKREFS_FBODY)//trim(params%ext))
                            ! nice
                            call mrc2jpeg_tiled(trim(PICKREFS_FBODY)//trim(params%ext), trim(PICKREFS_FBODY)//".jpeg", scale=pickrefs_thumbnail_scale, ntiles=n_pickrefs, n_xtiles=xtiles, n_ytiles=ytiles)
                            call random_number(rnd)
                            call nice_communicator%update_pick(pickrefs_thumbnail=trim(cwd_job)//'/'//trim(PICKREFS_FBODY)//".jpeg", pickrefs_thumbnail_id=floor(100000 * rnd), pickrefs_thumbnail_static_id=1000, pickrefs_thumbnail_n_tiles=n_pickrefs, pickrefs_thumbnail_scale=pickrefs_thumbnail_scale)
                            ! http stats
                            xtile = 0
                            ytile = 0
                            do i=0, n_pickrefs - 1
                                call communicator_add_picking_template(trim(cwd_job)//'/'//trim(PICKREFS_FBODY)//".jpeg", xtile * 100, ytile * 100, 100 * ytiles, 100 * xtiles)
                                xtile = xtile + 1
                                if(xtile .eq. xtiles) then
                                    xtile = 0
                                    ytile = ytile + 1
                                endif
                            end do
                            write(logfhandle,'(A)')'>>> PREPARED PICKING TEMPLATES'
                            nice_communicator%stat_root%stage = "picking references prepared"
                            call nice_communicator%cycle()
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
              !  call gui_stats%set('micrographs', 'micrographs_imported', int2commastr(project_buff%n_history * NMOVS_SET), primary=.true.)
               ! call gui_stats%set_now('micrographs', 'last_micrograph_imported')
                ! nice
               ! call nice_communicator%update_pick(micrographs_imported=project_buff%n_history * NMOVS_SET, last_micrograph_imported=.true.)
              !  call nice_communicator%cycle()
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "micrographs_imported",     project_buff%n_history * NMOVS_SET, found)
                call http_communicator%json%update(http_communicator%job_json, "last_micrograph_imported", stream_datestr(), found)
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
                call update_projects_list( projrecords, n_imported )
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
                    call histogram_moldiams%plot('moldiams', 'Diameters')
                    if( l_multipick .and. l_interactive)       call starproj_stream%stream_export_pick_diameters(params%outdir, histogram_moldiams, filename="pick_init.star")
                    if( l_multipick .and. .not. l_interactive) call starproj_stream%stream_export_pick_diameters(params%outdir, histogram_moldiams)
                    write(logfhandle,'(A,F8.2)') '>>> ESTIMATED MOLECULAR DIAMETER        : ',histogram_moldiams%mean()
                    call gui_stats%set('current_search', 'estimated_diameter', int2str(floor(histogram_moldiams%mean())) // 'Å')
                    call nice_communicator%update_pick(suggested_diameter=floor(histogram_moldiams%mean()))
                endif
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
              !  call gui_stats%set('micrographs', 'micrographs_picked', int2commastr(n_imported) // ' (' // int2str(ceiling(100.0 * real(n_imported) / real(project_buff%n_history * NMOVS_SET))) // '%)', primary=.true.)
              !  call gui_stats%set('micrographs', 'micrographs_rejected', int2commastr(n_failed_jobs + nmics_rejected_glob) // ' (' // int2str(floor(100.0 * real(n_failed_jobs + nmics_rejected_glob) / real(project_buff%n_history * NMOVS_SET))) // '%)', primary=.true.)
                ! nice
              !  call nice_communicator%update_pick(micrographs_picked=n_imported, micrographs_rejected=n_failed_jobs + nmics_rejected_glob)
              !  if(spproj_glob%os_mic%isthere("nptcls") .and. .not. l_multipick) then
              !      call gui_stats%set('micrographs', 'avg_number_picks', ceiling(spproj_glob%os_mic%get_avg("nptcls")), primary=.true.)
              !  end if
            !    if(.not. l_multipick .and. .not. l_interactive) then
                if(.not. l_interactive) then
                    ! http stats   
                    call http_communicator%json%update(http_communicator%job_json, "stage",               "finding, picking and extracting micrographs", found)  
                    call http_communicator%json%update(http_communicator%job_json, "box_size",              params%box,                                  found)
                    call http_communicator%json%update(http_communicator%job_json, "particles_extracted",   nptcls_glob,                                 found)
                    call http_communicator%json%update(http_communicator%job_json, "micrographs_processed", n_imported,                                  found)
                    call http_communicator%json%update(http_communicator%job_json, "micrographs_rejected",  n_failed_jobs + nmics_rejected_glob,         found)
                    if(spproj_glob%os_mic%isthere('thumb') .and. spproj_glob%os_mic%isthere('xdim') .and. spproj_glob%os_mic%isthere('ydim') \
                    .and. spproj_glob%os_mic%isthere('smpd') .and. spproj_glob%os_mic%isthere('boxfile')) then
                        i_max = 10
                        if(spproj_glob%os_mic%get_noris() < i_max) i_max = spproj_glob%os_mic%get_noris()
                        ! create an empty latest_picked_micrographs json array
                        call http_communicator%json%remove(latest_picked_micrographs, destroy=.true.)
                        call http_communicator%json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                        call http_communicator%json%add(http_communicator%job_json, latest_picked_micrographs)
                        do i_thumb=0, i_max - 1
                            call communicator_add_micrograph(trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris() - i_thumb, "thumb"))),\
                            100, 0, 100, 200,\
                            nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "xdim")),\
                            nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "ydim")),\
                            trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris() - i_thumb, "boxfile"))))
                        end do
                    endif
                end if
              !  if(spproj_glob%os_mic%isthere('intg') .and. spproj_glob%os_mic%isthere('boxfile')) then
               !     latest_boxfile = trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'boxfile'))
                 !   if(file_exists(latest_boxfile)) call gui_stats%set('latest', '', trim(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris(), 'intg')), thumbnail=.true., boxfile=trim(latest_boxfile))
               ! end if
                !nice
                if(spproj_glob%os_mic%isthere('intg') .and. spproj_glob%os_mic%isthere('xdim') .and. spproj_glob%os_mic%isthere('ydim') \
                    .and. spproj_glob%os_mic%isthere('smpd') .and. spproj_glob%os_mic%isthere('boxfile')) then
                    i_max = 20
                    if(spproj_glob%os_mic%get_noris() < i_max) then
                        i_max = spproj_glob%os_mic%get_noris()
                        thumbcount = 0
                    end if
                    if(spproj_glob%os_mic%get_noris() >= thumbcount + i_max) then
                        do i_thumb=0, i_max - 1
                            ! Generate jpegs 
                            i = spproj_glob%os_mic%get_noris() - i_thumb
                            if(trim(spproj_glob%os_mic%get_static(i, 'boxfile')) .eq. "") cycle
                            ldim_mic(1) = int(spproj_glob%os_mic%get(i, 'xdim'))
                            ldim_mic(2) = int(spproj_glob%os_mic%get(i, 'ydim'))
                            ldim_mic(3) = 1
                            call mic_thumb%new(ldim_mic, spproj_glob%os_mic%get(i, 'smpd'))
                            call mic_thumb%read(trim(spproj_glob%os_mic%get_static(i, 'intg')))
                            ! resize
                            jpg_scale = 1000.0 / ldim_mic(2)
                            call mic_thumb%fft
                            call mic_thumb%clip_inplace([2 * (nint(0.5 * jpg_scale * ldim_mic(1))), 2 * (nint(0.5 * jpg_scale * ldim_mic(2))), 1]) ! ensures new img has even dims
                            call mic_thumb%ifft
                            call mic_thumb%write_jpg(trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT)
                            call mic_thumb%kill
                            if(i_thumb == 0) then
                                call nice_communicator%update_pick(carousel=.true., clear_carousel=.true., thumbnail=trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT,&
                                thumbnail_id=i + thumbid_offset, thumbnail_static_id=i_thumb + 1, boxfile=trim(spproj_glob%os_mic%get_static(i, 'boxfile')), scale=jpg_scale)  
                            else
                                call nice_communicator%update_pick(carousel=.true., thumbnail=trim(cwd_job) // '/' // int2str(i_thumb) // THUMBNAIL_SUFFIX // JPG_EXT,&
                                thumbnail_id=i + thumbid_offset, thumbnail_static_id=i_thumb + 1, boxfile=trim(spproj_glob%os_mic%get_static(i, 'boxfile')), scale=jpg_scale)
                            end if
                        end do
                        thumbcount = thumbcount + i_max
                    end if
                end if
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                last_injection = simple_gettime()
                l_haschanged = .true.
                ! always write micrographs snapshot if less than 1000 mics, else every 100
                if( n_imported < 1000 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_micrographs_starfile
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
                    call update_user_params(cline)
                    call write_micrographs_starfile
                    nmic_star = n_imported
                endif
            else
                ! write snapshot
                if( .not.l_projects_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive
                        call write_project
                        call update_user_params(cline)
                        call write_micrographs_starfile
                        l_haschanged = .false.
                    endif
                endif
            endif
           ! if(l_multipick .and. l_interactive .and. pause_import .and. n_imported + n_failed_jobs >= params%ninit) then
           if(l_multipick .and. l_interactive .and. pause_import .and. spproj_glob%os_mic%get_noris() >= params%ninit) then
                if(.not. interactive_waiting) write(logfhandle,'(A)')'>>> WAITING FOR USER INPUT'
           !     nice_communicator%stat_root%user_input = .true.
          !      call gui_stats%set('current_search', 'status', 'waiting for refinement diameter')
                ! nice stats
         !       nice_communicator%stat_root%stage = "waiting for user input"
         !       call nice_communicator%cycle()
         !       nice_communicator%stat_root%user_input = .true.
                ! search diameters
                if(allocated(active_search_diameters)) then
                    if(.not. allocated(complete_search_diameters)) allocate(complete_search_diameters(0))
                    if(allocated(refined_search_diameters))        deallocate(refined_search_diameters)
                    complete_search_diameters = [complete_search_diameters, active_search_diameters]
                    allocate(refined_search_diameters(size(active_search_diameters)))
                    refined_search_diameters(:) = active_search_diameters(:)
                    deallocate(active_search_diameters)
                end if
                ! http stats 
                call http_communicator%json%update(http_communicator%job_json, "stage", "waiting for user input", found)
                call http_communicator%json%update(http_communicator%job_json, "user_input", .true., found)
                call http_communicator%json%update(http_communicator%job_json, "micrographs_processed", n_imported,                                  found)
                call http_communicator%json%update(http_communicator%job_json, "micrographs_rejected",  n_failed_jobs + nmics_rejected_glob,         found)
                ! create an empty picking_diameters json array
                call http_communicator%json%remove(picking_diameters, destroy=.true.)
                call http_communicator%json%create_array(picking_diameters, "picking_diameters")
                call http_communicator%json%add(http_communicator%job_json, picking_diameters)
                call hpsort(complete_search_diameters)
                do i=1, size(complete_search_diameters)
                    call http_communicator%json%add(picking_diameters, "", int(complete_search_diameters(i)))
                enddo
                if(spproj_glob%os_mic%isthere('thumb') .and. spproj_glob%os_mic%isthere('xdim') .and. spproj_glob%os_mic%isthere('ydim') \
                    .and. spproj_glob%os_mic%isthere('smpd') .and. spproj_glob%os_mic%isthere('boxfile')) then
                    i_max = 10
                    if(spproj_glob%os_mic%get_noris() < i_max) i_max = spproj_glob%os_mic%get_noris()
                    ! create an empty latest_picked_micrographs json array
                    call http_communicator%json%remove(latest_picked_micrographs, destroy=.true.)
                    call http_communicator%json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                    call http_communicator%json%add(http_communicator%job_json, latest_picked_micrographs)
                    do i_thumb=0, i_max - 1
                        call communicator_add_micrograph(trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris() - i_thumb, "thumb"))),\
                        100, 0, 100, 200,\
                        nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "xdim")),\
                        nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "ydim")),\
                        trim(adjustl(spproj_glob%os_mic%get_static(spproj_glob%os_mic%get_noris() - i_thumb, "boxfile"))))
                    end do
                endif
                interactive_waiting = .true.
            endif
          !  call update_user_params(cline, nice_communicator%update_arguments)
            ! guistats
            !call gui_stats%write_json
            if(interactive_waiting) then
                call sleep(1)
            else
                call sleep(WAITTIME)
            end if
          !  call nice_communicator%cycle()
            ! http stats send
            call http_communicator%send_jobstats()
            call update_user_params(cline, http_communicator%update_arguments)
        end do
        ! termination
        call write_project
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        call update_user_params(cline)
        call starproj_stream%copy_micrographs_optics(spproj_glob, write=.true., verbose=DEBUG_HERE)
        call write_micrographs_starfile(optics_set=.true.)
        call starproj_stream%stream_export_particles_2D(spproj_glob, params%outdir, optics_set=.true., verbose=.true.)
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
        call kill_projrecords(projrecords)
        call kill_projrecords(projrecords_main)
        ! end gracefully
        call nice_communicator%terminate()
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT NORMAL STOP ****')
        contains

            !>  write starfile snapshot
            subroutine write_micrographs_starfile( optics_set )
                logical, optional, intent(in) :: optics_set
                integer(timer_int_kind)       :: ms0
                real(timer_int_kind)          :: ms_export
                logical                       :: l_optics_set
                l_optics_set = .false.
                if( present(optics_set) ) l_optics_set = optics_set
                if (spproj_glob%os_mic%get_noris() > 0) then
                    if( DEBUG_HERE ) ms0 = tic()
                    if(l_multipick .and. l_interactive) then
                        call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set, filename="micrographs_init.star")
                    else
                        call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set)
                    endif
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                end if
            end subroutine write_micrographs_starfile

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
                        if( trim(projrecords(imic)%projname) /= trim(prev_projname) )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('stk', projrecords(imic)%projname)
                            prev_projname = trim(projrecords(imic)%projname)
                        endif
                        micind       = projrecords(imic)%micind
                        fromps(imic) = stream_spproj%os_stk%get_fromp(micind) ! fromp from individual project
                        fromp        = nptcls + 1
                        nptcls       = nptcls + projrecords(imic)%nptcls
                        top          = nptcls
                        call spproj_glob%os_stk%transfer_ori(imic, stream_spproj%os_stk, micind)
                        call spproj_glob%os_stk%set(imic, 'fromp',fromp)
                        call spproj_glob%os_stk%set(imic, 'top',  top)
                    enddo
                    call spproj_glob%write_segment_inside('stk', params%projfile)
                    ! particles
                    call spproj_glob%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                    iptcl         = 0
                    prev_projname = ''
                    do imic = 1,nmics
                        if( trim(projrecords(imic)%projname) /= prev_projname )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('ptcl2D', projrecords(imic)%projname)
                            prev_projname = trim(projrecords(imic)%projname)
                        endif
                        fromp = fromps(imic)
                        top   = fromp + projrecords(imic)%nptcls - 1
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
                type(projrecord), allocatable, intent(inout) :: records(:)
                integer,                       intent(inout) :: nimported
                type(sp_project), allocatable :: spprojs(:)
                type(projrecord), allocatable :: old_records(:)
                character(len=:), allocatable :: fname, abs_fname, new_fname
                type(sp_project)              :: tmpproj
                integer :: n_spprojs, n_old, j, nprev_imports, n_completed, nptcls, nmics, imic, iproj
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
                        allocate(projrecords(nmics))
                    else
                        call spproj_glob%os_mic%reallocate(n_completed)
                        call move_alloc(projrecords, old_records)
                        allocate(projrecords(n_completed))
                        if( n_old > 0 ) projrecords(1:n_old) = old_records(:)
                        call kill_projrecords(old_records)
                    endif
                    ! update records and global project
                    j = n_old
                    do iproj = 1,n_spprojs
                        if( spprojs(iproj)%os_mic%get_noris() == 0 ) cycle
                        ! move project to appropriate directory
                        fname = filepath(odir, completed_jobs_clines(iproj)%get_carg('projfile'))
                        new_fname = filepath(odir_completed, completed_jobs_clines(iproj)%get_carg('projfile'))
                        call simple_rename(fname, new_fname)
                        abs_fname = simple_abspath(new_fname, errmsg='stream pick_extract :: update_projects_list 1')
                        ! records & project
                        do imic = 1,spprojs(iproj)%os_mic%get_noris()
                            j = j + 1
                            projrecords(j)%projname = trim(abs_fname)
                            projrecords(j)%micind   = imic
                            if( l_multipick )then
                                projrecords(j)%nptcls = 0
                            else
                                nptcls      = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                                nptcls_glob = nptcls_glob + nptcls ! global update
                                projrecords(j)%nptcls = nptcls
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
                if(l_multipick .and. l_interactive) call spproj_here%projinfo%set(1,'init', 1.0)
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
                type(projrecord), allocatable, intent(inout) :: records(:)
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
                            records(irec)%nptcls   = spprojs(iproj)%os_mic%get_int(imic, 'nptcls')
                            call spproj_glob%os_mic%transfer_ori(irec, spprojs(iproj)%os_mic, imic)
                        endif
                    enddo
                    call spprojs(iproj)%kill
                enddo
                ! update global set counter
                do iproj = 1,n_spprojs
                    fname = basename_safe(completed_fnames(iproj))
                    fname = trim(get_fbody(trim(fname),trim(METADATA_EXT),separator=.false.))
                    id    = str2int(fname, iostat)
                    if( iostat==0 ) pick_extract_set_counter = max(pick_extract_set_counter, id)
                enddo
                nmics_rejected_glob = nmics - nsel_mics
                ! add previous projects to history
                do iproj = 1,n_spprojs
                    call project_buff%add2history(completed_fnames(iproj))
                enddo
                write(logfhandle,'(A,I6,A)')'>>> IMPORTED ',nsel_mics,' PREVIOUSLY PROCESSED MICROGRAPHS'
            end subroutine import_previous_mics

            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",               "initialising")
                call http_communicator%json%add(http_communicator%job_json, "micrographs_imported",     0)
                call http_communicator%json%add(http_communicator%job_json, "micrographs_processed",    0)
                call http_communicator%json%add(http_communicator%job_json, "micrographs_rejected",     0)
                call http_communicator%json%add(http_communicator%job_json, "particles_extracted",      dble(0.0))
                call http_communicator%json%add(http_communicator%job_json, "box_size",                 dble(0.0))
                call http_communicator%json%add(http_communicator%job_json, "user_input",               .false.)
                call http_communicator%json%add(http_communicator%job_json, "last_micrograph_imported", "")
                call http_communicator%json%create_array(latest_extracted_particles, "latest_extracted_particles")
                call http_communicator%json%add(http_communicator%job_json, latest_extracted_particles)
                call http_communicator%json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                call http_communicator%json%add(http_communicator%job_json, latest_picked_micrographs)
                call http_communicator%json%create_array(picking_templates, "picking_templates")
                call http_communicator%json%add(http_communicator%job_json, picking_templates)
                call http_communicator%json%create_array(picking_diameters, "picking_diameters")
                call http_communicator%json%add(http_communicator%job_json, picking_diameters)
            end subroutine communicator_init

            subroutine communicator_add_micrograph(path, spritex, spritey, spriteh, spritew, xdim, ydim, boxfile_path)
                character(*),     intent(in)  :: path, boxfile_path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew, xdim, ydim
                type(nrtxtfile)               :: boxfile
                type(json_value), pointer     :: micrograph, boxes, box
                real,             allocatable :: boxdata(:,:)
                integer                       :: i, x, y, diameter, type
                call http_communicator%json%create_object(micrograph, "")
                call http_communicator%json%add(micrograph, "path",    path)
                call http_communicator%json%add(micrograph, "spritex", spritex)
                call http_communicator%json%add(micrograph, "spritey", spritey)
                call http_communicator%json%add(micrograph, "spriteh", spriteh)
                call http_communicator%json%add(micrograph, "spritew", spritew)
                call http_communicator%json%add(micrograph, "xdim"   , xdim)
                call http_communicator%json%add(micrograph, "ydim",    ydim)
                call http_communicator%json%create_array(boxes, "boxes")
                call boxfile%new(boxfile_path, 1)
                allocate(boxdata(boxfile%get_nrecs_per_line(), boxfile%get_ndatalines()))
                if(boxfile%get_nrecs_per_line() == 5) then
                    ! standard boxfile
                    do i=1, boxfile%get_ndatalines()
                        call boxfile%readNextDataLine(boxdata(:,i))
                        call http_communicator%json%create_object(box, "")
                        x = nint(boxdata(1,i) + boxdata(3,i)/2)
                        y = nint(boxdata(2,i) + boxdata(4,i)/2)
                        call http_communicator%json%add(box, "x",    x)
                        call http_communicator%json%add(box, "y",    y)
                        call http_communicator%json%add(boxes, box)
                    enddo
                else if(boxfile%get_nrecs_per_line() == 6) then
                    ! multipick boxfile
                    do i=1, boxfile%get_ndatalines()
                        call boxfile%readNextDataLine(boxdata(:,i))
                        call http_communicator%json%create_object(box, "")
                        x = nint(boxdata(1,i) + boxdata(3,i)/2)
                        y = nint(boxdata(2,i) + boxdata(3,i)/2)
                        diameter = floor(boxdata(4,i))
                        type     = nint(boxdata(6,i))
                        call http_communicator%json%add(box, "x",        x)
                        call http_communicator%json%add(box, "y",        y)
                        call http_communicator%json%add(box, "diameter", diameter)
                        call http_communicator%json%add(box, "type",     type)
                        call http_communicator%json%add(boxes, box)
                    enddo
                endif
                call boxfile%kill()
                deallocate(boxdata)
                call http_communicator%json%add(micrograph, boxes)
                call http_communicator%json%add(latest_picked_micrographs, micrograph)
            end subroutine communicator_add_micrograph

            subroutine communicator_add_picking_template(path, spritex, spritey, spriteh, spritew)
                character(*),     intent(in)  :: path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew
                type(json_value), pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", spritex)
                call http_communicator%json%add(template, "spritey", spritey)
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(picking_templates, template)
            end subroutine communicator_add_picking_template

    end subroutine exec_stream_pick_extract

    subroutine exec_stream_gen_picking_refs( self, cline )
        class(commander_stream_gen_picking_refs), intent(inout) :: self
        class(cmdline),                           intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        integer,                   parameter   :: MAXPOP_DEFAULT  = 200000 ! # of particles after wich picking is stopped
        logical,                   parameter   :: DEBUG_HERE      = .FALSE.
        type(parameters)                       :: params
        type(stream_http_communicator)         :: http_communicator
        type(json_value),          pointer     :: latest_cls2D, selected_references
        type(simple_nice_communicator)         :: nice_communicator
        type(guistats)                         :: gui_stats
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(projrecord),          allocatable :: projrecords(:)
        type(qsys_env)                         :: qenv
        type(cmdline)                          :: cline_extract, cline_pick_extract
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(starproject_stream)               :: starproj_stream
        type(json_core)                        :: json
        type(nrtxtfile)                        :: boxsize_file
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=:),          allocatable :: output_dir, output_dir_extract, output_dir_picker, final_selection_source
        character(len=LONGSTRLEN)              :: cwd_job
        character(len=STDLEN)                  :: refgen_nthr_env, refgen_part_env
        integer                                :: extract_set_counter    ! Internal counter of projects to be processed
        integer                                :: i, nmics_sel, nmics_rej, nmics_rejected_glob, final_selection_boxsize
        integer                                :: nmics, nprojects, stacksz, prev_stacksz, iter, iproj, envlen
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, nmic_star
        integer                                :: xtile, ytile, nxtiles, nytiles
        integer,                   allocatable :: final_selection(:), final_boxsize(:)
        logical                                :: l_haschanged, l_once, l_pick_extract, l_user, json_found
        logical(LK)                            :: found 
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
        ! 2D analysis
        if( .not. cline%defined('ncls')            ) call cline%set('ncls',           30)
        if( .not. cline%defined('nptcls_per_cls')  ) call cline%set('nptcls_per_cls', 500)
        if( .not. cline%defined('ml_reg')          ) call cline%set('ml_reg',         'no')
        if( .not. cline%defined('refine')          ) call cline%set('refine',         'snhc_smpl')
        if( .not. cline%defined('maxpop')          ) call cline%set('maxpop',         MAXPOP_DEFAULT)
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_REFGEN_NTHR, refgen_nthr_env, envlen)
        if(envlen > 0) then
            call cline%set('nthr',   str2int(refgen_nthr_env))
            call cline%set('nthr2D', str2int(refgen_nthr_env))
        end if
        ! sanity check for restart
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                THROW_HARD('Previous directory does not exists: '//trim(cline%get_carg('dir_exec')))
            endif
        endif
        ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'gen_pick_refs')
            call cline%set('projfile', 'gen_pick_refs.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
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
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "generate_picking_refs")
        call communicator_init()
        call http_communicator%send_jobstats()
         ! write cmdline for GUI
        call cline%writeline(".cline")
        if( params%nptcls_per_cls*params%ncls > params%maxpop )then
            THROW_HARD('Incorrect inputs: MAXPOP < NPTCLS_PER_CLS x NCLS')
        endif
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
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
        ! nice communicator init
        call nice_communicator%init(1,'')
        call nice_communicator%cycle()
        ! master project file
        call spproj_glob%read( params%projfile )
        ! force local if nparts < 2
        if(params%nparts < 2) then
            call spproj_glob%compenv%set(1, 'qsys_name', 'local')
            call spproj_glob%write_non_data_segments(trim(params%projfile))
        endif
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! 2D parameters
        params%nchunks     = 0
        params%ncls_start  = params%ncls
        params%ml_reg_pool = params%ml_reg
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('particles', 'particles_extracted',         0,            primary=.true.)
        call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        call gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! nice
        call nice_communicator%update_cls2D(particles_extracted=0)
        call nice_communicator%cycle()
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
            nice_communicator%stat_root%stage = "erasing previous run data"
            call nice_communicator%cycle()
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
        else
            ! extraction only
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
        endif
        ! Infinite loop
        prev_stacksz          = 0
        iter                  = 0
        n_imported            = 0   ! global number of imported processed micrographs
        n_failed_jobs         = 0
        n_added               = 0   ! global number of micrographs added to processing stack
        l_haschanged          = .false.
        l_once                = .true.
        l_user                = .false.
        do
            if( file_exists(trim(TERM_STREAM)) .or. file_exists(STREAM_REJECT_CLS) .or. http_communicator%exit)then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
            !    call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM GEN_PICKING_REFERENCES USER STOP ****')
                call EXIT(0)
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
                !call gui_stats%set_now('particles', 'last_particles_imported')
                !call nice_communicator%update_cls2D(last_particles_imported=.true.)
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "last_particles_imported", stream_datestr(), found)    
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
                call update_records_with_project( projrecords, n_imported )
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
            ! project update
            if( n_imported > 0 )then
                n_imported = spproj_glob%os_mic%get_noris()
                write(logfhandle,'(A,I8)')       '>>> # MICROGRAPHS PROCESSED & IMPORTED  : ',n_imported
                write(logfhandle,'(A,I8)')       '>>> # PARTICLES EXTRACTED               : ',nptcls_glob
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                !call gui_stats%set('particles', 'particles_extracted', int2commastr(nptcls_glob), primary=.true.)
                !call nice_communicator%update_cls2D(particles_extracted=nptcls_glob)
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "particles_imported", nptcls_glob, found)
                ! update progress monitor
                !call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                l_haschanged = .true.
                ! always write micrographs snapshot if less than 1000 mics, else every 100
                if( n_imported < 1000 .and. l_haschanged )then
                    call update_user_params(cline)
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir)
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
                    call update_user_params(cline)
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir)
                    nmic_star = n_imported
                endif
                ! nice

            endif
            ! 2D section
            if( l_once .and. (nptcls_glob > params%nptcls_per_cls*params%ncls) )then
                if( .not.cline%defined('mskdiam') ) params%mskdiam = 0.85 * real(params%box) * params%smpd
                call init_cluster2D_stream( cline, spproj_glob, micspproj_fname, reference_generation=.true. )
                l_once = .false.
            endif
            call update_pool_status
            call update_pool
            call import_records_into_pool( projrecords )
            call analyze2D_pool
            if(  params%nparts > 1 ) then
                call sleep(WAITTIME)
            else
                call sleep(SHORTWAIT)
            endif
            ! nice
            if(l_once) then
                !nice_communicator%stat_root%stage = "importing and extracting >" // int2str(params%nptcls_per_cls*params%ncls) // " particles"
                call http_communicator%json%update(http_communicator%job_json, "stage", "importing and extracting >" // int2str(params%nptcls_per_cls*params%ncls) // " particles", found)
            else 
               ! call nice_communicator%update_cls2D(iteration=get_pool_iter()-1, number_classes=get_pool_n_classes(), number_classes_rejected=get_pool_n_classes_rejected(),&
               ! number_particles_assigned=get_pool_assigned(), number_particles_rejected=get_pool_rejected(), maximum_resolution=get_pool_res(), &
               ! thumbnail=trim(get_pool_cavgs_jpeg()), thumbnail_id=get_pool_iter(), thumbnail_static_id=1, last_iteration=get_pool_iter_time(), thumbnail_n_tiles=get_pool_cavgs_jpeg_ntiles(), scale=get_pool_cavgs_jpeg_scale())
               ! call nice_communicator%update_cls2D(stats_mask=get_pool_cavgs_mask(), stats_resolution=get_pool_cavgs_res(), stats_population=get_pool_cavgs_pop())
                if(get_pool_iter() > 2) then ! put back to 10!!!!
                !    nice_communicator%stat_root%stage = "iterating whilst waiting for user input"
                   ! if(.not. l_user) then
              !          nice_communicator%stat_root%user_input = .true.
                   !     l_user = .true.
                  ! end if
                  ! http stats
                    call http_communicator%json%update(http_communicator%job_json, "stage",      "iterating whilst waiting for user input", found)
                    call http_communicator%json%update(http_communicator%job_json, "user_input", .true.,                                    found)
                    xtile = 0
                    ytile = 0
                    call http_communicator%json%remove(latest_cls2D, destroy=.true.)
                    call http_communicator%json%create_array(latest_cls2D, "latest_cls2D")
                    call http_communicator%json%add(http_communicator%job_json, latest_cls2D)
                    if(allocated(pool_jpeg_map)) then
                        do i=0, size(pool_jpeg_map) - 1
                            call communicator_add_cls2D(trim(get_pool_cavgs_jpeg()), trim(cwd_glob) // '/' // trim(get_pool_cavgs_mrc()), pool_jpeg_map(i + 1), xtile * (100 / (get_pool_cavgs_jpeg_ntilesx() - 1)), ytile * (100 / (get_pool_cavgs_jpeg_ntilesy() - 1)), 100 * get_pool_cavgs_jpeg_ntilesy(), 100 * get_pool_cavgs_jpeg_ntilesx())
                            xtile = xtile + 1
                            if(xtile .eq. get_pool_cavgs_jpeg_ntilesx()) then
                                xtile = 0
                                ytile = ytile + 1
                            endif
                        end do
                    endif
                else
               !     nice_communicator%stat_root%stage = "performing 10 initial 2D iterations"
                    call http_communicator%json%update(http_communicator%job_json, "stage", "performing 10 initial 2D iterations", found)
                end if
                call http_communicator%json%update(http_communicator%job_json, "particles_processed", get_pool_assigned(), found)
                call http_communicator%json%update(http_communicator%job_json, "particles_rejected",  get_pool_rejected(), found)
                call http_communicator%json%update(http_communicator%job_json, "iteration",           get_pool_iter(),     found)
                call http_communicator%json%update(http_communicator%job_json, "box_size",            get_box(),           found)
            end if
            !call nice_communicator%cycle()
            ! guistats
            !if(file_exists(POOLSTATS_FILE)) then
            !    call gui_stats%merge(POOLSTATS_FILE)
            !end if
            !call gui_stats%write_json
            call http_communicator%send_jobstats()
        end do
        ! termination
        !nice_communicator%stat_root%stage = "terminating"
        !call nice_communicator%cycle()
        call http_communicator%json%update(http_communicator%job_json, "stage", "terminating", found)
        call http_communicator%json%update(http_communicator%job_json, "user_input", .false., found)
        call http_communicator%json%remove(latest_cls2D, destroy=.true.)
        if( l_once )then
            ! nothing to write
        else
            call terminate_stream2D(projrecords)
        endif
        if(file_exists(STREAM_REJECT_CLS)) call write_pool_cls_selected_user
       ! if(associated(nice_communicator%update_arguments)) then
       !     call json%get(nice_communicator%update_arguments, 'final_selection', final_selection, json_found)
       !     call write_pool_cls_selected_nice(final_selection)
       !     call nice_communicator%update_cls2D(thumbnail=trim(get_pool_cavgs_jpeg()), thumbnail_id=get_pool_iter(), thumbnail_static_id=1, thumbnail_n_tiles=get_pool_cavgs_jpeg_ntiles())
       !     call nice_communicator%cycle()
       !     if(allocated(final_selection)) deallocate(final_selection)
       ! end if
        call json%get(http_communicator%update_arguments, 'final_selection', final_selection, found)
        if(found) then
            call json%get(http_communicator%update_arguments, 'final_selection_source', final_selection_source, found)
            if(found) then
                call write_selected_references(final_selection_source, final_selection, nxtiles, nytiles)
                xtile = 0
                ytile = 0
                do i=0, size(final_selection) - 1
                    call communicator_add_selected_reference(trim(cwd_glob) // '/' // STREAM_SELECTED_REFS // JPG_EXT, xtile * (100 / (nxtiles - 1)), ytile * (100 / (nytiles - 1)), 100 * nytiles, 100 * nxtiles)
                    xtile = xtile + 1
                    if(xtile .eq. nxtiles) then
                        xtile = 0
                        ytile = ytile + 1
                    endif
                end do
            endif
        endif
        if(allocated(final_selection)) deallocate(final_selection)
        call json%get(http_communicator%update_arguments, 'final_selection_boxsize', final_selection_boxsize, found)
        if(found) then
            call http_communicator%json%update(http_communicator%job_json, "selected_boxsize", final_selection_boxsize, found)
            allocate(final_boxsize(1))
            final_boxsize(1) = final_selection_boxsize
            call boxsize_file%new(trim(cwd_glob) // '/' // STREAM_SELECTED_REFS // TXT_EXT, 2, 1)
            call boxsize_file%write(final_boxsize)
            call boxsize_file%kill()
            deallocate(final_boxsize)
        endif
      !  call gui_stats%delete('latest', '')
      !  call gui_stats%set('selected references', '', trim(cwd_glob) // '/' // STREAM_SELECTED_REFS // trim(JPG_EXT), thumbnail = .true., smpd = params%smpd, box = real(params%box))
        !call write_project
        !call update_user_params(cline)
        ! call starproj_stream%copy_micrographs_optics(spproj_glob, write=.true., verbose=DEBUG_HERE)
        ! call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=.true.)
        ! call starproj_stream%stream_export_particles_2D(spproj_glob, params%outdir, optics_set=.true., verbose=.true.)
        ! final stats
      !  call gui_stats%hide('compute', 'compute_in_use')
      !  call gui_stats%deactivate_section('compute')
      !  call gui_stats%write_json
      !  call gui_stats%kill
        ! cleanup
        call kill_projrecords(projrecords)
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
       ! call nice_communicator%terminate()
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_GENERATE_PICKING_REFERENCES NORMAL STOP ****')

        contains

            ! updates global project, returns records of processed micrographs
            subroutine update_records_with_project( records, nimported )
                type(projrecord), allocatable, intent(inout) :: records(:)
                integer,                       intent(inout) :: nimported
                type(sp_project), allocatable :: spprojs(:)
                type(projrecord), allocatable :: old_records(:)
                character(len=:), allocatable :: fname, abs_fname, new_fname
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
                        allocate(projrecords(nmics))
                    else
                        call spproj_glob%os_mic%reallocate(n_completed)
                        call move_alloc(projrecords, old_records)
                        allocate(projrecords(n_completed))
                        if( n_old > 0 ) projrecords(1:n_old) = old_records(:)
                    call kill_projrecords(old_records)
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
                            projrecords(j)%projname = trim(abs_fname)
                            projrecords(j)%micind   = imic
                            nptcls      = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                            nptcls_glob = nptcls_glob + nptcls ! global update
                            projrecords(j)%nptcls     = nptcls
                            projrecords(j)%nptcls_sel = nptcls
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
                        if( tmp_proj%os_mic%isthere(imic, 'icefrac') )then
                            if( tmp_proj%os_mic%get(imic,'icefrac') > (params%icefracthreshold-0.001) ) states(imic) = 0
                        end if
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
                        if(params%moldiam == 0.0) then
                            mask = .true.
                        else if(params%box == 0) then
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
                call spproj_here%projinfo%set(1,'projname', projname)
                call spproj_here%projinfo%set(1,'projfile', projfile)
                call spproj_here%projinfo%set(1,'cwd',      path)
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

            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",                   "initialising")
                call http_communicator%json%add(http_communicator%job_json, "particles_imported",      0)
                call http_communicator%json%add(http_communicator%job_json, "particles_processed",     0)
                call http_communicator%json%add(http_communicator%job_json, "particles_rejected",      0)
                call http_communicator%json%add(http_communicator%job_json, "iteration",               0)
                call http_communicator%json%add(http_communicator%job_json, "box_size",                0)
                call http_communicator%json%add(http_communicator%job_json, "selected_boxsize",        0)
                call http_communicator%json%add(http_communicator%job_json, "last_particles_imported", "")
                call http_communicator%json%add(http_communicator%job_json, "user_input",              .false.)
                call http_communicator%json%create_array(latest_cls2D, "latest_cls2D")
                call http_communicator%json%add(http_communicator%job_json, latest_cls2D)
                call http_communicator%json%create_array(selected_references, "selected_references")
                call http_communicator%json%add(http_communicator%job_json, selected_references)
            end subroutine communicator_init

            subroutine communicator_add_cls2D(path, mrcpath, mrc_idx, spritex, spritey, spriteh, spritew)
                character(*),     intent(in)  :: path, mrcpath
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew, mrc_idx
                type(json_value), pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "mrcpath", mrcpath)
                call http_communicator%json%add(template, "mrcidx",  mrc_idx)
                call http_communicator%json%add(template, "spritex", spritex)
                call http_communicator%json%add(template, "spritey", spritey)
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(latest_cls2D,        template)
            end subroutine communicator_add_cls2D

            subroutine communicator_add_selected_reference(path, spritex, spritey, spriteh, spritew)
                character(*),     intent(in)  :: path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew
                type(json_value), pointer     :: template
                call http_communicator%json%create_object(template, "")
                call http_communicator%json%add(template, "path",    path)
                call http_communicator%json%add(template, "spritex", spritex)
                call http_communicator%json%add(template, "spritey", spritey)
                call http_communicator%json%add(template, "spriteh", spriteh)
                call http_communicator%json%add(template, "spritew", spritew)
                call http_communicator%json%add(selected_references, template)
            end subroutine communicator_add_selected_reference

    end subroutine exec_stream_gen_picking_refs

    subroutine exec_stream_assign_optics( self, cline )
        class(commander_stream_assign_optics), intent(inout) :: self
        class(cmdline),                        intent(inout) :: cline
        type(parameters)                       :: params
        type(stream_http_communicator)         :: http_communicator
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj, spproj_part
        type(starproject_stream)               :: starproj_stream
        type(json_value),          pointer     :: optics_assignments, optics_group     
        type(json_value),          pointer     :: optics_group_coordinates,  optics_group_coordinate        
        character(len=LONGSTRLEN), allocatable :: projects(:)
        integer                                :: nprojects, iproj, iori, new_oris, nimported, i, j
        logical                                :: found
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
        ! generate own project file if projfile isnt set
        if(cline%get_carg('projfile') .eq. '') then 
            call cline%set('projname', 'assign_optics')
            call cline%set('projfile', 'assign_optics.simple')
            call spproj%update_projinfo(cline)
            call spproj%update_compenv(cline)
            call spproj%write
        endif
        ! master parameters
        call params%new(cline)
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver, "optics_assignment")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! master project file
        call spproj%read( params%projfile )
        if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(0, .false.)
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
        endif
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! initialise progress monitor
        call progressfile_init()
        ! Infinite loop
        nimported = 0
        do
            if( file_exists(trim(TERM_STREAM)) .or. http_communicator%exit) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj%kill
                call qsys_cleanup
             !   call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS USER STOP ****')
                call EXIT(0)
            endif
            ! http stats
            call http_communicator%json%update(http_communicator%job_json, "stage", "finding and processing new micrographs", found) 
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
                ! http stats
                call http_communicator%json%update(http_communicator%job_json, "micrographs_assigned",     spproj%os_mic%get_noris(),    found)
                call http_communicator%json%update(http_communicator%job_json, "optics_groups_assigned",   spproj%os_optics%get_noris(), found)
                call http_communicator%json%update(http_communicator%job_json, "last_micrograph_imported", stream_datestr(),             found)
                call http_communicator%json%remove(optics_assignments, destroy=.true.)
                call http_communicator%json%create_array(optics_assignments, "optics_assignments")
                call http_communicator%json%add(http_communicator%job_json, optics_assignments)
                do i = 1, spproj%os_optics%get_noris()
                    call http_communicator%json%create_array(optics_group_coordinates, "coordinates")
                    do j = 1, spproj%os_mic%get_noris()
                        if (spproj%os_mic%get(j, 'ogid') .eq. real(i)) then
                            call http_communicator%json%create_object(optics_group_coordinate, "")
                            call http_communicator%json%add(optics_group_coordinate, "x", dble(spproj%os_mic%get(i, 'shiftx')))
                            call http_communicator%json%add(optics_group_coordinate, "y", dble(spproj%os_mic%get(i, 'shifty')))
                            call http_communicator%json%add(optics_group_coordinates, optics_group_coordinate)
                        endif
                    enddo
                    call http_communicator%json%create_object(optics_group, "")
                    call http_communicator%json%add(optics_group, optics_group_coordinates)
                    call http_communicator%json%add(optics_group, "id", i)
                    call http_communicator%json%add(optics_assignments, optics_group)
                enddo
            else
                call sleep(WAITTIME) ! may want to increase as 3s default
            endif
            call update_user_params(cline)
            if(params%updated .eq. 'yes') then
                call starproj_stream%stream_export_optics(spproj, params%outdir)
                params%updated = 'no'
            end if
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        call http_communicator%json%update(http_communicator%job_json, "stage", "terminating", found) 
        call http_communicator%send_jobstats()
        if(allocated(projects)) deallocate(projects)
        ! cleanup
        call spproj%write
        call spproj%kill
        ! end gracefully
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS NORMAL STOP ****')

        contains
        
            subroutine communicator_init()
                call http_communicator%json%add(http_communicator%job_json, "stage",                    "initialising")
                call http_communicator%json%add(http_communicator%job_json, "micrographs_assigned ",    0)
                call http_communicator%json%add(http_communicator%job_json, "optics_groups_assigned",   0)
                call http_communicator%json%add(http_communicator%job_json, "last_micrograph_imported", "")
                call http_communicator%json%create_array(optics_assignments, "optics_assignments")
                call http_communicator%json%add(http_communicator%job_json, optics_assignments)
            end subroutine communicator_init

    end subroutine exec_stream_assign_optics

    subroutine exec_stream_cluster2D( self, cline )
        class(commander_stream_cluster2D), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        character(len=STDLEN),     parameter   :: micspproj_fname = './streamdata.simple'
        integer,                   parameter   :: PAUSE_NITERS    = 5   ! # of iterations after which 2D analysis is paused
        integer,                   parameter   :: PAUSE_TIMELIMIT = 600 ! time (secs) after which 2D analysis is paused
        integer(kind=dp),          parameter   :: FLUSH_TIMELIMIT = 900 ! time (secs) after which leftover particles join the pool IF the 2D analysis is paused
        type(projrecord),          allocatable :: projrecords(:)
        type(parameters)                       :: params
        type(simple_nice_communicator)         :: nice_communicator
        type(guistats)                         :: gui_stats
        type(oris)                             :: moldiamori, chunksizeori
        type(moviewatcher)                     :: project_buff
        type(sp_project)                       :: spproj_glob
        type(json_core)                        :: json
        type(json_value),          pointer     :: snapshot_json
        character(len=LONGSTRLEN), allocatable :: projects(:)
        character(len=LONGSTRLEN)              :: cwd_job
        integer(kind=dp)                       :: pool_time_last_chunk, time_last_import
        integer                                :: nmics_rejected_glob, nchunks_imported_glob, nchunks_imported, nprojects, iter
        integer                                :: n_imported, n_imported_prev, n_added, nptcls_glob, n_failed_jobs, ncls_in, pool_iter
        integer                                :: pool_iter_last_chunk, pool_iter_max_chunk_imported, i
        integer,                   allocatable :: snapshot_selection(:)
        logical                                :: l_nchunks_maxed, l_pause, l_params_updated
        real                                   :: nptcls_pool, moldiam
        nullify(snapshot_json)
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
        if( .not. cline%defined('ndev')         ) call cline%set('ndev',         CLS_REJECT_STD)
        if( .not. cline%defined('reject_cls')   ) call cline%set('reject_cls',   'yes')
        if( .not. cline%defined('objfun')       ) call cline%set('objfun',       'euclid')
        if( .not. cline%defined('ml_reg')       ) call cline%set('ml_reg',       'no')
        if( .not. cline%defined('tau')          ) call cline%set('tau',          5)
        if( .not. cline%defined('remove_chunks')) call cline%set('remove_chunks','yes')
        if( .not. cline%defined('refine')       ) call cline%set('refine',       'snhc_smpl')
        if( .not. cline%defined('dynreslim')    ) call cline%set('dynreslim',    'no')
        if( .not. cline%defined('cavg_ini')     ) call cline%set('cavg_ini',     'no')
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
            ncls_in = cline%get_iarg('ncls')
            call cline%delete('ncls')
        endif
        ! generate own project file
        call cline%set('projname', 'cluster2D')
        call cline%set('projfile', 'cluster2D.simple')
        call spproj_glob%update_projinfo(cline)
        call spproj_glob%update_compenv(cline)
        call spproj_glob%write 
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
        ! nice communicator init
        call nice_communicator%init(params%niceprocid, params%niceserver)
        call nice_communicator%cycle()
        ! restart
        if( cline%defined('dir_exec') )then
            call cline%delete('dir_exec')
            call del_file(micspproj_fname)
            call cleanup_root_folder
        endif
        ! needed for stream3d
        call simple_mkdir(DIR_SNAPSHOT)
        ! mskdiam
        if( .not. cline%defined('mskdiam') ) then
            ! nice communicator status
            nice_communicator%stat_root%stage = "waiting for mask diameter"
            call nice_communicator%cycle()
            write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 60 MINUTES FOR '//trim(STREAM_MOLDIAM)
            do i=1, 360
                if(file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) exit
                call sleep(10)
            end do
            if( .not. file_exists(trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM))) THROW_HARD('either mskdiam must be given or '// trim(STREAM_MOLDIAM) // ' exists in target_dir')
            ! read mskdiam from file
            call moldiamori%new(1, .false.)
            call moldiamori%read( trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            if( .not. moldiamori%isthere(1, "moldiam") ) THROW_HARD( 'moldiam missing from ' // trim(params%dir_target)//'/'//trim(STREAM_MOLDIAM) )
            moldiam = moldiamori%get(1, "moldiam")
           
            ! write acopy for stream 3d
            call moldiamori%write(1, trim(STREAM_MOLDIAM))
            call moldiamori%kill
            params%mskdiam = moldiam * 1.2
            call cline%set('mskdiam', params%mskdiam)
            write(logfhandle,'(A,F8.2)')'>>> MASK DIAMETER SET TO', params%mskdiam
        endif
        ! Resolution based class rejection
        if( .not.cline%defined('lpthres') )then
            call mskdiam2streamresthreshold(params%mskdiam, params%lpthres)
            call cline%set('lpthres', params%lpthres)
            call set_lpthres_type("auto")
        else if(params%lpthres .gt. LOWRES_REJECT_THRESHOLD) then
            call set_lpthres_type("off")
        else
            call set_lpthres_type("manual")
        endif
        ! Number of particles per class
        if(params%nptcls_per_cls == 0) write(logfhandle,'(A)')   '>>> # PARTICLES PER CLASS WILL BE AUTO DETERMINED AFTER 100 IMPORTED MICROGRAPHS'
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! wait if dir_target doesn't exist yet
        if(.not. dir_exists(trim(params%dir_target))) then
            write(logfhandle, *) ">>> WAITING FOR ", trim(params%dir_target), " TO BE GENERATED"
            do i=1, 360
                if(dir_exists(trim(params%dir_target))) then
                    write(logfhandle, *) ">>> ", trim(params%dir_target), " FOUND"
                    exit
                endif
                call sleep(10)
            end do
        endif
        ! movie watcher init
        project_buff = moviewatcher(LONGTIME, trim(params%dir_target)//'/'//trim(DIR_STREAM_COMPLETED), spproj=.true., nretries=10)
        ! Infinite loop
        nptcls_glob           = 0       ! global number of particles
        nchunks_imported_glob = 0       ! global number of completed chunks
        n_imported            = 0       ! global number of imported processed micrographs
        n_imported_prev       = 0
        nmics_rejected_glob   = 0       ! global number of micrographs rejected
        nprojects             = 0
        iter                  = 0       ! global number of infinite loop iterations
        n_failed_jobs         = 0
        l_nchunks_maxed       = .false. ! Whether a minimum number of chunks as been met
        l_pause               = .false. ! whether pool 2D analysis is skipped
        pool_iter_last_chunk  = -1                          ! used for pausing 2D analysis
        pool_time_last_chunk  = huge(pool_time_last_chunk)  ! used for pausing 2D analysis
        time_last_import      = huge(time_last_import)      ! used for flushing unprocessed particles
        ! guistats init
        call gui_stats%init(.true.)
        call gui_stats%set('particles', 'particles_imported',          0,            primary=.true.)
        call gui_stats%set('2D',        'iteration',                   0,            primary=.true.)
        ! Joe: nparts is not an input, also see project_buff below
        call gui_stats%set('compute',   'compute_in_use',      int2str(0) // '/' // int2str(params%nparts), primary=.true.)
        ! nice
        nice_communicator%stat_root%stage = "importing particles"
        call nice_communicator%update_cls2D(particles_imported=0)
        call nice_communicator%cycle()
        do
            if( file_exists(trim(TERM_STREAM)) .or. nice_communicator%exit)then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( nice_communicator%stop .or. test_repick() ) then
                if(test_repick()) call write_repick_refs("../repick_refs.mrc")
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
                call nice_communicator%terminate(stop=.true.)
                call simple_end('**** SIMPLE_STREAM_CLUSTER2D USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
            ! detection of new projects
            if( l_nchunks_maxed )then
                call project_buff%kill
                nprojects = 0
            else
                call project_buff%watch(nprojects, projects, max_nmovies=10*params%nparts)
            endif
            ! update global records
            if( nprojects > 0 )then
                call update_records_with_project(projects, n_imported )
                call project_buff%add2history(projects)
            endif
            ! project update
            if( nprojects > 0 )then
                n_imported = size(projrecords)
                write(logfhandle,'(A,I6,I8)') '>>> # MICROGRAPHS / PARTICLES IMPORTED : ',n_imported, nptcls_glob
                ! guistats
                call gui_stats%set('particles', 'particles_imported', int2commastr(nptcls_glob), primary=.true.)
                call gui_stats%set_now('particles', 'last_particles_imported')
                call nice_communicator%update_cls2D(particles_imported=nptcls_glob, last_particles_imported=.true.)

                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                time_last_import = time8()
                if( n_imported < 1000 )then
                    call update_user_params(cline)
                else if( n_imported > n_imported_prev + 100 )then
                    call update_user_params(cline)
                    n_imported_prev = n_imported
                endif
            endif
            ! 2D analysis section
            nice_communicator%view_cls2D%mskdiam  = params%mskdiam
            nice_communicator%view_cls2D%boxsizea = get_boxa()
            call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
            if( l_params_updated ) l_pause = .false.
            call update_chunks
            if( l_pause )then
                call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
                if( l_params_updated ) l_pause = .false.    ! resuming 2D analysis
            else
                call update_pool_status
                call update_pool
                call update_user_params2D(cline, l_params_updated, nice_communicator%update_arguments)
                call reject_from_pool
            endif
            call reject_from_pool_user
            if( l_nchunks_maxed )then
                ! # of chunks is above desired number
                if( is_pool_available() .and. (get_pool_iter()>=pool_iter_max_chunk_imported+10) ) exit
                call analyze2D_pool
            else
                call import_chunks_into_pool( nchunks_imported )
                if( nchunks_imported > 0 )then
                    nchunks_imported_glob = nchunks_imported_glob + nchunks_imported
                    pool_iter_last_chunk  = get_pool_iter()
                    pool_time_last_chunk  = time8()
                    l_pause = .false.   ! resuming 2D analysis
                endif
                if( nchunks_imported_glob >= params%maxnchunks )then
                    if( .not.l_nchunks_maxed ) pool_iter_max_chunk_imported = get_pool_iter()
                    l_nchunks_maxed = .true.
                endif
                if( l_pause )then
                    ! Whether to flush particles
                    if( (time8()-time_last_import > FLUSH_TIMELIMIT) .and. all_chunks_available() )then
                        ! Remaining unclassified particles will join the pool directly if already paused,
                        ! and all chunks are inactive and no new particles are being imported
                        call flush_remaining_particles(projrecords)
                    endif
                else
                    ! whether to pause
                    if( pool_iter_last_chunk > 0 )then
                        if( (get_pool_iter() >= pool_iter_last_chunk+PAUSE_NITERS) .or.&
                            &(time8()-pool_time_last_chunk > PAUSE_TIMELIMIT) )then
                            ! pause pool 2D analysis & rejection in absence of new chunks, resumes
                            ! when new chunks are added or 2D analysis parameters have been updated
                            l_pause = is_pool_available()
                            if( l_pause ) write(logfhandle,'(A)')'>>> PAUSING 2D ANALYSIS'
                        endif
                    endif
                endif
                ! 2D analyses
                if( l_pause )then
                    call generate_pool_stats
                else
                    pool_iter = get_pool_iter()
                    call analyze2D_pool
                    if(get_pool_iter() > pool_iter) nice_communicator%stat_root%user_input = .true.
                endif
                call analyze2D_new_chunks(projrecords)
                ! Optionally generates snapshot project for abinitio3D
                call generate_snapshot_for_abinitio
            endif
            call sleep(WAITTIME)
            ! guistats
            if(file_exists(POOLSTATS_FILE)) then
                call gui_stats%merge(POOLSTATS_FILE)
                call gui_stats%get('particles', 'particles_processed', nptcls_pool)
                if(nptcls_pool > 0.0) call gui_stats%set('particles', 'particles_processed', int2commastr(floor(nptcls_pool)) // ' (' // int2str(ceiling(100.0 * real(nptcls_pool) / real(nptcls_glob))) // '%)')
            end if
            call gui_stats%write_json
            ! nice
            if(l_pause) then
                nice_communicator%stat_root%stage = "paused pool 2D analysis"
                nice_communicator%stat_root%user_input = .true.
            else if(get_nchunks() > 0) then
                if(get_pool_iter() < 1) then
                    nice_communicator%stat_root%stage = "classifying initial chunks"
                    nice_communicator%stat_root%user_input = .false.
                else
                    nice_communicator%stat_root%stage = "classifying pool and new chunks"
                    nice_communicator%stat_root%user_input = .true.
                end if
            end if
            call nice_communicator%update_cls2D(iteration=get_pool_iter()-1 , number_classes=get_pool_n_classes(), number_classes_rejected=get_pool_n_classes_rejected(),&
            number_particles_assigned=get_pool_assigned(), number_particles_rejected=get_pool_rejected(), maximum_resolution=get_pool_res(), &
            thumbnail=trim(get_pool_cavgs_jpeg()), thumbnail_id=get_pool_iter(), thumbnail_static_id=1, last_iteration=get_pool_iter_time(), thumbnail_n_tiles=get_pool_cavgs_jpeg_ntiles())
            if(get_pool_rejected_jpeg_ntiles() .gt. 0) then
                call nice_communicator%update_cls2D(pool_rejected_thumbnail=trim(get_pool_rejected_jpeg()), pool_rejected_thumbnail_id=1000 + get_pool_rejected_thumbnail_id(), pool_rejected_thumbnail_static_id=2, &
                pool_rejected_thumbnail_n_tiles=get_pool_rejected_jpeg_ntiles(), pool_rejected_scale=get_pool_rejected_jpeg_scale())
            end if
            if(get_chunk_rejected_jpeg_ntiles() .gt. 0) then
                call nice_communicator%update_cls2D(chunk_rejected_thumbnail=trim(get_chunk_rejected_jpeg()), chunk_rejected_thumbnail_id=2000 + get_chunk_rejected_thumbnail_id(), chunk_rejected_thumbnail_static_id=3, &
                chunk_rejected_thumbnail_n_tiles=get_chunk_rejected_jpeg_ntiles(), chunk_rejected_scale=get_chunk_rejected_jpeg_scale())
            end if
            call nice_communicator%update_cls2D(stats_mask=get_pool_cavgs_mask(), stats_resolution=get_pool_cavgs_res(), stats_population=get_pool_cavgs_pop(), rejection_params=get_rejection_params())
            call nice_communicator%update_cls2D(snapshot_id=get_last_snapshot_id(), snapshot_time=get_last_snapshot())
            call nice_communicator%view_cls2D%res_histogram%zero()
            do i = 1, get_pool_n_classes()
                call nice_communicator%view_cls2D%res_histogram%update(get_pool_cavgs_res_at(i))
            enddo
            nice_communicator%view_cls2D%cutoff_res = params%lpthres
            nice_communicator%view_cls2D%cutoff_type = get_lpthres_type()
            ! project snapshot
            if(get_pool_iter() .gt. 0 .and. trim(params%snapshot) .ne. "") then
                call nice_communicator%update_cls2D(snapshot_json_clear=.true.)
                call write_project_stream2D(snapshot_projfile=trim(params%snapshot))
                params%snapshot = ""
                params%updated = "no"
                if(associated(snapshot_json)) then
                    call json%destroy(snapshot_json)
                    nullify(snapshot_json)
                end if
                call get_snapshot_json(snapshot_json)
                call nice_communicator%update_cls2D(snapshot_id=get_last_snapshot_id(), snapshot_time=get_last_snapshot(), snapshot_json=snapshot_json)
            else
                call nice_communicator%update_cls2D(snapshot_id=get_last_snapshot_id(), snapshot_time=get_last_snapshot())
            end if
            call nice_communicator%cycle()
        end do
        ! termination
        nice_communicator%stat_root%stage = "terminating"
        call nice_communicator%cycle()
        call terminate_stream2D( projrecords )
        call update_user_params(cline) ! Joe: bit late for this?
        ! final stats
        if(file_exists(POOLSTATS_FILE)) call gui_stats%merge(POOLSTATS_FILE, delete = .true.)
        call gui_stats%hide('compute', 'compute_in_use')
        call gui_stats%deactivate_section('compute')
        call gui_stats%write_json
        call gui_stats%kill
        ! cleanup
        call kill_projrecords(projrecords)
        call spproj_glob%kill
        call qsys_cleanup
        ! end gracefully
        call nice_communicator%terminate()
        call simple_end('**** SIMPLE_STREAM_CLUSTER2D NORMAL STOP ****')
        contains

            ! updates global records
            subroutine update_records_with_project( projectnames, n_imported )
                character(len=LONGSTRLEN), allocatable, intent(in)  :: projectnames(:)
                integer,                                intent(out) :: n_imported
                type(sp_project),     allocatable :: spprojs(:)
                type(projrecord),     allocatable :: old_records(:)
                character(len=:),     allocatable :: fname, abs_fname
                real    :: avgmicptcls, nptcls_per_cls
                integer :: iproj, n_spprojs, n_old, irec, n_completed, nptcls, nmics, imic, n_ptcls, first
                n_imported = 0
                n_ptcls    = 0
                if( .not.allocated(projectnames) ) return
                n_spprojs  = size(projectnames)
                if( n_spprojs == 0 )return
                n_old = 0 ! on first import
                if( allocated(projrecords) ) n_old = size(projrecords)
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
                ! import micrographs
                n_completed = n_old + nmics
                n_imported  = nmics
                ! reallocate records
                if( n_old == 0 )then
                    allocate(projrecords(nmics))
                else
                    call move_alloc(projrecords, old_records)
                    allocate(projrecords(n_completed))
                    projrecords(1:n_old) = old_records(:)
                    deallocate(old_records)
                endif
                ! update global records and some global variables
                irec = n_old
                do iproj = 1,n_spprojs
                    do imic = 1,spprojs(iproj)%os_mic%get_noris()
                        irec      = irec + 1
                        nptcls    = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                        n_ptcls   = n_ptcls + nptcls ! global update
                        fname     = trim(projectnames(iproj))
                        abs_fname = simple_abspath(fname, errmsg='stream_cluster2D :: update_projects_list 1')
                        projrecords(irec)%projname   = trim(abs_fname)
                        projrecords(irec)%micind     = imic
                        projrecords(irec)%nptcls     = nptcls
                        projrecords(irec)%nptcls_sel = nptcls
                        projrecords(irec)%included   = .false.
                    enddo
                enddo
                nptcls_glob = nptcls_glob + n_ptcls ! global update
                ! Updates global parameters once and init 2D
                if(params%nptcls_per_cls == 0) then
                    if(size(projrecords) .gt. 100) then
                        avgmicptcls = nptcls_glob / size(projrecords)
                        avgmicptcls = ceiling(avgmicptcls / 10) * 10.0
                        ! these parameters may need tweaking
                        nptcls_per_cls = 1000 * (20 + (0.15 * avgmicptcls))
                        nptcls_per_cls = nptcls_per_cls / real(params%ncls_start)
                        nptcls_per_cls = ceiling(nptcls_per_cls / 100) * 100.0
                        write(logfhandle,'(A,I6)')   '>>> AVERAGE # PARTICLES PER MICROGRAPH : ', int(avgmicptcls)
                        write(logfhandle,'(A,I6,A)') '>>> USING ', int(nptcls_per_cls), ' PARTICLES PER CLASS'
                        params%nptcls_per_cls = int(nptcls_per_cls)
                        call cline%set('nptcls_per_cls', nptcls_per_cls)
                        params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                        call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                        params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                        call init_cluster2D_stream(cline, spproj_glob, micspproj_fname)
                        call cline%delete('ncls')
                        ! write out for stream3d to pick up
                        call chunksizeori%new(1, .false.)
                        call chunksizeori%set(1, 'nptcls_per_cls', params%nptcls_per_cls)
                        call chunksizeori%write(1, trim(STREAM_CHUNKSIZE))
                        call chunksizeori%kill
                    end if
                else if( n_old == 0 )then
                    params%smpd = spprojs(first)%os_mic%get(1,'smpd')
                    call spprojs(first)%read_segment('stk', trim(projectnames(first)))
                    params%box  = nint(spprojs(first)%os_stk%get(1,'box'))
                    call init_cluster2D_stream(cline, spproj_glob, micspproj_fname)
                    call cline%delete('ncls')
                endif
                ! cleanup
                do iproj = 1,n_spprojs
                    call spprojs(iproj)%kill
                enddo
                deallocate(spprojs)
            end subroutine update_records_with_project

    end subroutine exec_stream_cluster2D

end module simple_commander_stream
