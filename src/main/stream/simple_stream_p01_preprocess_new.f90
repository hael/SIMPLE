!@descr: task 1 in the stream pipeline: pre-processing (movie registration, CTF estimation, segmentation-based picking)
module simple_stream_p01_preprocess_new
use unix,                        only: SIGTERM
use simple_stream_api
use simple_stream_mq_defs
use simple_gui_metadata_api
use simple_motion_correct_utils, only: flip_gain
use simple_histogram,            only: histogram
implicit none

public :: stream_p01_preprocess
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p01_preprocess
  contains
    procedure :: execute => exec_stream_p01_preprocess
end type stream_p01_preprocess

contains

    subroutine exec_stream_p01_preprocess( self, cline )
        class(stream_p01_preprocess), intent(inout) :: self
        class(cmdline),               intent(inout) :: cline
        type(parameters)               :: params
        logical,          parameter    :: DEBUG_HERE      = .false.
        class(cmdline),  allocatable   :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(cmdline)                  :: cline_exec
        type(qsys_env)                 :: qenv
        type(stream_watcher)           :: movie_buff
        type(sp_project)               :: spproj_glob    ! global project
        type(starproject_stream)       :: starproj_stream
        type(string), allocatable      :: movies(:), dir_movies(:)
        type(string)                   :: output_dir, output_dir_ctf_estimate, output_dir_motion_correct, projfile
        type(gui_metadata_stream_update)     :: meta_update
        type(gui_metadata_stream_preprocess) :: meta_preprocess
        type(gui_metadata_micrograph)        :: meta_preprocess_micrograph
        type(gui_metadata_histogram)         :: meta_preprocess_histogram_ctfres, meta_preprocess_histogram_icefrac, meta_preprocess_histogram_astig
        type(gui_metadata_timeplot)          :: meta_preprocess_timeplot_ctfres, meta_preprocess_timeplot_astig, meta_preprocess_timeplot_df, meta_preprocess_timeplot_rate
        character(len=:),        allocatable :: meta_buffer
        type(histogram)                      :: key_histogram
        real,                    allocatable :: histogram_rvec(:), histogram_rvec2(:), histogram_rvec3(:)
        integer,                 allocatable :: histogram_ivec(:)
        integer                              :: ikey, ibin, nbins, binsize, iori
        real :: average_ctfres, average_astig, average_icefrac, ave, sdev, var
        logical                      :: err
        integer                      :: fromto(2)

        character(len=STDLEN)          :: preproc_nthr_env, preproc_part_env, preproc_nparts_env
        real        :: stat_dfx_threshold, stat_dfy_threshold
        real        :: stat_astig_threshold, stat_icefrac_threshold, stat_ctfres_threshold
        integer     :: movies_set_counter, import_counter, nwaits, nmovs2importperiter
        integer     :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection, nsets, i, j, i_thumb, i_max
        integer     :: cnt, n_imported, n_added, n_failed_jobs, n_fail_iter, nmic_star, iset, envlen
        logical     :: l_movies_left, l_haschanged, l_restart, SJ_directory_structure, l_dir_found, l_terminate=.false.
        call signal(SIGTERM, sigterm_handler)
        call cline%set('oritype',     'mic')
        call cline%set('mkdir',       'yes')
        call cline%set('reject_mics', 'no')
        if( .not. cline%defined('walltime')         ) call cline%set('walltime',   29.0*60.0) ! 29 minutes
        ! motion correction
        if( .not. cline%defined('trs')              ) call cline%set('trs',              20.)
        if( .not. cline%defined('lpstart')          ) call cline%set('lpstart',          8.)
        if( .not. cline%defined('lpstop')           ) call cline%set('lpstop',           5.)
        if( .not. cline%defined('bfac')             ) call cline%set('bfac',             50.)
        if( .not. cline%defined('mcconvention')     ) call cline%set('mcconvention',     'simple')
        if( .not. cline%defined('eer_upsampling')   ) call cline%set('eer_upsampling',   1)
        if( .not. cline%defined('algorithm')        ) call cline%set('algorithm',        'patch')
        if( .not. cline%defined('mcpatch')          ) call cline%set('mcpatch',          'yes')
        if( .not. cline%defined('mcpatch_thres')    ) call cline%set('mcpatch_thres',    'yes')
        if( .not. cline%defined('tilt_thres')       ) call cline%set('tilt_thres',       0.05)
        if( .not. cline%defined('beamtilt')         ) call cline%set('beamtilt',         'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512)
        if( .not. cline%defined('hp_ctf_estimate')  ) call cline%set('hp_ctf_estimate',  HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate')  ) call cline%set('lp_ctf_estimate',  LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',            DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',            DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',        'yes')
        if( .not. cline%defined('ctfresthreshold') )  call cline%set('ctfresthreshold',  STREAM_CTFRES_THRESHOLD)
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NTHR, preproc_nthr_env, envlen)
        if( envlen > 0)  call cline%set('nthr', str2int(preproc_nthr_env))
        call get_environment_variable(SIMPLE_STREAM_PREPROC_NPARTS, preproc_nparts_env, envlen)
        if( envlen > 0 ) call cline%set('nparts', str2int(preproc_nparts_env))
        ! sanity check for restart
        l_restart = .false.
        if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
            l_restart = .true.
        endif
        ! generate own project file if it doesnt already exist
        projfile = cline%get_carg('projfile')
        if( .not.file_exists(projfile) )then
            call cline%set('projfile', projfile)
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
        ! initialise metadata
        call meta_update%new(                                            GUI_METADATA_STREAM_UPDATE_TYPE)
        call meta_preprocess%new(                                    GUI_METADATA_STREAM_PREPROCESS_TYPE)
        call meta_preprocess_micrograph%new(              GUI_METADATA_STREAM_PREPROCESS_MICROGRAPH_TYPE)
        call meta_preprocess_histogram_ctfres%new(  GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_CTFRES_TYPE)
        call meta_preprocess_histogram_icefrac%new(GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ICEFRAC_TYPE)
        call meta_preprocess_histogram_astig%new(    GUI_METADATA_STREAM_PREPROCESS_HISTOGRAM_ASTIG_TYPE)
        call meta_preprocess_timeplot_ctfres%new(    GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_CTFRES_TYPE)
        call meta_preprocess_timeplot_astig%new(      GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_ASTIG_TYPE)
        call meta_preprocess_timeplot_df%new(            GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_DF_TYPE)
        call meta_preprocess_timeplot_rate%new(        GUI_METADATA_STREAM_PREPROCESS_TIMEPLOT_RATE_TYPE)
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        ! gain reference
        call flip_gain(cline, params%gainref, params%flipgain)
        ! Sniffing for movies & subfolders: we have to wait for first movie/subfolder
        ! to know which directory structure we are working with
        l_dir_found = .false.
        nwaits      = 0
        do while( .not.l_dir_found )
            call workout_directory_structure(params%dir_movies, l_dir_found, SJ_directory_structure)
            if( .not.l_dir_found)then
                call sleep(SHORTWAIT)
                nwaits = nwaits + 1
                if( mod(nwaits*SHORTWAIT, 60)==0 )then
                    write(*,"(A,I3,A)")'>>> NO MOVIE HAS BEEN DETECTED FOR ',&
                    &nint(real(nwaits*SHORTWAIT)/60.),' MINS'
                endif
            endif
        enddo
        ! movie watcher init
        if( SJ_directory_structure )then
            ! add a suffix & multiple folders
            call sniff_folders_SJ(params%dir_movies, l_dir_found, dir_movies )
            if( .not.l_dir_found ) THROW_HARD('Fatal error directory structure')
            movie_buff = stream_watcher(LONGTIME, dir_movies(1), suffix_filter=string('_fractions'))
            call movie_buff%detect_and_add_dirs(params%dir_movies, SJ_directory_structure)
            deallocate(dir_movies)
        else
            ! no suffix, one folder
            movie_buff = stream_watcher(LONGTIME, params%dir_movies)
        endif
        ! restart
        movies_set_counter = 0  ! global number of movies set
        import_counter     = 0  ! global import id
        nmic_star          = 0
        if( l_restart )then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
            if( cline%defined('dir_exec') ) call cline%delete('dir_exec')
            call import_previous_projects
            nmic_star = spproj_glob%os_mic%get_noris()
            call write_mic_star_and_field(write_field=.true.)
        endif
        ! output directories
        call simple_mkdir(PATH_HERE//DIR_STREAM_COMPLETED)
        output_dir = PATH_HERE//DIR_STREAM
        call simple_mkdir(output_dir)
        call simple_mkdir(output_dir//STDERROUT_DIR)
        output_dir_ctf_estimate   = filepath(PATH_HERE, DIR_CTF_ESTIMATE)
        output_dir_motion_correct = filepath(PATH_HERE, DIR_MOTION_CORRECT)
        call simple_mkdir(output_dir_ctf_estimate)
        call simple_mkdir(output_dir_motion_correct)
        call cline%set('dir','../')
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PREPROC_PARTITION, preproc_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(params, 1,stream=.true.,qsys_partition=string(trim(preproc_part_env)))
        else
            call qenv%new(params, 1,stream=.true.)
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
        nmovs2importperiter    = 2*params%nparts*STREAM_NMOVS_SET
        stat_dfx_threshold     = 0.
        stat_dfy_threshold     = 0.
        stat_astig_threshold   = 0.
        stat_icefrac_threshold = 0.
        stat_ctfres_threshold  = 0.
        call cline_exec%set('fromp',1)
        call cline_exec%set('top',  STREAM_NMOVS_SET)
        do
            if( file_exists(TERM_STREAM) .or. l_terminate )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            iter = iter + 1
            ! detection of new folders & movies
            call movie_buff%detect_and_add_dirs(params%dir_movies, SJ_directory_structure)
            call movie_buff%watch( nmovies, movies, max_nmovies=nmovs2importperiter )
            ! append movies to processing stack
            if( nmovies >= STREAM_NMOVS_SET )then
                nsets = floor(real(nmovies) / real(STREAM_NMOVS_SET))
                cnt   = 0
                do iset = 1,nsets
                    i = (iset-1)*STREAM_NMOVS_SET+1
                    j = iset*STREAM_NMOVS_SET
                    call create_movies_set_project(movies(i:j))
                    call qenv%qscripts%add_to_streaming( cline_exec )
                    do imovie = i,j
                        call movie_buff%add2history( movies(imovie) )
                        cnt     = cnt     + 1
                        n_added = n_added + 1 ! global number of movie sets
                    enddo
                    if( cnt == min(nmovs2importperiter,nmovies) ) exit
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
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz*STREAM_NMOVS_SET
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
                ! ctfres histogram + timeplot
                if(spproj_glob%os_mic%isthere("ctfres")) then
                    allocate(histogram_rvec(size(CTFRES_BINS)), histogram_ivec(size(CTFRES_BINS)))
                    histogram_rvec = CTFRES_BINS
                    call key_histogram%new(histogram_rvec)
                    call key_histogram%zero()
                    do ikey = 1, spproj_glob%os_mic%get_noris()
                        call key_histogram%update(spproj_glob%os_mic%get(ikey, 'ctfres'))
                    enddo
                    do ikey=1, size(histogram_rvec)
                        histogram_ivec(ikey) = int(key_histogram%get(ikey))
                    end do
                    call meta_preprocess_histogram_ctfres%set(name=string('ctfres'), labels=histogram_rvec, data=histogram_ivec)
                    if( meta_preprocess%assigned() ) then
                        if( mq_stream_master_in%is_active() ) then
                            call meta_preprocess_histogram_ctfres%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                    deallocate(histogram_rvec, histogram_ivec)
                    binsize = 500
                    nbins = ceiling(spproj_glob%os_mic%get_noris() / real(binsize))
                    allocate(histogram_rvec(nbins), histogram_rvec2(nbins), histogram_rvec3(nbins))
                    do ibin=1, nbins
                        fromto(1) = 1 + ((ibin - 1)  * binsize)
                        fromto(2) = ibin * binsize
                        if(fromto(2) .gt. spproj_glob%os_mic%get_noris()) fromto(2) = spproj_glob%os_mic%get_noris()
                        call spproj_glob%os_mic%stats('ctfres', ave, sdev, var, err, fromto)
                        histogram_rvec(ibin)  = ibin
                        histogram_rvec2(ibin) = ave
                        histogram_rvec3(ibin) = sdev
                    enddo
                    call meta_preprocess_timeplot_ctfres%set(name=string('ctfres'), labels=histogram_rvec, data=histogram_rvec2, data2=histogram_rvec3)
                    if( meta_preprocess%assigned() ) then
                        if( mq_stream_master_in%is_active() ) then
                            call meta_preprocess_timeplot_ctfres%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                    deallocate(histogram_rvec, histogram_rvec2, histogram_rvec3)
                endif
                ! icefrac histogram
                if(spproj_glob%os_mic%isthere("icefrac")) then
                    allocate(histogram_rvec(size(ICESCORE_BINS)), histogram_ivec(size(ICESCORE_BINS)))
                    histogram_rvec = ICESCORE_BINS
                    call key_histogram%new(histogram_rvec)
                    call key_histogram%zero()
                    do ikey = 1, spproj_glob%os_mic%get_noris()
                        call key_histogram%update(spproj_glob%os_mic%get(ikey, 'icefrac'))
                    enddo
                    do ikey=1, size(histogram_rvec)
                        histogram_ivec(ikey) = int(key_histogram%get(ikey))
                    end do
                    call meta_preprocess_histogram_icefrac%set(name=string('icefrac'), labels=histogram_rvec, data=histogram_ivec)
                    if( meta_preprocess%assigned() ) then
                        if( mq_stream_master_in%is_active() ) then
                            call meta_preprocess_histogram_icefrac%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                    deallocate(histogram_rvec, histogram_ivec)
                endif
                ! astigmatism histogram + timeplot
                if(spproj_glob%os_mic%isthere("astig")) then
                    allocate(histogram_rvec(size(ASTIG_BINS)), histogram_ivec(size(ASTIG_BINS)))
                    histogram_rvec = ASTIG_BINS
                    call key_histogram%new(histogram_rvec)
                    call key_histogram%zero()
                    do ikey = 1, spproj_glob%os_mic%get_noris()
                        call key_histogram%update(spproj_glob%os_mic%get(ikey, 'astig'))
                    enddo
                    do ikey=1, size(histogram_rvec)
                        histogram_ivec(ikey) = int(key_histogram%get(ikey))
                    end do
                    call meta_preprocess_histogram_astig%set(name=string('astig'), labels=histogram_rvec, data=histogram_ivec)
                    if( meta_preprocess%assigned() ) then
                        if( mq_stream_master_in%is_active() ) then
                            call meta_preprocess_histogram_astig%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                    deallocate(histogram_rvec, histogram_ivec)
                    binsize = 500
                    nbins = ceiling(spproj_glob%os_mic%get_noris() / real(binsize))
                    allocate(histogram_rvec(nbins), histogram_rvec2(nbins), histogram_rvec3(nbins))
                    do ibin=1, nbins
                        fromto(1) = 1 + ((ibin - 1)  * binsize)
                        fromto(2) = ibin * binsize
                        if(fromto(2) .gt. spproj_glob%os_mic%get_noris()) fromto(2) = spproj_glob%os_mic%get_noris()
                        call spproj_glob%os_mic%stats('astig', ave, sdev, var, err, fromto)
                        histogram_rvec(ibin)  = ibin
                        histogram_rvec2(ibin) = ave
                        histogram_rvec3(ibin) = sdev
                    enddo
                    call meta_preprocess_timeplot_astig%set(name=string('astig'), labels=histogram_rvec, data=histogram_rvec2, data2=histogram_rvec3)
                    if( meta_preprocess%assigned() ) then
                        if( mq_stream_master_in%is_active() ) then
                            call meta_preprocess_timeplot_astig%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    endif
                    deallocate(histogram_rvec, histogram_rvec2, histogram_rvec3)
                endif
                ! df timeplot
                binsize = 500
                nbins = ceiling(spproj_glob%os_mic%get_noris() / real(binsize))
                allocate(histogram_rvec(nbins), histogram_rvec2(nbins), histogram_rvec3(nbins))
                do ibin=1, nbins
                    fromto(1) = 1 + ((ibin - 1)  * binsize)
                    fromto(2) = ibin * binsize
                    if(fromto(2) .gt. spproj_glob%os_mic%get_noris()) fromto(2) = spproj_glob%os_mic%get_noris()
                    call spproj_glob%os_mic%stats('df', ave, sdev, var, err, fromto)
                    histogram_rvec(ibin)  = ibin
                    histogram_rvec2(ibin) = ave
                    histogram_rvec3(ibin) = sdev
                enddo
                call meta_preprocess_timeplot_df%set(name=string('df'), labels=histogram_rvec, data=histogram_rvec2, data2=histogram_rvec3)
                if( meta_preprocess%assigned() ) then
                    if( mq_stream_master_in%is_active() ) then
                        call meta_preprocess_timeplot_df%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                endif
                deallocate(histogram_rvec, histogram_rvec2, histogram_rvec3)
                ! rate timeplot
                nbins = size(movie_buff%ratehistory)
                allocate(histogram_rvec(nbins), histogram_rvec2(nbins))
                do ibin=1, nbins
                    histogram_rvec(ibin)  = ibin
                    histogram_rvec2(ibin) = movie_buff%ratehistory(ibin)
                enddo
                call meta_preprocess_timeplot_rate%set(name=string('rate'), labels=histogram_rvec, data=histogram_rvec2, data2=histogram_rvec3)
                if( meta_preprocess%assigned() ) then
                    if( mq_stream_master_in%is_active() ) then
                        call meta_preprocess_timeplot_rate%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                endif
                deallocate(histogram_rvec, histogram_rvec2)

                if(spproj_glob%os_mic%isthere('thumb')) then
                    i_max = min(spproj_glob%os_mic%get_noris(), 10)
                    do iori=1, i_max
                        i_thumb = spproj_glob%os_mic%get_noris() - i_max + iori
                        call meta_preprocess_micrograph%set(path  =spproj_glob%os_mic%get_str(i_thumb, "thumb"),&
                                                            dfx   =spproj_glob%os_mic%get(i_thumb,       "dfx"),&
                                                            dfy   =spproj_glob%os_mic%get(i_thumb,       "dfy"),&
                                                            ctfres=spproj_glob%os_mic%get(i_thumb,    "ctfres"),&
                                                            i_max =i_max                                       ,&
                                                            i     =iori                                       )
                        if( meta_preprocess_micrograph%assigned() .and. mq_stream_master_in%is_active() ) then
                            call meta_preprocess_micrograph%serialise(meta_buffer)
                            call mq_stream_master_in%send(meta_buffer)
                        endif
                    end do
                end if
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
                    call write_mic_star_and_field
                else if( n_imported > nmic_star + 100 .and. l_haschanged )then
                    call write_mic_star_and_field
                    nmic_star = n_imported
                endif
            else
                ! wait & write snapshot
                if( .not.l_movies_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive...
                        call write_mic_star_and_field
                        l_haschanged = .false.
                    else
                        ! ...or wait
                        call sleep(WAITTIME)
                    endif
                endif
            endif
            ! send metadata
            average_ctfres  = 0
            average_astig   = 0
            average_icefrac = 0
            if(spproj_glob%os_mic%isthere('ctfres'))  average_ctfres  = spproj_glob%os_mic%get_avg("ctfres")
            if(spproj_glob%os_mic%isthere('icefrac')) average_icefrac = spproj_glob%os_mic%get_avg("icefrac")
            if(spproj_glob%os_mic%isthere('astig'))   average_astig   = spproj_glob%os_mic%get_avg("astig")
            call meta_preprocess%set(stage=string("finding and processing new movies"),          &
                                     movies_imported     = movie_buff%n_history,                 &
                                     movies_processed    = spproj_glob%os_mic%get_noris(),       & 
                                     movies_rejected     = n_failed_jobs,                        &
                                     movies_rate         = movie_buff%rate,                      &
                                     average_ctf_res     = average_ctfres,                       &
                                     average_ice_score   = average_icefrac,                      &
                                     average_astigmatism = average_astig,                        &
                                     cutoff_ctf_res      = params%ctfresthreshold,               &
                                     cutoff_ice_score    = params%icefracthreshold,              & 
                                 !   last_movie_imported = c_time(0_c_long),                     &
                                     cutoff_astigmatism  = params%astigthreshold)
            if( meta_preprocess%assigned() ) then
                if( mq_stream_master_in%is_active() ) then
                    write(logfhandle,*) "SENDING"
                    call meta_preprocess%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            endif
        !    call update_user_params(params, cline, http_communicator)
            ! update params
            if( mq_stream_master_out%is_active() ) then
                if( mq_stream_master_out%receive(meta_buffer) ) then
                    if( allocated(meta_buffer) ) then
                        ! add message back to queue for other processes
                        call mq_stream_master_out%send(meta_buffer)
                        ! deserialise buffer into meta_update
                        meta_update = transfer(meta_buffer, meta_update)
                        write(logfhandle, *) "GOT UPDATE", meta_update%get_ctfres_update()
                        if( abs(meta_update%get_ctfres_update()) > 0.001 .and. meta_update%get_ctfres_update() /= params%ctfresthreshold) then
                            params%ctfresthreshold = meta_update%get_ctfres_update()
                            params%updated         = 'yes' ! not sure this is still neccessary
                            call cline_exec%set('ctfresthreshold', params%ctfresthreshold)
                            write(logfhandle,'(A,F8.2)')'>>> CTF RESOLUTION THRESHOLD UPDATED TO: ', params%ctfresthreshold
                        endif
                        if( abs(meta_update%get_astigmatism_update()) > 0.001 .and. meta_update%get_astigmatism_update() /= params%astigthreshold) then
                            params%astigthreshold = meta_update%get_astigmatism_update()
                            params%updated        = 'yes' ! not sure this is still neccessary
                            call cline_exec%set('astigthreshold', params%astigthreshold)
                            write(logfhandle,'(A,F8.2)')'>>> ASTIGMATISM THRESHOLD UPDATED TO: ', params%astigthreshold
                        endif
                        if( abs(meta_update%get_icescore_update()) > 0.001 .and. meta_update%get_icescore_update() /= params%icefracthreshold) then
                            params%icefracthreshold = meta_update%get_icescore_update()
                            params%updated          = 'yes' ! not sure this is still neccessary
                            call cline_exec%set('icefracthreshold', params%icefracthreshold)
                            write(logfhandle,'(A,F8.2)')'>>> ICE SCORE THRESHOLD UPDATED TO: ', params%icefracthreshold
                        endif
                    endif
                endif
            endif
            call flush(logfhandle)
        end do
        ! termination
        call update_user_params(params, cline)
        if( spproj_glob%os_mic%get_noris() > 0) call write_mic_star_and_field(write_field=.true., copy_optics=.true.)
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup(params)
        ! end gracefully
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
                    call starproj_stream%stream_export_micrographs(params, spproj_glob, params%outdir, optics_set=.true.)
                else
                    call starproj_stream%stream_export_micrographs(params, spproj_glob, params%outdir)
                end if
                if( l_wfield )then
                    call spproj_glob%write_segment_inside('mic', params%projfile)
                    call spproj_glob%write_non_data_segments(params%projfile)
                endif
            end subroutine write_mic_star_and_field

            ! returns list of completed jobs
            subroutine update_projects_list( nimported )
                integer,          intent(out) :: nimported
                type(sp_project), allocatable :: streamspprojs(:)
                type(string),     allocatable :: completed_fnames(:)
                logical,          allocatable :: mics_mask(:)
                type(string) :: fname, abs_fname
                integer      :: i, n_spprojs, n_old, j, n2import, n_completed, iproj, nmics, imic, cnt
                n_completed = 0
                nimported   = 0
                n_spprojs = size(completed_jobs_clines) ! projects to import
                if( n_spprojs == 0 )return
                n_old = spproj_glob%os_mic%get_noris()  ! previously processed mmovies
                nmics = STREAM_NMOVS_SET * n_spprojs    ! incoming number of processed movies
                allocate(streamspprojs(n_spprojs), completed_fnames(n_spprojs), mics_mask(nmics))
                ! read all
                do iproj = 1,n_spprojs
                    fname     = output_dir//completed_jobs_clines(iproj)%get_carg('projfile')
                    abs_fname = simple_abspath(fname)
                    completed_fnames(iproj) = abs_fname
                    call streamspprojs(iproj)%read_segment('mic', completed_fnames(iproj))
                    cnt = 0
                    do imic = (iproj-1)*STREAM_NMOVS_SET+1, iproj*STREAM_NMOVS_SET
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
                        do i = 1,STREAM_NMOVS_SET
                            imic = imic+1
                            if( mics_mask(imic) )then
                                j = j + 1
                                call spproj_glob%os_mic%transfer_ori(j, streamspprojs(iproj)%os_mic, i)
                            endif
                        enddo
                        call streamspprojs(iproj)%kill
                    enddo
                endif
                ! finally we move the completed projects to appropriate directory
                do iproj = 1,n_spprojs
                    imic = (iproj-1)*STREAM_NMOVS_SET+1
                    if( any(mics_mask(imic:imic+STREAM_NMOVS_SET-1)) )then
                        fname = string(DIR_STREAM_COMPLETED)//basename(completed_fnames(iproj))
                        call simple_rename(completed_fnames(iproj), fname)
                    endif
                enddo
                ! cleanup
                call completed_jobs_clines(:)%kill
                deallocate(completed_jobs_clines,streamspprojs,mics_mask,completed_fnames)
            end subroutine update_projects_list

            subroutine create_movies_set_project( movie_names )
                type(string),  intent(in) :: movie_names(STREAM_NMOVS_SET)
                type(sp_project)          :: spproj_here
                type(ctfparams)           :: ctfvars
                type(string)              :: projname, projfile, xmlfile, xmldir, cwd, cwd_old
                integer :: imov
                cwd_old = CWD_GLOB
                call simple_chdir(output_dir)
                call simple_getcwd(cwd)
                CWD_GLOB = cwd%to_char()
                ! movies set
                movies_set_counter = movies_set_counter + 1
                projname   = int2str_pad(movies_set_counter,params%numlen)
                projfile   = projname//METADATA_EXT
                call cline_exec%set('projname', projname)
                call cline_exec%set('projfile', projfile)
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
                call spproj_here%add_movies(movie_names(1:STREAM_NMOVS_SET), ctfvars, verbose = .false.)
                do imov = 1,STREAM_NMOVS_SET
                    import_counter = import_counter + 1
                    call spproj_here%os_mic%set(imov, "importind", real(import_counter))
                    call spproj_here%os_mic%set(imov, "tiltgrp",   0.0)
                    call spproj_here%os_mic%set(imov, "shiftx",    0.0)
                    call spproj_here%os_mic%set(imov, "shifty",    0.0)
                    call spproj_here%os_mic%set(imov, "flsht",     0.0)
                    if(cline%defined('dir_meta')) then
                        if( SJ_directory_structure )then
                            ! multiple folders
                            xmldir = stemname(movie_names(imov))
                        else
                            ! single folder
                            xmldir = params%dir_meta
                        endif
                        xmlfile = basename(movie_names(imov))
                        if(xmlfile%substr_ind('_fractions') > 0) xmlfile = xmlfile%to_char([1,xmlfile%substr_ind('_fractions') - 1])
                        if(xmlfile%substr_ind('_EER')       > 0) xmlfile = xmlfile%to_char([1,xmlfile%substr_ind('_EER')       - 1])
                        xmlfile = xmldir//'/'//xmlfile//'.xml'
                        call spproj_here%os_mic%set(imov, "meta", xmlfile)
                    end if
                enddo
                call spproj_here%write
                call simple_chdir(cwd_old)
                CWD_GLOB = cwd_old%to_char()
                call spproj_here%kill
            end subroutine create_movies_set_project

            !>  import previous movies and updates global project & variables
            subroutine import_previous_projects
                type(sp_project),          allocatable :: spprojs(:)
                type(string), allocatable :: completed_fnames(:)
                logical,      allocatable :: mics_mask(:)
                type(string) :: fname, moviename
                integer :: n_spprojs, iproj, nmics, imic, jmic, cnt, iostat,id
                ! previously completed projects
                call simple_list_files_regexp(string(DIR_STREAM_COMPLETED), '\.simple$', completed_fnames)
                if( .not.allocated(completed_fnames) ) then
                    return ! nothing was previously completed
                endif
                if( size(completed_fnames)==0 ) then
                    deallocate(completed_fnames)    
                    return ! nothing was previously completed
                endif
                n_spprojs = size(completed_fnames)
                ! import into global project
                allocate(spprojs(n_spprojs), mics_mask(n_spprojs*STREAM_NMOVS_SET))
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
                        fname = string(DIR_STREAM_COMPLETED)//completed_fnames(iproj)
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
                    fname = basename(completed_fnames(iproj))
                    fname = get_fbody(fname,METADATA_EXT,separator=.false.)
                    id    = str2int(fname)
                    if( iostat==0 ) movies_set_counter = max(movies_set_counter, id)
                enddo
                ! update import id counter
                import_counter = spproj_glob%os_mic%get_noris()
                ! add previous movies to history
                do imic = 1,spproj_glob%os_mic%get_noris()
                    moviename = spproj_glob%os_mic%get_str(imic,'movie')
                    call movie_buff%add2history(moviename)
                enddo
                ! tidy files
                call simple_rmdir(DIR_STREAM)
                write(logfhandle,'(A,I6,A)')'>>> IMPORTED ',nmics,' PREVIOUSLY PROCESSED MOVIES'
            end subroutine import_previous_projects

            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED (PREPROCESSING)'
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p01_preprocess

end module simple_stream_p01_preprocess_new