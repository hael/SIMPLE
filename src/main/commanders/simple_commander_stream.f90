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
use FoX_dom
implicit none

public :: commander_stream

private
#include "simple_local_flags.inc"

character(len=STDLEN), parameter :: dir_preprocess  = trim(PATH_HERE)//'spprojs/'
character(len=STDLEN), parameter :: micspproj_fname = './streamdata.simple'
character(len=STDLEN), parameter :: USER_PARAMS     = 'stream_user_params.txt'


type, extends(commander_base) :: commander_stream
  contains
    procedure :: execute => exec_simple_stream
end type commander_stream

contains

    subroutine exec_simple_stream( self, cline )
        use simple_moviewatcher, only: moviewatcher
        use simple_timer
        class(commander_stream), intent(inout) :: self
        class(cmdline),          intent(inout) :: cline
        type(parameters)                       :: params
        type(guistats)                         :: gui_stats
        integer,                   parameter   :: WAITTIME        = 3    ! movie folder watched every WAITTIME seconds
        integer,                   parameter   :: LONGTIME        = 300  ! time lag after which a movie is processed
        integer,                   parameter   :: INACTIVE_TIME   = 900  ! inactive time trigger for writing project file
        logical,                   parameter   :: DEBUG_HERE      = .false.
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(qsys_env)                         :: qenv
        type(moviewatcher)                     :: movie_buff
        type(sp_project)                       :: spproj, stream_spproj, tmp_spproj
        type(starproject)                      :: starproj
        character(len=LONGSTRLEN), allocatable :: movies(:), completed_fnames(:)
        character(len=:),          allocatable :: output_dir, output_dir_ctf_estimate, output_dir_motion_correct
        character(len=LONGSTRLEN)              :: movie
        integer                                :: nmovies, imovie, stacksz, prev_stacksz, iter, last_injection, iproj
        integer                                :: cnt, n_imported, n_added, n_failed_jobs, n_fail_iter, nmic_star
        logical                                :: l_movies_left, l_haschanged
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
        if( .not. cline%defined('mcconvention')     ) call cline%set('mcconvention','simple')
        if( .not. cline%defined('eer_upsampling')   ) call cline%set('eer_upsampling',    1.)
        if( .not. cline%defined('algorithm')        ) call cline%set('algorithm',    'patch')
        if( .not. cline%defined('mcpatch')          ) call cline%set('mcpatch',        'yes')
        if( .not. cline%defined('mcpatch_thres')    ) call cline%set('mcpatch_thres','  yes')
        if( .not. cline%defined('tilt_thres')       ) call cline%set('tilt_thres',      0.05)
        if( .not. cline%defined('beamtilt')         ) call cline%set('beamtilt',        'no')
        call cline%set('groupframes', 'no')
        ! ctf estimation
        if( .not. cline%defined('pspecsz')          ) call cline%set('pspecsz',          512.)
        if( .not. cline%defined('hp_ctf_estimate')  ) call cline%set('hp_ctf_estimate',  HP_CTF_ESTIMATE)
        if( .not. cline%defined('lp_ctf_estimate')  ) call cline%set('lp_ctf_estimate',  LP_CTF_ESTIMATE)
        if( .not. cline%defined('dfmin')            ) call cline%set('dfmin',            DFMIN_DEFAULT)
        if( .not. cline%defined('dfmax')            ) call cline%set('dfmax',            DFMAX_DEFAULT)
        if( .not. cline%defined('ctfpatch')         ) call cline%set('ctfpatch',         'yes')
        if( .not. cline%defined('ctfresthreshold')  ) call cline%set('ctfresthreshold',  CTFRES_THRESHOLD_STREAM)
        if( .not. cline%defined('icefracthreshold') ) call cline%set('icefracthreshold', ICEFRAC_THRESHOLD_STREAM)
        ! write cmdline for GUI
        call cline%writeline(".cline")
        ! master parameters
        call cline%set('numlen', 5.)
        call cline%set('stream','yes')
        call params%new(cline)
        params_glob%split_mode = 'stream'
        params_glob%ncunits    = params%nparts
        call cline%set('mkdir', 'no')
        call cline%set('prg',   'preprocess')
        if( cline%defined('dir_prev') .and. .not.file_exists(params%dir_prev) )then
            THROW_HARD('Directory '//trim(params%dir_prev)//' does not exist!')
        endif
        ! initialise progress monitor
        call progressfile_init()
        ! master project file
        call spproj%read( params%projfile )
        call spproj%update_projinfo(cline)
        if( spproj%os_mic%get_noris() /= 0 ) THROW_HARD('PREPROCESS_STREAM must start from an empty project (eg from root project folder)')
        call spproj%write(micspproj_fname)
        ! output directories
        output_dir = trim(PATH_HERE)//trim(dir_preprocess)
        call simple_mkdir(output_dir)
        call simple_mkdir(trim(output_dir)//trim(STDERROUT_DIR))
        output_dir_ctf_estimate   = filepath(trim(PATH_HERE), trim(DIR_CTF_ESTIMATE))
        output_dir_motion_correct = filepath(trim(PATH_HERE), trim(DIR_MOTION_CORRECT))
        call simple_mkdir(output_dir_ctf_estimate,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call simple_mkdir(output_dir_motion_correct,errmsg="commander_stream_wflows :: exec_preprocess_stream;  ")
        call cline%set('dir','../')
        ! setup the environment for distributed execution
        call qenv%new(1,stream=.true.)
        ! movie watcher init
        movie_buff = moviewatcher(LONGTIME)
        ! import previous runs
        call import_prev_streams
        ! start watching
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
        do
            ! guistats init each loop
            call gui_stats%init(nlines=2)
            call gui_stats%set(1, "title", "micrographs")
            if( file_exists(trim(TERM_STREAM)) )then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PREPROCESS STREAM'
                exit
            endif
            iter = iter + 1
            ! detection of new movies
            call movie_buff%watch( nmovies, movies, max_nmovies=params%nparts )
            ! append movies to processing stack
            if( nmovies > 0 )then
                cnt = 0
                do imovie = 1, nmovies
                    movie = trim(adjustl(movies(imovie)))
                    call create_individual_project(movie)
                    call qenv%qscripts%add_to_streaming( cline )
                    call qenv%qscripts%schedule_streaming( qenv%qdescr, path=output_dir )
                    call movie_buff%add2history( movies(imovie) )
                    cnt = cnt+1
                    n_added = n_added+1
                    if( cnt == min(params%nparts,nmovies) ) exit
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
                write(logfhandle,'(A,I6)')'>>> MOVIES TO PROCESS:                ', stacksz
                ! guistats
               call gui_stats%set(1, 'primary_movies', int2str(spproj%os_mic%get_noris()) // '/' // int2str(stacksz + spproj%os_mic%get_noris()))
            endif
            ! fetch completed jobs list
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
                write(logfhandle,'(A,I8)')                         '>>> # MOVIES PROCESSED & IMPORTED       : ',n_imported
                write(logfhandle,'(A,I3,A2,I3)')                   '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                ! guistats
                call gui_stats%set(1, 'primary_movies',  int2str(n_imported) // '/' // int2str(stacksz + spproj%os_mic%get_noris()))
                call gui_stats%set(1, 'secondary_compute', int2str(qenv%get_navail_computing_units()) // '/' // int2str(params%nparts))
                if( n_failed_jobs > 0 ) call gui_stats%set(1, 'primary_rejected', n_failed_jobs)
                if(spproj%os_mic%isthere("ctfres"))  call gui_stats%set(1, 'primary_avg_ctf_res', spproj%os_mic%get_avg("ctfres"))
                ! update progress monitor
                call progressfile_update(progress_estimate_preprocess_stream(n_imported, n_added))
                ! write project for gui, micrographs field only
                call spproj%write_segment_inside('mic',micspproj_fname)
                last_injection = simple_gettime()
                ! guistats
                call gui_stats%set_now(1, 'secondary_last_injection')
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
            else
                ! wait & write snapshot
                if( .not.l_movies_left )then
                    if( (simple_gettime()-last_injection > INACTIVE_TIME) .and. l_haschanged )then
                        ! write project when inactive...
                        call spproj%write_segment_inside('mic', params%projfile)
                        call spproj%write_non_data_segments(params%projfile)
                        call update_user_params(cline)
                        call write_migrographs_starfile
                        l_haschanged = .false.
                    else
                        ! ...or wait
                        call sleep(WAITTIME)
                    endif
                endif
            endif
            ! read beamtilts
            if( cline%defined('dir_meta')) call read_xml_beamtilts()
            ! guistats
            call gui_stats%merge(POOLSTATS_FILE, 2)
            call gui_stats%write
        end do
        ! termination
        call spproj%write_segment_inside('mic', params%projfile)
        call spproj%write_non_data_segments(params%projfile)
        call update_user_params(cline)
        call write_migrographs_starfile
        ! final stats
        call gui_stats%merge(POOLSTATS_FILE, 2, delete = .true.)
        call gui_stats%remove(1, 'secondary_compute')
        call gui_stats%write
        call gui_stats%kill
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
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%assign_optics(cline, spproj)
                    if( DEBUG_HERE )then
                        ms_assign = toc(ms0)
                        print *,'ms_assign  : ', ms_assign; call flush(6)
                    endif
                    if( DEBUG_HERE ) ms0 = tic()
                    call starproj%export_mics(cline, spproj)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                    if(allocated(starproj%tiltinfo)) deallocate(starproj%tiltinfo)
                end if
            end subroutine write_migrographs_starfile

            ! returns list of completed jobs
            subroutine update_projects_list( completedfnames, nimported )
                character(len=LONGSTRLEN), allocatable, intent(inout) :: completedfnames(:)
                integer,                                intent(inout) :: nimported
                type(sp_project)                       :: streamspproj
                character(len=:),          allocatable :: fname, abs_fname
                character(len=LONGSTRLEN), allocatable :: old_fnames(:)
                logical,                   allocatable :: spproj_mask(:)
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
                    call check_nptcls(fname, nptcls_here, state)
                    spproj_mask(i) = state > 0
                enddo
                n2import      = count(spproj_mask)
                n_failed_jobs = n_failed_jobs + (n_spprojs-n2import)
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
                        deallocate(completedfnames)
                        allocate(completedfnames(n_completed))
                        if( n_old > 0 )then
                            completedfnames(1:n_old) = old_fnames(:)
                        endif
                        deallocate(old_fnames)
                    endif
                    j = 0
                    do i=1,n_spprojs
                        if( .not.spproj_mask(i) ) cycle
                        j = j+1
                        fname     = trim(output_dir)//trim(completed_jobs_clines(i)%get_carg('projfile'))
                        abs_fname = simple_abspath(fname, errmsg='preprocess_stream :: update_projects_list 1')
                        completedfnames(n_old+j) = trim(abs_fname)
                        call streamspproj%read_segment('mic', abs_fname)
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
                    enddo
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
                integer :: iproj,nprojs,icnt
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
                do iproj = 1,nprojs
                    call streamspproj%read_segment('mic', sp_files(iproj) )
                    if( streamspproj%os_mic%get_noris() /= 1 )then
                        THROW_WARN('Ignoring previous project'//trim(sp_files(iproj)))
                        cycle
                    endif
                    if( .not. streamspproj%os_mic%isthere(1,'intg') )cycle
                    spproj_mask(iproj) = .true.
                enddo
                if( count(spproj_mask) == 0 ) return
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
                    ! import mic segment
                    call streamspproj%os_mic%set_ori(1, o)
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
                    ! updating STREAM_SPPROJFILES
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

    end subroutine exec_simple_stream

        !> updates current parameters with user input
    subroutine update_user_params( cline_here )
        type(cmdline), intent(inout) :: cline_here
        type(oris) :: os
        real       :: lpthres, ndev, tilt_thres
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

    subroutine check_nptcls( fname, nptcls, state )
        character(len=*), intent(in)  :: fname
        integer,          intent(out) :: nptcls, state
        type(sp_project) :: spproj_here
        integer :: nmics, nstks
        state  = 0
        call spproj_here%read_data_info(fname, nmics, nstks, nptcls)
        if( nmics /= 1 )then
            THROW_WARN('No micrograph for: '//trim(fname)//'. Skipping')
        else
            call spproj_here%read_segment('mic',fname)
            state = spproj_here%os_mic%get_state(1)
            call spproj_here%kill
        endif
    end subroutine check_nptcls

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
