module simple_stream_p04_refpick_extract
include 'simple_lib.f08'
use simple_cmdline,               only: cmdline
use simple_commander_base,        only: commander_base
use simple_commanders_preprocess, only: commander_make_pickrefs
use simple_parameters,            only: parameters
use simple_qsys_env,              only: qsys_env
use simple_sp_project,            only: sp_project
use simple_starproject_stream,    only: starproject_stream
use simple_stream_watcher,        only: stream_watcher
use simple_gui_utils
use simple_progress
use simple_qsys_funs
use simple_rec_list
use simple_stream_communicator
use simple_stream_utils
use simple_timer
use json_kinds
implicit none

public :: stream_p04_refpick_extract
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p04_refpick_extract
  contains
    procedure :: execute => exec_stream_pick_extract
end type stream_p04_refpick_extract

contains

    subroutine exec_stream_pick_extract( self, cline )
        class(stream_p04_refpick_extract), intent(inout) :: self
        class(cmdline),                    intent(inout) :: cline
        type(commander_make_pickrefs)          :: xmake_pickrefs
        type(parameters)                       :: params
        logical,                   parameter   :: DEBUG_HERE    = .false.
        class(cmdline),            allocatable :: completed_jobs_clines(:), failed_jobs_clines(:)
        type(rec_list)                         :: project_list, project_list_main
        type(qsys_env),            pointer     :: qenv
        type(json_value),          pointer     :: latest_picked_micrographs, latest_extracted_particles, picking_templates!, picking_diameters
        type(json_core)                        :: json
        type(qsys_env),            target      :: qenv_main, qenv_interactive
        type(cmdline)                          :: cline_make_pickrefs, cline_pick_extract
        type(stream_watcher)                   :: project_buff
        type(sp_project)                       :: spproj_glob, stream_spproj
        type(starproject_stream)               :: starproj_stream
        type(stream_http_communicator)         :: http_communicator
        type(string),              allocatable :: projects(:)
        type(string)                           :: odir, odir_extract, odir_picker, odir_completed, str_mic, str_box
        type(string)                           :: odir_interactive, odir_interactive_picker, odir_interactive_completed, cwd_job, str_dir
        character(len=STDLEN)                  :: pick_nthr_env, pick_part_env
        real                                   :: pickrefs_thumbnail_scale
        integer                                :: nmics_sel, nmics_rej, nmics_rejected_glob, pick_extract_set_counter, i_max, i_thumb, i
        integer                                :: nmics, nprojects, stacksz, prev_stacksz, iter, last_injection, iproj, envlen
        integer                                :: cnt, n_imported, n_added, nptcls_glob, n_failed_jobs, n_fail_iter, nmic_star, thumbid_offset
        integer                                :: n_pickrefs, thumbcount, xtile, ytile, xtiles, ytiles
        logical                                :: l_templates_provided, l_projects_left, l_haschanged, l_extract, l_once
        logical                                :: pause_import, l_interactive, interactive_waiting, found, l_restart
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
        if( .not. cline%defined('pick_roi')        ) call cline%set('pick_roi',         'yes')
        if( .not. cline%defined('backgr_subtr')    ) call cline%set('backgr_subtr',     'no')
        ! extraction
        if( .not. cline%defined('pcontrast')       ) call cline%set('pcontrast',        'black')
        if( .not. cline%defined('extractfrommov')  ) call cline%set('extractfrommov',   'no')
        ! ev overrides
        call get_environment_variable(SIMPLE_STREAM_PICK_NTHR, pick_nthr_env, envlen)
        if(envlen > 0)  call cline%set('nthr', str2int(pick_nthr_env))
        ! sanity check for restart
        l_restart = .false.
        if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
            l_restart = .true.
        endif
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                str_dir = cline%get_carg('dir_exec')
                THROW_HARD('Previous directory does not exists: '//str_dir%to_char())
                call str_dir%kill
            endif
            l_restart = .true.
        endif
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            if(cline%get_carg('interactive') .eq. 'yes') then
                call cline%set('projname', 'initial_picking')
                call cline%set('projfile', 'initial_picking.simple')
            else
                call cline%set('projname', 'pick_extract')
                call cline%set('projfile', 'pick_extract.simple')
            endif
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
        ! picking
        l_interactive = params%interactive == 'yes'
        ! http communicator init
        if(l_interactive) then
            call http_communicator%create(params%niceprocid, params%niceserver%to_char(), "initial_picking")
        else
            call http_communicator%create(params%niceprocid, params%niceserver%to_char(), "pick_extract")   
        endif
        call communicator_init()
        call http_communicator%send_jobstats()
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs', '**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs_completed', '**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
        l_extract            = .true.
        l_templates_provided = cline%defined('pickrefs')
        if( l_templates_provided )then
            if( .not.file_exists(params%pickrefs) ) then
                write(logfhandle,'(A,F8.2)')'>>> WAITING UP TO 24 HOURS FOR '//params%pickrefs%to_char()
                do i=1, 8640
                    if(file_exists(params%pickrefs)) exit
                    call sleep(10)
                    call http_communicator%send_jobstats()
                    if( http_communicator%exit_status() )then
                        ! termination
                        write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                        call http_communicator%term()
                        call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                        call EXIT(0)
                    endif
                end do
            endif
            write(logfhandle,'(A)')'>>> PERFORMING REFERENCE-BASED PICKING'
            if( cline%defined('moldiam') )then
                call cline%delete('moldiam')
                write(logfhandle,'(A)')'>>> MOLDIAM IGNORED'
            endif
            ! test for STREAM_NMICS file in pickrefs folder and copy here if it exists
            if(file_exists(stemname(params%pickrefs) // '/' // STREAM_NMICS)) then
                call simple_copy_file(stemname(params%pickrefs) // '/' // STREAM_NMICS, string(STREAM_NMICS))
            endif
        else if( .not.cline%defined('moldiam') )then
            THROW_HARD('MOLDIAM required for picker=new reference-free picking')
            write(logfhandle,'(A)')'>>> PERFORMING SINGLE DIAMETER PICKING'
        endif
        ! endif
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! movie watcher init
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! directories structure & restart
        odir                       = DIR_STREAM
        odir_completed             = DIR_STREAM_COMPLETED
        odir_picker                = PATH_HERE//DIR_PICKER
        odir_interactive           = 'interactive/'//DIR_STREAM
        odir_interactive_picker    = 'interactive/'//DIR_PICKER
        odir_interactive_completed = 'interactive/'//DIR_STREAM_COMPLETED
        odir_extract               = PATH_HERE//DIR_EXTRACT
        pick_extract_set_counter = 0    ! global counter of projects to be processed
        nptcls_glob              = 0    ! global number of particles
        nmics_rejected_glob      = 0    ! global number of micrographs rejected
        nmic_star                = 0
        thumbid_offset           = 0
        thumbcount               = 0
        if( l_restart )then
            call del_file(TERM_STREAM)
            if(cline%defined('dir_exec')) call cline%delete('dir_exec')
            call simple_rmdir(odir)
            if( params%clear .eq. "yes")then
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
                call http_communicator%update_json("stage", "importing previously processed data", found)
                call http_communicator%send_jobstats()
                call import_previous_mics( project_list )
                if( project_list%size() > 0 )then
                    nptcls_glob = project_list%get_nptcls_tot()
                    nmic_star   = spproj_glob%os_mic%get_noris()
                    ! http stats
                    call http_communicator%update_json("micrographs_imported",  spproj_glob%os_mic%get_noris(), found)
                    call http_communicator%update_json("micrographs_processed", spproj_glob%os_mic%get_noris(), found)
                    call http_communicator%update_json("movies_rejected",       0,                              found)
                endif
            endif
        endif
        ! make directories structure
        call simple_mkdir(odir)
        call simple_mkdir(odir//STDERROUT_DIR)
        call simple_mkdir(odir_completed)
        call simple_mkdir(odir_picker)
        if( l_extract ) call simple_mkdir(odir_extract)
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PICK_PARTITION, pick_part_env, envlen)
        if(envlen > 0) then
            call qenv_main%new(1,stream=.true.,qsys_partition=string(trim(pick_part_env)))
            call qenv_interactive%new(1,stream=.true.,qsys_partition=string(trim(pick_part_env)))
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
        call cline_pick_extract%set('dir', PATH_PARENT)
        if( l_extract )then
            call cline_pick_extract%set('extract','yes')
        else
            call cline_pick_extract%set('extract','no')
        endif
        ! Infinite loop
        last_injection  = simple_gettime()
        prev_stacksz    = 0
        iter            = 0
        n_imported      = 0   ! global number of imported processed micrographs
        n_failed_jobs   = 0
        n_added         = 0   ! global number of micrographs added to processing stack
        l_projects_left = .false.
        l_haschanged    = .false.
        l_once          = .true.
        pause_import    = .false.
        do
            if( file_exists(TERM_STREAM) .or. http_communicator%exit_status()) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop_status() )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj_glob%kill
                call qsys_cleanup
                call simple_end('**** SIMPLE_STREAM_PICK_EXTRACT USER STOP ****')
                call EXIT(0)
            endif
            iter = iter + 1
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
                            call xmake_pickrefs%execute_safe(cline_make_pickrefs)
                            call cline_pick_extract%set('pickrefs', '../'//PICKREFS_FBODY//params%ext%to_char())
                            call mrc2jpeg_tiled(string(PICKREFS_FBODY)//params%ext, string(PICKREFS_FBODY)//".jpeg",&
                            &scale=pickrefs_thumbnail_scale, ntiles=n_pickrefs, n_xtiles=xtiles, n_ytiles=ytiles)
                            ! http stats
                            xtile = 0
                            ytile = 0
                            do i=0, n_pickrefs - 1
                                call json_add_picking_template(cwd_job%to_char()//'/'//PICKREFS_FBODY//".jpeg",&
                                &xtile * 100, ytile * 100, 100 * ytiles, 100 * xtiles)
                                xtile = xtile + 1
                                if(xtile .eq. xtiles) then
                                    xtile = 0
                                    ytile = ytile + 1
                                endif
                            end do
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
                ! http stats
                call http_communicator%update_json("micrographs_imported",     project_buff%n_history * STREAM_NMOVS_SET, found)
                call http_communicator%update_json("last_micrograph_imported", stream_datestr(), found)
            else
                l_projects_left = .false.
            endif
            ! submit jobs
            call qenv%qscripts%schedule_streaming( qenv%qdescr, path=odir )
            stacksz = qenv%qscripts%get_stacksz()
            if( stacksz .ne. prev_stacksz )then
                prev_stacksz = stacksz                          ! # of projects
                stacksz      = qenv%qscripts%get_stack_range()  ! # of micrographs
                write(logfhandle,'(A,I6)')'>>> MICROGRAPHS TO PROCESS:                 ', stacksz
            endif
            ! fetch completed jobs list & updates
            if( qenv%qscripts%get_done_stacksz() > 0 )then
                call qenv%qscripts%get_stream_done_stack( completed_jobs_clines )
                call update_projects_list( project_list, n_imported )
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
                write(logfhandle,'(A,I3,A2,I3)') '>>> # OF COMPUTING UNITS IN USE/TOTAL   : ',qenv%get_navail_computing_units(),'/ ',params%nparts
                if( n_failed_jobs > 0 ) write(logfhandle,'(A,I8)') '>>> # DESELECTED MICROGRAPHS/FAILED JOBS: ',n_failed_jobs
                if(.not. l_interactive) then
                    ! http stats   
                    call http_communicator%update_json("stage",               "finding, picking and extracting micrographs",  found)  
                    call http_communicator%update_json("box_size",              params%box,                                   found)
                    call http_communicator%update_json("particles_extracted",   nptcls_glob,                                  found)
                    call http_communicator%update_json("particles_per_mic",     nint(float(nptcls_glob) / float(n_imported)), found)
                    call http_communicator%update_json("micrographs_processed", n_imported,                                   found)
                    call http_communicator%update_json("micrographs_rejected",  n_failed_jobs + nmics_rejected_glob,          found)
                    if(spproj_glob%os_mic%isthere('thumb_den') .and. spproj_glob%os_mic%isthere('xdim') .and. spproj_glob%os_mic%isthere('ydim') \
                    .and. spproj_glob%os_mic%isthere('smpd') .and. spproj_glob%os_mic%isthere('boxfile')) then
                        i_max = 10
                        if(spproj_glob%os_mic%get_noris() < i_max) i_max = spproj_glob%os_mic%get_noris()
                        ! create an empty latest_picked_micrographs json array
                        call json%remove(latest_picked_micrographs, destroy=.true.)
                        call json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                        call http_communicator%add_to_json(latest_picked_micrographs)
                        do i_thumb=0, i_max - 1
                            str_mic = spproj_glob%os_mic%get_str(spproj_glob%os_mic%get_noris() - i_thumb, "thumb_den")
                            str_box = spproj_glob%os_mic%get_str(spproj_glob%os_mic%get_noris() - i_thumb, "boxfile")
                            call json_add_micrograph(str_mic, 100, 0, 100, 200,\
                            nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "xdim")),\
                            nint(spproj_glob%os_mic%get(spproj_glob%os_mic%get_noris() - i_thumb, "ydim")), str_box)
                            call str_mic%kill
                            call str_box%kill
                        end do
                    endif
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
            if(interactive_waiting) then
                call sleep(1)
            else
                call sleep(WAITTIME)
            end if
            ! http stats send
            call http_communicator%send_jobstats()
            call update_user_params(cline, http_communicator)
        end do
        ! termination
        call http_communicator%update_json("stage", "terminating", found)
        call http_communicator%send_jobstats()
        call write_project
        ! write star files (just in case you want ot import these particles/micrographs elsewhere)
        call spproj_glob%write_mics_star(string("micrographs.star"))
        call spproj_glob%write_ptcl2D_star(string("particles.star"))
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup
        call project_list%kill
        call project_list_main%kill
        ! end gracefully
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
                    call starproj_stream%stream_export_micrographs(spproj_glob, params%outdir, optics_set=l_optics_set)
                    if( DEBUG_HERE )then
                        ms_export = toc(ms0)
                        print *,'ms_export  : ', ms_export; call flush(6)
                    endif
                end if
            end subroutine write_micrographs_starfile

            subroutine write_project()
                integer, allocatable :: fromps(:)
                integer              :: nptcls,fromp,top,i,iptcl,nmics,imic,micind,optics_map_id
                type(string)         :: prev_projname, mapfileprefix
                type(project_rec)    :: prec
                type(rec_iterator)   :: it
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
                    it            = project_list%begin()
                    do imic = 1,nmics
                        ! retrieve one record from the list with the iterator
                        call it%get(prec)
                        if( prec%projname /= prev_projname )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('stk', prec%projname)
                            prev_projname = prec%projname
                        endif
                        micind       = prec%micind
                        fromps(imic) = stream_spproj%os_stk%get_fromp(micind) ! fromp from individual project
                        fromp        = nptcls + 1
                        nptcls       = nptcls + prec%nptcls
                        top          = nptcls
                        call spproj_glob%os_stk%transfer_ori(imic, stream_spproj%os_stk, micind)
                        call spproj_glob%os_stk%set(imic, 'fromp',fromp)
                        call spproj_glob%os_stk%set(imic, 'top',  top)
                        ! move the iterator
                        call it%next()
                    enddo
                    call spproj_glob%write_segment_inside('stk', params%projfile)
                    ! particles
                    call spproj_glob%os_ptcl2D%new(nptcls, is_ptcl=.true.)
                    iptcl         = 0
                    prev_projname = ''
                    it            = project_list%begin()
                    do imic = 1,nmics
                        ! retrieve one record from the list with the iterator
                        call it%get(prec)
                        if( prec%projname /= prev_projname )then
                            call stream_spproj%kill
                            call stream_spproj%read_segment('ptcl2D', prec%projname)
                            prev_projname = prec%projname
                        endif
                        fromp = fromps(imic)
                        top   = fromp + prec%nptcls - 1
                        do i = fromp,top
                            iptcl = iptcl + 1
                            call spproj_glob%os_ptcl2D%transfer_ori(iptcl, stream_spproj%os_ptcl2D, i)
                            call spproj_glob%os_ptcl2D%set_stkind(iptcl, imic)
                        enddo
                        ! move the iterator
                        call it%next()
                    enddo
                    call stream_spproj%kill
                    write(logfhandle,'(A,I8)')'>>> # PARTICLES EXTRACTED:          ',spproj_glob%os_ptcl2D%get_noris()
                    call spproj_glob%write_segment_inside('ptcl2D', params%projfile)
                    spproj_glob%os_ptcl3D = spproj_glob%os_ptcl2D
                    call spproj_glob%os_ptcl3D%delete_2Dclustering
                    call spproj_glob%write_segment_inside('ptcl3D', params%projfile)
                    call spproj_glob%os_ptcl3D%kill
                endif
                ! add optics
                if(cline%defined('optics_dir')) then
                    optics_map_id = get_latest_optics_map_id(params%optics_dir)
                    if(optics_map_id .gt. 0) then
                        mapfileprefix = params%optics_dir // '/' // OPTICS_MAP_PREFIX // int2str(optics_map_id)
                        write(logfhandle,'(A,I8)')'>>> IMPORTING OPTICS FROM : ' // mapfileprefix%to_char()
                        call spproj_glob%import_optics_map(mapfileprefix)
                        call spproj_glob%write_segment_inside('optics', params%projfile)
                    endif
                endif
                call spproj_glob%write_non_data_segments(params%projfile)
                ! benchmark
                if( DEBUG_HERE )then
                    rt_write = toc(t0)
                    print *,'rt_write  : ', rt_write; call flush(6)
                endif
            end subroutine write_project

            ! updates global project, returns records of processed micrographs
            subroutine update_projects_list( project_list, nimported )
                class(rec_list), intent(inout) :: project_list
                integer,         intent(inout) :: nimported
                type(sp_project), allocatable :: spprojs(:)
                type(string)       :: fname, abs_fname, new_fname
                type(sp_project)  :: tmpproj
                type(project_rec) :: prec
                integer :: n_spprojs, n_old, j, nprev_imports, n_completed, nptcls, nmics, imic, iproj
                n_completed = 0
                nimported   = 0
                ! previously imported
                n_old = project_list%size()
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
                    else
                        call spproj_glob%os_mic%reallocate(n_completed)
                    endif
                    ! update records and global project
                    j = n_old
                    do iproj = 1,n_spprojs
                        if( spprojs(iproj)%os_mic%get_noris() == 0 ) cycle
                        ! move project to appropriate directory
                        fname = filepath(odir, completed_jobs_clines(iproj)%get_carg('projfile'))
                        new_fname = filepath(odir_completed, completed_jobs_clines(iproj)%get_carg('projfile'))
                        call simple_rename(fname, new_fname)
                        abs_fname = simple_abspath(new_fname)
                        ! records & project
                        do imic = 1,spprojs(iproj)%os_mic%get_noris()
                            j = j + 1
                            prec%projname = abs_fname
                            prec%micind   = imic
                            nptcls        = spprojs(iproj)%os_mic%get_int(imic,'nptcls')
                            nptcls_glob   = nptcls_glob + nptcls ! global update
                            prec%nptcls   = nptcls
                            call spproj_glob%os_mic%transfer_ori(j, spprojs(iproj)%os_mic, imic)
                            call project_list%push_back(prec)
                        enddo
                    enddo
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
                class(string), intent(in)  :: project_fname
                integer,       intent(out) :: nselected, nrejected
                type(sp_project)           :: tmp_proj, spproj_here
                integer,      allocatable  :: states(:)
                type(string) :: proj_fname, projname, projfile, path
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
                path       = trim(CWD_GLOB)//'/'//odir%to_char()
                proj_fname = basename(project_fname)
                projname   = get_fbody(proj_fname, METADATA_EXT, separator=.false.)
                projfile   = projname//METADATA_EXT
                call spproj_here%projinfo%new(1, is_ptcl=.false.)
                call spproj_here%projinfo%set(1,'projname', projname)
                call spproj_here%projinfo%set(1,'projfile', projfile)
                call spproj_here%projinfo%set(1,'cwd',      path)
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
                ! cleanup from inipick preprocessing
                do imic = 1,spproj_here%os_mic%get_noris()
                    if(spproj_here%os_mic%isthere(imic, 'mic_den'))  call spproj_here%os_mic%delete_entry(imic, 'mic_den')
                    if(spproj_here%os_mic%isthere(imic, 'mic_topo')) call spproj_here%os_mic%delete_entry(imic, 'mic_topo')
                    if(spproj_here%os_mic%isthere(imic, 'mic_bin'))  call spproj_here%os_mic%delete_entry(imic, 'mic_bin')
                    if(spproj_here%os_mic%isthere(imic, 'mic_diam')) call spproj_here%os_mic%delete_entry(imic, 'mic_diam')
                enddo
                ! update for execution
                pick_extract_set_counter = pick_extract_set_counter + 1
                projname = int2str_pad(pick_extract_set_counter,params%numlen)
                projfile = projname//METADATA_EXT
                call cline_pick_extract%set('projname', projname)
                call cline_pick_extract%set('projfile', projfile)
                call cline_pick_extract%set('fromp',    1)
                call cline_pick_extract%set('top',      nselected)
                call spproj_here%write(path//'/'//projfile)
                call spproj_here%kill
                call tmp_proj%kill
            end subroutine create_individual_project

            !>  import previous run to the current project and reselect micrographs
            subroutine import_previous_mics( project_list )
                type(rec_list), intent(inout) :: project_list
                type(sp_project), allocatable :: spprojs(:)
                type(string),     allocatable :: completed_fnames(:)
                logical, allocatable :: mics_mask(:)
                type(project_rec) :: prec
                type(string)      :: fname
                integer           :: n_spprojs, iproj, nmics, imic, jmic, iostat,id, nsel_mics, irec
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
                call spproj_glob%os_mic%new(nsel_mics, is_ptcl=.false.)
                irec = 0
                jmic = 0
                do iproj = 1,n_spprojs
                    do imic = 1, spprojs(iproj)%os_mic%get_noris()
                        jmic = jmic+1
                        if( mics_mask(jmic) )then
                            irec = irec + 1
                            prec%projname = completed_fnames(iproj)
                            prec%micind   = imic
                            prec%nptcls   = spprojs(iproj)%os_mic%get_int(imic, 'nptcls')
                            call spproj_glob%os_mic%transfer_ori(irec, spprojs(iproj)%os_mic, imic)
                            call project_list%push_back(prec)
                        endif
                    enddo
                    call spprojs(iproj)%kill
                enddo
                ! update global set counter
                do iproj = 1,n_spprojs
                    fname = basename(completed_fnames(iproj))
                    fname = get_fbody(fname,METADATA_EXT,separator=.false.)
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
                call http_communicator%add_to_json("stage",               "initialising")
                call http_communicator%add_to_json("micrographs_imported",     0)
                call http_communicator%add_to_json("micrographs_processed",    0)
                call http_communicator%add_to_json("micrographs_rejected",     0)
                call http_communicator%add_to_json("particles_per_mic",        0)
                call http_communicator%add_to_json("particles_extracted",      dble(0.0))
                call http_communicator%add_to_json("box_size",                 dble(0.0))
                ! call http_communicator%add_to_json("best_diam",                dble(0.0))
                call http_communicator%add_to_json("user_input",               .false.)
                call http_communicator%add_to_json("last_micrograph_imported", "")
                call json%create_array(latest_extracted_particles, "latest_extracted_particles")
                call http_communicator%add_to_json(latest_extracted_particles)
                call json%create_array(latest_picked_micrographs, "latest_picked_micrographs")
                call http_communicator%add_to_json(latest_picked_micrographs)
                call json%create_array(picking_templates, "picking_templates")
                call http_communicator%add_to_json(picking_templates)
                ! call json%create_array(picking_diameters, "picking_diameters")
                ! call http_communicator%add_to_json(picking_diameters)
                ! call json%create_array(refinement_diameters, "refinement_diameters")
                ! call http_communicator%add_to_json(refinement_diameters)
            end subroutine communicator_init

            subroutine json_add_micrograph( path, spritex, spritey, spriteh, spritew, xdim, ydim, boxfile_path )
                class(string),     intent(in)  :: path, boxfile_path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew, xdim, ydim
                type(nrtxtfile)               :: boxfile
                type(json_value), pointer     :: micrograph, boxes, box
                real,             allocatable :: boxdata(:,:)
                integer                       :: i, x, y
                call json%create_object(micrograph, "")
                call json%add(micrograph, "path",    path%to_char())
                call json%add(micrograph, "spritex", spritex)
                call json%add(micrograph, "spritey", spritey)
                call json%add(micrograph, "spriteh", spriteh)
                call json%add(micrograph, "spritew", spritew)
                call json%add(micrograph, "xdim"   , xdim)
                call json%add(micrograph, "ydim",    ydim)
                call json%create_array(boxes, "boxes")
                call boxfile%new(boxfile_path, 1)
                allocate(boxdata(boxfile%get_nrecs_per_line(), boxfile%get_ndatalines()))
                if(boxfile%get_nrecs_per_line() == 5) then
                    ! standard boxfile
                    do i=1, boxfile%get_ndatalines()
                        call boxfile%readNextDataLine(boxdata(:,i))
                        call json%create_object(box, "")
                        x = nint(boxdata(1,i) + boxdata(3,i)/2)
                        y = nint(boxdata(2,i) + boxdata(4,i)/2)
                        call json%add(box, "x",    x)
                        call json%add(box, "y",    y)
                        call json%add(boxes, box)
                    enddo
                else if(boxfile%get_nrecs_per_line() == 6) then
                    THROW_HARD('Invalid boxfile format!')
                endif
                call boxfile%kill()
                deallocate(boxdata)
                call json%add(micrograph, boxes)
                call json%add(latest_picked_micrographs, micrograph)
            end subroutine json_add_micrograph

            subroutine json_add_picking_template( path, spritex, spritey, spriteh, spritew )
                character(*),     intent(in)  :: path
                integer,          intent(in)  :: spritex, spritey, spriteh, spritew
                type(json_value), pointer     :: template
                call json%create_object(template, "")
                call json%add(template, "path",    path)
                call json%add(template, "spritex", spritex)
                call json%add(template, "spritey", spritey)
                call json%add(template, "spriteh", spriteh)
                call json%add(template, "spritew", spritew)
                call json%add(picking_templates, template)
            end subroutine json_add_picking_template

    end subroutine exec_stream_pick_extract

end module simple_stream_p04_refpick_extract