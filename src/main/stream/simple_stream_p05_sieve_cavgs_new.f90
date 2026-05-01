!@descr: task 5 in the stream pipeline: multi-pass chunk-based 2D classification with automatic sieving of low-quality class averages
!==============================================================================
! MODULE: simple_stream_p05_sieve_cavgs_new
!
! PURPOSE:
!   Implements stream pipeline task 5: continuously ingests incoming project
!   files and drives a four-stage chunked2D classification pipeline that
!   produces progressively refined class averages, automatically rejecting
!   low-quality averages after each stage (sieving).
!
! TYPES:
!   stream_p05_sieve_cavgs - commander_base extension; entry point for the
!                            sieve-cavgs stream task.
!
! WORKFLOW:
!   1. Initialise parameters, queue environment, and chunked2D object.
!   2. Restore previously imported project history (if present).
!   3. Enter main loop (runs until termination signal):
!      a. Watch dir_target for newly completed project files.
!      b. Import new projects into the rec_list.
!      c. Call chunked_2D%cycle(), which performs per-cycle:
!           i.  collect_and_reject   — harvest completed chunks, sieve cavgs
!           ii. generate_microchunks_pass_1 — seed pass-1 chunks from new records
!           iii.generate_microchunks_pass_2 — promote pass-1 results to pass-2
!           iv. generate_refchunk          — build reference chunk from pass-2
!           v.  generate_microchunks_match — match-refine against reference
!           vi. submit                     — dispatch pending chunks to queue
!      d. Sleep for WAITTIME before the next cycle.
!
! PARAMETERS (hard-coded):
!   MAX_MOVIE_IMPORT    — maximum movies imported per loop cycle   (20)
!
! ENVIRONMENT:
!   SIMPLE_STREAM_CHUNK_PARTITION — queue partition for chunk jobs
!
! DEPENDENCIES:
!   simple_stream_api, simple_microchunked2D, simple_stream_pool2D_utils
!==============================================================================
module simple_stream_p05_sieve_cavgs_new
use unix,                        only: SIGTERM
use simple_stream_api
use simple_fileio,               only: read_filetable
use simple_stream_mq_defs,       only: mq_stream_master_in, mq_stream_master_out
use simple_microchunked2D,       only: microchunked2D
use simple_stream_pool2D_utils,  only: set_lpthres_type
use simple_gui_metadata_api,     only: gui_metadata_cavg2D,                                    &
                                       gui_metadata_stream_particle_sieving,                   &
                                       gui_metadata_stream_update,                             &
                                       sprite_sheet_pos,                                       &
                                       GUI_METADATA_STREAM_UPDATE_TYPE,                        &
                                       GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE,              &
                                       GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE,        &
                                       GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_REF_TYPE

implicit none
public  :: stream_p05_sieve_cavgs
private
#include "simple_local_flags.inc"
    
type, extends(commander_base) :: stream_p05_sieve_cavgs
    contains
    procedure :: execute => exec_stream_p05_sieve_cavgs
end type stream_p05_sieve_cavgs
    
contains
    
    ! Entry point for stream task 5. Continuously watches dir_target for
    ! completed project files, imports their micrographs and particles into a
    ! growing rec_list, and drives the four-stage chunked2D pipeline
    ! (collect/reject → pass-1 → pass-2 → refchunk → match → submit) on each
    ! loop cycle until a termination signal is detected.
    subroutine exec_stream_p05_sieve_cavgs( self, cline )
        class(stream_p05_sieve_cavgs), intent(inout) :: self
        class(cmdline),                intent(inout) :: cline
        
        ! Hard-coded classification and import parameters
        integer,            parameter   :: MAX_MOVIE_IMPORT = 20   ! max movies imported per cycle
        type(string),       allocatable :: projects(:)
        character(len=:),   allocatable :: meta_buffer
        integer,            allocatable :: jpeg_inds(:), jpeg_pops(:), ref_selection(:)
        real,               allocatable :: jpeg_res(:)
        type(string)                    :: refs_jpeg, refs_stk
        type(string)                    :: match_jpeg, match_stk
        type(rec_list)                  :: project_list
        type(qsys_env)                  :: qenv
        type(parameters)                :: params
        type(sp_project)                :: spproj_glob, spproj_tmp
        type(project_rec)               :: prec
        type(rec_iterator)              :: it
        type(microchunked2D)            :: chunked_2D
        type(stream_watcher)            :: project_buff
        type(gui_metadata_cavg2D)       :: meta_cavg2D
        type(gui_metadata_stream_particle_sieving) :: meta_particle_sieving
        integer :: nprojects, n_mics_imported, n_ptcls_imported, i, xtiles, ytiles
        logical :: l_refs_complete, l_terminate, l_once
        l_once      = .true.
        l_terminate = .false.
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        call cline%set('mkdir',        'yes')
        ! Apply user-overridable defaults
        if( .not. cline%defined('walltime') ) call cline%set('walltime', 29 * 60)
        if( .not. cline%defined('outdir')   ) call cline%set('outdir',        '')
        if( .not. cline%defined('nmics')    ) call cline%set('nmics',        100)
        ! initialise counters
        n_ptcls_imported = 0
        n_mics_imported  = 0
        l_refs_complete  = .false.
        ! initialise metadata
        call meta_particle_sieving%new(GUI_METADATA_STREAM_PARTICLE_SIEVING_TYPE)
        call send_meta(string('initialising'))
        ! Create project file, folder structure, and initialise parameters
        call create_stream_project(spproj_glob, cline, string('sieve_cavgs'))
        call params%new(cline)
        call simple_mkdir(PATH_HERE // DIR_STREAM_COMPLETED)
        ! initialise the queue environment and worker pool
        params%workers     = params%nchunks
        params%worker_nthr = params%nthr * params%nparts
        call init_stream_qenv(params, qenv, string(SIMPLE_STREAM_CHUNK_PARTITION))
        ! Sanity-check: must start from an empty project
        call spproj_glob%read(params%projfile)
        if( spproj_glob%os_mic%get_noris() /= 0 ) &
        THROW_HARD('stream_p05_sieve_cavgs must start from an empty project (e.g. from root project folder)')
        ! wait if dir_target doesn't exist yet
        call send_meta(string('waiting on reference picking'))
        call wait_for_folder2(params%dir_target)
        call wait_for_folder2(params%dir_target//'/spprojs_completed')
        ! Initialise the project watcher 
        project_buff = stream_watcher(LONGTIME, params%dir_target // '/' // DIR_STREAM_COMPLETED, &
        spproj=.true., nretries=10)
        ! Restore previously imported project history to avoid re-importing on restart
        if( file_exists('imported_projects.txt') ) then
            call send_meta(string('importing previous run'))
            call read_filetable(string('imported_projects.txt'), projects)
            if( allocated(projects) ) then
                do i=1, size(projects)
                    call project_buff%add2history(projects(i))
                end do
                deallocate(projects)
            endif
        endif      
        ! Main processing loop — runs until a termination signal is detected
        do
            if( file_exists(TERM_STREAM) .or. l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            ! Detect and import newly completed project files
            call project_buff%watch(nprojects, projects, max_nmovies=MAX_MOVIE_IMPORT)
            if( nprojects > 0 ) then
                call import_new_projects(project_list, projects, n_mics_imported, n_ptcls_imported)
                call project_buff%add2history(projects)
                write(logfhandle,'(A,I6,I9)') &
                '>>> # MICROGRAPHS / PARTICLES IMPORTED : ', n_mics_imported, n_ptcls_imported
            end if
            if( l_once ) then
                if( project_list%size() > 0 ) then
                    it = project_list%begin()
                    call it%get(prec)
                    call spproj_tmp%read_segment('out', prec%projname)
                    call spproj_tmp%get_mskdiam('pickrefs', params%mskdiam)
                    write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', params%mskdiam
                    call spproj_tmp%kill()
                    ! Initialise and drive the chunk pipeline for this first cycle
                    call chunked_2D%new(params, string(PATH_HERE // DIR_STREAM_COMPLETED))
                    call chunked_2D%cycle(project_list)
                    l_once = .false.
                end if 
            else
                ! Drive the chunk pipeline for this cycle
                call chunked_2D%cycle(project_list)
            end if
            if( n_ptcls_imported > 0) then
                call send_meta(string('importing and sieving particles'))
            else
                call send_meta(string('waiting on reference picking'))
            endif
            if( l_refs_complete ) then
                call meta_particle_sieving%set_user_input(.true.)
                if( chunked_2D%get_latest_match(jpeg_inds, jpeg_pops, jpeg_res, match_jpeg, match_stk, xtiles, ytiles) ) then
                    call send_matched_cavgs2D(match_jpeg, size(jpeg_inds))
                end if
            else
                if( chunked_2D%get_references(jpeg_inds, jpeg_pops, jpeg_res, refs_jpeg, refs_stk, xtiles, ytiles) ) then
                    if( chunked_2D%get_reference_selection(ref_selection) ) then
                        do i=1, size(ref_selection)
                            call meta_particle_sieving%set_initial_ref_selection(ref_selection(i))
                        end do
                    end if
                    call send_reference_cavgs2D(refs_jpeg, size(jpeg_inds))
                    l_refs_complete = .true.
                end if
            end if
            call sleep(WAITTIME)
        end do
        call meta_particle_sieving%set_user_input(.false.)
        call send_meta(string('terminating'))
        call chunked_2D%kill()
        call project_buff%kill()
        call project_list%kill()
        call qenv%kill()
        call spproj_glob%kill()
        call meta_cavg2D%kill()
        call meta_particle_sieving%kill()
        
        contains
        
        ! Broadcast progress to the GUI.
        subroutine send_meta( my_stage )
            type(string), intent(in) :: my_stage
            call meta_particle_sieving%set(                         &
            stage              = my_stage,                          &
            particles_imported = n_ptcls_imported,                  &
            particles_accepted = chunked_2D%get_n_accepted_ptcls(), &
            particles_rejected = chunked_2D%get_n_rejected_ptcls())
            if( meta_particle_sieving%assigned() .and. mq_stream_master_in%is_active() ) then
                call meta_particle_sieving%serialise(meta_buffer)
                call mq_stream_master_in%send(meta_buffer)
            endif
        end subroutine send_meta

        ! Serialise and broadcast a single class-average tile's metadata to the GUI.
        ! Precondition: jpeg_inds, jpeg_pops, jpeg_res, xtiles, ytiles must have been
        ! populated by a prior call to get_references() or get_latest_match().
        subroutine send_cavg2D_meta( my_path, my_stk, my_i, my_i_max, my_xtile, my_ytile )
            type(string), intent(in) :: my_path, my_stk
            integer,      intent(in) :: my_i, my_i_max, my_xtile, my_ytile
            integer                  :: my_idx
            real                     :: x_sprite, y_sprite
            my_idx = jpeg_inds(my_i)
            if( xtiles > 1 ) then
                x_sprite = my_xtile * (100.0 / real(xtiles - 1))
            else
                x_sprite = 0.0
            endif
            if( ytiles > 1 ) then
                y_sprite = my_ytile * (100.0 / real(ytiles - 1))
            else
                y_sprite = 0.0
            endif
            call meta_cavg2D%set(                                    &
                path    = my_path,                                   &
                mrcpath = my_stk,                                    &
                i       = my_i,                                      &
                i_max   = my_i_max,                                  &
                res     = jpeg_res(my_i),                            &
                pop     = jpeg_pops(my_i),                           &
                idx     = my_idx,                                    &
                sprite  = sprite_sheet_pos(                          &
                                x = x_sprite,                        &
                                y = y_sprite,                        &
                                h = 100 * ytiles,                    &
                                w = 100 * xtiles))
            if( meta_cavg2D%assigned() .and. mq_stream_master_in%is_active() ) then
                call meta_cavg2D%serialise(meta_buffer)
                call mq_stream_master_in%send(meta_buffer)
            endif
        end subroutine send_cavg2D_meta

        ! Send a batch of cavg2D metadata to the GUI, resetting tile counters.
        ! my_meta_type selects the GUI metadata subtype (ref or match).
        subroutine send_reference_cavgs2D( my_path, n )
            type(string), intent(in) :: my_path
            integer,      intent(in) :: n
            call send_cavgs2D_batch(my_path, refs_stk, n, GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_REF_TYPE)
        end subroutine send_reference_cavgs2D

        subroutine send_matched_cavgs2D( my_path, n )
            type(string), intent(in) :: my_path
            integer,      intent(in) :: n
            call send_cavgs2D_batch(my_path, match_stk, n, GUI_METADATA_STREAM_PARTICLE_SIEVING_CLS2D_TYPE)
        end subroutine send_matched_cavgs2D

        subroutine send_cavgs2D_batch( my_path, my_stk, n, meta_type )
            type(string), intent(in) :: my_path, my_stk
            integer,      intent(in) :: n, meta_type
            integer                  :: my_xtile, my_ytile, my_i
            call meta_cavg2D%kill()
            call meta_cavg2D%new(meta_type)
            my_xtile = 0
            my_ytile = 0
            do my_i = 1, n
                call send_cavg2D_meta(my_path, my_stk, my_i, n, my_xtile, my_ytile)
                my_xtile = my_xtile + 1
                if( my_xtile == xtiles ) then
                    my_xtile = 0
                    my_ytile = my_ytile + 1
                endif
            end do
        end subroutine send_cavgs2D_batch
        
        ! Called asynchronously on SIGTERM. Exits immediately after logging.
        subroutine sigterm_handler()
            write(logfhandle, '(A)') 'SIGTERM RECEIVED'
            l_terminate = .true.
        end subroutine sigterm_handler
        
    end subroutine exec_stream_p05_sieve_cavgs
    
end module simple_stream_p05_sieve_cavgs_new