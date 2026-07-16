!@descr: stream pipeline stage 6 — global 2D classification of pooled particles from sieving
!==============================================================================
! MODULE: simple_stream_p06_pool2D_new
!
! PURPOSE:
!   Drives the continuous global 2D classification loop for the streaming
!   pipeline.  Watches for completed particle-sieve sets, imports them into
!   a growing pool, runs iterative 2D clustering, and broadcasts progress
!   and class-average metadata to the GUI via ipc_pipe_pool2D_in.
!   Also responds to live GUI updates: mask-diameter changes and snapshot-2D
!   write requests received on ipc_pipe_pool2D_out.
!
! ENTRY POINT:
!   stream_p06_pool2D%execute(cline) — called by the stream master
!
! INTERNAL SUBROUTINES:
!   unpause_pool          — clear the pause flag and log resumption
!   import_sets_into_pool — read new sieve sets into the pool; initialise
!                           clustering parameters on first import
!   cleanup4restart       — remove stale files when restarting an existing job
!   send_meta             — broadcast pool-2D progress metadata to the GUI
!   send_meta_snapshot2D  — broadcast snapshot metadata to the GUI
!   send_cavg2D_meta      — serialise and send one class-average sprite tile
!   send_cavgs2D          — iterate all current class averages and send each
!   sigterm_handler       — SIGTERM handler: sets l_terminate for graceful exit
!
! DEPENDENCIES:
!   simple_stream_api, simple_stream2D_state, simple_stream_pool2D_utils,
!   simple_stream_state, simple_gui_metadata_api, unix
!==============================================================================
module simple_stream_p06_pool2D_new
use unix,                        only: SIGTERM, c_write, c_usleep, EAGAIN, EWOULDBLOCK, c_read
use, intrinsic :: iso_c_binding, only: c_char, c_size_t, c_int, c_loc
use simple_stream_api
use simple_stream2D_state,       only: snapshot_iteration, snapshot_selection, snapshot_last_nptcls
use simple_stream_pool2D_utils
use simple_stream_state,         only: ipc_pipe_pool2D_in, ipc_pipe_pool2D_out
use simple_gui_metadata_utils,   only: max_metadata_size
use simple_gui_metadata_api,     only: gui_metadata_cavg2D,                      &
                                       gui_metadata_stream_pool2D,               &
                                       gui_metadata_stream_pool2D_snapshot,       &
                                       gui_metadata_stream_update,               &
                                       sprite_sheet_pos,                         &
                                       GUI_METADATA_STREAM_POOL2D_TYPE,          &
                                       GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE,    &
                                       GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE, &
                                       GUI_METADATA_STREAM_UPDATE_TYPE
implicit none

public :: stream_p06_pool2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p06_pool2D
  contains
    procedure :: execute => exec_stream_p06_pool2D
end type stream_p06_pool2D

contains

    ! Manages Global 2D Clustering
    ! //TODO- handling of un-classified particles
    subroutine exec_stream_p06_pool2D( self, cline )
        class(stream_p06_pool2D), intent(inout) :: self
        class(cmdline),           intent(inout) :: cline
        integer,                  parameter     :: OPTICS_ID_DELTA = 500
        character(len=:),           allocatable   :: meta_buffer
        type(parameters)                          :: params
        type(rec_list)                            :: setslist
        type(stream_watcher)                      :: project_buff
        type(sp_project)                          :: spproj_glob, spproj_tmp
        type(chunk_rec)                           :: crec_mskdiam
        type(rec_iterator)                        :: it_mskdiam
        type(gui_metadata_stream_update)          :: meta_update
        type(gui_metadata_cavg2D)                 :: meta_cavg2D
        type(gui_metadata_stream_pool2D)          :: meta_pool2D
        type(gui_metadata_stream_pool2D_snapshot) :: meta_snapshot
        type(string),               allocatable   :: projects(:)
        character(len=:),           allocatable   :: update_pending
        logical,                    allocatable   :: l_imported(:)
        type(string)                              :: snapshot_filename, snapshot_dir
        integer(kind=dp)                          :: time_last_import
        integer                                   :: i, nprojects, nimported, nptcls_glob, pool_iter, iter_last_import
        integer                                   :: mskdiam_update, extra_pause_iters, last_sent_iter
        integer                                   :: snapshot_id, last_snapshot_id, nptcls_glob_state_1, nmics
        integer                                   :: nptcls_threshold, nptcls_max_threshold, nptcls_dynamic_threshold
        integer                                   :: state_1_particle_rate, optics_id_offset
        integer                                   :: update_expected_len
        logical                                   :: l_pause, l_terminate, l_once, l_changed, l_sieve_final
        real                                      :: final_mskdiam
        update_expected_len   = -1
        l_once               = .true.
        l_terminate          = .false.
        l_sieve_final        = .false. ! set when a sieved set flagged sieve_final=yes is imported
        nptcls_glob          = 0
        nptcls_glob_state_1  = 0
        nmics                = 0
        nptcls_threshold     = 0
        nptcls_max_threshold = 0
        nptcls_dynamic_threshold = 0
        final_mskdiam        = 0.0
        state_1_particle_rate = 0
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        call cline%set('oritype',      'mic')
        call cline%set('mkdir',        'yes')
        call cline%set('autoscale',    'yes')
        call cline%set('reject_mics',  'no')
        call cline%set('refine',       'snhc_smpl')
        call cline%set('ml_reg',       'no')
        call cline%set('objfun',       'euclid')
        call cline%set('sigma_est',    'global')
        call cline%set('cls_init',     'rand')
        call cline%set('numlen',       5)
        if( .not.cline%defined('dynreslim') ) call cline%set('dynreslim', 'yes')
        if( .not.cline%defined('center')    ) call cline%set('center',    'yes')
        if( .not.cline%defined('ncls')      ) call cline%set('ncls',       200)
        if( .not.cline%defined('projfile_optics') ) call cline%set('projfile_optics', OPTICS_JOB_NAME//METADATA_EXT)
        ! restart
        call cleanup4restart
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            call cline%set('projname', 'stream_abinitio2D')
            call cline%set('projfile', 'stream_abinitio2D.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        ! initialise metadata
        call meta_update%new(GUI_METADATA_STREAM_UPDATE_TYPE)
        call meta_pool2D%new(GUI_METADATA_STREAM_POOL2D_TYPE)
        call meta_cavg2D%new(GUI_METADATA_STREAM_POOL2D_CLS2D_TYPE)
        call meta_snapshot%new(GUI_METADATA_STREAM_POOL2D_SNAPSHOT_TYPE)
        call send_meta(string('initialising'))
        call create_stream_project(spproj_glob, cline, string('pool2D'))
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        optics_id_offset = max(params%nicedispid - 1, 0) * OPTICS_ID_DELTA
        write(logfhandle,'(A, I8)')'>>> OPTICS ID OFFSET', optics_id_offset
        ! wait if dir_target doesn't exist yet
        call send_meta(string('waiting on particle sieving'))
        call wait_for_folder2(params%dir_target)
        call wait_for_folder2(params%dir_target//'/spprojs_completed')
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_cluster2D must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! Infinite loop
        last_sent_iter    = 0
        nprojects         = 0        ! # of projects per iteration
        nimported         = 0        ! # of sets per iteration
        nptcls_glob       = 0        ! global # of particles present in the pool
        l_pause           = .false.  ! pause clustering
        pool_iter         = 0
        time_last_import  = huge(time_last_import)
        iter_last_import  = -1       ! Pool iteration # when set(s) last imported
        extra_pause_iters = 0        ! extra iters before pause
        last_snapshot_id  = 0
        do
            if( file_exists(TERM_STREAM) .or. l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            ! detection of new projects
            call project_buff%watch(nprojects, projects)
            if( nprojects > 0 )then
                ! memoize detected projects
                call project_buff%add2history(projects)
                do i = 1,nprojects
                    call setslist%push2chunk_list(projects(i), setslist%size() + 1, .true.)
                enddo
            endif
            if( l_once ) then
                if( setslist%size() > 0 ) then
                    it_mskdiam = setslist%begin()
                    call it_mskdiam%get(crec_mskdiam)
                    call spproj_tmp%read_segment('out', crec_mskdiam%projfile)
                    call spproj_tmp%get_mskdiam('cavg', final_mskdiam)
                    write(logfhandle,'(A,F8.2)') '>>> FINAL MASK DIAMETER SET TO : ', final_mskdiam
                    call spproj_tmp%kill()
                    l_once = .false.
                end if
            end if
            ! check on progress, updates particles & alignment parameters
            ! TODO: class remapping
            if( l_pause )then
                call generate_pool_stats(params)
                ! TODO Flush raw particles
            else
                ! progress
                call update_pool_status
                ! updates particles if iteration completes
                call update_pool(params)
            endif
            ! Import new particles & clustering init
            call import_sets_into_pool( nimported )
            if( nimported > 0 )then
                time_last_import = time8()
                iter_last_import = get_pool_iter()
                if( nptcls_glob_state_1 > nptcls_dynamic_threshold .or. l_sieve_final ) call unpause_pool()
            endif
            l_imported            = setslist%get_included_flags()
            ! Adaptive pause policy:
            ! - iterations 11..20: pause only if no imports for >1 iterations
            ! - iterations 21+:    pause after 1 iteration without imports
            ! When the final sieve set has been imported, run uninterrupted to
            ! iteration 25 (no pausing).
            if( l_sieve_final .and. pool_iter < 25 ) then
                call unpause_pool()
            else if( pool_iter > 1 .and. pool_iter <= 20 ) then
                if( pool_iter > iter_last_import + 1 ) then
                    if( .not.l_pause ) then
                        l_pause = is_pool_available()
                        extra_pause_iters        = 0
                        nptcls_dynamic_threshold = nptcls_glob_state_1 + max(params%ncls * 20, state_1_particle_rate * 50)
                        if( l_pause ) write(logfhandle,'(A,I8)')'>>> PAUSING 2D ANALYSIS UNTIL #PTCLS IS ', nptcls_dynamic_threshold
                    endif
                end if
            else if( pool_iter > 20 ) then
                if( pool_iter > iter_last_import ) then
                    if( .not.l_pause )then
                        l_pause = is_pool_available()
                        extra_pause_iters = 0
                        nptcls_dynamic_threshold = nptcls_glob_state_1 + max(params%ncls * 20, state_1_particle_rate * 500)
                        if( l_pause ) write(logfhandle,'(A,I8)')'>>> PAUSING 2D ANALYSIS UNTIL #PTCLS IS ', nptcls_dynamic_threshold
                    endif
                endif
            endif
            
            ! pause?
            ! if( (pool_iter >= iter_last_import+PAUSE_NITERS+extra_pause_iters) .or.&
            !     & (time8()-time_last_import>PAUSE_TIMELIMIT) )then
            !     if( .not.l_pause )then
            !         l_pause = is_pool_available()
            !         extra_pause_iters = 0
            !         if( l_pause ) write(logfhandle,'(A)')'>>> PAUSING 2D ANALYSIS'
            !     endif
            ! endif
            
            ! Performs clustering iteration
            if( get_pool_iter() == 0 ) then
                nptcls_threshold = max(params%ncls * 20, state_1_particle_rate * 500) ! minimum threshold to trigger clustering; scales with particle rate but has a floor to avoid stalling when rates are very low
                write(logfhandle,'(A,I8)') '>>> INITIAL PARTICLE THRESHOLD: ', nptcls_threshold
            end if
            if( l_pause )then
                ! skip iteration
                call send_meta(string('paused whilst awaiting new particles'))
            else if( nptcls_glob == 0 .or. nptcls_threshold == 0) then
                ! skip iteration but update metadata
                call send_meta(string('waiting for sieved particles'))
            else if( nptcls_glob_state_1 < nptcls_threshold .and. .not. l_sieve_final ) then
                ! skip iteration but update metadata
                call send_meta(string('waiting for minimum number sieved particles ... '// int2str(ceiling(100.0 * float(nptcls_glob_state_1) / float(nptcls_threshold) )) // '%'))
            else
                ! optionally updates alignment parameters
                call update_pool_aln_params
                ! initiates new iteration
                pool_iter = get_pool_iter()
                call iterate_pool(params)
                call send_meta(string('finding and classifying particles'))
            endif
            ! Adaptive mask policy
            pool_iter = get_pool_iter()
            if( pool_iter >= 10 .and. final_mskdiam > 0.0 ) then
               call update_mskdiam(params, nint(final_mskdiam))
               final_mskdiam = 0.0
            end if
            ! http stats
            if(get_pool_iter() > 1) then
                call meta_pool2D%set_user_input(.true.)
                if( last_sent_iter /= last_complete_iter ) then
                    call send_cavgs2D()
                    last_sent_iter = last_complete_iter
                endif
            endif
            ! update params
            if( receive_from_pool2D_out_pipe(meta_buffer) ) then
                if( allocated(meta_buffer) ) then
                    ! deserialise buffer into meta_update
                    meta_update    = transfer(meta_buffer, meta_update)
                    mskdiam_update = nint(meta_update%get_mskdiam2D_update())
                    if( mskdiam_update > 0 .and. mskdiam_update /= nint(params%mskdiam) ) then
                        call update_mskdiam(params, mskdiam_update)
                        if( pool_iter > iter_last_import) extra_pause_iters = PAUSE_NITERS
                        time_last_import = time8()
                        call unpause_pool()
                    endif
                    if( meta_update%get_sieverefs_selection_length() > 0 ) then
                        call update_match_class_states(meta_update%get_sieverefs_selection(), l_changed)
                        if( l_changed ) then
                            if( pool_iter > iter_last_import) extra_pause_iters = PAUSE_NITERS
                            time_last_import = time8()
                            call unpause_pool()
                        end if
                    end if
                    if( meta_update%has_snapshot2D_update() ) then
                        call meta_update%get_snapshot2D_update(snapshot_id, snapshot_iteration, &
                                                               snapshot_selection, snapshot_filename)
                        if( snapshot_id > last_snapshot_id ) then
                            ! build the snapshot directory path once
                            snapshot_dir = string(CWD_GLOB) // '/' // DIR_SNAPSHOT // '/' // &
                                          swap_suffix(snapshot_filename, "", ".simple")
                            call write_project_stream2D(params,                                                                &
                                snapshot_projfile      = snapshot_dir // '/' // snapshot_filename,                             &
                                snapshot_starfile_base = snapshot_dir // '/' // swap_suffix(snapshot_filename, "", ".simple"), &
                                optics_dir             = params%optics_dir,                                                    &
                                optics_offset          = optics_id_offset)
                            last_snapshot_id = snapshot_id
                            call send_meta_snapshot2D()
                        end if
                        if( allocated(snapshot_selection) ) deallocate(snapshot_selection)
                    endif
                    deallocate(meta_buffer)
                endif
            endif
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        if( allocated(l_imported) ) deallocate(l_imported)
        call terminate_stream2D(params, optics_dir=params%optics_dir)
        ! cleanup
        call spproj_glob%kill
        call qsys_cleanup(params)
        ! end gracefully 
        call simple_end('**** SIMPLE_STREAM_ABINITIO2D NORMAL STOP ****')
        contains

            ! resumes with new sets or analysis parameters have been updated
            subroutine unpause_pool
                if( l_pause ) write(logfhandle,'(A)')'>>> RESUMING 2D ANALYSIS'
                l_pause = .false.
            end subroutine unpause_pool

            ! imports new sets of pre-classified particles into the pool
            ! and initialize the clustering module
            subroutine import_sets_into_pool( nimported )
                integer,           intent(out) :: nimported
                type(sp_project),  allocatable :: spprojs(:)
                type(string),      allocatable :: sigmas(:)
                class(sp_project), pointer     :: pool
                type(rec_iterator)             :: it
                type(chunk_rec)                :: crec
                logical, allocatable :: l_processed(:), l_imported(:)
                integer :: nsets2import, iset, nptcls2import, nmics2import, pool_nmics, nptcls
                integer :: i, fromp, imic, ind, iptcl, jptcl, jmic, nptcls_sel, nptcls_sel_tot, ptcl_match_class
                nimported = 0
                if( setslist%size()== 0 ) return
                ! at other times only import when the pool is free
                if( nptcls_glob > 0 .and. .not.is_pool_available() ) return
                l_processed = setslist%get_processed_flags()
                l_imported  = setslist%get_included_flags()
                nsets2import = count(l_processed(:).and.(.not.l_imported(:)))
                if( nsets2import == 0 ) return
                ! read sets in
                allocate(spprojs(setslist%size()))
                nptcls2import = 0
                nmics2import  = 0
                it            = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    call spprojs(iset)%read_non_data_segments(crec%projfile)
                    call spprojs(iset)%read_segment('mic',    crec%projfile)
                    call spprojs(iset)%read_segment('stk',    crec%projfile)
                    call spprojs(iset)%read_segment('ptcl2D', crec%projfile)
                    call spprojs(iset)%read_segment('out',    crec%projfile)
                    ! detect the final sieve set
                    if( spprojs(iset)%os_out%get_noris() >= 1 )then
                        if( spprojs(iset)%os_out%isthere(1, 'sieve_final') )then
                            if( spprojs(iset)%os_out%get_str(1, 'sieve_final') == 'yes' )then
                                l_sieve_final = .true.
                                write(logfhandle,'(A)')'>>> FINAL SIEVE SET DETECTED - RUNNING UNINTERRUPTED TO ITERATION 25'
                            endif
                        endif
                    endif
                    nmics2import  = nmics2import  + spprojs(iset)%os_mic%get_noris()
                    nptcls2import = nptcls2import + spprojs(iset)%os_ptcl2D%get_noris()
                    call it%next()
                enddo
                ! reallocations
                call get_pool_ptr(pool)
                pool_nmics  = pool%os_mic%get_noris()
                nptcls = pool%os_ptcl2D%get_noris()
                if( pool_nmics == 0 )then
                    call pool%os_mic%new(nmics2import, is_ptcl=.false.)
                    call pool%os_stk%new(nmics2import, is_ptcl=.false.)
                    call pool%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
                    fromp = 1
                else
                    call pool%os_mic%reallocate(pool_nmics+nmics2import)
                    call pool%os_stk%reallocate(pool_nmics+nmics2import)
                    call pool%os_ptcl2D%reallocate(nptcls+nptcls2import)
                    fromp = pool%os_stk%get_top(pool_nmics)+1
                endif
                ! parameters transfer
                imic           = pool_nmics
                nptcls_sel_tot = 0
                it             = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    ind = 1
                    do jmic = 1,spprojs(iset)%os_mic%get_noris()
                        imic = imic + 1
                        ! micrograph
                        call pool%os_mic%transfer_ori(imic, spprojs(iset)%os_mic, jmic)
                        ! stack
                        call pool%os_stk%transfer_ori(imic, spprojs(iset)%os_stk, jmic)
                        nptcls = spprojs(iset)%os_stk%get_int(jmic,'nptcls')
                        call pool%os_stk%set(imic, 'fromp', fromp)
                        call pool%os_stk%set(imic, 'top',   fromp+nptcls-1)
                        ! particles
                        !$omp parallel do private(i,iptcl,jptcl,ptcl_match_class) default(shared) proc_bind(close)
                        do i = 1,nptcls
                            iptcl = fromp + i - 1
                            jptcl = ind   + i - 1
                            call pool%os_ptcl2D%transfer_ori(iptcl, spprojs(iset)%os_ptcl2D, jptcl)
                            call pool%os_ptcl2D%set_stkind(iptcl, imic)
                            ! new particle
                            call pool%os_ptcl2D%set(iptcl, 'updatecnt', 0)
                            call pool%os_ptcl2D%set(iptcl, 'frac',      0.)
                            call pool%os_ptcl2D%delete_2Dclustering(iptcl, keepshifts=.true.)
                            ! apply sieverefs-based selection state
                            if( allocated(match_selection) ) then
                                ptcl_match_class = pool%os_ptcl2D%get_int(iptcl, 'class_match')
                                if( .not. any(match_selection == ptcl_match_class) ) then
                                    call pool%os_ptcl2D%set(iptcl, 'state', 0)
                                end if
                            end if
                        enddo
                        !$omp end parallel do
                        ind   = ind   + nptcls
                        fromp = fromp + nptcls
                    enddo
                    nptcls_sel     = spprojs(iset)%os_ptcl2D%get_noris(consider_state=.true.)
                    nptcls_sel_tot = nptcls_sel_tot + nptcls_sel
                    ! display
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',nptcls_sel,' PARTICLES FROM SET ',crec%id
                    call flush(logfhandle)
                    ! global list update
                    crec%included = .true.
                    call setslist%replace_iterator(it, crec)
                    ! move iterator
                    call it%next()
                enddo
                nimported = nsets2import
                ! average all previously imported sigmas
                l_imported  = setslist%get_included_flags()
                allocate(sigmas(count(l_imported)))
                i  = 0
                it = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( .not.l_imported(iset) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    i         = i+1
                    ind       = crec%id
                    call spprojs(iset)%read_segment('out', crec%projfile)
                    call spprojs(iset)%get_sigma2(sigmas(i))
                    ! move iterator
                    call it%next()
                enddo
                call average_sigma2_groups(sigma2_star_from_iter(get_pool_iter()+1), sigmas)
                deallocate(sigmas)
                ! initialize fundamental parameters & clustering
                if( nptcls_glob == 0 )then
                    params%smpd = pool%get_smpd()
                    params%box  = pool%get_box()
                    if(params%mskdiam <= 0.0) then
                        !params%mskdiam = 0.9 * ceiling(params%box * params%smpd)
                        params%mskdiam = 2 * nint(real(params%box/2)-COSMSKHALFWIDTH-1.)
                        call cline%set('mskdiam', params%mskdiam)
                        write(logfhandle,'(A,F8.2)')'>>> INITIAL MASK DIAMETER SET TO', params%mskdiam
                    endif
                    call init_pool_clustering(params, cline, spproj_glob, string(MICSPPROJ_FNAME), reference_generation=.false.)
                endif
                ! global count
                nptcls_glob           = nptcls_glob + nptcls_sel_tot
                nptcls_glob_state_1   = pool%os_ptcl2D%count_state_gt_zero()
                nmics                 = pool%os_mic%get_noris()
                state_1_particle_rate = ceiling(real(nptcls_glob_state_1) / max(1, nmics))
                ! cleanup
                do iset = 1,setslist%size()
                    call spprojs(iset)%kill
                enddo
                if( allocated(l_imported)  ) deallocate(l_imported)
                if( allocated(l_processed) ) deallocate(l_processed)
                deallocate(spprojs)
                nullify(pool)
            end subroutine import_sets_into_pool

            ! Remove previous files from folder to restart
            subroutine cleanup4restart
                type(string) :: cwd_restart, str_dir
                logical      :: l_restart
                call simple_getcwd(cwd_restart)
                l_restart = .false.
                if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
                    l_restart = .true.
                    call simple_chdir(cline%get_carg('outdir'))
                endif
                if(cline%defined('dir_exec')) then
                    if( .not.file_exists(cline%get_carg('dir_exec')) )then
                        str_dir = cline%get_carg('dir_exec')
                        THROW_HARD('Previous directory does not exists: '//str_dir%to_char())
                    endif
                    l_restart = .true.
                endif
                if( l_restart ) then
                    write(logfhandle,'(A,A)') ">>> RESTARTING EXISTING JOB", cwd_restart%to_char()
                    if(cline%defined('dir_exec')) call cline%delete('dir_exec')
                    call del_file(micspproj_fname)
                    call cleanup_root_folder
                endif
                call simple_chdir(cwd_restart)
            end subroutine cleanup4restart

            ! Broadcast progress to the GUI.
            subroutine send_meta( my_stage )
                type(string), intent(in) :: my_stage
                call meta_pool2D%set(                             &
                stage              = my_stage,                    &
                iteration          = last_complete_iter,          &
                particles_imported = nptcls_glob,                 &
                particles_accepted = get_pool_assigned(),         &
                particles_rejected = get_pool_rejected(),         &
                mskdiam            = nint(params%mskdiam),        &
                mskscale           = params%box * params%smpd,    &
                resolution         = get_pool_resolution()        )            
                if( meta_pool2D%assigned() ) then
                    call meta_pool2D%serialise(meta_buffer)
                    call send_to_pool2D_in_pipe(meta_buffer)
                endif
            end subroutine send_meta

            ! Broadcast snapshot metadata to the GUI after writing a classification snapshot.
            subroutine send_meta_snapshot2D()
                call meta_snapshot%set(                                           &
                    id                = last_snapshot_id,                         &
                    snapshot_filename = snapshot_dir // '/' // snapshot_filename, &
                    snapshot_nptcls   = snapshot_last_nptcls)
                if( meta_snapshot%assigned() ) then
                    call meta_snapshot%serialise(meta_buffer)
                    call send_to_pool2D_in_pipe(meta_buffer)
                endif
            end subroutine send_meta_snapshot2D

            ! Serialise one class-average sprite tile and send it to the GUI.
            subroutine send_cavg2D_meta( my_path, my_stk, my_i, my_i_max, my_xtile, my_ytile )
                type(string), intent(in) :: my_path, my_stk
                integer,      intent(in) :: my_i, my_i_max, my_xtile, my_ytile
                integer                  :: my_idx, my_xtiles, my_ytiles
                my_idx    = pool_jpeg_map(my_i)
                my_xtiles = get_pool_cavgs_jpeg_ntilesx()
                my_ytiles = get_pool_cavgs_jpeg_ntilesy()
                call meta_cavg2D%set(                                    &
                    path    = my_path,                                   &
                    mrcpath = my_stk,                                    &
                    i       = my_i,                                      &
                    i_max   = my_i_max,                                  &
                    res     = pool_jpeg_res(my_i),                       &
                    pop     = pool_jpeg_pop(my_i),                       &
                    idx     = my_idx,                                    &
                    sprite  = sprite_sheet_pos(                                                           &
                                    x = my_xtile * (100.0 / max(1, my_xtiles - 1)),                       &
                                    y = my_ytile * (100.0 / max(1, my_ytiles - 1)),                       &
                                    h = 100 * my_ytiles,                                                  &
                                    w = 100 * my_xtiles))
                if( meta_cavg2D%assigned() ) then
                    call meta_cavg2D%serialise(meta_buffer)
                    call send_to_pool2D_in_pipe(meta_buffer)
                endif
            end subroutine send_cavg2D_meta

            ! Send metadata for every current class average to the GUI.
            subroutine send_cavgs2D()
                type(string) :: my_path, my_stk, my_cwd
                integer      :: my_i, my_xtile, my_ytile, my_n, my_xtiles
                call simple_getcwd(my_cwd)
                my_xtiles = get_pool_cavgs_jpeg_ntilesx()
                if(allocated(pool_jpeg_map)) then
                    my_path  = my_cwd // '/' // get_pool_cavgs_jpeg()
                    my_stk   = my_cwd // '/' // get_pool_cavgs_mrc()
                    my_xtile = 0
                    my_ytile = 0
                    my_n     = size(pool_jpeg_map)
                    do my_i = 1, my_n
                        call send_cavg2D_meta(my_path, my_stk, my_i, my_n, my_xtile, my_ytile)
                        my_xtile = my_xtile + 1
                        if( my_xtile == my_xtiles ) then
                            my_xtile = 0
                            my_ytile = my_ytile + 1
                        endif
                    end do
                end if
            end subroutine send_cavgs2D

            subroutine send_to_pool2D_in_pipe(buffer)
                character(len=*), intent(in)                :: buffer
                character(len=:), allocatable               :: framed
                character(kind=c_char), allocatable, target :: cbuf(:)
                integer(c_int)                              :: nwritten
                integer(c_int), target                      :: msg_len
                integer                                     :: err_no
                integer                                     :: sent, nbytes, header_bytes, framed_nbytes, retry_count, ich, rc_sleep
                integer, parameter                          :: MAX_RETRIES = 200
                integer, parameter                          :: RETRY_SLEEP_US = 10000

                if( ipc_pipe_pool2D_in(2) < 0 ) return
                nbytes = len(buffer)
                if( nbytes <= 0 ) return

                msg_len = int(nbytes, c_int)
                header_bytes = sizeof(msg_len)
                framed_nbytes = header_bytes + nbytes
                allocate(character(len=framed_nbytes) :: framed)
                framed(1:header_bytes) = transfer(msg_len, framed(1:header_bytes))
                framed(header_bytes + 1:) = buffer

                allocate(cbuf(framed_nbytes))
                do ich = 1, framed_nbytes
                    cbuf(ich) = transfer(framed(ich:ich), cbuf(ich))
                end do

                sent = 0
                retry_count = 0
                do while( sent < framed_nbytes )
                    nwritten = c_write(ipc_pipe_pool2D_in(2), c_loc(cbuf(sent + 1)), int(framed_nbytes - sent, c_size_t))
                    if( nwritten > 0 ) then
                        sent = sent + int(nwritten)
                        retry_count = 0
                        cycle
                    end if

                    err_no = ierrno()
                    if( err_no == int(EAGAIN) .or. err_no == int(EWOULDBLOCK) ) then
                        retry_count = retry_count + 1
                        if( retry_count > MAX_RETRIES ) then
                            THROW_HARD('failed to write pool2D metadata to ipc_pipe_pool2D_in: retry limit exceeded')
                        end if
                        rc_sleep = c_usleep(RETRY_SLEEP_US)
                        cycle
                    end if

                    THROW_HARD('failed to write pool2D metadata to ipc_pipe_pool2D_in')
                end do

                if( allocated(cbuf) ) deallocate(cbuf)
            end subroutine send_to_pool2D_in_pipe

            logical function receive_from_pool2D_out_pipe(buffer)
                character(len=:), allocatable, intent(inout) :: buffer
                character(kind=c_char), allocatable, target  :: raw(:)
                character(len=:), allocatable                :: chunk
                integer(c_int), target                       :: msg_len_c
                integer                                      :: header_bytes
                integer                                      :: nread, ibyte, err_no, max_meta_bytes

                receive_from_pool2D_out_pipe = .false.
                if( allocated(buffer) ) deallocate(buffer)
                header_bytes = sizeof(msg_len_c)
                max_meta_bytes = max_metadata_size()
                allocate(raw(max_meta_bytes))

                nread = c_read(ipc_pipe_pool2D_out(1), c_loc(raw(1)), int(size(raw), c_size_t))
                if( nread < 0 ) then
                    err_no = ierrno()
                    if( err_no == int(EAGAIN) .or. err_no == int(EWOULDBLOCK) ) return
                    return
                endif
                if( nread > 0 ) then
                    allocate(character(len=nread) :: chunk)
                    do ibyte = 1, nread
                        chunk(ibyte:ibyte) = transfer(raw(ibyte), 'a')
                    end do
                    if( allocated(update_pending) ) then
                        update_pending = update_pending // chunk
                    else
                        allocate(character(len=nread) :: update_pending)
                        update_pending = chunk
                    endif
                endif

                if( update_expected_len < 0 ) then
                    if( .not.allocated(update_pending) ) return
                    if( len(update_pending) < header_bytes ) return
                    msg_len_c = transfer(update_pending(1:header_bytes), msg_len_c)
                    update_expected_len = int(msg_len_c)
                    if( update_expected_len <= 0 .or. update_expected_len > max_meta_bytes ) then
                        THROW_HARD('invalid framed metadata length read from pool2D_out pipe')
                    endif
                    if( len(update_pending) == header_bytes ) then
                        deallocate(update_pending)
                    else
                        update_pending = update_pending(header_bytes + 1:)
                    endif
                endif

                if( .not.allocated(update_pending) ) return
                if( len(update_pending) < update_expected_len ) return

                allocate(character(len=update_expected_len) :: buffer)
                buffer = update_pending(1:update_expected_len)
                if( len(update_pending) == update_expected_len ) then
                    deallocate(update_pending)
                else
                    update_pending = update_pending(update_expected_len + 1:)
                endif
                update_expected_len = -1
                receive_from_pool2D_out_pipe = .true.
            end function receive_from_pool2D_out_pipe

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p06_pool2D

end module simple_stream_p06_pool2D_new
