!@descr: stream pipeline stage 7 — multistate 3D reconstruction/refinement of pooled particles
!==============================================================================
! MODULE: simple_stream_p07_abinitio3D_multistate
!
! PURPOSE:
!   Drives the continuous multistate 3D reconstruction/refinement loop for the
!   streaming pipeline. Watches for completed 3D-export sets (written by
!   stream_p06_pool2D) and imports them into a growing pool project.
!
!   Broadcasts stage/progress metadata and per-state gui_metadata_vol3D
!   (volume paths, population, FSC-derived resolution/curve, orientation
!   distribution histogram) to the GUI via ipc_pipe_abinitio3D_multstate_in.
!   refine3D is not yet wired in — only the abinitio3D stage currently
!   produces per-state volumes/FSCs.
!
! ENTRY POINT:
!   stream_p07_abinitio3D_multistate%execute(cline) — called by the stream master
!
! INTERNAL SUBROUTINES:
!   import_sets_into_pool          — read new exported sets into the pool
!   send_meta_abinitio3D_multistate — broadcast stage/progress metadata to the GUI
!   build_and_send_vol3D_states     — build/send per-state gui_metadata_vol3D once
!                                      abinitio3D volumes/FSCs are available
!   compute_oridist_for_state       — bin one state's particle orientations into
!                                      the 72x36 azimuth/elevation histogram
!   locate_state_jpeg               — locate a per-state output jpeg (reprojections
!                                      or orientation-distribution heatmap) alongside its volume
!   send_state_reprojtiles           — send one gui_metadata_cavg2D entry per
!                                      orthogonal reprojection tile in a state's sprite sheet
!   send_to_abinitio3D_multstate_in_pipe — frame and write a metadata buffer
!   sigterm_handler       — SIGTERM handler: sets l_terminate for graceful exit
!
! DEPENDENCIES:
!   simple_stream_api, simple_stream_state, simple_gui_metadata_api,
!   simple_refine3D_fnames, unix
!==============================================================================
module simple_stream_p07_abinitio3D_multistate
use unix,                        only: SIGTERM, c_write, c_usleep, EAGAIN, EWOULDBLOCK, EINTR
use, intrinsic :: iso_c_binding, only: c_char, c_size_t, c_int, c_loc
use simple_commanders_cavgs,     only: commander_model_cavgs_rejection
use simple_gui_utils,            only: mrc2jpeg_tiled
use simple_imghead,              only: get_mrc_minmax
use simple_qsys_env,             only: qsys_env
use simple_refine3D_fnames,      only: refine3D_oris_heatmap_fname
use simple_stream_state,         only: ipc_pipe_abinitio3D_multstate_in
use simple_gui_metadata_api,     only: gui_metadata_stream_abinitio3D_multistate, gui_metadata_vol3D, &
                                       gui_metadata_cavg2D, sprite_sheet_pos,                          &
                                       GUI_METADATA_STREAM_ABINITIO3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE, &
                                       GUI_METADATA_STREAM_ABINITIO3D_MULTISTATE_REPROJ_TYPE
use simple_stream_api
implicit none

public :: stream_p07_abinitio3D_multistate
private
#include "simple_local_flags.inc"

integer, parameter       :: NSTATES3D  = 3                 ! number of classes for abinitio3D
integer, parameter       :: NSTAGES3D  = 1!5                 ! number of stages for abinitio3D

type, extends(commander_base) :: stream_p07_abinitio3D_multistate
  contains
    procedure :: execute => exec_stream_p07_abinitio3D_multistate
end type stream_p07_abinitio3D_multistate

contains

    ! Manages multistate 3D reconstruction/refinement
    subroutine exec_stream_p07_abinitio3D_multistate( self, cline )
        class(stream_p07_abinitio3D_multistate), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)          :: params
        type(rec_list)            :: setslist
        type(stream_watcher)      :: project_buff
        type(sp_project)          :: spproj_glob
        type(qsys_env)            :: qenv
        type(string), allocatable :: projects(:)
        type(gui_metadata_stream_abinitio3D_multistate) :: meta_abinitio3D_multistate
        type(gui_metadata_vol3D), allocatable            :: meta_states_vol3D(:)
        type(gui_metadata_cavg2D), allocatable            :: meta_reprojtiles(:)
        character(len=:),         allocatable            :: meta_buffer
        integer                   :: i, nprojects, nimported, nptcls_glob, abinitio_stage, refine_stage
        integer                   :: envlen, refine_it, nptcls_at_last_refine
        character(len=STDLEN)     :: preproc_part_env
        logical                   :: l_terminate, l_pause_ingestion
        volatile :: l_terminate ! set asynchronously by sigterm_handler
        l_terminate           = .false.
        abinitio_stage        = 0
        refine_stage          = 0
        nptcls_glob           = 0
        refine_it             = 0
        nptcls_at_last_refine = 0
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        call cline%set('oritype', 'mic')
        call cline%set('mkdir',   'yes')
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            call cline%set('projname', 'stream_3Dmultistate')
            call cline%set('projfile', 'stream_3Dmultistate.simple')
            call spproj_glob%update_projinfo(cline)
            call spproj_glob%update_compenv(cline)
            call spproj_glob%write
        endif
        call create_stream_project(spproj_glob, cline, string('3Dmultistate'))
        ! master parameters
        call params%new(cline)
        call cline%set('mkdir', 'no')
        ! GUI metadata
        call meta_abinitio3D_multistate%new(GUI_METADATA_STREAM_ABINITIO3D_MULTISTATE_TYPE)
        ! setup the environment for distributed execution
        call get_environment_variable(SIMPLE_STREAM_PREPROC_PARTITION, preproc_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(params, 1,stream=.true.,qsys_partition=string(trim(preproc_part_env)))
        else
            call qenv%new(params, 1,stream=.true.)
        end if
        ! wait if dir_target doesn't exist yet
        call wait_for_folder2(params%dir_target)
        call wait_for_folder2(params%dir_target//'/spprojs_completed')
        ! master project file
        call spproj_glob%read( params%projfile )
        if( spproj_glob%os_mic%get_noris() /= 0 ) THROW_HARD('stream_abinitio3D_multistate must start from an empty project (eg from root project folder)')
        ! project watcher
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! Infinite loop
        nprojects = 0 ! # of projects per iteration
        nimported = 0 ! # of sets per iteration
        l_pause_ingestion = .false.
        do
            if( file_exists(TERM_STREAM) .or. l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            ! detection of new projects
            call project_buff%watch(nprojects, projects)
            if( nprojects > 0 )then
                ! zero-padded names sort lexicographically in export order
                call lex_sort(projects)
                ! memoize detected projects
                call project_buff%add2history(projects)
                do i = 1,nprojects
                    call setslist%push2chunk_list(projects(i), setslist%size() + 1, .true.)
                enddo
            endif
            ! Import new particles, paused while abinitio3D or refine3D is running
            if( .not. l_pause_ingestion ) call import_sets_into_pool( nimported )
            ! abinitio stage
            if( abinitio_stage < 2 .and. spproj_glob%os_ptcl2D%get_noris() /= 0 ) then
                if( abinitio_stage == 0 ) then
                    l_pause_ingestion = .true.
                    ! start abinitio 3D
                    call start_abinitio3D(spproj_glob, string('abinitio3D'), 5000)
                    abinitio_stage = 1
                else if( abinitio_stage == 1 ) then
                  ! Test for abinitio 3D completion
                  if( file_exists(string('abinitio3D')//'/'//TASK_FINISHED) ) then
                      ! stage complete
                      call finish_abinitio3D(spproj_glob, string('abinitio3D'))
                      abinitio_stage    = 2
                      l_pause_ingestion = .false.
                      call build_and_send_vol3D_states
                  end if
                end if
            end if
            ! refine stage
            if( abinitio_stage == 2 .and. spproj_glob%os_ptcl2D%get_noris() /= 0 ) then
                if( refine_stage == 0 ) then
                    ! only enter if the particle count has grown since the last refine stage
                    if( spproj_glob%os_ptcl2D%get_noris() > nptcls_at_last_refine ) then
                        refine_it             = refine_it + 1
                        nptcls_at_last_refine = spproj_glob%os_ptcl2D%get_noris()
                        write(logfhandle,'(A,I0)')'>>> ENTERING REFINE STAGE ', refine_it
                        l_pause_ingestion = .true.
                        ! start refine 3D
                     !   call start_refine3D(spproj_glob, string('refine3D/it_')//int2str(refine_it), 5000)
                        call spproj_glob%write(string('refine_') // int2str(refine_it) // METADATA_EXT)
                        refine_stage = 1
                    end if
                else if( refine_stage == 1 ) then
                    ! Test for refine stage completion
                !    if( file_exists(string('refine3D/it_')//int2str(refine_it)//'/'//TASK_FINISHED) ) then
                !        ! stage complete
                !        call finish_refine3D(spproj_glob, string('refine3D')//int2str(refine_it))
                !        refine_stage = 0
                !    end if
                    l_pause_ingestion = .false.
                    refine_stage = 0 ! for testing purposes
                end if
            end if
            ! broadcast progress to the GUI
            call send_meta_abinitio3D_multistate
            ! Wait
            call sleep(WAITTIME)
        enddo
        ! Cleanup and final project
        call spproj_glob%kill
        call qsys_cleanup(params)
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_ABINITIO3D_MULTISTATE NORMAL STOP ****')
        contains

            ! imports new sets of exported particles into the pool
            subroutine import_sets_into_pool( nimported )
                integer,           intent(out) :: nimported
                type(sp_project),  allocatable :: spprojs(:)
                type(rec_iterator)             :: it
                type(chunk_rec)                :: crec
                logical, allocatable :: l_processed(:), l_imported(:), l_include_now(:)
                integer :: nsets2import, iset, nptcls2import, nmics2import, pool_nmics, nptcls
                integer :: i, fromp, imic, ind, iptcl, jptcl, jmic, nptcls_sel_tot, first_iset
                type(cmdline)                         :: cline_quality
                type(commander_model_cavgs_rejection) :: xmodel_cavgs_rejection
                type(string)                          :: quality_dir, quality_projfile, cwd_before, selected_jpeg, rejected_jpeg
                nimported = 0
                if( setslist%size() == 0 ) return
                l_processed = setslist%get_processed_flags()
                l_imported  = setslist%get_included_flags()
                nsets2import = count(l_processed(:).and.(.not.l_imported(:)))
                if( nsets2import == 0 ) return
                allocate(l_include_now(setslist%size()), source=.false.)
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
                    quality_dir      = string('quality_selection/') // get_fbody(basename(crec%projfile), 'simple')
                    quality_projfile = quality_dir // '/' // basename(crec%projfile)
                    call simple_mkdir('quality_selection')
                    call simple_mkdir(quality_dir)
                    call simple_copy_file(crec%projfile, quality_projfile)
                    ! class-average quality selection, applied to the copy before import
                    call cline_quality%set('projfile',       basename(crec%projfile))
                    call cline_quality%set('mskdiam',        params%mskdiam)
                    call cline_quality%set('quality_mode',   'apply')
                    call cline_quality%set('rejection_type', 'pool')
                    call cline_quality%set('mkdir',          'no')
                    call simple_getcwd(cwd_before)
                    call simple_chdir(quality_dir)
                    call xmodel_cavgs_rejection%execute(cline_quality)
                    call mrc2jpeg_tiled(string('quality_selected_cavgs')//MRC_EXT, string('quality_selected_cavgs')//JPG_EXT)
                    call mrc2jpeg_tiled(string('quality_rejected_cavgs')//MRC_EXT, string('quality_rejected_cavgs')//JPG_EXT)
                    selected_jpeg = cwd_before//'/'//quality_dir//'/quality_selected_cavgs'//JPG_EXT
                    if( file_exists(selected_jpeg) ) then
                        write(logfhandle, '(A)')   '>>> QUALITY SELECTED CLASS AVERAGES'
                        write(logfhandle, '(A,A)') '>>> JPEG ', selected_jpeg%to_char()
                    endif
                    rejected_jpeg = cwd_before//'/'//quality_dir//'/quality_rejected_cavgs'//JPG_EXT
                    if( file_exists(rejected_jpeg) ) then
                        write(logfhandle, '(A)')   '>>> QUALITY REJECTED CLASS AVERAGES'
                        write(logfhandle, '(A,A)') '>>> JPEG ', rejected_jpeg%to_char()
                    endif
                    call simple_chdir(cwd_before)
                    ! read the entire project file for the current set
                    call spprojs(iset)%read(quality_projfile)
                    nmics2import  = nmics2import  + spprojs(iset)%os_mic%get_noris()
                    nptcls2import = nptcls2import + spprojs(iset)%os_ptcl3D%get_noris()
                    l_include_now(iset) = .true.
                    call it%next()
                enddo
                ! reallocations
                pool_nmics = spproj_glob%os_mic%get_noris()
                nptcls     = spproj_glob%os_ptcl3D%get_noris()
                if( pool_nmics == 0 )then
                    call spproj_glob%os_mic%new(nmics2import, is_ptcl=.false.)
                    call spproj_glob%os_stk%new(nmics2import, is_ptcl=.false.)
                    call spproj_glob%os_ptcl2D%new(nptcls2import, is_ptcl=.true.)
                    call spproj_glob%os_ptcl3D%new(nptcls2import, is_ptcl=.true.)
                    ! cls2D is only transferred once, to seed the pool project on the initial import
                    first_iset = findloc(l_include_now, .true., dim=1)
                    if( first_iset > 0 ) spproj_glob%os_cls2D = spprojs(first_iset)%os_cls2D
                    fromp = 1
                else
                    call spproj_glob%os_mic%reallocate(pool_nmics+nmics2import)
                    call spproj_glob%os_stk%reallocate(pool_nmics+nmics2import)
                    call spproj_glob%os_ptcl2D%reallocate(nptcls+nptcls2import)
                    call spproj_glob%os_ptcl3D%reallocate(nptcls+nptcls2import)
                    fromp = spproj_glob%os_stk%get_top(pool_nmics)+1
                endif
                ! parameters transfer
                imic           = pool_nmics
                nptcls_sel_tot = 0
                it             = setslist%begin()
                do iset = 1,setslist%size()
                    call it%get(crec)
                    if( crec%included .or. (.not.crec%processed .or. crec%busy) .or. .not.l_include_now(iset) )then
                        ! move iterator
                        call it%next()
                        cycle
                    endif
                    ind = 1
                    do jmic = 1,spprojs(iset)%os_mic%get_noris()
                        imic = imic + 1
                        ! micrograph
                        call spproj_glob%os_mic%transfer_ori(imic, spprojs(iset)%os_mic, jmic)
                        ! stack
                        call spproj_glob%os_stk%transfer_ori(imic, spprojs(iset)%os_stk, jmic)
                        nptcls = spprojs(iset)%os_stk%get_int(jmic,'nptcls')
                        call spproj_glob%os_stk%set(imic, 'fromp', fromp)
                        call spproj_glob%os_stk%set(imic, 'top',   fromp+nptcls-1)
                        ! particles
                        do i = 1,nptcls
                            iptcl = fromp + i - 1
                            jptcl = ind   + i - 1
                            call spproj_glob%os_ptcl2D%transfer_ori(iptcl, spprojs(iset)%os_ptcl2D, jptcl)
                            call spproj_glob%os_ptcl2D%set_stkind(iptcl, imic)
                            call spproj_glob%os_ptcl3D%transfer_ori(iptcl, spprojs(iset)%os_ptcl3D, jptcl)
                            call spproj_glob%os_ptcl3D%set_stkind(iptcl, imic)
                        enddo
                        ind   = ind   + nptcls
                        fromp = fromp + nptcls
                    enddo
                    nptcls_sel_tot = nptcls_sel_tot + spprojs(iset)%os_ptcl3D%get_noris(consider_state=.true.)
                    ! display
                    write(logfhandle,'(A,I6,A,I6)')'>>> TRANSFERRED ',spprojs(iset)%os_ptcl3D%get_noris(consider_state=.true.),' PARTICLES FROM SET ',crec%id
                    call flush(logfhandle)
                    ! global list update
                    crec%included = .true.
                    call setslist%replace_iterator(it, crec)
                    ! move iterator
                    call it%next()
                enddo
                nimported   = count(l_include_now)
                nptcls_glob = nptcls_glob + nptcls_sel_tot
                ! cleanup
                do iset = 1,setslist%size()
                    call spprojs(iset)%kill
                enddo
                if( allocated(l_imported)     ) deallocate(l_imported)
                if( allocated(l_processed)    ) deallocate(l_processed)
                if( allocated(l_include_now)  ) deallocate(l_include_now)
                deallocate(spprojs)
            end subroutine import_sets_into_pool

            ! Broadcast pipeline stage/progress and per-state population/resolution to the GUI.
            subroutine send_meta_abinitio3D_multistate
                type(string)      :: my_stage, fsc_fname
                integer           :: istate, my_pop, fsc_box
                real              :: my_res, res0143, res05
                real, allocatable :: fsc_arr(:), res_arr(:)
                if( abinitio_stage == 0 ) then
                    my_stage = string('importing particles')
                else if( abinitio_stage == 1 ) then
                    my_stage = string('running abinitio3D')
                else if( refine_stage == 1 ) then
                    my_stage = string('running refine3D')
                else
                    my_stage = string('idle')
                endif
                call meta_abinitio3D_multistate%set(                                  &
                    stage                    = my_stage,                              &
                    abinitio3D_stage         = abinitio_stage,                        &
                    refine_iteration         = refine_it,                             &
                    nstates                  = NSTATES3D,                             &
                    particles_imported       = spproj_glob%os_ptcl3D%get_noris(),     &
                    particles_at_last_refine = nptcls_at_last_refine,                 &
                    resolution               = 0.0)
                if( abinitio_stage == 2 ) then
                    do istate = 1, NSTATES3D
                        my_pop = spproj_glob%os_ptcl3D%get_pop(istate, 'state')
                        my_res = 0.0
                        if( spproj_glob%isthere_in_osout('fsc', istate) ) then
                            call spproj_glob%get_fsc(istate, fsc_fname, fsc_box)
                            if( fsc_fname%strlen() > 0 ) then
                                if( file_exists(fsc_fname) ) then
                                    fsc_arr = file2rarr(fsc_fname)
                                    res_arr = get_resarr(fsc_box, spproj_glob%get_smpd())
                                    call get_resolution(fsc_arr, res_arr, res05, res0143)
                                    my_res = res0143
                                endif
                            endif
                        endif
                        call meta_abinitio3D_multistate%set_state_stats(istate, my_pop, my_res)
                    enddo
                endif
                if( meta_abinitio3D_multistate%assigned() ) then
                    call meta_abinitio3D_multistate%serialise(meta_buffer)
                    call send_to_abinitio3D_multstate_in_pipe(meta_buffer)
                endif
            end subroutine send_meta_abinitio3D_multistate

            ! Build (once volumes/FSCs are available) and send one gui_metadata_vol3D
            ! entry per state: paths, population, FSC-derived resolution/curve, and
            ! the orientation-distribution histogram. Called once abinitio3D completes;
            ! refine3D is not yet wired in, so this reflects the abinitio3D output only.
            subroutine build_and_send_vol3D_states
                integer                      :: istate, my_pop, my_box, n_fsc_pts, k, fsc_box
                real                         :: my_smpd, res0143, res05
                real                         :: minval3D, maxval3D
                type(string)                 :: volpath, fsc_fname, pprocpath, lppath, pprocmirrpath, reprojpath, oridistpath
                real,          allocatable   :: fsc_arr(:), res_arr(:), invres_arr(:)
                integer                      :: hist(72, 36) ! matches gui_metadata_vol3D's ORIDIST_NBINS_X x ORIDIST_NBINS_Y (5-degree bins)
                logical                      :: l_have_fsc
                if( .not.allocated(meta_states_vol3D) ) then
                    allocate(meta_states_vol3D(NSTATES3D))
                    do istate = 1, NSTATES3D
                        call meta_states_vol3D(istate)%new(GUI_METADATA_VOL3D_TYPE)
                    enddo
                endif
                do istate = 1, NSTATES3D
                    if( .not.spproj_glob%isthere_in_osout('vol', istate) ) cycle
                    call spproj_glob%get_vol('vol', istate, volpath, my_smpd, my_box)
                    if( volpath%strlen() == 0 ) cycle
                    my_pop = spproj_glob%os_ptcl3D%get_pop(istate, 'state')
                    ! optional postprocessed/low-pass/mirrored products, if already present on disk
                    ! (no postprocessing step runs in this stage yet, so these are typically absent)
                    pprocpath = add2fbody(volpath, MRC_EXT, PPROC_SUFFIX)
                    if( .not.file_exists(pprocpath) ) pprocpath = string('')
                    lppath = add2fbody(volpath, MRC_EXT, LP_SUFFIX)
                    if( .not.file_exists(lppath) ) lppath = string('')
                    if( pprocpath%strlen() > 0 ) then
                        pprocmirrpath = add2fbody(pprocpath, MRC_EXT, MIRR_SUFFIX)
                        if( .not.file_exists(pprocmirrpath) ) pprocmirrpath = string('')
                    else
                        pprocmirrpath = string('')
                    endif
                    call locate_state_jpeg(volpath, string('orthogonal_reprojs_state')//int2str_pad(istate,2)//JPG_EXT, reprojpath)
                    call locate_state_jpeg(volpath, refine3D_oris_heatmap_fname(istate), oridistpath)
                    if( reprojpath%strlen() > 0 ) call send_state_reprojtiles(istate, reprojpath, volpath, my_pop)
                    ! FSC curve + resolution
                    l_have_fsc = .false.
                    res0143    = 0.0
                    res05      = 0.0
                    if( spproj_glob%isthere_in_osout('fsc', istate) ) then
                        call spproj_glob%get_fsc(istate, fsc_fname, fsc_box)
                        if( fsc_fname%strlen() > 0 ) then
                            if( file_exists(fsc_fname) ) then
                                fsc_arr    = file2rarr(fsc_fname)
                                res_arr    = get_resarr(fsc_box, my_smpd)
                                call get_resolution(fsc_arr, res_arr, res05, res0143)
                                l_have_fsc = .true.
                            endif
                        endif
                    endif
                    if( l_have_fsc ) then
                        call meta_states_vol3D(istate)%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, &
                            &istate, my_box, my_smpd, istate, NSTATES3D, res0143=res0143, res05=res05, pop=my_pop, &
                            &oridistpath=oridistpath)
                    else
                        call meta_states_vol3D(istate)%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, &
                            &istate, my_box, my_smpd, istate, NSTATES3D, pop=my_pop, oridistpath=oridistpath)
                    endif
                    ! MRC header min/max, read once here so GUI consumers don't need to
                    ! reopen each volume file per request
                    call get_mrc_minmax(volpath, minval3D, maxval3D)
                    call meta_states_vol3D(istate)%set_minmax('volpath', minval3D, maxval3D)
                    if( pprocpath%strlen() > 0 ) then
                        call get_mrc_minmax(pprocpath, minval3D, maxval3D)
                        call meta_states_vol3D(istate)%set_minmax('pprocpath', minval3D, maxval3D)
                    end if
                    if( lppath%strlen() > 0 ) then
                        call get_mrc_minmax(lppath, minval3D, maxval3D)
                        call meta_states_vol3D(istate)%set_minmax('lppath', minval3D, maxval3D)
                    end if
                    if( pprocmirrpath%strlen() > 0 ) then
                        call get_mrc_minmax(pprocmirrpath, minval3D, maxval3D)
                        call meta_states_vol3D(istate)%set_minmax('pprocmirrpath', minval3D, maxval3D)
                    end if
                    if( l_have_fsc ) then
                        n_fsc_pts = min(size(fsc_arr), 1000) ! matches gui_metadata_vol3D's MAX_FSC_VOL3D
                        allocate(invres_arr(n_fsc_pts))
                        do k = 1, n_fsc_pts
                            invres_arr(k) = 1.0 / res_arr(k)
                        enddo
                        call meta_states_vol3D(istate)%set_fsc(invres_arr(1:n_fsc_pts), fsc_arr(1:n_fsc_pts))
                        deallocate(invres_arr)
                    endif
                    call compute_oridist_for_state(istate, hist)
                    call meta_states_vol3D(istate)%set_oridist(hist)
                    if( meta_states_vol3D(istate)%assigned() ) then
                        call meta_states_vol3D(istate)%serialise(meta_buffer)
                        call send_to_abinitio3D_multstate_in_pipe(meta_buffer)
                    endif
                enddo
            end subroutine build_and_send_vol3D_states

            ! Bin one state's particle orientations (os_ptcl3D projection directions)
            ! into a 72x36 azimuth (-180..180) x elevation (-90..90) histogram, 5 degree bins.
            subroutine compute_oridist_for_state( istate, hist )
                integer, intent(in)  :: istate
                integer, intent(out) :: hist(72, 36)
                real    :: normal(3), azimuth, elevation
                integer :: iptcl, nptcls, ix, iy
                hist   = 0
                nptcls = spproj_glob%os_ptcl3D%get_noris()
                do iptcl = 1, nptcls
                    if( spproj_glob%os_ptcl3D%get_state(iptcl) /= istate ) cycle
                    normal    = spproj_glob%os_ptcl3D%get_normal(iptcl)
                    azimuth   = rad2deg(atan2(normal(2), normal(1)))
                    elevation = rad2deg(asin(max(-1.0, min(1.0, normal(3)))))
                    ix = min(72, max(1, floor((azimuth   + 180.0) / 5.0) + 1))
                    iy = min(36, max(1, floor((elevation + 90.0)  / 5.0) + 1))
                    hist(ix, iy) = hist(ix, iy) + 1
                enddo
            end subroutine compute_oridist_for_state

            ! Locate a per-state output jpeg already produced alongside the volume
            ! (by abinitio3D's calc_final_rec/gen_ortho_reprojs4viz and refine3D's
            ! orientation-distribution heatmap); '' if not present on disk.
            subroutine locate_state_jpeg( volpath, fname, jpegpath )
                type(string), intent(in)  :: volpath, fname
                type(string), intent(out) :: jpegpath
                jpegpath = get_fpath(volpath) // fname
                if( .not. file_exists(jpegpath) ) jpegpath = string('')
            end subroutine locate_state_jpeg

            ! Send the 3 orthogonal reprojection tiles of one state's sprite-sheet
            ! jpeg as individual gui_metadata_cavg2D entries (idx=state, sprite
            ! position selects the tile); flat-indexed i/i_max across all states.
            subroutine send_state_reprojtiles( istate, reprojpath, volpath, pop )
                integer,      intent(in) :: istate, pop
                type(string), intent(in) :: reprojpath, volpath
                integer, parameter :: NTILES = 3
                integer            :: itile, i_flat
                if( .not.allocated(meta_reprojtiles) ) then
                    allocate(meta_reprojtiles(NSTATES3D * NTILES))
                    do i_flat = 1, size(meta_reprojtiles)
                        call meta_reprojtiles(i_flat)%new(GUI_METADATA_STREAM_ABINITIO3D_MULTISTATE_REPROJ_TYPE)
                    enddo
                endif
                do itile = 1, NTILES
                    i_flat = (istate - 1) * NTILES + itile
                    call meta_reprojtiles(i_flat)%set(path=reprojpath, mrcpath=volpath, idx=istate, &
                        &sprite=sprite_sheet_pos(x=real(itile-1)*(100.0/real(NTILES-1)), y=0.0, h=100, w=100*NTILES), &
                        &i=i_flat, i_max=size(meta_reprojtiles), pop=pop)
                    call meta_reprojtiles(i_flat)%serialise(meta_buffer)
                    call send_to_abinitio3D_multstate_in_pipe(meta_buffer)
                enddo
            end subroutine send_state_reprojtiles

            ! Frame (length-prefix) and write a serialised metadata buffer to the
            ! master-facing abinitio3D_multistate IPC pipe, with EAGAIN/EINTR retry.
            subroutine send_to_abinitio3D_multstate_in_pipe(buffer)
                character(len=*), intent(in)                :: buffer
                character(len=:), allocatable               :: framed
                character(kind=c_char), allocatable, target :: cbuf(:)
                integer(c_int)                              :: nwritten
                integer(c_int), target                      :: msg_len
                integer                                     :: err_no
                integer                                     :: sent, nbytes, header_bytes, framed_nbytes, retry_count, ich, rc_sleep
                integer, parameter                          :: MAX_RETRIES    = 3000
                integer, parameter                          :: RETRY_SLEEP_US = 10000

                if( ipc_pipe_abinitio3D_multstate_in(2) < 0 ) return
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
                    nwritten = int(c_write(ipc_pipe_abinitio3D_multstate_in(2), c_loc(cbuf(sent + 1)), &
                        &int(framed_nbytes - sent, c_size_t)))
                    if( nwritten > 0 ) then
                        sent = sent + int(nwritten)
                        retry_count = 0
                        cycle
                    end if

                    err_no = ierrno()
                    if( err_no == int(EINTR) ) then
                        ! interrupted system call; not backpressure, just retry immediately
                        ! without counting against the retry budget or touching sent
                        cycle
                    end if

                    if( err_no == int(EAGAIN) .or. err_no == int(EWOULDBLOCK) ) then
                        retry_count = retry_count + 1
                        if( retry_count > MAX_RETRIES ) then
                            ! Bail out unconditionally once the retry budget is exhausted,
                            ! even if part of the frame already reached the pipe. Blocking
                            ! forever on a stalled/dead reader would silently hang the whole
                            ! polling loop (no further metadata, no reprojections, no log
                            ! output). A partial frame may desync the reader's length-prefixed
                            ! framing, but that is preferable to an indefinite hang; the next
                            ! reconnect/restart of the reader will resynchronise.
                            THROW_WARN('failed to write abinitio3D_multistate metadata to ipc_pipe_abinitio3D_multstate_in: retry limit exceeded')
                            exit
                        end if
                        rc_sleep = c_usleep(RETRY_SLEEP_US)
                        cycle
                    end if

                    THROW_WARN('failed to write abinitio3D_multistate metadata to ipc_pipe_abinitio3D_multstate_in')
                    exit
                end do

                if( allocated(cbuf) ) deallocate(cbuf)
            end subroutine send_to_abinitio3D_multstate_in_pipe

            ! Run ab-initio 3D classification 
            subroutine start_abinitio3D( spproj_stage, outdir, mskdiam_in )
                type(sp_project),   intent(inout) :: spproj_stage
                type(string),          intent(in) :: outdir
                integer,               intent(in) :: mskdiam_in
                type(cmdline)                     :: cline_abinitio3D
                type(string)                     :: cwd, cwd_abinitio3D, server_address
                call simple_getcwd(cwd)
                call simple_mkdir('abinitio3D')
                call simple_mkdir(outdir)
                call simple_chdir(outdir)
                call simple_getcwd(cwd_abinitio3D)
                CWD_GLOB       = cwd_abinitio3D%to_char()
                server_address = qenv%get_persistent_worker_server_address()
                call spproj_stage%write(string('abinitio3D.simple'))
                call cline_abinitio3D%kill()
                call cline_abinitio3D%set('prg',              'abinitio3D')
                call cline_abinitio3D%set('mkdir',                    'no')
                call cline_abinitio3D%set('pgrp',                     'c1')
                call cline_abinitio3D%set('nstates',             NSTATES3D)
                call cline_abinitio3D%set('lpstart',                    50) ! 20
                call cline_abinitio3D%set('lpstop',                     10) ! 6
                call cline_abinitio3D%set('force_lp_range',          'yes')
                call cline_abinitio3D%set('mskdiam',            mskdiam_in)
                call cline_abinitio3D%set('nparts',                      4)
                call cline_abinitio3D%set('nthr',                       16)
                call cline_abinitio3D%set('nstages',             NSTAGES3D)
                call cline_abinitio3D%set('projfile',  'abinitio3D.simple')
                call cline_abinitio3D%set('worker_priority',        'high')
                if( server_address%strlen() > 0 ) call cline_abinitio3D%set('worker_server', server_address)
                call cline_abinitio3D%printline()
                call qenv%exec_simple_prg_in_queue_async( cline_abinitio3D, string('./distr_abinitio3D'), string('simple_log_abinitio3D'), exec_bin=string('simple_exec') )
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
            end subroutine start_abinitio3D

            subroutine finish_abinitio3D( spproj_stage, outdir )
                type(sp_project), intent(inout) :: spproj_stage
                type(string),        intent(in) :: outdir
                type(string)                    :: cwd
                if( .not. file_exists(outdir) ) THROW_HARD('Output directory does not exist: ')
                call simple_getcwd(cwd)
                call simple_chdir(outdir)
                call spproj_stage%kill()
                call spproj_stage%read(string('abinitio3D.simple')) ! read the project with abinitio3D output
                call simple_chdir(cwd)
            end subroutine finish_abinitio3D

            ! Run refine 3D classification 
            subroutine start_refine3D( spproj_stage, outdir, mskdiam_in )
                type(sp_project),   intent(inout) :: spproj_stage
                type(string),          intent(in) :: outdir
                integer,               intent(in) :: mskdiam_in
                type(cmdline)                     :: cline_refine3D
                type(string)                      :: cwd, cwd_refine3D, server_address
                type(string)                      :: vol_fname
                real                               :: vol_smpd
                integer                           :: vol_box, state, nstates_stage
                call simple_getcwd(cwd)
                call simple_mkdir('refine3D')
                call simple_mkdir(outdir)
                call simple_chdir(outdir)
                call simple_getcwd(cwd_refine3D)
                CWD_GLOB       = cwd_refine3D%to_char()
                server_address = qenv%get_persistent_worker_server_address()
                nstates_stage  = spproj_stage%os_ptcl3D%get_n('state')
                call spproj_stage%write(string('refine3D.simple'))
                call cline_refine3D%kill()
                call cline_refine3D%set('prg',                'refine3D')
                call cline_refine3D%set('mkdir',                    'no')
                call cline_refine3D%set('balance',                  'no')
                call cline_refine3D%set('frac_best',                 1.0)
                call cline_refine3D%set('fillin',                   'no')
                call cline_refine3D%set('update_frac',               1.0)
                call cline_refine3D%set('trail_rec',                'np')
                call cline_refine3D%set('maxits',                      1)
                call cline_refine3D%set('refine',               'greedy')
                call cline_refine3D%set('greedy_sampling',         'yes')
                call cline_refine3D%set('update_missing',          'yes')
                call cline_refine3D%set('pgrp',                     'c1')
                call cline_refine3D%set('nstates',         nstates_stage)
                call cline_refine3D%set('lpstart',                    20)
                call cline_refine3D%set('lpstop',                      6)
                call cline_refine3D%set('force_lp_range',          'yes')
                call cline_refine3D%set('mskdiam',            mskdiam_in)
                call cline_refine3D%set('nparts',                      4)
                call cline_refine3D%set('nthr',                       16)
                call cline_refine3D%set('nstages',                     1)
                call cline_refine3D%set('sigma_est',            'global')
                call cline_refine3D%set('projfile',    'refine3D.simple')
                ! state volumes from the abinitio3D output segment
                do state = 1, nstates_stage
                    if( .not. spproj_stage%isthere_in_osout('vol', state) )then
                        THROW_HARD('missing state volume in project out segment; start_refine3D')
                    endif
                    call spproj_stage%get_vol('vol', state, vol_fname, vol_smpd, vol_box)
                    call cline_refine3D%set('vol'//int2str(state), vol_fname)
                enddo
                call cline_refine3D%printline()
                call qenv%exec_simple_prg_in_queue_async( cline_refine3D, string('./distr_refine3D'), string('simple_log_refine3D'), exec_bin=string('simple_exec') )
                call simple_chdir(cwd)
                CWD_GLOB = cwd%to_char()
                call vol_fname%kill
            end subroutine start_refine3D

            subroutine finish_refine3D( spproj_stage, outdir )
                type(sp_project), intent(inout) :: spproj_stage
                type(string),        intent(in) :: outdir
                type(string)                    :: cwd
                if( .not. file_exists(outdir) ) THROW_HARD('Output directory does not exist: ')
                call simple_getcwd(cwd)
                call simple_chdir(outdir)
                call spproj_stage%kill()
                call spproj_stage%read(string('refine3D.simple')) ! read the project with refine3D output
                call simple_chdir(cwd)
            end subroutine finish_refine3D

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p07_abinitio3D_multistate

end module simple_stream_p07_abinitio3D_multistate
