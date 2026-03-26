!@descr: task 3 in the stream pipeline: the first 2D analysis from segmentation picked particles used for initial screening and generation of picking references
!==============================================================================
! MODULE: simple_stream_p03_opening2D_new
!
! PURPOSE:
!   Stream pipeline stage 3 — opening 2D / generate picking references.
!   Imports micrographs from the upstream pipeline until a quality threshold
!   is reached, runs segmentation-based picking and particle extraction,
!   performs ab-initio 2D classification and shape ranking, then waits for
!   the user to select references via the GUI before handing off to the
!   reference-based picking stage.
!
! FLOW:
!   1. Accumulate micrographs (micimporter) until params%nmics accepted mics.
!   2. Segmentation-based picking  →  particle extraction.
!   3. Ab-initio 2D classification (abinitio2D)  →  shape ranking.
!   4. Send sprite-sheet of ranked classes to GUI; wait for user selection.
!   5a. User requests more micrographs → increment nmics, restart from step 1.
!   5b. User confirms reference selection → process_selected_refs, continue.
!
! RESTARTABILITY:
!   The outer `do` loop (lines 1–5 above) acts as a GOTO-free restart block.
!   `restart_requested = .true.` with `cycle` repeats from step 1;
!   `exit` after user confirmation proceeds to optics assignment and write-out.
!
! GUI MESSAGING:
!   Progress and image metadata are broadcast to the GUI via the
!   mq_stream_master_in queue using the gui_metadata_* types.
!   User updates (threshold changes, reference selections) arrive on
!   mq_stream_master_out.
!
! DEPENDENCIES:
!   simple_stream_api, simple_gui_metadata_api, simple_mini_stream_utils,
!   commander_extract, commander_abinitio2D, commander_shape_rank_cavgs
!==============================================================================
module simple_stream_p03_opening2D_new
use unix,                         only: SIGTERM
use simple_stream_api
use simple_stream_mq_defs,        only: mq_stream_master_in, mq_stream_master_out
use simple_commanders_pick,       only: commander_extract
use simple_commanders_cavgs,      only: commander_shape_rank_cavgs
use simple_commanders_abinitio2D, only: commander_abinitio2D
use simple_mini_stream_utils,     only: segdiampick_mics
use simple_gui_metadata_api

implicit none

public :: stream_p03_opening2D
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p03_opening2D
  contains
    procedure :: execute => exec_stream_p03_opening2D
end type stream_p03_opening2D

contains

    subroutine exec_stream_p03_opening2D( self, cline )
        implicit none
        class(stream_p03_opening2D), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,                  parameter       :: NCLS_MIN = 10, NCLS_MAX = 100, NPARTS2D = 4, NTHUMB_MAX = 10
        real,                     parameter       :: LPSTOP = 8.            ! low-pass stop resolution (A) for abinitio2D
        character(len=:),         allocatable     :: meta_buffer            ! serialised GUI metadata message
        type(string),             allocatable     :: projects(:)            ! batch of new project paths from the watcher
        integer,                  allocatable     :: cavg_inds(:)           ! shape-ranked class indices into os_cls2D
        type(oris)                                :: nmics_ori              ! single-ori container for writing STREAM_NMICS
        type(string)                              :: cavgsstk, mapfileprefix, projfile
        type(cmdline)                             :: cline_extract, cline_abinitio2D, cline_shape_rank
        type(parameters)                          :: params
        type(sp_project)                          :: spproj, spproj_part
        type(stream_watcher)                      :: project_buff           ! monitors dir_target for new partial projects
        type(commander_extract)                   :: xextract
        type(gui_metadata_cavg2D)                 :: meta_cavg2D
        type(commander_abinitio2D)                :: xabinitio2D
        type(gui_metadata_micrograph)             :: meta_micrograph
        type(commander_shape_rank_cavgs)          :: xshape_rank
        type(gui_metadata_stream_update)          :: meta_update            ! inbound: user selections and threshold updates
        type(gui_metadata_stream_opening2D)       :: meta_opening2D         ! outbound: 2D stage progress
        type(gui_metadata_stream_initial_picking) :: meta_initial_picking   ! outbound: micrograph / picking progress
        integer                                   :: nprojects, i
        integer                                   :: box_in_pix=0           ! box size (px) set by segdiampick_mics; 0 until known
        integer                                   :: box_for_pick           ! box size used for reference picking
        integer                                   :: box_for_extract        ! box size used for particle extraction
        integer                                   :: ithumb, xtiles, ytiles ! sprite-sheet thumbnail index and grid dims
        integer                                   :: ncls_stk, iori, i_max, optics_map_id
        integer                                   :: nmics                  ! local copy of params%nmics passed to callee
        integer                                   :: n_increase_cycles      ! tracks how many "more mics" requests have been applied
        integer                                   :: increase_nmics_gui     ! cached value of meta_update%get_increase_nmics()
        integer                                   :: ncavgs                 ! number of reference classes selected by user
        logical                                   :: restart_requested
        real                                      :: mskdiam_estimate=0     ! mask diameter (A) estimated by segdiampick_mics
        real                                      :: smpd_stk               ! pixel size of the cavgs stack
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        ! validate required args and apply defaults
        if( .not. cline%defined('dir_target')     ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls',     200)
        if( .not. cline%defined('pick_roi')       ) call cline%set('pick_roi',         'yes')
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',              '')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',              100)
        ! sanity check for restart
        if( cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir')) ) then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
        endif
        ! generate own project file
        call create_stream_project(spproj, cline, string('opening_2D'))
        ! master parameters
        call params%new(cline)
        projfile          = params%projfile
        n_increase_cycles = 0
        ! initialise metadata
        call meta_update%new(                        GUI_METADATA_STREAM_UPDATE_TYPE)
        call meta_initial_picking%new(      GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
        call meta_opening2D%new(                  GUI_METADATA_STREAM_OPENING2D_TYPE)
        call meta_micrograph%new(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
        call meta_cavg2D%new(               GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)

        ! ---------- Main pipeline loop: import → pick → extract → classify → user selection ----------
        ! cycles if the user requests more micrographs; exits on reference selection.
        do
            restart_requested = .false.
            if( file_exists(STREAM_NMICS) ) call del_file(STREAM_NMICS)
            call nmics_ori%new(1, .false.)
            call nmics_ori%set(1, 'nmics', params%nmics)
            call nmics_ori%write(1, string(STREAM_NMICS))
            call nmics_ori%kill
            ! read project fresh at the start of each restart iteration
            call spproj%read( params%projfile )
            if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(1, .false.) ! reset mic oris before fresh import
            ! wait if dir_target doesn't exist yet
            call wait_for_folder2(params%dir_target)
            call wait_for_folder2(params%dir_target//'/spprojs')
            call wait_for_folder2(params%dir_target//'/spprojs_completed')
            ! movie watcher init
            project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
            ! import at least params%nmics micrographs
            write(logfhandle, '(A,I6,A)') '>>> IMPORTING AT LEAST', params%nmics, ' MICROGRAPHS'
            call send_meta(string('waiting for micrographs'))
            call send_meta2D(string('waiting for particles'), box_in_pix)
            call micimporter( params%nmics )
            call send_meta(string('picking particles'))
            ! segmentation-based picking
            nmics = params%nmics ! local copy: segdiampick_mics must not modify params%nmics
            call segdiampick_mics(spproj, params%pcontrast, nmics, params%moldiam_max, box_in_pix, mskdiam_estimate)
            call send_meta(string('extracting particles'))
            ! send the NTHUMB_MAX most recent micrograph thumbnails to the GUI
            if( spproj%os_mic%isthere('thumb_den') .and. spproj%os_mic%isthere('xdim') .and. spproj%os_mic%isthere('ydim') &
            .and. spproj%os_mic%isthere('smpd') .and. spproj%os_mic%isthere('boxfile') ) then
                ithumb = spproj%os_mic%get_noris()   ! cache; repurposed as offset below
                i_max  = min(ithumb, NTHUMB_MAX)
                ithumb = ithumb - i_max              ! index just before the first thumbnail
                do iori = 1, i_max
                    call send_micrograph_meta(iori, i_max, ithumb + iori)
                end do
            endif
            call run_extract()
            call spproj%read(projfile)
            call send_meta(string('complete'))
            call send_meta2D(string('classifying particles'), box_in_pix)
            call run_abinitio2D()
            call send_meta2D(string('shape ranking particles'), box_in_pix)
            call run_shape_rank()
            call spproj%read(projfile)
            call spproj%shape_ranked_cavgs2jpg(cavg_inds, string("shape_ranked_")//int2str(params%nmics)//JPG_EXT,&
            &xtiles, ytiles, mskdiam_px=ceiling(mskdiam_estimate * spproj%get_smpd()))
            call spproj%get_cavgs_stk(cavgsstk, ncls_stk, smpd_stk)
            ! signal the GUI that user input is required, then send the cavgs
            call meta_opening2D%set_user_input(.true.)
            call send_meta2D(string('waiting for user selection'), box_in_pix)
            if( allocated(cavg_inds) ) call send_available_cavgs2D(&
                &string(trim(CWD_GLOB)//'/shape_ranked_'//int2str(params%nmics)//JPG_EXT), size(cavg_inds))
            ! wait for user interaction
            write(logfhandle, '(A)') ">>> WAITING FOR USER TO SELECT REFERENCES"
            do
                ! poll the outbound queue for a GUI update message
                if( mq_stream_master_out%is_active() ) then
                    if( mq_stream_master_out%receive(meta_buffer) ) then
                        if( allocated(meta_buffer) ) then
                            ! add message back to queue for other processes
                            call mq_stream_master_out%send(meta_buffer)
                            ! deserialise buffer into meta_update
                            meta_update = transfer(meta_buffer, meta_update)
                            increase_nmics_gui = meta_update%get_increase_nmics()
                            if( increase_nmics_gui > n_increase_cycles ) then
                                n_increase_cycles = increase_nmics_gui
                                params%nmics = params%nmics + NMICS_DELTA
                                call cline%set('nmics', params%nmics)
                                call cleanup_previous()
                                call meta_opening2D%set_user_input(.false.)
                                call send_meta2D(string('waiting for particles'), box_for_extract)
                                restart_requested = .true.
                                exit  ! exit inner wait-loop to restart outer loop
                            end if
                            ncavgs = meta_update%get_pickrefs_selection_length()
                            if( ncavgs > 0 ) then
                                call meta_opening2D%set_user_input(.false.)
                                call process_selected_refs(params, cavgsstk, spproj%get_smpd(),  &
                                    &meta_update%get_pickrefs_selection(),                        &
                                    &mskdiam_estimate, box_for_pick, box_for_extract, xtiles, ytiles)
                                call meta_cavg2D%kill()
                                call meta_cavg2D%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE)
                                call send_meta2D(string('applying user selection'), box_for_extract)
                                call send_available_cavgs2D(&
                                    &string(trim(CWD_GLOB)//'/'//STREAM_SELECTED_REFS//JPG_EXT), ncavgs)
                                exit  ! exit inner wait-loop; restart_requested already .false.
                            end if
                        end if
                    end if
                end if
                call sleep(WAITTIME)
            enddo
            if( restart_requested ) cycle  ! repeat the restartable block
            exit                           ! normal completion of restart loop
        end do ! end restart loop
        call send_meta2D(string('terminating'), box_for_extract)
        if( allocated(projects)  ) deallocate(projects)
        if( allocated(cavg_inds) ) deallocate(cavg_inds)
        ! add optics
        if( cline%defined('optics_dir') ) then
            optics_map_id = get_latest_optics_map_id(params%optics_dir)
            if( optics_map_id > 0 ) then
                mapfileprefix = params%optics_dir//'/'//OPTICS_MAP_PREFIX//int2str(optics_map_id)
                call spproj%import_optics_map(mapfileprefix)
            endif
        endif
        ! write project and star files (just in case you want to import these particles/micrographs elsewhere)
        call spproj%write
        call spproj%write_mics_star(string("micrographs.star"))
        call spproj%write_ptcl2D_star(string("particles.star"))
        call spproj%kill
        call simple_end('**** SIMPLE_GEN_PICKREFS NORMAL STOP ****')

        contains

            ! Delete all files in the working directory and re-create a clean project,
            ! called when the user requests more micrographs before picking references.
            subroutine cleanup_previous()
                type(string), allocatable :: list(:)
                integer :: i, nfiles
                call spproj%kill
                call cline%set('projfile', '')
                call simple_list_files('*', list)
                nfiles = size(list)
                if( nfiles > 0 ) write(logfhandle,'(A)')'>>> CLEANING UP'
                do i = 1, nfiles
                    call del_file(list(i))
                enddo
                call cline%set('projname', params%projname)
                call cline%set('projfile', params%projfile)
                call spproj%update_projinfo(cline)
                call spproj%update_compenv(cline)
                call spproj%write
            end subroutine cleanup_previous

            ! Extract particles using one thread per partition (nthr partitions).
            subroutine run_extract()
                call cline_extract%set('prg',        'extract')
                call cline_extract%set('mkdir',           'no')
                call cline_extract%set('nparts',   params%nthr)
                call cline_extract%set('nthr',               1)
                call cline_extract%set('projfile',    projfile)
                call xextract%execute(cline_extract)
            end subroutine run_extract

            ! Run ab-initio 2D classification on the extracted particles.
            ! ncls is clamped to [NCLS_MIN, NCLS_MAX] based on particle count.
            subroutine run_abinitio2D()
                integer :: nptcls, ncls, nthr2D
                nptcls = spproj%os_ptcl2D%get_noris()
                ncls   = min(NCLS_MAX, max(NCLS_MIN, nptcls/params%nptcls_per_cls))
                nthr2D = max(1,floor(real(params%nthr)/4.))
                call cline_abinitio2D%set('prg',               'abinitio2D')
                call cline_abinitio2D%set('mkdir',                     'no')
                call cline_abinitio2D%set('ncls',                      ncls)
                call cline_abinitio2D%set('sigma_est',             'global')
                call cline_abinitio2D%set('center',                   'yes')
                call cline_abinitio2D%set('autoscale',                'yes')
                call cline_abinitio2D%set('lpstop',                  LPSTOP)
                call cline_abinitio2D%set('mskdiam',       mskdiam_estimate)
                call cline_abinitio2D%set('nthr',                    nthr2D)
                call cline_abinitio2D%set('nparts',                NPARTS2D)
                call cline_abinitio2D%set('projfile',              projfile)
                call xabinitio2D%execute(cline_abinitio2D)
            end subroutine run_abinitio2D

            ! Shape-rank the 2D class averages by particle quality.
            subroutine run_shape_rank()
                call cline_shape_rank%set('nthr',               params%nthr)
                call cline_shape_rank%set('projfile',              projfile)
                call xshape_rank%execute(cline_shape_rank)
            end subroutine run_shape_rank

            ! Block until at least nmics micrographs pass quality thresholds
            ! (ctfres, icefrac, astig).  Polls project_buff for new partial
            ! projects and applies rejection after each batch.
            subroutine micimporter( nmics )
                integer, intent(in) :: nmics
                integer :: n_imported, n_new_oris, n_oris, iproj, iori, imic
                n_imported = 0
                do
                    if( file_exists(TERM_STREAM) ) then
                        ! termination
                        write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                        call spproj%kill
                        call qsys_cleanup(params)
                        call simple_end('**** SIMPLE_GEN_PICKREFS USER STOP ****')
                        call EXIT(0)
                    endif
                    ! detection of new projects
                    call project_buff%watch( nprojects, projects, max_nmovies=50 )
                    ! append projects to processing stack
                    if( nprojects > 0 ) then
                        n_imported = spproj%os_mic%get_noris()
                        if( n_imported > 0 ) then
                            n_new_oris = n_imported + nprojects * STREAM_NMOVS_SET
                            call spproj%os_mic%reallocate(n_new_oris)
                        else
                            n_new_oris = nprojects * STREAM_NMOVS_SET
                            call spproj%os_mic%new(n_new_oris, .false.)
                        end if
                        do iproj = 1, nprojects
                            call project_buff%add2history(projects(iproj)) ! mark as seen so the watcher won't deliver it again
                            call spproj_part%read(projects(iproj))
                            do iori = 1, STREAM_NMOVS_SET
                                n_imported = n_imported + 1
                                ! set state=0 mics to state=-1
                                if( spproj_part%os_mic%get(iori, 'state') < 1 ) call spproj_part%os_mic%set(iori, 'state', -1)
                                call spproj%os_mic%transfer_ori(n_imported, spproj_part%os_mic, iori)
                            end do
                            call spproj_part%kill()
                        enddo
                        write(logfhandle,'(A,I4,A,A)')'>>> ', nprojects * STREAM_NMOVS_SET, ' NEW MICROGRAPHS DETECTED; ', cast_time_char(simple_gettime())
                    else
                        call sleep(WAITTIME)
                    endif
                    ! apply quality thresholds (ctfres, icefrac, astig) to all accepted mics
                    n_oris = spproj%os_mic%get_noris()
                    do imic = 1, n_oris
                        if( spproj%os_mic%get(imic, 'state') < 0 ) cycle
                        call spproj%os_mic%set(imic, 'state', 1) ! reset to accepted before threshold test
                        if( spproj%os_mic%isthere(imic, 'ctfres') ) then
                            if( spproj%os_mic%get(imic,'ctfres') > (params%ctfresthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                        if( spproj%os_mic%isthere(imic, 'icefrac') ) then
                            if( spproj%os_mic%get(imic,'icefrac') > (params%icefracthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                        if( spproj%os_mic%isthere(imic, 'astig') ) then
                            if( spproj%os_mic%get(imic,'astig') > (params%astigthreshold-0.001) ) call spproj%os_mic%set(imic, 'state', 0)
                        end if
                    enddo
                    ! enough accepted micrographs — expose dead mics as state=0 and return
                    if( spproj%os_mic%count_state_gt_zero() >= nmics ) then
                        do imic = 1, n_oris
                            if( spproj%os_mic%get(imic, 'state') < 0 ) call spproj%os_mic%set(imic, 'state', 0)
                        enddo
                        return
                    endif
                end do
            end subroutine micimporter

            ! Broadcast initial-picking progress to the GUI.
            subroutine send_meta( my_stage )
                type(string), intent(in) :: my_stage
                call meta_initial_picking%set(                                   &
                    stage                = my_stage,                             &
                    micrographs_imported = spproj%os_mic%get_noris(),            &
                    micrographs_accepted = spproj%os_mic%count_state_gt_zero(),  &
                    particles_extracted  = spproj%os_ptcl2D%get_noris())
                if( meta_initial_picking%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_initial_picking%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_meta

            ! Broadcast 2D-classification stage progress to the GUI.
            ! Always sends — stage-only updates (e.g. 'waiting for particles') must
            ! reach the GUI even before any particles exist; smpd is 0 in that case.
            subroutine send_meta2D( my_stage, box_size )
                type(string), intent(in) :: my_stage
                integer,      intent(in) :: box_size
                real                     :: smpd
                integer                  :: nptcls
                nptcls = spproj%os_ptcl2D%get_noris()
                smpd   = 0.0
                if( nptcls > 0 ) smpd = spproj%get_smpd()
                call meta_opening2D%set(                                             &
                    stage              = my_stage,                                   &
                    particles_imported = nptcls,                                     &
                    particles_accepted = spproj%os_ptcl2D%count_state_gt_zero(),     &
                    mask_diam          = nint(mskdiam_estimate),                     &
                    mask_scale         = box_size * smpd,                            &
                    box_size           = box_size)
                if( meta_opening2D%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_opening2D%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_meta2D

            ! Send thumbnail, CTF values, and particle coordinates for one
            ! micrograph (index my_ithumb) as message my_iori of my_i_max.
            subroutine send_micrograph_meta( my_iori, my_i_max, my_ithumb )
                integer, intent(in) :: my_iori, my_i_max, my_ithumb
                type(nrtxtfile)     :: my_boxfile
                type(string)        :: my_boxpath
                real,  allocatable  :: boxdata(:)
                integer             :: my_i, my_x, my_y, my_nrecs, my_nlines, my_xdim, my_ydim
                call meta_micrograph%set(                                     &
                    path   = spproj%os_mic%get_str(my_ithumb, "thumb"),       &
                    dfx    = spproj%os_mic%get(my_ithumb,       "dfx"),       &
                    dfy    = spproj%os_mic%get(my_ithumb,       "dfy"),       &
                    ctfres = spproj%os_mic%get(my_ithumb,    "ctfres"),       &
                    i      = my_iori,                                         &
                    i_max  = my_i_max)
                call meta_micrograph%clear_coordinates()
                my_boxpath = spproj%os_mic%get_str(my_ithumb, "boxfile")
                if( my_boxpath%strlen() > 0 .and. file_exists(my_boxpath) ) then
                    call my_boxfile%new(my_boxpath, 1)
                    my_nrecs  = my_boxfile%get_nrecs_per_line()
                    my_nlines = my_boxfile%get_ndatalines()
                    my_xdim   = spproj%os_mic%get_int(my_ithumb, "xdim")
                    my_ydim   = spproj%os_mic%get_int(my_ithumb, "ydim")
                    if( my_nrecs >= 4 ) then
                        allocate(boxdata(my_nrecs))
                        do my_i = 1, my_nlines
                            call my_boxfile%readNextDataLine(boxdata)
                            my_x = nint(boxdata(1) + boxdata(3)/2)
                            my_y = nint(boxdata(2) + boxdata(4)/2)
                            call meta_micrograph%set_coordinate(my_i, my_x, my_y, my_xdim, my_ydim)
                        enddo
                        deallocate(boxdata)
                    endif
                    call my_boxfile%kill()
                    if( meta_micrograph%assigned() .and. mq_stream_master_in%is_active() ) then
                        call meta_micrograph%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                endif
            end subroutine send_micrograph_meta

            ! Send metadata for one 2D class average to the GUI.
            ! my_xtile/my_ytile are the 0-based grid position within the sprite sheet.
            subroutine send_cavg2D_meta( my_path, my_i, my_i_max, my_xtile, my_ytile )
                type(string), intent(in) :: my_path
                integer,      intent(in) :: my_i, my_i_max, my_xtile, my_ytile
                integer                  :: my_idx
                my_idx = cavg_inds(my_i)
                call meta_cavg2D%set(                                    &
                    path    = my_path,                                   &
                    mrcpath = cavgsstk,                                  &
                    i       = my_i,                                      &
                    i_max   = my_i_max,                                  &
                    res     = spproj%os_cls2D%get(my_idx,     'res'),    &
                    pop     = spproj%os_cls2D%get_int(my_idx, 'pop'),    &
                    idx     = my_idx,                                    &
                    sprite  = sprite_sheet_pos(                                        &
                                  x = merge(0.0, my_xtile * (100.0 / (xtiles - 1)), xtiles == 1), &
                                  y = merge(0.0, my_ytile * (100.0 / (ytiles - 1)), ytiles == 1), &
                                  h = 100 * ytiles,                      &
                                  w = 100 * xtiles))
                if( meta_cavg2D%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_cavg2D%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_cavg2D_meta

            ! Send a batch of cavg2D metadata to the GUI, resetting tile counters.
            subroutine send_available_cavgs2D( my_path, n )
                type(string), intent(in) :: my_path
                integer,      intent(in) :: n
                integer                  :: my_xtile, my_ytile
                my_xtile = 0
                my_ytile = 0
                do i = 1, n
                    write(*,*) 'send_available_cavgs2D', my_path%to_char() , i, n, my_xtile, my_ytile
                    call send_cavg2D_meta(my_path, i, n, my_xtile, my_ytile)
                    my_xtile = my_xtile + 1
                    if( my_xtile == xtiles ) then
                        my_xtile = 0
                        my_ytile = my_ytile + 1
                    endif
                end do
            end subroutine send_available_cavgs2D

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                call exit(1)
            end subroutine sigterm_handler

    end subroutine exec_stream_p03_opening2D

end module simple_stream_p03_opening2D_new
