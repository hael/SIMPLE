!@descr: task 3 in the stream pipeline: the first 2D analysis from segmentation picked particles used for initial screening and generation of picking references
!==============================================================================
! MODULE: simple_stream_p03_initial_analysis
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
module simple_stream_p03_initial_analysis
use unix,                         only: SIGTERM
use simple_stream_api
use simple_stream_mq_defs,        only: mq_stream_master_in, mq_stream_master_out
use simple_commanders_pick,       only: commander_extract, commander_reextract
use simple_commanders_cavgs,      only: commander_shape_rank_cavgs
use simple_commanders_abinitio2D, only: commander_abinitio2D
use simple_commanders_abinitio,   only: commander_abinitio3D_cavgs_reject, commander_abinitio3D_cavgs
use simple_commanders_reproject,  only: commander_reproject
use simple_commanders_cluster2D,  only: commander_cls_split
use simple_commanders_mkcavgs,    only: commander_make_cavgs
use simple_stream_cluster2D_microchunked, only: stream_cluster2D_microchunked
use simple_mini_stream_utils,     only: segdiampick_mics_multi
use simple_qsys_env,              only: qsys_env
use simple_cavg_quality_analysis, only: evaluate_cavg_quality
use simple_cavg_quality_model,    only: cavg_quality_model, CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_cavg_quality_types,    only: cavg_quality_result
use simple_imgarr_utils,          only: dealloc_imgarr, read_cavgs_into_imgarr, read_stk_into_imgarr
use simple_image_msk,             only: automask2D
use simple_projfile_utils,        only: merge_chunk_projfiles, merge_selected_project_files
use simple_procimgstk,            only: scale_imgfile
use simple_fileio,                only: swap_suffix
use simple_defs,                  only: MSK_EXP_FAC, BOX_EXP_FAC, COSMSKHALFWIDTH
use simple_abinitio_utils,        only: abinitio_rec_fbody
use simple_gui_metadata_api

implicit none

public :: stream_p03_initial_analysis
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p03_initial_analysis
  contains
    procedure :: execute => exec_stream_p03_initial_analysis
end type stream_p03_initial_analysis

contains

    subroutine exec_stream_p03_initial_analysis( self, cline )
        implicit none
        class(stream_p03_initial_analysis), intent(inout) :: self
        class(cmdline),              intent(inout) :: cline
        integer,                   parameter       :: NCLS_MIN = 10, NCLS_MAX = 100, NPARTS2D = 8, NTHUMB_MAX = 10
    !    integer,                   parameter       :: NMICS_PLAN(3) = [10, 20, 50]  ! number of micrographs to import for each cycle of the opening2D plan; must have at least 2 steps and no more than 9 (for IPC routing via single-digit cluster counts)
        integer,                   parameter       :: NMICS_PLAN(1) = [100]
        real,                      parameter       :: LPSTOP2D = 8.          ! low-pass stop resolution (A) for abinitio2D/3D setup
        integer,                   parameter       :: NSTATES3D = 3                 ! number of classes for abinitio3D
        integer,                   parameter       :: NSTAGES3D = 4                 ! number of stages for abinitio3D
        character(len=:),          allocatable     :: meta_buffer            ! serialised GUI metadata message
        type(string),              allocatable     :: projects(:)            ! batch of new project paths from the watcher
        type(string),              allocatable     :: imgfiles(:)            ! cache of cavgs stack paths for each cluster, for use in process_selected_refs
        type(oris)                                 :: nmics_ori              ! single-ori container for writing STREAM_NMICS
        type(string)                               :: projfile, cwd_master, cwd_cycle, cycle_dir, cycle_projfile, boxfile
        type(qsys_env)                             :: qsys
        type(cmdline)                              :: cline_extract, cline_abinitio2D, cline_shape_rank, cline_reextract, cline_abinitio3D, cline_cls_split
        type(parameters)                           :: params
        type(sp_project)                           :: spproj, spproj_part
        type(stream_watcher)                       :: project_buff           ! monitors dir_target for new partial projects
        type(commander_extract)                    :: xextract
        type(commander_reextract)                  :: xreextract
        type(gui_metadata_cavg2D)                  :: meta_cavg2D
        type(commander_abinitio2D)                 :: xabinitio2D
        type(commander_abinitio3D_cavgs_reject)    :: xabinitio3D_cavgs_reject
        type(commander_reproject)                  :: xreproject
        type(commander_cls_split)                  :: xcls_split
        type(commander_make_cavgs)                 :: xmake_cavgs
        type(gui_metadata_micrograph)              :: meta_micrograph
        type(commander_shape_rank_cavgs)           :: xshape_rank
        type(gui_metadata_stream_update)           :: meta_update            ! inbound: user selections and threshold updates
        type(gui_metadata_stream_opening2D)        :: meta_opening2D         ! outbound: 2D stage progress
        type(gui_metadata_stream_picking)          :: meta_initial_picking   ! outbound: micrograph / picking progress
        integer                                    :: nprojects, i
        integer                                    :: box_in_pix      =0      ! box size (px) set by segdiampick_mics; 0 until known
        integer                                    :: box_for_pick    =0      ! box size used for reference picking
        integer                                    :: box_for_extract =0      ! box size used for particle extraction
        integer                                    :: ithumb, xtiles, ytiles ! sprite-sheet thumbnail index and grid dims
        integer                                    :: iori, i_max, i_start
        integer                                    :: nmics                  ! local copy of params%nmics passed to callee
        integer                                    :: n_selected_cavgs       ! number of quality-selected cavgs in current cycle
        integer                                    :: n_cycles      ! tracks how many "more mics" requests have been applied
        integer                                    :: ncavgs                 ! number of reference classes selected by user
        logical                                    :: restart_requested
        real                                       :: mskdiam                ! mask diameter (A) for reference-based picking, set by segdiampick_mics; 0 until known
        real                                       :: mskdiam_estimate=0     ! mask diameter (A) estimated by segdiampick_mics
        real                                       :: smpd_stk               ! pixel size of the cavgs stack
        call signal(SIGTERM, sigterm_handler)   ! graceful shutdown on SIGTERM
        ! validate required args and apply defaults
        if( .not. cline%defined('dir_target')     ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('mkdir')          ) call cline%set('mkdir',            'yes')
        if( .not. cline%defined('nptcls_per_cls') ) call cline%set('nptcls_per_cls',     100)
        if( .not. cline%defined('pick_roi')       ) call cline%set('pick_roi',         'yes')
        if( .not. cline%defined('outdir')         ) call cline%set('outdir',              '')
        if( .not. cline%defined('nmics')          ) call cline%set('nmics',    NMICS_PLAN(1))
        ! arguments for automask2D
        if( .not. cline%defined('ngrow')  ) call cline%set('ngrow',    3)
        if( .not. cline%defined('winsz')  ) call cline%set('winsz',   5.)
        if( .not. cline%defined('amsklp') ) call cline%set('amsklp', 20.)
        if( .not. cline%defined('edge')   ) call cline%set('edge',     6)
        ! sanity check for restart
        if( cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir')) ) then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
        endif
        ! generate own project file
        call create_stream_project(spproj, cline, string('opening_2D'))
        ! master parameters
        call params%new(cline)
        projfile           = params%projfile
        params%workers     = NPARTS2D
        params%worker_nthr = max(1,floor(real(params%nthr)/4.))
        params%nmics       = NMICS_PLAN(1)
        call cline%set('nmics', params%nmics)
        n_cycles  = 0
        ! initialise metadata
        call meta_update%new(                        GUI_METADATA_STREAM_UPDATE_TYPE)
        call meta_initial_picking%new(      GUI_METADATA_STREAM_INITIAL_PICKING_TYPE)
        call meta_opening2D%new(                  GUI_METADATA_STREAM_OPENING2D_TYPE)
        call meta_micrograph%new(GUI_METADATA_STREAM_INITIAL_PICKING_MICROGRAPH_TYPE)
        call meta_cavg2D%new(               GUI_METADATA_STREAM_OPENING2D_CLS2D_TYPE)
        ! stash cwd for constructing absolute paths to send to the GUI
        call simple_getcwd(cwd_master)

        ! ---------- Main pipeline loop: import → pick → extract → classify → user selection ----------
        ! cycles if the user requests more micrographs; exits on reference selection.
        do
            restart_requested = .false.
            call simple_chdir(cwd_master) ! ensure we start each cycle in a known location
            if( file_exists(STREAM_NMICS) ) call del_file(STREAM_NMICS)
            call nmics_ori%new(1, .false.)
            call nmics_ori%set(1, 'nmics', params%nmics)
            call nmics_ori%write(1, string(STREAM_NMICS))
            call nmics_ori%kill
            ! read project fresh at the start of each restart iteration
            call spproj%read( params%projfile )
            if( spproj%os_mic%get_noris() /= 0    ) call spproj%os_mic%new(1, .false.) ! reset mic oris before fresh import
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
            cycle_dir      = string('cycle_'//int2str(n_cycles + 1))
            cycle_projfile = string(cycle_dir%to_char()//METADATA_EXT)
            call simple_mkdir(cycle_dir)
            call simple_chdir(cycle_dir)
            call spproj%update_projinfo(cycle_projfile)
            call spproj%write()
            call send_meta(string('picking particles'))
            ! segmentation-based picking
            nmics = params%nmics ! local copy: segdiampick_mics must not modify params%nmics
            call segdiampick_mics_multi(spproj, params%pcontrast, nmics, params%moldiam_max, box_in_pix, mskdiam)
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
            !
            call qsys%new(params, NPARTS2D)
            call qsys%kill()
            call run_extract(spproj, cycle_projfile, string('extract'), box_in_pix)
            do i = spproj%os_mic%get_noris(), 1, -1
                if( spproj%os_mic%get(i, 'state') < 1.0 )then
                    call spproj%os_mic%delete(i)
                    cycle
                endif
                if( .not. spproj%os_mic%isthere(i, 'nptcls') )then
                    call spproj%os_mic%delete(i)
                    cycle
                endif
                if( spproj%os_mic%get(i, 'nptcls') <= 0.0 )then
                    call spproj%os_mic%delete(i)
                    cycle
                endif
                if( .not. spproj%os_mic%isthere(i, 'boxfile') )then
                    call spproj%os_mic%delete(i)
                    cycle
                endif
                boxfile = spproj%os_mic%get_str(i, 'boxfile')
                if( boxfile%strlen() == 0 .or. .not. file_exists(boxfile) )then
                    call boxfile%kill
                    call spproj%os_mic%delete(i)
                    cycle
                endif
                call boxfile%kill
            enddo
            call spproj%write(cycle_projfile)
            call send_meta(string('complete'))
            call send_meta2D(string('classifying particles'), box_in_pix)
            call run_abinitio2D_microchunked(spproj, cycle_projfile, string('abinitio2D_microchunked'), nint(mskdiam))
            call run_abinitio2D(spproj, cycle_projfile, string('abinitio2D'), nint(mskdiam))
            call send_meta2D(string('evaluating class average quality'), box_in_pix)
            i_max   = 0
            i_start = 1
            call simple_getcwd(cwd_cycle) ! cache master CWD for constructing absolute paths to send to the GUI
            call simple_copy_file(string('abinitio2D/')//cycle_projfile, cycle_projfile)
            call run_cavg_quality_selection(spproj, cycle_projfile, string('quality_selection_1'),cwd_cycle, mskdiam, imgfiles, smpd_stk, n_selected_cavgs)
            call reestimate_box_size_from_selected_cavgs(spproj, params, mskdiam, box_in_pix)
           ! call run_reextract(spproj, cycle_projfile, string('reextract'), box_in_pix)
          !  call run_make_cavgs(spproj, cycle_projfile, string('make_cavgs_1'))
            if( n_cycles >= size(NMICS_PLAN) - 1 ) then 
                call balance_classes(spproj, cycle_projfile, string('balance_classes')) ! balance class populations before final abinitio2D, to improve quality of top classes and thus picking references
              !  call run_cls_split(spproj, cycle_projfile, string('cls_split'), nint(mskdiam))
              !  call run_make_cavgs(spproj, cycle_projfile, string('make_cavgs_2'))
             !   call run_cavg_quality_selection(spproj, cycle_projfile, string('quality_selection_2'), cwd_cycle, mskdiam, imgfiles, smpd_stk, n_selected_cavgs)
                call run_abinitio3D_and_reproject(spproj, cycle_projfile, string('abinitio3D'), nint(mskdiam))
            else
                write(logfhandle, '(A,I0,A)') '>>> SKIPPING ABINITIO3D UNTIL PLAN STEP 3 (CURRENT STEP=', &
                    n_cycles + 1, ')'
            end if
            i_max = i_max + n_selected_cavgs   ! accumulate count of selected classes for IPC routing
            if( n_cycles < size(NMICS_PLAN) - 1 ) then
                n_cycles = n_cycles + 1
                params%nmics      = NMICS_PLAN(n_cycles + 1)
                call cline%set('nmics', params%nmics)
                write(logfhandle, '(A,I0,A,I0)') '>>> ADVANCING TO PLAN STEP ', n_cycles + 1, ' WITH NMICS=', params%nmics
             !   call cleanup_previous()
                call restart_cleanup_allocatables()
                call meta_opening2D%set_user_input(.false.)
                call spproj%kill()
                cycle
            end if
            call meta_opening2D%set_user_input(.true.)
            call send_meta2D(string('waiting for user selection'), box_in_pix)
            i_start = 1
            call spproj%kill()
            call spproj%read(cycle_projfile) ! read the full project with all clusters merged in, to get the complete set of class averages for user selection
       !     call spproj%cavgs2jpg(cavg_inds, string("diameter_cluster_")//int2str(i_cluster)//"_"//int2str(params%nmics)//JPG_EXT, xtiles, ytiles, ignore_states=.false.)
       !     if( allocated(cavg_inds) ) then
       !         cavg_inds = pack(cavg_inds, cavg_inds > 0) ! filter out zero indices (unselected classes)
       !         if( size(cavg_inds) > 0 ) then
       !             call send_available_cavgs2D(cwd_cycle//'/diameter_cluster_'//int2str(i_cluster)//"_"//int2str(params%nmics)//JPG_EXT, size(cavg_inds), i_cluster, my_i_max=i_max, my_i_start=i_start)
       !             i_start = i_start + size(cavg_inds)
       !         end if
       !     end if
           ! end do
            ! wait for user interaction
            write(logfhandle, '(A)') ">>> WAITING FOR USER TO SELECT REFERENCES"

            exit

            do
                ! poll the outbound queue for a GUI update message
                if( mq_stream_master_out%is_active() ) then
                    if( mq_stream_master_out%receive(meta_buffer) ) then
                        if( allocated(meta_buffer) ) then
                            ! add message back to queue for other processes
                            call mq_stream_master_out%send(meta_buffer)
                            ! deserialise buffer into meta_update
                            meta_update = transfer(meta_buffer, meta_update)
                            ncavgs = meta_update%get_pickrefs_selection_length()
                            if( ncavgs > 0 ) then
                                call meta_opening2D%set_user_input(.false.)
                                call process_selected_refs_2(params, imgfiles, smpd_stk,            &
                                    &meta_update%get_pickrefs_selection(),                        &
                                    &meta_update%get_pickrefs_clusters(),                         &
                                    &mskdiam_estimate, box_for_pick, box_for_extract, xtiles, ytiles)
                                call meta_cavg2D%kill()
                                call meta_cavg2D%new(GUI_METADATA_STREAM_OPENING2D_CLS2D_FINAL_TYPE)
                                call send_meta2D(string('applying user selection'), box_for_extract)
                                !call send_available_cavgs2D(&
                                !    &cwd_cycle//'/'//STREAM_SELECTED_REFS//JPG_EXT, ncavgs)
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
        if( allocated(imgfiles)  ) deallocate(imgfiles)
        ! add optics
        ! if( cline%defined('optics_dir') ) then
        !     optics_map_id = get_latest_optics_map_id(params%optics_dir)
        !     if( optics_map_id > 0 ) then
        !         mapfileprefix = params%optics_dir//'/'//OPTICS_MAP_PREFIX//int2str(optics_map_id)
        !         call spproj%import_optics_map(mapfileprefix)
        !     endif
        ! endif
        ! ! write project and star files (just in case you want to import these particles/micrographs elsewhere)
        ! call spproj%write
        ! call spproj%write_mics_star(string("micrographs.star"))
        ! call spproj%write_ptcl2D_star(string("particles.star"))
        call spproj%kill
        call simple_end('**** SIMPLE_GEN_PICKREFS NORMAL STOP ****')

        contains

            ! Delete all files in the working directory and re-create a clean project,
            ! called when the user requests more micrographs before picking references.
            subroutine cleanup_previous()
                type(string), allocatable :: list(:), folders(:)
                integer :: i_file, nfiles
                call spproj%kill
                call cline%set('projfile', '')
                call simple_list_files('*', list)
                nfiles = size(list)
                if( nfiles > 0 ) write(logfhandle,'(A)')'>>> CLEANING UP'
                do i_file = 1, nfiles
                    call del_file(list(i_file))
                enddo
                folders = simple_list_dirs('.')
                if( allocated(folders) )then
                    do i_file = 1,size(folders)
                        if( folders(i_file)%has_substr('extract_cluster_')    ) call simple_rmdir(folders(i_file))
                        if( folders(i_file)%has_substr('abinitio2D_cluster_') ) call simple_rmdir(folders(i_file))
                    enddo
                endif
                call cline%set('projname', params%projname)
                call cline%set('projfile', params%projfile)
                call create_stream_project(spproj, cline, string('opening_2D'))
            end subroutine cleanup_previous

            ! Clear loop-carried allocatables before restarting the outer cycle.
            subroutine restart_cleanup_allocatables()
                if( allocated(imgfiles)       ) deallocate(imgfiles)
            end subroutine restart_cleanup_allocatables

            ! Extract particles using one thread per partition (nthr partitions).
            subroutine run_extract( spproj_inout, cluster_projfile, outdir, box )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in) :: cluster_projfile
                type(string), intent(in) :: outdir
                integer,      intent(in) :: box
                type(string)             :: cwd
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call cline_extract%kill()
                call cline_extract%set('prg',              'extract')
                call cline_extract%set('box',                    box)
                call cline_extract%set('nparts',         params%nthr)
                call cline_extract%set('nthr',                     1)
                call cline_extract%set('mkdir',                 'no')
                call cline_extract%set('projfile',  cluster_projfile)
                call cline_extract%printline()
                call xextract%execute(cline_extract)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy extracted particles back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_extract

            ! Re-extract particles at the updated box size after user selection.
            subroutine run_reextract( spproj_inout, cluster_projfile, outdir, box )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in) :: cluster_projfile
                type(string), intent(in) :: outdir
                integer,      intent(in) :: box
                type(string)             :: cwd
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to reextract dir for partitioning
                call simple_chdir(outdir)
                call cline_reextract%kill()
                call cline_reextract%set('prg',            'reextract')
                call cline_reextract%set('box',                    box)
                call cline_reextract%set('oritype',            'ptcl2D')
                call cline_reextract%set('nparts',         params%nthr)
                call cline_reextract%set('nthr',                    1)
                call cline_reextract%set('mkdir',                'no')
                call cline_reextract%set('projfile', cluster_projfile)
                call cline_reextract%printline()
                call xreextract%execute(cline_reextract)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy reextracted particles back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_reextract

            ! Run ab-initio 2D classification on the extracted particles.
            ! ncls is clamped to [NCLS_MIN, NCLS_MAX] based on particle count.
            subroutine run_abinitio2D( spproj_inout, cluster_projfile, outdir, mskdiam_in )
                type(sp_project), intent(inout) :: spproj_inout
                type(string),        intent(in) :: cluster_projfile
                type(string),        intent(in) :: outdir
                integer,             intent(in) :: mskdiam_in
                type(sp_project)                :: spproj_cluster
                type(string)                    :: cwd
                integer :: nptcls, ncls_job, nthr2D
                call spproj_cluster%read(cluster_projfile)
                nptcls = spproj_cluster%os_ptcl2D%get_noris()
                ncls_job = min(NCLS_MAX, max(NCLS_MIN, nptcls/params%nptcls_per_cls))
                nthr2D = max(1,floor(real(params%nthr)/4.))
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call spproj_cluster%kill()
                call cline_abinitio2D%kill()
                call cline_abinitio2D%set('prg',               'abinitio2D')
                call cline_abinitio2D%set('mkdir',                     'no')
                call cline_abinitio2D%set('ncls',                  ncls_job)
                call cline_abinitio2D%set('sigma_est',             'global')
                call cline_abinitio2D%set('center',                   'yes')
                call cline_abinitio2D%set('autoscale',                'yes')
                call cline_abinitio2D%set('lpstop',                LPSTOP2D)
                call cline_abinitio2D%set('mskdiam',             mskdiam_in)
                call cline_abinitio2D%set('nthr',                    nthr2D)
                call cline_abinitio2D%set('nparts',                NPARTS2D)
                call cline_abinitio2D%set('projfile',      cluster_projfile)
                call cline_abinitio2D%printline()
                call xabinitio2D%execute(cline_abinitio2D)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy abinitio2D output back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_abinitio2D
            
            subroutine run_abinitio2D_microchunked( spproj_inout, cluster_projfile, outdir, mskdiam_in )
                type(sp_project), intent(inout) :: spproj_inout
                type(string),        intent(in) :: cluster_projfile
                type(string),        intent(in) :: outdir
                integer,             intent(in) :: mskdiam_in
                type(stream_cluster2D_microchunked) :: xcluster2D_microchunked
                type(cmdline)                   :: cline_abinitio2D_microchunked
                type(string)                    :: cwd
                integer :: nptcls, ncls_job, nthr2D
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call cline_abinitio2D_microchunked%kill()
                call cline_abinitio2D_microchunked%set('prg',               'abinitio2D_microchunked')
                call cline_abinitio2D_microchunked%set('mkdir',                     'no')
                call cline_abinitio2D_microchunked%set('nmics',                       50)
                call cline_abinitio2D_microchunked%set('maxnptcls',                 5000)
                call cline_abinitio2D_microchunked%set('mskdiam',             mskdiam_in)
                call cline_abinitio2D_microchunked%set('nthr',                        16)
                call cline_abinitio2D_microchunked%set('nchunks',                      2)
                call cline_abinitio2D_microchunked%set('projfile',      cluster_projfile)
                call cline_abinitio2D_microchunked%printline()
                call xcluster2D_microchunked%execute(cline_abinitio2D_microchunked)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy abinitio2D output back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_abinitio2D_microchunked

            subroutine run_cls_split( spproj_inout, cluster_projfile, outdir, mskdiam_in )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in) :: cluster_projfile
                type(string), intent(in) :: outdir
                integer,      intent(in) :: mskdiam_in
                type(string)             :: cwd
                integer                  :: nthr2D
                nthr2D = max(1,floor(real(params%nthr)/4.))
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call cline_cls_split%kill()
                call cline_cls_split%set('prg',                'cls_split')
                call cline_cls_split%set('mkdir',                     'no')
                call cline_cls_split%set('mskdiam',             mskdiam_in)
                call cline_cls_split%set('nthr',                    nthr2D)
                call cline_cls_split%set('nparts',                NPARTS2D)
                call cline_cls_split%set('projfile',      cluster_projfile)
                call cline_cls_split%printline()
                call xcls_split%execute(cline_cls_split)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy cls_split output back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_cls_split

            subroutine run_make_cavgs( spproj_inout, cluster_projfile, outdir )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in) :: cluster_projfile
                type(string), intent(in) :: outdir
                type(cmdline)            :: cline_make_cavgs
                type(string)             :: cwd
                integer                  :: nthr2D
                nthr2D = max(1,floor(real(params%nthr)/4.))
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call cline_make_cavgs%kill()
                call cline_make_cavgs%set('prg',               'make_cavgs')
                call cline_make_cavgs%set('mkdir',                     'no')
                call cline_make_cavgs%set('nthr',                    nthr2D)
                call cline_make_cavgs%set('nparts',                NPARTS2D)
                call cline_make_cavgs%set('projfile',      cluster_projfile)
                call cline_make_cavgs%printline()
                call xmake_cavgs%execute(cline_make_cavgs)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy make_cavgs output back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
                
            end subroutine run_make_cavgs

            ! Run ab-initio 3D classification on the extracted particles.
            subroutine run_abinitio3D_and_reproject( spproj_inout, cluster_projfile, outdir, mskdiam_in )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in) :: cluster_projfile      ! cycle-local project file to process
                type(string), intent(in) :: outdir                ! output directory for abinitio3D run products
                integer,      intent(in) :: mskdiam_in            ! mask diameter (A) used by 3D/ref-projection steps
                integer,     allocatable :: states(:), projs(:), nunique_proj(:), uniqbuf(:)
                type(commander_abinitio3D_cavgs) :: xabinitio3D_cavgs ! commander for abinitio3D with class-average-based rejection
                type(cmdline)            :: cline_reproject       ! command line builder for the reproject commander
                type(string)             :: cwd, volpath          ! saved working directory and selected volume path
                type(oris)               :: voloris               ! volume metadata container from project
                integer                  :: ldim(3)               ! selected volume box dimensions
                integer                  :: ldim_clip(3)          ! particle-stack box dimensions for clipping/padding
                integer                  :: ldim_new(3)           ! reprojection box dimensions after Fourier rescaling
                integer                  :: nvols, nuniq, i_cls3d, proj_here
                integer :: ivol, bestvol                          ! volume loop index and best-population volume id
                real    :: vol_smpd
                real    :: smpd_part                              ! particle-stack sampling distance (target)
                real    :: smpd_new                               ! sampling distance after rescaling
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to extract dir for partitioning
                call simple_chdir(outdir)
                call cline_abinitio3D%kill()
            !     call cline_abinitio3D%set('prg',  'abinitio3D_cavgs_reject')
            !     call cline_abinitio3D%set('mkdir',                     'no')
            !     !call cline_abinitio3D%set('pgrp',                      'c1')
            !     !call cline_abinitio3D%set('outdir',                  outdir)
            !     call cline_abinitio3D%set('nstates',              NSTATES3D)
            !     call cline_abinitio3D%set('lpstop',                LPSTOP2D)
            !  !   call cline_abinitio3D%set('nstages',              NSTAGES3D)
            !     call cline_abinitio3D%set('mskdiam',             mskdiam_in)
            ! !    call cline_abinitio3D%set('nsample_start',              200)
            !     call cline_abinitio3D%set('nthr',                         16)
            !  !   call cline_abinitio3D%set('nparts',                       4)
            !     call cline_abinitio3D%set('projfile',      cluster_projfile)

                call cline_abinitio3D%set('prg',  'abinitio3D_cavgs')
                call cline_abinitio3D%set('mkdir',                     'no')
                call cline_abinitio3D%set('pgrp',                      'c1')
                call cline_abinitio3D%set('nstates',              NSTATES3D)
                call cline_abinitio3D%set('lpstop',                8)
                call cline_abinitio3D%set('mskdiam',             mskdiam_in)
                call cline_abinitio3D%set('lpstart',              20)
                call cline_abinitio3D%set('prune',              'no')
                call cline_abinitio3D%set('nthr',                         16)
                call cline_abinitio3D%set('nstages',                       4)
                call cline_abinitio3D%set('projfile',      cluster_projfile)

                call cline_abinitio3D%printline()
                call xabinitio3D_cavgs%execute(cline_abinitio3D)

                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile) ! read the project with abinitio3D output
                if( spproj_inout%os_cls3D%isthere('state') .and. spproj_inout%os_cls3D%isthere('proj') ) then
                    states = spproj_inout%os_cls3D%get_all_asint('state')
                    projs  = spproj_inout%os_cls3D%get_all_asint('proj')
                    nvols = maxval(states)
                    if( size(states) == size(projs) )then
                        allocate(nunique_proj(nvols), source=0)
                        allocate(uniqbuf(max(1, size(projs))), source=0)
                        do ivol = 1, nvols
                            nuniq = 0
                            do i_cls3d = 1, size(states)
                                if( states(i_cls3d) /= ivol ) cycle
                                proj_here = projs(i_cls3d)
                                if( nuniq == 0 )then
                                    nuniq = 1
                                    uniqbuf(1) = proj_here
                                else if( .not. any(uniqbuf(1:nuniq) == proj_here) )then
                                    nuniq = nuniq + 1
                                    uniqbuf(nuniq) = proj_here
                                endif
                            enddo
                            nunique_proj(ivol) = nuniq
                        enddo
                        bestvol = maxloc(nunique_proj, 1)
                        write(logfhandle,'(A,I0,A,I0,A)') '>>> BEST VOLUME BY UNIQUE PROJ IN CLS3D: STATE=', bestvol, &
                            ' (NUNIQUE_PROJ=', nunique_proj(bestvol), ')'
                        deallocate(nunique_proj, uniqbuf)
                    else
                        write(logfhandle,'(A)') '>>> WARNING: cannot rank volumes by unique proj, nonconforming cls3D arrays'
                        write(logfhandle,'(A,I0)') '>>> FALLING BACK TO VOLUME STATE=', bestvol
                    endif
                else
                    write(logfhandle,'(A)') '>>> WARNING: missing cls3D state/proj or empty volume set; unique-proj ranking skipped'
                    write(logfhandle,'(A,I0)') '>>> FALLING BACK TO VOLUME STATE=', bestvol
                endif
                
                volpath = string('recvol_state'//int2str_pad(bestvol,2)//MRC_EXT)
                if( .not. file_exists(volpath) ) THROW_HARD('Expected abinitio3D output volume not found: '//volpath%to_char())
                call find_ldim_nptcls(volpath, ldim, nuniq)
                vol_smpd = find_img_smpd(volpath)
                ldim_clip(1) = spproj_inout%os_stk%get_int(1, 'box') ! particle box size
                ldim_clip(2) = ldim_clip(1)
                ldim_clip(3) = 1
                call cline_reproject%kill()
                call cline_reproject%set('prg',                   'reproject')
                call cline_reproject%set('vol1',                      volpath)
                call cline_reproject%set('smpd',                     vol_smpd)
                call cline_reproject%set('nspace',                         10)
                call cline_reproject%set('pgrp',                         'c1')
                call cline_reproject%set('nthr',                            8)
                call cline_reproject%set('mskdiam',              mskdiam_in)
                call cline_reproject%printline()
                call xreproject%execute(cline_reproject)
                call mrc2jpeg_tiled(string('reprojs.mrcs'), string('reprojs.jpeg'))
                ! rescale reprojections to the particle sampling/box and write as selected references
                smpd_part   = spproj_inout%os_stk%get(1, 'smpd')
                ldim_new(1) = round2even(real(ldim(1)) * vol_smpd / smpd_part)
                ldim_new(2) = ldim_new(1)
                ldim_new(3) = 1
                write(logfhandle,'(A,I0,A,I0,A)') '>>> RESCALING AND CLIPPING REPROJECTIONS TO ', ldim_new(1), ' PIXEL BOX (', ldim_clip(1), ' A) FOR PICKING REFERENCES'
                call scale_imgfile( string('reprojs.mrcs'), string('selected_references.mrcs'), vol_smpd, ldim_new, smpd_part)
                call voloris%kill()
                if( allocated(states) ) deallocate(states)
                if( allocated(projs)  ) deallocate(projs)
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile) ! copy abinitio3D output back to main projfile for downstream steps
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_abinitio3D_and_reproject

            ! Block until at least nmics micrographs pass quality thresholds
            ! (ctfres, icefrac, astig).  Polls project_buff for new partial
            ! projects and applies rejection after each batch.
            subroutine micimporter( nmics_target )
                integer, intent(in) :: nmics_target
                integer :: n_imported, n_new_oris, n_oris, iproj, iori_loc, imic
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
                            do iori_loc = 1, STREAM_NMOVS_SET
                                n_imported = n_imported + 1
                                ! set state=0 mics to state=-1
                                if( spproj_part%os_mic%get(iori_loc, 'state') < 1 ) call spproj_part%os_mic%set(iori_loc, 'state', -1)
                                call spproj%os_mic%transfer_ori(n_imported, spproj_part%os_mic, iori_loc)
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
                    !send progress update to GUI
                    if( n_oris > 0 ) call send_meta(string('importing micrographs'))
                    ! enough accepted micrographs — expose dead mics as state=0 and return
                    if( spproj%os_mic%count_state_gt_zero() >= nmics_target ) then
                        do imic = 1, n_oris
                            if( spproj%os_mic%get(imic, 'state') < 0 ) call spproj%os_mic%set(imic, 'state', 0)
                        enddo
                        return
                    endif
                end do
            end subroutine micimporter

            ! Evaluate class-average quality, map accepted classes, and emit GUI-ready cavgs.
            subroutine run_cavg_quality_selection(spproj_inout, cluster_projfile, outdir, cwd_cycle_in, mskdiam_inout, imgfiles_inout, smpd_stk_out, n_selected )
                type(sp_project),            intent(inout) :: spproj_inout
                type(string),                intent(in)    :: cluster_projfile
                type(string),                intent(in)    :: outdir
                type(string),                intent(in)    :: cwd_cycle_in
                real,                        intent(inout) :: mskdiam_inout
                type(string), allocatable,   intent(inout) :: imgfiles_inout(:)
                real,                        intent(out)   :: smpd_stk_out
                integer,                     intent(out)   :: n_selected
                type(string)                :: cavgsstk_local
                type(string)                :: cwd
                type(image), allocatable    :: cavg_imgs_local(:)
                type(cavg_quality_model)    :: model_local
                type(cavg_quality_result)   :: quality_local
                integer, allocatable        :: cavg_inds_local(:)
                integer                     :: ncls_local, i_mic_local, xtiles_local, ytiles_local
                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_copy_file(cluster_projfile, outdir//'/'//cluster_projfile) ! copy projfile to quality selection dir for processing
                call simple_chdir(outdir)
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
                call spproj_inout%get_cavgs_stk(cavgsstk_local, ncls_local, smpd_stk_out)
                if( allocated(imgfiles_inout) ) deallocate(imgfiles_inout)
                allocate(imgfiles_inout(1))
                imgfiles_inout(1) = cavgsstk_local
                cavg_imgs_local = read_cavgs_into_imgarr(spproj_inout)
                smpd_stk_out    = spproj_inout%get_smpd()
                call model_local%init_preset(CAVG_QUALITY_MODEL_POOL_DEFAULT)
                call evaluate_cavg_quality(cavg_imgs_local, spproj_inout%os_cls2D, mskdiam_inout, quality_local, model_local)
                call model_local%kill()
                n_selected = count(quality_local%states > 0)
                call write_quality_stack(string('quality_selected_cavgs'//MRC_EXT), cavg_imgs_local, quality_local%states, ncls_local, selected=.true.)
                call write_quality_stack(string('quality_rejected_cavgs'//MRC_EXT), cavg_imgs_local, quality_local%states, ncls_local, selected=.false.)
                call spproj_inout%map_cavgs_selection(quality_local%states)
                do i_mic_local = spproj_inout%os_mic%get_noris(), 1, -1
                    if( spproj_inout%os_mic%get(i_mic_local, 'state') == 0.0 .or. spproj_inout%os_mic%get(i_mic_local, 'nptcls') == 0.0 ) then
                        call spproj_inout%os_mic%delete(i_mic_local)
                    endif
                end do
                call spproj_inout%cavgs2jpg(cavg_inds_local, string('quality_cavgs')//JPG_EXT, xtiles_local, ytiles_local, ignore_states=.false.)
                if( allocated(cavg_inds_local) ) then
                    cavg_inds_local = pack(cavg_inds_local, cavg_inds_local > 0) ! filter out zero indices (unselected classes)
                    if( size(cavg_inds_local) > 0 ) then
                        call send_available_cavgs2D(cwd_cycle_in//'/quality_cavgs'//JPG_EXT, size(cavg_inds_local), &
                            cavg_inds_local, cavgsstk_local, xtiles_local, ytiles_local, spproj_inout)
                    endif
                    deallocate(cavg_inds_local)
                end if
                call dealloc_imgarr(cavg_imgs_local)
                call spproj_inout%write(cluster_projfile) ! write project with quality-selected classes mapped in, for downstream steps
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile)
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
            end subroutine run_cavg_quality_selection

            subroutine balance_classes( spproj_inout, cluster_projfile, outdir )
                type(sp_project), intent(inout) :: spproj_inout
                type(string), intent(in)         :: cluster_projfile
                type(string), intent(in)         :: outdir
                integer, parameter :: TARGET_NCLS = 501
                type(string)                 :: cavgsstk, balanced_stk, odd_stk, even_stk, sigma2_stk
                type(string)                 :: odd_balanced_stk, even_balanced_stk, sigma2_balanced_stk
                type(image), allocatable     :: cavg_imgs(:)
                type(oris)                   :: os_cls2D_src
                integer, allocatable         :: src_inds(:), pops(:), reps(:), extra(:)
                real,    allocatable         :: frac(:)
                logical, allocatable         :: picked(:)
                integer :: ncls_all, nsrc, total_pop, icls, j, out_ind, istk, rem, imax, iout
                integer :: n_balanced
                real    :: smpd_dummy
                type(string) :: cwd

                call simple_getcwd(cwd)
                call simple_mkdir(outdir)
                call simple_chdir(outdir)

                call spproj_inout%get_cavgs_stk(cavgsstk, ncls_all, smpd_dummy, out_ind=out_ind)
                if( ncls_all <= 0 ) then
                    call simple_chdir(cwd)
                    return
                endif

                allocate(src_inds(spproj_inout%os_cls2D%count_state_gt_zero()))
                nsrc = 0
                do icls = 1, spproj_inout%os_cls2D%get_noris()
                    if( nint(spproj_inout%os_cls2D%get(icls, 'state')) > 0 ) then
                        nsrc = nsrc + 1
                        src_inds(nsrc) = icls
                    endif
                enddo
                if( nsrc <= 0 ) then
                    deallocate(src_inds)
                    call simple_chdir(cwd)
                    return
                endif
                if( nsrc >= TARGET_NCLS ) then
                    deallocate(src_inds)
                    call simple_chdir(cwd)
                    return
                endif

                cavg_imgs = read_cavgs_into_imgarr(spproj_inout)
                call os_cls2D_src%new(ncls_all, is_ptcl=.false.)
                do icls = 1, ncls_all
                    call os_cls2D_src%transfer_ori(icls, spproj_inout%os_cls2D, icls)
                enddo

                allocate(pops(nsrc), reps(nsrc), extra(nsrc), frac(nsrc), picked(nsrc))
                pops = 1
                do icls = 1, nsrc
                    if( spproj_inout%os_cls2D%isthere(src_inds(icls), 'pop') ) then
                        pops(icls) = max(1, nint(spproj_inout%os_cls2D%get(src_inds(icls), 'pop')))
                    endif
                enddo
                total_pop = sum(pops)
                rem = TARGET_NCLS - nsrc
                frac = real(rem) * real(pops) / real(total_pop)
                extra = int(frac)
                reps = 1 + extra
                rem = rem - sum(extra)
                if( rem > 0 ) then
                    frac = frac - real(extra)
                    picked = .false.
                    do j = 1, rem
                        imax = maxloc(frac, dim=1, mask=.not. picked)
                        reps(imax) = reps(imax) + 1
                        picked(imax) = .true.
                    enddo
                endif

                n_balanced = sum(reps)
                balanced_stk = stemname(cavgsstk)//'_balanced'//MRC_EXT
                if( file_exists(balanced_stk) ) call del_file(balanced_stk)

                call spproj_inout%os_cls2D%new(n_balanced, is_ptcl=.false.)
                istk = 0
                do icls = 1, nsrc
                    do j = 1, reps(icls)
                        istk = istk + 1
                        call cavg_imgs(src_inds(icls))%write(balanced_stk, istk)
                        call spproj_inout%os_cls2D%transfer_ori(istk, os_cls2D_src, src_inds(icls))
                        call spproj_inout%os_cls2D%set(istk, 'indstk', istk)
                        call spproj_inout%os_cls2D%set(istk, 'state', 1.)
                    enddo
                enddo
                if( spproj_inout%os_cls3D%get_noris() /= n_balanced ) then
                    call spproj_inout%os_cls3D%new(n_balanced, is_ptcl=.false.)
                    do icls = 1, n_balanced
                        call spproj_inout%os_cls3D%transfer_ori(icls, spproj_inout%os_cls2D, icls)
                    enddo
                endif

                call spproj_inout%os_out%set(out_ind, 'stk',        simple_abspath(balanced_stk))
                call spproj_inout%os_out%set(out_ind, 'nptcls',     n_balanced)
                call spproj_inout%os_out%set(out_ind, 'nptcls_stk', n_balanced)
                call spproj_inout%os_out%set(out_ind, 'fromp',      1)
                call spproj_inout%os_out%set(out_ind, 'top',        n_balanced)

                ! Keep odd/even class-average stacks synchronized when present.
                odd_stk = swap_suffix(cavgsstk, string('_odd'//MRC_EXT), string(MRC_EXT))
                if( file_exists(odd_stk) ) then
                    odd_balanced_stk = stemname(odd_stk)//'_balanced_odd'//MRC_EXT
                    call duplicate_balanced_stack(odd_stk, odd_balanced_stk, ncls_all, nsrc, src_inds, reps)
                    call update_os_out_stk(spproj_inout, odd_stk, odd_balanced_stk, n_balanced)
                endif
                even_stk = swap_suffix(cavgsstk, string('_even'//MRC_EXT), string(MRC_EXT))
                if( file_exists(even_stk) ) then
                    even_balanced_stk = stemname(even_stk)//'_balanced_even'//MRC_EXT
                    call duplicate_balanced_stack(even_stk, even_balanced_stk, ncls_all, nsrc, src_inds, reps)
                    call update_os_out_stk(spproj_inout, even_stk, even_balanced_stk, n_balanced)
                endif

                ! Carry sigma2 output forward to a balanced companion file.
                do iout = 1, spproj_inout%os_out%get_noris()
                    if( .not. spproj_inout%os_out%isthere(iout, 'imgkind') ) cycle
                    if( spproj_inout%os_out%get_str(iout, 'imgkind') /= 'sigma2' ) cycle
                    if( .not. spproj_inout%os_out%isthere(iout, 'sigma2') ) cycle
                    sigma2_stk = spproj_inout%os_out%get_str(iout, 'sigma2')
                    if( .not. file_exists(sigma2_stk) ) cycle
                    sigma2_balanced_stk = stemname(sigma2_stk)//'_balanced'//STAR_EXT
                    call simple_copy_file(sigma2_stk, sigma2_balanced_stk)
                    call spproj_inout%os_out%set(iout, 'sigma2', simple_abspath(sigma2_balanced_stk))
                enddo

                call spproj_inout%write(cluster_projfile) ! write project with quality-selected classes mapped in, for downstream steps
                call simple_chdir(cwd)
                call simple_copy_file(outdir//'/'//cluster_projfile, cluster_projfile)
                call spproj_inout%kill()
                call spproj_inout%read(cluster_projfile)
                write(logfhandle, '(A,I0,A)') '>>> BALANCED CLASS AVERAGES TO ', n_balanced, ' ENTRIES'
                call os_cls2D_src%kill()
                call dealloc_imgarr(cavg_imgs)
                deallocate(src_inds, pops, reps, extra, frac, picked)
            end subroutine balance_classes

            subroutine duplicate_balanced_stack( stk_in, stk_out, ncls_all, nsrc, src_inds, reps )
                type(string), intent(in) :: stk_in, stk_out
                integer,      intent(in) :: ncls_all, nsrc
                integer,      intent(in) :: src_inds(:), reps(:)
                type(image), allocatable :: imgs_in(:)
                integer :: ii, jj, kk
                imgs_in = read_stk_into_imgarr(stk_in)
                if( size(imgs_in) /= ncls_all ) then
                    call dealloc_imgarr(imgs_in)
                    return
                endif
                if( file_exists(stk_out) ) call del_file(stk_out)
                kk = 0
                do ii = 1, nsrc
                    do jj = 1, reps(ii)
                        kk = kk + 1
                        call imgs_in(src_inds(ii))%write(stk_out, kk)
                    enddo
                enddo
                call dealloc_imgarr(imgs_in)
            end subroutine duplicate_balanced_stack

            subroutine update_os_out_stk( spproj_inout, old_stk, new_stk, nstk )
                type(sp_project), intent(inout) :: spproj_inout
                type(string),     intent(in)    :: old_stk, new_stk
                integer,          intent(in)    :: nstk
                type(string) :: stk_here
                integer      :: io
                do io = 1, spproj_inout%os_out%get_noris()
                    if( .not. spproj_inout%os_out%isthere(io, 'stk') ) cycle
                    stk_here = spproj_inout%os_out%get_str(io, 'stk')
                    if( simple_abspath(stk_here) /= simple_abspath(old_stk) ) cycle
                    call spproj_inout%os_out%set(io, 'stk',        simple_abspath(new_stk))
                    call spproj_inout%os_out%set(io, 'nptcls',     nstk)
                    call spproj_inout%os_out%set(io, 'nptcls_stk', nstk)
                    call spproj_inout%os_out%set(io, 'fromp',      1)
                    call spproj_inout%os_out%set(io, 'top',        nstk)
                enddo
            end subroutine update_os_out_stk

            ! Re-estimate mask diameter and extraction box from quality-selected cavgs.
            subroutine reestimate_box_size_from_selected_cavgs( spproj_inout, params_in, mskdiam_inout, box_in_pix_inout )
                type(sp_project), intent(inout) :: spproj_inout
                real,             intent(inout) :: mskdiam_inout
                integer,          intent(inout) :: box_in_pix_inout
                type(parameters),    intent(in) :: params_in
                type(image),        allocatable :: masks_local(:)
                type(image),        allocatable :: cavg_imgs_local(:)
                integer,            allocatable :: selected_cavg_inds_local(:)
                integer,            allocatable :: cls_states_local(:)
                real,               allocatable :: diams_inliers_local(:)
                real,               allocatable :: diams_local(:), shifts_local(:,:)
                integer                         :: i_mask_local, mskdiam_pix, ncls_local
                real                            :: mad_local, diam_max_local
                type(stats_struct)              :: diam_stats_local
                ncls_local      = spproj_inout%os_cls2D%get_noris()
                cavg_imgs_local = read_cavgs_into_imgarr(spproj_inout)
                cls_states_local = nint(spproj_inout%os_cls2D%get_all('state'))
                selected_cavg_inds_local = pack((/(i_mask_local, i_mask_local = 1, ncls_local)/), cls_states_local > 0)
                if( size(selected_cavg_inds_local) == 0 ) THROW_HARD('No quality-selected class averages for box-size re-estimation')
                allocate(masks_local(size(selected_cavg_inds_local)))
                do i_mask_local = 1, size(selected_cavg_inds_local)
                    call masks_local(i_mask_local)%copy(cavg_imgs_local(selected_cavg_inds_local(i_mask_local)))
                end do
                call automask2D(params_in, masks_local, params_in%ngrow, nint(params_in%winsz), params_in%edge, diams_local, shifts_local)

                call calc_stats(diams_local, diam_stats_local)
                mad_local = mad_gau(diams_local, diam_stats_local%med)
                if( mad_local > 0.0 ) then
                    diams_inliers_local = pack(diams_local, abs((diams_local - diam_stats_local%med) / mad_local) < 3.0)
                    if( size(diams_inliers_local) > 0 ) then
                        diam_max_local = maxval(diams_inliers_local)
                    else
                        diam_max_local = maxval(diams_local)
                    endif
                    deallocate(diams_inliers_local)
                else
                    diam_max_local = maxval(diams_local)
                endif

                mskdiam_pix      = round2even((diam_max_local / params_in%smpd + 2. * COSMSKHALFWIDTH) * MSK_EXP_FAC)
                mskdiam_inout    = params_in%smpd * mskdiam_pix
               ! box_in_pix_inout = find_larger_magic_box(round2even((mskdiam_pix + mskdiam_pix * BOX_EXP_FAC) / params_in%smpd))
                box_in_pix_inout = find_larger_magic_box(round2even((mskdiam_pix * 1.2)))
                write(logfhandle, '(A,I0,A)') '>>> UPDATED MASK DIAMETER TO ', round2even(mskdiam_inout), ' A'
                write(logfhandle, '(A,I0,A)') '>>> UPDATED BOX SIZE TO ', box_in_pix_inout, ' PIXELS'
                call dealloc_imgarr(cavg_imgs_local)
                deallocate(cls_states_local)
                deallocate(selected_cavg_inds_local)
                deallocate(masks_local)
                deallocate(diams_local, shifts_local)
            end subroutine reestimate_box_size_from_selected_cavgs

            subroutine write_quality_stack( fname, cavg_imgs_in, quality_states, ncls_in, selected )
                type(string), intent(in) :: fname
                type(image),  intent(inout) :: cavg_imgs_in(:)
                integer,      intent(in) :: quality_states(:)
                integer,      intent(in) :: ncls_in
                logical,      intent(in) :: selected
                integer :: icls, istk
                if( file_exists(fname) ) call del_file(fname)
                istk = 0
                do icls = 1, ncls_in
                    if( selected .eqv. (quality_states(icls) > 0) )then
                        istk = istk + 1
                        call cavg_imgs_in(icls)%write(fname, istk)
                    endif
                enddo
                write(logfhandle,'(A,A,A,I6)') '>>> WROTE ', fname%to_char(), ' #CAVGS: ', istk
            end subroutine write_quality_stack

            ! Broadcast initial-picking progress to the GUI.
            subroutine send_meta( my_stage )
                type(string), intent(in) :: my_stage
                integer                  :: my_ntarget
                my_ntarget = max(params%nmics, spproj%os_mic%count_state_gt_zero())
                call meta_initial_picking%set(                                   &
                    stage                = my_stage,                             &
                    micrographs_imported = my_ntarget,                           &
                    micrographs_accepted = spproj%os_mic%count_state_gt_zero(),  &
                    particles_extracted  = spproj%os_ptcl2D%get_noris(),         &
                    box_size             = box_in_pix)
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
                    path   = spproj%os_mic%get_str(my_ithumb, "thumb"),       & !// TODO - update to use thumb_den. also change caller test
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
            subroutine send_cavg2D_meta( my_path, my_i, my_i_delta, my_i_max, my_xtile, my_ytile, cavg_inds_in, cavgsstk_in, xtiles_in, ytiles_in, spproj_in )
                type(string), intent(in) :: my_path
                integer,      intent(in) :: my_i, my_i_delta, my_i_max, my_xtile, my_ytile
                integer,      intent(in) :: cavg_inds_in(:)
                type(string), intent(in) :: cavgsstk_in
                integer,      intent(in) :: xtiles_in, ytiles_in
                type(sp_project), intent(in) :: spproj_in
                integer                  :: my_idx
                my_idx = cavg_inds_in(my_i)
                call meta_cavg2D%set(                                       &
                    path       = my_path,                                   &
                    mrcpath    = cavgsstk_in,                               &
                    i          = my_i + my_i_delta,                         &
                    i_max      = my_i_max,                                  &
                    res        = spproj_in%os_cls2D%get(my_idx,     'res'), &
                    pop        = spproj_in%os_cls2D%get_int(my_idx, 'pop'), &
                    idx        = my_idx,                                    &
                    sprite     = sprite_sheet_pos(                          &
                                  x = merge(0.0, my_xtile * (100.0 / (xtiles_in - 1)), xtiles_in == 1), &
                                  y = merge(0.0, my_ytile * (100.0 / (ytiles_in - 1)), ytiles_in == 1), &
                                  h = 100 * ytiles_in,                                                  &
                                  w = 100 * xtiles_in))
                if( meta_cavg2D%assigned() .and. mq_stream_master_in%is_active() ) then
                    call meta_cavg2D%serialise(meta_buffer)
                    call mq_stream_master_in%send(meta_buffer)
                endif
            end subroutine send_cavg2D_meta

            ! Send a batch of cavg2D metadata to the GUI, resetting tile counters.
            subroutine send_available_cavgs2D( my_path, n, cavg_inds_in, cavgsstk_in, xtiles_in, ytiles_in, spproj_in, my_i_max, my_i_start )
                type(string),      intent(in) :: my_path
                integer,           intent(in) :: n
                integer,           intent(in) :: cavg_inds_in(:)
                type(string),      intent(in) :: cavgsstk_in
                integer,           intent(in) :: xtiles_in, ytiles_in
                type(sp_project),  intent(in) :: spproj_in
                integer, optional, intent(in) :: my_i_max, my_i_start
                integer                       :: my_xtile, my_ytile, my_i, my_i_delta
                my_xtile   = 0
                my_ytile   = 0
                my_i_delta = 0
                if( present(my_i_start) ) my_i_delta = my_i_start - 1
                if( present(my_i_max) ) then
                    write(logfhandle,*) '>>> SENDING', n, ' CLASS AVERAGES TO GUI', my_i_max
                else
                    write(logfhandle,*) '>>> SENDING', n, ' CLASS AVERAGES TO GUI'
                endif
                do my_i = 1, n
                    if (present(my_i_max)) then
                        call send_cavg2D_meta(my_path, my_i, my_i_delta, my_i_max, my_xtile, my_ytile, &
                            cavg_inds_in, cavgsstk_in, xtiles_in, ytiles_in, spproj_in)
                    else
                        call send_cavg2D_meta(my_path, my_i, my_i_delta, n, my_xtile, my_ytile, &
                            cavg_inds_in, cavgsstk_in, xtiles_in, ytiles_in, spproj_in)
                    endif
                    my_xtile = my_xtile + 1
                    if( my_xtile == xtiles_in ) then
                        my_xtile = 0
                        my_ytile = my_ytile + 1
                    endif
                end do
            end subroutine send_available_cavgs2D

            ! Called asynchronously on SIGTERM. Exits immediately after logging.
            subroutine sigterm_handler()
                write(logfhandle, '(A)') 'SIGTERM RECEIVED'
                call exit(0)
            end subroutine sigterm_handler

    end subroutine exec_stream_p03_initial_analysis

end module simple_stream_p03_initial_analysis
