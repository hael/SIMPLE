!@descr: state and steps of stream task 7 (multistate 3D): import the particles pool 2D exports, run solve3D and then solve3D_addon as they grow, report to the GUI
!==============================================================================
! MODULE: simple_stream_stage_solve3D
!
! PURPOSE:
!   The body of stream p07 as a type. Pool 2D publishes its classified state
!   after each completed iteration (doc/policies/stream/stream_3D_ingestion_policy.md).
!   Each pass in which no 3D job runs takes the newest publication: its class
!   averages are selected once, and the selection and 2D parameters are merged
!   into the stage's rows. The first particles start solve3D; once it is
!   done, every growth of the rows starts an solve3D_addon run from the
!   latest result. After each run the stage's project holds the result, and
!   the GUI gets each state's volume, resolution, reprojections and orientation
!   distribution. The commander (simple_commanders_stream_p07_solve3D_multistate)
!   only normalises the command line and loops over iterate() until finished().
!
!   What happens is delegated:
!     - class-average selection -> simple_cavg_quality_selection
!     - 3D jobs                 -> simple_qsys_async_job (solve3D, solve3D_addon)
!     - GUI                     -> simple_stream_pipe, simple_stream_gui_senders,
!                                  gui_metadata_vol3D%set_oridist_from_oris
!
!   The stage's rows only ever grow, as solve3D_addon requires (rows are
!   read by index, never renumbered): a publication's stacks are matched to
!   the rows by stack name, so the order of a restarted pool does not matter.
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! RESTART:
!   A restart removes a leftover TERM_STREAM and starts again from the newest
!   publication; solve3D runs again in its folder.
!==============================================================================
module simple_stream_stage_solve3D
use simple_defs,                                      only: logfhandle, COSMSKHALFWIDTH
use simple_defs_fname,                                only: TERM_STREAM, METADATA_EXT, MRC_EXT, JPG_EXT, PPROC_SUFFIX,&
                                                           &LP_SUFFIX, MIRR_SUFFIX
use simple_defs_stream,                               only: DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME
use simple_defs_environment,                          only: SIMPLE_STREAM_PREPROC_PARTITION
use simple_error,                                     only: simple_exception
use simple_string,                                    only: string
use simple_string_utils,                              only: int2str, int2str_pad, lex_sort
use simple_fileio,                                    only: add2fbody, basename, del_file, file2rarr, file_exists, get_fbody,&
                                                           &get_fpath, simple_abspath, simple_getcwd
use simple_syslib,                                    only: dir_exists, simple_mkdir
use simple_math,                                      only: round2even
use simple_math_ft,                                   only: get_resarr
use simple_estimate_ssnr,                             only: get_resolution
use simple_imghead,                                   only: get_mrc_minmax
use simple_refine3D_fnames,                           only: refine3D_oris_heatmap_fname
use simple_cmdline,                                   only: cmdline
use simple_parameters,                                only: parameters
use simple_sp_project,                                only: sp_project
use simple_image,                                     only: image
use simple_imgarr_utils,                              only: dealloc_imgarr
use simple_qsys_env,                                  only: qsys_env
use simple_qsys_funs,                                 only: qsys_cleanup
use simple_qsys_async_job,                            only: qsys_async_job, ASYNC_JOB_RUNNING, ASYNC_JOB_DONE, ASYNC_JOB_FAILED
use simple_rec_list,                                  only: rec_list, rec_iterator, chunk_rec
use simple_stream_watcher,                            only: stream_watcher
use simple_stream_state,                              only: ipc_pipe_solve3D_multistate_in
use simple_stream_utils,                              only: create_stream_project, init_stream_qenv
use simple_cavg_quality_model,                        only: cavg_quality_model, CAVG_QUALITY_MODEL_POOL_DEFAULT
use simple_cavg_quality_types,                        only: cavg_quality_result
use simple_cavg_quality_selection,                    only: score_project_cavgs, write_cavg_selection_stacks
use simple_gui_utils,                                 only: mrc2jpeg_tiled
use simple_gui_metadata_utils,                        only: max_metadata_size
use simple_gui_metadata_types,                        only: GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE, GUI_METADATA_VOL3D_TYPE,&
                                                           &GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE
use simple_gui_metadata_cavg2D,                       only: gui_metadata_cavg2D
use simple_gui_metadata_vol3D,                        only: gui_metadata_vol3D
use simple_gui_metadata_stream_solve3D_multistate, only: gui_metadata_stream_solve3D_multistate
use simple_stream_pipe,                               only: stream_pipe
use simple_stream_gui_senders,                        only: send_reproj_tiles
implicit none

public :: stream_stage_solve3D
public :: PHASE_IMPORTING, PHASE_SOLVE3D, PHASE_IDLE, PHASE_ADDON
public :: JOB_NONE, JOB_SOLVE3D, JOB_ADDON
private
#include "simple_local_flags.inc"

! phases of the stage
integer, parameter :: PHASE_IMPORTING = 0 ! no 3D yet
integer, parameter :: PHASE_SOLVE3D  = 1 ! solve3D runs
integer, parameter :: PHASE_IDLE      = 2 ! a result exists; waiting for more particles
integer, parameter :: PHASE_ADDON     = 3 ! solve3D_addon runs
! the job to start next (next_job)
integer, parameter :: JOB_NONE        = 0
integer, parameter :: JOB_SOLVE3D    = 1
integer, parameter :: JOB_ADDON       = 2

character(len=*), parameter :: SOLVE3D_DIR     = 'solve3D'
character(len=*), parameter :: ADDON_DIR        = 'solve3D_addon'
character(len=*), parameter :: QUALITY_DIR      = 'quality_selection'
character(len=*), parameter :: SELECTED_CAVGS   = 'quality_selected_cavgs'
character(len=*), parameter :: REJECTED_CAVGS   = 'quality_rejected_cavgs'

! Components and steps are public so simple_stream_stage_solve3D_tester can assemble a stage
! and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_solve3D
    type(parameters)                                :: params
    type(sp_project)                                :: spproj        ! the pool: imported particles, then each run's result
    type(qsys_env)                                  :: qenv
    type(qsys_async_job)                            :: job           ! the running 3D job
    type(stream_watcher)                            :: project_buff  ! the publications of pool 2D
    type(rec_list)                                  :: setslist      ! one record per publication; included once taken or passed over
    type(stream_pipe)                               :: pipe          ! to the master
    type(gui_metadata_stream_solve3D_multistate) :: meta_status
    type(gui_metadata_cavg2D)                       :: meta_reproj
    type(string),   allocatable :: stk_names(:)    ! the stacks in the pool
    real,           allocatable :: state_res(:)    ! FSC=0.143 resolution per state of the latest run (0: none)
    type(string)                :: frozen_projfile ! the latest run's project, which the next addon run builds on
    integer :: phase              = PHASE_IMPORTING
    integer :: naddon_runs        = 0
    integer :: nptcls_at_last_run = 0 ! particles in the pool when the latest run started
    integer :: nptcls_selected    = 0 ! selected particles in the rows
    real    :: mskdiam            = 0. ! pool 2D's mask diameter (A), from its first export
    logical :: l_mskdiam_read     = .false.
    logical :: l_restart          = .false.
    logical :: l_attached         = .false. ! pool 2D's completed folder exists and is watched
    logical :: l_waiting_logged   = .false.
    logical :: l_exists           = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s           = SHORTWAIT ! an export is taken once untouched longer than this;
                                              ! pool 2D writes exports with a rename
    integer :: wait_s             = WAITTIME  ! pause at the end of a pass
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    ! the steps of new(), in order
    procedure :: init_params
    procedure :: init_queue
    procedure :: init_gui
    ! the steps of iterate() and their helpers
    procedure :: attach_upstream
    procedure :: watch_sets
    procedure :: read_mskdiam
    procedure :: import_sets
    procedure :: select_cavgs
    procedure :: merge_publication
    procedure :: stack_index
    procedure :: advance_jobs
    procedure :: start_solve3D
    procedure :: start_addon
    procedure :: finish_run
    procedure :: write_stage_project
    ! GUI
    procedure :: send_status
    procedure :: send_volumes
    procedure :: read_state_fsc
    ! rules that use no stage state
    procedure, nopass :: next_job
    procedure, nopass :: fit_mskdiam
end type stream_stage_solve3D

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline, the queue environment of the 3D jobs
    !! and the GUI.
    subroutine new( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        call self%kill()
        call self%init_params(cline)
        self%l_exists = .true. ! from here kill() releases what has been built
        call self%init_queue()
        call self%init_gui(ipc_pipe_solve3D_multistate_in(2))
        call self%send_status()
    end subroutine new

    !> The stage's project file (made and given a computing environment) and its parameters; the
    !! project must start without micrographs. A restart's leftover TERM_STREAM is removed in the
    !! stage's folder, where the loop looks for it.
    subroutine init_params( self, cline )
        class(stream_stage_solve3D), intent(inout) :: self
        class(cmdline),                 intent(inout) :: cline
        type(string) :: outdir, dir_exec
        ! a restart is recognised by its output directory (or dir_exec), before params%new makes one
        self%l_restart = .false.
        if( cline%defined('outdir') )then
            outdir = cline%get_carg('outdir')
            if( outdir%strlen() > 0 ) self%l_restart = dir_exists(outdir)
        endif
        if( cline%defined('dir_exec') )then
            dir_exec = cline%get_carg('dir_exec')
            if( .not. file_exists(dir_exec) ) THROW_HARD('Previous directory does not exist: '//dir_exec%to_char())
            self%l_restart = .true.
        endif
        call create_stream_project(self%spproj, cline, string('3Dmultistate'))
        call self%params%new(cline)
        if( self%l_restart )then
            write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
            call del_file(TERM_STREAM)
        endif
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 )then
            THROW_HARD('commander_stream_p07_solve3D_multistate must start from an empty project (e.g. from root project folder)')
        endif
        allocate(self%stk_names(0))
        allocate(self%state_res(self%params%nstates), source=0.)
    end subroutine init_params

    !> The queue environment the 3D jobs are submitted through; on the preprocessing partition,
    !! as before.
    subroutine init_queue( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call init_stream_qenv(self%params, self%qenv, string(SIMPLE_STREAM_PREPROC_PARTITION))
    end subroutine init_queue

    !> The GUI metadata objects and the pipe end to the master (-1: none).
    subroutine init_gui( self, fd_write )
        class(stream_stage_solve3D), intent(inout) :: self
        integer,                        intent(in)    :: fd_write
        call self%meta_status%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_TYPE)
        call self%meta_reproj%new(GUI_METADATA_STREAM_SOLVE3D_MULTISTATE_REPROJ_TYPE)
        call self%pipe%new(-1, fd_write, max_metadata_size(), 'solve3D_multistate')
    end subroutine init_gui

    !> One pass: wait for pool 2D's folder; take new exports (unless a job runs); advance the 3D
    !! jobs; report.
    subroutine iterate( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call self%send_status()
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%watch_sets()
        if( self%phase /= PHASE_SOLVE3D .and. self%phase /= PHASE_ADDON ) call self%import_sets()
        call self%advance_jobs()
        call self%send_status()
        call sleep(self%wait_s)
    end subroutine iterate

    !> .true. once the stream is told to stop.
    logical function finished( self )
        class(stream_stage_solve3D), intent(in) :: self
        finished = file_exists(TERM_STREAM)
    end function finished

    !> The last status, and the stage's project with the latest result and the particles imported
    !! since. A running job is left to finish.
    subroutine finalize( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call self%meta_status%set_user_input(.false.)
        call self%send_status(string('terminating'))
        if( self%spproj%os_ptcl3D%get_noris() > 0 ) call self%write_stage_project()
        call qsys_cleanup(self%params)
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_solve3D), intent(inout) :: self
        if( .not. self%l_exists ) return
        call self%spproj%kill
        call self%qenv%kill
        call self%job%kill
        call self%project_buff%kill
        call self%setslist%kill
        call self%pipe%kill
        call self%meta_status%kill
        call self%meta_reproj%kill
        call self%frozen_projfile%kill
        if( allocated(self%stk_names) ) deallocate(self%stk_names)
        if( allocated(self%state_res) ) deallocate(self%state_res)
        self%phase              = PHASE_IMPORTING
        self%naddon_runs        = 0
        self%nptcls_at_last_run = 0
        self%nptcls_selected    = 0
        self%mskdiam            = 0.
        self%l_mskdiam_read     = .false.
        self%l_restart          = .false.
        self%l_attached         = .false.
        self%l_waiting_logged   = .false.
        self%l_exists           = .false.
    end subroutine kill

    !---------------- steps ----------------

    ! Starts watching pool 2D's completed folder once it exists.
    subroutine attach_upstream( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: completed
        logical      :: l_ready
        completed = self%params%dir_target//'/'//DIR_STREAM_COMPLETED
        l_ready   = dir_exists(self%params%dir_target)
        if( l_ready ) l_ready = dir_exists(completed)
        if( .not. l_ready )then
            if( .not. self%l_waiting_logged )then
                write(logfhandle,'(A)') '>>> WAITING FOR '//completed%to_char()//' TO BE GENERATED'
                self%l_waiting_logged = .true.
            endif
            return
        endif
        write(logfhandle,'(A)') '>>> '//completed%to_char()//' FOUND'
        self%project_buff     = stream_watcher(self%settle_s, simple_abspath(completed), spproj=.true., nretries=10)
        self%l_attached       = .true.
        self%l_waiting_logged = .false.
    end subroutine attach_upstream

    ! One record per new export, in export order (the zero-padded names sort that way); the first
    ! gives pool 2D's mask diameter.
    subroutine watch_sets( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), allocatable :: projects(:)
        integer :: nprojects, i
        call self%project_buff%watch(nprojects, projects)
        if( nprojects == 0 ) return
        call lex_sort(projects)
        call self%project_buff%add2history(projects)
        do i = 1,nprojects
            call self%setslist%push2chunk_list(projects(i), self%setslist%size() + 1, .true.)
        enddo
        if( .not. self%l_mskdiam_read ) call self%read_mskdiam(projects(1))
    end subroutine watch_sets

    ! The mask diameter pool 2D records with the class averages of its exports.
    subroutine read_mskdiam( self, projfile )
        class(stream_stage_solve3D), intent(inout) :: self
        class(string),                  intent(in)    :: projfile
        type(sp_project) :: spproj
        type(string)     :: stk
        real             :: smpd
        integer          :: ncls
        call spproj%read_segment('out', projfile)
        call spproj%get_cavgs_stk(stk, ncls, smpd, fail=.false.)
        if( ncls > 0 ) call spproj%get_mskdiam('cavg', self%mskdiam)
        call spproj%kill
        self%l_mskdiam_read = .true.
        if( self%mskdiam > 0. )then
            write(logfhandle,'(A,F8.2)') '>>> MASK DIAMETER SET TO : ', self%mskdiam
        else
            write(logfhandle,'(A)') '>>> NO MASK DIAMETER IN THE EXPORTS; THE BOX DEFAULT IS USED'
        endif
    end subroutine read_mskdiam

    ! The newest publication not yet taken: its class averages are selected and it is merged into
    ! the rows. Older publications not taken are passed over: the newest holds what they held.
    subroutine import_sets( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(sp_project)   :: set
        type(rec_iterator) :: it
        type(chunk_rec)    :: crec, newest
        type(string)       :: stem
        integer            :: irec, inewest
        if( self%setslist%size() == 0 ) return
        ! the records are in publication order: the last one not taken is the newest
        inewest = 0
        it = self%setslist%begin()
        do irec = 1,self%setslist%size()
            call it%get(crec)
            if( .not. crec%included )then
                inewest = irec
                newest  = crec
            endif
            call it%next()
        enddo
        if( inewest == 0 ) return
        call set%read(newest%projfile)
        stem = get_fbody(basename(newest%projfile), METADATA_EXT, separator=.false.)
        call self%select_cavgs(set, stem)
        call self%merge_publication(set, newest%id)
        call set%kill
        ! the newest and every older publication are taken
        it = self%setslist%begin()
        do irec = 1,inewest
            call it%get(crec)
            if( .not. crec%included )then
                crec%included = .true.
                call self%setslist%replace_iterator(it, crec)
            endif
            call it%next()
        enddo
    end subroutine import_sets

    ! The class averages of @p set scored with the pool model, the selection mapped to its
    ! particles, and the selected and rejected class averages written (with JPEGs) in
    ! quality_selection/<stem>/. The relation parameters get the class averages' pixel size and a
    ! mask diameter that fits their box, as the model_cavgs_rejection commander gives them.
    subroutine select_cavgs( self, set, stem )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),              intent(inout) :: set
        class(string),                  intent(in)    :: stem
        type(cavg_quality_model)  :: model
        type(cavg_quality_result) :: quality
        type(image), allocatable  :: cavg_imgs(:)
        type(string)              :: stk, dir, fname
        real                      :: smpd, mskdiam, box_stk
        integer                   :: ncls
        box_stk = 0.
        call set%get_cavgs_stk(stk, ncls, smpd, box=box_stk) ! the os_out entry's box, as a real
        if( set%os_cls2D%get_noris() == 0 ) THROW_HARD('no cls2D entries in the export '//stem%to_char())
        mskdiam = fit_mskdiam(self%mskdiam, nint(box_stk), smpd)
        call model%init_preset(CAVG_QUALITY_MODEL_POOL_DEFAULT)
        call score_project_cavgs(set, model, mskdiam, cavg_imgs, quality, smpd=smpd)
        if( size(cavg_imgs) /= set%os_cls2D%get_noris() ) THROW_HARD('# class averages /= # cls2D entries in '//stem%to_char())
        dir = string(QUALITY_DIR//'/')//stem
        call simple_mkdir(QUALITY_DIR)
        call simple_mkdir(dir)
        call write_cavg_selection_stacks(cavg_imgs, quality%states, dir//'/'//SELECTED_CAVGS//MRC_EXT,&
            &dir//'/'//REJECTED_CAVGS//MRC_EXT)
        call mrc2jpeg_tiled(dir//'/'//SELECTED_CAVGS//MRC_EXT, dir//'/'//SELECTED_CAVGS//JPG_EXT)
        call mrc2jpeg_tiled(dir//'/'//REJECTED_CAVGS//MRC_EXT, dir//'/'//REJECTED_CAVGS//JPG_EXT)
        fname = simple_abspath(dir//'/'//SELECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY SELECTED CLASS AVERAGES JPEG ', fname%to_char()
        fname = simple_abspath(dir//'/'//REJECTED_CAVGS//JPG_EXT, check_exists=.false.)
        if( file_exists(fname) ) write(logfhandle,'(A,A)') '>>> QUALITY REJECTED CLASS AVERAGES JPEG ', fname%to_char()
        write(logfhandle,'(A,I6,A,I6)') '>>> CAVG QUALITY SELECTED / REJECTED : ', count(quality%states > 0), ' / ',&
            &count(quality%states <= 0)
        call set%map_cavgs_selection(quality%states)
        call dealloc_imgarr(cavg_imgs)
        call quality%kill
        call model%kill
    end subroutine select_cavgs

    ! Merges the publication @p set (number @p id) into the stage's rows, which only ever grow
    ! (solve3D_addon reads them by index). A stack the stage holds is matched by name: its
    ! particles take the publication's 2D parameters and selection in place and keep their 3D
    ! parameters, CTF and optics group. A new stack is appended with its micrograph and particles
    ! (2D and 3D). A stack the publication lacks keeps its rows, deselected. The classes are the
    ! publication's.
    subroutine merge_publication( self, set, id )
        class(stream_stage_solve3D), intent(inout) :: self
        class(sp_project),              intent(inout) :: set
        integer,                        intent(in)    :: id
        integer,      allocatable :: kstk(:)   ! the stage's stack of each publication stack, 0 when new
        logical,      allocatable :: l_seen(:) ! the stage's stacks the publication holds
        type(string), allocatable :: new_names(:)
        type(string) :: name
        integer :: nstks, nstks_pool, jstk, k, hint, nnew, nptcls_new, nptcls, i, iptcl, jptcl
        integer :: fromp, fromp_set, pool_nmics, pool_nptcls, imic, nupdated, ndeselected, s
        nstks = set%os_stk%get_noris()
        if( set%os_mic%get_noris() /= nstks ) THROW_HARD('# micrographs /= # stacks in publication '//int2str(id))
        nstks_pool = size(self%stk_names)
        allocate(kstk(nstks), source=0)
        allocate(l_seen(nstks_pool), source=.false.)
        ! match the stacks; a restarted pool may list them in another order
        hint       = 0
        nnew       = 0
        nptcls_new = 0
        do jstk = 1,nstks
            name = set%os_stk%get_str(jstk, 'stk')
            k    = self%stack_index(name, hint + 1)
            kstk(jstk) = k
            if( k > 0 )then
                l_seen(k) = .true.
                hint      = k
            else
                nnew       = nnew + 1
                nptcls_new = nptcls_new + set%os_stk%get_top(jstk) - set%os_stk%get_fromp(jstk) + 1
            endif
        enddo
        ! the stacks the stage holds: 2D parameters and selection in place
        nupdated = 0
        do jstk = 1,nstks
            k = kstk(jstk)
            if( k == 0 ) cycle
            fromp_set = set%os_stk%get_fromp(jstk)
            nptcls    = set%os_stk%get_top(jstk) - fromp_set + 1
            fromp     = self%spproj%os_stk%get_fromp(k)
            if( self%spproj%os_stk%get_top(k) - fromp + 1 /= nptcls ) THROW_HARD('a stack changed size: '//self%stk_names(k)%to_char())
            do i = 0,nptcls - 1
                iptcl = fromp     + i
                jptcl = fromp_set + i
                s     = set%os_ptcl2D%get_state(jptcl)
                call self%spproj%os_ptcl2D%transfer_2Dparams(iptcl, set%os_ptcl2D, jptcl)
                call self%spproj%os_ptcl2D%set_state(iptcl, s)
                call self%spproj%os_ptcl3D%set_state(iptcl, s)
            enddo
            nupdated = nupdated + nptcls
        enddo
        ! the stacks the publication lacks: kept, deselected
        ndeselected = 0
        do k = 1,nstks_pool
            if( l_seen(k) ) cycle
            fromp = self%spproj%os_stk%get_fromp(k)
            do iptcl = fromp,self%spproj%os_stk%get_top(k)
                call self%spproj%os_ptcl2D%set_state(iptcl, 0)
                call self%spproj%os_ptcl3D%set_state(iptcl, 0)
                ndeselected = ndeselected + 1
            enddo
        enddo
        ! the new stacks: appended, with their micrographs and particles
        if( nnew > 0 )then
            pool_nmics  = self%spproj%os_mic%get_noris()
            pool_nptcls = self%spproj%os_ptcl3D%get_noris()
            if( pool_nmics == 0 )then
                call self%spproj%os_mic%new(nnew,          is_ptcl=.false.)
                call self%spproj%os_stk%new(nnew,          is_ptcl=.false.)
                call self%spproj%os_ptcl2D%new(nptcls_new, is_ptcl=.true.)
                call self%spproj%os_ptcl3D%new(nptcls_new, is_ptcl=.true.)
                fromp = 1
            else
                call self%spproj%os_mic%reallocate(pool_nmics + nnew)
                call self%spproj%os_stk%reallocate(pool_nmics + nnew)
                call self%spproj%os_ptcl2D%reallocate(pool_nptcls + nptcls_new)
                call self%spproj%os_ptcl3D%reallocate(pool_nptcls + nptcls_new)
                fromp = self%spproj%os_stk%get_top(pool_nmics) + 1
            endif
            allocate(new_names(nnew))
            imic = pool_nmics
            do jstk = 1,nstks
                if( kstk(jstk) /= 0 ) cycle
                fromp_set = set%os_stk%get_fromp(jstk)
                nptcls    = set%os_stk%get_top(jstk) - fromp_set + 1
                imic      = imic + 1
                new_names(imic - pool_nmics) = set%os_stk%get_str(jstk, 'stk')
                call self%spproj%os_mic%transfer_ori(imic, set%os_mic, jstk)
                call self%spproj%os_stk%transfer_ori(imic, set%os_stk, jstk)
                call self%spproj%os_stk%set(imic, 'fromp', fromp)
                call self%spproj%os_stk%set(imic, 'top',   fromp + nptcls - 1)
                do i = 0,nptcls - 1
                    iptcl = fromp     + i
                    jptcl = fromp_set + i
                    call self%spproj%os_ptcl2D%transfer_ori(iptcl, set%os_ptcl2D, jptcl)
                    call self%spproj%os_ptcl2D%set_stkind(iptcl, imic)
                    call self%spproj%os_ptcl3D%transfer_ori(iptcl, set%os_ptcl3D, jptcl)
                    call self%spproj%os_ptcl3D%set_stkind(iptcl, imic)
                enddo
                fromp = fromp + nptcls
            enddo
            self%stk_names = [self%stk_names, new_names]
            deallocate(new_names)
        endif
        ! the publication's classes
        self%spproj%os_cls2D = set%os_cls2D
        self%nptcls_selected = self%spproj%os_ptcl3D%count_state_gt_zero()
        write(logfhandle,'(A,I6,A,I8,A,I8,A,I8,A,I8,A)') '>>> PUBLICATION ', id, ': ', nupdated, ' PARTICLES UPDATED, ',&
            &nptcls_new, ' APPENDED, ', ndeselected, ' DESELECTED (STACK NOT PUBLISHED); ', self%nptcls_selected, ' SELECTED'
        deallocate(kstk, l_seen)
    end subroutine merge_publication

    ! The stage's stack named @p name, 0 when it has none. The search starts at @p hint and wraps
    ! around: a publication lists the stacks in the order the stage holds them, so the next match
    ! is usually the next stack.
    integer function stack_index( self, name, hint )
        class(stream_stage_solve3D), intent(in) :: self
        class(string),                  intent(in) :: name
        integer,                        intent(in) :: hint
        integer :: n, i, k
        stack_index = 0
        n = size(self%stk_names)
        if( n == 0 ) return
        do i = 0,n - 1
            k = mod(max(hint, 1) - 1 + i, n) + 1
            if( self%stk_names(k) == name )then
                stack_index = k
                return
            endif
        enddo
    end function stack_index

    ! A running job is checked: done, its result becomes the pool; failed, the stage stops. Without
    ! a job, next_job decides whether to start solve3D or an addon run.
    subroutine advance_jobs( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: logfile
        integer      :: nptcls
        if( self%phase == PHASE_SOLVE3D .or. self%phase == PHASE_ADDON )then
            select case(self%job%status())
                case(ASYNC_JOB_RUNNING)
                    return
                case(ASYNC_JOB_FAILED)
                    logfile = self%job%get_log()
                    THROW_HARD('the 3D job failed; see '//logfile%to_char())
                case(ASYNC_JOB_DONE)
                    call self%finish_run()
            end select
            return
        endif
        nptcls = self%spproj%os_ptcl3D%get_noris()
        select case(next_job(self%phase, nptcls, self%nptcls_at_last_run))
            case(JOB_SOLVE3D)
                call self%start_solve3D()
            case(JOB_ADDON)
                call self%start_addon()
        end select
    end subroutine advance_jobs

    ! solve3D on the pool, in solve3D/.
    subroutine start_solve3D( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address
        call simple_getcwd(cwd)
        dir = cwd//'/'//SOLVE3D_DIR
        call simple_mkdir(dir)
        call self%spproj%write(dir//'/'//SOLVE3D_DIR//METADATA_EXT)
        call cline_job%set('prg',             'solve3D')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        SOLVE3D_DIR//METADATA_EXT)
        call cline_job%set('pgrp',            'c1')
        call cline_job%set('nstates',         self%params%nstates)
        call cline_job%set('nstages',         self%params%nstages)
        call cline_job%set('lpstart',         self%params%lpstart)
        call cline_job%set('lpstop',          self%params%lpstop)
        call cline_job%set('force_lp_range',  'yes')
        if( self%mskdiam > 0. ) call cline_job%set('mskdiam', self%mskdiam)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        call cline_job%printline()
        self%nptcls_at_last_run = self%spproj%os_ptcl3D%get_noris()
        call self%job%start(self%qenv, cline_job, dir, SOLVE3D_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_SOLVE3D
        call cline_job%kill
    end subroutine start_solve3D

    ! solve3D_addon on the grown pool, from the latest result, in solve3D_addon/it_<n>/.
    subroutine start_addon( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(cmdline) :: cline_job
        type(string)  :: cwd, dir, server_address
        self%naddon_runs = self%naddon_runs + 1
        write(logfhandle,'(A,I0)') '>>> ENTERING ADDON STAGE ', self%naddon_runs
        call simple_getcwd(cwd)
        call simple_mkdir(cwd//'/'//ADDON_DIR)
        dir = cwd//'/'//ADDON_DIR//'/it_'//int2str(self%naddon_runs)
        call simple_mkdir(dir)
        call self%spproj%write(dir//'/'//ADDON_DIR//METADATA_EXT)
        call cline_job%set('prg',             'solve3D_addon')
        call cline_job%set('mkdir',           'no')
        call cline_job%set('projfile',        ADDON_DIR//METADATA_EXT)
        call cline_job%set('projfile_frozen', self%frozen_projfile)
        call cline_job%set('nparts',          self%params%nparts3D)
        call cline_job%set('nthr',            self%params%nthr3D)
        call cline_job%set('worker_priority', 'high')
        server_address = self%qenv%get_persistent_worker_server_address()
        if( server_address%strlen() > 0 ) call cline_job%set('worker_server', server_address)
        call cline_job%printline()
        self%nptcls_at_last_run = self%spproj%os_ptcl3D%get_noris()
        call self%job%start(self%qenv, cline_job, dir, ADDON_DIR, exec_bin=string('simple_exec'))
        self%phase = PHASE_ADDON
        call cline_job%kill
    end subroutine start_addon

    ! The finished run's project becomes the pool and the stage's project, the base of the next
    ! addon run, and goes to the GUI state by state.
    subroutine finish_run( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string) :: dir
        dir = self%job%get_dir()
        if( self%phase == PHASE_SOLVE3D )then
            self%frozen_projfile = dir//'/'//SOLVE3D_DIR//METADATA_EXT
        else
            self%frozen_projfile = dir//'/'//ADDON_DIR//METADATA_EXT
        endif
        if( .not. file_exists(self%frozen_projfile) ) THROW_HARD('no project from the 3D job: '//self%frozen_projfile%to_char())
        call self%spproj%kill
        call self%spproj%read(self%frozen_projfile)
        call self%job%kill
        self%phase = PHASE_IDLE
        call self%write_stage_project()
        call self%send_volumes()
    end subroutine finish_run

    ! The pool as the stage's project, written as a temporary file and renamed.
    subroutine write_stage_project( self )
        class(stream_stage_solve3D), intent(inout) :: self
        call self%spproj%write(self%params%projfile, tempfile=.true.)
    end subroutine write_stage_project

    !---------------- GUI ----------------

    ! The phase, the run count, the particle counts and, once a run is done, each state's
    ! population and resolution.
    subroutine send_status( self, stage )
        class(stream_stage_solve3D), intent(inout) :: self
        type(string), optional,         intent(in)    :: stage
        type(string) :: stage_here
        integer      :: istate
        select case(self%phase)
            case(PHASE_IMPORTING)
                stage_here = 'importing particles'
                if( .not. self%l_attached ) stage_here = 'waiting on pool 2D'
            case(PHASE_SOLVE3D)
                stage_here = 'running solve3D'
            case(PHASE_ADDON)
                stage_here = 'running refine3D'
            case default
                stage_here = 'idle'
        end select
        if( present(stage) ) stage_here = stage
        call self%meta_status%set(stage=stage_here, solve3D_stage=min(self%phase, PHASE_IDLE),&
            &refine_iteration=self%naddon_runs, nstates=self%params%nstates,&
            &particles_imported=self%spproj%os_ptcl3D%get_noris(), particles_at_last_refine=self%nptcls_at_last_run,&
            &resolution=0.)
        if( self%frozen_projfile%strlen() > 0 .and. allocated(self%state_res) )then
            do istate = 1,self%params%nstates
                call self%meta_status%set_state_stats(istate, self%spproj%os_ptcl3D%get_pop(istate, 'state'),&
                    &self%state_res(istate))
            enddo
        endif
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

    ! Each state with a volume in the pool's project: its products, population, resolution and FSC
    ! curve, the minimum and maximum of each volume, the orientation distribution, and its three
    ! reprojections as tiles. The resolutions are kept for the status.
    subroutine send_volumes( self )
        class(stream_stage_solve3D), intent(inout) :: self
        type(gui_metadata_vol3D) :: meta, fresh_meta
        type(string)             :: volpath, pprocpath, lppath, pprocmirrpath, reprojpath, oridistpath
        real, allocatable        :: fsc(:), res(:)
        real                     :: smpd, res05, res0143, vmin, vmax
        integer                  :: istate, box, pop, n
        logical                  :: l_fsc
        self%state_res = 0.
        do istate = 1,self%params%nstates
            if( .not. self%spproj%isthere_in_osout('vol', istate) ) cycle
            call self%spproj%get_vol('vol', istate, volpath, smpd, box)
            if( volpath%strlen() == 0 ) cycle
            pop = self%spproj%os_ptcl3D%get_pop(istate, 'state')
            ! products of a postprocessing, when on disk (none runs in this stage)
            pprocpath = add2fbody(volpath, MRC_EXT, PPROC_SUFFIX)
            if( .not. file_exists(pprocpath) ) pprocpath = ''
            lppath = add2fbody(volpath, MRC_EXT, LP_SUFFIX)
            if( .not. file_exists(lppath) ) lppath = ''
            pprocmirrpath = ''
            if( pprocpath%strlen() > 0 )then
                pprocmirrpath = add2fbody(pprocpath, MRC_EXT, MIRR_SUFFIX)
                if( .not. file_exists(pprocmirrpath) ) pprocmirrpath = ''
            endif
            ! JPEGs the 3D writes beside the volume
            reprojpath = get_fpath(volpath)//'orthogonal_reprojs_state'//int2str_pad(istate, 2)//JPG_EXT
            if( .not. file_exists(reprojpath) ) reprojpath = ''
            oridistpath = get_fpath(volpath)//refine3D_oris_heatmap_fname(istate)
            if( .not. file_exists(oridistpath) ) oridistpath = ''
            if( reprojpath%strlen() > 0 ) call send_reproj_tiles(self%pipe, self%meta_reproj, reprojpath, volpath,&
                &istate, self%params%nstates, pop)
            l_fsc = self%read_state_fsc(istate, smpd, fsc, res, res05, res0143)
            ! every field back to its default: new and kill reset only the flags, and a state
            ! without an FSC curve or a product must not carry the previous state's
            meta = fresh_meta
            call meta%new(GUI_METADATA_VOL3D_TYPE)
            if( l_fsc )then
                self%state_res(istate) = res0143
                call meta%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, istate, box, smpd, istate,&
                    &self%params%nstates, res0143=res0143, res05=res05, pop=pop, oridistpath=oridistpath)
                n = min(size(fsc), size(res), 1000) ! gui_metadata_vol3D holds 1000 points
                call meta%set_fsc(1. / res(:n), fsc(:n))
            else
                call meta%set(reprojpath, volpath, pprocpath, lppath, pprocmirrpath, istate, box, smpd, istate,&
                    &self%params%nstates, pop=pop, oridistpath=oridistpath)
            endif
            ! header minimum and maximum, so the GUI need not open the volumes
            call get_mrc_minmax(volpath, vmin, vmax)
            call meta%set_minmax('volpath', vmin, vmax)
            if( pprocpath%strlen() > 0 )then
                call get_mrc_minmax(pprocpath, vmin, vmax)
                call meta%set_minmax('pprocpath', vmin, vmax)
            endif
            if( lppath%strlen() > 0 )then
                call get_mrc_minmax(lppath, vmin, vmax)
                call meta%set_minmax('lppath', vmin, vmax)
            endif
            if( pprocmirrpath%strlen() > 0 )then
                call get_mrc_minmax(pprocmirrpath, vmin, vmax)
                call meta%set_minmax('pprocmirrpath', vmin, vmax)
            endif
            call meta%set_oridist_from_oris(self%spproj%os_ptcl3D, istate)
            call self%pipe%send_meta(meta)
            call meta%kill
        enddo
    end subroutine send_volumes

    ! The FSC of state @p istate in the pool's project, its resolutions at 1/(Fourier index) for
    ! pixel size @p smpd, and the FSC=0.5 and FSC=0.143 resolutions; .false. without an FSC file.
    logical function read_state_fsc( self, istate, smpd, fsc, res, res05, res0143 ) result( l_fsc )
        class(stream_stage_solve3D), intent(in)    :: self
        integer,                        intent(in)    :: istate
        real,                           intent(in)    :: smpd
        real, allocatable,              intent(inout) :: fsc(:), res(:)
        real,                           intent(out)   :: res05, res0143
        type(string) :: fsc_fname
        integer      :: fsc_box, n
        l_fsc   = .false.
        res05   = 0.
        res0143 = 0.
        if( .not. self%spproj%isthere_in_osout('fsc', istate) ) return
        call self%spproj%get_fsc(istate, fsc_fname, fsc_box)
        if( fsc_fname%strlen() == 0 ) return
        if( .not. file_exists(fsc_fname) ) return
        fsc = file2rarr(fsc_fname)
        res = get_resarr(fsc_box, smpd)
        n   = min(size(fsc), size(res))
        if( n == 0 ) return
        fsc = fsc(:n)
        res = res(:n)
        call get_resolution(fsc, res, res05, res0143)
        l_fsc = .true.
    end function read_state_fsc

    !---------------- rules ----------------

    !> The job to start in @p phase with @p nptcls particles in the pool: solve3D on the first
    !! particles, an addon run once the pool has grown past the @p nptcls_last_run of the latest
    !! run; none while a job runs.
    pure integer function next_job( phase, nptcls, nptcls_last_run )
        integer, intent(in) :: phase, nptcls, nptcls_last_run
        next_job = JOB_NONE
        select case(phase)
            case(PHASE_IMPORTING)
                if( nptcls > 0 ) next_job = JOB_SOLVE3D
            case(PHASE_IDLE)
                if( nptcls > nptcls_last_run ) next_job = JOB_ADDON
        end select
    end function next_job

    !> The mask diameter (A) for class averages of @p box pixels at @p smpd: @p mskdiam, or the
    !! box default when it is not positive or too large for the box, as parameters falls back.
    pure real function fit_mskdiam( mskdiam, box, smpd )
        real,    intent(in) :: mskdiam, smpd
        integer, intent(in) :: box
        real :: msk, msk_default
        fit_mskdiam = mskdiam
        msk_default = round2even((real(box) - COSMSKHALFWIDTH) / 2.)
        if( msk_default <= 0. ) return
        msk = 0.
        if( smpd > 0. ) msk = round2even((mskdiam / smpd) / 2.)
        if( msk < 0.1 .or. msk > msk_default ) fit_mskdiam = (real(box) - COSMSKHALFWIDTH) * smpd
    end function fit_mskdiam

end module simple_stream_stage_solve3D
