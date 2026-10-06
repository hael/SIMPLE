!@descr: state and steps of stream task 2 (optics assignment): import preprocessed micrographs, group them, publish the groups
!==============================================================================
! MODULE: simple_stream_stage_optics
!
! PURPOSE:
!   The body of stream p02 as a type. Each pass imports the micrographs of
!   newly completed preprocessing projects, groups all micrographs again
!   (single linkage of their beam shifts, simple_optics_groups) keeping each
!   group's id from the previous grouping, writes the STAR files, publishes
!   the next optics map for the later stages, and reports to the GUI. The commander
!   (simple_commanders_stream_p02_assign_optics) only normalises the command line and
!   loops over iterate() until finished().
!
!   What happens is delegated:
!     - micrograph import       -> simple_mic_import
!     - optics-group assignment -> simple_optics_groups
!     - STAR files              -> starproject_stream (write only)
!     - optics-map files        -> simple_optics_maps
!     - GUI plots / pipe        -> simple_stream_meta_plots, simple_stream_pipe
!
! LIFECYCLE:
!   new(cline) -> { iterate() } until finished() -> finalize() -> kill()
!
! RESTART:
!   The micrograph segment starts empty and every completed upstream project
!   is imported again (the watcher history starts empty); the optics-map ids
!   continue from the newest map in the stage directory. The micrographs that
!   map lists come back with its group ids, and no map is published until all
!   of them are imported again (or a pass finds nothing more to import), so
!   the later stages keep applying a complete map meanwhile.
!==============================================================================
module simple_stream_stage_optics
use simple_defs,                                  only: logfhandle
use simple_defs_fname,                            only: TERM_STREAM
use simple_defs_stream,                           only: DIR_STREAM, DIR_STREAM_COMPLETED, SHORTWAIT, WAITTIME
use simple_string,                                only: string
use simple_fileio,                                only: del_file, file_exists, simple_getcwd
use simple_syslib,                                only: dir_exists
use simple_timer,                                 only: simple_gettime, cast_time_char
use simple_cmdline,                               only: cmdline
use simple_parameters,                            only: parameters
use simple_sp_project,                            only: sp_project
use simple_starproject_stream,                    only: starproject_stream
use simple_stream_watcher,                        only: stream_watcher
use simple_stream_state,                          only: ipc_pipe_assign_optics_in
use simple_gui_metadata_utils,                    only: max_metadata_size
use simple_gui_metadata_types,                    only: GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE,&
                                                       &GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE
use simple_gui_metadata_optics_group,             only: gui_metadata_optics_group, MAX_OPTICS_SHIFTS
use simple_gui_metadata_stream_optics_assignment, only: gui_metadata_stream_optics_assignment
use simple_stream_pipe,                           only: stream_pipe
use simple_mic_import,                            only: append_mics_from_projects
use simple_optics_groups,                         only: assign_optics_groups
use simple_optics_maps,                           only: publish_optics_map, latest_optics_map_id, latest_optics_map_table
use simple_stream_meta_plots,                     only: recent_shifts_by_optics_group
implicit none

public :: stream_stage_optics
private

integer, parameter :: MAX_PROJECTS_IMPORT = 50 ! completed upstream projects taken per pass
integer, parameter :: NMAPS_KEPT          = 5  ! optics maps left on disk for the readers

! Components and steps are public so simple_stream_stage_optics_tester can assemble a stage
! and run one step at a time; production code uses new/iterate/finished/finalize/kill.
type :: stream_stage_optics
    type(parameters), allocatable               :: params
    type(sp_project), allocatable               :: spproj           ! every imported micrograph and the optics groups (new..kill)
    type(stream_watcher)                        :: project_buff     ! completed preprocessing projects
    type(starproject_stream)                    :: starproj_stream
    type(stream_pipe)                           :: pipe             ! to the master
    type(gui_metadata_stream_optics_assignment) :: meta_status
    type(gui_metadata_optics_group)             :: meta_group
    type(string)                                :: map_dir          ! absolute; where the optics maps are published
    integer :: map_id            = 0
    integer :: last_ogid         = 0        ! the highest optics group id given; ids are kept across passes
    ! on a restart: the newest map's group id of each import index (0: not listed), the micrographs
    ! it lists not imported again yet, and whether no map is published until they are
    integer, allocatable :: restored_ogid(:)
    integer :: n_restore_left    = 0
    logical :: l_restoring       = .false.
    logical :: l_ungrouped       = .false.  ! micrographs imported while restoring, not grouped yet
    logical :: l_attached        = .false.  ! the upstream completed-projects folder exists and is watched
    logical :: l_waiting_logged  = .false.
    logical :: l_nmics_reached   = .false.
    logical :: l_exists          = .false.
    ! waits (s); tests set them to 0, and settle_s to -1 to take files written in the same second
    integer :: settle_s          = SHORTWAIT ! an upstream project is taken once untouched longer than this;
                                             ! preprocessing moves finished projects in with a rename
    integer :: wait_s            = WAITTIME  ! pause in a pass that imports nothing
contains
    procedure :: new
    procedure :: iterate
    procedure :: finished
    procedure :: finalize
    procedure :: kill
    procedure :: attach_upstream
    procedure :: import_new_projects
    procedure :: assign_and_publish
    procedure :: restore_groups
    procedure :: end_restore
    procedure :: send_group_shifts
    procedure :: send_status
end type stream_stage_optics

contains

    !---------------- lifecycle ----------------

    !> Builds the parameters from the normalised @p cline and prepares an empty micrograph segment.
    subroutine new( self, cline )
        class(stream_stage_optics), intent(inout) :: self
        class(cmdline),             intent(inout) :: cline
        type(string) :: outdir, projfile
        logical      :: l_restart
        call self%kill()
        allocate(self%spproj)
        ! a restart is recognised by its output directory, before params%new creates one
        l_restart = .false.
        outdir    = cline%get_carg('outdir')
        if( .not. (outdir == '') ) l_restart = dir_exists(outdir)
        ! own project file when none exists yet: params%new copies it into the stage's folder (the
        ! program requires a project), and the stage works on that copy (params%projfile); nothing
        ! reads the one left in the launch folder
        projfile = cline%get_carg('projfile')
        if( .not. file_exists(projfile) )then
            call self%spproj%update_projinfo(cline)
            call self%spproj%update_compenv(cline)
            call self%spproj%write
        endif
        if( .not. allocated(self%params) ) allocate(self%params)
        call self%params%new(cline)
        ! params%new has moved into the stage directory, where the maps are written and read
        ! and where finished() looks for the termination file: one left by the previous run
        ! would stop this one at once
        if( l_restart )then
            write(logfhandle,'(A)') '>>> RESTARTING EXISTING JOB'
            call del_file(TERM_STREAM)
        endif
        call simple_getcwd(self%map_dir)
        self%map_id = latest_optics_map_id(self%map_dir)
        if( self%map_id > 0 ) call self%restore_groups()
        call self%spproj%read(self%params%projfile)
        if( self%spproj%os_mic%get_noris() /= 0 ) call self%spproj%os_mic%new(0, is_ptcl=.false.)
        call self%meta_status%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
        call self%meta_group%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE)
        ! the master sends this stage nothing, so only the write end is used
        call self%pipe%new(-1, ipc_pipe_assign_optics_in(2), max_metadata_size(), 'assign_optics')
        self%l_exists = .true.
    end subroutine new

    !> One pass: attach to the upstream folder until it exists, then import, group, publish, report.
    subroutine iterate( self )
        class(stream_stage_optics), intent(inout) :: self
        integer :: nimported
        if( .not. self%l_attached )then
            call self%attach_upstream()
            if( .not. self%l_attached )then
                call sleep(self%wait_s)
                return
            endif
        endif
        call self%import_new_projects(nimported)
        if( self%l_restoring )then
            if( self%n_restore_left <= 0 )then
                call self%end_restore()
            else if( nimported == 0 .and. self%spproj%os_mic%get_noris() > 0 )then
                ! every completed upstream project is taken: those micrographs are gone
                write(logfhandle,'(A,I8,A)') '>>> WARNING: ', self%n_restore_left,&
                    &' MICROGRAPHS OF THE NEWEST OPTICS MAP DID NOT COME BACK; PUBLISHING WITHOUT THEM'
                call self%end_restore()
            endif
        endif
        if( self%l_restoring )then
            if( nimported > 0 ) self%l_ungrouped = .true.
            if( nimported == 0 ) call sleep(self%wait_s)
        else if( nimported > 0 .or. self%l_ungrouped )then
            call self%assign_and_publish()
            self%l_ungrouped = .false.
        else
            call sleep(self%wait_s)
        endif
        call self%send_status()
        if( self%params%nmics > 0 )then
            if( self%spproj%os_mic%get_noris() >= self%params%nmics .and. .not. self%l_nmics_reached )then
                write(logfhandle,'(A,I8)') '>>> TERMINATING PROCESS: requested number of micrographs reached: ', self%params%nmics
                self%l_nmics_reached = .true.
            endif
        endif
        call flush(logfhandle)
    end subroutine iterate

    !> .true. once the stream is told to stop or the requested number of micrographs is imported.
    logical function finished( self )
        class(stream_stage_optics), intent(in) :: self
        finished = self%l_nmics_reached .or. file_exists(TERM_STREAM)
    end function finished

    !> Writes the project with every micrograph and its optics group.
    subroutine finalize( self )
        class(stream_stage_optics), intent(inout) :: self
        call self%spproj%write(self%params%projfile)
    end subroutine finalize

    subroutine kill( self )
        class(stream_stage_optics), intent(inout) :: self
        if( allocated(self%spproj) )then
            call self%spproj%kill
            deallocate(self%spproj)
        endif
        if( .not. self%l_exists )then
            if( allocated(self%params) ) deallocate(self%params)
            return
        endif
        call self%project_buff%kill
        call self%pipe%kill
        call self%meta_status%kill
        call self%meta_group%kill
        call self%map_dir%kill
        if( allocated(self%params) ) deallocate(self%params)
        if( allocated(self%restored_ogid) ) deallocate(self%restored_ogid)
        self%map_id           = 0
        self%last_ogid        = 0
        self%n_restore_left   = 0
        self%l_restoring      = .false.
        self%l_ungrouped      = .false.
        self%l_attached       = .false.
        self%l_waiting_logged = .false.
        self%l_nmics_reached  = .false.
        self%l_exists         = .false.
    end subroutine kill

    !---------------- steps ----------------

    ! Starts watching the upstream completed-projects folder once preprocessing has created it.
    subroutine attach_upstream( self )
        class(stream_stage_optics), intent(inout) :: self
        type(string) :: completed
        logical      :: l_ready
        completed = self%params%dir_target//'/'//DIR_STREAM_COMPLETED
        l_ready   = dir_exists(self%params%dir_target)
        if( l_ready ) l_ready = dir_exists(self%params%dir_target//'/'//DIR_STREAM)
        if( l_ready ) l_ready = dir_exists(completed)
        if( .not. l_ready )then
            if( .not. self%l_waiting_logged )then
                write(logfhandle,'(A)') '>>> WAITING FOR '//completed%to_char()//' TO BE GENERATED'
                self%l_waiting_logged = .true.
            endif
            return
        endif
        write(logfhandle,'(A)') '>>> '//completed%to_char()//' FOUND'
        self%project_buff = stream_watcher(self%settle_s, completed, spproj=.true., nretries=10)
        self%l_attached   = .true.
    end subroutine attach_upstream

    ! Appends every micrograph, accepted or not, of the newly completed upstream projects. A new
    ! micrograph has no optics group, unless a restart restores the one the newest map gave it.
    subroutine import_new_projects( self, nimported )
        class(stream_stage_optics), intent(inout) :: self
        integer,                    intent(out)   :: nimported
        type(string), allocatable :: projects(:)
        integer :: nprojects, nbefore, imic, importind, ogid
        nimported = 0
        call self%project_buff%watch(nprojects, projects, max_nmovies=MAX_PROJECTS_IMPORT)
        if( nprojects == 0 ) return
        call self%project_buff%add2history(projects)
        nbefore = self%spproj%os_mic%get_noris()
        call append_mics_from_projects(self%spproj%os_mic, projects, .false., nimported)
        do imic = nbefore + 1,self%spproj%os_mic%get_noris()
            ogid = 0
            if( allocated(self%restored_ogid) )then
                importind = self%spproj%os_mic%get_int(imic, 'importind')
                if( importind >= 1 .and. importind <= size(self%restored_ogid) )then
                    ogid = self%restored_ogid(importind)
                    if( ogid > 0 ) self%n_restore_left = self%n_restore_left - 1
                endif
            endif
            call self%spproj%os_mic%set(imic, 'ogid', real(ogid))
        enddo
        write(logfhandle,'(A,I6,A,A)') '>>> ', nimported, ' NEW MICROGRAPHS IMPORTED; ', cast_time_char(simple_gettime())
    end subroutine import_new_projects

    ! A restart: the newest map's group id of each micrograph it lists, which the micrograph takes
    ! back when it is imported again, and the highest id given; no map until they are all back.
    subroutine restore_groups( self )
        class(stream_stage_optics), intent(inout) :: self
        integer, allocatable :: importinds(:), ogids(:)
        integer :: id, i
        id = latest_optics_map_table(self%map_dir, importinds, ogids)
        if( size(importinds) == 0 ) return
        allocate(self%restored_ogid(max(1, maxval(importinds))), source=0)
        do i = 1,size(importinds)
            if( importinds(i) >= 1 ) self%restored_ogid(importinds(i)) = max(0, ogids(i))
        enddo
        self%n_restore_left = count(self%restored_ogid > 0)
        self%last_ogid      = max(0, maxval(ogids))
        self%l_restoring    = self%n_restore_left > 0
        write(logfhandle,'(A,I6,A,I8,A)') '>>> RESTORING THE OPTICS GROUPS OF MAP ', id, ' FOR ', self%n_restore_left,&
            &' MICROGRAPHS; NO MAP IS PUBLISHED UNTIL THEY ARE IMPORTED AGAIN'
    end subroutine restore_groups

    ! The restored micrographs are back (or will not come): maps are published again.
    subroutine end_restore( self )
        class(stream_stage_optics), intent(inout) :: self
        if( self%n_restore_left <= 0 ) write(logfhandle,'(A)') '>>> THE MICROGRAPHS OF THE NEWEST OPTICS MAP ARE BACK'
        self%l_restoring    = .false.
        self%n_restore_left = 0
        if( allocated(self%restored_ogid) ) deallocate(self%restored_ogid)
    end subroutine end_restore

    ! Regroups every micrograph, keeping the groups' ids, writes the STAR files and the project,
    ! reports the groups and publishes the next optics map.
    subroutine assign_and_publish( self )
        class(stream_stage_optics), intent(inout) :: self
        call assign_optics_groups(self%spproj, self%params%tilt_thres, self%params%beamtilt == 'yes', 0,&
            &last_ogid=self%last_ogid)
        call self%starproj_stream%stream_write_optics(self%params, self%spproj, self%params%cwd)
        call self%starproj_stream%stream_export_micrographs(self%params, self%spproj, self%params%cwd, optics_set=.true.)
        ! the STAR exporter used to rewrite the project here as a side effect; kept, now explicit,
        ! to the stage's project file
        call self%spproj%write(self%params%projfile)
        call self%send_group_shifts()
        self%map_id = self%map_id + 1
        call publish_optics_map(self%spproj, self%map_dir, self%map_id, NMAPS_KEPT)
    end subroutine assign_and_publish

    ! One GUI message per optics group with the shifts of its newest micrographs.
    subroutine send_group_shifts( self )
        class(stream_stage_optics), intent(inout) :: self
        real,    allocatable :: xshifts(:,:), yshifts(:,:), xs(:), ys(:)
        integer, allocatable :: npoints(:)
        integer :: ngroups, igroup
        call recent_shifts_by_optics_group(self%spproj%os_mic, self%spproj%os_optics, MAX_OPTICS_SHIFTS,&
            &xshifts, yshifts, npoints)
        ngroups = size(npoints)
        do igroup = 1,ngroups
            xs = xshifts(:,igroup)
            ys = yshifts(:,igroup)
            call self%meta_group%set(igroup, ngroups, xs, ys, npoints(igroup))
            call self%pipe%send_meta(self%meta_group)
        enddo
    end subroutine send_group_shifts

    ! Imported: every micrograph taken in, accepted or not; assigned: the accepted ones.
    subroutine send_status( self )
        class(stream_stage_optics), intent(inout) :: self
        integer :: naccepted
        naccepted = self%spproj%os_mic%count_state_gt_zero()
        call self%meta_status%set(stage=string('finding and processing new micrographs'),&
            micrographs_assigned   = naccepted,                                         &
            optics_groups_assigned = self%spproj%os_optics%get_noris(),                &
            micrographs_imported   = self%spproj%os_mic%get_noris())
        call self%pipe%send_meta(self%meta_status)
    end subroutine send_status

end module simple_stream_stage_optics
