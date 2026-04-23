!@descr: task 2 in the stream pipeline: assign optics groups to streamed micrographs
!==============================================================================
! MODULE: simple_stream_p02_assign_optics_new
!
! PURPOSE:
!   Watches for completed preprocessing projects, imports micrograph metadata,
!   updates optics-group assignments/maps, exports STAR products, and streams
!   optics-assignment metadata to the GUI.
!
! ENTRY POINT:
!   stream_p02_assign_optics%execute(cline)
!
! NOTES:
!   - Runs continuously until TERM_STREAM or SIGTERM.
!   - Maintains rolling optics maps and removes stale map files.
!==============================================================================
module simple_stream_p02_assign_optics_new
use unix, only: c_time, SIGTERM
use simple_stream_api
use simple_stream_mq_defs    
use simple_gui_metadata_api
implicit none

public :: stream_p02_assign_optics
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p02_assign_optics
  contains
    procedure :: execute => exec_stream_p02_assign_optics
end type stream_p02_assign_optics

contains

    ! Main loop for stream task 2. Imports micrograph metadata from newly
    ! completed projects, updates optics outputs, and publishes GUI metadata.
    subroutine exec_stream_p02_assign_optics( self, cline )
        class(stream_p02_assign_optics), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        integer, parameter                         :: MAX_MOVIES_IMPORT = 50
        type(parameters)                           :: params
        type(stream_watcher)                       :: project_buff
        type(sp_project)                           :: spproj, spproj_part
        type(gui_metadata_stream_optics_assignment) :: meta_optics_assignment
        type(gui_metadata_optics_group)             :: meta_optics_group
        type(starproject_stream)                   :: starproj_stream
        type(string),                  allocatable :: projects(:)
        type(string)                               :: projfile
        integer,                       allocatable :: nmics_per_proj(:)
        character(len=:),        allocatable :: meta_buffer
        real,                    allocatable :: xshifts(:), yshifts(:)
        integer                                    :: nprojects, iproj, iori, nmics_here, nimported, nimported_cycle
        integer                                    :: i, j, map_count, imap, last_micrograph_imported, i_point, n_points, max_points
        logical                                    :: l_terminate=.false.
        call signal(SIGTERM, sigterm_handler)
        call cline%printline()
        call flush(logfhandle)
        map_count = 0
        last_micrograph_imported = 0
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('dir_target') ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('outdir')     ) call cline%set('outdir', '')
        ! sanity check for restart
        if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
        endif
        ! generate own project file if projfile doesn't exist
        projfile = cline%get_carg('projfile')
        if( .not.file_exists(projfile) )then
            call cline%set('projfile', projfile)
            call spproj%update_projinfo(cline)
            call spproj%update_compenv(cline)
            call spproj%write
        endif
        ! master parameters
        call params%new(cline)
        map_count = get_latest_optics_map_id(params%outdir)
        ! initialise metadata
        call meta_optics_assignment%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_TYPE)
        call meta_optics_group%new(GUI_METADATA_STREAM_OPTICS_ASSIGNMENT_OPTICS_GROUP_TYPE)
        ! master project file
        call spproj%read( params%projfile )
        if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(0, .false.)
        ! wait if dir_target doesn't exist yet
        call wait_for_stream_folder(params%dir_target)
        call wait_for_stream_folder(params%dir_target//'/spprojs')
        call wait_for_stream_folder(params%dir_target//'/spprojs_completed')
        ! movie watcher init
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! Infinite loop
        nimported = 0
        do
            if( file_exists(TERM_STREAM) .or.  l_terminate ) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            ! detection of new projects
            call project_buff%watch( nprojects, projects, max_nmovies=MAX_MOVIES_IMPORT )
            ! append projects to processing stack
            if( nprojects > 0 )then
                nimported = spproj%os_mic%get_noris()
                nimported_cycle = 0
                if( allocated(nmics_per_proj) ) deallocate(nmics_per_proj)
                allocate(nmics_per_proj(nprojects), source=0)
                ! Pre-scan projects to compute the exact import size.
                do iproj = 1, nprojects
                    call project_buff%add2history(projects(iproj))
                    call spproj_part%read_segment('mic', projects(iproj))
                    nmics_here = min(STREAM_NMOVS_SET, spproj_part%os_mic%get_noris())
                    nmics_per_proj(iproj) = nmics_here
                    nimported_cycle = nimported_cycle + nmics_here
                    call spproj_part%kill()
                enddo
                if( nimported_cycle > 0 )then
                    if( nimported > 0 )then
                        call spproj%os_mic%reallocate(nimported + nimported_cycle)
                    else
                        call spproj%os_mic%new(nimported_cycle, .false.)
                    endif
                endif
                do iproj = 1, nprojects
                    if( nmics_per_proj(iproj) == 0 ) cycle
                    call spproj_part%read_segment('mic', projects(iproj))
                    nmics_here = nmics_per_proj(iproj)
                    do iori = 1, nmics_here
                        nimported = nimported + 1
                        call spproj%os_mic%transfer_ori(nimported, spproj_part%os_mic, iori)
                        last_micrograph_imported = int(c_time(0_c_long))
                    end do
                    call spproj_part%kill()
                enddo
                if( allocated(nmics_per_proj) ) deallocate(nmics_per_proj)
                write(logfhandle,'(A,I6,A,A)')'>>> ', nimported_cycle, ' NEW MICROGRAPHS IMPORTED; ', cast_time_char(simple_gettime())
                call starproj_stream%stream_export_optics(params, spproj, params%outdir)
                call starproj_stream%stream_export_micrographs(params, spproj, params%outdir, optics_set=.true.)
                max_points = meta_optics_group%get_max_points()
                if( allocated(xshifts) ) deallocate(xshifts)
                if( allocated(yshifts) ) deallocate(yshifts)
                allocate(xshifts(max_points), yshifts(max_points))
                do i = 1, spproj%os_optics%get_noris()
                    i_point = 1
                    xshifts = 0.0
                    yshifts = 0.0
                    do j = spproj%os_mic%get_noris(), 1, -1 
                        if (spproj%os_mic%get(j, 'ogid') .eq. real(i)) then
                            if( i_point > max_points ) exit
                            xshifts(i_point) = spproj%os_mic%get(j, 'shiftx')
                            yshifts(i_point) = spproj%os_mic%get(j, 'shifty')
                            i_point = i_point + 1
                        endif
                    enddo
                    n_points = max(0, i_point - 1)
                    call meta_optics_group%set(i, spproj%os_optics%get_noris(), xshifts, yshifts, n_points)
                    if( meta_optics_group%assigned() .and. mq_stream_master_in%is_active() ) then
                        call meta_optics_group%serialise(meta_buffer)
                        call mq_stream_master_in%send(meta_buffer)
                    endif
                enddo
                if( allocated(xshifts) ) deallocate(xshifts)
                if( allocated(yshifts) ) deallocate(yshifts)
                map_count = map_count + 1
                call spproj%write_optics_map(OPTICS_MAP_PREFIX // int2str(map_count))
                ! remove old maps > 5 iterations ago
                if(map_count > 5) then
                    do imap=1, map_count - 5
                        if(file_exists(OPTICS_MAP_PREFIX // int2str(imap) // TXT_EXT))      call del_file(OPTICS_MAP_PREFIX // int2str(imap) // TXT_EXT)
                        if(file_exists(OPTICS_MAP_PREFIX // int2str(imap) // METADATA_EXT)) call del_file(OPTICS_MAP_PREFIX // int2str(imap) // METADATA_EXT)
                    enddo
                endif 
            else
                call sleep(WAITTIME) ! may want to increase as 3s default
            endif
            call update_user_params(params, cline)
            call meta_optics_assignment%set(stage=string('finding and processing new micrographs'),&
                                            micrographs_assigned=spproj%os_mic%get_noris(),        & 
                                            optics_groups_assigned=spproj%os_optics%get_noris(),   &
                                            last_micrograph_imported=last_micrograph_imported)
            if(params%updated .eq. 'yes') then
                call starproj_stream%stream_export_optics(params, spproj, params%outdir)
                params%updated = 'no'
            end if
            if( meta_optics_assignment%assigned() .and. mq_stream_master_in%is_active() ) then
                call meta_optics_assignment%serialise(meta_buffer)
                call mq_stream_master_in%send(meta_buffer)
            endif
            call flush(logfhandle)
        end do
        ! termination
        if(allocated(projects)) deallocate(projects)
        if(allocated(nmics_per_proj)) deallocate(nmics_per_proj)
        ! cleanup
        call spproj%write
        call spproj%kill
        ! end gracefully
        call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS NORMAL STOP ****')

        contains
        
            ! Poll for a required directory. Returns early on termination.
            subroutine wait_for_stream_folder( folder )
                class(string),                   intent(in)    :: folder
                integer :: iwait
                if(.not. dir_exists(folder)) then
                    write(logfhandle, *) ">>> WAITING FOR ", folder%to_char(), " TO BE GENERATED"
                    ! wait up to 24 hours
                    do iwait = 1, 8640
                        if( file_exists(TERM_STREAM) .or. l_terminate ) return
                        if(dir_exists(folder)) then
                            write(logfhandle, *) ">>> ", folder%to_char(), " FOUND"
                            exit
                        endif
                        call sleep(10)
                    end do
                endif
            end subroutine wait_for_stream_folder

            subroutine sigterm_handler()
                l_terminate = .true.
            end subroutine sigterm_handler

    end subroutine exec_stream_p02_assign_optics

end module simple_stream_p02_assign_optics_new
