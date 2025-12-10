module simple_stream_p02_assign_optics
include 'simple_lib.f08'
use simple_cmdline,            only: cmdline
use simple_commander_base,     only: commander_base
use simple_parameters,         only: parameters
use simple_sp_project,         only: sp_project
use simple_starproject_stream, only: starproject_stream
use simple_stream_watcher,     only: stream_watcher
use simple_gui_utils
use simple_progress
use simple_qsys_funs
use simple_stream_communicator
use simple_stream_utils
use json_kinds
implicit none

public :: stream_p02_assign_optics
private
#include "simple_local_flags.inc"

type, extends(commander_base) :: stream_p02_assign_optics
  contains
    procedure :: execute => exec_stream_p02_assign_optics
end type stream_p02_assign_optics

contains

    subroutine exec_stream_p02_assign_optics( self, cline )
        class(stream_p02_assign_optics), intent(inout) :: self
        class(cmdline),                  intent(inout) :: cline
        type(parameters)               :: params
        type(stream_http_communicator) :: http_communicator
        type(stream_watcher)           :: project_buff
        type(sp_project)               :: spproj, spproj_part
        type(starproject_stream)       :: starproj_stream
        type(json_value),  pointer     :: optics_assignments, optics_group     
        type(json_value),  pointer     :: optics_group_coordinates,  optics_group_coordinate
        type(json_core)                :: json   
        type(string),      allocatable :: projects(:)
        type(string)                   :: str_dir
        integer                        :: nprojects, iproj, iori, new_oris, nimported, i, j, map_count, imap
        logical                        :: found
        map_count = 0
        call cline%set('mkdir', 'yes')
        if( .not. cline%defined('dir_target') ) THROW_HARD('DIR_TARGET must be defined!')
        if( .not. cline%defined('outdir')     ) call cline%set('outdir', '')
        ! sanity check for restart
        if(cline%defined('outdir') .and. dir_exists(cline%get_carg('outdir'))) then
            write(logfhandle,'(A)') ">>> RESTARTING EXISTING JOB"
            call del_file(TERM_STREAM)
        endif
        ! below may be redundant
        if( cline%defined('dir_exec') )then
            if( .not.file_exists(cline%get_carg('dir_exec')) )then
                str_dir = cline%get_carg('dir_exec')
                THROW_HARD('Previous directory does not exist: '//str_dir%to_char())
                call str_dir%kill
            endif
            call del_file(TERM_STREAM)
        endif
        ! generate own project file if projfile isnt set
        if( .not.cline%defined('projfile') )then
            call cline%set('projname', 'assign_optics')
            call cline%set('projfile', 'assign_optics.simple')
            call spproj%update_projinfo(cline)
            call spproj%update_compenv(cline)
            call spproj%write
        endif
        ! master parameters
        call params%new(cline)
        ! http communicator init
        call http_communicator%create(params%niceprocid, params%niceserver%to_char(), "optics_assignment")
        call communicator_init()
        call http_communicator%send_jobstats()
        ! master project file
        call spproj%read( params%projfile )
        if( spproj%os_mic%get_noris() /= 0 ) call spproj%os_mic%new(0, .false.)
        ! wait if dir_target doesn't exist yet
        call wait_for_folder(http_communicator, params%dir_target, '**** SIMPLE_STREAM_ASSIGN_OPTICS USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs', '**** SIMPLE_STREAM_ASSIGN_OPTICS USER STOP ****')
        call wait_for_folder(http_communicator, params%dir_target//'/spprojs_completed', '**** SIMPLE_STREAM_ASSIGN_OPTICS USER STOP ****')
        ! movie watcher init
        project_buff = stream_watcher(LONGTIME, params%dir_target//'/'//DIR_STREAM_COMPLETED, spproj=.true., nretries=10)
        ! initialise progress monitor
        call progressfile_init()
        ! Infinite loop
        nimported = 0
        do
            if( file_exists(TERM_STREAM) .or. http_communicator%exit_status()) then
                ! termination
                write(logfhandle,'(A)')'>>> TERMINATING PROCESS'
                exit
            endif
            if( http_communicator%stop_status() )then
                ! termination
                write(logfhandle,'(A)')'>>> USER COMMANDED STOP'
                call spproj%kill
                call qsys_cleanup
                call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS USER STOP ****')
                call EXIT(0)
            endif
            ! http stats
            call http_communicator%update_json("stage", "finding and processing new micrographs", found) 
            ! detection of new projects
            call project_buff%watch( nprojects, projects, max_nmovies=50 )
            ! append projects to processing stack
            if( nprojects > 0 )then
                nimported = spproj%os_mic%get_noris()
                if(nimported > 0) then
                    new_oris  =  nimported + nprojects * STREAM_NMOVS_SET
                    call spproj%os_mic%reallocate(new_oris)
                else
                    new_oris = nprojects * STREAM_NMOVS_SET
                    call spproj%os_mic%new(new_oris, .false.)
                end if
                do iproj = 1, nprojects
                    call project_buff%add2history(projects(iproj))
                    call spproj_part%read(projects(iproj))
                    do iori = 1, STREAM_NMOVS_SET
                        nimported = nimported + 1
                        call spproj%os_mic%transfer_ori(nimported, spproj_part%os_mic, iori)
                    end do
                    call spproj_part%kill()
                enddo
                write(logfhandle,'(A,I4,A,A)')'>>> ' , nprojects * STREAM_NMOVS_SET, ' NEW MICROGRAPHS IMPORTED; ',cast_time_char(simple_gettime())
                call starproj_stream%stream_export_optics(spproj, params%outdir)
                call starproj_stream%stream_export_micrographs(spproj, params%outdir, optics_set=.true.)
                ! http stats
                call http_communicator%update_json("micrographs_assigned",     spproj%os_mic%get_noris(),    found)
                call http_communicator%update_json("optics_groups_assigned",   spproj%os_optics%get_noris(), found)
                call http_communicator%update_json("last_micrograph_imported", stream_datestr(),             found)
                call json%remove(optics_assignments, destroy=.true.)
                call json%create_array(optics_assignments, "optics_assignments")
                call http_communicator%add_to_json(optics_assignments)
                do i = 1, spproj%os_optics%get_noris()
                    call json%create_array(optics_group_coordinates, "coordinates")
                    do j = 1, spproj%os_mic%get_noris()
                        if (spproj%os_mic%get(j, 'ogid') .eq. real(i)) then
                            call json%create_object(optics_group_coordinate, "")
                            call json%add(optics_group_coordinate, "x", dble(spproj%os_mic%get(i, 'shiftx')))
                            call json%add(optics_group_coordinate, "y", dble(spproj%os_mic%get(i, 'shifty')))
                            call json%add(optics_group_coordinates, optics_group_coordinate)
                        endif
                    enddo
                    call json%create_object(optics_group, "")
                    call json%add(optics_group, optics_group_coordinates)
                    call json%add(optics_group, "id", i)
                    call json%add(optics_assignments, optics_group)
                enddo
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
            call update_user_params(cline)
            if(params%updated .eq. 'yes') then
                call starproj_stream%stream_export_optics(spproj, params%outdir)
                params%updated = 'no'
            end if
            ! http stats send
            call http_communicator%send_jobstats()
        end do
        ! termination
        call http_communicator%update_json("stage", "terminating", found) 
        call http_communicator%send_jobstats()
        if(allocated(projects)) deallocate(projects)
        ! cleanup
        call spproj%write
        call spproj%kill
        ! end gracefully
        call http_communicator%term()
        call simple_end('**** SIMPLE_STREAM_ASSIGN_OPTICS NORMAL STOP ****')

        contains
        
            subroutine communicator_init()
                call http_communicator%add_to_json("stage",                    "initialising")
                call http_communicator%add_to_json("micrographs_assigned ",    0)
                call http_communicator%add_to_json("optics_groups_assigned",   0)
                call http_communicator%add_to_json("last_micrograph_imported", "")
                call json%create_array(optics_assignments, "optics_assignments")
                call http_communicator%add_to_json(optics_assignments)
            end subroutine communicator_init

    end subroutine exec_stream_p02_assign_optics

end module simple_stream_p02_assign_optics
