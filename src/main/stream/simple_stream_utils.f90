!@descr: various stream utilities
module simple_stream_utils
use simple_core_module_api
use simple_cmdline,             only: cmdline
use simple_qsys_env,            only: qsys_env
use simple_rec_list,            only: rec_list, project_rec
use simple_sp_project,          only: sp_project
use simple_parameters,          only: parameters
use simple_syslib,              only: get_current_rss_bytes, get_peak_rss_bytes
implicit none
#include "simple_local_flags.inc"

contains

    ! The process's resident memory, logged at @p phase with the time, so a log shows both what a
    ! step holds and how long it took.
    subroutine log_rss( phase )
        use, intrinsic :: iso_c_binding, only: c_int64_t, c_double
        character(len=*), intent(in) :: phase
        integer(c_int64_t) :: current_rss, peak_rss
        real(c_double)     :: current_mib, peak_mib
        current_rss = get_current_rss_bytes()
        peak_rss    = get_peak_rss_bytes()
        if( current_rss < 0_c_int64_t .or. peak_rss < 0_c_int64_t ) return
        current_mib = real(current_rss, c_double) / 1048576.0_c_double
        peak_mib    = real(peak_rss,    c_double) / 1048576.0_c_double
        write(logfhandle,'(A,A,A,F10.1,A,F10.1,A,A)') '>>> RSS ', trim(phase), ': current=', current_mib, ' MiB peak=',&
            &peak_mib, ' MiB AT: ', cast_time_char(simple_gettime())
        call flush(logfhandle)
    end subroutine log_rss

    subroutine terminate_stream( params, msg )
        class(parameters), intent(in) :: params
        character(len=*), intent(in) :: msg
        if(trim(params%async).eq.'yes') then
            if( file_exists(TERM_STREAM) ) call simple_end('**** '//trim(msg)//' ****', print_simple=.false.)
        endif
    end subroutine terminate_stream

    subroutine create_stream_project( spproj, cline, projname )
        type(sp_project),  intent(inout) :: spproj
        type(cmdline),     intent(inout) :: cline
        type(string),         intent(in) :: projname
        if( .not. cline%defined('projfile') ) then
            call cline%set('projname', projname)
            call cline%set('projfile', string(projname%to_char() // METADATA_EXT))
        endif
        call spproj%update_projinfo(cline)
        call spproj%update_compenv(cline)
        call spproj%write()
    end subroutine create_stream_project

    !> .true. when the stage in @p dir_upstream has gone idle or stopped (its STREAM_IDLE or
    !! STREAM_FINISHED marker): it hands on nothing more until it removes its idle marker.
    logical function upstream_done( dir_upstream )
        class(string), intent(in) :: dir_upstream
        upstream_done = file_exists(dir_upstream//'/'//STREAM_IDLE_MARKER)
        if( .not. upstream_done ) upstream_done = file_exists(dir_upstream//'/'//STREAM_FINISHED_MARKER)
    end function upstream_done

    subroutine init_stream_qenv( params, qenv, envvar )
        type(parameters), intent(inout) :: params
        type(qsys_env),   intent(inout) :: qenv
        type(string),        intent(in) :: envvar
        character(len=STDLEN)           :: chunk_part_env
        integer                         :: envlen
        call get_environment_variable(envvar%to_char(), chunk_part_env, envlen)
        if(envlen > 0) then
            call qenv%new(params, 1, stream=.true., exec_bin=string('simple_exec'), qsys_partition=string(trim(chunk_part_env)))
        else
            call qenv%new(params, 1, stream=.true., exec_bin=string('simple_exec'))
        end if
    end subroutine init_stream_qenv 

    subroutine import_new_projects( project_list, projects, n_mics_imported, n_ptcls_imported, ignore_ptcls, check_state )
        type(rec_list),              intent(inout) :: project_list
        type(string),   allocatable, intent(inout) :: projects(:)
        integer,                     intent(inout) :: n_mics_imported, n_ptcls_imported
        logical,        optional,       intent(in) :: ignore_ptcls, check_state
        type(project_rec)                          :: prec
        type(sp_project)                           :: spproj
        type(string)                               :: projabspath
        integer :: iproj, imic
        logical :: l_check_state, l_ignore_ptcls
        l_check_state  = .false.
        l_ignore_ptcls = .false.
        if( present(check_state)  ) l_check_state  = check_state
        if( present(ignore_ptcls) ) l_ignore_ptcls = ignore_ptcls
        if( .not.allocated(projects) ) return
        if( size(projects) == 0      ) return
        do iproj=1, size(projects)
            call spproj%read(projects(iproj))
            ! because pick_extract purges state=0 and nptcls=0 mics,
            ! all mics can be assumed associated with particles
            if( spproj%os_mic%get_noris() == 0) then
                write(logfhandle, *) "ERROR: mic noris 0", projects(iproj)%to_char()
                call spproj%kill()
                cycle
            end if
            if( .not. l_ignore_ptcls .and. spproj%os_stk%get_noris() == 0) then
                write(logfhandle, *) "ERROR: stk noris 0", projects(iproj)%to_char()
                call spproj%kill()
                cycle
            end if
            projabspath = simple_abspath(projects(iproj))
            do imic = 1, spproj%os_mic%get_noris()
                if( l_check_state .and. spproj%os_mic%get_int(imic, 'state') <= 0 ) cycle
                prec%id          = project_list%size() + 1
                prec%projname    = projabspath
                prec%micind      = imic
                prec%nptcls      = spproj%os_mic%get_int(imic,'nptcls')
                prec%nptcls_sel  = prec%nptcls
                prec%included    = .false.
                n_mics_imported  = n_mics_imported + 1
                n_ptcls_imported = n_ptcls_imported + prec%nptcls
                call project_list%push_back(prec)
            enddo
            ! cleanup
            call spproj%kill()
        end do
    end subroutine import_new_projects

    function get_latest_optics_map_id(optics_dir) result (lastmap)
        class(string), optional, intent(in)   :: optics_dir
        type(string), allocatable :: map_list(:)
        type(string) :: map_str, map_i_str
        integer      :: imap, prefix_len, testmap, lastmap
        lastmap = 0
        if(optics_dir .ne. "") then
            if(dir_exists(optics_dir)) call simple_list_files(optics_dir%to_char()//'/'// OPTICS_MAP_PREFIX //'*'//TXT_EXT, map_list)
        endif
        if(allocated(map_list)) then
            prefix_len = len(optics_dir%to_char() // '/' // OPTICS_MAP_PREFIX) + 1
            do imap=1, size(map_list)
                map_str   = map_list(imap)%to_char([prefix_len,map_list(imap)%strlen_trim()])
                map_i_str = swap_suffix(map_str, "", TXT_EXT)
                testmap   = map_i_str%to_int()
                if(testmap > lastmap) lastmap = testmap
                call map_str%kill
                call map_i_str%kill
            enddo
            deallocate(map_list)
        endif
    end function get_latest_optics_map_id

end module simple_stream_utils
