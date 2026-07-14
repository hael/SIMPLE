!@descr: opt-in process memory telemetry for all SIMPLE commanders and processing phases
module simple_memory_monitor
use, intrinsic :: iso_c_binding, only: c_double, c_int, c_int64_t
use simple_cmdline, only: cmdline
use simple_defs,   only: STDLEN, logfhandle
use simple_string, only: string
use simple_syslib, only: get_current_rss_bytes, get_peak_rss_bytes, get_process_id, &
    &memory_monitor_start_c, memory_monitor_report_c, memory_monitor_stop_c
implicit none
private

public :: mem_monitor_init, mem_monitor_finish, mem_monitor_is_enabled
public :: mem_report

integer, parameter :: DEFAULT_REPORT_INTERVAL_SECONDS = 30

logical, save :: monitor_enabled = .false.
integer, save :: report_interval_seconds = DEFAULT_REPORT_INTERVAL_SECONDS
integer, save :: process_id = 0
integer, save :: run_id = 0
character(len=STDLEN), save :: program_label = ''
character(len=STDLEN), save :: telemetry_filename = ''

contains

    !> Enable automatic process-memory telemetry for every SIMPLE commander.
    !! memreport=yes starts a native POSIX sampler. memreport_interval is the
    !! interval in seconds. Explicit mem_report calls add named phase samples.
    subroutine mem_monitor_init( cline, label )
        class(cmdline),   intent(in) :: cline
        character(len=*), intent(in) :: label
        type(string) :: mode_arg
        character(len=STDLEN) :: mode
        integer :: requested_interval
        if( monitor_enabled ) call mem_monitor_finish
        if( .not. cline%defined('memreport') ) return
        mode_arg = cline%get_carg('memreport')
        mode = trim(adjustl(mode_arg%to_char()))
        select case(trim(mode))
            case('no')
                return
            case('yes')
                ! continue below
            case DEFAULT
                write(logfhandle,'(A,A,A)') '>>> WARNING: memreport=', trim(mode), &
                    &' is unsupported; expected yes or no, memory reporting disabled'
                return
        end select
        requested_interval = DEFAULT_REPORT_INTERVAL_SECONDS
        if( cline%defined('memreport_interval') ) requested_interval = cline%get_iarg('memreport_interval')
        if( requested_interval < 1 )then
            write(logfhandle,'(A,I0,A)') '>>> WARNING: memreport_interval=', requested_interval, &
                &' is invalid; using 1 second'
            requested_interval = 1
        endif
        call start_monitor(trim(label), requested_interval)
    end subroutine mem_monitor_init

    !> Emit a named memory sample immediately. Automatic periodic reporting
    !! continues independently, so this is only needed for phase attribution.
    subroutine mem_report( phase )
        character(len=*), intent(in) :: phase
        integer(c_int64_t) :: current_rss, peak_rss
        integer(c_int) :: rc
        real(c_double) :: elapsed
        if( .not. monitor_enabled ) return
        rc = memory_monitor_report_c(trim(phase), int(len_trim(phase),c_int), &
            &elapsed, current_rss, peak_rss)
        if( rc == 0_c_int ) call log_sample(trim(phase), elapsed, current_rss, peak_rss)
    end subroutine mem_report

    !> Stop periodic sampling, emit the final sample, and close telemetry.
    subroutine mem_monitor_finish
        integer(c_int64_t) :: current_rss, peak_rss
        integer(c_int) :: rc
        real(c_double) :: elapsed
        if( .not. monitor_enabled ) return
        rc = memory_monitor_stop_c(elapsed, current_rss, peak_rss)
        if( rc == 0_c_int ) call log_sample('finish', elapsed, current_rss, peak_rss)
        monitor_enabled = .false.
    end subroutine mem_monitor_finish

    logical function mem_monitor_is_enabled()
        mem_monitor_is_enabled = monitor_enabled
    end function mem_monitor_is_enabled

    subroutine start_monitor( label, interval_seconds )
        character(len=*), intent(in) :: label
        integer,          intent(in) :: interval_seconds
        integer(c_int64_t) :: current_rss, peak_rss
        integer(c_int) :: rc
        current_rss = get_current_rss_bytes()
        peak_rss    = get_peak_rss_bytes()
        if( current_rss < 0_c_int64_t .or. peak_rss < 0_c_int64_t )then
            write(logfhandle,'(A)') '>>> WARNING: process RSS reporting is unsupported on this platform'
            return
        endif
        process_id = get_process_id()
        run_id = run_id + 1
        program_label = ''
        if( len_trim(label) > 0 )then
            program_label(1:min(len_trim(label),len(program_label))) = &
                &label(1:min(len_trim(label),len(program_label)))
        endif
        report_interval_seconds = interval_seconds
        write(telemetry_filename,'(A,I0,A)') 'memory_usage_', process_id, '.csv'
        rc = memory_monitor_start_c(trim(telemetry_filename), int(len_trim(telemetry_filename),c_int), &
            &trim(program_label), int(len_trim(program_label),c_int), int(report_interval_seconds,c_int), &
            &int(run_id,c_int))
        if( rc /= 0_c_int )then
            write(logfhandle,'(A,A,A,I0)') '>>> WARNING: could not start memory telemetry file ', &
                &trim(telemetry_filename), '; status=', rc
            return
        endif
        monitor_enabled = .true.
        write(logfhandle,'(A,A,A,I0,A)') '>>> MEMORY telemetry file: ', trim(telemetry_filename), &
            &' (interval=', report_interval_seconds, ' seconds)'
        call mem_report('start')
    end subroutine start_monitor

    subroutine log_sample( phase, elapsed, current_rss, peak_rss )
        character(len=*), intent(in) :: phase
        real(c_double),   intent(in) :: elapsed
        integer(c_int64_t), intent(in) :: current_rss, peak_rss
        real(c_double) :: current_mib, peak_mib
        current_mib = real(current_rss, c_double) / 1048576.0_c_double
        peak_mib    = real(peak_rss,    c_double) / 1048576.0_c_double
        write(logfhandle,'(A,F10.3,A,I0,A,A,A,A,A,F10.1,A,F10.1,A)') &
            &'>>> MEMORY elapsed=', elapsed, ' s pid=', process_id, ' program=', trim(program_label), &
            &' phase=', trim(phase), ' current=', current_mib, ' MiB peak=', peak_mib, ' MiB'
        flush(logfhandle)
    end subroutine log_sample

end module simple_memory_monitor
