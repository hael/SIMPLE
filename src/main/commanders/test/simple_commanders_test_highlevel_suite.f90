!@descr: the suite runners of simple_commanders_test_highlevel: test=nightly, test=profile and test=identity
! nightly and profile run test programs as children of the test dispatcher (simple_test_dispatch),
! each in its own directory, and collect the gate tables (metrics.tsv) the tests write; identity
! compares two run directories of one entry over the outputs the registry declares. The children run
! through the test routers, which is why this is a submodule: the routers use the parent module.
submodule (simple_commanders_test_highlevel) simple_commanders_test_highlevel_suite
use simple_test_registry,      only: test_entry, get_test_registry, find_test_entry, ENTRY_GATE
use simple_test_dispatch,      only: test_dispatcher, test_execution, set_test_environment, exec_state_name, &
                                    &EXEC_PASSED, EXEC_INTERRUPTED, EXEC_NOT_RUN
use simple_test_profile_stats, only: gate_row, metric_dist, read_gate_table, profile_table, write_profile, &
                                    &read_profile, compare_profiles
implicit none
#include "simple_local_flags.inc"

integer, parameter :: GRACE_S = 60 !< seconds between SIGTERM and SIGKILL of a child past its deadline

contains

!> simple_test_exec test=nightly [ncunits=<cores>] [dir=<output directory>]
module subroutine exec_test_nightly( self, cline )
    class(commander_test_nightly), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(test_entry), allocatable :: entries(:)
    type(test_dispatcher) :: dispatcher
    type(test_execution)  :: exec
    character(len=:), allocatable :: root, simple_path
    integer :: i
    logical :: interrupted, all_ok
    call set_test_environment(simple_path)
    root    = output_root(cline, 'nightly')
    entries = get_test_registry()
    call dispatcher%new(core_budget(cline), GRACE_S)
    do i = 1, size(entries)
        if( entries(i)%kind /= ENTRY_GATE ) cycle
        exec          = test_execution()
        exec%label    = entries(i)%name
        exec%prg      = entries(i)%prg
        exec%args     = entries(i)%args
        exec%dir      = root//'/'//trim(entries(i)%name)
        exec%nthr     = entries(i)%nthr
        exec%ncpu     = entries(i)%ncpu
        exec%mem_gb   = entries(i)%mem_gb
        exec%disk_gb  = entries(i)%disk_gb
        exec%deadline = entries(i)%deadline
        call dispatcher%add(exec)
    enddo
    call log_start('nightly', root, simple_path, dispatcher%get_nexecs())
    call dispatcher%run(route_test, interrupted)
    call write_report(dispatcher, root//'/nightly_report.tsv', all_ok)
    call dispatcher%kill
    if( interrupted ) THROW_HARD('test=nightly was interrupted; report in '//root)
    if( .not. all_ok ) THROW_HARD('test=nightly: an entry or a gate failed; report in '//root)
    write(logfhandle,'(a)') 'PASS: every entry of the nightly suite passed; report in '//root
    call simple_end('**** SIMPLE_TEST_NIGHTLY NORMAL STOP ****')
end subroutine exec_test_nightly

!> simple_test_exec test=profile target=<entry[,entry...]|all> [suite=<sub-suite>] nseeds=<n> [nthr=<list>]
!!     [ncunits=<cores>] [dir=<output directory>]
!! or, to compare two profiles, simple_test_exec test=profile dir1=<profile> dir2=<profile>
module subroutine exec_test_profile( self, cline )
    class(commander_test_profile), intent(inout) :: self
    class(cmdline),                intent(inout) :: cline
    type(test_entry), allocatable :: entries(:)
    type(test_dispatcher) :: dispatcher
    type(test_execution)  :: exec
    character(len=:), allocatable :: root, simple_path, target, suite
    logical,          allocatable :: selected(:)
    character(len=64)  :: tag
    character(len=len(exec%dir)+12), allocatable :: tables(:)
    logical, allocatable :: oks(:)
    integer, allocatable :: nthrs(:)
    integer :: i, ient, iseed, ithr, nseeds, nexec, ntables
    real    :: scale
    logical :: interrupted, all_ok
    type(metric_dist), allocatable :: dists(:)
    type(test_execution) :: done
    type(string) :: val
    if( cline%defined('dir1') .and. cline%defined('dir2') )then
        call compare_profile_dirs(cline%get_carg('dir1'), cline%get_carg('dir2'))
        return
    endif
    if( .not. cline%defined('target') ) THROW_HARD('test=profile needs target=<entry[,entry...]|all>')
    ! the interface's default
    if( .not. cline%defined('nseeds') ) call cline%set('nseeds', 20)
    val    = cline%get_carg('target')
    target = val%to_char()
    nseeds = cline%get_iarg('nseeds')
    if( nseeds < 1 ) THROW_HARD('nseeds must be at least 1')
    suite = ''
    if( cline%defined('suite') )then
        val   = cline%get_carg('suite')
        suite = val%to_char()
    endif
    entries = get_test_registry()
    selected = selected_entries(entries, target)
    ! a sub-suite names one entry's sub-suites
    if( len(suite) > 0 .and. count(selected) /= 1 ) THROW_HARD('suite= applies to one entry: give target=<entry>')
    call set_test_environment(simple_path)
    root = output_root(cline, 'profile')
    call dispatcher%new(core_budget(cline), GRACE_S)
    do ient = 1, size(entries)
        if( .not. selected(ient) ) cycle
        nthrs = thread_list(cline, entries(ient)%nthr)
        do ithr = 1, size(nthrs)
            do iseed = 0, nseeds - 1
                write(tag,'(a,i0,a,i0)') '_t', nthrs(ithr), '_s', iseed
                exec             = test_execution()
                exec%label       = trim(entries(ient)%name)//trim(tag)
                exec%prg         = entries(ient)%prg
                exec%args        = thread_args(entries(ient), nthrs(ithr), suite)
                if( len(root) + 1 + len_trim(exec%label) > len(exec%dir) )&
                    &THROW_HARD('test=profile: the run directory name is too long')
                exec%dir         = root//'/'//trim(exec%label)
                exec%seed_offset = iseed
                exec%nthr        = nthrs(ithr)
                ! the reservation follows the threads of the master and of the workers, the deadline
                ! the threads inversely; the memory and the disk do not scale
                scale            = real(nthrs(ithr)) / real(entries(ient)%nthr)
                exec%ncpu        = scaled_ncpu(entries(ient), nthrs(ithr))
                exec%mem_gb      = entries(ient)%mem_gb
                exec%disk_gb     = entries(ient)%disk_gb
                exec%deadline    = nint(real(entries(ient)%deadline) * max(1., 1. / scale))
                call dispatcher%add(exec)
            enddo
        enddo
    enddo
    nexec = dispatcher%get_nexecs()
    allocate(tables(nexec), oks(nexec))
    call log_start('profile', root, simple_path, nexec)
    call dispatcher%run(route_test, interrupted)
    call write_report(dispatcher, root//'/executions.tsv', all_ok)
    ! one distribution table per entry and thread count
    do ient = 1, size(entries)
        if( .not. selected(ient) ) cycle
        nthrs = thread_list(cline, entries(ient)%nthr)
        do ithr = 1, size(nthrs)
            write(tag,'(a,i0)') '_t', nthrs(ithr)
            ntables = 0
            do i = 1, nexec
                done = dispatcher%get_execution(i)
                if( index(done%label, trim(entries(ient)%name)//trim(tag)//'_s') /= 1 ) cycle
                ! an execution the interrupt stopped, or never started, has no verdict
                if( done%state == EXEC_INTERRUPTED .or. done%state == EXEC_NOT_RUN ) cycle
                ntables = ntables + 1
                tables(ntables) = trim(done%dir)//'/metrics.tsv'
                oks(ntables)    = done%state == EXEC_PASSED
            enddo
            call profile_table(tables(1:ntables), oks(1:ntables), dists)
            call write_profile(root//'/profile_'//trim(entries(ient)%name)//trim(tag)//'.tsv', dists)
        enddo
    enddo
    call dispatcher%kill
    if( interrupted ) THROW_HARD('test=profile was interrupted; tables in '//root)
    write(logfhandle,'(a)') 'test=profile done; tables in '//root
    call simple_end('**** SIMPLE_TEST_PROFILE NORMAL STOP ****')

  contains

    !> the profile tables of two profile directories, compared where both have one
    subroutine compare_profile_dirs( dir1, dir2 )
        use simple_fileio, only: basename
        type(string), intent(in) :: dir1, dir2
        type(string), allocatable :: list(:)
        type(metric_dist), allocatable :: a(:), b(:)
        type(string) :: bname
        integer :: j, ncompared
        logical :: l_exists
        call simple_list_files(dir1%to_char()//'/profile_*.tsv', list)
        ncompared = 0
        do j = 1, size(list)
            bname = basename(list(j))
            inquire(file=dir2%to_char()//'/'//bname%to_char(), exist=l_exists)
            if( .not. l_exists ) cycle
            call read_profile(list(j)%to_char(), a)
            call read_profile(dir2%to_char()//'/'//bname%to_char(), b)
            call compare_profiles(a, b, 'comparison_'//bname%to_char())
            write(logfhandle,'(a)') '>>> compared '//bname%to_char()//': comparison_'//bname%to_char()
            ncompared = ncompared + 1
        enddo
        if( ncompared == 0 ) THROW_HARD('the two directories hold no profile table in common')
        call simple_end('**** SIMPLE_TEST_PROFILE COMPARISON NORMAL STOP ****')
    end subroutine compare_profile_dirs

end subroutine exec_test_profile

!> simple_test_exec test=identity target=<entry> dir1=<run directory> dir2=<run directory>
module subroutine exec_test_identity( self, cline )
    use simple_test_identity, only: compare_run_dirs
    class(commander_test_identity), intent(inout) :: self
    class(cmdline),                 intent(inout) :: cline
    type(test_entry), allocatable :: entries(:)
    type(string) :: dir1, dir2, target
    integer :: ient, nsame, ndiff
    if( .not. (cline%defined('target') .and. cline%defined('dir1') .and. cline%defined('dir2')) )&
        &THROW_HARD('test=identity needs target=<entry> dir1=<run directory> dir2=<run directory>')
    entries = get_test_registry()
    target  = cline%get_carg('target')
    ient    = find_test_entry(entries, target%to_char())
    if( ient == 0 ) THROW_HARD('no registry entry '//target%to_char())
    if( len_trim(entries(ient)%outputs) == 0 ) THROW_HARD('the registry declares no outputs of '//trim(entries(ient)%name))
    dir1 = simple_abspath(cline%get_carg('dir1'))
    dir2 = simple_abspath(cline%get_carg('dir2'))
    call compare_run_dirs(dir1%to_char(), dir2%to_char(), trim(entries(ient)%outputs), 'identity_report.txt', nsame, ndiff)
    write(logfhandle,'(a,i0,a,i0)') '>>> identity of '//trim(entries(ient)%name)//': files the same ', nsame, &
        &', differing or missing ', ndiff
    if( ndiff > 0 .or. nsame == 0 ) THROW_HARD('the run directories are not identical; see identity_report.txt')
    write(logfhandle,'(a)') 'PASS: the declared outputs are identical'
    call simple_end('**** SIMPLE_TEST_IDENTITY NORMAL STOP ****')
end subroutine exec_test_identity

!> runs one test program in this process, as simple_test_exec does
subroutine route_test( prg, cline )
    use simple_test_exec_class,     only: exec_test_class_commander
    use simple_test_exec_highlevel, only: exec_test_highlevel_commander
    use simple_test_exec_single,    only: exec_test_single_commander
    use simple_test_exec_stream,    only: exec_test_stream_commander
    character(len=*), intent(in)    :: prg
    class(cmdline),   intent(inout) :: cline
    logical :: l_silent, l_did_execute
    l_did_execute = .false.
    call exec_test_class_commander(    prg, cline, l_silent, l_did_execute)
    call exec_test_highlevel_commander(prg, cline, l_silent, l_did_execute)
    call exec_test_single_commander(   prg, cline, l_silent, l_did_execute)
    call exec_test_stream_commander(   prg, cline, l_silent, l_did_execute)
    if( .not. l_did_execute ) THROW_HARD('no test program '//prg)
end subroutine route_test

!> One row per execution and gate metric (an execution without a gate table has one row), with
!! its outcome, wall time, peak CPU use and memory. all_ok: every execution passed and every
!! metric with a floor passed.
subroutine write_report( dispatcher, fname, all_ok )
    class(test_dispatcher), intent(in)  :: dispatcher
    character(len=*),       intent(in)  :: fname
    logical,                intent(out) :: all_ok
    character(len=1), parameter :: T = achar(9)
    type(test_execution) :: exec
    type(gate_row), allocatable :: rows(:)
    character(len=512) :: head
    integer :: funit, i, j
    all_ok = .true.
    open(newunit=funit, file=fname, action='write', status='replace')
    write(funit,'(a)') 'execution'//T//'outcome'//T//'exit_code'//T//'wall_s'//T//'peak_cpu'//T//'reserved_cpu'//T//&
        &'peak_rss_GiB'//T//'seed_offset'//T//'nthr'//T//'metric'//T//'value'//T//'floor'//T//'pass'//T//'sense'
    do i = 1, dispatcher%get_nexecs()
        exec = dispatcher%get_execution(i)
        if( exec%state /= EXEC_PASSED ) all_ok = .false.
        write(head,'(a,a,a,a,i0,a,a,a,a,a,i0,a,a,a,i0,a,i0)') trim(exec%label), T, exec_state_name(exec%state), T, &
            &exec%exit_code, T, fixed(exec%wall, 1), T, fixed(exec%peak_cpu, 2), T, exec%ncpu, T, &
            &fixed(exec%peak_rss_gb, 2), T, exec%seed_offset, T, exec%nthr
        call read_gate_table(trim(exec%dir)//'/metrics.tsv', rows)
        if( size(rows) == 0 ) write(funit,'(a)') trim(head)//T//'-'//T//'-'//T//'-'//T//'-'//T//'-'
        do j = 1, size(rows)
            if( .not. rows(j)%pass ) all_ok = .false.
            write(funit,'(a,a,a,a,es13.6,a,es13.6,a,a,a,a)') trim(head), T, trim(rows(j)%name), T, rows(j)%value, T, &
                &rows(j)%floor, T, trim(merge('yes', 'no ', rows(j)%pass)), T, trim(rows(j)%sense)
        enddo
        write(logfhandle,'(a,a,f10.1,a,f6.2,a,i0,a)') trim(exec%label)//': '//exec_state_name(exec%state), ', ', &
            &exec%wall, ' s, peak CPU ', exec%peak_cpu, ' (reserved ', exec%ncpu, ')'
    enddo
    close(funit)
end subroutine write_report

!> x with ndec decimals and its leading zero
function fixed( x, ndec ) result( str )
    real,    intent(in) :: x
    integer, intent(in) :: ndec
    character(len=:), allocatable :: str
    character(len=32) :: buf
    character(len=16) :: fmt
    write(fmt,'(a,i0,a)') '(f32.', ndec, ')'
    write(buf,fmt) x
    str = trim(adjustl(buf))
end function fixed

!> the directory of a run: dir=, or <kind>_<date>_<time> in the working directory
function output_root( cline, kind ) result( root )
    class(cmdline),   intent(in)  :: cline
    character(len=*), intent(in)  :: kind
    character(len=:), allocatable :: root
    character(len=8)  :: date
    character(len=10) :: time
    type(string) :: dir, cwd
    if( cline%defined('dir') )then
        dir = cline%get_carg('dir')
    else
        call date_and_time(date=date, time=time)
        dir = kind//'_'//date//'_'//time(1:6)
    endif
    call simple_mkdir(dir)
    cwd  = simple_abspath(dir)
    root = cwd%to_char()
end function output_root

!> ncunits=, or the processors of the machine
integer function core_budget( cline )
    !$ use omp_lib, only: omp_get_num_procs
    class(cmdline), intent(in) :: cline
    core_budget = 1
    !$ core_budget = omp_get_num_procs()
    if( cline%defined('ncunits') ) core_budget = cline%get_iarg('ncunits')
end function core_budget

!> nthr=<n[,n...]>, or the entry's own thread count
function thread_list( cline, nthr_entry ) result( nthrs )
    use simple_string_utils, only: parsestr
    class(cmdline), intent(in) :: cline
    integer,        intent(in) :: nthr_entry
    integer, allocatable :: nthrs(:)
    character(len=256) :: list
    character(len=16)  :: items(32)
    type(string) :: val
    integer :: n, i, ios
    if( .not. cline%defined('nthr') )then
        nthrs = [nthr_entry]
        return
    endif
    ! one count is a number on the command line, a list a text
    list = ''
    val  = cline%get_carg('nthr')
    if( val%is_allocated() ) list = val%to_char()
    if( len_trim(list) == 0 ) list = int2str(cline%get_iarg('nthr'))
    call parsestr(list, ',', items, n)
    allocate(nthrs(n))
    do i = 1, n
        read(items(i),*,iostat=ios) nthrs(i)
        if( ios /= 0 .or. nthrs(i) < 1 ) THROW_HARD('nthr= is a list of thread counts: '//trim(list))
    enddo
end function thread_list

!> the registry entries target names: all, or a comma-separated list of entry names
function selected_entries( entries, target ) result( selected )
    use simple_string_utils, only: parsestr
    type(test_entry), intent(in) :: entries(:)
    character(len=*), intent(in) :: target
    logical, allocatable :: selected(:)
    character(len=1024) :: list
    character(len=64)   :: names(64)
    integer :: n, i, ient
    allocate(selected(size(entries)), source=trim(target) == 'all')
    if( trim(target) == 'all' ) return
    list = target
    call parsestr(list, ',', names, n)
    do i = 1, n
        ient = find_test_entry(entries, trim(names(i)))
        if( ient == 0 ) THROW_HARD('no registry entry '//trim(names(i)))
        selected(ient) = .true.
    enddo
end function selected_entries

!> The command line of an entry at nthr threads: a thread count the entry passes as nthr= scales
!! with nthr / its CTest threads; suite= becomes suite, when one is given.
function thread_args( entry, nthr, suite ) result( args )
    type(test_entry), intent(in) :: entry
    integer,          intent(in) :: nthr
    character(len=*), intent(in) :: suite
    character(len=256) :: args
    character(len=16)  :: val
    integer :: nthr_arg, ios
    args = entry%args
    val  = key_value(args, 'nthr')
    if( len_trim(val) > 0 )then
        read(val,*,iostat=ios) nthr_arg
        if( ios /= 0 ) THROW_HARD('the registry row of '//trim(entry%name)//' has a malformed nthr=')
        args = with_key(args, 'nthr', int2str(max(1, nint(real(nthr_arg * nthr) / real(entry%nthr)))))
    endif
    if( len_trim(suite) > 0 ) args = with_key(args, 'suite', suite)
end function thread_args

!> The CPU reservation of an entry run at nthr threads: its master's share (the reservation less
!! its workers' threads) scales with nthr / its CTest threads when the master's threads follow the
!! thread count, and so do the threads of each worker; an entry whose threads are its own keeps its
!! reservation.
integer function scaled_ncpu( entry, nthr )
    type(test_entry), intent(in) :: entry
    integer,          intent(in) :: nthr
    real    :: scale
    integer :: master, wthr
    scaled_ncpu = entry%ncpu
    if( nthr == entry%nthr .or. .not. entry%scales ) return
    scale       = real(nthr) / real(entry%nthr)
    master      = max(1, entry%ncpu - entry%nworkers * entry%wthr)
    wthr        = 0
    if( entry%nworkers > 0 ) wthr = max(1, nint(real(entry%wthr) * scale))
    scaled_ncpu = max(1, ceiling(real(master) * scale - 0.1) + entry%nworkers * wthr)
end function scaled_ncpu

!> the value of key= in a blank-separated command line, blank when absent
function key_value( args, key ) result( val )
    character(len=*), intent(in) :: args, key
    character(len=16) :: val
    integer :: pos, pend
    val = ''
    pos = index(' '//args, ' '//key//'=')
    if( pos == 0 ) return
    pos  = pos + len(key) + 1
    pend = index(args(pos:)//' ', ' ')
    val  = args(pos:pos+pend-2)
end function key_value

!> args with key= set to val, replaced where it is, appended otherwise
function with_key( args, key, val ) result( out )
    character(len=*), intent(in) :: args, key, val
    character(len=256) :: out
    integer :: pos, pend
    pos = index(' '//args, ' '//key//'=')
    if( pos == 0 )then
        out = trim(adjustl(trim(args)//' '//key//'='//val))
    else
        pend = index(args(pos:)//' ', ' ')
        out  = args(:pos-1)//key//'='//val//args(pos+pend-1:)
    endif
end function with_key

subroutine log_start( kind, root, simple_path, nexec )
    character(len=*), intent(in) :: kind, root, simple_path
    integer,          intent(in) :: nexec
    write(logfhandle,'(a,i0,a)') '>>> test='//kind//': ', nexec, ' executions in '//root
    write(logfhandle,'(a)')      '>>> SIMPLE_PATH of the children and their workers: '//simple_path
end subroutine log_start

end submodule simple_commanders_test_highlevel_suite
