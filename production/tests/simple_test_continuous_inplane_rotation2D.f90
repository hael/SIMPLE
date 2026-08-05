program simple_test_continuous_inplane_rotation2D
    implicit none

    character(len=4096) :: build_dir, vol1, src_proj, logdir, only_tests
    character(len=64)   :: threads, box_crop, diagnostic
    character(len=4096) :: off_proj, on_proj
    integer :: failures
    logical :: volume_exists, project_exists

    build_dir  = '.'
    vol1       = ''
    src_proj   = ''
    logdir     = '.'
    threads    = '24'
    box_crop   = '88'
    diagnostic = 'no'
    only_tests = 'all'
    failures   = 0

    call parse_arguments(build_dir, vol1, src_proj, logdir, threads, box_crop, &
        &diagnostic, only_tests)
    call make_directory(logdir)

    if( len_trim(vol1) == 0 )then
        call find_newest_file(build_dir, '1JYX_v4.mrc', logdir, 'input_vol1.txt', vol1)
        if( len_trim(vol1) == 0 )then
            call find_newest_file(build_dir, '1JYX_v3.mrc', logdir, 'input_vol1.txt', vol1)
        endif
    endif
    if( len_trim(src_proj) == 0 )then
        call find_newest_file(build_dir, 'onejyx_abinitio_v3.simple', logdir, 'input_proj.txt', src_proj)
    endif

    if( len_trim(vol1) == 0 )then
        write(*,'(a)') 'NOTICE: no 1JYX volume was found under '//trim(build_dir)//'.'
        write(*,'(a)') 'NOTICE: run simple_test_sgd_base_suite first; its V3/V4 fixture tests create the volume.'
        error stop 1
    else
        write(*,'(a)') 'USING VOLUME: '//trim(vol1)
        inquire(file=trim(vol1), exist=volume_exists)
        if( .not. volume_exists )then
            write(*,'(a)') 'NOTICE: volume was not found: '//trim(vol1)
            write(*,'(a)') 'NOTICE: run simple_test_sgd_base_suite first, or run the shorter V3/V4 fixture test.'
            error stop 1
        endif
    endif
    if( len_trim(src_proj) == 0 )then
        write(*,'(a)') 'NOTICE: no 1JYX V3 project was found under '//trim(build_dir)//'.'
        write(*,'(a)') 'NOTICE: run simple_test_sgd_base_suite first to create the 1JYX V3 project.'
        error stop 1
    else
        write(*,'(a)') 'USING SOURCE PROJECT: '//trim(src_proj)
        inquire(file=trim(src_proj), exist=project_exists)
        if( .not. project_exists )then
            write(*,'(a)') 'NOTICE: source project was not found: '//trim(src_proj)
            write(*,'(a)') 'NOTICE: run simple_test_sgd_base_suite first to create the '// &
                &'1JYX V3 project.'
            error stop 1
        endif
    endif

    if( test_selected(only_tests, 'route_identity') )then
        call run_test('route_identity', &
            &'simple_test_continuous_inplane_rotation2D_route_identity vol1='//trim(quote(vol1))// &
            &' mskdiam=120 lp=8 smpd=1.3', logdir, failures)
    endif

    if( test_selected(only_tests, 'stage1_validation') )then
        call run_test('stage1_validation', &
            &'simple_test_continuous_inplane_rotation2D_stage1_validation vol1='//trim(quote(vol1))// &
            &' mskdiam=120 smpd=1.3 lp=8 lpstop=2.7 angerr=37 xsh=2 ysh=-1.5', &
            logdir, failures)
    endif

    if( test_selected(only_tests, 'euclid_off_workflow') .or. &
        &test_selected(only_tests, 'euclid_off_metadata') )then
        call run_workflow('off', 'no', src_proj, threads, diagnostic, logdir, failures)
        call newest_project(build_dir, logdir, 'off_project.txt', off_proj, failures)
        if( test_selected(only_tests, 'euclid_off_metadata') .and. len_trim(off_proj) > 0 )then
            call run_test('euclid_off_metadata', &
                &'simple_test_continuous_inplane_rotation2D_metadata projfile='//trim(quote(off_proj))// &
                &' mskdiam=120 lp=9.853 smpd=1.3 box_crop='//trim(box_crop)// &
                &' inpl_cont=no', logdir, failures)
        endif
    endif

    if( test_selected(only_tests, 'euclid_on_workflow') .or. &
        &test_selected(only_tests, 'euclid_on_metadata') )then
        call run_workflow('on', 'yes', src_proj, threads, diagnostic, logdir, failures)
        call newest_project(build_dir, logdir, 'on_project.txt', on_proj, failures)
        if( test_selected(only_tests, 'euclid_on_metadata') .and. len_trim(on_proj) > 0 )then
            call run_test('euclid_on_metadata', &
                &'simple_test_continuous_inplane_rotation2D_metadata projfile='//trim(quote(on_proj))// &
                &' mskdiam=120 lp=9.853 smpd=1.3 box_crop='//trim(box_crop)// &
                &' inpl_cont=yes', logdir, failures)
        endif
    endif

    ! CC regression is intentionally disabled until CC derivatives support
    ! continuous in-plane refinement.  Re-enable it with a CC-specific fixture.

    if( failures /= 0 )then
        write(*,'(a,i0,a)') 'Continuous in-plane rotation 2D regression completed with ', &
            failures, ' failure(s).'
        error stop 1
    endif
    write(*,'(a)') 'Continuous in-plane rotation 2D regression completed successfully.'

contains

    subroutine parse_arguments(build_dir, vol1, src_proj, logdir, threads, box_crop, &
            &diagnostic, only_tests)
        character(len=*), intent(inout) :: build_dir, vol1, src_proj, logdir
        character(len=*), intent(inout) :: threads, box_crop, diagnostic, only_tests
        character(len=4096) :: arg, key, value
        integer :: i, eq
        do i = 1, command_argument_count()
            call get_command_argument(i, arg)
            eq = index(arg, '=')
            if( eq <= 1 ) cycle
            key = adjustl(arg(:eq-1))
            value = adjustl(arg(eq+1:))
            select case(trim(key))
            case('build');      build_dir = trim(value)
            case('vol1');       vol1 = trim(value)
            case('projfile');   src_proj = trim(value)
            case('logdir');     logdir = trim(value)
            case('threads');    threads = trim(value)
            case('box_crop');   box_crop = trim(value)
            case('diagnostic'); diagnostic = trim(value)
            case('only_tests'); only_tests = trim(value)
            case default
                write(*,'(a)') 'WARNING: ignored argument '//trim(arg)
            end select
        enddo
    end subroutine parse_arguments

    logical function test_selected(only, label)
        character(len=*), intent(in) :: only, label
        character(len=8192) :: list
        list = ','//trim(only)//','
        test_selected = trim(only) == 'all' .or. index(list, ','//trim(label)//',') > 0
    end function test_selected

    subroutine make_directory(path)
        character(len=*), intent(in) :: path
        integer :: status, cmdstat
        call execute_command_line('mkdir -p '//trim(quote(path)), wait=.true., &
            &exitstat=status, cmdstat=cmdstat)
        if( cmdstat /= 0 .or. status /= 0 ) error stop 'could not create regression log directory'
    end subroutine make_directory

    subroutine run_workflow(label, inpl_mode, src_proj, threads, diagnostic, logdir, failures)
        character(len=*), intent(in) :: label, inpl_mode, src_proj, threads, diagnostic, logdir
        integer,          intent(inout) :: failures
        character(len=8192) :: command
        command = 'simple_exec prg=abinitio2D projfile='//trim(quote(src_proj))// &
            &' mkdir=yes mskdiam=120 ncls=4 nstages=6 nits_per_stage=5 nthr='//trim(threads)// &
            &' objfun=euclid inpl_cont='//trim(inpl_mode)//' refine=snhc_smpl'
        if( trim(diagnostic) == 'yes' ) command = 'SIMPLE_INPL_DIAGNOSTIC=yes '//trim(command)
        call run_test('euclid_'//trim(label)//'_workflow', command, logdir, failures)
    end subroutine run_workflow

    subroutine newest_project(build_dir, logdir, scratch_name, project, failures)
        character(len=*), intent(in) :: build_dir, logdir, scratch_name
        character(len=*), intent(out) :: project
        integer, intent(inout) :: failures
        character(len=8192) :: scratch, command
        integer :: unit, ios, status, cmdstat
        project = ''
        scratch = trim(logdir)//'/'//trim(scratch_name)
        command = 'find '//trim(quote(build_dir))//' -type f -name onejyx_abinitio_v3.simple '// &
            '-printf ''%T@ %p\n'' | sort -nr | head -n 1 | cut -d'' '' -f2- > '//trim(quote(scratch))
        call execute_command_line(command, wait=.true., exitstat=status, cmdstat=cmdstat)
        if( cmdstat /= 0 .or. status /= 0 )then
            call record_failure(failures, 'could not locate newest output project')
            return
        endif
        open(newunit=unit, file=scratch, status='old', action='read', iostat=ios)
        if( ios == 0 )then
            read(unit,'(a)',iostat=ios) project
            close(unit, status='delete')
        endif
        if( ios /= 0 .or. len_trim(project) == 0 )then
            call record_failure(failures, 'newest output project not found')
            project = ''
        endif
    end subroutine newest_project

    subroutine find_newest_file(search_root, filename, logdir, scratch_name, selected)
        character(len=*), intent(in) :: search_root, filename, logdir, scratch_name
        character(len=*), intent(out) :: selected
        character(len=8192) :: scratch, command
        integer :: unit, ios, status, cmdstat

        selected = ''
        scratch = trim(logdir)//'/'//trim(scratch_name)
        command = 'find '//trim(quote(search_root))//' -type f -name '//trim(quote(filename))// &
            ' -printf ''%T@ %p\n'' | sort -nr | head -n 1 | cut -d'' '' -f2- > '//trim(quote(scratch))
        call execute_command_line(command, wait=.true., exitstat=status, cmdstat=cmdstat)
        if( cmdstat /= 0 .or. status /= 0 ) return
        open(newunit=unit, file=scratch, status='old', action='read', iostat=ios)
        if( ios == 0 )then
            read(unit,'(a)',iostat=ios) selected
            close(unit, status='delete')
        endif
        if( ios /= 0 ) selected = ''
    end subroutine find_newest_file

    subroutine run_test(label, command, logdir, failures)
        character(len=*), intent(in) :: label, command, logdir
        integer, intent(inout) :: failures
        character(len=8192) :: logfile, full_command
        integer :: status, cmdstat, start_count, end_count, rate
        call system_clock(start_count, rate)
        logfile = trim(logdir)//'/continuous_'//trim(label)//'.log'
        write(*,'(a)') '=== TEST ['//trim(label)//'] START ==='
        write(*,'(a)') '+ '//trim(command)
        full_command = trim(command)//' > '//trim(quote(logfile))//' 2>&1'
        call execute_command_line(full_command, wait=.true., exitstat=status, cmdstat=cmdstat)
        call system_clock(end_count)
        write(*,'(a,i0)') '=== TEST ['//trim(label)//'] ELAPSED_SECONDS: ', &
            (end_count-start_count)/max(1,rate)
        if( cmdstat /= 0 .or. status /= 0 )then
            call record_failure(failures, trim(label)//' exited with status '//int_to_string(status))
            write(*,'(a)') '=== TEST ['//trim(label)//'] FAIL ==='
        else
            write(*,'(a)') '=== TEST ['//trim(label)//'] PASS ==='
        endif
    end subroutine run_test

    subroutine record_failure(failures, message)
        integer, intent(inout) :: failures
        character(len=*), intent(in) :: message
        failures = failures + 1
        write(*,'(a)') 'ERROR: '//trim(message)
    end subroutine record_failure

    function quote(value) result(out)
        character(len=*), intent(in) :: value
        character(len=8192) :: out
        out = ''''//trim(value)//''''
    end function quote

    function int_to_string(value) result(out)
        integer, intent(in) :: value
        character(len=32) :: out
        write(out,'(i0)') value
    end function int_to_string

end program simple_test_continuous_inplane_rotation2D
