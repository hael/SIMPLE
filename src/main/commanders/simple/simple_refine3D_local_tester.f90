!@descr: local-refinement startup rejection checks
module simple_refine3D_local_tester
use, intrinsic :: iso_fortran_env, only: int8
use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
use simple_string, only: string
use simple_syslib, only: simple_mkdir, simple_chdir, simple_list_dirs, file_exists
use simple_cmdline, only: cmdline
use simple_image, only: image
use simple_sp_project, only: sp_project
use simple_commanders_refine3D, only: commander_refine3D_local
use simple_test_dispatch, only: test_dispatcher, test_execution, EXEC_FAILED, proc_available
use simple_test_utils, only: assert_true, assert_false, assert_int, enter_fixture, leave_fixture, tests_failed
implicit none
private
public :: run_all_refine3D_local_tests

integer, parameter :: BOX = 16
real, parameter :: SMPD = 1.5
character(len=24), parameter :: CASES(16) = [character(len=24) :: &
    &'missing_even', 'missing_odd', 'bad_dimensions', 'bad_sampling', 'bad_extent', 'bad_projection', &
    &'bad_angle', 'bad_shift', 'fractional_half', 'unassigned_half', 'one_half', 'no_active', &
    &'merged_lp', 'merged_nu', 'reference_halves', 'reference_state']
character(len=20), parameter :: INPUTS(4) = [character(len=20) :: &
    &'source.simple', 'reference.mrc', 'reference_even.mrc', 'reference_odd.mrc']
type :: file_snapshot
    integer(int8), allocatable :: bytes(:)
    logical :: exists = .false.
end type file_snapshot

contains

    ! Fatal checks run in isolated children; the parent checks exit, reason and source preservation.
    subroutine run_all_refine3D_local_tests()
        type(test_dispatcher) :: dispatcher
        type(test_execution) :: execution
        type(file_snapshot) :: before(4), after
        type(string) :: saved, root, case_dir, name
        type(string), allocatable :: dirs(:)
        integer :: icase, i, nfail0
        logical :: interrupted
#if !defined(__linux__)
        write(*,'(A)') 'SKIP: local startup rejection tests require the Linux test dispatcher'
        return
#endif
        if( .not. proc_available() )then
            write(*,'(A)') 'SKIP: local startup rejection tests require /proc'
            return
        endif
        nfail0 = tests_failed
        call enter_fixture('local_startup', saved, root)
        do icase = 1, size(CASES)
            case_dir = root%to_char()//'/'//trim(CASES(icase))
            call simple_mkdir(case_dir)
            call simple_chdir(case_dir)
            call write_inputs(trim(CASES(icase)))
            do i = 1, size(INPUTS)
                call snapshot(string(trim(INPUTS(i))), before(i))
            enddo
            call simple_chdir(root)
            call dispatcher%new(1, 1)
            execution = test_execution()
            execution%prg = trim(CASES(icase))
            execution%label = trim(CASES(icase))
            execution%dir = case_dir%to_char()
            execution%deadline = 10
            call dispatcher%add(execution)
            call dispatcher%run(rejection_route, interrupted)
            execution = dispatcher%get_execution(1)
            call assert_false(interrupted, trim(CASES(icase))//': child was not interrupted')
            call assert_int(EXEC_FAILED, execution%state, trim(CASES(icase))//': input was rejected')
            call assert_true(execution%exit_code /= 0, trim(CASES(icase))//': rejection has nonzero exit')
            name = case_dir%to_char()//'/run.log'
            call assert_true(log_has(name, expected_reason(trim(CASES(icase)))), &
                &trim(CASES(icase))//': intended preflight diagnostic')
            dirs = simple_list_dirs(case_dir)
            if( allocated(dirs) )then
                call assert_int(0, size(dirs), trim(CASES(icase))//': no parameter-created run directory')
                deallocate(dirs)
            else
                call assert_true(.true., trim(CASES(icase))//': no parameter-created run directory')
            endif
            do i = 1, size(INPUTS)
                name = case_dir%to_char()//'/'//trim(INPUTS(i))
                call snapshot(name, after)
                call assert_true(after%exists .eqv. before(i)%exists, &
                    &trim(CASES(icase))//': input existence unchanged '//trim(INPUTS(i)))
                if( .not. before(i)%exists .or. .not. after%exists ) cycle
                call assert_int(size(before(i)%bytes), size(after%bytes), 'source byte count unchanged')
                if( size(before(i)%bytes) == size(after%bytes) )then
                    call assert_true(all(before(i)%bytes == after%bytes), &
                        &trim(CASES(icase))//': source bytes unchanged '//trim(INPUTS(i)))
                endif
            enddo
            call dispatcher%kill()
        enddo
        call leave_fixture(saved, root, nfail0)
    end subroutine run_all_refine3D_local_tests

    ! The real public commander must reject before params%new, not a test-only copy of its checks.
    subroutine rejection_route( which, cline )
        character(len=*), intent(in) :: which
        class(cmdline), intent(inout) :: cline
        class(commander_refine3D_local), allocatable :: command
        allocate(command)
        call cline%delete('test')
        call cline%set('prg', 'refine3D_local')
        call cline%set('projfile', 'source.simple')
        call cline%set('vol1', 'reference.mrc')
        call cline%set('mskdiam', 12.)
        call cline%set('nthr', 1)
        call cline%set('maxits', 1)
        call cline%set('mkdir', 'yes')
        select case(which)
            case('merged_lp')
                call cline%set('filt_mode', 'none')
                call cline%set('lp', 8.)
            case('merged_nu')
                call cline%set('filt_mode', 'nonuniform_lpset')
        end select
        call command%execute(cline)
    end subroutine rejection_route

    subroutine write_inputs( which )
        use simple_parameters, only: parameters
        use simple_refine3D_provenance, only: begin_reference_set, publish_reference_set
        use simple_fileio, only: arr2file
        character(len=*), intent(in) :: which
        class(sp_project), allocatable :: project
        class(image), allocatable :: volume
        class(parameters), allocatable :: params
        type(string) :: halves(2)
        integer :: i, dims(3), producer_state
        real :: sampling, nan
        allocate(project, volume)
        call project%os_stk%new(1, is_ptcl=.false.)
        call project%os_stk%set(1, 'box', BOX)
        call project%os_stk%set(1, 'smpd', SMPD)
        call project%os_stk%set(1, 'stk', 'particles.mrc')
        call project%os_ptcl3D%new(2, is_ptcl=.true.)
        do i = 1,2
            call project%os_ptcl3D%set_euler(i, [0., 0., 0.])
            call project%os_ptcl3D%set_shift(i, [0., 0.])
            call project%os_ptcl3D%set_state(i, 1)
            call project%os_ptcl3D%set(i, 'proj', 1.)
            call project%os_ptcl3D%set(i, 'eo', real(i-1))
            call project%os_ptcl3D%set(i, 'stkind', 1)
            call project%os_ptcl3D%set(i, 'indstk', i)
        enddo
        nan = ieee_value(0., ieee_quiet_nan)
        select case(which)
            case('bad_projection')
                call project%os_ptcl3D%set(1, 'proj', 0.)
            case('bad_angle')
                call project%os_ptcl3D%set_euler(1, [nan, 0., 0.])
            case('bad_shift')
                call project%os_ptcl3D%set_shift(1, [nan, 0.])
            case('fractional_half')
                call project%os_ptcl3D%set(1, 'eo', 0.5)
            case('unassigned_half')
                call project%os_ptcl3D%set(1, 'eo', -1.)
            case('one_half')
                call project%os_ptcl3D%set(1, 'eo', 1.)
            case('no_active')
                call project%os_ptcl3D%set_all2single('state', 0.)
        end select
        call project%write(string('source.simple'))
        call project%kill()
        do i = 2, size(INPUTS)
            if( which == 'missing_even' .and. i == 3 ) cycle
            if( which == 'missing_odd' .and. i == 4 ) cycle
            dims = BOX
            sampling = SMPD
            if( which == 'bad_dimensions' .and. i == 3 ) dims = BOX/2
            if( which == 'bad_sampling' .and. i == 4 ) sampling = 2.
            if( which == 'bad_extent' ) sampling = 2.
            call volume%new(dims, sampling)
            call volume%write(string(trim(INPUTS(i))))
            call volume%kill()
        enddo
        if( which == 'reference_halves' .or. which == 'reference_state' )then
            allocate(params)
            call project%read(string('source.simple'))
            call project%update_projinfo(string('source.simple'))
            params%oritype = 'ptcl3D'
            params%projfile = 'source.simple'
            params%box_crop = BOX
            params%smpd_crop = SMPD
            params%l_ml_reg = .false.
            params%l_nonuniform = .false.
            halves = [string('reference_even.mrc'), string('reference_odd.mrc')]
            call arr2file([1.,0.5], string('reference_fsc.bin'))
            producer_state = 1
            if( which == 'reference_state' ) producer_state = 2
            call begin_reference_set(string('reference.mrc'))
            call publish_reference_set(params, project, project%os_ptcl3D, producer_state, &
                &string('reference.mrc'), halves, string('reference_fsc.bin'), 0.)
            if( which == 'reference_halves' )then
                call project%os_ptcl3D%set(1, 'eo', 1)
                call project%os_ptcl3D%set(2, 'eo', 0)
                call project%write(string('source.simple'))
            endif
            call project%kill()
        endif
    end subroutine write_inputs

    function expected_reason( which ) result(reason)
        character(len=*), intent(in) :: which
        character(len=96) :: reason
        select case(which)
            case('missing_even')
                reason = 'Continuous startup requires an existing reference: reference_even.mrc'
            case('missing_odd')
                reason = 'Continuous startup requires an existing reference: reference_odd.mrc'
            case('bad_dimensions', 'bad_sampling', 'bad_extent')
                reason = 'Continuous startup reference dimensions, sampling, or physical extent are incompatible'
            case('merged_lp')
                reason = 'refine3D_local requires half-specific references; use lpstop instead of lp'
            case('merged_nu')
                reason = 'refine3D_local does not support filt_mode=nonuniform_lpset'
            case('reference_halves')
                reason = 'Continuous startup reference provenance failed: reference half assignments'
            case('reference_state')
                reason = 'Continuous startup reference provenance failed: reference producer state'
            case('one_half')
                reason = 'Continuous startup requires particles in both half-sets'
            case('no_active')
                reason = 'refine=cont found no active particle'
            case default
                reason = 'refine=cont continues from the 3D poses'
        end select
    end function expected_reason

    logical function log_has( name, expected ) result(found)
        class(string), intent(in) :: name
        character(len=*), intent(in) :: expected
        character(len=2048) :: line
        integer :: unit, ios
        found = .false.
        open(newunit=unit, file=name%to_char(), status='old', action='read', iostat=ios)
        if( ios /= 0 ) return
        do
            read(unit,'(A)',iostat=ios) line
            if( ios /= 0 ) exit
            if( index(line, trim(expected)) > 0 ) found = .true.
        enddo
        close(unit)
    end function log_has

    subroutine snapshot( name, saved )
        class(string), intent(in) :: name
        type(file_snapshot), intent(inout) :: saved
        integer :: unit, nbytes
        if( allocated(saved%bytes) ) deallocate(saved%bytes)
        saved%exists = file_exists(name)
        if( .not. saved%exists ) return
        inquire(file=name%to_char(), size=nbytes)
        allocate(saved%bytes(nbytes))
        open(newunit=unit, file=name%to_char(), access='stream', form='unformatted', status='old', action='read')
        read(unit) saved%bytes
        close(unit)
    end subroutine snapshot

end module simple_refine3D_local_tester
