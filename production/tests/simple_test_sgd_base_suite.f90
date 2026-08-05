program simple_test_sgd_base_suite
    use base_sgd_direct_test, only: run_direct_shift
    use base_sgd_sampling_restore_test, only: run_sampling_and_restore_policy
    use base_sgd_v2_test, only: run_abinitio_v2
    use base_sgd_v3_test, only: run_abinitio_v3
    use base_sgd_v4_test, only: run_abinitio_v4
    use base_sgd_baseline_test, only: run_abinitio_baseline
    implicit none

    character(len=256) :: arg
    integer :: arg_len, arg_status

    call get_command_argument(1, arg, arg_len, arg_status)
    if (arg_status == 0 .and. arg_len > 5) then
        if (arg(:5) == 'case=') then
            call run_case(arg(6:arg_len))
        else
            call run_suite()
        end if
    else
        call run_suite()
    end if

contains

    subroutine run_suite()
        character(len=*), parameter :: labels(*) = [ character(len=20) :: &
            'direct_shift', 'sampling_restore', 'abinitio_v2', 'abinitio_v3', 'abinitio_v4', 'baseline' ]
        integer :: failures, i

        failures = 0
        write (*, '(a)') 'SIMPLE base SGD regression suite: START'
        do i = 1, size(labels)
            call run_case_in_child(trim(labels(i)), failures)
        end do

        write (*, '(/,a)') 'SIMPLE base SGD regression suite: SUMMARY'
        write (*, '(a,i0)') '  Tests run:    ', size(labels)
        write (*, '(a,i0)') '  Tests passed: ', size(labels) - failures
        write (*, '(a,i0)') '  Tests failed: ', failures
        if (failures /= 0) then
            write (*, '(a)') 'SIMPLE base SGD regression suite: FAIL'
            error stop 1
        end if
        write (*, '(a)') 'SIMPLE base SGD regression suite: PASS'
    end subroutine run_suite

    subroutine run_case_in_child(label, failures)
        character(len=*), intent(in) :: label
        integer, intent(inout) :: failures
        character(len=1024) :: executable, command
        integer :: executable_len, command_status, exit_status

        call get_command_argument(0, executable, executable_len, command_status)
        if (command_status /= 0 .or. executable_len < 1) error stop 'could not determine suite executable path'
        command = '"'//trim(executable(:executable_len))//'" case='//trim(label)
        write (*, '(/,a)') '>>> TEST ['//trim(label)//'] START'
        call execute_command_line(trim(command), wait=.true., exitstat=exit_status, cmdstat=command_status)
        if (command_status == 0 .and. exit_status == 0) then
            write (*, '(a)') '>>> TEST ['//trim(label)//'] PASS'
        else
            failures = failures + 1
            write (*, '(a,i0,a,i0)') '>>> TEST ['//trim(label)//'] FAIL (cmdstat=', command_status, &
                ', exitstat=', exit_status
        end if
    end subroutine run_case_in_child

    subroutine run_case(label)
        character(len=*), intent(in) :: label

        write (*, '(a)') '>>> SIMPLE base SGD regression child case: '//trim(label)
        select case (trim(label))
        case ('direct_shift')
            call run_direct_shift()
        case ('sampling_restore')
            call run_sampling_and_restore_policy()
        case ('abinitio_v2')
            call run_abinitio_v2()
        case ('abinitio_v3')
            call run_abinitio_v3()
        case ('abinitio_v4')
            call run_abinitio_v4()
        case ('baseline')
            call run_abinitio_baseline()
        case default
            error stop 'unknown SIMPLE base SGD regression case: '//trim(label)
        end select
    end subroutine run_case

end program simple_test_sgd_base_suite
